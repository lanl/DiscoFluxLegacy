"""
some helper functions for the postprocessing script
© 2026. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform
publicly and display publicly, and to permit others to do so.
"""
import lzma
import random
import numpy as np
if np.__version__ > "1.27":
    np.set_printoptions(legacy='1.25') ## numpy >= 2.0 prints its floats as np.float64(), this option avoids that

Nchar = 2 # number of character angles (only needed to convert .inc to .dat input files), needs to match Nchar in Modules.F90
## if Nchar>1, some variables will be converted to arrays with repeated entries
## (because the old .inc files did not know about char-dependence)
    
def readfortraninput(fname):
    """reads an input file for the fortran code and returns a dictionary containing the data"""
    inputparams = {}
    with open(fname,"r",encoding='utf8') as infile:
        lines = infile.readlines()
        for line in lines:
            line = line.split('!')[0].split('#')[0] ## remove comments
            currentline = line.lstrip().rstrip().split()
            if len(currentline)>2:
                for i in range(len(currentline)-1):
                    currentline[i+1] = currentline[i+1].replace('d','e')
                try:
                    value = list(np.asarray(currentline[1:],dtype='float'))
                except ValueError:
                    value = currentline[1:]
            elif len(currentline)==2:
                value = currentline[1]
                try:
                    value = float(value.replace('d','e'))
                except ValueError:
                    pass
            if len(currentline)>1:
                if currentline[0]=='euler_angle' and currentline[0] in inputparams:
                    if np.asarray(inputparams[currentline[0]]).shape==(3,):
                        inputparams[currentline[0]] = [inputparams[currentline[0]]]
                    inputparams[currentline[0]].append(value)
                else:
                    inputparams[currentline[0]] = value
    return inputparams
    
def writefortraninput(inputparams,fname,Nchar=Nchar):
    """
    writes dictionary 'inputparams' to an inputfile for the fortran code.
    Nchar denotes the number of character angles to consider for the dislocations.
    """
    if Nchar>1:
        for key in ['B0', 'wave_vel']:
            if key in inputparams and isinstance(inputparams[key],float):
                tmp = inputparams[key]
                inputparams[key] = [tmp]
                for ch in range(Nchar-1):
                    inputparams[key].append(tmp)
    with open(fname,"w", encoding='utf8') as outfile:
        for k in sorted(inputparams.keys()):
            if isinstance(inputparams[k],float):
                outfile.write(f"{k} {(str(inputparams[k])).replace('e','d')}\n")
            elif isinstance(inputparams[k],bool):
                outfile.write(f"{k} {('.'+str(inputparams[k]))+'.'}\n")
            elif isinstance(inputparams[k],str):
                outfile.write(f"{k} {inputparams[k]!s}\n")
            elif isinstance(inputparams[k],list):
                outfile.write(f"{k} {(str(inputparams[k]).replace('e','d').replace('[','').replace(']','').replace(',',''))}\n")
            else:
                print(f"skipping {k} {inputparams[k]}")

### some additional functions to read the output data of the fortran program:
def F90float(arg):
    """some quirks for reading floats from a discoflux legacy output file"""
    if len(arg)>10:
        char = arg[-4]
    else:
        char = ''
    if char == '-':
        out = arg[:-4] + 'E' + arg[-4:]
    else:
        out = arg
    return float(out)
    
def read_field_output(filename):
    """helper function for reading a discoflux legacy outout file"""
    data = {}
    increment = []
    time = {}
    prev_inc = None
    try:
        with open(filename,"r", encoding='utf8') as file1:
            lines = file1.readlines()
    except FileNotFoundError:
        # print("reading compressed file: ",filename+".xz")
        with lzma.open(filename+".xz","rt") as file1:
            lines = file1.readlines()
    for line in lines:
        if "#" != line[0]:
            col = line.strip().split()
            if len(col)>0:
                inc = int(col[0])
                if inc != prev_inc:
                    prev_inc = inc
                    data[inc] = []
                    time[inc] = float(col[1])
                    increment.append(inc)
                data[inc].append([F90float(v) for v in col[2:]])
    
    for k, v in data.items():
        data[k] = np.array(v)
        
    return increment, time, data

def read_altfield_output(filename):
    """helper function for reading a discoflux legacy outout file"""
    data = []
    time = []
    try:
        with open(filename,"r", encoding='utf8') as file1:
            lines = file1.readlines()
    except FileNotFoundError:
        # print("reading compressed file: ",filename+".xz")
        with lzma.open(filename+".xz","rt") as file1:
            lines = file1.readlines()
    for line in lines:
        if "#" != line[0]:
            col = line.strip().split()
            if len(col)>0:
                time.append(float(col[0]))
                data.append([F90float(v) for v in col[1:]])
  
      
    return time, np.array(data)

def random_euler(N,nrange=np.pi/2):
    '''generates an array of shape Nx3 filled with random numbers between 0 and nrange (np.pi/2 by default, don't need more for cubic symmetry)'''
    out = np.zeros((N,3))
    for i in range(N):
        for j in range(3):
            out[i,j] = random.random()
        print(f"euler_angle {nrange*out[i,0]:.8f} {nrange*out[i,1]:.8f} {nrange*out[i,2]:.8f}")
    return nrange*out

        
        
