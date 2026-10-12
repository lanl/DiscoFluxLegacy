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
import pandas as pd

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
    
def writefortraninput(inputparams,fname,Nchar=2):
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

def read_field_output(filename,multiindex=1,return_pd=True,Nslip=12):
    """
    helper function for reading a discoflux legacy outout file
    if option 'multiindex'=n>1, the first n columns will be used for a multi-index in the padas output.
    (multiindex>2 implies return_pd=True because the legacy output only supports multiindex=1 or 2.)
    Nslip is only needed to amend missing column names read from the output header (not necessary for new calculations)
    """
    indcol = 0
    if multiindex>1:
        indcol = [i for i in range(int(multiindex))]
    try:
        raw = pd.read_csv(filename,skip_blank_lines=True,sep=r'\s+',comment='#',header=None,index_col=indcol,low_memory=False)
        with open(filename,"r",encoding='utf8') as input_file:
            head = [next(input_file) for _ in range(2)]
    except FileNotFoundError:
        raw = pd.read_csv(filename+".xz",skip_blank_lines=True,sep=r'\s+',comment='#',header=None,index_col=indcol,low_memory=False)
        with lzma.open(filename+".xz","rt") as input_file:
            head = [next(input_file) for _ in range(2)]
    # now cleanup columns of dtype string because fortran sometime omits the "E" when the exponent is >99 or <-99 and python's float() trips over this:
    str_cols = raw.select_dtypes(include=['string','object']).columns
    # if len(str_cols)>0:
    #     print(f"{str_cols=}, cleaning data for ",filename)
    raw[str_cols] = raw[str_cols].map(F90float)
    col_names = head[1].split()[1:] # read column names from file
    # append missing column names:
    if col_names[-1] == "..." and "-3" in col_names[-2]:
        nm = col_names[-2][:-2]
        col_names[-1] = f"{nm}-4"
        for i in range(Nslip-4):
            col_names.append(f"{nm}-{i+5}")
    ln_col_nm = len(col_names)
    ln_cols = len(raw.columns)
    if return_pd or multiindex>2:
        if ln_col_nm>ln_cols+multiindex:
            col_names.remove("id")
            ln_col_nm -= 1
        if multiindex>1:
            raw.index.names = col_names[:multiindex]
        else:
            raw.index.name = col_names[0]
        if ln_cols+multiindex==ln_col_nm:
            raw.columns = col_names[multiindex:]
        return raw
    # convert to legacy format (TODO: change later code to use pandas data instead):
    elif multiindex==1:
        time = raw.index.to_list()
        data = raw.to_numpy(dtype=float)
        return time, data
    elif multiindex==2:
        increment = raw.index.get_level_values(0).unique().to_list()
        time_values = raw.index.get_level_values(1).unique().to_list()
        time = {inc:time_values[i] for i,inc in enumerate((increment))}
        data = {inc:raw.xs(inc).to_numpy(dtype=float) for inc in increment}
        return increment, time, data

def random_euler(N,nrange=np.pi/2):
    '''generates an array of shape Nx3 filled with random numbers between 0 and nrange (np.pi/2 by default, don't need more for cubic symmetry)'''
    out = np.zeros((N,3))
    for i in range(N):
        for j in range(3):
            out[i,j] = random.random()
        print(f"euler_angle {nrange*out[i,0]:.8f} {nrange*out[i,1]:.8f} {nrange*out[i,2]:.8f}")
    return nrange*out
