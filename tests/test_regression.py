#!/usr/bin/env python3
"""
Regression testing for DiscoFlux Legacy
© 2026. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform
publicly and display publicly, and to permit others to do so.
"""
import os
import sys
import pathlib
import subprocess
import shutil
import math
# import pytest
# import numpy as np
root_path = pathlib.Path(__file__).resolve().parents[1]
dir_path = str(root_path / "postprocessing")
if dir_path not in sys.path:
    sys.path.append(dir_path)
test_path = pathlib.Path(__file__).resolve().parents[0]
from helperfunctions import read_field_output
cwd =pathlib.Path.cwd()
fpm = shutil.which('fpm')
compiler = os.getenv("FC")
if compiler is None:
    basecommand = [fpm,"run","--"]
else:
    basecommand = [fpm,"run","--compiler",compiler,"--"]

def runscript(scriptname,args,logfname):
    '''Run script "scriptname" as a subprocess passing a list of command line arguments "args" and saving its stdout to a file "logfname"'''
    out = -1
    with open(logfname, 'w', encoding="utf8") as logfile:
        command = [pathlib.Path(dir_path,scriptname)]
        if sys.platform=='win32':
            command = ["python",pathlib.Path(dir_path,scriptname)]
        with subprocess.Popen(command+args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as subproc:
            for line in subproc.stdout:
                sys.stdout.write(line)
                logfile.write(line)
            subproc.wait()
            out = subproc.returncode
    return out

print(f"\nrunning regression tests using {compiler=}")

def test_impactCu(newbaseline=False):
    """regression test for 1D impact simulation on copper"""
    command = basecommand.copy()
    fname = "Cu1.dat"
    print(f"\nrunning a regression test for {fname}")
    baseline = {'T':246855.3369005, "Ei":0.0749191082212, "Fe":2460.001066108, "Fp":2460.002187059, 
                "tau":80914.152052, "tau_back":8643.48995368, "stress":562124.0340608, 
                "rho_pos.1":2659742105.104, "rho_neg.1":2659742105.104, "dis_vel_x.1":391654883.4191, 
                "dis_vel_y.1":195088655.52, "dis_vel.1":552489174.0105, "const_soln":800, 
                "dis_acc.1":1.623897325981e17, "dmb_vel":4208128.086478, "dmb_acc":395569933011000.0, 
                "dmb_pos":165.6705104228}
    command.append(test_path / fname)
    with open("test_impactCu.log", 'w', encoding="utf8") as logfile, \
      subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as subproc:
        for line in subproc.stdout:
            sys.stdout.write(line)
            logfile.write(line)
        subproc.wait()

    for k,v in baseline.items():
        if k in ['dmb_vel', 'dmb_acc', 'dmb_pos']:
            outfile = test_path /f"{fname}.node.{k}.F90txt"
        else:
            outfile = test_path /f"{fname}.ip.{k}.F90txt"
        f1 = read_field_output(outfile)
        f1_thinned = f1.iloc[::50].iloc[:,2:]
        x = f1_thinned.abs().sum().sum()
        # print(k,v,abs(f1_thinned.abs().sum().sum()/v - 1))
        if newbaseline:
            baseline[k] = v= round(x, 12 - math.floor(math.log10(abs(x))))

        acc = 1e11
        if k in ('Ei', 'tau', 'tau_back', 'dis_acc'):
            acc = 1e-4
        elif k in ('stress', 'rho_pos.1', 'rho_neg.1', 'dis_vel.1', 'dis_vel_x.1', 
                   'dis_vel_y.1', 'dmb_vel', 'dmb_acc'):
            acc = 1e-6
        elif k in ('T', 'Fe', 'Fp'):
            acc = 1e-9
        assert abs(x/v - 1) < acc, \
            f"{k} output changed: got {f1_thinned.abs().sum().sum()}, expected {v}"
    if newbaseline:
        return baseline

def test_plots():
    """make sure the plotting script runs without failing"""
    commandargs = [test_path / "Cu1.dat"]
    print("\nrunning the simple_postprocess.py script")
    assert runscript(root_path / "postprocessing" / "simple_postprocess.py",commandargs,'plots.log')==0
