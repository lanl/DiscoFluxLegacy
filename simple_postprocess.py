#!/usr/bin/env python3
### python 3.8 or higher required
### also: plot_model_stress=True requires PyDislocDyn >=1.2.9
import sys
import os
import numpy as np
from scipy.signal import argrelextrema
from scipy.optimize import fmin

import matplotlib as mpl
mpl.use('Agg', force=False) # don't need X-window, allow running in a remote terminal session
import matplotlib.pyplot as pyl
from matplotlib.ticker import AutoMinorLocator
## workaround for spyder's runfile() command when cwd is somewhere else:
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path)
##
from mkinput import readinputdata, read_field_output, read_altfield_output, random_euler

vel_direction = 1 ## decide which component of dmb velocities to plot (typically 1 for impact, 2 for shear; 'all' plots for 1,2, and 3)
slipsystem = 1 ## decide for which slip system to plot disloc. velocities (1-12 for fcc)
fs_direction = 1 ## decide which component of the fs velocity to plot
stress_comp = 'all' ## which component of the stress tensor to plot, 'all' plots all of [1='11',2='22',3='33',4='23',5='31,6='12']; ('uniaxial stress' in old version referred to '11')
fs_position = 'auto' ## node index at which to plot fs_stress over time, i.e. -1 is the free surface but with -n<-1 one might simulate a gauge; 'auto'= ipOut if given in input file and -1 otherwise
max_snapshot_time = 1/2 ## relative to total simulation time, i.e. must be <=1

## choose which plots to generate:
plot_rho_pos = True
plot_rho_neg = True
plot_dis_vel = True ## disloc. velocity |v|
plot_dis_vel_x = False ## global x component v_x
plot_dis_vel_y = False ## v_y
plot_dis_vel_z = False ## v_z
plot_av_dis_vel_x = False
plot_tau = True
plot_tau_back = True
plot_dmb_vel = True
plot_stress = True
plot_T = True
plot_fs_vel = True
plot_dev_stress = True
plot_model_stress = False ## additional plot showing comparison with a simpler model for precursor decay, requires PyDislocDyn
plot_fs_stress = True
plot_pressure = True
plot_dis_acc = True ## acceleration of a dislocation as a function of x at various time snapshots (compute time-derivatives of dis_vel data)

crystalstruct='fcc'
Nslip = 12 # default number of slip systems (used to normalize disloc. densities in plots)

if len(sys.argv) == 2:
    jobname = sys.argv[1]
else:
    print(f"USAGE: {os.path.basename(__file__)} 'jobname'")
    sys.exit()
    
## determine Nchar from number of .ip.dis_vel.{}.F90txt files for jobname (hence no longer necessary to edit default value of mkinput.py)
Nchar = len([X for X in os.listdir('.') if jobname + ".ip.dis_vel_x" in X and "F90txt" in X])
fname = 'input_parameters.'+jobname+'.inc'
inputdata = readinputdata(fname)
rho0 = inputdata['rho0']
L0 = inputdata['L0']
if 'Nchar' in inputdata.keys(): 
    Nchar = int(inputdata['Nchar'])
if Nchar==0:
    raise ValueError("cannot find required files for job ",jobname)
if 'crystalstruct' in inputdata.keys():
    crystalstruct = inputdata['crystalstruct']
    if crystalstruct=='fcc':
        Nslip=12
    elif crystalstruct=='bcc':
        Nslip=48
if 'Nslip' in inputdata.keys(): 
    Nslip = int(inputdata['Nslip'])
if 'ipOut' in inputdata.keys():
    fs_position = int(inputdata['ipOut'])
    ipOut = int(inputdata['ipOut'])
elif fs_position == 'auto':
    fs_position = -1
    ipOut = -1
elif fs_position == 'oldauto': # use this for output created by program version <= 2021.01.11
    fs_position = -1
    ipOut = 0
else:
    ipOut = 0
print(f"{crystalstruct=}, {Nslip=}, {Nchar=}")
    
try:
    from pydislocdyn import metal_props
    mpl.use('Agg') # ensure we are not using the LaTeX backend (even if pydislocdyn selected it)
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "text.usetex": False,
        "pgf.rcfonts": True,
    })
    have_poly = True
except ImportError:
    if plot_model_stress:
        raise ImportError("ERROR: Cannot find PyDislocDyn, need this packge to generate requested plot 'plot_model_stress': \
              either set plot_model_stress=False or download 'https://github.com/dblaschke-LANL/PyDislocDyn' and copy its contents to the same folder as this script.\n")
    have_poly = False

rho_pos = {}
rho_neg = {}
dis_vel = {}
dis_acc = {}
dis_vel_x = {}
dis_vel_y = {}
dis_vel_z = {}
characters = list(range(1,Nchar+1))
Nnode=0
inc_list = []
len_time = np.inf
Nveldir=1
Nslipsystems=1
Nstresscomp=1
if vel_direction == 'all':
    Nveldir=3
if slipsystem == 'all':
    Nslipsystems=Nslip
if stress_comp == 'all':
    Nstresscomp=6
number_of_plots = Nchar*Nslipsystems*(plot_rho_pos + plot_rho_neg + plot_dis_vel + plot_dis_vel_x + plot_dis_vel_y + plot_dis_vel_z + plot_dis_acc) \
                + Nchar*plot_av_dis_vel_x + Nslipsystems*(plot_tau + plot_tau_back) + Nveldir*plot_dmb_vel + Nstresscomp*plot_stress \
                + plot_T + plot_fs_vel + plot_dev_stress + plot_fs_stress + plot_pressure + plot_model_stress
print(f"reading data for {number_of_plots} plots")
for ich in characters:
    if plot_rho_pos or plot_av_dis_vel_x:
        inc_list, time_values, rho_pos[ich]   = read_field_output(jobname + ".ip.rho_pos.{}.F90txt".format(ich))
        Nnode = len(rho_pos[ich][inc_list[0]])+1
        len_time = min(len(inc_list),len_time)
    if plot_rho_neg or plot_av_dis_vel_x:
        inc_list, time_values, rho_neg[ich]   = read_field_output(jobname + ".ip.rho_neg.{}.F90txt".format(ich))
        Nnode = len(rho_neg[ich][inc_list[0]])+1
        len_time = min(len(inc_list),len_time)
    if plot_dis_vel:
        inc_list, time_values, dis_vel[ich]   = read_field_output(jobname + ".ip.dis_vel.{}.F90txt".format(ich))
        Nnode = len(dis_vel[ich][inc_list[0]])+1
        len_time = min(len(inc_list),len_time)
    if plot_dis_acc:
        inc_list, time_values, dis_acc[ich]   = read_field_output(jobname + ".ip.dis_acc.{}.F90txt".format(ich))
        Nnode = len(dis_acc[ich][inc_list[0]])+1
        len_time = min(len(inc_list),len_time)
    if plot_dis_vel_x or plot_av_dis_vel_x:
        inc_list, time_values, dis_vel_x[ich]   = read_field_output(jobname + ".ip.dis_vel_x.{}.F90txt".format(ich))
        Nnode = len(dis_vel_x[ich][inc_list[0]])+1
        len_time = min(len(inc_list),len_time)
    if plot_dis_vel_y:
        inc_list, time_values, dis_vel_y[ich]   = read_field_output(jobname + ".ip.dis_vel_y.{}.F90txt".format(ich))
        Nnode = len(dis_vel_y[ich][inc_list[0]])+1
        len_time = min(len(inc_list),len_time)
    if plot_dis_vel_z: ## TODO: determine from sqrt(vel^2 - vel_x^2 - vel_y^2), dont need to write this file
        inc_list, time_values, dis_vel_z[ich]   = read_field_output(jobname + ".ip.dis_vel_z.{}.F90txt".format(ich))
        Nnode = len(dis_vel_z[ich][inc_list[0]])+1
        len_time = min(len(inc_list),len_time)
###
if plot_tau:
    inc_list, time_values, tau  = read_field_output(jobname + ".ip.tau.F90txt")
    Nnode = len(tau[inc_list[0]])+1
    len_time = min(len(inc_list),len_time)
if plot_tau_back:
    inc_list, time_values, tau_back  = read_field_output(jobname + ".ip.tau_back.F90txt")
    Nnode = len(tau_back[inc_list[0]])+1
    len_time = min(len(inc_list),len_time)
if plot_dmb_vel or plot_fs_vel:
    inc_list, time_values, dmb_vel   = read_field_output(jobname + ".node.dmb_vel.F90txt")
    Nnode = len(dmb_vel[inc_list[0]])
    len_time = min(len(inc_list),len_time)
if plot_stress or plot_dev_stress or plot_fs_stress or plot_pressure or plot_model_stress:
    inc_list, time_values, stress    = read_field_output(jobname + ".ip.stress.F90txt")
    Nnode = len(stress[inc_list[0]])+1
    len_time = min(len(inc_list),len_time)
if plot_T:
    inc_list, time_values, tmp  = read_field_output(jobname + ".ip.T.F90txt")
    Nnode = len(tmp[inc_list[0]])+1
    len_time = min(len(inc_list),len_time)
if plot_fs_vel:
    fs_time_values, fs_vel = read_altfield_output(jobname + ".th.vel.F90txt")
if plot_fs_stress:
    fs_time_values, fs_stress = read_altfield_output(jobname + ".th.stress.F90txt")

## if fortran was still running, the files above may have different len(inc_list), take the shortest in this case:
inc_list = inc_list[:len_time]

Xref = np.linspace(0., L0, Nnode)
Xip  = 0.5*(Xref[:-1]+Xref[1:])

def generate_dof_history(jobname, node_index = -1, dof='dmb_vel'):
  # inc_list, time_values, dmb_vel   = read_field_output(jobname + ".node.dmb_vel.F90txt") ## no need to read this twice
  if dof=='dmb_vel':
      if ipOut == 0:
          time_history = np.empty((len(inc_list), 4 ))
          data = dmb_vel
      else:
          time_history = np.empty((len(fs_time_values), 4 ))
          data = fs_vel
  elif dof=='stress':
      if ipOut == 0:
          time_history = np.empty((len(inc_list), 7 ))
          data=stress
      else:
          time_history = np.empty((len(fs_time_values), 7 ))
          data=fs_stress
  elif dof=='devstress':
      time_history = np.empty((len(inc_list), 7 ))
      data = stress.copy()
      for ind in range(len(inc_list)):
        cauchy     = stress[inc_list[ind]][:,1:]
        pressure   = -np.sum(cauchy[:,:3], axis=1) / 3.
        # dev        = cauchy
        # dev[:,:3] += np.tile(pressure, (3,1)).T
        data[inc_list[ind]][:,1:4] += np.tile(pressure, (3,1)).T
  if ipOut == 0:
      for ii, inc in enumerate(inc_list):
        time_history[ii,:] = np.hstack((time_values[inc], data[inc][node_index,1:] ))
  else:
      for ii in range(len(fs_time_values)):
        time_history[ii,:] = np.hstack((fs_time_values[ii], data[ii,1:] ))

  return time_history
  

def compute_ssd_gnd_profile(jobname,ich):
### for a given character index ich (where ich=Nchar=edge always, and if Nchar>1, ich=1 encodes a screw)
  # inc_list, time_values, rho_pos[ich]   = read_field_output(jobname + ".ip.rho_pos.{}.F90txt".format(ich))
  # inc_list, time_values, rho_neg[ich]   = read_field_output(jobname + ".ip.rho_neg.{}.F90txt".format(ich))

  kappa     = {}
  rho_total = {}
  for inc in inc_list:
    kappa[inc]      = rho_pos[ich][inc] - rho_neg[ich][inc]
    kappa[inc][:,0] = rho_pos[ich][inc][:,0]

    rho_total[inc]      = rho_pos[ich][inc] + rho_neg[ich][inc]
    rho_total[inc][:,0] = rho_pos[ich][inc][:,0]
    
  return rho_total, kappa
  

#--------------------------------------------------------------------------------------------------------------------#
#                                           Plotting Parameters                                                      #
#--------------------------------------------------------------------------------------------------------------------#
#        suggested font size 10 for paper, 14 for presentation
#        set your figure size to the final desired size, especially for presentations so that you don't resize in powerpoint
fw          = 3.3    # width of figure in inches
fh          = 3.0     # height of figure in inches
txt_size    = 10      # size for text in axis labels in points
nmrks       = 0      # number of 'markers' to place along data line
markers     = ['d','^','o','*','s','v','+','h']
linecolors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd', '#8c564b'] ## mimic matplotlibs default color cycle (of mpl >=2)
# timestamps  = [15,30,45,60,75] ## old hard coded indices for impact problem
maxtimestamp = max_snapshot_time*len(inc_list) ## in impact problem this was 0.375*len(inc_list) assuming the input file ran up to 0.4 microsec, and index 75 corresponded to 0.15 microsec.
timestamps = list(np.linspace(int(maxtimestamp/5),maxtimestamp,5).astype(int)) ## perhaps better: automatic distribution of time stamp indices
timestamplist = np.copy(timestamps)
## if calculation is ongoing, some or all timestamps may not be available yet
for i in range(len(timestamps)):
    if len(inc_list)-2<timestamplist[-i]:
        timestamps.pop(-1)
if len(timestamps)<len(timestamplist): ## if none are left, plot only the last time step computed so far
    timestamps.append(-2)
##########################
if have_poly:
    metal = metal_props('fcc')
    metal.rho, metal.c11, metal.c12, metal.c44 = (1e12*inputdata['rhobar0'],1e6*inputdata['C11'],1e6*inputdata['C12'],1e6*inputdata['C44'])
    metal.init_all()
###########

def modelprec(pos,modelstress,rhom,G_over_cl='auto_edge'):
    '''required parameters: position array 'pos', initial stress at position 0 'modelstress', mobile dislocation density 'modelrhom',
       and shear-modulus-over-longitudinal-soundspeed 'G_over_cl'. The latter can either be an explicit float or one of two keywords:
       G_over_cl='auto_lame' will compute an average shear modulus and an  average long. sound speed from c11, c12, c44, whereas
       G_over_cl='auto_edge' will compute the effective shear modulus and longitudinal sound speed in the direction of edge dislocation glide motion in an fcc crystal.
       For both computations, we require the Python package 'PyDislocDyn'', which is available on Github.'''
    burgers = inputdata['burger']/1e3 ## m ~3e-10 for copper
    B0 = 1e6*min(inputdata['B0']) ## Pas, 24.0e-6 for edge dislocations in copper
    vcrit = min(inputdata['wave_vel'])/1e3 ## m/s, 1620 for edge dislocations in copper
    rho = 1e12*inputdata['rhobar0'] ## material density in kg/m^3
    shear = rho*vcrit**2 ## effective shear modulus corresponding to limiting velocity
    c44 = 1e6*inputdata['C44']
    if G_over_cl=='auto_lame':
        G_over_cl = metal.mu/metal.cl
        # print("autolame:",metal.mu/1e9,metal.cl,G_over_cl)
    elif G_over_cl=='auto_edge':
        cl110 = max(metal.computesound([1,1,0]))
        G_over_cl = shear / cl110
        # print("autoedge:",metal.c44/1e9,cl110,G_over_cl)
    ### don't really know what G/cl should be in the anisotropic case, 
    ### but looking at the assumptions Taylor made in his 1965 paper, the relevant quantities in the anisotropic case
    ### should be the critical velocity (c44 in the case of fcc edge) and the longitudinal sound speed in the direction
    ### of dislocation motion (here cl110 for fcc edge); this means we have different values for fcc screw e.g., but since edge disloc.
    ### are slower (smaller vcrit), we assume here that their values dominate, concluding that 'auto_edge' should work best
    ## could it be G computed from faster shear wave speed and longitudinal wave speed in 110 direction?
    
    def dragB(sigma,B0=B0,vcrit=vcrit,b=burgers):
        return B0*np.sqrt(1+(b*sigma/(vcrit*B0))**2)
    
    def vel(sigma,B0=B0,vcrit=vcrit,b=burgers):
        return b*sigma/dragB(sigma,B0,vcrit,b)
    
    def eps(sigma,B0=B0,vcrit=vcrit,b=burgers,rhom=rhom):
        return rhom*b*vel(sigma,B0,vcrit,b)
    
    prec = np.empty(pos.shape)
    for i,x in enumerate(pos):
        if i==0:
            prec[0] = modelstress
        else:
            epsprec = eps(prec[i-1])
            if isinstance(epsprec, np.ndarray):
                epsprec = epsprec[0]
            prec[i] = prec[i-1]-(pos[i]-pos[i-1])*(2/3)*(2*G_over_cl)*epsprec
    return prec


def plot_snapshots(data,pl_typ,vel_dir=vel_direction,char=0,annotate=True,setlimits=True,adjustticks=False,slipsystem=slipsystem,stresscomp=stress_comp,weights=None,minmax=False):
    '''pl_typ= one of the following keywords: 'T' for temperature, 'dmb_vel' for material velocity, 'stress',
       'rho_pos'/'rho_neg' for positive/negative dislocation density, 'dis_vel', for dislocation velocity, 'tau' for resolved shear stress,
       and 'tau_back' for back stress.
       Types 'dmb_vel', and 'dis_vel' take the additional argument 'vel_dir' to determine which component of the 3-vector to plot.
       Types 'rho_pos', 'rho_neg', and 'dis_vel' take the additional argument 'char' (for dislocation character index) for their output filenames,
       assuming 'data' was passed just for that character angle.'''
    legendops={'loc':'best'}
    figname = ".fig."+pl_typ
    if pl_typ=='T':
        ylabel=r'$\mathregular{T: \/ K }$'
    elif pl_typ=='dmb_vel':
        figname+="_x{}".format(vel_dir)
        ylabel=r'$\mathregular{V_{m}: \/ m/s }$'
        legendops={'loc':'upper left','bbox_to_anchor':(1.01,1)}
        locmin = np.zeros((2,len(timestamps)))
        locmax = np.zeros((2,len(timestamps)))
    elif pl_typ=='stress':
        figname+="_{}".format(stresscomp)
        ylabel=r'$\mathregular{Stress:  \/GPa }$'
        legendops={'loc':'upper left','bbox_to_anchor':(1.01,1)}
    elif pl_typ=='pressure':
        ylabel=r'$\mathregular{pressure:  \/GPa }$'
    elif pl_typ=='rho_pos':
        figname+="_s{0}.{1}".format(slipsystem,char)
        ylabel=r'$(\varrho_+-\varrho_0)/1e3$'
    elif pl_typ=='rho_neg':
        figname+="_s{0}.{1}".format(slipsystem,char)
        ylabel=r'$(\varrho_--\varrho_0)/1e3$'
    elif pl_typ=='dis_vel' or pl_typ=='dis_vel_x' or pl_typ=='dis_vel_y' or pl_typ=='dis_vel_z':
        figname+="_s{0}.{1}".format(slipsystem,char)
        ylabel=r'$\mathregular{v}_{dis}: \/ \mathregular{m/s}$'
    elif pl_typ=='av_dis_vel_x':
        figname+="_av.{}".format(char)
        ylabel=r'$\mathregular{v}_{dis}: \/ \mathregular{m/s}$'
    elif pl_typ=='dis_acc':
        figname+="_s{0}.{1}".format(slipsystem,char)
        ylabel=r'$\dot{\mathregular{v}}_{dis}: \/ \mathregular{m/s^2}$'
    elif pl_typ=='tau':
        figname+="_s{}".format(slipsystem)
        ylabel=r'$\tau: \/ \mathregular{MPa}$'
    elif pl_typ=='tau_back':
        figname+="_s{}".format(slipsystem)
        ylabel=r'$\tau_b: \/ \mathregular{MPa}$'
        
    fig = pyl.figure(figsize=(fw, fh))
    ax  = pyl.subplot(111)
    filename = "./" + jobname + figname
    # print(pl_typ,data[inc_list[0]].shape)
    lendata = data[inc_list[0]].shape[0]
    if nmrks==1:
        markevery = slice( 0, lendata, lendata )
    else:
        markevery = slice( 0, lendata, int(lendata / (nmrks-1))-1 )
    
    Xdata=Xip
    for ii, ind in enumerate(timestamps):
        if pl_typ=='T':
            currentdata=data[inc_list[ind]][:,1]
        elif pl_typ=='dmb_vel':
            currentdata=data[inc_list[ind]][:,vel_dir]/1.e3
            Xdata=Xref
        elif pl_typ=='stress':
            currentdata=data[inc_list[ind]][:,stresscomp]/1.e3
        elif pl_typ=='pressure':
            cauchy     = data[inc_list[ind]][:,1:]
            currentdata=-np.sum(cauchy[:,:3], axis=1) / 3./1.e3 ## pressure
        elif pl_typ=='rho_pos' or pl_typ=='rho_neg':
            currentdata=(data[inc_list[ind]][:,slipsystem]-rho0/(2*Nchar*Nslip))/1e3
        elif pl_typ=='dis_vel' or pl_typ=='dis_vel_x' or pl_typ=='dis_vel_y' or pl_typ=='dis_vel_z':
            currentdata=(data[inc_list[ind]][:,slipsystem])/1e3
        elif pl_typ=='av_dis_vel_x':
            if weights is None:
                currentdata=np.average(data[inc_list[ind]][:,1:],axis=-1)/1e3
            else:
                currentdata=np.average(data[inc_list[ind]][:,1:],axis=-1,weights=weights[inc_list[ind]][:,1:])/1e3
        elif pl_typ=='dis_acc':
            # deltat = (time_values[inc_list[ind]] - time_values[inc_list[ind-1]])  ## output-timestep too large to resolve disloc-acceleration
            # currentdata=(data[inc_list[ind]][:,slipsystem] - data[inc_list[ind-1]][:,slipsystem])/deltat/1e3
            currentdata=data[inc_list[ind]][:,slipsystem]/1e3
        elif pl_typ=='tau':
            currentdata=data[inc_list[ind]][:,slipsystem]
        elif pl_typ=='tau_back':
            currentdata=data[inc_list[ind]][:,slipsystem]
        
        if nmrks==0:
            markers[ii]=None
        if np.any(np.isnan(currentdata)):
            print(f"warning: nan encountered in {figname[1:]} at timestamp {ind}")
        else:
            pyl.plot(Xdata, currentdata, '-', color = linecolors[ii], marker=markers[ii],
                markevery = markevery, label = r'$\mathregular{t=%5.1f ns}$' % (time_values[inc_list[ind]]/1.e-9) )
            if minmax and pl_typ=='dmb_vel':
                local_min = argrelextrema(currentdata, np.less)
                local_max = argrelextrema(currentdata, np.greater)
                mindata = currentdata[local_min]
                mindata = min(mindata[abs(mindata)>1e-2])
                minind = np.where(currentdata==mindata)
                maxdata = max(currentdata[local_max])
                maxind = np.where(currentdata==maxdata)
                locmin[0,ii] = Xdata[minind]
                locmin[1,ii] = currentdata[minind]
                locmax[0,ii] = Xdata[maxind]
                locmax[1,ii] = currentdata[maxind]
                
    if minmax and pl_typ=='dmb_vel':
        pyl.plot(locmin[0],locmin[1], '--', color = 'gray')
        pyl.plot(locmax[0],locmax[1], '--', color = 'gray')
        pyl.plot(Xdata,np.repeat(min(locmin[1]),len(Xdata)), '-.', color = 'lightgray')
        pyl.plot(Xdata,np.repeat(min(locmax[1]),len(Xdata)), '-.', color = 'lightgray')
        pyl.plot(np.repeat(locmin[0][np.where(locmin[1]==min(locmin[1]))],2),[min(currentdata),max(currentdata)], ':', color = 'lightgray')
        pyl.plot(np.repeat(locmax[0][np.where(locmax[1]==min(locmax[1]))],2),[min(currentdata),max(currentdata)], ':', color = 'lightgray')
    
    if annotate:
        # annotate plot axes
        ax.set_xlabel(r'$\mathregular{Position:\/mm  }$',  fontsize=txt_size, family='serif')
        ax.set_ylabel(ylabel,  fontsize=txt_size, family='serif')
        pyl.setp(ax.get_yticklabels(), fontsize=txt_size, family='serif')
        pyl.setp(ax.get_xticklabels(), fontsize=txt_size, family='serif') 
        # Create legend (only if there is more than one curve)
        ax.legend(prop={'size':0.8*txt_size, 'family':'serif'}, **legendops, borderpad=0.3, borderaxespad=0.75, labelspacing=0.005,
                    markerscale=0.7, numpoints=1, handlelength=1,shadow=True)             
    
    # Optional stuff (often unecessary)
    if setlimits:
        # set axis limits (not necessary, default is often good enough
        ax.set_xlim(0.0, L0)
        # ax.set_xlim(0.0, 0.5)
        # ax.set_ylim(250,450)
        if pl_typ=='tau_back':
            ax.set_ylim(-0.4, 0.4)
            # ax.set_ylim(-5., 0.1)
        if pl_typ=='rho_pos' or pl_typ=='rho_neg':
            ax.set_ylim(np.min(currentdata[10:])-10,np.max(currentdata[10:])+10)
        # Grid - if needed
        # pyl.grid()

    if adjustticks: #adjust tick marks as desired
        # ax.set_yticks([1.e-12, 1.e-10, 1.e-8, 1.e-6, 1.e-4, 1.e-2, 1.])
        # ax.xaxis.set_minor_locator(AutoMinorLocator())
        # ax.yaxis.set_minor_locator(AutoMinorLocator())  
        ax.tick_params(which='both',  width  = 1)
        ax.tick_params(which='major', length = 4)
        ax.tick_params(which='minor', length = 2)

    # Adjust these values to shift the plot region inside the figure. The values reflect the location of
    #  the left, right, top, and bottom of the plot axes expressed as a fraction of the total figure size
    pyl.subplots_adjust(left=0.22, right=0.98, top=0.975, bottom=0.15,wspace=0.0,hspace=0.0)
    
    # save the figure in a bitmap and/or vector graphics form
    # pyl.savefig(filename + '.png', dpi=300)
    # pyl.savefig(filename + '.eps'         )
    pyl.savefig(filename + '.pdf', format='pdf', bbox_inches='tight')
    pyl.close()


def plot_history(pl_typ, node_index=-1, fs_direction=fs_direction,stresscomp=1,annotate=True,setlimits=False,adjustticks=False):
    fig = pyl.figure(figsize=(fw, fh))
    ax  = pyl.subplot(111) 
    if pl_typ=='dmb_vel':
        filename = "./" + jobname + ".fig.fs_velocity_x{}".format(fs_direction)
        ylabel = r'$\mathregular{Velocity:\/ m/s }$'
    elif pl_typ=='stress':
        filename = "./" + jobname + ".fig.fs_stress_{}".format(stresscomp)
        ylabel = r'$\mathregular{Stress:\/ GPa }$'
    elif pl_typ=='devstress':
        filename = "./" + jobname + ".fig.fs_devstress_{}".format(stresscomp)
        ylabel = r'$\mathregular{Deviatoric \/Stress:\/ GPa }$'
  
    # markevery = stress[inc_list[0]].shape[0] / nmrks
    
    free_surface = generate_dof_history(jobname, node_index=node_index, dof=pl_typ)
    if pl_typ=='dmb_vel':
        pyl.plot(free_surface[:,0] / 1.e-6, free_surface[:,fs_direction] /1.e3, '-', color = linecolors[0] , label="position={:.2f}mm".format(Xref[node_index]))
    if pl_typ=='stress' or pl_typ=='devstress':
        pyl.plot(free_surface[:,0] / 1.e-6, np.abs(free_surface[:,stresscomp]) /1.e3, '-', color = linecolors[0] , label="position={:.2f}mm".format(Xref[node_index]))
    
    if annotate:
        # annotate plot axes
        ax.set_xlabel(r'$\mathregular{Time:    \/ \mu s  }$',  fontsize=txt_size, family='serif')
        ax.set_ylabel(ylabel,  fontsize=txt_size, family='serif')
        pyl.setp(ax.get_yticklabels(), fontsize=txt_size, family='serif')
        pyl.setp(ax.get_xticklabels(), fontsize=txt_size, family='serif') 
    
    # Create legend (only if there is more than one curve)
    ax.legend(prop={'size':0.8*txt_size, 'family':'serif'}, loc=0, borderpad=0.3, borderaxespad=0.75, labelspacing=0.005,
                  markerscale=1.0, numpoints=1, handlelength=1,shadow=True)             
    
    # Optional stuff (often unecessary)
    if setlimits:
        # set axis limits (not necessary, default is often good enough
        ax.set_xlim(0.05, 0.15)
        # ax.set_ylim(-5., 0.1)
        # ax.set_xlim(0.15, 0.2)
        # ax.set_ylim(-0.3, 0.1)
        
    # Grid - if needed
    #pyl.grid()
    
    if adjustticks: #adjust tick marks as desired
        ax.set_xticks(np.linspace(0,max(time_values.values())*1e6,6))
        # ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4])
        #ax.xaxis.set_minor_locator(AutoMinorLocator())
        #ax.yaxis.set_minor_locator(AutoMinorLocator())  
        ax.tick_params(which='both',  width  = 1)
        ax.tick_params(which='major', length = 4)
        ax.tick_params(which='minor', length = 2)

    # Adjust these values to shift the plot region inside the figure. The values reflect the location of
    #  the left, right, top, and bottom of the plot axes expressed as a fraction of the total figure size
    pyl.subplots_adjust(left=0.20, right=0.98, top=0.975, bottom=0.15,wspace=0.0,hspace=0.0)
    
    # save the figure in a bitmap and vector graphics form
    # pyl.savefig(filename + '.png', dpi=300)
    # pyl.savefig(filename + '.eps'         )
    pyl.savefig(filename + '.pdf', format='pdf', bbox_inches='tight')
    pyl.close()


####################################
print("generating plots")

if plot_T: # Plot temperature
    plot_snapshots(tmp,'T')

if plot_dmb_vel: # Plot velocity
    if vel_direction=='all':
        for i in [1,2,3]:
            plot_snapshots(dmb_vel,'dmb_vel',vel_dir=i)
    else:
        plot_snapshots(dmb_vel,'dmb_vel',vel_dir=vel_direction,minmax=False)

if plot_stress: # Plot stress
    if stress_comp=='all':
        for i in [1,2,3,4,5,6]:
            plot_snapshots(stress,'stress',stresscomp=i)
    else:
        plot_snapshots(stress,'stress',stresscomp=stress_comp)

if plot_pressure:
    plot_snapshots(stress,'pressure')

if plot_rho_pos: # Plot positive dislocation density
    for ich in characters:
        if slipsystem=='all':
            for i in range(Nslip):
                plot_snapshots(rho_pos[ich],'rho_pos',char=ich,slipsystem=i+1)
        else:
            plot_snapshots(rho_pos[ich],'rho_pos',char=ich,slipsystem=slipsystem)

if plot_rho_neg: # Plot negative dislocation density
    for ich in characters:
        if slipsystem=='all':
            for i in range(Nslip):
                plot_snapshots(rho_neg[ich],'rho_neg',char=ich,slipsystem=i+1)
        else:
            plot_snapshots(rho_neg[ich],'rho_neg',char=ich,slipsystem=slipsystem)

if plot_tau: # Plot Resolved shear stress
    if slipsystem=='all':
        for i in range(Nslip):
            plot_snapshots(tau,'tau',slipsystem=i+1)
    else:
        plot_snapshots(tau,'tau',slipsystem=slipsystem)

if plot_tau_back: # Plot back stress
    if slipsystem=='all':
        for i in range(Nslip):
            plot_snapshots(tau_back,'tau_back',slipsystem=i+1)
    else:
        plot_snapshots(tau_back,'tau_back',slipsystem=slipsystem)
      
if plot_dis_vel: # Plot dislocation velocity magnitudes per slip system / character
    for ich in characters:
        if slipsystem=='all':
            for i in range(Nslip):
                plot_snapshots(dis_vel[ich],'dis_vel',char=ich,slipsystem=i+1)
        else:
            plot_snapshots(dis_vel[ich],'dis_vel',char=ich,slipsystem=slipsystem)

if plot_dis_acc: # Plot dislocation acceleration magnitudes per slip system / character
    for ich in characters:
        if slipsystem=='all':
            for i in range(Nslip):
                # plot_snapshots(dis_vel[ich],'dis_acc',char=ich,slipsystem=i+1)
                plot_snapshots(dis_acc[ich],'dis_acc',char=ich,slipsystem=i+1)
        else:
            # plot_snapshots(dis_vel[ich],'dis_acc',char=ich,slipsystem=slipsystem)
            plot_snapshots(dis_acc[ich],'dis_acc',char=ich,slipsystem=slipsystem)

if plot_dis_vel_x: # Plot dislocation velocities projected onto the x-direction
    for ich in characters:
        if slipsystem=='all':
            for i in range(Nslip):
                plot_snapshots(dis_vel_x[ich],'dis_vel_x',char=ich,slipsystem=i+1)
        else:
            plot_snapshots(dis_vel_x[ich],'dis_vel_x',char=ich,slipsystem=slipsystem)

if plot_dis_vel_y: # Plot dislocation velocities projected onto the x-direction
    for ich in characters:
        if slipsystem=='all':
            for i in range(Nslip):
                plot_snapshots(dis_vel_y[ich],'dis_vel_y',char=ich,slipsystem=i+1)
        else:
            plot_snapshots(dis_vel_y[ich],'dis_vel_y',char=ich,slipsystem=slipsystem)

if plot_dis_vel_z: # Plot dislocation velocities projected onto the x-direction
    for ich in characters:
        if slipsystem=='all':
            for i in range(Nslip):
                plot_snapshots(dis_vel_z[ich],'dis_vel_z',char=ich,slipsystem=i+1)
        else:
            plot_snapshots(dis_vel_z[ich],'dis_vel_z',char=ich,slipsystem=slipsystem)

if plot_av_dis_vel_x: # Plot average dislocation velocity in the x-direction ## TODO: plot also y,z directions for shear problem (requires additional OUTPUT.F90 code)
    for ich in characters:
        # plot_snapshots(dis_vel_x[ich],'av_dis_vel_x',char=ich)
        weights = rho_pos[ich].copy() # use rho_pos as weights below
        for k in weights.keys(): ## uncomments to use sum of rho_pos/neg as weihgts below
            weights[k] += rho_neg[ich][k]
        plot_snapshots(dis_vel_x[ich],'av_dis_vel_x',char=ich,weights=weights) ## weighted average

if plot_fs_vel: # Plot Free Surface Velocity
    plot_history('dmb_vel',fs_direction=fs_direction,setlimits=False,node_index=fs_position)
  
if plot_fs_stress: # Plot Free Surface Stress
    plot_history('stress',stresscomp=1,setlimits=False,node_index=fs_position)

if plot_dev_stress: # begin plot functions

  # ----------------------- #
  # Plot Deviatoric Stress  # <----------------------------------------------------------------------
  # ----------------------- #

  fig = pyl.figure(figsize=(fw, fh))
  ax  = pyl.subplot(111) 
  filename = "./" + jobname + ".fig.deviatoric_stress"
  
  # markevery = stress[inc_list[0]].shape[0] / nmrks
  
  for ii, ind in enumerate(timestamps):
      
    cauchy     = stress[inc_list[ind]][:,1:]
    pressure   = -np.sum(cauchy[:,:3], axis=1) / 3.
    dev        = cauchy
    dev[:,:3] += np.tile(pressure, (3,1)).T
    
    if np.any(np.isnan(dev[:,:3])):
        print(f"warning: nan encountered in fig.deviatoric_stress at timestamp {ind}")
    pyl.plot(Xip, dev[:,0], '-', color=linecolors[ii], label = r'$\sigma^{\prime}_{11}$, ' + r'$\mathregular{t=%5.1f ns}$' % (time_values[inc_list[ind]]/1.e-9) )
    pyl.plot(Xip, dev[:,1], '--', color=linecolors[ii], label = r'$\sigma^{\prime}_{22}$')
    pyl.plot(Xip, dev[:,2], ':', color=linecolors[ii], label = r'$\sigma^{\prime}_{33}$')
    # elif ii == 0:
    #   pyl.plot(Xip, dev[:,0], 'k-', label = r'$\sigma^{\prime}_{11}$' )
    #   pyl.plot(Xip, dev[:,1], 'b-', label = r'$\sigma^{\prime}_{22}$' )
    #   pyl.plot(Xip, dev[:,2], 'r-', label = r'$\sigma^{\prime}_{33}$' )
    # else:
    #   pyl.plot(Xip, dev[:,0], 'k--' )
    #   pyl.plot(Xip, dev[:,1], 'b--' )
    #   pyl.plot(Xip, dev[:,2], 'r--' )
      

  if 1:
    # annotate plot axes
    ax.set_xlabel(r'$\mathregular{Position:\/mm  }$',  fontsize=txt_size, family='serif')
    ax.set_ylabel(r'$\mathregular{Stress:  \/MPa }$',  fontsize=txt_size, family='serif')
    pyl.setp(ax.get_yticklabels(), fontsize=txt_size, family='serif')
    pyl.setp(ax.get_xticklabels(), fontsize=txt_size, family='serif') 
  
    # Create legend (only if there is more than one curve)
    ax.legend(prop={'size':0.8*txt_size, 'family':'serif'}, loc='upper left', bbox_to_anchor=(1.01,1), borderpad=0.3, borderaxespad=0.75, labelspacing=0.005,
                  markerscale=1.0, numpoints=1, handlelength=1.5,shadow=True)             
  
  # Optional stuff (often unecessary)
  if 1:
    # set axis limits (not necessary, default is often good enough
    ax.set_xlim(0.0, L0*1.05)
##    ax.set_ylim(-200., 120)

    # Grid - if needed
    #pyl.grid()

  if 0: #adjust tick marks as desired
      
    ax.set_yticks([1.e-12, 1.e-10, 1.e-8, 1.e-6, 1.e-4, 1.e-2, 1.])
    
    #ax.xaxis.set_minor_locator(AutoMinorLocator())
    #ax.yaxis.set_minor_locator(AutoMinorLocator())  

    ax.tick_params(which='both',  width  = 1)
    ax.tick_params(which='major', length = 4)
    ax.tick_params(which='minor', length = 2)

  # Adjust these values to shift the plot region inside the figure. The values reflect the location of
  #  the left, right, top, and bottom of the plot axes expressed as a fraction of the total figure size
  pyl.subplots_adjust(left=0.22, right=0.98, top=0.975, bottom=0.15,wspace=0.0,hspace=0.0)

  # save the figure in a bitmap and vector graphics form
  # pyl.savefig(filename + '.png', dpi=300)
  # pyl.savefig(filename + '.eps'         )
  pyl.savefig(filename + '.pdf', format='pdf', bbox_inches='tight')
  pyl.close()


        
  # ----------------------- #
  # Plot Deviatoric Stress vs model  # <----------------------------------------------------------------------
  # ----------------------- #
if plot_model_stress:
  fig = pyl.figure(figsize=(fw*1.1, fh))
  ax  = pyl.subplot(111) 
  filename = "./" + jobname + ".fig.deviatoric_stress_Xmodel"
    
  for ii, ind in enumerate(timestamps):
      
    cauchy     = stress[inc_list[ind]][:,1:]
    pressure   = -np.sum(cauchy[:,:3], axis=1) / 3.
    dev        = cauchy
    dev[:,:3] += np.tile(pressure, (3,1)).T
    devx = np.abs(dev[:,0])
    xmaxind = np.where(max(devx)==devx)[0][0]
    xmax = Xip[xmaxind]
    
    if np.any(np.isnan(dev[:,:3])):
        print(f"warning: nan encountered in fig.deviatoric_stress at timestamp {ind}")
    pyl.plot(Xip, np.abs(dev[:,0]), '-', color=linecolors[ii], label = r'$\sigma^{\prime}_{11}$, ' + r'$\mathregular{t=%5.1f ns}$' % (time_values[inc_list[ind]]/1.e-9) )
  ii=0 ## skip the first nanosecond, then read the highest deviatoric stress in the 11 direction
  while time_values[inc_list[ii]]<1e-9:
    ii+=1
  modelcauchy     = stress[inc_list[ii]][:,1:] ## look at stress at left boundary at first time-step that the code wrote to a file
  modelpressure   = -np.sum(modelcauchy[:,:3], axis=1) / 3.
  modeldev        = modelcauchy
  modeldev[:,:3] += np.tile(modelpressure, (3,1)).T
  modelstress = 1e6*np.max(np.abs(modeldev[:,0]))
  # print("modelstress: {:.3e}, time={:.2e}, time-index={}".format(modelstress,time_values[inc_list[ii]],ii))
  pos = np.linspace(0,L0*1e-3,1000) ## m
  posdevxi = np.min(np.where(pos>int(1000*xmax)/1e6))
  posdevx = pos[posdevxi]
  maxdevx = 1e6*max(devx)
  # print(maxdevx/1e6,xmax,1e3*posdevx)
  def maxmodel(x):
      return abs(float(modelprec(np.asarray([pos[0],posdevx]), modelstress, x)[-1]) - maxdevx)
  modelrhom = float(fmin(maxmodel,1e11,disp=False)[0])
  # print("initial total disloc density:{:.2e}={:.1f}".format(1e6*inputdata['rho0'],1e6*inputdata['rho0']/modelrhom),"times modelrhom")
  prec = modelprec(pos,modelstress,modelrhom)
  pyl.plot(1e3*pos,prec/1e6,'-.',color='gray', label=r"Model: $\rho_m$="+"{:.2e}".format(modelrhom))
      

  if 1:
    # annotate plot axes
    ax.set_xlabel(r'$\mathregular{Position:\/mm  }$',  fontsize=txt_size, family='serif')
    ax.set_ylabel(r'$\sigma^\prime_{11}\mathregular{:  \/MPa }$',  fontsize=txt_size, family='serif')
    pyl.setp(ax.get_yticklabels(), fontsize=txt_size, family='serif')
    pyl.setp(ax.get_xticklabels(), fontsize=txt_size, family='serif') 
  
    # Create legend (only if there is more than one curve)
    ax.legend(prop={'size':0.8*txt_size, 'family':'serif'}, loc='upper left', bbox_to_anchor=(1.01,1), borderpad=0.3, borderaxespad=0.75, labelspacing=0.005,
                  markerscale=1.0, numpoints=1, handlelength=1.5,shadow=True)             
  
  # Optional stuff (often unecessary)
  if 1:
    # set axis limits (not necessary, default is often good enough
    ax.set_xlim(0.0, L0*1.05)
##    ax.set_ylim(-200., 120)

  if 0: #adjust tick marks as desired
      
    ax.set_yticks([1.e-12, 1.e-10, 1.e-8, 1.e-6, 1.e-4, 1.e-2, 1.])
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())  

    ax.tick_params(which='both',  width  = 1)
    ax.tick_params(which='major', length = 4)
    ax.tick_params(which='minor', length = 2)

  # Adjust these values to shift the plot region inside the figure. The values reflect the location of
  #  the left, right, top, and bottom of the plot axes expressed as a fraction of the total figure size
  pyl.subplots_adjust(left=0.22, right=0.98, top=0.975, bottom=0.15,wspace=0.0,hspace=0.0)

  # save the figure in a bitmap and vector graphics form
  # pyl.savefig(filename + '.png', dpi=300)
  # pyl.savefig(filename + '.eps'         )
  pyl.savefig(filename + '.pdf', format='pdf', bbox_inches='tight')
  pyl.close()


