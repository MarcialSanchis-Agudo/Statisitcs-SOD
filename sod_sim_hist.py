import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys
import os

fileName = sys.argv[1]

file_1 = 'analysis_'+fileName+'.dat';
file_2 = 'surf_code_1-'+fileName+'.dat';
file_3 = 'monitor_'+fileName+'.dat'

data = np.loadtxt(file_1,skiprows=0)
data_2 = np.loadtxt(file_2,delimiter=",",skiprows=1)

#visco   = 5.3566e-05;
Lz = 0.2;
visco   = 1.0;
lw = 1.0

startInd = 100;
time    = data[startInd::,0];
Ek      = data[startInd::,1];
eps_S   = data[startInd::,2];
eps_D   = data[startInd::,3];
eps_T   = data[startInd::,4];
Mmax    = data[startInd::,5];

time_2    = data_2[startInd::,1];
area      = data_2[startInd::,2]
Fpx       = data_2[startInd::,3];
Fpy       = data_2[startInd::,4];
Fpz       = data_2[startInd::,5];
Fvx       = data_2[startInd::,6];
Fvy       = data_2[startInd::,7];
Fvz       = data_2[startInd::,8];

if(os.path.isfile(file_3)):
  data_3 = np.loadtxt(file_3,skiprows=0)
  time_3    = data_3[startInd::,0];
  maxRho    = data_3[startInd::,1]
  maxMue    = data_3[startInd::,2]/visco;
  maxVel    = data_3[startInd::,3];

if('pipe' in fileName):
 print("--||SOD :: Calculating utau for pipe")
 utau      = np.sqrt(abs(Fvz/area))
 print("--||SOD :: time-average utau =%f",np.mean(utau,axis=None))
elif('channel' in fileName):
 print("--||SOD :: Calculating utau for channel")
 utau      = np.sqrt(abs(Fvx/area))
 print("--||SOD :: time-average utau = ",np.mean(utau,axis=None))
elif('naca' in fileName):
 print("--||SOD :: Calculating Cl/Cd for NACA")
 print("----||Fvx,Fpx,Fvy,Fpy = ",np.mean(Fvx,axis=None),np.mean(Fpx,axis=None),np.mean(Fvy,axis=None),np.mean(Fpy,axis=None))
 cdp     = 2.0*(Fpx)/Lz
 cdv     = 2.0*(Fvx)/Lz
 cd      = 2.0*(Fvx+Fpx)/Lz
 cl      = 2.0*(Fvy+Fpy)/Lz
 print("--||SOD :: time-average Cdp = ",np.mean(cdp,axis=None))
 print("--||SOD :: time-average Cdv = ",np.mean(cdv,axis=None))
 print("--||SOD :: time-average Cd  = ",np.mean(cd,axis=None))
 print("--||SOD :: time-average Cl  = ",np.mean(cl,axis=None))
else:
 print("--||SOD :: Error. fileName handling not defined!")

if(os.path.isfile(file_3)):
  #ylab = [r'$E_k$',r'$\epsilon_{S}$',r'$\epsilon_{D}$',r'$\epsilon_{T}$',r'$M_{max}$','']
  ylab = [r'$E_k$',r'$eps_{T}$',r'$\max(M)$',r'$\max(U)$',r'$\max(\mu_e)/\mu_f$',r'$\max(\rho)$']
  
  #fig, axs = plt.subplots(3, 2, sharex=True, dpi=300)
  fig = plt.figure(dpi=300)
  grid = GridSpec(4, 2, figure=fig)
  
  #------------------#
  axs = fig.add_subplot(grid[0,0])
  axs.plot(time,Ek,linewidth=lw,label='Ek')
  plt.xlabel(r'$t$')
  plt.ylabel(ylab[0])
  plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  #plt.gca().set_xticklabels([])
  #------------------#
  #------------------#
  axs = fig.add_subplot(grid[0,1])
  axs.plot(time_3,maxVel,linewidth=lw,label='')
  plt.xlabel(r'$t$')
  plt.ylabel(ylab[3])
  #plt.gca().set_xticklabels([])
  #------------------#
  axs = fig.add_subplot(grid[1,0])
  #axs.plot(time,eps_S,linewidth=lw,label=r'$eps_S$')
  axs.plot(time,eps_T,linewidth=lw,label=r'$eps_T$')
  plt.xlabel(r'$t$')
  plt.ylabel(ylab[1])
  plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  #plt.gca().set_xticklabels([])
  #------------------#
  #------------------#
  axs = fig.add_subplot(grid[1,1])
  axs.plot(time_3,maxMue,linewidth=lw,label='')
  plt.xlabel(r'$t$')
  plt.ylabel(ylab[4])
  plt.gca().set_ylim(bottom=0)
  #plt.gca().set_xticklabels([])
  #------------------#
  axs = fig.add_subplot(grid[2,0])
  #axs.plot(time,eps_D,linewidth=lw,label=r'$eps_D$')
  axs.plot(time,Mmax,linewidth=lw,label=r'$M_{max}$')
  plt.xlabel(r'$t$')
  plt.ylabel(ylab[2])
  plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  #------------------#
  #------------------#
  axs = fig.add_subplot(grid[2,1])
  axs.plot(time_3,maxRho,linewidth=lw,label='')
  plt.xlabel(r'$t$')
  plt.ylabel(ylab[5])
  #------------------#
  if('naca' in fileName):
    axs = fig.add_subplot(grid[3,0])
    axs.plot(time_2,cl,linewidth=lw,label=r'$c_{l}$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$c_{l}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))

    axs = fig.add_subplot(grid[3,1])
    axs.plot(time_2,cd,linewidth=lw,label=r'$c_{d}$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$c_{d}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  else:
    axs = fig.add_subplot(grid[3,:])
    axs.plot(time_2,utau,linewidth=lw,label=r'$u_{\tau}$')
    #axs[2,1].plot(time,Fy,linewidth=lw,label=r'$F_y$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$u_{\tau}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  #------------------#
  
  fig.tight_layout()
  plt.savefig('sod_analysis.png')
  #plt.show()
else:
  ylab = [r'$E_k$',r'$eps_{S}$',r'$eps_{D}$',r'$eps_{T}$',r'$M_{max}$','']
  
  #fig, axs = plt.subplots(3, 2, sharex=True, dpi=300)
  fig = plt.figure(dpi=300)
  grid = GridSpec(3, 2, figure=fig)
  
  axs = fig.add_subplot(grid[0,0])
  axs.plot(time,Ek,linewidth=1.5,label='Ek')
  plt.xlabel(r'$t$')
  plt.ylabel(ylab[0])
  axs = fig.add_subplot(grid[0,1])
  axs.plot(time,eps_S,linewidth=1.5,label=r'$eps_S$')
  plt.xlabel(r'$t$')
  plt.ylabel(ylab[1])
  axs = fig.add_subplot(grid[1,0])
  axs.plot(time,eps_D,linewidth=1.5,label=r'$eps_D$')
  plt.xlabel(r'$t$')
  plt.ylabel(ylab[2])
  axs = fig.add_subplot(grid[1,1])
  #axs[1,1].plot(time,eps_T,linewidth=1.5,label=r'$eps_T$')
  axs.plot(time,Mmax,linewidth=1.5,label=r'$M_{max}$')
  plt.xlabel(r'$t$')
  plt.ylabel(ylab[4])
  if('naca' in fileName):
    axs = fig.add_subplot(grid[2,0])
    axs.plot(time_2,cl,linewidth=lw,label=r'$c_{l}$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$c_{l}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))

    axs = fig.add_subplot(grid[2,1])
    axs.plot(time_2,cd,linewidth=lw,label=r'$c_{d}$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$c_{d}$')
    plt.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
  else:
    axs = fig.add_subplot(grid[2,:])
    axs.plot(time_2,utau,linewidth=1.5,label=r'$u_{\tau}$')
    #axs[2,1].plot(time,Fy,linewidth=1.5,label=r'$F_y$')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$u_{\tau}$')

  plt.tight_layout()
  plt.savefig('sod_analysis.png')
  #plt.show()

