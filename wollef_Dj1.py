import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
import os, sys
try:
   import cPickle as pickle
except:
   import pickle
#connect echoRD Tools
pathdir='../echoRD/' #path to echoRD
lib_path = os.path.abspath(pathdir)
sys.path.append(lib_path)
import vG_conv as vG
from hydro_tools import plotparticles_t,hydroprofile,plotparticles_specht



# Prepare echoRD

#connect to echoRD
import run_echoRD as rE
#connect and load project
[dr,mc,mcp,pdyn,cinf,vG]=rE.loadconnect(pathdir='../',mcinif='mcini_wollef_Dj1',experimental=True)
mc = mcp.mcpick_out(mc,'wollef_Dj1.pickle')

runname='Wollef_Dj1'

mc.advectref='Shipitalo'
mc.soilmatrix=pd.read_csv(mc.matrixbf, sep=' ')
mc.soilmatrix['m'] = np.fmax(1-1/mc.soilmatrix.n,0.1)
mc.md_macdepth=mc.md_depth[np.fmax(2,np.sum(np.ceil(mc.md_contact),axis=1).astype(int))]
mc.md_macdepth[mc.md_macdepth<=0.]=0.065


precTS=pd.read_csv(mc.precf, sep=',',skiprows=3)

precTS.tstart=60
precTS.tend=60+3600
precTS.total=0.04
precTS.intense=precTS.total/(precTS.tend-precTS.tstart)


#use modified routines for binned retention definitions
mc.part_sizefac=500
mc.gridcellA=mc.mgrid.vertfac*mc.mgrid.latfac
mc.particleA=abs(mc.gridcellA.values)/(2*mc.part_sizefac) #assume average ks at about 0.5 as reference of particle size
mc.particleD=2.*np.sqrt(mc.particleA/np.pi)
mc.particleV=3./4.*np.pi*(mc.particleD/2.)**3.
mc.particleV/=np.sqrt(abs(mc.gridcellA.values)) #assume grid size as 3rd dimension
mc.particleD/=np.sqrt(abs(mc.gridcellA.values))
mc.particlemass=dr.waterdensity(np.array(20),np.array(-9999))*mc.particleV #assume 20C as reference for particle mass
                                                                        #DEBUG: a) we assume 2D=3D; b) change 20C to annual mean T?
mc=dr.ini_bins(mc)
mc=dr.mc_diffs(mc,np.max(np.max(mc.mxbin)))

[mc,particles,npart]=dr.particle_setup(mc)

#define bin assignment mode for infiltration particles
mc.LTEdef='instant'#'ks' #'instant' #'random'
mc.LTEmemory=mc.soilgrid.ravel()*0.

#new reference
mc.maccon=np.where(mc.macconnect.ravel()>0)[0] #index of all connected cells
mc.md_macdepth=np.abs(mc.md_macdepth)
mc.prects=False
#theta=mc.zgrid[:,1]*0.+0.273
#[mc,particles,npart]=rE.particle_setup_obs(theta,mc,vG,dr,pdyn)
[thS,npart]=pdyn.gridupdate_thS(particles.lat,particles.z,mc)
#[A,B]=plotparticles_t(particles,thS/100.,mc,vG,store=True)



# Run Model

mc.LTEpercentile=70 #new parameter


t_end=24.*3600.
saveDT=True

#1: MDA
#2: MED
#3: rand
infiltmeth='MDA'
#3: RWdiff
#4: Ediss
#exfiltmeth='RWdiff'
exfiltmeth='Ediss'
#5: film_uconst
#6: dynamic u
film=True
#7: maccoat1
#8: maccoat10
#9: maccoat100
macscale=1. #scale the macropore coating 
clogswitch=False
infiltscale=False

#mc.dt=0.11
#mc.splitfac=5
#pdyn.part_diffusion_binned_pd(particles,npart,thS,mc)

#import profile
#%prun -D diff_pd_prof.prof pdyn.part_diffusion_binned_pd(particles,npart,thS,mc)

wdir='/work/kit/iwg/oj4748/wollefN'
try:
    #unpickle:
    with open(''.join([wdir,'/results/Z',runname,'_Mstat.pick']),'rb') as handle:
        pickle_l = pickle.load(handle)
        dummyx = pickle.loads(pickle_l)
        particles = pickle.loads(dummyx[0])
        [leftover,drained,t] = pickle.loads(dummyx[1])
    print('resuming into stored run at t='+str(t)+'...')
except:
    print('starting new run...')

drained=pd.DataFrame(np.array([]))
leftover=0
output=60. #mind to set also in TXstore.index definition

dummy=np.floor(t_end/output)
t=0.
#loop through plot cycles
for i in np.arange(dummy.astype(int)):
    plotparticles_specht(particles,mc,pdyn,vG,runname,t,i,saving=True,relative=False,wdir=wdir)
    [particles,npart,thS,leftover,drained,t]=rE.CAOSpy_rundx1(i*output,(i+1)*output,mc,pdyn,cinf,precTS,particles,leftover,drained,6.,splitfac=4,prec_2D=False,maccoat=macscale,saveDT=saveDT,clogswitch=clogswitch,infilt_method=infiltmeth,exfilt_method=exfiltmeth,film=film,infiltscale=infiltscale)
    
    if i/5.==np.round(i/5.):
        with open(''.join([wdir,'/results/Z',runname,'_Mstat.pick']),'wb') as handle:
        	pickle.dump(pickle.dumps([pickle.dumps(particles),pickle.dumps([leftover,drained,t])]), handle, protocol=2)

    #pickle at reference states
    if i==29:
        with open(''.join([wdir,'/results/30_',runname,'_Mstat.pick']),'wb') as handle:
            pickle.dump(pickle.dumps([pickle.dumps(particles),pickle.dumps([leftover,drained,t])]), handle, protocol=2)        
    if i==59:
        with open(''.join([wdir,'/results/60_',runname,'_Mstat.pick']),'wb') as handle:
            pickle.dump(pickle.dumps([pickle.dumps(particles),pickle.dumps([leftover,drained,t])]), handle, protocol=2)
    if i==119:
        with open(''.join([wdir,'/results/120_',runname,'_Mstat.pick']),'wb') as handle:
            pickle.dump(pickle.dumps([pickle.dumps(particles),pickle.dumps([leftover,drained,t])]), handle, protocol=2)
    if i==239:
        with open(''.join([wdir,'/results/240_',runname,'_Mstat.pick']),'wb') as handle:
            pickle.dump(pickle.dumps([pickle.dumps(particles),pickle.dumps([leftover,drained,t])]), handle, protocol=2)
