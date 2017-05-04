import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
import scipy.constants as const
import os, sys

def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

def loadconnect(pathdir='./', mcinif='mcini', oldvers=False, experimental=False):
    lib_path = os.path.abspath(pathdir)
    sys.path.append(lib_path)

    import dataread as dr
    if oldvers:
        import mcpickle as mcp
    else:
        if sys.version_info[0]==3:
            import mcpickle3 as mcp
        else:
            import mcpickle2 as mcp
    import infilt as cinf
    if experimental:
        import partdyn_d5 as pdyn
    else:
        import partdyn_d2 as pdyn
    mc = __import__(mcinif)
    import vG_conv as vG

    return(dr,mc,mcp,pdyn,cinf,vG)

def preproc_echoRD(mc, dr, mcp, pickfile='test.pickle'):
    mc=dr.dataread_caos(mc)
    mcp.mcpick_in(mc,pickfile)
    return mc

def pickup_echoRD(mc, mcp, dr, pickfile='test.pickle'):
    mcp.mcpick_out(mc,pickfile)
    [mc,particles,npart]=dr.particle_setup(mc)
    precTS=pd.read_csv(mc.precf, sep=',',skiprows=3)

    return(mc,particles,npart,precTS)

def particle_setup_obs(theta_obs,mc,vG,dr,pdyn):
    moistdomain=np.tile(theta_obs,int(mc.mgrid.latgrid)).reshape((mc.mgrid.latgrid,mc.mgrid.vertgrid)).T
        
    # define particle size
    # WARNING: as in any model so far, we have a volume problem here. 
    #          we consider all parts of the domain as static in volume at this stage. 
    #          however, we will work on a revision of this soon.
    mc.gridcellA=mc.mgrid.vertfac*mc.mgrid.latfac
    mc.particleA=abs(mc.gridcellA.values)/(2*mc.part_sizefac) #assume average ks at about 0.5 as reference of particle size
    mc.particleD=2.*np.sqrt(mc.particleA)/np.pi
    mc.particleV=3./4.*np.pi*(mc.particleD/2.)**3.
    mc.particlemass=dr.waterdensity(np.array(20),np.array(-9999))*mc.particleV #assume 20C as reference for particle mass
                                                                            #DEBUG: a) we assume 2D=3D; b) change 20C to annual mean T?

    # define macropore capacity based on particle size
    # we introduce a scale factor for converting macropore space and particle size
    mc.maccap=np.round(mc.md_area/((mc.particleD**2)*np.pi*mc.macscalefac)).astype(int)

    # convert theta to particles
    # npart=moistdomain*(2*mc.part_sizefac)
    npart=np.floor(mc.part_sizefac*vG.thst_theta(moistdomain,mc.soilmatrix.ts[mc.soilgrid.ravel()-1].reshape(np.shape(mc.soilgrid)), mc.soilmatrix.tr[mc.soilgrid.ravel()-1].reshape(np.shape(mc.soilgrid)))).astype(int)

    # setup particle domain
    particles=pd.DataFrame(np.zeros(int(np.sum(npart))*8).reshape(int(np.sum(npart)),8),columns=['lat', 'z', 'conc', 'temp', 'age', 'flag', 'fastlane', 'advect'])
    particles['cell']=pd.Series(np.zeros(int(np.sum(npart)),dtype=int), index=particles.index)
    # distribute particles
    k=0
    npartr=npart.ravel()
    cells=len(npartr)

    for i in np.arange(cells):
        j=int(npartr[i])
        particles.cell[k:(k+j)]=i
        rw,cl=np.unravel_index(i,(mc.mgrid.vertgrid,mc.mgrid.latgrid))
        particles.lat[k:(k+j)]=(cl+np.random.rand(j))*mc.mgrid.latfac.values
        particles.z[k:(k+j)]=(rw+np.random.rand(j))*mc.mgrid.vertfac.values
        k+=j

    particles.fastlane=np.random.randint(len(mc.t_cdf_fast.T), size=len(particles))
    particles.advect=pdyn.assignadvect(int(np.sum(npart)),mc,particles.fastlane.values,True)

    mc.mgrid['cells']=cells
    return [mc,particles.iloc[0:k,:],npart]


def CAOSpy_rundx1(tstart,tstop,mc,pdyn,cinf,precTS,particles,leftover=0,drained=0.1,dt_max=1.,splitfac=10,prec_2D=False,maccoat=10.,exfilt_method='Ediss',saveDT=True,vertcalfac=1.,latcalfac=1.,clogswitch=False,infilt_method='MDA',film=True,infiltscale=False,dynamic_pedo=False,counteractpow=3):
    if run_from_ipython():
        from IPython import display

    timenow=tstart
    prec_part=0. #precipitation which is less than one particle to accumulate
    acc_mxinf=0. #matrix infiltration may become very small - this shall handle that some particles accumulate to infiltrate
    exfilt_p=0. #exfiltration from the macropores
    s_red=0.
    mc.splitfac=splitfac
    mc.counteract_pow=counteractpow

    if type(drained)==float:
        drained=pd.DataFrame(np.array([]))
    #loop through time
    while timenow < tstop:
        [thS,npart]=pdyn.gridupdate_thS(particles.lat,particles.z,mc)
        if saveDT==True:
            #define dt as Courant/Neumann criterion
            dt_D=(mc.mgrid.vertfac.values[0])**2 / (6*np.amax(mc.D[particles.LTEbin.max(),:]))
            dt_ku=-mc.mgrid.vertfac.values[0]/np.amax(mc.ku[particles.LTEbin.max(),:])
            dt=np.amin([dt_D,dt_ku,dt_max,tstop-timenow])
            dt=np.amax([dt,0.5])
            mc.dt=dt
        else:
            if type(saveDT)==float:
                #define dt as pre-defined
                dt=np.amin([saveDT,tstop-timenow])
                mc.dt=dt
            elif type(saveDT)==int:
                #define dt as modified  Corant/Neumann criterion
                dt_D=(mc.mgrid.vertfac.values[0])**2 / (6*np.nanmax(mc.D[np.amax(thS),:]))*saveDT
                dt_ku=-mc.mgrid.vertfac.values[0]/np.nanmax(mc.ku[np.amax(thS),:])*saveDT
                dt=np.amin([dt_D,dt_ku,dt_max,tstop-timenow])
                mc.dt=dt
        
        pcount=len(particles)
        #INFILTRATION (new version with binned)
        [p_inf,prec_part,acc_mxinf]=cinf.pmx_infilt(timenow,precTS,prec_part,acc_mxinf,thS,mc,pdyn,mc.dt,leftover,prec_2D,particles.index[-1],infilt_method,infiltscale,npart) #drain all ponding // leftover <-> 0.
        particles=pd.concat([particles,p_inf])
        
        #DEBUG:
        #if len(particles)<100000:
        #    particles.to_pickle('DEBUG_particlesX.pic')
        #    assert False
        #particles.to_pickle('DEBUG_particles.pic')
        #particles=particles.groupby('cell').apply(pdyn.binupdate_pd,mc)
        #particles.loc[particles.LTEbin>(len(mc.ku)-1),'LTEbin']=len(mc.ku)-1
        #plot_mac(particles,'./results/Mtest'+str(np.round(timenow,2)))

        #DIFFUSION
        #[particles,thS,npart,phi_mx]=pdyn.part_diffusion_split(particles,npart,thS,mc,dt,False,splitfac,vertcalfac,latcalfac)
        #[particles,thS,npart,phi_mx]=pdyn.part_diffusion_binned(particles,npart,thS,mc,dt)
        
        [particles,thS,npart,phi_mx]=pdyn.part_diffusion_binned_pd(particles,npart,thS,mc)
        #ADVECTION
        if not particles.loc[(particles.flag>0) & (particles.flag<len(mc.maccols)+1)].empty:
            [particles,s_red,exfilt_p,mc]=pdyn.mac_advectionX(particles,mc,thS,mc.dt,clogswitch,maccoat,exfilt_method,film=film,dynamic_pedo=dynamic_pedo,npart=npart,fc_check=False)
        #INTERACT
        particles=pdyn.mx_mp_interact_nobulk(particles,npart,thS,mc,mc.dt,dynamic_pedo)

        #CLEAN UP DATAFRAME
        drained=drained.append(particles[particles.flag==len(mc.maccols)+1])
        particles=particles.loc[particles.flag!=len(mc.maccols)+1]
        #leftover+=np.intp((particles.z>=0.).sum())
        leftover=particles.z[particles.z>=0.].count()
        particles=particles.loc[particles.z<0.] #drop all leftover back to infilt routine

        if run_from_ipython():
            display.clear_output()
            display.display_pretty(''.join(['time: ',str(timenow),'s  |  n_particles: ',str(pcount),'  |  precip: ',str(len(p_inf)),'  |  advect: ',str(np.sum(((particles.flag>0) & (particles.flag<=len(mc.maccols))))),'  |  exfilt: ',str(int(exfilt_p)),'  |  mean v(adv): ',str(particles.loc[((particles.flag>0) & (particles.flag<=len(mc.maccols))),'advect'].mean())[:7]+str(particles.loc[((particles.flag>0) & (particles.flag<=len(mc.maccols))),'advect'].mean())[-4:],' m/s']))
            display.display_pretty(''.join(['time: ',str(timenow),'s  |  n_particles now: ',str(len(particles)),'  |  leftover: ',str(leftover)]))
        else:
            print('time: ',timenow,'s')

        #particles.loc[particles.cell<0,'cell']=mc.mgrid.cells.values
        #particles=particles[pondparts]
        timenow=timenow+mc.dt

    return(particles,npart,thS,leftover,drained,timenow)

def part_store(particles,mc):
    dummy_n2=pd.Series(np.zeros(mc.mgrid.cells),index=np.arange(mc.mgrid.cells))
    dummy_o2=pd.Series(np.zeros(mc.mgrid.cells),index=np.arange(mc.mgrid.cells))
    dummy_n=particles.loc[((particles.age>0.)),'cell'].value_counts().sort_index()
    dummy_o=particles.loc[((particles.age<=0.)),'cell'].value_counts().sort_index()
    dummy_n2.loc[dummy_n.index]+=dummy_n
    dummy_o2.loc[dummy_o.index]+=dummy_o
    dummy=pd.concat([dummy_n2,dummy_o2],axis=1)
    dummy.columns=['new','old']
    return dummy.values

#define plot routine for macropore states
def plot_mac(particles,savef=False):
    import seaborn as sns
    sns.set(font='Roboto',style='white', palette='deep',font_scale=1.2)
    zres=50.
    tmacs=int(particles.flag.max())
    for i in particles.flag.unique().astype(np.intp):
        if any(particles.flag==(i+1)):
            plt.subplot(1,tmacs,i+1)
            dummy=np.zeros((np.intp(zres),2))
            dummy[:,0]=np.linspace(0., -1., num=np.intp(zres), endpoint=False)
            filling=np.bincount(np.fmax(0,np.floor(-particles.loc[particles.flag==(i+1),'z']*zres).astype(np.intp)))
            dummy[:len(filling),1]=filling
            plt.plot([0.,0.],[0.,-1.],'grey')
            plt.fill_betweenx(dummy[:,0],0,dummy[:,1],alpha=0.4)
            plt.plot(dummy[:,1],dummy[:,0],'-')
            plt.xlim(0,50)
            plt.ylim(-1.,0.)
            plt.title(str(i))
            plt.xticks([])
            plt.yticks([])
            sns.despine(left=True, bottom=True)
    if savef:
        plt.savefig(savef+'.pdf')


