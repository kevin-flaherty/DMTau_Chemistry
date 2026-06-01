from astropy.io import fits
import numpy as np
import disk
from vis_sample import vis_sample
from single_model import *
from scipy.optimize import minimize

def run_simplex(model='cofid'):
    '''Run a downhill simplex fitting of data to a model for N2H+. Returns the best fit abundances (more specifically, the base-10 logarithms of the abundances) and the chi-squared of the model.

    REQUIRES: alma.n2hdata.vis.fits: The visibility fits files for N2H+.
    
    param model (default = 'cofid'): A string that specifies the model that is considered. The options are:

    'cofid': The fiducial model, based on the CO-derived temperature strucure. There is one ring of N2H+, corresponding to the region of the disk where ?? < Tgas < 20 K with no additional modifications. OUTPUTS: Logarithm of the N2H+ abundance.
    
    
    'smallwarmout': Similar to 'cofid', but it assumes that the midplane temperature beyond ?? au is increased by 30% above its fiducial value. OUTPUTS: Logarithm of the N2H+ abundance.
    
    'midflat': Similar to 'cofid', but it assumes that the midplane temperature profile is flat with radius (i.e., a constant). OUTPUTS: Logarithm of the N2H+ abundance.
    
    'highphotod': Similar to 'cofid', but it adds an extra ring of N2H+ in the outer disk. The log(abundance) of this ring, and the central radius of this ring, are additional parameters in the model that are fit using the downhill simplex method. The extra ring is assumed to have a width of 50 au, and is vertically located between .01 < Sigma_21 < .79, which is above the region in the outer disk with photodesorbed CO. OUTPUTS: Logarithm of the N2H+ abundance in the inner ring and the outer ring, and the central radius (in au) of the outer ring.
    
    'lowphotod': Similar to 'photod', but puts the outer ring of N2H+ closer to the midplane, at 3 < Sigma_21 < 100. OUTPUTS: Logarithm of the N2H+ abundance in the inner ring and the outer ring, and the central radius (in au) of the outer ring.
    
    '''

    if model == 'cofid':
        q=[-.371,-.371,0] #fiducial model, photod
        initial_guess = [-11.1,] #log(abundance of ring)


    if model =='smallwarmout':
        q=[-.371,1.3,2] #fiducial model, photod
        initial_guess = [-11.1,] #log(abundance of ring)

    
    if model == 'midflat':
        q=[0,-.371,3] #fiducial model, photod
        initial_guess = [-11.1,] #log(abundance of ring)

    if model == 'highphotod':
        q=[-.371,-.371,0] #fiducial model, photod

        initial_guess = [-11.1,-10.0,325] #log(abundance of ring)

    if model == 'lowphotod':
        q=[-.371,-.371,0] #fiducial model, photod
        initial_guess = [-11.1,-10.0,325] #log(abundance of ring)


    
    datfile = 'alma.n2hdata'
    hdr=fits.getheader(datfile+'.vis.fits')
    nu = 2*hdr['naxis4']*hdr['gcount']-1-2086120 #227478
    freq = (np.arange(hdr['naxis4'])+1-hdr['crpix4'])*hdr['cdelt4']+hdr['crval4']
    obsv = (hdr['restfreq']-freq)/hdr['restfreq']*2.99e5
    vsys=5.95#5.76 from grid_search
    chanstep = np.abs(obsv[1]-obsv[0])#-0.208337
    nchans = 2*np.ceil(np.abs(obsv-vsys).max()/chanstep)-2
    chanmin = -(nchans/2.-.5)*chanstep
    xoff = 0.
    yoff = 0.
    resolution = 0.0325

    obj = fits.open(datfile+'.vis.fits')
    vis_obj = (obj[0].data['data']).squeeze()
    real_obj = (vis_obj[:,:,0,0]+vis_obj[:,:,1,0])/2.
    imag_obj = (vis_obj[:,:,0,1]+vis_obj[:,:,1,1])/2.
    weight_real = vis_obj[:,:,0,2]
    weight_imag = vis_obj[:,:,1,2]
    weight_real[real_obj==0] = 0.
    weight_imag[imag_obj==0] = 0.

    obj.close()

    

    def objective_function(x,model):
        params = [q, #qq
              0.04, #Mdisk
            1.,   #pp
            1.,  #Rin
            1000.,#Rout
            10**(2.444), #Rc
            -36.0,  #inclination, which will be negative for IM Lup.
            .54,  #Mstar
            10**(x[0]),  #Xco
            0,  #vturb
            70., #Zq0
            14.3, #Tmid0
            24.68, #Tatm0
            19., #Tco
            [.79,1000], #upper and lower boundaries in column density
            [9.,800.], #inner and outer boundaries for abundance
            -1, 200]   #handed

        d=disk.Disk(params,rtg=False)
        obs = [150,101,300,170] #150,101,280,170 rt grid nr,nphi,nz,zmax

        d.set_obs(obs)
        d.set_rt_grid()
        d.set_line('n2h32')
        
        if model =='lowphotod':
            d.add_mol_ring(x[2],x[2]+50,3,100,10**(x[1]),just_frozen=False) 
        if model == 'highphotod':
            d.add_mol_ring(x[2],x[2]+50,.01,.79,10**(x[1]),just_frozen=False) 
        total_model(disk=d,chanmin=chanmin,nchans=nchans,chanstep=chanstep,offs=[xoff,yoff],modfile='alma',imres=resolution,obsv=obsv,vsys=vsys,freq0=279.51170100,Jnum=2,distance=144.5,hanning=True,PA=154.8,bin=4)
        
        # - Generate model visibilities
        model_vis = vis_sample(uvfile=datfile+'.vis.fits',imagefile='alma.fits',mod_interp=False)
        real_model = model_vis.real
        imag_model = model_vis.imag
        chi = ((real_model-real_obj)**2*weight_real).sum() + ((imag_model-imag_obj)**2*weight_imag).sum()
        return chi/nu
    

    print('starting minimization')
    
    result = minimize(objective_function,initial_guess,method='nelder-mead',tol=.1,args=(model))
    if result.success:
        print('Minimum found at ',result.x)
        print('Minimum value: ',result.fun)
    else:
        print('Optimization failed: ',result.message)
