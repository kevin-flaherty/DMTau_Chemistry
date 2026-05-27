from astropy.io import fits
import numpy as np
import disk_dco
from vis_sample import vis_sample
from single_model_dco import *
from scipy.optimize import minimize
def run_simplex_dco(model='cofid'):
    '''Run a downhill simplex fitting of data to a model for DCO+. Returns the best fit abundances (more specifically, the base-10 logarithms of the abundances) and the chi-squared of the model.

    REQUIRES: alma.dcodata.vis.fits: The visibility fits files for DCO+.
    
    param model (default = 'cofid'): A string that specifies the model that is considered. The options are:

    'cofid': The fiducial model, based on the CO-derived temperature strucure. There are two rings of DCO+, corresponding to the warm and cold pathway, with no additional modifications. OUTPUTS: Logarithm of the DCO+ abundance in the warm and cold pathway.
    
    'cofid_coldonly': Similar to 'cofid', but it only models the cold pathways. The abundance in the warm pathway is assumed to be 10^(-30). OUTPUTS: Logarithm of the DCO+ abundance in the cold pathway.
    
    'smallwarmout': Similar to 'cofid', but it assumes that the midplane temperature beyond ?? au is increased by 30% above its fiducial value. OUTPUTS: Logarithm of the DCO+ abundance in the warm and cold pathway.
    
    'midflat': Similar to 'cofid', but it assumes that the midplane temperature profile is flat with radius (i.e., a constant). OUTPUTS: Logarithm of the DCO+ abundance in the warm and cold pathway.
    
    'photod': Similar to 'cofid', but it adds an extra ring of DCO+ in the outer disk. The log(abundance) of this ring, and the central radius of this ring, are additional parameters in the model that are fit using the downhill simplex method. The extra ring is assumed to have a width of 50 au, and is vertically located between .79 < Sigma_21 < 3. OUTPUTS: Logarithm of the DCO+ abundance in the warm pathway, cold pathway, and the outer ring, and the central radius (in au) of the outer ring.
    
    'photod_coldonly': Similar to 'photod', but excludes the warm pathway. The warm pathway is assumed to have a DCO+ abundance of 10^(-30), while the rest of the parameters are maintained. OUTPUTS: Logarithm of the DCO+ abundance in the cold pathway, and the outer ring, and the central radius (in au) of the outer ring.
    
    '''

    if model == 'cofid':
        q=[-.371,-.371,0] #fiducial model, photod
        initial_guess = [-13.5,-12.5] #log(abundance of ring)

    if model == 'cofid_coldonly':
        q=[-.371,-.371,0] #fiducial model, photod, only the cold component
        initial_guess = [-12.5] #log(abundance of ring)

    if model =='smallwarmout':
        q=[-.371,1.3,2] #fiducial model, photod
        initial_guess = [-13.2,-12.8]  #log(abundance of rings)


    if model == 'midflat':
        q=[0,-.371,3] #fiducial model, photod
        initial_guess = [-13.5,-12.5] #[-11.1,] #log(abundance of ring)

    
    if model == 'photod':
        
        initial_guess = [-13.5,-12.5,-10.0,325]  #log(abundance of ring)
        

    if model == 'photod_coldonly':
        #photod, but only the cold pathway
        q = [-.371,-.371,0]

        initial_guess = [-12.5,-10,325]
       

    datfile = 'alma.dcodata'
    hdr=fits.getheader(datfile+'.vis.fits')
    nu = 2*hdr['naxis4']*hdr['gcount']-1-2295568 
    freq = (np.arange(hdr['naxis4'])+1-hdr['crpix4'])*hdr['cdelt4']+hdr['crval4']
    obsv = (hdr['restfreq']-freq)/hdr['restfreq']*2.99e5
    vsys=5.95
    chanstep = np.abs(obsv[1]-obsv[0])
    nchans = 2*np.ceil(np.abs(obsv-vsys).max()/chanstep)-2
    chanmin = -(nchans/2.-.5)*chanstep
    xoff = 0.
    yoff = 0.
    resolution = 0.03


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
              0.04,#10**(p[1]), #Mdisk
            1.,   #pp
            1.,  #Rin
            1000.,#Rout
            10**(2.444),# #Rc
            -36.0,  #inclination, which will be negative for IM Lup.
            .54,  #Mstar
            10**(-4),  #Xco
            0,  #vturb
            70., #Zq0
            14.3,#19.,  #Tmid0
            24.68, #Tatm0
            #19., #Tco
            [.79,1000], #upper and lower boundaries in column density
            [9.,800.], #inner and outer boundaries for abundance
            -1, 200]   #handed

        if model =='cofid_coldonly' or model =='photod_coldonly':
            params[8] = [10**(-30.),10**(x[0])]
        elif model=='smallwarmout':
            params[8] = [10**(x[0]),10**(-12.3)] 
        else:
            params[8] = [10**(x[0]),10**(x[1])]

        d=disk_dco.Disk(params,rtg=False)
        obs = [150,101,300,170] #150,101,280,170 rt grid nr,nphi,nz,zmax
        #disk_dco.Tco = 18. # Un-comment in order to change the CO freeze-out temperature to a different value (default = 19 K).

        d.set_obs(obs)
        d.set_rt_grid()
        d.set_line('dco')
        
        if model =='photod':
            d.add_mol_ring(x[3],x[3]+50,.79,3.,10**(x[2]),just_frozen=False) 
        if model =='photod_coldonly':
            d.add_mol_ring(x[2],x[2]+50,.79,3.,10**(x[1]),just_frozen=False) 
        total_model(disk=d,chanmin=chanmin,nchans=nchans,chanstep=chanstep,offs=[xoff,yoff],modfile='alma',imres=resolution,obsv=obsv,vsys=vsys,freq0=288.143858,Jnum=3,distance=144.5,hanning=True,PA=154.8,bin=2)

        
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

