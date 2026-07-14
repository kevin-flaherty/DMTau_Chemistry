
from astropy.io import fits
import os
from disk_dco import *
from raytrace_dco import *
from scipy.optimize import curve_fit
import scipy.interpolate
#from scipy.integrate import cumtrapz,trapz
#from galario import double as gdouble
import time
import uuid
from vis_sample import vis_sample


##############################################################################
def make_model_vis(datfile='alma.dcodata',modfile='testpy_alma_dco',isgas=True,freq0=288.143858):
    if isgas:
        cmd = ' ./sample_alma.csh '+modfile+' '+datfile+' '+str(freq0)
    else:
        cmd = ' ./sample_cont.csh '+modfile+' '+datfile+' '+str(freq0)
    os.system(cmd)

def compare_vis_sample(datfile='alma.n2hdata',modfile='testpy_alma_n2h',new_weight=[1,],systematic=False,isgas=True,plot_resid=False):
    '''Calculate the raw chi-squared based on the difference between the model and data visibilities. This function uses vis_sample to perform this calculation.'''

    # - Read in object visibilities
    obj = fits.open(datfile+'.vis.fits')
    vis_obj = (obj[0].data['data']).squeeze()
    if isgas:
        if obj[0].header['telescop'] == 'ALMA':
            if obj[0].header['naxis3'] == 2:
                real_obj = (vis_obj[:,:,0,0]+vis_obj[:,:,1,0])/2.
                imag_obj = (vis_obj[:,:,0,1]+vis_obj[:,:,1,1])/2.
                weight_real = vis_obj[:,:,0,2]
                weight_imag = vis_obj[:,:,1,2]
            else:
                real_obj = vis_obj[::2,:,0]
                imag_obj = vis_obj[::2,:,1]
    else:
        if obj[0].header['telescop'] == 'ALMA':
            if obj[0].header['naxis3'] == 2:
                real_obj = (vis_obj[:,0,0]+vis_obj[:,1,0])/2.
                imag_obj = (vis_obj[:,0,1]+vis_obj[:,1,1])/2.
                weight_real = vis_obj[:,0,2]
                weight_imag = vis_obj[:,1,2]

    obj.close()

    # - Generate model visibilities
    model_vis = vis_sample(uvfile=datfile+'.vis.fits',imagefile=modfile+'.fits',mod_interp=False)
    real_model = model_vis.real
    imag_model = model_vis.imag

    if systematic:
        real_model = real_model/systematic
        imag_model = imag_model/systematic

    if len(new_weight)>1:
        weight_real = new_weight
        weight_imag = new_weight

    weight_real[real_obj==0] = 0.
    weight_imag[imag_obj==0] = 0.
    print('Removed data %i' % ((weight_real ==0).sum()+(weight_imag==0).sum()))

    if plot_resid:
        #Code to plot, and fit, residuals
        #If errors are Gaussian, then residuals should have gaussian shape
        #If error size is correct, residuals will have std=1
        obj = fits.open(datfile+'.vis.fits')
        freq0 = obj[0].header['crval4']
        u_obj,v_obj = (obj[0].data['UU']*freq0).astype(np.float64),(obj[0].data['VV']*freq0).astype(np.float64)
        vis_obj = (obj[0].data['data']).squeeze()
        obj.close()
        uv = np.sqrt(u_obj**2+v_obj**2)
        use = (weight_real > .05) & (weight_imag>.05)
        diff = np.concatenate((((real_model[use]-real_obj[use])*np.sqrt(weight_real[use])),((imag_model[use]-imag_obj[use])*np.sqrt(weight_imag[use]))))
        diff = diff.flatten()
        n,bins,patches = plt.hist(diff,10000,density=1,histtype='step',color='k',label='Data',lw=3)
        popt,pcov = curve_fit(gaussian,bins[1:],n)
        y=gaussian(bins,popt[0],popt[1],popt[2])
        print('Gaussian fit parameters (amp,width,center): ',popt)
        print('If errors are properly scaled, then width should be close to 1')
        plt.plot(bins,y,'r--',lw=6,label='gaussuian')
        #slight deviations from gaussian, but gaussian is still the best...
        plt.xlabel('(Model-Data)/$\sigma$',fontweight='bold',fontsize=20)
        ax=plt.gca()
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20)
            tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20)
            tick.label1.set_fontweight('bold')

        plt.show()



    chi = ((real_model-real_obj)**2*weight_real).sum() + ((imag_model-imag_obj)**2*weight_imag).sum()
    return chi


def gaussian(x,amp,width,center):
    return amp/(width*np.sqrt(2*np.pi))*np.exp(-.5*((x-center)**2)/width**2)

#def lorentzian(x,amp,width,center):
#    return 1/np.pi*(width/((x-center)**2.+width**2))*amp

#def laplace(x,amp,width,center):
#    return amp/(2*width)*np.exp(-np.abs(x-center)/width)

def lnlike(p,massprior=False,cleanup=False,systematic=False,line='dco',vcs=True,exp_temp=False,add_ring=False,gs25=False,Tco=18.):
    '''Calculate the log-likelihood (=-0.5*chi-squared) for a given model.
    
    REQUIRES: alma.dcodata.vis.fits: The visibility fits files for DCO+.
    
    PARAMS: 
    
    p (default = None): Required list of parameter values to generate the model. p = [q, log(abund), log(abund2), log(abund3), Rring, turbulence, x-offset, y-offset]. For example, p=((-.371,-.371,0),-20,-10.4,-20,200,0.,0.,0.)
    
    massprior (default = False): Include a prior on the disk mass. Note that this functionality is not used in analyzing the DCO+ data, but is included here as an example of how to include a prior on one of the input parameters.

    cleanup (default = False): Delete the model file from the current directory when done. This is used when you want the chi-squared value, but don't need to model image itself. If set to TRUE, the model is given a generic name, otherwise the file is named 'alma.fits'.

    systematic (default = False): Include a systematic error in the amplitude calibration of the data. If set to True, then the systematic error is assumed to be the last value in the p array. The systematic uncertainty is a multiplicative factor that is applied to the model. For example, a systematic value of 1.2 would mean that the true flux of the data is 20% brighter than what has been observed, with this scaling applied to the model instead of changing the data.

    line (default = 'dco'): The molecule to model. For now, only DCO+ is implemented.

    vcs (default = True): If true, then the turbulence is treated as proportional to the local thermal broadening in the disk. If false, then the turbulence is treated as a constant value, in units of km/s.

    exp_temp (default = False): If true, then the temperature is treated as an exponential function of height, as opposed to a Dartois et al. 2013 type II profile.

    add_ring (default = False): If true, then a ring is added to the model. NOTE: This is not used in analyzing the DCO+ data.

    OUTPUTS:
    The code will print out the p array, the reduced chi-squared value, and the time it took to run the code. It will return the log-likelihood value. 

    If cleanup=False then the code will retain the model image, in the file 'alma.fits'. If cleanup=True then the code will delete the model image. The units for the model image are included in the fits file header; the flux units are Jy/pixel while the pixel size is matched to that of the data. 

    '''

    start=time.time()
    all_params = {
    'q':p[0], #q
    'Mdisk':0.04, #solar masses
    'p':1, #gamma
    'Rin':1., #Model inner domain - generally no need to change
    'Rout':1000., #Model outer domain - generally no need to change
    'Rc':10**(2.444), #10^(p[1]), #Rc
    'incl':-36,#-36, #inclination, degrees
    'Mstar':0.54, #solar masses
    'Xdco':[10**(p[1]),10**(p[2]),10**(p[3])], #Abundance
    'vturb':p[5], #turbulence, as a fraction of the thermal broadening for this line
    'Zq0':70., #Zq0
    'Tmid0':14.3, #p[6], #K
    'Tatm0':24.68,#K
    'Zabund':[.79,1000], #upper and lower boundaries in column density
    'Rabund':[10.,1000.], #inner and outer boundaries for abundance
    'handed':-1, #handed
    'vsys':5.95,#6.06, #systemic velocity, km/s
    'offs':[p[6],p[7]], #position offset, arcseconds
    'PA':154.8, #position angle, degrees
    'distance':144.5, #distance
    'Rbreak':200} #Radius of temperature break
    if gs25 == True:
        ### Adjust the underlying temperature structure to match that derived by Galloway-Sprietsma et al. 2025
        all_params['Zq0'] = 56.25
        all_params['Tmid0'] = 18.8
        all_params['Tatm0'] = 34.2
    if line.lower() =='co21' or line.lower()=='co32' or line.lower()=='svco21' or line.lower()=='dco':
        params = [all_params['q'],all_params['Mdisk'],all_params['p'],all_params['Rin'],all_params['Rout'],all_params['Rc'],all_params['incl'],all_params['Mstar'],all_params['Xdco'],all_params['vturb'],all_params['Zq0'],all_params['Tmid0'],all_params['Tatm0'],all_params['Zabund'],all_params['Rabund'],all_params['handed'],all_params['Rbreak']]
    if all_params['Mdisk'] <0 or all_params['Mdisk']>all_params['Mstar'] or all_params['Rin']<0 or all_params['Rin']>all_params['Rout'] or all_params['Rout']<0 or all_params['Rc']<0 or all_params['Mstar']<0 or all_params['vturb']<0 or all_params['Zq0']<0 or all_params['Tmid0']<0 or all_params['Tmid0']>all_params['Tatm0'] or all_params['Tatm0']<0 or all_params['Zabund'][0]<0 or all_params['Zabund'][1]<0 or all_params['Zabund'][1]<all_params['Zabund'][0] or all_params['Rabund'][0]<0 or all_params['Rabund'][0]<all_params['Rin'] or all_params['Rabund'][0]>all_params['Rabund'][1] or all_params['Rabund'][1]<0 or all_params['Rabund'][1]>all_params['Rout']:
        print(all_params)
        chi = np.inf
        nu = 1
    else:
        if add_ring:
            if systematic:
                if p[-3]<0:
                    print('Bad ring parameters ',p[-3],p[-2])
                    return -np.inf
                else:
                    disk_structure=Disk(params,rtg=False,exp_temp=exp_temp,ring=[(params[3]+p[-3])/2.,p[-3]-params[3],p[-2]])
            else:
                if p[-2]<0:
                    print('Bad Ring parameters',p[-2:])
                    return -np.inf
                else:
                #disk_structure = Disk(params,rtg=False,exp_temp=exp_temp,ring=[p[-3],p[-2],p[-1]])
                    disk_structure=Disk(params,rtg=False,exp_temp=exp_temp,ring=[(params[3]+p[-2])/2.,p[-2]-params[3],p[-1]])
        else:
            obs = [180,131,300,270]
            disk_structure=Disk(params,rtg=False,exp_temp=exp_temp,obs=obs)
            disk_structure.Tco=Tco
        if cleanup:
            tf = tempfile.NamedTemporaryFile()
            modfile = tf.name[-9:]
            tf.close()
        else:
            modfile = 'alma'

        # The next series of lines are very specific to the CO3-2 ALMA data for HD 163296, and would need to be modified to use another data set. Basically you need to set the keyword datfile, read in the weights, set the degrees of freedom, set the chanmin,nchans,chanstep keywords (specific to this spectra you are trying to simulate) as well as the image offset.
        if line.lower() == 'dco':
            datfile = '~/alma.dcodata'
            hdr=fits.getheader(datfile+'.vis.fits')
            nu = 2*hdr['naxis4']*hdr['gcount']-len(p)-2295568 #227478
            freq = (np.arange(hdr['naxis4'])+1-hdr['crpix4'])*hdr['cdelt4']+hdr['crval4']
            obsv = (hdr['restfreq']-freq)/hdr['restfreq']*2.99e5
            vsys=all_params['vsys']#5.76 from grid_search
            chanstep = np.abs(obsv[1]-obsv[0])#-0.208337
            nchans = 2*np.ceil(np.abs(obsv-vsys).max()/chanstep)
            chanmin = -(nchans/2.-.5)*chanstep
            offs = all_params['offs']
            resolution = 0.03
            obs = [180,131,300,170] #180,131,300,170 rt grid nr,nphi,nz,zmax

            disk_structure.set_obs(obs)
            disk_structure.set_rt_grid(vcs=vcs)
            disk_structure.set_line(line)
            disk_structure.add_mol_ring(p[4],p[4]+50,.79,3.,all_params['Xdco'][2],just_frozen=True) #.79,3, 325,375
            total_model(disk=disk_structure,chanmin=chanmin,nchans=nchans,chanstep=chanstep,offs=offs,modfile=modfile,imres=resolution,obsv=obsv,vsys=vsys,freq0=288.143858,Jnum=3,distance=all_params['distance'],hanning=True,PA=all_params['PA'],bin=2)



        if systematic:
            sys = p[-1]#1+0.2*np.random.randn()
        else:
            sys = None

        chi = compare_vis_sample(datfile=datfile,modfile=modfile,systematic=sys)

        if cleanup:
            # Clean up files
            files = [modfile+'p.fits',modfile+'p.han.im',modfile+'p.im',modfile+'p.model.vis',modfile+'p.model.vis.fits']
            for file in files:
                os.system('rm -r '+file)

    if np.isnan(chi):
        chi = 100*nu

    print(p)
    print(chi/nu)


    if massprior:
        # Include a prior on mass with the likelihood estimate
        # Assume a gaussian prior with given mean and
        mean_mdisk = 0.09 #mean Mdisk
        sig_mdisk = 0.01 # standard deviation on prior
        lnp = -np.log(sig_mdisk*np.sqrt(2*np.pi))-(p[1]-mean_mdisk)**2/(2*sig_mdisk**2)
    else:
        lnp = 0.0

    #if systematic:
        #if np.abs(p[-1]-1)>.2:
        #    lnp = -np.inf
        #else:
        #    lnp=0.
        #lnp -= (p[-1]-1.)**2/(2*.2**2) #prior on the gain, centered at 1 with dispersion of 0.2

    print('%r minutes' % ((time.time()-start)/60.))
    return -0.5*chi+lnp
    #return chi/nu


