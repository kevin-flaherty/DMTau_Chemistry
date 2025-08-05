from astropy.io import fits
import os
import numpy as np
import matplotlib.pylab as plt
#from galario import double as gdouble
from gofish import imagecube
import disk
from vis_sample import vis_sample
from scipy.optimize import minimize

def plot_profile(file,overplot=False,norm_peak=False):
    '''file is the moment zero map created from bettermoments.
    Before calling plot_moments, run bettermoments from within python as:
        !bettermoments alma.model.fits -method zeroth -clip 3
    where alma.model.fits is the cleaned model image. The other inputs specify that you want to make a moment zero map, and that you want to clip anything below 3 sigma.
    This will create a file called alma.model_I0.fits, which is the moment zero map, and is the file that you want to use with this code.

    The shaded area represents the area within one beam of the center of the disk. This is the region that I have found is heavily dependent on the sigma clipping used when creating the moment map.
    '''
    #DM Tau distance = 144.5 pc
    cube = imagecube(file)
    x,y,dy = cube.radial_profile(inc=36.0,PA=154.8)
    y, dy = y*1e3,dy*1e3 # for continuum, convert from Jy/beam to mJy/beam
    if file[-8:] == 'Fnu.fits':
        scale=.002
        label = 'Flux (mJy/bm)'
        if norm_peak:
            dy/=y.max()
            scale/=y.max()
            y/=y.max()
            label='Normalized Flux'
    else:
        scale = 1
        label = 'Integrated Intensity (mJy bm$^{-1}$ km s$^{-1}$)'
        if norm_peak:
            dy/=y.max()
            scale/=y.max()
            y/=y.max()
            label='Normalized Intensity'

    if not overplot:
        plt.figure()
        plt.errorbar(x,y,dy,fmt=' ',capsize=1.25,capthick=1.25,color='k',lw=1.0)
        plt.step(x,y,where='mid',color='k',lw=1.0)
        plt.xlim(0,3)
        plt.xlabel('Radius (arcsec)',fontsize=14)
        plt.ylabel(label,fontsize=14)
        #plt.ylabel('Intensity (mJy bm$^{-1}$)') #for continuum
        plt.axhline(0.,ls=':',color='k')
        ax = plt.gca()
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(12)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(12)

        #plot a Gaussian beam
        xgauss = np.linspace(0.5,1.5,50)
        #xgauss = np.linspace(0.7,1.3,50) # for continuum
        sigmaj = cube.beam[0]/(2*np.sqrt(2*np.log(2))) #convert from FWHM to sigma
        #ygauss = 1*np.exp(-(xgauss-1)**2./(2*sigmaj**2.))+2 #for continuum
        ygauss = 3e3*scale*np.exp(-(xgauss-1)**2./(2*sigmaj**2.))
        plt.plot(xgauss,ygauss,color='k')

        #shade in the area within one beam of the center
        ylim = plt.ylim()
        plt.fill_between([0,cube.beam[0]],[ylim[0],ylim[0]],[ylim[1],ylim[1]],color='k',alpha=.3)
        plt.ylim(ylim)

        ax = plt.gca()
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        #ax2.set_xticks([0.0,0.5,1.0,1.5,2.0,2.5,3.0])
        #ax2.set_xticklabels(['{:0.1f}'.format(x) for x in 144.5*np.array([0.,.5,1.,1.5,2.,2.5,3.])])
        ax2.set_xticks([x for x in np.arange(0,450,100)/144.5])
        ax2.set_xticklabels(['0','100','200','300','400'],fontsize=12)
        ax2.set_xlabel('Radius (au)',fontsize=14)
    else:
        plt.errorbar(x,y,dy,fmt=' ',capsize=1.25,capthick=1.25,color='r',lw=1.0)
        plt.step(x,y,where='mid',color='r',lw=1.0,ls=':')


def write_model_vis_sample(file,modfile,model_vis_file = 'alma.model.vis.fits'):
    '''Use vis_sample and the model image to generate model visibilities, which are then written to a file.'''

    os.system('rm -rf '+model_vis_file)
    vis_sample(imagefile=modfile,uvfile=file,outfile=model_vis_file,mod_interp=False)

def im_plot_spec(file='data/HD163296.CO32.regridded.cen15.cm.fits',size=10,fwhm=False,gain=1.,threshold=None,norm_peak=False,mirror=False,line='co21',show_cloud=False,skip2 = False,**kwargs):
    'Given an image, plot the spectrum based on averaging over a given image area'



    # - Read in the data
    im,data,hdr,ra,dec,noise = get_cropped_image(file,size)
    #im = fits.open(file)
    #data = im[0].data.squeeze()
    #hdr = im[0].header
    ##noise = np.std(data[:,:,:10])
    #noise = calc_noise(data,np.round(size/(3600.*hdr['cdelt1'])))
    #ra = 3600*hdr['cdelt1']*(np.arange(hdr['naxis1'])-hdr['naxis1']/2.-0.5)
    #dec = 3600*hdr['cdelt2']*(np.arange(hdr['naxis2'])-hdr['naxis2']/2.-0.5)
    #offs = [.0,.0]
    #offs=[3.,0.]
    #ira = np.abs(ra-offs[0]) < size/2.
    #ide = np.abs(dec-offs[1]) < size/2.
    #data = data[:,:,ira]
    #data = data[:,ide,:]
    #npix = ira.sum()*ide.sum()
    data *= gain
    noise *= gain

    #If displaying a model, convert flux from Jy/beam to Jy/pixel
    if (hdr['object'])[:5] != 'model':
        try:
            bmaj,bmin = hdr['bmaj'],hdr['bmin']
        except KeyError:
            #multiple beams are present, and listed in a table
            #convert from arceconds to degrees
            bmaj,bmin = im[1].data[0][0]/3600.,im[1].data[0][1]/3600.
        beam = np.pi*bmaj*bmin/(4*np.log(2))
        pix = np.abs(hdr['cdelt1'])*np.abs(hdr['cdelt2'])
        #data *= (pix/beam)
        #noise *= (pix/beam)
        sigmaj = 3600*bmaj/(2*np.sqrt(2*np.log(2)))
        sigmin = 3600*bmin/(2*np.sqrt(2*np.log(2)))
        ram,dem = np.meshgrid(ra,dec)
        area = np.exp(-(ram**2/(2*sigmaj**2)+dem**2/(2*sigmin**2))).sum()
    else:
        area = 1.

    if (hdr['ctype3'] == 'VELO-LSR') or (hdr['ctype3']=='VELO-OBS') or (hdr['ctype3']=='VRAD'):
        xaxis = ((np.arange(hdr['naxis3'])+1-hdr['crpix3'])*hdr['cdelt3']+hdr['crval3'])/1e3
        #xaxis -= np.median(xaxis)
    else:
        nfreq = hdr['naxis3']
        freq = (np.arange(hdr['naxis3'])+1-hdr['crpix3'])*hdr['cdelt3']+hdr['crval3']
        try:
            xaxis = (hdr['restfrq']-freq)/hdr['restfrq']*2.99e5
        except KeyError:
            xaxis = (np.median(freq)-freq)/np.median(freq)*2.99e5

    #im_dec = data.sum(axis=2)
    #spec = im_dec.sum(axis=1)/area
    spec = np.zeros(hdr['naxis3'])
    if threshold is None:
        threshold = 3*noise
    print('threshold: ',threshold)
    for i in range(hdr['naxis3']):
        bright = data[i,:,:]>threshold
        spec[i] = data[i,:,:][bright].sum()/area

    if norm_peak:
        spec = spec/spec.max()
    if mirror:
        xaxis = xaxis[::-1]-.25#-.1

    plt.rc('axes',lw=2)
    if skip2:
        #skip the last two channels
        plt.plot(xaxis[:-2],spec[:-2],lw=4,**kwargs)
    else:
        plt.plot(xaxis,spec,lw=4,**kwargs)
    ax = plt.gca() #apply_aspect,set_adjustable,set_aspect,get_adjustable,get_aspect
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(14)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(14)
        tick.label1.set_fontweight('bold')
    plt.xlabel('Velocity (km/sec)',fontweight='bold',fontsize=14)
    if norm_peak:
        plt.ylabel('Normalized Flux',fontweight='bold',fontsize=14)
    else:
        plt.ylabel('Flux (Jy)',fontweight='bold',fontsize=14)
    if show_cloud:
        #Highlight the region affected by foreground absorption
        plt.fill_between([4,6],[plt.ylim()[0],plt.ylim()[0]],[plt.ylim()[1],plt.ylim()[1]],alpha=.2,color='k')
    im.close()

    #print(xaxis)
    #print(xaxis[:49],xaxis[98:])
    #print(xaxis[:14],xaxis[27:])

    if fwhm:
        #Calculate the fwhm of the line
        m=spec.max()
        dv1 = np.interp(m/2,(spec[xaxis<0].squeeze())[::-1],(xaxis[xaxis<0].squeeze())[::-1])
        dv2 = np.interp(m/2,(spec[xaxis>0].squeeze()),(xaxis[xaxis>0].squeeze()))
        return np.abs(dv1)+np.abs(dv2)

    print('dv,v_centroid:',xaxis[1]-xaxis[0],(spec*xaxis).sum()/spec.sum())
    return spec.sum()*(xaxis[1]-xaxis[0])
    #return xaxis,spec

def get_cropped_image(file,imx,offs=[0.,0.]):
    '''Given a fits file name and a size, return a cropped image.'''

    alma = fits.open(file)
    image = alma[0].data.squeeze()
    hdr = alma[0].header
    ra = 3600*hdr['cdelt1']*(np.arange(hdr['naxis1'])-hdr['naxis1']/2.-0.5)

    noise = calc_noise(image,np.round(imx/(3600.*hdr['cdelt1'])))
    de = -1*ra
    ira = np.abs(ra-offs[0]) < imx/2.
    ide = np.abs(de-offs[1]) < imx/2.
    cm_tmp = image[:,:,ira]
    cm_tmp = cm_tmp[:,ide,:]
    ra = ra[ira]
    de = de[ide]
    #alma.close()
    return alma,cm_tmp,hdr,ra,de,noise

def plot_temp_av2(model='cofid'):
    '''For the disk structure derived by various models, plot the radial profiles of N2H+ and DCO+, as well as the temperature and surface density structure. This can help show how the distribution of N2H+ and DCO+ changes with different factors.
    alma.n2h_cofid_I0.fits: model='cofid'
    alma.n2h_warmout_I0.fits: model='warmout'
    alma.n2h_midflat_I0.fits: model='midflat'
    alma.n2h_smallwarmout_I0.fits: model='smallwarmout'
    alma.n2h.photod.fits: model='photod'
    '''

    if model == 'cofid':
        q = [-.371,-.371,0] #fiducial model, photod
        n2hfile = 'cofid/alma.n2h_cofid.model_M0.fits'
        dcofile = 'cofid/alma.dco_cofid.model_M0.fits'
        label = 'CO-derived'
        Tmid=14.3
    if model == 'warmout':
        q = [-.371,1.7,2] #fiducial model, photod
        n2hfile = 'alma.n2h_warmout.model_M0.fits'
        dcofile = 'DCO_code/alma.dco_warmout.model_M0.fits'
        label = '70% warmer outer disk'
        Tmid=14.3
    if model =='smallwarmout':
        q = [-.371,1.3,2] #fiducial model, photod
        n2hfile = 'smallwarmout/alma.n2h_smallwarmout.model_M0.fits'
        dcofile = 'smallwarmout/alma.dco_smallwarmout.model_M0.fits'
        label = '30% warmer outer disk'
        Tmid=14.3
    if model =='warmout':
        q = [-.371,1.5,2] #fiducial model, photod
        n2hfile = 'warmout/alma.n2h_warmout.model_M0.fits'
        dcofile = 'warmout/alma.dco_warmout.model_M0.fits'
        label = '50% warmer outer disk'
        Tmid=14.3
    if model =='photodwarmout':
        q = [-.371,1.5,2] #fiducial model, photod
        n2hfile = 'alma.n2h.warmout.photodlow_M0.fits'
        dcofile = 'DCO_code/alma.dco.warmout.photod_M0.fits'
        label = 'Model'
        Tmid=14.3
    if model =='photodwarmoutmid0p15':
        q = [-0.15,-.371,4,1.5] #fiducial model, photod
        n2hfile = 'alma.n2h.warmout.mid0p15_M0.fits'
        dcofile = 'DCO_code/alma.dco.warmout.mid0p15.photod_M0.fits'
        label = '50% warmer outer disk, Photodesorption'
        Tmid=14.3
    if model =='photodwarmoutmid0p5':
        q = [-0.5,-.371,4,1.5] #fiducial model, photod
        n2hfile = 'alma.n2h.warmout.mid0p5_M0.fits'
        dcofile = 'DCO_code/alma.dco.warmout.mid0p5.photod_M0.fits'
        label = '50% warmer outer disk, Photodesorption'
        Tmid=14.3
    if model == 'photod':
        q = [-.371,-.371,0] #fiducial model, photod
        n2hfile = 'alma.n2h.photodlow_M0.fits'
        dcofile = 'DCO_code/alma.dco.photod_M0.fits'
        label = 'Photodesorbtion'
        Tmid=14.3
    if model == 'lowphotod':
        q = [-.371,-.371,0] #fiducial model, photod
        n2hfile = 'lowphotod/alma.n2h_lowphotod.model_M0.fits'
        dcofile = 'lowphotod/alma.dco_photod.model_M0.fits'
        label = 'Photodesorbtion'
        Tmid=14.3
    if model == 'highphotod':
        q = [-.371,-.371,0] #fiducial model, photod
        n2hfile = 'highphotod/alma.n2h_highphotod.model_M0.fits'
        dcofile = 'lowphotod/alma.dco_photod.model_M0.fits'
        label = 'Photodesorbtion'
        Tmid=14.3
    if model == 'midflat':
        q = [0,-.371,3] #fiducial model, photod
        n2hfile = 'midflat/alma.n2h_midflat.model_M0.fits'
        dcofile = 'midflat/alma.dco_midflat.model_M0.fits'
        label = 'q$_{mid}$=0'
        Tmid=14.3
    if model =='smallwarmoutmidflat':
        q = [0.,-.371,4,1.3] #fiducial model, photod
        n2hfile = 'alma.n2h.smallwarmout.midflat_M0.fits'
        dcofile = 'DCO_code/alma.dco.smallwarmout.midflat_M0.fits'
        label = '30% warmer outer disk'
        Tmid=14.3
    if model =='smallwarmoutmid0p15':
        q = [-0.15,-.371,4,1.3] #fiducial model, photod
        n2hfile = 'alma.n2h.smallwarmout.mid0p15_M0.fits'
        dcofile = 'DCO_code/alma.dco.smallwarmout.mid0p15_M0.fits'
        label = '30% warmer outer disk'
        Tmid=14.3
    if model =='midcold':
        q = [-.371,-0.371,0,1.3] #fiducial model, photod
        n2hfile = 'alma.n2h.midcold_M0.fits'
        dcofile = 'DCO_code/alma.dco.midcold_M0.fits'
        label = 'T$_{mid}$=10 K'
        Tmid = 10.
    if model =='midcoldsmallwarmout':
        q = [-.371,-0.371,0,1.3] #fiducial model, photod
        n2hfile = 'alma.n2h.midcold.smallwarmout_M0.fits'
        dcofile = 'DCO_code/alma.dco.midcold.smallwarmout_M0.fits'
        label = 'T$_{mid}$=10 K'
        Tmid = 10.


    #q = [-.371,-.371,0] #fiducial model, photod
    #q = [-.371,1.3,2] #smallwarmout
    #q = [-.371,1.7,2] #warmout
    #q = [0,-.371,3] #midflat

    params = [q, #qq
              0.04,#10**(p[1]), #Mdisk
              1.,#1.181526,#p[6],#1.,   #pp
              1.,#10.,  #Rin
              1000.,#Rout
              10**(2.444),#10**(p[2]), #Rc
              -36.0, #51.5 #inclination, which will be negative for IM Lup.
              .54,  #Mstar
              2*10**(-5.),#1e-4, #Xco
              .279*3.438,#p[2],#10**(p[3]), #vturb
              70.,#2*3.4707*np.sqrt(p[4]),#****** CHANGE BACK WHEN DONE ***** 70.,#p[4], #Zq0
              Tmid,#19.,  #Tmid0
              24.68, #Tatm0
              19., #Tco
              [.79,1000], #upper and lower boundaries in column density
              [9.,800.], #inner and outer boundaries for abundance
              -1]   #handed

    d=disk.Disk(params,rtg=False)
    obs = [150,101,300,170] #150,101,280,170 rt grid nr,nphi,nz,zmax

    d.set_obs(obs)
    d.set_rt_grid()
    d.set_line('CO21')
    #d.add_mol_ring(params[14][0],params[14][1],.79,3.,params[8],just_frozen=True)

    #co21_hpcc3
    #CO emitting region
    r = np.arange(10,500,20)
    zup_co21 = np.array([3.1,8.7,15.5,23.5,31.4,41.6,51.7,61.9,72.1,82.3,93.6,103.8,113.9,123.0,130.9,136.6,139.9,143.4,145.6,147.9,149.0,150.2,151.3,151.3,151.3])
    zlow_co21 = np.array([-1.4,-3.7,-7.1,-10.5,-12.7,-15.0,-17.3,-20.7,-22.9,-26.3,-28.6,-30.8,-33.1,-35.4,-38.8,-41.0,-44.4,-48.9,-51.2,-55.7,-58.0,-60.3,-62.5,-63.7,-64.8])

    #Plot N2H+ and DCO+ radial profiles
    cube_n2h = imagecube('alma.n2hdata_I0.fits')
    #cube_dco = imagecube('radialprofilesandmomentmaps/dcolb.contsub_clean.pbcor_I0.fits')
    cube_dco = imagecube('DCO_code/dco_M0.fits')
    xn2h,yn2h,dyn2h = cube_n2h.radial_profile(inc=36.0,PA=154.8)
    xdco,ydco,dydco = cube_dco.radial_profile(inc=36.0,PA=154.8)

    cube_n2hmodel = imagecube('radialprofilesandmomentmaps/radprofile_models/'+n2hfile)
    xn2hmodel,yn2hmodel,dyn2hmodel = cube_n2hmodel.radial_profile(inc=36.0,PA=154.8)

    cube_dcomodel = imagecube('radialprofilesandmomentmaps/radprofile_models/'+dcofile)
    xdcomodel,ydcomodel,ddcomodel = cube_dcomodel.radial_profile(inc=36.0,PA=154.8)

    #Pull in CO and continuum data
    cube_co = imagecube('/Volumes/Untitled2/data/DMTau/dmtau_co21sblb_new_I0.fits')
    xco,yco,dco = cube_co.radial_profile(inc=36.0,PA=154.8)
    cube_cont = imagecube('/Volumes/Untitled2/data/DMTau/dmtau_cont.cm.fits')
    xcont,ycont,dcont = cube_cont.radial_profile(inc=36.0,PA=154.8)

    #Read in the beam size
    xgauss = np.linspace(400,500,50)
    #xgauss = np.linspace(0.7,1.3,50) # for continuum
    sigmaj = cube_n2h.beam[0]/(2*np.sqrt(2*np.log(2)))*158 #convert from FWHM to sigma
    #ygauss = 1*np.exp(-(xgauss-1)**2./(2*sigmaj**2.))+2 #for continuum
    ygauss = .15*np.exp(-(xgauss-450)**2./(2*sigmaj**2.))+.5



    plt.figure()

    plt.subplot(221)
    plt.errorbar(xn2h*144.5,yn2h/yn2h.max(),dyn2h/yn2h.max(),fmt='.',capsize=1.25,capthick=1.25,color='k')
    plt.step(xn2h*144.5,yn2h/yn2h.max(),where='mid',color='k')
    plt.step(xn2hmodel*144.5,yn2hmodel/yn2hmodel.max(),where='mid',color='b')
    plt.plot(xgauss,ygauss,color='k',alpha=.8)
    plt.axhline(0.,color='k',alpha=.2)
    plt.text(420,.45,'beam',alpha=.8)
    plt.legend(('Data',label),frameon=False)
    plt.title('N$_2$H$^+$')
    plt.xlim(0,500)
    plt.ylabel('Normalized Flux')

    plt.subplot(222)
    plt.errorbar(xdco*144.5,ydco/ydco.max(),dydco/ydco.max(),fmt='.',capsize=1.25,capthick=1.25,color='k')
    plt.step(xdco*144.5,ydco/ydco.max(),where='mid',color='k')
    plt.step(xdcomodel*144.5,ydcomodel/ydcomodel.max(),where='mid',color='b')
    #plt.step(xco*144.5,yco/yco.max(),where='mid',color='r',alpha=.8,ls='--')
    #plt.step(xcont*144.5,ycont/ycont.max(),where='mid',color='g',alpha=.8,ls=':')
    plt.axhline(0.,color='k',alpha=.2)
    #plt.legend(('Data',label,'CO','Dust'),frameon=False)
    plt.legend(('Data',label),frameon=False)
    plt.title('DCO$^+$')
    plt.xlim(0,500)
    plt.ylabel('Normalized Flux')

    plt.subplot(223)
    cs = plt.contourf(d.r[0,:,:]/d.AU,d.Z[0,:,:]/d.AU,d.T[0,:,:],(9,19),cmap=plt.cm.Set1)
    cs2 = plt.contour(d.r[0,:,:]/d.AU,d.Z[0,:,:]/d.AU,d.T[0,:,:],np.linspace(9,19,11),colors='k',alpha=.2)
    if model == 'lowphotod':
        #self.sig_col*Disk.Hnuctog/Disk.m0>Sig0*Disk.sc
        cs = plt.contour(d.r[0,:,:]/d.AU,d.Z[0,:,:]/d.AU,np.abs(d.sig_col[0,:,:]*d.Hnuctog/d.m0/d.sc),(3,100.),linestyles=':',colors='k')
    if model == 'highphotod':
        #self.sig_col*Disk.Hnuctog/Disk.m0>Sig0*Disk.sc
        cs = plt.contour(d.r[0,:,:]/d.AU,d.Z[0,:,:]/d.AU,np.abs(d.sig_col[0,:,:]*d.Hnuctog/d.m0/d.sc),(.01,.79),linestyles=':',colors='k')
    #plt.fill_between(r,-zlow_co21,zup_co21,color='k',alpha=.4)
    plt.text(120,60,'CO snowline',rotation='vertical',fontsize=8)
    if model =='cofid':
        plt.text(520,5,'N$_2$\n snowline',rotation='vertical',fontsize=8)
    if model =='midcold':
        plt.text(290,5,'N$_2$\n snowline',rotation='vertical',fontsize=8)
    #plt.text(300,140,'N$_2$H$^+$')
    plt.xlim(0,500)
    plt.xlabel('R (au)')
    plt.ylabel('Z (au)')
    plt.ylim(0,170)

    plt.subplot(224)
    cs = plt.contourf(d.r[0,:,:]/d.AU,d.Z[0,:,:]/d.AU,d.T[0,:,:],(19,30),cmap=plt.cm.Set1)
    cs2 = plt.contour(d.r[0,:,:]/d.AU,d.Z[0,:,:]/d.AU,d.T[0,:,:],np.linspace(19,30,12),colors='k',alpha=.2)
    if model == 'lowphotod' or model=='highphotod':
        #cs = plt.contour(d.r[0,:,:]/d.AU,d.Z[0,:,:]/d.AU,np.abs(d.sig_col[0,:,:]),(.01,.1),linestyles=':',colors='r')
        cs = plt.contour(d.r[0,:,:]/d.AU,d.Z[0,:,:]/d.AU,np.abs(d.sig_col[0,:,:]*d.Hnuctog/d.m0/d.sc),(.79,3.),linestyles=':',colors='k')
    #plt.fill_between(r,-zlow_co21,zup_co21,color='k',alpha=.4)
    plt.xlim(0,500)
    plt.text(160,70,'CO snowline',rotation='vertical',fontsize=8)
    plt.xlabel('R (au)')
    plt.ylabel('Z (au)')
    plt.ylim(0,170)


def plot_chan_n2models(model='cofid',include_fluxscale=False):
    '''For the disk structure derived by various models, plot channel maps for the data and the models to provide a comparison.
    alma.n2h_cofid_I0.fits: model='cofid'
    alma.n2h_warmout_I0.fits: model='warmout'
    alma.n2h_midflat_I0.fits: model='midflat'
    alma.n2h_smallwarmout_I0.fits: model='smallwarmout'
    alma.n2h.photod.fits: model='photod'
    '''

    if model == 'cofid':
        q = [-.371,-.371,0] #fiducial model, photod
        n2hfile = 'cofid/alma.n2h_cofid.model.fits'
        dcofile = 'cofid/alma.dco_cofid.model.fits'
        label = 'CO-derived'
        Tmid=14.3
    if model =='smallwarmout':
        q = [-.371,1.3,2] #fiducial model, photod
        n2hfile = 'smallwarmout/alma.n2h_smallwarmout.model.fits'
        dcofile = 'smallwarmout/alma.dco_smallwarmout.model.fits'
        label = '30% warmer outer disk'
        Tmid=14.3
    if model =='warmout':
        q = [-.371,1.5,2] #fiducial model, photod
        n2hfile = 'warmout/alma.n2h_warmout.model.fits'
        dcofile = 'warmout/alma.dco_warmout.model.fits'
        label = '50% warmer outer disk'
        Tmid=14.3
    if model =='photodwarmout':
        q = [-.371,1.5,2] #fiducial model, photod
        n2hfile = 'alma.n2h.warmout.photodlow.fits'
        dcofile = 'DCO_code/alma.dco.warmout.photod.fits'
        label = 'Model'
        Tmid=14.3
    if model =='photodwarmoutmid0p15':
        q = [-0.15,-.371,4,1.5] #fiducial model, photod
        n2hfile = 'alma.n2h.warmout.mid0p15.fits'
        dcofile = 'DCO_code/alma.dco.warmout.mid0p15.photod.fits'
        label = '50% warmer outer disk, Photodesorption'
        Tmid=14.3
    if model =='photodwarmoutmid0p5':
        q = [-0.5,-.371,4,1.5] #fiducial model, photod
        n2hfile = 'alma.n2h.warmout.mid0p5.fits'
        dcofile = 'DCO_code/alma.dco.warmout.mid0p5.photod.fits'
        label = '50% warmer outer disk, Photodesorption'
        Tmid=14.3
    if model == 'photod':
        q = [-.371,-.371,0] #fiducial model, photod
        n2hfile = 'alma.n2h.photodlow.fits'
        dcofile = 'DCO_code/alma.dco.photod.fits'
        label = 'Photodesorbtion'
        Tmid=14.3
    if model == 'lowphotod':
        q = [-.371,-.371,0] #fiducial model, photod
        n2hfile = 'lowphotod/alma.n2h_lowphotod.model.fits'
        dcofile = 'lowphotod/alma.dco_photod.model.fits'
        label = 'Photodesorbtion'
        Tmid=14.3
    if model == 'highphotod':
        q = [-.371,-.371,0] #fiducial model, photod
        n2hfile = 'highphotod/alma.n2h_highphotod.model.fits'
        dcofile = 'lowphotod/alma.dco_photod.model.fits'
        label = 'Photodesorbtion'
        Tmid=14.3
    if model == 'midflat':
        q = [0,-.371,3] #fiducial model, photod
        n2hfile = 'midflat/alma.n2h_midflat.model.fits'
        dcofile = 'midflat/alma.dco_midflat.model.fits'
        label = 'q$_{mid}$=0'
        Tmid=14.3
    if model =='smallwarmoutmidflat':
        q = [0.,-.371,4,1.3] #fiducial model, photod
        n2hfile = 'alma.n2h.smallwarmout.midflat.fits'
        dcofile = 'DCO_code/alma.dco.smallwarmout.midflat.fits'
        label = '30% warmer outer disk'
        Tmid=14.3
    if model =='smallwarmoutmid0p15':
        q = [-0.15,-.371,4,1.3] #fiducial model, photod
        n2hfile = 'alma.n2h.smallwarmout.mid0p15.fits'
        dcofile = 'DCO_code/alma.dco.smallwarmout.mid0p15.fits'
        label = '30% warmer outer disk'
        Tmid=14.3
    if model =='midcold':
        q = [-.371,-0.371,0,1.3] #fiducial model, photod
        n2hfile = 'alma.n2h.midcold.fits'
        dcofile = 'DCO_code/alma.dco.midcold.fits'
        label = 'T$_{mid}$=10 K'
        Tmid = 10.
    if model =='midcoldsmallwarmout':
        q = [-.371,-0.371,0,1.3] #fiducial model, photod
        n2hfile = 'alma.n2h.midcold.smallwarmout.fits'
        dcofile = 'DCO_code/alma.dco.midcold.smallwarmout.fits'
        label = 'T$_{mid}$=10 K'
        Tmid = 10.


    data = fits.open('alma.n2hdata.fits')
    model = fits.open('radialprofilesandmomentmaps/radprofile_models/'+n2hfile)
    mod_hdr,data_hdr = model[0].header,data[0].header
    model_im,data_im = model[0].data.squeeze(),data[0].data.squeeze()
    data2 = fits.open('DCO_code/dco.fits')
    model2 = fits.open('radialprofilesandmomentmaps/radprofile_models/'+dcofile)
    mod2_hdr,data2_hdr = model2[0].header,data2[0].header
    model2_im,data2_im = model2[0].data.squeeze(),data2[0].data.squeeze()
    ra = 3600*mod_hdr['cdelt1']*(np.arange(mod_hdr['naxis1'])-mod_hdr['naxis1']/2.-0.5)
    ra2 = 3600*mod2_hdr['cdelt1']*(np.arange(mod2_hdr['naxis1'])-mod2_hdr['naxis1']/2.-0.5)

    # Crop image
    imx = 8.
    noise = calc_noise(data_im,np.round(imx/(3600.*data_hdr['cdelt1'])))
    de = -1*ra
    ira = np.abs(ra) < imx/2.
    ide = np.abs(de) < imx/2.
    data_im = data_im[:,:,ira]
    data_im = data_im[:,ide,:]
    model_im = model_im[:,:,ira]
    model_im = model_im[:,ide,:]
    ra = ra[ira]
    de = de[ide]

    noise2 = calc_noise(data2_im,np.round(imx/(3600.*data2_hdr['cdelt1'])))
    de2 = -1*ra2
    ira = np.abs(ra2) < imx/2.
    ide = np.abs(de2) < imx/2.
    data2_im = data2_im[:,:,ira]
    data2_im = data2_im[:,ide,:]
    model2_im = model2_im[:,:,ira]
    model2_im = model2_im[:,ide,:]
    ra2 = ra2[ira]
    de2 = de2[ide]

    try:
        bmaj,bmin,bpa = data_hdr['bmaj'],data_hdr['bmin'],data_hdr['bpa']
    except KeyError:
        bmaj,bmin,bpa = data[1].data[0][0]/3600.,data[1].data[0][1]/3600.,data[1].data[0][1]
    #Beam size
    a,b,phi = bmaj/2*3600.,bmin/2*3600.,np.radians(bpa)
    t=np.linspace(0,2*np.pi,100)
    x = -(imx/2.-1.)+a*np.cos(t)*np.cos(phi)-b*np.sin(t)*np.sin(phi)
    y = -(imx/2.-1.)+a*np.cos(t)*np.sin(phi)+b*np.sin(t)*np.cos(phi)

    try:
        bmaj2,bmin2,bpa2 = data2_hdr['bmaj'],data2_hdr['bmin'],data2_hdr['bpa']
    except KeyError:
        bmaj2,bmin2,bpa2 = data2[1].data[0][0]/3600.,data2[1].data[0][1]/3600.,data2[1].data[0][1]
    #Beam size
    a,b,phi = bmaj2/2*3600.,bmin2/2*3600.,np.radians(bpa2)
    t=np.linspace(0,2*np.pi,100)
    x2 = -(imx/2.-1.)+a*np.cos(t)*np.cos(phi)-b*np.sin(t)*np.sin(phi)
    y2 = -(imx/2.-1.)+a*np.cos(t)*np.sin(phi)+b*np.sin(t)*np.cos(phi)

    if (data_hdr['ctype3'] == 'VELO-LSR') or (data_hdr['ctype3'] == 'VELO-OBS'):
        velo = ((np.arange(hdr['naxis3'])+1-data_hdr['crpix3'])*data_hdr['cdelt3']+data_hdr['crval3'])/1e3
        #velo -= 5.79#np.median(velo)
    else:
        nfreq = data_hdr['naxis3']
        freq = (np.arange(data_hdr['naxis3'])+1-data_hdr['crpix3'])*data_hdr['cdelt3']+data_hdr['crval3']
        try:
            velo = (data_hdr['restfrq']-freq)/data_hdr['restfrq']*2.99e5
        except KeyError:
            velo = (np.median(freq)-freq)/np.median(freq)*2.99e5

    levels = (np.arange(20))/20.*(np.max(data_im)-np.min(data_im))+np.min(data_im)
    print('Min level, noise:',levels.min(),noise)

    levels2 = (np.arange(20))/20.*(np.max(data2_im)-np.min(data2_im))+np.min(data2_im)
    print('Min level, noise:',levels2.min(),noise2)

    plt.figure()
    ax=plt.subplot(221)
    ax.set(aspect='equal',adjustable='box')
    cs = plt.contourf(ra,de,data_im[9,:,:],40,cmap=plt.cm.afmhot)
    plt.contour(ra,de,data_im[9,:,:],[3*noise,],colors='w',linewidths=2)
    plt.gca().invert_xaxis()
    plt.text(3.5,2.5,'N$_2$H$^+$',color='w',fontsize=16)
    #plt.xlabel(r'$\Delta\alpha$ (")',fontsize=18)
    #plt.ylabel('$\Delta\delta$ (")',fontsize=18)
    plt.fill(x-.5,y-.5,lw=2,color='w')
    #ax = plt.gca() #apply_aspect,set_adjustable,set_aspect,get_adjustable,get_aspect
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(14)
    if include_fluxscale:
        fig.subplots_adjust(top=.8)
        cbar_ax = fig.add_axes([.125,.805,.775,.03])
        cb=fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',format='%0.3f')
        cb.ax.xaxis.set_ticks_position('top')
        cb.ax.xaxis.set_label_position('top')
        #cb.ax.invert_xaxis()
        cb.set_label(label='Jy/beam',size=16)#,weight='bold')
    plt.plot((-(imx-1.5)/2.+1,-(imx-1.5)/2.+1-100./145.1),(imx/2.-1,imx/2.-1),lw=3,color='w')


    ax=plt.subplot(223)
    ax.set(aspect='equal',adjustable='box')
    cs = plt.contourf(ra,de,model_im[9,:,:],40,cmap=plt.cm.afmhot)
    plt.contour(ra,de,model_im[9,:,:],[3*noise,],colors='w',linewidths=2)
    plt.gca().invert_xaxis()
    plt.xlabel(r'$\Delta\alpha$ (")',fontsize=18)
    plt.ylabel('$\Delta\delta$ (")',fontsize=18)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(14)

    ax=plt.subplot(222)
    ax.set(aspect='equal',adjustable='box')
    cs = plt.contourf(ra2,de2,data2_im[11,:,:],40,cmap=plt.cm.afmhot)
    plt.contour(ra2,de2,data2_im[11,:,:],[3*noise2,],colors='w',linewidths=2)
    plt.gca().invert_xaxis()
    plt.text(3.5,2.5,'DCO$^+$',color='w',fontsize=16)
    #plt.xlabel(r'$\Delta\alpha$ (")',fontsize=18)
    #plt.ylabel('$\Delta\delta$ (")',fontsize=18)
    plt.fill(x2-.5,y2-.5,lw=2,color='w')
    if include_fluxscale:
        fig.subplots_adjust(top=.8)
        cbar_ax = fig.add_axes([.125,.805,.775,.03])
        cb=fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',format='%0.3f')
        cb.ax.xaxis.set_ticks_position('top')
        cb.ax.xaxis.set_label_position('top')
        #cb.ax.invert_xaxis()
        cb.set_label(label='Jy/beam',size=16)#,weight='bold')
    plt.plot((-(imx-1.5)/2.+1,-(imx-1.5)/2.+1-100./145.1),(imx/2.-1,imx/2.-1),lw=3,color='w')
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(14)


    ax=plt.subplot(224)
    ax.set(aspect='equal',adjustable='box')
    cs = plt.contourf(ra2,de2,model2_im[11,:,:],40,cmap=plt.cm.afmhot)
    plt.contour(ra2,de2,model2_im[11,:,:],[3*noise2,],colors='w',linewidths=2)
    plt.gca().invert_xaxis()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(14)

def run_simplex(model='cofid'):
    '''Run a downhill simplex fitting of data to a model. This function is used for N2H+.'''
    import disk
    from single_model import *


    if model == 'cofid':
        q=[-.371,-.371,0] #fiducial model, photod
        initial_guess = [-11.1,] #log(abundance of ring)
        #With tol=.1, min=-11.308125, chi^2=0.99385442

    if model =='smallwarmout':
        q=[-.371,1.3,2] #fiducial model, photod
        initial_guess = [-11.1,] #log(abundance of ring)
        #With tol=0.1, min=-11.3775, chi2=0.993766365679
        #  Significantly better than cofid, with Delta_AIC=Delta_BIC = 1400 (prob=8.9e-305)

    if model == 'midflat':
        q=[0,-.371,3] #fiducial model, photod
        initial_guess = [-11.1,] #log(abundance of ring)
        #With tol=0.1, min=-11.308125, chi2=0.99376027955
        #  Significantly better than smallwarmout, with Delta_AIC=Delta_BIC=96.8, prob=9e-22

    if model == 'highphotod':
        q=[-.371,-.371,0] #fiducial model, photod

        ### To fit only the abundance of the two rings, uncomment the next line, comment out the other initial_guess, and replace x[2] in the call to add_mol_ring with 325.
        #initial_guess = [-11.1,-10.0] #log(abundance of ring)
        #With tol=0.1, min=-11.36483093, -9.56600952, chi2=0.993715526489
        #  Significantly better than midflat, with Delta_BIC=696, Delta_AIC=711 (prob=2.9e-155)
        initial_guess = [-11.1,-10.0,325] #log(abundance of ring)
        #With tol=0.1, min=-11.39164, -9.4834387, 277.89, chi2=0.99369215
        #  Significantly better than 2 param highphotod, with Delta_AIC=371 (prob=1.9e-81), Delta_BIC=356

    #models.append('lowphotod')
    if model == 'lowphotod':
        q=[-.371,-.371,0] #fiducial model, photod
        initial_guess = [-11.1,-10.0,325] #log(abundance of ring)
        #with tol=0.1, min=-11.38986, -10.814, 314.52, chi2 = 0.993631575
        #  Significantly better than 3 param high photod with Delta_AIC=Delta_BIC=963.2(prob=6.8e-210)
        #  Significantly better than 2 param low photod, with Delta_AIC=155.9 (prob=1e-34), Delta_BIC=140.3

        ### To fit only the abundance of the two rings, uncomment the next line, comment out the other initial_guess, and replace x[2] in the call to add_mol_ring with 325.
        #initial_guess = [-11.1,-10.0]
        #with tol=0.1, min=-11.3775, -10.875, chi2=0.9936413787060615
        #  Significantly better than 2 param highphotod, with Delta_AIC=1179 (prob=9.3e-257), Delta_BIC=1179
        #  Significantly better than 3 param highphotod, with Delta_AIC=807 (prob=4.9e-176), Delta_BIC=822

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
              0.04,#10**(p[1]), #Mdisk
            1.,#1.181526,#p[6],#1.,   #pp
            1.,#10.,  #Rin
            1000.,#Rout
            10**(2.444),#10**(p[2]), #Rc
            -36.0, #51.5 #inclination, which will be negative for IM Lup.
            .54,  #Mstar
            10**(x[0]), #2*10**(-5.),#1e-4, #Xco
            .279*3.438,#p[2],#10**(p[3]), #vturb
            70.,#2*3.4707*np.sqrt(p[4]),#****** CHANGE BACK WHEN DONE ***** 70.,#p[4], #Zq0
            14.3,#19.,  #Tmid0
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
            d.add_mol_ring(x[2],x[2]+50,3,100,10**(x[1]),just_frozen=False) #3,100 #.01, .79
        if model == 'highphotod':
            d.add_mol_ring(x[2],x[2]+50,.01,.79,10**(x[1]),just_frozen=False) #3,100 #.01, .79
        total_model(disk=d,chanmin=chanmin,nchans=nchans,chanstep=chanstep,offs=[xoff,yoff],modfile='alma',imres=resolution,obsv=obsv,vsys=vsys,freq0=279.51170100,Jnum=2,distance=144.5,hanning=True,PA=154.8,bin=4)

        # - Generate model visibilities
        model_vis = vis_sample(uvfile=datfile+'.vis.fits',imagefile='alma.fits',mod_interp=False)
        real_model = model_vis.real
        imag_model = model_vis.imag
        chi = ((real_model-real_obj)**2*weight_real).sum() + ((imag_model-imag_obj)**2*weight_imag).sum()
        return chi/nu


    print('starting minimization')

    ### Minimize the reduced chi-squared, to within a tolerance of 0.1
    result = minimize(objective_function,initial_guess,method='nelder-mead',tol=.1,args=(model))
    if result.success:
        print('Minimum found at ',result.x)
        print('Minimum value: ',result.fun)
    else:
        print('Optimization failed: ',result.message)

def run_simplex_dco(model='cofid'):
    '''Run a downhill simplex fitting of data to a model for DCO+.'''
    import disk_dco
    from single_model_dco import *

    if model == 'cofid':
        #models.append('cofid')
        q=[-.371,-.371,0] #fiducial model, photod
        initial_guess = [-11.1,] #log(abundance of ring)
        #With tol=.1, min=-11.655, chi^2=0.9186662777

    if model =='smallwarmout':
        #models.append('smallwarmout')
        q=[-.371,1.3,2] #fiducial model, photod
        initial_guess = [-11.1,] #log(abundance of ring)
        #With tol=0.1, min=-11.5856, chi2=0.9186099438 -> Significantly better than cofid, with Delta_AIC=Delta_BIC=997.1, prob=3e-217

    if model == 'midflat':
        #models.append('midflat')
        q=[0,-.371,3] #fiducial model, photod
        initial_guess = [-11.1,] #log(abundance of ring)
        #With tol=0.1, min=-11.030625, chi2=0.91839062838
        #  Significantly better than cofid, with Delta_AIC=Delta_BIC=4879, prob=0.
        #  Significantly better than cofid, with Delta_AIC=Delta_BIC=3881, prob=0

    if model == 'photod':
        #models.append('highphotod')
        q=[-.371,-.371,0] #fiducial model, photod

        ### To fit only the abundance of the two rings, uncomment the next line, comment out the other initial_guess, and replace x[2] in the call to add_mol_ring with 325.
        #initial_guess = [-11.1,-10.0] #log(abundance of ring)
        #With tol=0.1, min=-11.63061, -9.604248 chi2=0.918313619
        #  Significantly better than midflat, with Delta_AIC=1363 (prob=1e-216), Delata_BIC=1347
        initial_guess = [-11.1,-10.0,325] #log(abundance of ring)
        #With tol=0.1, min=-11.698575, -9.47846, 230.693 chi2=0.9180635
        #  Significantly better than two param photod, with Delta_AIC=4427 (prob=0.), Delta_BIC=4411


    #disks = []

    #for q in qs:
    datfile = 'alma.dcodata'
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
    resolution = 0.03
    #obs = [180,131,300,170] #150,101,280,170 rt grid nr,nphi,nz,zmax

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
            1.,#1.181526,#p[6],#1.,   #pp
            1.,#10.,  #Rin
            1000.,#Rout
            10**(2.444),#10**(p[2]), #Rc
            -36.0, #51.5 #inclination, which will be negative for IM Lup.
            .54,  #Mstar
            [10**(-20.),10**(x[0])], #2*10**(-5.),#1e-4, #Xco
            .279*3.438,#p[2],#10**(p[3]), #vturb
            70.,#2*3.4707*np.sqrt(p[4]),#****** CHANGE BACK WHEN DONE ***** 70.,#p[4], #Zq0
            14.3,#19.,  #Tmid0
            24.68, #Tatm0
            #19., #Tco
            [.79,1000], #upper and lower boundaries in column density
            [9.,800.], #inner and outer boundaries for abundance
            -1, 200]   #handed

        d=disk_dco.Disk(params,rtg=False)
        obs = [150,101,300,170] #150,101,280,170 rt grid nr,nphi,nz,zmax

        d.set_obs(obs)
        d.set_rt_grid()
        d.set_line('dco')

        if model =='photod':
            d.add_mol_ring(x[2],x[2]+50,.79,3.,10**(x[1]),just_frozen=False) #3,100 #.01, .79
        total_model(disk=d,chanmin=chanmin,nchans=nchans,chanstep=chanstep,offs=[xoff,yoff],modfile='alma',imres=resolution,obsv=obsv,vsys=vsys,freq0=288.143858,Jnum=3,distance=144.5,hanning=True,PA=154.8,bin=2)

        # - Generate model visibilities
        model_vis = vis_sample(uvfile=datfile+'.vis.fits',imagefile='alma.fits',mod_interp=False)
        real_model = model_vis.real
        imag_model = model_vis.imag
        chi = ((real_model-real_obj)**2*weight_real).sum() + ((imag_model-imag_obj)**2*weight_imag).sum()
        return chi/nu


    print('starting minimization')

    ### Minimize the reduced chi-squared, to within a tolerance of 0.1
    result = minimize(objective_function,initial_guess,method='nelder-mead',tol=.1,args=(model))
    if result.success:
        print('Minimum found at ',result.x)
        print('Minimum value: ',result.fun)
    else:
        print('Optimization failed: ',result.message)
