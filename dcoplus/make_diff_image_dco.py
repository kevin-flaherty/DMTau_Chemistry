#Import visibility fits file into CASA format.
os.system('rm -rf alma.diff.ms')
importuvfits('alma.diff.vis.fits','alma.diff.ms')

#Clean the model visibilities
os.system('rm -rf diff_image.*')
tclean(vis='alma.diff.ms',imagename='diff_image',weighting='briggs', robust = 1.5,interactive=False, cell='0.03arcsec',specmode='cubedata',imsize=[600,600],  mask = 'box [ [170pix , 173pix],[403pix , 411pix] ]', threshold ='0.002mJy',niter=5000,outframe='LSRK')

#Output the cleaned image into fits format
os.system('rm -rf alma.diff.fits')
exportfits('diff_image.image','alma.diff.fits')

#Clean up
os.system('rm -rf alma.diff.ms')
os.system('rm -rf diff_image.*')
