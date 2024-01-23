#Import visibility fits file into CASA format.
os.system('rm -rf alma.model.ms')
importuvfits('alma.model.vis.fits','alma.model.ms')

#Clean the model visibilities
os.system('rm -rf model_image.*')
tclean(vis='alma.model.ms',imagename='model_image',outframe='LSRK',weighting='briggs', robust = 1.5,cell='0.03arcsec',specmode='cubedata',imsize=[600,600], mask = 'box [ [170pix , 173pix],[403pix , 411pix] ]', interactive = False, threshold = '0.002mJy',niter=5000)
#imlup_co21_line_image.mask
#Output the cleaned image into fits format
os.system('rm -rf alma.model.fits')
exportfits('model_image.image','alma.model.fits')

#Clean up
os.system('rm -rf alma.model.ms')
#os.system('rm -rf model_image.*')
