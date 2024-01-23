#Import visibility fits file into CASA format.
os.system('rm -rf alma.diff.ms')
importuvfits('alma.diff.vis.fits','alma.diff.ms')

#Clean the model visibilities
os.system('rm -rf diff_image.*')
tclean(vis='alma.diff.ms',imagename='diff_image',weighting='briggs', robust = 1.5,interactive=False, cell='0.0325arcsec',specmode='cubedata',imsize=[576,576],  mask = 'box [ [155pix , 177pix],[422pix , 396pix] ]', threshold ='0.0mJy',niter=50000,outframe='LSRK',pbcor=True)

#Output the cleaned image into fits format
os.system('rm -rf alma.diff.fits')
exportfits('diff_image.image','alma.diff.fits')

#Clean up
os.system('rm -rf alma.diff.ms')
os.system('rm -rf diff_image.*')
