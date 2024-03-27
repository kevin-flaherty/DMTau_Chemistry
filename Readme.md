# DM Tau chemistry code

## Contents
This repository contains the code used by Diop et al. to model and analyze ALMA observations of the N2H+(4-3) and DCO+(3-2) emission from around DM Tau obtained as part of 2018.1.01119.S (PI: K. Flaherty).

In particular the code models the emission from N2H+ and DCO+ under various conditions (e.g. a warm outer midplane, CO photo-desorption). The underlying disk model is based on a parametric structure described in more detail in Diop et al., as well as in Flaherty et al. 2020 (which in turn is based on earlier work by Dartois et al. 2003 and Rosenfeld et al. 2013).

This repository contains three folders. The folders *n2hplus* and *dcoplus* contain the code used to model N2H+ and DCO+ respectively. The folder models contains documents describing the exact function calls used to generate the models discussed in the paper.

**n2hplus**:
- *disk.py*: Code for calculating density, temperature, and velocity at r,z throughout the disk, based on the input parametric structure.
- *raytrace.py*: Takes as input a disk object created by disk.py, and performs radiative transfer to create an image at a specified velocity/spatial resolution, distance, position angle, etc.
- *single_model.py*: Wrapper for raytrace.py and disk.py. Given a particular set of parameters, functions within single_model_dco.py will create a model image and compare the visibilities from this image to the data.
- *make_model_image.py*: CASA commands for creating cleaned model channel maps from model visibilities.
- *make_diff_image.py*: CASA commands for created cleaned residual channel maps from the difference in visibilities between the model and data.
- *n2h.dat*: Molecular constants
- *alma.n2hdata.vis.fits*: Visibility fits file for the N2H+ data.

**dcoplus**:
- *disk.py*: Code for calculating density, temperature, and velocity at r,z throughout the disk, based on the input parametric structure.
- *raytrace.py*: Takes as input a disk object created by disk.py, and performs radiative transfer to create an image at a specified velocity/spatial resolution, distance, position angle, etc.
- *single_model_dco.py*: Wrapper for raytrace.py and disk.py. Given a particular set of parameters, functions within single_model_dco.py will create a model image and compare the visibilities from this image to the data.
- *make_model_image_dco.py*: CASA commands for creating cleaned model channel maps from model visibilities.
- *make_diff_image_dco.py*: CASA commands for created cleaned residual channel maps from the difference in visibilities between the model and data.
- *dco.dat*: Molecular constants
- *alma.dcodata.vis.fits*: Visibility fits file for the DCO+ data.

**Models**:
- *cofid.txt*: Code used to generate the CO-fiducial model, i.e. the model based the structure derived from CO, with no further modifications.
- *midflat.txt*: Code used to generate the model with a midplane radial temperature profile that is shallower than at the surface layers.
- *smallwarmout.txt*: Code used to generate models in which the midplane temperature beyond 200 au is increased by 30%.
- *warmout.txt*: Code used to generate models in which the midplane temperature beyond 200 au is increased by 50%.
- *highphotod.txt*: Code used to generate models in which an extra ring of DCO+ and N2H+ is placed in the outer disk. The extra N2H+ is placed high in the disk, mimicking the effect of CO photo-dissociation.
- *lowphotod.txt*: Code used to generate models in which an extra ring of DCO+ and N2H+ is placed in the outer disk. The extra N2H+ is placed close to the midplane, mimicking the effect of additional ionization from X-rays or cosmic rays.
- *photod_warmout.txt*: Same as above, but with an increase in the midplane temperature beyond 200 au.

## Usage

The most direct way to generate models is using the *lnlike* function within *single_model.py or *single_model_dco.py*. This function takes as an input a subset of model parameters, and calls functiosn within *disk.py* and *raytrace.py* to create a model image at the same spatial and spectral resolution as the data.

For N2H+ an example usage of the *lnlike* function is:

```
lnlike(((-.371,-.371,0),-11.1,325,-100,0.,0.,0.),cleanup=False)
```

The first element is a list of model parameters:

```
[q,abund,Rc,abund2,turbulence, x-offset, y-offset, Rbreak]
```
The abundance is the logarithm of the abundance relative to H2, Rc is in units of au, turbulence is in units of the thermal broadening of N2H+, and the x and y offsets are in units of arc-seconds. Rbreak refers to the radii beyond which e.g. the temperature jumps (see discussion of q).

The start of the *lnlike* function includes the full list of model parameters beyond these inputs. `q` can take multiple forms, but is always a list. In this case the first element of the list is the `q` value for the midplane, the second is the `q` value for the disk atmosphere, and the third is a flag (which can take values of 0, 1, 2, 3, or 4) which in this case specifies that the midplane and the atmosphere have the same temperature profile. The other values of this flag specify:
- `q=1`: a double power-law shape for the temperature profile (i.e. T~r^(q1) inside of Rbreak and T~r^(q2) outside of Rbreak)
- `q=2`: a jump in the temperature beyond Rbreak. In this case the second element of this list specifies the multiplicative factor applied to the temperature beyond Rbreak. For example (-.25, 1.3,2) specifies that T~r^(-.25), with a 30% increase in the temperature beyond Rbreak.
- `q=3`: separate temperature profiles for the midplane and atmosphere, with the first element specifying the temperature profile of the midplane and the second element specifying the temperature profile of the atmosphere.
- `q=4`: separate temperature profiles for the midplane and atmosphere, along with a jump in the temperature beyond Rbreak. In this case the list must include four elements, i.e. (-.15,-.5,4.1.3), where the midplane follows T~r^(-.15), the atmosphere follows T~r^(-0.5), and there is a 30% jump in temperature beyond Rbreak. 


## Attribution
If you make use of this package in your research, please cite Diop et al. in prep.
