# ADIC
Adaptive Digital Image Correlation for biological materials in medicine and material science

# Description
Computes deformation fields out of image natural texture by performing digital image correlation between test and reference images. 

The ADIC_deformation algorithm is a texture correlation (TC) algorithm, which is optimized and expanded to robustly extract
deformation information from neutron radiographies while ignoring "shadow regions", where only poor correlation statistics are available.
The inner TC core performs a zero-order search (rigid subset) with a zero-normalized cross-correlation calculated at integer positions and
refined with bicubic interpolation. The center position and size of the search window are adaptively adjusted. Error control ensure correlation,
uniqueness and continuity of deformation estimates.

The ADIC_strain algorithm computes strain fields by fitting discontinuous deformation fields to
continuous cubic smoothing splines and differentiating along the spatial coordinates
Strains are defined as: eps_xx = dux/dx; eps_yy = duy/dy; eps_xy = 0.5*(dux/dy + duy/dx).

## Getting started
Run EXAMPLE1.m, EXAMPLE2.m, and EXAMPLE3.m in /examples to get started
Each of the examples are explained in detail in the literature reference below.

Input data is included in /data

Results are stored in /results

Source code functions are available in /src

## About the examples

EXAMPLE1.m: CORRELATION OF OPTICAL IMAGES OF SPECKLED SAMPLES

Moisture-induced swelling gradients in growth rings at hygroscopic equilibrium
are determined from optical images of speckled softwood samples
Reference image: RH = 0%, Test image: RH = 95%


EXAMPLE2.m : CORRELATION OF NEUTRON RADIOGRAPHIES OF SOFTWOOD GROWTH RINGS

Moisture-induced swelling gradients in growth rings at hygroscopic equilibrium
are determined from neutron radiographies
Reference image: RH = 0%, Test image: RH = 95%


EXAMPLE3.m: CORRELATION OF NEUTRON RADIOGRAPHIES OF WOOD FIBER COMPOSITES

One-dimensional moisture diffusion in three-layer wood fiber composites
is determined from neutron radiographies
Reference image: t = 0 h (dry state), Test image: t= 17h


# Authors and acknowledgment
The author list for this code repository is:

Sanabria SJ^{a,*},

Lanvermann C^{b},

Franco M^{a},

Mannes D^{a},

Niemz P^{a},


The participating institutions are: 

{a} Institute for Building Materials ETH Zurich, Stefano-Franscini-Platz 6, CH-8093, Zurich, Switzerland

{b} Neutron Imaging and Activation Group, Paul Scherrer Institute, CH-5232 PSI, Villigen, Switzerland

*Corresponding author: E-mail: sanse@stanford.edu

# License
[MIT](https://choosealicense.com/licenses/mit/)

Please cite the original work in:
1. Sanabria Martin, Sergio & Lanvermann, Christian & Michel, Franco & Mannes, David & Niemz, Peter. (2014). Adaptive Neutron Radiography Correlation for Simultaneous Imaging of Moisture Transport and Deformation in Hygroscopic Materials. Experimental Mechanics. 55. 403-415. 10.1007/s11340-014-9955-2. 
