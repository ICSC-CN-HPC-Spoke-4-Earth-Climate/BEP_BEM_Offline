BEP-BEM Offline
-------------

These files form part of the ICSC-CN-HPC-Spoke-4-Earth-Climate project. 

This work presents a novel approach employing the Building Effect Parameterization
(BEP, Martilli et al., 2002) coupled with a Building Energy Model (BEM, Salamanca et al., 2009) 
in an offline configuration to simulate urban climates.

These scripts are provided to run the BEP+BEM Offline model and input/output file generation 
for multiple urban sites.

Contents
--------
These are fortran-90 scripts. 

The following example scripts will work for any urban site:

- main_program.f90: atmospheric column model that reads site geometric data and meteorological forcings 
- module_sf_bep_bem.f90: building effect parameterisation scheme    
- module_sf_bem.f90: simple building energy model 
 
The following folders are an example test site (Melbourne, Australia) with annual data.

- **/input_file/Input_2d/d_ **: text files with site characteristic information
- **/input_file/Input_2d/w_ **: meteorological forcing data in .txt format (1 year data)

Instructions
------------
Meteorological forcings
You can download data from https://cds.climate.copernicus.eu/ to force the model, or from 
observation data or other climate models. 
The details of the forcings variables are listed below, and we recommend taking them from the 
following ERA5 datasets: 

ERA5 hourly data on pressure levels from 1940 to present (1000hPa)
- 1- Specific humidity (kg kg-1)
- 2- Temperature (K)
- 3- U-component of wind (m s-1)
- 4- V-component of wind (m s-1)

ERA5 hourly data on single levels from 1940 to present
- 5- Surface direct short-wave radiation flux  (W m-2)
- 6- Surface downward long-wave radiation flux  (W m-2)
- 7- Surface downward short-wave radiation flux (W m-2)
- 8- Surface pressure (Pa)
- 9- Total precipitation (m)

Urban Canopy Parameters (UCPs)
If available, the following urban morphological parameters can be used to appropriately 
describe the geometry of the urban fabric. If not available, tabulated values may be used. 

- λu - Fraction of urban land use (1)
- Ps - Pervious surface fraction (1)
- λp - Roof area fraction (1)
- Hm - Building mean height (m)
- λb - Wall to plan area ratio (m)
- hi – Distribution of building height (m) 

To produce a simulation with BEP+BEM Offline, run the code: 

main_program.f90

The outputs of the main atmospheric variables and urban surfaces variables are saved 
within the folder named Output

More detailed instructions are contained in each script.

These files are distributed freely in the hope that it will be useful, but with NO WARRANTY, 
implied or otherwise. Users are responsible for ensuring they are fit for their purpose. 

Reference
------------
Martilli, A., Clappier, A. and Rotach, M. W. (2002), 'An urban surface exchange parameterisation
for mesoscale models', Boundary-Layer Meteorology 104(2), 261(2)304.
URL: https://doi.org/10.1023/A:1016099921195

Salamanca, F., Krpo, A., Martilli, A. and Clappier, A. (2009), 'A new building energy model
coupled with an urban canopy parameterization for urban climate simulations-part I. formulation,
verification, and sensitivity analysis of the model', Theoretical and Applied Climatology
99(3), 331.
URL: https://doi.org/10.1007/s00704-009-0142-9
