MLUCM BEP-BEM Offline
-------------

These files form part of the ICSC-CN-HPC-Spoke-4-Earth-Climate project. 

MLUCM BEP+BEM, incorporates the vertical turbulent diffusion scheme of Santiago and Martilli (2010) 
with the Building Effect Parameterization (BEP, Martilli et al., 2002) and the 
Building Energy Model (BEM, Salamanca et al., 2009). 

MLUCM BEP+BEM bridges mesoscale and microscale phenomena within the urban canopy, capturing scale interactions and feedbacks. 
The accuracy and low computational cost of this 1-dimensional model makes it ideal for offline climate projections to assess 
urban climate impacts under different emission scenarios. The model's features allow analyzing urban overheating, 
energy demands, and evaluating the efficiency of strategies like green roofs, cool roofs, and photovoltaic panels

These scripts are provided to run the MLUCM BEP+BEM Offline model and input/output file generation 
for case study at Preston, Melborune (Australia), reported at Pappaccogli et al., [preprint]

Contents
--------
These are fortran-90 scripts. 

The following example scripts will work for any urban site:

- main_program.f90: atmospheric column model that reads site geometric data and meteorological forcings 
- module_sf_bep_bem.f90: building effect parameterisation scheme    
- module_sf_bem.f90: simple building energy model 
 
The following folders are an example test site (Melbourne, Australia) with annual data.

- **/input_file/Input_2d/d_ **: text files with site characteristic information
- **/input_file/Input_ERA5/w_ **: meteorological forcing data in .txt format (case study at Preston, Melbourne - Australia)

Instructions
------------
Meteorological forcings - ERA5 
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

- Î»u - Fraction of urban land use (1)
- Ps - Pervious surface fraction (1)
- Î»p - Roof area fraction (1)
- Hm - Building mean height (m)
- Î»b - Wall to plan area ratio (m)
- hi â€“ Distribution of building height (m) 

To produce a simulation with BEP+BEM Offline, run the code: 

main_program.f90

The outputs of the main atmospheric variables and urban surfaces variables are saved 
within the folder 

More detailed instructions are contained in each script.

Authors
-----------
Gianluca Pappaccogli, University of Salento
Dipartimento di Scienze e Tecnologie Biologiche ed Ambientali, University of Salento, Lecce, 73100, Italy
Correspondence to: Gianluca Pappaccogli (gianluca.pappaccogli@unisalento.it)

Reference
------------
Pappaccogli, G., Zonato, A., Martilli, A., Buccolieri, R. and Lionello, P. 
'MLUCM BEP+BEM: An offline one-dimensional Multi-Layer Urban Canopy Model based on the BEP+BEM Scheme',
Geoscientific Model Development, EGUsphere [preprint]

Pappaccogli, G., Zonato, A., Martilli, A., Buccolieri, R., and Lionello, P.: The Implementation of the 
BEP+BEM Offline Parameterization Scheme: Exploring Urban Dynamics through Climatic Projections, 
EGU General Assembly 2024, Vienna, Austria, 14â€“19 Apr 2024, EGU24-6308, https://doi.org/10.5194/egusphere-egu24-6308, 2024.

Martilli, A., Clappier, A. and Rotach, M. W. (2002), 'An urban surface exchange parameterisation
for mesoscale models', Boundary-Layer Meteorology 104(2), 261(2)304.
URL: https://doi.org/10.1023/A:1016099921195

Salamanca, F., Krpo, A., Martilli, A. and Clappier, A. (2009), 'A new building energy model
coupled with an urban canopy parameterization for urban climate simulations-part I. formulation,
verification, and sensitivity analysis of the model', Theoretical and Applied Climatology
99(3), 331.
URL: https://doi.org/10.1007/s00704-009-0142-9

Santiago, J. L., & Martilli, A. (2010). A Dynamic Urban Canopy Parameterization for Mesoscale Models 
Based on Computational Fluid Dynamics Reynolds-Averaged Navier-Stokes Simulations. Boundary-Layer Meteorology, 
137(2), 231-245.
