PROGRAM main_program

!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +    by A. Martilli,     CIEMAT  SP 2040 MADRID                  +
!     +                        phone: ++34-9-13-46-62-99               +
!     +                        email:alberto.martilli@ciemat.es        +
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +  Modified by G. Pappaccogli, University of Salento             +
!     +  Modified by A. Zonato, KNMI 				       +
!     +				              Mar 2023 - May 2024      +
!     +                    email:gianluca.pappaccogli@unisalento.it    +
!     +                    email:Andrea.zonato@knmi.nl                 +
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  USE module_sf_bep_bem
  IMPLICIT NONE
! ----------------------------------------------------------------------- 
!  Dimension for the array used in the BEP+BEM module
! -----------------------------------------------------------------------  
integer nurbmax_mp ! Maximum number of urban classes
      parameter (nurbmax_mp=3)
 integer ndm_mp ! Maximum number of street directions
      parameter (ndm_mp=2)
  integer nz_um_mp ! Maximum number of vertical levels in the urban grid 
     parameter(nz_um_mp=41)
  integer ng_u_mp ! Number of grid levels in the ground
      parameter (ng_u_mp=10)
  integer ngr_u_mp ! Number of grid levels in green roof
      parameter (ngr_u_mp=10)
  integer nwr_u_mp ! Number of grid levels in the walls or roofs (layer esterno 10)
      parameter (nwr_u_mp=10)
  integer nf_u_mp !Number of grid levels in the floors (BEM)
      parameter (nf_u_mp=10)
  integer ngb_u_mp !Number of grid levels in the ground below building (BEM)
      parameter (ngb_u_mp=10)
  real dz_u_mp ! Urban grid resolution
      parameter (dz_u_mp=1.)
  integer nbui_max_mp !maximum number of types of buildings in an urban class
      parameter (nbui_max_mp=9)   !must be less or equal than nz_um
  integer history_interval
      parameter(history_interval=1)

  real, dimension(nbui_max_mp) :: hi_urb2d = [0.0,0.0,0.494,0.377,0.095,0.017,0.012,0.005,0.000] 
!==============================================================================================================
! Define the number of rows to be read from the weather model files/ERA5 - 1-minute time resolution
  integer, parameter :: n_loop = ((9 * 365.0 + 2*366.0 + 332) * 24.0 * 60.0)+13*60
    
  real :: dt=60					!Temporal resolution (seconds)
  real :: gmt=10				!Greenwich Mean Time simulation start
  integer :: itimestep=1
  integer :: julday=1				!julian day start simulation
  integer :: doy=1	                	!day of the year start simulation
  real, parameter :: DeltaGMT=10                !difference with Greenwich Meridian Time [h]

  INTEGER, PARAMETER :: ix=1, iy=1		!indices for printing
  INTEGER, PARAMETER :: ims = 1, ime = 1	!number of elements along x-axis
  INTEGER, PARAMETER :: kms = 1, kme = nz_um_mp	!number of elements along z-axis
  INTEGER, PARAMETER :: jms = 1, jme = 1	!number of elements along y-axis
  INTEGER, PARAMETER :: start_year=1993, final_year=2004
  INTEGER :: count_year(final_year-start_year+1) 
 !===================================================================================
  real, parameter :: deg2rad=3.14159265358979323846 / 180.0
  
 ! Define the number of files to read from weather model/ERA5
  integer, parameter :: num_files = 10

  ! Define arrays to store the data
  real :: input_vector(num_files, n_loop)
  real :: vector(n_loop)
  real :: z_bou(kms:kme+1),dz_bou(kms:kme+1)
  real :: dld(kms:kme),dlu(kms:kme),th0(kms:kme) 
  real :: mean_dl
  
  ! Define file-related variables
  character(256) :: file_path_vector 
  character(256) :: file_list_vector(num_files)
  character(256) :: filename_vector
  integer :: i_wm=1
  real :: declin_urb=0.0

!------------------------------------------------------------------------------------------------
  INTEGER :: i, k, j, z, i_o, file_num, file_num_2d, file_num_vector, iii, status
  INTEGER :: ibui, id, iz_um, iz_u, iw, ig, iz, ir, kkk

! arrays used to collapse indexes  
  integer :: ind_zwd_mp(nbui_max_mp,nz_um_mp,nwr_u_mp,ndm_mp)
  integer :: ind_gd_mp(ng_u_mp,ndm_mp)
  integer :: ind_zd_mp(nbui_max_mp,nz_um_mp,ndm_mp)
  integer :: ind_zdf_mp(nz_um_mp,ndm_mp)
  integer :: ind_zrd_mp(nz_um_mp,nwr_u_mp,ndm_mp)
  integer :: ind_grd_mp(nz_um_mp,ngr_u_mp,ndm_mp)
  integer :: ind_bd_mp(nbui_max_mp,nz_um_mp)
  integer :: ind_wd_mp(nbui_max_mp,nz_um_mp,ndm_mp)
  integer :: ind_gbd_mp(nbui_max_mp,ngb_u_mp,ndm_mp)
  integer :: ind_fbd_mp(nbui_max_mp,nf_u_mp,nz_um_mp-1,ndm_mp)

  INTEGER, dimension(ims:ime,jms:jme) :: UTYPE_URB2D
!WEATHER MODEL INPUT
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme )::   u_phy 	!%"U_PHY"   "x-wind component at mass point"         "m s-1"
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme )::   v_phy 	!%"V_PHY"   "y-wind component at mass point"         "m s-1"
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme )::   th_phy	!% th_phy ?? air temperature
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme )::   p_phy        !% p_phy  ??
  
  real a_u(ims:ime,kms:kme,jms:jme)         ! Implicit component for the momemtum in X-direction (center)
  real a_v(ims:ime,kms:kme,jms:jme)         ! Implicit component for the momemtum in Y-direction (center)
  real a_t(ims:ime,kms:kme,jms:jme)         ! Implicit component for the temperature
  real a_e(ims:ime,kms:kme,jms:jme)         ! Implicit component for the TKE
  real b_u(ims:ime,kms:kme,jms:jme)         ! Explicit component for the momemtum in X-direction (center)
  real b_v(ims:ime,kms:kme,jms:jme)         ! Explicit component for the momemtum in Y-direction (center)
  real b_t(ims:ime,kms:kme,jms:jme)         ! Explicit component for the temperature
  real b_e(ims:ime,kms:kme,jms:jme)         ! Explicit component for the TKE
  real b_q(ims:ime,kms:kme,jms:jme)         ! Explicit component for the Humidity
  real dlg(ims:ime,kms:kme,jms:jme)         ! Height above ground (L_ground in formula (24) of the BLM paper).
  real dl_u(ims:ime,kms:kme,jms:jme)        ! Length scale (lb in formula (22) ofthe BLM paper).
  real sf(ims:ime,kms:kme,jms:jme)           ! surface of the urban grid cells
  real vl(ims:ime,kms:kme,jms:jme)             ! volume of the urban grid cells

  real, dimension(ims:ime, jms:jme) :: arr_2d				! 2D array to be read from file
  real, dimension(27,ims:ime, jms:jme) :: input_2d			! 2D arrat to be assign for each input file
  character(len=256) :: file_path, file_path_2d, filename, filename_2d  ! file path
  character(len=256) :: path_output, filename_output			! path output
  character(len=10), dimension(47) :: file_list		 		! list of filenames 3D
  character(len=10), dimension(28) :: file_list_2d			! list of filenames 2D
  character(len=256) :: file_path_wm

  real :: value, value_0, value_1, value_1b, value_2, value_3, value_4, value_4b, value_5, value_6, value_7, value_99

!===============================================================================================================
  INTEGER, PARAMETER :: ids=1, ide=1, jds=1, jde=1, kds=1, kde=1, its=1, ite=1, jts=1, jte=1, kts=1, kte=kme, &
  	urban_map_zrd=nz_um_mp*nwr_u_mp*ndm_mp,  urban_map_zwd=nbui_max_mp*nz_um_mp*nwr_u_mp*ndm_mp, &
  	urban_map_gd=ng_u_mp*ndm_mp, urban_map_zd=nbui_max_mp*nz_um_mp*ndm_mp, urban_map_zdf=nz_um_mp*ndm_mp, &
  	urban_map_bd=nbui_max_mp*nz_um_mp, urban_map_wd=nbui_max_mp*nz_um_mp*ndm_mp, urban_map_gbd=nbui_max_mp*ngb_u_mp*ndm_mp, &
  	urban_map_fbd=nbui_max_mp*nf_u_mp*nz_um_mp*ndm_mp, urban_map_zgrd=nz_um_mp*nwr_u_mp*ndm_mp, &
  	num_urban_hi=nbui_max_mp, num_urban_ndm=2, num_urban_ncomf=10

!======= variabili solar position ======================================
  real :: Lat_rad, Delta_TSL, h_S_mean, theta_Z, t, sinh_S, h_S 

  REAL, DIMENSION( ims:ime,jms:jme )::  FRC_URB2D
  REAL, DIMENSION( ims:ime,jms:jme )::  swdown
  REAL, DIMENSION( ims:ime,jms:jme )::  GLW
  REAL, DIMENSION(ims:ime,jms:jme ) :: XLAT, XLONG
  REAL, DIMENSION(ims:ime,jms:jme ) :: COSZ_URB2D
  REAL, DIMENSION(ims:ime,jms:jme ) :: OMG_URB2D
  REAL, DIMENSION(ims:ime,jms:jme) ::  lp_urb2d
  REAL, DIMENSION(ims:ime,jms:jme) ::  lb_urb2d
  REAL, DIMENSION(ims:ime,jms:jme) ::  hgt_urb2d
  REAL, DIMENSION(ims:ime,jms:jme ) :: sf_ac_urb3d
  REAL, DIMENSION(ims:ime,jms:jme ) :: lf_ac_urb3d
  REAL, DIMENSION(ims:ime,jms:jme ) :: cm_ac_urb3d
  REAL, DIMENSION(ims:ime,jms:jme ) :: sfvent_urb3d
  REAL, DIMENSION(ims:ime,jms:jme ) :: lfvent_urb3d
  REAL, DIMENSION(ims:ime,jms:jme ) :: ep_pv_urb3d
  REAL, DIMENSION(ims:ime,jms:jme ) :: qgr_urb3d
  REAL, DIMENSION(ims:ime,jms:jme ) :: tgr_urb3d
  REAL, DIMENSION(ims:ime,jms:jme ) :: draingr_urb3d
  REAL, DIMENSION(ims:ime,jms:jme ) :: rainbl
  REAL, DIMENSION(ims:ime,jms:jme ) :: swddir
  REAL, DIMENSION(ims:ime,jms:jme ) :: swddif
  REAL  PBLH
  REAL tsk_av 
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme )::   DZ8W
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme )::   RHO
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme )::   qv_phy
 
  real rl_up(its:ite,jts:jte) ! upward long wave radiation
  real rs_abs(its:ite,jts:jte) ! absorbed short wave radiation
  real emiss(its:ite,jts:jte)  ! emissivity averaged for urban surfaces
  real grdflx_urb(its:ite,jts:jte)  ! ground heat flux for urban areas
  real rs_up(its:ite,jts:jte)       ! upward radiation based on swdown-rs_up
  real H_plumber, LE_plumber

  REAL, DIMENSION(ims:ime, 1:urban_map_zrd, jms:jme) :: trb_urb4d	!%"TRB_URB4D" "ROOF LAYER TEMPERATURE"          "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_zwd, jms:jme) :: tw1_urb4d	!"TW1_URB4D" "WALL LAYER TEMPERATURE"          "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_zwd, jms:jme) :: tw2_urb4d	!"TW2_URB4D" "WALL LAYER TEMPERATURE"          "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_gd , jms:jme) :: tgb_urb4d	!"TGB_URB4D" "ROAD LAYER TEMPERATURE"          "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_zgrd, jms:jme) :: trv_urb4d	!"TRV_URB4D" "GREEN ROOF LAYER TEMPERATURE"          "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_zgrd, jms:jme) :: qr_urb4d	!"QR_URB4D" "GREEN ROOF LAYER MOISTURE"          "dimensionless"
  
  REAL, DIMENSION(ims:ime, 1:urban_map_zdf, jms:jme) :: drain_urb4d	!"DRAIN_URB4D" "GREEN ROOF DRAINAGE"          "mm"
  REAL, DIMENSION(ims:ime, 1:urban_map_bd, jms:jme ) :: tlev_urb3d	!"TFLEV_URB3D" "FLOOR TEMPERATURE"                       "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_bd , jms:jme) :: qlev_urb3d	!"QLEV_URB3D" "SPECIFIC HUMIDITY"              "dimensionless"
  REAL, DIMENSION(ims:ime, 1:urban_map_wd , jms:jme) :: tw1lev_urb3d	!"TW1LEV_URB3D" "WINDOW TEMPERATURE"           "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_wd , jms:jme) :: tw2lev_urb3d	!"TW2LEV_URB3D" "WINDOW TEMPERATURE"           "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_gbd, jms:jme) :: tglev_urb3d	!"TGLEV_URB3D" "GROUND TEMPERATURE BELOW A BUILDING"     "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_fbd, jms:jme) :: tflev_urb3d	!"TFLEV_URB3D" "FLOOR TEMPERATURE"                       "K"

  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme) :: t_pv_urb3d	!"T_PV_URB3D"  "PHOTOVOLTAIC PANELS TEMPERATURE " "K" !PVP
  REAL, DIMENSION( ims:ime, 1:urban_map_wd , jms:jme) :: sfwin1_urb3d	!"SFWIN1_URB3D" "SENSIBLE HEAT FLUX FROM URBAN SFC WINDOW"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:urban_map_wd , jms:jme) :: sfwin2_urb3d	!"SFWIN2_URB3D" "SENSIBLE HEAT FLUX FROM URBAN SFC WINDOW"  "W m{-2}"

  REAL, DIMENSION( ims:ime, 1:urban_map_zd , jms:jme) :: sfw1_urb3d	!"SFW1_URB3D"  "SENSIBLE HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:urban_map_zd , jms:jme) :: sfw2_urb3d	!"SFW2_URB3D"  "SENSIBLE HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme) :: sfr_urb3d	!"SFR_URB3D"  "SENSIBLE HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:num_urban_ndm, jms:jme) :: sfg_urb3d	!"SFG_URB3D"  "SENSIBLE HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme) :: sfrv_urb3d	!"SFRV_URB3D"  "SENSIBLE HEAT FLUX FROM GREEN ROOF"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme) :: lfrv_urb3d
  REAL,  DIMENSION( ims:ime, 1:urban_map_gd, jms:jme) :: tgv_urb4d !GARDEN
  REAL, DIMENSION( ims:ime, 1:urban_map_gd, jms:jme)  :: qg_urb4d !GARDEN
  REAL, DIMENSION( ims:ime, 1:num_urban_ndm, jms:jme) :: sfgv_urb3d !GARDEN
  REAL, DIMENSION( ims:ime, 1:num_urban_ndm, jms:jme) :: lfgv_urb3d !GARDEN
  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme) :: dgr_urb3d !GRZ"	DGR_URB3D" "ROOF LAYER DEPTH WATER RETENTION"          "mm"
  REAL, DIMENSION( ims:ime, 1:num_urban_ndm, jms:jme) :: dg_urb3d !GRZ		"DG_URB4D" "ROOF LAYER DEPTH WATER RETENTION"          "mm"
  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme) :: lfr_urb3d !GRZ	"LFR_URB3D"  "LATENT HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:num_urban_ndm, jms:jme) :: lfg_urb3d !G	"LFG_URB3D"  "LATENT HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  
! WRF COMFORT
 REAL, DIMENSION( ims:ime, jms:jme )  :: tmr_11,tmr_12,tmr_13,tmr_21,tmr_22,tmr_23
 REAL, DIMENSION( ims:ime, jms:jme )  :: comf_10,comf_50,comf_90
 REAL, DIMENSION( ims:ime, 1:num_urban_ncomf,jms:jme ) :: hist_comf
! WRF COMFORT

!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +    by A. Martilli,     CIEMAT  SP 2040 MADRID                  +
!     +                        phone: ++34-9-13-46-62-99               +
!     +                        email:alberto.martilli@ciemat.es        +
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           integer nzm
           parameter (nzm=nz_um_mp) 
           integer jz
           real, dimension(1) :: lambda_p		!from .txt input files
           real, dimension(nzm+1) ::  hgt_column
           real, dimension(nzm+1) ::  prb_column
! flow variables
           real vx(nzm) ! x component of the wind
           real vy(nzm) ! y component of the wind
           real sstot(nzm+1)
           real pt(nzm) ! potential temperature 
           real qv(nzm) ! specific humidity
           real tke(nzm) ! tke
           real trc(nzm) ! tracer
           real dpdx,dpdy ! external pressure gradient in x and y direction
           real utau_x,utau_y ! utau 
! variable needed in the calculation of the TKE and vertical diffusion  
           real dlk(nzm) ! mixing length
           real dls(nzm) ! dissipation length scale
           real dlk_u(nzm) ! mixing length form the morphology
           real dls_u(nzm) ! dissipation length scale from the morphology
           real cdz(nzm+1) ! vertical diffusion coefficients for momentum
	   real cdt(nzm+1) ! vertical diffusion coefficients for scalars
           real sh(nzm) ! shear production term
           real bu(nzm) ! buoyancy production term !AM to be added
           real ceps,ck,cmu,g,temin ! constants for the k-l scheme
	   real pr ! prandtl number
           parameter (ceps=1/1.4,ck=0.4,cmu=0.09,g=9.81,temin=0.001)
! variables needed in the resolution of the diffusion equation: here source/sinks for a generic variable c are defined as c*srim_c+srex_c
           real srim_vx(nzm) ! implicit part of the source/sinks terms of v1
           real srim_vy(nzm) ! implicit part of the source/sinks terms of v2
           real srim_tke(nzm) ! implicit part of the source/sinks terms of tke
           real srim_trc(nzm) ! implicit part of the source/sinks terms of tracer
           real srim_pt(nzm) ! implicit part of the source/sinks terms of potential temperature
           real srex_vx(nzm) ! explicit part of the source/sinks terms of v1
           real srex_vy(nzm) ! explicit part of the source/sinks terms of v2
           real srex_tke(nzm) ! explicit part of the source/sinks terms of tke
           real srex_trc(nzm) ! explicit part of the source/sinks terms of tracer
           real srex_pt(nzm) ! explicit part of the source/sinks terms of potential temperature
           real srex_qv(nzm) !Explicit component for the Humidity   !_gl
           real srim_qv(nzm) 
           real vl_c(nzm) ! fraction of air in each cell
           real sf_c(nzm+1) ! fraction of air at the interface between cells
! dispersive flux
           real dwtrc(nzm+1)
	   real duw(nzm+1)
           real aaa(4,4),bbb(4)

! drag coefficient
          ! real cdragx(nzm),cdragy(nzm) ! drag coefficient in x and y direction           
! building information
           real ss(nzm+1) ! ss(iz)=probability to have a building of height z(iz)
           real pb(nzm+1) ! pb(iz)=probbaility to have a building taller or equal to z(iz)
           real wx,wy ! distance between buildings at street level in the x and y direction respectively
           real bx,by ! building dimension in the x and y direction respectively
           real hmean ! mean building height
           real hgt(nzm+1),prb(nzm+1)
! emission at surface
           real emi
! numerical parameters
           real lmo
           real zlmo(nzm+1),z_ground(nzm+1),prandtl(nzm+1)
           real, parameter :: dz=dz_u_mp ! vertical resolution
           integer nz
           parameter (nz=nz_um_mp) ! number of vertical levels
           real z_c(nzm+1) ! height of the faces of the numerical grid
           integer iz_c ! loop index for vertical
           real time ! time from the start
           real time_pr ! time from the last print

! output variables 
           real uw(nzm+1) ! turbuelnt flux of vx
           real uwm
           real duwdz(nzm) ! vertical derivative of the turbulent flux of vx
           real vw(nzm+1) ! turbuelnt flux of vy
           real dvwdz(nzm) ! vertical derivative of the turbulent flux of vy
           real wtke(nzm+1) ! turbuelnt flux of tke
           real dwtkedz(nzm) ! vertical derivative of the turbulent flux of tke
           real wtrc(nzm+1) ! turbuelnt flux of tracer
           real dwtrcdz(nzm) ! vertical derivative of the turbulent flux of tracer
           real wt(nzm+1) ! turbuelnt flux of pt
           real dwtdz(nzm) ! vertical derivative of the turbulent flux of pt
           real dwq(nzm+1) !_gl
           real dwqdz(nzm)  

! flag to distinguish between alligned and staggered arrays (1: staggered, 2: alligned)
           integer iconfig
! working
           real lambdaf_x,lambdaf_y
           real tke_old(nzm) ! tke
           real dtot,wrk(8)
           real trc_an(nzm),trc_2(nzm)
           integer i_col
! convert th to real temperature
           real rcp,cp,p0_rt,r_rt
           REAL :: pr_t(nz)
           REAL :: tr(nz)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
write(*,*) "n_loop", n_loop

    ! Initialize the vector with years from 1993 to 2024
    do i = start_year, final_year
        count_year(i-(start_year-1)) = i
    end do

    write(*,*) "count_year", count_year
    do i = kms, kme+1
     z_bou(i) = i*dz_u_mp
     dz_bou(i)= dz_u_mp
    end do

! define dlg
DO i=ims,ime
  DO j=kms,kme
    DO  k=jms,jme
        dlg(i, j, k) = z_bou(j)
    end do
  end do
end do

! ASSIGN A VALUE TO ALL ELEMENTS OF THE 3D ARRAY
DO i=ims,ime
  DO j=kms,kme
    DO  k=jms,jme
       DZ8W(i, j, k) = dz_u_mp
       RHO(i, j, k) = 1.293
       qv_phy(i, j ,k) = 0.011
    END DO
  END DO
END DO

value = 290.00;
value_0 = 0.0;
value_99 = 0.0;

DO i=ims,ime
  DO j=kms,kme
    DO  k=jms,jme
        a_u(i, j, k) = value_99
        a_v(i, j, k) = value_99
        a_t(i, j, k) = value_99
        a_e(i, j, k) = value_99
        b_u(i, j, k) = value_99
        b_v(i, j, k) = value_99
        b_t(i, j, k) = value_99
        b_e(i, j, k) = value_99
        b_q(i, j, k) = value_99
        dl_u(i, j, k) = value_99
        sf(i, j, k) = value_99
        vl(i, j, k) = value_99
       END DO
     END DO
  END DO

DO i = ims, ime
    DO j = 1, urban_map_zrd   
      DO k = jms,jme 
        trb_urb4d(i, j, k) = 293.0
        qr_urb4d(i, j, k) = value_0
        trv_urb4d(i, j, k) = value
      END DO
    END DO
  END DO
 
 value_1 = 290.00;
 value_1b = 0.0;
DO i = ims, ime
    DO j = 1, urban_map_zd    
      DO k = jms,jme
        tw1lev_urb3d(i, j, k) = value_1
        tw2lev_urb3d(i, j, k) = value_1
        sfwin1_urb3d(i, j, k) = value_1b
        sfwin2_urb3d(i, j, k) = value_1b
        sfw1_urb3d(i, j, k) = value_1b
        sfw2_urb3d(i, j, k) = value_1b
      END DO
    END DO
  END DO  

 value_2 = 290.00;

DO i = ims, ime
    DO j = 1, urban_map_zwd   
      DO k = jms,jme
        tw1_urb4d(i, j, k) = value_2
        tw2_urb4d(i, j, k) = value_2
      END DO
    END DO
  END DO

 value_3 = 290.00;

DO i = ims, ime
    DO j = 1, urban_map_gd    
      DO k = jms,jme
        tgb_urb4d(i, j, k) = value_3
        tgv_urb4d(i, j, k) = value_3 
        qg_urb4d(i, j, k)  = 0.8
     END DO
    END DO
  END DO

 value_4 = 290.00;
 value_4b = 0.0001 
DO i = ims, ime
    DO j = 1, urban_map_zdf   
      DO k = jms,jme
        drain_urb4d(i, j, k) = value_4b
        t_pv_urb3d(i, j, k) = value_4
        sfr_urb3d(i, j, k) = value_1b
        sfrv_urb3d(i, j, k) = value_1b
        lfrv_urb3d(i, j, k) = value_1b
        dgr_urb3d(i, j, k) = value_0
        lfr_urb3d(i, j, k) = value_1b
      END DO
    END DO
  END DO

 value_5 = 290.00;

DO i = ims, ime
    DO j = 1, urban_map_bd    
      DO k = jms,jme
        tlev_urb3d(i, j, k) = value_5	!"INDOOR TEMPERATURE"
        qlev_urb3d(i, j, k) = value_0	!"SPECIFIC HUMIDITY"
      END DO
    END DO
  END DO
 
  value_6 = 290.00;

DO i = ims, ime
    DO j = 1, urban_map_gbd     
      DO k = jms,jme
        tglev_urb3d(i, j, k) = value_6	!GROUND TEMPERATURE BELOW A BUILDING
      END DO
    END DO
  END DO

  value_7 = 290.00;

DO i = ims, ime
    DO j = 1, urban_map_fbd     
      DO k = jms,jme
        tflev_urb3d(i, j, k) = value_7	!"FLOOR TEMPERATURE"
      END DO
    END DO
  END DO

 

  DO i = ims, ime
    DO j = 1, num_urban_ndm     
      DO k = jms,jme
        sfg_urb3d(i, j, k) = 0.  	!SENSIBLE HEAT FLUX FROM URBAN SFC
        dg_urb3d(i, j, k) = 0.		!"ROOF LAYER DEPTH WATER RETENTION"          "mm"
        lfg_urb3d(i, j, k) = 0.		!LATENT HEAT FLUX FROM URBAN SFC"
        sfgv_urb3d(i, j, k) = 0.
        lfgv_urb3d(i, j, k) = 0. 
     END DO
    END DO
  END DO

! Assign a value of 1 to all elements of the vector
   ! Initialize the array with values of 1
    do i = ims, ime
        do j = jms, jme
            UTYPE_URB2D(i, j) = 1
        end do
    end do
    
tmr_11=0.
tmr_12=0.
tmr_13=0.
tmr_21=0.
tmr_22=0.
tmr_23=0.
comf_10=0.
comf_50=0.
comf_90=0.
hist_comf=0.

! Define the common path to the input files
  file_path_2d = "input_file/Input_2d/"
  file_path_vector = "input_file/Input_ERA5/"
  path_output = "Output/"

  ! Define the list of filenames
  file_list_2d = (/ "d_01", "d_02", "d_03", "d_04", "d_05", "d_06", "d_07", "d_08", "d_09", &
                    "d_10", "d_11", "d_12", "d_13", "d_14", "d_15", "d_16", "d_17", &
                    "d_18", "d_19", "d_20", "d_21", "d_22", "d_23", "d_24", "d_25", &
                    "d_26", "d_27", "d_28" /)
 
 ! Define the list for ERA5
  file_list_vector = (/ "w_01", "w_02", "w_03", "w_04", "w_05", "w_06", "w_07", "w_08", &
                        "w_09", "w_10"/)

!-------------------------------------------------------------------------------------------
  ! Read the 2D array from each file with the same dimension
!-------------------------------------------------------------------------------------------
  do file_num_2d = 1, 28
    ! Define the filename for the current file
    filename_2d = trim(file_path_2d)//trim(file_list_2d(file_num_2d))//".txt"

    ! Open the current input file
    open(unit=10, file=filename_2d, status="old")
    
    ! Read the 2D array from the file
   do i = 1, ime
     read(10,*) (arr_2d(i,k), k=1,jme)
   end do
   
   ! Assign the values of the current input array to the appropriate element of input array
   input_2d(file_num_2d,:,:) = arr_2d(:,:)

  ! Close the current input file
  close(10)
end do

!-------------------------------------------------------------------------------------------
! Assign the 2D variables - costant in time 1-row 1-column
!------------------------------------------------------------------------------------------
   do i = 1, ime
        do j = 1, jme
            FRC_URB2D(i, j) = input_2d(1,i,j)
            xlong(i,j) = input_2d(5,i,j)
            xlat(i,j) = input_2d(6,i,j)
            cosz_urb2d(i,j) = input_2d(7,i,j)
            omg_urb2d(i,j) = input_2d(8,i,j)
            sf_ac_urb3d(i,j) = input_2d(9,i,j)
            lf_ac_urb3d(i,j) = input_2d(10,i,j)
            cm_ac_urb3d(i,j) = input_2d(11,i,j)
            sfvent_urb3d(i,j) = input_2d(12,i,j)
            lfvent_urb3d(i,j) = input_2d(11,i,j)
            ep_pv_urb3d(i,j) = input_2d(14,i,j)
            qgr_urb3d(i,j) = input_2d(15,i,j)
            tgr_urb3d(i,j) = input_2d(16,i,j)
            draingr_urb3d(i,j) = input_2d(17,i,j)
            rainbl(i,j) = 0.!input_2d(18,i,j)
            lp_urb2d(i,j) =input_2d(21,i,j)
            lb_urb2d(i,j) =input_2d(22,i,j)
            hgt_urb2d(i,j) = input_2d(23,i,j)
            rl_up(i,j) = input_2d(24,i,j)
            rs_abs(i,j) = input_2d(25,i,j)
            emiss(i,j) = input_2d(26,i,j)
            grdflx_urb(i,j) = input_2d(27,i,j)
        end do
    end do


 bx=2.*hgt_urb2d(1,1)*lp_urb2d(1,1)/(lb_urb2d(1,1)-lp_urb2d(1,1))               !_gl equation from BEP+BEM
 wx=2.*hgt_urb2d(1,1)*lp_urb2d(1,1)*((frc_urb2d(1,1)/lp_urb2d(1,1))-1.)/(lb_urb2d(1,1)-lp_urb2d(1,1))  !_gl equation from BEP+BEM

write(*,*) "urban morphology parameters UCP lp_urb2d", lp_urb2d
write(*,*) "urban morphology parameters UCP lb_urb2d", lb_urb2d
write(*,*) "urban morphology parameters UCP frc_urb2d", frc_urb2d
write(*,*) "urban morphology parameters UCP hgt_urb2d", hgt_urb2d
write(*,*) "urban morphology parameters UCP hi_urb2d", hi_urb2d
write(*,*) "latitude-longitude coordinate of cell", xlat, xlong
write(*,*) "bx,wx",bx,wx
!--------------------------------------------------------------------------------------------
! Read the variable from weather model 
!-------------------------------------------------------------------------------------------
do file_num_vector = 1, num_files
    ! Define the filename for the current file
    filename_vector = trim(file_path_vector)//trim(file_list_vector(file_num_vector))//".txt"
    
     ! Open the current input file
    open(unit=10, file=filename_vector, status="old")

    ! Read the vector from the file
    read(10,*) vector(1:n_loop)

    ! Assign the values of the current input vector to the appropriate element of input_vector
    input_vector(file_num_vector, :) = vector(1:n_loop)

    ! Close the current input file
    close(10)
end do


!----------------------------------------------------------------------------
!APRO I FILE PER SALVARE VARIBILI 2D & 3D
!----------------------------------------------------------------------------
OPEN(19, FILE='PRANDTL.txt')
OPEN(20, FILE='TW1_URB4D.txt')
OPEN(21, FILE='TW2_URB4D.txt')

OPEN(32, FILE='TRB_URB4D.txt')
OPEN(33, FILE='TRV_URB4D.txt')
OPEN(34, FILE='QR_URB4D.txt')
OPEN(35, FILE='TGB_URB4D.txt')

OPEN(61, FILE='SF_AC_URB3D.txt')
OPEN(62, FILE='LF_AC_URB3D.txt')
OPEN(63, FILE='CM_AC_URB3D.txt')

OPEN(104, FILE='TMR_11.txt')
OPEN(105, FILE='TMR_22.txt')
OPEN(106, FILE='TMR_13.txt')
OPEN(107, FILE='TMR_23.txt')
OPEN(108, FILE='COMF_10.txt')
OPEN(109, FILE='COMF_50.txt')
OPEN(110, FILE='COMF_90.txt')
!------------------------------------------------------------------------------------------
! CREZIONE DEGLI INDICI PER RILEGGERE I FILE IN OUTPUT
!-----------------------------------------------------------------------------------------
ind_zwd_mp=0 
ind_gd_mp=0
ind_zd_mp=0
ind_zdf_mp=0
ind_zrd_mp=0
ind_grd_mp=0
ind_bd_mp=0
ind_wd_mp=0
ind_gbd_mp=0
ind_fbd_mp=0

!IND_ZWD
 iii=0
 do ibui=1,nbui_max_mp
  do iz_u=1,nz_um_mp
   do iw=1,nwr_u_mp
      do id=1,ndm_mp
      iii=iii+1
      ind_zwd_mp(ibui,iz_u,iw,id)=iii
    enddo
   enddo
  enddo
enddo

! IND_GD 
iii=0
  do ig=1,ng_u_mp
   do id=1,ndm_mp
      iii=iii+1
      ind_gd_mp(ig,id)=iii
    enddo
   enddo

!IND_ZD
iii=0
  do ibui=1,nbui_max_mp
    do iz_u=1,nz_um_mp
      do id=1,ndm_mp
	iii=iii+1
	ind_zd_mp(ibui,iz_u,id)=iii
	enddo
     enddo
  enddo

!IND_ZRD
iii=0
  do iz_u=1,nz_um_mp
    do iw=1,nwr_u_mp
    do id=1,ndm_mp
       iii=iii+1
	ind_zrd_mp(iz_u,iw,id)=iii
     enddo
    enddo
  enddo

!IND_GRD
iii=0
  do iz_u=1,nz_um_mp
   do iw=1,ngr_u_mp
   do id=1,ndm_mp
     iii=iii+1
     ind_grd_mp(iz_u,iw,id)=iii
   enddo
  enddo
 enddo

!IND_ZDF
iii=0
  do iz_u=1,nz_um_mp
    do id=1,ndm_mp
       iii=iii+1
       ind_zdf_mp(iz_u,id)=iii
    enddo
  enddo

!IND_BD
iii=0
 do ibui=1,nbui_max_mp  !Type of building
   do iz_u=1,nz_um_mp     !vertical levels
      iii=iii+1
      ind_bd_mp(ibui,iz_u)=iii
   enddo
 enddo

!IND_WD
iii=0
do ibui=1,nbui_max_mp 	!type of building
  do iz_u=1,nz_um_mp 	!vertical levels
    do id=1,ndm_mp 	!direction
    iii=iii+1
    ind_wd_mp(ibui,iz_u,id)=iii
    enddo
  enddo
enddo

!IND_GBD
iii=0
 do ibui=1,nbui_max_mp	!type of building
    do iw=1,ngb_u_mp 	!layers in the wall (ground below a building)
     do id=1,ndm_mp 	!direction
        iii=iii+1
        ind_gbd_mp(ibui,iw,id)=iii  
     enddo
   enddo
 enddo

!IND_FBD
iii=0
  do ibui=1,nbui_max_mp !type of building
    do iw=1,nf_u_mp !layers in the wall (floor)
    do iz_u=1,nz_um_mp-1 !vertical levels
    do id=1,ndm_mp  !direction
       iii=iii+1
       ind_fbd_mp(ibui,iw,iz_u,id)=iii
    enddo
   enddo
   enddo
  enddo

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! COLUMN PROGRAM
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         !------------------------------------------------------------------------
         ! REPLACE READ_INPUT SUBROUTINE  !_gp
         !------------------------------------------------------------------------ 
         !AM below is just to have an initial value different than zero for the
         !different variables. In your case I would put vx,vy,pt, and tke equal to
         !those of the weather model at initialization.
           trc=0.
           cdz=0.0001
           emi=1.
           vx=input_vector(1,1)
           tke=temin
           vy=input_vector(2,1)
           pt=input_vector(3,1)
           time=0.
           time_pr=0.
           qv=input_vector(9,1)

           open(unit=9999,file='solar_position')
           open(unit=9998,file='surface_temperatures')
           open(unit=9997,file='time_check')
           write(9999,'(a3,11(2x,a15))')'doy','gmt','declin_urb','omg_urb2d','cosz_urb2d','sinh_S','h_S','swddir','swdown'
           write(9998,'(a3,11(2x,a15))')'doy','gmt','TGB_URB4D','temp','TW1_URB4D','TW2_URB4D','TRB_URB4D'
           write(*,*) "vx before loop time", vx

!=================== variabili solar position =====================
 Lat_rad = 0.
 Delta_TSL = 0. 
 h_S_mean = 0.
 theta_Z= 0.
 sinh_S = 0.  
 h_S=0.
 kkk=1

do i=kms,kme
   th0(i)=300.
enddo

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INIZIO CICLO NEL TEMPO
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
do itimestep =1,n_loop	
!AM here the loop over time starts. I would introduce in the loop a
!reading of the weather model variables (wind components, temperature,
!etc., every hour or every time they are available). 

   !=============================================================
   ! progress on the hour 
   if (mod(itimestep, 60) == 0) then
    gmt=gmt+1.
   end if
   
   ! progress on the day
   if (mod(gmt, 24.0) == 0 .and. mod(itimestep, 60)==0) then
      doy=doy+1
      gmt=0.0
   end if

    ! Check for leap years from start_year to current_year
        if (mod(count_year(kkk), 4) == 0 .and. (mod(count_year(kkk), 100) /= 0 .or. mod(count_year(kkk), 400) == 0)) then
          if (mod(doy, 367) == 0) then
              doy=1
              gmt=0.0
              kkk=kkk+1
          end if        
        else
           if (mod(doy, 366) == 0) then
              doy=1
              gmt=0.0
              kkk=kkk+1
           end if
        end if

 ! Check if it's the 60th step (every hour)
    if (itimestep == 1 .or. mod(itimestep,60) == 0) then
    	declin_urb= 23.45 * 3.14159265358979323846 / 180.0 * COS(2.0 * 3.14159265358979323846 / 365.0 * (172.0 - doy))

        ! Compute time difference between standard and local meridian
   	IF (xlong(1,1) < 0) THEN
      	  Delta_TSL = -1.0 / 15.0 * (15.0 * ABS(DeltaGMT) - ABS(xlong(1,1)))
        ELSE
      	  Delta_TSL = 1.0 / 15.0 * (15.0 * ABS(DeltaGMT) - ABS(xlong(1,1)))
   	END IF

   	Lat_rad = xlat(1,1) * deg2rad

	omg_urb2d = 15.0 * (gmt - 12) 
  	omg_urb2d = (omg_urb2d* deg2rad)

   ! Compute solar altitude
   	sinh_S = SIN(Lat_rad) * SIN(declin_urb) + COS(Lat_rad) * COS(declin_urb) * COS(omg_urb2d(1,1))
    	h_S = ASIN(sinh_S)

   ! Calculate solar zenith angle
   	theta_Z = pi / 2.0 - h_S
  
   ! Check and adjust theta_Z if it's outside the valid range
   	IF (theta_Z <= -pi / 2.0 .OR. theta_Z >= pi / 2.0) THEN
      	theta_Z = pi / 2.0
   	END IF
      	cosz_urb2d = cos(theta_Z)
     end if

!=====================================================
! weather model variables or reanalysis - every 1min
!=====================================================
   if (itimestep >= 1) then 
    DO j=nz-1,nz
      vx(j)=input_vector(1, i_wm)
      vy(j)=input_vector(2, i_wm)
      pt(j)=input_vector(3, i_wm)
      qv(j)=input_vector(9, i_wm)                       
    end do
end if
    
      !variables from WM 3D 
      DO i=ims,ime
         DO j=kms,kme
           DO  k=jms,jme
            u_phy(i, j, k) = vx(j)	
            v_phy(i, j, k) = vy(j)	
            th_phy(i, j, k) = pt(j)	
            qv_phy(i, j, k) = qv(j) 
           end do
        end do
       end do   

if (itimestep >= 1) then
   DO i=ims,ime
      DO j=kms,kme
           DO  k=jms,jme
               p_phy(i, j, k) = input_vector(7, i_wm)     !atmospheric pressure 
        end do
     end do
  end do
      !variables from WD 2D
      DO i=ims,ime
        DO  k=jms,jme
           swddir(i, k) = input_vector(5, i_wm)		
           glw(i, k) = input_vector(4, i_wm)
           swdown(i, k) = input_vector(6, i_wm)
           swddif(i, k) =  input_vector(10, i_wm)
           rainbl(i, k) = input_vector(8, i_wm)
       end do
      end do
!      PBLH = input_vector(11, i_wm)

    !condition on the zenith angle
    if(cosz_urb2d(1,1) <= 0.0001) then    
      swddir(1, 1) = 0. 
      swdown(1, 1) = 0.
      swddif(1, 1) =  0.
    end if
      i_wm = i_wm+1
  end if 

!-----------------------------------------------------------------------
           srex_vx=0.
           srim_vx=0.
           srex_vy=0.
           srim_vy=0.
           srex_pt=0.
           srim_pt=0.
           srex_tke=0.
           srim_tke=0.
           srex_trc=0.
           srim_trc=0.
           
           srex_qv=0.
           srim_qv=0.  

           time=time+dt
           time_pr=time_pr+dt
        
!introduce the call to BEP and BEM.
    DO i=ims,ime
         DO j=kms,kme
           DO  k=jms,jme
            b_u(i,j,k) = srex_vx(j)
            a_u(i,j,k) = srim_vx(j)
            b_v(i,j,k) = srex_vy(j)
            a_v(i,j,k) = srim_vy(j)
            b_t(i,j,k) = srex_pt(j)
            a_t(i,j,k) = srim_pt(j)
            b_e(i,j,k) = srex_tke(j)
            a_e(i,j,k) = srim_tke(j)
            b_q(i,j,k) = srex_qv(j)
           END DO
         END DO
     END DO

   CALL BEP_BEM(FRC_URB2D,UTYPE_URB2D,itimestep,dz8w,dt,u_phy,v_phy,    & !1
              th_phy,rho,p_phy,swdown,glw, 	 			& !2
              gmt,julday,xlong,xlat,       	   			& !3
              declin_urb,cosz_urb2d,omg_urb2d,		   		& !4
              num_urban_ndm,urban_map_zrd,urban_map_zwd,urban_map_gd, 	& !5
              urban_map_zd,urban_map_zdf,urban_map_bd,urban_map_wd,  	& !6
	      urban_map_gbd,urban_map_fbd,                              & !7
	      urban_map_zgrd,num_urban_hi,				& !8
	      trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,		        & !9
              tlev_urb3d,qlev_urb3d,tw1lev_urb3d,tw2lev_urb3d,	        & !10
	      tglev_urb3d,tflev_urb3d,sf_ac_urb3d,lf_ac_urb3d,    	& !11
              cm_ac_urb3d, 		 & !12
	      sfvent_urb3d,lfvent_urb3d, & !13
              sfwin1_urb3d,sfwin2_urb3d, & !14
              sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d, & !15
              ep_pv_urb3d,t_pv_urb3d, &  !16
              trv_urb4d,qr_urb4d,qgr_urb3d,tgr_urb3d, & !17
              drain_urb4d,draingr_urb3d,&  !18
              sfrv_urb3d,lfrv_urb3d, 	&  !19
              dgr_urb3d,dg_urb3d, 	&  !20
              tgv_urb4d,sfgv_urb3d,lfgv_urb3d,qg_urb4d,                    & !GARDEN
              tmr_11,tmr_12,tmr_13,tmr_21,tmr_22,tmr_23,                   & !COMFORT
              comf_10,comf_50,comf_90,hist_comf,num_urban_ncomf,           & !COMFORT
              lfr_urb3d,lfg_urb3d,rainbl,swddir,swddif, & !21
              lp_urb2d,hi_urb2d,lb_urb2d,hgt_urb2d, 	& !22
              a_u,a_v,a_t,a_e,b_u,b_v, & !23
              b_t,b_e,b_q,dlg,dl_u,sf_c,vl_c, 	& !24
	      rl_up,rs_abs,emiss,grdflx_urb,qv_phy, 	& !25
              history_interval,tsk_av,         &
              ids,ide,jds,jde,kds,kde, 	& !26
              ims,ime,jms,jme,kms,kme, 	& !27
              its,ite,jts,jte,kts,kte) 	  !28

           call dissip_bougeault(ix,iy,g,kms,kme,kts,kte,z_bou,dz_bou,tke,dlu,dld,pt,th0)

           call length_bougeault(ix,iy,kms,kme,kts,kte,dld,dlu,dlg,dl_u,dls,dlk)
                   
           call cdtur_bougeault(ix,iy,kms,kme,kts,kte,tke,z_bou,dz_bou,cdz,dlk)
 

    DO i=ims,ime
         DO j=kms,kme
           DO  k=jms,jme 
           srex_vx(j)= b_u(i,j,k)
           srim_vx(j)= a_u(i,j,k)
           srex_vy(j)= b_v(i,j,k)
           srim_vy(j)= a_v(i,j,k)
           srex_pt(j)= b_t(i,j,k)
           srim_pt(j)= a_t(i,j,k)
           srex_tke(j)= b_e(i,j,k) 
           srim_tke(j)= a_e(i,j,k)
           srex_qv(j) = b_q(i,j,k)
           END DO
         END DO
     END DO


  !Monin–Obukhov length (L) is estimated using formulation of Louis (Louis, 1979)

       call lmo_calc(dz,2.0,vx(1),vy(1),pt(1),300.,tsk_av,lmo)     
            z_ground(1)=0.
            zlmo(1)=0.

	DO k=2,nzm+1
         z_ground(k)=z_ground(k-1)+dz
         zlmo(k)=z_ground(k)/lmo

     	 if(zlmo(k)< -2)then         
          zlmo(k)=-2
         end if

        if(zlmo(k)> 2)then         
         zlmo(k)=2 
       end if
       END DO

!Vertical diffusion coefficient for scalar considering atmospheric stability Businger et al. (1971) and Dyer (1974)
       call calculate_prandtl(nzm+1,zlmo,prandtl)

       cdt=cdz/prandtl

       call tke_bougeault(nzm,nz,ceps,dz,vx,vy,pt, &
                   tke,cdz,cdt,dls,sf_c, &                      
                   srim_tke,srex_tke)

             srex_trc(1)=emi/dz/vl_c(1)
             srex_trc(nz)=-emi/dz
 
! compute the vertical diffusion       
           call diff(nzm,nz,1,2,dt,vx,cdz,srim_vx,srex_vx,sf_c,vl_c,dz,uw,duwdz)  
           call diff(nzm,nz,1,2,dt,vy,cdz,srim_vy,srex_vy,sf_c,vl_c,dz,vw,dvwdz)
           call diff(nzm,nz,1,2,dt,pt,cdt,srim_pt,srex_pt,sf_c,vl_c,dz,wt,dwtdz)
           call diff(nzm,nz,1,2,dt,qv,cdt,srim_qv,srex_qv,sf_c,vl_c,dz,dwq,dwqdz)  
           call diff(nzm,nz,1,1,dt,tke,cdz*3.5,srim_tke,srex_tke,sf_c,vl_c,dz,wtke,dwtkedz)  
           
           call diff(nzm,nz,1,2,dt,trc,cdt,srim_trc,srex_trc,sf_c,vl_c,dz,wtrc,dwtrcdz)

! avoid negative values for tke
             do i=1,nzm
                if(tke(i).lt.temin)then
                tke(i)=temin
               end if
             end do
 
!convert pt to real temperature
          r_rt=287.
          cp=1004.
          rcp= r_rt/cp
 
          pr_t(1)=p0_rt
          do iz=2,nz
           pr_t(iz) = (pr_t(iz-1)**rcp - g/cp * (p0_rt**rcp) * &
                (1/pt(iz) + 1/pt(iz-1)) * 0.5 * dz)**(1./rcp)
          end do
          do iz=1,nz
             tr(iz)=pt(iz)*(pr_t(iz)/p0_rt)**rcp;
          end do
        
         !printing files containing solar position and surfaces temperatures for different surfaces
         if (itimestep == 1 .or. mod(itimestep,30) == 0) then
         !if (mod(itimestep,30) == 0) then
            write(*,*)'doy=',doy
            write(9999,*)doy,gmt,declin_urb,omg_urb2d,cosz_urb2d,sinh_S*180/3.14,h_S*180/3.14,swddir,swdown
            write(9998,*)doy,gmt,tgb_urb4d(1,ind_gd_mp(10,1),1)-273.15,pt(1)-273.15,tr(1)-273.15, & 
                         tw1_urb4d(1,ind_zwd_mp(1,1,10,1),1)-273.15, &
                         tw2_urb4d(1,ind_zwd_mp(1,1,10,1),1)-273.15,trb_urb4d(1,ind_zrd_mp(6,10,1),1)-273.15
            write(9997,*)count_year(kkk),doy,gmt,itimestep,i_wm

            write(*,*) "ENTRA NEL CICLO DI STAMPA", itimestep,gmt,i_wm
            time_pr=0.
                    
!------------------------------------------------------------------------------
! SCRITTURA DEI DATI IN BASE AGLI INDICI
!-----------------------------------------------------------------------------
OPEN(502, FILE='TGR_urb3D.txt')

OPEN(601, FILE='U_component.txt')
OPEN(602, FILE='V_component.txt')
OPEN(603, FILE='Temperature.txt')
OPEN(604, FILE='Tke.txt')
OPEN(605, FILE='Sensible_flux.txt')
OPEN(606, FILE='Latent_flux.txt')
OPEN(607, FILE='Longwave_upward.txt')
OPEN(608, FILE='Shortwave_upward.txt')
OPEN(609, FILE='Momentum_flux_u.txt')
OPEN(610, FILE='Momentum_flux_v.txt')
OPEN(611, FILE='H.txt')
OPEN(612, FILE='LE.txt')
OPEN(616, FILE='QG.txt')
OPEN(617, FILE='TGARDEN.txt')
OPEN(618, FILE='TROAD.txt')
OPEN(613, FILE='SWDOWN.txt')
OPEN(614, FILE='GLW.txt')
OPEN(615, FILE='SW_abs.txt')


!STAMPA
write(601,'(24(1x,f10.2))') vx(40),vx(41) 	!impostare il livello verticale considerandone l'altezza in metri 
write(602,'(24(1x,f10.2))') vy(40),vy(41)
write(603,'(24(1x,f10.2))') pt(2),pt(41)
write(604,'(24(1x,f10.2))') tke(40)

write(61,'(24(1x,f10.2))') sf_ac_urb3d 
write(62,'(24(1x,f10.2))') lf_ac_urb3d
write(63,'(24(1x,f10.2))') cm_ac_urb3d


write(104,'(24(1x,f10.2))') tmr_11
write(105,'(24(1x,f10.2))') tmr_22
write(106,'(24(1x,f10.2))') tmr_13
write(107,'(24(1x,f10.2))') tmr_23
write(108,'(24(1x,f10.2))') comf_10
write(109,'(24(1x,f10.2))') comf_50
write(110,'(24(1x,f10.2))') comf_90

! faccio l'integrale del flussi sull'intera colonna
 H_plumber=0.0
 LE_plumber=0.0

  do i = 1, kme-2
    H_plumber = H_plumber + srex_pt(i)
    LE_plumber = LE_plumber + srex_qv(i)
  end do

LE_plumber = LE_plumber*1000


! stampa flussi per validazione progetto pu
write(605,'(24(1x,f10.3))') H_plumber
write(606,'(24(1x,f10.3))') LE_plumber
write(607,'(24(1x,f10.3))') rl_up
write(609,'(24(1x,f10.3))') uw(40)
write(610,'(24(1x,f10.3))') vw(40)
write(611,'(24(1x,f10.3))') wt(40)
write(612,'(24(1x,f10.3))') dwq(40)*1000
write(613,'(24(1x,f10.3))') swdown
write(19,'(24(1x,f10.3))') prandtl(5)
write(614,'(24(1x,f10.3))') glw
write(615,'(24(1x,f10.3))') rs_abs
write(616,'(24(2x,f10.3))') qg_urb4d(1,size(qg_urb4d),1),rainbl
write(617,'(24(2x,f10.3))') tgv_urb4d(1,size(tgv_urb4d),1)
write(618,'(24(2x,f10.3))') tgb_urb4d(1,size(tgb_urb4d),1)
id=1        !strada

!do iz=1,nz_um_mp
!   if(iz.eq.39)then !un livello verticale superiore a quello indicato 
!     write(42,'(24(1x,f10.3))') t_pv_urb3d(ix,ind_zdf_mp(iz,id),iy), t_pv_urb3d(ix,ind_zdf_mp(iz,id+1),iy) !"PHOTOVOLTAIC PANELS TEMPERATURE " "K"
!     write(43,'(24(1x,f10.3))')	sfr_urb3d(ix,ind_zdf_mp(iz,id),iy), sfr_urb3d(ix,ind_zdf_mp(iz,id+1),iy)   !"SENSIBLE HEAT FLUX FROM URBAN SFC"
!     write(46,'(24(1x,f10.3))')	dgr_urb3d(ix,ind_zdf_mp(iz,id),iy), dgr_urb3d(ix,ind_zdf_mp(iz,id+1),iy)   !"ROOF LAYER DEPTH WATER RETENTION" mm  
!     write(47,'(24(1x,f10.3))')	lfr_urb3d(ix,ind_zdf_mp(iz,id),iy), lfr_urb3d(ix,ind_zdf_mp(iz,id+1),iy)   !"LATENT HEAT FLUX FROM URBAN SFC"
!!	if(gr_flag_u.eq.1)then
!     write(44,'(24(1x,f10.3))')	sfrv_urb3d(ix,ind_zdf_mp(iz,id),iy), sfrv_urb3d(ix,ind_zdf_mp(iz,id+1),iy) !SENSIBLE HEAT FLUX FROM GREEN ROOF
!     write(45,'(24(1x,f10.3))')	lfrv_urb3d(ix,ind_zdf_mp(iz,id),iy), lfrv_urb3d(ix,ind_zdf_mp(iz,id+1),iy) !"LATENT HEAT FLUX FROM GREEN ROOF"
!     write(41,'(24(1x,f10.3))')	drain_urb4d(ix,ind_zdf_mp(iz,id),iy), drain_urb4d(ix,ind_zdf_mp(iz,id+1),iy) !"GREEN ROOF DRAINAGE"   "mm"
!      endif
!   enddo

do ibui=1,nbui_max_mp
   do iz_u=1,nz_um_mp
     do iw=1,nwr_u_mp
        if(ibui.eq.1.and.iz_u.eq.1.and.iw.eq.10)then
        write(20,'(24(1x,f10.3))') tw1_urb4d(ix,ind_zwd_mp(ibui,iz_u,iw,id),iy), tw1_urb4d(ix,ind_zwd_mp(ibui,iz_u,iw,id+1),iy) !"WALL LAYER TEMPERATURE"
        write(21,'(24(1x,f10.3))') tw2_urb4d(ix,ind_zwd_mp(ibui,iz_u,iw,id),iy), tw2_urb4d(ix,ind_zwd_mp(ibui,iz_u,iw,id+1),iy)
       endif
      enddo
    enddo
  enddo

do ig=1,ng_u_mp
   if(ig.eq.10)then
   write(35,'(24(1x,f10.3))') tgb_urb4d(ix,ind_gd_mp(ig,id),iy), tgb_urb4d(ix,ind_gd_mp(ig,id+1),iy) !"ROAD LAYER TEMPERATURE"
   endif
enddo

!do id=1,ndm   !gl
do iz_u=1,nz_um_mp
 do ir=1,nwr_u_mp
 !   if(iz_u.eq.1+1.and.ir.eq.10)then !un livello nz_um  superiore a quello indagato
    if(iz_u.eq.11.and.ir.eq.10)then !un livello nz_um  superiore a quello indagato
!    write(*,*) "write index of trb_urb4d", iz_u, ir, id
    write(32,'(24(1x,f10.3))') trb_urb4d(ix,ind_zrd_mp(iz_u,ir,id),iy), trb_urb4d(ix,ind_zrd_mp(iz_u,ir,id+1),iy) !"ROOF LAYER TEMPERATURE"
    write(33,'(24(1x,f10.3))') trv_urb4d(ix,ind_grd_mp(iz_u,ir,id),iy), trv_urb4d(ix,ind_grd_mp(iz_u,ir,id+1),iy) !"GREEN ROOF LAYER TEMPERATURE"
    write(34,'(24(1x,f10.3))') qr_urb4d(ix,ind_grd_mp(iz_u,ir,id),iy), qr_urb4d(ix,ind_grd_mp(iz_u,ir,id+1),iy) !"GREEN ROOF LAYER MOISTURE"
    endif
   enddo
 enddo

!Outputs of BEM
!do ibui=1,nbui_max_mp !type of building
!  ! do iz_u=1,nz_um_mp !vertical levels
!      if(ibui.eq.1.and.iz_u.gt.1)then     !extract data of five  nz_um and ibui=1 
!      write(50,'(24(1x,f10.3))') tlev_urb3d(ix,ind_bd_mp(ibui,1),iy),tlev_urb3d(ix,ind_bd_mp(ibui,2),iy), &
!   tlev_urb3d(ix,ind_bd_mp(ibui,3),iy),tlev_urb3d(ix,ind_bd_mp(ibui,4),iy),tlev_urb3d(ix,ind_bd_mp(ibui,5),iy)
!    write(51, '(24(1x,f10.3))') qlev_urb3d(ix,ind_bd_mp(ibui,1),iy),qlev_urb3d(ix,ind_bd_mp(ibui,2),iy), &
!   qlev_urb3d(ix,ind_bd_mp(ibui,3),iy),qlev_urb3d(ix,ind_bd_mp(ibui,4),iy),qlev_urb3d(ix,ind_bd_mp(ibui,5),iy)
!      endif  
!  !  enddo
!enddo

!do ibui=1,nbui_max_mp !type of building
!  do iz_u=1,nz_um_mp !vertical levels
!    if(ibui.eq.1.and.iz_u.eq.1)then  	!extract data from nz_um=1 ibui=1 on both street directions
! 	!write(*,*) "nuovo loop", ibui,iz_u,id 
!  	!write(*,*) "ultimo indice", ind_wd_mp(ibui,iz_u,id)
!  	write(112,'(24(1x,f10.3))') tw1lev_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), tw1lev_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1) !"WINDOW TEMPERATURE"
!        write(113,'(24(1x,f10.3))') tw2lev_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), tw2lev_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1)
!        write(114,'(24(1x,f10.3))') sfwin1_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), sfwin1_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1) !"SENSIBLE HEAT FLUX FROM URBAN SFC WINDOW"
!        write(115,'(24(1x,f10.3))') sfwin2_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), sfwin2_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1)
!        write(116,'(24(1x,f10.3))') sfw1_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), sfw1_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1) !SENSIBLE HEAT FLUX FROM URBAN SFC
!        write(117,'(24(1x,f10.3))') sfw2_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), sfw2_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1)
!  endif
! enddo
!enddo

!do ibui=1,nbui_max_mp  !type of building
!  do iw=1,ngb_u_mp !layers in the walls
!     if(ibui.eq.1.and.iw.eq.10)then   !Ground temperature below a building in BEM
!       write(200,'(24(1x,f10.3))') tglev_urb3d(ix,ind_gbd_mp(ibui,iw,id),iy), tglev_urb3d(ix,ind_gbd_mp(ibui,iw,id+1),iy) 
!     endif
!  enddo
!enddo

!do ibui=1,nbui_max_mp !type of building
!  do iw=1,nf_u_mp !layer in the walls
!   do iz_u=1,nz_um_mp-1 !verticals levels
!    if(ibui.eq.1.and.iw.eq.10.and.iz_u.eq.1)then   !FLOOR TEMPERATURE
!    write(201,'(24(1x,f10.3))') tflev_urb3d(ix,ind_fbd_mp(ibui,iw,iz_u,id),iy), tflev_urb3d(ix,ind_fbd_mp(ibui,iw,iz_u,id+1),iy)
!    endif
!    enddo
!   enddo
!  enddo

end if  !condizione sul ciclo stampa
enddo 	!loop su subroutine bep_bem (loop sul tempo del modello colonna)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ULTIMA PARTE PROGRAM COLUMN ALBERTO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           close(66)
           trc=trc-trc(nz)
           trc_2(nz)=0.
           do iz_c=nz-1,1,-1
            trc_an(iz_c)=trc_an(iz_c+1)+dz/cdz(iz_c+1)*(emi-dwtrc(iz_c+1))
            trc_2(iz_c)=trc_2(iz_c+1)+dz/cdz(iz_c+1)*(emi)
           enddo
           if(iconfig.eq.1)open(unit=68,file='output_staggered_r')
           if(iconfig.eq.2)open(unit=68,file='output_aligned_r')
            write(68,'(a3,11(2x,a15))')'iz_c','vx/utau','tke/utau**2.','uwm/utau**2.','trc','wtrc','cdz','dls'
            do iz_c=1,nz
              uwm=(uw(iz_c)+uw(iz_c+1))/2.
            write(68,'(i3,11(2x,f15.8))')iz_c,vx(iz_c),tke(iz_c), &
                    uwm,trc(iz_c),wtrc(iz_c),cdz(iz_c),dls(iz_c)
            enddo

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
close(20)
close(21)
close(32)
close(33)
close(34)
close(35)
close(41)
close(42)
close(43)
close(44)
close(45)
close(46)
close(47)
close(104)
close(105)
close(106)
close(107)
close(108)
close(109)
close(110)

close(112)
close(113)
close(114)
close(115)
close(116)
close(117)
close(200)
close(201)

END PROGRAM main_program

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!SUBROUTINE PROGRAM COLUMN
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      
      subroutine cdtur(nzm,nz,dz,ck,tke,dlk,cdz)                                  


      !implicit none

! Input
! -----
      integer nzm                              ! maximum number of vertical levels
      integer nz                               ! number of vertical levels
      real dz                                  ! levels size [m]
      real ck                                  ! von Karman constant
      real tke(nzm)                            ! turbulent kinetic energy
      real dlk(nzm)                            ! lenth scale

! Output
! ------
      real cdz(nzm+1)                            ! diffusion coefficient

! Local
! -----
      integer iz_c
      real tke_m
      real dlk_m

! ----------------------------------------------------------------------

       cdz(1)=0.

!       do iz=2,nz-1
       do iz_c=2,nz
        tke_m=(tke(iz_c-1)+tke(iz_c))/2.
        dlk_m=(dlk(iz_c-1)+dlk(iz_c))/2.
        cdz(iz_c)=ck*dlk_m*sqrt(tke_m)
       enddo
       cdz(nz+1)=cdz(nz)

       return
       end  

! ===6================================================================72

       subroutine tke_bougeault (nzm,nz,ceps,dz,vx,vy,pt, &
                                tke,cdz,cdt,dls,sf_c, &                       ! Input
                                srim_tke,srex_tke)                            ! Output


! ----------------------------------------------------------------------
!   Calculation of the sources (shear) and the dissipation
!    terms for the TKE equation.
! ----------------------------------------------------------------------

      !implicit none

! Input
! -----
      integer nz,nzm                               ! number of vertical levels
    
      real ceps
      real dz                                  ! levels size [m]
      real vx(nzm)                              ! wind velocity
      real vy(nzm)                              ! wind velocity
      real pt(nzm)     				! potential temperature _gl
      real tke(nzm)                             ! turbulent kinetic energy
      real cdz(nzm+1),cdt(nzm+1)                           ! turbulent diffusion coefficient
      real dls(nzm)                             ! lrngth scale leps
      real sf_c(nzm+1)                            ! ?????
   
! Ouput
! -----
      real srim_tke(nzm)                           ! inplicit term in the tke equation
      real srex_tke(nzm)                           ! explicit term in the tke equation

! Local
! -----
      integer iz_c
      real sh(nzm)                             ! shear turbulent kinetic energy source
      real td(nz)                              ! dissipation
      real bu(nzm)                             !buoyancy _gl
! ----------------------------------------------------------------------
 
      call shear_2(nzm,nz,dz,vx,vy,cdz, &                               ! Input
                  sh)                                                   ! Ouput

     call buoy(nzm,nz,dz,pt,g,cdt,bu)
          
       do iz_c=1,nz
       if (dls(iz_c).ne.0.) then
        td(iz_c)=-ceps*sqrt(tke(iz_c))/dls(iz_c)
       else
        td(iz_c)=0.
       end if 
       sh(iz_c)=sh(iz_c)*sf_c(iz_c)
       bu(iz_c)=bu(iz_c)*sf_c(iz_c)
       srim_tke(iz_c)=td(iz_c)
       srex_tke(iz_c)=srex_tke(iz_c)+sh(iz_c)+bu(iz_c) 
      end do 
!    write(*,*) "buoyancy", bu
      return
      end  ! 

! ===6================================================================72
! ===6================================================================72

      subroutine shear(nzm,nz,dz,vx,vy,cdz, &                               ! Input
                       sh)                                              ! Ouput



! ----------------------------------------------------------------------
!     Calculation of the shear source for the TKE equation.
! ----------------------------------------------------------------------

      !implicit none

! Input
! -----
      integer nzm,nz                            ! number of vertical levels
      real dz                              	! levels size [m]
      real vx(nzm)                              ! wind velocity
      real vy(nzm)                              ! wind velocity
      real cdz(nzm+1)                           ! turbulent diffusion coefficient

! Ouput
! -----
      real sh(nzm)                              ! shear turbulent kinetic energy source

! Local
! -----
      integer iz_c
      real dudz1
      real dvdz1
      real dudz2
      real dvdz2
      real cdm
      real dumdz
      real cdmmin
      parameter (cdmmin=0.01)

! ----------------------------------------------------------------------

      sh(1)=0.

      do iz_c=2,nz-1      
        dudz1=(vx(iz_c)-vx(iz_c-1))/dz
        dvdz1=(vy(iz_c)-vy(iz_c-1))/dz
        dudz2=(vx(iz_c+1)-vx(iz_c))/dz
        dvdz2=(vy(iz_c+1)-vy(iz_c))/dz

        cdm=max(0.5*(cdz(iz_c)+cdz(iz_c+1)),cdmmin)

        dumdz=0.5*((dudz1**2.+dvdz1**2)+(dudz2**2.+dvdz2**2))

        sh(iz_c)=cdm*dumdz
      enddo  ! iz

      sh(nz)=0.     

      return ! shear
      end

! ===6================================================================72 
! ===6================================================================72

      subroutine shear_2(nzm,nz,dz,vx,vy,cdz, &                               ! Input
                       sh)                                              ! Ouput



! ----------------------------------------------------------------------
!     Calculation of the shear source for the TKE equation.
! ----------------------------------------------------------------------

      !implicit none

! Input
! -----
      integer nzm,nz                               ! number of vertical levels
      real dz                              ! levels size [m]
      real vx(nzm)                              ! wind velocity
      real vy(nzm)                              ! wind velocity
      real cdz(nzm+1)                           ! turbulent diffusion coefficient

! Ouput
! -----
      real sh(nzm)                              ! shear turbulent kinetic energy source

! Local
! -----
      integer iz_c
      real dudz1
      real dvdz1
      real dudz2
      real dvdz2
      real cd1,cd2

     

! ----------------------------------------------------------------------

      sh(1)=0.

      do iz_c=2,nz-1    
        dudz1=(vx(iz_c)-vx(iz_c-1))/dz
        dvdz1=(vy(iz_c)-vy(iz_c-1))/dz
        dudz2=(vx(iz_c+1)-vx(iz_c))/dz
        dvdz2=(vy(iz_c+1)-vy(iz_c))/dz
     
        sh(iz_c)=0.5*(cdz(iz_c)*(dudz1**2.+dvdz1**2)+cdz(iz_c+1)*(dudz2**2.+dvdz2**2))
      enddo  ! iz

      sh(nz)=0.    

      return ! shear_2
      end

! ===6================================================================72 
! ===6=8===============================================================72

       subroutine diff(nzm,nz,iz1,izf,dt,co,cd,aa,bb,sf_c,vl_c,dz,fc,df)

!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +           Subroutine prepared by A.Martilli                     +
!     +        Ecole Polytechnique Federale de Lausanne                 +
!     +   DGR - IGE - Laboratoire de Pollution de l'Air et des Sols     +
!     +            tel.: (021)-693-61-60                                +
!     +            Email: alberto.martilli@dgr.epfl.ch                  + 
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!------------------------------------------------------------------------
!           Calculation of the diffusion in 1D        
!------------------------------------------------------------------------
!  - Input:
!       nz    : number of points
!       iz1   : first calculated point
!       co    : concentration of the variable of interest
!       dz    : vertical levels
!       cd    : diffusion coefficients
!       dtext : external time step
!       itest : if itest eq 1 then update co, else store in a flux array
!  - Output:
!       co :concentration of the variable of interest

!  - Internal:
!       cddz  : constant terms in the equations 
!       dt    : diffusion time step
!       nt    : number of the diffusion time steps
!       cstab : ratio of the stability condition for the time step
!---------------------------------------------------------------------

         !implicit none
                
         integer nz,nzm,iz_c,iz1,izf
         real co(nzm),cd(nzm+1),dz,dt,dzv

         real cddz(nzm+2),fc(nzm+1),df(nzm)
         real a(nzm,3),c(nzm)
         real sf_c(nzm+1),vl_c(nzm)
         real aa(nzm),bb(nzm)
        
! Compute cddz=2*cd/dz  
        
        cddz(1)=sf_c(1)*cd(1)/dz
        do iz_c=2,nz
         cddz(iz_c)=2.*sf_c(iz_c)*cd(iz_c)/(2.*dz)
        enddo
        if(izf.gt.0)then
         cddz(nz+1)=sf_c(nz+1)*cd(nz+1)/dz
        else
         cddz(nz+1)=0.
        endif
         do iz_c=1,iz1-1
          a(iz_c,1)=0.
          a(iz_c,2)=1.
          a(iz_c,3)=0.
          c(iz_c)=co(iz_c)
         enddo
          
          do iz_c=iz1,nz-izf   !_gl
           dzv=vl_c(iz_c)*dz
           a(iz_c,1)=-cddz(iz_c)*dt/dzv         
           a(iz_c,2)=1+dt*(cddz(iz_c)+cddz(iz_c+1))/dzv-aa(iz_c)*dt
           a(iz_c,3)=-cddz(iz_c+1)*dt/dzv            
           c(iz_c)=co(iz_c)+bb(iz_c)*dt                     
          enddo

          if(izf.eq.1)then
           dzv=vl_c(nz)*dz
           a(nz,1)=-cddz(nz)*dt/dzv        
           a(nz,2)=1+dt*(cddz(nz))/dzv-aa(nz)*dt
           a(nz,3)=0.            
           c(nz)=co(nz)+bb(nz)*dt   
          else
           do iz_c=nz-izf+1,nz
            a(iz_c,1)=0.
            a(iz_c,2)=1.
            a(iz_c,3)=0.
            c(iz_c)=co(iz_c)
           enddo
          endif

          call invert (nzm,nz,a,c,co)
                  
          do iz_c=1,iz1 
           fc(iz_c)=0.
          enddo
          
          do iz_c=iz1+1,nz 
           fc(iz_c)=-(cddz(iz_c)*(co(iz_c)-co(iz_c-1)))
          enddo
        
          do iz_c=1,iz1
           df(iz_c)=0.
          enddo
          
          do iz_c=iz1+1,nz-izf   
           dzv=vl_c(iz_c)*dz
           if(iz.lt.nz)then
            df(iz_c)=+(co(iz_c-1)*cddz(iz_c)&
                 -co(iz_c)*(cddz(iz_c)+cddz(iz_c+1))&
                 +co(iz_c+1)*cddz(iz_c+1))/dzv
           else
            df(iz_c)=+(co(iz_c-1)*cddz(iz_c)&
                 -co(iz_c)*(cddz(iz_c)+cddz(iz_c+1)))/dzv
           endif
          enddo
          
          do iz_c=nz-izf,nz    
           df(iz_c)=0.
         enddo
                                        
       return
       end
!----------------------------------------------------------------------------

       subroutine invert(nzm,nn,a,c,x)
       !implicit none
       
!ccccccccccccccccccccccccccccccc       
! Aim: INversion and resolution of a tridiagonal matrix
!          A X = C
! Input:
!  a(*,1) lower diagonal (Ai,i-1)
!  a(*,2) principal diagonal (Ai,i)
!  a(*,3) upper diagonal (Ai,i+1)
!  c      
! Output
!  x     results
!ccccccccccccccccccccccccccccccc 
       integer nzm,nn,in
       real a(nzm,3),c(nzm),x(nzm)                       
        
        do in=nn-1,1,-1                 
         c(in)=c(in)-a(in,3)*c(in+1)/a(in+1,2)
         a(in,2)=a(in,2)-a(in,3)*a(in+1,1)/a(in+1,2)
        enddo
        
        do in=2,nn        
         c(in)=c(in)-a(in,1)*c(in-1)/a(in-1,2)
        enddo
        
        do in=1,nn
         x(in)=c(in)/a(in,2)
        enddo

        return
        end
          
         

! ===6=8===============================================================72

      subroutine gaussj(a,n,b,np)

! ----------------------------------------------------------------------
! This routine solve a linear system of n equations of the form
!              A X = B
!  where  A is a matrix a(i,j)
!         B a vector and X the solution
! In output b is replaced by the solution
! ----------------------------------------------------------------------

!      include 'wind.h'

! ----------------------------------------------------------------------
! INPUT:
! ----------------------------------------------------------------------
      integer np
      real a(np,np)

! ----------------------------------------------------------------------
! OUTPUT:
! ----------------------------------------------------------------------
      real b(np)

! ----------------------------------------------------------------------
! LOCAL:
! ----------------------------------------------------------------------
      integer nm
      parameter (nm=150)

      real*8 big,dum
      integer i_gauss,icol,irow
      integer j,k,l,ll,n
      integer ii,jj
      integer ipiv(nm)
      real*8 pivinv

! ----------------------------------------------------------------------
! END VARIABLES DEFINITIONS
! ----------------------------------------------------------------------


      do j=1,n
         ipiv(j)=0.
      enddo

      do i_gauss=1,n
         big=0.
         do j=1,n
            if(ipiv(j).ne.1)then
               do k=1,n
                  if(ipiv(k).eq.0)then
                     if(abs(a(j,k)).ge.big)then
                        big=abs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  elseif(ipiv(k).gt.1)then
                     pause 'singular matrix in gaussj'
                  endif
               enddo
            endif
         enddo

         ipiv(icol)=ipiv(icol)+1

         if(irow.ne.icol)then
            do l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
            enddo

            dum=b(irow)
            b(irow)=b(icol)
            b(icol)=dum

         endif

         if(a(icol,icol).eq.0)pause 'singular matrix in gaussj'

         pivinv=1./a(icol,icol)
         a(icol,icol)=1

         do l=1,n
            a(icol,l)=a(icol,l)*pivinv
         enddo

         b(icol)=b(icol)*pivinv

         do ll=1,n
            if(ll.ne.icol)then
               dum=a(ll,icol)
               a(ll,icol)=0.
               do l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
               enddo

               b(ll)=b(ll)-b(icol)*dum

            endif
         enddo


      enddo

      return
      end !subroutine gaussj

!===========================================================================
! SUBROUTINE TEMP-BUOYANCY
!==========================================================================

subroutine buoy(nzm,nz,dz,th,g,cdz,bu)

! compute buoyancy term 
integer nzm, nz 
integer iz 
real dz
real dtdz1,dtdz2,cdm,dtmdz,g 
real th(nzm),cdz(nzm+1),bu(nzm) 

       bu(1)=0.	!_gl
       g=9.81
   
    do iz=2,nzm-1
          dtdz1=(th(iz)-th(iz-1))/dz
          dtdz2=(th(iz+1)-th(iz))/dz

          dtmdz=0.5*(dtdz1+dtdz2) 
          cdm=0.5*(cdz(iz+1)+cdz(iz)) 
          
          bu(iz)=-cdm*dtmdz*g/300.
     !  enddo
if(iz>=nzm-2)then
 !      write(*,*) "this is value of buoyancy inside the subrotune", bu(iz), cdm, dtmdz, dz          
 !      write(*,*) "this is value of th(iz)-th(iz-1) inside the subrotune", th(iz), th(iz+1), th(iz)-th(iz+1)
 !      write(*,*) "this is value of dtmdz inside the subrotune", dtmdz
 !      write(*,*) "this is value of g inside the subrotune", g
   end if 
      enddo       


              
       return
       end

SUBROUTINE calculate_prandtl(nzm,zlmo,prandtl)


       IMPLICIT NONE
       INTEGER :: nzm
       REAL, INTENT(IN)  :: zlmo(nzm)
       REAL, INTENT(OUT) :: prandtl(nzm) 
       INTEGER :: k
       DO k=1,nzm
       IF(zlmo(k).le.0)THEN
       prandtl(k)=(0.74*(1-9.*zlmo(k))**(-0.5))/((1-15.*zlmo(k))**(-0.25))
       ELSE
!       prandtl(k)=0.74 + (0.26 / (1.0 + EXP(-10.0 * (zlmo(k) - 0.5))))
       prandtl(k)=(0.74 + 4.7*zlmo(k))/(1+4.7*zlmo(k))
       ENDIF
       ENDDO

      
    END SUBROUTINE calculate_prandtl


     subroutine lmo_calc(dz,z0,ua,va,pt,pt0,ptg,lmo)

!----------------------------------------------------------------------
!           Calculation of the flux at the ground
!           Formulation of Louis (Louis, 1979)
!----------------------------------------------------------------------

      implicit none



! ----------------------------------------------------------------------
! INPUT:
! ----------------------------------------------------------------------
      real dz                   ! first vertical level
      real pt                   ! potential temperature
      real pt0                  ! reference potential temperature
      real ptg                  ! ground potential temperature
      real ua                   ! wind speed
      real va                   ! wind speed
      real z0                   ! Roughness length

! ----------------------------------------------------------------------
! OUTPUT:
! ----------------------------------------------------------------------
! Explicit component of the momentum, temperature and TKE sources or sinks on horizontal
!  surfaces (roofs and street)
! The fluxes can be computed as follow: Fluxes of X = B
!  Example: Momentum fluxes on horizontal surfaces =  uhb_u


! ----------------------------------------------------------------------
! LOCAL:
! ----------------------------------------------------------------------
      real aa
      real al
      real buu
      real c
      real fbuw
      real fbpt
      real fh
      real fm
      real ric
      real tstar
      real ustar
      real utot
      real wstar
      real zz
      real lmo
      real b,cm,ch,rr,tol
      parameter(b=9.4,cm=7.4,ch=5.3,rr=0.74,tol=.001)

! ----------------------------------------------------------------------
! END VARIABLES DEFINITIONS
! ----------------------------------------------------------------------


! computation of the ground temperature

      utot=(ua**2+va**2)**.5


!!!! Louis formulation
!
! compute the bulk Richardson Number

      zz=dz/2.

!        if(tstar.lt.0.)then
!         wstar=(-ustar*tstar*g*hii/pt)**(1./3.)
!        else
!         wstar=0.
!        endif
!
!      if (utot.le.0.7*wstar) utot=max(0.7*wstar,0.00001)

      utot=max(utot,0.01)

      ric=2.*9.81*zz*(pt-ptg)/((pt+ptg)*(utot**2))

      aa=0.4/log(zz/z0)

! determine the parameters fm and fh for stable, neutral and unstable conditions

      if(ric.gt.0)then
         fm=1/(1+0.5*b*ric)**2
         fh=fm
      else
         c=b*cm*aa*aa*(zz/z0)**.5
         fm=1-b*ric/(1+c*(-ric)**.5)
         c=c*ch/cm
         fh=1-b*ric/(1+c*(-ric)**.5)
      endif

      fbuw=-aa*aa*utot*utot*fm
      fbpt=-aa*aa*utot*(pt-ptg)*fh/rr

      ustar=(-fbuw)**.5
      tstar=-fbpt/ustar

      lmo=(0.4*9.81*tstar)/(pt*ustar*ustar)

!!!!!!!!!!!!!!!

      return
      end subroutine lmo_calc


! ===6=8===============================================================72

! ===6=8===============================================================72
         subroutine dissip_bougeault(ix,iy,g,kms,kme,kts,kte,z,dz,te,dlu,dld,th,th0)

! compute the length scales up and down
         implicit none
         integer kms,kme,kts,kte,iz,izz,ix,iy
         real dzt,zup,beta,zup_inf,bbb,tl,zdo,zdo_sup,zzz,g
         real te(kms:kme),dlu(kms:kme),dld(kms:kme),dz(kms:kme)
         real th(kms:kme),th0(kms:kme),z(kms:kme)

        do iz=kts,kte
          zup=0.
          dlu(iz)=z(kte+1)-z(iz)-dz(iz)/2.
          zzz=0.
          zup_inf=0.
          beta=g/th0(iz)      !Buoyancy coefficient
          do izz=iz,kte-1
           dzt=(dz(izz+1)+dz(izz))/2.
           zup=zup-beta*th(iz)*dzt
           zup=zup+beta*(th(izz+1)+th(izz))*dzt/2.
           zzz=zzz+dzt
           if(te(iz).lt.zup.and.te(iz).ge.zup_inf)then
            bbb=(th(izz+1)-th(izz))/dzt
            if(bbb.ne.0)then
             tl=(-beta*(th(izz)-th(iz))+sqrt( max(0.,(beta*(th(izz)-th(iz)))**2.+2.*bbb*beta*(te(iz)-zup_inf))))/bbb/beta
            else
             if(th(izz).ne.th(iz))then
              tl=(te(iz)-zup_inf)/(beta*(th(izz)-th(iz)))
             else
              tl=0.
             endif
            endif            
            dlu(iz)=zzz-dzt+tl
           endif
           zup_inf=zup
          enddo
                  
          zdo=0.
          zdo_sup=0.
          dld(iz)=z(iz)+dz(iz)/2.
          zzz=0.
          do izz=iz,kts+1,-1
           dzt=(dz(izz-1)+dz(izz))/2.
           zdo=zdo+beta*th(iz)*dzt
           zdo=zdo-beta*(th(izz-1)+th(izz))*dzt/2.
           zzz=zzz+dzt
           if(te(iz).lt.zdo.and.te(iz).ge.zdo_sup)then
            bbb=(th(izz)-th(izz-1))/dzt
            if(bbb.ne.0.)then
             tl=(beta*(th(izz)-th(iz))+sqrt( max(0.,(beta*(th(izz)-th(iz)))**2.+2.*bbb*beta*(te(iz)-zdo_sup))))/bbb/beta
            else
             if(th(izz).ne.th(iz))then
              tl=(te(iz)-zdo_sup)/(beta*(th(izz)-th(iz)))
             else
              tl=0.
             endif
            endif
            
            dld(iz)=zzz-dzt+tl
           endif
           zdo_sup=zdo
          enddo
         enddo
        ! write(*,*)  "value at end of g inside dissip_bouge", g    
                   
         end subroutine dissip_bougeault

!=================================================================================
!=================================================================================

      subroutine length_bougeault(ix,iy,kms,kme,kts,kte,dld,dlu,dlg,dl_u,dls,dlk) 
  ! compute the length scales for dissipation and turbulent coefficients implicit none 
  
         integer kms,kme,kts,kte,iz,ix,iy 
         real dlu(kms:kme),dld(kms:kme),dl_u(kms:kme) 
         real dls(kms:kme),dlk(kms:kme),dlg(kms:kme)
         
         do iz=kts,kte 
           dld(iz)=min(dld(iz),dlg(iz)) 
           dls(iz)=sqrt(dlu(iz)*dld(iz)) 
           dlk(iz)=min(dlu(iz),dld(iz)) 
           if(dl_u(iz).gt.0.)then
           dls(iz)=1./(1./dls(iz)+1./dl_u(iz)) 
           dlk(iz)=1./(1./dlk(iz)+1./dl_u(iz))
           endif 
         enddo
                   
         return
         end subroutine length_bougeault

! ===6=8===============================================================72 
! ===6=8===============================================================72 
	subroutine cdtur_bougeault(ix,iy,kms,kme,kts,kte,te,z,dz,exch,dlk)
	! compute turbulent coefficients
	implicit none
	integer iz,kms,kme,kts,kte,ix,iy
	real te_m,dlk_m,ck_b
        parameter(ck_b=0.4)
	real te(kms:kme),exch(kms:kme),z(kms:kme)
	real dz(kms:kme)
	real dlk(kms:kme)
	real fact

	 exch(kts)=0.
	do iz=kts+1,kte
	   te_m=(te(iz-1)*dz(iz)+te(iz)*dz(iz-1))/(dz(iz)+dz(iz-1))
	   dlk_m=(dlk(iz-1)*dz(iz)+dlk(iz)*dz(iz-1))/(dz(iz)+dz(iz-1))
	   exch(iz)=ck_b*dlk_m*sqrt(te_m)
           exch(iz)=max(exch(iz),0.1) 
	enddo
	exch(kte+1)=0.1
	return
	end subroutine cdtur_bougeault




