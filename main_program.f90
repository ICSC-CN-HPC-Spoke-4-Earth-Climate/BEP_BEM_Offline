PROGRAM main_program

  USE module_sf_bep_bem
!  USE module_column
  
  IMPLICIT NONE
! ----------------------------------------------------------------------- 
!  Dimension for the array used in the BEP module
! -----------------------------------------------------------------------  
integer nurbmax_mp ! Maximum number of urban classes
      parameter (nurbmax_mp=3)
 integer ndm_mp ! Maximum number of street directions
      parameter (ndm_mp=2)
  integer nz_um_mp ! Maximum number of vertical levels in the urban grid 
     parameter(nz_um_mp=5)
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
      parameter (dz_u_mp=5.)
  integer nbui_max_mp !maximum number of types of buildings in an urban class
      parameter (nbui_max_mp=1)   !must be less or equal than nz_um
!-----------------------------------------------------------------------------------------------  
  INTEGER, PARAMETER :: iz_u_target=1, iw_target=1, nz_um_target=1, nbui_max_target=1   !usati???
  INTEGER, PARAMETER :: ix=1, iy=1

  INTEGER, PARAMETER :: ims = 1, ime = 1	!number of elements along x-axis
  INTEGER, PARAMETER :: kms = 1, kme = 18	!number of elements along y-axis
  INTEGER, PARAMETER :: jms = 1, jme = 1	!number of elements along z-axis
!------------------------------------------------------------------------------------------------
  INTEGER :: i, k, j, z, i_o, file_num, file_num_2d, iii, status
  INTEGER :: ibui, id, iz_um, iz_u, iw, ig, iz, ir

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

  real, dimension(ims:ime, kms:kme, jms:jme) :: arr_3d 			! 3D array to be read from file
  real, dimension(ims:ime, jms:jme) :: arr_2d				! 2D array to be read from file
  real, dimension(47,ims:ime, kms:kme, jms:jme) :: input_3d	        ! 3D array to be assign for each input file
  real, dimension(27,ims:ime, jms:jme) :: input_2d			! 2D arrat to be assign for each input file
  character(len=256) :: file_path, file_path_2d, filename, filename_2d  ! file path
  character(len=256) :: path_output, filename_output			! path output
  character(len=10), dimension(47) :: file_list		 		! list of filenames 3D
  character(len=10), dimension(27) :: file_list_2d			! list of filenames 2D
  
  real :: value, value_0, value_1, value_1b, value_2, value_3, value_4, value_4b, value_5, value_6, value_7, value_99
  real :: value_hi

  ! Inizializzo tutte le variabili scalari contenute nella subroutine BEP_BEM (tot num 35 compresi indici sopra)
  real :: dt=60			!"TEMPORAL RESOLUTION"      "SECONDS"
  real :: gmt=0.0 		
  real :: declin_urb=45.0
  integer :: n_loop = 1
  real :: time_max = 60   !value for module_column
  real :: prtime = 60     !value for module_column

  INTEGER, PARAMETER :: julday=150, ids=1, ide=1, jds=1, jde=1, kds=1, kde=18, its=1, ite=1, jts=1, jte=1, kts=1, kte=18, & 
             itimestep=1, &
  urban_map_zrd=nz_um_mp*nwr_u_mp*ndm_mp,  urban_map_zwd=nbui_max_mp*nz_um_mp*nwr_u_mp*ndm_mp, &
  urban_map_gd=ng_u_mp*ndm_mp, urban_map_zd=nbui_max_mp*nz_um_mp*ndm_mp, urban_map_zdf=nz_um_mp*ndm_mp, &
  urban_map_bd=nbui_max_mp*nz_um_mp, urban_map_wd=nbui_max_mp*nz_um_mp*ndm_mp, urban_map_gbd=nbui_max_mp*ngb_u_mp*ndm_mp, &
  urban_map_fbd=nbui_max_mp*nf_u_mp*nz_um_mp*ndm_mp, urban_map_zgrd=nz_um_mp*nwr_u_mp*ndm_mp, &
  num_urban_hi=1, num_urban_ndm=2
  
  REAL, DIMENSION(ims:ime, 1:urban_map_zrd, jms:jme ) :: trb_urb4d	!%"TRB_URB4D" "ROOF LAYER TEMPERATURE"          "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_zwd, jms:jme ) :: tw1_urb4d	!"TW1_URB4D" "WALL LAYER TEMPERATURE"          "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_zwd, jms:jme ) :: tw2_urb4d	!"TW2_URB4D" "WALL LAYER TEMPERATURE"          "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_gd , jms:jme ) :: tgb_urb4d	!"TGB_URB4D" "ROAD LAYER TEMPERATURE"          "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_zgrd, jms:jme ) :: trv_urb4d	!"TRV_URB4D" "GREEN ROOF LAYER TEMPERATURE"          "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_zgrd, jms:jme ) :: qr_urb4d	!"QR_URB4D" "GREEN ROOF LAYER MOISTURE"          "dimensionless"
  
  REAL, DIMENSION(ims:ime, 1:urban_map_zdf, jms:jme ) :: drain_urb4d	!"DRAIN_URB4D" "GREEN ROOF DRAINAGE"          "mm"
  REAL, DIMENSION(ims:ime, 1:urban_map_bd, jms:jme  ) :: tlev_urb3d	!"TFLEV_URB3D" "FLOOR TEMPERATURE"                       "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_bd , jms:jme ) :: qlev_urb3d	!"QLEV_URB3D" "SPECIFIC HUMIDITY"              "dimensionless"
  REAL, DIMENSION(ims:ime, 1:urban_map_wd , jms:jme ) :: tw1lev_urb3d	!"TW1LEV_URB3D" "WINDOW TEMPERATURE"           "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_wd , jms:jme ) :: tw2lev_urb3d	!"TW2LEV_URB3D" "WINDOW TEMPERATURE"           "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_gbd, jms:jme ) :: tglev_urb3d	!"TGLEV_URB3D" "GROUND TEMPERATURE BELOW A BUILDING"     "K"
  REAL, DIMENSION(ims:ime, 1:urban_map_fbd, jms:jme ) :: tflev_urb3d	!"TFLEV_URB3D" "FLOOR TEMPERATURE"                       "K"

  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ) :: t_pv_urb3d	!"T_PV_URB3D"  "PHOTOVOLTAIC PANELS TEMPERATURE " "K" !PVP
  REAL, DIMENSION( ims:ime, 1:urban_map_wd , jms:jme ) :: sfwin1_urb3d	!"SFWIN1_URB3D" "SENSIBLE HEAT FLUX FROM URBAN SFC WINDOW"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:urban_map_wd , jms:jme ) :: sfwin2_urb3d	!"SFWIN2_URB3D" "SENSIBLE HEAT FLUX FROM URBAN SFC WINDOW"  "W m{-2}"

  REAL, DIMENSION( ims:ime, 1:urban_map_zd , jms:jme ) :: sfw1_urb3d	!"SFW1_URB3D"  "SENSIBLE HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:urban_map_zd , jms:jme ) :: sfw2_urb3d	!"SFW2_URB3D"  "SENSIBLE HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ) :: sfr_urb3d	!"SFR_URB3D"  "SENSIBLE HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:num_urban_ndm, jms:jme ) :: sfg_urb3d	!"SFG_URB3D"  "SENSIBLE HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ) :: sfrv_urb3d	!"SFRV_URB3D"  "SENSIBLE HEAT FLUX FROM GREEN ROOF"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ) :: lfrv_urb3d

  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ) :: dgr_urb3d !GRZ"	DGR_URB3D" "ROOF LAYER DEPTH WATER RETENTION"          "mm"
  REAL, DIMENSION( ims:ime, 1:num_urban_ndm, jms:jme ) :: dg_urb3d !GRZ		"DG_URB4D" "ROOF LAYER DEPTH WATER RETENTION"          "mm"
  REAL, DIMENSION( ims:ime, 1:urban_map_zdf, jms:jme ) :: lfr_urb3d !GRZ	"LFR_URB3D"  "LATENT HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  REAL, DIMENSION( ims:ime, 1:num_urban_ndm, jms:jme ) :: lfg_urb3d !G	"LFG_URB3D"  "LATENT HEAT FLUX FROM URBAN SFC"  "W m{-2}"
  
  REAL, DIMENSION( ims:ime, 1:num_urban_hi, jms:jme ) :: hi_urb2d	!"HEIGHT_HISTOGRAMS" "DISTRIBUTION OF BUILDING HEIGHTS" "dimensionless"

! we will solve the drag and vertical diffusion in a column for neutral case
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +    by A. Martilli,     CIEMAT  SP 2040 MADRID                  +
!     +                        phone: ++34-9-13-46-62-99               +
!     +                        email:alberto.martilli@ciemat.es        +
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           integer nzm
           parameter (nzm=18)
           integer jz
           real, dimension(1) :: lambda_p	!da cambiare con dimensioni
           real, dimension(nzm+1) ::  hgt_column
           real, dimension(nzm+1) ::  prb_column

! flow variables
           real vx(nzm) ! x component of the wind
           real vy(nzm) ! y component of the wind
           real sstot
!AM: this is for wind only becasue we have a neutral case. I add below
!potential temperature, even if no sources are present in this column
!model. You will have to add them based on BEP_BEM. 
           real pt(nzm) ! potential temperature 
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
           real ceps,ck,cmu ! constants for the k-l scheme
	   real pr ! prandtl number
           parameter (ceps=1/1.4,ck=0.4,cmu=0.09)
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
           real vl_c(nzm) ! fraction of air in each cell
           real sf_c(nzm+1) ! fraction of air at the interface between cells
! dispersive flux
           real dwtrc(nzm+1)
	   real duw(nzm+1)
           real aaa(4,4),bbb(4)

! drag coefficient
           real cdragx(nzm),cdragy(nzm) ! drag coefficient in x and y direction           
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
           real dz ! vertical resolution
           integer nz ! number of vertical levels
           real z_c(nzm+1) ! height of the faces of the numerical grid
           !real dt ! time step
           integer iz_c ! loop index for vertical
           integer it ! loop index for time
           real time ! time from the start
           real time_pr ! time from the last print
          ! real time_max ! total time of the simulation (definiti all'inizio)
           integer ntime ! total number of timesteps
          ! real prtime ! frequency of output            (definiti all'inzio)
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
           
! flag to distinguish between alligned and staggered arrays (1: staggered, 2: alligned)
           integer iconfig
! working
           real lambdaf_x,lambdaf_y
           real tke_old(nzm) ! tke
           real dtot,wrk(8)
           real trc_an(nzm),trc_2(nzm)
           integer i_col

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ASSIGN A VALUE TO ALL ELEMENTS OF THE 3D ARRAY
value = 290.00;
value_0 = 0.1;
value_99 = 0.0;
value_hi = 1;

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
        dlg(i, j, k) = value_99
        dl_u(i, j, k) = value_99
        sf(i, j, k) = value_99
        vl(i, j, k) = value_99
       END DO
     END DO
  END DO

DO i=ims,ime
 DO j = 1, num_urban_hi
    DO k = jms,jme
      hi_urb2d(i,j,k)=value_hi
    END DO
  END DO
END DO

DO i = ims, ime
    DO j = 1, urban_map_zrd   !360
      DO k = jms,jme 
        trb_urb4d(i, j, k) = value
        qr_urb4d(i, j, k) = value_0
        trv_urb4d(i, j, k) = value
      END DO
    END DO
  END DO
 
 value_1 = 290.00;
 value_1b = 100.1;
DO i = ims, ime
    DO j = 1, urban_map_zd	!540
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
    DO j = 1, urban_map_zwd    !5400
      DO k = jms,jme
        tw1_urb4d(i, j, k) = value_2
        tw2_urb4d(i, j, k) = value_2
      END DO
    END DO
  END DO

 value_3 = 290.00;

DO i = ims, ime
    DO j = 1, urban_map_gd	!20
      DO k = jms,jme
        tgb_urb4d(i, j, k) = value_3
      END DO
    END DO
  END DO

 value_4 = 290.00;
 value_4b = 0.01 
DO i = ims, ime
    DO j = 1, urban_map_zdf      !36
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
    DO j = 1, urban_map_bd      !270
      DO k = jms,jme
        tlev_urb3d(i, j, k) = value_5	!"INDOOR TEMPERATURE"
        qlev_urb3d(i, j, k) = value_0	!"SPECIFIC HUMIDITY"
      END DO
    END DO
  END DO
 
  value_6 = 290.00;

DO i = ims, ime
    DO j = 1, urban_map_gbd      !300
      DO k = jms,jme
        tglev_urb3d(i, j, k) = value_6	!GROUND TEMPERATURE BELOW A BUILDING
      END DO
    END DO
  END DO

  value_7 = 290.00;

DO i = ims, ime
    DO j = 1, urban_map_fbd      !5100
      DO k = jms,jme
        tflev_urb3d(i, j, k) = value_7	!"FLOOR TEMPERATURE"
      END DO
    END DO
  END DO

DO i = ims, ime
    DO j = 1, num_urban_ndm      !2
      DO k = jms,jme
        sfg_urb3d(i, j, k) = 100.1  	!SENSIBLE HEAT FLUX FROM URBAN SFC
        dg_urb3d(i, j, k) = 0.01	!"ROOF LAYER DEPTH WATER RETENTION"          "mm"
        lfg_urb3d(i, j, k) = 100.1	!LATENT HEAT FLUX FROM URBAN SFC"
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
    
! Define the common path to the input files
  file_path = "/home/gpappaccogli/BEP_BEM/bep_bem_may/input_file/Output_3d/"
  file_path_2d = "/home/gpappaccogli/BEP_BEM/bep_bem_may/input_file/Output_2d/"
  path_output = "/home/gpappaccogli/BEP_BEM/Output/"

  ! Define the list of filenames
  file_list = (/"t_01", "t_02", "t_03", "t_04", "t_05", "t_06", "t_07", "t_08", "t_09", & 
                "t_10", "t_11", "t_12", "t_13", "t_14", "t_15", "t_16", "t_17", & 
		"t_18", "t_19", "t_20", "t_21", "t_22", "t_23", "t_24", "t_25", & 
		"t_26", "t_27", "t_28", "t_29", "t_30", "t_31", "t_32", "t_33", & 
		"t_34", "t_35", "t_36", "t_37", "t_38", "t_39", "t_40", "t_41", & 
		"t_42", "t_43", "t_44", "t_45", "t_46", "t_47" /)

  file_list_2d = (/ "d_01", "d_02", "d_03", "d_04", "d_05", "d_06", "d_07", "d_08", "d_09", &
                    "d_10", "d_11", "d_12", "d_13", "d_14", "d_15", "d_16", "d_17", &
                    "d_18", "d_19", "d_20", "d_21", "d_22", "d_23", "d_24", "d_25", &
                    "d_26", "d_27" /)

  ! Read the 3D array from each file with the same dimension
  do file_num = 1, 47
    ! Define the filename for the current file
    filename = trim(file_path)//trim(file_list(file_num))//".txt"
    
    ! Open the current input file
    open(unit=10, file=filename, status="old")
    
    ! Read the 3D array from the current file
  do j = 1, jme
    do k = 1, kme
       read(10,*) (arr_3d(i,k,j), i=1,ime)
   end do
 end do
    
   ! Assign the values of the current input array to the appropriate element of input array
   input_3d(file_num,:,:,:) = arr_3d(:,:,:)

  ! Close the current input file
  close(10)
end do

  ! Read the 2D array from each file with the same dimension
  do file_num_2d = 1, 27
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

!----------------------------------------------------------------------------
!APRO I FILE PER SALVARE VARIBILI 2D & 3D
!----------------------------------------------------------------------------
OPEN(501, FILE='QGR_URB3D.txt')
OPEN(502, FILE='TGR_URB3D.txt')
OPEN(503, FILE='DRAINGR_URB3D.txt')
OPEN(504, FILE='LF_AC_URB3D.txt')
OPEN(505, FILE='SF_AC_URB3D.txt')
OPEN(506, FILE='CM_AC_URB3D.txt')
OPEN(507, FILE='EP_PV_URB3D.txt')
OPEN(508, FILE='SFVENT_URB3D.txt')
OPEN(509, FILE='LFVENT_URB3D.txt')

OPEN(20, FILE='TW1_URB4D.txt')
OPEN(21, FILE='TW2_URB4D.txt')

OPEN(32, FILE='TRB_URB4D.txt')
OPEN(33, FILE='TRV_URB4D.txt')
OPEN(34, FILE='QR_URB4D.txt')
OPEN(35, FILE='TGB_URB4D.txt')

OPEN(41, FILE='DRAIN_URB4D.txt')
OPEN(42, FILE= 'T_PV_URB3D.txt')
OPEN(43, FILE= 'SFR_URB3D.txt')
OPEN(44, FILE= 'SFRV_URB3D.txt')
OPEN(45, FILE='LFRV_URB3D.txt')
OPEN(46, FILE='DGR_URB3D.txt')
OPEN(47, FILE='LFR_URB3D.txt')

OPEN(50, FILE='TLEV_URB3D.txt')
OPEN(51, FILE='QLEV_URB3D.txt')

OPEN(112, FILE='TW1LEV_URB3D.txt')
OPEN(113, FILE='TW2LEV_URB3D.txt')
OPEN(114, FILE='SFWIN1_URB3D.txt')
OPEN(115, FILE='SFWIN2_URB3D.txt')
OPEN(116, FILE='SFW1_URB3D.txt')
OPEN(117, FILE='SFW2_URB3D.txt')

OPEN(200, FILE='TGLEV_URB3D.txt')
OPEN(201, FILE='TFLEV_URB3D.txt')
!WRITE(*,*) "a_u", input_3d(4,:,:,:)

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
!    WRITE(*,*) "esempio ind_wd loop", ind_wd_mp
!    WRITE(*,*) "id", id
!    WRITE(*,*) "iz_u", iz_u
!    WRITE(*,*) "ibui", ibui
!    WRITE(*,*) "iii", iii
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
          jz=1
          ss=0.
          pb=0.
          hgt = -100
          prb = -100
          hgt_column = -100
          prb_column = -100. 
          lambda_p = 0.
           write(*,*) "hgt_column", hgt_column

          ! convert urban parameter to column
          do i = 1, jme
    		hgt_column(i) = input_2d(23,1,i) ! Replace '1' with the desired column index (modificare in seguito)
                write(*,*) "hgt_column  in loop", hgt_column
          end do
          do i = 1, num_urban_hi      
                prb_column(i) = hi_urb2d(1,i,1) ! Replace '1' with the desired column index (modificare in seguito)
                write(*,*) "prb_column in loop", prb_column
          end do
 !!!         do i = 1, jme
 !!!               hmean(i) = hmean(,1,i) ! Replace '1' with the desired column index (modificare in seguito)
 !!!               write(*,*) "hgt_column", hgt_column
 !!!         end do

          ! define hgt and prb with hgt_column and prb_column
          hgt = hgt_column    !_gl dato da hgt_urb2d
          prb = prb_column    !_gl dato da hi_urb2d

          ! define the mesh
          do iz_c=1,nzm+1
             z_c(iz)=real(iz_c-1)*dz
          enddo

          ! define probabilities
           do jz=1,nzm+1
               if(hgt(jz).ge.0.)then
                 do iz_c=1,nzm+1 
                   if(z_c(iz_c).eq.hgt(jz))then
                     ss(iz_c)=prb(jz)
                   endif
                 enddo
             endif
            enddo
            do iz_c=1,nzm+1
               sstot=0.
               do jz=iz_c,nzm+1
                  sstot=sstot+ss(jz)
               enddo
               pb(iz_c)=sstot
            enddo

          ! compute fraction of air for each grid cell
          do i = 1, jme
               lambda_p(i) = input_2d(21,1,i) ! Replace '1' with the desired column index (modificare in seguito)
                write(*,*) "lambda_p", lambda_p
          end do
!!             write(*,*) "size of pb", size(pb)
             !write(*,*) "size of lambda_p", size(lambda_p)   !non funziona dimensione
         !!!    write(*,*) "size of vl_c", size(vl_c)
         !!!    write(*,*) "size of sf_c", size(sf_c)
         !!!    write(*,*) "size of ss", size(ss)
         !!!    write(*,*) "write value of nzm", nzm
         !!!    write(*,*) "write value of iz_c", iz_c

             !lambda_p=bx*by/((wx+bx)*(wy+by))   !AM  (LO PRENDIAMO COME INPUT)
             do iz_c=1,nzm		               !AM
               vl_c(iz_c)=1.-lambda_p(1)*pb(iz_c+1)    !AM
               sf_c(iz_c)=1.-lambda_p(1)*ss(iz_c)      !AM
             enddo                                     !AM
            sf_c(nzm+1)=1.
          
            !compute mean height
            hmean=0.
            do iz_c=1,nz
               hmean=hmean+z_c(iz_c)*ss(iz_c)
            enddo

! pressure gradient
            dpdx=(utau_x**2.)/z_c(nz)
            dpdy=(utau_y**2.)/z_c(nz)
            ntime=time_max/dt

       !------------------------------------------------------------
       ! end ex subroutine read_input
       !----------------------------------------------------------
          iconfig=1
!AM the two options below refers to the paper of Andres SImon-Moral
!(BLM). For your case, I would keep staggered, but you can still make
!some tests with alligned if you like to see the sensitivity of the
!results.
	 if(iconfig.eq.1)write(*,*)'staggered'
	 if(iconfig.eq.2)write(*,*)'alligned'
! compute the drag coefficent and the length scales
           if(iconfig.eq.1)then
            call drag_length_stag(nzm,nz,z_c,wx,wy,bx,by,hmean,ceps,ck,cmu,cdragx,cdragy,dls_u,dlk_u,lambda_p(1))
           elseif(iconfig.eq.2)then
            call drag_length_alligned(nzm,nz,z_c,wx,wy,bx,by,hmean,ceps,ck,cmu,cdragx,cdragy,dls_u,dlk_u,lambda_p(1))
	 endif
	 write(*,*)'cdragx',cdragx(1)
! initialize vx,vy,tke - just put some value, the code should forget it since we impose a pressure gradient
           lambdaf_x=hmean*by/((wx+bx)*(wy+by))
           lambdaf_y=hmean*bx/((wx+bx)*(wy+by))
!AM below is just to have an initial value different than zero for the
!different variables. In your case I would put vx,vy,pt, and tke equal to
!those of the weather model at initialization.
           trc=0.
           cdz=0.1
           emi=1.
           vx=10
           tke=2.
           vy=0.
 
           time=0.
           time_pr=0.
           write(*,*)'prtime',prtime
           open(unit=66,file='output')

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INIZIO CICLO NEL TEMPO
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

do i_o =1, n_loop		
!AM here the loop over time starts. I would introduce in the loop a
!reading of the weather model variables (wind components, temperature,
!etc., every hour or every time they are available).
!put the values of the weather model at the top of the column, like:
!vx(nz)=U(weather model) 
!vy(nz)=V(weather model) 
!pt(nz)=PT(weather model) 
!tke(nz)=TKE(weather model) if you have TKE form the weather model
!in tbis exampel all these variables are kept equal to the initial value
!at the top of the colum during the whole simulation.

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

             time=time+dt
             time_pr=time_pr+dt
!AM how to estimate the lenght scales will be the main challenge here
!I called dlk_u and dls_u those determined just by the morphology, that
!are unchanged during the run, but we should adapt them based on the
!stability and PBLH.
            dlk=dlk_u
            dls=dls_u
!as first guess we may try the following:
!
!            dlk=dlk_u
!            dls=sqrt(dls_u*PBLH)
!with PBLH fomr the weather model
!but further adaptation  may be done based on atmospheric stability, as
!suggested by Andrea (phy functions?). 

! first compute the diffusion coefficients

             call cdtur(nzm,nz,dz,ck,tke,dlk,cdz)

            write(*,*) "srex_vx", srex_vx, size(srex_vx)
            write(*,*) "srim_vx", srim_vx, size(srim_vx)
            write(*,*) "srex_vy", srex_vy, size(srex_vy)
            write(*,*) "srim_vy", srim_vy, size(srim_vy)

            write(*,*) "srex_pt", srex_pt, size(srex_pt)
            write(*,*) "srim_pt", srim_pt, size(srim_pt)
            write(*,*) "srex_tke", srex_tke, size(srex_tke)
            write(*,*) "srim_tke", srim_tke, size(srim_tke)

! compute the implicit and explicit terms of the equation for vx,vy and tke

!AM Here I would intorduce the call to BEP and BEM.
! mind that the correspondence is as follow:
!            b_u = srex_vx
!            a_u = srim_vx
!            b_v = srex_vy
!            a_v = srim_vy
!            b_t = srex_pt
!            a_t = srim_pt
!            b_e = srex_tke
!            a_e = srim_tke
do i = 1, nzm
    b_u(1,i,1) = srex_vx(i)
end do


! Since you will have BEP_BEM you should not call building below

   CALL BEP_BEM(input_2d(1,:,:),UTYPE_URB2D,i,input_3d(1,:,:,:),dt,input_3d(2,:,:,:),input_3d(3,:,:,:),  &	!1
              input_3d(4,:,:,:),input_3d(5,:,:,:),input_3d(6,:,:,:),input_2d(3,:,:),input_2d(4,:,:), 	 & !2
              gmt,julday,input_2d(5,:,:),input_2d(6,:,:),       	   & !3
              declin_urb,input_2d(7,:,:),input_2d(8,:,:),		   & !4
              num_urban_ndm,  urban_map_zrd,  urban_map_zwd, urban_map_gd, & !5
              urban_map_zd,  urban_map_zdf,   urban_map_bd, urban_map_wd,  & !6
	      urban_map_gbd,  urban_map_fbd,                               & !7
	      urban_map_zgrd,  num_urban_hi,				   & !8
	      trb_urb4d,tw1_urb4d,tw2_urb4d,tgb_urb4d,		           & !9
              tlev_urb3d,qlev_urb3d,tw1lev_urb3d,tw2lev_urb3d,	           & !10
	      tglev_urb3d,tflev_urb3d,input_2d(9,:,:),input_2d(10,:,:),    & !11
              input_2d(11,:,:), & !12
	      input_2d(12,:,:),input_2d(13,:,:), & !13
              sfwin1_urb3d,sfwin2_urb3d, &  !14
              sfw1_urb3d,sfw2_urb3d,sfr_urb3d,sfg_urb3d, & !15
              input_2d(14,:,:),t_pv_urb3d, &  !16
              trv_urb4d,qr_urb4d,input_2d(15,:,:),input_2d(16,:,:), & !17
              drain_urb4d,input_2d(17,:,:), &  !18
              sfrv_urb3d,lfrv_urb3d, &  !19
              dgr_urb3d,dg_urb3d, &   !20
              lfr_urb3d,lfg_urb3d,input_2d(18,:,:),input_2d(19,:,:),input_2d(20,:,:), &  !21
              input_2d(21,:,:),hi_urb2d,input_2d(22,:,:),input_2d(23,:,:), & !22
              srim_vx,srim_vy,srim_pt,srim_tke,srex_vx,srex_vy, & !23
              srex_pt,srex_tke,b_q,dlg,dl_u,sf,vl, & !24
	      input_2d(24,:,:),input_2d(25,:,:),input_2d(26,:,:),input_2d(27,:,:),input_3d(47,:,:,:), & !25
              ids,ide,jds,jde,kds,kde, & !26
              ims,ime,jms,jme,kms,kme, & !27
              its,ite,jts,jte,kts,kte) !28

            write(*,*) "main_program call tke_bougeault"

            write(*,*) "srex_vx dp_bep", srex_vx, size(srex_vx)
            write(*,*) "srim_vx dp_bep", srim_vx, size(srim_vx)
            write(*,*) "srex_vy dp_bep", srex_vy, size(srex_vy)
            write(*,*) "srim_vy dp_bep", srim_vy, size(srim_vy)

            write(*,*) "srex_pt dp_bep", srex_pt, size(srex_pt)
            write(*,*) "srim_pt dp_bep", srim_pt, size(srim_pt)
            write(*,*) "srex_tke dp_bep", srex_tke, size(srex_tke)
            write(*,*) "srim_tke dp_bep", srim_tke, size(srim_tke)

call tke_bougeault(nzm,nz,ceps,dz,vx,vy, &
                   tke,cdz,dls,sf_c, &                      
                   srim_tke,srex_tke)

 write(*,*) "main_program dopo call tke_bougeault"

!    write(*,*)'srex_vx',srex_vx

! add the pressure gradient
!AM no need to have the lines below, since we will force with the 
             srex_vx=srex_vx+dpdx !AM comemnt this
             srex_vy=srex_vy+dpdy !AM comment this
           !  srex_vx(nz)=utau_x*utau_x

             srex_trc(1)=emi/dz/vl_c(1)
             srex_trc(nz)=-emi/dz
 write(*,*) "main_program level 1"
     !         write(*,*)'dtot=',dtot
             pr=1. !AM this is the Prandtl number. Currently equal to 1, but we may adapt it based on stabilility as suggested by Andrea.
             cdt=pr*cdz
! compute the vertical diffusion !AM below I change to force at the top
           ! call diff(nzm,nz,1,1,dt,vx,cdz,srim_vx,srex_vx,sf,vl,dz,uw,duwdz)
           ! call diff(nzm,nz,1,1,dt,vy,cdz,srim_vy,srex_vy,sf,vl,dz,vw,dvwdz)
           ! call diff(nzm,nz,1,1,dt,tke,cdz,srim_tke,srex_tke,sf,vl,dz,wtke,dwtkedz)
           ! call diff(nzm,nz,1,2,dt,trc,cdt,srim_trc,srex_trc,sf,vl,dz,wtrc,dwtrcdz)

           write(*,*) "main_program call diff subroutine"            
            call diff(nzm,nz,1,2,dt,vx,cdz,srim_vx,srex_vx,sf_c,vl_c,dz,uw,duwdz)
            call diff(nzm,nz,1,2,dt,vy,cdz,srim_vy,srex_vy,sf_c,vl_c,dz,vw,dvwdz)
            call diff(nzm,nz,1,2,dt,pt,cdt,srim_pt,srex_pt,sf_c,vl_c,dz,wt,dwtdz)
            call diff(nzm,nz,1,2,dt,tke,cdz,srim_tke,srex_tke,sf_c,vl_c,dz,wtke,dwtkedz)

!AM here abvoe I made the assumption that you will have TKE fomr the
!weather model. If not, you can keep the top b.c. open, and call as
!below:
!            call diff(nzm,nz,1,1,dt,tke,cdz,srim_tke,srex_tke,sf,vl,dz,wtke,dwtkedz)
 
            call diff(nzm,nz,1,2,dt,trc,cdt,srim_trc,srex_trc,sf_c,vl_c,dz,wtrc,dwtrcdz)
        !    write(*,*)time/3600.,trc(1)

              if(time_pr.gt.prtime)then
              !adapt your prinitng as you like
               time_pr=dt
               write(*,*)'time=',time,'time_pr=',time_pr
               write(66,*)'time=',time,'time_pr=',time_pr
               do iz_c=1,nz
                uwm=(uw(iz_c)+uw(iz_c+1))/2.
                write(66,*)iz_c,vx(iz_c)/utau_x,tke(iz_c)/(utau_x**2.),uwm/(utau_x**2.)
               enddo
          
              endif
!!!           enddo  !_gl (chiusura loop ntime column_alberto)

           
!           end   
                                       
!   end do
                    
!--------------------------------------------------------------------------
! SEZIONE PER STAMPARE LE VARIABILI - 2D
!-------------------------------------------------------------------------
! WRITE(*,*) "qgr_urb3d:", input_2d(15,:,:)
! WRITE(*,*) "tgr_urb3d:", input_2d(16,:,:)
! WRITE(*,*) "draingr_urb3d:", input_2d(17,:,:)
! WRITE(*,*) "lf_ac_urb3d:", input_2d(10,:,:)
! WRITE(*,*) "sf_ac_urb3d:", input_2d(9,:,:)
! WRITE(*,*) "cm_ac_urb3d:", input_2d(11,:,:)
! WRITE(*,*) "ep_pv_urb3d:", input_2d(14,:,:)
! WRITE(*,*) "sfvent_urb3d:", input_2d(12,:,:)
! WRITE(*,*) "lfvent_urb3d:", input_2d(13,:,:)
! WRITE(*,*) "a_u:", input_3d(34,:,:,:)   

!------------------------------------------------------------------------------
! SCRITTURA DEI DATI IN BASE AGLI INDICI
!-----------------------------------------------------------------------------
OPEN(501, FILE='QGR_urb3D.txt')
OPEN(502, FILE='TGR_urb3D.txt')
OPEN(503, FILE='DRAINGR_urb3D.txt')
OPEN(504, FILE='LF_AC_urb3D.txt')
OPEN(505, FILE='SF_AC_urb3D.txt')
OPEN(506, FILE='CM_AC_urb3D.txt')
OPEN(507, FILE='EP_PV_urb3D.txt')
OPEN(508, FILE='SFVENT_urb3D.txt')
OPEN(509, FILE='LFVENT_urb3D.txt')

write(501,'(24(1x,f10.3))') input_2d(15,:,:)
write(502,'(24(1x,f10.3))') input_2d(16,:,:)
write(503,'(24(1x,f10.3))') input_2d(17,:,:)
write(504,'(24(1x,f10.3))') input_2d(10,:,:)
write(505,'(24(1x,f10.3))') input_2d(9,:,:)
write(506,'(24(1x,f10.3))') input_2d(11,:,:)
write(507,'(24(1x,f10.3))') input_2d(14,:,:)
write(508,'(24(1x,f10.3))') input_2d(12,:,:)
write(509,'(24(1x,f10.3))') input_2d(13,:,:)

id=1

do iz=1,nz_um_mp
   if(iz.gt.0)then !un livello verticale superiore a quello indicato 
     write(42,'(24(1x,f10.3))') t_pv_urb3d(ix,ind_zdf_mp(iz,id),iy), t_pv_urb3d(ix,ind_zdf_mp(iz,id+1),iy) !"PHOTOVOLTAIC PANELS TEMPERATURE " "K"
     write(43,'(24(1x,f10.3))')	sfr_urb3d(ix,ind_zdf_mp(iz,id),iy), sfr_urb3d(ix,ind_zdf_mp(iz,id+1),iy)   !"SENSIBLE HEAT FLUX FROM URBAN SFC"
     write(46,'(24(1x,f10.3))')	dgr_urb3d(ix,ind_zdf_mp(iz,id),iy), dgr_urb3d(ix,ind_zdf_mp(iz,id+1),iy)   !"ROOF LAYER DEPTH WATER RETENTION" mm  
     write(47,'(24(1x,f10.3))')	lfr_urb3d(ix,ind_zdf_mp(iz,id),iy), lfr_urb3d(ix,ind_zdf_mp(iz,id+1),iy)   !"LATENT HEAT FLUX FROM URBAN SFC"
!	if(gr_flag_u.eq.1)then
     write(44,'(24(1x,f10.3))')	sfrv_urb3d(ix,ind_zdf_mp(iz,id),iy), sfrv_urb3d(ix,ind_zdf_mp(iz,id+1),iy) !SENSIBLE HEAT FLUX FROM GREEN ROOF
     write(45,'(24(1x,f10.3))')	lfrv_urb3d(ix,ind_zdf_mp(iz,id),iy), lfrv_urb3d(ix,ind_zdf_mp(iz,id+1),iy) !"LATENT HEAT FLUX FROM GREEN ROOF"
     write(41,'(24(1x,f10.3))')	drain_urb4d(ix,ind_zdf_mp(iz,id),iy), drain_urb4d(ix,ind_zdf_mp(iz,id+1),iy) !"GREEN ROOF DRAINAGE"   "mm"
      endif
   enddo

do ibui=1,nbui_max_mp
   do iz_u=1,nz_um_mp
     do iw=1,nwr_u_mp
        if(ibui.eq.1.and.iz_u.eq.1.and.iw.eq.1)then
        write(20,'(24(1x,f10.3))') tw1_urb4d(ix,ind_zwd_mp(ibui,iz_u,iw,id),iy), tw1_urb4d(ix,ind_zwd_mp(ibui,iz_u,iw,id+1),iy) !"WALL LAYER TEMPERATURE"
        write(21,'(24(1x,f10.3))') tw2_urb4d(ix,ind_zwd_mp(ibui,iz_u,iw,id),iy), tw2_urb4d(ix,ind_zwd_mp(ibui,iz_u,iw,id+1),iy)
       endif
      enddo
    enddo
  enddo

do ig=1,ng_u_mp
   if(ig.eq.1)then
   write(35,'(24(1x,f10.3))') tgb_urb4d(ix,ind_gd_mp(ig,id),iy), tgb_urb4d(ix,ind_gd_mp(ig,id+1),iy) !"ROAD LAYER TEMPERATURE"
   endif
enddo

!do id=1,ndm   !gl
do iz_u=1,nz_um_mp
 do ir=1,nwr_u_mp
    if(iz_u.eq.1+1.and.ir.eq.1)then !un livello nz_um  superiore a quello indagato
!    write(*,*) "write index of trb_urb4d", iz_u, ir, id
    write(32,'(24(1x,f10.3))') trb_urb4d(ix,ind_zrd_mp(iz_u,ir,id),iy), trb_urb4d(ix,ind_zrd_mp(iz_u,ir,id+1),iy) !"ROOF LAYER TEMPERATURE"
    write(33,'(24(1x,f10.3))') trv_urb4d(ix,ind_grd_mp(iz_u,ir,id),iy), trv_urb4d(ix,ind_grd_mp(iz_u,ir,id+1),iy) !"GREEN ROOF LAYER TEMPERATURE"
    write(34,'(24(1x,f10.3))') qr_urb4d(ix,ind_grd_mp(iz_u,ir,id),iy), qr_urb4d(ix,ind_grd_mp(iz_u,ir,id+1),iy) !"GREEN ROOF LAYER MOISTURE"
    endif
   enddo
 enddo

!Outputs of BEM
do ibui=1,nbui_max_mp !type of building
  ! do iz_u=1,nz_um_mp !vertical levels
      if(ibui.eq.1.and.iz_u.gt.1)then     !extract data of five  nz_um and ibui=1 
      write(50,'(24(1x,f10.3))') tlev_urb3d(ix,ind_bd_mp(ibui,1),iy),tlev_urb3d(ix,ind_bd_mp(ibui,2),iy), &
   tlev_urb3d(ix,ind_bd_mp(ibui,3),iy),tlev_urb3d(ix,ind_bd_mp(ibui,4),iy),tlev_urb3d(ix,ind_bd_mp(ibui,5),iy)
    write(51, '(24(1x,f10.3))') qlev_urb3d(ix,ind_bd_mp(ibui,1),iy),qlev_urb3d(ix,ind_bd_mp(ibui,2),iy), &
   qlev_urb3d(ix,ind_bd_mp(ibui,3),iy),qlev_urb3d(ix,ind_bd_mp(ibui,4),iy),qlev_urb3d(ix,ind_bd_mp(ibui,5),iy)
      endif  
  !  enddo
enddo

do ibui=1,nbui_max_mp !type of building
  do iz_u=1,nz_um_mp !vertical levels
    if(ibui.eq.1.and.iz_u.eq.1)then  	!extract data from nz_um=1 ibui=1 on both street directions
 	!write(*,*) "nuovo loop", ibui,iz_u,id 
  	!write(*,*) "ultimo indice", ind_wd_mp(ibui,iz_u,id)
  	write(112,'(24(1x,f10.3))') tw1lev_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), tw1lev_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1) !"WINDOW TEMPERATURE"
        write(113,'(24(1x,f10.3))') tw2lev_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), tw2lev_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1)
        write(114,'(24(1x,f10.3))') sfwin1_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), sfwin1_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1) !"SENSIBLE HEAT FLUX FROM URBAN SFC WINDOW"
        write(115,'(24(1x,f10.3))') sfwin2_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), sfwin2_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1)
        write(116,'(24(1x,f10.3))') sfw1_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), sfw1_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1) !SENSIBLE HEAT FLUX FROM URBAN SFC
        write(117,'(24(1x,f10.3))') sfw2_urb3d(1,ind_wd_mp(ibui,iz_u,id),1), sfw2_urb3d(1,ind_wd_mp(ibui,iz_u,id+1),1)
  endif
 enddo
enddo

do ibui=1,nbui_max_mp  !type of building
  do iw=1,ngb_u_mp !layers in the walls
     if(ibui.eq.1.and.iw.eq.1)then   !Ground temperature below a building in BEM
       write(200,'(24(1x,f10.3))') tglev_urb3d(ix,ind_gbd_mp(ibui,iw,id),iy), tglev_urb3d(ix,ind_gbd_mp(ibui,iw,id+1),iy) 
     endif
  enddo
enddo

do ibui=1,nbui_max_mp !type of building
  do iw=1,nf_u_mp !layer in the walls
   do iz_u=1,nz_um_mp-1 !verticals levels
    if(ibui.eq.1.and.iw.eq.10.and.iz_u.eq.1)then   !FLOOR TEMPERATURE
    write(201,'(24(1x,f10.3))') tflev_urb3d(ix,ind_fbd_mp(ibui,iw,iz_u,id),iy), tflev_urb3d(ix,ind_fbd_mp(ibui,iw,iz_u,id+1),iy)
    endif
    enddo
   enddo
  enddo

enddo 	!loop su subroutine bep_bem

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
              write(68,'(i3,11(2x,f15.8))')iz_c,vx(iz_c)/utau_x,tke(iz_c)/(utau_x**2.), & 
                    uwm/(utau_x**2.),trc(iz_c),wtrc(iz_c),cdz(iz_c),dls(iz_c)
            enddo

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!fine ciclo subroutine BEP_BEM
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

      subroutine drag_length_stag(nzm,nz,z,wx,wy,bx,by,hmean,ceps,ck,cmu,cdragx,cdragy,dls,dlk,lambda_p_dls)
      !implicit none
     ! input
      integer :: nzm,nz,iz_c
      real :: z(nzm+1)
      real dls(nzm),dlk(nzm)      
      real wx,wy,bx,by,hmean
      real ceps,ck,cmu
      real, intent(in) :: lambda_p_dls
! output
      real cdragx(nzm),cdragy(nzm)
!local
      real disp,a1,a2,zc,d2      
! Here I follow the formualtion for staggered arrays based on BLM2010. To be modified for aligned

! compute the displacement height
      !lambda_p_dls=bx*by/((wx+bx)*(wy+by))
      disp=hmean*(lambda_p_dls)**(0.13)

! Drag coefficient
      if(lambda_p_dls.le.0.29)then
       do iz_c=1,nz
        cdragx(iz_c)=3.32*(lambda_p_dls)**(0.47)
        cdragy(iz_c)=3.32*(lambda_p_dls)**(0.47)
       enddo
      else
       do iz_c=1,nz
        cdragx(iz_c)=1.85
        cdragy(iz_c)=1.85
       enddo
      endif
!
      cdragx=cdragx*1.7
!

!	  write(*,*)cdragx(1),(lambda_p)**(0.47)
! length scales
      a1=2.19
      a2=1.2
      do iz_c=1,nz
       zc=(z(iz_c)+z(iz_c+1))/2.
       if((zc/hmean).le.1.)then
        dls(iz_c)=ceps*a1*(hmean-disp)
       elseif((zc/hmean).gt.1..and.(zc/hmean).le.1.5)then
        dls(iz_c)=ceps*a1*(zc-disp)
        
       elseif((zc/hmean).gt.1.5)then
        d2=(1.-a1/a2)*1.5*hmean+a1/a2*disp
        dls(iz_c)=ceps*a2*(zc-d2)
       endif
       dlk(iz_c)=cmu/(ceps*ck)*dls(iz_c)
      enddo
      
      return
      end

!----------------------------------------
! ===6================================================================72
!-----------------------------------------------------------------
      subroutine drag_length_alligned(nzm,nz,z,wx,wy,bx,by,hmean,ceps,ck,cmu,cdragx,cdragy,dls,dlk,lambda_p_dla)
      !implicit none
! input
      integer nzm,nz,iz_c
      real z(nzm+1)
      real dls(nzm),dlk(nzm)      
      real wx,wy,bx,by,hmean
      real ceps,ck,cmu
! output
      real cdragx(nzm),cdragy(nzm)
!local
      real disp,lambda_p_dla,a1,a2,zc,d2,lambda_c,lambda_s,a,b,c,f,i_all,j,fs,fc,fcs      
! Here I follow the formualtion for alligned arrays 

! compute the displacement height
      lambda_p_dla=bx*by/((wx+bx)*(wy+by))
      disp=hmean*(lambda_p_dla)**(0.13)
! define the sheltering and channeling paramter
          lambda_c=wy/by
	  lambda_s=wx/hmean
	  write(*,*)'lambda_c=',lambda_c,'lambda_s',lambda_s
          a=0.24
	  b=1.67
	  c=2.07
	  f=0.6
	  i_all=1.4
	  j=4.
          fs=1.-exp(-a*(lambda_s**b))
	  fc=c/lambda_c
	  fcs=1.+f/(lambda_s**i_all)/(lambda_c**j)
	  
! Drag coefficient
      ! put to 0 cdragy for the moment
       do iz_c=1,nz
        cdragx(iz_c)=fc*fs*fcs
        cdragy(iz_c)=cdragx(iz_c)
       enddo
  !    write(*,*)cdragx(1),fc,fs,fcs
! length scales
      a1=2.19
      a2=1.2
      do iz_c=1,nz
       zc=(z(iz_c)+z(iz_c+1))/2.
       if((zc/hmean).le.1.)then
        !dls(iz)=ceps*a1*(hmean-disp)*2.
		dls(iz_c)=ceps*a1*(hmean-disp)
       elseif((zc/hmean).gt.1..and.(zc/hmean).le.1.5)then
        dls(iz_c)=ceps*a1*(zc-disp)
        
       elseif((zc/hmean).gt.1.5)then
        d2=(1.-a1/a2)*1.5*hmean+a1/a2*disp
        dls(iz_c)=ceps*a2*(zc-d2)
       endif
       dlk(iz_c)=cmu/(ceps*ck)*dls(iz_c)
      enddo
      
      return
      end

!----------------------------------------
! ===6================================================================72
! ===6================================================================72

      subroutine cdtur(nzm,nz,dz,ck,tke,dlk,cdz)                                  


      !implicit none

! Input
! -----
      integer nzm                              ! maximum number of vertical levels
      integer nz                               ! number of vertical levels
      real dz                              ! levels size [m]
      real ck                                  ! von Karman constant
      real tke(nzm)                             ! turbulent kinetic energy
      real dlk(nzm)                             ! lenth scale

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

      subroutine building(nzm,nz,z,dz,pb,vl_c,vx,vy,cdragx,cdragy,wx,wy,bx,by,srim_vx,srim_vy,srex_tke)
      !implicit none
! input      
      integer nzm
      integer nz,iz_c
      real z(nzm+1),dz,pb(nzm+1),ss(nzm+1),vl(nzm)
      real cdragx(nzm),cdragy(nzm),wx,wy,bx,by
      real vx(nzm),vy(nzm)
! output      
      real srim_vx(nzm),srim_vy(nzm),srex_tke(nzm)
! local
      real lfx,lfy

      do iz_c=1,nz
       lfx=dz*by/((wx+bx)*(wy+by))*pb(iz_c+1)
       lfy=dz*bx/((wx+bx)*(wy+by))*pb(iz_c+1)
       srim_vx(iz_c)=-lfx/vl_c(iz_c)/dz*cdragx(iz_c)*abs(vx(iz_c))
       srim_vy(iz_c)=-lfy/vl_c(iz_c)/dz*cdragy(iz_c)*abs(vy(iz_c))
! not sure it is like this for tke, but for wind orthotogonal to the face of the cube, 
! vy=0.
       srex_tke(iz_c)=cdragx(iz_c)*(lfx/vl_c(iz_c)/dz*(abs(vx(iz_c))**3.))+cdragy(iz_c)*(lfy/vl_c(iz_c)/dz*(abs(vy(iz_c))**3.))
      enddo

      return
      end
      
!------------------------------------------------     

      subroutine tke_bougeault (nzm,nz,ceps,dz,vx,vy, &
                                tke,cdz,dls,sf_c, &                       ! Input
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
      real tke(nzm)                             ! turbulent kinetic energy
      real cdz(nzm+1)                           ! turbulent diffusion coefficient
      real dls(nzm)                             ! lrngth scale leps
      real sf_c(nzm+1)                            ! ?????

! Ouput
! -----
      real srim_tke(nzm)                           ! inplicit term in the tke equation
      real srex_tke(nzm)                           ! explicit term in the tke equation

! Local
! -----
      integer iz_c
      real sh(nzm)                              ! shear turbulent kinetic energy source
      real td(nz)                              ! dissipation

! ----------------------------------------------------------------------

      call shear_2(nzm,nz,dz,vx,vy,cdz, &                               ! Input
                  sh)                                              ! Ouput
!AM: here you should add the calculation of the buoyancy. You can take
!the one from boulac
      do iz_c=1,nz

       if (dls(iz_c).ne.0.) then
        td(iz_c)=-ceps*sqrt(tke(iz_c))/dls(iz_c)
       else
        td(iz_c)=0.
       end if ! dls
       sh(iz_c)=sh(iz_c)*sf_c(iz_c)
       srim_tke(iz_c)=td(iz_c)
       srex_tke(iz_c)=srex_tke(iz_c)+sh(iz_c)!+bu(iz_c) ! Am add here the buoyancy
      end do ! iz
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
   !    write(68,*)iz,sh(iz),dudz2,vx(iz+1)

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
    !   write(68,*)iz,sh(iz),dudz2,vx(iz+1)

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
          
          do iz_c=iz1,nz-izf
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

!END PROGRAM main_program

