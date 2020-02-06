program test3

use gw_convect, only : BeresSourceDesc
use gw_rdg
use gw_rdg_phunit_mod, only : gw_rdg_phunit
use gw_oro_phunit_mod, only : gw_oro_phunit
use gw_oro
use gw_utils
use gw_common, only: pver,gw_prof,gw_drag_prof,gw_common_init, GWBand
use physconst
use coords_1d,  only: Coords1D
use interpolate_data, only: lininterp
use physics_types,  only: physics_ptend,physics_ptend_init
use constituent, only: pcnst

implicit none
#include <netcdf.inc>

!!integer, parameter :: r8 = selected_real_kind(12)

!integer, parameter :: ncol =16
!integer, parameter :: nrdgs=16
!integer, parameter :: xpver=30
!integer, parameter :: xpgwv=5
      
REAL, PARAMETER :: PIc=3.1415927

integer :: ncol,ntim,nlat,nlon
integer :: nrdgs,Lchnk
integer :: xpver
integer :: pcols !!!!! WHY WHY WHY do we need both pcols and ncol ??????!!!!!

  ! Top level for gravity waves.
  integer, parameter :: ktop = 1

type physstate
   real(r8) , allocatable, dimension(:) :: lat
end type physstate

type(physstate) :: state

type(physics_ptend) :: ptend

character(len=128)  :: longchar$,err$
  ! Allow reporting of error messages.
  character(len=128) :: errstring

   integer :: SEVERITY=1

  ! Midpoint zonal/meridional winds.
  real(r8)  , allocatable :: u(:,:,:), v(:,:,:),ps(:,:),ai(:),bi(:)
  ! Midpoint temperatures.
  real(r8)  , allocatable :: t(:,:,:),q(:,:,:)
  ! Layers .
  real(r8)  , allocatable :: clh(:)
  ! Midpoint and interface pressures.
  real(r8)  , allocatable :: pmid(:,:), pint(:,:), piln(:,:),rhom(:,:)
  ! Layer thickness (pressure delta).
  real(r8)  , allocatable :: dpm(:,:),rdpm(:,:)
  ! Altitudes.
  real(r8)  , allocatable :: zm(:,:),zi(:,:)
  ! Mi
  real(r8)  , allocatable :: nm(:,:),rhoi(:,:),ti(:,:),ni(:,:)
  real(r8)  , allocatable :: dse(:,:),lats(:),lons(:),area(:)
  real(r8)  :: dt,p00


  ! Ridge parameters
  real(r8)  , allocatable :: mxdis(:,:), angll(:,:), aniso(:,:),wghts(:,:),effgw_x(:)
  real(r8)  , allocatable :: anixy(:,:), kwvrdg(:), ilnrdg(:),mxdisX(:),angllX(:),sgh(:),sgh_scaled(:)
  real(r8)  , allocatable :: hwdth(:,:), clngt(:,:), gbxar(:),leewv(:),tau0x(:),tau0y(:)
  integer   , allocatable :: isoflag(:)

  integer :: itime


  real(r8) effgw_rdg, effgw_rdg_max , rdg_cd_llb, effgw_oro
  real(r8) effgw_rdg_beta, effgw_rdg_beta_max , rdg_beta_cd_llb


  ! Newtonian cooling
  real(r8)  , allocatable :: alphi(:)


  type(Coords1D) :: p               ! Pressure coordinates



  real(r8) , allocatable :: kvtt(:,:)

  real(r8), allocatable  :: flx_heat(:),landfrac(:)

  
  real(r8) :: Rgas=287.058_r8
  real(r8) :: SHR_CONST_REARTH  = 6.37122e6_R8 

 
  integer :: nlevs,i,j,l,k,n,irdg,nwvs,nrdg_beta,nn,n_rdg_beta

  integer  :: pver_in
  integer  :: pgwv_in,xpgwv
  real(r8)  :: dc_in
  logical  :: do_molec_diff_in
  logical  :: tau_0_ubc_in
  integer  :: nbot_molec_in
  integer  :: ktop_in
  integer  :: kbotbg_in
  real(r8)  :: fcrit2_in
  real(r8)  :: kwv_in
  real(r8)  :: gravit_in
  real(r8)  :: rair_in,gscale,aniso_outer_scale, aniso_base_gridln
  real(r8)  :: qbo_hdepth_scaling, prndl

  logical   , allocatable :: trappedwave(:),lqptend(:)
  logical   :: tms_like,do_trappedwaves,do_old_gwdrag

  logical   :: trpd_leewv_rdg_beta

   integer :: nrdg_mtn, icol

  integer ncid,status, dimid,varid  ! for netCDF data file

  ! Whether or not to enforce an upper boundary condition of tau = 0.
  ! (Like many variables, this is only here to hold the value between
  ! the readnl phase and the init phase of the CAM physics; only gw_common
  ! should actually use it.)
  logical :: tau_0_ubc = .false.

  ! Maximum wave number and width of spectrum bins.
  real(r8) :: gw_dc = 2.5D0

  ! Horzontal wavelengths [m].
  real(r8), parameter :: wavelength_mid = 1.e5_r8
  real(r8), parameter :: wavelength_long = 1.e6_r8

  ! fcrit2 for the mid-scale waves has been made a namelist variable to
  ! facilitate backwards compatibility with the CAM3 version of this
  ! parameterization.  In CAM3, fcrit2=0.5.
  real(r8) :: fcrit2 = 1.0   ! critical froude number squared

  ! A mid-scale "band" with only stationary waves (l = 0).
  type(GWBand) :: band_oro
  ! A mid-scale "band" with 2*pgwv+1 spectrum of waves.
  type(GWBand) :: band_mid


  ! Efficiency for a gravity wave source.
  real(r8) , allocatable :: effgw(:)  ! (state%ncol)


  ! Stuff for Beres convective gravity wave source.
  real(r8), allocatable :: mfcc(:,:,:), hdcc(:)
  integer  :: hd_mfcc , mw_mfcc, ps_mfcc, ngwv_file
  type(BeresSourceDesc)  ::  beres_dp_desc

  
  ! Interpolated Newtonian cooling coefficients.
  real(r8) , allocatable :: alpha(:),pref_edge(:)

  ! Levels of pre-calculated Newtonian cooling (1/day).
  ! The following profile is digitized from:
  ! Wehrbein and Leovy (JAS, 39, 1532-1544, 1982) figure 5

  integer, parameter :: nalph = 71
  real(r8) :: alpha0(nalph) = [ &
       0.1_r8,         0.1_r8,         0.1_r8,         0.1_r8,         &
       0.1_r8,         0.1_r8,         0.1_r8,         0.1_r8,         &
       0.1_r8,         0.1_r8,         0.10133333_r8,  0.104_r8,       &
       0.108_r8,       0.112_r8,       0.116_r8,       0.12066667_r8,  &
       0.126_r8,       0.132_r8,       0.138_r8,       0.144_r8,       &
       0.15133333_r8,  0.16_r8,        0.17_r8,        0.18_r8,        &
       0.19_r8,        0.19933333_r8,  0.208_r8,       0.216_r8,       &
       0.224_r8,       0.232_r8,       0.23466667_r8,  0.232_r8,       &
       0.224_r8,       0.216_r8,       0.208_r8,       0.20133333_r8,  &
       0.196_r8,       0.192_r8,       0.188_r8,       0.184_r8,       &
       0.18266667_r8,  0.184_r8,       0.188_r8,       0.192_r8,       &
       0.196_r8,       0.19333333_r8,  0.184_r8,       0.168_r8,       &
       0.152_r8,       0.136_r8,       0.12133333_r8,  0.108_r8,       &
       0.096_r8,       0.084_r8,       0.072_r8,       0.061_r8,       &
       0.051_r8,       0.042_r8,       0.033_r8,       0.024_r8,       &
       0.017666667_r8, 0.014_r8,       0.013_r8,       0.012_r8,       &
       0.011_r8,       0.010333333_r8, 0.01_r8,        0.01_r8,        &
       0.01_r8,        0.01_r8,        0.01_r8                         &
       ]

  ! Pressure levels that were used to calculate alpha0 (hPa).
  real(r8) :: palph(nalph) = [ &
       2.06115E-06_r8, 2.74280E-06_r8, 3.64988E-06_r8, 4.85694E-06_r8, &
       6.46319E-06_r8, 8.60065E-06_r8, 1.14450E-05_r8, 1.52300E-05_r8, &
       2.02667E-05_r8, 2.69692E-05_r8, 3.58882E-05_r8, 4.77568E-05_r8, &
       6.35507E-05_r8, 8.45676E-05_r8, 0.000112535_r8, 0.000149752_r8, &
       0.000199277_r8, 0.000265180_r8, 0.000352878_r8, 0.000469579_r8, &
       0.000624875_r8, 0.000831529_r8, 0.00110653_r8,  0.00147247_r8,  &
       0.00195943_r8,  0.00260744_r8,  0.00346975_r8,  0.00461724_r8,  &
       0.00614421_r8,  0.00817618_r8,  0.0108801_r8,   0.0144783_r8,   &
       0.0192665_r8,   0.0256382_r8,   0.0341170_r8,   0.0453999_r8,   &
       0.0604142_r8,   0.0803939_r8,   0.106981_r8,    0.142361_r8,    &
       0.189442_r8,    0.252093_r8,    0.335463_r8,    0.446404_r8,    &
       0.594036_r8,    0.790490_r8,    1.05192_r8,     1.39980_r8,     &
       1.86273_r8,     2.47875_r8,     3.29851_r8,     4.38936_r8,     &
       5.84098_r8,     7.77266_r8,     10.3432_r8,     13.7638_r8,     &
       18.3156_r8,     24.3728_r8,     32.4332_r8,     43.1593_r8,     &
       57.4326_r8,     76.4263_r8,     101.701_r8,     135.335_r8,     &
       180.092_r8,     239.651_r8,     318.907_r8,     424.373_r8,     &
       564.718_r8,     751.477_r8,     1000._r8                        &
       ]

integer :: UNIT
logical :: use_ridge_param=.false.
logical :: use_isotropic_param=.false.
logical :: FV_dycore, SE_dycore

namelist /cntrls/  use_ridge_param , use_isotropic_param
     


! ---- Begin --------------------------------

UNIT=222
OPEN( UNIT=UNIT, FILE="control.nml" ) !, NML =  cntrls )
READ( UNIT=UNIT, NML=cntrls)
CLOSE(UNIT=UNIT)

  call  gw_rdg_readnl("control.nml")

  nrdg_mtn=1
  do_old_gwdrag=.false.
  tms_like = .false.

  longchar$='profiles_L30.dat'

  pver_in=xpver
  pgwv_in=10
  dc_in=10.
  Lchnk = 1 

  FV_dycore = .TRUE.
  SE_dycore = .FALSE.

!=================================================================
!   Odd bits of initialization
!================================================================  

pcols= -9999999  ! force model to croak if pcols is used

  band_oro = GWBand(0, gw_dc, fcrit2, wavelength_mid)

  !if (use_ridge_param)     call set_band_rdg(band_oro)
  !!if (use_isotropic_param) call set_band_oro(band_oro)

! Need to call GWBand for convective waves
!-----------------
! band_mid = GWBand(pgwv, gw_dc, 1.0_r8, wavelength_mid)
!------------------------------------------
!!! WHERE is gw_dc set? Where is pgwv set?
!!! ==> They are set in atm_in. Out of the box values seem to be
!!!  pgwv = 32
!!!  gw_dc = 2.5D0  

write(*,*) " Made band_oro !!!!!!!  "

write(*,*) " NGWV = ", band_oro%ngwv

xpgwv= band_oro%ngwv

do_trappedwaves=.true.


!===========================================================================

write(*,*) "Get Ridge data  ....  "

!=============================================================================
  status = nf_open('ridgedata.nc', 0, ncid)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
 
  status = NF_INQ_DIMID(ncid, 'nrdg', dimid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimid, nrdgs)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)


  if (SE_dycore) then
    status = NF_INQ_DIMID(ncid, 'ncol', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, ncol)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  end if

  if (FV_dycore) then
    status = NF_INQ_DIMID(ncid, 'lon', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, nlon)
    status = NF_INQ_DIMID(ncid, 'lat', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, nlat)
    ncol=nlon*nlat
  end if

  

     allocate(  mxdis(ncol,nrdgs), angll(ncol,nrdgs), anixy(ncol,nrdgs), aniso(ncol,nrdgs), &
                kwvrdg(ncol), ilnrdg(ncol), sgh(ncol),  sgh_scaled(ncol), & 
                hwdth(ncol,nrdgs), clngt(ncol,nrdgs), gbxar(ncol), & 
                wghts(ncol,nrdgs ) , mxdisX(ncol), angllX(ncol), effgw(ncol)  )
     
     ! Input needed by gw_oro
     allocate( landfrac(ncol) )


  status = NF_INQ_VARID(ncid, 'MXDIS', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, mxdis )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " MXDIS ", shape(mxdis)
     write(*,*) " MXDIS ", minval(mxdis),maxval(mxdis)

  status = NF_INQ_VARID(ncid, 'WGHTS', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, wghts )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " WGHTS ", shape(wghts)
     write(*,*) " WGHTS ", minval(wghts),maxval(wghts)

  status = NF_INQ_VARID(ncid, 'ANISO', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, aniso )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " ANISO ", shape(aniso)
     write(*,*) " ANISO ", minval(aniso),maxval(aniso)

  status = NF_INQ_VARID(ncid, 'ANIXY', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, anixy )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " ANIXY ", shape(anixy)
     write(*,*) " ANIXY ", minval(anixy),maxval(anixy)

  status = NF_INQ_VARID(ncid, 'HWDTH', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, hwdth )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " HWDTH ", shape(hwdth)
     write(*,*) " HWDTH ", minval(hwdth),maxval(hwdth)

  status = NF_INQ_VARID(ncid, 'CLNGT', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, clngt )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " CLNGT ", shape(clngt)
     write(*,*) " CLNGT ", minval(clngt),maxval(clngt)

  status = NF_INQ_VARID(ncid, 'ANGLL', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, angll )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " ANGLL ", shape(angll)
     write(*,*) " ANGLL ", minval(angll),maxval(angll)
 
  status = NF_INQ_VARID(ncid, 'SGH', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, sgh )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " SGH ", shape(sgh)
     write(*,*) " SGH ", minval(sgh),maxval(sgh)

  status = NF_INQ_VARID(ncid, 'GBXAR', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, gbxar )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " GBXAR ", shape(gbxar)
     write(*,*) " GBXAR ", minval(gbxar),maxval(gbxar)

  status = nf_close (ncid)
  if (status .ne. NF_NOERR) call handle_err(status)          



  ! Convert grid box areas from Steradians to km+2
  gbxar = gbxar*( SHR_CONST_REARTH/1000._r8)*( SHR_CONST_REARTH/1000._r8)

     write(*,*) " ... "
     write(*,*) " scaled HWDTH ", minval(hwdth),maxval(hwdth)
 
  landfrac(:) = 1.0
 
  sgh_scaled = sgh 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  status = nf_open('/Users/juliob/cesm_inputdata/newmfspectra40_dc25.nc', 0, ncid)


    status = NF_INQ_DIMID(ncid, 'PS', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, ps_mfcc )

    status = NF_INQ_DIMID(ncid, 'MW', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, mw_mfcc )

    status = NF_INQ_DIMID(ncid, 'HD', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, hd_mfcc )

 !  allocate( mfcc(ps_mfcc , mw_mfcc, hd_mfcc ) )
  allocate( mfcc(hd_mfcc , mw_mfcc, ps_mfcc ) )
  allocate( hdcc(hd_mfcc) )
   
  status = NF_INQ_VARID(ncid, 'HD', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, hdcc )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = NF_INQ_VARID(ncid, 'mfcc', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, mfcc )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)



     write(*,*) " MFCC ", shape(mfcc)
     write(*,*) " MFCC ", minval(mfcc),maxval(mfcc)


     write(*,*) " HD_MFCC=",hd_mfcc
     write(*,*) " MW_MFCC=",mw_mfcc
     write(*,*) " PS_MFCC=",ps_mfcc

  status = nf_close (ncid)


  ! Get HD (heating depth) dimension.
  beres_dp_desc%maxh = HD_MFCC  ! get_pio_dimlen(gw_file_desc, "HD", file_path)

  ! Get MW (mean wind) dimension.
  beres_dp_desc%maxuh = MW_MFCC ! get_pio_dimlen(gw_file_desc, "MW", file_path)

  ! Get PS (phase speed) dimension.
  ngwv_file = ps_mfcc  ! get_pio_dimlen(gw_file_desc, "PS", file_path)

  ! Number in each direction is half of total (and minus phase speed of 0).
  beres_dp_desc%maxuh = (beres_dp_desc%maxuh-1)/2
  ngwv_file = (ngwv_file-1)/2



!===========================================================================
    write(*,*) " Read in MFCC "
!===========================================================================
    write(*,*) " NOw will read in Metdata "
!===========================================================================


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  status = nf_open('metdata.nc', 0, ncid)
  IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
  
  if (SE_dycore) then
    status = NF_INQ_DIMID(ncid, 'ncol', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, ncol)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  end if

  if (FV_dycore) then
    status = NF_INQ_DIMID(ncid, 'lon', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, nlon)
    status = NF_INQ_DIMID(ncid, 'lat', dimid)
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    status = NF_INQ_DIMLEN(ncid, dimid, nlat)
    ncol=nlon*nlat
  end if

  status = NF_INQ_DIMID(ncid, 'lev', dimid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimid, xpver)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
 
  status = NF_INQ_DIMID(ncid, 'time', dimid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimid, ntim)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

   
      if (FV_dycore)  allocate( lats(nlat) , lons(nlon) )
      if (SE_dycore)  allocate( lats(ncol) , lons(ncol) )
      allocate(  state%lat(ncol)  , area(ncol)  )

     ! Dynamics import state 
     allocate(  u(ncol,xpver,ntim), v(ncol,xpver,ntim), t(ncol,xpver,ntim), q(ncol,xpver,1), ps(ncol,ntim) )
     !------
     ! A's and B's     
     allocate(  ai(xpver+1), bi(xpver+1)  )

  status = NF_INQ_VARID(ncid, 'P0', varid)
  IF (status .NE. NF_NOERR) then 
  !   CALL HANDLE_ERR(status)
      p00=100000._r8
      write(*,*) " Manually setting P00 "
  ELSE
    status = NF_GET_VAR_DOUBLE(ncid, varid, p00 )
    IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  ENDIF
     write(*,*) " P_0 ",p00

  status = NF_INQ_VARID(ncid, 'hyai', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, ai )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
    

  status = NF_INQ_VARID(ncid, 'hybi', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, bi )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = NF_INQ_VARID(ncid, 'lon', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, lons )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " LONS ", shape(lons)
     write(*,*) " LONS ", minval(lons),maxval(lons)

  
  
  status = NF_INQ_VARID(ncid, 'lat', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, lats )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " LATS ", shape(lats)
     write(*,*) " LATS ", minval(lats),maxval(lats)

  do j=1,nlat
  do i=1,nlon
     icol = nlon*(j-1)+i
     state%lat(icol) = lats(j)
  end do
  end do



!-----------------------------
! Read in UVT etc

  status = NF_INQ_VARID(ncid, 'U', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, u )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " U ", shape(u)
     write(*,*) " U ", minval(u),maxval(u)

  status = NF_INQ_VARID(ncid, 'V', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, v )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " V ", shape(v)
     write(*,*) " V ", minval(v),maxval(v)

  status = NF_INQ_VARID(ncid, 'T', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, t )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " T ", shape(t)
     write(*,*) " T ", minval(t),maxval(t)

  status = NF_INQ_VARID(ncid, 'PS', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, ps )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " PS ", shape(ps)
     write(*,*) " PS ", minval(ps),maxval(ps)

  status = nf_close (ncid)
  if (status .ne. NF_NOERR) call handle_err(status)          

!===========================================================================
!  Finished specifying U,V,T,PS  
!===========================================================================
    write(*,*) " NOw creating extra fields needed by OGWD "
!===========================================================================



     ! Variables for 3D pressure, geopht manipulations etc..
     allocate(  zm(ncol,xpver),zi(ncol,xpver+1) )
     allocate(  pmid(ncol,xpver), pint(ncol,xpver+1), piln(ncol,xpver+1),dpm(ncol,xpver),rdpm(ncol,xpver)  )
     !------

     ! Other Thermo profiles
     allocate(  rhom(ncol,xpver), rhoi(ncol,xpver+1) )
     allocate(  nm(ncol,xpver), ti(ncol,xpver+1), ni(ncol,xpver+1), dse(ncol,xpver)  )
     !------
     
     ! Needed outputs:
     allocate( kvtt(ncol,xpver+1) )
     allocate(  flx_heat(ncol) )
     
     ! Need for pre-calculated Newtonian cooling
     allocate(alphi(xpver+1))
     allocate(alpha(xpver+1))
     allocate(pref_edge(xpver+1))


    
    !-----------------------------------------------------
    ! calculate 3D pressure
    ! Much of this is also done within Coords1D

  DO l=1,xpver+1
    pint(:,l) = ai(l)* p00 + bi(l)*ps(:,1)
  end do
  DO l=1,xpver+1
     pref_edge(l) = ai(l)* p00 + bi(l)*p00
  end do
  
     write(*,*) " PINT ", shape(pint)
     write(*,*) " PINT ", minval(pint),maxval(pint)

  pmid = 0.5*( pint(:,2:xpver+1) +  pint(:,1:xpver) )
  dpm  = pint(:,2:xpver+1) -  pint(:,1:xpver)
     write(*,*) " PMID ", shape(pmid)
     write(*,*) " PMID ", minval(pmid),maxval(pmid)
     write(*,*) " DPM  ", shape(dpm)
     write(*,*) " DPM  ", minval(dpm),maxval(dpm)

  dt=1800.

  rdpm=1./dpm

  piln=log( pint )

  rhom = pmid / (Rgas * T(:,:,1) )

     write(*,*) " RHOM ", shape(rhom)
     write(*,*) " RHOM ", minval(rhom),maxval(rhom)

  zi =0.
  DO l=xpver,1,-1
    zi(:,l) = zi(:,l+1) +dpm(:,l)/rhom(:,l) 
  end do

  zi = zi / 9.80616_r8
 
    write(*,*) " ZI ", shape(zi)
     write(*,*) " ZI ", minval(zi),maxval(zi)
  zm = 0.5*( zi(:,2:xpver+1) +  zi(:,1:xpver) )

    write(*,*) " ZM ", shape(zm)
     write(*,*) " ZM ", minval(zm),maxval(zm)
     

  ! pre-calculated newtonian damping:
  !     * convert to 1/s
  !     * ensure it is not smaller than 1e-6
  !     * convert palph from hpa to pa
  do k=1,nalph
     alpha0(k) = alpha0(k) / 86400._r8
     alpha0(k) = max(alpha0(k), 1.e-6_r8)
     palph(k) = palph(k)*1.e2_r8
  end do

  ! interpolate to current vertical grid and obtain alpha
  call lininterp (alpha0  ,palph, nalph , alpha  , pref_edge , xpver+1)


     q=0.
     !c=0.

     !!call gw_common_init(xpver,xpgwv,rair,gravit)

     !! Create all 3D P-fields from intfc P's
     !! This is a nice idea. For some reason GW codes need "piln" - log of pressure - also
  p = Coords1D(pint(:ncol,:))

write(*,*) " Pressure derived type from Coords1D "
write(*,*) shape(p%ifc)
write(*,*) shape(p%del)
write(*,*) minval( p%ifc(:,31) ),maxval( p%ifc(:,31) )



!!!  IMHO stupidly, this call to gw_common_init is needed
!!!  to tell every GW routine how many levels (pver/xpver) there are.
!!!
!!!  Newtonian cooling coeffs alpha are set here for use by gw_drag_prof.
!!!  These alpha are declared in pre-amble of gw_common, and therefore are
!!!  available to gw_drag_prof in same module.
!!!--------------------------------------------
  prndl		= 0.5d0  
  qbo_hdepth_scaling		= 0.25D0
  alphi=0.0 ! 0.01
  call gw_common_init(xpver,&
       tau_0_ubc, ktop, gravit, rair, alpha, & 
       prndl, qbo_hdepth_scaling, & 
       errstring)


   allocate( lqptend(pcnst))
   lqptend(:)=.TRUE.
   call physics_ptend_init(ptend, ncol , 'OGWD',.TRUE.,.TRUE.,.TRUE., lqptend )

write(*,*) " Met specification complete  "
write(*,*) "   ...   "

    !eff2(:,:)=0. ! ++jtb Init for safety
     !---------------------------------------------------------------------
     ! Anisotropic orographic gravity waves
     ! Should transition to old GW if do_old_gwdrag=.true.
     !---------------------------------------------------------------------


write(*,*) " All preps done  ..."
write(*,*) " Ready to do gw drag calculations "

write(*,*) "P before gw prof"
write(*,*) minval( p%ifc(:,31) ),maxval( p%ifc(:,31) )

do itime=1,ntim
write(*,*) "in time loop "

  !-----------------------------------------------
  ! Profiles of Thermo. background state variables
  ! needed by gravity wave calculations.
  !-----------------------------------------------
  call gw_prof(ncol, p, cpair, t(:,:,itime), rhoi, nm, ni)


write(*,*) itime," out of ",ntim
write(*,*) "P after gw prof"
write(*,*) minval( p%ifc(:,31) ),maxval( p%ifc(:,31) )


    write(*,*) " NM ", shape(nm)
     write(*,*) " NM ", minval(nm),maxval(nm)

write(*,*) " Got through gw_prof .... "



write(*,*) "P before oro src"
write(*,*) minval( p%ifc(:,31) ),maxval( p%ifc(:,31) )


where(mxdis < 0.)
  mxdis=0.
end where


     ! Stupid double setting of parameters
     ! Not sure anymore why this here.
     !------------------------------------
     effgw_rdg = 1.0_r8
     effgw_rdg_max = 1.0_r8
     rdg_cd_llb = 1.0_r8
     rdg_beta_cd_llb = 1.0_r8
     nrdg_beta =10

     trpd_leewv_rdg_beta = .FALSE.
     effgw_rdg_beta =  effgw_rdg
     effgw_rdg_beta_max =  effgw_rdg_max
     if (use_ridge_param) then
         n_rdg_beta = nrdg_beta
     else if (use_isotropic_param) then
         n_rdg_beta=1
     end if

     ! nml settings for old GW scheme
     !-------------------------------
     effgw_oro = 0.125D0

     kvtt = 0.000001_r8
    
    write(511) nlon,nlat,pver,n_rdg_beta,itime,ntim
    write(511) lons,lats,state%lat

    if (use_ridge_param) then
     call gw_rdg_phunit(                           &
       'BETA ', band_oro,                         &
        ncol, lchnk, n_rdg_beta, dt,              &
        u(:,:,itime), v(:,:,itime), t(:,:,itime), &
        p, piln, zm, zi,                          &
        nm, ni, rhoi, kvtt, q, dse,               &
        effgw_rdg_beta, effgw_rdg_beta_max,       &
        hwdth, clngt, gbxar, mxdis, angll, anixy, &
        rdg_beta_cd_llb, trpd_leewv_rdg_beta,     &
        ptend, flx_heat)
    endif

    if (use_isotropic_param) then
       call gw_oro_phunit( &
      band_oro, &
      ncol, lchnk, dt, effgw_oro, &
      u(:,:,itime), v(:,:,itime), t(:,:,itime), &
      p, piln, zm, zi, &
      nm, ni, rhoi, kvtt, q, dse, &
      sgh, landfrac,lats, &
      ptend, flx_heat)
    endif


end do

end program test3


!************************************************************************
!!handle_err
!************************************************************************
!
!!ROUTINE:      handle_err
!!DESCRIPTION:  error handler
!--------------------------------------------------------------------------

subroutine handle_err(status)
  
  implicit         none
  
#include <netcdf.inc>
  
  integer          status
  
  if (status .ne. nf_noerr) then
    print *, nf_strerror(status)
    stop 'Stopped'
  endif
  
end subroutine handle_err

