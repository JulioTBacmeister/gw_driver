program xesm
  
use dycore_mod, only: SE_dycore,FV_dycore,read_dynamics_state,read_convective_heating
use gw_convect, only : BeresSourceDesc,gw_beres_init,gw_beres_ifc
use gw_rdg, only : gw_rdg_init, gw_rdg_readnl,gw_rdg_ifc
!use gw_rdg_phunit_mod, only  : gw_rdg_phunit
!use gw_oro_phunit_mod, only  : gw_oro_phunit
!use gw_dpcu_phunit_mod, only : gw_dpcu_phunit
!use gw_oro
use gw_utils
!use gw_newtonian, only: gw_newtonian_set
use gw_common, only: pver,gw_prof,gw_drag_prof,gw_common_init, GWBand, gw_newtonian_set
use physconst
use coords_1d,  only: Coords1D
use interpolate_data, only: lininterp
use physics_types,  only: physics_ptend,physics_ptend_init
use constituent, only: pcnst

implicit none
#include <netcdf.inc>

      
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
  ! Convective heating
  real(r8)  , allocatable :: netdt(:,:,:)


  ! Ridge parameters
  real(r8)  , allocatable :: mxdis(:,:), angll(:,:), aniso(:,:),wghts(:,:),effgw_x(:)
  real(r8)  , allocatable :: anixy(:,:), kwvrdg(:), ilnrdg(:),mxdisX(:),angllX(:),sgh(:),sgh_scaled(:)
  real(r8)  , allocatable :: hwdth(:,:), clngt(:,:), gbxar(:),leewv(:),tau0x(:),tau0y(:)
  integer   , allocatable :: isoflag(:)

  integer :: itime


  real(r8) effgw_rdg, effgw_rdg_max , rdg_cd_llb, effgw_oro, effgw_dpcu
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

integer :: UNIT
logical :: use_ridge_param=.false.
logical :: use_isotropic_param=.false.

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


!=================================================================
!   Odd bits of initialization
!================================================================  

pcols= -9999999  ! force model to croak if pcols is used


do_trappedwaves=.true.


!===========================================================================

write(*,*) "Get Ridge data  ....  "

!=============================================================================
  !!!status = nf_open('ridgedata.nc', 0, ncid)

  call gw_rdg_init( '/Users/juliob/cesm_inputdata/ridgedata.nc' , hwdth, clngt, gbxar, mxdis, angll, anixy, sgh, band_oro )

  ncol = SIZE( SGH, 1)
  allocate(  landfrac(ncol) )

  ! Convert grid box areas from Steradians to km+2
  gbxar = gbxar*( SHR_CONST_REARTH/1000._r8)*( SHR_CONST_REARTH/1000._r8)

  landfrac(:) = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !status = nf_open('/Users/juliob/cesm_inputdata/newmfspectra40_dc25.nc', 0, ncid)

  call gw_beres_init ( '/Users/juliob/cesm_inputdata/newmfspectra40_dc25.nc' , band_mid,  beres_dp_desc )


!===========================================================================
    write(*,*) " Read in MFCC "

  
!===========================================================================
    write(*,*) " NOw will read in Metdata "
!===========================================================================


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  

  call read_dynamics_state ( '/Users/juliob/cesm_inputdata/Nudge_FWsc_e20b07_BSLN_Conc_05.cam.h1.2011-01-01-00000.nc' , ai, bi, p00, u, v, t, q, ps, lons, lats  )
  nlon=size( lons, 1)
  nlat=size(lats, 1)
  ntim=1
  xpver=size( ai , 1)-1
  if (FV_dycore) then
      ncol=nlat*nlon
      allocate(  state%lat(ncol)   )

      do j=1,nlat
      do i=1,nlon
         icol = nlon*(j-1)+i
         state%lat(icol) = lats(j)
      end do
      end do
   endif
  if (SE_dycore) then
      ncol=nlat
      allocate(  state%lat(ncol)   )
      state%lat(:) = lats(:)
   endif

!         10        20        30        40        50        60        70        80        90       100       110       120       130
!         |         |         |         |         |         |         |         |         |         |         |         |         |   
!1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
   
   call read_convective_heating ( '/Users/juliob/cesm_inputdata/Nudge_FWsc_e20b07_BSLN_Conc_05.cam.h3.2011-01-01-00000.nc'  , netdt  )

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
     allocate(alpha(xpver+1))
     allocate(pref_edge(xpver+1))

    ! allocate(  netdt(ncol,xpver)  )

    
    !-----------------------------------------------------
    ! calculate 3D pressure
    ! Much of this is also done within Coords1D

  DO l=1,xpver+1
    pint(:,l) = ai(l)* p00 + bi(l)*ps(:,1)
  end do
  DO l=1,xpver+1
     pref_edge(l) = ai(l)* p00 + bi(l)*p00
  end do
  

  pmid = 0.5*( pint(:,2:xpver+1) +  pint(:,1:xpver) )
  dpm  = pint(:,2:xpver+1) -  pint(:,1:xpver)

  dt=1800.

  rdpm=1./dpm

  piln=log( pint )

  rhom = pmid / (Rgas * T(:,:,1) )


  zi =0.
  DO l=xpver,1,-1
    zi(:,l) = zi(:,l+1) +dpm(:,l)/rhom(:,l) 
  end do

  zi = zi / 9.80616_r8

  
    zm = 0.5*( zi(:,2:xpver+1) +  zi(:,1:xpver) )


     call gw_newtonian_set( xpver, pref_edge , alpha )
  
     q=0.

     !! Create all 3D P-fields from intfc P's
     !! This is a nice idea. For some reason GW codes need "piln" - log of pressure - also
  p = Coords1D(pint(:ncol,:))


     !! This should be in beres_init.  Pref_edge should be available immediately shouldnt it??
     do k = 0, pver
        ! 700 hPa index
        if (pref_edge(k+1) < 70000._r8) beres_dp_desc%k = k+1
     end do



!!!  IMHO stupidly, this call to gw_common_init is needed
!!!  to tell every GW routine how many levels (pver/xpver) there are.
!!!
!!!  Newtonian cooling coeffs alpha are set here for use by gw_drag_prof.
!!!  These alpha are declared in pre-amble of gw_common, and therefore are
!!!  available to gw_drag_prof in same module.
!!!--------------------------------------------
  prndl		= 0.5d0  
  qbo_hdepth_scaling		= 0.25D0
  call gw_common_init(xpver,&
       tau_0_ubc, ktop, gravit, rair, alpha, & 
       prndl, qbo_hdepth_scaling, & 
       errstring)


   allocate( lqptend(pcnst))
   lqptend(:)=.TRUE.
   call physics_ptend_init(ptend, ncol , 'OGWD',.TRUE.,.TRUE.,.TRUE., lqptend )



do itime=1,ntim
write(*,*) "in time loop "

  !-----------------------------------------------
  ! Profiles of Thermo. background state variables
  ! needed by gravity wave calculations.
  ! This should probably go into gw_*_ifc routines.
  !-----------------------------------------------
  call gw_prof(ncol, p, cpair, t(:,:,itime), rhoi, nm, ni)



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

    use_ridge_param=.TRUE.
    if (use_ridge_param) then

     write(*,*) " gw_rdg_ifc "
     call gw_rdg_ifc(                           &
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
     write(*,*) "After gw_rdg_ifc:"
     write(*,*) "Min Max ptend%U "
     write(*,*) minval( ptend%U ), maxval(ptend%U)

     

       effgw_dpcu=1.0_r8

     
!, band_mid,  beres_dp_desc )
     call gw_beres_ifc( band_mid, &
          ncol, lchnk, dt, effgw_dpcu,  &
        u(:,:,itime), v(:,:,itime), t(:,:,itime), &
        p, piln, zm, zi,                          &
        nm, ni, rhoi, kvtt, q, dse,               &
          netdt,beres_dp_desc,lats, &
          ptend, flx_heat )

     write(*,*) "gw_beres_ifc:"
     write(*,*) "Min Max ptend%U "
     write(*,*) minval( ptend%U ), maxval(ptend%U)



end do

end program xesm


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

