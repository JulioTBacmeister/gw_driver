program xesm
  
use dycore_mod, only: SE_dycore,FV_dycore,read_dynamics_state,read_convective_heating,read_ridge_file
use gw_convect, only : gw_beres_init,gw_beres_run
use gw_rdg,     only : gw_rdg_init, gw_rdg_run
use gw_utils
use gw_common, only: gw_prof,gw_drag_prof,gw_common_init, GWBand
use physconst
use coords_1d,  only: Coords1D
use interpolate_data, only: lininterp
!use physics_types,  only: physics_ptend,physics_ptend_init

implicit none
#include <netcdf.inc>

      
REAL, PARAMETER :: PIc=3.1415927

integer :: ncol,ntim,nlat,nlon
integer :: nrdgs
integer :: xpver , pver, pcnst, pverp
integer :: pcols !!!!! WHY WHY WHY do we need both pcols and ncol ??????!!!!!

  ! Top level for gravity waves.
  integer, parameter :: ktop = 1

type physstate
   real(r8) , allocatable, dimension(:) :: lat
end type physstate

type(physstate) :: state

!!type(physics_ptend) :: ptend

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
  real(r8)  , allocatable :: delp(:,:),rdelp(:,:)
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



  real(r8) , allocatable :: kvtt(:,:) , utgw(:,:) , vtgw(:,:), ttgw(:,:), qtgw(:,:,:)

  real(r8), allocatable  :: flx_heat(:),landfrac(:)

  
  real(r8) :: Rgas=287.058_r8
  real(r8) :: SHR_CONST_REARTH  = 6.37122e6_R8 

 
  integer :: nlevs,i,j,l,k,n,irdg,nwvs,nrdg_beta,nn,n_rdg_beta

  !!integer  :: pver_in
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

  !!! call  gw_rdg_readnl("control.nml")

  nrdg_mtn=1
  do_old_gwdrag=.false.
  tms_like = .false.

  longchar$='profiles_L30.dat'

  !pver_in=xpver
  pgwv_in=10
  dc_in=10.
  pcnst = 1


!=================================================================
!   Odd bits of initialization
!================================================================  

pcols= -9999999  ! force model to croak if pcols is used


do_trappedwaves=.true.


!===========================================================================

write(*,*) "Get Ridge data  ....  "

!=============================================================================
  !!!status = nf_open('ridgedata.nc', 0, ncid)

  call read_ridge_file( '/Users/juliob/cesm_inputdata/ridgedata.nc' , hwdth, clngt, gbxar, mxdis, angll, anixy, sgh  )





!===========================================================================
    write(*,*) " Read in MFCC "

  
!===========================================================================
    write(*,*) " NOw will read in Metdata "
!===========================================================================


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  

  call read_dynamics_state ( '/Users/juliob/cesm_inputdata/Nudge_FWsc_e20b07_BSLN_Conc_05.cam.h1.2011-01-01-00000.nc' , ai, bi, p00, u, v, t, q, ps, lons, lats, landfrac  )
  nlon=size( lons, 1)
  nlat=size(lats, 1)
  ntim=1
  xpver=size( ai , 1)-1
  pver  = xpver
  pverp = pver+1
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
   
   call read_convective_heating ( '/Users/juliob/cesm_inputdata/Nudge_FWsc_e20b07_BSLN_Conc_05.cam.h3.NETDT.2011-01-01-00000.nc'  , netdt  )

!===========================================================================
!  Finished specifying U,V,T,PS  
!===========================================================================
    write(*,*) " NOw creating extra fields needed by OGWD "
!===========================================================================



     ! Variables for 3D pressure, geopht manipulations etc..
     allocate(  zm(ncol,pver),zi(ncol,pver+1) )
     allocate(  pmid(ncol,pver), pint(ncol,pver+1), piln(ncol,pver+1),delp(ncol,pver),rdelp(ncol,pver)  )
     !------

     ! Other Thermo profiles
     allocate(  rhom(ncol,pver), rhoi(ncol,pver+1) )
     allocate(  nm(ncol,pver), ti(ncol,pver+1), ni(ncol,pver+1), dse(ncol,pver)  )
     !------
     ! Molecular diffusivity input
     allocate( kvtt(ncol,pver+1)  )

     ! Needed outputs:
     allocate( utgw(ncol,pver),  vtgw(ncol,pver),  ttgw(ncol,pver)  )
     allocate( qtgw(ncol,pver,pcnst)  )
     allocate(  flx_heat(ncol) )
     
     ! Need for pre-calculated Newtonian cooling
     ! allocate(alpha(pver+1))
     allocate(pref_edge(pver+1))

    ! allocate(  netdt(ncol,pver)  )

   dt=1800.
   
    !-----------------------------------------------------
    ! calculate 3D pressure
 
  DO l=1,pver+1
    pint(:,l) = ai(l)* p00 + bi(l)*ps(:,1)
  end do
  DO l=1,pver+1
     pref_edge(l) = ai(l)* p00 + bi(l)*p00
  end do
  pmid = 0.5*( pint(:,2:pver+1) +  pint(:,1:pver) )
  delp  = pint(:,2:pver+1) -  pint(:,1:pver)


  rdelp=1./delp
  piln=log( pint )

  rhom = pmid / (Rgas * T(:,:,1) )

  ! Calculate geopotential height
  !---------------------------------
  zi =0.
  DO l=pver,1,-1
    zi(:,l) = zi(:,l+1) +delp(:,l)/rhom(:,l) 
  end do
  zi = zi / 9.80616_r8
  zm = 0.5*( zi(:,2:pver+1) +  zi(:,1:pver) )

!!!  IMHO stupidly, this call to gw_common_init is needed
!!!  to tell every GW routine how many levels (pver/xpver) there are.
!!!
!!!  Newtonian cooling coeffs alpha are set here for use by gw_drag_prof.
!!!  These alpha are declared in pre-amble of gw_common, and therefore are
!!!  available to gw_drag_prof in same module.
!!!--------------------------------------------

  call gw_common_init(pverp, pref_edge )

  call gw_rdg_init( )

  call gw_beres_init ( '/Users/juliob/cesm_inputdata/newmfspectra40_dc25.nc' , pver, pverp, pref_edge  )



   !allocate( lqptend(pcnst))
   !lqptend(:)=.TRUE.
   !call physics_ptend_init(ptend, ncol, pver, pcnst, 'OGWD',.TRUE.,.TRUE.,.TRUE., lqptend )



do itime=1,ntim
write(*,*) "in time loop "



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


       utgw=0._r8
       vtgw=0._r8
       ttgw=0._r8
       qtgw=0._r8


    
    use_ridge_param=.TRUE.
    if (use_ridge_param) then

     write(*,*) " gw_rdg_run "
     call gw_rdg_run(                           &
       'BETA ',                                  &
        ncol, pver, pverp, pcnst, n_rdg_beta, dt,              &
        u(:,:,itime), v(:,:,itime), t(:,:,itime), &
        pint,pmid,delp, piln, zm, zi,                          &
        kvtt, q , dse,               &
        effgw_rdg_beta, effgw_rdg_beta_max,       &
        hwdth, clngt, gbxar, mxdis, angll, anixy, &
        rdg_beta_cd_llb, trpd_leewv_rdg_beta,     &
        flx_heat, utgw, vtgw, ttgw, qtgw )
     
     
     endif
     write(*,*) "After gw_rdg_run:"
     write(*,*) "Min Max ptend%U SUM ABS "
     write(*,*) minval( utgw ), maxval(utgw), sum( abs(utgw) )

     

       effgw_dpcu=1.0_r8


       utgw=0._r8
       vtgw=0._r8
       ttgw=0._r8
       qtgw=0._r8


       write(*,*) " shape LATS ",shape(lats)
       
     call gw_beres_run(  &
          ncol, pver,pverp, pcnst, dt, effgw_dpcu,  &
        u(:,:,itime), v(:,:,itime), t(:,:,itime), &
        pint, pmid, delp, piln, zm, zi,                          &
        kvtt, q, dse,               &
          netdt , state%lat, &
          flx_heat, utgw, vtgw, ttgw, qtgw  )

     write(*,*) "gw_beres_run:"
     write(*,*) "Min Max ptend%U "
     write(*,*) minval( utgw ), maxval(utgw), sum( abs(utgw) )



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

