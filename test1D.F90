program test1

use gw_rdg
use gw_rdg_phunit_mod, only : gw_rdg_phunit,set_band_rdg,pcols
use gw_oro_calc_mod, only : gw_oro_calc,set_band_oro
use gw_oro
use gw_utils
use gw_common, only: gw_prof,gw_drag_prof,gw_common_init, GWBand
use physconst
use coords_1d,  only: Coords1D
use interpolate_data, only: lininterp
use physics_types,  only: physics_ptend,physics_ptend_init
!!use constituent, only: pcnst

implicit none
      

   integer, parameter :: pver=30
   !integer, parameter :: pcols=100
   integer, parameter :: ncol=1000
   integer, parameter :: prdg=16
   integer, parameter :: pcnst=1

   character(len=5) :: type         ! BETA or GAMMA
   !!integer    :: ncol         ! number of atmospheric columns
   integer    :: lchnk        ! chunk identifier
   integer    :: n_rdg
   real(r8)  :: dt           ! Time step.

   real(r8)  :: u(ncol,pver)    ! Midpoint zonal winds. ( m s-1)
   real(r8)  :: v(ncol,pver)    ! Midpoint meridional winds. ( m s-1)
   real(r8)  :: t(ncol,pver)    ! Midpoint temperatures. (K)
   real(r8)  :: th(ncol,pver)    ! Midpoint temperatures. (K)
   type(Coords1D)  :: p               ! Pressure coordinates.
   real(r8)  :: piln(ncol,pver+1)  ! Log of interface pressures.
   real(r8)  :: pint(ncol,pver+1)  ! interface pressures.
   real(r8)  :: zm(ncol,pver)   ! Midpoint altitudes above ground (m).
   real(r8)  :: zi(ncol,pver+1) ! Interface altitudes above ground (m).
   real(r8)  :: nm(ncol,pver)   ! Midpoint Brunt-Vaisalla frequencies (s-1).
   real(r8)  :: ni(ncol,pver+1) ! Interface Brunt-Vaisalla frequencies (s-1).
   real(r8)  :: rhoi(ncol,pver+1) ! Interface density (kg m-3).
   real(r8)  :: rhom(ncol,pver) ! Midpoint density (kg m-3).
   real(r8)  :: kvtt(ncol,pver+1) ! Molecular thermal diffusivity.
   real(r8)  :: q(ncol,pver,pcnst)        ! Constituent array.
   real(r8)  :: dse(ncol,pver)  ! Dry static energy.


   real(r8)  :: effgw_rdg ,  effgw_rdg_beta       ! Tendency efficiency.
   real(r8)  :: effgw_rdg_max, effgw_rdg_beta_max
   real(r8)  :: hwdth(ncol,prdg) ! width of ridges.
   real(r8)  :: clngt(ncol,prdg) ! length of ridges.
   real(r8)  :: gbxar(ncol)      ! gridbox area
   real(r8)  :: sgh(ncol)        ! 
   real(r8)  :: landfrac(ncol)   ! 

   real(r8)  :: mxdis(ncol,prdg) ! Height estimate for ridge (m).
   real(r8)  :: angll(ncol,prdg) ! orientation of ridges.
   real(r8)  :: anixy(ncol,prdg) ! Anisotropy parameter.

   real(r8)  :: rdg_beta_cd_llb,rdg_cd_llb      ! Drag coefficient for low-level flow
   logical    :: trpd_leewv

   type(physics_ptend)  :: ptend   ! Parameterization net tendencies.

   real(r8)   :: flx_heat(ncol)
   real(r8)   :: effgw_oro

   integer :: nrdg_beta, n_rdg_beta, nlon, nlat, ntim, itime, ktop, i,j,k, L

     logical :: trpd_leewv_rdg_beta 

  ! A mid-scale "band" with only stationary waves (l = 0).
  type(GWBand) :: band_oro

   real(r8) :: lons(ncol),lats(ncol)

  ! Maximum wave number and width of spectrum bins.
  real(r8) :: gw_dc = 5.0

  ! Horzontal wavelengths [m].
  real(r8), parameter :: wavelength_mid = 1.e5_r8
  ! fcrit2 for the mid-scale waves has been made a namelist variable to
  ! facilitate backwards compatibility with the CAM3 version of this
  ! parameterization.  In CAM3, fcrit2=0.5.
  real(r8) :: fcrit2 = 1.0   ! critical froude number squared
  logical :: tau_0_ubc = .false.


  !real(r8)  :: rair,gravit
  real(r8)  :: qbo_hdepth_scaling, prndl

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

  ! Interpolated Newtonian cooling coefficients.
  real(r8) :: alpha(pver+1), alphi(pver+1),pref_edge(pver+1), p00,ai(pver+1),bi(pver+1)
  character(len=128) :: errstring
  real(r8) :: Rgas=287.058_r8
  real(r8) :: SHR_CONST_REARTH  = 6.37122e6_R8 
  logical   , allocatable :: lqptend(:)
  logical :: use_ridge_param=.false.
  logical :: use_isotropic_param=.false.
  integer :: UNIT
  namelist /cntrls/  use_ridge_param , use_isotropic_param

UNIT=222
OPEN( UNIT=UNIT, FILE="control.nml" ) !, NML =  cntrls )
READ( UNIT=UNIT, NML=cntrls)
CLOSE(UNIT=UNIT)

  call  gw_rdg_readnl("control.nml")

  open(unit=111,file='ab31edge.dat',form='FORMATTED')
  DO L=1,pver+1
     read(111,911) ai(l),bi(l)
     !!write(6,911)  ai(l),bi(l)
  end do
  close(unit=111)
911 format( 2(f20.15,1x) )

p00 = 100000._r8
  dt=1800.


  DO L=1,pver+1
     pref_edge(l) = ai(l)* p00 + bi(l)*p00
     !!write(6,*)  ai(l),bi(l),pref_edge(L)
  end do


  DO i=1,ncol
     pint(i,:)  = pref_edge(:)
  end do


  p = Coords1D(pint(:ncol,:))
write(*,*) " Pressure derived type from Coords1D "
write(*,*) shape(p%ifc)
write(*,*) shape(p%del)
write(*,*) minval( p%ifc(:,31) ),maxval( p%ifc(:,31) )


  band_oro = GWBand(0, gw_dc, fcrit2, wavelength_mid)
  if (use_ridge_param)      call set_band_rdg(band_oro)
  if (use_isotropic_param)  call set_band_oro(band_oro)


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
  call lininterp (alpha0  ,palph, nalph , alpha  , pref_edge , pver+1)


  pcols=ncol
  ktop=1
  prndl		= 0.5d0  
  qbo_hdepth_scaling		= 0.25D0
  alphi=0.0 ! 0.01
  call gw_common_init(pver,&
       tau_0_ubc, ktop, gravit, rair, alpha, & 
       prndl, qbo_hdepth_scaling, & 
       errstring)

   allocate( lqptend(pcnst))
   lqptend(:)=.TRUE.
   call physics_ptend_init(ptend, ncol , 'OGWD',.TRUE.,.TRUE.,.TRUE., lqptend )

  write(*,*) " finished init phase "

  th = 275._r8 *((p00/p%mid)**0.15_r8 )
  !! th = 275._r8 *((p00/p%mid)**0.40_r8 )
  !! th = 300._r8

  t  = th*(( p%mid/p00 )**(rair/cpair))
  
  !DO i=1,ncol
  !   T(i,1:10)  = 260. !K
  !   T(i,11:)  = 260.+ 3.*7.5*( log(pref_edge(11:pver)/p00)-log(pref_edge(11)/p00) ) !K
  !end do

    ! Profiles of background state variables
  call gw_prof(ncol, p, cpair, T, rhoi, nm, ni)



  piln=log( pint )
  rhom = p%mid / (Rgas * T )


! Kluge in constant profiles of rho and N
!----------------------------------------
    rhom=1._r8
    rhoi=1._r8
    nm = 0.0175_r8
    ni = 0.0175_r8
!-----------------------------------------

     write(*,*) " RHOM ", shape(rhom)
     write(*,*) " RHOM ", minval(rhom),maxval(rhom)

  zi =0.
  DO l=pver,1,-1
    zi(:,l) = zi(:,l+1) +p%del(:,l)/rhom(:,l) 
  end do

  zi = zi / 9.80616_r8
  zm = 0.5*( zi(:,2:pver+1) +  zi(:,1:pver) )


 do L=1,pver
    write(*,*) nm(1,L), T(1,L), th(1,L) , zm(1,L)
 end do

 u = 10.
 v = 0.
 dse = 0.


#if 0 
do i=1,ncol
    where(zm(i,:)<2500. )
       u(i,:) = 10.*zm(i,:)/2500.
    end where
 end do
#endif

!** nonhydrostatic
!  Add shear to test nonhydro
!-----------------------------
#if 0
do i=1,ncol
    where(zm(i,:)<5000. )
       u(i,:) = 20.*zm(i,:)/2500.
    end where
    where(zm(i,:)>=5000. )
       u(i,:) = 40.
    end where
 end do
#endif

 gbxar = 0.000275_r8 
! Convert grid box areas from Steradians to km+2
gbxar = gbxar*( SHR_CONST_REARTH/1000._r8)*( SHR_CONST_REARTH/1000._r8)


! -- Ridge width and length (km)
hwdth = 80._r8

!** nonhydrostatic
!  narrow ridge to test nonhydro
!-----------------------------
!hwdth = 2._r8

clngt = 160._r8

angll = 0._r8
anixy = 1.0_r8 !0.8_r8

do i=1,ncol
   mxdis(i,:) = 3000.*( 1._r8 - 1._r8 * (i-1)/ncol )
end do
sgh=mxdis(:,1)
landfrac(:)=1._r8

     ! Stupid double setting of parameters
     !------------------------------------
     effgw_rdg = 1.0_r8
     effgw_rdg_max = 1.0_r8
     rdg_cd_llb = 1.0_r8
     rdg_beta_cd_llb = 3.0_r8 !1.0_r8
     nrdg_beta =1


     trpd_leewv_rdg_beta = .FALSE.
!** nonhydrostatic
     !trpd_leewv_rdg_beta = .TRUE.
     effgw_rdg_beta =  effgw_rdg
     effgw_rdg_beta_max =  effgw_rdg_max
     n_rdg_beta = nrdg_beta

     ! nml settings for old GW scheme
     !-------------------------------
     effgw_oro = 0.125D0

     kvtt = 0.000001_r8
    
    lchnk = 1
    nlon = ncol
    nlat = 1
    lons = 100._r8 
    lats=45._r8

    write(511) nlon,nlat,pver,n_rdg_beta,itime,ntim
    write(511) lons,lats,lats


    !use_ridge_param = .true.
    if (use_ridge_param) then
     call gw_rdg_calc(&
        'BETA ', ncol , lchnk, n_rdg_beta, dt,     &
        u , v , t, &
        p, piln, zm, zi,                          &
        nm, ni, rhoi, kvtt, q, dse,               &
        effgw_rdg_beta, effgw_rdg_beta_max,       &
        hwdth, clngt, gbxar, mxdis, angll, anixy, &
        rdg_beta_cd_llb, trpd_leewv_rdg_beta,     &
        ptend, flx_heat)
    endif
#if 1
    if (use_isotropic_param) then
     call gw_oro_calc( &
      ncol, lchnk, dt, effgw_oro, &
      u, v, t, &
      p, piln, zm, zi, &
      nm, ni, rhoi, kvtt, q, dse, &
      sgh, landfrac,lats, &
      ptend, flx_heat)
    endif
#endif


end program test1

