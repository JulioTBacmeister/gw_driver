!-------------------------------------------------------------------------------
!physics data types module
!-------------------------------------------------------------------------------
module physics_types

  use shr_kind_mod,      only: r8 => shr_kind_r8
  use gw_common,         only: pver
  use constituent,       only: pcnst

  implicit none
  private          ! Make default type private to the module

  logical, parameter :: adjust_te = .FALSE.

! Public types:

  public physics_state
  public physics_tend
  public physics_ptend
  public physics_ptend_init

  

!-------------------------------------------------------------------------------
  type physics_state
     integer                                     :: &
          lchnk,                &! chunk index
          ngrdcol,              &! -- Grid        -- number of active columns (on the grid)
          psetcols=0,           &! --             -- max number of columns set - if subcols = pcols*psubcols, else = pcols
          ncol=0                 ! --             -- sum of nsubcol for all ngrdcols - number of active columns
     real(r8), dimension(:), allocatable         :: &
          lat,     &! latitude (radians)
          lon,     &! longitude (radians)
          ps,      &! surface pressure
          psdry,   &! dry surface pressure
          phis,    &! surface geopotential
          ulat,    &! unique latitudes  (radians)
          ulon      ! unique longitudes (radians)
     real(r8), dimension(:,:),allocatable        :: &
          t,       &! temperature (K)
          u,       &! zonal wind (m/s)
          v,       &! meridional wind (m/s)
          s,       &! dry static energy
          omega,   &! vertical pressure velocity (Pa/s) 
          pmid,    &! midpoint pressure (Pa) 
          pmiddry, &! midpoint pressure dry (Pa) 
          pdel,    &! layer thickness (Pa)
          pdeldry, &! layer thickness dry (Pa)
          rpdel,   &! reciprocal of layer thickness (Pa)
          rpdeldry,&! recipricol layer thickness dry (Pa)
          lnpmid,  &! ln(pmid)
          lnpmiddry,&! log midpoint pressure dry (Pa) 
          exner,   &! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
          zm        ! geopotential height above surface at midpoints (m)

     real(r8), dimension(:,:,:),allocatable      :: &
          q         ! constituent mixing ratio (kg/kg moist or dry air depending on type)

     real(r8), dimension(:,:),allocatable        :: &
          pint,    &! interface pressure (Pa)
          pintdry, &! interface pressure dry (Pa) 
          lnpint,  &! ln(pint)
          lnpintdry,&! log interface pressure dry (Pa) 
          zi        ! geopotential height above surface at interfaces (m)

     real(r8), dimension(:),allocatable          :: &
          te_ini,  &! vertically integrated total (kinetic + static) energy of initial state
          te_cur,  &! vertically integrated total (kinetic + static) energy of current state
          tw_ini,  &! vertically integrated total water of initial state
          tw_cur    ! vertically integrated total water of new state
     integer :: count ! count of values with significant energy or water imbalances
     integer, dimension(:),allocatable           :: &
          latmapback, &! map from column to unique lat for that column
          lonmapback, &! map from column to unique lon for that column
          cid        ! unique column id
     integer :: ulatcnt, &! number of unique lats in chunk
                uloncnt   ! number of unique lons in chunk

  end type physics_state

!-------------------------------------------------------------------------------
  type physics_tend

     integer   ::   psetcols=0 ! max number of columns set- if subcols = pcols*psubcols, else = pcols

     real(r8), dimension(:,:),allocatable        :: dtdt, dudt, dvdt
     real(r8), dimension(:),  allocatable        :: flx_net
     real(r8), dimension(:),  allocatable        :: &
          te_tnd,  &! cumulative boundary flux of total energy
          tw_tnd    ! cumulative boundary flux of total water
  end type physics_tend

!-------------------------------------------------------------------------------
! This is for tendencies returned from individual parameterizations
  type physics_ptend

     integer   ::   psetcols=0 ! max number of columns set- if subcols = pcols*psubcols, else = pcols

     character*24 :: name    ! name of parameterization which produced tendencies.

     logical ::             &
          ls = .false.,               &! true if dsdt is returned
          lu = .false.,               &! true if dudt is returned
          lv = .false.                 ! true if dvdt is returned

     logical,dimension(pcnst) ::  lq = .false.  ! true if dqdt() is returned
     !logical,dimension(1) ::  lq = .false.  ! true if dqdt() is returned

     integer ::             &
          top_level,        &! top level index for which nonzero tendencies have been set
          bot_level          ! bottom level index for which nonzero tendencies have been set

     real(r8), dimension(:,:),allocatable   :: &
          s,                &! heating rate (J/kg/s)
          u,                &! u momentum tendency (m/s/s)
          v                  ! v momentum tendency (m/s/s)
     real(r8), dimension(:,:,:),allocatable :: &
          q                  ! consituent tendencies (kg/kg/s)

! boundary fluxes
     real(r8), dimension(:),allocatable     ::&
          hflux_srf,     &! net heat flux at surface (W/m2)
          hflux_top,     &! net heat flux at top of model (W/m2)
          taux_srf,      &! net zonal stress at surface (Pa)
          taux_top,      &! net zonal stress at top of model (Pa)
          tauy_srf,      &! net meridional stress at surface (Pa)
          tauy_top        ! net meridional stress at top of model (Pa)
     real(r8), dimension(:,:),allocatable   ::&
          cflx_srf,      &! constituent flux at surface (kg/m2/s)
          cflx_top        ! constituent flux top of model (kg/m2/s)

  end type physics_ptend


contains
!===============================================================================
  subroutine physics_ptend_init(ptend, psetcols, name, ls, lu, lv, lq)
!-----------------------------------------------------------------------
! Allocate the fields in the structure which are specified.
! Initialize the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(out)    :: ptend    ! Parameterization tendencies
    integer, intent(in)                 :: psetcols ! maximum number of columns
    character(len=*)                    :: name     ! optional name of parameterization which produced tendencies.
    logical, optional                   :: ls       ! if true, then fields to support dsdt are allocated
    logical, optional                   :: lu       ! if true, then fields to support dudt are allocated
    logical, optional                   :: lv       ! if true, then fields to support dvdt are allocated
    logical, dimension(pcnst),optional  :: lq       ! if true, then fields to support dqdt are allocated
    
!-----------------------------------------------------------------------
    if (allocated(ptend%s)) then 
       call endrun(' physics_ptend_init: ptend should not be allocated before calling this routine')
    end if

    ptend%name     = name
    ptend%psetcols =  psetcols

    ! If no fields being stored, initialize all values to appropriate nulls and return
    if (.not. present(ls) .and. .not. present(lu) .and. .not. present(lv) .and. .not. present(lq) ) then
       ptend%ls       = .false.
       ptend%lu       = .false.
       ptend%lv       = .false.
       ptend%lq(:)    = .false.
       ptend%top_level = 1
       ptend%bot_level = pver
       return
    end if

    if (present(ls)) then
       ptend%ls = ls
    else
       ptend%ls = .false.
    end if

    if (present(lu)) then
       ptend%lu = lu
    else
       ptend%lu = .false.
    end if

    if (present(lv)) then
       ptend%lv = lv
    else
       ptend%lv = .false.
    end if

    if (present(lq)) then
       ptend%lq(:) = lq(:)
    else
       ptend%lq(:) = .false.
    end if

    call physics_ptend_alloc(ptend, psetcols)

    call physics_ptend_reset(ptend)

    return
  end subroutine physics_ptend_init

!===============================================================================


!===============================================================================

  subroutine physics_ptend_reset(ptend)
!-----------------------------------------------------------------------
! Reset the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies
!-----------------------------------------------------------------------
    integer :: m             ! Index for constiuent
!-----------------------------------------------------------------------

    if(ptend%ls) then
       ptend%s = 0._r8
       ptend%hflux_srf = 0._r8
       ptend%hflux_top = 0._r8
    endif
    if(ptend%lu) then
       ptend%u = 0._r8
       ptend%taux_srf = 0._r8
       ptend%taux_top = 0._r8
    endif
    if(ptend%lv) then
       ptend%v = 0._r8
       ptend%tauy_srf = 0._r8
       ptend%tauy_top = 0._r8
    endif
    if(any (ptend%lq(:))) then
       ptend%q = 0._r8
       ptend%cflx_srf = 0._r8
       ptend%cflx_top = 0._r8
    end if

    ptend%top_level = 1
    ptend%bot_level = pver

    return
  end subroutine physics_ptend_reset

!===============================================================================

!===============================================================================

subroutine physics_ptend_alloc(ptend,psetcols)

! allocate the individual ptend components

  type(physics_ptend), intent(inout) :: ptend

  integer, intent(in)                :: psetcols

  integer :: ierr = 0

  ptend%psetcols = psetcols

  if (ptend%ls) then
     allocate(ptend%s(psetcols,pver), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%s')

     allocate(ptend%hflux_srf(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%hflux_srf')

     allocate(ptend%hflux_top(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%hflux_top')
  end if

  if (ptend%lu) then 
     allocate(ptend%u(psetcols,pver), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%u')

     allocate(ptend%taux_srf(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%taux_srf')

     allocate(ptend%taux_top(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%taux_top')
  end if

  if (ptend%lv) then 
     allocate(ptend%v(psetcols,pver), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%v')

     allocate(ptend%tauy_srf(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%tauy_srf')

     allocate(ptend%tauy_top(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%tauy_top')
  end if

  if (any(ptend%lq)) then 
     allocate(ptend%q(psetcols,pver,pcnst), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%q')

     allocate(ptend%cflx_srf(psetcols,pcnst), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%cflx_srf')

     allocate(ptend%cflx_top(psetcols,pcnst), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%cflx_top')
  end if

end subroutine physics_ptend_alloc

!===============================================================================


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine endrun(msg)

   integer :: iulog

   character(len=*), intent(in), optional :: msg    ! string to be printed

    iulog=6

   if (present (msg)) then
      write(iulog,*)'ENDRUN:', msg
   else
      write(iulog,*)'ENDRUN: called without a message string'
   end if

   stop

end subroutine endrun

end module physics_types
