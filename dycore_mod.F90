!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fake dycore module 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module dycore_mod
  implicit none
  private

  public :: SE_dycore
  public :: FV_dycore
  public :: read_dynamics_state
  public :: read_convective_heating
  public :: read_ridge_file
  

  logical, parameter :: SE_dycore=.FALSE.
  logical, parameter :: FV_dycore=.TRUE.

contains
  
!------------------------------------
subroutine read_dynamics_state (file_name , ai, bi, p00, u, v, t, q, ps, lons, lats, landfrac )
use shr_kind_mod,   only: r8=>shr_kind_r8
#include <netcdf.inc>

  character(len=*), intent(in) :: file_name

  real(r8)               , intent(out) :: p00
  real(r8)  , allocatable, intent(out) :: ai(:)
  real(r8)  , allocatable, intent(out) :: bi(:)
  real(r8)  , allocatable, intent(out) :: u(:,:,:)
  real(r8)  , allocatable, intent(out) :: v(:,:,:)
  real(r8)  , allocatable, intent(out) :: ps(:,:)
  real(r8)  , allocatable, intent(out) :: t(:,:,:)
  real(r8)  , allocatable, intent(out) :: q(:,:,:)
  real(r8)  , allocatable, intent(out) :: lons(:)
  real(r8)  , allocatable, intent(out) :: lats(:)
  real(r8)  , allocatable, intent(out) :: landfrac(:)

  integer  :: xpver, ncol, nlat, nlon, ntim

  ! Vars needed by NetCDF operators
  integer  :: ncid, dimid, varid, status


  status = nf_open(file_name, 0, ncid)
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

  if (FV_dycore)  ncol=nlat*nlon
  
  allocate(  landfrac(ncol) )


  landfrac(:) = 1.0

  write(*,*) "SHAPE Q ", shape( q )
  q=0.


end subroutine read_dynamics_state

!------------------------------------
subroutine read_convective_heating (file_name , netdt   )
use shr_kind_mod,   only: r8=>shr_kind_r8
#include <netcdf.inc>

  character(len=*), intent(in) :: file_name


  real(r8)  , allocatable, intent(out) :: netdt(:,:,:)

  integer  :: xpver, ncol, nlat, nlon, ntim


  ! Vars needed by NetCDF operators
  integer  :: ncid, dimid, varid, status



  status = nf_open(file_name, 0, ncid)
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

   

     ! Convective heating
     allocate(  netdt(ncol,xpver,ntim) )

  status = NF_INQ_VARID(ncid, 'NETDT', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, netdt )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
     write(*,*) " NETDT ", shape(netdt)
     write(*,*) " NETDT ", minval(netdt),maxval(netdt)


  status = nf_close (ncid)
  if (status .ne. NF_NOERR) call handle_err(status)          


end subroutine read_convective_heating

!=================================================

subroutine read_ridge_file( file_name , hwdth, clngt, gbxar, mxdis, angll, anixy, sgh ) ! , band )
use shr_kind_mod,   only: r8=>shr_kind_r8
use shr_const_mod,  only: shr_const_rearth
#include <netcdf.inc>

  character(len=*), intent(in) :: file_name

  ! Ridge parameters
  real(r8)  , allocatable, intent(out) :: mxdis(:,:)
  real(r8)  , allocatable, intent(out) :: hwdth(:,:)
  real(r8)  , allocatable, intent(out) :: clngt(:,:)
  real(r8)  , allocatable, intent(out) :: angll(:,:)
  real(r8)  , allocatable, intent(out) :: anixy(:,:)
  real(r8)  , allocatable, intent(out) :: gbxar(:)
  real(r8)  , allocatable, intent(out) :: sgh(:)

  integer  :: nlon,nlat,ncol,nrdgs
  
  ! Vars needed by NetCDF operators
  integer  :: ncid, dimid, varid, status


!===========================================================================


!=============================================================================
  status = nf_open( file_name , 0, ncid)
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

      allocate( mxdis(ncol,nrdgs), angll(ncol,nrdgs), anixy(ncol,nrdgs), &
                hwdth(ncol,nrdgs), clngt(ncol,nrdgs), gbxar(ncol), sgh(ncol)  )
     



  status = NF_INQ_VARID(ncid, 'MXDIS', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, mxdis )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

#if 0
  status = NF_INQ_VARID(ncid, 'WGHTS', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, wghts )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = NF_INQ_VARID(ncid, 'ANISO', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, aniso )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
#endif
     
  status = NF_INQ_VARID(ncid, 'ANIXY', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, anixy )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = NF_INQ_VARID(ncid, 'HWDTH', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, hwdth )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = NF_INQ_VARID(ncid, 'CLNGT', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, clngt )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = NF_INQ_VARID(ncid, 'ANGLL', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, angll )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
 
  status = NF_INQ_VARID(ncid, 'SGH', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, sgh )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = NF_INQ_VARID(ncid, 'GBXAR', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, gbxar )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = nf_close (ncid)
  if (status .ne. NF_NOERR) call handle_err(status)          

  
  ! Convert grid box areas from Steradians to km+2
  gbxar = gbxar*( SHR_CONST_REARTH/1000._r8)*( SHR_CONST_REARTH/1000._r8)

end subroutine read_ridge_file

end module dycore_mod
