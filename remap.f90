module remap_mod

  use MOM_remapping, only : remapping_CS
  use MOM_remapping, only : initialize_remapping
  use MOM_remapping, only : remapping_core_h
  use MOM_remapping, only : remapping_set_param
  use MOM_error_handler, only : MOM_error, FATAL

  implicit none

  private

  public :: remap

contains

  function remap(u_in,zi,zo,method,bndy_extrapolation,missing)
    real(kind=8), dimension(:,:,:), intent(in) :: u_in
    real(kind=8), dimension(:,:,:), intent(in) :: zi,zo
    logical, intent(in) :: bndy_extrapolation
    real(kind=8), intent(in), optional :: missing
    character(len=*), intent(in) :: method
    real(kind=8), dimension(size(u_in,1),size(u_in,2),size(zo,3)-1) :: remap


    real, parameter :: epsln=1.e-10
    real, parameter :: min_thickness=1.e-9

    type(remapping_CS) :: CS
    real(kind=8), dimension(size(zi,3)-1) :: h0
    real(kind=8), dimension(size(zi,3)-1) :: u0
    real(kind=8), dimension(size(zo,3)-1) :: h1
    real(kind=8), dimension(size(zo,3)-1) :: u1


    integer :: nz,ni,nj,i,j,k,imethod, degree,nz2, n1, n2
    integer :: ierr

    ni=size(u_in,1);nj=size(u_in,2);nz=size(u_in,3);nz2=size(zo,3)-1

    !print *,'u_in size= ',size(u_in)
    !print *,'zi size= ',size(zi)
    !print *,'zo size= ',size(zo)    
    
    if (ni /= size(zi,1) .or. nj .ne. size(zi,2) .or. nz /= size(zi,3)-1) call MOM_error(FATAL,'size mismatch u_in/zi')
    if (ni /= size(zo,1) .or. nj .ne. size(zo,2)) call MOM_error(FATAL,'size mismatch u_in/zo')

    call initialize_remapping(CS, method, boundary_extrapolation=bndy_extrapolation )

    do j=1,nj
      do i=1,ni
        do k=1,nz
          h0(k)=zi(i,j,k)-zi(i,j,k+1)
          u0(k)=u_in(i,j,k)
        enddo
        do k=1,nz2
          h1(k)=zo(i,j,k)-zo(i,j,k+1)
        enddo
        call remapping_core_h(CS, nz,h0,u0,nz2,h1,u1)
        do k=1,nz2
          remap(i,j,k)=u1(k)
        enddo
      enddo
    enddo

    return

  end function remap



end module remap_mod
