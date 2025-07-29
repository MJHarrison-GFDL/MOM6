!> Regrid columns for the HyCOM coordinate
module coord_hycom_2d

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_remapping,     only : remapping_CS, remapping_core_h
use MOM_EOS,           only : EOS_type, calculate_density
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid, regridding_set_ppolys
use regrid_interp,     only : DEGREE_MAX

implicit none ; private

!> Control structure containing required parameters for the HyCOM coordinate
type, public :: hycom_2D_CS !; private

  !> Number of layers/levels in generated grid
  integer :: nk

  !> Number of Hycom1 grid types
  integer :: ng

  !> Nominal near-surface resolution [Z ~> m]
  real, allocatable, dimension(:,:) :: coordinateResolution

  !> Nominal density of interfaces [R ~> kg m-3]
  real, allocatable, dimension(:,:) :: target_density

  !> A 2-dimensonal (i,j) identifier for the grid type (1..ng) [nondim]
  real, allocatable, dimension(:,:) :: rmask

  !> Maximum depths of interfaces [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: max_interface_depths

  !> Maximum thicknesses of layers [H ~> m or kg m-2]
  real, allocatable, dimension(:) :: max_layer_thickness

  !> If true, an interface only moves if it improves the density fit
  logical :: only_improves = .false.

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS
end type hycom_2D_CS

public init_coord_hycom_2d, set_hycom_2d_params, build_hycom_2d_column, end_coord_hycom_2d

contains

!> Initialise a hycom_2d_CS with pointers to parameters
subroutine init_coord_hycom_2d(CS, nk, ng, coordinateResolution, target_density, interp_CS)
  type(hycom_2d_CS),       pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,              intent(in) :: nk !< Number of layers in generated grid
  integer,              intent(in) :: ng !< Number of unique grid constructors
  real, dimension(nk,ng),  intent(in) :: coordinateResolution !< Nominal near-surface resolution [Z ~> m]
  real, dimension(nk+1,ng),intent(in) :: target_density !< Interface target densities [R ~> kg m-3]
  type(interp_CS_type), intent(in) :: interp_CS !< Controls for interpolation

  if (associated(CS)) call MOM_error(FATAL, "init_coord_hycom: CS already associated!")
  allocate(CS)
  allocate(CS%coordinateResolution(nk,ng))
  allocate(CS%target_density(nk+1,ng))

  CS%nk                      = nk
  CS%ng                      = ng
  if (size(coordinateResolution,1) /= CS%nk) &
       call MOM_error(FATAL, "set_hycom_2d_params: coordinateResolution inconsistent size")
  if (size(coordinateResolution,2) /= CS%ng) &
       call MOM_error(FATAL, "set_hycom_2d_params: coordinateResolution inconsistent size")
  CS%coordinateResolution(:,:) = coordinateResolution(:,:)
  if (size(target_density,1) /= CS%nk+1) &
       call MOM_error(FATAL, "set_hycom_2d_params: target_density inconsistent size")
  if (size(target_density,2) /= CS%ng) &
       call MOM_error(FATAL, "set_hycom_2d_params: target_density inconsistent size")
  CS%target_density(:,:)       = target_density(:,:)
  CS%interp_CS               = interp_CS

end subroutine init_coord_hycom_2d

!> This subroutine deallocates memory in the control structure for the coord_hycom module
subroutine end_coord_hycom_2d(CS)
  type(hycom_2d_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return
  deallocate(CS%coordinateResolution)
  deallocate(CS%target_density)
  if (allocated(CS%max_interface_depths)) deallocate(CS%max_interface_depths)
  if (allocated(CS%max_layer_thickness)) deallocate(CS%max_layer_thickness)
  deallocate(CS)
end subroutine end_coord_hycom_2d

!> This subroutine can be used to set the parameters for the coord_hycom module
subroutine set_hycom_2d_params(CS, max_interface_depths, max_layer_thickness, only_improves, interp_CS)
  type(hycom_2d_CS),                 pointer    :: CS !< Coordinate control structure
  real, dimension(:),   optional, intent(in) :: max_interface_depths !< Maximum depths of interfaces [H ~> m or kg m-2]
  real, dimension(:),   optional, intent(in) :: max_layer_thickness  !< Maximum thicknesses of layers [H ~> m or kg m-2]
  logical, optional, intent(in) :: only_improves !< If true, an interface only moves if it improves the density fit
  type(interp_CS_type), optional, intent(in) :: interp_CS !< Controls for interpolation

  if (.not. associated(CS)) call MOM_error(FATAL, "set_hycom_params: CS not associated")

  if (present(max_interface_depths)) then
    if (size(max_interface_depths,1) /= CS%nk+1) &
         call MOM_error(FATAL, "set_hycom_params: max_interface_depths inconsistent size")
    allocate(CS%max_interface_depths(CS%nk+1))
    CS%max_interface_depths(:) = max_interface_depths(:)
  endif

  if (present(max_layer_thickness)) then
    if (size(max_layer_thickness,1) /= CS%nk) &
         call MOM_error(FATAL, "set_hycom_params: max_layer_thickness inconsistent size")
    allocate(CS%max_layer_thickness(CS%nk))
    CS%max_layer_thickness(:) = max_layer_thickness(:)
  endif

  if (present(only_improves)) CS%only_improves = only_improves

  if (present(interp_CS)) CS%interp_CS = interp_CS
end subroutine set_hycom_2d_params

!> Build a HyCOM coordinate column
subroutine build_hycom_2d_column(CS, rmask, remapCS, eqn_of_state, nz, depth, h, T, S, p_col, &
                               z_col, z_col_new, zScale, h_neglect, h_neglect_edge)
  type(hycom_2d_CS),        intent(in)    :: CS    !< Coordinate control structure
  real, dimension(CS%ng),   intent(in)    :: rmask  !< region mask
  type(remapping_CS),    intent(in)    :: remapCS !< Remapping parameters and options
  type(EOS_type),        intent(in)    :: eqn_of_state !< Equation of state structure
  integer,               intent(in)    :: nz    !< Number of levels
  real,                  intent(in)    :: depth !< Depth of ocean bottom (positive [H ~> m or kg m-2])
  real, dimension(nz),   intent(in)    :: T     !< Temperature of column [C ~> degC]
  real, dimension(nz),   intent(in)    :: S     !< Salinity of column [S ~> ppt]
  real, dimension(nz),   intent(in)    :: h     !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz),   intent(in)    :: p_col !< Layer pressure [R L2 T-2 ~> Pa]
  real, dimension(nz+1), intent(in)    :: z_col !< Interface positions relative to the surface [H ~> m or kg m-2]
  real, dimension(CS%nk+1), intent(inout) :: z_col_new !< Absolute positions of interfaces [H ~> m or kg m-2]
  real, optional,        intent(in)    :: zScale !< Scaling factor from the input coordinate thicknesses in [Z ~> m]
                                                !! to desired units for zInterface, perhaps GV%Z_to_H in which
                                                !! case this has units of [H Z-1 ~> nondim or kg m-3]
  real,                  intent(in)    :: h_neglect !< A negligibly small width for the purpose of
                                                !! cell reconstruction [H ~> m or kg m-2]
  real,        optional, intent(in)    :: h_neglect_edge !< A negligibly small width for the purpose of
                                                !! edge value calculation [H ~> m or kg m-2]

  ! Local variables
  integer   :: k, i, j
  real, dimension(nz)      :: rho_col   ! Layer densities in a column [R ~> kg m-3]
  real, dimension(CS%nk)   :: h_col_new ! New layer thicknesses [H ~> m or kg m-2]
  real, dimension(CS%nk)   :: h1_col_new ! New layer thicknesses [H ~> m or kg m-2]
  real, dimension(CS%nk+1)   :: z1_col_new ! New layer interfacel positions relative to the surface [H ~> m or kg m-2]
  real, dimension(CS%nk)   :: r_col_new ! New layer densities [R ~> kg m-3]
  real, dimension(CS%nk)   :: T_col_new ! New layer temperatures [C ~> degC]
  real, dimension(CS%nk)   :: S_col_new ! New layer salinities [S ~> ppt]
  real, dimension(CS%nk)   :: p_col_new ! New layer pressure [R L2 T-2 ~> Pa]
  real, dimension(CS%nk+1) :: RiA_ini   ! Initial nk+1 interface density anomaly w.r.t. the
                                        ! interface target densities [R ~> kg m-3]
  real, dimension(CS%nk+1) :: RiA_new   ! New interface density anomaly w.r.t. the
                                        ! interface target densities [R ~> kg m-3]
  real :: z_1, z_nz  ! mid point of 1st and last layers [H ~> m or kg m-2]
  real :: z_scale    ! A scaling factor from the input thicknesses to the target thicknesses,
                     ! perhaps 1 or a factor in [H Z-1 ~> 1 or kg m-3]
  real :: stretching ! z* stretching, converts z* to z [nondim].
  real :: nominal_z ! Nominal depth of interface when using z* [H ~> m or kg m-2]
  logical :: maximum_depths_set ! If true, the maximum depths of interface have been set.
  logical :: maximum_h_set      ! If true, the maximum layer thicknesses have been set.
  real :: wt, coord_res

  maximum_depths_set = allocated(CS%max_interface_depths)
  maximum_h_set = allocated(CS%max_layer_thickness)

  z_scale = 1.0 ; if (present(zScale)) z_scale = zScale

  ! Work bottom recording potential density
  call calculate_density(T, S, p_col, rho_col, eqn_of_state)
  ! This ensures the potential density profile is monotonic
  ! although not necessarily single valued.
  do k = nz-1, 1, -1
    rho_col(k) = min( rho_col(k), rho_col(k+1) )
  enddo


  ! Interpolates for the target interface position with the rho_col profile
  ! Based on global density profile, interpolate to generate a new grid
  h_col_new(:)=0.0;z_col_new(:)=0.0
  do k=1,CS%ng
    wt=rmask(k)
    if (wt>0.0) then
      call build_and_interpolate_grid(CS%interp_CS, rho_col, nz, h(:), z_col, &
           CS%target_density(:,k), CS%nk, h1_col_new, z1_col_new, h_neglect, h_neglect_edge)
      h_col_new(:)=wt*h1_col_new(:)+h_col_new(:)
      z_col_new(:)=wt*z1_col_new(:)+z_col_new(:)
    endif
  enddo


  ! Sweep down the interfaces and make sure that the interface is at least
  ! as deep as a nominal target z* grid
  nominal_z = 0.
  stretching = z_col(nz+1) / depth ! Stretches z* to z
  do k = 2, CS%nk+1
    coord_res=0.0
    do i=1,CS%ng
      wt=rmask(i)
      coord_res = coord_res+wt*CS%coordinateResolution(k-1,i)
    enddo
    nominal_z = nominal_z + (z_scale * (coord_res * stretching))
    z_col_new(k) = max( z_col_new(k), nominal_z )
    z_col_new(k) = min( z_col_new(k), z_col(nz+1) )
  enddo

  if (maximum_depths_set .and. maximum_h_set) then ; do k=2,CS%nk
    ! The loop bounds are 2 & nz so the top and bottom interfaces do not move.
    ! Recall that z_col_new is positive downward.
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K), &
                       z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; elseif (maximum_depths_set) then ; do K=2,CS%nk
    z_col_new(K) = min(z_col_new(K), CS%max_interface_depths(K))
  enddo ; elseif (maximum_h_set) then ; do k=2,CS%nk
    z_col_new(K) = min(z_col_new(K), z_col_new(K-1) + CS%max_layer_thickness(k-1))
  enddo ; endif
end subroutine build_hycom_2d_column


end module coord_hycom_2d
