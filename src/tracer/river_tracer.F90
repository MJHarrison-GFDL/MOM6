!> A tracer package to track river inputs.
module river_tracer

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,            only : EFP_type
use MOM_coupler_types,   only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator,   only : diag_ctrl
use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_grid,            only : ocean_grid_type
use MOM_hor_index,       only : hor_index_type
use MOM_io,              only : file_exists, MOM_read_data, slasher
use MOM_io,              only : vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : query_initialized, set_initialized, MOM_restart_CS
use MOM_spatial_means,   only : global_mass_int_EFP
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type, time_type_to_real
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_river_tracer, initialize_river_tracer
public river_tracer_column_physics
public river_stock, river_tracer_end

type :: stringType
   character(len=:), allocatable :: string
end type stringType

type :: maskType2D
   integer, allocatable, dimension(:,:) :: mask
end type maskType2D

!> The control structure for the river tracer package
type, public :: river_tracer_CS ; private
   logical :: rivers_are_3d = .false. !< If true, river discharge can occur below the surface.
                                     !! (Option currently not suppported)
  integer :: num_rivers = 0  !< The number of rivers in the model domain
  type(maskType2D)                    :: mask2D  !< A 2-dimensional map of river identifiers
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM tracer registry
  real, pointer :: tr(:,:,:,:) => NULL() !< The array of tracers used in this subroutine, [R L-3] ~> [kg m-3]
  logical :: river_may_reinit  !< If true, river tracers may be reset by the initialization code
 !! if they are not found in the restart files.
  type(vardesc), pointer :: tr_desc(:) => NULL() !< Descriptions and metadata for the tracer
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure
  character(len=128) :: map_file !< The filename of the river map
  character(len=32), allocatable, dimension(:) :: names !< A list of river names
end type river_tracer_CS

contains

!> Register river tracer fields and subroutines to be used with MOM.
function register_river_tracer(G, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(ocean_grid_type),       intent(in) :: G   !< A horizontal index type structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(river_tracer_CS),        pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                 !! structure for the tracer advection and
                                                 !! diffusion module
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct
  ! Local variables
  character(len=40)  :: mdl = "river_tracer" ! This module's name.
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=3)   :: name_tag ! String for creating identifying rivers
  real, pointer :: tr_ptr(:,:,:) => NULL() ! The river tracer concentration [kg m-3]
  logical :: register_river_tracer
  integer :: isd, ied, jsd, jed, nz, m
  character(len=64) :: filename, inputdir, flux_units
  real, allocatable, dimension(:,:) :: rmask ! A temporary array for storing the river mask

  isd = G%HI%isd ; ied = G%HI%ied ; jsd = G%HI%jsd ; jed = G%HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_river_tracer called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NUM_RIVER_PASSIVE_TRACERS", CS%num_rivers, &
                 "The number of rivers in the model domain.", units="none", &
                 fail_if_missing=.false.)
  if (CS%num_rivers > 0) then
    allocate( CS%names(CS%num_rivers))
    allocate( CS%tr_desc(CS%num_rivers))
    allocate(CS%mask2d%mask(isd:ied,jsd:jed))
    allocate(rmask(isd:ied,jsd:jed))
    CS%mask2d%mask=0.0
    call get_param(param_file,mdl, "RIVER_MAP_FILENAME", CS%map_file,fail_if_missing=.true.)
    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
    filename=trim(slasher(inputdir))//trim(CS%map_file)
    call MOM_read_data(filename,'river_tracer_map', rmask, G%Domain)
    CS%mask2d%mask(:,:)=anint(rmask)
    deallocate(rmask)
    call get_param(param_file,mdl, "RIVER_NAMES", CS%names,fail_if_missing=.true.)
    do m=1, CS%num_rivers
      CS%tr_desc(m) = var_desc(trim(CS%names(m)), "kg m-3", trim(CS%names(m))//" River Tracer", caller=mdl)
    enddo
  else
     return
  endif

  ! This needs to be changed if the units of tracer are changed above.
  if (GV%Boussinesq) then ; flux_units = "kg s-1"
  else ; flux_units = "kg m-3 kg s-1" ; endif

  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%num_rivers), source=0.0)

  do m=1,CS%num_rivers
    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=CS%names(m), caller="register_river_tracer")
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, G%HI, GV, tr_desc=CS%tr_desc(m), &
                         registry_diags=.true., flux_units=flux_units, restart_CS=restart_CS, &
                         mandatory=.not.CS%river_may_reinit)
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_river_tracer = .true.

end function register_river_tracer

!> Initialize the river tracers and set up tracer output
subroutine initialize_river_tracer(restart, G, GV, US, h, diag, OBC, CS, &
                                  sponge_CSp)
  logical,                            intent(in) :: restart !< .true. if the fields have already
                                                         !! been read from a restart file.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US   !< The dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl),            target, intent(in) :: diag !< A structure that is used to regulate
                                                         !! diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(river_tracer_CS),                pointer    :: CS !< The control structure returned by a previous
                                                       !! call to register_river_tracer.
  type(sponge_CS),                    pointer    :: sponge_CSp !< Pointer to the control structure for the sponges.

  ! Local variables
  character(len=16) :: name
  logical :: OK
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB

  if (.not.associated(CS)) return
  if (CS%num_rivers < 1) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%diag => diag

  do m=1,CS%num_rivers
    call query_vardesc(CS%tr_desc(m), name=name, caller="initialize_river_tracer")
    if ((.not.restart) .or. (CS%river_may_reinit .and. .not. &
        query_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp))) then
       CS%tr(:,:,:,m)=0.0
       call set_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp)
    endif
  enddo


  if (associated(OBC)) then
  ! Put something here...
  endif

end subroutine initialize_river_tracer

!> Apply sources, sinks, diapycnal mixing and rising motions to the river tracers
subroutine river_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, tv, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(river_tracer_CS),     pointer    :: CS   !< The control structure returned by a previous
                                              !! call to register_river_tracer.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]
!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  real :: vol_scale ! A conversion factor for volumes into m3 [m3 H-1 L-2 ~> 1 or m3 kg-1]
  real :: h_total   ! A running sum of thicknesses [H ~> m or kg m-2]

  integer :: i, j, k, is, ie, js, je, nz, m, k_max
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (CS%num_rivers < 1) return

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,CS%num_rivers
      do k=1,nz ;do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  else
    do m=1,CS%num_rivers
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

  do m=2,CS%num_rivers
    do j=js,je ; do i=is,ie
       CS%tr(i,j,1,m) = CS%tr(i,j,1,m)+(fluxes%lrunoff(i,j)+fluxes%frunoff(i,j))*dt/ &
                         (vol_scale * (h_new(i,j,1)+GV%H_subroundoff) * G%areaT(i,j) )
    enddo ; enddo
  enddo


end subroutine river_tracer_column_physics

!> Calculate the mass-weighted integral of the river tracer stocks, returning the number of stocks it
!! has calculated.  If the stock_index is present, only the stock corresponding to that coded index is returned.
function river_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(EFP_type), dimension(:),       intent(out)   :: stocks !< The mass-weighted integrated amount of each
                                                            !! tracer, in kg times concentration units [kg conc]
  type(river_tracer_CS),                pointer       :: CS   !< The control structure returned by a previous
                                                            !! call to register_river_tracer.
  character(len=*), dimension(:),     intent(out)   :: names  !< the names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units  !< the units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< the coded index of a specific stock
                                                                   !! being sought.
  integer                                           :: river_stock !< The number of stocks calculated here.

  ! Local variables
  integer :: m

  river_stock = 0
  if (.not.associated(CS)) return
  if (CS%num_rivers < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%num_rivers
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="river_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = global_mass_int_EFP(h, G, GV, CS%tr(:,:,:,m), on_PE_only=.true.)
  enddo
  river_stock = CS%num_rivers

end function river_stock

!> Deallocate memory associated with this tracer package
subroutine river_tracer_end(CS)
  type(river_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                     !! call to register_river_tracer.

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine river_tracer_end

!> \namespace river_tracer
!!
!!  By Matthew Harrison, Dmitry Dhukovsky and Theresa Cordero

end module river_tracer
