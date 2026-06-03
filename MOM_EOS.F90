!> Provides subroutines for quantities specific to the equation of state
module MOM_EOS

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_string_functions, only : uppercase
use MOM_unit_scaling, only : unit_scale_type
use MOM_EOS_linear, only : calculate_density_linear
use MOM_EOS_Wright, only : calculate_density_wright, calculate_density_wright_full
use MOM_file_parser, only : param_file_type, get_param

implicit none ; private





!> A control structure for the equation of state
type, public :: EOS_type ; private
  integer :: form_of_EOS = 0 !< The equation of state to use.
  integer :: form_of_TFreeze = 0 !< The expression for the potential temperature
                             !! of the freezing point.
  logical :: EOS_quadrature  !< If true, always use the generic (quadrature)
                             !! code for the integrals of density.
  logical :: Compressible = .true. !< If true, in situ density is a function of pressure.
! The following parameters are used with the linear equation of state only.
  real :: Rho_T0_S0 !< The density at T=0, S=0 [kg m-3]
  real :: dRho_dT   !< The partial derivative of density with temperature [kg m-3 degC-1]
  real :: dRho_dS   !< The partial derivative of density with salinity [kg m-3 ppt-1]
! The following parameters are use with the linear expression for the freezing
! point only.
  real :: TFr_S0_P0 !< The freezing potential temperature at S=0, P=0 [degC]
  real :: dTFr_dS   !< The derivative of freezing point with salinity [degC ppt-1]
  real :: dTFr_dp   !< The derivative of freezing point with pressure [degC Pa-1]

! Unit conversion factors (normally used for dimensional testing but could also allow for
! change of units of arguments to functions)
  real :: m_to_Z = 1.      !< A constant that translates distances in meters to the units of depth.
  real :: kg_m3_to_R = 1.  !< A constant that translates kilograms per meter cubed to the units of density.
  real :: R_to_kg_m3 = 1.  !< A constant that translates the units of density to kilograms per meter cubed.
  real :: RL2_T2_to_Pa = 1.!< Convert pressures from R L2 T-2 to Pa.
  real :: L_T_to_m_s = 1.  !< Convert lateral velocities from L T-1 to m s-1.

!  logical :: test_EOS = .true. ! If true, test the equation of state
end type EOS_type

! The named integers that might be stored in eqn_of_state_type%form_of_EOS.
integer, parameter, public :: EOS_LINEAR = 1 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_UNESCO = 2 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_WRIGHT = 3 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_WRIGHT_FULL = 4 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_TEOS10 = 5 !< A named integer specifying an equation of state
integer, parameter, public :: EOS_NEMO   = 6 !< A named integer specifying an equation of state

character*(10), parameter :: EOS_LINEAR_STRING = "LINEAR" !< A string for specifying the equation of state
character*(10), parameter :: EOS_UNESCO_STRING = "UNESCO" !< A string for specifying the equation of state
character*(10), parameter :: EOS_WRIGHT_STRING = "WRIGHT" !< A string for specifying the equation of state
character*(12), parameter :: EOS_WRIGHT_FULL_STRING = "WRIGHT_FULL" !< A string for specifying the equation of state
character*(10), parameter :: EOS_TEOS10_STRING = "TEOS10" !< A string for specifying the equation of state
character*(10), parameter :: EOS_NEMO_STRING   = "NEMO"   !< A string for specifying the equation of state
character*(12), parameter :: EOS_DEFAULT = EOS_WRIGHT_FULL_STRING !< The default equation of state

integer, parameter :: TFREEZE_LINEAR = 1  !< A named integer specifying a freezing point expression
integer, parameter :: TFREEZE_MILLERO = 2 !< A named integer specifying a freezing point expression
integer, parameter :: TFREEZE_TEOS10 = 3  !< A named integer specifying a freezing point expression
character*(10), parameter :: TFREEZE_LINEAR_STRING = "LINEAR" !< A string for specifying the freezing point expression
character*(10), parameter :: TFREEZE_MILLERO_STRING = "MILLERO_78" !< A string for specifying
                                                              !! freezing point expression
character*(10), parameter :: TFREEZE_TEOS10_STRING = "TEOS10" !< A string for specifying the freezing point expression
character*(10), parameter :: TFREEZE_DEFAULT = TFREEZE_LINEAR_STRING !< The default freezing point expression



public calculate_density, EOS_init, EOS_allocate, EOS_end


interface calculate_density
  module procedure calculate_density_scalar, calculate_density_array, calculate_density_1d
end interface calculate_density

contains

  !> Initializes EOS_type by allocating and reading parameters
subroutine EOS_init(param_file, EOS, US)
  type(param_file_type), intent(in) :: param_file !< Parameter file structure
  type(EOS_type)    :: EOS !< Equation of state structure
  type(unit_scale_type), intent(in) :: US  !< A dimensional unit scaling type
  optional :: US
  ! Local variables
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_EOS" ! This module's name.
  character(len=40)  :: tmpstr

!  if (.not.associated(EOS)) call EOS_allocate(EOS)

  ! Read all relevant parameters and write them to the model log.
!  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "EQN_OF_STATE", tmpstr, &
                 "EQN_OF_STATE determines which ocean equation of state "//&
                 "should be used.  Currently, the valid choices are "//&
                 '"LINEAR", "UNESCO", "WRIGHT", "NEMO" and "TEOS10". '//&
                 "This is only used if USE_EOS is true.", default=EOS_DEFAULT)
  select case (uppercase(tmpstr))
    case (EOS_LINEAR_STRING)
      EOS%form_of_EOS = EOS_LINEAR
    case (EOS_UNESCO_STRING)
      EOS%form_of_EOS = EOS_UNESCO
    case (EOS_WRIGHT_STRING)
       EOS%form_of_EOS = EOS_WRIGHT
    case (EOS_WRIGHT_FULL_STRING)
      EOS%form_of_EOS = EOS_WRIGHT_FULL       
    case (EOS_TEOS10_STRING)
      EOS%form_of_EOS = EOS_TEOS10
    case (EOS_NEMO_STRING)
      EOS%form_of_EOS = EOS_NEMO
    case default
      call MOM_error(FATAL, "interpret_eos_selection: EQN_OF_STATE "//&
                              trim(tmpstr) // "in input file is invalid.")
  end select
!  call MOM_mesg('interpret_eos_selection: equation of state set to "' // &
!                trim(tmpstr)//'"', 5)

  if (EOS%form_of_EOS == EOS_LINEAR) then
    EOS%Compressible = .false.
    call get_param(param_file, mdl, "RHO_T0_S0", EOS%Rho_T0_S0, &
                 "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", "//&
                 "this is the density at T=0, S=0.", units="kg m-3", &
                 default=1000.0)
    call get_param(param_file, mdl, "DRHO_DT", EOS%dRho_dT, &
                 "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", "//&
                 "this is the partial derivative of density with "//&
                 "temperature.", units="kg m-3 K-1", default=-0.2)
    call get_param(param_file, mdl, "DRHO_DS", EOS%dRho_dS, &
                 "When EQN_OF_STATE="//trim(EOS_LINEAR_STRING)//", "//&
                 "this is the partial derivative of density with "//&
                 "salinity.", units="kg m-3 PSU-1", default=0.8)
  endif

  call get_param(param_file, mdl, "EOS_QUADRATURE", EOS%EOS_quadrature, &
                 "If true, always use the generic (quadrature) code "//&
                 "code for the integrals of density.", default=.false.)

  call get_param(param_file, mdl, "TFREEZE_FORM", tmpstr, &
                 "TFREEZE_FORM determines which expression should be "//&
                 "used for the freezing point.  Currently, the valid "//&
                 'choices are "LINEAR", "MILLERO_78", "TEOS10"', &
                 default=TFREEZE_DEFAULT)
  select case (uppercase(tmpstr))
    case (TFREEZE_LINEAR_STRING)
      EOS%form_of_TFreeze = TFREEZE_LINEAR
    case (TFREEZE_MILLERO_STRING)
      EOS%form_of_TFreeze = TFREEZE_MILLERO
    case (TFREEZE_TEOS10_STRING)
      EOS%form_of_TFreeze = TFREEZE_TEOS10
    case default
      call MOM_error(FATAL, "interpret_eos_selection:  TFREEZE_FORM "//&
                              trim(tmpstr) // "in input file is invalid.")
  end select

  if (EOS%form_of_TFreeze == TFREEZE_LINEAR) then
    call get_param(param_file, mdl, "TFREEZE_S0_P0",EOS%TFr_S0_P0, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", "//&
                 "this is the freezing potential temperature at "//&
                 "S=0, P=0.", units="deg C", default=0.0)
    call get_param(param_file, mdl, "DTFREEZE_DS",EOS%dTFr_dS, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", "//&
                 "this is the derivative of the freezing potential "//&
                 "temperature with salinity.", &
                 units="deg C PSU-1", default=-0.054)
    call get_param(param_file, mdl, "DTFREEZE_DP",EOS%dTFr_dP, &
                 "When TFREEZE_FORM="//trim(TFREEZE_LINEAR_STRING)//", "//&
                 "this is the derivative of the freezing potential "//&
                 "temperature with pressure.", &
                 units="deg C Pa-1", default=0.0)
  endif

  if ((EOS%form_of_EOS == EOS_TEOS10 .OR. EOS%form_of_EOS == EOS_NEMO) .AND. &
      EOS%form_of_TFreeze /= TFREEZE_TEOS10) then
      call MOM_error(FATAL, "interpret_eos_selection:  EOS_TEOS10 or EOS_NEMO \n" //&
      "should only be used along with TFREEZE_FORM = TFREEZE_TEOS10 .")
  endif

  ! Unit conversions
  EOS%m_to_Z = 1. ; if (present(US)) EOS%m_to_Z = US%m_to_Z
  EOS%kg_m3_to_R = 1. ; if (present(US)) EOS%kg_m3_to_R = US%kg_m3_to_R
  EOS%R_to_kg_m3 = 1. ; if (present(US)) EOS%R_to_kg_m3 = US%R_to_kg_m3
  EOS%RL2_T2_to_Pa = 1. ; if (present(US)) EOS%RL2_T2_to_Pa = US%RL2_T2_to_Pa
  EOS%L_T_to_m_s = 1. ; if (present(US)) EOS%L_T_to_m_s = US%L_T_to_m_s

end subroutine EOS_init

subroutine extract_member_EOS(EOS, form_of_EOS, form_of_TFreeze, EOS_quadrature, Compressible, &
                              Rho_T0_S0, drho_dT, dRho_dS, TFr_S0_P0, dTFr_dS, dTFr_dp)
  type(EOS_type)     :: EOS !< Equation of state structure
  integer, optional, intent(out) :: form_of_EOS !< A coded integer indicating the equation of state to use.
  integer, optional, intent(out) :: form_of_TFreeze !< A coded integer indicating the expression for
                                       !! the potential temperature of the freezing point.
  logical, optional, intent(out) :: EOS_quadrature !< If true, always use the generic (quadrature)
                                       !! code for the integrals of density.
  logical, optional, intent(out) :: Compressible !< If true, in situ density is a function of pressure.
  real   , optional, intent(out) :: Rho_T0_S0 !< Density at T=0 degC and S=0 ppt [kg m-3]
  real   , optional, intent(out) :: drho_dT   !< Partial derivative of density with temperature
                                              !! in [kg m-3 degC-1]
  real   , optional, intent(out) :: dRho_dS   !< Partial derivative of density with salinity
                                              !! in [kg m-3 ppt-1]
  real   , optional, intent(out) :: TFr_S0_P0 !< The freezing potential temperature at S=0, P=0 [degC]
  real   , optional, intent(out) :: dTFr_dS   !< The derivative of freezing point with salinity
                                              !! [degC PSU-1]
  real   , optional, intent(out) :: dTFr_dp   !< The derivative of freezing point with pressure
                                              !! [degC Pa-1]

  if (present(form_of_EOS    ))  form_of_EOS     = EOS%form_of_EOS
  if (present(form_of_TFreeze))  form_of_TFreeze = EOS%form_of_TFreeze
  if (present(EOS_quadrature ))  EOS_quadrature  = EOS%EOS_quadrature
  if (present(Compressible   ))  Compressible    = EOS%Compressible
  if (present(Rho_T0_S0      ))  Rho_T0_S0       = EOS%Rho_T0_S0
  if (present(drho_dT        ))  drho_dT         = EOS%drho_dT
  if (present(dRho_dS        ))  dRho_dS         = EOS%dRho_dS
  if (present(TFr_S0_P0      ))  TFr_S0_P0       = EOS%TFr_S0_P0
  if (present(dTFr_dS        ))  dTFr_dS         = EOS%dTFr_dS
  if (present(dTFr_dp        ))  dTFr_dp         = EOS%dTFr_dp

end subroutine extract_member_EOS

subroutine calculate_density_scalar(T, S, pressure, rho, EOS, rho_ref, scale)
  real,           intent(in)  :: T        !< Potential temperature referenced to the surface [degC]
  real,           intent(in)  :: S        !< Salinity [ppt]
  real,           intent(in)  :: pressure !< Pressure [Pa] or [R L2 T-2 ~> Pa]
  real,           intent(out) :: rho      !< Density (in-situ if pressure is local) [kg m-3] or [R ~> kg m-3]
  type(EOS_type)              :: EOS      !< Equation of state structure
  real, optional, intent(in)  :: rho_ref  !< A reference density [kg m-3]
  real, optional, intent(in)  :: scale    !< A multiplicative factor by which to scale density in
                                          !! combination with scaling given by US [various]

  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: p_scale   ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]

!  if (.not.associated(EOS)) call MOM_error(FATAL, &
!    "calculate_density_scalar called with an unassociated EOS_type EOS.")

  p_scale = EOS%RL2_T2_to_Pa

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_linear(T, S, p_scale*pressure, rho, &
                                    EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS, rho_ref)
!    case (EOS_UNESCO)
!      call MPP_ERROR(FATAL,'UNESCO not implememted here')
    case (EOS_WRIGHT)
       call calculate_density_wright(T, S, p_scale*pressure, rho, rho_ref)
    case (EOS_WRIGHT_FULL)
      call calculate_density_wright_full(T, S, p_scale*pressure, rho)       
!    case (EOS_TEOS10)
!      call MPP_ERROR(FATAL,'TEOS10 not implememted here')
!      call calculate_density_teos10(T, S, p_scale*pressure, rho, rho_ref)
!    case (EOS_NEMO)
!      call MPP_ERROR(FATAL,'NEMO EOS not implememted here')
!      call calculate_density_nemo(T, S, p_scale*pressure, rho, rho_ref)
    case default
      call MOM_error(FATAL, "calculate_density_scalar: EOS is not valid.")
  end select

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  rho = rho_scale * rho

end subroutine calculate_density_scalar


!> Calls the appropriate subroutine to calculate the density of sea water for 1-D array inputs,
!! potentially limiting the domain of indices that are worked on.
!! If rho_ref is present, the anomaly with respect to rho_ref is returned.
subroutine calculate_density_1d(T, S, pressure, rho, EOS, dom, rho_ref, scale)
  real, dimension(:),    intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:),    intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:),    intent(in)    :: pressure !< Pressure [R L2 T-2 ~> Pa]
  real, dimension(:),    intent(inout) :: rho      !< Density (in-situ if pressure is local) [R ~> kg m-3]
  type(EOS_type)                      :: EOS      !< Equation of state structure
  integer, dimension(2), optional, intent(in) :: dom   !< The domain of indices to work on, taking
                                                       !! into account that arrays start at 1.
  real,                  optional, intent(in) :: rho_ref !< A reference density [kg m-3]
  real,                  optional, intent(in) :: scale !< A multiplicative factor by which to scale density
                                                   !! in combination with scaling given by US [various]
  ! Local variables
  real :: p_scale   ! A factor to convert pressure to units of Pa [Pa T2 R-1 L-2 ~> 1]
  real :: rho_scale ! A factor to convert density from kg m-3 to the desired units [R m3 kg-1 ~> 1]
  real :: rho_unscale ! A factor to convert density from R to kg m-3 [kg m-3 R-1 ~> 1]
  real :: rho_reference ! rho_ref converted to [kg m-3]
  real, dimension(size(rho)) :: pres  ! Pressure converted to [Pa]
  integer :: i, is, ie, npts

!  if (.not.associated(EOS)) call MOM_error(FATAL, &
!    "calculate_density_1d called with an unassociated EOS_type EOS.")

  if (present(dom)) then
    is = dom(1) ; ie = dom(2) ; npts = 1 + ie - is
  else
    is = 1 ; ie = size(rho) ; npts = 1 + ie - is
  endif

  p_scale = EOS%RL2_T2_to_Pa
  rho_unscale = EOS%R_to_kg_m3

  if ((p_scale == 1.0) .and. (rho_unscale == 1.0)) then
    call calculate_density_array(T, S, pressure, rho, is, npts, EOS, rho_ref=rho_ref)
  elseif (present(rho_ref)) then ! This is the same as above, but with some extra work to rescale variables.
    do i=is,ie ; pres(i) = p_scale * pressure(i) ; enddo
    rho_reference = rho_unscale*rho_ref
    call calculate_density_array(T, S, pres, rho, is, npts, EOS, rho_ref=rho_reference)
  else  ! There is rescaling of variables, but rho_ref is not present. Passing a 0 value of rho_ref
        ! changes answers at roundoff for some equations of state, like Wright and UNESCO.
    do i=is,ie ; pres(i) = p_scale * pressure(i) ; enddo
    call calculate_density_array(T, S, pres, rho, is, npts, EOS)
  endif

  rho_scale = EOS%kg_m3_to_R
  if (present(scale)) rho_scale = rho_scale * scale
  if (rho_scale /= 1.0) then ; do i=is,ie
    rho(i) = rho_scale * rho(i)
  enddo ; endif

end subroutine calculate_density_1d

!> Calls the appropriate subroutine to calculate the density of sea water for 1-D array inputs.
!! If rho_ref is present, the anomaly with respect to rho_ref is returned.
subroutine calculate_density_array(T, S, pressure, rho, start, npts, EOS, rho_ref, scale)
  real, dimension(:), intent(in)    :: T        !< Potential temperature referenced to the surface [degC]
  real, dimension(:), intent(in)    :: S        !< Salinity [ppt]
  real, dimension(:), intent(in)    :: pressure !< Pressure [Pa] or [R L2 T-2 ~> Pa]
  real, dimension(:), intent(inout) :: rho      !< Density (in-situ if pressure is local) [kg m-3] or [R ~> kg m-3]
  integer,            intent(in)    :: start    !< Start index for computation
  integer,            intent(in)    :: npts     !< Number of point to compute
  type(EOS_type)       :: EOS      !< Equation of state structure
  real,                  optional, intent(in) :: rho_ref  !< A reference density [kg m-3]
  real,                  optional, intent(in) :: scale    !< A multiplicative factor by which to scale density
                                                !! in combination with scaling given by US [various]
  integer :: j

!  if (.not.associated(EOS)) call MOM_error(FATAL, &
!    "calculate_density_array called with an unassociated EOS_type EOS.")

  select case (EOS%form_of_EOS)
    case (EOS_LINEAR)
      call calculate_density_linear(T, S, pressure, rho, start, npts, &
                                    EOS%Rho_T0_S0, EOS%dRho_dT, EOS%dRho_dS, rho_ref)
!    case (EOS_UNESCO)
!      call calculate_density_unesco(T, S, pressure, rho, start, npts, rho_ref)
    case (EOS_WRIGHT)
       call calculate_density_wright(T, S, pressure, rho, start, npts, rho_ref)
    case (EOS_WRIGHT_FULL)
      call calculate_density_wright_full(T, S, pressure, rho, start, npts)       
!    case (EOS_TEOS10)
!      call calculate_density_teos10(T, S, pressure, rho, start, npts, rho_ref)
!    case (EOS_NEMO)
!    call calculate_density_nemo(T, S, pressure, rho, start, npts, rho_ref)
    case default
      call MOM_error(FATAL, "calculate_density_array: EOS%form_of_EOS is not valid.")
  end select

  if (present(scale)) then ; if (scale /= 1.0) then ; do j=start,start+npts-1
    rho(j) = scale * rho(j)
  enddo ; endif ; endif

end subroutine calculate_density_array

!> Allocates EOS_type
subroutine EOS_allocate(EOS)
  type(EOS_type) :: EOS !< Equation of state structure

!  if (.not.associated(EOS)) allocate(EOS)
end subroutine EOS_allocate

!> Deallocates EOS_type
subroutine EOS_end(EOS)
  type(EOS_type) :: EOS !< Equation of state structure

!  if (associated(EOS)) deallocate(EOS)
end subroutine EOS_end

end module MOM_EOS

!> \namespace mom_eos
!!
!! The MOM_EOS module is a wrapper for various equations of state (e.g. Linear,
!! Wright, UNESCO) and provides a uniform interface to the rest of the model
!! independent of which equation of state is being used.
