!> Provides transparent structures with groups of MOM6 variables and supporting routines
module MOM_variables

! This file is part of MOM6. See LICENSE.md for the license.

  use MOM_eos, only : EOS_type
  
public ocean_grid_type, thermo_var_ptrs


!> Pointers to an assortment of thermodynamic fields that may be available, including
!! potential temperature, salinity, heat capacity, and the equation of state control structure.
!type, public :: thermo_var_ptrs
!   logical :: used
!end type thermo_var_ptrs

type :: thermo_var_ptrs
   real, dimension(:,:,:), pointer :: T=>NULL()
   real, dimension(:,:,:), pointer :: S=>NULL()
   type(EOS_type) :: eqn_of_state
   real :: p_ref
end type thermo_var_ptrs

!> Pointers to an assortment of thermodynamic fields that may be available, including
!! potential temperature, salinity, heat capacity, and the equation of state control structure.
type, public :: ocean_grid_type
   integer :: isc=1, iec=1, jsc=1, jec=1
   integer :: isd=1, ied=1, jsd=1, jed=1
   real, pointer, dimension(:,:) :: mask2dT
   real, pointer, dimension(:,:) :: bathyT
end type ocean_grid_type


end module MOM_variables
