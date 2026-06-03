!> Routines for error handling and I/O management
module MOM_error_handler

implicit none ; private

! These routines are found in this module.
public :: MOM_error, MOM_mesg, is_root_pe
!> Integer parameters encoding the severity of an error message
public :: NOTE, WARNING, FATAL
public :: stdlog, stdout

integer :: NOTE = 0
integer :: WARNING = 1
integer :: FATAL = 2

contains


integer function stdout()

  stdout = 6

end function stdout


integer function stdlog()

  stdlog = -1

end function stdlog

logical function is_root_pe()

  is_root_pe = .true.

end function is_root_pe

!> This provides a convenient interface for writing an informative comment, depending
!! on the model's current verbosity setting and the verbosity level for this message.
subroutine MOM_mesg(message, verb, all_print)
  character(len=*), intent(in)  :: message !< A message to write out
  integer, optional, intent(in) :: verb !< A level of verbosity for this message
  logical, optional, intent(in) :: all_print !< If present and true, any PEs are
                                             !! able to write this message.
  ! This provides a convenient interface for writing an informative comment.
  integer :: verb_msg
  logical :: write_msg

  print *,message

end subroutine MOM_mesg

!> This provides a convenient interface for writing an error message
!! with run-time filter based on a verbosity and the severity of the error.
subroutine MOM_error(level, message, all_print)
  integer,           intent(in) :: level !< The severity level of this message
  character(len=*),  intent(in) :: message !< A message to write out
  logical, optional, intent(in) :: all_print !< If present and true, any PEs are
                                             !! able to write this message.
  ! This provides a convenient interface for writing an error message
  ! with run-time filter based on a verbosity.
  logical :: write_msg

  print *,message
end subroutine MOM_error

end module MOM_error_handler
