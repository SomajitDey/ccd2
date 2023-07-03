! Just a sample unit test for now
! Usage: ccd_unit_test --opt=testarg dummyargs
program unit_test
    use utilities, only: cmd_line_opt
    implicit none

    character(len=:), allocatable :: argument
    integer :: arglen

    call cmd_line_opt('--opt', length=arglen)
    allocate (character(len=arglen) :: argument)
    call cmd_line_opt('--opt', argument)
    print"(a,1x,'length=',i0)", argument, arglen
end program unit_test
