!TODO: Check if memory can be freed by deallocating
!TODO: Check if shared module is really required. All 'use module' statement must have 'only:'
module shared
    use utilities
    use parameters
    use state_vars
    use files
    implicit none
    public
    double precision, parameter::  pi = dacos(-1.0d0)
end module shared
