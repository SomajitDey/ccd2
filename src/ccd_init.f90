! Help:Begin
!Brief: Initializes state. Positions are initialized randomly. Motility vectors initialized isotropically by default.
!Usage: ccd_init [--no-check] [--box=<value> | --box=<scale>%] [--pl0]
! -n | --no-check : Do not check parameters for consistency
! --box : Provide desired box length. E.g. --box=5.5. Can also provide relative length in percentage. E.g. --box=200%
! --pl0 : Use p*l0 for pressure force when computing equilibrium radius from spring-pressure force balance
! Help:End

program ccd_init
    use parameters
    use state_vars
    use init
    use files
    use utilities, only: cmd_line_flag, help_handler
    implicit none

    integer :: clock_tick
    real, dimension(:), allocatable :: rands

    call help_handler()

    call assign_params(params_fname, nocheck=cmd_line_flag('-n') .or. cmd_line_flag('--no-check'))

    if (allocate_array_stat(m, n) /= 0) error stop 'Fatal: Array allocation in heap encountered some problem.'
    allocate (rands(size(prng_seeds)))

    call system_clock(count=clock_tick)
    call random_number(rands)
    rands = rands + rands - 1.0  ! random numbers in [-1, 1]
    prng_seeds = nint(rands*clock_tick)
    call random_seed(put=prng_seeds)

    call initial()

    call cpt_write(timepoint, recnum, 0, 1)

    write (err_fd, '(a,1x,"(",es23.16,",",es23.16,")")') 'Average motility vector:', sum(mx)/(m*n), sum(my)/(m*n)
    write (err_fd, '(a,1x,es23.16)') 'Box:', box
end program ccd_init
