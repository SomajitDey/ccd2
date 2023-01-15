!Brief: Initializes state. Positions are initialized randomly. Motility vectors initialized isotropically by default.
!Usage: ccd_init [-r | --random]
! -r | --random : Initialize motility unit vectors randomly. Vicsek order parameter may not be zero.

program ccd_init
    use parameters
    use state_vars
    use init
    use files
    use utilities, only: cmd_line_flag
    implicit none
    
    integer :: clock_tick
    real, dimension(:), allocatable :: rands
    
    call assign_params(params_fname)

    if(allocate_array_stat(m, n) /= 0) error stop 'Array allocation in heap encountered some problem.'
    allocate(rands(size(prng_seeds)))

    call system_clock(count=clock_tick)
    call random_number(rands)
    rands = rands + rands - 1.0  ! random numbers in [-1, 1] 
    prng_seeds = nint(rands*clock_tick)
    call random_seed(put = prng_seeds)

    call initial()
    if (cmd_line_flag('-r') .or. cmd_line_flag('--random')) call initial_angle()
    
    call cpt_write(timepoint, recnum, 0, 1)
    
    write(*,'(a,1x,"(",es23.16,",",es23.16,")")') 'Average motility vector:', sum(mx)/(m*n), sum(my)/(m*n)
end program ccd_init
