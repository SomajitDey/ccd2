program ccd_rinit
    use parameters
    use state_vars
    use init
    use files
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
    call initial_angle()
    
    call cpt_write(timepoint, recnum, 0)
end program ccd_rinit
