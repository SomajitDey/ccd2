! Brief: The main run engine. Produces trajectory starting from an initial state provided by a checkpoint.
! Prerequisites: checkpoint, parameter file, trajectory to append to if --append flag is on.
! Synopsis: ccd_run [--append | -a] [--force | -f] [-n | --no-status-dump]
! --append : append to an existing trajectory thus extending a previous run
! --force : ignore lockfile left behind by a previous incomplete run
! --no-status-dump : Doesn't show live progess. This makes the run faster. 

program ccd_run
	use shared
    use prerun
    use grid_linked_list
    use forces
    use integrator
    use ring_nb, only: init_ring_nb

    implicit none

	integer:: j1,jf
    real:: cpusec,wcsec

    ! Initialize/Pre-run setup
    call prerun_setup(jf)
    
    call log_this('Setting up neighbor list grids')
	call gridmaps

    call log_this('Starting the main run')

    call timestamp()

	!$omp parallel default(shared) private(j1)
    timeseries: do j1=1,jf

	call force()

    !TODO: Turning the 2 omp singles below into omp sections gives either seg fault or divide by 0 error. Why? 
    
    !$omp master
        ! Master makes sure random_number is called from a particular thread alone, i.e. the master thread.
        ! This is important as gfortran provides different prng sequences (different seeds) for different threads.
        call random_seed(get = prng_seeds)
             CALL gasdev(noise,mean,var) ! Updates prng_seeds as side-effect        
    !$omp end master

    !$omp single
    ! links() couldn't be parallelized as most of it needs to run sequentially. 
    ! omp ordered would just increase overhead.
    call links()
    !$omp end single nowait
    
    !$omp single
        if(mod(j1,traj_dump_int).eq.0) call init_ring_nb()
    !$omp end single

	call interaction(store_ring_nb = mod(j1,traj_dump_int).eq.0)

    !$omp sections
    !$omp section
    traj_dump: if(mod(j1,traj_dump_int).eq.0) then
            call traj_write(recnum, timepoint)

            cpt_dump: if(mod(j1,cpt_dump_int).eq.0) then
                call cpt_write(timepoint, recnum, jf-j1)
                call log_this('Created checkpoint @ timestep = '//int_to_char(j1))
             end if cpt_dump

            recnum = recnum + 1 ! Update record number must be after cpt_dump, if any
        end if traj_dump

    !$omp section
        for_pv: if(mod(j1,status_dump_int).eq.0) then
            call status_dump()
        end if for_pv
    !$omp end sections
    
	call move_noise() ! Update state

    !$omp single
	timepoint = timepoint + dt ! Update timepoint
    !$omp end single

    end do timeseries
    !$omp end parallel

    call timestamp(cpusec,wcsec)

    call log_this('Run complete. Writing final checkpoint')
    call cpt_write(timepoint-dt, recnum-1, 0)

    call close_traj()

    call metadata_dump()

    call perf_dump(cpusec, wcsec, jf)

    call log_this('Done')

    call release_lock
end program ccd_run
