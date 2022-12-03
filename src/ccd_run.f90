! Brief: The main run engine. Produces trajectory starting from an initial state provided by a checkpoint.
! Prerequisites: checkpoint, parameter file, trajectory to append to if --append flag is on.
! Synopsis: ccd_run [--append | -a] [--force | -f]
! --append : append to an existing trajectory thus extending a previous run
! --force : ignore lockfile left behind by a previous incomplete run

program ccd_run
	use shared
    use prerun
    use grid_linked_list
    use forces
    use integrator

    implicit none

	integer:: j1,jf
    real:: cpusec,wcsec

    ! Initialize/Pre-run setup
    call prerun_setup(jf)
    
    call log_this('Setting up neighbor list grids')
	call maps
   
    call log_this('Starting the main run')

    call timestamp()

	!$omp parallel default(shared) private(j1)
    timeseries: do j1=1,jf

    !$omp single
    !TODO: Parallelize links()
    call links()
    !$omp end single
    
	call force()

	call interaction()

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
    
    !$omp single
    !TODO: Parallelize move_noise()
	call move_noise() ! Update state
    !$omp end single

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
