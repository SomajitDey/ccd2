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

	timeseries: do j1=1,jf

    call links
	call force
	call interaction()

        traj_dump: if(mod(j1,traj_dump_int).eq.0) then
            call traj_write(recnum, timepoint)

            cpt_dump: if(mod(j1,cpt_dump_int).eq.0) then
                call cpt_write(timepoint, recnum, jf-j1)
                call log_this('Created checkpoint @ timestep = '//int_to_char(j1))
             end if cpt_dump

            recnum = recnum + 1 ! Update record number must be after cpt_dump, if any
        end if traj_dump

        for_pv: if(mod(j1,status_dump_int).eq.0) then
            call status_dump()
        end if for_pv	

	call move_noise ! Update state

	timepoint = timepoint + dt ! Update timepoint
    
    end do timeseries

    call timestamp(cpusec,wcsec)

    call log_this('Run complete. Writing final checkpoint')
    call cpt_write(timepoint-dt, recnum-1, 0)

    call close_traj()

    call metadata_dump()

    call perf_dump(cpusec, wcsec, jf)

    call log_this('Done')

    call release_lock
end program ccd_run
