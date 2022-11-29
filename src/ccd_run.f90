program ccd_run
	use shared
    use init
    use grid_linked_list
    use forces
    use integrator

    implicit none

	integer:: j1,recnum=0
    real:: cpusec,wcsec
    logical:: another_run_is_live

    call assign_params(params_fname)
    call log_this('Run parameters read in')


    call log_this('Creating all necessary arrays (such as positions) @ heap')
    if(allocate_array_stat(m,n) /= 0) error stop 'Array allocation in heap encountered some problem.'

    ! Check if another run is live and Open status file to manifest live status.
    ! This has to be done before the run changes any state / writes anything to memory.
    inquire(file=status_fname, exist=another_run_is_live)
    if(another_run_is_live) error stop 'Another run is going on'
    open(newunit=status_fd,file=status_fname, access='sequential', form='formatted',status='new', &
        asynchronous='yes', action='write')

    write(status_fd,*,asynchronous='yes') jf/status_dump_int
    flush(status_fd) 
    ! This doesn't guarantee data would be written to disk, but warrants it would be available to other processes

    ! Open trajectory file
    call log_this('Initiating trajectory file: '//traj_fname)
    call open_traj('write', 'replace')
    ! Initialize/Pre-run setup
    call log_this('Pre-run setups...')
    call initial      
	call maps
    call initial_angle
   
    call timestamp()
    
    call log_this('Starting the main run')
	timeseries: do j1=1,jf

    call links
	call force
	call interaction()

        traj_dump: if(mod(j1,traj_dump_int).eq.0) then
             recnum = recnum + 1
            call traj_write(recnum, j1)

            cpt_dump: if(mod(j1,cpt_dump_int).eq.0) then
                call cpt_write(recnum, jf-j1)
                call log_this('Created checkpoint @ timestep = '//int_to_char(j1))
             end if cpt_dump
        end if traj_dump

        for_pv: if(mod(j1,status_dump_int).eq.0) then
            write(status_fd,*,asynchronous='yes')
            flush(status_fd)
        end if for_pv	

	call move_noise

	end do timeseries

    call log_this('Run complete. Writing final checkpoint')
    call traj_write(recnum+1, jf)
    call cpt_write(recnum+1, 0)

    call timestamp(cpusec,wcsec)

    call close_traj()

    call log_this('Writing final config: '//final_fname)
    call xy_dump(final_fname)

        close(status_fd, status='delete')

    call metadata_dump()

    call log_this('Done')

    call perf_dump(cpusec, wcsec, jf)
end program ccd_run
