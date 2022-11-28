program ccd_run
	use shared
    use grid_linked_list
    use forces
    use integrator

    implicit none

	integer:: l,i,j1,reclen,recnum
    real:: cpusec,wcsec
    logical:: another_run_is_live
    character(len=10):: buffer  ! Internal file

    call assign_params(params_fname)
    call log_this('Run parameters read in')

    write(*,'(/,a,/)') 'Metadata:'
    write(*,nml=params) ! Dump all params on STDOUT

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
	inquire(iolength=reclen) j1, x, y, mx, my, fx, fy, f_rpx, f_rpy, f_adx, f_ady
    open(newunit=traj_fd,file=traj_fname, access='direct', recl=reclen, form='unformatted', status='replace', &
        asynchronous='yes', action='write')

    ! Initialize/Pre-run setup
    call log_this('Pre-run setups...')
    call initial      
	call maps
    call initial_angle
   
    call timestamp()
    
    call log_this('Starting the main run')
	timeseries: do j1=0,jf

    call links
	call force
	call interaction()

        traj_dump: if(mod(j1,traj_dump_int).eq.0) then
             recnum = j1/traj_dump_int + 1
             write(traj_fd, asynchronous='yes', rec=recnum) &
                j1, x, y, mx, my, fx, fy, f_rpx, f_rpy, f_adx, f_ady

            cpt_dump: if(mod(j1,cpt_dump_int).eq.0) then
                ! Complete trajectory-file dumps so far and flush output buffer
                wait(traj_fd)
                flush(traj_fd)

                ! Dumping to .cpt.tmp instead of *.cpt for now
                open(newunit=cpt_fd,file='.cpt.tmp', access='sequential', form='unformatted', &
                    status='replace', action='write')
                    call random_seed(get = prng_seeds)
                    write(cpt_fd) prng_seeds ! Saves current state of the PRNG. To be `put=` in `random_seeds` call...
                    write(cpt_fd) x,y,mx,my ! Saves current state of the physical system
                    write(cpt_fd) recnum ! Last record number of trajectory file, as of now
                close(cpt_fd)
                ! Atomically moving .cpt.tmp to *.cpt now
                call execute_command_line('mv '//'.cpt.tmp '//cpt_fname, wait=.false.)

                write(buffer,'(i0)') j1
                call log_this('Created checkpoint @ timestep = '//trim(buffer))
             end if cpt_dump
        end if traj_dump

        status_dump: if(mod(j1,status_dump_int).eq.0) then
            write(status_fd,*,asynchronous='yes')
            flush(status_fd)
        end if status_dump	

	call move_noise

	end do timeseries

    call timestamp(cpusec,wcsec)

    call log_this('Run complete...writing final config: '//final_fname)
        ! Open final config file
        open(newunit=final_fd, file=final_fname, access='sequential', form='formatted',status='replace', &
            action='write')
          final_dump: do l=1,M
            write(final_fd,'(a,1x,i0)') '#Cell:', l
            do i=1,N        
				   x(l,i) = x(l,i) - box*floor(x(l,i)/box)
				   y(l,i) = y(l,i) - box*floor(y(l,i)/box)
                write(final_fd,*) x(l,i),y(l,i)
            end do
            write(final_fd,'(a,1x,i0,/)') '#End_Cell:', l
          end do final_dump

    call close_files()

    final_cpt_hash = sha1(cpt_fname)
    traj_hash = sha1(traj_fname)
    write(*,nml=checksums)

    call log_this('Done')

    write(err_fd,'(/,a)') 'Performance:'
    write(err_fd,'(a)') 'CPU = '//dhms(cpusec)
    write(err_fd,'(a)') 'Wallclock = '//dhms(wcsec)
    write(err_fd,'(a,1x,i0,/)') '# Threads = ', nint(cpusec/wcsec) 

end program ccd_run
