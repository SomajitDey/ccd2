module files
    use iso_fortran_env, only: err_fd => error_unit
    use state_vars
    use utilities
    implicit none
    public
    ! File Names
    character(len=*), parameter :: traj_fname='traj.bin'
    character(len=*), parameter :: params_fname='params.in'
    character(len=*), parameter :: status_fname='status.lock'
    character(len=*), parameter :: cpt_fname='state.cpt'
    ! File Descriptors
    integer :: traj_fd, final_fd, status_fd, cpt_fd
    character(len=40) :: init_cpt_hash, final_cpt_hash, traj_hash
    
    namelist /checksums/ init_cpt_hash, final_cpt_hash, traj_hash

    logical :: do_status_dump = .true. ! This flag may be set and unset using signals
    
    contains
    
    subroutine open_traj(read_or_write, old_or_replace)
        character(len=*), intent(in) :: read_or_write, old_or_replace
        integer :: reclen, io_stat

        inquire(iolength=reclen) timepoint, x, y, mx, my, fx, fy, f_rpx, f_rpy, f_adx, f_ady
        open(newunit=traj_fd,file=traj_fname, access='direct', recl=reclen, form='unformatted', &
            status=old_or_replace, asynchronous='yes', action=read_or_write, iostat=io_stat)
        if(io_stat /= 0) error stop 'Problem with opening '//traj_fname
    end subroutine open_traj
    
    subroutine traj_read(recnum, timepoint)
        integer, intent(in) :: recnum
        double precision, intent(out) :: timepoint
        integer :: io_stat
        read(traj_fd, asynchronous='no', rec=recnum, iostat=io_stat) &
            timepoint, x, y, mx, my, fx, fy, f_rpx, f_rpy, f_adx, f_ady
        if(io_stat /= 0) error stop 'Problem with reading from '//traj_fname//' @ record= '//int_to_char(recnum)
    end subroutine traj_read

    subroutine traj_write(recnum, timepoint)
        integer, intent(in) :: recnum
        double precision, intent(in) :: timepoint
        integer :: io_stat
        write(traj_fd, asynchronous='yes', rec=recnum, iostat=io_stat) &
            timepoint, x, y, mx, my, fx, fy, f_rpx, f_rpy, f_adx, f_ady
        if(io_stat /= 0) error stop 'Problem with writing to '//traj_fname//' @ record= '//int_to_char(recnum)
    end subroutine traj_write

    subroutine close_traj()
        close(traj_fd, status='keep')
        traj_hash = sha1(traj_fname)
    end subroutine close_traj

    subroutine cpt_read(timepoint, recnum, pending_steps, params_hash)
                double precision, intent(out) :: timepoint
                integer, intent(out) :: recnum, pending_steps
                character(len=40), intent(out) :: params_hash
                integer :: ncells, nbeads_per_cell, io_stat

                open(newunit=cpt_fd,file=cpt_fname, access='sequential', form='unformatted', &
                    status='old', action='read', iostat=io_stat)
                if(io_stat /= 0) error stop 'Problem with opening '//cpt_fname
                    read(cpt_fd) ncells, nbeads_per_cell, box
                    if(.not.(allocated(prng_seeds))) then
                        if(allocate_array_stat(ncells, nbeads_per_cell) /= 0) &
                            error stop 'Array allocation in heap encountered some problem.'
                    end if
                    read(cpt_fd) prng_seeds ! Saves current state of the PRNG. To be `put=` in `random_seeds` call...
                    read(cpt_fd) x,y,mx,my ! Saves current state of the physical system
                    read(cpt_fd) timepoint ! Saves current time instant for the timeseries
                    read(cpt_fd) recnum ! Last record number of trajectory file, as of now
                    read(cpt_fd) pending_steps ! How many steps are still pending for the current run
                    read(cpt_fd) params_hash
                close(cpt_fd)

                call random_seed(put = prng_seeds)

                init_cpt_hash = sha1(cpt_fname)
    end subroutine cpt_read
    
    subroutine cpt_write(timepoint, recnum, pending_steps)
                use parameters, only: m,n
                double precision, intent(in) :: timepoint
                integer, intent(in) :: recnum, pending_steps

                ! Complete trajectory-file dumps so far and flush output buffer
                wait(traj_fd)
                flush(traj_fd)

                call random_seed(get = prng_seeds)

                ! Dumping to .cpt.tmp instead of *.cpt for now
                open(newunit=cpt_fd,file='.cpt.tmp', access='sequential', form='unformatted', &
                    status='replace', action='write')
                    write(cpt_fd) m,n,box
                    write(cpt_fd) prng_seeds ! Saves current state of the PRNG. To be `put=` in `random_seeds` call...
                    write(cpt_fd) x,y,mx,my ! Saves current state of the physical system
                    write(cpt_fd) timepoint ! Saves current time instant for the timeseries
                    write(cpt_fd) recnum ! Last record number of trajectory file, as of now
                    write(cpt_fd) pending_steps ! How many steps are still pending for the current run
                    write(cpt_fd) sha1(params_fname)
                close(cpt_fd)
                ! Atomically moving .cpt.tmp to *.cpt now
                call execute_command_line('mv '//'.cpt.tmp '//cpt_fname)

                final_cpt_hash = sha1(cpt_fname)
    end subroutine cpt_write

    ! Dumps xy file for any frame/timestep to be consumed by third party apps like gnuplot
    subroutine xy_dump(fname, boxlen)
        double precision, intent(in) :: boxlen
        character(len=*), intent(in) :: fname
        integer :: fd, l, i

        open(newunit=fd, file=fname, access='sequential', form='formatted',status='replace', &
            action='write')
            write(fd,'(a,1x,es23.16)') '#Box:', boxlen
          do l=1,size(x,1)
            write(fd,'(a,1x,i0)') '#Cell:', l
            do i=1,size(x,2)     
				   x(l,i) = x(l,i) - boxlen*floor(x(l,i)/boxlen)
				   y(l,i) = y(l,i) - boxlen*floor(y(l,i)/boxlen)
                write(fd,*) x(l,i),y(l,i)
            end do
            write(fd,'(a,1x,i0,/)') '#End_Cell:', l
          end do
        close(fd, status='keep')
    end subroutine xy_dump
    
    ! Status file is a simple text file that serves to show live progress with pv as well as acts as a basic lockfile
    logical function acquire_lock(force)
        integer :: io_stat
        logical, intent(in), optional :: force
        character(len=len('replace')) :: lock_status

        if(present(force).and.force)then
            lock_status='replace'
        else
            lock_status='new'
        end if

        open(newunit=status_fd,file=status_fname, access='sequential', form='formatted',status=lock_status, &
            action='write', iostat=io_stat)
        acquire_lock = io_stat == 0
    end function acquire_lock
    
    subroutine release_lock
        close(status_fd, status='delete')
    end subroutine release_lock
    
    subroutine status_dump()
        integer :: pending_dumps = 0 ! Holds number of dumps pending because of do_status_dump=.false.

        if(do_status_dump)then
            if(pending_dumps /= 0) then
                write(status_fd, '('//int_to_char(pending_dumps)//'/)')
            else
                write(status_fd,*)
            end if
            flush(status_fd)
            ! Flush doesn't guarantee data would be written to disk, but warrants it would be available to other processes
            pending_dumps=0
        else
            pending_dumps = pending_dumps+1
        end if
    end subroutine status_dump
    
    subroutine metadata_dump()
        use parameters, only: params
        use iso_fortran_env, only: compiler_version, compiler_options

        write(*,'(/,a,/)') 'Metadata:'
        write(*,nml=params)
        write(*,*)
        write(*,nml=checksums)
        write(*,'(/,a,/,a,/)') 'Compiler version: ', compiler_version()
        write(*,'(a,/,a,/)') 'Compiler Options: ', compiler_options()
    end subroutine metadata_dump
    
    subroutine perf_dump(cpusec, wcsec, steps)
        real, intent(in) :: cpusec, wcsec
        integer, intent(in) :: steps

        write(err_fd,'(/,a)') 'Performance:'
        write(err_fd,'(a,1x,i0)') '# Steps =', steps
        write(err_fd,'(a)') 'CPU = '//dhms(cpusec)
        write(err_fd,'(a)') 'Wallclock = '//dhms(wcsec)
        write(err_fd,'(a,1x,i0)') '# Threads = ', nint(cpusec/wcsec)
        write(err_fd,'(a,1x,i0,1x,a,/)')'Rate =', nint(steps*3600/wcsec), 'steps/hour'
    end subroutine perf_dump
   
    subroutine log_this(msg)
        character(len=*), intent(in) :: msg
        character(len=8) :: curr_date
        character(len=10) :: curr_time
        call date_and_time(date=curr_date, time=curr_time)
        write(err_fd,'(/,a,1x,a,1x,a,1x,a)') curr_date(7:8) // '-' // curr_date(5:6) // '-' // curr_date(1:4), &
            curr_time(1:2) // ':' // curr_time(3:4) // ':' // curr_time(5:6), ' => ', msg
    end subroutine log_this
    
end module files
