module files
    use iso_fortran_env, only: err_fd => error_unit
    use parameters, only: traj_dump_int, status_dump_int, cpt_dump_int
    use state_vars
    use utilities
    implicit none
    public
    ! File Names
    character(len=*), parameter :: traj_fname = 'traj.bin'
    character(len=*), parameter :: params_fname = '.params.in'
    character(len=*), parameter :: status_fname = 'status.lock'
    character(len=*), parameter :: cpt_fname = 'state.cpt'
    ! File Descriptors
    integer, protected :: traj_fd, status_fd, cpt_fd

    ! File hashes(checksums)
    character(len=40), protected :: init_cpt_hash, final_cpt_hash, traj_hash
    namelist /checksums/ init_cpt_hash, final_cpt_hash, traj_hash

    logical :: do_status_dump = .true. ! This flag may be set and unset using signals

    real, dimension(:, :, :), allocatable, private :: compressed_fp_for_io
    ! Double to single precision compression. Multiple arrays compressed into one big array
    ! Declaring as permanent to avoid overhead of repetitive allocation and deallocation as temporary array
    ! Single precision storage (trajectory only, not checkpoint) and double precision computation following:
    ! https://github.com/FortRun/resources/blob/main/Fortran%20Handbook.md
    ! #real8-or-real16double-precision-why-store-in-real8-and-compute-in-real16

contains

    subroutine open_traj(read_or_write, old_or_replace)
        use ring_nb, only: ring_nb_io
        character(len=*), intent(in) :: read_or_write, old_or_replace
        integer :: reclen, io_stat, alloc_stat

        allocate (compressed_fp_for_io(size(x, 1), size(x, 2), 10), stat=alloc_stat)
        if (alloc_stat /= 0) error stop 'Fatal: Problem in allocating compressed_fp_for_io'
        !TODO: All error stops should contain the subroutine/module name and variable name

        inquire (iolength=reclen) timepoint, compressed_fp_for_io, ring_nb_io
        open (newunit=traj_fd, file=traj_fname, access='direct', recl=reclen, form='unformatted', &
              status=old_or_replace, asynchronous='yes', action=read_or_write, iostat=io_stat)
        if (io_stat /= 0) error stop 'Fatal: Problem with opening '//traj_fname
    end subroutine open_traj

    ! The optional boolean arg `cmframe` asks traj_read to pass coordinates in com frame. Default: .false.
    subroutine traj_read(recnum, timepoint, cmframe)
        use ring_nb, only: ring_nb_io, unpack_ring_nb
        integer, intent(in) :: recnum
        logical, intent(in), optional :: cmframe
        real, intent(out) :: timepoint
        integer :: io_stat, ring, nbeads_per_cell, ncells
        double precision, dimension(size(mx, 1), size(mx, 2)) :: m_norm
        double precision :: gcmx, gcmy ! Global centre of mass

        read (traj_fd, asynchronous='no', rec=recnum, iostat=io_stat) &
            timepoint, compressed_fp_for_io, ring_nb_io
        if (io_stat /= 0) error stop &
            'Fatal: Problem with reading from '//traj_fname//' @ record= '//int_to_char(recnum)
        call unpack_ring_nb()

        x = dble(compressed_fp_for_io(:, :, 1))
        y = dble(compressed_fp_for_io(:, :, 2))
        mx = dble(compressed_fp_for_io(:, :, 3))
        my = dble(compressed_fp_for_io(:, :, 4))
        fx = dble(compressed_fp_for_io(:, :, 5))
        fy = dble(compressed_fp_for_io(:, :, 6))
        f_rpx = dble(compressed_fp_for_io(:, :, 7))
        f_rpy = dble(compressed_fp_for_io(:, :, 8))
        f_adx = dble(compressed_fp_for_io(:, :, 9))
        f_ady = dble(compressed_fp_for_io(:, :, 10))

        nbeads_per_cell = size(x, 1)
        ncells = size(x, 2)

        do ring = 1, ncells
            cmx(ring) = sum(x(:, ring))/nbeads_per_cell
            cmy(ring) = sum(y(:, ring))/nbeads_per_cell
        end do

        m_norm = hypot(mx, my)
        mx = mx/m_norm
        my = my/m_norm

        if (present(cmframe) .and. cmframe) then
            gcmx = sum(cmx)/ncells
            gcmy = sum(cmy)/ncells
            x = x - gcmx
            y = y - gcmy
            cmx = cmx - gcmx
            cmy = cmy - gcmy
        end if
    end subroutine traj_read

    ! Reads only x,y info from trajectory in threadsafe manner
    ! It is threadsafe because it doesn't access any global variable
    subroutine threadsafe_traj_read_xy_only(recnum, timepoint, x, y)
        integer, intent(in) :: recnum
        real, intent(out) :: timepoint
        double precision, dimension(:, :), intent(out) :: x, y
        real, dimension(size(x, 1), size(x, 2), 10) :: compressed_fp_for_io
        integer :: io_stat

        read (traj_fd, asynchronous='no', rec=recnum, iostat=io_stat) timepoint, compressed_fp_for_io
        if (io_stat /= 0) error stop &
            'Fatal: Problem with reading from '//traj_fname//' @ record= '//int_to_char(recnum)

        x = dble(compressed_fp_for_io(:, :, 1))
        y = dble(compressed_fp_for_io(:, :, 2))
    end subroutine threadsafe_traj_read_xy_only

    subroutine traj_write(recnum, timepoint)
        use ring_nb, only: pack_ring_nb, ring_nb_io
        integer, intent(in) :: recnum
        real, intent(in) :: timepoint
        integer :: io_stat

        compressed_fp_for_io(:, :, 1) = real(x)
        compressed_fp_for_io(:, :, 2) = real(y)
        compressed_fp_for_io(:, :, 3) = real(mx)
        compressed_fp_for_io(:, :, 4) = real(my)
        compressed_fp_for_io(:, :, 5) = real(fx)
        compressed_fp_for_io(:, :, 6) = real(fy)
        compressed_fp_for_io(:, :, 7) = real(f_rpx)
        compressed_fp_for_io(:, :, 8) = real(f_rpy)
        compressed_fp_for_io(:, :, 9) = real(f_adx)
        compressed_fp_for_io(:, :, 10) = real(f_ady)

        call pack_ring_nb()
        write (traj_fd, asynchronous='no', rec=recnum, iostat=io_stat) timepoint, compressed_fp_for_io, ring_nb_io
        !TODO: asynchronous='yes' above with introduction of wait(traj_fd) and asynchronous attributes for x,y
        if (io_stat /= 0) error stop 'Fatal: Problem with writing to '//traj_fname//' @ record= '//int_to_char(recnum)
    end subroutine traj_write

    subroutine close_traj()
        close (traj_fd, status='keep')
        traj_hash = sha1(traj_fname)
        deallocate (compressed_fp_for_io)
    end subroutine close_traj

    subroutine cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)
        real, intent(out) :: timepoint
        integer, intent(out) :: recnum, pending_steps, current_step
        character(len=40), intent(out) :: params_hash
        integer :: ncells, nbeads_per_cell, io_stat

        open (newunit=cpt_fd, file=cpt_fname, access='sequential', form='unformatted', &
              status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) error stop 'Fatal: Problem with opening '//cpt_fname
        read (cpt_fd) traj_dump_int, status_dump_int, cpt_dump_int
        read (cpt_fd) ncells, nbeads_per_cell, box
        if (.not. (allocated(prng_seeds))) then
            if (allocate_array_stat(ncells, nbeads_per_cell) /= 0) &
                error stop 'Fatal: Array allocation in heap encountered some problem.'
        end if
        read (cpt_fd) prng_seeds ! Saves current state of the PRNG. To be `put=` in `random_seeds` call...
        read (cpt_fd) x, y, mx, my ! Saves current state of the physical system
        read (cpt_fd) timepoint ! Saves current time instant for the timeseries
        read (cpt_fd) recnum ! Last record number of trajectory file, as of now
        read (cpt_fd) current_step, pending_steps ! How many steps are still pending for the current run
        read (cpt_fd) params_hash
        close (cpt_fd)

        call random_seed(put=prng_seeds)

        init_cpt_hash = sha1(cpt_fname)
    end subroutine cpt_read

    subroutine cpt_write(timepoint, recnum, pending_steps, current_step)
        use parameters, only: m, n
        real, intent(in) :: timepoint
        integer, intent(in) :: recnum, pending_steps, current_step

        ! Complete trajectory-file dumps so far and flush output buffer
        wait(traj_fd)
        flush (traj_fd)

        ! Dumping to .cpt.tmp instead of *.cpt for now
        open (newunit=cpt_fd, file='.cpt.tmp', access='sequential', form='unformatted', &
              status='replace', action='write')
        write (cpt_fd) traj_dump_int, status_dump_int, cpt_dump_int
        write (cpt_fd) m, n, box, box !TODO: Future proof against box_x,box_y for rectangular sim box
        write (cpt_fd) prng_seeds ! Saves current state of the PRNG. To be `put=` in `random_seeds` call...
        write (cpt_fd) x, y, mx, my ! Saves current state of the physical system
        write (cpt_fd) timepoint ! Saves current time instant for the timeseries
        write (cpt_fd) recnum ! Last record number of trajectory file, as of now
        write (cpt_fd) current_step, pending_steps ! How many steps are still pending for the current run
        write (cpt_fd) sha1(params_fname)
        close (cpt_fd)
        ! Atomically moving .cpt.tmp to *.cpt now
        call execute_command_line('mv '//'.cpt.tmp '//cpt_fname)

        final_cpt_hash = sha1(cpt_fname)
    end subroutine cpt_write

    ! Dumps xy file for any frame/timestep to be consumed by third party apps like gnuplot
    ! This routine is threadsafe provided different threads use different `fname`s
    ! x and y are passed as arguments to aid threadsafety
    subroutine xy_dump(fname, boxlen, x, y, title)
        character(len=*), intent(in) :: fname
        double precision, intent(in) :: boxlen
        double precision, dimension(:, :), intent(in) :: x, y
        character(len=*), intent(in), optional :: title
        integer :: fd, l, i

        open (newunit=fd, file=fname, access='sequential', form='formatted', status='replace', &
              action='write')
        if (present(title)) write (fd, '(a,1x,a)') '#Title:', title
        write (fd, '(a,1x,es23.16)') '#Box:', boxlen
        write (fd, '(a)') '#Column headers:'
        write (fd, '(a,4x,a)') 'x', 'y'

        do l = 1, size(x, 2)
            write (fd, '(/)') ! Two consecutive blank records for separating datasets, each containing single cell info
            write (fd, '(a,1x,i0)') '#Cell:', l
            do i = 1, size(x, 1)
                ! Output folded coordinates
                write (fd, '(es23.16,1x,es23.16)') modulo(x(i, l), boxlen), modulo(y(i, l), boxlen)
            end do
            write (fd, '(a,1x,i0)') '#End_Cell:', l
        end do
        close (fd, status='keep')
    end subroutine xy_dump

    ! Status file is a simple text file that serves to show live progress with pv as well as acts as a basic lockfile
    logical function acquire_lock(force)
        integer :: io_stat
        logical, intent(in), optional :: force
        character(len=len('replace')) :: lock_status

        if (present(force) .and. force) then
            lock_status = 'replace'
        else
            lock_status = 'new'
        end if

        open (newunit=status_fd, file=status_fname, access='sequential', form='formatted', status=lock_status, &
              action='write', iostat=io_stat)
        acquire_lock = io_stat == 0
    end function acquire_lock

    subroutine release_lock
        close (status_fd, status='delete')
    end subroutine release_lock

    subroutine status_dump()
        integer :: pending_dumps = 0 ! Holds number of dumps pending because of do_status_dump=.false.

        if (do_status_dump) then
            if (pending_dumps /= 0) then
                write (status_fd, '('//int_to_char(pending_dumps)//'/)')
            else
                write (status_fd, *)
            end if
            flush (status_fd)
            ! Flush doesn't guarantee data would be written to disk,
            ! but warrants it would be available to other processes
            pending_dumps = 0
        else
            pending_dumps = pending_dumps + 1
        end if
    end subroutine status_dump

    subroutine metadata_dump()
        use parameters, only: params
        use iso_fortran_env, only: compiler_version, compiler_options

        write (*, '(/,a,/)') 'Metadata:'
        write (*, nml=params)
        write (*, *)
        write (*, nml=checksums)
        write (*, '(/,a,/,a,/)') 'Compiler version: ', compiler_version()
        write (*, '(a,/,a,/)') 'Compiler Options: ', compiler_options()
    end subroutine metadata_dump

    subroutine perf_dump(cpusec, wcsec, steps)
        real, intent(in) :: cpusec, wcsec
        integer, intent(in) :: steps

        write (err_fd, '(/,a)') 'Performance:'
        write (err_fd, '(a,1x,i0)') '# Steps =', steps
        write (err_fd, '(a)') 'CPU = '//dhms(cpusec)
        write (err_fd, '(a)') 'Wallclock = '//dhms(wcsec)
        write (err_fd, '(a,1x,i0)') '# Threads = ', nint(cpusec/wcsec)
        write (err_fd, '(a,1x,i0,1x,a,/)') 'Rate =', nint((steps/wcsec)*3600), 'steps/hour'
    end subroutine perf_dump

    subroutine log_this(msg)
        character(len=*), intent(in) :: msg
        character(len=8) :: curr_date
        character(len=10) :: curr_time
        call date_and_time(date=curr_date, time=curr_time)
        write (err_fd, '(/,a,1x,a,1x,a,1x,a)') curr_date(7:8)//'-'//curr_date(5:6)//'-'//curr_date(1:4), &
            curr_time(1:2)//':'//curr_time(3:4)//':'//curr_time(5:6), ' => ', msg
    end subroutine log_this

end module files
