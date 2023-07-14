! Help:Begin
! NOTE: This program requires the last checkpoint too.
! Usage: ccd_traj_to_xy [--records=<begin>:<end>] <dump directory path> ! Creates the dump directory if non-existent
! --records : Pass range of records to work with. Omit either <begin> or <end> to assume default. E.g. --records=3:4
! Help:End

program ccd_traj_to_xy
    use files
    use utilities, only: cmd_line_opt, cmd_line_arg, int_to_char, help_handler
!$  use omp_lib, only: omp_get_max_threads
    implicit none
    integer :: pending_steps, current_step, rec_index, begin_rec, end_rec
    character(len=40) :: params_hash
    integer :: frame
    character(len=*), parameter :: dump_fname_prefix = 'frame_', dump_fname_suffix = '.xy'
    character(len=:), allocatable :: dump_dir, opt_arg
    integer :: dump_dir_str_length, exitcode, opt_arg_len

    ! Following variables with trailing _ would be threadprivate
    double precision, dimension(:, :), allocatable :: x_, y_
    real :: timepoint_

    call help_handler()

    ! Get (and create, if needed) the dump directory
    call cmd_line_arg(1, length=dump_dir_str_length)
    if (dump_dir_str_length == 0) error stop 'Fatal: Pass a directory path as argument'
    allocate (character(len=dump_dir_str_length) :: dump_dir)
    call cmd_line_arg(1, dump_dir)
    dump_dir = dump_dir//'/'
    call execute_command_line('mkdir -p '//dump_dir, exitstat=exitcode)
    if (exitcode /= 0) error stop 'Fatal: Failed to create directory '//dump_dir

    call cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)
    allocate (x_(size(x, 1), size(x, 2)), y_(size(y, 1), size(y, 2)))

    call open_traj('read', 'old')

    ! Sort out the begin and end record number from --records=<from>:<to> cmd line option, if any
    begin_rec = 1 ! default
    end_rec = recnum ! default
    call cmd_line_opt('--records', length=opt_arg_len)
    if (opt_arg_len /= 0) then
        allocate (character(len=opt_arg_len) :: opt_arg)
        call cmd_line_opt('--records', opt_arg)
        read (opt_arg(:scan(opt_arg, ':') - 1), *, iostat=exitcode) begin_rec
        read (opt_arg(scan(opt_arg, ':') + 1:), *, iostat=exitcode) end_rec
        if (begin_rec > end_rec) error stop 'Fatal: Records range provided must be in ascending order'
        if ((begin_rec < 1) .or. (end_rec > recnum)) error stop 'Fatal: Provided records range out of bounds'
        deallocate (opt_arg)
    end if

!$  call log_this('Using '//int_to_char(int(omp_get_max_threads(), kind=kind(rec_index)))//' OpenMP threads')

!$omp parallel do default(shared) private(rec_index, frame, x_, y_, timepoint_)
    do rec_index = begin_rec, end_rec
        call threadsafe_traj_read_xy_only(rec_index, timepoint_, x_, y_)
        frame = rec_index*traj_dump_int
        call xy_dump(fname=dump_dir//dump_fname_prefix//int_to_char(frame)//dump_fname_suffix, &
                     boxlen=box, x=x_, y=y_, title='Frame: '//int_to_char(frame))
        call log_this('Dumped #frame= '//int_to_char(frame))
    end do
!$omp end parallel do

    call close_traj()
    call log_this('Done')
end program ccd_traj_to_xy
