! NOTE: This program requires the last checkpoint too.
! Usage: ccd_traj_to_xy <dump directory path> ! Creates the dump directory if non-existent

program ccd_traj_to_xy
    use files
    use utilities, only: int_to_char
    implicit none
    integer :: pending_steps, rec_index
    character(len=40) :: params_hash
    integer :: frame
    character(len=*), parameter :: dump_fname_prefix = 'frame_', dump_fname_suffix='.xy'
    character(len=:), allocatable :: dump_dir
    integer :: dump_dir_str_length, exitcode
    
    ! Get (and create, if needed) the dump directory
    call get_command_argument(1, length=dump_dir_str_length, status=exitcode)
    if(exitcode /= 0) error stop 'Pass a directory path as argument'
    allocate(character(len=dump_dir_str_length) :: dump_dir)
    call get_command_argument(1, dump_dir)
    dump_dir=dump_dir//'/'
    call execute_command_line('mkdir -p '//dump_dir, exitstat=exitcode)
    if(exitcode /= 0) error stop 'Failed to create directory '//dump_dir
    
    call cpt_read(timepoint, recnum, pending_steps, params_hash)
    
    call open_traj('read', 'old')

    !TODO: Parallelize
    do rec_index = 1, recnum
        call traj_read(rec_index, timepoint)        
        frame=rec_index*traj_dump_int
        call xy_dump(fname=dump_dir//dump_fname_prefix//int_to_char(frame)//dump_fname_suffix, &
                boxlen=box, title='Frame: '//int_to_char(frame))
        call log_this('Dumped #frame= '//int_to_char(frame))
    end do
    
    call close_traj()
    call log_this('Done')
end program ccd_traj_to_xy