! Help:Begin
! NOTE: This program requires the metadata, last checkpoint and trajectory. Outputs datafile heatmap.xy containing
! x y coordinates along with some data z for each coordinate. One can visualize z as heatmap using `ccd visual -z:`
! Usage: ccd_heatmap --record=<recnum> <metadata file path>
! --record : Pass record number to work with. Assume last proper record if not provided.
! Help:End

program ccd_heatmap
    use utilities, only: cmd_line_opt, cmd_line_arg, help_handler
    use parameters, only: m, n, assign_params
    use state_vars
    use files, only: traj_read, cpt_read, open_traj, err_fd
    use gnuplot, only: gp_xy_dump
    use analysis, only: nbeads_per_ring, cell_vicsekop
    implicit none

    character(len=:), allocatable :: metadata_fname, opt_arg
    integer :: metadata_fname_length, opt_arg_len, exitcode
    integer:: rec_index, ring
    double precision, dimension(:, :), allocatable :: z ! Stores the variable for heatmap
    double precision :: vopx, vopy
    integer :: pending_steps, current_step
    character(len=40) :: params_hash

    call help_handler()

    ! Get the metadata file path
    call cmd_line_arg(1, length=metadata_fname_length)
    if (metadata_fname_length == 0) error stop 'Fatal: Pass metadata path as argument'
    allocate (character(len=metadata_fname_length) :: metadata_fname)
    call cmd_line_arg(1, metadata_fname)
    
    ! Assign parameters from metadata file
    call assign_params(fname=metadata_fname, nocheck=.true.)

    call cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)
    
    allocate(z(n,m))
    
    call open_traj('read', 'old')

    ! Get the record number from --record option. Assume last record if --record is not provided.
    call cmd_line_opt('--record', length=opt_arg_len)
    if (opt_arg_len /= 0) then
        allocate (character(len=opt_arg_len) :: opt_arg)
        call cmd_line_opt('--record', opt_arg)
        read (opt_arg, *, iostat=exitcode) rec_index
        if ((rec_index < 1) .or. (rec_index > recnum)) error stop 'Fatal: Provided record range out of bounds'
        deallocate (opt_arg)
    else
        rec_index = recnum
    end if
    write(err_fd,'(a,1x,i0)') 'Using record number:', rec_index

    call traj_read(rec_index, timepoint)
    
    nbeads_per_ring = n
    do ring = 1,m
        call cell_vicsekop(ring, vopx, vopy)
        z(:, ring) = vopx
    end do

    call gp_xy_dump('heatmap.xy', box, x, y, z)
end program ccd_heatmap
