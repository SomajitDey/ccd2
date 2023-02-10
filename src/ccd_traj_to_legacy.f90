! Help:Begin
! NOTE: This program requires the last checkpoint too.
! Usage: ccd_traj_to_legacy
! Help:End

program ccd_traj_to_legacy
    use files
    use utilities, only: help_handler
    implicit none
    integer :: pending_steps, current_step, rec_index, ncell, nbeads_per_cell, l, i, legacy_fd
    character(len=40) :: params_hash
    character(len=*), parameter :: legacy_fname = 'legacy.traj.txt'
    
    call help_handler()
    
    call cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)
    
    call open_traj('read', 'old')

    open(newunit=legacy_fd, file=legacy_fname, status='replace')

    ncell = size(x,2)
    nbeads_per_cell = size(x,1)
    
    do rec_index = 1, recnum
        call traj_read(rec_index, timepoint)
        do l=1,ncell
            do i=1,nbeads_per_cell
                write(legacy_fd,'(3(i0,1x), es23.16, 1x, es23.16)') rec_index*traj_dump_int,l,i,x(i,l),y(i,l)
            end do
        end do
    end do

end program ccd_traj_to_legacy
