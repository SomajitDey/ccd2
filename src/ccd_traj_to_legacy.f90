! NOTE: This program requires the last checkpoint too.

program ccd_traj_to_legacy
    use files
    use parameters, only: traj_dump_int
    implicit none
    integer :: pending_steps, rec_index, ncell, nbeads_per_cell, l, i, legacy_fd
    character(len=40) :: params_hash
    character(len=*), parameter :: legacy_fname = 'legacy.traj.txt'
    
    call cpt_read(timepoint, recnum, pending_steps, params_hash)
    
    call open_traj('read', 'old')

    open(newunit=legacy_fd, file=legacy_fname, status='replace')

    ncell = size(x,1)
    nbeads_per_cell = size(x,2)

    write(legacy_fd,*) 'Ignore this line'

    
    do rec_index = 1, recnum
        call traj_read(rec_index, timepoint)
        do l=1,ncell
            do i=1,nbeads_per_cell
                write(legacy_fd,*)rec_index*traj_dump_int,l,i,x(l,i),y(l,i)
            end do
        end do
    end do

end program ccd_traj_to_legacy
