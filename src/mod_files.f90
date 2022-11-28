module files
    use iso_fortran_env, only: err_fd => error_unit
    implicit none
    public
    ! File Names
    character(len=*), parameter :: traj_fname='traj.bin'
    character(len=*), parameter :: final_fname='final.xy'
    character(len=*), parameter :: params_fname='params.in'
    character(len=*), parameter :: status_fname='status.live'
    character(len=*), parameter :: cpt_fname='state.cpt'
    ! File Descriptors
    integer :: traj_fd, final_fd, status_fd, cpt_fd
    character(len=40) :: init_cpt_hash, final_cpt_hash, traj_hash
    
    namelist /checksums/ init_cpt_hash, final_cpt_hash, traj_hash
    
    contains
    
    subroutine close_files()
        close(traj_fd, status='keep')
        close(final_fd, status='keep')
        close(status_fd, status='delete')
    end subroutine close_files
    
    subroutine log_this(msg)
        character(len=*), intent(in) :: msg
        character(len=8) :: curr_date
        character(len=10) :: curr_time
        call date_and_time(date=curr_date, time=curr_time)
        write(err_fd,'(/,a,1x,a,1x,a,1x,a)') curr_date(7:8) // '-' // curr_date(5:6) // '-' // curr_date(1:4), &
            curr_time(1:2) // ':' // curr_time(3:4) // ':' // curr_time(5:6), ' => ', msg
    end subroutine log_this
    
end module files
