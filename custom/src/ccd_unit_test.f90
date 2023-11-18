! Just a sample unit test for now
! Usage: ccd_unit_test --opt=testarg
program unit_test
    use utilities, only: cmd_line_opt, mktemp
    use files, only: cpt_read, open_traj, traj_read
    use state_vars
    use ring_nb, only: ring_nb_io
    use voronoi
    implicit none

    character(len=:), allocatable :: argument
    integer :: arglen
    integer :: record = 100
    integer :: pending_steps, current_step
    character(len=40) :: params_hash

    call cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)
    call open_traj('read', 'old')
    call traj_read(record, timepoint)
    call periodic_voronoi(cmx, cmy, box, ring_nb_io, 'vor.xy')
end program unit_test
