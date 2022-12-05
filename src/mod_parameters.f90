module parameters
    implicit none
    public
    double precision, protected:: k=120.0d0      !  Single cell spring constant
    double precision, protected:: p=25.0d0       !  Single cell internal hydrostatic pressure coefficient
    double precision, protected:: l0=0.1d0       !  Single cell natural spring-length
    double precision, protected:: rc_adh=0.28d0  ! Adhesion interaction cut-off
    double precision, protected:: rc_rep=0.18d0  ! Repulsion interaction cut-off
    double precision, protected:: k_adh=0.001d0   !  Adhesion interaction strength
    double precision, protected:: k_rep=1000.0d0 !  Adhesion interaction strength
    double precision, parameter:: mean=0.0d0     !  Mean of the gaussian white noise
    double precision, protected:: var=0.05d0      !  Variance of the gaussian white noise
    double precision, protected:: Vo=0.05d0       !  Self propulsion of the beads
    double precision, protected:: c = 0.5d0         ! c is coeff. of viscous damping      
    double precision, protected:: dt=0.001d0   ! Integration timestep
    integer, protected:: tau_align=10 ! Tau for Vicsek alignment in multiples of dt
    integer, protected:: nsamples=100500  !! No. of Iterations in terms of traj_dump_int
    integer,protected:: n = 50    ! No. of beads
    integer,protected:: m = 256   ! No. of cell

    namelist /params/ dt, nsamples, Vo, k_adh, tau_align, var,m,n

    integer, parameter:: traj_dump_int=100 ! Trajectory file dump interval
    integer, parameter:: status_dump_int=100 ! Status file dump interval
    integer, parameter:: cpt_dump_int=50*traj_dump_int ! Checkpoint file dump interval

    contains
    
    subroutine assign_params(fname)
        character(len=*), intent(in) :: fname
        integer:: fd

        open(newunit=fd,file=fname, access='sequential', form='formatted', status='old', action='read', err=100)
            read(fd, nml=params, err=100, end=100)
        close(fd)
        
        return
100  error stop 'Problem with parameter input file : '//fname
    end subroutine assign_params
end module parameters
