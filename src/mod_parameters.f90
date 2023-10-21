module parameters
    implicit none
    public

    !! USER PARAMETERS and THEIR DEFAULT VALUES

    !!!! MODEL PARAMETERS
    double precision, protected :: c = 1.0d0         ! c is coeff. of viscous damping
    double precision, protected :: k = 1.0d0      !  Single cell spring constant
    double precision, protected :: l0 = 0.02d0       !  Single cell natural spring-length
    double precision, protected :: p = 0.25d0       !  Single cell internal hydrostatic pressure coefficient
    double precision, protected :: rc_adh = 0.06d0  ! Adhesion interaction cut-off
    double precision, protected :: rc_rep = 0.04d0  ! Repulsion interaction cut-off
    double precision, protected :: ovrlp_trshld = dcos(atan(1.d0)) ! Cell-Cell overlap metric threshold
    double precision, protected :: k_adh = 0.0025d0   !  Adhesion interaction strength
    double precision, protected :: k_rep = 10.0d0 !  Adhesion interaction strength
    double precision, protected :: Vo = 0.0004d0       !  Self propulsion of the beads
    double precision, protected :: dt = 0.1d0   ! Integration timestep
    integer, protected :: tau_align = 1000 ! Timescale for Vicsek alignment in multiples of dt
    integer, protected :: tau_noise = 100 ! Timescale for rotational diffusion in multiples of dt

    !!!! RUN LENGTH
    integer, protected :: nsamples = 100  !! No. of Iterations in terms of traj_dump_int

    !!!! SYSTEM SIZE
    integer, protected :: n = 50    ! No. of beads
    integer, protected :: m = 256   ! No. of cell

    !!!! DUMP INTERVALS
    integer :: traj_dump_int = 100 ! Trajectory file dump interval
    integer :: status_dump_int = 100 ! Status file dump interval
    integer :: cpt_dump_int = 5000 ! Checkpoint file dump interval

    !! END OF USER PARAMETERS

    namelist /params/ c, k, p, l0, rc_adh, rc_rep, k_adh, k_rep, tau_noise, Vo, dt, tau_align, nsamples, n, m
    namelist /params/ traj_dump_int, status_dump_int, cpt_dump_int
    namelist /params/ ovrlp_trshld

    double precision, protected :: noise_strength = 0.0d0 ! Constant coeff. (rot. diff. related) in noise term
    double precision, protected :: align_strength = 0.0d0 ! Constant coeff. in the Vicsek term
    private :: check_params

contains

    subroutine assign_params(fname, nocheck)
        character(len=*), intent(in) :: fname
        logical, intent(in), optional :: nocheck
        integer :: fd

        open (newunit=fd, file=fname, access='sequential', form='formatted', status='old', action='read', err=100)
        read (fd, nml=params, err=100, end=100)
        close (fd)

        if (tau_align /= 0) align_strength = 1.0d0/(tau_align*dt)
        if (tau_noise /= 0) noise_strength = 1.0d0/(tau_noise*dt)

        if (.not. present(nocheck) .or. .not. nocheck) call check_params() !i.e. checking is the default behavior
        return
100     error stop 'Fatal: Problem with parameter input file : '//fname
    end subroutine assign_params

    subroutine check_params()
        use iso_fortran_env, only: err_fd => error_unit
        double precision, parameter :: pi = dacos(-1.0d0)
        double precision :: factor, radius_circular_cell, k_adh_estimate
        logical :: fatal

        fatal = .false.

        write (err_fd, '(a,/,29("="))') 'PARAMETER CONSISTENCY REPORT:'

        write (err_fd, '(/,a)') 'POSSIBLE UNITS:'
        write (err_fd, '(a,1x,f0.3)') 'c =', c
        write (err_fd, '(a,1x,f0.3)') 'k =', k
        write (err_fd, '(a,1x,f0.3)') 'nl0 =', n*l0

        write (err_fd, '(/,a)') 'TIMESCALES:'

        factor = (c/k)/dt
        write (err_fd, '(a,1x,f0.3,1x,a)') 'c/k =', factor, 'dt'
        if (factor < 10.0d0) then
            write (err_fd, '(a)') '**Fatal:  c/k < 10 dt. Upsets assumption of slowly varying force'
            fatal = .true.
        end if

        factor = (c/k_rep)/dt
        write (err_fd, '(a,1x,f0.3,1x,a)') 'c/k_rep =', factor, 'dt'
        if (factor < 10.0d0) then
            write (err_fd, '(a)') '**Fatal:  c/k_rep < 10 dt. Upsets assumption of slowly varying force'
            fatal = .true.
        end if

        factor = (c/k_adh)/dt
        write (err_fd, '(a,1x,f0.3,1x,a)') 'c/k_adh =', factor, 'dt'
        if (factor < 10.0d0) then
            write (err_fd, '(a)') '**Fatal:  c/k_adh < 10 dt. Upsets assumption of slowly varying force'
            fatal = .true.
        end if

        factor = (c/p)/dt
        write (err_fd, '(a,1x,f0.3,1x,a)') 'c/p =', factor, 'dt'
        if (factor < 10.0d0) then
            write (err_fd, '(a)') '**Fatal:  c/p < 10 dt. Upsets assumption of slowly varying force'
            fatal = .true.
        end if

        write (err_fd, '(a,1x,i0,1x,a)') 'tau_align =', tau_align, 'dt'

        write (err_fd, '(a,1x,i0,1x,a)') 'tau_rotational_diffusion =', tau_noise, 'dt'

        write (err_fd, '(/,a)') 'LENGTHSCALES:'

        write (err_fd, '(a,1x,f0.3)') 'l0 =', l0

        if (Vo > epsilon(0.0d0)) then !Only if motility is not insignificant

            factor = l0/(Vo*dt)
            write (err_fd, '(a,1x,f0.3,1x,a)') 'l0 =', factor, 'Vo*dt'
            if (factor < 10.0d0) then
                write (err_fd, '(a)') '**Fatal:  l0 < 10 Vo*dt. Upsets assumption of slowly varying force'
                fatal = .true.
            end if

            factor = rc_rep/(Vo*dt)
            write (err_fd, '(a,1x,f0.3,1x,a)') 'rc_rep =', factor, 'Vo*dt'
            if (factor < 10.0d0) then
                write (err_fd, '(a)') '**Fatal:  rc_rep < 10 Vo*dt. Upsets assumption of slowly varying force'
                fatal = .true.
            end if

            factor = rc_adh/(Vo*dt)
            write (err_fd, '(a,1x,f0.3,1x,a)') 'rc_adh =', factor, 'Vo*dt'
            if (factor < 10.0d0) then
                write (err_fd, '(a)') '**Fatal:  rc_adh < 10 Vo*dt. Upsets assumption of slowly varying force'
                fatal = .true.
            end if

        end if

        factor = rc_rep/l0
        write (err_fd, '(a,1x,f0.3,1x,a)') 'rc_rep =', factor, 'l0'
        if (factor < 1.0d0) write (err_fd, '(a)') '**Warning: l0 > rc_rep &
            &Beads from neighboring cells may interlock or the cells may penetrate each other'

        factor = rc_adh/l0
        write (err_fd, '(a,1x,f0.3,1x,a)') 'rc_adh =', factor, 'l0'
        if (factor < 1.0d0) write (err_fd, '(a)') '**Warning: l0 > rc_adh'

        if (rc_rep > rc_adh) then
            write (err_fd, '(a)') '**Fatal: rc_rep > rc_adh. Defeats the purpose of steric and attractive forces'
            fatal = .true.
        end if

        write (err_fd, '(/,a)') 'SPRING CONSTANT AND PRESSURE SCALES:'

        factor = k_rep/(k*k_adh*p) ! Because k_rep works against k, k_adh and p
        write (err_fd, '(a,1x,f0.3,1x,a)') 'k_rep =', factor, 'k*k_adh*p'
        if (factor < 1.0d0) write (err_fd, '(a)') '**Warning: k_rep < k*k_adh*p'

        factor = p/(2*k*dtan(pi/n)) ! Derived from Free Body Diagram of regular n-gon. factor equals dl/l
        ! Also this measure of p agrees with P = dEnergy/dArea for a regular n-gon that is stretched uniformly
        write (err_fd, '(a,1x,f0.3,1x,a)') 'p =', factor, '2k*tan(pi/n)'
        if (factor > 1.0d0) write (err_fd, '(a)') '**Warning: p > 2k*tan(pi/n)'

        radius_circular_cell = 0.5d0*l0/dsin(pi/n) ! circumcirle of a regular n-gon with side l0
        k_adh_estimate = 0.1d0*k*(radius_circular_cell**2)/(n*(rc_adh - rc_rep))**2
        ! Derived from logic of Eq. 5 in A. Mkrtchyan et al., Soft Matter, 2014, 10, 4332
        !TODO: But where is the effect of p taken into account in the above estimate?
        factor = k_adh/k_adh_estimate
        write (err_fd, '(a,1x,f0.3,1x,a)') 'k_adh =', factor, 'k_adh_estimate [Eq. 5 in Mkrtchyan(Soft Matter, 2014)]'
        if (factor > 10.0d0 .or. factor < 10.0d0) write (err_fd, '(a)') '**Warning: k_adh does not respect estimate'

        if (ovrlp_trshld > 1.d0) then
            write (err_fd, '(a)') '**Fatal:  ovrlp_trshld must be within [0,1]'
            fatal = .true.
        end if

        write (err_fd, '(/,a,/,29("="),/)') 'END OF CONSISTENCY REPORT:' ! Demarcates End of Report

        if (fatal) error stop 'Encountered Fatal conditions above. Exiting...'

    end subroutine check_params
end module parameters
