!TODO: parallelize (OMP) later
! Many subroutines/functions are prefixed with cell_. This means they compute cell properties only.
! To get average over all cells/rings, use them inside a loop over cells. This is to eliminate
! multiple loops.

module analysis
    use files
    use parameters
    implicit none
    integer :: nrings, nbeads_per_ring
    double precision, dimension(:), allocatable :: init_xcm, init_ycm
    double precision :: init_sys_xcm, init_sys_ycm
    character(len=*), parameter :: analysis_dump_fname = 'analysis.txt'
    integer :: analysis_dump_fd

contains

    ! Initialization/Setup: Reading in parameters, Opening traj file, analysis dump file etc.
    ! Provide metadata file for reading in parameters as the original parameter file may not be present
    subroutine init(metadata_fname)
        character(len=*), intent(in) :: metadata_fname
        integer :: pending_steps, current_step, ring
        character(len=40) :: params_hash

        call assign_params(fname=metadata_fname, nocheck=.true.)

        call cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)
        call open_traj('read', 'old')

        ! Set nrings and nbeads_per_ring info for all other analysis routines to use
        nrings = size(x, 2)
        nbeads_per_ring = size(x, 1)

        allocate (init_xcm(nrings), init_ycm(nrings))
        call traj_read(1, timepoint)

        do ring = 1, nrings
            call cell_cm(ring, init_xcm(ring), init_ycm(ring))
        end do

        init_sys_xcm = sum(x)/size(x)
        init_sys_ycm = sum(y)/size(y)

        open (newunit=analysis_dump_fd, file=analysis_dump_fname, status='replace')

        ! Write the column headers in analysis dump file
        write (analysis_dump_fd, '(11(a,2x))') &
            'rec', 'time', 'msd', 'alpha2', 'shapeind', 'hexop1', 'hexop2', 'vicsekop', 'areafrac', 'tension', 'nemop'
    end subroutine init

    ! Dump analysis results
    subroutine dump(rec_index, time, msd, alpha2, shapeind, hexop1, hexop2, vicsekop, areafrac, tension, nemop)
        integer, intent(in) :: rec_index
        real, intent(in) :: time
        double precision, intent(in) :: msd, alpha2, shapeind, hexop1, hexop2, vicsekop, areafrac, tension, nemop
        write (analysis_dump_fd, '(i0,2x,10(es23.16,2x))') &
            rec_index, time, msd, alpha2, shapeind, hexop1, hexop2, vicsekop, areafrac, tension, nemop
    end subroutine dump

    ! Outputs global/system centre of mass square deviation
    double precision function sys_sd(init_sys_xcm, init_sys_ycm)
        double precision, intent(in) :: init_sys_xcm, init_sys_ycm
        double precision :: sys_xcm_dev, sys_ycm_dev

        sys_xcm_dev = sum(x)/size(x) - init_sys_xcm
        sys_ycm_dev = sum(y)/size(y) - init_sys_ycm

        sys_sd = (sys_xcm_dev*sys_xcm_dev) + (sys_ycm_dev*sys_ycm_dev)
    end function sys_sd

    ! Takes cell/ring index and outputs cm coordinates
    subroutine cell_cm(ring, xcm, ycm)
        integer, intent(in) :: ring
        double precision, intent(out) :: xcm, ycm

        xcm = sum(x(:, ring))/nbeads_per_ring
        ycm = sum(y(:, ring))/nbeads_per_ring
    end subroutine cell_cm

    ! Takes cell/ring index and outputs square deviation (of cm) for that cell, for the current state
    double precision function cell_sd(ring)
        integer, intent(in) :: ring
        double precision :: xcm, ycm, dxcm, dycm

        call cell_cm(ring, xcm, ycm)
        dxcm = init_xcm(ring) - xcm
        dycm = init_ycm(ring) - ycm
        cell_sd = dxcm*dxcm + dycm*dycm
    end function cell_sd

    ! Takes cell/ring index and outputs perimetry info
    subroutine cell_perimetry(ring, area, perimeter, tension)
        use utilities, only: circular_next
        integer, intent(in) :: ring
        double precision, intent(out) :: area, perimeter, tension
        integer :: this_bead, next_bead
        double precision :: x_this_bead, y_this_bead, x_next_bead, y_next_bead, l

        area = 0.0d0
        perimeter = 0.0d0
        tension = 0.0d0

        do this_bead = 1, nbeads_per_ring
            x_this_bead = x(this_bead, ring)
            y_this_bead = y(this_bead, ring)

            next_bead = circular_next(this_bead, +1, nbeads_per_ring)
            x_next_bead = x(next_bead, ring)
            y_next_bead = y(next_bead, ring)

            area = area + x_this_bead*y_next_bead - x_next_bead*y_this_bead

            l = hypot(x_this_bead - x_next_bead, y_this_bead - y_next_bead)
            perimeter = perimeter + l

            tension = tension + dabs(l - l0)
        end do
        area = dabs(0.5d0*area)
        tension = k*tension/nbeads_per_ring
    end subroutine cell_perimetry

    ! Takes cell/ring index and outputs its major and minor axis unit vectors
    subroutine cell_shape(ring, major_axis, minor_axis)
        use utilities, only: eigen_2x2_mat
        integer, intent(in) :: ring
        double precision, dimension(:), intent(out) :: major_axis, minor_axis
        double precision :: xcm, ycm, shape_tensor(2, 2), eval1, eval2, evec1(2), evec2(2)

        call cell_cm(ring, xcm, ycm)
        shape_tensor(1, 1) = sum((x(:, ring) - xcm)*(x(:, ring) - xcm))/nbeads_per_ring
        shape_tensor(2, 2) = sum((y(:, ring) - ycm)*(y(:, ring) - ycm))/nbeads_per_ring
        shape_tensor(2, 1) = sum((y(:, ring) - ycm)*(x(:, ring) - xcm))/nbeads_per_ring
        shape_tensor(1, 2) = shape_tensor(2, 1)

        call eigen_2x2_mat(shape_tensor, eval1, eval2, evec1, evec2)

        if (eval1 > eval2) then
            major_axis = evec1
            minor_axis = evec2
        else
            major_axis = evec2
            minor_axis = evec1
        end if
    end subroutine cell_shape

    ! Takes cell/ring index and outputs its vicsek order parameter contribution
    subroutine cell_vicsekop(ring, vopx, vopy)
        integer, intent(in) :: ring
        double precision, intent(out) :: vopx, vopy
        double precision :: norm

        vopx = sum((fx(:, ring) + f_adx(:, ring) + f_rpx(:, ring))/c + Vo*mx(:, ring))/nbeads_per_ring
        vopy = sum((fy(:, ring) + f_ady(:, ring) + f_rpy(:, ring))/c + Vo*my(:, ring))/nbeads_per_ring
        norm = hypot(vopx, vopy)

        vopx = vopx/norm
        vopy = vopy/norm
    end subroutine cell_vicsekop

    ! Hexatic/Bond-orientational order parameter: h.o.p. It is defined in two ways.
    ! hexop1 is from Revalee et al., J. Chem. Phys. 128, 035102 (2008); https://doi.org/10.1063/1.2825300
    ! hexop2 is from Loewe et al., Phy. Rev. Lett. 125(3):038003, 2020
    ! Due to triangle law of complex numbers hexop1 and hexop2 may differ a lot.
    ! hexop1 seems more acceptable to us.
    subroutine psi_6(nrings, hexop1, hexop2)
        use ring_nb, only: are_nb_rings
        integer, intent(in) :: nrings ! Number of rings/cells
        double precision, intent(out) :: hexop1, hexop2
        integer :: ring1, ring2 ! ring/cell index
        double precision :: re, im, xcm_ring1, ycm_ring1, xcm_ring2, ycm_ring2
        complex, dimension(nrings) :: hop_z_sum ! Stores the complex sum for every cell/ring
        integer, dimension(nrings) :: num_nb ! Number of nearest neighbors
        complex :: z ! Just to store any complex value

        hop_z_sum = (0.0, 0.0)
        num_nb = 0

        do ring1 = 1, nrings - 1
            call cell_cm(ring1, xcm_ring1, ycm_ring1)
            do ring2 = ring1 + 1, nrings
                if (are_nb_rings(ring1, ring2)) then
                    call cell_cm(ring2, xcm_ring2, ycm_ring2)
                    re = xcm_ring2 - xcm_ring1
                    im = ycm_ring2 - ycm_ring1
                    z = cmplx(re, im)**6; z = z/abs(z) ! gives e(i6theta)
                    hop_z_sum(ring1) = hop_z_sum(ring1) + z
                    num_nb(ring1) = num_nb(ring1) + 1
                    hop_z_sum(ring2) = hop_z_sum(ring2) + conjg(z)
                    num_nb(ring2) = num_nb(ring2) + 1
                end if
            end do
        end do

        hexop1 = sum(abs(hop_z_sum/num_nb))/nrings
        hexop2 = abs(sum(hop_z_sum/num_nb)/nrings)
    end subroutine psi_6

end module analysis
