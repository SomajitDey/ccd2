! This module contains procedures that you can use for your custom (analysis) codes. Consequently,
!! they are independent, pure/elemental as far as possible. Procedures prefixed with cell_ compute
!! properties of the given cell/ring only. To get average over all cells, call cell_* procedures
!! inside a loop over cells. This eliminates unnecessary use of multiple loops.

!NOTE: When adding new procedures to this module try to make them pure/elemental and independent.
!! To retain independence, quantities such as beads per cell may be derived from size of x.

module analysis
    use files
    use parameters
    implicit none
    integer :: nrings, nbeads_per_ring
    double precision, dimension(:), allocatable :: init_xcm, init_ycm
    double precision :: bead_area
    character(len=*), parameter :: analysis_dump_fname = 'analysis.txt'
    integer :: analysis_dump_fd

contains

    ! Reading in parameters, Opening traj file, analysis dump file etc.
    ! Provide metadata file for reading in parameters as the original parameter file may not be present
    subroutine setup(metadata_fname)
        character(len=*), intent(in) :: metadata_fname
        integer :: pending_steps, current_step
        character(len=40) :: params_hash

        call assign_params(fname=metadata_fname, nocheck=.true.)

        call cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)
        call open_traj('read', 'old')

        ! Set nrings and nbeads_per_ring info for all other analysis routines to use
        nrings = size(x, 2)
        nbeads_per_ring = size(x, 1)

        ! Area of all the beads (excluded volume) not taken into account by cell_perimetry().
        bead_area = (nrings*nbeads_per_ring)*0.5d0*(dacos(-1.d0)*rc_rep*rc_rep/4)

        allocate (init_xcm(nrings), init_ycm(nrings))

        open (newunit=analysis_dump_fd, file=analysis_dump_fname, status='replace')

        ! Write the column headers in analysis dump file
        write (analysis_dump_fd, '(11(a,2x))') 'frame', 'time', 'msd', 'alpha2', &
            'shapeind', 'hexop1', 'hexop2', 'vicsekop', 'areafrac', 'tension', 'nemop'
    end subroutine setup

    ! Initialize
    subroutine init(begin_rec)
        integer, intent(in) :: begin_rec

        call traj_read(begin_rec, timepoint, cmframe=.true.)

        init_xcm = cmx
        init_ycm = cmy
    end subroutine init

    ! Dump analysis results
    subroutine dump(rec_index, time, msd, alpha2, shapeind, hexop1, hexop2, vicsekop, areafrac, tension, nemop)
        integer, intent(in) :: rec_index
        real, intent(in) :: time
        double precision, intent(in) :: msd, alpha2, shapeind, hexop1, hexop2, vicsekop, areafrac, tension, nemop
        write (analysis_dump_fd, '(i0,2x,10(es23.16,2x))') &
            rec_index, time, msd, alpha2, shapeind, hexop1, hexop2, vicsekop, areafrac, tension, nemop
    end subroutine dump

    ! Takes cell/ring index and outputs square deviation (of cm) for that cell, for the current state
    elemental double precision function cell_sd(ring)
        integer, intent(in) :: ring
        double precision :: dxcm, dycm

        dxcm = init_xcm(ring) - cmx(ring)
        dycm = init_ycm(ring) - cmy(ring)
        cell_sd = dxcm*dxcm + dycm*dycm
    end function cell_sd

    ! Takes cell/ring index as well as k and l0. Outputs perimetry info
    pure subroutine cell_perimetry(ring, k, l0, area, perimeter, tension)
        use utilities, only: circular_next
        integer, intent(in) :: ring
        double precision, intent(in) :: k, l0 ! Explicit intent in for retaining independence of this procedure
        double precision, intent(out) :: area, perimeter, tension
        integer :: this_bead, next_bead, nbeads_per_ring
        double precision :: x_this_bead, y_this_bead, x_next_bead, y_next_bead, l

        nbeads_per_ring = size(x, 1) ! Independent derivation based on x without depending on setup() call

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
    pure subroutine cell_shape(ring, major_axis, minor_axis)
        use utilities, only: eigen_2x2_mat
        integer, intent(in) :: ring
        double precision, dimension(:), intent(out) :: major_axis, minor_axis
        double precision :: xcm, ycm, shape_tensor(2, 2), eval1, eval2, evec1(2), evec2(2)
        integer :: nbeads_per_ring

        nbeads_per_ring = size(x, 1) ! Independent derivation based on x without depending on setup() call

        xcm = cmx(ring)
        ycm = cmy(ring)
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
    pure subroutine cell_vicsekop(ring, c, Vo, vopx, vopy)
        integer, intent(in) :: ring
        double precision, intent(in) :: c, Vo ! Explicit intent in for retaining independence of this procedure
        double precision, intent(out) :: vopx, vopy
        double precision :: norm
        integer :: nbeads_per_ring

        nbeads_per_ring = size(x, 1) ! Independent derivation based on x without depending on setup() call

        vopx = sum((fx(:, ring) + f_adx(:, ring) + f_rpx(:, ring))/c + Vo*mx(:, ring))/nbeads_per_ring
        vopy = sum((fy(:, ring) + f_ady(:, ring) + f_rpy(:, ring))/c + Vo*my(:, ring))/nbeads_per_ring
        norm = hypot(vopx, vopy)

        vopx = vopx/norm
        vopy = vopy/norm
    end subroutine cell_vicsekop

    ! Hexatic/Bond-orientational order parameter: PSI_6.
    ! Ref: Jami et al., arXiv:2307.09327v2 (Look for definition of capital \PSI_6)
    function psi_6(nrings)
        use ring_nb, only: index_to_pair
        integer, intent(in) :: nrings ! Number of rings/cells
        integer :: i, valu, ring1, ring2
        double precision :: re, im, xcm_ring1, ycm_ring1, xcm_ring2, ycm_ring2
        complex, dimension(nrings) :: psi_6
        complex :: z ! Just to store any complex value

        psi_6 = (0.0, 0.0)

        ! Looping over neighboring pairs, i.e. elements of ring_nb_io, is faster than loop over all possible pairs
        !! and asking if they are are_nb_rings()
        do i = 1, size(ring_nb_io)
            valu = ring_nb_io(i)
            if (valu == 0) exit
            call index_to_pair(valu, ring1, ring2)
            xcm_ring1 = cmx(ring1)
            ycm_ring1 = cmy(ring1)
            xcm_ring2 = cmx(ring2)
            ycm_ring2 = cmy(ring2)
            re = xcm_ring2 - xcm_ring1
            im = ycm_ring2 - ycm_ring1
            z = cmplx(re, im)**6; z = z/abs(z) ! gives e(i6theta)
            psi_6(ring1) = psi_6(ring1) + z
            psi_6(ring2) = psi_6(ring2) + conjg(z)
        end do

        psi_6 = psi_6/coord_num
    end function psi_6

end module analysis
