! Help:Begin
! NOTE: This program requires the metadata, last checkpoint and trajectory. Outputs analysis dump.
! Usage: ccd_analysis <metadata file path>
! Help:End

! To minimize the number of local variables, we shall be using the same variable for storing sum and the desired avg
!! as much as possible.
program ccd_analysis
    use analysis
    use utilities, only: help_handler
    implicit none

    character(len=:), allocatable :: metadata_fname
    integer :: metadata_fname_length, exitcode
    integer :: ring, rec_index
    double precision :: msd, alpha2, shapeind, hexop1, hexop2, vicsekop, areafrac, tension, nemop
    double precision :: cell_sd_, cell_area, cell_perim, sum_area, cell_tension
    double precision :: cell_vicsekop_x, cell_vicsekop_y, vicsekop_x, vicsekop_y
    double precision :: cell_major_axis(2), cell_minor_axis(2), cell_nemop_cos_theta

    call help_handler()

    ! Get the metadata file path
    call get_command_argument(1, length=metadata_fname_length, status=exitcode)
    if (exitcode /= 0) error stop 'Fatal: Pass metadata path as argument'
    allocate (character(len=metadata_fname_length) :: metadata_fname)
    call get_command_argument(1, metadata_fname)

    call init(metadata_fname)

    traj_records: do rec_index = 1, recnum
        call traj_read(rec_index, timepoint)

        msd = 0.d0
        alpha2 = 0.d0
        shapeind = 0.d0
        vicsekop_x = 0.d0
        vicsekop_y = 0.d0
        sum_area = 0.d0
        tension = 0.d0
        nemop = 0.d0

        cells: do ring = 1, nrings
            cell_sd_ = cell_sd(ring)
            msd = msd + cell_sd_
            alpha2 = alpha2 + cell_sd_*cell_sd_

            call cell_perimetry(ring, cell_area, cell_perim, cell_tension)
            sum_area = sum_area + cell_area
            tension = tension + cell_tension
            shapeind = shapeind + cell_perim/dsqrt(cell_area)

            call cell_vicsekop(ring, cell_vicsekop_x, cell_vicsekop_y)
            vicsekop_x = vicsekop_x + cell_vicsekop_x
            vicsekop_y = vicsekop_y + cell_vicsekop_y

            ! Nematic-like order parameter from Giavazzi et al. Soft Matter, 2018, 14, Sec.3.4.1 : nemop
            call cell_shape(ring, cell_major_axis, cell_minor_axis)
            cell_nemop_cos_theta = dot_product([cell_vicsekop_x, cell_vicsekop_y], cell_major_axis)
            nemop = nemop + cell_nemop_cos_theta*cell_nemop_cos_theta
        end do cells

        ! Non-gaussian parameter: M. Chiang and D. Marenduzzo, EPL, 116, 10 2016
        alpha2 = alpha2*nrings/(2*msd*msd) - 1.d0
        msd = msd/nrings
        shapeind = shapeind/nrings
        areafrac = sum_area/(box*box)
        tension = tension/nrings
        vicsekop = hypot(vicsekop_x/nrings, vicsekop_y/nrings)
        nemop = 2*(nemop/nrings) - 1.d0

        call psi_6(nrings, hexop1, hexop2)

        call dump(rec_index*traj_dump_int, timepoint, &
                  msd, alpha2, shapeind, hexop1, hexop2, vicsekop, areafrac, tension, nemop)
    end do traj_records

end program ccd_analysis
