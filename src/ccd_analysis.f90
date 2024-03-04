! Help:Begin
! NOTE: This program requires the metadata, last checkpoint and trajectory. Outputs analysis dump.
! Usage: ccd_analysis [--records=<begin>:<end>] [--voronoi] <metadata file path>
! --records : Pass range of records to work with. Omit either <begin> or <end> to assume default. E.g. --records=3:4
! --voronoi : Compute hexatic order of periodic Voronoi tesselations derived from the ring/cell centres.
! Help:End

! To minimize the number of local variables, we shall be using the same variable for storing sum and the desired avg
!! as much as possible.
program ccd_analysis
    use analysis
    use utilities, only: cmd_line_opt, cmd_line_arg, help_handler
    use voronoi, only: periodic_voronoi
    implicit none

    character(len=:), allocatable :: metadata_fname, opt_arg
    integer :: metadata_fname_length, opt_arg_len, exitcode
    integer :: ring, rec_index, begin_rec, end_rec
    double precision :: msd, alpha2, shapeind, hexop1, hexop2, vicsekop, areafrac, tension, nemop
    double precision :: cell_sd_, cell_area, cell_perim, sum_area, cell_tension
    double precision :: cell_vicsekop_x, cell_vicsekop_y, vicsekop_x, vicsekop_y
    double precision :: cell_major_axis(2), cell_minor_axis(2), cell_nemop_cos_theta
    complex, dimension(:), allocatable :: hop ! To hold psi6 for every cell/ring

    call help_handler()

    ! Get the metadata file path
    call cmd_line_arg(1, length=metadata_fname_length)
    if (metadata_fname_length == 0) error stop 'Fatal: Pass metadata path as argument'
    allocate (character(len=metadata_fname_length) :: metadata_fname)
    call cmd_line_arg(1, metadata_fname)

    call setup(metadata_fname)

    ! Sort out the begin and end record number from --records=<from>:<to> cmd line option, if any
    begin_rec = 1 ! default
    end_rec = recnum ! default
    call cmd_line_opt('--records', length=opt_arg_len)
    if (opt_arg_len /= 0) then
        allocate (character(len=opt_arg_len) :: opt_arg)
        call cmd_line_opt('--records', opt_arg)
        read (opt_arg(:scan(opt_arg, ':') - 1), *, iostat=exitcode) begin_rec
        read (opt_arg(scan(opt_arg, ':') + 1:), *, iostat=exitcode) end_rec
        if (begin_rec > end_rec) error stop 'Fatal: Records range provided must be in ascending order'
        if ((begin_rec < 1) .or. (end_rec > recnum)) error stop 'Fatal: Provided records range out of bounds'
        deallocate (opt_arg)
    end if

    call init(begin_rec)

    traj_records: do rec_index = begin_rec, end_rec
        call traj_read(rec_index, timepoint, cmframe=.true.)

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

            call cell_perimetry(ring, k, l0, cell_area, cell_perim, cell_tension)
            sum_area = sum_area + cell_area
            tension = tension + cell_tension
            shapeind = shapeind + cell_perim/dsqrt(cell_area)

            call cell_vicsekop(ring, c, Vo, cell_vicsekop_x, cell_vicsekop_y)
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
        areafrac = min((sum_area + bead_area)/(box*box), 1.d0)
        tension = tension/nrings
        vicsekop = hypot(vicsekop_x/nrings, vicsekop_y/nrings)
        nemop = 2*(nemop/nrings) - 1.d0

        ! hexop1 is from Revalee et al., J. Chem. Phys. 128, 035102 (2008); https://doi.org/10.1063/1.2825300
        ! hexop2 is from Loewe et al., Phy. Rev. Lett. 125(3):038003, 2020
        ! Due to triangle law of complex numbers hexop1 and hexop2 may differ a lot.
        ! hexop1 seems more acceptable to us.
        if (cmd_line_flag('--voronoi')) call periodic_voronoi(cmx, cmy, box)
        hop = psi_6(nrings)
        hexop1 = sum(abs(hop))/nrings
        hexop2 = abs(sum(hop))/nrings

        call dump(rec_index*traj_dump_int, timepoint, &
                  msd, alpha2, shapeind, hexop1, hexop2, vicsekop, areafrac, tension, nemop, poten/nrings)
    end do traj_records

end program ccd_analysis
