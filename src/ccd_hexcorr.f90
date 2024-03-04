! Help:Begin
! Computes orientation correlation (g6) histogram. Outputs hexcorr.xy.
! User can optionally provide number of bins (default: 1000).
! NOTE: This program requires the last checkpoint too.
! Usage: ccd_hexcorr [--records=<begin>:<end>] [--voronoi] [<nbins>]
! --records : Pass range of records to work with. Omit either <begin> or <end> to assume default. E.g. --records=3:4
! --voronoi : Compute hexatic order of periodic Voronoi tesselations derived from the ring/cell centres.
! Help:End

! Ref: Code 8.1 in Allen and Tildesley's book "Computer Simulation of Liquids".
! Ref: Jami et al., arXiv:2307.09327v2
! Note: Read-in trajectory contains single precision data only. Hence, analysed output is dumped as single precision.

program ccd_hexcorr
    use files
    use utilities, only: cmd_line_opt, cmd_line_arg, int_to_char, help_handler
    use analysis, only: psi_6
    use voronoi, only: periodic_voronoi
    implicit none
    integer :: pending_steps, current_step, rec_index, begin_rec, end_rec
    character(len=40) :: params_hash
    character(len=*), parameter :: fname = 'hexcorr.xy'
    character(len=:), allocatable :: nbins_arg, opt_arg
    integer :: nbins_arg_len, opt_arg_len, exitcode
    integer :: fd, nbeads_per_ring, nrings, nbins, l, q, bin
    integer, dimension(:), allocatable :: h
    double precision :: cm_dx, cm_dy, r, dr
    double precision, dimension(:), allocatable :: g6
    complex, dimension(:), allocatable :: hop ! To hold psi6 for every cell/ring

    call help_handler()

    ! Get nbins from cmd line arg. Assume default if not provided as arg.
    call cmd_line_arg(1, length=nbins_arg_len)
    if (nbins_arg_len == 0) then
        nbins = 1000 ! Default
    else
        allocate (character(len=nbins_arg_len) :: nbins_arg)
        call cmd_line_arg(1, nbins_arg)
        read (nbins_arg, *) nbins
    end if

    ! Setup and initialize
    call cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)
    nbeads_per_ring = size(x, 1)
    nrings = size(x, 2)
    dr = (box/2)/nbins
    allocate (h(nbins), g6(nbins))
    h = 0
    g6 = 0.d0
    open (newunit=fd, file=fname, access='sequential', form='formatted', status='replace', &
          action='write')
    call open_traj('read', 'old')

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

    timesteps: do rec_index = begin_rec, end_rec
        call traj_read(rec_index, timepoint)
        if (cmd_line_flag('--voronoi')) call periodic_voronoi(cmx, cmy, box)
        hop = psi_6(nrings)

        pairs: do l = 1, nrings - 1
            do q = l + 1, nrings
                cm_dx = cmx(q) - cmx(l)
                cm_dy = cmy(q) - cmy(l)
                cm_dx = cm_dx - box*nint(cm_dx/box)
                cm_dy = cm_dy - box*nint(cm_dy/box)
                r = hypot(cm_dx, cm_dy)
                bin = floor(r/dr) + 1
                if (bin <= nbins) then
                    g6(bin) = g6(bin) + hop(l)*conjg(hop(q)) + hop(q)*conjg(hop(l))
                    h(bin) = h(bin) + 2
                end if
            end do
        end do pairs
    end do timesteps

    bins: do bin = 1, nbins
        r = (bin - 1)*dr + dr/2 ! Mid point of bin
        if (h(bin) /= 0) then
            write (fd, '(es15.8,1x,es15.8)') r, g6(bin)/h(bin)
        else
            write (fd, '(es15.8,1x,es15.8)') r, 0.d0
        end if
    end do bins

    call close_traj()
    close (fd)
    write (err_fd, '(a,1x,es23.16,1x,a,1x,i0,1x,a,1x,i0,":",i0)') &
        'Used: dr =', dr, 'nbins =', nbins, 'records =', begin_rec, end_rec
    write (err_fd, '(a)') 'Orientational correlation g6(r) histogram dumped at '//fname
end program ccd_hexcorr
