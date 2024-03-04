! Help:Begin
! Computes pair correlation histogram. Outputs paircorr.xy. User can optionally provide number of bins (default: 1000).
! NOTE: This program requires the last checkpoint too.
! Usage: ccd_paircorr [--records=<begin>:<end>] [<nbins>]
! --records : Pass range of records to work with. Omit either <begin> or <end> to assume default. E.g. --records=3:4
! Help:End

! Ref: Code 8.1 in Allen and Tildesley's book "Computer Simulation of Liquids".
! Note: Read-in trajectory contains single precision data only. Hence, analysed output is dumped as single precision.

program ccd_paircorr
    use files
    use utilities, only: cmd_line_opt, cmd_line_arg, int_to_char, help_handler
!$  use omp_lib, only: omp_get_max_threads
    implicit none
    integer :: pending_steps, current_step, rec_index, begin_rec, end_rec
    character(len=40) :: params_hash
    character(len=*), parameter :: fname = 'paircorr.xy'
    character(len=:), allocatable :: nbins_arg, opt_arg
    integer :: nbins_arg_len, opt_arg_len, exitcode
    integer :: fd, nbeads_per_ring, nrings, nbins, l, q, bin
    integer, dimension(:), allocatable :: h
    double precision :: cm_dx, cm_dy, r, dr, density, factor

    ! Following variables with trailing _ would be threadprivate
    double precision, dimension(:, :), allocatable :: x_, y_
    real :: timepoint_
    double precision, dimension(:), allocatable :: cmx_, cmy_

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
    allocate (h(nbins))
    h = 0
    density = nrings/(box*box)
    open (newunit=fd, file=fname, access='sequential', form='formatted', status='replace', &
          action='write')
    allocate (x_(nbeads_per_ring, nrings), y_(nbeads_per_ring, nrings), cmx_(nrings), cmy_(nrings))
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

    ! factor below includes all the constants involved in scaling h(k) to get g(r)
    factor = nrings*(end_rec - begin_rec + 1)*dacos(-1.d0)*density*dr*2

!$  write (err_fd, '(a,1x,i0,1x,a)') 'Using', omp_get_max_threads(), 'OpenMP threads'

!$omp parallel do default(shared) private(rec_index, x_, y_, timepoint_, cmx_, cmy_, cm_dx, cm_dy, r, bin, l, q) &
!$omp reduction(+: h)
    do rec_index = begin_rec, end_rec
        call threadsafe_traj_read_xy_only(rec_index, timepoint_, x_, y_)

        do l = 1, nrings
            cmx_(l) = sum(x_(:, l))/nbeads_per_ring
            cmy_(l) = sum(y_(:, l))/nbeads_per_ring
        end do

        do l = 1, nrings - 1
            do q = l + 1, nrings
                cm_dx = cmx_(q) - cmx_(l)
                cm_dy = cmy_(q) - cmy_(l)
                cm_dx = cm_dx - box*nint(cm_dx/box)
                cm_dy = cm_dy - box*nint(cm_dy/box)
                r = hypot(cm_dx, cm_dy)
                bin = floor(r/dr) + 1
                if (bin <= nbins) h(bin) = h(bin) + 2
            end do
        end do
    end do
!$omp end parallel do

    do bin = 1, nbins
        r = (bin - 1)*dr + dr/2 ! Mid point of bin
        write (fd, '(es15.8,1x,es15.8)') r, h(bin)/(factor*r)
    end do

    call close_traj()
    close (fd)
    write (err_fd, '(a,1x,es23.16,1x,a,1x,i0,1x,a,1x,i0,":",i0)') &
        'Used: dr =', dr, 'nbins =', nbins, 'records =', begin_rec, end_rec
    write (err_fd, '(a)') 'Pair correlation g(r) histogram dumped at '//fname
end program ccd_paircorr
