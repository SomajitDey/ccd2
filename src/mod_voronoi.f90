module voronoi
    implicit none

contains

    ! Constructs periodic Voronoi tesselation for given 2D square system. Takes points and boxsize.
    ! Outputs neighborlist in a packed form similar to ring_nb_io from module ring_nb.
    ! Optionally, dumps xyfile at provided path for visualization of the Voronoi tesselation using `ccd visual`.
    subroutine periodic_voronoi(x, y, box, nblist, xyfile)
        use utilities, only: mktemp
        use gnuplot, only: dump_cell_xy
        use ring_nb, only: pair_to_index
        double precision, dimension(:), intent(in) :: x, y
        double precision, intent(in) :: box
        integer, dimension(:), intent(out) :: nblist
        character(len=*), intent(in), optional :: xyfile
        character(len=len('/tmp/tmp.XXXXXXXXXX')) :: pts_in, cells_out, nblist_out
        integer :: point, fd_pts, fd_cells, fd_nblist, fd_xy
        character(len=240) :: cmd
        double precision, dimension(:), allocatable :: cell_vx, cell_vy ! Vertices of any single Voronoi cell
        double precision :: vx, vy
        integer :: ios, nsides, i, cell_id, nb_cell_id

        pts_in = mktemp('/tmp')
        cells_out = mktemp('/tmp')
        nblist_out = mktemp('/tmp')

        ! Dump folded x y to pts_in. Also dump z=0 column for compatibility with ccd_voronoi.py interface.
        open (newunit=fd_pts, file=pts_in, status='replace', action='write')
        do point = 1, size(x)
            write (fd_pts, '(3(es23.16,1x))') modulo(x(point), box), modulo(y(point), box), 0.d0
        end do
        flush (fd_pts) ! Make it available for other programs to read without closing file

        write (cmd, '(a,1x,es23.16,4(1x,a))') 'ccd_voronoi.py', box, pts_in, cells_out, '>', nblist_out
        call execute_command_line(cmd)
        close (fd_pts, status='delete')

        ! Prepare the neighborlist array
        open (newunit=fd_nblist, file=nblist_out, status='old', action='read')
        i = 0
        do
            if (i > size(nblist)) exit
            read (fd_nblist, *, iostat=ios) cell_id, nb_cell_id
            if (is_iostat_end(ios)) then
                nblist(i) = 0
            else if (cell_id > nb_cell_id) then
                cycle ! Avoids double counting
            else
                nblist(i) = pair_to_index(cell_id, nb_cell_id)
            end if
            i = i + 1
        end do
        close (fd_nblist, status='delete')

        ! Prepare xyfile to be used by `ccd visual`, if opted for
        open (newunit=fd_cells, file=cells_out, status='old', action='read')
        if (present(xyfile)) then
            allocate (cell_vx(0), cell_vy(0))
            cell_id = 0
            open (newunit=fd_xy, file=xyfile, status='replace', action='write')
            write (fd_xy, '(a,1x,es23.16)') '#Box:', box
            write (fd_xy, '(a)') '#Column headers:'
            write (fd_xy, '(4(a,4x))') 'x', 'y', 'z', 'bead'

            do
                read (fd_cells, *, iostat=ios) vx, vy

                if (ios == 0) then
                    cell_vx = [cell_vx, vx]
                    cell_vy = [cell_vy, vy]
                else if (is_iostat_end(ios)) then
                    exit
                else
                    cell_id = cell_id + 1
                    nsides = size(cell_vx)
                    write (fd_xy, '(//)') ! 2-blank lines to separate each Vornoi cell
                    call dump_cell_xy(fd_xy, box, cell_vx, cell_vy, [(dble(nsides), i=1, nsides)])
                    ! Print comment line containing: cm_x cm_y cell_id
                    write (fd_xy, '(a,1x,es23.16,1x,es23.16,1x,i0)') '#ID@COM', &
                        modulo(sum(cell_vx)/nsides, box), modulo(sum(cell_vy)/nsides, box), cell_id
                    deallocate (cell_vx, cell_vy)
                    allocate (cell_vx(0), cell_vy(0))
                end if

            end do

            close (fd_xy, status='keep')

        end if
        close (fd_cells, status='delete')

    end subroutine periodic_voronoi
end module voronoi
