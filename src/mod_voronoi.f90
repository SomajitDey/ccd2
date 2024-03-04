module voronoi
    implicit none

contains

    ! Constructs Periodic Voronoi Tesselation (PVT) for given 2D square system. Takes points and boxsize.
    ! Prepares cell-cell neighborlist using module ring_nb.
    ! Optionally, dumps xyfile at provided path for visualization of the PVT using `ccd visual`.
    ! If logical optional arg `centroidal` is present and is true, transforms x and y to the centroids of the PVT.
    !! This is useful for getting iterative approximations to centroidal Voronoi tesselation using Lloyd's algorithm.
    subroutine periodic_voronoi(x, y, box, xyfile, centroidal)
        use utilities, only: mktemp
        use gnuplot, only: dump_cell_xy
        use ring_nb, only: reset_ring_nb, assert_are_nb_rings, pack_ring_nb
        double precision, dimension(:), intent(inout) :: x, y
        double precision, intent(in) :: box
        character(len=*), intent(in), optional :: xyfile
        logical, intent(in), optional :: centroidal
        character(len=:), allocatable :: pts_in, cells_out, nblist_out
        logical :: read_cells_out
        integer :: point, fd_pts, fd_cells, fd_nblist, fd_xy
        character(len=240) :: cmd
        double precision, dimension(:), allocatable :: cell_vx, cell_vy ! Vertices of any single Voronoi cell
        double precision :: vx, vy
        integer :: ios, nsides, i, cell_id, nb_cell_id

        pts_in = mktemp('/tmp') ! Input xy(z) file for ccd_voronoy.py
        cells_out = mktemp('/tmp') ! Contains voronoi polytope vertices output by ccd_voronoy.py
        nblist_out = mktemp('/tmp') ! Holds stdout (nblist info) by ccd_voronoy.py
        ! Decide whether to read voronoi polytope vertices coordinates or not
        read_cells_out = present(xyfile) .or. present(centroidal)

        ! Dump folded x y to pts_in. Also dump z=0 column for compatibility with ccd_voronoi.py interface.
        open (newunit=fd_pts, file=pts_in, status='replace', action='write')
        do point = 1, size(x)
            write (fd_pts, '(3(es23.16,1x))') modulo(x(point), box), modulo(y(point), box), 0.d0
        end do
        flush (fd_pts) ! Make it available for other programs to read without closing file

        ! Invoke ccd_voronoy.py to do the main magic
        write (cmd, '(a,1x,es23.16,4(1x,a))') 'ccd_voronoi.py', box, pts_in, cells_out, '>', nblist_out
        call execute_command_line(cmd)
        close (fd_pts, status='delete')

        ! Prepare the neighborlist array by reading nblist_out dumped by ccd_voronoy.py
        call reset_ring_nb()
        open (newunit=fd_nblist, file=nblist_out, status='old', action='read')
        do
            read (fd_nblist, *, iostat=ios) cell_id, nb_cell_id
            if (is_iostat_end(ios)) then
                exit
            else if (cell_id >= nb_cell_id) then
                ! = makes sure initialization with small number of cells work where a cell's image may be its neighbor
                cycle ! Avoids double counting
            else
                call assert_are_nb_rings(cell_id, nb_cell_id)
            end if
        end do
        call pack_ring_nb() ! Preps ring_nb_io and coord_num arrays
        close (fd_nblist, status='delete')

        ! Read cells_out dumped by ccd_voronoi.py if necessary
        open (newunit=fd_cells, file=cells_out, status='old', action='read')
        if (read_cells_out) then
            allocate (cell_vx(0), cell_vy(0))
            cell_id = 0
            if (present(xyfile)) then
                open (newunit=fd_xy, file=xyfile, status='replace', action='write')
                write (fd_xy, '(a,1x,es23.16)') '#Box:', box
                write (fd_xy, '(a)') '#Column headers:'
                write (fd_xy, '(4(a,4x))') 'x', 'y', 'z', 'bead'
            end if

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

                    ! Do the needful with the voronoi polytope vertex coordinates
                    ! E.g. Prepare xyfile to be used by `ccd visual`, if opted for
                    ! E.g. Transform x and y into centroids of the voronoi tesselation, if opted for
                    process_vertices: if (present(xyfile)) then
                        write (fd_xy, '(//)') ! 2-blank lines to separate each Vornoi cell
                        call dump_cell_xy(fd_xy, box, cell_vx, cell_vy, [(dble(nsides), i=1, nsides)])
                        ! Print comment line containing: cm_x cm_y cell_id
                        write (fd_xy, '(a,1x,es23.16,1x,es23.16,1x,i0)') '#ID@COM', &
                            modulo(sum(cell_vx)/nsides, box), modulo(sum(cell_vy)/nsides, box), cell_id
                    else if (present(centroidal) .and. centroidal) then
                        x(cell_id) = modulo(sum(cell_vx)/nsides, box)
                        y(cell_id) = modulo(sum(cell_vy)/nsides, box)
                    end if process_vertices

                    deallocate (cell_vx, cell_vy)
                    allocate (cell_vx(0), cell_vy(0))
                end if

            end do

            if (present(xyfile)) close (fd_xy, status='keep')

        end if
        close (fd_cells, status='delete')

    end subroutine periodic_voronoi
end module voronoi
