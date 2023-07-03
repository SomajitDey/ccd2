module init
    use shared
    implicit none

contains

!!! Subroutine for random initial configurations
    subroutine initial
        double precision :: radius_l0   ! initial/seed cell radius
        double precision :: mindist, mindist2 ! initial separation between cell centres to ensure non-overlap
        double precision :: radius_eq ! equilibrium radius, achieved by force balance between spring and pressure
        double precision :: dist_eq ! equilibrium distance (centre-centre)

        character(len=:), allocatable :: argument
        integer :: arglen
        double precision :: box_scale, box_scale_percentage

        real :: rands(2)
        integer :: this, other, iter_count, fail_count
        double precision, dimension(size(x, dim=2)) :: xcell, ycell ! centre coordinates of any circular cell/ring
        double precision, dimension(size(x, dim=2)) :: xdisp, ydisp ! displacements during overlap elimination
        double precision :: xcell_this, ycell_this, dx, dy, dr, dr2, disp_rate

! Max trials while seeding, and max iterations while eliminating overlaps
        integer, parameter :: max_trials_seed = 1000, max_iters_emin = 20000

! For stability/convergence during overlap elimination, number of cycles needed for overlap elimination of any pair
        integer, parameter :: cycles_pair_nonoverlap = 100

        logical :: no_overlap_found

        integer :: i, l
        double precision :: angle, m_tmp, mx_tmp, my_tmp, mx_cell, my_cell

! Seed cell properties
        radius_l0 = 0.5d0*l0/dsin(pi/n) ! circumcirle of a regular n-gon with side l0
        mindist = 2*radius_l0 + rc_rep ! rc_rep is the minimum distance between cell peripheries
        mindist2 = mindist*mindist

! If box length is provided by user with --box=<val>, use that. Estimate otherwise.
! To estimate, use zero/negative stress condition where box can accommodate all equilibrium
!! cells without overlap.
! Instead of box, user may also provide a scaling factor in percentage using --box=<scale>%
! In that case, scale the initial box estimate accordingly to get the final box.
        box_scale_percentage = 100.d0
        call cmd_line_opt('--box', length=arglen)
        allocate (character(len=arglen) :: argument)
        call cmd_line_opt('--box', argument)
        if ((argument(arglen:arglen) == '%') .or. arglen == 0) then
            if (argument(arglen:arglen) == '%') read (argument(1:arglen - 1), *, err=100, end=100) box_scale_percentage
100         box_scale = box_scale_percentage/100
            radius_eq = (2*k*dtan(pi/n)*radius_l0)/(2*k*dtan(pi/n) - p)
            if (radius_eq < 0.d0) error stop &
                'Spring-Pressure force balance non-existent. Cannot estimate boxlength. Please provide with --box='
            dist_eq = 2*radius_eq + rc_rep
            box = dsqrt(m*pi*dist_eq*dist_eq/4)*box_scale ! box should accommodate m equilibrium cells
        else
            read (argument, *) box
        end if
        deallocate (argument)

! Seeding the first cell centre at origin
        xcell(1) = 0.d0
        ycell(1) = 0.d0

! Find non-overlapping centres for circular cells by random trial and error
! In case overlap is unavoidable even after so many trials, accept it
        seed_cell_centres: do this = 2, m
            fail_count = 0

            trial: do
                fail_count = fail_count + 1
                if (fail_count > max_trials_seed) exit trial

                call random_number(rands)
                xcell_this = box*rands(1)
                ycell_this = box*rands(2)

                distance_from_other_cells: do other = 1, this - 1
                    dx = xcell_this - xcell(other)
                    dy = ycell_this - ycell(other)
                    dx = dx - box*nint(dx/box)
                    dy = dy - box*nint(dy/box)
                    if ((dx*dx + dy*dy) .lt. mindist2) cycle trial
                end do distance_from_other_cells
                exit trial
            end do trial

            xcell(this) = xcell_this
            ycell(this) = ycell_this
        end do seed_cell_centres

! Eliminating overlaps
        disp_rate = mindist/cycles_pair_nonoverlap ! constant displacement per cycle
        iter_count = 0 ! initial count of iterations
        overlap_elim: do
            no_overlap_found = .true.
            if (iter_count > max_iters_emin) error stop 'Fatal: Took too many cycles for overlap elimination'
            xdisp = 0.d0
            ydisp = 0.d0

            ! Loop over all pairs
            do this = 1, m - 1
                xcell_this = xcell(this)
                ycell_this = ycell(this)
                do other = this + 1, m
                    dx = xcell_this - xcell(other)
                    dy = ycell_this - ycell(other)
                    dx = dx - box*nint(dx/box)
                    dy = dy - box*nint(dy/box)
                    dr2 = dx*dx + dy*dy

                    ! Detect overlap
                    if (dr2 .lt. mindist2) then
                        no_overlap_found = .false.
                        dr = dsqrt(dr2)
                        xdisp(this) = disp_rate*dx/dr ! dx/dr provides direction cosine of unit vector along dr
                        xdisp(other) = -xdisp(this)
                        ydisp(this) = disp_rate*dy/dr ! dy/dr provides direction cosine of unit vector along dr
                        ydisp(other) = -ydisp(this)
                    end if
                end do
            end do

            if (no_overlap_found) exit overlap_elim

            ! Move cell centres towards non-overlap
            xcell = xcell + xdisp
            ycell = ycell + ydisp

            iter_count = iter_count + 1
        end do overlap_elim
        write (err_fd, '(a,1x,i0,1x,a)') 'Overlap eliminated after', iter_count, 'iterations'

! Construct circular cells from the seeded centres
! Also initialize the motility unit vectors randomly for each cell
        do l = 1, m
            call random_number(rands)
            mx_tmp = 2.0d0*rands(1) - 1.0d0
            my_tmp = 2.0d0*rands(2) - 1.0d0
            m_tmp = hypot(mx_tmp, my_tmp)
            mx_cell = mx_tmp/m_tmp
            my_cell = my_tmp/m_tmp
            do i = 1, n
                angle = (i - 1)*2.0d0*pi/n
                x(i, l) = radius_l0*dcos(angle) + xcell(l)
                y(i, l) = radius_l0*dsin(angle) + ycell(l)
                mx(i, l) = mx_cell
                my(i, l) = my_cell
            end do
        end do

        return
    end subroutine initial

end module init
