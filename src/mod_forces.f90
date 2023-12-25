module forces
    implicit none
    procedure(force_pl), pointer :: force => null() ! alias for the subroutine computing intracellular forces

contains

!! Subroutine for intracellular forces: pressure force is pl
!! It can be proved that the sum of these forces for each cell is null from the fact that
!! all the side vectors (r12, r23, ...) sum upto zero for a polygon
    subroutine force_pl

        use shared

        integer :: i, l
        double precision :: l1, l2, dx1, dx2, dy1, dy2
        integer :: i_minus_1, i_plus_1

!$omp do private(i,l, l1,l2,dx1,dx2,dy1,dy2, i_minus_1,i_plus_1)
        do l = 1, m
            do i = 1, n

                i_minus_1 = circular_next(i, -1, n)
                i_plus_1 = circular_next(i, +1, n)

                dx1 = x(i_minus_1, l) - x(i, l)
                dy1 = y(i_minus_1, l) - y(i, l)
                dx2 = x(i, l) - x(i_plus_1, l)
                dy2 = y(i, l) - y(i_plus_1, l)

                ! TODO: Is the following necessary?
                ! dx1 = dx1 - box*nint(dx1/box)
                ! dx2 = dx2 - box*nint(dx2/box)
                ! dy1 = dy1 - box*nint(dy1/box)
                ! dy2 = dy2 - box*nint(dy2/box)

                l1 = hypot(dx1, dy1)

                l2 = hypot(dx2, dy2)

                fx(i, l) = k*((l1 - l0)*dx1/l1 - (l2 - l0)*dx2/l2) &
                           - 0.5d0*p*(dy1 + dy2)

                fy(i, l) = k*((l1 - l0)*dy1/l1 - (l2 - l0)*dy2/l2) &
                           + 0.5d0*p*(dx1 + dx2)
            end do

            cmx(l) = sum(x(:, l))/n
            cmy(l) = sum(y(:, l))/n
        end do
!$omp end do nowait

        return
    end subroutine force_pl

!! Subroutine for intracellular forces: pressure force is pl0
!! Sum of intracellular forces for each cell must be null
    subroutine force_pl0

        use shared

        integer :: i, l
        double precision :: l1, l2, dx1, dx2, dy1, dy2
        integer :: i_minus_1, i_plus_1
        double precision, dimension(size(x, dim=1)) :: f_bead_x, f_bead_y
        double precision :: f_bead_x_avg, f_bead_y_avg

!$omp do private(i,l, l1,l2,dx1,dx2,dy1,dy2, i_minus_1,i_plus_1, f_bead_x,f_bead_y,f_bead_x_avg,f_bead_y_avg)
        do l = 1, m
            f_bead_x_avg = 0.d0
            f_bead_y_avg = 0.d0

            do i = 1, n

                i_minus_1 = circular_next(i, -1, n)
                i_plus_1 = circular_next(i, +1, n)

                dx1 = x(i_minus_1, l) - x(i, l)
                dy1 = y(i_minus_1, l) - y(i, l)
                dx2 = x(i, l) - x(i_plus_1, l)
                dy2 = y(i, l) - y(i_plus_1, l)

                ! TODO: Is the following necessary?
                ! dx1 = dx1 - box*nint(dx1/box)
                ! dx2 = dx2 - box*nint(dx2/box)
                ! dy1 = dy1 - box*nint(dy1/box)
                ! dy2 = dy2 - box*nint(dy2/box)

                l1 = hypot(dx1, dy1)

                l2 = hypot(dx2, dy2)

                f_bead_x(i) = k*((l1 - l0)*dx1/l1 - (l2 - l0)*dx2/l2) &
                              - 0.5d0*p*l0*(dy1/l1 + dy2/l2)

                f_bead_x_avg = f_bead_x_avg + f_bead_x(i)

                f_bead_y(i) = k*((l1 - l0)*dy1/l1 - (l2 - l0)*dy2/l2) &
                              + 0.5d0*p*l0*(dx1/l1 + dx2/l2)

                f_bead_y_avg = f_bead_y_avg + f_bead_y(i)
            end do

            f_bead_x_avg = f_bead_x_avg/n
            f_bead_y_avg = f_bead_y_avg/n

            fx(:, l) = f_bead_x(:) - f_bead_x_avg ! Because total intra force for any cell/ring must be zero
            fy(:, l) = f_bead_y(:) - f_bead_y_avg ! e.g. a live amoeba in vacuum can have no net movement

            cmx(l) = sum(x(:, l))/n
            cmy(l) = sum(y(:, l))/n
        end do
!$omp end do nowait

        return
    end subroutine force_pl0

         !!*** Subroutine for intercellular forces of interaction ***!!
    subroutine interaction(store_ring_nb)

        use shared
        use grid_linked_list
        use ring_nb, only: assert_are_nb_rings

        logical, intent(in) :: store_ring_nb ! flag to store ring-ring neighborhood info
        integer :: i, j, l, q
        double precision :: r, frepx, frepy, dx, dy, fadhx, fadhy, factor
        double precision :: cm_d, cm_dx, cm_dy, overlap
        integer :: icell, jcell, nabor
        integer :: bead_index, other_bead_index

!$omp do private(l)
        do l = 1, m
            f_rpx(:, l) = 0.0d0
            f_rpy(:, l) = 0.0d0
            f_adx(:, l) = 0.0d0
            f_ady(:, l) = 0.0d0
        end do
!$omp end do nowait

        !! Loop Over All Cells !!

!$omp do private(i,j,l,q, r,frepx,frepy,dx,dy,fadhx,fadhy,factor, icell,jcell,nabor, bead_index) &
!$omp private(other_bead_index) &
!$omp private(cm_d, cm_dx, cm_dy, overlap) &
!$omp reduction(+: f_rpx, f_rpy) &
!$omp reduction(+: f_adx, f_ady)
        grids: do icell = 1, ncell
            bead_index = bead_nl_head(icell)

            beads_downlist: do
                if (bead_index == 0) exit beads_downlist

                intra_or_intergrid: do nabor = 0, 4

                    if (nabor == 0) then ! intragrid
                        other_bead_index = bead_nl_body(bead_index)
                    else !intergrid
                        jcell = gridmap(4*(icell - 1) + nabor)
                        other_bead_index = bead_nl_head(jcell)
                    end if

                    other_beads_downlist: do
                        if (other_bead_index == 0) exit other_beads_downlist

                        not_within_same_ring: if ((bead_index - 1)/n /= (other_bead_index - 1)/n) then

                            l = (bead_index - 1)/n + 1 ! Ring index of bead
                            i = mod((bead_index - 1), n) + 1 ! Intraring serial number of bead
                            q = (other_bead_index - 1)/n + 1 ! Ring index of other bead
                            j = mod((other_bead_index - 1), n) + 1 ! Intraring serial number of other bead

                            ! Get position vector joining images of i and j that are nearest to each other
                            dx = x(j, q) - x(i, l)
                            dy = y(j, q) - y(i, l)
                            dx = dx - box*nint(dx/box) ! Folding w.r.t PBC (minimum image convention)
                            dy = dy - box*nint(dy/box)
                            r = hypot(dx, dy)

                            within_cutoff: if (r .lt. rc_adh) then

                                if (store_ring_nb) call assert_are_nb_rings(l, q)

                                if (r .lt. rc_rep) then ! Repulsion

                                    ! Get position vector joining nearest images of c.o.ms of cells q and l
                                    cm_dx = cmx(q) - cmx(l)
                                    cm_dy = cmy(q) - cmy(l)
                                    cm_dx = cm_dx - box*nint(cm_dx/box)
                                    cm_dy = cm_dy - box*nint(cm_dy/box)
                                    cm_d = hypot(cm_dx, cm_dy)

                                    ! Overlap metric is the projection of bead-bead vector along centre-centre
                                    !! joining unit vector, scaled by rc_rep to make the metric dimensionless
                                    overlap = (cm_dx*dx + cm_dy*dy)/(rc_rep*cm_d)
                                    ! i-j joining position vector [dx,dy] makes obtuse angle with the position vector
                                    !! joining the com's [cm_dx,cm_dy] only when i and j penetrate each other's cell.
                                    !! Negativity of the above dot product thus indicates cell-cell overlap.

                                    factor = k_rep*(r - rc_rep)/r
                                    frepx = factor*dx
                                    frepy = factor*dy

                                    ! In case of imminent overlap, push cells away with extra force
                                    if (overlap < ovrlp_trshld) then
                                        factor = -k_rep*rc_rep/cm_d
                                        frepx = frepx + factor*cm_dx
                                        frepy = frepy + factor*cm_dy
                                    end if

                                    f_rpx(i, l) = f_rpx(i, l) + frepx
                                    f_rpx(j, q) = f_rpx(j, q) - frepx

                                    f_rpy(i, l) = f_rpy(i, l) + frepy
                                    f_rpy(j, q) = f_rpy(j, q) - frepy

                                else ! Adhesion

                                    factor = k_adh*(rc_adh - r)/r
                                    fadhx = factor*dx
                                    fadhy = factor*dy

                                    f_adx(i, l) = f_adx(i, l) + fadhx
                                    f_adx(j, q) = f_adx(j, q) - fadhx

                                    f_ady(i, l) = f_ady(i, l) + fadhy
                                    f_ady(j, q) = f_ady(j, q) - fadhy

                                end if
                            end if within_cutoff

                        end if not_within_same_ring

                        other_bead_index = bead_nl_body(other_bead_index)
                    end do other_beads_downlist

                end do intra_or_intergrid

                bead_index = bead_nl_body(bead_index)
            end do beads_downlist

        end do grids
!$omp end do

        return
    end subroutine interaction

end module forces
