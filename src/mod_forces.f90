module forces
implicit none

contains

!! Subroutine for intracellular forces
!! It can be proved that the sum of these forces for each cell is null from the fact that
!! all the side vectors (r12, r23, ...) sum upto zero for a polygon
	subroutine force

	use shared	
	

	integer:: i,l
        double precision::l1,l2,dx1,dx2,dy1,dy2
        integer :: i_minus_1, i_plus_1          

         
     
    !$omp do private(i,l, l1,l2,dx1,dx2,dy1,dy2, i_minus_1,i_plus_1)
  	do l=1,m
        do i=1,n

                   i_minus_1 = circular_next(i,-1,n)
                   i_plus_1 = circular_next(i,+1,n)
                   
                   dx1 = x(i_minus_1,l)-x(i,l)
                   dy1 = y(i_minus_1,l)-y(i,l)
                   dx2 = x(i,l)-x(i_plus_1,l)
                   dy2 = y(i,l)-y(i_plus_1,l)

                   ! TODO: Is the following necessary?
                   ! dx1 = dx1 - box*nint(dx1/box)
                   ! dx2 = dx2 - box*nint(dx2/box)
                   ! dy1 = dy1 - box*nint(dy1/box)             
                   ! dy2 = dy2 - box*nint(dy2/box)

                   l1 = hypot(dx1,dy1)

                   l2 = hypot(dx2,dy2)

                   fx(i,l)=k*((l1-l0)*dx1/l1 - (l2-l0)*dx2/l2) &
                                - 0.5d0*p*(dy1 + dy2)

                   fy(i,l)=k*((l1-l0)*dy1/l1 - (l2-l0)*dy2/l2) &
                                + 0.5d0*p*(dx1 + dx2)
	  end do
    end do
    !$omp end do nowait
		 

	return
	end subroutine force

 	!!*** Subroutine for intercellular forces of interaction ***!!
	subroutine interaction(store_ring_nb)

	use shared
    use grid_linked_list
    use ring_nb, only: assert_are_nb_rings, init_ring_nb
          
    logical, intent(in) :: store_ring_nb ! flag to store ring-ring neighborhood info
    integer:: i,j,l,q
    double precision:: r,frepx,frepy,dx,dy,fadhx,fadhy,factor
	integer:: icell,jcell,nabor
    integer :: bead_index, other_bead_index
             
    !$omp do private(l)
        do l=1,m
			f_rpx(:,l)=0.0d0
			f_rpy(:,l)=0.0d0
			f_adx(:,l)=0.0d0
			f_ady(:,l)=0.0d0
        end do
    !$omp end do nowait

    if(store_ring_nb) call init_ring_nb()

        !! Loop Over All Cells !!

	!$omp do private(i,j,l,q, r,frepx,frepy,dx,dy,fadhx,fadhy,factor, icell,jcell,nabor, bead_index) &
    !$omp private(other_bead_index) &
    !$omp reduction(+: f_rpx, f_rpy) &
    !$omp reduction(+: f_adx, f_ady)
    grids: do icell=1,ncell
        bead_index = bead_nl_head(icell)

        beads_downlist: do
            if(bead_index == 0) exit beads_downlist

            intra_or_intergrid: do nabor=0,4

                if(nabor == 0) then ! intragrid
                    other_bead_index = bead_nl_body(bead_index)
                else !intergrid
                    jcell = gridmap(4*(icell-1) + nabor)
                    other_bead_index = bead_nl_head(jcell)
                end if
                
                other_beads_downlist: do
                    if(other_bead_index == 0) exit other_beads_downlist

                    not_within_same_ring: if((bead_index-1)/n /= (other_bead_index-1)/n) then

                        l = (bead_index-1)/n + 1 ! Ring index of bead
                        i = mod((bead_index-1), n) + 1 ! Intraring serial number of bead
                        q = (other_bead_index-1)/n + 1 ! Ring index of other bead
                        j = mod((other_bead_index-1), n) + 1 ! Intraring serial number of other bead
                    
							dx = x(j,q)-x(i,l)
							dy = y(j,q)-y(i,l)
							dx = dx - box*nint(dx/box)
							dy = dy - box*nint(dy/box)					        
							r = hypot(dx,dy)
							

                      					within_cutoff: if(r.lt.rc_adh) then

				          		if(store_ring_nb) call assert_are_nb_rings(l, q)

                                if(r.lt.rc_rep) then ! Repulsion
                                factor = k_rep*(r-rc_rep)/r
				          		frepx = factor*dx
					  			frepy = factor*dy

				          			f_rpx(i,l) = f_rpx(i,l) + frepx 
				          			f_rpx(j,q) = f_rpx(j,q) - frepx 

				          			f_rpy(i,l) = f_rpy(i,l) + frepy 
				          			f_rpy(j,q) = f_rpy(j,q) - frepy 

                        				else ! Adhesion
                                
                                factor = k_adh*(rc_adh-r)/r
								fadhx = factor*dx
								fadhy = factor*dy

                                f_adx(i,l) = f_adx(i,l) + fadhx
								f_adx(j,q) = f_adx(j,q) - fadhx

						        f_ady(i,l) = f_ady(i,l) + fadhy 
								f_ady(j,q) = f_ady(j,q) - fadhy

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
