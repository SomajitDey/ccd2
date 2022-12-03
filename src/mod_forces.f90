module forces
implicit none

contains

!! Subroutine for intracellular forces 
	subroutine force

	use shared	
	

	integer:: i,l
        double precision::l1,l2,dx1,dx2,dy1,dy2
   
          

         
     
    !$omp do private(i,l, l1,l2,dx1,dx2,dy1,dy2)
  	do l=1,m
        ! boundary condition
	
          x(l,0) = x(l,n)
          y(l,0) = y(l,n)
          x(l,n+1) = x(l,1)
          y(l,n+1) = y(l,1)

        do i=1,n

                   dx1 = x(l,i-1)-x(l,i)
                   dy1 = y(l,i-1)-y(l,i)
                   dx2 = x(l,i)-x(l,i+1)
                   dy2 = y(l,i)-y(l,i+1)

                   dx1 = dx1 - box*nint(dx1/box)
                   dx2 = dx2 - box*nint(dx2/box)
                   dy1 = dy1 - box*nint(dy1/box)             
                   dy2 = dy2 - box*nint(dy2/box)

                   l1 = dsqrt(dx1*dx1 + dy1*dy1)

                   l2 = dsqrt(dx2*dx2 + dy2*dy2)

                   fx(l,i)=k*(l1-l0)*dx1/l1 - k*(l2-l0)*dx2/l2 - & 
                           0.5d0*p*l0*(dy1/l1 + dy2/l2) 
                     

                   fy(l,i)=k*(l1-l0)*dy1/l1 - k*(l2-l0)*dy2/l2 + &  
                            0.5d0*p*l0*(dx1/l1 + dx2/l2)

	  end do
    end do
    !$omp end do
		 

	return
	end subroutine force

 	!!*** Subroutine for intercellular forces of interaction ***!!
    ! Below `ring_a` and `ring_b` denote any ring pair within the same cell or grid
    ! `ring_a` and `ring_c` denote any ring pair within neighbouring cells/grids
	subroutine interaction()

	use shared
    use grid_linked_list
          
    integer:: i,j,l,q
    double precision:: r,frepx,frepy,dx,dy,fadhx,fadhy,rlist,factor
	integer:: icell,jcell,jcell0,nabor 
             
    !$omp do private(l)
        do l=1,m
			f_rpx(l,:)=0.0d0
			f_rpy(l,:)=0.0d0
			f_adx(l,:)=0.0d0
			f_ady(l,:)=0.0d0
        end do
    !$omp end do nowait

        !! Loop Over All Cells !!

	!$omp do private(i,j,l,q, r,frepx,frepy,dx,dy,fadhx,fadhy,rlist,factor, icell,jcell,jcell0,nabor)
    grids: do icell=1,ncell

			l=headcell(icell)  !! Highest ring index in a cell(icell)
			
        ring_a: do	 ! For explanation of ring_a see above
            if(l.eq.0) exit ring_a

				i=headbead(icell,l) !! Highest bead index in the part of the ring(l) in a cell(icell)

			ring_a_beads: do	
                if(i.eq.0) exit ring_a_beads

					q=listcell(icell,l)  !! The next lower ring index after l-th ring in a cell(icell) 

            ring_b: do
					if((q.eq.0)) exit ring_b
					       
						j=headbead(icell,q) !! The highest bead index in the q-th ring in a cell(icell) 

                ring_b_beads: do
						if(j.eq.0) exit ring_b_beads
											
							dx = x(q,j)-x(l,i)
							dy = y(q,j)-y(l,i)
							dx = dx - box*nint(dx/box)
							dy = dy - box*nint(dy/box)					        
							r = dsqrt(dx*dx + dy*dy)
							

                      					if(r.lt.rc_rep) then

				          		factor = (r-rc_rep)/r
				          		frepx = factor*dx
					  			frepy = factor*dy

                                !$omp flush (f_rpx, f_rpy)

								!$omp atomic
				          			f_rpx(l,i) = f_rpx(l,i) + frepx 
								!$omp atomic
				          			f_rpx(q,j) = f_rpx(q,j) - frepx 

								!$omp atomic
				          			f_rpy(l,i) = f_rpy(l,i) + frepy 
								!$omp atomic
				          			f_rpy(q,j) = f_rpy(q,j) - frepy 

                        				else if((r.le.rc_adh).and.(r.ge.rc_rep)) then

                                factor = (rc_adh-r)/r
								fadhx = factor*dx
								fadhy = factor*dy

                                !$omp flush (f_adx, f_ady)

								!$omp atomic
                                f_adx(l,i) = f_adx(l,i) + fadhx
								!$omp atomic
								f_adx(q,j) = f_adx(q,j) - fadhx

								!$omp atomic
						        f_ady(l,i) = f_ady(l,i) + fadhy 
								!$omp atomic
								f_ady(q,j) = f_ady(q,j) - fadhy

       							end if

							j=listbead(q,j) !! the next lower bead index in the q-th ring for the same cell(icell)
						end do ring_b_beads

						q=listcell(icell,q) !! the next lower ring index in the same cell(icell)
					end do ring_b



	                      !! Loop Over Neighbouring Cells !!

					jcell0 = 4*(icell-1)     
					
					do nabor=1,4
						jcell = map(jcell0 + nabor)  !!

	!! Loop Over All Rings & Beads Of The Neighbouring Cells!! 

						q=headcell(jcell)
																	
                ring_c: do
						if(q.eq.0) exit ring_c
						    
                                                     if(q.eq.l) then
                                                        !!Considering the next ring in the neighbour cell
                                                        q=listcell(jcell,q)
                                                        cycle ring_c
                                                     end if
	
							j=headbead(jcell,q)

                    ring_c_beads: do
							if(j.eq.0) exit ring_c_beads

							dx = x(q,j)-x(l,i)
							dy = y(q,j)-y(l,i)
							dx = dx - box*nint(dx/box)
							dy = dy - box*nint(dy/box)					        
							r = dsqrt(dx*dx + dy*dy)
							

                      					if(r.lt.rc_rep) then

				          		factor = (r-rc_rep)/r
				          		frepx = factor*dx
					  			frepy = factor*dy

                                !$omp flush (f_rpx, f_rpy)
								!$omp atomic
				          			f_rpx(l,i) = f_rpx(l,i) + frepx 
								!$omp atomic
				          			f_rpx(q,j) = f_rpx(q,j) - frepx 

								!$omp atomic
				          			f_rpy(l,i) = f_rpy(l,i) + frepy 
								!$omp atomic
				          			f_rpy(q,j) = f_rpy(q,j) - frepy 

                        				else if((r.le.rc_adh).and.(r.ge.rc_rep)) then

                                factor = (rc_adh-r)/r
								fadhx = factor*dx
								fadhy = factor*dy

                                !$omp flush (f_adx, f_ady)
								!$omp atomic
								f_adx(l,i) = f_adx(l,i) + fadhx
								!$omp atomic
								f_adx(q,j) = f_adx(q,j) - fadhx

								!$omp atomic
						        f_ady(l,i) = f_ady(l,i) + fadhy 
								!$omp atomic
								f_ady(q,j) = f_ady(q,j) - fadhy
								
       							end if


								j=listbead(q,j)!! Considering the next bead of the current ring in the neighbour cell
							end do ring_c_beads

							q=listcell(jcell,q) !!Considering the next ring in the neighbour cell


						end do ring_c
					end do


					i=listbead(l,i)  !! Considering the next bead of the current ring in the current cell
				end do ring_a_beads

				l=listcell(icell,l)  !!Considering the next ring in the current cell  	
			end do ring_a
	end do grids		
    !$omp end do
	
    !$omp do private(l)
    do l=1,m
    f_adx(l,:) = k_adh*f_adx(l,:)
    f_ady(l,:) = k_adh*f_ady(l,:)
    f_rpx(l,:) = k_rep*f_rpx(l,:)
    f_rpy(l,:) = k_rep*f_rpy(l,:)
    end do
    !$omp end do

	return
        end subroutine interaction            

end module forces
