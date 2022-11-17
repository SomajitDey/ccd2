

       module parameters

       implicit none
	integer,parameter:: n = 50    ! No. of beads
	integer,parameter:: m = 256   ! No. of cell
    integer,parameter:: s = n+1
    double precision,parameter:: box = 46.0d0   !  Box length
	double precision,parameter:: rcut= 1.50d0 
	double precision,parameter:: k=120.0d0      !  Single cell spring constant
	double precision,parameter:: p=25.0d0       !  Single cell internal hydrostatic pressure coefficient
	double precision,parameter:: l0=0.1d0       !  Single cell natural spring-length
	double precision,parameter:: rc_adh=0.28d0  ! Adhesion interaction cut-off
	double precision,parameter:: rc_rep=0.18d0  ! Repulsion interaction cut-off
	double precision,parameter:: k_adh=0.001d0   !  Adhesion interaction strength
	double precision,parameter:: k_rep=1000.0d0 !  Adhesion interaction strength
	double precision,parameter:: mean=0.0d0     !  Mean of the gaussian white noise
	double precision,parameter:: var=0.05d0      !  Variance of the gaussian white noise
	double precision,parameter:: Vo=0.05d0       !  Self propulsion of the beads
    double precision:: tau_align ! Tau for Vicsek alignment
       end module parameters


       module position_arrays

       use parameters
 
       implicit none
       double precision:: x1(m),y1(m),x(m,0:s),y(m,0:s)
       double precision:: mx(m,0:s),my(m,0:s) ! denotes motility unit vector for every bead
	
       end module position_arrays


       module forces
         
       use parameters

       implicit none
       double precision:: fx(m,n),fy(m,n),f_intx(m,n),f_inty(m,n),frpx(m,n),frpy(m,n),fadx(m,n),fady(m,n)
     
       end module forces


       module map_arrays
         
       use parameters

       implicit none
	integer,parameter:: w=int(box/rcut), ncell=w*w
	integer,parameter:: mapsiz=4*ncell
	integer:: map(mapsiz)
 
       end module map_arrays


       module list_arrays
         
       use map_arrays

       implicit none
	integer:: headcell(ncell),headbead(ncell,m),listcell(ncell,m),listbead(m,n),icell_memory(m,n)
 
       end module list_arrays


        program many_cell       ! Main Program Starts

	use parameters
	use position_arrays
	use forces
	use map_arrays
	use list_arrays

	integer:: icell,jcell,l,i,j1,jf,q,j
	double precision:: r,cell,celli,dt,x2(m,n),y2(m,n),ti,tf,frepx,frepy,dx,dy,rint,h,sys_xcm,sys_ycm
	Character(len=len('With_Noise/Coords_shift_wrt_SysCom/Adhesions/Diff_noise_var/Var0.2/Kadh20/Ensmbls_for_alpha2/')) &
	& ::filepath

	filepath = 'With_Noise/Coords_shift_wrt_SysCom/Adhesions/Diff_noise_var/vo0.05_var0.05/'
	
	! open(14,file='Grid1.5_P2.5/Fin_5_10_4_highadh9_2.0_dt0.005_krep100_cnf_5.dat',status='unknown')
	 open(16,file=trim(adjustl(filepath))//'adh0.001_1.5_Noise(var0.05_v_0.05)P2.5_Kspr120_ite_1E7_cnf_2.dat',&
                   & status='unknown',position='append')
         !open(18,file=trim(adjustl(filepath))//'adh0.01_1.5_Noise(var0.5_v_0.2)P4.0_Kspr120_ite_1E7_cnf_1_Fin.dat',&
        ! & status='unknown')


        call cpu_time(ti)
	r=1.0d0
	dt=0.001d0
    tau_align=dt*10
	jf=10050000  !! No. of Iterations
	!jf=1
	call initial(r)   
          do l=1,M  
            do i=1,N
           
               write(16,*)l,i,x(l,i),y(l,i)
                
            end do
          end do
      
	call maps
        call initial_angle
	!write(*,*)'No. of boxes',ncell
	do j1=1,jf

	if(mod(j1,1).eq.0) then  !! Step to calculate all the coordinates w.r.t. System COM.
	sys_xcm = 0.0d0
	sys_ycm = 0.0d0
	do l=1,M
		do i=1,N
			sys_xcm = sys_xcm + x(l,i)
			sys_ycm = sys_ycm + y(l,i)
		end do
	end do
	sys_xcm = sys_xcm/(M*N)
	sys_ycm = sys_ycm/(M*N)
	x = x - sys_xcm
 	y = y - sys_ycm
	end if

        call links(j1,jf)
	call force(j1)
	call interaction(j1,frepx,frepy,dx,dy,rint,icell,jcell,l)

	if(j1.gt.50000) then
							
	call move_noise(dt)

        if(mod(j1,5000).eq.0) then

          do l=1,M  
            do i=1,N
           
               write(16,*)j1-50000,l,i,x(l,i),y(l,i)
                
            end do
          end do
        end if


	else

	call move_deterministic(dt)

        if(j1.eq.50000) then

          do l=1,M  
            do i=1,N
           
               write(16,*)1,l,i,x(l,i),y(l,i)
             
            end do
          end do
        end if


	end if


       if(j1.eq.jf) then
          do l=1,M  
            do i=1,N        
				   !x(l,i) = x(l,i) - box*floor(x(l,i)/box)
				   !y(l,i) = y(l,i) - box*floor(y(l,i)/box)
                !write(18,*)j1,l,i,x(l,i),y(l,i)	              
            end do
          end do
       end if
	
	end do

	call cpu_time(tf)

	write(*,*)'time=',tf-ti

        end program many_cell       ! Main Program Ends


	!!*** Subroutine to set up the cells and to make the list of neighbouring cells through PBC ***!!
	subroutine maps

        use map_arrays

	implicit none
	integer:: ix,iy,imap,icell

	!! Statement Function To Give Cell Index !!

	icell(ix,iy) = 1 + mod(ix-1+w , w) + mod(iy-1+w , w) * w

	!! Find Half The Nearest Neighbours Of Each Cell !!

	do iy=1,w
        	do ix=1,w

			
			imap = (icell(ix,iy) - 1) * 4
		
			map(imap+1) = icell(ix+1,iy)
			map(imap+2) = icell(ix+1,iy+1)
			map(imap+3) = icell(ix,iy+1)
			map(imap+4) = icell(ix-1,iy+1)

			!write(11,*)icell(ix,iy),map(imap+1),map(imap+2),map(imap+3),map(imap+4)
		end do
	end do	

	return
	end


	!!*** Subroutine to make linked lists & head of chain arrays ***!!
	subroutine links(j1,jf) 

	use position_arrays
	use list_arrays
	use map_arrays	

	implicit none
	integer:: j1,jf,icell,l,i,headcell_dummy(ncell)
	double precision:: cell,celli,x3(m,n),y3(m,n),delta
	
   	!! Zero Head Of Chain Array & List Array !!


	headcell      =0 !! head of chain array for the rings in a cell
	headcell_dummy=0
	listcell      =0 !! List array for the rings in a cell
	headbead      =0 !! head of chain array for the beads in a ring in the cell
	listbead      =0

	delta = 0.0d0

	celli = dble(w/box) !! celli is the inverse of cell length
	cell  = 1.0d0/celli !! cell is the cell length

	!if(j1.eq.1) write(*,*) 'cell length=',cell,'No.of cells=',ncell,'inverse length=',celli
	!write(*,*) 'cell length=',cell,'No.of cells=',ncell,'inverse length=',celli

	if(cell.lt.rcut) then
		stop 'cell size to small for cut off'
	end if

	!! Sort All Beads !!

	do l=1,m	 !! Loop over rings
		do i=1,n !! Loop over beads

		        x3(l,i) = x(l,i) - box*floor(x(l,i)/box)
                        y3(l,i) = y(l,i) - box*floor(y(l,i)/box)


8			icell = 1 + int(x3(l,i) * celli) + int(y3(l,i) * celli) * w  !! determining the cell index for a particular bead

			icell_memory(l,i) = icell
			
			
10			listcell(icell,l)     = headcell(icell)     !! The next lower ring index in that cell(icell)
			listbead(l,i)         = headbead(icell,l)   !! The next lower bead index of a part of a ring(l) in a cell(icell) 
			headbead(icell,l)     = i	            !! Highest bead index in the part of the ring(l) in a cell(icell)
			headcell_dummy(icell) = l

		end do	

		headcell = headcell_dummy                   !! Highest ring index in a cell(icell) 

	end do


	return
	end
	
	
 	!!*** Subroutine for the forces of interaction ***!!
	subroutine interaction(j1,frepx,frepy,dx,dy,r,icell,jcell,l)

	use parameters
        use position_arrays
	use forces
	use map_arrays
	use list_arrays
       
       
        implicit none
        integer:: i,j,j1,l,q
        double precision:: r,frepx,frepy,dx,dy,fadhx,fadhy,rlist,ti,tf
	integer:: icell,jcell,jcell0,nabor 
	
        	        f_intx=0.0d0
        	        f_inty=0.0d0
			frpx=0.0d0
			frpy=0.0d0
			fadx=0.0d0
			fady=0.0d0
              
        
        !! Loop Over All Cells !!

	do icell=1,ncell
			!write(*,*)'icell=',icell,'ncell=',ncell

			l=headcell(icell)  !! Highest ring index in a cell(icell)
			
50			if(l.ne.0) then
				!write(*,*)'icell=',icell,'l(headcell(icell)=',l

				i=headbead(icell,l) !! Highest bead index in the part of the ring(l) in a cell(icell)

100				if(i.ne.0) then
					!write(14,*) 'icell=',icell,'l(headcell(icell)=',l,'i(bead(icell,l)=',i

					q=listcell(icell,l)  !! The next lower ring index after l-th ring in a cell(icell) 

200					if((q.ne.0)) then !.and.(q.ne.l)) then
					       
						j=headbead(icell,q) !! The highest bead index in the q-th ring in a cell(icell) 

300						if(j.ne.0) then
											
							dx = x(q,j)-x(l,i)
							dy = y(q,j)-y(l,i)
							dx = dx - box*nint(dx/box)
							dy = dy - box*nint(dy/box)					        
							r = dsqrt(dx*dx + dy*dy)
							
							frepx=0.0d0
							frepy=0.0d0
							fadhx=0.0d0
							fadhy=0.0d0

                      					if(r.lt.rc_rep) then

				          			frepx = -k_rep*(rc_rep-r)*(dx)/r
				          			!fadhx = k_adh*(rc_adh-r)*(dx)/r
					  			frepy = -k_rep*(rc_rep-r)*(dy)/r
				          			!fadhy = k_adh*(rc_adh-r)*(dy)/r

				          			f_intx(l,i) = f_intx(l,i) + frepx !+ fadhx
				          			f_intx(q,j) = f_intx(q,j) - frepx !- fadhx

				          			f_inty(l,i) = f_inty(l,i) + frepy !+ fadhy 
				          			f_inty(q,j) = f_inty(q,j) - frepy !- fadhy

							! 	frpx(l,i) = frpx(l,i) + frepx
							!	frpx(q,j) = frpx(q,j) - frepx
							!	frpy(l,i) = frpy(l,i) + frepy
							!	frpy(q,j) = frpy(q,j) - frepy

							!	fadx(l,i) = fadx(l,i) + fadhx
							!	fadx(q,j) = fadx(q,j) - fadhx
							!	fady(l,i) = fady(l,i) + fadhy
							!	fady(q,j) = fady(q,j) - fadhy


                        				else if((r.le.rc_adh).and.(r.ge.rc_rep)) then

								fadhx = k_adh*(rc_adh-r)*(dx)/r
								fadhy = k_adh*(rc_adh-r)*(dy)/r

								f_intx(l,i) = f_intx(l,i) + fadhx
								f_intx(q,j) = f_intx(q,j) - fadhx

						                f_inty(l,i) = f_inty(l,i) + fadhy 
								f_inty(q,j) = f_inty(q,j) - fadhy

								
								!fadx(l,i) = fadx(l,i) + fadhx
								!fadx(q,j) = fadx(q,j) - fadhx
								!fady(l,i) = fady(l,i) + fadhy
								!fady(q,j) = fady(q,j) - fadhy

       							end if

							j=listbead(q,j) !! the next lower bead index in the q-th ring for the same cell(icell)
							goto 300
						end if

						q=listcell(icell,q) !! the next lower ring index in the same cell(icell)
						goto 200
					end if
	                      !! Loop Over Neighbouring Cells !!

				!	if(ncell.lt.6) goto 700

					jcell0 = 4*(icell-1)     
					
					do nabor=1,4
						jcell = map(jcell0 + nabor)  !!

	!! Loop Over All Rings & Beads Of The Neighbouring Cells!! 

						q=headcell(jcell)
																	
400						if(q.ne.0)then
						    
                                                     if(q.eq.l) goto 600
	
							j=headbead(jcell,q)

500							if(j.ne.0) then

							dx = x(q,j)-x(l,i)
							dy = y(q,j)-y(l,i)
							dx = dx - box*nint(dx/box)
							dy = dy - box*nint(dy/box)					        
							r = dsqrt(dx*dx + dy*dy)
							
							frepx=0.0d0
							frepy=0.0d0
							fadhx=0.0d0
							fadhy=0.0d0

                      					if(r.lt.rc_rep) then

				          			frepx = -k_rep*(rc_rep-r)*(dx)/r
				          			!fadhx = k_adh*(rc_adh-r)*(dx)/r
					  			frepy = -k_rep*(rc_rep-r)*(dy)/r
				          			!fadhy = k_adh*(rc_adh-r)*(dy)/r

				          			f_intx(l,i) = f_intx(l,i) + frepx !+ fadhx
				          			f_intx(q,j) = f_intx(q,j) - frepx !- fadhx

				          			f_inty(l,i) = f_inty(l,i) + frepy !+ fadhy 
				          			f_inty(q,j) = f_inty(q,j) - frepy !- fadhy

							! 	frpx(l,i) = frpx(l,i) + frepx
							!	frpx(q,j) = frpx(q,j) - frepx
							!	frpy(l,i) = frpy(l,i) + frepy
							!	frpy(q,j) = frpy(q,j) - frepy

							!	fadx(l,i) = fadx(l,i) + fadhx
							!	fadx(q,j) = fadx(q,j) - fadhx
							!	fady(l,i) = fady(l,i) + fadhy
							!	fady(q,j) = fady(q,j) - fadhy


                        				else if((r.le.rc_adh).and.(r.ge.rc_rep)) then

								fadhx = k_adh*(rc_adh-r)*(dx)/r
								fadhy = k_adh*(rc_adh-r)*(dy)/r

								f_intx(l,i) = f_intx(l,i) + fadhx
								f_intx(q,j) = f_intx(q,j) - fadhx

						                f_inty(l,i) = f_inty(l,i) + fadhy 
								f_inty(q,j) = f_inty(q,j) - fadhy

								
								!fadx(l,i) = fadx(l,i) + fadhx
								!fadx(q,j) = fadx(q,j) - fadhx
								!fady(l,i) = fady(l,i) + fadhy
								!fady(q,j) = fady(q,j) - fadhy

       							end if


								j=listbead(q,j)!! Considering the next bead of the current ring in the neighbour cell
								goto 500
							end if

600							q=listcell(jcell,q) !!Considering the next ring in the neighbour cell


							goto 400
						end if
					end do


700					i=listbead(l,i)  !! Considering the next bead of the current ring in the current cell
					goto 100
				end if

				l=listcell(icell,l)  !!Considering the next ring in the current cell  	
				goto 50
			end if
	end do		
					
	return
        end                

	

!!! Subroutine for random initial configurations
subroutine initial(r)

use position_arrays
 
implicit none
double precision:: a,b,a1,b1,g,r,dr,q,pi,theta1,ran2
integer:: l,k1,i,t,idum

q = 0.5d0*dsqrt(m*7.3d0)
dr = 0.2d0
t = 0

CALL SYSTEM_CLOCK(COUNT=idum)

!call srand(653214)

4   a = q*2.0d0*ran2(idum)
    b = q*2.0d0*ran2(idum)
!4   a = q*2.0d0*rand(0)
!    b = q*2.0d0*rand(0)

   if((a.lt.(r+0.05d0)).or.(b.lt.(r+0.05d0))) goto 4

   x1(1) = a
   y1(1) = b


do l=2,m

5   a1 = q*2.0d0*ran2(idum)
    b1 = q*2.0d0*ran2(idum)
!5   a1 = q*2.0d0*rand(0)
!    b1 = q*2.0d0*rand(0)

    if((a1.lt.(r+0.05d0)).or.(b1.lt.(r+0.05d0))) goto 5

   x1(l) = a1
   y1(l) = b1

    do k1=1,l-1

       g = dsqrt((x1(l)-x1(k1))**2 +(y1(l)-y1(k1))**2)

       if (g.lt.(2*r+dr)) then
      !  t=t+1
      !  write(*,*)t
        goto 5
       end if
    end do

end do

       pi = 4.0d0*datan(1.0d0)
       theta1 = 0.0d0
       
  
     do l=1,m  
       do i=1,n
        
           x(l,i) = r*cos(theta1) + 0.01d0*(2.0d0*ran2(idum)-1.0d0) + x1(l)
           y(l,i) = r*sin(theta1) + 0.01d0*(2.0d0*ran2(idum)-1.0d0) + y1(l)
           !x(l,i) = r*cos(theta1) + x1(l) !+ 0.01d0*(2.0d0*rand(0)-1.0d0) 
           !y(l,i) = r*sin(theta1) + y1(l) !+ 0.01d0*(2.0d0*rand(0)-1.0d0)  
           
           theta1 = theta1 + 2.0d0*pi/n
                                                         
       end do
     end do

	return
	end     


!! Subroutine for force calculation in a single cell 
	subroutine force(j1)

	use parameters	
	use position_arrays
	use forces

	implicit none
	integer:: i,j,l,j1
        double precision::l1,l2,dx1,dx2,dy1,dy2
   
        ! boundary condition
	
	do l=1,m
          x(l,0) = x(l,n)
          y(l,0) = y(l,n)
          x(l,n+1) = x(l,1)
          y(l,n+1) = y(l,1)
        end do
          

         
     
  	  do l=1,m             ! loop for cell no.
              do i=1,n         ! loop for beads in each cell

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
		 

	return
	end


!! Subroutine for the movement step in time using Euler-Maruyama algo ( including noise)
       subroutine move_noise(dt)

       use position_arrays
       use forces

       implicit none
       integer:: i,j,l,idum,j1
       double precision:: c,dt,noise,tau
       double precision:: vx,vy,wz
       double precision :: theta_x, theta_y, theta_sq_by_4 ! as in Dey arxiv: 1811.06450
          
        c = 0.5d0       ! c is coeff. of viscous damping      
 
       do l=1,m


         do i=1,n

	     
           vx = (fx(l,i) + f_intx(l,i))*dt/c + Vo*mx(l,i)*dt
           vy = (fy(l,i) + f_inty(l,i))*dt/c + Vo*my(l,i)*dt       
           x(l,i) = x(l,i) + vx
           y(l,i) = y(l,i) + vy
            
            
            CALL SYSTEM_CLOCK(COUNT=idum)

            CALL gasdev(noise,idum,mean,var)

            wz = (mx(l,i)*vy - my(l,i)*vx)/tau_align + noise
            theta_x = -my(l,i)*wz*dt
            theta_y = mx(l,i)*wz*dt
            theta_sq_by_4 = (theta_x*theta_x + theta_y*theta_y)/4.0d0
            
            ! Norm preserving rotation of m with ang vel w -> ang dispacement wz*dt
            mx(l,i) = ((1.0d0 - theta_sq_by_4)*mx(l,i) + theta_x)/(1.0d0 + theta_sq_by_4)
            my(l,i) = ((1.0d0 - theta_sq_by_4)*my(l,i) + theta_y)/(1.0d0 + theta_sq_by_4)
            
            

         end do

      end do

      return
      end


!! Subroutine for the movement step in time using Euler algo (No noise)
      subroutine move_deterministic(dt)

       use position_arrays
       use forces

       implicit none
       integer:: i,j,l,idum,j1
       double precision:: c,dt  
       
       c = 0.5d0       ! c is coeff. of viscous damping      
      
       do l=1,m


         do i=1,n

            x(l,i) = x(l,i) + (fx(l,i) + f_intx(l,i))*dt/c
            y(l,i) = y(l,i) + (fy(l,i) + f_inty(l,i))*dt/c
      

         end do

      end do

      return
      end


!! Subroutine for initializaion of noise-angle
       subroutine initial_angle           

        use position_arrays

        implicit none

        double precision :: pi,ran2,m_tmp,mx_tmp,my_tmp
        integer :: idum,i,l

         pi = 4.0d0*datan(1.0d0)
       !  idum = 564326
        CALL SYSTEM_CLOCK(COUNT=idum)
         
        do l=1,m
            do i=1,n
                mx_tmp = 2.0d0*ran2(idum) - 1.0d0
                my_tmp = 2.0d0*ran2(idum) - 1.0d0
                m_tmp = dsqrt(mx_tmp*mx_tmp + my_tmp*my_tmp)
                mx(l,i)=mx_tmp/m_tmp
                my(l,i)=my_tmp/m_tmp     
            end do
        end do

        return
        end

 
!!!!!!!!/////// Gaussian random no. generator \\\\\\\\\\\\\\\\\\\\\\\\\\\

 Subroutine gasdev(g2,idum,mean,variance) 
  !use numz
  Implicit none
      INTEGER,INTENT(INOUT)::idum
      DOUBLE PRECISION,INTENT(IN)::variance,mean
      DOUBLE PRECISION,INTENT(OUT)::g2
      DOUBLE PRECISION::ran2,g1
      INTEGER:: iset
      DOUBLE PRECISION:: fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
    !  if (iset.eq.0) then
        DO
        v1=2.0d0*ran2(idum)-1.0d0
        v2=2.0d0*ran2(idum)-1.0d0
        rsq=v1**2+v2**2
        if((rsq<1.0d0).AND.(rsq/=0.0d0))EXIT
        ENDDO
        fac=variance*SQRT(-2.0d0*log(rsq)/(rsq))
        g1=v1*fac+mean
        g2=v2*fac+mean      
      END subroutine gasdev


!!!!/////// Uniform Random number generators////////////////////////////////////

FUNCTION ran2(idum)
 ! USE numz
  IMPLICIT NONE
  DOUBLE PRECISION:: ran2
  !INTEGER,INTENT(inout),OPTIONAL::idum
  INTEGER,INTENT(inout)::idum
  !INTEGER :: idum
  INTEGER,PARAMETER::IM1=2147483563,IM2=2147483399,IMM1=IM1-1
  INTEGER,PARAMETER::IA1=40014,IA2=40692,IQ1=53668
  INTEGER,PARAMETER::IQ2=52774,IR1=12211,IR2=3791   
  INTEGER,PARAMETER::NTAB=32,NDIV=1+IMM1/NTAB
  DOUBLE PRECISION,PARAMETER::AM=1.0d0/IM1,EPS=1.2e-7,RNMX=1.0d0-EPS
  INTEGER::idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  IF (idum<0) THEN
     idum=MAX(-idum,1)
     idum2=idum
      DO j=NTAB+8,1,-1
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         IF (idum<0) idum=idum+IM1
         IF (j.LE.NTAB) iv(j)=idum
      ENDDO
      iy=iv(1)
   ENDIF
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   IF (idum<0) idum=idum+IM1
   k=idum2/IQ2
   idum2=IA2*(idum2-k*IQ2)-k*IR2
   IF (idum2<0) idum2=idum2+IM2
   j=1+iy/NDIV
   iy=iv(j)-idum2
   iv(j)=idum
   IF(iy.LT.1)iy=iy+IMM1
   ran2=MIN(AM*iy,RNMX)
   RETURN
 END FUNCTION ran2

