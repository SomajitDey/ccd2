module init
use shared
implicit none

contains

!!! Subroutine for random initial configurations
subroutine initial
double precision:: radius   ! Initial/seed cell radius
double precision:: acell,bcell,dr,mindist,dx,dy, theta1
integer:: l,kcell,i,fail_count
real:: rands(2)
double precision, dimension(size(x,dim=2)) :: xcell, ycell ! centre coordinates of any circular cell/ring

radius = 0.5d0*l0/dsin(pi/n) ! circumcirle of a regular n-gon with side l0
dr = rc_rep ! initial separation between cell peripheries
mindist = 2*radius+dr
box = dsqrt(m*pi*mindist*mindist/4) ! box should accommodate m circumcircles with dr separation

! Seeding the first cell centre at origin
xcell(1) = 0.d0
ycell(1) = 0.d0

! Find a suitable centre for each circular cell by random trial and error
seed_cell_centres: do l=2,m
    fail_count=0

    trial: do
        fail_count = fail_count+1 
        if(fail_count > 1000) then
            ! Too many failed seed_cell_centress means box is too small to accomodate yet another cell
            ! Let's increase the box size
            box = box*1.1d0
            fail_count = 0
        end if 
        
        CALL RANDOM_NUMBER(rands)
        acell = box*rands(1)
        bcell = box*rands(2)

        distance_from_other_cells: do kcell=1,l-1
            dx = acell-xcell(kcell)
            dy = bcell-ycell(kcell)
            dx = dx - box*nint(dx/box)
            dy = dy - box*nint(dy/box)
            if ((dx*dx + dy*dy).lt.(mindist*mindist)) cycle trial
        end do distance_from_other_cells

        xcell(l) = acell
        ycell(l) = bcell

        exit trial
    end do trial

end do seed_cell_centres

! Construct circular cells from the seeded centres
! Also initialize the motility unit vectors symmetrically (radially outward) such that the total for each cell is null
do l=1,m
    do i=1,n        
        theta1 = (i-1)*2.0d0*pi/n
        x(i,l) = radius*dcos(theta1) + xcell(l)
        y(i,l) = radius*dsin(theta1) + ycell(l)
        mx(i,l)= dcos(theta1)
        my(i,l)= dsin(theta1)
     end do
end do

	return
	end subroutine initial

!! Subroutine for random initializaion of motility unit vectors
       subroutine initial_angle()
        double precision :: m_tmp,mx_tmp,my_tmp
        integer :: i,l
        DOUBLE PRECISION:: rands(2)
        
        do l=1,m
            do i=1,n
               CALL RANDOM_NUMBER(rands)

                mx_tmp = 2.0d0*rands(1) - 1.0d0
                my_tmp = 2.0d0*rands(2) - 1.0d0
                m_tmp = dsqrt(mx_tmp*mx_tmp + my_tmp*my_tmp)
                mx(i,l)=mx_tmp/m_tmp
                my(i,l)=my_tmp/m_tmp
            end do
        end do

        return
        end subroutine initial_angle


end module init
