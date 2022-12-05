module init
use shared
implicit none

contains

!!! Subroutine for random initial configurations
subroutine initial
double precision:: radius   ! Initial/seed cell radius
double precision:: acell,bcell,dr,mindist,dx,dy, theta1
integer:: l,kcell,i,counter
DOUBLE PRECISION:: rands(2)
double precision, dimension(size(x,dim=1)) :: xcell, ycell

radius = 0.5d0*l0/dsin(pi/n) ! circumcirle of a regular n-gon with side l0
dr = rc_rep ! initial separation between cell peripheries
mindist = 2*radius+dr
!box = dsqrt(m*dsqrt(3.0d0)/24.0d0) * n*l0
! Derived from m*(area of hexagon with perimeter n*l0 or 2*pi*radius) = box^2
box = dsqrt(m*pi*mindist*mindist/4) ! box should accommodate m circumcircles with dr separation

x=0.d0 !DEBUG:
y=0.d0 !DEBUG:
xcell(1) = 0.d0
ycell(1) = 0.d0

trial: do l=2,m
    counter=0

    lcell: do
        counter = counter+1 
        if(counter > 1000) then
            !DEBUG: if(tradius .lt. radius/2) then 
                !DEBUG: print*, 'I am cell num: ', l 
                !DEBUG: exit lcell 
            !DEBUG: end if 
            box = box*1.1d0
            counter = 0
        end if 
        
        CALL RANDOM_NUMBER(rands)
        acell = box*rands(1)
        bcell = box*rands(2)

        do kcell=1,l-1
            dx = acell-xcell(kcell)
            dy = bcell-ycell(kcell)
            dx = dx - box*nint(dx/box)
            dy = dy - box*nint(dy/box)
            if ((dx*dx + dy*dy).lt.(mindist*mindist)) cycle lcell
        end do

        xcell(l) = acell
        ycell(l) = bcell

        exit lcell
    end do lcell

    do i=1,n        
        theta1 = (i-1)*2.0d0*pi/n
        x(l,i) = radius*dcos(theta1) + xcell(l)
        y(l,i) = radius*dsin(theta1) + ycell(l)
     end do
end do trial


	return
	end subroutine initial

!! Subroutine for initializaion of noise-angle
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
                mx(l,i)=mx_tmp/m_tmp
                my(l,i)=my_tmp/m_tmp     
            end do
        end do

        return
        end subroutine initial_angle


end module init
