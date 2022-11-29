module init
implicit none

contains

!!! Subroutine for random initial configurations
subroutine initial

use shared

double precision:: a,b,a1,b1,g,dr,q,theta1
integer:: l,k1,i,t
      DOUBLE PRECISION:: rands(2)
      double precision, dimension(size(x,dim=1)) :: x1, y1


q = 0.5d0*dsqrt(m*7.3d0)
dr = 0.2d0
t = 0


    do
        CALL RANDOM_NUMBER(rands)
    a = q*2.0d0*rands(1)
    b = q*2.0d0*rands(2)
    if((a.ge.(radius+0.05d0)).and.(b.ge.(radius+0.05d0))) exit
    end do


   x1(1) = a
   y1(1) = b


do l=2,m
        lcell: do
    do
        CALL RANDOM_NUMBER(rands)
    a1 = q*2.0d0*rands(1)
    b1 = q*2.0d0*rands(2)
    if((a1.ge.(radius+0.05d0)).and.(b1.ge.(radius+0.05d0))) exit
    end do

   x1(l) = a1
   y1(l) = b1

    do k1=1,l-1

       g = dsqrt((x1(l)-x1(k1))**2 +(y1(l)-y1(k1))**2)

       if (g.lt.(2*radius+dr)) cycle lcell
    end do
        exit lcell
        end do lcell
end do

       theta1 = 0.0d0
       
  
     do l=1,m  
       do i=1,n
        CALL RANDOM_NUMBER(rands)
        
           x(l,i) = radius*cos(theta1) + 0.01d0*(2.0d0*rands(1)-1.0d0) + x1(l)
           y(l,i) = radius*sin(theta1) + 0.01d0*(2.0d0*rands(2)-1.0d0) + y1(l)
           theta1 = theta1 + 2.0d0*pi/n
                                                         
       end do
     end do

	return
	end subroutine initial

!! Subroutine for initializaion of noise-angle
       subroutine initial_angle           

        use shared

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
