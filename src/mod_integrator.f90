module integrator

implicit none

contains

!! Subroutine for the movement step in time using Euler-Maruyama algo ( including noise)
       subroutine move_noise

       use shared
       
       integer:: i,j,l
       double precision:: vx,vy,wz
       double precision :: theta_x, theta_y, theta_sq_by_4 ! as in Dey arxiv: 1811.06450
       double precision :: noise(m,n), g(m*n)
          
 
             CALL gasdev(g,mean,var)
             noise=reshape(g, [m,n]) ! reshapes 1D array g into 2D

       do concurrent (l=1:m, i=1:n)
           vx = (fx(l,i) + f_adx(l,i) + f_rpx(l,i))*dt/c + Vo*mx(l,i)*dt
           vy = (fy(l,i) + f_ady(l,i) + f_rpy(l,i))*dt/c + Vo*my(l,i)*dt
           x(l,i) = x(l,i) + vx
           y(l,i) = y(l,i) + vy
            
      
            wz = (mx(l,i)*vy - my(l,i)*vx)/(tau_align*dt) + noise(l,i)
            theta_x = -my(l,i)*wz*dt
            theta_y = mx(l,i)*wz*dt
            theta_sq_by_4 = (theta_x*theta_x + theta_y*theta_y)/4.0d0
            
            ! Norm preserving rotation of m with ang vel w -> ang dispacement wz*dt
            mx(l,i) = ((1.0d0 - theta_sq_by_4)*mx(l,i) + theta_x)/(1.0d0 + theta_sq_by_4)
            my(l,i) = ((1.0d0 - theta_sq_by_4)*my(l,i) + theta_y)/(1.0d0 + theta_sq_by_4)
      end do

      return
      end subroutine move_noise

end module integrator
