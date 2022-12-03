module integrator

implicit none

contains

!! Subroutine for the movement step in time using Euler-Maruyama algo ( including noise)
       subroutine move_noise

       use shared
       
       integer:: i,l, tmp
       double precision:: vx,vy,wz
       double precision :: theta_x, theta_y, theta_sq_by_4 ! as in Dey arxiv: 1811.06450
       double precision :: factor1, factor2

        !$omp do private(i,l, vx,vy,wz, theta_x, theta_y, theta_sq_by_4, tmp)
       do l=1,m
        tmp = (l-1)*n
        do i=1,n
           vx = (fx(l,i) + f_adx(l,i) + f_rpx(l,i))*dt/c + Vo*mx(l,i)*dt
           vy = (fy(l,i) + f_ady(l,i) + f_rpy(l,i))*dt/c + Vo*my(l,i)*dt
           x(l,i) = x(l,i) + vx
           y(l,i) = y(l,i) + vy
            
      
            wz = (mx(l,i)*vy - my(l,i)*vx)/(tau_align*dt) + noise(tmp+i)
            theta_x = -my(l,i)*wz*dt
            theta_y = mx(l,i)*wz*dt
            theta_sq_by_4 = (theta_x*theta_x + theta_y*theta_y)/4.0d0
            
            ! Norm preserving rotation of m with ang vel w -> ang dispacement wz*dt
            factor1 = 1.0d0 - theta_sq_by_4
            factor2 = 1.0d0 + theta_sq_by_4
            
            mx(l,i) = (factor1*mx(l,i) + theta_x)/factor2
            my(l,i) = (factor1*my(l,i) + theta_y)/factor2
        end do
      end do
      !$omp end do nowait

      return
      end subroutine move_noise

end module integrator
