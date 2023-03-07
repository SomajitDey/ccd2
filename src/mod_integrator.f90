module integrator

implicit none
double precision :: evolve_motility_bool=1.0d0

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
           vx = (fx(i,l) + f_adx(i,l) + f_rpx(i,l))/c + Vo*mx(i,l)
           vy = (fy(i,l) + f_ady(i,l) + f_rpy(i,l))/c + Vo*my(i,l)
           x(i,l) = x(i,l) + vx*dt
           y(i,l) = y(i,l) + vy*dt
            
      
            wz = ((mx(i,l)*vy - my(i,l)*vx)/(tau_align*dt) + noise(tmp+i))*evolve_motility_bool
            theta_x = -my(i,l)*wz*dt
            theta_y = mx(i,l)*wz*dt
            theta_sq_by_4 = (theta_x*theta_x + theta_y*theta_y)/4.0d0
            
            ! Norm preserving rotation of m with ang vel w -> ang dispacement wz*dt
            factor1 = 1.0d0 - theta_sq_by_4
            factor2 = 1.0d0 + theta_sq_by_4
            
            mx(i,l) = (factor1*mx(i,l) + theta_x)/factor2
            my(i,l) = (factor1*my(i,l) + theta_y)/factor2
        end do
      end do
      !$omp end do

      return
      end subroutine move_noise

end module integrator
