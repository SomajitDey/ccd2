module integrator

    implicit none
    double precision :: evolve_motility_bool = 1.0d0

contains

!! Subroutine for the movement step in time using Euler-Maruyama algo ( including noise)
    subroutine move_noise

        use shared

        integer :: i, l
        double precision :: vx, vy, vx_sum, vy_sum, vnorm, wz, mx_cell, my_cell
        double precision :: theta_x, theta_y, theta_sq_by_4 ! as in Dey arxiv: 1811.06450
        double precision :: factor1, factor2

!$omp do private(i,l,vx,vy,vx_sum,vy_sum,vnorm,wz,mx_cell,my_cell,theta_x,theta_y,theta_sq_by_4,factor1,factor2)
        do l = 1, m
            vx_sum = 0.d0
            vy_sum = 0.d0

            do i = 1, n
                vx = (fx(i, l) + f_adx(i, l) + f_rpx(i, l))/c + Vo*mx(i, l)
                vy = (fy(i, l) + f_ady(i, l) + f_rpy(i, l))/c + Vo*my(i, l)
                vx_sum = vx + vx_sum
                vy_sum = vy + vy_sum
                x(i, l) = x(i, l) + vx*dt
                y(i, l) = y(i, l) + vy*dt
            end do
            vnorm = hypot(vx_sum, vy_sum)
            if (vnorm < epsilon(0.0d0)) then
                ! Null vector : protection against division by vnorm=0
                vx_sum = 0.d0
                vy_sum = 0.d0
                vnorm = 1.d0
            end if
            mx_cell = mx(1, l)
            my_cell = my(1, l)
            wz = (align_strength*(mx_cell*vy_sum - my_cell*vx_sum)/vnorm + noise_strength*noise(l)) &
                 *evolve_motility_bool
            theta_x = -my_cell*wz*dt
            theta_y = mx_cell*wz*dt
            theta_sq_by_4 = (theta_x*theta_x + theta_y*theta_y)*0.25d0

            ! Norm preserving rotation of m with ang vel w -> ang dispacement wz*dt
            factor1 = 1.0d0 - theta_sq_by_4
            factor2 = 1.0d0/(1.0d0 + theta_sq_by_4)

            mx_cell = (factor1*mx_cell + theta_x)*factor2
            my_cell = (factor1*my_cell + theta_y)*factor2

            mx(:, l) = mx_cell
            my(:, l) = my_cell
        end do
!$omp end do

        return
    end subroutine move_noise

end module integrator
