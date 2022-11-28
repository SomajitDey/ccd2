module state_vars
    implicit none
    public
    double precision, dimension(:,:), allocatable :: x,y ! Position coordinates of every bead
    double precision, dimension(:,:), allocatable :: mx,my ! Motility unit vector for every bead
    double precision, dimension(:,:), allocatable :: fx,fy ! Intracellular force
    double precision, dimension(:,:), allocatable :: f_rpx,f_rpy ! Intercellular steric repulsion / volume exclusion
    double precision, dimension(:,:), allocatable :: f_adx,f_ady ! Intercellular attraction / adhesion
    integer, dimension(:), allocatable :: prng_seeds ! Stores the state of the P(seudo) R(andom) N(um) G(enerator)

    contains

    integer function allocate_array_stat(cellnum, beadnum)
        integer, intent(in) :: cellnum, beadnum
        integer :: seeds_size

        call random_seed(size = seeds_size)
        
        allocate( &
            x(1:cellnum, 0:beadnum+1), y(1:cellnum, 0:beadnum+1), &
            mx(1:cellnum, 1:beadnum), my(1:cellnum, 1:beadnum), &
            fx(1:cellnum, 1:beadnum), fy(1:cellnum, 1:beadnum), &
            f_rpx(1:cellnum, 1:beadnum), f_rpy(1:cellnum, 1:beadnum), &
            f_adx(1:cellnum, 1:beadnum), f_ady(1:cellnum, 1:beadnum), &
            prng_seeds(seeds_size), &
            stat=allocate_array_stat )
        ! 0:beadnum+1 is to expand a circular array into linear one later in force
    end function allocate_array_stat
end module state_vars
