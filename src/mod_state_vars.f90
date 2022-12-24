module state_vars
    use ring_nb, only: ring_nb_io, ring_nb_yesno_packed
    implicit none
    public
    double precision :: box
    double precision, dimension(:,:), allocatable :: x,y ! Position coordinates of every bead
    double precision, dimension(:,:), allocatable :: mx,my ! Motility unit vector for every bead
    double precision, dimension(:,:), allocatable :: fx,fy ! Intracellular force
    double precision, dimension(:,:), allocatable :: f_rpx,f_rpy ! Intercellular steric repulsion / volume exclusion
    double precision, dimension(:,:), allocatable :: f_adx,f_ady ! Intercellular attraction / adhesion
    double precision, dimension(:), allocatable :: noise
    integer, dimension(:), allocatable :: prng_seeds ! Stores the state of the P(seudo) R(andom) N(um) G(enerator)
    integer :: recnum=1 ! Record number the trajectory file is currently positioned at
    real :: timepoint = 0.0 ! Time instant (#steps x dt). No need to store in double precision

    contains

    integer function allocate_array_stat(cellnum, beadnum)
        use ring_nb, only: allocate_ring_nb
        integer, intent(in) :: cellnum, beadnum
        integer :: seeds_size

        call random_seed(size = seeds_size)
        
        allocate( &
            x(1:beadnum, 1:cellnum), y(1:beadnum, 1:cellnum), &
            mx(1:beadnum, 1:cellnum), my(1:beadnum, 1:cellnum), &
            fx(1:beadnum, 1:cellnum), fy(1:beadnum, 1:cellnum), &
            f_rpx(1:beadnum, 1:cellnum), f_rpy(1:beadnum, 1:cellnum), &
            f_adx(1:beadnum, 1:cellnum), f_ady(1:beadnum, 1:cellnum), &
            noise(1:cellnum*beadnum), &
            prng_seeds(seeds_size), &
            stat=allocate_array_stat )
        
        call allocate_ring_nb(cellnum)
    end function allocate_array_stat
end module state_vars
