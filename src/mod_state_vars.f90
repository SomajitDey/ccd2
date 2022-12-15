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
    double precision :: timepoint = 0.0d0 ! Time instant (#steps x dt)

    contains

    integer function allocate_array_stat(cellnum, beadnum)
        use ring_nb, only: allocate_ring_nb
        integer, intent(in) :: cellnum, beadnum
        integer :: seeds_size

        call random_seed(size = seeds_size)
        
        allocate( &
            x(1:cellnum, 1:beadnum), y(1:cellnum, 1:beadnum), &
            mx(1:cellnum, 1:beadnum), my(1:cellnum, 1:beadnum), &
            fx(1:cellnum, 1:beadnum), fy(1:cellnum, 1:beadnum), &
            f_rpx(1:cellnum, 1:beadnum), f_rpy(1:cellnum, 1:beadnum), &
            f_adx(1:cellnum, 1:beadnum), f_ady(1:cellnum, 1:beadnum), &
            noise(1:cellnum*beadnum), &
            prng_seeds(seeds_size), &
            stat=allocate_array_stat )
        
        call allocate_ring_nb(cellnum)
    end function allocate_array_stat
end module state_vars
