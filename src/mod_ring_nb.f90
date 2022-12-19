!Brief: Bookkeeping for ring-ring neighborhood
module ring_nb
    implicit none

    integer, dimension(:), allocatable, protected :: ring_nb_yesno_packed ! holds ring ring neighborhood info
    ! the above is a 1D packing of the strictly upper triangular matrix - whether_neighbors(i,j) for i<j, j=1,m
    private :: index_packed_strictly_upper_triang_mat
    integer, dimension(:), allocatable :: ring_nb_io
    ! compressed ring_nb info : stores only non-zero elements of ring_nb_yesno_packed

contains

    subroutine allocate_ring_nb(m)
        integer, intent(in) :: m
        integer :: alloc_stat

        allocate(ring_nb_yesno_packed((m-1)*(m-2)/2 + (m-1)), ring_nb_io(5*m), stat=alloc_stat)
        ! Assuming a max of 10 neighbors per ring. Avoiding double counting : 10*m/2 = 5*m elements for ring_nb_io
        if(alloc_stat /= 0) error stop 'Problem while allocating ring_nb'
    end subroutine allocate_ring_nb
    
    subroutine init_ring_nb()
        integer :: i
        !$omp do private(i)
        do i=1,size(ring_nb_yesno_packed)
            ring_nb_yesno_packed(i) = 0
        end do
        !$omp end do
    end subroutine init_ring_nb
    
    pure integer function index_packed_strictly_upper_triang_mat(l, q) result(array_1d_index)
        integer, intent(in) :: l, q
        integer :: row, col
        
        row = min(l,q)
        col = max(l,q)
        array_1d_index = (col-1)*(col-2)/2 + row
    end function index_packed_strictly_upper_triang_mat

    subroutine assert_are_nb_rings(l, q)
        integer, intent(in) :: l, q
        integer :: i
        
        i = index_packed_strictly_upper_triang_mat(l, q)
        ring_nb_yesno_packed(i) = i
    end subroutine assert_are_nb_rings
    
    pure logical function are_nb_rings(l, q) result(check)
        integer, intent(in) :: l, q

        check = ( ring_nb_yesno_packed( index_packed_strictly_upper_triang_mat(l, q) ) /= 0 )
    end function are_nb_rings
    
    subroutine pack_ring_nb()
        integer :: i
        
        if(size(ring_nb_io) < count(ring_nb_yesno_packed /= 0)) error stop &
            'Too small output buffer for holding the compressed ring-ring neighborhood information'
        
        ring_nb_io = pack( array=ring_nb_yesno_packed, mask=ring_nb_yesno_packed/=0, &
                        vector=[ (0, i=1,size(ring_nb_io)) ] )
    end subroutine pack_ring_nb
    
    subroutine unpack_ring_nb()
        integer :: i, valu
        
        call init_ring_nb()
        i=1
        do
            valu = ring_nb_io(i)
            if(valu == 0) exit
            ring_nb_yesno_packed(valu) = valu
            i = i+1
        end do
    end subroutine unpack_ring_nb

end module ring_nb
