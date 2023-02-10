module grid_linked_list
       use shared

       implicit none

	   integer, protected:: w, ncell
       double precision, protected :: celli
	   integer, protected:: gridmapsiz
	   integer, dimension(:), allocatable:: gridmap
       integer, dimension(:), allocatable :: bead_nl_head, bead_nl_body ! bead-only neighborlist: head and body
       private :: grid_index
       
    contains
    
    !! 2D to 1D representation of the grid points !!
    ! Assuming that within the central simulation box, ix and iy are within [1,w] 
    pure integer function grid_index(ix, iy)
        integer, intent(in):: ix,iy

      	grid_index = 1 + modulo(ix-1, w) + modulo(iy-1, w)*w
    end function grid_index

	!!*** Subroutine to set up the cells and to make the list of neighbouring cells through PBC ***!!
	subroutine gridmaps()
	integer:: ix,iy,igridmap,alloc_stat
    double precision:: rcut

    rcut = max(rc_rep, rc_adh)  ! Grid-size estimate
    ! Note: l0 isn't in max() because the intra-force computation where l0 is required doesnt use the neighborlist
    w=floor(box/rcut)
	celli = dble(w/box) ! celli is the inverse of cell length
	if((1.d0/celli).lt.rcut) error stop 'Fatal: Grid size smaller than interaction cutoff'
    ncell=w*w
    gridmapsiz=4*ncell
    
    allocate(gridmap(gridmapsiz), bead_nl_head(ncell), bead_nl_body(m*n), stat=alloc_stat)
        
     if(alloc_stat /= 0) error stop 'Fatal: Problem while allocating grid_linked_list'
    
	!! Find Half The Nearest Neighbours Of Each Cell !!

	do iy=1,w
        	do ix=1,w

			
			igridmap = (grid_index(ix,iy) - 1)*4
		
			gridmap(igridmap+1) = grid_index(ix+1,iy)
			gridmap(igridmap+2) = grid_index(ix+1,iy+1)
			gridmap(igridmap+3) = grid_index(ix,iy+1)
			gridmap(igridmap+4) = grid_index(ix-1,iy+1)
		end do
	end do
	end subroutine gridmaps


	!!*** Subroutine to make linked lists & head of chain arrays ***!!
	subroutine links()
	integer:: icell,l,i, bead_index
	
   	!! Zero Head Of Chain Array & List Array !!
    bead_nl_head = 0
    bead_nl_body = 0

	!! Sort All Beads !!

    do l=1,m	 !! Loop over rings
        do i=1,n !! Loop over beads

            ! determining the grid index for a particular bead
			icell = grid_index(ceiling(x(i,l)*celli), ceiling(y(i,l)*celli))
			
			bead_index = (l-1)*n+i ! Global serial of the bead : [1,mn]

			bead_nl_body(bead_index) = bead_nl_head(icell)
			bead_nl_head(icell) = bead_index
		end do	
	end do
	end subroutine links

end module grid_linked_list
