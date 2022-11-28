module grid_linked_list
       use shared

       implicit none

	   integer, protected:: w, ncell
	   integer, protected:: mapsiz
	   integer, dimension(:), allocatable:: map
	   integer, dimension(:), allocatable:: headcell
       integer, dimension(:,:), allocatable:: headbead,listcell,listbead,icell_memory
       private :: grid_index
       
    contains
    
    !! Function To Give Cell Index !!
    pure integer function grid_index(ix, iy)
        integer, intent(in):: ix,iy

      	grid_index = 1 + mod(ix-1+w , w) + mod(iy-1+w , w) * w
    end function grid_index

	!!*** Subroutine to set up the cells and to make the list of neighbouring cells through PBC ***!!
	subroutine maps()
	integer:: ix,iy,imap,alloc_stat

    w=int(box/rcut)
    ncell=w*w
    mapsiz=4*ncell
    
    allocate(map(mapsiz), headcell(ncell), headbead(ncell,m),listcell(ncell,m),listbead(m,n),icell_memory(m,n), &
        stat=alloc_stat)
        
     if(alloc_stat /= 0) error stop 'Problem while allocating grid_linked_list'
    
	!! Find Half The Nearest Neighbours Of Each Cell !!

	do iy=1,w
        	do ix=1,w

			
			imap = (grid_index(ix,iy) - 1) * 4
		
			map(imap+1) = grid_index(ix+1,iy)
			map(imap+2) = grid_index(ix+1,iy+1)
			map(imap+3) = grid_index(ix,iy+1)
			map(imap+4) = grid_index(ix-1,iy+1)
		end do
	end do
	end subroutine maps


	!!*** Subroutine to make linked lists & head of chain arrays ***!!
	subroutine links()
	integer:: icell,l,i
    integer, dimension(size(headcell)):: headcell_dummy
	double precision:: cell,celli,delta
    double precision, dimension(size(x,1),size(x,2)):: x3,y3
	
   	!! Zero Head Of Chain Array & List Array !!


	headcell      =0 !! head of chain array for the rings in a cell
	headcell_dummy=0
	listcell      =0 !! List array for the rings in a cell
	headbead      =0 !! head of chain array for the beads in a ring in the cell
	listbead      =0

	delta = 0.0d0

	celli = dble(w/box) !! celli is the inverse of cell length
	cell  = 1.0d0/celli !! cell is the cell length

	if(cell.lt.rcut) error stop 'Grid size too small compared to interaction cut-off'

	!! Sort All Beads !!

	do l=1,m	 !! Loop over rings
		do i=1,n !! Loop over beads

		        x3(l,i) = x(l,i) - box*floor(x(l,i)/box)
                        y3(l,i) = y(l,i) - box*floor(y(l,i)/box)


			icell = 1 + int(x3(l,i) * celli) + int(y3(l,i) * celli) * w  !! determining the cell index for a particular bead

			icell_memory(l,i) = icell
			
			
			listcell(icell,l)     = headcell(icell)     !! The next lower ring index in that cell(icell)
			listbead(l,i)         = headbead(icell,l)   !! The next lower bead index of a part of a ring(l) in a cell(icell) 
			headbead(icell,l)     = i	            !! Highest bead index in the part of the ring(l) in a cell(icell)
			headcell_dummy(icell) = l

		end do	

		headcell = headcell_dummy                   !! Highest ring index in a cell(icell) 

	end do
	end subroutine links

end module grid_linked_list
