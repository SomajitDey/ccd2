program ccd_cpt_to_xy
    use files
    implicit none
    integer :: pending_steps
    character(len=40) :: params_hash
    double precision :: boxlen
    
    call cpt_read(timepoint, recnum, pending_steps, params_hash, boxlen)
    call xy_dump('config.xy', boxlen)
end program ccd_cpt_to_xy
