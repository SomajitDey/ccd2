program ccd_cpt_to_xy
    use files
    implicit none
    integer :: pending_steps, current_step
    character(len=40) :: params_hash
    
    call cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)
    call xy_dump('config.xy', box, x, y)
end program ccd_cpt_to_xy
