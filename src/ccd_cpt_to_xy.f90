! Help:Begin
! Usage: ccd_cpt_to_xy
! Help:End

program ccd_cpt_to_xy
    use files
    use utilities, only: help_handler
    implicit none
    integer :: pending_steps, current_step
    character(len=40) :: params_hash

    call help_handler()

    call cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)
    call xy_dump('config.xy', box, x, y)
end program ccd_cpt_to_xy
