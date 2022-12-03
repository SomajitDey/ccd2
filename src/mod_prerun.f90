module prerun
    implicit none
    
    contains
    
    subroutine prerun_setup(jf)
        use state_vars
        use parameters
        use files
        
        integer, intent(out) :: jf
        integer :: pending_steps
        character(len=40) :: params_hash
        character(len=len('replace')) :: traj_status
        ! Holds either 'old' or 'replace', same as status= specifier in an open statement
        logical :: append_flag_present, finish_prev_run

        if(.not. acquire_lock(force = (cmd_line_flag('-f') .or. cmd_line_flag('--force')))) error stop &
            'Uh-oh...seems like another run is going on in the current working directory. I better stop than mess up'
        
        call log_this('Reading initial state from '//cpt_fname)
        call cpt_read(timepoint, recnum, pending_steps, params_hash)

        call log_this('Reading run parameters from '//params_fname)
        call assign_params(params_fname)
        if((size(x,1) /= m).or.(size(x,2) /= n+2)) error stop &
            'System size as read in from checkpoint: '//cpt_fname// &
                ' does not match that given in parameter file: '//params_fname

        append_flag_present = cmd_line_flag('-a') .or. cmd_line_flag('--append')
        finish_prev_run = (pending_steps /= 0) .and. (params_hash == sha1(params_fname))
        
        if(finish_prev_run .or. append_flag_present) then
            traj_status='old'
        else
            traj_status='replace'
            recnum = 1
            timepoint = 0.0d0
        end if
        if(.not. finish_prev_run) pending_steps = 0
        jf = nsamples*traj_dump_int + pending_steps

        write(status_fd,'(i0)') jf/status_dump_int
        flush(status_fd)
        
        call log_this('Opening trajectory file: '//traj_fname)
        call open_traj('write', traj_status)

        if(finish_prev_run) call log_this('Continuing with the incomplete previous run')
        if(append_flag_present) call log_this('Extending the existing trajectory @ '//traj_fname)
    end subroutine prerun_setup
end module prerun
