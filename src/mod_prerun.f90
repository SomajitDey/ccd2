module prerun
    implicit none
    
    contains
    
    subroutine prerun_setup(ji,jf)
        use state_vars
        use parameters
        use files
        !$ use omp_lib, only: omp_get_max_threads
        
        integer, intent(out) :: ji,jf
        integer :: pending_steps, current_step
        character(len=40) :: params_hash
        character(len=len('replace')) :: traj_status
        ! Holds either 'old' or 'replace', same as status= specifier in an open statement
        logical :: append_flag_present, finish_prev_run

        if(.not. acquire_lock(force = (cmd_line_flag('-f') .or. cmd_line_flag('--force')))) error stop &
            'Uh-oh...seems like another run is going on in the current working directory. I better stop than mess up'
        
        call log_this('Reading initial state from '//cpt_fname)
        call cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)

        call log_this('Reading run parameters from '//params_fname)
        call assign_params(params_fname)
        if((size(x,2) /= m).or.(size(x,1) /= n)) error stop &
            'System size as read in from checkpoint: '//cpt_fname// &
                ' does not match that given in parameter file: '//params_fname

        !TODO: Include the following tests in assign_params by calling check_params() in mod_parameters
        if(l0 .gt. rc_adh) call log_this('Warning: Spring length is bigger than adhesion cutoff.')
        if(l0 .gt. rc_rep) call log_this('Warning: Spring length is bigger than repulsion cutoff.')
        !TODO: Should the above be error stops instead of warnings ?
        if(rc_rep .gt. rc_adh) error stop 'Repulsion cutoff is bigger than adhesion cutoff.'
        if(rc_adh .gt. box/2) error stop &
            'Minimum image convention is at stake. Make box bigger than 2 x interaction-cutoff.'
        !TODO: The above should also check box/2 > max possible spring length
        ! estimate involves pressure pl0 and spring constant K - force balance FBD of regular n-gon

        append_flag_present = cmd_line_flag('-a') .or. cmd_line_flag('--append')
        finish_prev_run = (pending_steps /= 0) .and. (params_hash == sha1(params_fname))
        
        if(finish_prev_run .or. append_flag_present) then
            traj_status='old'
        else
            traj_status='replace'
            recnum = 0
            timepoint = 0.0
        end if

        if(.not. finish_prev_run) then
            pending_steps = nsamples*traj_dump_int - 1
            current_step = 1
        end if

        ji = current_step
        jf = ji + pending_steps
        if(append_flag_present .and. finish_prev_run) jf = jf + nsamples*traj_dump_int

        write(status_fd,'(i0, 1x, a)') jf/status_dump_int, 'new lines to follow'
        if(ji/status_dump_int == 1) then
            write(status_fd,*)
        else if(ji/status_dump_int > 1) then
            write(status_fd, '('//int_to_char(ji/status_dump_int - 1)//'/)')
        end if
        flush(status_fd)
        
        call log_this('Opening trajectory file: '//traj_fname)
        call open_traj('write', traj_status)

        if(finish_prev_run) call log_this('Continuing with the incomplete previous run')
        if(append_flag_present) call log_this('Extending the existing trajectory @ '//traj_fname)

        !$ call log_this('Using '//int_to_char(int(omp_get_max_threads(), kind=kind(jf)))//' OpenMP threads')
        
        do_status_dump = .not. (cmd_line_flag('-n') .or. cmd_line_flag('--no-status-dump'))
    end subroutine prerun_setup
end module prerun
