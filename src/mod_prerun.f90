module prerun
    implicit none

contains

    subroutine prerun_setup(ji, jf)
        use state_vars
        use parameters
        use files
        use integrator, only: evolve_motility_bool
        use forces, only: force, force_pl, force_pl0
!$      use omp_lib, only: omp_get_max_threads

        integer, intent(out) :: ji, jf
        integer :: pending_steps, current_step
        character(len=40) :: params_hash
        character(len=len('replace')) :: traj_status
        ! Holds either 'old' or 'replace', same as status= specifier in an open statement
        logical :: append_flag_present, finish_prev_run

        if (.not. acquire_lock(force=(cmd_line_flag('-f') .or. cmd_line_flag('--force')))) error stop &
            'Fatal: Seems like another run is going on in the current working directory. Exiting...'

        call log_this('Reading initial state from '//cpt_fname)
        call cpt_read(timepoint, recnum, pending_steps, current_step, params_hash)

        call log_this('Reading run parameters from '//params_fname)
        call assign_params(params_fname, nocheck=cmd_line_flag('--no-check'))
        if ((size(x, 2) /= m) .or. (size(x, 1) /= n)) error stop &
            'Fatal: System size as read in from checkpoint: '//cpt_fname// &
            ' does not match that given in parameter file: '//params_fname
        if (rc_adh .gt. box/2) error stop &
            'Fatal: Minimum image convention is at stake. Make box bigger than 2 x interaction-cutoff.'
        !TODO: The above should also check box/2 > max possible spring length
        ! estimate involves pressure pl0 and spring constant K - force balance FBD of regular n-gon

        append_flag_present = cmd_line_flag('-a') .or. cmd_line_flag('--append')
        finish_prev_run = (pending_steps /= 0) .and. (params_hash == sha1(params_fname))

        if (finish_prev_run .or. append_flag_present) then
            traj_status = 'old'
        else
            traj_status = 'replace'
            recnum = 0
            timepoint = 0.0
        end if

        if (.not. finish_prev_run) then
            pending_steps = nsamples*traj_dump_int - 1
            current_step = 1
        end if

        ji = current_step
        jf = ji + pending_steps
        if (append_flag_present .and. finish_prev_run) jf = jf + nsamples*traj_dump_int

        write (status_fd, '(i0, 1x, a)') jf/status_dump_int, 'new lines to follow'
        if (ji/status_dump_int == 1) then
            write (status_fd, *)
        else if (ji/status_dump_int > 1) then
            write (status_fd, '('//int_to_char(ji/status_dump_int - 1)//'/)')
        end if
        flush (status_fd)

        call log_this('Opening trajectory file: '//traj_fname)
        call open_traj('write', traj_status)

        if (finish_prev_run) call log_this('Continuing with the incomplete previous run')
        if (append_flag_present) call log_this('Extending the existing trajectory @ '//traj_fname)

!$      call log_this('Using '//int_to_char(int(omp_get_max_threads(), kind=kind(jf)))//' OpenMP threads')

        do_status_dump = .not. (cmd_line_flag('-n') .or. cmd_line_flag('--no-status-dump'))

        if (cmd_line_flag('--no-evolve-motility')) then
            evolve_motility_bool = 0.0d0
            call log_this('Motility unit vectors are not being evolved')
        end if

        ! Associate intracellular force call with the user-chosen subroutine
        if (cmd_line_flag('--pl0')) then
            force => force_pl0
        else
            force => force_pl
        end if
    end subroutine prerun_setup
end module prerun
