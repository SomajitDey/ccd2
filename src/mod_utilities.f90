module utilities
implicit none
private :: div_rem 
contains

!!!!!!!!/////// Gaussian random no. generator \\\\\\\\\\\\\\\\\\\\\\\\\\\
! Polar rejection method (Knop[1969]) related to Box-Muller transform
 Subroutine gasdev(g,mean,variance) 
      DOUBLE PRECISION,INTENT(IN)::variance,mean
      DOUBLE PRECISION,INTENT(OUT)::g(:)
      DOUBLE PRECISION:: fac,rsq,v1,v2
      DOUBLE PRECISION:: rands(2)
      INTEGER:: i, size_g
      double precision, parameter :: small=epsilon(0.0d0)
      
      size_g=size(g)
     
     harvest_g_array: DO i=1,size_g,2
        DO
        CALL RANDOM_NUMBER(rands)
        v1=2.0d0*rands(1)-1.0d0
        v2=2.0d0*rands(2)-1.0d0
        rsq=v1**2+v2**2
        if((rsq<1.0d0).AND.(rsq>small))EXIT
        ENDDO
        fac=variance*DSQRT(-2.0d0*dlog(rsq)/(rsq))
        g(i)=v1*fac+mean
        if(i/=size_g) g(i+1)=v2*fac+mean 
     END DO harvest_g_array
      END subroutine gasdev



! Brief: This subroutine outputs the time spent in secs (cpu & wall clock) since its previous invocation
! Note: All arguments are optional and of type real
! Note: CPU usage = cpu x 100 % / wclock
! Note: #Threads = nint(cpu/wclock)

!TODO: Rename timestamp -> stopwatch. Add split capture mechanism: split_cpu, split_wclock.
! cpu and wclock would give time spent since first/resetting invocation
! Add logical reset flag as well
subroutine timestamp(cpu, wclock)
    real, intent(out), optional :: cpu, wclock
    integer, save :: sys_clock_count_prev
    integer :: sys_clock_count, sys_clock_max, sys_clock_rate, sys_clock_diff
    real, save :: cpu_sec_prev
    real :: cpu_sec
    integer :: call_count = 1 ! implicit save attribute

    call cpu_time(cpu_sec)
    
    if (present(cpu)) then
        if (call_count == 1) then
            cpu = 0.0
        else
            cpu = cpu_sec - cpu_sec_prev
         end if
    end if
    cpu_sec_prev = cpu_sec

    call system_clock(sys_clock_count, sys_clock_rate, sys_clock_max)

    if (present(wclock)) then
        if (call_count == 1) then
            sys_clock_diff = 0
        else
            sys_clock_diff = sys_clock_count - sys_clock_count_prev
        end if
        
        if (sys_clock_diff < 0) then
            wclock = real(sys_clock_diff + sys_clock_max) / sys_clock_rate
        else
            wclock = real(sys_clock_diff) / sys_clock_rate
        end if
    end if
    sys_clock_count_prev = sys_clock_count
     
    call_count = call_count + 1
end subroutine timestamp

    ! Transforms real seconds into human readable format
    pure character(len=15) function dhms(sec)
        real, intent(in) :: sec
        integer :: days,hrs,mins,secs
     
        secs=int(sec)
        call div_rem(secs, 3600*24, days)
        call div_rem(secs, 3600, hrs)
        call div_rem(secs, 60, mins)   
     
        write(dhms,"(i0,'d',1x,i0,'h',1x,i0,'m',1x,i0,'s')") days, hrs, mins, secs
    end function dhms
     
    pure subroutine div_rem(divided_becomes_remainder, divisor, quotient)
        integer, intent(inout) :: divided_becomes_remainder
        integer, intent(in) :: divisor
        integer, intent(out) :: quotient
    
        quotient = divided_becomes_remainder / divisor
        divided_becomes_remainder = mod(divided_becomes_remainder, divisor)
    end subroutine div_rem

    ! Outputs the sha1 hash of any given file
    character(len=40) function sha1(fname)
        character(len=*), intent(in) :: fname
        character(len=*), parameter :: tmpfile = '.sha1.tmp'
        integer :: tmpunit, ios

        call execute_command_line("sha1sum "//fname//" 2>/dev/null > "//tmpfile)
        open(newunit=tmpunit, file=tmpfile, status='old', action='read')
            read(tmpunit, '(a)', iostat=ios) sha1
            if (ios /= 0) sha1=''
        close(tmpunit, status='delete')
    end function sha1

    function int_to_char(intarg)
        integer, intent(in) :: intarg
        character(len=floor(log10(real(intarg)))+1) :: int_to_char
        write(int_to_char,'(i0)') intarg
    end function int_to_char

    ! Returns true if `flag` is present as a command line flag/option/argument
    logical function cmd_line_flag(flag)
        character(len=*), intent(in) :: flag
        character(len=:), allocatable :: cmd_line
        integer :: cmd_line_length
        character(len=*), parameter :: delimiter=' ' ! Flags in a command line are always delimited
        
        call get_command(length=cmd_line_length)
        allocate(character(len=cmd_line_length) :: cmd_line)
        call get_command(command=cmd_line)
        
        cmd_line_flag = index(cmd_line//delimiter, delimiter//flag//delimiter) /= 0
    end function cmd_line_flag
end module utilities
