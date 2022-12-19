program ccd_show_params
    use parameters
    use files
    
    call assign_params(params_fname)
    write(*,nml=params)
end program ccd_show_params
