    type(c_ptr), intent(inout) :: D
    
    call c_f_pointer(D, f_D)
    call Finalize(f_D)
    deallocate(f_D)
    D = c_null_ptr;
