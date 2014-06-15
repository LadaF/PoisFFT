#if dims==1
#define colons :
#elif dims==2
#define colons :,:
#else
#define colons :,:,:
#endif

    type(c_ptr), value :: D
    type(c_ptr), value :: Phi, RHS
    real(rp), pointer :: f_Phi(colons), f_RHS(colons)
    integer(c_int), optional :: ngPhi(dims), ngRHS(dims)
    integer :: i
    
    call c_f_pointer(D, f_D)
    
    if (present(ngPhi)) then
      call c_f_pointer(Phi, f_Phi, [(f_D%nxyz(i)+2*ngPhi(dims+1-i), i=1,dims)])
    else
      call c_f_pointer(Phi, f_Phi, f_D%nxyz)
    end if
    
    if (present(ngRHS)) then
      call c_f_pointer(RHS, f_RHS, [(f_D%nxyz(i)+2*ngRHS(dims+1-i), i=1,dims)] )
    else
      call c_f_pointer(RHS, f_RHS, f_D%nxyz)
    end if

    call Execute(f_D, f_Phi, f_RHS)
#undef colons