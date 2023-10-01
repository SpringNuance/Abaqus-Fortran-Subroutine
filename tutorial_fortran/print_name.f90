program print_name

    implicit none
    character*20 :: name
    character (len = 20) :: f_name, l_name
    print *, "What's your first name "
    read *, f_name
    print *, "What's your last name "
    read *, l_name
    print *, "Hello ", trim(f_name), " ", trim(l_name)

end program print_name

! How to run the code
! gfortran print_name.f90 -o program 
! ./program.exe
