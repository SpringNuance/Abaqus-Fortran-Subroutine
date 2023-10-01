program print_array
  implicit none
  ! There are two ways to define an array
  ! integer, dimension(10) :: array = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
  integer :: array(10) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /)
  integer :: i

  do i = 1, size(array)
    if (mod(array(i), 2) == 0) then
      print *, array(i), "is divisible by 2"
    else
      print *, array(i), "is not divisible by 2"
    end if
  end do

end program print_array

! How to run the code
! gfortran print_array.f90 -o program 
! ./program.exe
