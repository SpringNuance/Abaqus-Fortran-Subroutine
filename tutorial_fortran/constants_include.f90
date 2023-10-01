PROGRAM MainProgram
    INCLUDE 'constants_include.inc'
    
    REAL :: radius, area
    
    PRINT *, "Enter the radius of a circle:"
    READ *, radius
    
    area = pi * radius**2
    PRINT *, "The area of the circle is:", area
END PROGRAM MainProgram

! How to run the code
! gfortran constants_include.f90 -o program 
! ./program.exe