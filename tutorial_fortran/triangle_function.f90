PROGRAM Triangle
    IMPLICIT NONE
    REAL :: a, b, c, Area
    PRINT *, 'Welcome, please enter the &
            &lengths of the 3 sides.'
    PRINT *, 'Side a: '
    READ *, a
    PRINT *, 'Side b: '
    READ *, b
    PRINT *, 'Side c: '
    READ *, c
    PRINT *, 'Triangle''s area: ', Area(a,b,c)
END PROGRAM Triangle

FUNCTION Area(x,y,z)
    IMPLICIT NONE
    REAL :: Area            ! function type
    REAL, INTENT( IN ) :: x, y, z
    REAL :: theta, height
    theta = ACOS((x**2+y**2-z**2)/(2.0*x*y))
    height = x*SIN(theta)
    Area = 0.5*y*height
END FUNCTION Area

! How to run the code
! gfortran triangle_function.f90 -o program 
! ./program.exe
