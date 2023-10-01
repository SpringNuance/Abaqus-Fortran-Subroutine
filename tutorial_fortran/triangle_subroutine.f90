PROGRAM Triangle
    IMPLICIT NONE
    REAL :: a, b, c, TriangleArea
    PRINT *, 'Welcome, please enter the &
            &lengths of the 3 sides.'
    PRINT *, 'Side a: '
    READ *, a
    PRINT *, 'Side b: '
    READ *, b
    PRINT *, 'Side c: '
    READ *, c
    CALL Area(a,b,c, TriangleArea)
    PRINT *, 'Triangle''s area: ', TriangleArea
END PROGRAM Triangle

SUBROUTINE Area(x,y,z, TriangleArea)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x, y, z
    REAL, INTENT(OUT) :: TriangleArea
    REAL :: theta, height
    theta = ACOS((x**2+y**2-z**2)/(2.0*x*y))
    height = x*SIN(theta) 
    TriangleArea = 0.5*y*height
END SUBROUTINE Area

! How to run the code
! gfortran triangle_subroutine.f90 -o program 
! ./program.exe
