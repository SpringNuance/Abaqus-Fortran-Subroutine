MODULE CircleModule
  IMPLICIT NONE

  TYPE :: Circle
     REAL :: radius = 0.0
  CONTAINS
     PROCEDURE :: set_radius
     PROCEDURE :: get_radius
     PROCEDURE :: compute_area
  END TYPE Circle

CONTAINS

  SUBROUTINE set_radius(self, r)
    CLASS(Circle), INTENT(INOUT) :: self
    REAL, INTENT(IN) :: r
    self%radius = r
  END SUBROUTINE set_radius

  FUNCTION get_radius(self) RESULT(r)
    CLASS(Circle), INTENT(IN) :: self
    REAL :: r
    r = self%radius
  END FUNCTION get_radius

  FUNCTION compute_area(self) RESULT(area)
    CLASS(Circle), INTENT(IN) :: self
    REAL :: area
    area = 3.14159 * self%radius * self%radius
  END FUNCTION compute_area

END MODULE CircleModule

PROGRAM TestCircle
  USE CircleModule
  IMPLICIT NONE

  TYPE(Circle) :: myCircle
  REAL :: area

  ! Set radius
  CALL myCircle%set_radius(5.0)

  ! Get radius
  PRINT *, "Radius of the circle:", myCircle%get_radius()

  ! Compute area
  area = myCircle%compute_area()
  PRINT *, "Area of the circle:", area

END PROGRAM TestCircle

! How to run the code
! gfortran circle_class.f90 -o program 
! ./program.exe
