Fortran:

      INTEGER*8 SMALocalIntArrayCreate(ID,SIZE,INITVAL)
      INTEGER*8 SMALocalFloatArrayCreate(ID,SIZE,INITVAL)

Example:

#include <SMAAspUserSubroutines.hdr>

      integer a(100)
      pointer(ptra,a)

      real*8 b(*)
      pointer(ptrb,b)

      ! create a local array with ID=1 and SIZE=100
      ptra = SMALocalIntArrayCreate(1,100) 
      a(1) = 11  ! use as a native Fortran array
      a(2) = 22  ! use as a native Fortran array

      ! create a local float array with ID=1 and SIZE=100, and
      ! initial value = -1.0
      ptrb = SMALocalFloatArrayCreate(1,100,-1.0) 

C++: 

       #include <SMAAspUserSubroutines.h>

       // Create a local integer array of with ID=1 and size=100
       int* a = SMALocalIntArrayCreate(1,100); 

       // Create a local float array of with ID=1, size=20, and 
       // initial value = -1.0
       real* b = SMALocalFloatArrayCreate(1,100,-1.0);   

NOTE: Float Arrays can store both SINGLE PRECISION 
      and DOUBLE PRECISION numbers. Internally, 
      memory is allocated in units of 64 bits (double/real*8).

NOTE: To resize an array, simply call Create() with the same ID, 
      but give it a new SIZE parameter. If the new size is larger, 
      the old data are copied over to the new array. No data are lost 
      during resizing.

For example:

      ! resize array with ID=1 to 300 integers
      ptra = SMALocalIntArrayCreate(1,300)  

NOTE: In Create() functions, there is an optional third
      argument -- initial value. If not supplied, all Int 
      arrays are initialized with INT_MAX ( 2,147,483,647 ). 
      All Float Arrays are initialized with Signaling NANs.  
      The values of INT_MAX and signaling NANs are accessible 
      via the 'SMAAspNumericLimits.h' and 'SMAAspNumericLimit.hdr' 
      header files.



      Fortran interface:

      INTEGER*8 SMALocalIntArrayAccess(ID)
      INTEGER*8 SMALocalFloatArrayAccess(ID)

Example:

#include <SMAAspUserSubroutines.hdr>

      integer a(100)
      pointer(ptra,a)

C Locate local Array(1) and associate a native array pointer with it

      ptra = SMALocalIntArrayAccess(1) 
      
      a(1) = 11 ! use as a native Fortran array
      a(2) = 22 ! use as a native Fortran array



C++ interface: 

      #include <SMAAspUserSubroutines.h>

      // Locate and open array with ID=1
      int* a =  SMALocalArrayIntAccess(1);
     
      a[1] = 11;  // use as a native array
      a[2] = 22;  // use as a native array

NOTE: If a request is made to access an array that has 
      not been created, the function will return 0.    


    Fortran interface:

      INTEGER*4 SMALocalIntArraySize(ID)
      INTEGER*4 SMALocalFloatArraySize(ID)

Example:

#include <SMAAspUserSubroutines.hdr>

      integer a_size, d_size

C Get the size of Array(1) as the number of INTEGERs
      a_size = SMALocalIntArraySize(1)   

! Get the size of Array(1) as the number of REALs

      d_size = SMALocalFloatArraySize(1) 

      do k=1,a_size
          ... 
      end do

C++:

      #include <SMAAspUserSubroutines.h>

      // Lookup the size of Array(1) as the number of ints

      int a_size = SMALocalIntArraySize(1);   

      // Lookup the size of Array(1) as the number of doubles
      
      int d_size = SMALocalFloatArraySize(1); 

      for(int i=1; i<=size; i++) {
          ...
      }  

    Fortran interface:

      INTEGER*4 SMALocalIntArraySize(ID)
      INTEGER*4 SMALocalFloatArraySize(ID)

Example:

#include <SMAAspUserSubroutines.hdr>

      integer a_size, d_size

C Get the size of Array(1) as the number of INTEGERs
      a_size = SMALocalIntArraySize(1)   

! Get the size of Array(1) as the number of REALs

      d_size = SMALocalFloatArraySize(1) 

      do k=1,a_size
          ... 
      end do

C++:

      #include <SMAAspUserSubroutines.h>

      // Lookup the size of Array(1) as the number of ints

      int a_size = SMALocalIntArraySize(1);   

      // Lookup the size of Array(1) as the number of doubles
      
      int d_size = SMALocalFloatArraySize(1); 

      for(int i=1; i<=size; i++) {
          ...
      }


      Fortran interface:

      INTEGER*8 SMAIntArrayCreate(ID,SIZE,INITVAL)
      INTEGER*8 SMAFloatArrayCreate(ID,SIZE,INITVAL)

Example:

#include <SMAAspUserSubroutines.hdr>

      integer a(100)
      pointer(ptra,a)
      double precision b(100)
      pointer(ptrb,b)

      ! create a global array with ID=1, SIZE=100, and
      ! INITVAL=-1.0
      ptra = SMAIntArrayCreate(1,100,-1) 

      a(1) = 11  ! use as a native Fortran array
      a(2) = 22  ! use as a native Fortran array

      ! create a global array with ID=2, SIZE=100, and
      ! INITVAL=-1.0
      ptrb = SMAFloatArrayCreate(2,100,-1.0) 

C++ interface: 

       #include <SMAAspUserSubroutines.h>

       // Create an integer array of with ID=1, size=100, 
       // and initial value=-1.0
       int* a = SMAIntArrayCreate(1,100,-1); 

       // Create a float array of with ID=2, size=20, 
       // and initial value=-1.0
       Real* b = SMAFloatArrayCreate(2,20,-1.0);   

NOTE: Float Arrays can store both SINGLE PRECISION and 
      DOUBLE PRECISION numbers. Internally, they
      allocate storage in 64-bit units (double/real*8).

NOTE: To resize an array, simply call Create() with the same ID, 
      but give it a new SIZE parameter. If the size has increased, 
      the old data will be copied over to the new array. 
      No data is lost during resizing.


For example:

      ! resize array with ID=1 to 300 integers
      ptra = SMAIntArrayCreate(1,300,-1)  

Fortran interface:

      INTEGER*8 SMAIntArrayAccess(ID)
      INTEGER*8 SMAFloatArrayAccess(ID)

Example:

#include <SMAAspUserSubroutines.hdr>

      integer a(100)
      pointer(ptra,a)

C Locate Array(1) and associate a native array pointer with it

      ptra = SMAIntArrayAccess(1) 
      
      a(1) = 11 ! use as a native Fortran array
      a(2) = 22 ! use as a native Fortran array

C++ interface: 

      #include <SMAAspUserSubroutines.h>

      // Locate and open array with ID=1
      int* a =  SMAIntArrayAccess(1);
     
      a[1] = 11;  // use as a native array
      a[2] = 22;  // use as a native array


NOTE: If a request is made to access an array which has 
      not been created, the function will return 0.



Fortran interface:

      INTEGER*8 SMAIntArrayAccess(ID)
      INTEGER*8 SMAFloatArrayAccess(ID)

Example:

#include <SMAAspUserSubroutines.hdr>

      integer a(100)
      pointer(ptra,a)

C Locate Array(1) and associate a native array pointer with it

      ptra = SMAIntArrayAccess(1) 
      
      a(1) = 11 ! use as a native Fortran array
      a(2) = 22 ! use as a native Fortran array

C++ interface: 

      #include <SMAAspUserSubroutines.h>

      // Locate and open array with ID=1
      int* a =  SMAIntArrayAccess(1);
     
      a[1] = 11;  // use as a native array
      a[2] = 22;  // use as a native array


NOTE: If a request is made to access an array which has 
      not been created, the function will return 0. 


Fortran:

#include <SMAAspUserSubroutines.hdr>

      call SMAIntArrayDelete(1) ! Delete global Array(1)

C++:

      #include <SMAAspUserSubroutines.h>

      SMAIntArrayDelete(1);  // Delete global Array(1)

NOTE: Deletion of arrays is optional. All storage allocated 
      for these arrays will be freed when Abaqus terminates
      (at the very end of the analysis). It is, however, a good 
      programming practice to delete all allocations explicitly, 
      especially when they are no longer needed, as this will 
      free up memory for use somewhere else. 


Fortran:

#include <SMAAspUserSubroutines.hdr>

      call SMAIntArrayDelete(1) ! Delete global Array(1)

C++:

      #include <SMAAspUserSubroutines.h>

      SMAIntArrayDelete(1);  // Delete global Array(1)

NOTE: Deletion of arrays is optional. All storage allocated 
      for these arrays will be freed when Abaqus terminates
      (at the very end of the analysis). It is, however, a good 
      programming practice to delete all allocations explicitly, 
      especially when they are no longer needed, as this will 
      free up memory for use somewhere else.
Fortran:

#include <SMAAspUserSubroutines.hdr>

      call SMAIntArrayDelete(1) ! Delete global Array(1)

C++:

      #include <SMAAspUserSubroutines.h>

      SMAIntArrayDelete(1);  // Delete global Array(1)

NOTE: Deletion of arrays is optional. All storage allocated 
      for these arrays will be freed when Abaqus terminates
      (at the very end of the analysis). It is, however, a good 
      programming practice to delete all allocations explicitly, 
      especially when they are no longer needed, as this will 
      free up memory for use somewhere else.