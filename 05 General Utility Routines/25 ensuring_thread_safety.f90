Fortran:
#include <SMAAspUserSubroutines.hdr>
            
       ! Initialization in UEXTERNALDB/VEXTERNALDB

       call MutexInit( 1 )      ! initialize Mutex #1

       ! Use in all other user subs after being initialized

        call MutexLock( 1 )      ! lock Mutex #1

        < critical section : update shared variables >

       call MutexUnlock( 1 )   ! unlock Mutex #1
C++:
       #include <SMAAspUserSubroutines.h>

       // Initialization in UEXTERNALDB/VEXTERNALD

       MutexInit( 1 );         // initialize Mutex #1

       // Use in all other user subs after being initialized

       MutexLock( 1 );         // lock Mutex #1

       < critical section : update shared variables >
       
       MutexUnlock( 1 );       // unlock Mutex #1
            
NOTE: IDs are arbitrary chosen by the user, from the pool of 1-100.
      Other threads, when encountering a locked mutex, will sleep.
      Once the first entering thread unlocks the mutex and leaves, 
      other threads will be able to come in and execute the critical 
      section (one at a
time).