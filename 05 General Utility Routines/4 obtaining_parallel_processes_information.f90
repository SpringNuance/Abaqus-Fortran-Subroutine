CALL GETNUMCPUS( NUMPROCESSES )
CALL VGETNUMCPUS( NUMPROCESSES )
...
CALL GETRANK( KPROCESSNUM )
CALL VGETRANK( KPROCESSNUM )
...

Fortran:
#include <SMAAspUserSubroutines.hdr>

         integer numThreads 
         numThreads = GETNUMTHREADS()

C++:
         #include <SMAAspUserSubroutines.h>

         int numThreads = GETNUMTHREADS();


FOTRAN:
#include <SMAAspUserSubroutines.hdr>    

         INTEGER myThreadID

         myThreadID = get_thread_id()

C++:
         #include <SMAAspUserSubroutines.h>

         int myThreadID = get_thread_id()


#include <SMAAspUserSubroutines.hdr>

         integer ABA_COMM_WORLD

         ABA_COMM_WORLD = GETCOMMUNICATOR()

         if (ABA_COMM_WORLD.ne.0) then
             ...do some parallel work, using MPI ...
         else
             ...do some work in a single process ...
         end if

#include <mpi.h>
#include <SMAAspUserSubroutines.h>

MPI_Comm ABA_COMM_WORLD = get_communicator()

if (ABA_COMM_WORLD) {
// ... do some parallel work, using MPI ... 
}
else{
// ...do some work in a single process ...
}