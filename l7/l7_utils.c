#include <signal.h>
#include <time.h>
#include "l7.h"
#include "l7p.h"

int L7_Get_Rank(void)
{
#ifdef HAVE_MPI
    return(l7.penum);
#else
    return(0);
#endif
}

int l7_get_rank_(void)
{
#ifdef HAVE_MPI
    return(l7.penum);
#else
    return(0);
#endif
}

int L7_Get_Numpes(void)
{
#ifdef HAVE_MPI
   return(l7.numpes);
#else
   return(1);
#endif
}

int l7_get_numpes_(void)
{
#ifdef HAVE_MPI
   return(l7.numpes);
#else
   return(1);
#endif
}

double L7_Wtime(void)
{
    double l7_wtime = 0.0;
#ifdef HAVE_MPI
    if (l7.initialized_mpi){
        l7_wtime = MPI_Wtime();
    }
    else {
        l7_wtime = (double)(clock()/CLOCKS_PER_SEC);
    }
#else
   l7_wtime = clock()/CLOCKS_PER_SEC;
#endif
    return(l7_wtime);
}

double l7_wtime_(void)
{
    double l7_wtime = 0.0;
#ifdef HAVE_MPI
    if (l7.initialized_mpi){
        l7_wtime = MPI_Wtime();
    }
    else {
        l7_wtime = (double)(clock()/CLOCKS_PER_SEC);
    }
#else
   l7_wtime = clock()/CLOCKS_PER_SEC;
#endif
    return(l7_wtime);
}

int L7_Get_Timeout_Signal(void)
{
    return(SIGURG);
}

int l7_get_timeout_signal_(void)
{
    return(SIGURG);
}

long long L7_Address(void *var)
{
   return((long long)var);
}
long long l7_address_(void *var)
{
   return((long long)var);
}

