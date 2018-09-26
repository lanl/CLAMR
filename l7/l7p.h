#ifndef L7P_H_
#define L7P_H_

#include <sys/types.h>
#include <sys/time.h>
#include <stdio.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mpi.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef HAVE_QUO
#include <quo.h>

typedef unsigned int uint;

typedef struct QUO_SubComm {
   // context for holding state of QUO
   QUO_context context;
   // the SubComm communicator
   MPI_Comm comm;
   // my rank in the communicator
   int rank;
   // the size of the communicator
   int size;
   // whether or not I'm a member of the communicator
   int member;
} QUO_SubComm;
#endif

/*
 * Some L7 parameters.
 */

   /* these are the functions we will profile.  they MUST start at 0, be
    * unique integers, and increase +1 up to the last function.  if you add
    * a new function in to profile, add it to the bottom and choose the 
    * next number! Enum does this for us so we use it to enforce these rules.
    * To remove a function from profiling, set it to a negative number.
    */
enum
{
   L7_IO_PROF_OPEN      = 0,
   L7_IO_PROF_CLOSE,
   L7_IO_PROF_READ,
   L7_IO_PROF_WRITE,
   L7_IO_PROF_SEEK,
   L7_IO_PROF_TELL,
   L7_IO_PROF_FILESIZE,
   L7_IO_PROF_INQUIRE,
   L7_IO_PROF_LINK,
   L7_IO_PROF_UNLINK,
   L7_IO_PROF_SYMLINK,

   /*
    *  L7_MAX_FILE_FUNC CALLS must be equal to the number of functions
    * we are profiling - so the value of the #define above +1. Putting
    * it at the end of the enum forces this to be the case.
    */
   L7_MAX_FILE_FUNC_CALLS,
   L7_MIN_FILE_FUNC_CALLS = L7_IO_PROF_OPEN
};

#if !defined L7_EXTERN
#define L7_EXTERN extern
#endif

#include "l7_assert.h"

#if ! defined (HAVE_MPI)

#define l7_id_database int

#define Comm_Datatype  int

#define Comm_Op        int  /* Reduction operation type. */

#else
   

/* Some parameters. */

#define Comm_Datatype MPI_Datatype

#define Comm_Op       MPI_Op  /* Reduction operation type.  */

//#define L7_MPI_INT  MPI_INT
//#define L7_MPI_REAL MPI_FLOAT

#define L7_MAX_NUM_DBS  50 /* May be more space, but exceeding this
                              probably indicates a leak. */

/*
 * Message tag management.
 */

#define L7_SETUP_SEND_COUNT_TAG      1000
#define L7_SETUP_INDICES_NEEDED_TAG  1001

#define L7_UPDATE_TAGS_MIN           2001
#define L7_UPDATE_TAGS_MAX           2999

#define L7_MIN_MPI_REQS                50 /* Number of outstanding
                                             MPI_Requests initially
                                             allocated, times "num_recvs". */

/*
 * Struct for data associated with specified L7 handle.
 */
   
typedef struct l7_id_database
{
   int
     *global_ids,              /* Global ids associated with local indices.
                                  After l7p_gid2index has been executed,
                                  contains gids for needed indices as well  */
     *index_for_gid,           /* Local indices associated with global ids. */
     *indices_needed,          /* As input to L7_SETUP.                     */
     *indices_global_to_send,  /* Array of global indices this pe sends,    */
     *indices_local_to_send,   /* Array of local indices this pe sends.     */
     indices_to_send_len,      /* Length (in int) of indices_global_to_send
                                  and indices_local_to_send.                */
     indices_needed_len,       /* Allocated space for indices_needed.       */
     l7_id,                    /* As input to L7_SETUP.                     */
     mpi_request_len,          /* Allocated number of mpi_requests          */
     mpi_status_len,           /* Allocated number of mpi_statuses          */
     my_start_index,           /* As input to L7_SETUP.                     */
     num_indices_needed,       /* As input to L7_SETUP.                     */
     num_indices_owned,        /* Number of on-process indices.             */
     num_recvs,                /* Number of processes this pe recvs from.   */
     num_reqs_outstanding,     /* Num MPI_Requests posted in packing model. */
     num_sends,                /* Number of processes this pe sends to.     */
     numpes,                   /* Size of MPI_COMM_WORLD                    */
     penum,                    /* Process rank in MPI_COMM_WORLD            */
     *recv_from,               /* Processes this pe receives from.          */
     recv_from_len,            /* Length (in int) of recv_from.             */
     *recv_counts,             /* Array of msg counts for recv_from pes.    */
     recv_counts_len,          /* Length (in int) of recv_counts.           */
     *send_to,                 /* Array of processes this pe sends to.      */
     send_to_len,              /* Length (in int) of send_to.               */
     *send_counts,             /* Msg counts for send_to pes.               */
     send_counts_len,          /* Length (in int) of send_counts_len.       */
     *starting_indices,        /* Array of my_start_index from each pe.     */
     this_tag_update;          /* Msg tag for updates.                      */
   
   /* MPI parameters */
   
   MPI_Request
     *mpi_request;
   
   MPI_Status
     *mpi_status;
   
#ifdef HAVE_OPENCL
   int
     num_indices_have,         /* Count of indices needed for send in update */
     *indices_have;            /* list of indices have on pe for send       */

   cl_mem dev_indices_have;    /* list of indices on the device             */
#endif

   struct l7_id_database
     *next_db;                 /* Link to next database.                    */
   
} l7_id_database;

typedef struct l7_push_id_database
{
   int
     num_comm_partners,        /* Number of processors to communicate with. */
     *comm_partner,            /* List of processors to communicate with.
                                  Length of this and next two arrays will 
                                  will be num_comm_partners                 */
     *send_buffer_count,       /* Count of data to send to each processor   */
     *recv_buffer_count,       /* Count of data to receive from each proc   */
     **send_database,          /* List of data indices to send to each
                                  processor. This is a 'ragged right' array
                                  with first index num_comm_partners and
                                  the second send_buffer_count              */
     receive_count_total,      /* Total receive size for each processor.
                                  This is a sum of the recv_buffer_count
                                  array                                     */
     **send_buffer,            /* A preallocated int datatype buffer to be
                                  used for packing data to send in
                                  L7_Push_Update. It is allocated in setup
                                  and deallocated in free.                  */
     l7_push_id;               /* As input to L7_PUSH_SETUP.                */
   
   struct l7_push_id_database
     *next_push_db;            /* Link to next database.                    */
   
} l7_push_id_database;

#endif /* HAVE_MPI */

/*
 * main structure for L7.
 */
       
typedef struct
{
   /* Some workspace, for message data. */
   
   void
     *send_buffer,
     *workspace;
   
   int
     sizeof_send_buffer,
     sizeof_workspace;
   
   l7_id_database
     *first_db,                /* For linked list of dbs.              */
     *last_db;
   
   l7_push_id_database
     *first_push_db,           /* For linked list of push dbs.         */
     *last_push_db;
   
   int
     data_check_len,           /* Number of bytes in array data_check. */
     initialized,              /* 1 if L7 initialized, else 0          */
     initialized_mpi,          /* 1 if L7 initialized MPI, else 0      */
     mpi_initialized,          /* 1 if L7_init sets use mpi, else 0    */
     num_dbs,                  /* Number of databases allocated.       */
     num_push_dbs,             /* Number of push databases allocated.  */
     numpes,                   /* Number of processors in mpi job      */
     penum;                    /* Process id for currently set db.     */

#ifdef HAVE_QUO
   QUO_SubComm subComm;
#endif
   
   void
     *data_check;              /* Workspace for use in l7_update_check */
   
   FILE
     *assert_out_file,         /* output file for MAYAP_ASSERT.
                                * if NULL, default is stderr           */
     *assert_out_file_none;    /* /dev/null file pointer.
                                * temporarily set assert_out_file
                                * to assert_out_file_none when you
                                * wish to disable error messages.
                                * When done, reset assert_outfile to
                                * NULL.
                                */
#ifdef HAVE_OPENCL
   cl_kernel
     kernel_pack_int_have_data,
     kernel_pack_float_have_data,
     kernel_pack_double_have_data,
     kernel_copy_ghost_int_data,
     kernel_copy_ghost_float_data,
     kernel_copy_ghost_double_data;
#endif
   
   int
     io_profiling_level;       /* L7_IO_PROF_OFF / SIMPLE / VERBOSE    */
} l7_globals;

L7_EXTERN l7_globals l7;

/* for simple IO profiling, this is our struct we use for each function name */

struct simple_timing {
   int num_calls;
   double total_time;
};

/* for the verbose IO profiling we have this much more, well, verbose and
 * full data structure to record everything */

struct verbose_timing {
   struct verbose_timing *next; /* linked list - the next guy in the chain */
   int fname_index;             /* the function id (L7_IO_PROF_OPEN ...) */
   int file_id;                 /* the file ID associated w/ this op */
   struct timeval time_in;
   struct timeval time_out;
   enum L7_DiskPatternType  l7_disk_pattern;  /* like distributed/replicated */
   /* these next three are left open to be used or not used by each timing.
    * some of the functions will want to store things in them while others
    * will not
    */
   long long int var_long;
   int           var_int;
   char          *var_string;
};

struct verbose_fd {
   struct verbose_fd *next;
   int file_id;
   char *file_name;
   /* fds are valid between an open and a close, after it's closed the FD can
    * be reused by the system so we need to record what times the open/close
    * was done
    */
   struct timeval valid_start;
   struct timeval valid_end;
};

union l7_file_handler {
   int  fd;
   FILE *fstream;
};

struct l7_file_record {
   struct l7_file_record    *next;
   struct l7_file_record    *prev;
   int                      file_id;
   enum L7_DiskPatternType  l7_disk_pattern;
   char                     *file_name;
   union l7_file_handler    file_handler;
};


/*
 * C Private Prototypes.
 */

MPI_Datatype l7p_mpi_type (
      const enum L7_Datatype  l7_datatype
      );

int l7p_sizeof (
      const enum L7_Datatype  l7_datatype
      );

l7_id_database *l7p_set_database(
      const int l7_id
      );

/*
 * L7 File Private Prototypes.
 */


void l7p_io_prof_vfdprint(void);

void l7p_io_prof_vfdrecord(
      const int            file_id,
      const int            func_idx,
      const char           *file_name,
      const struct timeval begin,
      const struct timeval end
      );

char *l7p_io_prof_get_filename(
      const struct verbose_timing *cur
      );

void l7p_io_prof_print_timings(void);

void l7p_io_prof_vrecord(
      const int                      fname_index,
      const struct timeval           begin_tv,
      const struct timeval           end_tv,
      const int                      file_id,
      const enum L7_DiskPatternType  l7_disk_pattern,
      const long long                var_long,
      const int                      var_int,
      const char                     *var_string
      );

void l7p_io_prof_record(
      const int            fname_index,
      const struct timeval begin_tv,
      const struct timeval end_tv
      );

void l7p_io_prof_init(
      const int new_io_profiling_level
      );

void l7p_io_prof_shutdown(void);

int l7p_file_record_new(
      const char  *fname,
      const int   disk_pattern,
      const union l7_file_handler new_file_handler
      );

struct l7_file_record* l7p_file_record_get(
      const int file_id
      );

void l7p_file_record_delete(
      const int file_id
      );

void l7p_file_shutdown(void);

void l7p_file_record_printall(void);

long long l7p_readfd(
      const int  fd,
      void       *buf,
      const      size_t num
      );

long long l7p_writefd(
      const int    fd,
      void         *buf,
      const size_t num
      );

long long l7p_fseek(
      const int       fid,
      const long long disk_loc
      );

long long l7p_ftell(
      const int fid
      );

/*
 * End prototypes.
 */

/*
 * remove typesafe linkage if compiling under c++
 */

#ifdef __cplusplus
}
#endif

#endif /*L7P_H_*/
