/* ---------------------------------------------------------------------
Author:     H. Carter Edwards 
            hcedwar@sandia.gov

Copyright:  Copyright (C) 1997   H. Carter Edwards
            Graduate Student
            University of Texas

Re-release: Copyright (C) 2011-2012   H. Carter Edwards

Purpose:    Domain paritioning based upon Hilbert Space-Filling Curve
            ordering.

License:    Re-release under the less-restrictive CLAMR software terms.
            Permitted by email with H. Carter Edwards on 9/13/2011

Disclaimer:

    These routines comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it
    under certain conditions. See License terms in file 'LICENSE'.
--------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>

//#include "cfsimple.h"
#include "hsfc.h"

/*--------------------------------------------------------------------*/
/* Make it callable from FORTRAN:
 *   Interface types: INTEGER and REAL*8
 */

//#define _hsfc2part  _F_NAME_(HSFC2PART,hsfc2part)
#define _hsfc2part  hsfc2part
//#define _fint       _F_INTEGER
#define _fint       int

#if   defined(_F_REAL_8)
#define _fdouble  _F_REAL_8
#elif defined(_F_DOUBLE)
#define _fdouble  _F_DOUBLE
#else
#define _fdouble  double
#endif

/*--------------------------------------------------------------------*/

#define MaxBits ( sizeof(unsigned) * CHAR_BIT )

#define SFC_NKEY  2

typedef struct {
  int      index ;         /* Original index */
  unsigned cell ;          /* Node is assigned to this cell */
  unsigned key[SFC_NKEY] ; /* SFC key of this node */
} sfc_node ;

typedef struct {
  int      part ;     /* Partition containing this cell */
  unsigned weight ;   /* Weight of nodes in the cell */
  unsigned addweight ;
  int      newcell ;
} sfc_grid ;

typedef void (*maptonorm)(
  const _fdouble * X , const _fdouble * Y , const _fdouble * Work );

/*--------------------------------------------------------------------*/

static int nodecellkeycomp( const void * X , const void * Y )
{
  const sfc_node * const x = (const sfc_node *) X ;
  const sfc_node * const y = (const sfc_node *) Y ;

  return ( x->cell   < y->cell   ) ? -1 :
       ( ( x->cell   > y->cell   ) ?  1 :
       ( ( x->key[0] < y->key[0] ) ? -1 :
       ( ( x->key[0] > y->key[0] ) ?  1 :
       ( ( x->key[1] < y->key[1] ) ? -1 :
       ( ( x->key[1] > y->key[1] ) ?  1 : 0 ) ) ) ) );
}

/*--------------------------------------------------------------------*/

void hsfc2part(
  const int      Level , /* IN: Background grid level of partitioning */
  const int      Limit , /* IN: Number of levels to consider for 'gaps' */
  const int      NPart , /* IN: Target number of partitions */
  const int      N ,     /* IN: Number of points */
  const double * X ,     /* IN: array of X-Coordinates */
  const double * Y ,     /* IN: array of Y-Coordinates */
  const int      ibase , /* IN: base - 0 for C, 1 for Fortran */
        int    * Info ,  /* IN:  Array of computational weights,
                              OUT: (1 <= LDInfo) [ Partitioning ]
                                   (2 <= LDInfo) [ Adjusted HSFC ordering ]
                                   (3 <= LDInfo) [ Original HSFC index, #1 ]
                                   (4 <= LDInfo) [ Original HSFC index, #2 ] */
        int      LDInfo )/* IN:  Leading dimension of Info */
{
  /*------------------------------------------------------------------*/

  const char name[] = "hsfc2part" ;

  const double   imax      = ((double) ~(0u)) ;
  const unsigned sfc_nkey  = SFC_NKEY ;

  const unsigned limit_gap = 01 << ( Limit + Limit );
  const unsigned level     = Level ; /* 2^level x 2^level grid */
  const unsigned ldinfo    = LDInfo ;
  const unsigned npt       = N ;
  const unsigned npart     = NPart ;
  const unsigned ngrid     = 01 << level ;
  const unsigned ncell     = 01 << ( level + level ); /* 2D */
  const unsigned shiftC    = MaxBits - level ;
  const unsigned shiftK    = MaxBits - ( level + level );

  /*------------------------------------------------------------------*/

  sfc_node *  const SFC  = (sfc_node *) malloc(sizeof(sfc_node) * npt  );
  sfc_grid *  const GRID = (sfc_grid *) malloc(sizeof(sfc_grid) * ncell);
  
  /*------------------------------------------------------------------*/

  double total_weight ;
  int i , ix , iy , ii ;

  /*------------------------------------------------------------------*/
  /* Initialize the GRID cells */

  if ( NULL == SFC || NULL == GRID ) {
    fprintf(stderr,"%s malloc failed, aborting\n",name);
    exit(-1);
  }

  for ( i = 0 ; (unsigned int)i < ncell ; ++i ) {
    sfc_grid * g = GRID + i ;
    g->part      = -1 ;
    g->newcell   = i ;
    g->weight    =
    g->addweight = 0 ;
  }

  /*------------------------------------------------------------------*/
  /* Fill SFC data structure and determine cell/total weight */

  total_weight = 0 ;

  for ( i = ii = ix = iy = 0 ; (unsigned int)i < npt ;
        ++i , ix++ , iy++ , ii += ldinfo ) {
    double xy[2] ;
    unsigned coord[2] ;

    xy[0] = X[ix] ;
    xy[1] = Y[iy] ;

    coord[0] = xy[0] * imax ;
    coord[1] = xy[1] * imax ;

    hsfc2d( coord , sfc_nkey , SFC[i].key );

    GRID[ SFC[i].cell = SFC[i].key[0] >> shiftK ].weight += Info[ii] ;

    total_weight += Info[ii] ;
  }

  /*------------------------------------------------------------------*/
  /* Identify SFC isolated cells */

  {
    int igap_beg[2] ;
    int igap_end[2] ;
    int jgap ;

    igap_beg[0] = igap_end[0] = 0 ;
    igap_beg[1] = - limit_gap ;

    for ( i = 0 ; (unsigned int)i < ncell && ! GRID[i].weight ; ++i );
    igap_end[1] = i ;
    jgap = 0 ;

    for ( ; (unsigned int)i < ncell ; ) {
      int gap ;

      for ( ; (unsigned int)i < ncell && GRID[i].weight ; ++i );

      for ( gap = 0 ; (unsigned int)i < ncell && ! GRID[i].weight ; ++i ) ++gap ;

      if ( limit_gap <= (unsigned int)gap ) {
        int j = jgap ^ 01 ;
        igap_beg[jgap] = i - gap ; /* Beginning of gap */
        igap_end[jgap] = i ;       /* End of gap       */

        if ( limit_gap <= (unsigned int)(igap_end[j] - igap_beg[j]) &&
             (unsigned int)(igap_beg[jgap] - igap_end[j]) <= limit_gap ) {
          int k ;
          for ( k = igap_end[j] ; k < igap_beg[jgap] ; ++k )
            if ( GRID[k].weight ) GRID[k].newcell = -1 ;
        }

        jgap = j ;
      }
    }
  }

  /*------------------------------------------------------------------*/
  /* Shift isolated cell's nodes to a neighbor */
  /* Precompute grid coordinate to SFC cell to reduce expense
     by a factor of four. */

  {
    sfc_grid ** const MAPG = (sfc_grid **) malloc(sizeof(sfc_grid *) * ncell);

    unsigned cell_shift_count = 0 ;
    unsigned fail_shift_count = 0 ;

    if ( NULL == MAPG ) {
      fprintf(stderr,"%s malloc failed, aborting\n",name);
      exit(-1);
    }

    for ( iy = 0 ; (unsigned int)iy < ngrid ; ++iy ) {
      const unsigned one  = 1 ;
      const unsigned icol = iy * ngrid ;
      unsigned coord[2] ;
      unsigned cell ;
      coord[1] = iy << shiftC ;
      for ( ix = 0 ; (unsigned int)ix < ngrid ; ++ix ) {
        coord[0] = ix << shiftC ;
        hsfc2d( coord , one , &cell );
        MAPG[ ix + icol ] = GRID + ( cell >> shiftK );
      }
    }

    for ( iy = 0 ; (unsigned int)iy < ngrid ; ++iy ) {
      const unsigned icol = iy * ngrid ;
      for ( ix = 0 ; (unsigned int)ix < ngrid ; ++ix ) {
        const unsigned ig = ix + icol ;
        sfc_grid * const g = MAPG[ ig ];

        if ( g->newcell == -1 ) { /* An SFC isolated cell */

          sfc_grid * neigh[4] ;
          unsigned max , imax=0 ;

          /* Edge neighbors */

          neigh[0] = (    0 < ix    ) ? ( MAPG[ig-1]     ) : NULL ;
          neigh[1] = (    0 < iy    ) ? ( MAPG[ig-ngrid] ) : NULL ;
          neigh[2] = ( (unsigned int)(ix+1) < ngrid ) ? ( MAPG[ig+1]     ) : NULL ;
          neigh[3] = ( (unsigned int)(iy+1) < ngrid ) ? ( MAPG[ig+ngrid] ) : NULL ;

          /* Find largest immediate neighbor for the cell shift */

          max = 0 ;
          for ( i = 0 ; i < 4 ; ++i )
            if ( neigh[i] && neigh[i]->weight &&
                 GRID + neigh[i]->newcell == neigh[i] &&
                 max < neigh[i]->weight ) {
              max = neigh[i]->weight ;
              imax = i ;
            }

          if ( max ) {
            g->newcell = neigh[imax] - GRID ;
            /* g->weight = 0 ; */
            neigh[imax]->addweight += g->weight ;
            cell_shift_count++ ;
          }
          else {
            fail_shift_count++ ;
          }
        }
      }
    }

    free( MAPG );

    if ( cell_shift_count || fail_shift_count ) {
      fprintf(stdout,
        "  Shifted %d cells, failed to shift %d cells, Total %d cells\n",
        cell_shift_count,fail_shift_count,ncell);
    }

    if ( fail_shift_count ) {
      /* Restore any failed shifts */

      for ( i = 0 ; (unsigned int)i < ncell ; ++i ) {
        if ( GRID[i].newcell == -1 ) GRID[i].newcell = i ;
      }
    }
  }

  /*------------------------------------------------------------------*/
  /* Partition the grid cells by weight,
     last-to-first so that Processor #0 has the least.
   */

  {
    const double target_weight = total_weight / ((double)npart) ;

    double current_weight = 0.0 ;
    int    current_part   = npart - 1 ;

    for ( i = ncell - 1 ; 0 <= i ; --i ) {
      if ( current_part && target_weight < current_weight ) {
        current_part-- ;
        current_weight = 0 ;
      }
      if ( GRID[i].newcell == i ) { /* Unmoved cells */
        current_weight += GRID[i].weight + GRID[i].addweight ;
      }
      GRID[i].part = current_part ;
    }
  }

  /*------------------------------------------------------------------*/
  /* Reassign nodes to cells, assign nodes to partitions */

  for ( ii = i = 0 ; (unsigned int)i < npt ; ++i , ii += ldinfo ) {
    SFC[i].cell = GRID[ SFC[i].cell ].newcell ;
    Info[ii] = GRID[ SFC[i].cell ].part ;
    SFC[i].index = i ;
  }

  /*------------------------------------------------------------------*/
  /* SFC Key output */
  
  if ( 3 < ldinfo && 1 < SFC_NKEY) {
    for ( ii = 2 , i = 0 ; (unsigned int)i < npt ; ++i , ii += ldinfo ) {
      Info[ii]   = SFC[i].key[0] ;
      Info[ii+1] = SFC[i].key[1] ;
    }
  }
  else if ( 2 < ldinfo ) {
    for ( ii = 2 , i = 0 ; (unsigned int)i < npt ; ++i , ii += ldinfo ) {
      Info[ii] = SFC[i].key[0] ;
    }
  }

  /*------------------------------------------------------------------*/
  /* Ordering output, sort by { cell , key } */

  if ( 1 < ldinfo ) {
    qsort( SFC , (size_t)npt , sizeof(sfc_node) , nodecellkeycomp );

    /* Output */

    for ( ii = 1 , i = 0 ; (unsigned int)i < npt ; ++i , ii += ldinfo ) {
      Info[ii] = SFC[i].index + ibase ; /* FORTRAN convention */
    }
  }

  /*------------------------------------------------------------------*/

  free( (void *) SFC );
  free( (void *) GRID );

  return ;
}

