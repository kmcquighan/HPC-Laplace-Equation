// PDLAPRNT2 Usage and code is exactly the same as PDLAPRNT
// except that the file argument is of type MPI_File.
//
// Written by Kelly McQuighan, copied exactly from PDLAPRNT
// last modified 8/2/13
//
// ************NOTES*****************
//
// assumes the matrix is column-major
// output is row-major
// 
//***********USAGE*************** 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>

#include <blacs_headers.h>
#include <pblas_headers.h>
#include <scalapack_tools_headers.h>
#include "pblasIOtools.h"
#include "stringTools.h"

static int imax( int a, int b ){
  if (a>b) return(a); else return(b);
}
static int imin( int a, int b ){
  if (a<b) return(a); else return(b);
}

//**********PDLAPRNT2*********************
void pdlaprnt2(int* msub, int* nsub, double* A, int* ia, int* ja, int* descA, int* irprnt, int* icprnt, char* mat_name, FILE* outfile, double* work, int len_mat_name)
{
int lwork, ictxt, mm, isIOprocessor, nprow, npcol, myrow, mycol;
int nn, mb, nb, rsrc, csrc, ldd, descWORK[9], info;
int istart, iend, isize, jstart, jend, jsize, i, j, ione=1;
int glob_i, glob_j, pad_i, pad_j, len_glob_i, len_glob_j, ip;
double alpha, beta;
char num[33], glob_ic[33], glob_jc[33], pad_ic[33], pad_jc[33];

  lwork = descA[4];
  ictxt = descA[1];
  rsrc = descA[6];
  csrc = descA[7];
  Cblacs_gridinfo ( ictxt, &nprow, &npcol, &myrow, &mycol );
  mm = imax( 1, imin(*msub, lwork) );
  nn = imax( 1, (int) lwork/mm );
  mb = mm;
  nb = nn;
  isIOprocessor = ( ( myrow==*irprnt ) && ( mycol==*icprnt ) );
  ldd = imax( 1, mm );
  descinit_( descWORK, &mm, &nn, &mb, &nb, &rsrc, &csrc, &ictxt, &ldd, &info );

  for (jstart=*ja; jstart < *ja+*nsub; jstart+=nn) { //cols
    jend = imin( *ja+*nsub-1, jstart+nn-1 );
    jsize = jend - jstart + 1;
    for (istart=*ia; istart < *ia+*msub; istart+=mm) { //rows
      iend = imin( *ia+*msub-1, istart+mm-1 );
      isize = iend - istart + 1;
      alpha = 1.0e0; beta = 0.0e0;
      pdgeadd_("N", &isize, &jsize, &alpha, A, &istart, &jstart, descA, &beta, work, &ione, &ione, descWORK );

      for (j=0; j<jsize; j++ ) { //cols
        for (i=0; i<isize; i++ ) { //rows
          if ( isIOprocessor ) {
            glob_i = istart+i; glob_j = jstart+j;
            len_glob_i = itoa_len(glob_i, glob_ic); len_glob_j = itoa_len(glob_j, glob_jc);
            pad_i = 6-len_glob_i; pad_j = 6-len_glob_j;
            if (pad_i<0) pad_i=0; if (pad_j<0) pad_j=0;
            for (ip=0; ip<pad_j; ip++) pad_jc[ip] = ' ';
            for (ip=0; ip<pad_i; ip++) pad_ic[ip] = ' ';
            pad_jc[pad_j] = '\0'; pad_ic[pad_i] = '\0';
            
            dtoa(*(work+j*ldd+i), num, "f", 6);
            fprintf ( outfile,"%s(%s%d, %s%d) =\t%s\n", mat_name, pad_ic, glob_i, pad_jc, glob_j, num );
          }
        }
      }
    }
  }

return;
}
//**********PSLAPRNT2*********************
void pslaprnt2(int* msub, int* nsub, float* A, int* ia, int* ja, int* descA, int* irprnt, int* icprnt, char* mat_name, FILE* outfile, float* work, int len_mat_name)
{
int lwork, ictxt, mm, isIOprocessor, nprow, npcol, myrow, mycol;
int nn, mb, nb, rsrc, csrc, ldd, descWORK[9], info;
int istart, iend, isize, jstart, jend, jsize, i, j, ione=1;
int glob_i, glob_j, pad_i, pad_j, len_glob_i, len_glob_j, ip;
float alpha, beta;
char num[33], glob_ic[33], glob_jc[33], pad_ic[33], pad_jc[33];

  lwork = descA[4];
  ictxt = descA[1];
  rsrc = descA[6];
  csrc = descA[7];
  Cblacs_gridinfo ( ictxt, &nprow, &npcol, &myrow, &mycol );
  mm = imax( 1, imin(*msub, lwork) );
  nn = imax( 1, (int) lwork/mm );
  mb = mm;
  nb = nn;
  isIOprocessor = ( ( myrow==*irprnt ) && ( mycol==*icprnt ) );
  ldd = imax( 1, mm );
  descinit_( descWORK, &mm, &nn, &mb, &nb, &rsrc, &csrc, &ictxt, &ldd, &info );

  for (jstart=*ja; jstart < *ja+*nsub; jstart+=nn) { //cols
    jend = imin( *ja+*nsub-1, jstart+nn-1 );
    jsize = jend - jstart + 1;
    for (istart=*ia; istart < *ia+*msub; istart+=mm) { //rows
      iend = imin( *ia+*msub-1, istart+mm-1 );
      isize = iend - istart + 1;
      alpha = 1.0e0; beta = 0.0e0;
      psgeadd_("N", &isize, &jsize, &alpha, A, &istart, &jstart, descA, &beta, work, &ione, &ione, descWORK );

      for (j=0; j<jsize; j++ ) { //cols
        for (i=0; i<isize; i++ ) { //rows
          if ( isIOprocessor ) {
            glob_i = istart+i; glob_j = jstart+j;
            len_glob_i = itoa_len(glob_i, glob_ic); len_glob_j = itoa_len(glob_j, glob_jc);
            pad_i = 6-len_glob_i; pad_j = 6-len_glob_j;
            if (pad_i<0) pad_i=0; if (pad_j<0) pad_j=0;
            for (ip=0; ip<pad_j; ip++) pad_jc[ip] = ' ';
            for (ip=0; ip<pad_i; ip++) pad_ic[ip] = ' ';
            pad_jc[pad_j] = '\0'; pad_ic[pad_i] = '\0';
            
            stoa(*(work+j*ldd+i), num, "f", 6);
            fprintf ( outfile,"%s(%s%d, %s%d) =\t%s\n", mat_name, pad_ic, glob_i, pad_jc, glob_j, num );
          }
        }
      }
    }
  }

return;
}
