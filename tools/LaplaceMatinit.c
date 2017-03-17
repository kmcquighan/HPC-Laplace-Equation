/*
desc	Name	Description	Value	Scope
0	DTYPE_	Descriptor type	DTYPE_ = 1	Global
1	CTXT_	BLACS context	CTXT_ = ictxt	Global
2	M_	    Number of rows in the global matrix	M_ = max(m, 0)	Global
3	N_	    Number of columns in the global matrix	N_ = max(n, 0)	Global
4	MB_	    Row block size	MB_ = max(mb, 1)	Global
5	NB_	    Column block size	NB_ = max(nb, 1)	Global
6	RSRC_	The process row of the p × q grid over which the first row of the global matrix is distributed	RSRC_ = max(0, min(irsrc, p-1))	Global
7	CSRC_	The process column of the p × q grid over which the first column of the global matrix is distributed	CSRC_ = max(0, min(icsrc, q-1))	Global
8	LLD_	The leading dimension of the local array	LLD_ = max(lld, max(1, LOCp(M_)))	Local
*/

/* Right now, this is supported for block sizes equal to the number of grid points. However, the structure is
   in place so that it is easily extendable to arbitrary block sizes.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <blacs_headers.h>
#include <pblas_headers.h>
#include <scalapack_tools_headers.h>
#include <scalapack_headers.h>
#include "LaplaceMatinit.h"

#define PI 3.14159265359

void Ablock(double* AA, int ia, int ja, int mloc, int mb, double fact);
void Iblock(double* AA, int ia, int ja, int mloc, int mb, double fact);
void Nblock(double* AA, int ia, int ja, int mloc, int mb, double fact);

//****************************************
// Kelly McQuighan 
//
//***************************************
void LaplaceMatinit( double* AA, int* descA, double* b, int* descB, int myrow, int mycol, int nprow, int npcol )
{
int m, ioffc, ioffr, mb, ic, ir, mloc, nloc;
int moff, noff, nend, mend, i, jk;
int izero=0, rblockglob, cblockglob, mendglob, nendglob;
double dmloc, dnloc, dmb, dm, h, x, y, fact;

  m = descA[2]; // n=m by assumption
  mb = descA[4]; // nb = mb by assumption
  
  mloc    = numroc_( &m   , &mb    , &myrow, &izero, &nprow );
  nloc    = numroc_( &m   , &mb    , &mycol, &izero, &npcol );
  dmb = (double) mb; dmloc = (double) mloc; dnloc = (double) nloc; dm = (double) m;

  moff   = 0; // my block number for row offset
  noff   = 0; // my block number for column offset
  mend   = ceil( dmloc / dmb ); // ending row block index, local # of blocks
  nend   = ceil( dnloc / dmb ); // end column block index, local # of blocks
  mendglob = ceil( dm / dmb ); nendglob = mendglob;

  h = 2.0*PI / dmb;
  fact = 1.0 / (h*h);

// Build Laplace matrix  
  for ( ic = noff; ic < nend; ic++ ){
    cblockglob = ic*npcol+mycol;
    ioffc = ic*mb;
    for ( ir = moff; ir < mend; ir++ ){
      rblockglob = ir*nprow+myrow;
      ioffr = ir*mb;
      if ( rblockglob == cblockglob ) Ablock(AA, ioffr, ioffc, mloc, mb, fact);
      else if ( ((rblockglob-1)==cblockglob) || ((rblockglob+1)==cblockglob) ) Iblock(AA, ioffr, ioffc, mloc, mb, fact);
      else if ( ((rblockglob==0)&&(cblockglob==(nendglob-1))) || ((rblockglob==(mendglob-1))&&(cblockglob==0)) ) Iblock(AA, ioffr, ioffc, mloc, mb, fact);
      else Nblock(AA, ioffr, ioffc, mloc, mb, fact);
    }
  }

 // Build RHS
  if (mycol ==0) {
    for ( ir = moff; ir < mend; ir++ ) {
      for ( i=0; i<mb; i++ ) {
        y = i*h;
        x = (ir*nprow+myrow)*h;
        b[ir*mb+i] = -2*sin(x)*sin(y);
      }
    }
  }

  // Ensure uniqueness
  if (myrow==0) {
    jk = 0;
    for (ic= noff; ic < nend; ic++ ) {
      for (i=0; i<mb; i++) { AA[jk*mb] = 0.0; jk++; }
    }
    if (mycol==0) { b[0] = 0.0; AA[0] = 1.0 * fact; }
  }

return;
}

//***************SUBFUNCTIONS**************************
//
//****************A block*****************************
void Ablock(double* AA, int ia, int ja, int mloc, int mb, double fact){
int i, j;
  for (j=0; j<mb; j++) {//col loop
    for (i=0; i<mb; i++) { // row loop
      if (i==j) AA[(ja+j)*mloc+(ia+i)] = -4.0e0 * fact;
      else if ( ((i-1)==j) || ((i+1)==j) ) AA[(ja+j)*mloc+(ia+i)] = 1.0e0 * fact;
      else if ( ((i==0)&&(j==(mb-1))) || ((i==(mb-1))&&(j==0)) ) AA[(ja+j)*mloc+(ia+i)] = 1.0e0 * fact;
      else AA[(ja+j)*mloc+(ia+i)] = 0.0e0;
    }
  }

return;	
}
//
//*************Identity block**********************
void Iblock(double* AA, int ia, int ja, int mloc, int mb, double fact){
int i, j;
  for (j=0; j<mb; j++) {//col loop
    for (i=0; i<mb; i++) { // row loop
      if (i == j) AA[(ja+j)*mloc+(ia+i)] = 1.0e0 * fact;
      else AA[(ja+j)*mloc+(ia+i)] = 0.0e0;
    }
  }

return;	
}
//
//**************Zero block*************************
void Nblock(double* AA, int ia, int ja, int mloc, int mb, double fact){
int i, j;
  for (j=0; j<mb; j++) {//col loop
    for (i=0; i<mb; i++) { // row loop
      AA[(ja+j)*mloc+(ia+i)] = 0.0e0;
    }
  }

return;
}
