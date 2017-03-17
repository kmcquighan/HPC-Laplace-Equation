// by Kelly McQuighan, 8/2/13
/*
This code is a template for how to use the ScaLaPACK high performance parallel linear
algebra libraries. This code solves the 2D Laplace equation u_xx + u_yy = -2*sin(x)*sin(y)
with periodic boundary conditions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "blacs_headers.h"
#include "pblas_headers.h"
#include "scalapack_tools_headers.h"
#include "scalapack_headers.h"
#include "scalapackinfo.h"
#include "LaplaceMatinit.h"
#include "pblasIOtools.h"
#include "stringTools.h"
#include "nameFile.h"

static int imax( int a, int b ){
  if (a>b) return(a); else return(b);
}

//------MAIN-----------------
int main(int argc, char *argv[])
{
  int i_am;
  int m, n, mn, nrhs, mb, nb, nb_rhs, nprow, npcol, num_procs, status, iamIOproc = 0;
  int descA[9], descB[9], info, izero=0, ione=1, ctxt, myrow, mycol;
  int maximn_rows, imn_rows, imn_cols, nrhsloc, maximn_cols, maxnrhsloc, maxall;
  int *ipiv;
  double eps, Anorm, Bnorm, xnorm, resid, alpha, beta; 
  char folder[]="data", suffix_out[]="out", infile_name[] = "data/ScalapackMatrixInfo.dat";
  FILE *outfile=NULL;
  char outfile_name[100], residStr[33], usr_info[100], file_prefix[100];
  double *A=NULL, *b=NULL, *work, *A0=NULL, *b0=NULL, *b1=NULL;

//****************************************************
//                INITIALIZE MPI & BLACS
//****************************************************

  Cblacs_pinfo( &i_am, &num_procs );
  if (i_am == 0) iamIOproc = 1;
  // read data from file BLAS.dat. Exit program if file not formatted properly
  scalapackinfo( &m, &n, &nrhs, &mb, &nb, &nb_rhs, &nprow, &npcol, iamIOproc, &num_procs, infile_name, file_prefix, usr_info, &status );
  if ( status < 0 ) { 
    return EXIT_FAILURE;
  }


  if (iamIOproc) {    
    // create all file names based on the prefix
    nameFile ( folder, file_prefix, "", suffix_out, outfile_name );
   
    outfile = fopen ( outfile_name, "w" );

    // output data
    fprintf ( outfile, "Solution for Laplace's Equation for ScaLAPACK example practice.\n" );
    fprintf ( outfile, "By Kelly McQuighan.\n" );
    fprintf ( outfile, "Solving Delta u = -2*sin(x)*sin(y) on [0,2pi]x[0,2pi] with boundary conditions\n");
    fprintf ( outfile, "- u(x,2pi) = u(x,0),\n- u(2pi,y) = u(0,y), and\n- u(0,0) = 0, and where\nthe domain is discretized into a ");
    fprintf ( outfile, "%dx%d grid with block size of %dx%d and \n", m, n, mb, nb);
    fprintf ( outfile, "f is a %dx%d vector with block size of %dx%d.\n", n, nrhs, nb, nb_rhs );
    fprintf ( outfile, "Running on %d processes, where the process grid is %dx%d.\n", nprow*npcol, nprow, npcol );
    fprintf ( outfile, "****************************************************************************************************\n\n");
  }


  if (nprow*npcol > num_procs){
    if (i_am==0) {
      fprintf ( stdout, "\n****************************************************************************************************\n\n");
      fprintf ( stdout, "ERROR : we do not have enough processes available to make a %d x %d process grid\n\n", nprow, npcol );
      fprintf ( stdout, "****************************************************************************************************\n");
    }

    return EXIT_FAILURE;
  }
  
  if (m != n) {
    if (i_am==0){
      fprintf ( stdout, "\n****************************************************************************************************\n\n");
      fprintf ( stdout, "ERROR : this code only supports M=N. Setting N=M\n\n" );
      fprintf ( stdout, "****************************************************************************************************\n");
    }	
  	
  	n=m;
  }
  
  if (mb != m) {
    if (i_am==0){
      fprintf ( stdout, "\n****************************************************************************************************\n\n");
      fprintf ( stdout, "ERROR : this code only supports MB=M. Setting MB=M\n\n" );
      fprintf ( stdout, "****************************************************************************************************\n");
    }	
  	
  	mb=m;
  }
  
  if (mb != m) {
    if (i_am==0){
      fprintf ( stdout, "\n****************************************************************************************************\n\n");
      fprintf ( stdout, "ERROR : this code only supports MB=M. Setting MB=M\n\n" );
      fprintf ( stdout, "****************************************************************************************************\n");
    }	
  	
  	mb=m;
  }
  
  if (nb != mb) {
    if (i_am==0){
      fprintf ( stdout, "\n****************************************************************************************************\n\n");
      fprintf ( stdout, "ERROR : this code only supports MB=NB. Setting NB=MB\n\n" );
      fprintf ( stdout, "****************************************************************************************************\n");
    }	
  	
  	nb=mb;
  }
  
  mn = m*n;

  Cblacs_get( -1, 0, &ctxt );
  Cblacs_gridinit( &ctxt, "Row", nprow, npcol );
  Cblacs_gridinfo( ctxt, &nprow, &npcol, &myrow, &mycol );

// only continue if I'm on the grid
  if ( ( myrow < nprow ) & ( mycol < npcol ) ) {

//****************************************************
//                INITIALIZE MATRICES
//****************************************************
    // allocate memory for the local part of matrices A, B, and C
    imn_rows       = numroc_( &mn   , &mb    , &myrow, &izero, &nprow );
    imn_cols       = numroc_( &mn   , &nb    , &mycol, &izero, &npcol );
    nrhsloc        = numroc_( &nrhs, &nb_rhs, &mycol, &izero, &npcol );
    maximn_rows    = imax(1, imn_rows   );
    maximn_cols    = imax(1, imn_cols   );
    maxnrhsloc     = imax(1, nrhsloc);
    maxall         = imax(imax(maximn_rows, maximn_cols), maxnrhsloc);

    // initialize descriptors for the matrices A, B
    descinit_( descA, &mn, &mn   , &mb, &nb    , &izero, &izero, &ctxt, &imn_rows, &info );
    descinit_( descB, &mn, &nrhs, &mb, &nb_rhs, &izero, &izero, &ctxt, &imn_rows, &info );

    A = (double*) malloc(maximn_rows*maximn_cols   *sizeof(double));
    b = (double*) malloc(maximn_rows*maxnrhsloc*sizeof(double));
    work = (double*) malloc(maxall*maxall*sizeof(double));

    // generate matrices A and B
    LaplaceMatinit( A, descA, b, descB, myrow, mycol, nprow, npcol );
//    pdlaprnt_( &mn, &mn, A, &ione, &ione, descA, &izero, &izero, "A", &isix, work, 1 );

    // make a copy of matrices A and B for reference later
    A0 = (double*) malloc(maximn_rows*maximn_cols   *sizeof(double));
    b0 = (double*) malloc(maximn_rows*maxnrhsloc*sizeof(double));
    b1 = (double*) malloc(maximn_rows*maxnrhsloc*sizeof(double));
    pdlacpy_( "All", &mn, &mn, A, &ione, &ione, descA, A0, &ione, &ione, descA );
    pdlacpy_( "All", &mn, &nrhs, b, &ione, &ione, descB, b0, &ione, &ione, descB );
    pdlacpy_( "All", &mn, &nrhs, b, &ione, &ione, descB, b1, &ione, &ione, descB );

//***************************************************
//	      CALL THE SCALAPACK ROUTINE
//***************************************************
 
    ipiv = (int*) malloc((imn_rows+mb)*sizeof(int));
    pdgesv_( &mn, &nrhs, A, &ione, &ione, descA, ipiv, b, &ione, &ione, descB, &info );

    if (iamIOproc) fprintf( outfile, "INFO code returned by PSGESV = %d.\n\n", info );

    // compute residual for verification || AX-B || / ( ||X||*||A||*eps*N )
    eps = pdlamch_( &ctxt, "Epsilon" );
    Anorm = pdlange_ ( "I", &mn, &mn, A, &ione, &ione, descA, work );
    Bnorm = pdlange_ ( "I", &mn, &nrhs, b, &ione, &ione, descB, work );
    alpha = 1.0; beta = -1.0;
    pdgemm_ ( "N", "N", &mn, &nrhs, &mn, &alpha, A0, &ione, &ione, descA, b, &ione, &ione, descB, &beta, b0, &ione, &ione, descB );
    xnorm = pdlange_ ( "I", &mn, &nrhs, b0, &ione, &ione, descB, work );
    resid = xnorm / ( Anorm * Bnorm * eps * (double) n );
    dtoa_len( resid, residStr, "e", 4 );
    if (iamIOproc) {
      fprintf( outfile, "According to the normalized residual:\n");
      fprintf( outfile, "||A*u-b||/( ||u||*||A||*eps*N ) = %s,\n", residStr );
      if (resid <= 10.0e0 ) fprintf( outfile, "the solution is correct.\n\n");
      else                  fprintf( outfile, "the solution is incorrect.\n\n");
    }

    alpha = -0.5;
    pdaxpy_ (&mn, &alpha, b0, &ione, &ione, descB, &ione, b1, &ione, &ione, descB, &ione );
    Bnorm = pdlange_ ( "I", &mn, &nrhs, b1, &ione, &ione, descB, work );
    resid = xnorm / ( Anorm * Bnorm * eps * (double) n );
    dtoa_len( resid, residStr, "e", 4 );
    if (iamIOproc) {
      fprintf( outfile, "Additionally, we expect u(x) = sin(x)*sin(y). According to the norm:\n");
      fprintf( outfile, "||sin(x)*sin(y)-u(x,y)||/( ||u||*||sin(x)*sin(y)||*eps*N ) = %s,\nthe solution is ", residStr );
      if (resid <= 10.0e0 ) fprintf ( outfile, "correct.\n\n");
      else fprintf ( outfile, "incorrect.\n\n");
    }

//****************************************************
//                FINALIZE ALL
//****************************************************

    free(A);
    free(b);
    free(work);
    free(A0);
    free(b0);
    free(ipiv);

    Cblacs_gridexit( ctxt );    
  }

  if (iamIOproc) {
    fprintf ( outfile, "****************************************************************************************************\n" );
    fprintf ( outfile, "End of Laplace's Equation Solver.");
    fclose  ( outfile );
  }

  Cblacs_exit( 0 );
  return 0;
}
