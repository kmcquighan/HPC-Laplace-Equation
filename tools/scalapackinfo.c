// -- PBLAS Example code --
// Original fortran code from http://acts.nersc.gov/scalapack/hands-on/
// Conversion into C by Kelly McQuighan,
// last modified 7/31/13
// 
// reads system parameters from the file BLAS.dat

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <blacs_headers.h>
#include <pblas_headers.h>
#include "scalapackinfo.h"
#include <string.h>

// returns -1 if unsuccessful, 0 if successful
void scalapackinfo( int* m, int* n, int* nrhs, int* mb, int* nb, int* nb_rhs, int* nprow, int* npcol, int iamIOproc, int* num_procs, char* infile_name, char* file_prefix, char* usr_info, int* status ){

//  char usr_info[100];
  char line[100], info[100], trash[100];
  int contxt, sys_param[10], i;
  FILE* infile;

  Cblacs_get( -1, 0, &contxt );
  Cblacs_gridinit( &contxt, "Row-major", 1, *num_procs );

// processor 0 reads data from PBLAS.dat and broadcasts to other processors.
  if ( iamIOproc ) {
    infile = fopen(infile_name, "r");
    if (infile == NULL) {
      fprintf(stdout, "***********************************************************\n");
      fprintf(stdout, "Error: input file %s expected!\n", infile_name);
      fprintf(stdout, "***********************************************************\n");
      *status = -1;
      Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
    }
// header line
    if (fgets(line, 100, infile)==NULL ) {
      fprintf(stdout, "\n***********************************************************\n\n");
      fprintf(stdout, "Error: input file does not contain any data!\n");
      fprintf(stdout, "\n***********************************************************\n\n");
      *status = -1;
      Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
    } 
// user info
    if (fgets(usr_info, 100, infile)==NULL ) {
      fprintf(stdout, "\n***********************************************************\n\n");
      fprintf(stdout, "Error: input file does not contain any data!\n");
      fprintf(stdout, "\n***********************************************************\n\n");
      *status = -1;
      Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
    }
// output file name
    if (fgets(line, 100, infile)==NULL )  {
      fprintf(stdout, "\n***********************************************************\n\n");
      fprintf(stdout, "Error: input file should contain the output file name on the 3rd line!\n");
      fprintf(stdout, "\n***********************************************************\n\n");
      *status = -1;
      Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
    }
    if ( sscanf(line,"%[^\t]%[^\n]", file_prefix, trash) < 1 ) {
      fprintf(stdout, "\n***********************************************************\n\n");
      fprintf(stdout, "Error: input file not formatted correctly! 3rd line should read (output file name) (tab) (additional text).\n");
      fprintf(stdout, "\n***********************************************************\n\n");
      *status = -1;
      Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
    }

//  system parameters
    for (i=0; i< 8; i++ ) {

      if (fgets(line, 100, infile) == NULL ) {
        fprintf(stdout, "\n***********************************************************\n\n");
        fprintf(stdout, "Error: input file not formatted correctly! Lines 4-13 should contain Scalapack matrix information in the format\n\n(value of m) (tab) (additional text)\n");
        fprintf(stdout, "in the following order: m, n, nrhs, mb, nb, nb_rhs, max_lldA, max_lldB, nprow, npcol\n");
        fprintf(stdout, "\n***********************************************************\n\n");
        *status = -1;
        Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
      return;
      }
      if ( sscanf(line,"%[^\t]%[^\n]", info, trash) < 2 ) {
        fprintf(stdout, "\n***********************************************************\n\n");
        fprintf(stdout, "Error: input file not formatted correctly on %dth line! It should read\n\n(value) (tab) (additional text)\n", i+4);
        fprintf(stdout, "\n***********************************************************\n\n");
        *status = -1;
        Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
        return;
      }
      sys_param[i] = atoi(info);
    }

// completed successfully. broadcast parameters to other processes
    *status = 0;
    Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );  
    Cigebs2d( contxt, "All", " ", 8, 1, sys_param, 8);

  } else {

// receive system variables and store in corresponding variables
    Cigebs2d( contxt, "All", " ", 1, 1, status, 1 );
    if ( *status < 0 ) return;
    Cigebs2d( contxt, "All", " ", 8, 1, sys_param, 8 );
  } // end if-then-else

// store system parameters in their appropriate values
  *m        = sys_param[0];
  *n        = sys_param[1];
  *nrhs     = sys_param[2];
  *mb       = sys_param[3];
  *nb       = sys_param[4];
  *nb_rhs   = sys_param[5];
  *nprow    = sys_param[6];
  *npcol    = sys_param[7];

  Cblacs_gridexit( contxt );

  return;
}
