// *************VARIABLE DESCRIPTIONS****************************
//
//  CMACH   (global input) CHARACTER*1
//          Specifies the value to be returned by PDLAMCH:
//          = 'E' or 'e',   PDLAMCH := eps
//          = 'S' or 's ,   PDLAMCH := sfmin
//          = 'B' or 'b',   PDLAMCH := base
//          = 'P' or 'p',   PDLAMCH := eps*base
//          = 'N' or 'n',   PDLAMCH := t
//          = 'R' or 'r',   PDLAMCH := rnd
//          = 'M' or 'm',   PDLAMCH := emin
//          = 'U' or 'u',   PDLAMCH := rmin
//          = 'L' or 'l',   PDLAMCH := emax
//          = 'O' or 'o',   PDLAMCH := rmax
//
//          where
//
//          eps   = relative machine precision
//          sfmin = safe minimum, such that 1/sfmin does not overflow
//          base  = base of the machine
//          prec  = eps*base
//          t     = number of (base) digits in the mantissa
//          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
//          emin  = minimum exponent before (gradual) underflow
//          rmin  = underflow threshold - base**(emin-1)
//          emax  = largest exponent before overflow
//          rmax  = overflow threshold  - (base**emax)*(1-eps)
//
// DESC    (output) INTEGER array of dimension DLEN_.
// 		The array descriptor of a distributed matrix to be set.
//  ICSRC   (global input) INTEGER
//	          The process column over which the first column of the
//	          matrix is distributed. 0 <= ICSRC < NPCOL.
//
//  ICTXT   (global input) INTEGER
//	          The BLACS context handle, indicating the global context of
//	          the operation on the matrix. The context itself is global.
//
//  INDXGLOB  (global input) INTEGER
//            The global index of the element.
//
//	          = 0: successful exit
//	         < 0: if INFO = -i, the i-th argument had an illegal value
//
// IPROC     (local input) INTEGER
//	           The coordinate of the process whose local array row or
//	           column is to be determined.
//
//  IRSRC   (global input) INTEGER
//	          The process row over which the first row of the matrix is
//	          distributed. 0 <= IRSRC < NPROW.
//
//  ISRCPROC  (global input) INTEGER
//	            The coordinate of the process that possesses the first
//	            row or column of the distributed matrix.
//
//  LLD     (local input)  INTEGER
//  	          The leading dimension of the local array storing the local
//	          blocks of the distributed matrix. LLD >= MAX(1,LOCr(M)).
//
// M       (global input) INTEGER
//		The number of rows in the distributed matrix. M >= 0.
//
// MB      (global input) INTEGER
//	         The blocking factor used to distribute the rows of the
//	         matrix. MB >= 1.
//
// N       (global input) INTEGER
//		The number of columns in the distributed matrix. N >= 0.
//
// NB      (global input) INTEGER
//	          The blocking factor used to distribute the columns of the
//	          matrix. NB >= 1.
//
// NORM    (global input) CHARACTER
//         Specifies the value to be returned in PDLANGE as described
//          above.
//
//  NPROCS    (global input) INTEGER
//	            The total number processes over which the matrix is
//	            distributed.
//
//  WORK      (local workspace) DOUBLE PRECISION. See documentation for
//	      individual functions for the meaning of the return values
//	      of work
//
//****************************************************************
//
// numroc 	: computes the numberof rows or columns of a distributed
//		  matrix owned by the process indicated by IPROC 
//
// descinit 	: initializes the descriptor array for use in pblas and
//                ScaLAPACK routines
//
// pdelset 	: sets the distributed matrix entry A( ia, ja ) to alpha
//
// pdlamch	: determines double precision machine parameters
//
// indxg2p	: computes the process coordinate which possesses the entry
//		  of a distributed matrix specified by global index indxglob
//
// indxg2l	: computes the local index of a distributed matrix entry
//		  pointed to by the global index INDXGLOB.
//
// pdlaset	: initializes an M-by-N distributed matrix sub( A ) denoting
//		  A(IA:IA+M-1,JA:JA+N-1) to BETA on the diagonal and ALPHA on the
//		  offdiagonals.
//
// pdlange	:  PDLANGE returns the value of the one norm, or the Frobenius norm,
//		  or the infinity norm, or the element of largest absolute value of a
//		  distributed matrix sub( A ) = A(IA:IA+M-1, JA:JA+N-1).
//
//		  PDLANGE returns the value
//
//		  max(abs(A(i,j))),  NORM = 'M' or 'm' with IA <= i <= IA+M-1,
//                                           and  JA <= j <= JA+N-1,
//    
//	         norm1( sub( A ) ), NORM = '1', 'O' or 'o'
//
//		 normI( sub( A ) ), NORM = 'I' or 'i'
//
//		 normF( sub( A ) ), NORM = 'F', 'f', 'E' or 'e'
//
//		where norm1 denotes the  one norm of a matrix (maximum column sum),
//		normI denotes the  infinity norm  of a matrix  (maximum row sum) and
//		normF denotes the  Frobenius norm of a matrix (square root of sum of
//		Note that  max(abs(A(i,j)))  is not a  matrix norm.
//
// pdlacpy	: copies all or part of a distributed matrix A to another
//		  distributed matrix B.  No communication is performed, PDLACPY
//		  performs a local copy sub( A ) := sub( B ), where sub( A ) denotes
//		  A(IA:IA+M-1,JA:JA+N-1) and sub( B ) denotes B(IB:IB+M-1,JB:JB+N-1).
//
// pdlaprnt	: prints distributed matrix A to the file POINTED TO by outfile. Note,
//		  this doesn't seem to work with outfiles of type FILE or MPI_File.
//		  This is why I wrote the utility pdlaprnt2. outfile=6 refers to stdout  
//
// pdlawrite	: writes distributed matrix to to file with NAME filename. Overwrites 
//		  an existing file
//
// pdlaread	: reads a distributed matrix from the file with NAME filename.
//*******************Function Descriptions**********************
//
// **********************************************************
extern int    numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void   descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);
extern int    indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern int    indxg2l_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
//
//------------------SINGLE PRECISION REAL--------------
extern void   pselset_( float *A, int *ia, int *ja, int *desca, float *alpha);
extern float  pslamch_( int *ictxt, char *cmach);
extern void   pslaset_( char *uplo, int *m, int *n, float *alpha, float *beta, float *A, int *ia, int *ja, int *descA );
extern float  pslange_( char *norm, int *m, int *n, float *A, int *ia, int *ja, int *desca, float *work);
extern void   pslaprnt_(int* M, int* N, float* A, int* ia, int* ja, int* descA, int* irprnt, int* icprnt, char* matrixname, int* outfile, float* WORK, int len_m);
extern void   pslawrite_(char* filename, int* m, int* n, float* A, int* ia, int* ja, int* descA, int* irwrit, int* icwrit, float* work, int name_len );
extern void   pslaread_(char* filename, float* A, int* descA, int* ia, int* ja, float* work, int name_len );
//extern void psmatgen_( int *ictxt, char *aform, char *diag, int *m, int *n, int *mb, int *nb, float *A, int *lda, int *iarow, int *iacol, int *iseed, int *iroff, int *irnum, int *icoff, int *icnum, int *myrow, int *mycol, int *nprow, int *npcol, int aform_len, int diag_len );
//
//------------DOUBLE PRECISION REAL------------------
extern void   pdelset_( double *A, int *ia, int *ja, int *desca, double *alpha);
extern double pdlamch_( int *ictxt, char *cmach);
extern void   pdlaset_( char *uplo, int *m, int *n, double *alpha, double *beta, double *A, int *ia, int *ja, int *descA );
extern double pdlange_( char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);
extern void   pdlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);
extern void   pdlaprnt_(int* M, int* N, double* A, int* ia, int* ja, int* descA, int* irprnt, int* icprnt, char* matrixname, int* outfile, double* WORK, int len_m);
extern void   pdlawrite_(char* filename, int* m, int* n, double* A, int* ia, int* ja, int* descA, int* irwrit, int* icwrit, double* work, int name_len );
extern void   pdlaread_(char* filename, double* A, int* descA, int* ia, int* ja, double* work, int name_len );
//extern void pdmatgen_( int *ictxt, char *aform, char *diag, int *m, int *n, int *mb, int *nb, double *A, int *lda, int *iarow, int *iacol, int *iseed, int *iroff, int *irnum, int *icoff, int *icnum, int *myrow, int *mycol, int *nprow, int *npcol, int aform_len, int diag_len );
//
//------------------SINGLE PRECISION COMPLEX--------------
extern void   pcelset_( float *A, int *ia, int *ja, int *desca, float *alpha);
extern float  pclamch_( int *ictxt, char *cmach);
extern void   pclaset_( char *uplo, int *m, int *n, float *alpha, float *beta, float *A, int *ia, int *ja, int *descA );
extern float  pclange_( char *norm, int *m, int *n, float *A, int *ia, int *ja, int *desca, float *work);
extern void   pclaprnt_(int* M, int* N, float* A, int* ia, int* ja, int* descA, int* iaprnt, int* icprnt, char* matrixname, int* outfile, float* WORK, int len_m);
extern void   pclawrite_(char* filename, int* m, int* n, float* A, int* ia, int* ja, int* descA, int* irwrit, int* icwrit, float* work, int name_len );
extern void   pclaread_(char* filename, float* A, int* descA, int* ia, int* ja, float* work, int name_len );
//
//------------DOUBLE PRECISION COMPLEX------------------
extern void   pzelset_( double *A, int *ia, int *ja, int *desca, double *alpha);
extern double pzlamch_( int *ictxt, char *cmach);
extern void   pzlaset_( char *uplo, int *m, int *n, double *alpha, double *beta, double *A, int *ia, int *ja, int *descA );
extern double pzlange_( char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);
extern void   pzlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);
extern void   pzlaprnt_(int* M, int* N, double* A, int* ia, int* ja, int* descA, int* iaprnt, int* icprnt, char* matrixname, int* outfile, double* WORK, int len_m);
extern void   pzlawrite_(char* filename, int* m, int* n, double* A, int* ia, int* ja, int* descA, int* irwrit, int* icwrit, double* work, int name_len );
extern void   pzlaread_(char* filename, double* A, int* descA, int* ia, int* ja, double* work, int name_len );


