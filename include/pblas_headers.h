//
//******************Variable descriptions*****************************
//
//  DESCX   (global and local input) INTEGER array
//          On entry, DESCX  is an integer array of dimension DLEN_. This
//          is the array descriptor for the matrix X.
//
//  DIAG    (global input) CHARACTER*1
//          On entry,  DIAG  specifies  whether or not  sub( A )  is unit
//          triangular as follows:
//
//             DIAG = 'U' or 'u'  sub( A )  is  assumed to be unit trian-
//                                gular,
//
//             DIAG = 'N' or 'n'  sub( A ) is not assumed to be unit tri-
//                                angular.
//
//  INCX    (global input) INTEGER
//          On entry,  INCX   specifies  the  global  increment  for  the
//          elements of  X.  Only two values of  INCX   are  supported in
//          this version, namely 1 and M_X. INCX  must not be zero.
//
//  IX      (global input) INTEGER
//          On entry, IX  specifies X's global row index, which points to
//          the beginning of the submatrix sub( X ).
//
//  JX      (global input) INTEGER
//          On entry, JX  specifies X's global column index, which points
//          to the beginning of the submatrix sub( X ).
//  N       (global input) INTEGER
//          On entry,  N  specifies the  length of the  subvectors to  be
//          swapped. N must be at least zero.
//
// TRANS   (global input) CHARACTER*1
//          On entry,  TRANS  specifies whether or not to transpose matrix A
//	    Options: 'N'=no transpose, 
//		     'T'=transpose, 
//		     'C'=conjugate transpose. 
//	    See individual function descriptions for the effect of each option.
//
// UPLO    (global input) CHARACTER*1
//          On  entry,   UPLO  specifies  whether  the  local  pieces  of
//          the array  A  containing the  upper or lower triangular  part
//          of the symmetric submatrix  sub( A )  are to be referenced as
//          follows:
//
//             UPLO = 'U' or 'u'   Only the local pieces corresponding to
//                                 the   upper  triangular  part  of  the
//                                 symmetric submatrix sub( A ) are to be
//                                 referenced,
//
//             UPLO = 'L' or 'l'   Only the local pieces corresponding to
//                                 the   lower  triangular  part  of  the
//                                 symmetric submatrix sub( A ) are to be
//                                 referenced.
//  X       (local input/local output) array
//          On entry, X is an array of dimension (LLD_X, Kx), where LLD_X
//          is   at  least  MAX( 1, Lr( 1, IX ) )  when  INCX = M_X   and
//          MAX( 1, Lr( 1, IX+N-1 ) )  otherwise,  and,  Kx  is  at least
//          Lc( 1, JX+N-1 )  when  INCX = M_X  and Lc( 1, JX ) otherwise.
//          Before  entry,  this array  contains the local entries of the
//          matrix X. On exit, sub( X ) is overwritten with sub( Y ).
//
//  Y       (local input/local output) array
//          On entry, Y is an array of dimension (LLD_Y, Ky), where LLD_Y
//          is   at  least  MAX( 1, Lr( 1, IY ) )  when  INCY = M_Y   and
//          MAX( 1, Lr( 1, IY+N-1 ) )  otherwise,  and,  Ky  is  at least
//          Lc( 1, JY+N-1 )  when  INCY = M_Y  and Lc( 1, JY ) otherwise.
//          Before  entry,  this array  contains the local entries of the
//          matrix Y. On exit, sub( Y ) is overwritten with sub( X ).
//
// -----------level 3 variables:----------------
//
//  K       (global input) INTEGER
//          On entry, K specifies the number of columns of the  submatrix
//          op( sub( A ) )  and  the  number of rows   of  the  submatrix
//          op( sub( B ) ). K must be at least  zero.
//
//  M       (global input) INTEGER
//          On entry,  M  specifies  the number of rows of the  submatrix
//          op( sub( A ) ) and of the submatrix sub( C ). M  must  be  at
//          least  zero.
//
//  N       (global input) INTEGER
//          On entry, N specifies the number of columns of the  submatrix
//          op( sub( B ) )  and  the  number of columns of the  submatrix
//          sub( C ). N must be at least zero.
//
// SIDE    (global input) CHARACTER*1
//          On entry,  SIDE  specifies  whether the  symmetric  submatrix
//          sub( A )  appears  on  the left or right in the operation  as
//          follows:
//
//             SIDE = 'L' or 'l'
//                   sub( C ) := alpha*sub( A )*sub( B ) + beta*sub( C ),
//
//             SIDE = 'R' or 'r'
//                   sub( C ) := alpha*sub( B )*sub( A ) + beta*sub( C ).
//
// ********************************************************************
//
// *******************LEVEL 1 ROUTINES*****************************
//
// ---------------single precision, real-------------------------
extern void psswap_(int* N, float* X, int* IX,	int* JX, int* DESCX, int* INCX,	float* Y, int* IY, int*	JY, int* DESCY,	int* INCY);
extern void psscal_(int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pscopy_(int* N, float* X, int* IX,	int* JX, int* DESCX, int* INCX,	float* Y, int* IY, int*	JY, int* DESCY,	int* INCY);
extern void psaxpy_(int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, int* INCX, float* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void psdot_(int* N, float* DOT, float* X, int* IX, int* JX, int* DESCX, int* INCX, float* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void psnrm2_(int* N, float* NORM2, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void psasum_(int* N, float* ASUM, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void psamax_(int* N, float* AMAX, int* INDX, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void psgeadd_( char* trans, int* m, int* n, float* alpha, float* A, int* ia, int* ja, int* descA, float* beta, float *C, int *ic, int *jc, int *descc);
//
// ---------------------double precision, real------------------
extern void pdswap_(int* N, double* X, int* IX,	int* JX, int* DESCX, int* INCX,	double* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pdscal_(int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pdcopy_(int* N, double* X, int* IX,	int* JX, int* DESCX, int* INCX,	double* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pdaxpy_(int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, int* INCX,	double* Y, int* IY, int* JY, int* DESCY,	int* INCY);
extern void pddot_(int* N, double* DOT, double* X, int* IX, int* JX, int* DESCX, int* INCX, double* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pdnrm2_(int* N, double* NORM2, double* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pdasum_(int* N, double* ASUM, double* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pdamax_(int* N, double* AMAX, int* INDX, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pdgeadd_( char* trans, int* m, int* n, double* alpha, double* A, int* ia, int* ja, int* descA, double* beta, double *C, int *ic, int *jc, int *descc);
//
// -------------------single precision, complex----------------------
extern void pcswap_(int* N, float* X, int* IX,	int* JX, int* DESCX, int* INCX,	float* Y, int* IY, int*	JY, int* DESCY,	int* INCY);
extern void pcscal_(int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pccopy_(int* N, float* X, int* IX,	int* JX, int* DESCX, int* INCX,	float* Y, int* IY, int*	JY, int* DESCY,	int* INCY);
extern void pcaxpy_(int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, int* INCX,	float* Y, int* IY, int*	JY, int* DESCY,	int* INCY);
extern void pcdotu_(int* N, float* DOTU, float* X, int* IX, int* JX, int* DESCX, int* INCX, float* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pcdotc_(int* N, float* DOTC, float* X, int* IX, int* JX, int* DESCX, int* INCX, float* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pscnrm2_(int* N, float* NORM2, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pscasum_(int* N, float* ASUM, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pcamax_(int* N, float* AMAX, int* INDX, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pcgeadd_( char* trans, int* m, int* n, float* alpha, float* A, int* ia, int* ja, int* descA, float* beta, float *C, int *ic, int *jc, int *descc);
//
// --------------------double precision, complex----------------------------
extern void pzswap_(int* N, double* X, int* IX,	int* JX, int* DESCX, int* INCX,	double* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pzscal_(int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pzcopy_(int* N, double* X, int* IX,	int* JX, int* DESCX, int* INCX,	double* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pzaxpy_(int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, int* INCX,	double* Y, int* IY, int* JY, int* DESCY,	int* INCY);
extern void pzdotu_(int* N, double* DOTU, double* X, int* IX, int* JX, int* DESCX, int* INCX, double* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pzdotc_(int* N, double* DOTC, double* X, int* IX, int* JX, int* DESCX, int* INCX, double* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pdznrm2_(int* N, double* NORM2, double* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pdzasum_(int* N, double* ASUM, double* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pzamax_(int* N, double* AMAX, int* INDX, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pzgeadd_( char* trans, int* m, int* n, double* alpha, double* A, int* ia, int* ja, int* descA, double* beta, double *C, int *ic, int *jc, int *descc);
//
//------------------------mixed type--------------------------------------
extern void pcsscal_(int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pdzscal_(int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pzlaprnt_(int* M, int* N, double* A, int* ia, int* ja, int* descA, int* irprnt, int* icprnt, char* matrixname, int* outfile, double* WORK, int len_m);
extern void pzlawrite_(char* filename, int* m, int* n, double* A, int* ia, int* ja, int* descA, int* irwrit, int* icwrit, double* work, int name_len );
extern void pzgeadd_( char* trans, int* m, int* n, double* alpha, double* A, int* ia, int* ja, int* descA, double* beta, double *C, int *ic, int *jc, int *descc);
//
// ************************LEVEL 2 ROUTINES***************************
//
// -------------------single precision, real------------------------------
extern void psgemv_(char* TRANS, int* M, int* N, float* ALPHA, float* A, int* IA, int* JA, int* DESCA, float* X, int* IX, int* JX, int* DESCX, int* INCX, float* BETA, float* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pssymv_( char* UPLO, int* N, float* ALPHA, float* A, int* IA, int* JA, int* DESCA, float* X, int* IX, int* JX, int* DESCX, int* INCX, float* BETA, float* Y, int* IY, int* JY, int* DESCY, int* INCY); 	
extern void pstrmv_( char* UPLO, char* TRANS, char* DIAG, int* N, float* A, int* IA, int* JA, int* DESCA, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pstrsv_( char* UPLO, char* TRANS, char* DIAG, int* N, float* A, int* IA, int* JA, int* DESCA, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void psger_(int* M, int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, float* Y, int* IY, int* JY, int* DESCY, int* INCY, float* A, int* IA, int* JA, int* DESCA);
extern void pssyr_(char* UPLO, int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, float* A, int* IA, int* JA, int* DESCA);
extern void pssyr2_(char* UPLO, int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, float* Y, int* IY, int* JY, int* DESCY, int* INCY, float* A, int* IA, int* JA, int* DESCA);
//
// -----------------------double precision, real---------------------------------
extern void pdgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* IA, int* JA, int* DESCA, double* X, int* IX, int* JX, int* DESCX, int* INCX, double* BETA, double* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pdsymv_( char* UPLO, int* N, double* ALPHA, double* A, int* IA, int* JA, int* DESCA, double* X, int* IX, int* JX, int* DESCX, int* INCX, double* BETA, double* Y, int* IY, int* JY, int* DESCY, int* INCY); 	
extern void pdtrmv_( char* UPLO, char* TRANS, char* DIAG, int* N, double* A, int* IA, int* JA, int* DESCA, double* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pdtrsv_( char* UPLO, char* TRANS, char* DIAG, int* N, double* A, int* IA, int* JA, int* DESCA, double* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pdger_(int* M, int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, double* Y, int* IY, int* JY, int* DESCY, int* INCY, double* A, int* IA, int* JA, int* DESCA);
extern void pdsyr_(char* UPLO, int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, double* A, int* IA, int* JA, int* DESCA);
extern void pdsyr2_(char* UPLO, int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, double* Y, int* IY, int* JY, int* DESCY, int* INCY, double* A, int* IA, int* JA, int* DESCA);
//
// ------------------------single precision, complex-----------------------------
extern void pcgemv_(char* TRANS, int* M, int* N, float* ALPHA, float* A, int* IA, int* JA, int* DESCA, float* X, int* IX, int* JX, int* DESCX, int* INCX, float* BETA, float* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pchemv_( char* UPLO, int* N, float* ALPHA, float* A, int* IA, int* JA, int* DESCA, float* X, int* IX, int* JX, int* DESCX, int* INCX, float* BETA, float* Y, int* IY, int* JY, int* DESCY, int* INCY); 	
extern void pctrmv_( char* UPLO, char* TRANS, char* DIAG, int* N, float* A, int* IA, int* JA, int* DESCA, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pctrsv_( char* UPLO, char* TRANS, char* DIAG, int* N, float* A, int* IA, int* JA, int* DESCA, float* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pcgeru_(int* M, int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, float* Y, int* IY, int* JY, int* DESCY, int* INCY, float* A, int* IA, int* JA, int* DESCA);
extern void pcgerc_(int* M, int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, float* Y, int* IY, int* JY, int* DESCY, int* INCY, float* A, int* IA, int* JA, int* DESCA);
extern void pcher_(char* UPLO, int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, float* A, int* IA, int* JA, int* DESCA);
extern void pcher2_(char* UPLO, int* N, float* alpha, float* X, int* IX, int* JX, int* DESCX, float* Y, int* IY, int* JY, int* DESCY, int* INCY, float* A, int* IA, int* JA, int* DESCA);
//
// -----------------------double precision, complex-------------------------------
extern void pzgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* IA, int* JA, int* DESCA, double* X, int* IX, int* JX, int* DESCX, int* INCX, double* BETA, double* Y, int* IY, int* JY, int* DESCY, int* INCY);
extern void pzhemv_( char* UPLO, int* N, double* ALPHA, double* A, int* IA, int* JA, int* DESCA, double* X, int* IX, int* JX, int* DESCX, int* INCX, double* BETA, double* Y, int* IY, int* JY, int* DESCY, int* INCY); 	
extern void pztrmv_( char* UPLO, char* TRANS, char* DIAG, int* N, double* A, int* IA, int* JA, int* DESCA, double* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pztrsv_( char* UPLO, char* TRANS, char* DIAG, int* N, double* A, int* IA, int* JA, int* DESCA, double* X, int* IX, int* JX, int* DESCX, int* INCX);
extern void pzgeru_(int* M, int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, double* Y, int* IY, int* JY, int* DESCY, int* INCY, double* A, int* IA, int* JA, int* DESCA);
extern void pzgerc_(int* M, int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, double* Y, int* IY, int* JY, int* DESCY, int* INCY, double* A, int* IA, int* JA, int* DESCA);
extern void pzher_(char* UPLO, int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, double* A, int* IA, int* JA, int* DESCA);
extern void pzher2_(char* UPLO, int* N, double* alpha, double* X, int* IX, int* JX, int* DESCX, double* Y, int* IY, int* JY, int* DESCY, int* INCY, double* A, int* IA, int* JA, int* DESCA);
//
// *********************LEVEL 3 ROUTINES*************************************
//
// ---------------------single precision, real-------------------------------------
extern void psgemm_( char* TRANSA, char* TRANSB, int* M, int* N, int* K, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pssymm_( char* SIDE, char* UPLO, int* M, int* N, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pssyrk_( char* UPLO, char* TRANS, int* N, int* K, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pssyr2k_( char* UPLO, char* TRANS, int* N, int* K, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pstran_( int* M, int* N, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pstrmm_( char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB);
extern void pstrsm_( char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB);
//
//----------------------double precision, real------------------------------------
extern void pdgemm_( char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB, double* BETA, double* C, int* IC, int* JC, int* DESCC);
extern void pdsymm_( char* SIDE, char* UPLO, int* M, int* N, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB, double* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pdsyrk_( char* UPLO, char* TRANS, int* N, int* K, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* BETA, double* C, int* IC, int* JC, int* DESCC);
extern void pdsyr2k_( char* UPLO, char* TRANS, int* N, int* K, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB, double* BETA, double* C, int* IC, int* JC, int* DESCC);
extern void pdtran_( int* M, int* N, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* BETA, double* C, int* IC, int* JC, int* DESCC);
extern void pdtrmm_( char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB);
extern void pdtrsm_( char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB);
//
// -------------------single precision, complex----------------------------------
extern void pcgemm_( char* TRANSA, char* TRANSB, int* M, int* N, int* K, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pcsymm_( char* SIDE, char* UPLO, int* M, int* N, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pchemm_( char* SIDE, char* UPLO, int* M, int* N, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pcsyrk_( char* UPLO, char* TRANS, int* N, int* K, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pcherk_( char* UPLO, char* TRANS, int* N, int* K, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pcher2k_( char* UPLO, char* TRANS, int* N, int* K, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pcsyr2k_( char* UPLO, char* TRANS, int* N, int* K, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pctranu_( int* M, int* N, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pctranc_( int* M, int* N, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pctrmm_( char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB);
extern void pctrsm_( char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, float* alpha, float* A, int* IA, int* JA, int* DESCA, float* B, int* IB, int* JB, int* DESCB);
//
// --------------------------double precision, complex-------------------------
extern void pzgemm_( char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB, double* BETA, double* C, int* IC, int* JC, int* DESCC);
extern void pzsymm_( char* SIDE, char* UPLO, int* M, int* N, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB, double* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pzhemm_( char* SIDE, char* UPLO, int* M, int* N, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB, double* BETA, float* C, int* IC, int* JC, int* DESCC);
extern void pzsyrk_( char* UPLO, char* TRANS, int* N, int* K, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* BETA, double* C, int* IC, int* JC, int* DESCC);
extern void pzherk_( char* UPLO, char* TRANS, int* N, int* K, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* BETA, double* C, int* IC, int* JC, int* DESCC);
extern void pzsyr2k_( char* UPLO, char* TRANS, int* N, int* K, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB, double* BETA, double* C, int* IC, int* JC, int* DESCC);
extern void pzher2k_( char* UPLO, char* TRANS, int* N, int* K, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB, double* BETA, double* C, int* IC, int* JC, int* DESCC);
extern void pztranu_( int* M, int* N, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* BETA, double* C, int* IC, int* JC, int* DESCC);
extern void pztranc_( int* M, int* N, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* BETA, double* C, int* IC, int* JC, int* DESCC);
extern void pztrmm_( char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB);
extern void pztrsm_( char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* alpha, double* A, int* IA, int* JA, int* DESCA, double* B, int* IB, int* JB, int* DESCB);
