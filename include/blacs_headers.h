//*************ABREVIATIONS ***********
//
//  A        DOUBLE PRECISION, dimension (LDA,N)
//            Pointer to the location into the local memory to an array of 
//	      dimension (LLD,N) containing the local pieces of the 
//	      distributed matrix A.
//
//  CA      Integer Array, dimension (RCFLAG,N)
//           Contains process column that the amx of each element
//           of A was found on: i.e., cA(1,2) contains the process
//           column that the max/min of A(1,2) was found on.
//           Values are left on process {rdest, cdest} only, others
//           may be modified, but not left with interesting data.
//           If rdest == -1, then result is left on all processes in scope.
//           If LDIA == -1, this array is not accessed, and need not exist.
//
// CONTEXT   INTEGER
//            This integer is passed to the BLACS to indicate a context.
//            A context is a universe where messages exist and do not
//            interact with other context's messages.  The context includes
//            the definition of a grid, and each process's coordinates in it.
//
//  CSRC      The process column over which the first column of the array A 
//	      is distributed.
//
//  DIAG    Ptr to char
//	     Specifies whether the matrix is unit diagonal or not.
//		 = 'U':      Matrix is unit diagonal, diagonal not 
// 			     communicated.
//		ELSE :      Matrix is not unit diagonal, diagonal is 
//			    communicated.
//
//  LDA      INTEGER LDA >= M
//            Leading Dimension of A.
//
//  M        INTEGER, M >= 0
//            Number of rows of the global matrix owned by this process.
//
//  N       INTEGER N >= 0
//            Number of columns of the global matrix owned by this process.
//
//  RA      Integer Array, dimension (RCFLAG, N)
//	     Contains process row that the amx of each element
//	     of A was found on: i.e., rA(1,2) contains the process
//	     row that the amx of A(1,2) was found on.
//	     Values are left on process {rdest, cdest} only, others
//	     may be modified, but not left with interesting data.
//	     If rdest == -1, then result is left on all processes in scope.
//	     If LDIA == -1, this array is not accessed, and need not exist.
//
// RCFLAG    Ptr to int
//          If (RCFLAG == -1), then the arrays RA and CA are not accessed.
//          ELSE leading dimension of the arrays RA and CA.  LDIA >= M.
//
//  RSRC      The process row over which the first row of the array A is 
//	      distributed.
//
//  SCOPE   CHARACTER*1
//          The BLACS scope of the operation.
//          If SCOPE = 'R', the operation is performed only in the process
//                          row containing A( IA, JA ),
//          If SCOPE = 'C', the operation is performed only in the process 
//                          column containing A( IA, JA ),
//          If SCOPE = 'A', the operation is performed in all the processes 
//                          of the grid,
//          otherwise the operation is performed only in the process 
//           containing A( IA, JA ).
//
//  TOP     CHARACTER*1
//          The topology to be used if broadcast is needed.
//
//  UPLO    Ptr to char
//	     Specifies the part of the matrix to be sent.
//		= 'U':      Upper trapezoidal part
//		ELSE :      Lower trapezoidal part
//
// ************************************************
//
// ***************INITIALIZATION********************
extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
extern void   Cblacs_setup( int* mypnum, int* numprocs);
extern void   Cblacs_get( int context, int request, int* value );
extern void   Clacs_set( int context, int what, int* value );
extern void   Cblacs_gridinit( int* context, const char* order, int np_row, int np_col );
extern void   Cblacs_gridmap( int* context, char* pmap, int ldpmap, int nprow, int npcol );
//
//******************DESTRUCTION*************************
extern void   Cblacs_freebuff( int context, int wait );
extern void   Cblacs_gridexit( int context );
extern void   Clacs_abort( int context, int errornum );
extern void   Cblacs_exit( int error_code );
//
//*****************INFORMATIONAL AND MISC.******************
extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col );
extern int    Cblacs_pnum( int context, int prow, int pcol );
extern void   Cblacs_pcoord( int context, int nodenum, int* prow, int* pcol );
extern void   Cblacs_barrier( int context, char *scope );
//
//****************NON-STANDARD**********************
extern void   Csetpvmtids( int ntasks, int* tids );
extern double Cdcputime00( );
extern double Cdwalltime00( );
extern int    Cksendid( int context, int rdest, int cdest );
extern int    Ckrecvid( int context, int rsrc, int csrc );
extern int    CKbsid( int context, char* scope );
extern int    Ckbrid( int context, char* scope, int rsrc, int csrc );
//
//************SENDING*************************
//
//---------------INT-------------------------------
extern void   Cigesd2d( int context, int m, int n, int* A, int lda, int rdest, int cdest );
extern void   Cigebs2d( int context, char *scope, char *top, int m, int n, int* A, int lda );
extern void   Citrsd2d( int context, char* uplo, char* diag, int m, int n, int* A, int lda, int rdest, int cdest );
extern void   Citrbs2d( int icontext, char* scope, char *top, char* uplo, char* diag, int m, int n, int* A, int lda );
//
//-----------SINGLE PRECISION REAL--------------------
extern void   Csgesd2d( int context, int m, int n, float* A, int lda, int rdest, int cdest );
extern void   Csgebs2d( int context, char* scope, char* top, int m, int n, float* A, int lda );
extern void   Cstrsd2d( int context, char* uplo, char* diag, int m, int n, float* A, int lda, int rdest, int cdest );
extern void   Cstrbs2d( int icontext, char* scope, char *top, char* uplo, char* diag, int m, int n, float* A, int lda );
//
//------------DOUBLE PRECISION REAL-------------------------
extern void   Cdgesd2d( int context, int m, int n, double* A, int lda, int rdest, int cdest );
extern void   Cdgebs2d( int context, char *scope, char *top, int m, int n, double* A, int lda );
extern void   Cdtrsd2d( int context, char* uplo, char* diag, int m, int n, double* A, int lda, int rdest, int cdest );
extern void   Cdtrbs2d( int icontext, char* scope, char *top, char* uplo, char* diag, int m, int n, double* A, int lda );
//
//-----------SINGLE PRECISION COMPLEX-------------------
extern void   Ccgesd2d( int context, int m, int n, float* A, int lda, int rdest, int cdest );
extern void   Ccgebs2d( int context, char *scope, char *top, int m, int n, float* A, int lda );
extern void   Cctrsd2d( int context, char* uplo, char* diag, int m, int n, float* A, int lda, int rdest, int cdest );
extern void   Cctrbs2d( int icontext, char* scope, char *top, char* uplo, char* diag, int m, int n, float* A, int lda );
//
//-----------------DOUBLE PRECISION COMPLEX-------------
extern void   Czgesd2d( int context, int m, int n, double* A, int lda, int rdest, int cdest );
extern void   Czgebs2d( int context, char *scope, char *top, int m, int n, double* A, int lda );
extern void   Cztrsd2d( int context, char* uplo, char* diag, int m, int n, double* A, int lda, int rdest, int cdest );
extern void   Cztrbs2d( int icontext, char* scope, char *top, char* uplo, char* diag, int m, int n, double* A, int lda );
//
//****************RECEIVING*************************
//
//-----------------INT-----------------------
extern void   Cigerv2d( int context, int m, int n, int* A, int lda, int rsrc, int csrc );
extern void   Cigebr2d( int context, char *scope, char *top, int m, int n, int* A, int lda, int rsrc, int csrc );
extern void   Citrrv2d( int context, char* uplo, char* diag, int m, int n, int* A, int lda, int rsrc, int csrc );
extern void   Citrbr2d( int icontext, char* scope, char *top, char* uplo, char* diag, int m, int n, int* A, int lda, int rsrc, int csrc );
//
//------------SINGLE PRECISION REAL----------
extern void   Csgerv2d( int context, int m, int n, float* A, int lda, int rsrc, int csrc );
extern void   Csgebr2d( int context, char *scope, char *top, int m, int n, float* A, int lda, int rsrc, int csrc );
extern void   Cstrrv2d( int context, char* uplo, char* diag, int m, int n, float* A, int lda, int rsrc, int csrc );
extern void   Cstrbr2d( int icontext, char* scope, char *top, char* uplo, char* diag, int m, int n, float* A, int lda, int rsrc, int csrc );
//
//----------DOUBLE PRECISION REAL---------------
extern void   Cdgerv2d( int context, int m, int n, double* A, int lda, int rsrc, int csrc );
extern void   Cdgebr2d( int context, char *scope, char *top, int m, int n, double* A, int lda, int rsrc, int csrc );
extern void   Cdtrrv2d( int context, char* uplo, char* diag, int m, int n, double* A, int lda, int rsrc, int csrc );
extern void   Cdtrbr2d( int icontext, char* scope, char *top, char* uplo, char* diag, int m, int n, double* A, int lda, int rsrc, int csrc );
//
//---------SINGLE PRECISION COMPLEX-----------------
extern void   Ccgerv2d( int context, int m, int n, float* A, int lda, int rsrc, int csrc );
extern void   Ccgebr2d( int context, char *scope, char *top, int m, int n, float* A, int lda, int rsrc, int csrc );
extern void   Cctrrv2d( int context, char* uplo, char* diag, int m, int n, float* A, int lda, int rsrc, int csrc );
extern void   Cctrbr2d( int icontext, char* scope, char *top, char* uplo, char* diag, int m, int n, float* A, int lda, int rsrc, int csrc );
//
//--------------DOUBLE PRECISION COMPLEX----------------
extern void   Czgerv2d( int context, int m, int n, double* A, int lda, int rsrc, int csrc );
extern void   Czgebr2d( int context, char *scope, char *top, int m, int n, double* A, int lda, int rsrc, int csrc );
extern void   Cztrrv2d( int context, char* uplo, char* diag, int m, int n, double* A, int lda, int rsrc, int csrc );
extern void   Cztrbr2d( int icontext, char* scope, char *top, char* uplo, char* diag, int m, int n, double* A, int lda, int rsrc, int csrc );
//
//**********************COMBINE********************
//
//-------------INT------------------------------------
extern void   Cigamx2d( int contxt, char* scope, char* top, int m, int n, int* A, int lda, int* RA, int* CA, int RCflag, int rdest, int cdest );
extern void   Cigamn2d( int contxt, char* scope, char* top, int m, int n, int* A, int lda, int* RA, int* CA, int RCflag, int rdest, int cdest );
extern void   Cigsum2d( int contxt, char* scope, char* top, int m, int n, int* A, int lda, int rdest, int cdest );
//
//----------SINGLE PRECISION REAL-----------------
extern void   Csgamx2d( int contxt, char* scope, char* top, int m, int n, float* A, int lda, int* RA, int* CA, int RCflag, int rdest, int cdest );
extern void   Csgamn2d( int contxt, char* scope, char* top, int m, int n, float* A, int lda, int* RA, int* CA, int RCflag, int rdest, int cdest );
extern void   Csgsum2d( int contxt, char* scope, char* top, int m, int n, float* A, int lda, int rdest, int cdest );
//
//------------DOUBLE PRECISION REAL-----------------
extern void   Cdgamx2d( int contxt, char* scope, char* top, int m, int n, double* A, int lda, int* RA, int* CA, int RCflag, int rdest, int cdest );
extern void   Cdgamn2d( int contxt, char* scope, char* top, int m, int n, double* A, int lda, int* RA, int* CA, int RCflag, int rdest, int cdest );
extern void   Cdgsum2d( int contxt, char* scope, char* top, int m, int n, double* A, int lda, int rdest, int cdest );
//
//---------SINGLE PRECISION COMPLEX----------------
extern void   Ccgamx2d( int contxt, char* scope, char* top, int m, int n, float* A, int lda, int* RA, int* CA, int RCflag, int rdest, int cdest );
extern void   Ccgamn2d( int contxt, char* scope, char* top, int m, int n, float* A, int lda, int* RA, int* CA, int RCflag, int rdest, int cdest );
extern void   Ccgsum2d( int contxt, char* scope, char* top, int m, int n, float* A, int lda, int rdest, int cdest );
//
//--------------DOUBLE PRECISION COMPLEX------------------
extern void   Czgamx2d( int contxt, char* scope, char* top, int m, int n, double* A, int lda, int* RA, int* CA, int RCflag, int rdest, int cdest );
extern void   Czgamn2d( int contxt, char* scope, char* top, int m, int n, double* A, int lda, int* RA, int* CA, int RCflag, int rdest, int cdest );
extern void   Czgsum2d( int contxt, char* scope, char* top, int m, int n, double* A, int lda, int rdest, int cdest );

