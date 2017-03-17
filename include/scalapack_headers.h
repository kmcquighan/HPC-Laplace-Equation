//******************** VARIABLE DESCRIPTIONS*************
//
// also see pblas_headers.h for variable descriptions
//
//  A       (local input/local output) REAL pointer into
//          local memory to an array with first dimension
//          LLD_A >=(bw+1) (stored in DESCA).
//          On entry, this array contains the local pieces of the
//          This local portion is stored in the packed banded format
//            used in LAPACK. Please see the Notes below and the
//            ScaLAPACK manual for more detail on the format of
//            distributed matrices.
//          On exit, this array contains information containing details
//            of the factorization.
//          Note that permutations are performed on the matrix, so that
//            the factors returned are different from those returned
//            by LAPACK.
//
//          On exit, if EQUED .ne. 'N', A(IA:IA+N-1,JA:JA+N-1) is scaled
//          as follows:
//          EQUED = 'R':  A(IA:IA+N-1,JA:JA+N-1) :=
//                                      diag(R) * A(IA:IA+N-1,JA:JA+N-1)
//          EQUED = 'C':  A(IA:IA+N-1,JA:JA+N-1) :=
//                                      A(IA:IA+N-1,JA:JA+N-1) * diag(C)
//          EQUED = 'B':  A(IA:IA+N-1,JA:JA+N-1) :=
//                            diag(R) * A(IA:IA+N-1,JA:JA+N-1) * diag(C).
//
//  ABSTOL  (global input) REAL
//          If JOBZ='V', setting ABSTOL to PSLAMCH( CONTEXT, 'U') yields
//          the most orthogonal eigenvectors.
//
//          The absolute error tolerance for the eigenvalues.
//          An approximate eigenvalue is accepted as converged
//          when it is determined to lie in an interval [a,b]
//          of width less than or equal to
//
//                  ABSTOL + EPS *   max( |a|,|b| ) ,
//
//          where EPS is the machine precision.  If ABSTOL is less than
//          or equal to zero, then EPS*norm(T) will be used in its place,
//          where norm(T) is the 1-norm of the tridiagonal matrix
//          obtained by reducing A to tridiagonal form.
//
//          Eigenvalues will be computed most accurately when ABSTOL is
//          set to twice the underflow threshold 2*PSLAMCH('S') not zero.
//          If this routine returns with ((MOD(INFO,2).NE.0) .OR.
//          (MOD(INFO/8,2).NE.0)), indicating that some eigenvalues or
//          eigenvectors did not converge, try setting ABSTOL to
//          2*PSLAMCH('S').
//
//          See "Computing Small Singular Values of Bidiagonal Matrices
//          with Guaranteed High Relative Accuracy," by Demmel and
//          Kahan, LAPACK Working Note #3.
//
//          See "On the correctness of Parallel Bisection in Floating
//          Point" by Demmel, Dhillon and Ren, LAPACK Working Note #70
//
//  AF      (local input or local output) REAL pointer
//          into the local memory to an array of local dimension
//          (LLD_AF,LOCc(JA+N-1)).  If FACT = 'F', then
//          AF(IAF:IAF+N-1,JAF:JAF+N-1) is an input argument and on
//          entry contains the factors L and U from the factorization
//          A(IA:IA+N-1,JA:JA+N-1) = P*L*U as computed by PSGETRF.
//          If EQUED .ne. 'N', then AF is the factored form of the
//          equilibrated matrix A(IA:IA+N-1,JA:JA+N-1).
//
//          If FACT = 'N', then AF(IAF:IAF+N-1,JAF:JAF+N-1) is an output
//          argument and on exit returns the factors L and U from the
//          factorization A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the original
//          matrix A(IA:IA+N-1,JA:JA+N-1).
//
//          If FACT = 'E', then AF(IAF:IAF+N-1,JAF:JAF+N-1) is an output
//          argument and on exit returns the factors L and U from the
//          factorization A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the equili-
//          brated matrix A(IA:IA+N-1,JA:JA+N-1) (see the description of
//          A(IA:IA+N-1,JA:JA+N-1) for the form of the equilibrated
//          matrix).
//
//  B       (local input/local output) REAL pointer into
//          local memory to an array of local lead dimension lld_b>=NB.
//          On entry, this array contains the
//          the local pieces of the right hand sides
//          B(IB:IB+N-1, 1:NRHS).
//          On exit, this contains the local piece of the solutions
//          distributed matrix X.
//
//  BERR    (local output) REAL array, dimension LOCc(N_B).
//          The componentwise relative backward error of each solution
//          vector X(j) (i.e., the smallest relative change in
//          any entry of A(IA:IA+N-1,JA:JA+N-1) or
//          B(IB:IB+N-1,JB:JB+NRHS-1) that makes X(j) an exact solution).
//          BERR is replicated in every process row, and is aligned
//          with the matrices B and X.
//
//  BW      (global input) INTEGER
//          Number of subdiagonals in L or U. 0 <= BW <= N-1
//
//  BWL     (global input) INTEGER
//          Number of subdiagonals. 0 <= BWL <= N-1
//
//  BWU     (global input) INTEGER
//          Number of superdiagonals. 0 <= BWU <= N-1
//
//  C       (local input or local output) REAL array,
//                                                  dimension LOCc(N_A).
//          The column scale factors for A(IA:IA+N-1,JA:JA+N-1).
//          If EQUED = 'C' or 'B', A(IA:IA+N-1,JA:JA+N-1) is multiplied
//          on the right by diag(C); if EQUED = 'N' or 'R', C is not
//          accessed.  C is an input variable if FACT = 'F'; otherwise,
//          C is an output variable.  If FACT = 'F' and EQUED = 'C' or
//          'B', each element of C must be positive.
//          C is replicated in every process row, and is aligned with
//          the distributed matrix A.
//
//  D       (local input/local output) REAL pointer to local
//          part of global vector storing the main diagonal of the
//          matrix.
//          On exit, this array contains information containing the
//            factors of the matrix.
//          Must be of size >= DESCA( NB_ ).
//
//  DL      (local input/local output) REAL pointer to local
//          part of global vector storing the lower diagonal of the
//          matrix. Globally, DL(1) is not referenced, and DL must be
//          aligned with D.
//          Must be of size >= DESCA( NB_ ).
//          On exit, this array contains information containing the
//            factors of the matrix.
//
//  DU       (local input/local output) REAL pointer to local
//          part of global vector storing the upper diagonal of the
//          matrix. Globally, DU(n) is not referenced, and DU must be
//          aligned with D.
//          On exit, this array contains information containing the
//            factors of the matrix.
//          Must be of size >= DESCA( NB_ ).
//
//  EQUED   (global input or global output) CHARACTER
//          Specifies the form of equilibration that was done.
//          = 'N':  No equilibration (always true if FACT = 'N').
//          = 'R':  Row equilibration, i.e., A(IA:IA+N-1,JA:JA+N-1) has
//                  been premultiplied by diag(R).
//          = 'C':  Column equilibration, i.e., A(IA:IA+N-1,JA:JA+N-1)
//                  has been postmultiplied by diag(C).
//          = 'B':  Both row and column equilibration, i.e.,
//                  A(IA:IA+N-1,JA:JA+N-1) has been replaced by
//                  diag(R) * A(IA:IA+N-1,JA:JA+N-1) * diag(C).
//          EQUED is an input variable if FACT = 'F'; otherwise, it is an
//          output variable.
//
//  FACT    (global input) CHARACTER
//          Specifies whether or not the factored form of the matrix
//          A(IA:IA+N-1,JA:JA+N-1) is supplied on entry, and if not,
//          whether the matrix A(IA:IA+N-1,JA:JA+N-1) should be
//          equilibrated before it is factored.
//          = 'F':  On entry, AF(IAF:IAF+N-1,JAF:JAF+N-1) and IPIV con-
//                  tain the factored form of A(IA:IA+N-1,JA:JA+N-1).
//                  If EQUED is not 'N', the matrix
//                  A(IA:IA+N-1,JA:JA+N-1) has been equilibrated with
//                  scaling factors given by R and C.
//                  A(IA:IA+N-1,JA:JA+N-1), AF(IAF:IAF+N-1,JAF:JAF+N-1),
//                  and IPIV are not modified.
//          = 'N':  The matrix A(IA:IA+N-1,JA:JA+N-1) will be copied to
//                  AF(IAF:IAF+N-1,JAF:JAF+N-1) and factored.
//          = 'E':  The matrix A(IA:IA+N-1,JA:JA+N-1) will be equili-
//                  brated if necessary, then copied to
//                  AF(IAF:IAF+N-1,JAF:JAF+N-1) and factored.
//
//  FERR    (local output) REAL array, dimension LOCc(N_B)
//          The estimated forward error bounds for each solution vector
//          X(j) (the j-th column of the solution matrix
//          X(IX:IX+N-1,JX:JX+NRHS-1). If XTRUE is the true solution,
//          FERR(j) bounds the magnitude of the largest entry in
//          (X(j) - XTRUE) divided by the magnitude of the largest entry
//          in X(j).  The estimate is as reliable as the estimate for
//          RCOND, and is almost always a slight overestimate of the
//          true error.  FERR is replicated in every process row, and is
//          aligned with the matrices B and X.
//
//  GAP     (global output) REAL array,
//             dimension (NPROW*NPCOL)
//          This array contains the gap between eigenvalues whose
//          eigenvectors could not be reorthogonalized. The output
//          values in this array correspond to the clusters indicated
//          by the array ICLUSTR. As a result, the dot product between
//          eigenvectors correspoding to the I^th cluster may be as high
//          as ( C * n ) / GAP(I) where C is a small constant.
//
// IBTYPE   (global input) INTEGER
//          Specifies the problem type to be solved:
//          = 1:  sub( A )*x = (lambda)*sub( B )*x
//          = 2:  sub( A )*sub( B )*x = (lambda)*x
//          = 3:  sub( B )*sub( A )*x = (lambda)*x
//
//  ICLUSTR (global output) integer array, dimension (2*NPROW*NPCOL)
//          This array contains indices of eigenvectors corresponding to
//          a cluster of eigenvalues that could not be reorthogonalized
//          due to insufficient workspace (see LWORK, ORFAC and INFO).
//          Eigenvectors corresponding to clusters of eigenvalues indexed
//          ICLUSTR(2*I-1) to ICLUSTR(2*I), could not be
//          reorthogonalized due to lack of workspace. Hence the
//          eigenvectors corresponding to these clusters may not be
//          orthogonal.  ICLUSTR() is a zero terminated array.
//          (ICLUSTR(2*K).NE.0 .AND. ICLUSTR(2*K+1).EQ.0) if and only if
//          K is the number of clusters
//          ICLUSTR is not referenced if JOBZ = 'N'
//
//  IFAIL   (global output) INTEGER array, dimension (N)
//          If JOBZ = 'V', then on normal exit, the first M elements of
//          IFAIL are zero.  If (MOD(INFO,2).NE.0) on exit, then
//          IFAIL contains the
//          indices of the eigenvectors that failed to converge.
//          If JOBZ = 'N', then IFAIL is not referenced.
//
//  IL      (global input) INTEGER
//          If RANGE='I', the index (from smallest to largest) of the
//          smallest eigenvalue to be returned.  IL >= 1.
//          Not referenced if RANGE = 'A' or 'V'.
//
//  INFO    (global output) INTEGER
//          = 0:  successful exit
//          < 0:  If the i-th argument is an array and the j-entry had
//                an illegal value, then INFO = -(i*100+j), if the i-th
//                argument is a scalar and had an illegal value, then
//                INFO = -i.
//          > 0:  If INFO = K, U(IA+K-1,JA+K-1) is exactly zero.
//                The factorization has been completed, but the factor U
//                is exactly singular, so the solution could not be
//                computed. (for psgesv. see documentation for indiviual
//                functions what information info > 0 gives)
//
//  IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
//          This array contains the pivoting information.
//          IPIV(i) -> The global row local row i was swapped with.
//          This array is tied to the distributed matrix A.
//
//          If FACT = 'F', then IPIV is an input argu-
//          ment and on entry contains the pivot indices from the fac-
//          torization A(IA:IA+N-1,JA:JA+N-1) = P*L*U as computed by
//          PSGETRF; IPIV(i) -> The global row local row i was
//          swapped with.  This array must be aligned with
//          A( IA:IA+N-1, * ).
//
//          If FACT = 'N', then IPIV is an output argument and on exit
//          contains the pivot indices from the factorization
//          A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the original matrix
//          A(IA:IA+N-1,JA:JA+N-1).
//
//          If FACT = 'E', then IPIV is an output argument and on exit
//          contains the pivot indices from the factorization
//          A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the equilibrated matrix
//          A(IA:IA+N-1,JA:JA+N-1).
//
//  IU      (for pssyevx) (global input) INTEGER
//          If RANGE='I', the index (from smallest to largest) of the
//          largest eigenvalue to be returned.  min(IL,N) <= IU <= N.
//          Not referenced if RANGE = 'A' or 'V'.
//
//  IWORK   (local workspace/local output) INTEGER array,
//                                                  dimension (LIWORK)
//          On exit, IWORK(1) returns the minimal and optimal LIWORK.
//
//  JOBU    (global input) CHARACTER*1
//          Specifies options for computing U:
//          = 'V':  the first SIZE columns of U (the left singular
//                  vectors) are returned in the array U;
//          = 'N':  no columns of U (no left singular vectors) are
//                  computed.
//
//  JOBVT   (global input) CHARACTER*1
//          Specifies options for computing V**T:
//          = 'V':  the first SIZE rows of V**T (the right singular
//                  vectors) are returned in the array VT;
//          = 'N':  no rows of V**T (no right singular vectors) are
//                  computed.
//
//  JOBZ    (global input) CHARACTER*1
//               Specifies whether or not to compute the eigenvectors:
//               = 'N':  Compute eigenvalues only.
//               = 'V':  Compute eigenvalues and eigenvectors.
//
//  LWORK   (local input or global input) INTEGER
//          Size of user-input workspace WORK.
//          If LWORK is too small, the minimal acceptable size will be
//          returned in WORK(1) and an error code is returned. LWORK>=
//          NB*(bwl+bwu)+6*max(bwl,bwu)*max(bwl,bwu)
//          +max((max(bwl,bwu)*NRHS), max(bwl,bwu)*max(bwl,bwu))
//          (description for psdbsv) see documentation for individual
//          functions for the meaning of specific return values of lwork
//
//  M       ( pssyevx) (global output) INTEGER
//          Total number of eigenvalues found.  0 <= M <= N.
//
//  NRHS    (global input) INTEGER
//          The number of right hand sides, i.e., the number of columns
//          of the distributed submatrix sub( B ). NRHS >= 0.
//
//  NZ      (global output) INTEGER
//          Total number of eigenvectors computed.  0 <= NZ <= M.
//          The number of columns of Z that are filled.
//          If JOBZ .NE. 'V', NZ is not referenced.
//          If JOBZ .EQ. 'V', NZ = M unless the user supplies
//          insufficient space and PSSYEVX is not able to detect this
//          before beginning computation.  To get all the eigenvectors
//          requested, the user must supply both sufficient
//          space to hold the eigenvectors in Z (M .LE. DESCZ(N_))
//          and sufficient workspace to compute them.  (See LWORK below.)
//          PSSYEVX is always able to detect insufficient space without
//          computation unless RANGE .EQ. 'V'.
//
//  ORFAC   (global input) REAL
//          Specifies which eigenvectors should be reorthogonalized.
//          Eigenvectors that correspond to eigenvalues which are within
//          tol=ORFAC*norm(A) of each other are to be reorthogonalized.
//          However, if the workspace is insufficient (see LWORK),
//          tol may be decreased until all eigenvectors to be
//          reorthogonalized can be stored in one process.
//          No reorthogonalization will be done if ORFAC equals zero.
//          A default value of 10^-3 is used if ORFAC is negative.
//          ORFAC should be identical on all processes.
//
//  R       (local input or local output) REAL array,
//                                                  dimension LOCr(M_A).
//          The row scale factors for A(IA:IA+N-1,JA:JA+N-1).
//          If EQUED = 'R' or 'B', A(IA:IA+N-1,JA:JA+N-1) is multiplied
//          on the left by diag(R); if EQUED='N' or 'C', R is not acces-
//          sed.  R is an input variable if FACT = 'F'; otherwise, R is
//          an output variable.
//          If FACT = 'F' and EQUED = 'R' or 'B', each element of R must
//          be positive.
//          R is replicated in every process column, and is aligned
//          with the distributed matrix A.
//
// RANGE   (global input) CHARACTER*1
//          = 'A': all eigenvalues will be found.
//          = 'V': all eigenvalues in the interval [VL,VU] will be found.
//          = 'I': the IL-th through IU-th eigenvalues will be found.
//
//  RCOND   (global output) REAL
//          The estimate of the reciprocal condition number of the matrix
//          A(IA:IA+N-1,JA:JA+N-1) after equilibration (if done).  If
//          RCOND is less than the machine precision (in particular, if
//          RCOND = 0), the matrix is singular to working precision.
//          This condition is indicated by a return code of INFO > 0.
//
//  S       (global output) REAL               array, dimension SIZE
//          The singular values of A, sorted so that S(i) >= S(i+1).
//
//  SC      (local input/local output) REAL array,
//                                              dimension (LOC(N_A))
//          The scale factors for A distributed across
//          process columns; not accessed if EQUED = 'N'. SC is an input
//          variable if FACT = 'F'; otherwise, SC is an output variable.
//          If FACT = 'F' and EQUED = 'Y', each element of SC must be
//          positive.
//
//  SR      (local input/local output) REAL array,
//                                                    dimension (LLD_A)
//          The scale factors for A distributed across process rows;
//          not accessed if EQUED = 'N'.  SR is an input variable if
//          FACT = 'F'; otherwise, SR is an output variable.
//          If FACT = 'F' and EQUED = 'Y', each element of SR must be
//          positive.
//
//  U       (local output) REAL               array, local dimension
//          (MP, SIZEQ), global dimension (M, SIZE)
//          if JOBU = 'V', U contains the first min(m,n) columns of U
//          if JOBU = 'N', U is not referenced.
//
//  VL      (global input) REAL
//          If RANGE='V', the lower bound of the interval to be searched
//          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
//
//  VT      (local output) REAL               array, local dimension
//          (SIZEP, NQ), global dimension (SIZE, N).
//          If JOBVT = 'V', VT contains the first SIZE rows of
//          V**T. If JOBVT = 'N', VT is not referenced.
//
//  VU      (global input) REAL
//          If RANGE='V', the upper bound of the interval to be searched
//          for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
//
//  W       (global output) DOUBLE PRECISION array, dimension (N)
//               On normal exit, the first M entries contain the selected eigen-
//               values in ascending order.
//
//  WORK    (local workspace/local output)
//          REAL temporary workspace. This space may
//          be overwritten in between calls to routines. WORK must be
//          the size given in LWORK.
//          On exit, WORK( 1 ) contains the minimal LWORK.
//
//  X       (local input/local output) REAL pointer
//          into the local memory to an array of local dimension
//          (LLD_X, LOCc(JX+NRHS-1)).  If INFO = 0, the N-by-NRHS
//          solution matrix X(IX:IX+N-1,JX:JX+NRHS-1) to the original
//          system of equations.  Note that A(IA:IA+N-1,JA:JA+N-1) and
//          B(IB:IB+N-1,JB:JB+NRHS-1) are modified on exit if
//          EQUED .ne. 'N', and the solution to the equilibrated system
//          is inv(diag(C))*X(IX:IX+N-1,JX:JX+NRHS-1) if TRANS = 'N'
//          and EQUED = 'C' or 'B', or
//          inv(diag(R))*X(IX:IX+N-1,JX:JX+NRHS-1) if TRANS = 'T' or 'C'
//          and EQUED = 'R' or 'B'.
//
//  Z       (local output) DOUBLE PRECISION array,
//               global  dimension (N, N), local dimension ( LLD_Z, LOCc(JZ+N-1)
//               ) If JOBZ = 'V', then on normal exit the first M columns  of  Z
//               contain    the   orthonormal   eigenvectors   of   the   matrix
//               corresponding to the selected eigenvalues.  If JOBZ = 'N', then
//               Z is not referenced.
//
//*******************************************************
//
//************SYSTEMS OF LINEAR EQUATIONS*****************
//
//----------------single precision----------------------------
extern void psgesv_( int* n, int* nrhs, float* A, int* ia, int* ja, float* descA, int* ipiv, float* B, int* ib, int* jb, float* descB, int* info );
extern void pcgesv_( int* n, int* nrhs, float* A, int* ia, int* ja, float* descA, int* ipiv, float* B, int* ib, int* jb, float* descB, int* info );
extern void psdbsv_( int* n, int* bwl, int* bwu, int* nrhs, float* A, int* ja, int* desca, float* B, int* ib, int* descb, float* work, int* lwork, int* info);
extern void pcdbsv_( int* n, int* bwl, int* bwu, int* nrhs, float* A, int* ja, int* desca, float* B, int* ib, int* descb, float* work, int* lwork, int* info);
extern void psgbsv_( int* n, int* bwl, int* bwu, int* nrhs, float* A, int* ja, int* desca, int* ipiv, float* B, int* ib, int* descb, float* work, int* lwork, int* info);
extern void pcgbsv_( int* n, int* bwl, int* bwu, int* nrhs, float* A, int* ja, int* desca, int* ipiv, float* B, int* ib, int* descb, float* work, int* lwork, int* info);
extern void psdtsv_( int* n, int* nrhs, float* DL, float* D, float* DU, int* ja, int* desca, float* B, int* ib, int* descb, float* work, int* lwork, int* info);
extern void pcdtsv_( int* n, int* nrhs, float* DL, float* D, float* DU, int* ja, int* desca, float* B, int* ib, int* descb, float* work, int* lwork, int* info);
extern void psposv_( char* uplo, int* n, int* nrhs, float* A, int* ia, int* ja, int* desca, float* B, int* ib, int* jb, int* descb, int* info);
extern void pcposv_( char* uplo, int* n, int* nrhs, float* A, int* ia, int* ja, int* desca, float* B, int* ib, int* jb, int* descb, int* info);
extern void pspbsv_( char* uplo, int* n, int* bw, int* nrhs, float* A, int* ja, int* desca, float* B, int* ib, int* descb, float* work, int* lwork, int* info);
extern void pcpbsv_( char* uplo, int* n, int* bw, int* nrhs, float* A, int* ja, int* desca, float* B, int* ib, int* descb, float* work, int* lwork, int* info);
extern void psptsv_( int* n, int* nrhs, float* D, float* DU, int* ja, int* desca, float* B, int* ib, int* descb, float* work, int* lwork, int* info );
extern void pcptsv_( int* n, int* nrhs, float* D, float* DU, int* ja, int* desca, float* B, int* ib, int* descb, float* work, int* lwork, int* info );
//
//----------------double precision----------------------------
extern void pdgesv_( int* n, int* nrhs, double* A, int* ia, int* ja, int* descA, int* ipiv, double* B, int* ib, int* jb, int* descB, int* info );
extern void pzgesv_( int* n, int* nrhs, double* A, int* ia, int* ja, int* descA, int* ipiv, double* B, int* ib, int* jb, int* descB, int* info );
extern void pddbsv_( int* n, int* bwl, int* bwu, int* nrhs, double* A, int* ja, int* desca, double* B, int* ib, int* descb, double* work, int* lwork, int* info);
extern void pzdbsv_( int* n, int* bwl, int* bwu, int* nrhs, double* A, int* ja, int* desca, double* B, int* ib, int* descb, double* work, int* lwork, int* info);
extern void pdgbsv_( int* n, int* bwl, int* bwu, int* nrhs, double* A, int* ja, int* desca, int* ipiv, double* B, int* ib, int* descb, double* work, int* lwork, int* info);
extern void pzgbsv_( int* n, int* bwl, int* bwu, int* nrhs, double* A, int* ja, int* desca, int* ipiv, double* B, int* ib, int* descb, double* work, int* lwork, int* info);
extern void pddtsv_( int* n, int* nrhs, double* DL, double* D, double* DU, int* ja, int* desca, double* B, int* ib, int* descb, double* work, int* lwork, int* info);
extern void pzdtsv_( int* n, int* nrhs, double* DL, double* D, double* DU, int* ja, int* desca, double* B, int* ib, int* descb, double* work, int* lwork, int* info);
extern void pdposv_( char* uplo, int* n, int* nrhs, double* A, int* ia, int* ja, int* desca, double* B, int* ib, int* jb, int* descb, int* info);
extern void pzposv_( char* uplo, int* n, int* nrhs, double* A, int* ia, int* ja, int* desca, double* B, int* ib, int* jb, int* descb, int* info);
extern void pdpbsv_( char* uplo, int* n, int* bw, int* nrhs, double* A, int* ja, int* desca, double* B, int* ib, int* descb, double* work, int* lwork, int* info);
extern void pzpbsv_( char* uplo, int* n, int* bw, int* nrhs, double* A, int* ja, int* desca, double* B, int* ib, int* descb, double* work, int* lwork, int* info);
extern void pdptsv_( int* n, int* nrhs, double* D, double* DU, int* ja, int* desca, double* B, int* ib, int* descb, double* work, int* lwork, int* info );
extern void pzptsv_( int* n, int* nrhs, double* D, double* DU, int* ja, int* desca, double* B, int* ib, int* descb, double* work, int* lwork, int* info );
//
//********************LEAST SQUARES PROBLEMS************
//
//----------------single precision----------------------------
extern void psgela_( char* trans, int* m, int* n, int* nrhs, float* A, int* ia, int* ja, int* desca, float* B, int* ib, int* jb, int* descb, float* work, int* lwork, int* info );
extern void pcgela_( char* trans, int* m, int* n, int* nrhs, float* A, int* ia, int* ja, int* desca, float* B, int* ib, int* jb, int* descb, float* work, int* lwork, int* info );
//
//----------------double precision----------------------------
extern void pdgela_( char* trans, int* m, int* n, int* nrhs, double* A, int* ia, int* ja, int* desca, double* B, int* ib, int* jb, int* descb, double* work, int* lwork, int* info );
extern void pzgela_( char* trans, int* m, int* n, int* nrhs, double* A, int* ia, int* ja, int* desca, float* B, int* ib, int* jb, int* descb, double* work, int* lwork, int* info );
//
//********EIGENVALUE AND SINGULAR VALUE PROBLEMS*****************
//
//------------------single precision---------------------------
extern void pssyev_ ( char *jobz, char *uplo, int *n, float *A, int *ia, int *ja, int *desca, float *w, float *z, int *iz, int *jz, int *descz, float *work, int *lwork, int *info );
extern void psgesvd_( char *jobu, char *jobvt, int* m, int *n, float *A, int *ia, int *ja, int *desca, float *s, float *u, int *iu, int *ju, int *descu, float* vt, int* ivt, int* jvt, int* descvt, float *work, int *lwork, int *info );
//
//------------------double precision---------------------------
extern void pdsyev_ ( char *jobz, char *uplo, int *n, double *A, int *ia, int *ja, int *desca, double *w, double *z, int *iz, int *jz, int *descz, double *work, int *lwork, int *info );
extern void pzgesvd_( char *jobu, char *jobvt, int* m, int *n, double *A, int *ia, int *ja, int *desca, double *s, double *u, int *iu, int *ju, int *descu, double* vt, int* ivt, int* jvt, int* descvt, double *work, int *lwork, int *info );
//
//**************EXPERT DRIVERS: SYSTEMS OF LINEAR EQUATIONS************
//
//----------------single precision------------------------------
extern void psgesvx_(char* fact, char* trans, int* n, int* nrhs, float* A, int* ia, int* ja, int* desca, float* AF, int* iaf, int* jaf, int* descaf, int* ipiv, char* equed, float* R, float* C, float* B, int* ib, int* jb, int* descb, float* X, int* ix, int* jx, int* descx, float* rcond, float* ferr, float* berr, float *work, int *lwork, int* iwork, int* liwork, int *info); 
extern void pcgesvx_(char* fact, char* trans, int* n, int* nrhs, float* A, int* ia, int* ja, int* desca, float* AF, int* iaf, int* jaf, int* descaf, int* ipiv, char* equed, float* R, float* C, float* B, int* ib, int* jb, int* descb, float* X, int* ix, int* jx, int* descx, float* rcond, float* ferr, float* berr, float *work, int *lwork, int* iwork, int* liwork, int *info);
extern void psposcx_(char* fact, char* trans, int* n, int* nrhs, float* A, int* ia, int* ja, int* desca, float* AF, int* iaf, int* jaf, int* descaf, char* equed, float* S, float* B, int* ib, int* jb, int* descb, float* X, int* ix, int* jx, int* descx, float* rcond, float* ferr, float* berr, float *work, int *lwork, int* iwork, int* liwork, int *info);
extern void pcposvx_(char* fact, char* trans, int* n, int* nrhs, float* A, int* ia, int* ja, int* desca, float* AF, int* iaf, int* jaf, int* descaf, char* equed, float* SR, float* SC, float* B, int* ib, int* jb, int* descb, float* X, int* ix, int* jx, int* descx, float* rcond, float* ferr, float* berr, float *work, int *lwork, int* iwork, int* liwork, int *info);
//
//----------------double precision------------------------------
extern void pdgesvx_(char* fact, char* trans, int* n, int* nrhs, double* A, int* ia, int* ja, int* desca, double* AF, int* iaf, int* jaf, int* descaf, int* ipiv, char* equed, double* R, double* C, double* B, int* ib, int* jb, int* descb, double* X, int* ix, int* jx, int* descx, double* rcond, double* ferr, double* berr, double *work, int *lwork, int* iwork, int* liwork, int *info); 
extern void pzgesvx_(char* fact, char* trans, int* n, int* nrhs, double* A, int* ia, int* ja, int* desca, double* AF, int* iaf, int* jaf, int* descaf, int* ipiv, char* equed, double* R, double* C, double* B, int* ib, int* jb, int* descb, double* X, int* ix, int* jx, int* descx, double* rcond, double* ferr, double* berr, double *work, int *lwork, int* iwork, int* liwork, int *info);
extern void pdposcx_(char* fact, char* trans, int* n, int* nrhs, double* A, int* ia, int* ja, int* desca, double* AF, int* iaf, int* jaf, int* descaf, char* equed, double* S, double* B, int* ib, int* jb, int* descb, double* X, int* ix, int* jx, int* descx, double* rcond, double* ferr, double* berr, double *work, int *lwork, int* iwork, int* liwork, int *info);
extern void pzposvx_(char* fact, char* trans, int* n, int* nrhs, double* A, int* ia, int* ja, int* desca, double* AF, int* iaf, int* jaf, int* descaf, char* equed, double* SR, double* SC, double* B, int* ib, int* jb, int* descb, double* X, int* ix, int* jx, int* descx, double* rcond, double* ferr, double* berr, double *work, int *lwork, int* iwork, int* liwork, int *info);
//
//*************EXPERT EIGENVALUE PROBLEMS***********************
//
//---------------single precision------------------------------
extern void pssyevx_( char* jobz, char* range, char* uplo, int* n, float* A, int* ia, int* ja, int* desca, float* vl, float* vu, int* il, int* iu, float* abstol, int* m, int* nz, float* W, float* orfac, float* z, int* iz, int* jz, int* descz, float* work, int* lwork, int* iwork, int* liwork, int* ifail, int* iclustr, float* gap, int* info );
extern void pcheevx_( char* jobz, char* range, char* uplo, int* n, float* A, int* ia, int* ja, int* desca, float* vl, float* vu, int* il, int* iu, float* abstol, int* m, int* nz, float* W, float* orfac, float* z, int* iz, int* jz, int* descz, float* work, int* lwork, int* iwork, int* liwork, int* ifail, int* iclustr, float* gap, int* info );
extern void pssygvx_( int* ibtype, char* jobz, char* range, char* uplo, int* n, float* A, int* ia, int* ja, int* desca, float* B, int*ib, int* jb, int* descB, float* vl, float* vu, int* il, int* iu, float* abstol, int* m, int* nz, float* W, float* orfac, float* z, int* iz, int* jz, int* descz, float* work, int* lwork, int* iwork, int* liwork, int* ifail, int* iclustr, float* gap, int* info );
extern void pchegvx_( int* ibtype, char* jobz, char* range, char* uplo, int* n, float* A, int* ia, int* ja, int* desca, float* B, int*ib, int* jb, int* descB, float* vl, float* vu, int* il, int* iu, float* abstol, int* m, int* nz, float* W, float* orfac, float* z, int* iz, int* jz, int* descz, float* work, int* lwork, int* iwork, int* liwork, int* ifail, int* iclustr, float* gap, int* info );
//
//---------------double precision------------------------------
extern void pdsyevx_( char* jobz, char* range, char* uplo, int* n, double* A, int* ia, int* ja, int* desca, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, int* nz, double* W, double* orfac, double* z, int* iz, int* jz, int* descz, double* work, int* lwork, int* iwork, int* liwork, int* ifail, int* iclustr, double* gap, int* info );
extern void pzheevx_( char* jobz, char* range, char* uplo, int* n, double* A, int* ia, int* ja, int* desca, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, int* nz, double* W, double* orfac, double* z, int* iz, int* jz, int* descz, double* work, int* lwork, int* iwork, int* liwork, int* ifail, int* iclustr, double* gap, int* info );
extern void pdsygvx_( int* ibtype, char* jobz, char* range, char* uplo, int* n, double* A, int* ia, int* ja, int* desca, double* B, int*ib, int* jb, int* descB, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, int* nz, double* W, double* orfac, double* z, int* iz, int* jz, int* descz, double* work, int* lwork, int* iwork, int* liwork, int* ifail, int* iclustr, double* gap, int* info );
extern void pzhegvx_( int* ibtype, char* jobz, char* range, char* uplo, int* n, double* A, int* ia, int* ja, int* desca, double* B, int*ib, int* jb, int* descB, double* vl, double* vu, int* il, int* iu, double* abstol, int* m, int* nz, double* W, double* orfac, double* z, int* iz, int* jz, int* descz, double* work, int* lwork, int* iwork, int* liwork, int* ifail, int* iclustr, double* gap, int* info );

/*#ifdef F77_WITH_NO_UNDERSCORE
#define   numroc_      numroc
#define   descinit_    descinit
#define   pdlamch_     pdlamch
#define   pdlange_     pdlange
#define   pdlacpy_     pdlacpy
#define   pdgesv_      pdgesv
#define   pdgemm_      pdgemm
#define   indxg2p_     indxg2p
#endif*/
