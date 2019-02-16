
* .. Parameters ..
      program ex01
      implicit none
      INTEGER          N
      PARAMETER        ( N = 4 )
      INTEGER          LDA
      PARAMETER        ( LDA = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK, LIWORK, LRWORK
*
*     .. Local Arrays ..
      INTEGER          IWORK( LWMAX )
      REAL             W( N ), RWORK( LWMAX )
      COMPLEX          A( LDA, N ), WORK( LWMAX )
      DATA             A/
     $ ( 3.40, 0.00),(-2.36, 1.93),(-4.68,-9.55),( 5.37, 1.23),
     $ ( 0.00, 0.00),( 6.94, 0.00),( 8.13, 1.47),( 2.07, 5.78),
     $ ( 0.00, 0.00),( 0.00, 0.00),(-2.14, 0.00),( 4.68,-7.44),
     $ ( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),(-7.42, 0.00)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         CHEEVD
      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'CHEEVD Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      LIWORK = -1
      LRWORK = -1
      CALL CHEEVD( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,
     $             LRWORK, IWORK, LIWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LRWORK = MIN( LWMAX, INT( RWORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK( 1 ) )
*
*     Solve eigenproblem.
*
      CALL CHEEVD( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,
     $             LRWORK, IWORK, LIWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
*
*     Print eigenvalues.
*
      CALL PRINT_RMATRIX( 'Eigenvalues', 1, N, W, 1 )
*
*     Print eigenvectors.
*
      CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', N, N, A,
     $                   LDA )
      STOP
      END
*
*     End of CHEEVD Example.
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      COMPLEX          A( LDA, * )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END
*
*     Auxiliary routine: printing a real matrix.
*
      SUBROUTINE PRINT_RMATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      REAL             A( LDA, * )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,F6.2) )
      RETURN
      END
