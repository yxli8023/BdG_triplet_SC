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
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
*     RWORK dimension should be at least MAX(1,3*N-2)
      REAL             W( N ), RWORK( 3*N-2 )
      COMPLEX          A( LDA, N ), WORK( LWMAX )
      DATA             A/
     $ ( 9.14, 0.00),(-4.37, 9.22),(-1.98, 1.72),(-8.96, 9.50),
     $ ( 0.00, 0.00),(-3.35, 0.00),( 2.25, 9.51),( 2.57,-2.40),
     $ ( 0.00, 0.00),( 0.00, 0.00),(-4.82, 0.00),(-3.24,-2.04),
     $ ( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),( 8.44, 0.00)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         CHEEV
      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
      write(*,*)"A:",A
*     .. Executable Statements ..
      WRITE(*,*)'CHEEV Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL CHEEV( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,
     $            INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve eigenproblem.
*
      CALL CHEEV( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,
     $            INFO )
*
      write(*,*)"A change:",A
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
*     End of CHEEV Example.
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
