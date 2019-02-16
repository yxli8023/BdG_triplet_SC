      program ex01
      implicit none
      INTEGER          N
      PARAMETER        ( N = 4 )
      INTEGER          LDA, LDVL, LDVR
      PARAMETER        ( LDA = N, LDVL = N, LDVR = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
*     RWORK dimension should be at least 2*N
      DOUBLE PRECISION RWORK( 2*N )
      COMPLEX*16       A( LDA, N ), VL( LDVL, N ), VR( LDVR, N ),
     $                 W( N ), WORK( LWMAX )
      DATA             A/
     $ (-3.84, 2.25),(-0.66, 0.83),(-3.99,-4.73),( 7.74, 4.18),
     $ (-8.94,-4.75),(-4.40,-3.82),(-5.88,-6.60),( 3.66,-7.53),
     $ ( 8.95,-6.53),(-3.50,-4.26),(-3.36,-0.40),( 2.58, 3.60),
     $ (-9.87, 4.82),(-3.15, 7.36),(-0.75, 5.23),( 4.59, 5.41)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         ZGEEV
      EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'ZGEEV Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL ZGEEV( 'Vectors', 'Vectors', N, A, LDA, W, VL, LDVL,
     $            VR, LDVR, WORK, LWORK, RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      write(*,*)lwork,work(1)
*
*     Solve eigenproblem.
*
      CALL ZGEEV( 'Vectors', 'Vectors', N, A, LDA, W, VL, LDVL,
     $            VR, LDVR, WORK, LWORK, RWORK, INFO )
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
      CALL PRINT_MATRIX( 'Eigenvalues', 1, N, W, 1 )
*
*     Print left eigenvectors.
*
      CALL PRINT_MATRIX( 'Left eigenvectors', N, N, VL, LDVL )
*
*     Print right eigenvectors.
*
      CALL PRINT_MATRIX( 'Right eigenvectors', N, N, VR, LDVR )
      STOP
      END
*
*     End of ZGEEV Example.
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      COMPLEX*16       A( LDA, * )
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
