        !COMPILER-GENERATED INTERFACE MODULE: Wed Feb 20 11:16:05 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DELTA__genmod
          INTERFACE 
            FUNCTION DELTA(X,Y,EIG_VAL,EIG_VEC)
              INTEGER(KIND=4) :: X
              INTEGER(KIND=4) :: Y
              REAL(KIND=4) :: EIG_VAL(36)
              COMPLEX(KIND=4) :: EIG_VEC(36,36)
              COMPLEX(KIND=4) :: DELTA
            END FUNCTION DELTA
          END INTERFACE 
        END MODULE DELTA__genmod
