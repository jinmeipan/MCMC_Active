! -------------------------------------------------------------------------
!
SUBROUTINE INTEGRMUO(XX,MUI,MINO,MAXO,STEPS,INTEGR,NUM,ARG_LENGTH)
!
! -------------------------------------------------------------------------
!
!     CODE ORIGINALLY OBTAINED IN MATLAB FROM MATZLER
!
!   CALCULATES INTEGRATION  OVER INCIDENT DIRECTIONS, MUI
!       FROM MINO TO MAXO IN STEPS INTERVALS OF THE FI
!       PHASE FUNCTION.
!
!
!   INTEGR = INTEGRMUO(XX,MUI,MINO,MAXO,STEPS)
!
!   VERSION HISTORY:
!      1.0     WI 27.05.98
!      2.0    MD 1 APR 05 TRANSLATED TO FORTRAN FROM MATLAB
!      2.1    MD 21 NOV 05 MADE ALL LOCALS ALLOCATABLE
!
!   USES: INTEGRFI
!
!   COPYRIGHT (C) 1998 BY THE INSTITUTE OF APPLIED PHYSICS,
!   UNIVERSITY OF BERN, SWITZERLAND

IMPLICIT NONE

INTEGER,INTENT(IN) :: NUM,STEPS,ARG_LENGTH(2)
REAL(8),INTENT(IN) :: MUI(NUM),MINO(ARG_LENGTH(1)),MAXO(ARG_LENGTH(2)),XX(NUM)
REAL(8),INTENT(OUT) :: INTEGR(NUM)
REAL(8) :: MINOS,MAXOS
REAL(8),DIMENSION(:), ALLOCATABLE ::  DMU,DELTA,F0,MUO,FUNC
INTEGER IMU,I

ALLOCATE( DMU(NUM),DELTA(NUM),F0(NUM),MUO(NUM),FUNC(NUM) )

IF (ARG_LENGTH(1)==1) THEN
! N.B. IF MINO IS SCALAR, MAXO IS VECTOR
MINOS=MINO(1)
DMU=MAXO-MINOS
ELSEIF (ARG_LENGTH(2)==1) THEN
! N.B. IF MAXO IS SCALAR, MINO IS VECTOR
MAXOS=MAXO(1)
DMU=MAXOS-MINO
ELSE
DMU=MAXO-MINO
END IF

DELTA=DMU/STEPS
F0=0.5d0*DELTA
INTEGR=0.0d0

DO IMU=1,STEPS
    IF (ARG_LENGTH(1)==1) THEN
        MUO=MINOS+F0+(IMU-1)*DELTA
    ELSE
        MUO=MINO+F0+(IMU-1)*DELTA
    END IF
    CALL INTEGRFI(XX,MUI,MUO,STEPS,FUNC,NUM)
    INTEGR=INTEGR+FUNC*DELTA
END DO

INTEGR=INTEGR/2d0

DEALLOCATE( DMU,DELTA,F0,MUO,FUNC )

END SUBROUTINE INTEGRMUO