SUBROUTINE EPSSOIL_WS80_TUNED_CALVET(MV,TD,FGHz,RHOKG,SAND,SILT,CLAY,EPS_SOIL)


!     PARAMETERS:
!     MV - SOIL VOLUMETRIC MOISTURE[FRAC]
!     TD - SOIL TEMPERATURE [degC]
!     FGHZ - FREQUENCY [GHZ]
!     RHOKG - SOIL BULK DENSITY [KG/M3]
!     SAND - SAND CONTENT[%]
!     SILT - SILT CONTENT[%]
!     CLAY - CLAY CONTENT[%]

IMPLICIT NONE

REAL(8), INTENT(IN) :: MV,TD,FGHZ,RHOKG,SAND,SILT,CLAY
COMPLEX(8),INTENT(OUT) :: EPS_SOIL

REAL(8) :: RHO_R,RHO_S,WP,GAMMA,WT
COMPLEX(8) :: EPS_ROCK, EPS_ICE, EPS_FW, EPS_X


!PARAMETER CONVERSION FIRST
RHO_S=RHOKG/1000.0d0

!CONSTANTS
!EPS_R: PERMITTIVITY OF ROCK
!EPS_I: PERMITTIVITY OF ICE
RHO_R=2.65d0
EPS_ROCK=COMPLEX(5.5d0,0.2d0)
EPS_ICE=COMPLEX(3.5d0,0.1d0)


!BEGIN
WP = 0.06774d0 - 0.064d0* SAND/100.0d0 + 0.478d0 *CLAY/100.0d0

IF (FGHZ.LE.15.0d0)then
    gamma = -0.57d0*WP + 0.481d0
    Wt = 0.49d0*WP + 0.165d0
ELSE
    gamma=0.60d0
    Wt=0.17d0

    IF(FGHZ.LT.37d0 .AND. FGHZ.GT.35d0)THEN
        GAMMA=0.81d0
    END IF

    IF(FGHZ.LT.91d0.AND.FGHZ.GT.88d0)THEN
        GAMMA=1.19d0
    END IF
END IF

!WATER PERMITTIVITY
CALL EPSW(FGHZ,TD+273.15d0,EPS_FW)

!CALCUALTE
!EPS_X IS THE EQUIVALENT BOUND WATER PERMITTIVITY THOUGHT BY WANG
IF(MV.LT.WT)THEN
    EPS_X=EPS_ICE + GAMMA *MV/WT * (EPS_FW-EPS_ICE)
    EPS_SOIL = 1.0d0 + RHO_S/RHO_R * (EPS_ROCK-1.0d0) + MV*EPS_X - MV
ELSE
    EPS_X=EPS_ICE + GAMMA * (EPS_FW - EPS_ICE)
    EPS_SOIL = 1.0d0 + RHO_S/RHO_R * (EPS_ROCK-1.0d0) + WT*EPS_X - MV + (MV-WT)*EPS_FW
END IF

END SUBROUTINE EPSSOIL_WS80_TUNED_CALVET