      SUBROUTINE NEWTON(TD4,TD5,TD6,IP,DTU)
*
*
*       Newton-Raphson iteration.
*       -------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      PARAMETER (ONE30=1.0/30.0D0,ONE120=1.0/120.0D0,ZZ=1.0/720.0D0)
*
*
*       Include possible slow-down factor.
      IMOD = KSLOW(IP)
*     IF (IMOD.GT.1) THEN
          ZMOD = FLOAT(ISLOW(IMOD))
*     END IF
*
*       Include improved tolerance near small pericentre (R/a < 1.0D-02).
      IF (R(IP)*ABS(H(IP))/BODY(N+IP).LT.5.0D-03) THEN
          TOL = 1.0D-10
      ELSE
          TOL = 1.0D-08
      END IF
*
*       Note TD4 is 1/2 fourth derivative and TD5 is 1/4 fifth derivative.
      DO 10 I = 1,20   ! above 4 iterations is quite rare.
          FX = (((((ZZ*TD6*DTU + ONE30*TD5)*DTU + ONE12*TD4)*DTU +
     &             ONE6*TDOT3(IP))*DTU + 0.5*TDOT2(IP))*DTU + R(IP))*DTU
          FD = ((((ONE120*TD6*DTU + ONE6*TD5)*DTU + ONE3*TD4)*DTU +
     &             0.5*TDOT3(IP))*DTU + TDOT2(IP))*DTU + R(IP)
          IF (IMOD.GT.1) THEN
              FX = ZMOD*FX
              FD = ZMOD*FD
          END IF
          FX = FX - STEP(2*IP-1)
*      Include fast safety measures to prevent large values of FX/FD.
          TMP = FX/FD
          IF (TMP.LT.-0.2*DTU) THEN
              TMP = -0.2*DTU
          ELSE IF (TMP.GT.0.5*DTU) THEN
              TMP = 0.5*DTU
          END IF
          DTU = DTU - TMP
          IF (DTU.LT.0.0D0) THEN
              DTU = 0.0
          END IF
*      NB! Last value of DTU ensures actual tolerance < 1D-12 (TOL = 1D-10).
          IF (ABS(FX).LT.TOL.AND.DTU.GT.0.0D0) GO TO 30
   10 CONTINUE
*
      WRITE (6,20)  NSTEPU, FX, FD, DTU, DTAU(IP)
   20 FORMAT (' DANGER NEWTON!    # FX FD DTU DTAU ',I11,1P,4E10.2)
      STOP
*
   30 RETURN
*
      END
