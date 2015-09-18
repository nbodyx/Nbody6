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
      DTU0 = DTAU(IP)
*       Note TD4 is 1/2 fourth derivative and TD5 is 1/4 fifth derivative.
      DO 10 I = 1,15   ! above 10 iterations is quite rare.
          FX = (((((ZZ*TD6*DTU + ONE30*TD5)*DTU + ONE12*TD4)*DTU +
     &             ONE6*TDOT3(IP))*DTU + 0.5*TDOT2(IP))*DTU + R(IP))*DTU
          FD = ((((ONE120*TD6*DTU + ONE6*TD5)*DTU + ONE3*TD4)*DTU +
     &             0.5*TDOT3(IP))*DTU + TDOT2(IP))*DTU + R(IP)
          IF (IMOD.GT.1) THEN
              FX = ZMOD*FX
              FD = ZMOD*FD
          END IF
          FX = FX - STEP(2*IP-1)
          DTU = DTU - FX/FD
          IF (DTU.LT.0.0D0) THEN
              DTU = 0.0
          ELSE IF (DTU.GT.DTU0) THEN
              DTU = DTU0
          END IF
*      Note last value of DTU ensures the actual tolerance is < 1D-12.
          IF (ABS(FX).LT.TOL.AND.DTU.GT.0.0D0) GO TO 30
   10 CONTINUE
*
      WRITE (6,20)  NSTEPU, FX, FD, DTU, DTU0
   20 FORMAT (' DANGER NEWTON!    # FX FD DTU DTU0 ',I11,1P,4E10.2)
      STOP
*
   30 RETURN
*
      END
