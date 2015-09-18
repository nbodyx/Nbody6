      SUBROUTINE CPERTX(JP,KCASE,GAMX,GAMY)
*
*
*       External chain perturbation.
*       ----------------------------
*
      INCLUDE 'common6.h'
         REAL*8  M,MASS,MC
         PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
         COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &   XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
      REAL*8  ACC(NMX3),CMX(3),XJ(3)
*
*
*       Initialize the relative perturbations.
      GAMX = 0.0
      GAMY = 0.0
*
*       Define body #JP in local coordinates.
      DO 2 K = 1,3
          XJ(K) = X(K,JP) - X(K,ICH)
    2 CONTINUE
*
*       Determine indices of closest separation (safer than RINV).
      RX2 = 10.0
      DO 5 I = 1,NN-1
          DO 4 J = I+1,NN
              RIJ2 = 0.0
              DO 3 K = 1,3
                  RIJ2 = RIJ2 + (X4(K,I) - X4(K,J))**2
    3         CONTINUE
              IF (RIJ2.LT.RX2) THEN
                  RX2 = RIJ2
                  I1 = I
                  I2 = J
              END IF
    4     CONTINUE
    5 CONTINUE
*
*       Form local coordinates of binary c.m.
      ZMB = BODYC(I1) + BODYC(I2)
      DO 6 K = 1,3
          CMX(K) = (BODYC(I1)*X4(K,I1) + BODYC(I2)*X4(K,I2))/ZMB
    6 CONTINUE
*
*       Check for treating second case only.
      IF (KCASE.EQ.2) GO TO 50
*
*       Consider the c.m. and each member.
      LX3 = 3*I1 - 3
      A1 = XJ(1) - CMX(1)
      A2 = XJ(2) - CMX(2)
      A3 = XJ(3) - CMX(3)
      RIJ2 = A1*A1 + A2*A2 + A3*A3
      A6 = BODY(JP)/(RIJ2*SQRT(RIJ2))
*       Use one of the free locations (#I1) for c.m. force.
      ACC(LX3+1) = A1*A6
      ACC(LX3+2) = A2*A6
      ACC(LX3+3) = A3*A6
*       Include special case of triple (no other contributions below).
      IF (NN.EQ.3) THEN
*       Employ tidal approximation (small binary size).
          GAMX = 2.0*A6*RX2*SQRT(RX2)/ZMB
      END IF
      IK = -3
*       Obtain perturber contributions from other components (skip I1 & I2).
      DO 10 I = 1,NN
          IK = IK + 3
          IF (I.EQ.I1.OR.I.EQ.I2) GO TO 10
          A1 = XJ(1) - X4(1,I)
          A2 = XJ(2) - X4(2,I)
          A3 = XJ(3) - X4(3,I)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(JP)/(RIJ2*SQRT(RIJ2))
          ACC(IK+1) = A1*A6
          ACC(IK+2) = A2*A6
          ACC(IK+3) = A3*A6
   10 CONTINUE
*
*       Evaluate the relative perturbation and save maximum.
      DO 30 I = 1,NN
          IF (I.EQ.I1.OR.I.EQ.I2) GO TO 30
          LI = 3*(I - 1)
          DF2 = 0.0
          RIJ2 = 0.0
          DO 25 K = 1,3
              DF = ACC(LI+K) - ACC(LX3+K)
              DF2 = DF2 + DF**2
              RIJ2 = RIJ2 + (X4(K,I) - CMX(K))**2
   25     CONTINUE
          GAM = SQRT(DF2)/((BODYC(I) + ZMB)/RIJ2)
          GAMX = MAX(GAM,GAMX)
   30 CONTINUE
*
*       Exit on control index of 1.
      IF (KCASE.EQ.1) GO TO 100
*
*       Evaluate maximum perturbation on c.m. and body #JP.
   50 NNB = LISTC(1) 
      DO 60 L = 2,NNB+1
          J = LISTC(L)
          IF (J.EQ.JP) GO TO 60
          A1 = (X(1,J) - X(1,ICH)) - CMX(1)
          A2 = (X(2,J) - X(2,ICH)) - CMX(2)
          A3 = (X(3,J) - X(3,ICH)) - CMX(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(J)/(RIJ2*SQRT(RIJ2))
          ACC(1) = A1*A6
          ACC(2) = A2*A6
          ACC(3) = A3*A6
*
*       Subtract the contribution from body #JP.
          A1 = X(1,J) - X(1,JP)
          A2 = X(2,J) - X(2,JP)
          A3 = X(3,J) - X(3,JP)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(J)/(RIJ2*SQRT(RIJ2))
*       Form relative force and dimensionless perturbation.
          ACC(1) = ACC(1) - A1*A6
          ACC(2) = ACC(2) - A2*A6
          ACC(3) = ACC(3) - A3*A6
          DF2 = 0.0
          RIJ2 = 0.0
          DO 55 K = 1,3
              DF2 = DF2 + ACC(K)**2
              RIJ2 = RIJ2 + (XJ(K) - CMX(K))**2
   55     CONTINUE
          GAM = SQRT(DF2)/((BODY(JP) + ZMB)/RIJ2)
          GAMY = MAX(GAM,GAMY)
   60 CONTINUE
*
  100 RETURN
*
      END
