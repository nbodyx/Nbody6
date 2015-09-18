      SUBROUTINE CHPERT(GAMX)
*
*
*       Chain perturber selection.
*       --------------------------
*
      INCLUDE 'common6.h'
         REAL*8  M,MASS,MC
         PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
         COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &   XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      REAL*8  ACC(NMX3),FP(NMX3),RIX2(NMX),CMX(3),FCM(3)
*
*
*       Skip on zero perturbers.
      GAMX = 0.0
      IF (LISTC(1).EQ.0) GO TO 40
*
*       Determine global values of XC (not known after ABSORB or REDUCE).
      CALL XCPRED(0)
*
*       Initialize the external perturbations.
      NK = 3*NN
      DO 1 K = 1,NK
          ACC(K) = 0.0D0
    1 CONTINUE
*
*       Copy provisional perturber list for evaluation loop.
      NPC = LISTC(1) + 1
      DO 2 L = 2,NPC
          JPERT(L) = LISTC(L)
    2 CONTINUE
*
*       Determine indices of closest separation (safer than RINV).
      RX2 = 10.0
      DO 5 I = 1,NN-1
          DO 4 J = I+1,NN
              RIJ2 = 0.0
              DO 3 K = 1,3
                  RIJ2 = RIJ2 + (XC(K,I) - XC(K,J))**2
    3         CONTINUE
              IF (RIJ2.LT.RX2) THEN
                  RX2 = RIJ2
                  I1 = I
                  I2 = J
              END IF
    4     CONTINUE
    5 CONTINUE
*
*       Form global coordinates of binary c.m.
      ZMB = BODYC(I1) + BODYC(I2)
      DO 6 K = 1,3
          CMX(K) = (BODYC(I1)*XC(K,I1) + BODYC(I2)*XC(K,I2))/ZMB
    6 CONTINUE
*
*       Consider each provisional perturber in turn.
      NP = 1
      LX3 = 3*I1 - 3
      DO 20 L = 2,NPC
          J = JPERT(L)
          A1 = X(1,J) - CMX(1)
          A2 = X(2,J) - CMX(2)
          A3 = X(3,J) - CMX(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(J)/(RIJ2*SQRT(RIJ2))
          FCM(1) = A1*A6
          FCM(2) = A2*A6
          FCM(3) = A3*A6
*       Use one of the free locations (#I1) for summing c.m. force.
          ACC(LX3+1) = ACC(LX3+1) + FCM(1)
          ACC(LX3+2) = ACC(LX3+2) + FCM(2)
          ACC(LX3+3) = ACC(LX3+3) + FCM(3)
*       Include special case of binary (no other contributions below).
          IF (NN.EQ.2) THEN
*       Employ tidal approximation (small binary size).
              GAM = 2.0*A6*RX2*SQRT(RX2)/ZMB
              GAMX = MAX(GAM,GAMX)
              IF (GAM.GT.0.1*GMIN) THEN
                  NP = NP + 1
                  LISTC(NP) = J
              END IF
          END IF
          IK = -3
          ITIME = 0
*       Sum perturber contributions over each chain component (skip I1 & I2).
          DO 10 I = 1,NN
              IK = IK + 3
              IF (I.EQ.I1.OR.I.EQ.I2) GO TO 10
              A1 = X(1,J) - XC(1,I)
              A2 = X(2,J) - XC(2,I)
              A3 = X(3,J) - XC(3,I)
              RIJ2 = A1*A1 + A2*A2 + A3*A3
              A6 = BODY(J)/(RIJ2*SQRT(RIJ2))
              FP(IK+1) = A1*A6
              FP(IK+2) = A2*A6
              FP(IK+3) = A3*A6
              ACC(IK+1) = ACC(IK+1) + FP(IK+1)
              ACC(IK+2) = ACC(IK+2) + FP(IK+2)
              ACC(IK+3) = ACC(IK+3) + FP(IK+3)
*       Use standard perturbation test (GMIN) wrt c.m. for acceptance.
              DF2 = 0.0
              RIJ2 = 0.0
              DO 8 K = 1,3
                  DF = FP(IK+K) - FCM(K)
                  DF2 = DF2 + DF**2
                  RIJ2 = RIJ2 + (XC(K,I) - CMX(K))**2
    8         CONTINUE
              RIX2(I) = RIJ2
              GAM = SQRT(DF2)/((BODYC(I) + ZMB)/RIJ2)
*       Add accepted perturber to LISTC first time only.
              IF (GAM.GT.0.1*GMIN.AND.ITIME.EQ.0) THEN
                  ITIME = 1
                  NP = NP + 1
                  LISTC(NP) = J
              END IF
   10     CONTINUE
   20 CONTINUE
*
*       Evaluate the total relative perturbations and save maximum.
      DO 30 I = 1,NN
          IF (I.EQ.I1.OR.I.EQ.I2) GO TO 30
          LI = 3*(I - 1)
          DF2 = 0.0
          DO 25 K = 1,3
              DF = ACC(LI+K) - ACC(LX3+K)
              DF2 = DF2 + DF**2
   25     CONTINUE
*       Note that rejected perturbers are included in each final value.
          GAM = SQRT(DF2)/((BODYC(I) + ZMB)/RIX2(I))
          GAMX = MAX(GAM,GAMX)
   30 CONTINUE
*
*       Save new membership.
      LISTC(1) = NP - 1
*
   40 RETURN
*
      END
