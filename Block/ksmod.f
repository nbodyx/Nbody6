      SUBROUTINE KSMOD(IPAIR,KMOD)
*
*
*       Modified KS motion.
*       -------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
*
*
*       Set current modification level.
      IMOD = KSLOW(IPAIR)
*
*       Check transition type (KMOD <= 1 denotes standard restart).
      IF (KMOD.GT.1) THEN
*       Determine provisional index for KS slow-down.
          DO 5 K = 2,10
              ISBIN = K - 1
              IF (ISLOW(K).GT.KMOD) GO TO 10
    5     CONTINUE
*       Restrict increase to two levels.
   10     ISBIN = MIN(ISBIN,IMOD+2)
*       See whether standard solution is called for.
          IF (ISBIN.EQ.1) GO TO 30
      ELSE
*       Include termination (IMOD > 1 & KMOD <= 1).
          ISBIN = 1
          GO TO 30
      END IF
*
*       Estimate time interval to reach largest permitted perturbation.
      GX = RANGE*GMIN
      IP = IPAIR
      CALL TPERT(IPAIR,GX,DT)
      IPAIR = IP
*
*       Evaluate the unmodified Kepler period.
      ICM = N + IPAIR
      SEMI = -0.5*BODY(ICM)/H(IPAIR)
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(ICM))
*
*       Reduce level if modification factor is too large.
      DO 20 K = 2,10
          IF (TK*FLOAT(ISLOW(ISBIN)).LT.DT.OR.ISBIN.EQ.1) GO TO 30
          ISBIN = ISBIN - 1
   20 CONTINUE
*
*       Exit if level is unchanged.
   30 IF (ISBIN.EQ.IMOD) GO TO 100
*
*       Ensure the Kepler period is defined for transition to standard case.
      IF (ISBIN.EQ.1) THEN
          ICM = N + IPAIR
          SEMI = -0.5*BODY(ICM)/H(IPAIR)
          TK = TWOPI*SEMI*SQRT(SEMI/BODY(ICM))
      END IF
*
*       Predict current coordinates & velocities of ICM.
      CALL XVPRED(ICM,0)
*
*       Set new KS level and increase restart counter (ISBIN > 1).
      ISBIN = MIN(ISBIN,8)
      KSLOW(IPAIR) = ISBIN
      IF (ISBIN.GT.1) NKSMOD = NKSMOD + 1
*
*       Perform perturbed restart of KS motion.
      CALL KSPOLY(IPAIR,ISBIN)
*
*       Ensure that apocentre criterion will fail after next step.
      TDOT2(IPAIR) = 0.0D0
*
*       Set indicator = -1 to skip perturber selection in routine KSINT.
      KMOD = -1
      IPHASE = 0
*
*       Enforce new KSLIST after small perturbation (CHAOS/SPIRAL case).
      IF (ISBIN.GE.6) KMOD = 0
*
  100 RETURN
*
      END
