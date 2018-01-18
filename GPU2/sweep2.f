      SUBROUTINE SWEEP2(I,IKS,DTCL,RCL)
*
*
*       Enforced KS search for massive wide binary.
*       -------------------------------------------
*
      INCLUDE 'common6.h'
      SAVE IT,IBIN
      DATA IT,IBIN /0,0/
*
*
*       Search neighbour list for KS candidates (STEP < DTCL).
      RX2 = 100.0
      JMIN = 0
      NNB1 = LIST(1,I) + 1
      DO 10 L = 2,NNB1
          J = LIST(L,I)
          IF (STEP(J).GT.DTCL) GO TO 10
          RIJ2 = 0.0
          DO 5 K = 1,3
              RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
    5     CONTINUE
          IF (RIJ2.LT.RX2) THEN
              RX2 = RIJ2
              JMIN = J
          END IF
   10 CONTINUE
*
*       Skip any close c.m./chain body (small STEP treated by IMPACT).
      IF (JMIN.EQ.0.OR.JMIN.GT.N) GO TO 30
      IF (NAME(I).LE.0.OR.NAME(JMIN).LE.0) GO TO 30
      IF (I.GT.N) GO TO 30
*
*       Form inverse semi-major axis.
      VIJ2 = 0.0
      DO 15 K = 1,3
          VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,JMIN))**2
   15 CONTINUE
      AINV = 2.0/SQRT(RX2) - VIJ2/(BODY(I) + BODY(JMIN))
*
*       Define KS binary candidates for augmented parameters (cf. sweep.f).
      IF (AINV.GT.0.0.AND.AINV.GT.1.0/RCL) THEN
          ICOMP = I
          JCOMP = JMIN
*       Skip rare case of chain & standard c.m. as binary component.
          IF (NAME(ICOMP).LE.0.OR.I.GT.N.OR.NSUB.GT.0) GO TO 30
          IKS = 1
          IT = IT + 1
          IF (IT.EQ.1) GO TO 30
          IBIN = IBIN + 1
          IT = 0
          IF (IBIN.LT.20) THEN
              WRITE (6,20)  NAME(I), NAME(JMIN), SQRT(RX2), 1.0/AINV
   20         FORMAT (' WIDE & MASSIVE KS    NMI NMJ RX A ',
     &                                       2I6,1P,2E10.2)
          END IF
      END IF
*
   30 RETURN
*
      END
