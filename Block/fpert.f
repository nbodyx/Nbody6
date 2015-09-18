      SUBROUTINE FPERT(I,J,NP,PERT)
*
*
*       Perturbing force on dominant two-body system.
*       ---------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  FP(3)
*
*
      DO 1 K = 1,3
          FP(K) = 0.0D0
    1 CONTINUE
*
*       Obtain perturbation on body #I & J due to NP members of JLIST.
      RCRIT2 = 1.0D+04*RMIN2
      ITER = 0
    2 DO 10 KK = 1,NP
          K = JLIST(KK)
          IF (K.EQ.I.OR.K.EQ.J) GO TO 10
          L = I
    5     A1 = X(1,K) - X(1,L)
          A2 = X(2,K) - X(2,L)
          A3 = X(3,K) - X(3,L)
          RIJ2 = A1**2 + A2**2 + A3**2
*       Restrict evaluation RCRIT2 (neighbour list may be used).
          IF (L.EQ.I.AND.RIJ2.GT.RCRIT2) GO TO 10
          A4 = BODY(K)/(RIJ2*SQRT(RIJ2))
          IF (L.EQ.J) A4 = -A4
          FP(1) = FP(1) + A1*A4
          FP(2) = FP(2) + A2*A4
          FP(3) = FP(3) + A3*A4
          IF (L.EQ.I) THEN
              L = J
              GO TO 5
          END IF
   10 CONTINUE
*
*       Increase search distance on zero contributions (< 10 times).
      IF (FP(1).EQ.0.0D0) THEN
          RCRIT2 = 4.0*RCRIT2
          ITER = ITER + 1
          IF (ITER.LT.10) GO TO 2
      END IF
*
      PERT = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)
*
      RETURN
*
      END
