      SUBROUTINE SWAP
*
*
*       Randomized particle swapping.
*       -----------------------------
*
      INCLUDE 'common6.h'
      REAL*8  RAN2,SAVE(7)
*
*       Fisher-Yates-Knuth shuffling algorithm.
*       ---------------------------------------
*
*       Decide on swapping here (strongly recommended).
*     IF (N.GT.0) RETURN
*
      KDUM = IDUM1
      JMIN = 2*NBIN0 + 1
      JMAX = N
      DO 20 J = JMAX,JMIN+1,-1
          XR1 = RAN2(KDUM)
          I = JMIN + FLOOR((J-JMIN+1)*XR1)
          I = MIN(I,J)
          I = MAX(I,JMIN)
          IF (I.EQ.J) GO TO 20
          NAMI = NAME(I)
          SAVE(1) = BODY(I)
          DO 15 K = 1,3
              SAVE(K+1) = X(K,I)
              SAVE(K+4) = XDOT(K,I)
   15     CONTINUE
          NAME(I) = NAME(J)
          BODY(I) = BODY(J)
          DO 16 K = 1,3
              X(K,I) = X(K,J)
              XDOT(K,I) = XDOT(K,J)
              X0DOT(K,I) = XDOT(K,J)
   16     CONTINUE
          NAME(J) = NAMI
          BODY(J) = SAVE(1)
          DO 18 K = 1,3
              X(K,J) = SAVE(K+1)
              XDOT(K,J) = SAVE(K+4)
              X0DOT(K,J) = SAVE(K+4)
   18     CONTINUE
   20 CONTINUE
*
      RETURN
*
      END
