      SUBROUTINE JPRED(I)
*
*
*       Neighbour prediction of single particle.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/PREDICT/ TPRED(NMAX)
*
*
*       Skip prediction if done at current time.
      IF (TPRED(I).EQ.TIME) GO TO 10
*
*       Adopt low order prediction for standard single particle.
      S = TIME - T0(I)
      S1 = 1.5*S
      S2 = 2.0*S
*       Note X may become corrected value X0 for zero interval.
      X(1,I) = ((FDOT(1,I)*S + F(1,I))*S + X0DOT(1,I))*S + X0(1,I)
      X(2,I) = ((FDOT(2,I)*S + F(2,I))*S + X0DOT(2,I))*S + X0(2,I)
      X(3,I) = ((FDOT(3,I)*S + F(3,I))*S + X0DOT(3,I))*S + X0(3,I)
      XDOT(1,I) = (FDOT(1,I)*S1 + F(1,I))*S2 + X0DOT(1,I)
      XDOT(2,I) = (FDOT(2,I)*S1 + F(2,I))*S2 + X0DOT(2,I)
      XDOT(3,I) = (FDOT(3,I)*S1 + F(3,I))*S2 + X0DOT(3,I)
*
*       Resolve the components of any perturbed pair.
      IF (I.GT.N) THEN
          JPAIR = I - N
          IF (LIST(1,2*JPAIR-1).GT.0) THEN
!$omp critical
              ZZ = 1.0
*       Distinguish between low and high-order prediction of U & UDOT.
              IF (GAMMA(JPAIR).GT.1.0D-04) ZZ = 0.0
              CALL KSRES2(JPAIR,J1,J2,ZZ)
!$omp end critical
          END IF
      END IF
*
      TPRED(I) = TIME
*
   10 RETURN
*
      END
