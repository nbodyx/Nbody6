      SUBROUTINE JPRED2(I,IPRED,CMX,CMV)
*
*
*       Neighbour prediction of single particle.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  CMX(3),CMV(3)
*
*
*       Choose between copy and standard prediction.
      IF (IPRED.EQ.0) THEN
          DO 10 K = 1,3
              CMX(K) = X(K,I)
              CMV(K) = XDOT(K,I)
   10     CONTINUE
      ELSE
          S = TIME - T0(I)
          S1 = 1.5*S
          S2 = 2.0*S
          CMX(1) = ((FDOT(1,I)*S + F(1,I))*S + X0DOT(1,I))*S + X0(1,I)
          CMX(2) = ((FDOT(2,I)*S + F(2,I))*S + X0DOT(2,I))*S + X0(2,I)
          CMX(3) = ((FDOT(3,I)*S + F(3,I))*S + X0DOT(3,I))*S + X0(3,I)
          CMV(1) = (FDOT(1,I)*S1 + F(1,I))*S2 + X0DOT(1,I)
          CMV(2) = (FDOT(2,I)*S1 + F(2,I))*S2 + X0DOT(2,I)
          CMV(3) = (FDOT(3,I)*S1 + F(3,I))*S2 + X0DOT(3,I)
      END IF
*
      RETURN
*
      END
