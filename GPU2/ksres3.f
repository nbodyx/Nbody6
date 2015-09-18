      SUBROUTINE KSRES3(J,J1,J2,IPRED,CMX,CMV,XK,VK)
*
*
*       Coordinates & velocities of KS pair (local arrays).
*       ---------------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      PARAMETER (ONE24=1.0/24.0D0)
      REAL*8  UI(4),RDOT(3),V(4),A1(3,4),XK(6),VK(6),CMX(3),CMV(3)
*
*
*       Resolve components of pair #J at curent time (according to GAMMA).
      J2 = J + J
      J1 = J2 - 1
      A2 = 1.0/R(J)
      A3 = A2*(TIME - T0(J1))
      IF (KSLOW(J).GT.1) THEN
          IMOD = KSLOW(J)
          A3 = A3/FLOAT(ISLOW(IMOD))
      END IF
*
*       Decide appropriate order for interpolation & prediction.
      IF (GAMMA(J).LT.1.0D-04) THEN
*       Convert physical interval to regularized time (second order only).
          DTU = (1.0 - 0.5D0*TDOT2(J)*A2*A3)*A3
          IF (ABS(DTU).GT.DTAU(J)) DTU = DTAU(J)
*
*       Predict regularized coordinates of distant pair to second order.
          DO 10 K = 1,4
              UI(K) = (FU(K,J)*DTU + UDOT(K,J))*DTU + U0(K,J)
              V(K) = (3.0*FUDOT(K,J)*DTU + 2.0*FU(K,J))*DTU + UDOT(K,J)
   10     CONTINUE
      ELSE
*        Expand regularized time interval to third order.
          A4 = 3.0D0*TDOT2(J)**2*A2 - TDOT3(J)
          DTU = ((ONE6*A4*A3 - 0.5D0*TDOT2(J))*A2*A3 + 1.0)*A3
*       Apply safety test near small pericentre.
          IF (DTU.GT.DTAU(J)) DTU = 0.8*DTAU(J)
*
*       Predict regularized coordinates to third order.
          DTU1 = ONE24*DTU
          DTU2 = ONE6*DTU
          DO 20 K = 1,4
              UI(K) = (((FUDOT2(K,J)*DTU1 + FUDOT(K,J))*DTU +
     &                           FU(K,J))*DTU + UDOT(K,J))*DTU + U0(K,J)
              V(K) = ((FUDOT2(K,J)*DTU2 + 3.0D0*FUDOT(K,J))*DTU +
     &                                    2.0D0*FU(K,J))*DTU + UDOT(K,J)
   20     CONTINUE
      END IF
*
*       Employ KS transformation.
      Q1 = UI(1)**2 - UI(2)**2 - UI(3)**2 + UI(4)**2
      Q2 = UI(1)*UI(2) - UI(3)*UI(4)
      Q3 = UI(1)*UI(3) + UI(2)*UI(4)
      Q2 = Q2 + Q2
      Q3 = Q3 + Q3
      I = N + J
*
*       Choose between copy c.m. or predict (cf. TPRED test in NBINTP).
      IF (IPRED.EQ.0) THEN
          DO 25 K = 1,3
              CMX(K) = X(K,I)
              CMV(K) = XDOT(K,I)
   25     CONTINUE
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
*       Set global coordinates of regularized components.
      A2 = BODY(J2)/BODY(I)
      XK(1) = CMX(1) + A2*Q1
      XK(2) = CMX(2) + A2*Q2
      XK(3) = CMX(3) + A2*Q3
      XK(4) = XK(1) - Q1
      XK(5) = XK(2) - Q2
      XK(6) = XK(3) - Q3
*
*       Set current transformation matrix and two-body separation.
      CALL MATRIX(UI,A1)
      RI = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
      RINV = 2.0/RI
*
*       Obtain relative velocities from KS transformation.
      DO 40 L = 1,3
          RDOT(L) = 0.0D0
          DO 30 K = 1,4
              RDOT(L) = RDOT(L) + A1(L,K)*V(K)*RINV
   30     CONTINUE
   40 CONTINUE
*
*       Set global velocities of KS components.
      VK(1) = CMV(1) + A2*RDOT(1)
      VK(2) = CMV(2) + A2*RDOT(2)
      VK(3) = CMV(3) + A2*RDOT(3)
      VK(4) = VK(1) - RDOT(1)
      VK(5) = VK(2) - RDOT(2)
      VK(6) = VK(3) - RDOT(3)
*
      RETURN
*
      END
