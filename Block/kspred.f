      SUBROUTINE KSPRED(IPAIR,I,BODYIN,DTU,UI,UIDOT,Q1,Q2,Q3,RDOT)
*
*
*       Prediction for KS regularization.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  UI(4),UIDOT(4),RDOT(3),A1(3,4)
      PARAMETER (ONE24=1.0/24.0D0,ONE120=1.0/120.0D0)
*
*
*       Check for stabilization of binaries (skip GAMMA > GMAX).
      IF (GAMMA(IPAIR).LT.GMAX) THEN
          A2 = 2.0D0*(UDOT(1,IPAIR)**2 + UDOT(2,IPAIR)**2 +
     &                UDOT(3,IPAIR)**2 + UDOT(4,IPAIR)**2) -
     &                BODY(I) - H(IPAIR)*R(IPAIR)
*       Include the stabilization term in the predicted force only.
          STAB = 0.2D0*A2*BODYIN/DTAU(IPAIR)
*       Note DTU may be small when called from KSPERT2 (hence full value).
      ELSE
          STAB = 0.0D0
      END IF
*
*       Predict U, UDOT & R to order FUDOT3.
      DO 20 K = 1,4
          FSTAB = FU(K,IPAIR) - STAB*UDOT(K,IPAIR)
          U4 = ONE24*FUDOT2(K,IPAIR)
          U5 = ONE120*FUDOT3(K,IPAIR)
          UI(K) = ((((U5*DTU + U4)*DTU + FUDOT(K,IPAIR))*DTU +
     &                FSTAB)*DTU + UDOT(K,IPAIR))*DTU + U0(K,IPAIR)
          UIDOT(K) = (((5.0*U5*DTU + 4.0*U4)*DTU +
     &          3.0*FUDOT(K,IPAIR))*DTU + 2.0*FSTAB)*DTU + UDOT(K,IPAIR)
   20 CONTINUE
*       Note predicted R as scalar for consistency with t' = R in newton.f.
      RI = UI(1)**2 + UI(2)**2 + UI(3)**2 + UI(4)**2
*
*       Form relative coordinates obtained from explicit KS transformation.
      Q1 = UI(1)**2 - UI(2)**2 - UI(3)**2 + UI(4)**2
      Q2 = UI(1)*UI(2) - UI(3)*UI(4)
      Q3 = UI(1)*UI(3) + UI(2)*UI(4)
      Q2 = Q2 + Q2
      Q3 = Q3 + Q3
*
*       Set current transformation matrix.
      CALL MATRIX(UI,A1)
*
*       Obtain relative velocities from KS transformation.
      RINV = 2.0D0/RI
      DO 30 L = 1,3
          RDOT(L) = 0.0D0
          DO 25 K = 1,4
              RDOT(L) = RDOT(L) + A1(L,K)*UIDOT(K)
   25     CONTINUE
          RDOT(L) = RDOT(L)*RINV
   30 CONTINUE
*
      RETURN
*
      END
