      SUBROUTINE SPINUP(IPAIR,ITERM)
*
*
*       Spin synchronization of hierarchical binary.
*       -------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      REAL*8  M1,M2,MX,MB
      DATA  rg2 /0.1/
      EXTERNAL RL
*
*
*       Set basic pointers.
      I = N + IPAIR
      I1 = 2*IPAIR - 1
      ITERM = 0
      IM = 0
      DO 1 K = 1,NMERGE
         IF (NAMEM(K).EQ.NAME(I)) IM = K
    1 CONTINUE
*
*       Find location of the secondary (i.e. ghost).
      CALL FINDJ(I,I2,IM)
      IF (I2.LT.0) GO TO 50
      M1 = CM(1,IM)
      M2 = CM(2,IM)
      MB = M1 + M2
      SEMI0 = -0.5*MB/HM(IM)
      SEMI1 = SEMI0
*
*       Obtain current separation and t'' = 2*U*U'.
      RB = 0.0
      TD2 = 0.0
      V20 = 0.0
      DO 5 K = 1,4
          RB = RB + UM(K,IM)**2
          TD2 = TD2 + 2.0*UM(K,IM)*UMDOT(K,IM)
          V20 = V20 + UMDOT(K,IM)**2
    5 CONTINUE
*
*       Evaluate inner eccentricity and exit non-circular orbit.
      ECC2 = (1.0 - RB/SEMI0)**2 + TD2**2/(MB*SEMI0)
      ECC = SQRT(ECC2)
      IF (ECC.GT.0.1) THEN
          GO TO 50
      END IF
*
*       Determine largest angular momentum and Roche radius.
      IF (SPIN(I2).GT.SPIN(I1)) THEN
          RX = RADIUS(I2)
          MX = M2
          J1 = I2
          Q = M2/M1
          RL1 = RL(Q)*SEMI0
      ELSE
          RX = RADIUS(I1)
          MX = M1
          J1 = I1
          Q = M1/M2
          RL1 = RL(Q)*SEMI0
      END IF
*
      IF (RX.GT.RL1) THEN
          ITERM = 1
          GO TO 50
      END IF
*
*       Solve for joint angular frequency for star and binary by iteration.
      ZMU = M1*M2/MB
      AMB0 = ZMU*SQRT(MB*SEMI0)
      SPIN0 = SPIN(J1)
      DO 10 L = 1,10
          W = SPIN(J1)/(rg2*MX*RX**2)
          WB = SQRT(MB/SEMI0**3)
          SPIN(J1) = SPIN(J1)*(WB/W)
          FAC = (AMB0 + SPIN0 - SPIN(J1))**2
*       Obtain new semi-major axis from angular momentum conservation.
          SEMI = FAC/(ZMU**2*MB)
*       Copy new value for continued iteration.
          DW = (WB - W)/WB
          ITER = L
          IF (ABS(DW).LT.1.0D-03) THEN
*             WRITE (6,8)  L, SEMI, WB, W, DW
              GO TO 15
          END IF
          SEMI0 = SEMI
*         WRITE (6,8)  L, SEMI, WB, W, DW
*   8     FORMAT (' ITERATE    L A WB W DW/W ',I4,1P,E12.4,3E10.2)
   10 CONTINUE
*
*       Check minimum size.
   15 IF (SEMI.LT.RX) THEN
          SEMI = RX
      END IF
*
*       Update binding energy.
      HI = HM(IM)
      HM(IM) = -0.5*MB/SEMI
*
*       Correct EMERGE & ECOLL (consistency; no net effect).
      DECORR = ZMU*(HI - HM(IM))
      EMERGE = EMERGE - DECORR
      ECOLL = ECOLL + DECORR
      EGRAV = EGRAV + DECORR
*
*       Specify KS coordinate & velocity scaling factors at general point.
      C2 = SQRT(SEMI/SEMI1)
      V2 = 0.5*(MB + HM(IM)*RB*(SEMI/SEMI1))
      C1 = SQRT(V2/V20)
*
*       Re-scale KS variables to new energy with constant eccentricity.
      RB = 0.0D0
      DO 20 K = 1,4
          UM(K,IM) = C2*UM(K,IM)
          UMDOT(K,IM) = C1*UMDOT(K,IM)
*         RB = RB + UM(K,IM)**2
   20 CONTINUE
*
*       Rectify the hierarchical KS variables.
      CALL HIRECT(IM)
*
*       Set new look-up time slightly ahead of components.
      TEV0(I) = TEV(I)
      TEV(I) = TEV(J1) + 0.01
*
      WRITE (6,25)  TIME+TOFF, NAME(I1), KSTAR(J1), ITER, SEMI, W, DW
   25 FORMAT (' SPINUP    T NM K* IT A W DW/W ',F9.3,I6,2I4,1P,3E10.2)
*
*     WRITE (6,30)  KSTAR(J1), KSTAR(I), MX*SMU, RX, RL1,
*    &              TEV(I1), TEV(I2)
*  30 FORMAT (' SPINCHECK   K*1 K*I MX RX RL TEV ',2I4,F7.2,1P,4E10.2)
*
   50 RETURN
*
      END
