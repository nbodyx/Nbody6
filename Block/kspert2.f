      SUBROUTINE KSPERT2(I1,I,NNB0,BODYIN,Q1,Q2,Q3,RDOT,FP,FD)
*
*
*       Perturbation on parallel KS pair.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)
      REAL*8  XG(3),XGDOT(3),FM(3),FMD(3),FS(3),FSD(3)
      REAL*8  XI(6),VI(6),FP(6),FD(6),TF(3),TD(3),XK(6),VK(6),RDOT(3)
      REAL*8  UI(4),UIDOT(4),XREL(3),VREL(3)
      SAVE IC,IX
      DATA IC,IX /0,0/
      SAVE HIP
      DATA HIP /0.0D0/
*
*
*       Restore < 0 KS index and use ISKIP > 0 to avoid PN energy correction.
      IF (I1.LT.0) THEN
          ISKIP = 1
          I1 = -I1
      ELSE
          ISKIP = 0
      END IF
*
*       Initialize the perturbing force & first derivative.
      DO 10 K = 1,6
          FP(K) = 0.0D0
          FD(K) = 0.0D0
   10 CONTINUE
*
*       Predict current c.m.
      S = TIME - T0(I)
      S1 = 1.5*S
      S2 = 2.0*S
      XI(1) = ((FDOT(1,I)*S + F(1,I))*S + X0DOT(1,I))*S + X0(1,I)
      XI(2) = ((FDOT(2,I)*S + F(2,I))*S + X0DOT(2,I))*S + X0(2,I)
      XI(3) = ((FDOT(3,I)*S + F(3,I))*S + X0DOT(3,I))*S + X0(3,I)
      VI(1) = (FDOT(1,I)*S1 + F(1,I))*S2 + X0DOT(1,I)
      VI(2) = (FDOT(2,I)*S1 + F(2,I))*S2 + X0DOT(2,I)
      VI(3) = (FDOT(3,I)*S1 + F(3,I))*S2 + X0DOT(3,I)
*       Set global coordinates of regularized components.
      A2 = BODY(I1+1)*BODYIN
      XI(1) = XI(1) + A2*Q1
      XI(2) = XI(2) + A2*Q2
      XI(3) = XI(3) + A2*Q3
      XI(4) = XI(1) - Q1
      XI(5) = XI(2) - Q2
      XI(6) = XI(3) - Q3
*       Form global velocities of KS components.
      VI(1) = VI(1) + A2*RDOT(1)
      VI(2) = VI(2) + A2*RDOT(2)
      VI(3) = VI(3) + A2*RDOT(3)
      VI(4) = VI(1) - RDOT(1)
      VI(5) = VI(2) - RDOT(2)
      VI(6) = VI(3) - RDOT(3)
*
*       Determine index of the last single perturber.
      NNB2 = NNB0 + 1
   15 IF (LIST(NNB2,I1).LE.N) GO TO 20
      NNB2 = NNB2 - 1
      IF (NNB2.GT.1) GO TO 15
*       Include special case of only c.m. perturbers.
      GO TO 30
*
*       Obtain the perturbation from single particles.
   20 DO 25 L = 2,NNB2
          J = LIST(L,I1)
          S = TIME - T0(J)
          S1 = 1.5*S
          S2 = 2.0*S
*       Predict each perturber in turn.
          XK(1) = ((FDOT(1,J)*S + F(1,J))*S + X0DOT(1,J))*S + X0(1,J)
          XK(2) = ((FDOT(2,J)*S + F(2,J))*S + X0DOT(2,J))*S + X0(2,J)
          XK(3) = ((FDOT(3,J)*S + F(3,J))*S + X0DOT(3,J))*S + X0(3,J)
          VK(1) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
          VK(2) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
          VK(3) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
*
*       Form perturbation on first component.
          A1 = XK(1) - XI(1)
          A2 = XK(2) - XI(2)
          A3 = XK(3) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(J)/(RIJ2*SQRT(RIJ2))
          FP(1) = FP(1) + A1*A6
          FP(2) = FP(2) + A2*A6
          FP(3) = FP(3) + A3*A6
          V1 = VK(1) - VI(1)
          V2 = VK(2) - VI(2)
          V3 = VK(3) - VI(3)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
          FD(1) = FD(1) + (V1 - A1*A9)*A6
          FD(2) = FD(2) + (V2 - A2*A9)*A6
          FD(3) = FD(3) + (V3 - A3*A9)*A6
*
*       Form perturbation on second component.
          A1 = XK(1) - XI(4)
          A2 = XK(2) - XI(5)
          A3 = XK(3) - XI(6)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(J)/(RIJ2*SQRT(RIJ2))
          FP(4) = FP(4) + A1*A6
          FP(5) = FP(5) + A2*A6
          FP(6) = FP(6) + A3*A6
          V1 = VK(1) - VI(4)
          V2 = VK(2) - VI(5)
          V3 = VK(3) - VI(6)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
          FD(4) = FD(4) + (V1 - A1*A9)*A6
          FD(5) = FD(5) + (V2 - A2*A9)*A6
          FD(6) = FD(6) + (V3 - A3*A9)*A6
   25 CONTINUE
*
*       See whether to include any remaining c.m. perturbers.
      IF (NNB2.GT.NNB0) GO TO 40
*
*       Loop over binary perturbers.
   30 DO 35 L = NNB2+1,NNB0+1
          K = LIST(L,I1)
          S = TIME - T0(K)
          S1 = 1.5*S
          S2 = 2.0*S
*       Predict c.m. coordinates and velocity on the fly (thread-safe).
          XK(1) = ((FDOT(1,K)*S + F(1,K))*S + X0DOT(1,K))*S + X0(1,K)
          XK(2) = ((FDOT(2,K)*S + F(2,K))*S + X0DOT(2,K))*S + X0(2,K)
          XK(3) = ((FDOT(3,K)*S + F(3,K))*S + X0DOT(3,K))*S + X0(3,K)
          VK(1) = (FDOT(1,K)*S1 + F(1,K))*S2 + X0DOT(1,K)
          VK(2) = (FDOT(2,K)*S1 + F(2,K))*S2 + X0DOT(2,K)
          VK(3) = (FDOT(3,K)*S1 + F(3,K))*S2 + X0DOT(3,K)
          J = K - N
          J1 = 2*J - 1
          IF (LIST(1,J1).EQ.0) GO TO 32
*
*       Resolve the binary components in a consistent way (not by RESOLV).
          BODYJN = 1.0/BODY(K)
          IF (T0(J1).EQ.TIME) THEN
              DTU = 0.0
          ELSE
              DTU = DTAU(J)*(TIME - T0(J1))/STEP(J1)
          END IF
          CALL KSPRED(J,K,BODYJN,DTU,UI,UIDOT,Q1,Q2,Q3,RDOT)
*
*       Set global coordinates of regularized components.
          A2 = BODY(J1+1)*BODYJN
          XK(1) = XK(1) + A2*Q1
          XK(2) = XK(2) + A2*Q2
          XK(3) = XK(3) + A2*Q3
          XK(4) = XK(1) - Q1
          XK(5) = XK(2) - Q2
          XK(6) = XK(3) - Q3
*       Form global velocities of KS components.
          VK(1) = VK(1) + A2*RDOT(1)
          VK(2) = VK(2) + A2*RDOT(2)
          VK(3) = VK(3) + A2*RDOT(3)
          VK(4) = VK(1) - RDOT(1)
          VK(5) = VK(2) - RDOT(2)
          VK(6) = VK(3) - RDOT(3)
*
*       Obtain perturbation on first component.
          K = J1
   32     A1 = XK(1) - XI(1)
          A2 = XK(2) - XI(2)
          A3 = XK(3) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          V1 = VK(1) - VI(1)
          V2 = VK(2) - VI(2)
          V3 = VK(3) - VI(3)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
          A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(1) = FP(1) + A1*A6
          FP(2) = FP(2) + A2*A6
          FP(3) = FP(3) + A3*A6
          FD(1) = FD(1) + (V1 - A1*A9)*A6
          FD(2) = FD(2) + (V2 - A2*A9)*A6
          FD(3) = FD(3) + (V3 - A3*A9)*A6
*
*       Form perturbation on second component.
          A1 = XK(1) - XI(4)
          A2 = XK(2) - XI(5)
          A3 = XK(3) - XI(6)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(4) = FP(4) + A1*A6
          FP(5) = FP(5) + A2*A6
          FP(6) = FP(6) + A3*A6
          V1 = VK(1) - VI(4)
          V2 = VK(2) - VI(5)
          V3 = VK(3) - VI(6)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
          FD(4) = FD(4) + (V1 - A1*A9)*A6
          FD(5) = FD(5) + (V2 - A2*A9)*A6
          FD(6) = FD(6) + (V3 - A3*A9)*A6
*
*       Check for individual component summation.
          IF (K.EQ.J1) THEN
*       Copy second KS component to first for simplicity.
              DO 34 KK = 1,3
                  XK(KK) = XK(KK+3)
                  VK(KK) = VK(KK+3)
   34         CONTINUE
              K = K + 1
              GO TO 32
          END IF
   35 CONTINUE
*
*       Check perturbation correction due to regularized chain.
   40 IF (NCH.GT.0) THEN
*!$omp critical
          DO 45 L = 2,NNB2
*       Note ICH <= N is ensured.
              J = LIST(L,I1)
              IF (J.GT.ICH) GO TO 50
              IF (J.EQ.ICH) THEN
                  J1 = I1
                  CALL FCHAIN(J1,0,XI(1),VI(1),FP(1),FD(1))
                  J1 = J1 + 1
                  CALL FCHAIN(J1,1,XI(4),VI(4),FP(4),FD(4))
                  GO TO 50
              END IF
   45     CONTINUE
*!$omp end critical
      END IF 
*
*       Set the relative perturbing force and first derivative.
   50 DO 55 K = 1,3
          FP(K) = FP(K) - FP(K+3)
          FD(K) = FD(K) - FD(K+3)
          TF(K) = 0.0D0
          TD(K) = 0.0D0
   55 CONTINUE
*
*       See whether the linearized perturbation should be included.
      IF (KZ(14).GT.0.AND.KZ(14).LT.3) THEN
          Q1 = XI(1) - XI(4)
          Q3 = XI(3) - XI(6)
          CALL XTRNLP(Q1,Q3,TF)
*
*       Use same formalism for the first derivative (omit Coriolis force).
          VX = VI(1) - VI(4)
          VZ = VI(3) - VI(6)
          CALL XTRNLP(VX,VZ,TD)
          DO 60 K = 1,3
             FP(K) = FP(K) + TF(K)
             FD(K) = FD(K) + TD(K)
   60     CONTINUE
      END IF
*
*       Check optional Plummer potential.
      IF (MP.GT.0.0) THEN
          RI2 = AP2
          RRDOT = 0.0
*       Form one central distance and scalar product of relative motion.
          DO 65 K = 1,3
              RI2 = RI2 + XI(K)**2
              RRDOT = RRDOT + (XI(K) - XI(K+3))*(VI(K) - VI(K+3))
   65     CONTINUE
          ZF = 1.0/RI2
*       Write current mass inside RI as MP*R3*ZF^{3/2} (Heggie & Hut p.73).
          FMP = MP*ZF*SQRT(ZF)
          DO 70 K = 1,3
              XREL(K) = XI(K) - XI(K+3)
              VREL(K) = VI(K) - VI(K+3)
              FP(K) = FP(K) - XREL(K)*FMP
              FD(K) = FD(K) - (VREL(K) - 3.0*RRDOT*ZF*XREL(K))*FMP
   70     CONTINUE
      END IF
*
*       Include Post-Newtonian terms in the perturbation.
      IF (KZ(49).GT.0.AND.MAX(KSTAR(I1),KSTAR(I1+1)).GE.13) THEN
          IP = KVEC(I1)
          SEMI = -0.5*BODY(I)/H(IP)
          ECC2 = (1.0 - R(IP)/SEMI)**2 + TDOT2(IP)**2/(BODY(I)*SEMI)
          DW = 3.0*TWOPI*BODY(I)/(SEMI*CLIGHT**2*(1.0 - ECC2))
*         IF (R(IP).LT.1000.0*RZ) THEN
*       Adopt Einstein shift with constant elements instead of a variable R.
          IF (DW.GT.1.0D-04) THEN
              DH = H(IP) - HIP
*       Skip the energy correction after significant change (new KS).
              IF (ABS(DH).GT.4.0D-04*ABS(H(IP)).AND.ISKIP.EQ.0) THEN
                  IC = 0
                  HIP = H(IP)
              ELSE
                  IC = IC + 1
              END IF
              DO 75 K = 1,3
                  XREL(K) = XI(K) - XI(K+3)
                  VREL(K) = VI(K) - VI(K+3)
                  TF(K) = 0.0
                  TD(K) = 0.0
   75         CONTINUE
              DT = STEP(I1)
*       Obtain the PN perturbation.
              CALL PNPERT2(BODY(I1),BODY(I1+1),XREL,VREL,TF,TD,DT,IC,
     &                     ISKIP,CLIGHT,DE)
*       Accumulate PN energy loss in ECOLL after call from KSINT and IC > 0.
              IF (IC.GT.0.AND.ISKIP.EQ.0) THEN
                  ECOLL = ECOLL - DE
              END IF
*       Combine the pertubations.
              DO 80 K = 1,3
                  FP(K) = FP(K) + TF(K)
                  FD(K) = FD(K) + TD(K)
   80         CONTINUE
*       Signal termination for pericentre < 10*RSCH (experimental but suppressed).
              PMIN = SEMI*(1.0 - SQRT(ECC2))
*             RSCH = 6.0*BODY(I)/CLIGHT**2
*             IF (PMIN.LT.10000.0*RSCH.AND.IX.GT.1000) IPHASE = -1
*             IF (DW.GT.1.0D-02.AND.IX.GT.1000) IPHASE = -1
              IX = IX + 1
              IF (MOD(IX,100).EQ.0) THEN
                  WRITE (50,90)  TIME+TOFF, SQRT(ECC2), SEMI, DW, DE,
     &                           ECOLL, PMIN
   90             FORMAT (' PNPERT    T ECC A DW DE ECOLL PM ',
     &                                F10.4,F8.4,1P,5E10.2)
                  CALL FLUSH(50)
              END IF
          END IF
      END IF
*
*       Include possible contribution from central point-mass.
      IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
          DO 95 K = 1,3
              XG(K) = RG(K) + XI(K)   ! Note XG used for each KS component.
              XGDOT(K) = VG(K) + VI(K)
   95     CONTINUE
          CALL FNUC(XG,XGDOT,FS,FSD)
          DO 100 K = 1,3
              XG(K) = RG(K) + XI(K+3)
              XGDOT(K) = VG(K) + VI(K+3)
  100     CONTINUE
          CALL FNUC(XG,XGDOT,FM,FMD)
          DO 105 K = 1,3
              FP(K) = FP(K) + (FS(K) - FM(K))
              FD(K) = FD(K) + (FSD(K) - FMD(K))
  105     CONTINUE
      END IF
*
      RETURN
*
      END
