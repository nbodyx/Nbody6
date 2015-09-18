      SUBROUTINE KSINIT
*
*
*       Initialization of KS regularization.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/FSAVE/ SAVEIT(6)
      REAL*8  Q(3),RDOT(3),UI(4),VI(4),A1(3,4)
      REAL*8  A(9),F1(3),F1DOT(3)
      INTEGER IPIPE(9)
      SAVE  IPIPE
      DATA  IPIPE /9*0/
*
*
*       Set new global indices of the components and current pair index.
      ICOMP = 2*NPAIRS - 1
      JCOMP = ICOMP + 1
      IPAIR = NPAIRS
*
*       Add body #N in case the only neighbour was removed in KSREG.
      IF (LIST(1,ICOMP).EQ.0) THEN
          LIST(2,ICOMP) = N
          LIST(2,JCOMP) = N
          LIST(1,ICOMP) = 1
          LIST(1,JCOMP) = 1
      END IF
*
*       Specify mass, neighbour radius, name, radius & type for new c.m.
      BODY(NTOT) = BODY(ICOMP) + BODY(JCOMP)
      RS(NTOT) = RS(ICOMP)
      NAME(NTOT) = NZERO + NAME(ICOMP)
      RADIUS(NTOT) = 0.0
      TEV(NTOT) = 1.0E+10
      TEV0(NTOT) = 1.0E+10
      BODY0(NTOT) = BODY(NTOT)
      EPOCH(NTOT) = TIME*TSTAR
      KSTAR(NTOT) = 0
      IMOD = 1
*
*       Define c.m. coordinates & velocities and set XDOT for components.
      DO 10 K = 1,3
          X(K,NTOT) = (BODY(ICOMP)*X(K,ICOMP) + BODY(JCOMP)*X(K,JCOMP))/
     &                                                        BODY(NTOT)
          X0DOT(K,NTOT) = (BODY(ICOMP)*X0DOT(K,ICOMP) + BODY(JCOMP)*
     &                                        X0DOT(K,JCOMP))/BODY(NTOT)
          XDOT(K,NTOT) = X0DOT(K,NTOT)
          XDOT(K,ICOMP) = X0DOT(K,ICOMP)
          XDOT(K,JCOMP) = X0DOT(K,JCOMP)
          X0(K,NTOT) = X(K,NTOT)
      X0(K,JCOMP) = X(K,JCOMP)
   10 CONTINUE
*
*       Obtain force polynomial for c.m. with components ICOMP & JCOMP.
      NNB = LIST(1,ICOMP)
*
*       Predict current coordinates & velocities for the neighbours.
      CALL XVPRED(ICOMP,NNB)
*
*       Choose between full FPOLY1/2 or just neighbours (skip non-standard).  
      IF (N.LT.5000.OR.IPHASE.NE.1) THEN
*
*       Obtain new polynomials & steps (first F & FDOT, then F2DOT & F3DOT).
          CALL FPOLY1(ICOMP,JCOMP,1)
          CALL FPOLY2(NTOT,NTOT,1)
*
      ELSE
*
*       Treat each component in turn (initialize, then evaluate F & FDOT).
          I = ICOMP
          ICM = NTOT
  110     DO 115 K = 1,3
              FI(K,I) = 0.0
              D1(K,I) = 0.0
  115     CONTINUE
*
*       Obtain irregular force & first derivative for body #I.
          KDUM = 0
          NNB = LIST(1,ICM)
*       Loop over neighbours only using c.m. list (cf. FPOLY1).
          DO 140 L = 2,NNB+1
              J = LIST(L,ICM)
              IF (J.GT.N) THEN
                  JPAIR = J - N
*       Use c.m. approximation for unperturbed binary.
                  IF (LIST(1,2*JPAIR-1).GT.0) THEN
                      KDUM = 2*JPAIR - 1
                      J = KDUM
                  END IF
              END IF
*
  120         DO 125 K = 1,3
                  A(K) = X(K,J) - X(K,I)
                  A(K+3) = XDOT(K,J) - XDOT(K,I)
  125         CONTINUE
*
              A(7) = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
              A(8) = BODY(J)*A(7)*SQRT(A(7))
              A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
*
*       Accumulate irregular force and first derivative.
              DO 130 K = 1,3
                  F1(K) = A(K)*A(8)
                  F1DOT(K) = (A(K+3) - A(K)*A(9))*A(8)
                  FI(K,I) = FI(K,I) + F1(K)
                  D1(K,I) = D1(K,I) + F1DOT(K)
  130         CONTINUE
*
*       Check for KS component.
              IF (J.EQ.KDUM) THEN
                  J = J + 1
                  GO TO 120
              END IF
  140     CONTINUE
*
*       Check option for external force (note FR already done).
          IF (KZ(14).GT.0) THEN
              CALL XTRNLD(I,I,1)
          END IF
*
*       Treat second component in the same way.
          IF (I.EQ.ICOMP) THEN
              I = JCOMP
              GO TO 110
          END IF
*
*       Form c.m. force and derivative from mass-weighted terms.
          FIRR = 0.0
          FDIRR = 0.0
          DO 150 K = 1,3
              FI(K,ICM) = (BODY(ICOMP)*FI(K,ICOMP) +
     &                     BODY(JCOMP)*FI(K,JCOMP))/BODY(ICM)
              D1(K,ICM) = (BODY(ICOMP)*D1(K,ICOMP) +
     &                     BODY(JCOMP)*D1(K,JCOMP))/BODY(ICM)
              FIRR = FIRR + FI(K,ICM)**2
              FDIRR = FDIRR + D1(K,ICM)**2
  150     CONTINUE
*
*       Obtain irregular time-step and check commensurability.
          DT = ETAI*SQRT(FIRR/FDIRR)
          CALL STEPK(DT,DTN)
          STEP(ICM) = DTN
          ITER = 0
  160     IF (DMOD(TIME,STEP(ICM)).NE.0.0D0) THEN
              STEP(ICM) = 0.5D0*STEP(ICM)
              ITER = ITER + 1
              IF (ITER.LT.16.OR.STEP(ICM).GT.DTK(40)) GO TO 160
              STEP(ICM) = DTK(40)
          END IF
          T0(ICM) = TIME
          T0R(ICM) = TIME
          TNEW(ICM) = TIME + STEP(ICM)
*
*       Form force components and first derivatives for c.m.
          DO 170 K = 1,3
              FR(K,ICM) = SAVEIT(K)
              D1R(K,ICM) = SAVEIT(K+3)
              FRDOT(K,ICM) = SAVEIT(K+3)
              F(K,ICM) = FI(K,ICM) + FR(K,ICM)
              FDOT(K,ICM) = D1(K,ICM) + D1R(K,ICM)
              F(K,ICM) = 0.5*F(K,ICM)
              FDOT(K,ICM) = ONE6*FDOT(K,ICM)
              D0(K,ICM) = FI(K,ICM)
              D2(K,ICM) = 0.0
              D3(K,ICM) = 0.0
              D2R(K,ICM) = 0.0
              D3R(K,ICM) = 0.0
  170     CONTINUE
*
*       Assign regular time-step and perform commensurability check.
          FIRR = 0.0
          FDIRR = 0.0
          DO 175 K = 1,3
              FIRR = FIRR + FR(K,ICM)**2
              FDIRR = FDIRR + D1R(K,ICM)**2
  175     CONTINUE
          DT = ETAR*SQRT(FIRR/FDIRR)
          CALL STEPK(DT,DTN)
          STEPR(ICM) = DTN
          ITER = 0
  180     IF (DMOD(TIME,STEPR(ICM)).NE.0.0D0) THEN
              STEPR(ICM) = 0.5D0*STEPR(ICM)
              ITER = ITER + 1
              IF (ITER.LT.16.OR.STEPR(ICM).GT.DTK(40)) GO TO 180
              STEPR(ICM) = DTK(40)
          END IF
      END IF
*
*       Skip KS initialization at merger termination (H, U & UDOT in RESET).
      IF (IPHASE.EQ.7) THEN
          EB = 2.0*EBH
          GO TO 50
      END IF
*
*       Define relative coordinates and velocities in physical scaled units.
      DO 20 K = 1,3
          Q(K) = X(K,ICOMP) - X(K,JCOMP)
          RDOT(K) = X0DOT(K,ICOMP) - X0DOT(K,JCOMP)
   20 CONTINUE
*
*       Introduce regularized variables using definition of Book.
      R(IPAIR) = SQRT(Q(1)**2 + Q(2)**2 + Q(3)**2)
*
*       Initialize the regularized coordinates according to sign of Q(1).
      IF (Q(1).LE.0.0D0) THEN
          UI(3) = 0.0D0
          UI(2) = SQRT(0.5D0*(R(IPAIR) - Q(1)))
          UI(1) = 0.5D0*Q(2)/UI(2)
          UI(4) = 0.5D0*Q(3)/UI(2)
      ELSE
          UI(4) = 0.0D0
          UI(1) = SQRT(0.5D0*(R(IPAIR) + Q(1)))
          UI(2) = 0.5D0*Q(2)/UI(1)
          UI(3) = 0.5D0*Q(3)/UI(1)
      END IF
*
*       Set current transformation matrix.
      CALL MATRIX(UI,A1)
*
*       Form regularized velocity and set initial KS coordinates.
      TDOT2(IPAIR) = 0.0
      DO 30 K = 1,4
          UDOT(K,IPAIR) = 0.50D0*(A1(1,K)*RDOT(1) + A1(2,K)*RDOT(2) +
     &                                              A1(3,K)*RDOT(3))
*       Note that A1(J,K) is the transpose of A1(K,J).
          U(K,IPAIR) = UI(K)
          U0(K,IPAIR) = U(K,IPAIR)
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0D0*UI(K)*UDOT(K,IPAIR)
   30 CONTINUE
*
*       Evaluate initial binding energy per unit mass (singular form).
      H(IPAIR) = (2.0D0*(UDOT(1,IPAIR)**2 + UDOT(2,IPAIR)**2 +
     &                   UDOT(3,IPAIR)**2 + UDOT(4,IPAIR)**2) -
     &                                              BODY(NTOT))/R(IPAIR)
      EB = H(IPAIR)*BODY(ICOMP)*BODY(JCOMP)/BODY(NTOT)
*
*       Form perturber list.
   50 CALL KSLIST(IPAIR)
*
*       Transform any unperturbed hard binary to apocentre and set time-step.
      SEMI = -0.5*BODY(NTOT)/H(IPAIR)
      IF (LIST(1,ICOMP).EQ.0.AND.EB.LT.EBH.AND.SEMI.LT.RMIN) THEN
          TK = TWOPI*SEMI*SQRT(SEMI/BODY(NTOT))
*       Note TIME is not commensurate after KSPERI (cf. CHTERM & STEPS).
          IF (IPHASE.NE.7.AND.IPHASE.NE.8) THEN
              DO 55 K = 1,4
                  VI(K) = UDOT(K,IPAIR)
   55         CONTINUE
*       Determine pericentre time (TP < 0 if TDOT2 < 0) and add TK/2.
              CALL TPERI(SEMI,UI,VI,BODY(NTOT),TP)
*       Note: apocentre to apocentre gives almost zero step.
              STEP(ICOMP) = 0.5*MIN(TK,STEP(NTOT)) - TP
*       Transform KS variables to peri and by pi/2 to apocentre (skip apo).
              IF (ABS(TDOT2(IPAIR)).GT.1.0E-12.OR.R(IPAIR).LT.SEMI) THEN
                  TIME0 = TIME
                  CALL KSPERI(IPAIR)
                  CALL KSAPO(IPAIR)
*       Reset TIME to quantized value (small > 0 or < 0 possible initially).
                  TIME = TIME0
                  TIME = MAX(TIME,0.0D0)
              ELSE IF (TDOT2(IPAIR).GT.0.0) THEN
                  TDOT2(IPAIR) = -1.0E-20
              END IF
          END IF
      END IF
*
*       Estimate an appropriate KS slow-down index for G < GMIN.
      IF (LIST(1,ICOMP).EQ.0.AND.SEMI.GT.0.0) THEN
          TK = TWOPI*SEMI*SQRT(SEMI/BODY(NTOT))
          IF (KZ(26).GT.0.AND.STEP(NTOT).GT.TK) THEN
              IMOD = 1 + LOG(STEP(NTOT)/TK)/0.69
              IMOD = MIN(IMOD,5)
          END IF
      END IF
*
*       Specify zero membership and set large steps for second component.
      LIST(1,JCOMP) = 0
*       Set large step for second component to avoid detection.
      STEP(JCOMP) = 1.0E+06
      STEPR(JCOMP) = 1.0E+06
*
*       Obtain polynomials for perturbed KS motion (standard case & merger).
      CALL KSPOLY(IPAIR,IMOD)
*
*       Obtain apocentre distance.
      SEMI = -0.5*BODY(NTOT)/H(IPAIR)
      ECC2 = (1.0-R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(NTOT)*SEMI)
      RAP = SEMI*(1.0 + SQRT(ECC2))
*
*       Include suggestion for monitoring hyperbolic encounters (suppressed).
*     IF (SEMI.LT.0.0) THEN
*         PMIN = SEMI*(1.0 - SQRT(ECC2))
*         WRITE (31,56)  TIME+TOFF, NAME(ICOMP), NAME(JCOMP), PMIN
*  56     FORMAT (' HYPERBOLIC    T NAM PMIN ',F8.2,2I6,1P,E10.2)
*     END IF
*
*       Set 2*SEMI as termination scale for hard binary if 2*SEMI < RS/20.
      IF (EB.LT.EBH.AND.RAP.LT.0.05*RS(NTOT)) THEN 
          R0(IPAIR) = MAX(RMIN,-BODY(NTOT)/H(IPAIR))
          R0(IPAIR) = MIN(R0(IPAIR),5.0*RMIN)
      ELSE
          R0(IPAIR) = R(IPAIR)
      END IF
*
*       Increase regularization counters (#9 & NKSHYP for hyperbolic orbits).
      NKSREG = NKSREG + 1
      IF (H(IPAIR).GT.0.0) NKSHYP = NKSHYP + 1
*
*       Check optional output for new KS.
      IF (KZ(10).GT.0) THEN
          RI = SQRT((X(1,NTOT) - RDENS(1))**2 +
     &              (X(2,NTOT) - RDENS(2))**2 +
     &              (X(3,NTOT) - RDENS(3))**2)
          WRITE (6,60)  TIME+TOFF, NAME(ICOMP), NAME(JCOMP),DTAU(IPAIR),
     &                  R(IPAIR), RI, H(IPAIR), IPAIR, GAMMA(IPAIR),
     &                  STEP(NTOT), LIST(1,ICOMP), LIST(1,NTOT)
   60     FORMAT (/,' NEW KSREG    TIME =',F8.2,2I6,F12.3,1P,E10.1,
     &                                0P,F7.2,F9.2,I5,F8.3,1P,E10.1,2I5)
      END IF
*
*       Include limited diagnostics for NS or BH hard binary formation.
      IF (MAX(KSTAR(ICOMP),KSTAR(JCOMP)).GE.13.AND.EB.LT.EBH.AND.
     &    IPHASE.NE.7) THEN
          ID = 0
          NP = IPIPE(1)
*       See whether at least one component was recorded in previous pair.
          DO 62 J = 2,NP+1
              IF (IPIPE(J).EQ.NAME(ICOMP).OR.
     &            IPIPE(J).EQ.NAME(JCOMP)) ID = ID + 1
   62     CONTINUE
*       Print diagnostics if NS/BH binary not listed among four last pairs.
          IF (ID.LE.1) THEN
              PD = DAYS*SEMI*SQRT(ABS(SEMI)/BODY(NTOT))
              WRITE (6,65)  TIME+TOFF, NAME(ICOMP), NAME(JCOMP),
     &                      KSTAR(ICOMP), KSTAR(JCOMP), KSTAR(NTOT),
     &                      SQRT(ECC2), PD, EB
   65         FORMAT (' NS/BH BINARY    T NM K* E P EB ',
     &                                  F8.1,2I6,3I4,F7.3,1P,E9.1,E11.2)
          END IF
*       Update table and membership by removing first two.
          IF (NP.GE.8) THEN
*       Note there are at most 4 pairs with entries in IPIPE(2:9).
              DO 66 J = 2,6
                  IPIPE(J) = IPIPE(J+2)
                  IPIPE(J+1) = IPIPE(J+3)
   66         CONTINUE
              NP = NP - 2
          END IF
*       Add NAME of each component in NP+2/3 and increase membership by 2.
          IPIPE(NP+2) = NAME(ICOMP)
          IPIPE(NP+3) = NAME(JCOMP)
          IPIPE(1) = NP + 2
      END IF
*
*       Modify the termination criterion according to value of NPAIRS.
      IF (NPAIRS.GT.KMAX - 3) GMAX = 0.8*GMAX
      IF (NPAIRS.LT.KMAX - 5.AND.GMAX.LT.0.001) GMAX = 1.2*GMAX
      IF (NPAIRS.EQ.KMAX) WRITE (6,70)  NPAIRS, TIME+TOFF
   70 FORMAT (5X,'WARNING!   MAXIMUM KS PAIRS   NPAIRS TIME',I5,F9.2)
*
*       Initialize prediction variables to avoid skipping KS components.
      DO 75 KCOMP = 1,2
          JDUM = 2*NPAIRS - 2 + KCOMP
          DO 72 K = 1,3
              X0(K,JDUM) = X(K,JDUM)
              X0DOT(K,JDUM) = 0.0D0
              F(K,JDUM) = 0.0D0
              FDOT(K,JDUM) = 0.0D0
              D2(K,JDUM) = 0.0D0
              D3(K,JDUM) = 0.0D0
              D2R(K,JDUM) = 0.0D0
              D3R(K,JDUM) = 0.0D0
   72     CONTINUE
   75 CONTINUE
*
*       See whether either component has been regularized recently.
      NNB = LISTD(1) + 1
      K = 0
*       Check case of initial binary and loop over disrupted pairs.
      IF (IABS(NAME(ICOMP) - NAME(JCOMP)).EQ.1) THEN
          IF (NAME(ICOMP).LE.2*NBIN0) K = -1
      END IF
      DO 80 L = 2,NNB
          IF (NAME(ICOMP).EQ.LISTD(L).OR.NAME(JCOMP).EQ.LISTD(L)) K = -1
   80 CONTINUE
      IF (H(IPAIR).GT.0.0) K = 0
*
*       Treat mergers as new binaries and ensure chain/hard KS as standard.
      IF (IPHASE.EQ.6) THEN
          K = 0
      ELSE IF (K.EQ.0) THEN
          IF (IPHASE.EQ.8.OR.EB.LT.EBH) K = -1
      END IF
*       Set flag to distinguish between existing and new binaries (K = -1/0).
      LIST(2,JCOMP) = K
*
*       Check diagnostic output of new hard binary.
      IF (KZ(8).GT.0.AND.K.EQ.0) THEN
          IF (EB.GT.EBH) GO TO 100
          SEMI = -0.5*BODY(NTOT)/H(IPAIR)
          RI = SQRT((X(1,NTOT) - RDENS(1))**2 +
     &              (X(2,NTOT) - RDENS(2))**2 +
     &              (X(3,NTOT) - RDENS(3))**2)
          WRITE (8,90)  TIME+TOFF, NAME(ICOMP), NAME(JCOMP), K,
     &                  BODY(ICOMP),BODY(JCOMP), EB, SEMI, R(IPAIR),
     &                  GAMMA(IPAIR), RI
   90     FORMAT (' NEW BINARY   T =',F7.1,'  NAME = ',2I6,I3,
     &                        '  M =',1P,2E9.1,'  EB =',E9.1,
     &                        '  A =',E9.1,'  R =',E9.1,'  G =',E9.1,
     &                        '  RI =',E9.1)
          CALL FLUSH(8)
      END IF
*
  100 RETURN
*
      END
