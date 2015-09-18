      SUBROUTINE KSINT(I1,KCASE)
*
*
*       Regularized integration.
*       ------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               ZP(NTMAX),ES(NTMAX),CZ(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      COMMON/SLOW0/  RANGE,ISLOW(10)
      COMMON/KSPAR/  ISTAT(KMAX)
      PARAMETER (ZZ=1.0/120.0D0)
      REAL*8  UI(4),UIDOT(4),XI(6),VI(6),FP(6),FD(6),RIDOT(3)
      REAL*8  EI0(3),EI(3)
      LOGICAL IQ
      SAVE ITERM
      SAVE  EI0
      DATA  EI0 /3*0/
      SAVE RI,HI
*
*
*       Set second component (= 2), pair index (= 1) & c.m. index.
      I2 = I1 + 1
      IPAIR = KVEC(I1)
      I = N + IPAIR
      JPHASE = 0
*
*       Define perturber membership & inverse c.m. mass.
      NNB0 = LIST(1,I1)
      BODYIN = 1.0/BODY(I)
*
*       Check for further unperturbed motion.
      IF (NNB0.EQ.0.AND.H(IPAIR).LT.0.0) THEN
          DT1 = STEP(I1)
          CALL UNPERT(IPAIR,KCASE)
          JPHASE = ISTAT(KCASE)
*
*       Update any unperturbed relativistic KS binary.
          IF (KZ(11).NE.0.AND.LIST(1,I1).EQ.0) THEN
              CALL BRAKE4(I1,I2,KCASE,DT1)
              IF (ISTAT(KCASE).LT.0) GO TO 100
          END IF
 
*      Try re-initialize chain WD/BH system after dormant KS (#11 only).
          IF (KZ(11).NE.0.AND.NCH.EQ.0) THEN   !  note absence of LIST(1,I1).
              IF (MIN(KSTAR(I1),KSTAR(I2)).GE.10) THEN
                  SEMI = -0.5*BODY(I)/H(IPAIR)
                  RIJ2 = 0.0
                  RD = 0.0
                  J = LIST(2,I1)
                  DO 220 K = 1,3
                      RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
                      RD = RD + (X(K,I)-X(K,J))*(XDOT(K,I)-XDOT(K,J))
  220             CONTINUE
                  RP = SQRT(RIJ2)
*       Skip outward motion or separation > 5*RMIN.
                  IF (RD.GT.0.0.OR.RP.GT.5.0*RMIN) GO TO 100
                  WRITE (6,222)  TIME+TOFF, JCLOSE, KSTAR(I1),
     &                           KSTAR(I2), LIST(1,I1), GAMMA(IPAIR),
     &                           SEMI, R(IPAIR)
  222             FORMAT (' ACTIVATE CHAIN    T JCL K* NP G A R ',
     &                                        F9.1,I7,3I4,1P,3E10.2)
                  KSPAIR = IPAIR
*       Restore unperturbed motion from BRAKE4 (NP = 0 fixes some problem).
                  IF (GAMMA(IPAIR).LT.1.0D-10) THEN
                      JCLOSE = 0
                      LIST(1,I1) = 0
                  END IF
                  KS2 = 0
*       Include case of binary as dominant perturber.
                  IF (JCLOSE.GT.N) THEN
                      KS2 = JCLOSE - N
                      IF (KS2.GT.KSPAIR) KS2 = KS2 - 1
                      JCOMP = JCLOSE
                      JP = JCLOSE - N
                      WRITE (6,223)  KSPAIR, KS2, JCLOSE, GAMMA(JP)
  223                 FORMAT (' BINARY PERT    KSP KS2 JCLOSE GJP ',
     &                                         2I4,I7,1P,E10.2)
                  ELSE
*       Avoid JCOMP > N & JCLOSE < N for spurious CALL KSTERM in DELAY.
                      JCOMP = JCLOSE
                  END IF
                  JPHASE = 8
                  CALL DELAY(JPHASE,-2)
*       Save KS parameters until end of block-step (JCMAX=0: no extra pert).
                  CALL DELAY(1,KS2)
                  JCMAX = 0
              END IF
          END IF
          GO TO 100
      END IF
*
*       Perform KS prediction of U & UDOT.
      DTU = DTAU(IPAIR)
      CALL KSPRED(IPAIR,I,BODYIN,DTU,UI,UIDOT,Q1,Q2,Q3,RIDOT)
*
*       Obtain the perturbing force & derivative.
      CALL KSPERT2(I1,I,NNB0,BODYIN,Q1,Q2,Q3,RIDOT,FP,FD)
*
*       Save old radial velocity & relative perturbation and set new GAMMA.
      RDOT = TDOT2(IPAIR)
      GI = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)*R(IPAIR)**2*BODYIN
      GAMMA(IPAIR) = GI
*
*       Apply the Hermite corrector.
      CALL KSCORR(IPAIR,UI,UIDOT,FP,FD,TD2,TDOT4,TDOT5,TDOT6)
      IF (IPHASE.LT.0) GO TO 90     ! KSPERT2 termination for PN.
*
*       Increase regularization time-step counter and update the time.
      NSTEPU = NSTEPU + 1
      T0(I1) = TIME
*
*       Define useful scalars.
      RI = R(IPAIR)
      HI = H(IPAIR)
*
*       Initialize termination indicator and check for large perturbation.
      IQ = .FALSE.
      IF (GI.GT.0.1) GO TO 2
      IF (GI.LT.0.03) THEN
          JCOMP = 0
          GO TO 20
      END IF
      CALL FLYBY(I,ITERM)
      IF (ITERM.EQ.0.AND.KSTAR(I).LT.0) THEN
          IQ = .TRUE.
      END IF
      IF (ITERM.EQ.1) THEN
*       Delay chain regularization search until near end of block-step.
          IF (TIME + STEP(I1).GT.TBLOCK) THEN
              CALL IMPACT(I,JPHASE)
              IF (IPHASE.GT.0) GO TO 100
          END IF
      ELSE IF (ITERM.EQ.2) THEN
          IQ = .TRUE.
          GO TO 20
      ELSE
          GO TO 20
      END IF
*
*       Find the dominant body for large perturbations.
    2 S = 4.0*STEP(I)
      FMAX = BODY(I)/RI**2
*       Initialize JCOMP for prediction and optional diagnostics in KSTERM.
      JCOMP = 0
      DO 10 L = 2,NNB0+1
          J = LIST(L,I1)
*       Only search bodies within twice the c.m. time-step.
          IF (STEP(J).GT.S) GO TO 10
*       Compare strong perturber and either component with current pair.
          DO 5 K = I1,I2
              RIJ2 = (X(1,J) - X(1,K))**2 + (X(2,J) - X(2,K))**2 +
     &                                      (X(3,J) - X(3,K))**2
              IF (BODY(J) + BODY(K).GT.RIJ2*FMAX) JCOMP = J
    5     CONTINUE
   10 CONTINUE
*
*       Try chain reg if close perturber <= N forms dominant pair.
      IF (JCOMP.GT.0.OR.GI.GT.0.1) THEN
*         IF (JCOMP.LE.N.OR.GI.GT.0.1) IQ = .TRUE.
          SEMI = -0.5*BODY(I)/HI
          IF (SEMI.LT.3.0*RMIN.AND.SEMI.GT.0.0) GO TO 84
      END IF
*
*       Check termination of hyperbolic encounter (R > R0 or R > 2*RMIN).
   20 IF (HI.GT.0.0D0.AND.NAME(I).GT.0) THEN
          IF ((RI.GT.R0(IPAIR).AND.GI.GT.GMAX).OR.RI.GT.2.0*RMIN.OR.
     &        (GI.GT.0.5.AND.TD2.GT.0.0)) THEN
*       Skip termination delay in case of velocity kick (cf. routine KSTERM).
              IF (HI.LT.100.0.OR.GI.GT.0.1.OR.RI.GT.5.0*RMIN) THEN
                  IQ = .TRUE.
              END IF
          END IF
      END IF
*
*       Choose basic regularized step using binding energy.
      W1 = 0.5/ABS(HI)
      W2 = R0(IPAIR)*BODYIN
      W1 = MIN(W1,W2)
      W2 = SQRT(W1)
*
*       Set new regularized step and convert to physical time units.
      DTU = 4.0*ETAU*W2
*
*       Include convergence criterion DH = H'*DTU + H''*DTU**2/2 = 0.001*|H|.
      IF (GI.GT.1.0D-05) THEN
          DH = 1.0E-03*MAX(ABS(HI),0.1D0)
          XF = 2.0*DH/ABS(HDOT2(IPAIR))
          YF = HDOT(IPAIR)/HDOT2(IPAIR)
          DTU2 = SQRT(XF + YF**2) - ABS(YF)
          DTU = MIN(DTU2,DTU)
      END IF
*
*       Convert to physical time units.
      IT = 0
   25 STEP(I1) = (((((ZZ*TDOT6*DTU + 0.2D0*TDOT5)*DTU + 0.5D0*TDOT4)*DTU
     &                     + TDOT3(IPAIR))*ONE6*DTU + TD2)*DTU + RI)*DTU
      IF (STEP(I1).LT.1.0D-12) THEN
          WRITE (6,30)  NAME(I1), KSLOW(IPAIR), HI, RI, DTU, STEP(I1),GI
   30     FORMAT (' NEGATIVE STEP    NM KSL H R DTU S1 G ',
     &                               I7,I4,1P,5E10.2)
      IF (IT.GE.3) THEN
      WRITE (6,31)  
   31 FORMAT (' KSINT! ')
      STOP
      END IF
          CALL FLUSH(6)
          DTU = 0.5*DTU
          IT = IT + 1
          IF (IT.LT.10) GO TO 25
      END IF
      DTAU(IPAIR) = DTU
*
*       See whether the KS slow-down procedure is activated.
      IMOD = KSLOW(IPAIR)
      IF (IMOD.GT.1) THEN
          ZMOD = FLOAT(ISLOW(IMOD))
          STEP(I1) = ZMOD*STEP(I1)
      END IF
*
*       Truncate STEP to quantized value.
      DT = STEP(I1)
      CALL STEPK(DT,DTN)
      STEP(I1) = DTN
*
*       Perform Newton-Raphson iteration.
      DTU = DTU*DTN/DT
      CALL NEWTON(TDOT4,TDOT5,TDOT6,IPAIR,DTU)
      DTAU(IPAIR) = DTU
*
*       Check diagnostics print option.
      IF (KZ(10).GE.3) THEN
          WRITE (6,40)  IPAIR, TIME, H(IPAIR), RI, DTAU(IPAIR), GI,
     &                  STEP(I1), IMOD, LIST(1,I1)
   40     FORMAT (3X,'KS MOTION',I6,2F10.4,1P,4E10.2,0P,2I4)
      END IF
*
*       Employ special termination criterion in merger case.
      IF (NAME(I).LT.0) THEN
*       Terminate if apocentre perturbation > 0.25 (R > SEMI) or GI > 0.25.
          IF (HI.LT.0.0) THEN
*             SEMI = -0.5*BODY(I)/HI
*             ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
*             A0 = SEMI*(1.0 + SQRT(ECC2))/RI
*       Replace eccentricity calculation with typical value.
*             A0 = 1.5*SEMI/RI
*             GA = GI*A0*A0*A0
*             IF (GA.GT.0.25.AND.RI.GT.SEMI) IQ = .TRUE.
*       Terminate wide orbits if perturber number is relatively large.
              IF (RI.GT.20*RMIN.AND.NNB0.GT.0.8*LIST(1,I)) IQ = .TRUE.
              IF (GI.GT.0.1.AND.RI.GT.RMIN) IQ = .TRUE.
              IF (GI.GT.0.01.AND.RI.GT.5.0*RMIN) IQ = .TRUE.
              IF (GI.GT.0.25) IQ = .TRUE.
*       Include extra condition for planet case.
              IF (MIN(BODY(I1),BODY(I2)).LT.0.05*BODYM) THEN
                  IF (GI.GT.2.0D-04) IQ = .TRUE.
              END IF
          ELSE
              IF (TD2.GT.0.0.AND.(GI.GT.GMAX.OR.RI.GT.RMIN)) IQ = .TRUE.
          END IF
          IF (.NOT.IQ) GO TO 60
      END IF
*
*       Delay termination until end of block for large perturbations.
      IF (IQ) THEN
          DTR = TBLOCK - TIME
*         WRITE (6,45)  IPAIR, TTOT, GI, RI, DTR, STEP(I1)
*  45     FORMAT (' TERM TEST    KS T G R DTR DT  ',
*    &                           I4,F10.4,F7.3,1P,E10.2,2E9.1)
          IF (DTR.LT.STEP(I1)) GO TO 90
      END IF
*
*       Check standard termination criterion (suppress on IQ = .true.).
      IF (RI.GT.R0(IPAIR).AND.RI.GT.2.0*RMIN.AND..NOT.IQ) THEN
*       Include termination for rare tidal capture starting at pericentre.
          IF (KSTAR(I).LT.0.AND.RI.GT.5.0*RMIN) GO TO 90
*       Impose a limit using size of neighbour sphere (100*R > 0.80*RS).
          IF (RI.GT.8.0D-03*RS(I).AND.GI.GT.1.0D-03) GO TO 90
*       See whether termination can be delayed for sufficient perturbers.
          IF (NNB0.LT.0.80*LIST(1,I).AND.GI.LT.0.1) GO TO 60
*       Check updating of R0 for newly hardened binary orbit.
          IF (HI.LT.-ECLOSE) THEN
              SEMI = -0.5*BODY(I)/HI
              R0(IPAIR) = MAX(RMIN,2.0D0*SEMI)
              R0(IPAIR) = MIN(R0(IPAIR),5.0*RMIN)
              GO TO 70
          END IF
*       Avoid repeated terminations for large distance (other limits above).
          IF (NNB0.GT.2.AND.LIST(1,I).GT.2) GO TO 90
      END IF
*
*       End integration cycle for hyperbolic motion.
   60 IF (HI.GE.0.0D0) THEN
          IF (RDOT*TD2.LT.0.0D0) THEN
*       Determine pericentre for hyperbolic two-body motion.
              SEMI = -0.5D0*BODY(I)/HI
              ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
              QPERI = SEMI*(1.0D0 - SQRT(ECC2))
*             DMIN2 = MIN(DMIN2,QPERI)
*
*       Check optional tidal interaction or stellar collision.
              IF (KZ(19).GE.3.AND.KSTAR(I).LT.10) THEN
                  VINF = SQRT(2.0*HI)*VSTAR
                  KS1 = KSTAR(I1)
                  KS2 = KSTAR(I2)
                  RX = MAX(RADIUS(I1),RADIUS(I2))
*       Determine maximum periastron factor for capture (VINF in km/sec).
                  IF (KZ(27).LE.2) THEN
                      RFAC = RPMAX2(RADIUS(I1),RADIUS(I2),BODY(I1),
     &                              BODY(I2),KS1,KS2,VINF)
                      RCAP = RFAC*RX
                  ELSE
                      DV = SQRT(2.0*HI)
*       Note that Quinlan & Shapiro function returns actual distance.
                      RCAP = RPMAX(BODY(I1),BODY(I2),VSTAR,DV,QPERI)
                  END IF
                  IF (QPERI.LT.5.0*RX) THEN
                      WRITE (54,54)  TTOT, NAME(I1), NAME(I2), KS1,
     &                               KS2, VINF, RCAP*SU, RX*SU, QPERI*SU
   54                 FORMAT (' CLOSE   T NAM K* VINF RCAP RX QP  ',
     &                                  F7.1,2I6,2I4,F6.2,3F6.1)
                  END IF
                  IF (QPERI.LT.RCAP) THEN
                      J1 = I1
                      IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                      FAC = 0.5*BODY(I)/BODY(J1)
*       Set possible BH index and check disruption condition (& #43) first.
                      J2 = 2*IPAIR + 1 - J1
                      IF (KZ(43).GE.2.AND.KSTAR(J2).EQ.14) THEN
                          RCOLL = (BODY(J2)/BODY(J1))**0.3333*RADIUS(J1)
                      ELSE IF (KZ(27).LE.2) THEN
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                          RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                      ELSE
                          RCOLL = 6.0*BODY(I)/CLIGHT**2
                      END IF
                      WRITE (55,58)  TTOT, IPAIR, NAME(I1), NAME(I2),
     &                               KS1, KS2, KSTAR(I), VINF
   58                 FORMAT (' RPMAX:    T KS NAM K* VINF ',
     &                                    F7.1,I5,2I6,3I4,F7.2)
                      WRITE (55,59)  SQRT(ECC2), HI, R(IPAIR), SEMI,
     &                               QPERI, BODY(I1), BODY(I2),
     &                               BODY(I)*ZMBAR
   59                 FORMAT (' E H R A QP BODY MT ',
     &                          F9.5,1P,6E10.2,0P,F6.1)
                      RI2 = 0.0
                      VI2 = 0.0
                      DO 61 K = 1,3
                          RI2 = RI2 + (X(K,I) - CMR(K))**2
                          VI2 = VI2 + XDOT(K,I)**2
   61                 CONTINUE
                      WRITE (55,62)  SQRT(RI2)/RC, SQRT(VI2)*VSTAR,
     &                               RHOD, RADIUS(I1)*SU, RADIUS(I2)*SU,
     &                               RCAP, RADIUS(J1)/QPERI, RCOLL/QPERI
   62                 FORMAT (' r/RC V* <C> R* RCAP R1/QP RCOLL/QP ',
     &                          2F5.1,3F6.1,3F5.1)
                      CALL FLUSH(55)
                      IF (QPERI.LT.RCOLL) THEN
*       Obtain KS variables at pericentre before merging into one body.
                          CALL KSPERI(IPAIR)
                          KSPAIR = IPAIR
                          JPHASE = -1
                          IQCOLL = -2
*                         CALL CMBODY(QPERI,2)
                      ELSE IF (KSTAR(I).GE.0.AND.KZ(27).GT.0) THEN
                          CALL KSTIDE(IPAIR,KCASE,QPERI)
                      END IF
                  END IF
*       Check options for artificial collisions.
              ELSE IF (KZ(27).EQ.-1.AND.KZ(13).LT.0) THEN
                  RFAC = 2.0
                  IF (QPERI.LT.RFAC*MAX(RADIUS(I1),RADIUS(I2))) THEN
                      J1 = I1
                      IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                      FAC = 0.5*BODY(I)/BODY(J1)
                      IF (KZ(27).LE.2) THEN
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                          RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                      ELSE
                          RCOLL = 6.0*BODY(I)/CLIGHT**2
                      END IF
                      IF (QPERI.LT.RCOLL) THEN
                          CALL TOUCH(IPAIR,I1,I2,RCOLL)
                      END IF
                  END IF
              END IF
          END IF
          GO TO 100
      END IF
*
*       Check optional output of periapse angle on unit #12.
      IF (KZ(29).GE.0) GO TO 70
*
*       Evaluate Runge-Lenz eccentricity vector from physical variables.
      CALL KSRES(IPAIR,I1,I2,0.0D0)
      VI2 = 0.0
      RD = 0.0
      RR = 0.0
      DO 63 K = 1,3
          XI(K) = X(K,ICOMP) - X(K,JCOMP)
          VI(K) = XDOT(K,ICOMP) - XDOT(K,JCOMP)
          RR = RR + XI(K)**2
          VI2 = VI2 + VI(K)**2
          RD = RD + XI(K)*VI(K)
   63 CONTINUE
      RR = SQRT(RR)
      EI2 = 0.0
      DO 64 K = 1,3
          EI(K) = (VI2*XI(K) - RD*VI(K))/BODY(NTOT) - XI(K)/RR
          EI2 = EI2 + EI(K)**2
   64 CONTINUE
*       Save initial values for error plot.
      IF (EI0(1).EQ.0.0D0) THEN
          DO 65 K = 1,3
              EI0(K) = EI(K)
   65     CONTINUE
      ELSE
*       Construct cross product to yield deviation (in degrees).
          ERR = (EI(1)*EI0(2) - EI(2)*EI0(1))/EI2
          ERR = 360.0*ASIN(ERR)/TWOPI
          WRITE (49,66)  TIME/TWOPI, ERR
   66     FORMAT (' ',F12.6,1P,E12.4)
      END IF
*
*       Determine new perturbers for binary at apocentre turning point.
   70 IF (RDOT*TD2.GE.0.0D0) GO TO 100
      SEMI = -0.5D0*BODY(I)/HI
*
*       Check minimum two-body separation just after pericentre.
      IF (RDOT.LT.0.0D0) THEN
*       Obtain pericentre by Mikkola's algorithm (GAMMA < 0.001).
          IF (GI.LT.0.001) THEN
              CALL PERI(UI,UIDOT,RI,BODY(I1),BODY(I2),QPERI)
          ELSE
              QPERI = RI
          END IF
*         DMIN2 = MIN(DMIN2,QPERI)
*
*       Check optional tidal interaction or stellar collision (skip merger).
          IF (KZ(19).GE.3.AND.KSTAR(I).LE.10.AND.NAME(I).GT.0) THEN
              RFAC = 5.0
              IF (KZ(27).LE.2) THEN
                  IF (KZ(27).EQ.1) RFAC = 4.0
                  RX = RFAC*MAX(RADIUS(I1),RADIUS(I2))
              ELSE
                  RX = RPMIN(BODY(I1),BODY(I2),VSTAR,HI,QPERI)
              END IF
              IF (QPERI.LT.RX) THEN
                  J1 = I1
                  IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                  FAC = 0.5*BODY(I)/BODY(J1)
*       Set possible BH index and check disruption condition (& #43) first.
                  J2 = 2*IPAIR + 1 - J1
                  IF (KZ(43).GE.2.AND.KSTAR(J2).EQ.14) THEN
                      RCOLL = (BODY(J2)/BODY(J1))**0.3333*RADIUS(J1)
                  ELSE IF (KZ(27).LE.2) THEN
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                      RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                  ELSE
                      RCOLL = 6.0*BODY(I)/CLIGHT**2
                  END IF
                  IF (QPERI.LT.RCOLL) THEN
*       Obtain KS variables at pericentre before merging into one body.
                      CALL KSPERI(IPAIR)
                      KSPAIR = IPAIR
                      JPHASE = -1
                      IQCOLL = -2
*                     CALL CMBODY(QPERI,2)
                  ELSE IF (KSTAR(I).GE.0) THEN
*       Distinguish between sequential, standard and GR circularization.
                      IF (KZ(27).EQ.1) THEN
                          ICIRC = 1
                          TC = 0.0
                      ELSE IF (KZ(27).EQ.2.AND.KSTAR(I).LT.10) THEN
                          ECC2 = (1.0 - RI/SEMI)**2 +
     &                                    TDOT2(IPAIR)**2/(BODY(I)*SEMI)
                          ECC = SQRT(ECC2)
                          ICIRC = 0
                          CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC)
                      ELSE
                          ICIRC = 1
                          TC = 0.0
                      END IF
                      IF (KSTAR(I).GE.10) ICIRC = 0
*       Skip tidal effects for circularization time above 100 Myr (07/08).
                      IF (ICIRC.GT.0.AND.KZ(27).GT.0.AND.
     &                    TC.LT.100.0) THEN
                          CALL KSTIDE(IPAIR,KCASE,QPERI)
                      END IF
                  END IF
              END IF
*       Check for perturbed spiral or chaos case (skip collision).
              IF (KSTAR(I).EQ.-2.AND.IPHASE.EQ.0) THEN
                  CALL SPIRAL(IPAIR)
              ELSE IF (KSTAR(I).EQ.-1.AND.IPHASE.EQ.0) THEN
                  CALL KSTIDE(IPAIR,KCASE,QPERI)
              END IF
*       Check options for artificial collisions.
          ELSE IF (KZ(27).EQ.-1.AND.KZ(13).LT.0) THEN
              RFAC = 2.0
              IF (QPERI.LT.RFAC*MAX(RADIUS(I1),RADIUS(I2))) THEN
                  J1 = I1
                  IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                  FAC = 0.5*BODY(I)/BODY(J1)
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                  RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                  IF (QPERI.LT.RCOLL) THEN
                      CALL TOUCH(IPAIR,I1,I2,RCOLL)
                  END IF
              END IF
          END IF
          GO TO 100
      END IF
*
*       Save maximum separation of persistent binary (skip parallel case).
*     RMAX = MAX(RMAX,RI)
*
*       Check binary reference radius or merger stability criterion.
      IF (NAME(I).GT.0) THEN
*       Update termination length scale in case of initial soft binary.
          EB = BODY(I1)*BODY(I2)*HI*BODYIN
          IF (EB.LT.EBH) R0(IPAIR) = MAX(RMIN,2.0*SEMI)
      ELSE
          ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          ECC = SQRT(ECC2)
          RP = SEMI*(1.0 - ECC)*(1.0 - 2.0*GI)
*       Find merger index.
          IM = 0
          DO 72 K = 1,NMERGE
              IF (NAMEM(K).EQ.NAME(I)) IM = K
   72     CONTINUE
*       Exclude inner planets from the general stability test.
          IF (MIN(BODY(I1),BODY(I2)).LT.0.05*BODYM) THEN
              IF (RP.LT.R0(IPAIR)) GO TO 90
          END IF
*       Include optional Kozai diagnostics (large EMAX) & circularization.
          IF (KZ(42).GT.0) THEN
              CALL KOZAI(IPAIR,IM,ECC,SEMI,ITERM)
*       Check merger termination and perform circularization or collision.
              IF (ITERM.GT.0) THEN
                  IPHASE = 7
                  KSPAIR = IPAIR
                  CALL RESET
                  IF (ITERM.GT.1) CALL CIRC(KCASE)
                  GO TO 100
              END IF
          END IF
*       Assess the stability inside critical pericentre (safety factor 1.04).
          IF (RP.LT.1.04*R0(IPAIR)) THEN
*       Note: assessment needs to use same eccentricity as for acceptance.
              CALL ASSESS(IPAIR,IM,ECC,SEMI,ITERM)
              IF (ITERM.GT.0) THEN
                  INSTAB = INSTAB + 1
                  GO TO 90
              END IF
          END IF
*       Check possible eccentricity modulation or t_circ update.
          IF (IM.GT.0.AND.(TIME.GT.TMDIS(IM).OR.
     &        TMDIS(IM).GT.1.0D+06)) THEN
              IF (KZ(27).EQ.2) THEN
                  CALL ECCMOD(I,ITERM)
                  IF (ITERM.GT.0) THEN
*                     WRITE (6,76)  RP, R0(IPAIR)
*  76                 FORMAT (' ECCMOD TERM    RP R0 ',1P,2E10.2)
                      GO TO 90
                  END IF
*       Consider both inner and possible outer circularizing binary.
                  DO 78 K = 1,NCHAOS
                      IF (NAMEC(K).EQ.NZERO - NAMEM(IM).AND.
     &                    KSTARM(IM).EQ.-2) THEN
*       Update unperturbed binary if T - TOSC > 10 Myr (cf. IMPACT & DECIDE).
                          IF ((TIME - TOSC(K))*TSTAR.GT.10.0) GO TO 90
                      END IF
                      IF (NAMEC(K).EQ.NAMEG(IM).AND.
     &                    KSTARM(IM).EQ.-2) THEN
                          IF ((TIME - TOSC(K))*TSTAR.GT.10.0) GO TO 90
                      END IF
*       Note: perturbed binary is treated if pericentre before next IMPACT.
   78             CONTINUE
              END IF
          END IF
      END IF
*
*       See whether KS slow-down procedure should be (re)-checked (no Chaos).
      IF (KZ(26).GT.0.AND.KSTAR(I).GE.0) THEN
          KMOD = RANGE*GMIN/MAX(GI,1.0D-10)
          IF (KMOD.GT.1.OR.IMOD.GT.1) THEN
              CALL KSMOD(IPAIR,KMOD)
              IF (KMOD.LT.0) GO TO 100
              GO TO 80
          END IF
      END IF
*
*       Set approximate value of next period.
      TK = TWOPI*SEMI*SQRT(SEMI*BODYIN)
      IF (IMOD.GT.1) THEN
          TK = ZMOD*TK
      END IF
*
*       Use old perturber list if next apocentre is before the c.m. step.
      IF (TIME + TK.LT.T0(I) + STEP(I)) THEN
          GO TO 100
      END IF
*
*       Select new perturbers (adopt unperturbed period if none found).
   80 CALL KSLIST(IPAIR)
*
*       Check rectification of chaotic spiral at start of unperturbed motion.
      IF (KSTAR(I).EQ.-2.AND.LIST(1,I1).EQ.0) THEN
          DMR = 0.D0
          CALL CHRECT(IPAIR,DMR)
          IF (IPHASE.LT.0) JPHASE = -1
          IF (IPHASE.LT.0) GO TO 100
      ELSE
          CALL KSRECT(IPAIR)
      END IF
*
*       See whether a massive BH subsystem can be selected.
*     IF (KZ(11).NE.0.AND.NCH.EQ.0.AND.BODY(I).GT.10.0*BODYM.AND.
  84  IF (NCH.EQ.0.AND.SEMI.LT.3.0*RMIN.AND.NAME(I).GT.0) THEN
*
*       Check optional BH condition (prevents mass-loss complications).
          IF (KZ(11).LE.-2) THEN
              IF (KSTAR(I1).NE.14.OR.KSTAR(I2).NE.14) GO TO 88
          END IF
*
          RCR2 = RMIN2*BODY(I)/BODYM
          RCR2 = 4.0*MAX(4.0*RMIN2,RCR2)
          NNB1 = LIST(1,I1) + 1
          RX2 = 1000.0
          JCLOSE = 0
          DO 85 L = 2,NNB1
              J = LIST(L,I1)
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
              RD = (X(1,I)-X(1,J))*(XDOT(1,I)-XDOT(1,J)) +
     &             (X(2,I)-X(2,J))*(XDOT(2,I)-XDOT(2,J)) +
     &             (X(3,I)-X(3,J))*(XDOT(3,I)-XDOT(3,J))
              IF (RIJ2.LT.RCR2.AND.RD.LT.0.0) THEN
                  IF (RIJ2.LT.RX2) THEN
                      RX2 = RIJ2
                      RRD = RD
                      JCLOSE = J
                      VIJ2 = (XDOT(1,I) - XDOT(1,J))**2 +
     &                       (XDOT(2,I) - XDOT(2,J))**2 +
     &                       (XDOT(3,I) - XDOT(3,J))**2
                  END IF
              END IF
   85     CONTINUE
          IF (JCLOSE.GT.0) THEN
              IF (NAME(JCLOSE).LE.0.OR.NCH.GT.0) GO TO 88
*       Skip chain selection if close perturber > 5*SEMI or impact > 3*SEMI.
              RX = SQRT(RX2)
              SEMI1 = 2.0/RX - VIJ2/(BODY(I) + BODY(JCLOSE))
              SEMI1 = 1.0/SEMI1
              ECC2 = (1.0-RX/SEMI1)**2 + RRD**2/(BODY(I)+BODY(JCLOSE))
              ECC1 = SQRT(ECC2)
              PMIN = SEMI1*(1.0 - ECC1)
              IF (RX.GT.5.0*SEMI.OR.PMIN.GT.3.0*SEMI) GO TO 88
*       Include distance criterion for B-B (avoids possible chain problem).
              IF (JCLOSE.GT.N) THEN
                  JP = JCLOSE - N
                  AJ = -0.5*BODY(JCLOSE)/H(JP) ! needs proper stability test.
                  IF (RX.GT.2.0*(SEMI + AJ)) GO TO 88
              END IF
*       Limit energy of triple system (< 5*EBH) using radial velocity.
*             ZMU = BODY(I)*BODY(JCLOSE)/(BODY(I) + BODY(JCLOSE))
*             EBT = EB + ZMU*(0.5*(RRD/RX)**2-(BODY(I)+BODY(JCLOSE)/RX))
*             IF (EBT.GT.5.0*EBH) GO TO 88
              WRITE (6,86)  TTOT, NAME(JCLOSE), LIST(1,I1), STEP(I),
     &                      STEP(JCLOSE), SEMI, RX, GI
   86         FORMAT (' NEW CHAIN   T NMJ NP STEPI STEPJ A RIJ GI ',
     &                              F9.3,I6,I4,1P,5E10.2)
*       Initiate chain regularization directly (B-B or B-S: see IMPACT).
              JCOMP = JCLOSE
              JCMAX = 0
              KSPAIR = IPAIR
              IPHASE = 8
              EBCH0 = EB
*       Distinguish between case of single and binary intruder.
              IF (JCLOSE.LE.N) THEN
                  KS2 = 0
              ELSE
                  KS2 = JCLOSE - N
*       Reduce c.m. body location and pair index if IPAIR terminated first.
                  IF (KS2.GT.IPAIR) JCLOSE = JCLOSE - 1
                  IF (KS2.GT.IPAIR) KS2 = KS2 - 1
              END IF
*       Specify index for ARchain initialization (cf. SUBINT).
              ISTAT(KCASE) = 8
*       Initialize new ARchain.
              CALL DELAY(1,KS2)
              GO TO 100
          END IF
      END IF
*
*       Include KS termination after failed chain test GO TO 84.
   88 IF (IQ) GO TO 90
*
*       Check optional search criterion for multiple encounter or merger.
      IF (KZ(15).GT.0.AND.STEP(I).LT.DTMIN) THEN
          CALL IMPACT(I,JPHASE)
      END IF
      GO TO 100
*
*       Terminate regularization of current pair (IPAIR set in KSPAIR).
   90 KSPAIR = IPAIR
*       Set indicator for calling KSTERM in MAIN (permits phase overlay).
      JPHASE = 2
*       Check case of hierarchical binary.
      IF (NAME(I).LT.0) JPHASE = 7
*
*       Save activity index unless already non-zero.
  100 IF (ISTAT(KCASE).EQ.0) ISTAT(KCASE) = JPHASE
*
      RETURN
*
      END
