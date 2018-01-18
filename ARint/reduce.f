      SUBROUTINE REDUCE(IESC,JESC,ISUB)
*
*
*       Reduction of chain.
*       -------------------
*
      INCLUDE 'common6.h'
        REAL*8  M,MASS,MC,MMIJ
        PARAMETER (NMX=10,NMX3=3*NMX,NMX4=4*NMX,NMXm=NMX*(NMX-1)/2)
         COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &   XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/ARCHAIN2/ MMIJ,CMX(3),CMV(3),ENERGY,EnerGR,CHTIME
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/ARZERO/  ISTAR0(NMX),SIZE0(NMX)
      COMMON/POSTN/  CVEL,TAUGR,RZ1,GAMMAZ,TKOZ,EMAX,TSP,KZ24,IGR,IPN
      COMMON/POSTN2/ SEMIGR,ECCGR,DEGR,ISPIN
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/ECHAIN/  ECH
      COMMON/KSAVE/ K1,K2
      DATA TCALL /0.0D0/
      REAL*8  XCM(3),VCM(3),DXC(3),DVC(3),FIRR(3),FD(3),XREL(3),
     &        VREL(3),XX(3,3),VV(3,3),X0S(3),V0S(3),XS(3),VS(3),G(3)
*     REAL*8  CG(6)
      SAVE   ! this was a bug fix 14/5/11 - JESC lost.
*
*
*       Include possible testing of energy budget (checked OK).
*     ID = 0
*     IF (JESC.GT.0) ID = 1
*
*       Copy chain variables to standard form (CHMOD may be skipped).
      LK = 0
      DO 4 L = 1,NCH
          DO 3 K = 1,3
              LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
    3     CONTINUE
    4 CONTINUE
*
*       Determine dominant chain force for elements of innermost binary.
      FX = 0.0
      DO 5 K = 1,NCH-1
         K1 = INAME(K)
         K2 = INAME(K+1)
         FF = (M(K1) + M(K2))*RINV(K)**2
         IF (FF.GT.FX) THEN
             I1 = K1
             I2 = K2
             FX = FF
         END IF
    5 CONTINUE
*
*       Employ dominant components instead of IESC & JESC (JESC may be = 0).
*     I1 = IESC
*     I2 = JESC
      ISING = 0
      SEMIX = 1.0
      KK = I1
      LX = I2    ! use LX for reference to possible binary component.
      RIJ2 = 0.0
      VIJ2 = 0.0
      RDOT = 0.0
*       Evaluate binary diagnostics before initializing escaper.
      DO 11 K = 1,3
          RIJ2 = RIJ2 + (X4(K,KK) - X4(K,LX))**2
          VIJ2 = VIJ2 + (XDOT4(K,KK) - XDOT4(K,LX))**2
          RDOT = RDOT + (X4(K,KK) - X4(K,LX))*
     &                  (XDOT4(K,KK) - XDOT4(K,LX))
          XREL(K) = X4(K,KK) - X4(K,LX)
          VREL(K) = XDOT4(K,KK) - XDOT4(K,LX)
          XX(K,1) = X4(K,KK)
          XX(K,2) = X4(K,LX)
          VV(K,1) = XDOT4(K,KK)
          VV(K,2) = XDOT4(K,LX)
          XCM(K) = (M(KK)*X4(K,KK) + M(LX)*X4(K,LX))/(M(KK) + M(LX))
          VCM(K) = (M(KK)*XDOT4(K,KK)+M(LX)*XDOT4(K,LX))/(M(KK)+M(LX))
   11 CONTINUE
      RIJ = SQRT(RIJ2)
      A2 = 2.0/RIJ - VIJ2/(M(KK) + M(LX))
      SEMIX = 1.0/A2
      ECC2 = (1.0 - RIJ/SEMIX)**2 + RDOT**2/(SEMIX*(M(KK) + M(LX)))
      ECC = SQRT(ECC2)
      EB = -0.5*M(I1)*M(I2)/SEMIX
*
*       Obtain the relativistic elements if relevant.
      IF (SEMIX.GT.0.0.AND.CVEL.GT.0.0D0) THEN
          CALL GRBIN(M(I1),M(I2),XREL,VREL,SEMI,ECC)
          SEMIX = SEMI
          EB = -0.5*M(I1)*M(I2)/SEMI
      END IF
      IF (SEMIX.GT.0.0) THEN
          SEMIGR = SEMIX
          ECCGR = ECC
      END IF
*
*       Determine smallest semi-major axis around the c.m. (exclude I1-I2).
      SEMI2 = 1.0
      ECC2 = 1.0
      ZI = 0.0
      KX = 1
      DO 13 KK = 1,NCH
          IF (KK.EQ.I1.OR.KK.EQ.I2) GO TO 13
          RIJ2 = 0.0
          VIJ2 = 0.0
          RDOT = 0.0
          DO 12 K = 1,3
              RIJ2 = RIJ2 + X4(K,KK)**2
              VIJ2 = VIJ2 + XDOT4(K,KK)**2
              RDOT = RDOT + X4(K,KK)*XDOT4(K,KK)
   12     CONTINUE
          RIJ = SQRT(RIJ2)
          A2 = 2.0/RIJ - VIJ2/MASS
          A2 = 1.0/A2
          IF (A2.GT.0.0.AND.A2.LT.SEMI2) THEN
              SEMI2 = A2
              ECC2 = (1.0 - RIJ/A2)**2 + RDOT**2/(A2*MASS)
              ECC2 = SQRT(ECC2)
              KX = KK
          END IF
   13 CONTINUE
      DO 7 K = 1,3
          XX(K,3) = X4(K,KX)
          VV(K,3) = XDOT4(K,KX)
    7 CONTINUE
      IF (ECC2.LT.1.0) THEN
          CALL INCLIN(XX,VV,XCM,VCM,ALPH)
          ZI = 360.0*ALPH/TWOPI
      ELSE
          ZI = 0.0D0
      END IF
      IF (SEMIX.LT.0.0) THEN
          SEMIX = SEMI2
          ECC = 0.0
          ZI = 0.0
      END IF
      K1 = I1
      IF (M(I1).GT.M(I2)) K1 = I2
*       Swap values in rare case of SEMI2 being smaller (wrong ID).
      IF (SEMIX.GT.SEMI2) THEN
          K1 = KX
          A2 = SEMIX
          SEMIX = SEMI2
          SEMI2 = A2
          ECC = ECC2
          EB = -0.5*M(K1)*MASS/SEMI2
          ZI = 0.0
      END IF
*       Skip large semi-major axis (dominant body in eccentric orbit).
      IF (ECC.GT.0.0.AND.SEMIX.LT.1.0D-03) THEN
          RI = SQRT((X(1,ICH)-RDENS(1))**2 + (X(2,ICH)-RDENS(2))**2 +
     &                                       (X(3,ICH)-RDENS(3))**2)
          WRITE (20,14)  TIME+TOFF, NAMEC(K1), ECC, SEMIX, EB,
     &                   NCH, NPERT, IPN, ZI, GPERT, TKOZ, RI
   14     FORMAT (' ',F9.3,I6,F8.4,1P,E12.4,0P,F9.4,3I4,F7.1,1P,2E9.1,
     &                0P,F7.3)
          CALL FLUSH(20)
      END IF
*
      IF (MOD(NCHAIN,10).EQ.0) THEN
          NM1 = LISTC(1) + 1
          JM = 0
          BM = 50.0*BODYM
          DO 105 L = 2,NM1
              J = LISTC(L)
              IF (BODY(J).GT.BM) THEN
                  BM = BODY(J)
                  JM = J
              END IF
  105     CONTINUE
          IF (JM.GT.0) THEN
              RIJ2 = 0.0
              VIJ2 = 0.0
              RDOT = 0.0
              DO 106 K = 1,3
                  RIJ2 = RIJ2 + (X(K,JM) - X(K,ICH))**2
                  VIJ2 = VIJ2 + (XDOT(K,JM) - XDOT(K,ICH))**2
                  RDOT = RDOT +
     &                  (X(K,JM) - X(K,ICH))*(XDOT(K,JM) - XDOT(K,ICH))
  106         CONTINUE
              RIJ = SQRT(RIJ2)
              SEMIJ = 2.0/RIJ - VIJ2/(BODY(JM) + BODY(ICH))
              SEMIJ = 1.0/SEMIJ
              ECC2 = (1.0 - RIJ/SEMIJ)**2 +
     &                RDOT**2/((BODY(JM) + BODY(ICH))*SEMIJ)
              ECC = SQRT(ECC2)
              IF (ZI.GT.0.0) THEN
                  DO 107 K = 1,3
                      XX(K,3) = X(K,JM) - X(K,ICH)
                      VV(K,3) = XDOT(K,JM) - XDOT(K,ICH)
  107             CONTINUE
                  CALL INCLIN(XX,VV,XCM,VCM,ALPH)
                  ZI = 360.0*ALPH/TWOPI
              END IF
              WRITE (40,108)  NAME(JM), TIME, ZI, BODY(JM), RIJ,
     &                        ECC, SEMIJ
  108         FORMAT (' HEAVY   NAM T INC M RIJ E A ',
     &                          I7,2F7.1,1P,2E10.2,0P,F7.3,1P,E10.2)
              CALL FLUSH(40)
          END IF
      END IF
*
*       Obtain statistics for any single BHs at regular intervals.
      IF (IABS(KZ(11)).GT.1.AND.TCALL.LT.TTOT) THEN
          TCALL = TTOT + 1.0
          CALL BHSTAT
      END IF
*
*       Define discrete time for polynomials (new attempt 12/2015).
      TIME0 = TIME
      TIME = TBLOCK
*
*       Check whether to treat two singles instead of KS (case of JESC > 0).
      IF (JESC.GT.0) THEN
          ZMB = BODYC(IESC) + BODYC(JESC)
*       Apply perturbation limit of 10 % for new KS (otherwise two singles).
          RCR = (0.1*ZMB/MASS)**0.3333*RSUM
          IB = 1
          IF (IESC.NE.INAME(1).OR.1.0/RINV(NCH-1).LT.RCR) IB = NCH - 1
          IF (1.0/RINV(1).LT.1.0/RINV(NCH-1)) IB = 1   ! Ensure the smallest.
          RB = 1.0/RINV(IB)
          IF (RB.GT.RMIN.OR.RB.GT.RCR) THEN
              ISING = 1
*       Switch to KS for small pericentre.
              IF (SEMIX*(1.0 - ECC).LT.0.1*RMIN) ISING = 0
          END IF
      ELSE
          ISING = 0
      END IF
*
*       Save global chain variables for correction procedure.
      ZM0S = BODY(ICH)
      DO 120 K = 1,3
          X0S(K) = X(K,ICH)
          V0S(K) = XDOT(K,ICH)
  120 CONTINUE
*
*       Make new chain from remaining members.
    2 LK = 0
      DPOT = 0.0
      DO 10 L = 1,NCH
          IF (L.EQ.IESC) GO TO 10
          RIJ2 = 0.0
          DO 8 K = 1,3
              LK = LK + 1
              XCH(LK) = X4(K,L)
              VCH(LK) = XDOT4(K,L)
              RIJ2 = RIJ2 + (XCH(LK) - X4(K,IESC))**2
    8     CONTINUE
          DPOT = DPOT + BODYC(L)*BODYC(IESC)/SQRT(RIJ2)
   10 CONTINUE
*
*       Reduce chain membership and mass (global & local COMMON).
      NCH = NCH - 1
      NN = NCH
      MASS = MASS - M(IESC)
*       Improve coordinates & velocities of c.m. body to order F3DOT.
      IF (ABS(TIME-T0(ICH)).LE.STEP(ICH)) THEN
          CALL XVPRED(ICH,-1)
      END IF
*
*       Set new c.m. for reduced system and save old c.m. variables.
      DO 20 K = 1,3
          DXC(K) = -BODYC(IESC)*X4(K,IESC)/(BODY(ICH) - BODYC(IESC))
          DVC(K) = -BODYC(IESC)*XDOT4(K,IESC)/(BODY(ICH) - BODYC(IESC))
          XCM(K) = X(K,ICH) + DXC(K)
          VCM(K) = XDOT(K,ICH) + DVC(K)
          CM(K) = X(K,ICH)
          CM(K+3) = XDOT(K,ICH)
   20 CONTINUE
*
*       Re-define new chain variables w.r. to modified c.m. (NB! retain X4).
      LK = 0
      DO 30 L = 1,NCH
*       Note opposite sign of DXC & DVC above (actually dx = m_p*x_p/(M-m_p).
          DO 25 K = 1,3
              LK = LK + 1
              XCH(LK) = XCH(LK) - DXC(K)
              VCH(LK) = VCH(LK) - DVC(K)
   25     CONTINUE
   30 CONTINUE
*
*       Save original mass and neighbour radius of c.m. body.
      BODYCH = BODY(ICH)
      RS0 = RS(ICH)
*
*       Search for global index of escaper.
      I = 0
      DO 40 J = IFIRST,NTOT
          IF (NAME(J).EQ.NAMEC(IESC)) THEN
              I = J
              IF (BODY(J).GT.0.0D0) WRITE (6,35)  I, IESC, NAMEC(IESC)
   35         FORMAT (' WARNING!   NON-ZERO GHOST    I IESC NAMEC ',3I5)
              GO TO 55
          END IF
   40 CONTINUE
*
*       Switch to another reference body since #ICH is escaping (NAME = 0).
      I = ICH
      IF (IESC.GT.1) THEN
          NEW = 1
      ELSE
          NEW = 2
      END IF
*
*       Identify global index of new reference body (including binary c.m.).
      DO 45 J = IFIRST,NTOT
          IF (NAME(J).EQ.NAMEC(NEW)) THEN
              ICH = J
              GO TO 50
          END IF
   45 CONTINUE
*
*       Include warning if no reference body (this should not occur).
      WRITE (6,48)  IESC, NAMEC(NEW)
   48 FORMAT (' REDUCE:   DANGER!   NO REFERENCE BODY    IESC NAME',2I5)
      NCH = NCH + 1
      NN = NCH
      MASS = MASS + M(IESC)
      GO TO 100
*
*       Restore ghost to lists of neighbours (body #I will be skipped).
   50 NNB = LIST(1,I)
      DO 52 L = 2,NNB+1
          JPERT(L-1) = LIST(L,I)
   52 CONTINUE
      JLIST(1) = ICH
      CALL NBREST(I,1,NNB)
*
*     IF (KZ(30).GT.1) THEN
*         WRITE (6,53)  NAME0, NAME(ICH), ICH
*  53     FORMAT (' REDUCE:    SWITCH C.M.    NAME0 NAMECH ICH ',3I7)
*     END IF
*
*       Exchange name of reference body and initialize new c.m. name.
      NAME(I) = NAME0
      NAME0 = NAME(ICH)
      NAME(ICH) = 0
*
*       Update total mass and initialize new c.m. body variables.
   55 BODY(ICH) = BODYCH - BODYC(IESC)
      CM(7) = BODY(ICH)
      T0(ICH) = TIME
      DO 60 K = 1,3
          X(K,ICH) = XCM(K)
          X0(K,ICH) = XCM(K)
          XDOT(K,ICH) = VCM(K)
          X0DOT(K,ICH) = VCM(K)
   60 CONTINUE
*
*       Save escaper in common6 for possible use as single in CHTERM2.
      JCMAX = I
*       Restore the mass and transform to global coordinates & velocities.
      BODY(I) = BODYC(IESC)
      T0(I) = TIME
      DKP = 0.0
      DECM = 0.0
      DKE = 0.0
      DO 65 K = 1,3
          X(K,I) = X4(K,IESC) + CM(K)
          XDOT(K,I) = XDOT4(K,IESC) + CM(K+3)
          X0(K,I) = X(K,I)
          X0DOT(K,I) = XDOT(K,I)
*       Form kinetic energy terms (definition of DVC needs negative sign).
          DKP = DKP - BODYC(IESC)*DVC(K)*XDOT4(K,IESC)
          DECM = DECM + 0.5*BODY(ICH)*DVC(K)**2
          DKE = DKE + 0.5*BODYC(IESC)*XDOT4(K,IESC)**2
   65 CONTINUE
*       Reduce look-up time to compensate for inactive ghost or merged star.
      IF (KZ(19).GE.3.AND.KZ(11).NE.0) TEV(I) = TIME
*
*       Remove chain (and clump) mass & reference name of escaper.
      DO 70 L = IESC,NCH
          M(L) = M(L+1)
          BODYC(L) = BODYC(L+1)
          NAMEC(L) = NAMEC(L+1)
          INAME(L) = INAME(L+1)
          SIZE(L) = SIZE(L+1)
          ISTAR(L) = ISTAR(L+1)
          SIZE0(L) = SIZE0(L+1)
          ISTAR0(L) = ISTAR0(L+1)
          BODYS(L,ISUB) = BODYS(L+1,ISUB)
          NAMES(L,ISUB) = NAMES(L+1,ISUB)
   70 CONTINUE
*
*       Ensure current coordinates & velocities for chain components.
      CALL XCPRED(0)
*
*     RR2 = 0.0
*     VR2 = 0.0
*     DO 72 K = 1,3
*     RR2 = RR2 + (X4(K,1) - X4(K,2))**2
*     VR2 = VR2 + (XDOT4(K,1) - XDOT4(K,2))**2
*  72 CONTINUE
*     IF (ID.GT.0) EREL = 0.5*VR2 - (BODYC(1)+BODYC(2))/SQRT(RR2)
*     IF (ID.GT.0) WRITE (6,74)  EREL
*  74 FORMAT (' EREL2  ',F12.7)
*       Re-initialize c.m. & perturber list (2nd time on binary escape).
      IF (JESC.LE.0) THEN
          CALL REINIT(ISUB)
      END IF
*
*       Copy new chain coordinates & velocities to standard variables.
      LK = 0
      DO 80 L = 1,NCH
          DO 75 K = 1,3
              LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
   75     CONTINUE
   80 CONTINUE
*
*     EN0 = ENERGY
*       Update total energy using regular expression (checked OK).
      ENERGY = ENERGY + DECM - DKP + DPOT - DKE
*     IF (ID.GT.0) THEN
*     WRITE (6,81)  DECM, DKP, DPOT, DKE, EN0-ENERGY, ENERGY
*  81 FORMAT (' CHECK!   DECM DKP DPOT DKE DE ENER  ',6F10.6)
*     CALL CONST(XCH,VCH,M,NCH,ENER1,G,ALAG)
*     END IF
*
*       Obtain consistent value of gravitational radius for decision-making.
      SUM = 0.0
      DO 83 J = 1,NCH-1
          DO 82 L = J+1,NCH
              SUM = SUM + M(J)*M(L)
   82     CONTINUE
   83 CONTINUE
      RGRAV = SUM/ABS(ENERGY)
*
*       Copy neighbour list of cm. body for routine NBREST.
      NNB = LIST(1,ICH)
      DO 84 L = 2,NNB+1
          JPERT(L-1) = LIST(L,ICH)
   84 CONTINUE
*
*       Restore ghost to lists of neighbours (body #ICH will be skipped).
      JLIST(1) = I
      CALL NBREST(ICH,1,NNB)
*
*       Make new neighbour list for escaper.
      CALL NBLIST(I,RS0)
*
*      Distinguish between single particle(s) and binary (JESC = 0 or > 0).
      IF (JESC.EQ.0.OR.ISING.GT.0) THEN
*       Initialize force polynomials & time-steps (add differential F & FD).
          IPHASE = 8
          CALL FPOLY1(I,I,0)
          DO 85 K = 1,3
              FIRR(K) = 0.0
              FD(K) = 0.0
   85     CONTINUE
          CALL FCHAIN(I,0,X(1,I),XDOT(1,I),FIRR,FD)
          RIJ2 = 0.0
          VIJ2 = 0.0
          RDOT = 0.0
          DO 86 K = 1,3
              F(K,I) = F(K,I) + FIRR(K)
              FDOT(K,I) = FDOT(K,I) + FD(K)
              RIJ2 = RIJ2 + (X(K,I) - X(K,ICH))**2
              VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,ICH))**2
              RDOT = RDOT + (X(K,I) - X(K,ICH))*(XDOT(K,I)-XDOT(K,ICH))
   86     CONTINUE
          CALL FPOLY2(I,I,0)
*       Check high-velocity ejection for outwards hyperbolic motion.
          HI = 0.5*VIJ2 - (BODY(I) + BODY(ICH))/SQRT(RIJ2)
          IF (HI.GT.0.0.AND.RDOT.GT.0.0.AND.KZ(37).GT.0) THEN
               CALL HIVEL(I)
          END IF
          IF (ISING.EQ.1) THEN
              IF (JESC.EQ.0) GO TO 100
              IF (JESC.GT.IESC) JESC = JESC - 1
              IESC = JESC
              JESC = -1
              ISING = 2
              GO TO 2
          END IF
*      Ensure single escape removal for original JESC = 0.
          IF (JESC.EQ.0) JESC = -1
          IPHASE = -1
*       Include case of escaping binary (set JESC < 0 for new KS).
      ELSE IF (JESC.GT.0) THEN
          ICLOSE = I
          IF (JESC.GT.IESC) JESC = JESC - 1
          IESC = JESC
          JESC = -1
          ISING = 0
          GO TO 2
      ELSE IF (ISING.EQ.0) THEN
*       Initialize KS regularization after second reduction.
          ICOMP = MIN(ICLOSE,I)
          JCOMP = MAX(ICLOSE,I)
*       Obtain new neighbour list to include current chain c.m. body.
          CALL NBLIST(ICOMP,RS0)
          IPHASE = 2
          CALL KSREG
          IPHASE = -1
          CALL CHLIST(ICH)
      END IF
*
*       Re-activate any dormant binary.
      IF (I.GT.N.AND.JESC.LE.0) THEN
          CALL RENEW(I)
          NAMEC(10) = 0
      END IF
*
*       Include experimental procedure to catch new perturbers.
      NP = LISTC(1)
      JP = I
      NNB = 0
      DO 230 L = 2,NNB+1
          I = LIST(L,ICH)
          IF (I.EQ.JP) GO TO 230
          FX1 = F(1,I)
          DO 210 LL = 2,NP+1
              IF (I.EQ.LISTC(LL)) GO TO 230
  210     CONTINUE
          SS = -1.0
          ZMS = ZM0S
          DO 212 K = 1,3
              XS(K) = X0S(K)
              VS(K) = V0S(K)
  212     CONTINUE
          CALL FFDOT(I,ZMS,XS,VS,SS)
          ZMS = BODY(ICH)
          DO 215 K = 1,3
              XS(K) = X(K,ICH)
              VS(K) = XDOT(K,ICH)
  215     CONTINUE
          SS = 1.0
          CALL FFDOT(I,ZMS,XS,VS,SS)
          ZMS = BODY(JP)
          DO 220 K = 1,3
              XS(K) = X(K,JP)
              VS(K) = XDOT(K,JP)
  220     CONTINUE
          CALL FFDOT(I,ZMS,XS,VS,SS)
      RIJ2 = 0.0
      DO K = 1,3
      RIJ2 = RIJ2 + (X(K,I) - XS(K))**2
      END DO
      RIJ = SQRT(RIJ2)
      WRITE (6,216)  NAME(I), FX1, F(1,I)-FX1, STEP(I), RIJ
  216 FORMAT (' FFDOT    NM FX0 DFX S RIJ ',I6,1P,4E10.2)
  230 CONTINUE
*
*       Check for possible escape of chain c.m. by strong recoil.
      IF (JESC.EQ.0.AND.HI.GT.0.0) THEN
          VX = SQRT(2.0*HI)*VSTAR
          V1 = VX*BODY(I)/BODYCH
          V2 = VX*BODY(ICH)/BODYCH
*       Evaluate relativistic time-scale (CVEL = 0 is OK).
          ECC2 = ECCGR**2
          FE = 1.0 + (73.0/24.0 + 37.0*ECC2/96.0)*ECC2
          GE = (1.0 - ECC2)**3.5/FE
          ZMX = MAX(M(I1),M(I2))
          RATIO = MIN(M(I1),M(I2))/ZMX
*       Replace physical time-scale by N-body units (cf. Lee 1993).
*         TZ = 5/64*TAUGR*GE*SEMIGR**4/(RATIO*(1.0 + RATIO)*ZMX**3)
          TZ = 5.0/64.0*GE*SEMIGR**4/(RATIO*(1.0 + RATIO)*ZMX**3)
*       Adopt real value of C for Newtonian case (19/8/11).
          IF (CVEL.GT.0.0) THEN
              TZ = CVEL**5*TZ
          ELSE
              CZ = 3.0D+05/VSTAR
              TZ = CZ**5*TZ
          END IF
*       Estimate escape velocity (or ejection) from cluster for chain c.m.
          IF (V1.GT.3.0*VRMS) THEN
              RIJ = SQRT(RIJ2)
              WRITE (6,240)  NCH, V1, V2, BODY(ICH)*SMU, BODY(I)*SMU,
     &                       TZ, NAME(I), (NAMEC(K),K=1,NCH)
  240         FORMAT (' CHAIN ESCAPE    NCH V1 V2 MC MI TZ NAMI NAMC ',
     &                                  I4,2F7.1,2F6.1,1P,E9.1,0P,6I7)
              CALL CONST(XCH,VCH,M,NCH,ENER1,G,ALAG)
*       Form the net energy change (note ECH = ENERGY - DEGR).
              DE = ENER1 - (ECH + DEGR)   ! opposite DEGR sign used 18/6/10.
              EB = HI*BODY(I)*BODY(ICH)/BODYCH
*       Evaluate tidal energy correction wrt single BH.
              DPHI = DPOT - BODY(I)*BODY(ICH)/RIJ
              WRITE (6,245)  ECH, EB, DE, DEGR, DPHI, ENER1-ENERGY,
     &                       SEMIGR, ECCGR
  245         FORMAT (' ENERGY CHECK    ECH EB DE DEGR DPHI ERROR A E ',
     &                                  4F12.6,1P,3E10.2,0P,F9.5)
*       Subtract the two-body part to compensate total energy budget.
              ECOLL = ECOLL + (ENER1 - DEGR) - DPHI
*       Clean current chain variables to allow another new case.
              ECH = 0.0
*       Set large time to prevent CALL CHAIN while c.m. is escaping.
              TS(ISUB) = 1.0D+10
*       Specify zero indicator for immediate exit from chain with NCH = 0.
*             ISUB = 0
*             NSUB = MAX(NSUB - 1,0)
*             NCH = 0             ! Elegance sacrificed for possible danger.
*       Redefine chain c.m. as single to allow new case.
              NAME(ICH) = 99999
              IF (KZ(37).GT.0) CALL HIVEL(ICH)
          ELSE IF (V2.GT.2.0*VRMS) THEN
              CALL CONST(XCH,VCH,M,NCH,ENER1,G,ALAG)
*       Set net energy gain (Newtonian value).
              DE = ENER1 - ECH  ! conceptually wrong to use - (ECH - DEGR).
              WRITE (6,250)  TTOT, NCH, V1, V2, BODY(I)*SMU, DE,
     &                       NAME(I), (NAMEC(K),K=1,NCH)
  250         FORMAT (' FAST ESCAPE    T NCH V1 V2 MI DE NAMI NAMC ',
     &                                 F8.1,I4,2F7.1,F6.1,F10.5,9I7)
              WRITE (6,252) VSTAR, VRMS,  ENER1, DEGR, ENER1-DEGR,
     &                      SEMIGR, ECCGR
  252         FORMAT (' VSTAR VRMS ENER1 DEGR EN-DEGR A E  ',
     &                                   2F7.2,2F10.5,1P,2E10.2,0P,F9.5)
              WRITE (6,254)  ENERGY, ECH, TZ
  254         FORMAT (' ENERGY ECH TZ  ',2F10.6,1P,E10.2)
          END IF
      END IF
*
*       Re-activate any dormant binary.
*     IF (I.GT.N.AND.JESC.EQ.0) THEN
*         CALL RENEW(I)
*     END IF
*
*       Prepare c.m. check.
*     DO 88 K = 1,6
*         CG(K) = 0.0
*  88 CONTINUE
*
*     LK = 0
*     DO 95 L = 1,NCH
*         DO 90 K = 1,3
*             LK = LK + 1
*             CG(K) = CG(K) + BODYC(L)*XCH(LK)
*             CG(K+3) = CG(K+3) + BODYC(L)*VCH(LK)
*  90     CONTINUE
*  95 CONTINUE
*
*     DO 96 K = 1,6
*         CG(K) = CG(K)/BODY(ICH)
*  96 CONTINUE
*
*     IF (KZ(30).GT.2) THEN
*         WRITE (6,97)  TIME+TOFF, (CG(K),K=1,6)
*  97     FORMAT (' REDUCE:   T CG ',F10.5,1P,6E9.1)
*     END IF
*
*       Ensure zero second indicator for next call of single escaper.
      IF (JESC.LT.0) JESC = 0
      TIME = TIME0
*
  100 RETURN
*
      END
