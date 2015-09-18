      SUBROUTINE CHAIN(ISUB)
*
*       Algorithmic Regularization.
*       ---------------------------
*
*       Method of Mikkola & Merritt (MN 372, 219, 2006)
*       ...............................................
*
*       Regularization routines coded by Seppo Mikkola
*       ..............................................
*
        INCLUDE 'ARCCOM2e2.CH'
        COMMON/DIAGNOSTICS/GAMMA,H,IWR
        common/justforfun/Tkin,Upot
        common/collision/icollision,IBH,JBH,iwarning
        COMMON/CLUMP/  BODYS(NMX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                 NAMES(NMX,5),ISYS(5)
        COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                  ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,JCOLL,NDISS1
        COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),
     &                 NSTEP1,KZ27,KZ30
        COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
        COMMON/POSTN/  CVEL,TAUGR,RZ,GAMMAZ,TKOZ,EMAX,TSP,KZ24,IGR,IPN
        COMMON/POSTN2/ SEMIGR,ECCGR,DEGR,ISPIN
        COMMON/ECHAIN/ ECH
        COMMON/SOFT/  EPS2
        COMMON/EXTRA2/  INJ  ! maybe for later
        EXTERNAL CHMOD
        REAL*8  G0(3),XREL(3),VREL(3),XCM(3),VCM(3),XX(3,3),VV(3,3),
     &          CMXX(3),CMVX(3)
        DATA IEND,ICHECK,IT /0,0,0/
        LOGICAL NEWREG
        SAVE
*
*
      ITERM = ISUB
      IF (ISUB.GT.0) THEN
*       Choose small step for termination (routine INTGRT).
          IF (STEPS(ISUB).EQ.0.0D0) THEN
***           STEP = 1.0D-06*STEP
              GO TO 100
          END IF
*       Update maximum prediction interval at start of every call.
          CALL TCHAIN(ISUB,TSMIN)
          STEPS(ISUB) = TSMIN
*       Synchronize next time interval with subsystem step.
          TMAX = TIMEC + STEPS(ISUB)
          GO TO 100
      END IF
*
*       Copy initial conditions from N-body COMMON and prepare chain.
      TIMEC = 0.0
      CALL CHINIT(ISUB)
*       Read velocity of light and disruption option (first time only).
      IF (IEND.EQ.0) THEN
          READ (5,*) Clight, NBH, IDIS
          CVEL = CLIGHT
          IEND = 1
      END IF
*
      INJ = 0
      ICOAL = 0
      ICOLL = 0
*       Skip star - BH collision search for pure BH treatment.
      IF (IDIS.EQ.0) ICOLL = -1
      JCOLL = 0
      CHTIME = 0.0
      ISYS(5) = ISUB
      ESUM = 0.0
      icollision = 0
      TIME = 0.0
*       Copy Clight into dummy of /POSTN/ (needed elsewhere).
      CVEL = Clight
      NSTEP1 = 0
*       Specify method coefficients (suggested by Seppo).
      cmethod(1) = 1.0
      cmethod(2) = 1.0D-20
      cmethod(3) = 0.0
*       Initialize the spin (gopu called but nothing happens if zero).
      DO K=1,3
      spin(K) = 0.0
      END DO
      spin(1) = 0.0
      spin(2) = 0.0
*     spin(3) = 0.999
      ISPIN = 1
      ISPIN = 0
      JGR = 0
      IBH = 0
      JBH = 0
      I2BH = 0
      J2BH = 0
      ISTAB = 0
      IESC = 0
      JESC = 0
      IPN = 0
      IGR = 0
      IEI = 0
      TZ = 1.0D+04
      TSTAB = 1.0D+06
      TKOZ = 0.0
      TWOPI = 8.0*ATAN(1.0D0)
      IF (CLIGHT.EQ.0.0D0) THEN
          NBH = 0
      END IF
      NBH2 = 0
      DEGR = 0.0
      WTTL = 0.0
      EnerGR = 0.0
      EPS = 1.0D-10
      tolerance = EPS
      RCOLL = 0.0
      RSUB = 0.0
      ESUB = 0.0
      ECOLL1 = 0.0
      IMOVE = 0
*
*       Prepare next step initially or after membership change.
   30 CONTINUE
      TMAX = TIMEC + STEPS(ISUB)
      IWR=-1 ! write some info (set -1 for no diagnostics)
*
      CALL CONST(X,V,M,N,ENER0,G0,AL)
      IF (TIME.EQ.0.0D0.OR.N.EQ.2) THEN
          SUM = 0.0
          MASS = M(N)
          DO 50 I = 1,N-1
              MASS = MASS + M(I)
              DO 45 L = I+1,N
                  SUM = SUM + M(I)*M(L)
   45         CONTINUE
   50     CONTINUE
          RGRAV = SUM/ABS(ENER0)
*       Set provisional dominant elements until routine REDUCE.
          SEMIGR = 0.5*RGRAV
          ECCGR = 0.0
*       Rename softening to sft to avoid clash with old variable.
*         sft=1.0D-03*RGRAV
          sft = 1.0D-20
          ee=sft**2 ! square of softening
          EPS2 = ee
          ENERGY = ENER0
      ENINIT = ENERGY
      WRITE (6,47)  ENERGY
   47 FORMAT (' INITIAL EN   ',F10.6)
      END IF
*
*     TCR=MASS**2.5/(2.*ABS(ENER0))**1.5
      KSMX=100000 ! only this many steps without return
c     Ixc=1 ! 1 for exact time, 0 for not exact time
      Ixc=1 ! activated for new version but no iteration
      NEWREG = .TRUE.
      KCASE = 0
      DO K = 1,3
          CMXX(K) = 0.0
          CMVX(K) = 0.0
      END DO
*
*       Begin the main loop for the block-step interval DELTAT.
  100 DELTAT = STEPS(ISUB)
*       Update energy budget on each NEWREG.
      IF (NEWREG.AND.Clight.GT.0.0) THEN
          IF (IPN.GT.0.AND.NSTEP1.LT.10) THEN
          WRITE (6,102)  N, NPERT, IPN, ENERGY, RGRAV
  102     FORMAT (' NEWREG    N NP IPN ENERGY RGRAV ',
     &                        3I3,2F10.6)
          END IF
      END IF
*       Check termination (positive energy possible without member change).
      IF (N.EQ.2.AND.ENERGY.GT.0.0) THEN
          CALL CONST(X,V,M,N,ENER0,G0,AL)
          ECH = ENERGY
          GO TO 250
      END IF
      DTREM = TMAX - TIMEC
      DELTAT = MIN(DTREM,DELTAT)
      EPREV = ENERGY
      IF (N.EQ.2) RSUM = SEMIGR
      IF (IGR.EQ.0) THEN
          CVEL = 0.0
      ELSE
          CVEL = Clight
      END IF
      IF (ICOAL.GT.0) THEN
          ICOAL = 0
          NEWREG = .TRUE.
      END IF
*
*       Re-determine active GR pointers after absorb or escape.
      IF (IGR.GT.0.AND.NEWREG) THEN
          FX = 0.0
          DO 110 I = 1,N-1
              LI = 3*(I - 1)
              DO 105 J = I+1,N
                  LJ = 3*(J - 1)
                  RIJ2 = (X(LI+1)-X(LJ+1))**2 + (X(LI+2)-X(LJ+2))**2
     &                                        + (X(LI+3)-X(LJ+3))**2
                  FF = (M(I) + M(J))/RIJ2
                  IF (FF.GT.FX) THEN
                      FX = FF
                      I1 = I
                      I2 = J
                  END IF
  105         CONTINUE
  110     CONTINUE
          IBH = MIN(I1,I2)
          JBH = MAX(I1,I2)
      END IF
*
*       Omit higher orders for nearly isolated binary during final stages.
      IF (N.EQ.2.AND.TZ.LT.1.0.AND.GPERT.LT.1.0D-07) THEN
          IPN = 1
      END IF
*
*       Perform the next integration step (note Clight changes with CVEL).
      COPYC = Clight
*       Ensure PN is active after unperturbed KS (IPN set in BRAKE4).
      IF (IPN.GT.0.AND.IPN.LE.3) THEN
          CVEL = CLIGHT
          IGR = 1
      END IF
      IF (JGR.GT.0) IBH = 0
      CALL ARC(N,X,V,M,TIME,DELTAT,EPS,NEWREG,KSMX,sft,cvel,Ixc,NBH,
     &         spin,CMXX,CMVX)
      Clight = COPYC
*
      IF (ISPIN.GT.0.AND.IGR.GT.0.AND.MOD(NSTEP1,1000).EQ.0) THEN
          CALL CONST(X,V,M,N,ENER0,G0,AL)
          ERR = (ENERGY - ENER0)/ENER0
          SS = SQRT(spin(1)**2 + spin(2)**2 + spin(3)**2)
          IF (IT.EQ.0) THEN
              SS0 = SS
              IT = 1
          ELSE
          DS = (SS - SS0)/SS0
          END IF
          WRITE (53,120)  TNOW, 1.0/RINV(1), ERR, spin, DS
  120     FORMAT (' SPIN    T R DE/E spin DS/S ',
     &                      F10.4,1P,3E10.2,3E12.4,E10.2)
          CALL FLUSH(53)
      END IF
      TIMEC = TIME
      TNOW = TSP + TIMEC
      NSTEP1 = NSTEP1 + 1
*
*       Check movie output.
*     IF (IMOVE.GT.0.AND.TMOVE.LT.100.0) THEN
*     IF (TIMEC.GT.TMOVE) THEN
*         CALL MOVIE_DATA(TNOW)
*         TMOVE = TIMEC + DTMOVE
*     END IF
*     END IF
*
*       Include extra BH information at late stages (when CALL REDUCE rare).
      IF (IPN.GE.2.AND.MOD(NSTEP1,100).EQ.0) THEN
          CALL BHSTAT
      END IF
*
*       Predict perturbers & XC, UC and form new LISTC every 10 steps.
      IF (MOD(NSTEP1,10).EQ.0) THEN
          JJ = 0
          CALL XCPRED(2)
          CALL CHLIST(JJ)
      ELSE
*       Perform fast prediction of XC & UC every step (#ICH in INTGRT).
          CALL XCPRED(0)
      END IF
*
*       Specify indicator for absorber check (only one more member).
      IF (GPERT.GT.0.05.AND.N.LT.4) INJ = -1
*       Check indicator for membership injection.
      IF (INJ.LT.0) THEN
          N0 = N
          CALL INJECT(ISUB)
          INJ = 0
          IBH = -1
          IF (N.GT.N0) GO TO 30
      END IF
      ESUM = ESUM + (ENERGY - EPREV)
      IF (NSTEP1.GT.2000000000) NSTEP1 = 0
*
*       Locate index of most massive body (save MX for later).
      MX = 0.0
      DO 130 L = 1,N
          IF (M(L).GT.MX) THEN
              MX = M(L)
              LX = L
          END IF
  130 CONTINUE
*
*       Set relevant coalescence (4*R_Sch) even with disruption.
      IF (CLIGHT.GT.0.0) THEN
          IF (IBH.GT.0) THEN
              RZ = 8.0*(M(IBH) + M(JBH))/CLIGHT**2
          ELSE
              RZ = 8.0*M(LX)/CLIGHT**2
          END IF
      ELSE
*       Do not allow coalescence if CLIGHT is inactive.
          RZ = 0.0
      END IF
*
*       Search for the closest binary BH-BH and/or star-BH pair.
      AINX = -1.0D+04
      RDIS2 = 1.0
      RPERT = 1.0
      RX = 0.0
      I1 = 0
      DO 135 I = 1,N-1
          RPERT = MAX(1.0/RINV(I),RPERT)
          IF (1.0/RINV(I).GT.RX) THEN
              RX = 1.0/RINV(I)
          END IF
          LI = 3*(I - 1)
          DO 134 J = I+1,N
              LJ = 3*(J - 1)
              RIJ2 = (X(LI+1)-X(LJ+1))**2 + (X(LI+2)-X(LJ+2))**2
     &                                    + (X(LI+3)-X(LJ+3))**2
              VIJ2 = (V(LI+1)-V(LJ+1))**2 + (V(LI+2)-V(LJ+2))**2
     &                                    + (V(LI+3)-V(LJ+3))**2
              AIN = 2.0/SQRT(RIJ2) - VIJ2/(M(I) + M(J))
              IF (AIN.GT.AINX.AND.
     &            ISTAR(I).EQ.14.AND.ISTAR(J).EQ.14) THEN
                  AINX = AIN
                  I1 = I
                  I2 = J
*       Consider black hole - single star encounter.
              ELSE IF (IDIS.GT.0.AND.RIJ2.LT.RDIS2.AND.
     &             AIN.GT.AINX.AND.
     &             ((ISTAR(I).EQ.14.AND.ISTAR(J).NE.14).OR.
     &              (ISTAR(J).EQ.14.AND.ISTAR(I).NE.14))) THEN
*       Note rare case of two stars inside RCOLL is skipped below.
                  RDIS2 = RIJ2
                  RDIS = SQRT(RDIS2)
                  RD = 0.0
                  VIJ2 = 0.0
                  DO 133 K = 1,3
                      VIJ2 = VIJ2 + (V(LI+K) - V(LJ+K))**2
                      RD = RD + (X(LI+K) - X(LJ+K))*(V(LI+K) - V(LJ+K))
  133             CONTINUE
                  ADIS = 2.0/RDIS - VIJ2/(M(I) + M(J))
                  ADIS = 1.0/ADIS
                  EDIS = (1.0-RDIS/ADIS)**2 + RD**2/((M(I)+M(J))*ADIS)
                  EDIS = SQRT(EDIS)
                  SZ = MAX(SIZE(I),SIZE(J))
                  RATIO = MAX(M(I),M(J))/MIN(M(I),M(J))
                  RCOLL = RATIO**0.3333*SZ
*       Include factor of 1000 to eliminate WD subsystem (too expensive).
*                 IF (MIN(ISTAR(I),ISTAR(J)).GE.10) RCOLL = 1000.*RCOLL
*       Check disruption distance for pericentre or actual separation.
                  PMDIS = ADIS*(1.0 - EDIS)
                  IF (RDIS.GT.0.1*RPERT) PMDIS = RDIS
                  IF (ICOLL.EQ.0.AND.PMDIS.LT.RCOLL.AND.RD.LT.0.0) THEN
                      ICOLL = I
                      JCOLL = J
                  END IF
                  I1 = I
                  I2 = J
                  AINX = AIN
*       Note final check condition for ordinary stars.
              ELSE IF (AIN.GT.AINX) THEN
                  AINX = AIN
                  RDIS2 = RIJ2
                  RDIS = SQRT(RDIS2)
                  I1 = I
                  I2 = J
              END IF
  134     CONTINUE
  135 CONTINUE
*       Include safety check (should not occur).
      IF (I1.EQ.0) THEN
          I1 = 1
          I2 = 2
      END IF
      IBH = I1
      JBH = I2
*
*       Form classical two-body elements for dominant interaction I1 & I2.
      RIJ2 = 0.0
      VIJ2 = 0.0
      RDOT = 0.0
      KI = 3*(I1 - 1)
      KN = 3*(I2 - 1)
      DO 140 K = 1,3
          XREL(K) = X(K+KI) - X(K+KN)
          VREL(K) = V(K+KI) - V(K+KN)
          RIJ2 = RIJ2 + XREL(K)**2
          VIJ2 = VIJ2 + VREL(K)**2
          RDOT = RDOT + XREL(K)*VREL(K)
  140 CONTINUE
      R12 = SQRT(RIJ2)
      SEMI = 2.0/R12 - VIJ2/(M(I1) + M(I2))
      SEMI = 1.0/SEMI
      ECC2 = (1.0 - R12/SEMI)**2 + RDOT**2/(SEMI*(M(I1) + M(I2)))
      ECC = SQRT(ECC2)
      IF (SEMI.LT.0.0) THEN
          SEMI = MIN(0.5*RSUM,SEMIGR)
          ECC = 0.5
      END IF
*
      PMIN = SEMI*(1.0 - ECC)
      SEMI0 = SEMI
      ECC0 = ECC
      IF (SEMI.GT.0.0.AND.IPN.EQ.0) THEN
          ECCGR = ECC
          SEMIGR = SEMI
      END IF
*
*       Obtain relativistic elements or velocity ratio.
      IF ((SEMI.GT.0.0.AND.IPN.GT.0)) THEN
*       Evaluate relativistic elements.
          CALL GRBIN(M(I1),M(I2),XREL,VREL,SEMI,ECC)
          PMIN = SEMI*(1.0 - ECC)
          ECC2 = ECC**2
*       Update GR elements for CHLIST (otherwise only in REDUCE; 23/8/11).
          IF (SEMI.GT.0.0) THEN
              ECCGR = ECC
              SEMIGR = SEMI
          END IF
      ELSE IF (CLIGHT.GT.0.0.AND.PMIN.LT.100.0*RZ) THEN
          VC2 = VIJ2/CLIGHT**2
          IF (VC2.GT.1.0D-06.AND.DW.LT.1.0D-03) THEN
              IGR = 1
              IPN = 1
          ELSE
              IGR = 0
          END IF
      ELSE
          IGR = 0
      END IF
      IF (IGR.EQ.0) IPN = 0
      JGR = 0        ! use JGR > 0 for switching off PN (experimental).
      IF (SEMI.LT.0.0) THEN
          IGR = 0
          IPN = 0
      END IF
*
*       Evaluate the Einstein shift per orbit and check IPN.
      IF (IPN.LE.1.AND.SEMI.GT.0.0) THEN
          DW = 3.0*TWOPI*(M(I1) + M(I2))/(SEMI*Clight**2*(1.0 - ECC2))
*       Ensure IPN activated if shift exceeds 1.0D-04 per orbit.
          IX = MAX(ISTAR(I1),ISTAR(I2))
          IF (DW.GT.1.0D-04.AND.IX.EQ.14) THEN
              IPN = 1
              IGR = 1
              IEI = 1   ! Einstein shift indicator (8/14).
*       Allow for 2nd order correction (Mikkola & Merritt ApJ 135, 2398).
              IF (DW.GT.1.0D-03) IPN = 2   ! 2nd order is 5 times bigger.
              IF (DW.GT.1.0D-03) THEN
                  IF (DW.GT.1.0D-02) IPN = 3
                  DW = DW*(1.0 + 5.0*DW)
                  WRITE (6,142)  IPN, IX, ECC, SEMI, DW
  142             FORMAT  (' EINSTEIN SHIFT    IPN IX E A DW ',
     &                                         2I4,F9.5,1P,2E10.2)
              END IF
          ELSE
              IEI = 0
          END IF
      END IF
*
      IF (IPN.GE.2.AND.MOD(NSTEP1,100).EQ.0) THEN
          CALL CONST(X,V,M,N,ENER0,G0,AL)
          WRITE (6,144)  NSTEP1, IPN, ECC, ENERGY, ENER0, SEMI, GPERT,
     &                   ESUM, (1.0/RINV(K),K=1,N-1)
  144     FORMAT (' WATCH    # IPN E EN EN0 A GP ES R ',
     &                       I8,I4,F8.4,2F12.7,1X,1P,E12.4,7E12.2)
      CALL FLUSH(6)
      END IF
*
*       Determine c.m. of dominant pair and closest chain member.
      IF (N.GT.2) THEN
          MB = M(I1) + M(I2)
          K1 = 3*(I1 - 1)
          K2 = 3*(I2 - 1)
          DO 200 K = 1,3
              XCM(K) = (M(I1)*X(K+K1) + M(I2)*X(K+K2))/MB
              VCM(K) = (M(I1)*V(K+K1) + M(I2)*V(K+K2))/MB
  200     CONTINUE
          RX2 = 1.0
          DO 205 I = 1,N
              IF (I.EQ.I1.OR.I.EQ.I2) GO TO 205
              RIJ2 = 0.0
              LI = 3*(I - 1)
              DO 202 K = 1,3
                  RIJ2 = RIJ2 + (X(K+LI) - XCM(K))**2
  202         CONTINUE
              IF (RIJ2.LT.RX2) THEN
                  RX2 = RIJ2
                  IM = I
              END IF
  205     CONTINUE
*       Form hierarchical elements and Kozai period.
          RIJ2 = 0.0
          VIJ2 = 0.0
          RRD = 0.0
          LX = 3*(IM - 1)
          DO 210 K = 1,3
              RIJ2 = RIJ2 + (X(K+LX) - XCM(K))**2
              VIJ2 = VIJ2 + (V(K+LX) - VCM(K))**2
              RRD = RRD + (X(K+LX) - XCM(K))*(V(K+LX) - VCM(K))
  210     CONTINUE
          RCJ = SQRT(RIJ2)
          AOUT = 2.0/RCJ - VIJ2/(MB + M(IM))
          IF (SEMI.GT.0.0.AND.AOUT.GT.0.0) THEN
              AOUT = 1.0/AOUT
              ECC1 = (1.0 - RCJ/AOUT)**2 + RRD**2/(AOUT*(MB + M(IM)))
              TIN = TWOPI*SEMI*SQRT(SEMI/MB)
              TOUT = TWOPI*AOUT*SQRT(AOUT/(MB + M(IM)))
*       Note small correction from (1 - e1) to (1 - e1**2) 6/2012.
              TKOZ = TOUT**2/TIN*(1.0 - ECC1)**1.5*MB/M(IM)
*       Include numerical factor quoted by Fabrycky & Tremaine 2007 (9/11).
              TKOZ = 2.0/(3.0*3.14)*TKOZ  ! Kiseleva et al MN 300, 292, 1998.
*       Check stability criterion near apocentre during GR or TOUT/TIN > 5.
              IF (TSTAB.GT.20000.0.AND.
     &            (IPN.GT.0.OR.TOUT.GT.5.0*TIN)) THEN
                  TSTAB = TNOW
              END IF
              IF (IPN.EQ.0.AND.TOUT.LT.5.0*TIN) TSTAB = 1.0D+06
              IF (N.EQ.3.AND.NPERT.EQ.0.AND.IPN.GT.0.AND.
     &            TNOW.GE.TSTAB.AND.R12.GT.0.9*SEMI*(1.0+ECC)) THEN
                  TSTAB = TSTAB + 100.0*TOUT
                  DO 215 K = 1,3
                      J1 = K1 + K
                      J2 = K2 + K
                      J3 = LX + K
                      XX(K,1) = X(J1)
                      XX(K,2) = X(J2)
                      XX(K,3) = X(J3)
                      VV(K,1) = V(J1)
                      VV(K,2) = V(J2)
                      VV(K,3) = V(J3)
  215             CONTINUE
*       Obtain the inclination & EMAX and perform stability test.
                  CALL INCLIN(XX,VV,XCM,VCM,ALPH)
                  QM = M(IM)/(MB + M(IM))
                  E1 = SQRT(ECC1)
                  NST = NSTAB(SEMI,AOUT,ECC,E1,ALPH,M(I1),M(I2),M(IM))
*       Note ECCGR replaced by ECC 3/7/11 because reduce.f called rarely.
                  FAC = (1.0 + QM)*(1.0 + E1)/SQRT(1.0 - E1)
                  RPC = 2.8*FAC**0.4*SEMI
                  ALPH = 180.0*ALPH/3.1415
*       Evaluate perturbation at mean separation.
                  GA = 2.0*M(IM)/MB*(SEMI/RCJ)**3
                  CALL EMAX1(MB,XX,VV,XCM,VCM,ECC2,EX,EM)
                  IF (NST.EQ.0.AND.MOD(NSTEP1,100).EQ.0) THEN
                      PM = AOUT*(1.0 - E1)
                      WRITE (97,220) TNOW, IPN, NAMEC(IM), ALPH, EX, EM,
     &                               PM, RPC, TKOZ, GA
  220                 FORMAT (' BHSTAB    T IPN NM IN EX EM PM RPC TK ',
     &                         'GA ',F9.2,I3,I6,F6.1,F9.5,F7.3,1P,4E9.1)
                      CALL FLUSH(97)
*       Check long-lived inclined triples for switching off PN using JGR > 0.
                      IF (IPN.EQ.1.AND.TZ.GT.20.0.AND.ALPH.GT.130.0.AND.
     &                    PM.GT.20.0*SEMI.AND.RRD.GT.0.0) THEN
                          JGR = 1
                          WRITE (6,222)  NSTEP1, ECC, ALPH, PM/SEMI, TZ
  222                     FORMAT (' PN SWITCH-OFF    # E IN PM/A TZ ',
     &                                          I10,F8.4,F7.1,1P,2E9.1)
*       Terminate chain by existing procedure.
                          GO TO 250
                      END IF
*       Include termination for weakly perturbed eccentric binary.
                      IF (IPN.EQ.1.AND.ECC.GT.0.99.AND.GA.LT.1.D-08.AND.
     &                    TZ.GT.10.0.AND.RRD.GT.0.0) THEN
                          WRITE (6,222)  NSTEP1, ECC, ALPH, PM/SEMI, TZ
                          GO TO 250
                      END IF
                  ELSE IF (AOUT*(1.0-E1).LT.RPC) THEN
                      WRITE (98,225) TNOW, IPN, NAMEC(IM), ALPH, EX, EM,
     &                               AOUT*(1.0-E1), RPC, GA
  225                 FORMAT (' UNSTAB    T IPN NM IN EX EM PM RPC GA ',
     &                               F9.2,I3,I6,F7.1,F9.5,F7.3,1P,3E9.1)
                      CALL FLUSH(98)
                  END IF
                  WRITE (6,230)  TNOW, IPN, ECC, EX, EM, ALPH, SEMI,
     &                           TZ, TKOZ
  230             FORMAT (' EMAX    T IPN E EX EM IN A TZ TKOZ ',
     &                              F9.2,I3,2F9.5,F8.4,F8.2,1P,3E9.1)
              IF (EX.GT.0.99998.AND.MIN(ISTAR(I1),ISTAR(I2)).GE.13.AND.
     &            TZ.LT.1.0) THEN
                  icollision = 1
                  ICOLL = I1
                  JCOLL = I2
              END IF
              END IF
          ELSE
              TKOZ = 1.0D+04
          END IF
      ELSE
          TKOZ = 1.0D+04
      END IF
*
*       Determine radiation time-scale and corresponding indicators.
      IF ((ECC.LT.1.0.AND.CLIGHT.GT.0.0.AND.ECC.GT.0.97).OR.
     &    (ECC.LT.1.0.AND.TZ.LT.100.0).OR.
     &    (DW.GT.1.0D-04.AND.ECC.LT.1.0)) THEN
          FE = 1.0 + (73.0/24.0 + 37.0*ECC2/96.0)*ECC2
          GE = (1.0 - ECC2)**3.5/FE
          MX = MAX(M(I1),M(I2))
          RATIO = MIN(M(I1),M(I2))/MX
*       Replace physical time-scale by N-body units (cf. Lee 1993).
*         TZ = TAUGR*GE*SEMI**4/(RATIO*(1.0 + RATIO)*MX**3)
          TZ = GE*SEMI**4/(RATIO*(1.0 + RATIO)*MX**3)
          TZ = 5.0/64.0*CLIGHT**5*TZ
          ZN = SQRT(MX/SEMI**3)
          PDOT = 3.0*ZN/(1.0 - ECC2)*MX/(SEMI*CLIGHT**2)
          TPOM = 6.283/PDOT
          IF (NSTEP1.EQ.1.AND.DW.GT.1.0D-04) THEN
              WRITE (6,145)  ECC, SEMI, PMIN, RZ, TZ, TPOM, DW
  145         FORMAT (' RELATIVISTIC    ECC AX PM RZ TZ TPOM DW ',
     &                                  F8.4,1P,6E10.2)
          END IF
*       Specify IGR & IPN according to time-scale (experimental).
          IGR = 1
          IF (TZ.LT.1.0) THEN
              IPN = 3
          ELSE IF (TZ.LT.50.0) THEN
              IPN = 2
          ELSE IF (TZ.LT.500.0.AND.DW.LT.1.0D-03) THEN
              IPN = 1
          ELSE IF (IEI.EQ.0) THEN
              IGR = 0
              IPN = 0
          END IF
*       Reduce GR indicator from 2 to 1 for small GPERT & MIN(TZ,TKOZ) > 10.
          TYZ = MIN(TZ,TKOZ)
*         IF (IPN.EQ.2.AND.GPERT.LT.1.0D-07.AND.TYZ.GT.10.0) THEN
*             IPN = 1
*         END IF
      ELSE
          IGR = 0
          IPN = 0
          CVEL = 0.0
      END IF
*
*       Perform occasional GR check for high eccentricity.
      IF ((IPN.GT.1.AND.ECC.GT.0.99.AND.MOD(NSTEP1,1000).EQ.0).OR.
     &    (IPN.GT.2.AND.ECC.GT.0.999.AND.ICHECK.LT.20000)) THEN
          ICHECK = ICHECK + 1
          WRITE (66,146)  TNOW, ECC, IGR, IPN, NPERT, R12/SEMI, TZ
  146     FORMAT (' GR CHECK    T E IGR IPN NP R/A TZ ',
     &                          F10.3,F9.5,2I3,I4,F7.3,1P,E9.1)
          CALL FLUSH(66)
      END IF
*       Define component indices for GR terms.
      IF (IGR.GT.0) THEN
          IBH = MIN(I1,I2)
          JBH = MAX(I1,I2)
      ELSE
          IBH = 0
          JBH = 0
          TZ = 1.0D+04
      END IF
*
*       Include extra diagnostics during late GR stages (CALL REDUCE rare).
      IF (IPN.GE.2.AND.MOD(NSTEP1,1000).EQ.0.AND.TZ.LT.50) THEN
          WRITE (6,147)  TNOW, NPERT, IPN, ECC, SEMI, TZ
  147     FORMAT (' INSPIRAL    T NP IP E A TZ ',
     &                          F10.3,2I4,F9.5,1P,2E10.2)
      END IF
*
*       Look for additional GR interaction terms (suppressed by IGR.LT.0).
      IF (IGR.LT.0.AND.N.GT.2.AND.IBH.GT.0) THEN
          RY = 1.0
          I = IBH
  148     LI = 3*(I - 1)
          DO 150 J = 1,N
              IF (J.EQ.IBH.OR.J.EQ.JBH) GO TO 150
              LJ = 3*(J - 1)
              RIJ2 = (X(LI+1)-X(LJ+1))**2 + (X(LI+2)-X(LJ+2))**2
     &                                    + (X(LI+3)-X(LJ+3))**2
              IF (RIJ2.LT.RY) THEN
                  RY = RIJ2
                  I2BH = I
                  J2BH = J
              END IF
  150     CONTINUE
          IF (I.EQ.IBH) THEN
              I = JBH
              GO TO 148
          END IF
          RY = SQRT(RY)
*       Consider stellar disruption and accretion of body #J2BH.
          IF (IDIS.GT.0) THEN
              RZ = (M(I2BH)/M(J2BH))**0.3333*SIZE(J2BH)
          END IF
          IF (RY.LT.RZ) THEN
              RX = RY
              JBH = J2BH
              WRITE (6,152)  INAME(J2BH), NAMEC(J2BH), RY, M(IBH),
     &                       M(JBH)
  152         FORMAT (' DISRUPT STAR    INM NM RIJ M12 ',
     &                                  I4,I6,1P,3E10.2)
          END IF
*     ELSE
          I2BH = 0
          J2BH = 0
      END IF
*
*       Exclude star-star collisions in disruption cases (also check for BH).
      IF (IDIS.GT.0.AND.ICOLL.GT.0) THEN
          RY2 = 1.0
          DO 158 I = 1,N-1
              LI = 3*(I - 1)
              DO 156 J = I+1,N
                  LJ = 3*(J - 1)
                  RIJ2 = (X(LI+1)-X(LJ+1))**2 + (X(LI+2)-X(LJ+2))**2
     &                                        + (X(LI+3)-X(LJ+3))**2
*       Avoid two stars being close for star-BH interaction..
                  IF (RIJ2.LT.RY2) THEN
                      RY2 = RIJ2
                      IY = I
                      IZ = J
                  END IF
  156         CONTINUE
  158     CONTINUE
*       Impose necessary condition for BH-star pair being closest.
          IF (IY + IZ.NE.ICOLL + JCOLL) THEN
              ICOLL = 0
              JCOLL = 0
          END IF
*       Ensure that ICOLL OR JCOLL represents the BH.
          IF (ICOLL.GT.0) THEN
              IF (ISTAR(ICOLL).NE.14.AND.ISTAR(JCOLL).NE.14) THEN
                  ICOLL = 0
                  JCOLL = 0
              ELSE
                  ECCGR = ECC
              END IF
          END IF
      END IF
*
*       Check collision for non-BH dominant stars (IBH = 0 if IGR = 0).
      IF (ISTAR(I1).LT.14.AND.ISTAR(I2).LT.14) THEN
          J1 = I1
          IF (SIZE(I1).LT.SIZE(I2)) J1 = I2
          RCOLL1 = 1.4*((M(I1) + M(I2))/M(J1))**0.3333*SIZE(J1)
          IF (PMIN.LT.RCOLL1) GO TO 258
*       Note: needs to be consistent with KSINT criterion (24/7/14).
      END IF
*
*       Perform collision test based on multiple criteria (BH-BH or BH-S).
      IF (icollision.gt.0.and.ICOLL.GT.0) icollision = 0
      IF (icollision.gt.0.OR.
     &   (IPN.GT.1.AND.PMIN.LT.RZ).OR.
     &   (IPN.GE.3.AND.TZ.LT.1.0.AND.NPERT.EQ.0.AND.N.EQ.2).OR.
*      Allow coalescence for wide outer orbit and short GR time-scale.
     &   (N.EQ.3.AND.TZ.LT.1.0.AND.TKOZ.GT.25.0.AND.
     &    AOUT*(1.0-SQRT(ECC1)).GT.100*SEMI).OR.
*      Note osculating orbit within block-step may be strongly perturbed.
     &    ICOLL.GT.0) THEN
          IF (ICOLL.GT.0) THEN
              IBH = ICOLL
              JBH = JCOLL
              ECC = EDIS
              SEMI = ADIS
              GO TO 165
          END IF
*       Ensure N > 2 before coalescence.
          IF (N.GE.2) THEN
*       Include injection for testing purposes to simulate N > 2.
*             CALL INJECT(ISUB)
              IBH = 1
              JBH = 2
*       Switch to termination for two stars (N = 3: CALL REDUCE first).
              IF (ISTAR(IBH).LT.14.AND.ISTAR(JBH).LT.14) GO TO 258
              GO TO 165
          END IF
*       Identify possible missing components.
          IF (IBH.EQ.0.OR.JBH.EQ.0) THEN
              RY = 0.0
              DO 160 K = 1,N-1
                  IF (RINV(K).GT.RY) THEN
                      RY = RINV(K)
                      KK = K
                   END IF
  160         CONTINUE
              IBH = INAME(KK)
              JBH = INAME(KK+1)
          END IF
          IF (IBH.EQ.0) THEN
              IBH = I1
              JBH = I2
          END IF
*       Check that the smallest mass will be absorbed by biggest.
  165     IF (M(JBH).LT.M(IBH)) THEN
              KK = JBH
              JBH = IBH
              IBH = KK
          END IF
          IF (IBH.GT.JBH) THEN
              WRITE (6,166)  IBH, JBH, ISTAR(IBH), ISTAR(JBH),
     &                       M(IBH), M(JBH), TZ
  166         FORMAT (' REVERSE INFALL    IBH JBH I* MI MJ TZ ',
     &                                    4I4,1P,3E10.2)
              KK = IBH
              IBH = MIN(IBH,JBH)
              JBH = MAX(KK,JBH)
          END IF
          EB = -0.5*M(IBH)*M(JBH)/SEMI
*       Define number of black holes in the binary for COAL & kick purposes.
          NBH2 = 0
          IF (ISTAR(IBH).EQ.14) NBH2 = NBH2 + 1
          IF (ISTAR(JBH).EQ.14) NBH2 = NBH2 + 1
*       Enforce termination for two colliding stars (IDIS = 0).
          IF (NBH2.EQ.0.AND.N.EQ.2) THEN
              WRITE (6,167)  ECC, SEMI, N, (ISTAR(K),K=1,N)
  167         FORMAT (' S + S CHAIN TERM    E SEMI N ISTAR ',
     &                                      F9.5,1P,E10.2,0P,5I4)
              CALL CHTERM2(NBH2)
              GO TO 290
          END IF
          IC = icollision
*       Note coalescence assumption of two BHs (decide on BH + NS later).
          IF (IDIS.EQ.0.OR.(IPN.GE.2.AND.NBH2.EQ.2)) THEN
              WRITE (6,168)  TNOW, IC, IPN, NAMEC(IBH), NAMEC(JBH),
     &                       ECC, EB, PMIN, TZ
  168         FORMAT (' COALESCENCE    T IC IPN NAM E EB PM TZ '
     &                                 F10.3,2I3,2I6,F9.5,1P,3E10.2)
*       Reduce NBH2 for single BH binary coalescence (CHTERM2 needs N = 2).
              IF (N.EQ.2) NBH2 = NBH2 - 1
              ICOAL = 1
          ELSE
              IF (N.GT.2) THEN
*       Determine likely perturber mass (RPERT may not be quite right).
                  DO 169 I = 1,N
                      IF (I.NE.IBH.AND.I.NE.JBH) MSTAR = M(I)
  169             CONTINUE
                  GX = 2.0*MSTAR/(M(IBH)+M(JBH))*(RDIS/RPERT)**3
              ELSE
                  GX = 0.0
              END IF
              KS = MIN(ISTAR(ICOLL),ISTAR(JCOLL))
              WRITE (6,170)  TNOW, N, IPN, NAMEC(IBH), NAMEC(JBH), KS,
     &                       ECC, PMDIS, RCOLL, SZ, EB, GX
  170         FORMAT (' DISRUPT    T N IPN NAM K* E PM RCOLL SZ EB GX '
     &                             F10.3,2I3,2I6,I4,F10.6,1P,5E10.2)
              ICOLL = 0
              NBH2 = 0
*       Ensure IBH is BH mass (otherwise segmentation error!). 
              IF (M(IBH).LT.M(JBH)) THEN
                  KK = JBH
                  JBH = IBH
                  IBH = KK
              END IF
          END IF
*
*       Combine components and make new chain with reduced membership.
          CALL INFALL(IBH,JBH,NBH2,ISUB)
          ENER0 = 0.0
          IF (N.GT.1) CALL CONST(X,V,M,N,ENER0,G0,AL)
          icollision = 0
          IBH = -1
          JBH = 0
*       Evaluate energy difference (old - new) for correction purpose.
          DE = ECH - ENER0
       DE = ENINIT - ENER0
*       Include energy loss via ECOLL.
          CALL DECORR(DE)
          WRITE (6,180)  ENERGY, ENER0, DE, ECH
  180     FORMAT (' CHAIN CHECK    ENERGY ENER0 DE ECH ',
     &                             1P,4E12.4)
*       Update current binding energy and initialize if no termination.
          ENERGY = ENER0
          ECH = 0.0
          TZ = 1.0D+04
          IPN = 0
          IGR = 0
*       Continue integration if N > 2.
          IF (N.GT.2) GO TO 30
*       Enforce termination on tidal disruption (avoids small steps).
          IF (IDIS.GT.0.AND.NBH2.LT.2) GO TO 250
*
*       Terminate for N=1 or wide binary (avoids perturber problems).
          IF (N.EQ.1.OR.N.GE.2.AND.ABS(ENER0).LT.0.1*ABS(EB)) THEN
              GO TO 250
          END IF
          DELTAT = TMAX - TIME
          IF (DELTAT.LT.0.0) THEN
              WRITE (6,185)  TSMIN, STEPS(ISUB), DELTAT
  185         FORMAT (' NEGATIVE!    TSMIN SS DT  ',1P,3E10.2)
              DELTAT = 0.001*TSMIN
          END IF
          STEPS(ISUB) = DELTAT
*       Continue reduced chain with new value of ECH.
          ECH = ENERGY
          GO TO 30
      END IF
*
*       Perform three-body stability test every 1000 steps (IPN = 0).
      IF (IPN.EQ.0.AND.N.EQ.3.AND.MOD(NSTEP1,1000).EQ.0) THEN
          CALL CHSTAB(ITERM)
          IF (ITERM.LT.0) GO TO 250
      END IF
*
*       Check enforced termination (N = 2, GPERT < 1D-03, DW < 1.0D-04).
      IF (N.EQ.2.AND.(GPERT.LT.1.0D-03.OR.GPERT.GT.0.5).AND.
     &    (IX.LT.13.OR.DW.LT.1.0D-04)) THEN   ! also include standard stars.
          WRITE (6,188)  IPN, 1.0/RINV(1), PMIN, RZ, GPERT, TZ, DW
  188     FORMAT (' ENFORCED CHAIN TERM    IPN R PM RZ G TZ DW ',
     &                                     I4,1P,6E10.2)
          GO TO 250
      END IF
*
*       Include four-body termination test after 10,000 steps.
      IF (IPN.LE.1.AND.NSTEP1.GT.10000.AND.N.EQ.4) THEN
          IF (MAX(1.0/RINV(1),1.0/RINV(N-1)).GT.50.0*RGRAV) GO TO 250
      END IF
*
      IF (NSTEP1.GT.100000.AND.GPERT.LT.1.0D-07) GO TO 250
*       Consider removal of distant member(s) (consistent with acceptance).
      KCASE = 0
***   IF (N.EQ.3.AND.IPN.GE.2.AND.RX.GT.50.0*SEMI) GO TO 300
*       Define perturbation-dependent scaling factor.
      XFAC = 5.0
      IF (GPERT.LT.1.0D-06) XFAC = 10.0
      IF (NPERT.EQ.0.OR.GPERT.LT.1.0D-08) XFAC = 20.0
      IF (IPN.GE.2) XFAC = 2.0*XFAC
      SX = MIN(5.0*SEMI,RGRAV)
*       Delay termination test for quadruple (NEW HIARCH may be possible).
      IF (N.EQ.4.AND.GPERT.LT.1.0D-06) XFAC = 40.0
*       Avoid termination near pericentre of dominant binary (accuracy loss).
      IF (N.EQ.3.AND.IPN.LE.1) THEN
          RB = MIN(1.0/RINV(1),1.0/RINV(2))
*       Note RINV is not quite latest value (differs from X1-X2 by 2-3 %).
          IF (RB.LT.0.5*SEMI) XFAC = 100.0
      END IF
*
*       Check termination (IESC, JESC) for different cases (IPN & GPERT).
      IF ((IPN.LE.1.AND.RX.GT.XFAC*SX.AND.GPERT.LT.1.0D-03).OR.
     &    (GPERT.GT.1.0D-03.AND.RX.GT.2.0*SX).OR.
     &    (IPN.GE.2.AND.RX.GT.XFAC*SX)) THEN
*
          CALL CHMOD(ISUB,KCASE,IESC,JESC)    ! Controls decision-making.
          IF (IESC.GT.0.AND.N.GT.3) THEN
              CALL REDUCE(IESC,JESC,ISUB)
              IESC = 0
              JESC = 0
              GO TO 30
          END IF
      END IF
*
*       Include termination condition for N = 2 or large perturbation.
      IF ((N.EQ.2.AND.DW.LT.1.0D-03.AND.GPERT.LT.1.0D-07).OR.
     &    (N.GT.2.AND.GPERT.GT.1.0)) THEN
          WRITE (6,190)  NSTEP1, IPN, DW, GPERT
  190     FORMAT (' QUIT CHAIN    # IPN DW GP  ',I6,I4,1P,2E10.2)
          GO TO 250 
      END IF
*
*       Decide between increased membership, escape removal or termination.
      IF (KCASE.EQ.0.OR.ISUB.EQ.0) GO TO 400
      IF (KCASE.EQ.1) THEN
          IF (IESC.EQ.0) GO TO 30
          ITERM = -1
          GO TO 258
      END IF
*
*       Note KCASE = 0 for standard return and = -1 and -1 for binary escape.
      IF (KCASE.EQ.-1) THEN
          CALL CONST(X,V,M,N,ENER0,G0,AL)
          ERR = (ENER0 - ENERGY)/ENER0
          IF (ABS(ERR).GT.1.0D-04) THEN
              WRITE (6,192)  N, KCASE, TNOW, NSTEP1, ERR
  192         FORMAT (' CHAIN CHANGE    N KCASE T # ERR ',
     &                                  2I4,F10.3,I9,1P,E11.2)
              WRITE (6,195)  (1.0/RINV(K),K=1,N-1)
  195         FORMAT (' DISTANCES   ',1P,5E10.2)
          END IF
          GO TO 250
      END IF
*
*       Note KCASE = -1 for binary escape with N = 4.
      IF (KCASE.NE.0) THEN
          NEWREG = .TRUE.
          KCASE = 0
      END IF
*
*       Exit normally.
      GO TO 300
*
*       Set zero step to define termination (just in case).
  250 STEPS(ISUB) = 0.0D0
      WRITE (6,255)  TNOW, NSTEP1, N, NBH2, ECC, SEMI, RX, ECH,
     &               GPERT
  255 FORMAT (' END CHAIN    T # N NBH ECC SEMI RX ECH G ',
     &                       F10.4,I9,2I4,F9.5,1P,4E10.2)
*       Distinguish between GR coalescence and standard termination.
      IF (ICOAL.EQ.0) NBH2 = 0
      ITERM = -1
*
*       Treat different termination cases separately: N = 1, 2 or > 2.
  258 IF (N.LE.2) THEN
*       Re-initialize N=1 or N=2 (KS energy contains c.m. & injected body).
          NBH2 = 0
          IF (ISTAR(1).EQ.14) NBH2 = 1
          IF (N.EQ.1) THEN
              CALL CHTERM(NBH2)
          ELSE
              IF (ICOAL.GT.0) NBH2 = NBH2 + 1
              CALL CHTERM2(NBH2)
          END IF
          GO TO 290
      ELSE
*
*       Determine most distant member for removal (binary may be central).
          RX = 0.0
          LK = 0
          DO 275 L = 1,N
              RI2 = 0.0
              DO 270 K = 1,3
                  RI2 = RI2 + X(K+LK)**2
  270         CONTINUE
              IF (RI2.GT.RX) THEN
                  RX = RI2
                  LX = L
              END IF
              LK = LK + 3
  275     CONTINUE
          RX = SQRT(RX)
*
*       Adopt 5*RGRAV as criterion to identify IESC since RMIN not available.
          IF (IESC+JESC.EQ.0.AND.
     &       (RX.GT.5.0*RGRAV.OR.ITERM.LT.0)) IESC = LX
*
*       Delay termination for small external perturbation (< 10000 steps).
          IF (IESC.GT.0.AND.GPERT.LT.1.0D-04) THEN
              XF = 10.0
              IF (GPERT.LT.1.0D-07) XF = 20.0
*       Note that stability condition requires termination (ITERM < 0).
              IF (RX.LT.XF*RGRAV.AND.NSTEP1.LT.9999.AND.ITERM.EQ.0) THEN
                  IESC = 0
                  JESC = 0
              END IF
          END IF
*
*       Check single versus binary escape in four-body case.
          IF (IESC.GT.0.AND.N.EQ.4.AND.R2.GT.0.5*RSUM) THEN
              R1 = 1.0/RINV(1)
              R2 = 1.0/RINV(2)
              R3 = 1.0/RINV(3)
*       Note four-body stability test still missing.
              IF (IESC.EQ.1.AND.R1.LT.3.0*R2) THEN
                  JESC = 2
              ELSE IF (IESC.EQ.4.AND.R3.LT.3.0*R2) THEN
                  IESC = 3
                  JESC = 4
              END IF
          END IF
*
*       Remove body #IESC > 0 using the standard procedure (repeat for N > 2).
          IF (IESC.GT.0) CALL REDUCE(IESC,JESC,ISUB)
          IF (N.GE.2) THEN
	      IF (N.EQ.2.AND.IPN.GE.2) GO TO 30
              IF (N.GE.3) THEN
                  IF (IESC.GT.0) THEN
                      IF (STEPS(ISUB).EQ.0.0D0) STEPS(ISUB) = 0.5*TSMIN
                      GO TO 30   ! Note looping on STEPS = 0.
                  END IF
                  IESC = 0
                  JESC = 0
                  GO TO 400
              END IF
              IF (NSTEP1.LT.10000) GO TO 250
          END IF
*       Terminate chain for two last members and exit after setting IGR = 0.
          CALL CHTERM2(NBH2)   ! may not be reached but no harm if N=2.
      END IF
  290 ISUB = -1
      IGR = 0
      GO TO 400
*
*       Save current global time for next CALL CHAIN.
  300 IF (ITERM.GE.0) TS(ISUB) = T0S(ISUB) + TIMEC
      ISUB = ITERM
*       Include GR energy for employing standard ADJUST (ETOT + ECH).
      ECH = ENERGY
*       Note that explicit energy check gives DE/E ~ 1D-10.
      CALL CONST(X,V,M,N,ENER1,G0,AL)
      SUM = 0.0
      DO 305 I = 1,N-1
          DO 304 L = I+1,N
              SUM = SUM + M(I)*M(L)
  304     CONTINUE
  305 CONTINUE
      RGRAV = SUM/ABS(ENER1)
      ERR = (ENERGY - ENER1)/ENERGY
      IF ((IPN.EQ.1.AND.TZ.GT.50.0.AND.MOD(NSTEP1,1000).EQ.0).OR.
     &    (IPN.EQ.1.AND.TZ.LT.50.0.AND.MOD(NSTEP1,100).EQ.0).OR.
     &    (IPN.GT.1.AND.MOD(NSTEP1,25).EQ.0.AND.TZ.LT.0.1).OR.
     &    (IPN.GT.2.AND.MOD(NSTEP1,10).EQ.0.AND.TZ.LT.0.01)) THEN
          ZN = SQRT(MX/SEMI**3)
          PDOT = 3.0*ZN/(1.0 - ECC2)*MX/(SEMI*CLIGHT**2)
          TPOM = 6.283/PDOT
*       Produce stability statistics for N = 2 (including EB & RSUB).
          ESUB = 0.0
          RSUB = 0.0
          IF (N.EQ.2) CALL BHSTAB(ESUB,RSUB,ISTAB)
          DECC = ECC - ECC0
          WRITE (6,306)  IPN, N, NPERT, TNOW, ECC, SEMI, TZ,
     &                   RGRAV, DECC
  306     FORMAT (' WATCH!   IPN N NP T E A TZ RG DECC ',
     &                       3I3,F11.4,F8.4,F10.6,1P,E10.2,E9.1,E10.2)
*       Perform extra disruption check.
          ITRY = 0
          IF (SEMI*(1.0 - ECC).LT.RCOLL) ITRY = 1
*       Include swallowing condition for close weakly perturbed WD binary.
          IF (N.EQ.2.AND.GPERT.LT.1.0D-08) THEN
              IF (MIN(ISTAR(1),ISTAR(2)).GE.10.AND.
     &            MAX(ISTAR(1),ISTAR(2)).EQ.14) THEN
                  IF (TZ.LT.200.0) THEN
                      ITRY = 1
                  END IF
              ELSE
*       Switch to unperturbed PN treatment unless GR time-scale is small.
                  IF (IPN.LT.2) GO TO 258
              END IF
          END IF
          IF (ITRY.GT.0) THEN
              IF (M(1).GE.M(2)) THEN
                  IBH = 1
              ELSE
                  IBH = 2
              END IF
              IESC = 3 - IBH
              NBH2 = 0
              CALL INFALL(IBH,IESC,NBH2,ISUB)
              IF (N.EQ.1) THEN
*       Correct for ECH in ECOLL before setting to zero and terminate.
                  CALL DECORR(ECH)
                  CALL CHTERM(NBH2)
              END IF
              ISUB = -1
              IGR = 0
              GO TO 400
          END IF
*       Refresh non-zero perturber list during significant shrinkage.
          IF (IPN.GE.2.AND.NPERT.GT.0) THEN
              KDUM = 1
              CALL CHLIST(KDUM)
          END IF
      ELSE IF (IPN.EQ.0.AND.MOD(NSTEP1,10000).EQ.0) THEN
          WRITE (23,308)  N, NPERT, TNOW, ECC, SEMI, GPERT, TZ
  308     FORMAT (' CHECK    N NP T E A G TZ ',2I3,F11.4,F8.4,1P,3E10.2)
          CALL FLUSH(23)
      END IF
***   IF (N.EQ.2) GO TO 250  ! tested OK
      IF (ICOLLISION.GT.0) THEN
          WRITE (6,310)  IPN, I1, I2, ECC0, ECC, SEMI, RX, RZ
  310     FORMAT (' MISSED COLLISION    IPN I1 I2 E0 E A RX RZ ',
     &                                  3I4,2F8.4,1P,3E10.2)
      END IF
*
*       Update time and energy (STEPS(ISUB) may change by TSMIN on entry).
  400 IF (ISUB.GT.0) TS(ISUB) = T0S(ISUB) + TIMEC
      ECH = ENERGY - ESUM
      IF (N.EQ.0) CALL DECORR(-ESUM)
*
      RETURN
*
      END
