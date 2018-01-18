      SUBROUTINE CHMOD(ISUB,KCASE,IESC,JESC)
*
*
*       Modification of chain member(s).
*       --------------------------------
*
      INCLUDE 'common6.h'
        REAL*8  M,MASS,MC
        PARAMETER (NMX=10,NMX3=3*NMX,NMX4=4*NMX,NMXm=NMX*(NMX-1)/2)
         COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &   XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
      INTEGER  ISORT(NMX)
      REAL*8  XCM(3),VCM(3),XJ(3),VJ(3)
      SAVE IT,IWARN,NAMESC,NAMESC2,NAMESC3,NEXT
      DATA IT,IWARN,NAMESC,NAMESC2,NAMESC3,NEXT /6*0/
*
*
      IESC = 0
      JESC = 0
      RI = 0.0
*       Copy chain variables to standard form.
      LK = 0
      DO 3 L = 1,NCH
          DO 2 K = 1,3
              LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
    2     CONTINUE
    3 CONTINUE
*
*       Identify the dominant perturber (skip if none, check on NN >= 7).
      ITRY = 0
      JCLOSE = 0
      NNB = LISTC(1)
      IF (NNB.EQ.0) GO TO 10
      PMAX = 0.0
      DO 5 L = 2,NNB+1
          J = LISTC(L)
          RIJ2 = (X(1,J) - X(1,ICH))**2 + (X(2,J) - X(2,ICH))**2 +
     &                                    (X(3,J) - X(3,ICH))**2
          PIJ = BODY(J)/(RIJ2*SQRT(RIJ2))
          IF (PIJ.GT.PMAX) THEN
              PMAX = PIJ
              RJMIN2 = RIJ2
              JCLOSE = J
          END IF
    5 CONTINUE
      RIJ = SQRT(RJMIN2)
*
*       Restrict membership increase above 6 for distant perturber.
      IF (NCH.GE.6.AND.RIJ.GT.30.0*RMIN) THEN
          GO TO 10
      END IF
*
*       Check enforced escape removal to allow a closer perturber.
      IF ((GPERT.GT.0.001.AND.NN.GE.6).OR.
     &    (GPERT.GT.0.01.AND.NN.GE.5)) THEN
          RX2 = 0.0
          IT = 0
*       Determine most distant outgoing member (if any).
          DO 105 I = 1,NN
              RIJ2 = 0.0
              RD = 0.0
              DO 102 K = 1,3
                  RIJ2 = RIJ2 + X4(K,I)**2
                  RD = RD + X4(K,I)*XDOT4(K,I)
  102         CONTINUE
              IF (RIJ2.GT.RX2.AND.RD.GT.0.0) THEN
                  IT = I
                  RX2 = RIJ2
                  RDOT = RD
              END IF
  105     CONTINUE
*
*       Include safety skip on no identification.
          IF (IT.EQ.0) THEN
              GO TO 10
          END IF
*
*       Find distance to nearest member for safety check (> 2*RMIN).
          RP2 = 1.0
          DO 110 I = 1,NN
              IF (I.EQ.IT) GO TO 110
              RIJ2 = 0.0
              DO 108 K = 1,3
                  RIJ2 = RIJ2 + (X4(K,I) - X4(K,IT))**2
  108         CONTINUE
              IF (RIJ2.LT.RP2) THEN
                  RP2 = RIJ2
              END IF
  110     CONTINUE
*
*       Check central distance and closest separation.
          IF (RX2.GT.RJMIN2.AND.RP2.GT.4.0*RMIN**2) THEN
              RI = SQRT(RX2)
*       Skip if escaper would increase the perturbation.
*             GX = 2.0*BODYC(IT)/(MASS - BODYC(IT))*((RSUM-RI)/RI)**3
              CALL CPERTJ(IT,GX)
              IF (GX.GT.MAX(GPERT,1.0D-05).AND.RI.LT.30.0*RMIN.AND.
     &            RSUM.LT.50.0*RMIN) THEN
                  KCASE = 0
                  GO TO 60
              END IF
              RDOT = RDOT/RI
              HI = 0.5*RDOT**2 - BODY(ICH)/RI
*             WRITE (6,112)  NAMEC(IT), NAME(JCLOSE), GPERT, RI,
*    &                       RDOT, RIJ, SQRT(RP2)
* 112         FORMAT (' ENFORCED ESCAPE    NAMC NAMJ GP RI RD RJ RP ',
*    &                                     2I7,1P,5E10.2)
              IESC = IT
              JESC = 0
              KCASE = 1
              VINF = 0.0
              GO TO 40
          END IF
      END IF
      IF (NN.GE.9) GO TO 10
*       Prevent absorbed particle escaping from CHAIN/REDUCE
      IF (RIJ.GT.2.0*RMIN) GO TO 10
*
*       Form the scalar product R*V for sign of radial velocity.
      RDOT = (X(1,JCLOSE) - X(1,ICH))*(XDOT(1,JCLOSE) - XDOT(1,ICH)) +
     &       (X(2,JCLOSE) - X(2,ICH))*(XDOT(2,JCLOSE) - XDOT(2,ICH)) +
     &       (X(3,JCLOSE) - X(3,ICH))*(XDOT(3,JCLOSE) - XDOT(3,ICH))
*
      IF (GPERT.GT.1.0.AND.IWARN.LT.10) THEN
          IWARN = IWARN + 1
          WRITE (6,115)  NN, NAME(JCLOSE), NPERT, GPERT, RGRAV, RIJ,
     &                   RDOT/RIJ
  115     FORMAT (' WARNING!   N NM NP G RG RIJ RD ',
     &                         I4,I6,I4,F7.2,1P,4E9.1)
      END IF
      IF (NN.GT.3.AND.GPERT.GT.1.0) IWARN = 0
*
*       Skip outgoing orbit with GPERT < 1D-04 or NCH > 4 & GPERT < 0.2.
      IF (RDOT.GT.0.0.AND.(NCH.GE.5.OR.GPERT.LT.0.2)) GO TO 10
      IF (RDOT.GT.0.0.AND.GPERT.LT.1.0D-04) GO TO 10
      IF (RDOT.GT.0.0) THEN
          HI = 0.5*(RDOT/RIJ)**2 - BODY(ICH)/RIJ
          RX = -BODY(ICH)/HI
          IF (RX.GT.10.0*RMIN.OR.HI.GT.0.0) GO TO 10
          IF (GPERT.LT.0.1) GO TO 10
      END IF
*       Quit on RIJ > 5*RMIN (or NCH = 2 & RDOT > 0) to prevent switching.
      IF (NCH.EQ.2.AND.(RIJ.GT.5.0*RMIN.OR.RDOT.GT.0.0)) GO TO 10
      IF (NN.GE.7.AND.JCLOSE.GT.N) GO TO 10
*
*       Note looping for NCH = 2 which gives GAMX = 0.
      IF (NCH.LE.5.AND.GPERT.LT.0.2.AND.NNB.GT.1.AND.
     &    RIJ.GT.RSUM) THEN
          CALL CPERTX(JCLOSE,3,GAMX,GAMY)
          IF (NCH.GT.2.AND.GAMY.LT.0.1.AND.GAMX.LT.GAMY) GO TO 10
      END IF
*
*       Include extra procedure for closest perturber.
      DO 120 K = 1,3
          XJ(K) = X(K,JCLOSE) - X(K,ICH)
  120 CONTINUE
      RX2 = 1.0
      DO 130 L = 1,NCH
          RIJ2 = 0.0
          DO 125 K = 1,3
              RIJ2 = RIJ2 + (X4(K,L) - XJ(K))**2
  125     CONTINUE
          IF (RIJ2.LT.RX2) THEN
              RX2 = RIJ2
          END IF
  130 CONTINUE
*
*       Check pericentre for three-body case with massive intruder.
      IF (NCH.EQ.2.AND.KSTAR(JCLOSE).EQ.14) THEN
          RD = 0.0
          RIJ2 = 0.0
          VI2 = 0.0
          DO 135 K = 1,3
              RIJ2 = RIJ2 + XJ(K)**2
              VJ(K) =  XDOT(K,JCLOSE) - XDOT(K,ICH)
              VI2 = VI2 + VJ(K)**2
              RD = RD + XJ(K)*VJ(K)
  135     CONTINUE
          RJ = SQRT(RIJ2)
          SEMI = 2.0/RJ - VI2/(BODY(ICH) + BODY(JCLOSE))
          SEMI = 1.0/SEMI
          ECC2 = (1.0 - RJ/SEMI)**2 +
     &                  RD**2/(SEMI*(BODY(ICH) + BODY(JCLOSE)))
          ECC = SQRT(ECC2)
          PMIN = SEMI*(1.0 - ECC)
          IF (PMIN.GT.5.0*RGRAV.AND.GPERT.LT.0.001.AND.NPERT.GT.4) THEN
              IF (MOD(NSTEP1,1000).EQ.0) THEN
              WRITE (6,140)  NSTEP1, NAME(JCLOSE), ECC, SEMI,
     &                       SEMI*(1.0-ECC), 5.0*RGRAV, GPERT
  140         FORMAT (' SKIP HEAVY   # NM E A PM 5*RG GP ',
     &                               I8,I4,F8.4,1P,4E10.2)
              CALL FLUSH(6)
              END IF
              GO TO 10
          END IF
      END IF
*
*       Delay acceptance of approaching third light body.
      IF (NCH.EQ.2.AND.KSTAR(JCLOSE).LT.14.AND.RDOT.LT.0.0) THEN
          IF (RIJ.GT.5.0*RGRAV.AND.GPERT.LT.0.0001) GO TO 10
      END IF
*
*       Evaluate relative perturbation on #LL due to #JCLOSE.
*     RI2 = X4(1,LL)**2 + X4(2,LL)**2 + X4(3,LL)**2
*     FF = BODY(ICH)/RI2
*     GP = BODY(JCLOSE)/(RX2*FF)
*
*       Adopt 4*RMIN criterion or combined larger distance & perturbation.
      RX2 = MIN(RX2,RJMIN2)
*       Exclude large contrast to prevent repeated absorb/escape.
      IF (RIJ.GT.10.0*RSUM) GO TO 10
      RABS = 2.0*RMIN
      GXX = 0.1
      IF (JCLOSE.GT.N) GXX = 0.03
      GSTAR = 0.001*BODY(JCLOSE)/BODY(ICH)
      IF ((RX2.LT.RABS**2.AND.NCH.LE.6).OR.
     &   (GPERT.GT.GSTAR.AND.RX2.LT.RABS**2.AND.NCH.LE.5).OR.
     &   (GPERT.GT.GXX.AND.NCH.LE.7).OR.(KSTAR(JCLOSE).EQ.14)) THEN
*
*       Include radial velocity condition for allowing membership increase.
          RD = RDOT/RIJ
          IF (N.GE.5.AND.RD**2.LT.0.25.AND.GPERT.LT.0.001) THEN
              GO TO 10
          END IF
*
*       Perform safety check.
          DO 145 L = 1,NCH
              IF (NAMEC(L).EQ.NAME(JCLOSE)) THEN
              WRITE (6,144)  NAME(JCLOSE)
  144         FORMAT (' DANGER!   DUPLICATED   NMJ  ',I6)
              STOP
              END IF
  145     CONTINUE
*
*       Treat extreme mixed or classical binary as inert (but not BH-BH).
          IF (JCLOSE.GT.N) THEN
              KX = MIN(KSTAR(2*(JCLOSE-N)-1),KSTAR(2*(JCLOSE-N)))
              SEMI = -0.5*BODY(JCLOSE)/H(JCLOSE-N)
              RM = 1.0
              DO 150 K = 1,NCH-1
                  RM = MIN(1.0/RINV(K),RM)
  150         CONTINUE
              SMALL = 0.001*(RSUM - RIJ)
              IF (RM.LT.SMALL.AND.SEMI.LT.SMALL.AND.KX.LT.14) THEN
                  WRITE (6,155)  NAME(JCLOSE), GPERT, SEMI, RIJ,
     &                                         (1.0/RINV(K),K=1,NCH-1)
  155             FORMAT (' ABSORB INERT    NM G A RIJ R ',
     &                                      I6,F6.2,1P,6E10.2)
*       Reduce option #26 and set #27 = 1 temporarily for routine SETSYS.
                  KZ26 = KZ(26)
                  KZ(26) = 1
                  KZ(27) = 1
                  CALL ABSORB(ISUB)
                  KZ(26) = KZ26
                  KZ(27) = KZ27
                  KCASE = 1
                  GO TO 50
              END IF
          END IF
*
          IF (JCLOSE.GT.N) THEN
              WRITE (6,160)  NAME(JCLOSE), R(JCLOSE-N), GAMMA(JCLOSE-N)
  160         FORMAT (' ABSORB BINARY    NM R G ',I6,1P,2E10.2)
          END IF
          IF (KSTAR(JCLOSE).EQ.14.AND.NCH.EQ.2) THEN
              WRITE (6,165)  NSTEP1, ECC, SEMI, SEMI*(1-ECC), GPERT, RIJ
  165         FORMAT (' ABSORB!    # E A PM GP RIJ ',I8,F8.4,1P,4E10.2)
          END IF
*
          IF (KZ(30).GT.1) THEN
              HI = 0.5*RD**2 - BODY(ICH)/RIJ
              WRITE (6,170) TIME+TOFF, NCH, NAME(JCLOSE), GPERT, RIJ,
     &                      RD, HI
  170         FORMAT (' ABSORB    T NCH NM GP RIJ RD HI ',
     &                            F9.3,I4,I7,1P,4E9.1)
          END IF
*
*       Absorb the perturber (single particle or binary).
          IF (KZ(30).GT.2) CALL ABSORB(ISUB)
*
*       Reduce time since new c.m. step may be very small.
          TIME = MIN(TIME,TBLOCK)
*
*       Activate indicator for new chain treatment.
          KCASE = 1
          ITRY = ITRY + 1
*         GO TO 10
          RETURN
      END IF
*
*       Check for rejection (RIJ > 2*MIN(RSUM,RMIN); RDOT > 0 & G < 0.05).
      IF (GPERT.LT.0.02*(BODY(JCLOSE)/BODY(ICH)).OR.
     &    NAME(JCLOSE).LT.0) GO TO 10
      IF (RDOT.GT.0.0.AND.GPERT.LT.0.05*(BODY(JCLOSE)/BODY(ICH))) THEN
          GO TO 10
      END IF
      IF (RSUM.GT.10.0*RMIN.AND.
     &    GPERT.LT.0.01*(BODY(JCLOSE)/BODY(ICH))) GO TO 10
      IF (RIJ.GT.2.0*RSUM.OR.RSUM.GT.10.0*RMIN) GO TO 10
*
*       Perform impact parameter test (abandoned 19/7/07).
      VR2 = (XDOT(1,ICH) - XDOT(1,JCLOSE))**2 +
     &      (XDOT(2,ICH) - XDOT(2,JCLOSE))**2 +
     &      (XDOT(3,ICH) - XDOT(3,JCLOSE))**2
      AINV = 2.0/RIJ - VR2/(MASS + BODY(JCLOSE))
      ECC2 = (1.0 - RIJ*AINV)**2 + RDOT**2*AINV/(MASS + BODY(JCLOSE))
      ECC = SQRT(ECC2)
      IF (NN.EQ.2.AND.GPERT.LT.0.02*(BODY(JCLOSE)/BODY(ICH))) THEN
          GO TO 10
      END IF
*
*       Include additional criteria (V^2 > VC^2/2; JCL > N & RIJ < RSUM).
      RDOT = RDOT/RIJ
      VC2 = (MASS + BODY(JCLOSE))/RIJ
      IF ((RSUM + RIJ.LT.RMIN).OR.(NN.EQ.2.AND.GPERT.GT.0.005).OR.
     &    (RIJ.LT.RSUM.AND.RDOT**2.GT.0.5*VC2).OR.
     &    (JCLOSE.GT.N.AND.RIJ.LT.RSUM).OR.RDOT.LT.0.0) THEN
*
*       Do not allow multiple absorptions (NAMES(NMX,5) = 0 each new chain).
          IF (NAME(JCLOSE).NE.NAMES(NMX,5).AND.NAMES(NMX,5).EQ.0) THEN
              NAMES(NMX,5) = NAME(JCLOSE)
*       Include possible return after ejection in bound orbit.
          ELSE IF (RDOT.GT.0.0) THEN
              GO TO 10
          END IF
*
          RSUM = RSUM + RIJ
          IF (KZ(30).GT.1) THEN
              SEMI1 = 1.0/AINV
              PMIN = SEMI1*(1.0 - ECC)
              WRITE (6,6)  NSTEPI, NAME(JCLOSE), GPERT, ECC, SEMI1,
     &                     RIJ, RSUM, PMIN
    6         FORMAT (' ABSORB:    # NMJ GPERT E A RIJ RSUM PMIN ',
     &                             I11,I7,F8.4,F9.5,1P,5E9.1)
              CALL FLUSH(6)
          END IF
*
          IF (JCLOSE.GT.N) THEN
              WRITE (6,155)  NAME(JCLOSE), R(JCLOSE-N), GAMMA(JCLOSE-N)
          END IF
*
*       Absorb the perturber (single particle or binary but note INJECT).
          IF (KZ(30).GT.2) CALL ABSORB(ISUB)
*
*       Reduce time since new c.m. step may be very small.
          TIME = MIN(TIME,TBLOCK)
*
*       Activate indicator for new chain treatment and try a second search.
          KCASE = 1
          ITRY = ITRY + 1
*         IF (ITRY.EQ.1) GO TO 1
      END IF
*
*       Exit if any particles have been absorbed.
   10 IF (ITRY.GT.0) GO TO 50
*
      KCASE = 0
      JESC = 0
      IF (NN.EQ.2) GO TO 60
*     IF (NN.EQ.2) THEN
*         R1 = 1.0/RINV(1)
*         IF (R1.GT.100.0*RMIN) THEN
*             KCASE = -1
*             GO TO 60
*         END IF
*     END IF
*
*       Place index of the smallest INVERSE distance in ISORT(1).
      CALL HPSORT(NN-1,RINV,ISORT)
*
*       Determine index of escaper candidate (single star or binary).
*       Distinguish two cases of each type (beginning & end of chain).
      IBIN = 1     ! Needs to be defined for KCASE = 2.
      IF (ISORT(1).EQ.1) THEN
          IESC = INAME(1)
          KCASE = 1
      ELSE IF (ISORT(1).EQ.NN-1) THEN
          IESC = INAME(NN)
          KCASE = 1
*       Check for possible binary escaper (NN = 3 implies single escaper).
      ELSE IF (ISORT(1).EQ.2) THEN
          IESC = INAME(1)
          JESC = INAME(2)
          IBIN = 1
*       Switch binary index in case last separation is smallest (NN >= 4).
          IF (RINV(1).LT.RINV(NN-1)) THEN
              IBIN = NN - 1
              IF (BODYC(IESC) + BODYC(JESC).GT.0.5*MASS) THEN
                  IBIN = NN - 1
              END IF
          END IF
          KCASE = 2
      ELSE IF (ISORT(1).EQ.NN-2) THEN
          IESC = INAME(NN-1)
          JESC = INAME(NN)
          IBIN = NN - 1
          KCASE = 2
      END IF
*
*       Skip if biggest chain distance is not dominant.
      IF (KCASE.EQ.1) THEN
          IS = ISORT(1)
          IF (1.0/RINV(IS).LT.0.8*RSUM) THEN
              KCASE = 0
              IESC = 0
              GO TO 60
          END IF
      END IF
*       Ensure massive binary is not removed (big single should be OK).
      IF (KCASE.EQ.2.AND.NN.GT.4) THEN   ! Avoid problems with NN >  4.
          IF (IESC.EQ.INAME(1)) THEN
              J1 = INAME(NN-1)
              J2 = INAME(NN)
          ELSE
              J1 = INAME(1)
              J2 = INAME(2)
          END IF
          ZMB = BODYC(IESC) + BODYC(JESC)
          ZMB2 = BODYC(J1) + BODYC(J2)
*       Retain the most massive binary for AR_Chain case.
          IF (ZMB.GT.10.0*ZMB2) THEN       ! Factor 10 included 1/16.
              IESC = J1
              JESC = J2
              IBIN = NN - IBIN
          END IF
      END IF
*
*       Exclude difficult cases (like small first binary).
      IF (KCASE.EQ.2.AND.IBIN.EQ.1.AND.1.0/RINV(NN-1).GT.2.0*RMIN) THEN
*       Switch to removal of last (single) member (R2>20*RMIN may take long).
          IESC = INAME(NN)
          JESC = 0
*       Skip case of last distance exceeding RMIN.
      ELSE IF (KCASE.EQ.2.AND.NCH.EQ.4.AND.1.0/RINV(IBIN).GT.RMIN) THEN
          IESC = 0
          JESC = 0
          KCASE = 0
          RI = 1.0/RINV(2)
          GO TO 60
      END IF
*
*       Include escape check if middle distance is largest.
      IF (KCASE.EQ.0.AND.NN.GE.5.AND.
     &    1.0/RINV(ISORT(1)).GT.2.0*RMIN) THEN
          R1 = 1.0/RINV(1)
          RN = 1.0/RINV(NN-1)
*       Set relevant indices for beginning or end of chain.
          IF (R1.GT.RN) THEN
              IESC = INAME(1)
              JX = INAME(2)
              IB = 1
              ISORT(1) = 1
              R2 = 1.0/RINV(2)
          ELSE
              IESC = INAME(NN)
              JX = INAME(NN-1)
              IB = NN - 1
              ISORT(1) = NN - 1
              R2 = 1.0/RINV(NN-2)
          END IF
*       Define binary indices for large second separation.
          IF (R2.GT.MAX(R1,RN)) THEN
              JESC = JX
              IBIN = IB
              KCASE = 2
          ELSE
              KCASE = 1
          END IF
      END IF
*
*       Skip on JESC > 0 not at beginning/end but check large chain.
      IF (JESC.GT.0.AND.(JESC.NE.1.AND.JESC.NE.NCH)) KCASE = 0
      IF (JESC.GT.0.AND.RSUM.GT.10.0*RMIN) KCASE = 2
*
*       Skip if no identified case (should be rare).
      IF (KCASE.EQ.0) GO TO 60
*
*       Ensure enforced escape of binary in wide four-body system.
      IF (KCASE.EQ.2.AND.NCH.EQ.4.AND.1.0/RINV(2).GT.20.0*RMIN) THEN
          IF (NAMEC(IESC).NE.NAMESC3.AND.JESC.GT.0) THEN    ! JESC bug 2/17.
              IF (NSTEP1.GT.NEXT) THEN
              WRITE (6,11)  NSTEP1, IESC, JESC, NPERT, NAMEC(IESC),
     &                      NAMEC(JESC), 1.0/RINV(IBIN), 1.0/RINV(2)
   11         FORMAT (' ENFORCED ESCAPE    # IESC JESC NP NAM RB R2 ',
     &                                     I6,3I4,2I7,1P,2E10.2)
              NAMESC3 = NAMEC(IESC)
              NEXT = NSTEP1 + 50
              END IF
          END IF
          RI = 1.0/RINV(2)
*       Mark the smallest (last) binary for removal even if IESC/JESC = 1/2.
          IF (IESC+JESC.EQ.3.AND.1.0/RINV(3).LT.1.0/RINV(1)) THEN
              IESC = 3
              JESC = 4
          END IF
          KCASE = -3
          GO TO  60    ! attempt Dec 2015.
      END IF
      IF (KCASE.EQ.2.AND.NCH.EQ.4) THEN
          RI = 0.0
          KCASE = 0
          GO TO 60
      END IF
*
*       Enforce escape for distant member with positive radial velocity.
      RD = X4(1,IESC)*XDOT4(1,IESC) + X4(2,IESC)*XDOT4(2,IESC)
     &                              + X4(3,IESC)*XDOT4(3,IESC)
      IF (NCH.EQ.4.AND.RSUM.GT.20.0*RMIN.AND.RD.GT.0.0D0) THEN
          KCASE = 1
          RI = 20.0*RMIN
          VINF = 0.0
          IF (NSTEP1.GT.NEXT) THEN
      WRITE (6,41)  NSTEP1, IESC, JESC, GPERT, (1.0/RINV(K),K=1,NCH-1)
   41 FORMAT (' WIDE ESCAPE    # IESC JESC G R  ',I8,2I4,F8.4,1P,4E10.2)
          NEXT = NSTEP1 + 20
          END IF
          GO TO 40
      END IF
*
      IF (KZ(30).GT.3) THEN
          WRITE (6,12)  IESC, JESC, NSTEP1, ISORT(1),
     &                  (1.0/RINV(K),K=1,NN-1)
   12     FORMAT (' CHMOD:    IESC JESC # ISORT1 R ',2I3,I5,I3,1P,5E9.1)
          CALL FLUSH(6)
      END IF
*
*       Distinguish between binary and single particle.
      RB = 0.0
      IF (JESC.EQ.0) GO TO 20
*
*       Switch to possible non-dominant mass binary.
      ZMB = BODYC(IESC) + BODYC(JESC)
      IF (ZMB.GT.0.9*MASS) THEN         ! Factor 0.5 changed to 0.9 (1/16).
          IF (IESC.EQ.INAME(1)) THEN
              IESC = INAME(NN-1)
              JESC = INAME(NN)
              IBIN = NN - 1
          ELSE
              IESC = INAME(1)
              JESC = INAME(2)
              IBIN = 1
          END IF
      END IF
*
*       First check case of escaping binary (JESC > 0 & RB < RJB/4).
      IF (JESC.GT.0) THEN
          RB = 1.0/RINV(IBIN)
          JBIN = IBIN + 1
          IF (IBIN.EQ.NN - 1) JBIN = IBIN - 1
          RJB = 1.0/RINV(JBIN)
*       Consider removal of outermost particle instead if binary is wide.
          IF (RB.GT.0.25*RJB) THEN
*       Change pointer to end of chain and redefine IESC (with JESC = 0).
              IF (ISORT(1).EQ.2) THEN
                  ISORT(1) = 1
              ELSE
                  ISORT(1) = NN - 1
              END IF
              IESC = JESC
              JESC = 0
              GO TO 20
          END IF
*       Skip if IESC & JESC are not consecutive.
          IF (IABS(IESC-JESC).NE.1) THEN
              KCASE = 0
              GO TO 60
          END IF
*       Ensure removal of widest binary (maybe not as KS) at RSUM > 10*RMIN.
          IF (RSUM.GT.10.0*RMIN.AND.RB.GT.0.25*RJB) THEN
              JBIN = NN - 1
              IF (JBIN.EQ.IBIN) IBIN = 1
              RB = 1.0/RINV(IBIN)
              RJB = 4.01*RB
              WRITE (6,13) IESC, JESC, IBIN, RB, RJB, ZMB*SMU
   13         FORMAT (' IMPOSED ESCAPE    IESC JESC IB RB RJB MB ',
     &                                    3I4,1P,3E10.2)
          END IF
      ELSE
          GO TO 20
      END IF
*
*       Form coordinates & velocities of distant binary (local c.m. frame).
      RB = 1.0/RINV(IBIN)
      BCM = BODYC(IESC) + BODYC(JESC)
      RI2 = 0.0
      RDOT = 0.0
      VREL2 = 0.0
      DO 15 K = 1,3
          XCM(K) = (BODYC(IESC)*X4(K,IESC) + BODYC(JESC)*X4(K,JESC))/BCM
          VCM(K) = (BODYC(IESC)*XDOT4(K,IESC) +
     &              BODYC(JESC)*XDOT4(K,JESC))/BCM
          RI2 = RI2 + XCM(K)**2
          RDOT = RDOT + XCM(K)*VCM(K)
          VREL2 = VREL2 + (XDOT4(K,IESC) - XDOT4(K,JESC))**2
   15 CONTINUE
*
*       Convert to relative distance & radial velocity w.r.t. inner part.
      FAC = BODY(ICH)/(BODY(ICH) - BCM)
      RI = SQRT(RI2)
      RDOT = FAC*RDOT/RI
      RI = FAC*RI
      RJB = 0.0
*       Adopt harmonic mean of RSUM and RMAXS for delaying escape.
*     DESC = 0.5*SQRT(RSUM*RMAXS(ISUB))
      DESC = SQRT(RMIN*RMAXS(ISUB))
*       Reduce delay test for large length contrasts.
      IM = ISORT(1)
      RM = 1.0/RINV(IM)
      CX = RSUM/(RSUM - RM)
      IF (CX.GT.25.0) DESC = 25.0*RGRAV
      RESC = MAX(3.0*RGRAV,DESC)
      AINV = 2.0/RB - VREL2/BCM
      EB = -0.5*BODYC(IESC)*BODYC(JESC)*AINV
      HI = 0.5*RDOT**2 - BODY(ICH)/RI
      IF (HI.GT.0.0) THEN
          VINF = SQRT(2.0*HI)*VSTAR
      ELSE
          VINF = 0.0
      END IF
*
*       Employ parabolic escape criterion (terminate if RI > RMIN & NCH < 5).
      IF ((RI.GT.RESC.AND.RDOT.GT.0.0.AND.RB*AINV.GT.0.99).OR.
     &    (HI.GT.0.0.AND.RI.GT.3.0*RMIN)) THEN
*       Ensure outward motion to prevent repeated acceptance (4/13).
          IF (RDOT.GT.0.0.AND.RDOT**2.LT.2.0*BODY(ICH)/RI) THEN
*       Define effective perturbation using remaining size.
              GB = 2.0*BCM/(BODY(ICH) - BCM)*((RSUM - RJB - RB)/RI)**3
*       Split into two KS solutions if binaries are well separated.
              IF (NCH.EQ.4.AND.GB.LT.0.001) THEN
*                 KCASE = -1
                  KCASE = 1
                  GO TO 40
              END IF
*       Enforce termination if RI > MAX(RMAXS/2,RMIN) and NCH <= 4.
              IF (RI.GT.MAX(0.5*RMAXS(ISUB),RMIN).AND.NCH.LE.4) THEN
*                 KCASE = -1
                  KCASE = 1
                  GO TO 40
*       Accept binary escape for small perturbation & V**2 > M/R (NCH > 4).
              ELSE IF (GB.LT.0.01.AND.NCH.GT.4.AND.
     &                 (RDOT**2.GT.BODY(ICH)/RI.OR.RI.GT.RMIN)) THEN
                  WRITE (6,18)  IESC, JESC, NAMEC(IESC), NAMEC(JESC),
     &                          RI, RDOT**2, 2.0*BODY(ICH)/RI, RB
                  CM(9) = CM(9) - EB
                  KCASE = 1
                  GO TO 40
              ELSE
*       Check splitting into two KS solutions (R1 + R3 < R2/5 & RDOT > VC).
                  IF (NCH.EQ.4) THEN
                      R13 = 1.0/RINV(1) + 1.0/RINV(3)
                      VC2 = RDOT**2*RI/BODY(ICH)
*       Ensure RSUM > 0.1*RMIN for reducing differential force corrections.
                      IF (R13.LT.0.2/RINV(2).AND.VC2.GT.1.0.AND.
     &                    RSUM.GT.0.1*RMIN) THEN
*                         KCASE = -1
                          KCASE = 1
                          GO TO 40
                      END IF
                  END IF
                  KCASE = 0
*       Try single escape for wide binary.
                  IF (1.0/RINV(IBIN).GT.4.0*RMIN) THEN
                      JESC = 0
                      GO TO 20
                  END IF
                  GO TO 60
              END IF
          ELSE
              IF (HI.GT.0.0) THEN
                  VINF = SQRT(2.0*HI)*VSTAR
              ELSE
                  VINF = 0.0
              END IF
              IF (KZ(30).GT.1.OR.VINF.GT.1.0) THEN
                  IF (NAMEC(IESC).NE.NAMESC) THEN
                  WRITE (6,18)  IESC, JESC, NAMEC(IESC), NAMEC(JESC),
     &                          RI, RDOT**2, 2.0*BODY(ICH)/RI, RB, VINF
   18             FORMAT (' CHAIN ESCAPE:    IESC JESC NM RI RDOT2 ',
     &                                      '2*M/R RB VINF ',
     &                                       2I3,2I6,1P,4E9.1,0P,F6.1)
                      NAMESC = NAMEC(IESC)
                  END IF
                  KCASE = -3
                  IF (VINF.GT.1.0) GO TO 60
              END IF
*       Enforce termination (KCASE < 0) if NCH <= 4 (final membership <= 2).
              IF (NCH.LE.4) THEN
                  KCASE = -1
                  KCASE = 2     ! Escape of binary enforced via REDUCE.
                  GO TO 40
              END IF
              CM(9) = CM(9) - EB
              KCASE = 1
              GO TO 40
          END IF
      ELSE
          KCASE = 0
          IF (1.0/RINV(IBIN).GT.4.0*RMIN) THEN
              JESC = 0
              GO TO 20
          END IF
          GO TO 60
      END IF
*
*       Include a basic escape criterion for distant single members.
   20 RX = 0.0
      DO 25 I = 1,NN
          RIJ2 = 0.0
          DO 22 K = 1,3
              RIJ2 = RIJ2 + X4(K,I)**2
   22     CONTINUE
          IF (RIJ2.GT.RX) THEN
              IT = I
              RX = RIJ2
          END IF
   25 CONTINUE
*
*       Employ safety check for heavy single closer to c.m. (NCH = 3).
      IF (NN.EQ.3) THEN
          IF (ISORT(1).EQ.1) THEN
              IT = INAME(1)
          ELSE IF (ISORT(1).EQ.NN-1) THEN
              IT = INAME(NN)
          END IF
      END IF
*
*       Bypass complications for distant member at end of chain.
      IF ((IESC.EQ.1.OR.IESC.EQ.NN).AND.
     &               1.0/RINV(IESC).GT.50.0*RMIN) GO TO 210
      RI = SQRT(RX)
      RDOT = X4(1,IT)*XDOT4(1,IT) + X4(2,IT)*XDOT4(2,IT) +
     &                              X4(3,IT)*XDOT4(3,IT)
*       Include safety termination (6/2014).
      IF (NN.EQ.3.AND.RI.GT.5.0*RMIN) THEN
          KCASE = 1
          HI = 0.5*RDOT**2 - BODY(ICH)/RI
          IF (HI.GT.0.0) GO TO 60
      END IF
      IF (RI.GT.20.0*RMIN.AND.RDOT.GT.0.0) THEN
*       Delay for BH and small perturbation (orbit may turn around).
          IF (ISTAR(IT).EQ.14.AND.GPERT.LT.1.0D-04.AND.
     &        RI.LT.50.0*RMIN) THEN
              KCASE = 0
              GO TO 60
          ELSE IF (IT.EQ.1.OR.IT.EQ.NN) THEN
              RDOT = RDOT/RI
              HI = 0.5*RDOT**2 - BODY(ICH)/RI
              RX = 1.0
              IF (HI.LT.0.0) RX = -BODY(ICH)/HI
              RY = MAX(RI,RX)
              DR = ABS(RI - RX)/RX
              IF ((RY.LT.20.0*RMIN.AND.GPERT.LT.1.0D-03).OR.
     &            (RY.LT.30.0*RMIN.AND.GPERT.LT.1.0D-04).OR.
     &            (DR.LT.0.2.AND.GPERT.LT.1.0D-04)) THEN
                  KCASE = 0
                  GO TO 60
              END IF
              CALL CPERTJ(IT,GX)
              IF (GX.GT.2.0*MAX(GPERT,1.D-05).AND.RI.LT.30.0*RMIN) THEN
                  KCASE = 0
                  GO TO 60
              END IF
*             WRITE (6,27)  IT, NAMEC(IT), GPERT, RI, RDOT, RX
*  27         FORMAT (' BASIC ESCAPE:    IESC NM GP RI RDOT RX ',
*    &                                   I3,I7,1P,4E10.2)
              IESC = IT
              JESC = 0
              KCASE = 1
              GO TO 40
          END IF
      END IF
*
          RP2 = 1.0
          DO 200 L = 2,NNB+1
              J = LISTC(L)
              RIJ2 = 0.0
              DO 195 K = 1,3
*       Include extra procedure for closest perturber.
                  XJ(K) = X4(K,IESC) + X(K,ICH)
                  RIJ2 = RIJ2 + (XJ(K) - X(K,J))**2
  195         CONTINUE
              IF (RIJ2.LT.RP2) RP2 = RIJ2
  200     CONTINUE
          IF (RP2.LT.4.0*RMIN**2) THEN
              KCASE = 0
              GO TO 60
          END IF
*
  210 CONTINUE
      RI = SQRT(X4(1,IESC)**2 + X4(2,IESC)**2 + X4(3,IESC)**2)
      IF (ISTAR(IESC).EQ.14.AND.NCH.LE.6) THEN
*         GX = 2.0*BODYC(IESC)/(MASS - BODYC(IESC))*((RSUM-RI)/RI)**3
          CALL CPERTJ(IESC,GX)
          IF (GX.GT.2.0*MAX(GPERT,1.D-05).AND.RI.LT.30.0*RMIN) THEN
              KCASE = 0
              GO TO 60
          END IF
      END IF
*
*       Form relative distance and radial velocity for single particle.
*       Skip massive body test because of large factor FAC.
      IF (RI.LT.RGRAV) THEN
          KCASE = 0
          GO TO 60
      END IF
      RDOT = X4(1,IESC)*XDOT4(1,IESC) + X4(2,IESC)*XDOT4(2,IESC) +
     &                                  X4(3,IESC)*XDOT4(3,IESC)
      FAC = BODY(ICH)/(BODY(ICH) - BODYC(IESC))
      RDOT = FAC*RDOT/RI
      RI = FAC*RI
      IM = ISORT(1)
*       Ensure that escaper is not close to another body.
      RM = MIN(1.0/RINV(IM),RI)
*
*       Check distance to closest chain member.
      RX2 = 1.0
      DO 30 I = 1,NN
          IF (I.EQ.IESC) GO TO 30
          RIJ2 = 0.0
          DO 28 K = 1,3
              RIJ2 = RIJ2 + (X4(K,I) - X4(K,IESC))**2
   28     CONTINUE
          IF (RIJ2.LT.RX2) THEN
              RX2 = RIJ2
          END IF
   30 CONTINUE
      RABS = MIN(RGRAV,2.0*RMIN)
      IF (RX2.LT.RABS**2.AND.RI.LT.10.0*RMIN) THEN
          KCASE = 0
          GO TO 60
      END IF
*
*       Impose termination of wide BH binary after collision (not ARchain).
*     IF (NN.EQ.2.AND.RI.GT.50.0*RMIN) THEN
*         KCASE = -1
*         GO TO 60
*     END IF
*
*       Skip escape removal if closest body is near pericentre.
*     IX = ISORT(NN-1)
*     IF (1.0/RINV(IX).LT.0.05*RGRAV) THEN
*         KCASE = 0
*         GO TO 60
*     END IF
*
*       Include additional tests for removal (RDOT > 0).
      HI = 0.5*RDOT**2 - BODY(ICH)/RI
      RX = 10.0
      IF (HI.LT.0.0) RX = -BODY(ICH)/HI
      DR = ABS(RI - RX)/RX
      IF ((HI.GT.0.0.AND.RI.GT.10.0*RGRAV.AND.RDOT.GT.0.0).OR.
     &    (RI.GT.20.0*RMIN.AND.RDOT.GT.0.0).OR.
     &    (DR.LT.0.2.AND.GPERT.LT.5.0D-03)) THEN
          IF (ISTAR(IESC).EQ.14.AND.GPERT.GT.0.001) THEN
              VINF = 0.0
              IF (HI.GT.0.0) VINF = SQRT(2.0*HI)*VSTAR
*             WRITE (6,32)  NAMEC(IESC), VINF, GPERT, RI, RX
*  32         FORMAT (' HEAVY ESCAPER   NAM VINF G R RX ',
*    &                                  I4,F6.1,1P,3E10.2)
              KCASE = 1
              GO TO 40
          ELSE
              RX = 10.0
              IF (HI.LT.0.0) RX = -BODY(ICH)/RI
              RY = MAX(RI,RX)
              IF (RY.LT.30.0*RMIN.AND.GPERT.LT.1.0D-04) THEN
                  KCASE = 0
                  GO TO 60
              ELSE
                  KCASE = 1
                  GO TO 40
              END IF
          END IF
      END IF
*
*       Check approximate escape criterion outside 3*RGRAV and DESC.
      DESC = 0.5*SQRT(RSUM*RMAXS(ISUB))
*       Note that harmonic mean tends to delay escape.
      RESC = MAX(5.0*RGRAV,DESC)
      RESC = MIN(RESC,3.0*RMIN)
      IF (NN.EQ.3.AND.RESC.LT.0.001*RMIN) RESC = 5.0*RESC
      IF (RM.GT.RESC.AND.RDOT.GT.0.0) THEN
          IF (RDOT**2.LT.2.0*BODY(ICH)/RI) THEN
              FAC2 = 2.0*BODYC(IESC)/(BODY(ICH) - BODYC(IESC))
              GI = FAC2*((RSUM - 1.0/RINV(IM))/RI)**3
*       Do not permit large length contrast even for bound orbit.
              IF (NN.GE.4.AND.RI.GT.5.0*RMIN.AND.GPERT.GT.0.001) THEN
                  KCASE = 1
                  GO TO 35
              END IF
              ER = 0.5*RDOT**2 - BODY(ICH)/RI
              RX = -BODY(ICH)/ER
              IF (ER.LT.0.0.AND.RX.LT.MIN(2.0*RSUM,2.0*RMIN)) THEN
                  KCASE = 0
                  GO TO 60
              END IF
              IF (RI.GT.MAX(0.5*RMAXS(ISUB),3.0*RMIN)) THEN
*       Decide between termination or removal for sub-parabolic motion.
                  IF (NN.LT.4) THEN
*       Delay termination for eccentric binary inside R = -m1*m2/2*E < a.
                      IB = 3 - IM
                      J1 = INAME(IB)
                      J2 = INAME(IB+1)
*       Obtain semi-major axis from ENERGY = -(m1*m2 + (m1+m2)*m3)/RGRAV.
                      ZMU = BODYC(J1)*BODYC(J2)/(BODYC(J1) + BODYC(J2))
                      ZF = BODYC(IESC)/ZMU
*       Write E = E_b + E_1 with E_b = -m1*m2/2*a and E_1 = MU3*ER.
                      SEMI = 0.5*RGRAV/(1.0 + ZF*(1.0 + RGRAV*ER/MASS))
                      RY = 1.0/RINV(IB)
                      IF (RY.LT.0.9*SEMI.AND.RI.LT.2.0*RMIN) THEN
                          WRITE (6,33)  NSTEP1, RY/SEMI, RI, RDOT**2,
     &                                  2.0*BODY(ICH)/RI, SEMI
                          KCASE = 0
                          GO TO 60
                      END IF
                      IF (RI.GT.20.0*SEMI) THEN
*       Delay single BH for small perturbation (orbit may turn around).
                          IF (ISTAR(IESC).EQ.14.AND.
     &                        GPERT.LT.1.0D-03) THEN
                              KCASE = 0
                              GO TO 60
                          ELSE
                              KCASE = 1
                              GO TO 35
                          END IF
                      ELSE
                          KCASE = 0
                          GO TO 60
                      END IF
                  ELSE
                      IF (ISTAR(IESC).EQ.14.AND.
     &                    GPERT.LT.1.0D-03) THEN
                          KCASE = 0
                          GO TO 60
                      ELSE
                          KCASE = 1
                          GO TO 35
                      END IF
                  END IF
              ELSE
                  KCASE = 0
                  GO TO 60
              END IF
          ELSE IF (NN.GE.4) THEN
*       Accept escape of hyperbolic body if separation > RMIN.
              RBIG = 1.0/RINV(IM)
              IF (RBIG.GT.RMIN) THEN
                  KCASE = 1
                  GO TO 35
              END IF
*       Include three-body stability test for distant fourth body.
*             IF (RBIG.GT.0.75*RSUM.AND.NN.EQ.4.AND.
*    &            RSUM.GT.0.5*RMIN) THEN
*                 CALL CHSTAB(ITERM)
*                 IF (ITERM.LT.0) THEN
*                     KCASE = -1
*                     GO TO 60
*                 END IF
*             END IF
*       Skip if largest separation is not dominant.
              IF (RBIG.LT.0.66*RSUM) THEN
                  KCASE = 0
                  GO TO 60
              END IF
          ELSE
*       Note possibility a << RGRAV/6 after strong interaction with NN = 3.
              IB = 3 - IM
              J1 = INAME(IB)
              J2 = INAME(IB+1)
*       Obtain semi-major axis from ENERGY = -(m1*m2 + (m1+m2)*m3)/RGRAV.
              ZMU = BODYC(J1)*BODYC(J2)/(BODYC(J1) + BODYC(J2))
              ZF = BODYC(IESC)/ZMU
              ER = 0.5*RDOT**2 - BODY(ICH)/RI
*       Write E = E_b + E_1 with E_b = -m1*m2/2*a and E_1 = MU3*ER.
              SEMI = 0.5*RGRAV/(1.0 + ZF*(1.0 + RGRAV*ER/MASS))
              RY = 1.0/RINV(IB)
*       Note RGRAV = MIN(0.5*RSUM,RGRAV) is OK and RGRAV adjusted on QPMOD.
              IF (RY.LT.0.8*SEMI.AND.RI.LT.2.0*RMIN) THEN
                  WRITE (6,33)  NSTEP1, RY/SEMI, RI, RDOT**2,
     &                          2.0*BODY(ICH)/RI, SEMI
   33             FORMAT (' CHAIN DELAY    # R/A RI RD2 VP2 A ',
     &                                     I8,F6.2,1P,4E9.1)
                  KCASE = 0
                  GO TO 60
              END IF
          END IF
*
*       Check that escaper is well separated (ratio > 2).
          RR = RSUM - 1.0/RINV(IM)
          RATIO = 1.0/(RINV(IM)*RR)
          IF (RATIO.LT.2.0.AND.RSUM.LT.RMIN) THEN
              KCASE = 0
              GO TO 60
          END IF
*
   35     HI = 0.5*RDOT**2 - BODY(ICH)/RI
          IF (HI.LT.0.0) THEN
*       Skip small radial velocity (0.01 circular value).
              RATIO = RI*RDOT**2/BODY(ICH)
              IF (RATIO.LT.0.0001) THEN
                  KCASE = 0
                  GO TO 60
              END IF
              VINF = 0.0
              RX = -BODY(ICH)/HI
              IF (RI.LT.5.0*RMIN.AND.GPERT.LT.1.0D-03) THEN
                  KCASE = 0
                  GO TO 60
              END IF
          ELSE
              VINF = SQRT(2.0*HI)*VSTAR
              RX = 10.0
          END IF
          IF ((RX.LT.5.0*RMIN.AND.GPERT.LT.1.0D-04).OR.
     &        (HI.GT.0.0.AND.RI.LT.5.0*RMIN)) THEN
              KCASE = 0
              GO TO 60
          END IF
          IF (KZ(30).GT.1.OR.VINF.GT.2.0) THEN
              IF (NAMEC(IESC).NE.NAMESC2) THEN
              WRITE (6,36)  IESC, NAMEC(IESC), RI, RDOT**2,
     &                      2.0*BODY(ICH)/RI, VINF, GPERT
   36         FORMAT (' CHAIN ESCAPE:    IESC NM RI RDOT2 2*M/R VF GP ',
     &                                   I3,I6,1P,3E9.1,0P,F6.1,1P,E9.1)
                  NAMESC2 = NAMEC(IESC)
              END IF
*       Ensure single body is removed in case of wide binary.
              JESC = 0
              IF (VINF.GT.2.0) GO TO 60
          END IF
      ELSE
          KCASE = 0
          GO TO 60
      END IF
*
*       Reduce chain membership (NCH > 3) or specify termination.
   40 IF (NCH.GE.3) THEN
          IF (VINF.GT.0.0.AND.RDOT.GT.0.0) THEN
              IF (GPERT.GT.0.05.AND.RI.GT.40.0*RMIN) THEN
                  GO TO 60
              END IF
          END IF
          IF (GPERT.GT.0.005.AND.RI.GE.30.0*RMIN) GO TO 60
*       Subtract largest chain distance from system size (also binary).
*         IM = ISORT(1)
*         RSUM = RSUM - 1.0/RINV(IM) - RB
          IF (JESC.EQ.0) THEN
              HI = 0.5*RDOT**2 - BODY(ICH)/RI
              RX = 10.0
              IF (HI.LT.0.0) RX = -BODY(ICH)/HI
              IF ((RX.LT.25.0*RMIN.AND.GPERT.LT.0.001).OR.
     &            (RDOT.LT.0.0.AND.NCH.LE.5.AND.GPERT.LT.0.01)) THEN
                  KCASE = 0
                  GO TO 60
              END IF
*       Skip if escaper would increase the perturbation by factor 1.5.
*             GX = 2.0*BODYC(IESC)/(MASS-BODYC(IESC))*((RSUM-RI)/RI)**3
              CALL CPERTJ(IESC,GX)
              IF ((GPERT.LT.1.5*GX.AND.RI.LT.20.0*RMIN).OR.
     &            (ISTAR(IESC).EQ.14.AND.GPERT.LT.0.01).OR.
     &            (RI.LT.10.0*RMIN.AND.GPERT.LT.0.01)) THEN
                  KCASE = 0
                  GO TO 60
              END IF
*       Check distance to nearest member.
              RP2 = 1.0
              DO 45 I = 1,NN
                  IF (I.EQ.IESC) GO TO 45
                  RIJ2 = 0.0
                  DO 44 K = 1,3
                      RIJ2 = RIJ2 + (X4(K,I) - X4(K,IESC))**2
   44             CONTINUE
                  IF (RIJ2.LT.RP2) THEN
                      RP2 = RIJ2
                  END IF
   45         CONTINUE
              IF (RP2.LT.4.0*RMIN**2.AND.RI.LT.10.0*RMIN) THEN
                  KCASE = 0
                  GO TO 60
              END IF
              RP = SQRT(RP2)
              IF (GPERT.GT.0.05.AND.RP.GT.40.0*RMIN) KCASE = -4
*             WRITE (6,46) NSTEP1, KCASE, NCH, NAMEC(IESC), HI, RI,
*    &                     GPERT, GX, RP, RDOT
*  46         FORMAT (' REDUCE!   # KCASE NCH NMC HI RI GP GX RP RD ',
*    &                            I8,2I4,I6,F8.3,1P,E10.2,4E9.1)
          END IF
*
          IF (JESC.GT.0) THEN
              WRITE (6,48)  NCH,NAMEC(IESC),NAMEC(JESC), 1.0/RINV(IBIN)
   48         FORMAT (' BINARY!!!   NCH NAM RB  ',I4,2I6,1P,E10.2)
              WRITE (6,49)  NSTEP1, IESC, JESC, GPERT,
     &                      (1.0/RINV(KK),KK=1,NN-1)
   49         FORMAT (' # GP R ',I8,2I4,1P,7E10.2)
              KCASE = -3         ! Set condition for enforced termination.
          END IF
          IF (KZ(30).GT.3) THEN
              WRITE (6,55)  TIME+TOFF, NN, NAMEC(IESC), RI, RDOT
   55         FORMAT (' REDUCE    T NCH NAM RI RD ',
     &                            F10.3,I4,I7,1P,2E10.2)
          END IF
*
*         CALL REDUCE(IESC,JESC,ISUB)
*     ELSE
*         KCASE = -1
      END IF
*
*       Set phase indicator < 0 (not clear if needed).
   50 IPHASE = -1
*
   60 CONTINUE
*       Ensure no action on zero indices (otherwise looping).
      IF (IESC.EQ.0.AND.JESC.EQ.0) KCASE = 0
*       Allow larger termination size for improved initialization.
      IF (KCASE.EQ.-3.OR.NCH.GT.3) THEN
          IF (GPERT.LT.1.0D-04.AND.RI.LT.10.0*RMIN) KCASE = 0
          IF (GPERT.LT.1.0D-02.AND.RI.LT.5.0*RMIN) KCASE = 0
          IF (NN.GT.3.AND.RI.LT.5.0*RMIN) KCASE = 0
*       Include safety limit to enforce termination.
          IF (RI.GT.10.0*RMIN) KCASE = -3
      END IF
*
*       Enforce binary removal for dominant middle distance.
      RX = MAX(1.0/RINV(1),1.0/RINV(NN-1))
      IF (NN.EQ.4.AND.((1.0/RINV(2).GT.5.0*RX.AND.GPERT.GT.0.01).OR.
     &    (GPERT.GT.0.1.AND.1.0/RINV(2).GT.2.0*RX))) THEN
          R1 = 1.0/RINV(1)
          R2 = 1.0/RINV(2)
          R3 = 1.0/RINV(3)
          KCASE = -3
*       Repeat specification of indices in case set to zero.
          IF (ISORT(1).EQ.2) THEN
              IESC = INAME(1)
              JESC = INAME(2)
          END IF
          WRITE (6,75)  KCASE, IESC, JESC, R1, R2, R3
   75     FORMAT (' FOUR-BODY TERM!    KC IE JE R1 R2 R3 ',3I4,1P,3E9.1)
          IF (R2.GT.1.0) STOP
      ELSE IF (NN.EQ.4.AND.((RX.GT.5.0*RMIN.AND.GPERT.GT.0.02).OR.
     &    (GPERT.GT.0.1.AND.RX.GT.3.0*RMIN))) THEN
*       Check for removal of first or last distant single. (SJA 12/2016)
          CALL HPSORT(NN-1,RINV,ISORT)
          IESC = 0
          IF (ISORT(1).EQ.1) THEN
              IESC = INAME(1)
          ELSE IF (ISORT(1).EQ.NN-1) THEN
              IESC = INAME(NN)
          END IF
          IF (IESC.GT.0) THEN     ! Refinement 1/2017 avoids ABSORB repeats.
              RD = X4(1,IESC)*XDOT4(1,IESC) + X4(2,IESC)*XDOT4(2,IESC)
     &                                      + X4(3,IESC)*XDOT4(3,IESC)
*       Ensure positive radial velocity for single escaper (include ABSORB).
              IF (RD.GT.0.0) THEN
                  KCASE = -3
                  WRITE (6,80)  KCASE, IESC, JESC, RX, RX/RD
   80             FORMAT (' FOUR-BODY TERM!    KC IE JE RX RDOT ',
     &                                         3I4,1P,E10.2,0P,F7.2)
                  JESC = 0
              ELSE
                  KCASE = 0
              END IF
          END IF
      ELSE IF (NN.EQ.4.AND.IESC.GT.0) THEN
          RD = X4(1,IESC)*XDOT4(1,IESC) + X4(2,IESC)*XDOT4(2,IESC)
     &                                  + X4(3,IESC)*XDOT4(3,IESC)
          IF (1.0/RINV(2).GT.20.0*RMIN.AND.RD.GT.0.0) THEN
              KCASE = -3
          ELSE
              KCASE = 0
          END IF
      END IF
*
*       Include provisional check of three-body system using perturbation.
      IF (NN.EQ.3) THEN
          IF (RX.GT.10.0*RMIN.AND.IESC.GT.0) THEN
              RD = X4(1,IESC)*XDOT4(1,IESC) + X4(2,IESC)*XDOT4(2,IESC)
     &                                      + X4(3,IESC)*XDOT4(3,IESC)
              IF (GPERT.GT.0.001.AND.RD.GT.0.0) THEN
                  KCASE = -2
              ELSE
                  KCASE = 0
              END IF
          ELSE IF (GPERT.LT.0.001) THEN
              KCASE = 0
          END IF
      END IF
*
*       Include algorithm for NN = 5 (JESC = 0 is easiest).
      IF (NN.EQ.5.AND.GPERT.GT.0.01) THEN
*       Set current chain distances and corresponding inverse sorting.
          R1 = 1.0/RINV(1)
          R2 = 1.0/RINV(2)
          R3 = 1.0/RINV(3)
          R4 = 1.0/RINV(4)
          CALL HPSORT(NN-1,RINV,ISORT)
          IF (ISORT(1).EQ.1) THEN
              IESC = INAME(1)
              JESC = INAME(2)
              KCASE = -3
              IF (R2.LT.0.1*R1) THEN
                  JESC = 0
                  RD = X4(1,IESC)*XDOT4(1,IESC)+X4(2,IESC)*XDOT4(2,IESC)
     &                                         +X4(3,IESC)*XDOT4(3,IESC)
                  IF (RD.LT.0.0) KCASE = 0
              END IF
          ELSE IF (ISORT(1).EQ.4) THEN
              IESC = INAME(NN-1)
              JESC = INAME(NN)
              KCASE = -3
              IF (R3.LT.0.1*R4) THEN
                  IESC = JESC
                  JESC = 0
                  RD = X4(1,IESC)*XDOT4(1,IESC)+X4(2,IESC)*XDOT4(2,IESC)
     &                                         +X4(3,IESC)*XDOT4(3,IESC)
                  IF (RD.LT.0.0) KCASE = 0
              END IF
          ELSE IF (ISORT(1).EQ.2) THEN
              IESC = INAME(1)
              JESC = INAME(2)
              KCASE = -3
          ELSE IF (ISORT(1).EQ.3) THEN
              IESC = INAME(NN-1)
              JESC = INAME(NN)
              KCASE = -3
          ELSE
              IESC = 0
          END IF
      ELSE IF (NN.EQ.5.AND.(GPERT.GT.0.1.OR.RSUM.GT.20.0*RMIN)) THEN
*       Ensure membership reduction for large perturbation (any IESC > 0).
          KCASE = -2
          IESC = 1
      ELSE IF (NN.EQ.5) THEN
          IESC = 0
      END IF
*
*       Include trigger for rare case of NN = 6 (it happened as 3*KS).
      IF (NN.EQ.6.AND.GPERT.GT.0.01) THEN
          CALL HPSORT(NN-1,RINV,ISORT)
          WRITE (6,85)  ISORT(1), GPERT, (1.0/RINV(K),K=1,NN-1)
   85     FORMAT (' SIX CHAIN    IS GP R ',I4,F7.3,1P,5E10.2)
          CALL FLUSH(6)
          KCASE = -2
          IESC = 1
      END IF
*
*       Avoid repeated ABSORB and ESCAPE by ensuring RD > 0.
      IF (KCASE.EQ.-3.AND.IESC.GT.0.AND.JESC.EQ.0) THEN
          RD = X4(1,IESC)*XDOT4(1,IESC) + X4(2,IESC)*XDOT4(2,IESC)
     &                                  + X4(3,IESC)*XDOT4(3,IESC)
          IF (RD.LT.0.0) KCASE = 0
      END IF
      IF (KCASE.EQ.0) IESC = 0  ! Otherwise different path in CHAIN.
*
*       Include inconsistency check.
      IF (IESC.EQ.0) KCASE = 0
*
*       Delay for weakly perturbed higher-order systems inside 30*RMIN.
      IF (KCASE.LT.0.AND.NN.GE.4) THEN
          IF (GPERT.LT.0.001.AND.RSUM.LT.30.0*RMIN) KCASE = 0
      END IF
*
      RETURN
*
      END
