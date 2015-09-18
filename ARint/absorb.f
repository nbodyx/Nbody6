      SUBROUTINE ABSORB(ISUB)
*
*
*       Absorption of chain member(s).
*       -----------------------------
*
      INCLUDE 'common6.h'
        REAL*8  M,MASS,MC,MMIJ
        PARAMETER (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &  NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
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
      COMMON/SOFT/  EPS2
      REAL*8  XCM(3),VCM(3)
      REAL*8  X0S(3),V0S(3),XS(3),VS(3),X20S(3),V20S(3)
*
*
*       Define discrete time for new polynomial (DT2 < 0 is OK).
      TIME0 = TIME
      TIMET = TIME
      TIME = T0S(ISUB) + TIMEC
      TIME = TBLOCK
*       Re-define initial epoch for consistency (ignore phase error).
***   T0S(ISUB) = TIME - TIMEC  !looks dangerous!!!
*
      IF (KZ(30).GT.2) THEN
          WRITE (6,1)  TIME0+TOFF, TIME+TOFF, TBLOCK
    1     FORMAT (' ABSORB:   TIME0 TIME TBLOCK ',3F10.6)
      END IF
*
      ZM0S = BODY(ICH)
      ZM20S = BODY(JCLOSE)
      DO 2 K = 1,3
          X0S(K) = X(K,ICH)
          V0S(K) = XDOT(K,ICH)
          X20S(K) = X(K,JCLOSE)
          V20S(K) = XDOT(K,JCLOSE)
    2 CONTINUE
*
*       Increase membership of chain (JCLOSE: single body or KS pair).
      NCH0 = NCH
      CALL SETSYS
      INAME(NCH) = NCH  ! what is this??
*
*       Improve coordinates & velocities of c.m. body to order F3DOT.
      CALL XVPRED(ICH,-1)
*     TIME = TIMET
*
      SUM = 0.0
      DO 5 K = 1,3
          XCM(K) = 0.0
          VCM(K) = 0.0
    5 CONTINUE
*
*       Accumulate mass-weighted moments of absorbed particle(s).
      DO 15 L = NCH0+1,NCH
          J = JLIST(L)
*       Ensure a large look-up value for AR_Chain ghost.
          IF (KZ(19).GE.3.AND.KZ(11).GT.0) THEN
              TEV(J) = MAX(TIME + 1.0D+04,TEV(J))
          END IF
          SUM = SUM + BODY(J)
          DO 10 K = 1,3
              XCM(K) = XCM(K) + BODY(J)*X(K,J)
              VCM(K) = VCM(K) + BODY(J)*XDOT(K,J)
   10     CONTINUE
   15 CONTINUE
*
*       Form combined c.m. of old chain and new perturber(s).
      DO 20 K = 1,3
          XCM(K) = (BODY(ICH)*X(K,ICH) + XCM(K))/(BODY(ICH) + SUM)
          VCM(K) = (BODY(ICH)*XDOT(K,ICH) + VCM(K))/(BODY(ICH) + SUM)
   20 CONTINUE
*
*       Define new relative coordinates & velocities and add to chain.
      LK = 3*NCH0
      DO 30 L = NCH0+1,NCH
          J = JLIST(L)
          SIZE(L) = RADIUS(J)
          ISTAR(L) = KSTAR(J)
          DO 25 K = 1,3
              LK = LK + 1
              XCH(LK) = X(K,J) - XCM(K)
              VCH(LK) = XDOT(K,J) - VCM(K)
   25     CONTINUE
   30 CONTINUE
*
*       Re-define old chain variables with respect to new c.m.
      LK = 0
      DECM = 0.0
      DO 40 L = 1,NCH0
          DO 35 K = 1,3
              LK = LK + 1
              XCH(LK) = XCH(LK) - (XCM(K) - X(K,ICH))
              VCH(LK) = VCH(LK) - (VCM(K) - XDOT(K,ICH))
              DECM = DECM + BODYC(L)*(VCM(K) - XDOT(K,ICH))**2
   35     CONTINUE
   40 CONTINUE
*
*       Create ghost particle(s) and remove from perturber lists.
      DO 50 L = NCH0+1,NCH
          J = JLIST(L)
          CALL GHOST(J)
   50 CONTINUE
*
*       Update total mass and initialize new c.m. body variables (also X0).
      BODY(ICH) = BODY(ICH) + SUM
      CM(7) = BODY(ICH)
      T0(ICH) = TIME
      LX = 3*NCH0
      DO 55 K = 1,3
          X(K,ICH) = XCM(K)
          X0(K,ICH) = XCM(K)
          XDOT(K,ICH) = VCM(K)
          X0DOT(K,ICH) = VCM(K)
   55 CONTINUE
*
*       Perform re-initialization of c.m. polynomials & perturber list.
      CALL REINIT(ISUB)
*
*       Update energies using regular expression (NB! include binary).
      DKE = 0.0
      DPOT = 0.0
      JX = NCH0 + 1
   60 LK = 0
      DO 62 K = 1,3
          DKE = DKE + 0.5*BODYC(JX)*VCH(LX+K)**2
   62 CONTINUE
      DO 70 L = 1,NCH0
          RIJ2 = 0.0
          DO 65 K = 1,3
              LK = LK + 1
              RIJ2 = RIJ2 + (XCH(LK) - XCH(LX+K))**2
   65     CONTINUE
          DPOT = DPOT + BODYC(L)*BODYC(JX)/SQRT(RIJ2)
   70 CONTINUE
*
*       Include mutual two-body interaction for second KS component.
      IF (JX.LT.NCH) THEN
      DP = DPOT
          RIJ2 = 0.0
          DO 75 K = 1,3
              RIJ2 = RIJ2 + (XCH(LX+K) - XCH(LX+3+K))**2
   75     CONTINUE
          DPOT = DPOT + BODYC(JX)*BODYC(NCH)/SQRT(RIJ2)
          JX = NCH
          LX = LX + 3
          GO TO 60
      END IF
*
*       Add kinetic energy corrections to current energy (checked OK).
      DECM = 0.5*DECM
      ENERGY = ENERGY + DKE - DPOT + DECM
*     DE = DKE - DPOT
*     CALL CONST(XCH,VCH,M,NCH,ENER1,G0,ALAG)
*     WRITE (6,72)  (ENERGY+EnerGR - ENER0), DE, DECM, DE+DECM
*  72 FORMAT (' ABS ERROR!   DE DEJ DEC DE+DEC  ',1P,4E16.6,2E10.2)
*     WRITE (6,74)  (ENER1 - ENERGY)/ENER1
*  74 FORMAT (' CHECK ABS   DE/E  ',1P,E10.2)
*
*       Include experimental procedure to catch new perturbers.
      NNB = LIST(1,ICH)
      NNB = 0
      NP = LISTC(1)
      DO 230 L = 2,NNB+1
          I = LIST(L,ICH)
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
          ZMS = ZM20S
          DO 215 K = 1,3
              XS(K) = X20S(K)
              VS(K) = V20S(K)
  215     CONTINUE
          CALL FFDOT(I,ZMS,XS,VS,SS)
          SS = 1.0
          ZMS = BODY(ICH)
          DO 220 K = 1,3
              XS(K) = X(K,ICH)
              VS(K) = XDOT(K,ICH)
  220     CONTINUE
          CALL FFDOT(I,ZMS,XS,VS,SS)
      RIJ2 = 0.0
      DO K = 1,3
      RIJ2 = RIJ2 + (X(K,I) - X(K,ICH))**2
      END DO
      RIJ = SQRT(RIJ2)
      WRITE (6,216)  NAME(I), FX1, F(1,I)-FX1, STEP(I),RIJ
  216 FORMAT (' ABSORB   NM FX0 DFX S RIJ ',I6,1P,4E10.2)
  230 CONTINUE
*       Check centre of mass condition (suppressed after testing).
*     DO 80 K = 1,6
*         CM(K) = 0.0
*  80 CONTINUE
*
*     LK = 0
*     DO 90 L = 1,NCH
*         DO 85 K = 1,3
*             LK = LK + 1
*             CM(K) = CM(K) + BODYC(L)*XCH(LK)
*             CM(K+3) = CM(K+3) + BODYC(L)*VCH(LK)
*  85     CONTINUE
*  90 CONTINUE
*
*     DO 95 K = 1,6
*         CM(K) = CM(K)/CM(7)
*  95 CONTINUE
*
*     WRITE (6,99)  (CM(K),K=1,6)
*  99 FORMAT (' ABSORB:   CM ',1P,6E9.1)
*     CALL FLUSH(6)
*
      RETURN
*
      END
