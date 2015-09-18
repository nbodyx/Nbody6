      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/PREDICT/ TPRED(NMAX)
      INTEGER ISTAT(NMAX)
      PARAMETER (NIMAX=1024,NPMAX=16)
      REAL*8   H2I(NIMAX),XI(3,NIMAX),VI(3,NIMAX),GPUACC(3,NIMAX),
     &         GPUJRK(3,NIMAX),GPUPHI(NIMAX),GF(3,NMAX),GFD(3,NMAX),
     &         DTR(NIMAX)
      INTEGER  NXTLST(NMAX),LISTQ(NMAX),NL(20)
      INTEGER  IRR(NMAX),IREG(NMAX),LISTGP(LMAX,NIMAX)
      LOGICAL LOOP,LSTEPM
      SAVE IQ,ICALL,NQ,LQ,LOOP,LSTEPM,STEPM,ISAVE,JSAVE,ISTART,TBH
      DATA IQ,ICALL,LQ,LOOP,LSTEPM,STEPM /0,2,11,.TRUE.,.FALSE.,0.03125/
      DATA ISAVE,JSAVE,ISTART /0,0,0/
      SAVE TLISTQ,TMIN,ICOMP0,JCOMP0,NPACT,GF,GFD
      SAVE NXTLST,LISTQ,IRR,IREG,LISTGP,ISTAT
*     SAVE CPRED,CPRED2,CPRED3
*     DATA CPRED,CPRED2,CPRED3 /0.0D0,0.0D0,0.0D0/
*
*
*       Update quantized value of STEPM for large N (first time only).
      IF (.NOT.LSTEPM.AND.NZERO.GT.1024) THEN
          K = INT((FLOAT(NZERO)/1024.0)**0.333333)
          STEPM = 0.03125D0/2**(K-1)
          LSTEPM = .TRUE.
      END IF
*
*       Open the GPU libraries on each new run (note nnbmax = NN is printed).
      IF (ISTART.EQ.0) THEN
          NN = NTOT + 100
          CALL GPUNB_OPEN(NN)
*       Set larger value for GPUIRR (note further possible increase of NTOT).
          NNN = NTOT + 100
*       Increase NNN further for cold start to deal with many KS.
          IF (QVIR.LT.0.01.AND.TIME.EQ.0.0D0) NNN = NTOT + 500
*       Allow for missing primordial binaries.
          NNN = NNN + NBIN0 - NPAIRS
          CALL GPUIRR_OPEN(NNN,LMAX)
          ISTART = 1
*       Define parameter for predicting active particles in library.
          IF (N.LE.5000) THEN
              NPACT = 75 
          ELSE IF (N.LE.50000) THEN
              NPACT = 100 
          ELSE IF (N.LE.100000) THEN
              NPACT = 150
          ELSE
              NPACT = 200
          END IF
          TBH = TTOT
      END IF
*
*       Search for high velocities after escape or KS/chain termination.
  999 IF (KZ(37).GT.0.AND.(IPHASE.EQ.-1.OR.IPHASE.GE.2)) THEN
          CALL HIVEL(0)
      END IF
*
*       Initialize new force times and predict to previous TIME.
      IF (IPHASE.EQ.1) THEN
*       Treat case of new KS with components from ICOMP0 & JCOMP0.
          I = ICOMP0
*       Avoid updating unchanged locations.
 1000     IF (I.GE.IFIRST) THEN
              TNEW(I) = T0(I) + STEP(I)
              TPRED(I) = -1.0
              CALL JPRED(I)
              CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
     &                                                BODY(I),T0(I))
          END IF
*       Select second component or new c.m. (not connected with ICOMP0).
          IF (I.EQ.ICOMP0) THEN
              I = JCOMP0
              GO TO 1000
          ELSE IF (I.EQ.JCOMP0) THEN
              I = NTOT
              GO TO 1000
          END IF
!$omp parallel do private(J)
          DO 1010 J = IFIRST,NTOT
              CALL GPUIRR_SET_LIST(J,LIST(1,J))
 1010     CONTINUE
!$omp end parallel do
*
*       Remove ICOMP0 & JCOMP0 from LISTQ and add NTOT.
          CALL REPAIR(ICOMP0,JCOMP0,NTOT,0,NQ,LISTQ,TMIN)
*       Check KS termination (locations IFIRST & IFIRST + 1).
      ELSE IF (IPHASE.EQ.2) THEN
*       Update relevant variables for the KS components and later c.m.
          I1 = IFIRST
          I2 = I1 + 1
 1015     DO 1030 I = I1,I2
              TNEW(I) = T0(I) + STEP(I)
              TPRED(I) = -1.0
              CALL JPRED(I)
              CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
     &                                                BODY(I),T0(I))
 1030     CONTINUE
*
*       Include all the more recent c.m. bodies.
          IF (I1.EQ.IFIRST) THEN
              I1 = N + KSPAIR
              I2 = NTOT
              GO TO 1015
*       Note loop is ignored for case of last c.m. (KSPAIR = NPAIRS + 1).
          END IF
*
*       Re-send all neighbour lists using parallel directive (cheap loop).
!$omp parallel do private(I)
          DO 1035 I = IFIRST,NTOT
              CALL GPUIRR_SET_LIST(I,LIST(1,I))
 1035     CONTINUE
!$omp end parallel do
*
*       Add two first single particles to LISTQ and remove terminated c.m.
          CALL REPAIR(N+KSPAIR,0,IFIRST,IFIRST+1,NQ,LISTQ,TMIN)
      END IF
*
*       Update all SET_JP & SET_LIST except for frequent IPHASE = 1 or 2.
      IF (IPHASE.LE.0.OR.IPHASE.GE.3) THEN
*       Include at TIME = 0, CHAIN (7, 8), COLL & COAL, escape & MERGE/RESET.
!$omp parallel do private(I)
*       Note IPHASE = 3 also after ENFORCED KS in ADJUST (if no escape).
          DO 1040 I = IFIRST,NTOT
              CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
     &                                                BODY(I),T0(I))
              CALL GPUIRR_SET_LIST(I,LIST(1,I))
*       Initialize end-point of new times and predict all at previous TIME.
              TNEW(I) = T0(I) + STEP(I)
              TPRED(I) = -1.0
              CALL JPRED(I)
*       Note that this prediction is strictly not necessary.
 1040     CONTINUE
!$omp end parallel do
*       Prescribe level search on return, except for new and terminated KS.
          LOOP = .TRUE.
          TLISTQ = TIME
      END IF
*
*       Reset control & regularization indicators.
      IPHASE = 0
      IKS = 0
      TPREV = TIME
*       Take close encounter step as typical lower value for time-step list.
      CALL STEPK(DTMIN,DTM)
*
*       Determine level for the smallest step (ignore extreme values).
      LQS = 20
      DO 1050 L = 6,20
          IF (DTM.EQ.DTK(L)) THEN
              LQS = L
          END IF
 1050 CONTINUE
*
*       Specify upper level for optimized membership and reset IQ.
      LQB = LQS - 8
      IF (IQ.LT.0) ICALL = 0
      IQ = 0
*
*       Check updating new list of block steps with T0 + STEP =< TLISTQ.
    1 ICALL = ICALL + 1
*       Reset TMIN second & third time after change to catch new chain step.
      IF (TIME.GE.TLISTQ.OR.ICALL.LE.3) THEN
*       Update interval by optimization at major times (sqrt of N-NPAIRS).
          IF (DMOD(TLISTQ,2.0D0).EQ.0.0D0.OR.LOOP) THEN
              LOOP = .FALSE.
              DO 10 L = 1,20
                  NL(L) = 0
   10         CONTINUE
              DO 14 I = IFIRST,NTOT
*       Count steps at five different levels for the smallest values.
                  DO 12 L = LQB,LQS
                      IF (STEP(I).LT.DTK(L)) NL(L) = NL(L) + 1
   12             CONTINUE
   14         CONTINUE
              NLSUM = 0
*       Determine interval by summing smallest steps until near sqrt(N-N_b).
              NSQ = INT(SQRT(FLOAT(N - NPAIRS)))
              LQ = LQS
              DO 15 L = LQS,LQB,-1
                  NLSUM = NLSUM + NL(L)
                  IF (NLSUM.LE.NSQ) LQ = L
   15         CONTINUE
*             WRITE (6,16)  TIME+TOFF,NQ,NLSUM,LQ,(NL(K),K=LQB,LQS)
*  16         FORMAT (' LEVEL CHECK:    T NQ NLSUM LQ NL  ',
*    &                                  F9.3,3I5,2X,7I4)
          END IF
*
*       Increase interval by optimized value.
          NQ = 0
          TMIN = 1.0D+10
   18     TLISTQ = TLISTQ + DTK(LQ)
          DO 20 I = IFIRST,NTOT
              IF (TNEW(I).LE.TLISTQ) THEN
                  NQ = NQ + 1
                  LISTQ(NQ) = I
                  TMIN = MIN(TNEW(I),TMIN)
              END IF
   20     CONTINUE
*       Increase interval in rare case of zero membership.
          IF (NQ.EQ.0) GO TO 18
*       Make a slight adjustment for high levels and small membership.
          IF (LQ.LE.15.AND.NQ.LE.2) LQ = LQ - 1
      END IF
*
*     TTT = TIME     ! Used for profiling (needs wtime.o).
*     TLAST = TIME   ! Activate on negative time test.
*
*       Find all particles in next block (TNEW = TMIN) and set TIME.
      CALL  INEXT(NQ,LISTQ,TMIN,NXTLEN,NXTLST)
*
*       Set new time and save block time (for regularization terminations).
      I = NXTLST(1)
      TIME = T0(I) + STEP(I)
      TBLOCK = TIME
*
*       Include diagnostics for negative TIME and block-step diagnostics.
*     IF (TIME.LT.TLAST.AND.TIME.LT.TLISTQ) THEN
*     WRITE (6,1300) NXTLEN, NQ, TIME, TIME-TLAST,STEP(I)
*1300 FORMAT (' NEGATIVE!    LEN NQ T T-TL S ',2I6,F10.5,1P,2E10.2)
*     CALL FLUSH(6)
*     END IF
*     WRITE (6,22)  I, NXTLEN, NSTEPU, NSTEPI, TIME, STEP(I), STEPR(I)
*  22 FORMAT (' INTGRT   I LEN #U #I T S SR  ',2I6,2I11,F10.4,1P,2E10.2)
*     CALL FLUSH(6)
*
*       Terminate on small irregular time-step (means problems).
      IF (STEP(I).LT.1.0D-11) THEN
          WRITE (6,24)  I, NAME(I), NXTLEN, NSTEPR, STEP(I), STEPR(I)
   24     FORMAT (' SMALL STEP!!    I NAME LEN #R SI SR ',
     &                              3I6,I11,1P,2E10.2)
          STOP
      END IF
*
*     TT2 = DBLTIM()
*     CPRED3 = CPRED3 + (TT2 - TT1)
*       Re-determine list if current time exceeds boundary.
      IF (TIME.GT.TLISTQ) GO TO 1
*
*       Check option for advancing interstellar clouds.
      IF (KZ(13).GT.0) THEN
          CALL CLINT
      END IF
*
*       Check optional integration of cluster guiding centre.
      IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
          IF (KZ(14).EQ.3.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
              CALL GCINT
          END IF
*       Include mass loss by gas expulsion (Kroupa et al. MN 321, 699).
          IF (MPDOT.GT.0.0D0.AND.TIME + TOFF.GT.TDELAY) THEN
              CALL PLPOT1(PHI1)
              MP = MP0/(1.0 + MPDOT*(TIME + TOFF - TDELAY))
*       Replace by exponential mass loss for faster decrease.
*             DT = TIME + TOFF - TDELAY
*             MP = MP0*EXP(-MPDOT*DT)
              CALL PLPOT1(PHI2)
*       Add differential correction for energy conservation.
              EMDOT = EMDOT + (PHI1 - PHI2)
          END IF
      END IF
*
*       Include commensurability test (may be suppressed if no problems).
*     IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
*         WRITE (6,25)  I, NAME(I), NSTEPI, TIME, STEP(I), TIME/STEP(I)
*  25     FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
*    &                        2I6,I11,F12.5,1P,E9.1,0P,F16.4)
*         STOP
*     END IF
*
*       Check for new regularization at end of previous block.
      IF (IKS.GT.0) THEN
          IPHASE = 1
*       Copy the saved KS component indices and time.
          ICOMP = ISAVE
          JCOMP = JSAVE
          TIME = TSAVE
          ICOMP0 = ICOMP
          JCOMP0 = JCOMP
          GO TO 100
      END IF
*
*     TT1 = DBLTIM()
*       Check next adjust time before beginning a new block.
      IF (TIME.GT.TADJ) THEN
          TIME = TADJ
          IPHASE = 3
          GO TO 100
      END IF
*
*       Check output time in case DTADJ & DELTAT not commensurate.
      IF (TIME.GT.TNEXT) THEN
          TIME = TNEXT
          CALL OUTPUT
          GO TO 1
      END IF
*
*       See whether to advance ARchain or KS at first new time.
      IF (TIME.GT.TPREV) THEN
          CALL SUBINT(IQ,I10)
*       Check collision/coalescence indicator.
          IF (IQ.LT.0) GO TO 999
      END IF
*
*       Form lists of candidates for new irregular and regular force.
      NFR = 0
      DO 28 L = 1,NXTLEN
          J = NXTLST(L)
          IF (TNEW(J).GE.T0R(J) + STEPR(J)) THEN
              NFR = NFR + 1
              IREG(NFR) = J
              IRR(L) = J
          ELSE
              IRR(L) = 0
          END IF
   28 CONTINUE
* 
*       Decide between predicting <= NPACT active (NFR = 0) or all particles.
      IF (NXTLEN.LE.NPACT.AND.NFR.EQ.0) THEN
*
*       Predict active particles in the GPUIRR library.
          CALL GPUIRR_PRED_ACT(NXTLEN,NXTLST,TIME)
      ELSE
          CALL GPUIRR_PRED_ALL(IFIRST,NTOT,TIME)
          NNPRED = NNPRED + 1
*!$omp parallel do private(S, S1, S2)
*         DO 40 J = IFIRST,NTOT
*             IF (BODY(J).EQ.0.0D0) GO TO 40
*             S = TIME - T0(J)
*             S1 = 1.5*S
*             S2 = 2.0*S
*             X(1,J) = ((FDOT(1,J)*S + F(1,J))*S +X0DOT(1,J))*S +X0(1,J)
*             X(2,J) = ((FDOT(2,J)*S + F(2,J))*S +X0DOT(2,J))*S +X0(2,J)
*             X(3,J) = ((FDOT(3,J)*S + F(3,J))*S +X0DOT(3,J))*S +X0(3,J)
*             XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
*             XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
*             XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
*  40     CONTINUE
*!$omp end parallel do
      END IF
*
*       Save new time (output time at TIME > TADJ) and increase # blocks.
      TPREV = TIME
      NBLOCK = NBLOCK + 1
      TMIN = 1.0D+10
      IKS0 = IKS
*
*       Predict chain variables and perturbers at new block-time.
      IF (NCH.GT.0) THEN
          CALL JPRED(ICH)
          CALL XCPRED(2)
      END IF
*
*       Evaluate all irregular forces & derivatives in the GPUIRR library.
*     DO 46 II = 1,NXTLEN
*         I = NXTLST(II)
*         CALL GPUIRR_FIRR(I,GF(1,II),GFD(1,II))
*  46 CONTINUE
      CALL GPUIRR_FIRR_VEC(NXTLEN,NXTLST,GF,GFD)
*
*       Choose between standard and parallel irregular integration.
      IF (NXTLEN.LE.NPMAX) THEN
*
*       Correct the irregular steps sequentially.
      DO 48 II = 1,NXTLEN
          I = NXTLST(II)
*       Advance the irregular step (no copy of X0 to GPUIRR needed here).
          CALL NBINT(I,IKS,IRR(II),GF(1,II),GFD(1,II))
*         CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
*    &                                            BODY(I),T0(I))
*
*       Save indices and TIME of first KS candidates in the block.
          IF (IKS0.EQ.0.AND.IKS.GT.0) THEN
              ISAVE = ICOMP
              JSAVE = JCOMP
              TSAVE = TIME
              IKS0 = IKS
          END IF
   48 CONTINUE
*
      ELSE
*       Initialize array for repeated cases to ensure thread-safety.
      DO 49 II = 1,NXTLEN
          ISTAT(II) = -1
   49 CONTINUE
*       Perform irregular correction in parallel (suppressed currently).
*!$omp parallel do private(II, I)
      DO 50 II = 1,NXTLEN
          I = NXTLST(II)
          CALL NBINTP(I,IRR(II),GF(1,II),GFD(1,II),ISTAT(II))
   50 CONTINUE
*!$omp end parallel do
*
*       Correct the irregular steps sequentially.
      IPHASE = 0
      IKS = 0
      DO 500 II = 1,NXTLEN
          IF (ISTAT(II).GT.0) THEN
              I = NXTLST(II)
*       Correct exceptional irregular steps in serial.
              CALL NBINT(I,IKS,IRR(II),GF(1,II),GFD(1,II))
*     WRITE (6,490)  TIME, IKS, NAME(I), STEP(I)
* 490 FORMAT (' INTGRT!     T IKS NM SI ',F12.4,I4,I7,1P,E10.2)
          NSTEPQ = NSTEPQ + 1
          END IF
  500 CONTINUE
      END IF
      NSTEPI = NSTEPI + NXTLEN
*
*       Check regular force updates (NFR members on block-step #NBLCKR).
      IF (NFR.GT.0) THEN
      NBLCKR = NBLCKR + 1
      NN = NTOT - IFIRST + 1
*     TT3 = DBLTIM()
*
*       Predict all particles (except TPRED=TIME) in C++ on host.
      CALL CXVPRED(IFIRST,NTOT,TIME,T0,X0,X0DOT,F,FDOT,X,XDOT,TPRED)
*       Note it would be wrong to predict TPRED=TIME by corrected X0 & X0DOT.
*
*     TT4 = DBLTIM()
*     CPRED = CPRED + (TT4 - TT3)
*     IF (DMOD(TIME,2.0D0).EQ.0.0D0) WRITE (6,540)  CPRED
* 540 FORMAT (' TIMING CXVPRED  ',1P,E10.2)
*       Send all single particles and c.m. bodies to the GPU.
      CALL GPUNB_SEND(NN,BODY(IFIRST),X(1,IFIRST),XDOT(1,IFIRST))
*
*       Perform regular force loop.
      NOFL(1) = 0
      NOFL2 = 0
      JNEXT = 0
      DO 55 II = 1,NFR,NIMAX
  550     NI = MIN(NFR-JNEXT,NIMAX)
*       Copy neighbour radius, STEPR and state vector for each block.
!$omp parallel do private(LL, I, K)
          DO 52 LL = 1,NI
              I = IREG(JNEXT+LL)
              H2I(LL) = RS(I)**2
              DTR(LL) = STEPR(I)
              DO 51 K = 1,3
                  XI(K,LL) = X(K,I)
                  VI(K,LL) = XDOT(K,I)
   51         CONTINUE
   52     CONTINUE
!$omp end parallel do
*
*       Evaluate forces, derivatives and neighbour lists for new block.
          CALL GPUNB_REGF(NI,H2I,DTR,XI,VI,GPUACC,GPUJRK,GPUPHI,LMAX,
     &                                                   NBMAX,LISTGP)
*       Copy neighbour lists from the GPU after overflow check.
          DO 54 LL = 1,NI
              I = IREG(JNEXT + LL)
              NNB = LISTGP(1,LL)
*       Repeat last block with reduced RS(I) on NNB < 0 (at end of loop).
              IF (NNB.LT.0) THEN
                  RI2 = (X(1,I)-RDENS(1))**2 + (X(2,I)-RDENS(2))**2 +
     &                                         (X(3,I)-RDENS(3))**2
                  WRITE (41,556)  NSTEPR, NAME(I), LIST(1,I), NNB,
     &                            RS(I), SQRT(RI2)
  556             FORMAT (' OVERFLOW!   #R NAME NB0 NB RS ri ',
     &                                  I11,I6,2I5,2F8.3)
                  CALL FLUSH(41)
                  NOFL(1) = NOFL(1) + 1
                  RS(I) = 0.9*RS(I)
                  H2I(LL) = RS(I)**2
                  NOFL2 = NOFL2 + 1
              END IF
   54     CONTINUE
          IF (NOFL2.GT.0) THEN
              NOFL2 = 0
              GO TO 550
           END IF
!$omp parallel do private(LL, I, ITEMP, NNB, L1, L)
          DO 56 LL = 1,NI
              I = IREG(JNEXT + LL)
              NNB = LISTGP(1,LL)
              L1 = 1
              DO 53 L = 2,NNB+1
*       Note GPU address starts from 0 (hence add IFIRST to neighbour list).
                  ITEMP = LISTGP(L,LL) + IFIRST
                  IF (ITEMP.NE.I) THEN
                      L1 = L1 + 1
                      LISTGP(L1,LL) = ITEMP
                  END IF
   53         CONTINUE
              LISTGP(1,LL) = L1 - 1
              CALL GPUIRR_SET_LIST(I,LISTGP(1,LL))
   56     CONTINUE
!$omp end parallel do
*
*       Evaluate current irregular forces by vector procedure.
          CALL GPUIRR_FIRR_VEC(NI,IREG(II),GF(1,1),GFD(1,1))
*       NB! Note suppression of parallel procedure here and for nbintp.f.
!$omp parallel do private(I, LX)
          DO 57 LL = 1,NI
              I = IREG(JNEXT+LL)
*       Send new irregular force and perform correction.
*             CALL GPUIRR_FIRR(I,GF(1,LL),GFD(1,LL))
              CALL GPUCOR(I,XI(1,LL),VI(1,LL),GPUACC(1,LL),GPUJRK(1,LL),
     &                               GF(1,LL),GFD(1,LL),LISTGP(1,LL),LX)
*       Update neighbour lists in GPUIRR library (only if changed: LX > 0).
              IF (LX.GT.0) THEN
                  CALL GPUIRR_SET_LIST(I,LIST(1,I))
                  NBPRED = NBPRED + 1
              END IF
*       Note possible errors in the potential when using predicted values.
*             POT = POT + BODY(I)*GPUPHI(LL)
   57     CONTINUE
!$omp end parallel do
          JNEXT = JNEXT + NI
   55 CONTINUE
*
*       Accumulate the sum of overflows (NOFL(1) holds current number).
      NOFL(2) = NOFL(2) + NOFL(1)
      NSTEPR = NSTEPR + NFR
*
      END IF
*
*       Determine next block time (note STEP may shrink in GPUCOR).
      DO 350 L = 1,NXTLEN
          I = NXTLST(L)
          IF (TNEW(I).LT.TMIN) THEN
              TMIN = TNEW(I)
          END IF
  350 CONTINUE
*
*       Copy current coordinates & velocities from corrected values, IQ.NE.0.
      IF (IKS.GT.0.OR.IQ.NE.0.OR.TIME.GE.TADJ.OR.TIME.GE.TNEXT.OR.
     &   (KZ(19).GE.3.AND.DMOD(TIME,STEPX).EQ.0.0D0)) THEN
*       Note: need copy X & XDOT at output (skip in XVPRED; also cf. #31 >0).
!$omp parallel do private(I, L, K)
      DO 360 L = 1,NXTLEN
          I = NXTLST(L)
          DO 58 K = 1,3
              X(K,I) = X0(K,I)
              XDOT(K,I) = X0DOT(K,I)
   58     CONTINUE
  360 CONTINUE
!$omp end parallel do
      END IF
*
*       Send corrected active particles to GPUIRR library.
!$omp parallel do private(I, L)
      DO 60 L = 1,NXTLEN
          I = NXTLST(L)
          CALL GPUIRR_SET_JP(I,X0(1,I),X0DOT(1,I),F(1,I),FDOT(1,I),
     &                                            BODY(I),T0(I))
   60 CONTINUE
!$omp end parallel do
*
*     TT2 = DBLTIM()
*     CPRED2 = CPRED2 + (TT2 - TT1)
*     IF (TIME.EQ.TCRIT) WRITE (6,380)  CPRED2, TT2, CPRED, CPRED3
* 380 FORMAT (' TIMING     CPRED2 TT2 CPRED CPRED3 ',1P,4E12.2)
*
*       Check integration of tidal tail members.
      IF (NTAIL.GT.0) THEN
*       Allow large quantized interval with internal iteration.
          IF (DMOD(TIME,0.25D0).EQ.0.0D0) THEN
              DO 65 J = ITAIL0,NTTOT
                  IF (TNEW(J).LE.TIME) THEN
                      CALL NTINT(J)
                  END IF
   65         CONTINUE
          END IF
      END IF
*
*       Exit on KS termination, new multiple regularization or merger.
      IF (IQ.NE.0) THEN
          NBPREV = 0
          IF (IQ.GE.4.AND.IQ.NE.7) THEN
              CALL DELAY(IQ,-1)
          ELSE
*       Ensure correct KS index (KSPAIR may denote second termination).
              KSPAIR = KVEC(I10)
              IPHASE = IQ
          END IF
          GO TO 100
      END IF
*
*       Perform optional check on high-velocity particles at major times.
      IF (KZ(37).GT.0.AND.LISTV(1).GT.0) THEN
          IF (DMOD(TIME,STEPM).EQ.0.0D0) THEN
              CALL SHRINK(TMIN)
              IF (LISTV(1).GT.0) THEN
                  CALL HIVEL(-1)
              END IF
          END IF
      END IF
*
*       Check optional mass loss time at end of block-step.
      IF (KZ(19).GT.0) THEN
*       Delay until time commensurate with 100-year step (new polynomials).
          IF (TIME.GT.TMDOT.AND.DMOD(TIME,STEPX).EQ.0.0D0) THEN
              IF (KZ(19).GE.3) THEN
                  CALL MDOT
              ELSE
                  CALL MLOSS
              END IF
              IF (IPHASE.LT.0) GO TO 999
          END IF
      END IF
*
*       Include optional plotting of single BH orbits (1 or 2).
      IF (KZ(45).GT.0.AND.TBH.LT.TIME+TOFF) THEN
          CALL BHPLOT
*       Update time interval (try 10 points per time unit).
          TBH = TBH + 1.0D-01
      END IF
*
*       Advance counters and check timer & optional COMMON save (NSUB = 0).
      NTIMER = NTIMER + NXTLEN
      IF (NTIMER.LT.NMAX) GO TO 1

      NTIMER = 0
      NSTEPS = NSTEPS + NMAX
*
      IF (NSTEPS.GE.10*N.AND.N.GT.1000.AND.NSUB.EQ.0) THEN
          NSTEPS = 0
          IF (KZ(1).EQ.3) CALL MYDUMP(1,1)
      END IF
*
*       Check option for general binary search.
      IF (KZ(4).GT.0.AND.TIME - TLASTS.GT.DELTAS) THEN  
          CALL EVOLVE(0,0)
      END IF
*
*       Include facility for termination of run (create dummy file STOP).
      OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) THEN
          CLOSE (99)
          IF (NSUB.EQ.0)  WRITE (6,70)
   70     FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
*
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      IF (TCOMP.LT.CPU) GO TO 1
*
*       Do not terminate during triple, quad or chain regularization.
      IF (NSUB.GT.0) THEN
*       Specify zero step to enforce termination.
          DO 75 L = 1,NSUB
              STEPS(L) = 0.0D0
   75     CONTINUE
          NTIMER = NMAX
          GO TO 1
      END IF
*
*       Obtain elapsed wall-clock time (hours, minutes & seconds).
      CALL WTIME(IHOUR,IMIN,ISEC)
      SECS = 3600.0*IHOUR + 60.0*IMIN + ISEC
      WTOT = WTOT + SECS - WTOT0
      WTOT0 = SECS
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).GT.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          WT = WTOT/3600.0
          CALL MYDUMP(1,1)
          WRITE (6,80)  TIME+TOFF, TCOMP, CPUTOT/60.0, ERRTOT, DETOT, WT
   80     FORMAT (/,9X,'COMMON SAVED AT TIME =',F8.2,'  TCOMP =',F7.1,
     &                 '  CPUTOT =',F6.1,'  ERRTOT =',F10.6,
     &                 '  DETOT =',F10.6,'  WTOT =',F7.1)
      END IF
*
*       Close GPU & GPUIRR library and stop.
      CALL GPUNB_CLOSE
      CALL GPUIRR_CLOSE
      STOP
*
*       Set current global time.
  100 TTOT = TIME + TOFF
*
      RETURN
*
      END
