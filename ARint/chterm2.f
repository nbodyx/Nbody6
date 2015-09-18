      SUBROUTINE CHTERM2(NBH2)
*
*
*       Termination of two-body ARC.
*       ----------------------------
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
      COMMON/BINARY/  ZZ(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
      COMMON/ECHAIN/  ECH
*
*
*       Define indices for two components and set commensurate time.
      I1 = 1
      I2 = 2
      JLIST(6) = NAMEC(I1)
      JLIST(7) = NAMEC(I2)
      NAME(ICH) = NAME0
      TIME = TBLOCK
*
*       Identify c.m. body and find global indices.
      ICM = 0
      DO 10 J = IFIRST,N
          DO 5 L = 1,NN
              IF (NAME(J).EQ.JLIST(L+5)) THEN
                  JLIST(L) = J
                  IF (BODY(J).GT.0.0D0) ICM = J
              END IF
    5     CONTINUE
   10 CONTINUE
      IF (ICM.EQ.0) ICM = ICH
*
*       Ensure ICOMP < JCOMP for KS regularization.
      ICOMP = MIN(JLIST(1),JLIST(2))
      JCOMP = MAX(JLIST(1),JLIST(2))
*
*       Copy final coordinates & velocities to standard variables.
      LK = 0
      DO 20 L = 1,NN     ! Note NCH may be zero after INFALL.
          DO 15 K = 1,3
              LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
   15     CONTINUE
   20 CONTINUE
*
*       Predict current coordinates & velocities to F3DOT before termination.
*     CALL XVPRED(ICM,-1)
*
*       Copy c.m. coordinates & velocities.
      DO 25 K = 1,3
          CM(K) = X(K,ICM)
          CM(K+3) = XDOT(K,ICM)
   25 CONTINUE
*
*       Set configuration pointers for KS candidates.
      JLIST(7) = I1
      JLIST(8) = I2
*
*       Place new global coordinates in the original locations.
      DO 30 L = 1,NN
          J = JLIST(L)
*       Copy the respective masses (BODY(ICM) holds the sum).
          IF (L.EQ.1) BODY(J) = BODYC(1)
          IF (L.EQ.2) BODY(J) = BODYC(2)
*       Transform to global coordinates & velocities using c.m. values.
          LL = JLIST(L+6)
          DO 28 K = 1,3
              X(K,J) = X4(K,LL) + CM(K)
              XDOT(K,J) = XDOT4(K,LL) + CM(K+3)
              X0(K,J) = X(K,J)
              X0DOT(K,J) = XDOT(K,J)
   28     CONTINUE
   30 CONTINUE
*
*       Ensure new neighbour list (may be zero for ICOMP).
      RS0 = RS(ICM)
      CALL NBLIST(ICOMP,RS0)
*
*       Save NAME of dominant perturbers for neighbour list check.
      NAM1 = NAME(ICOMP)
      NAM2 = NAME(JCOMP)
*
*       Perform KS regularization of dominant components.
      CALL KSREG
*
*       Search for missing dominant bodies in perturber list.
      J1 = 2*NPAIRS - 1
      NP = LIST(1,J1)
      DO 40 L = 2,NP+1
          J = LIST(L,J1)
          NB1 = LIST(1,J) + 1
          I0 = 0
*       Check whether the former IC1 & IC2 are members.
          DO 33 LL = 2,NB1
              JJ = LIST(LL,J)
              IF (NAME(JJ).EQ.NAM1.OR.NAME(JJ).EQ.NAM2) I0 = I0 + 1
   33     CONTINUE
*
*       Add new c.m. at the end after failed search (NBREST not used).
          IF (I0.EQ.0.AND.LIST(NB1,J).NE.NTOT) THEN
              LIST(NB1+1,J) = NTOT
              LIST(1,J) = LIST(1,J) + 1
          END IF
   40 CONTINUE
*
*       Include optional kick velocity of 3*VRMS km/s after GR coalescence.
      IF (KZ(43).GT.0.AND.NBH2.EQ.2) THEN
          VI20 = 0.0
          DO 42 K = 1,3
              VI20 = VI20 + XDOT(K,NTOT)**2
   42     CONTINUE
          VF = 3.0*(VRMS/VSTAR)/SQRT(VI20)
*       Note c.m. is assigned kick because both members would escape.
          DO 44 K = 1,3
              XDOT(K,NTOT) = VF*XDOT(K,NTOT)
              X0DOT(K,NTOT) = XDOT(K,NTOT)
   44     CONTINUE
          ECD0 = ECDOT
          VESC = 3.0*VRMS
          ECDOT = ECDOT + 0.5*BODY(NTOT)*VI20*(1.0 - VF**2)
          WRITE (6,45)  LIST(1,2*NPAIRS-1), VF, ECD0-ECDOT, VESC
   45     FORMAT (' COALESCENCE KICK    NP VF ECDOT VESC ',
     &                                  I4,F7.3,F10.6,F6.1)
*       Define BH index for c.m. diagnostic escape output.
          KSTAR(NTOT) = 14
      END IF
*
      IF (NSTEP1.GT.100.OR.NBH2.EQ.2) THEN
          NP = LIST(1,2*NPAIRS-1)
          ZMU = BODY(2*NPAIRS-1)*BODY(2*NPAIRS)/BODY(NTOT)
          EB = ZMU*H(NPAIRS)
          WRITE (6,46)  NSTEP1, NP, LIST(1,NTOT), EB, ECH, H(NPAIRS),
     &                  R(NPAIRS), STEP(NTOT)
   46     FORMAT (' TERMINATE ARC    # NP NNB EB ECH H R STEP ',
     &                                 I7,2I4,1P,5E10.2)
      END IF
*
*       Reduce subsystem counter and initialize membership & internal energy.
      NSUB = MAX(NSUB - 1,0)
      NCH = 0
      NN = 0
      ECH = 0.0
*       Note stellar collision energies are accumulated in ECOLL by DECORR.
      NSTEPC = NSTEPC + NSTEP1
*
*       Find reference pointer for ARC (should have ISYS(L)=4).
      ISUB = 1
      DO 47 L = 1,NSUB+1
          IF (ISYS(L).EQ.4) ISUB = L
   47 CONTINUE
*
*       Update subsystem pointer and relevant variables (ISUB <= NSUB).
      DO 50 L = ISUB,NSUB
          DO 48 K = 1,6
              BODYS(K,L) = BODYS(K,L+1)
              NAMES(K,L) = NAMES(K,L+1)
   48     CONTINUE
          T0S(L) = T0S(L+1)
          TS(L) = TS(L+1)
          STEPS(L) = STEPS(L+1)
          RMAXS(L) = RMAXS(L+1)
          ISYS(L) = ISYS(L+1)
   50 CONTINUE
*
*       Set escape condition for ghost after termination/coalescence.
      DO 60 I = IFIRST,N
          IF ((NAME(I).EQ.NAMEC(I1).AND.BODY(I).EQ.0.0D0).OR.
     &        (NAME(I).EQ.NAMEC(I2).AND.BODY(I).EQ.0.0D0)) THEN
*       Skip any ghosts associated with stable hierarchies.
              DO 52 L = 1,NMERGE
                  IF (NAME(I).EQ.NAMEG(L)) GO TO 60
   52         CONTINUE
*       Ensure that ghost will escape next output (far from fast escapers).
              DO 55 K = 1,3
                  X0(K,I) = 1.0d+04 + 0.001*FLOAT(K*I)
                  X(K,I) = X0(K,I)
                  X0DOT(K,I) = SQRT(0.04d0*ZMASS/RSCALE)+FLOAT(K*I)/100.
                  XDOT(K,I) = X0DOT(K,I2)
   55         CONTINUE
*       Note that T0 may be set large for ghost (F & D2 checked as zero).
              T0(I) = TIME + DTADJ
          END IF
   60 CONTINUE
*
      RETURN
*
      END
