      SUBROUTINE FPOLY0(RS0)
*
*
*       GPU initialization of forces and neighbour lists.
*       -------------------------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER (NIMAX=1024)
      REAL*8   H2I(NIMAX),GPUACC(3,NIMAX),GPUJRK(3,NIMAX),GPUPHI(NIMAX),
     &         GF(3,NIMAX),GFD(3,NIMAX),DTR(NIMAX)
      INTEGER  LISTGP(LMAX,NIMAX),IREG(NMAX)
      SAVE  GF,GFD,LISTGP,IREG
*
*
*       Open the GPU libraries on new run (note nnbmax = NN is printed).
      NN = NTOT + 100
      CALL GPUNB_OPEN(NN)
*       Set larger value for GPUIRR (note further possible increase of NTOT).
      NNN = NN + 10
      CALL GPUIRR_OPEN(NNN,LMAX)
*
*       Set provisional neighbour radius and regular time-step for GPUNB.
*!$omp parallel do private(I, NFI)
      NFI = 0
*       Note use of NFI & parallel is risky (standard loop is cheap).
      DO 10 I = IFIRST,NTOT
          NFI = NFI + 1
          IREG(NFI) = I
          T0(I) = 0.0
*       Modify neighbour radius according to optional NBLIST procedure.
          RI2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
          RS(I) = RS0*SQRT(1.0 + RI2)
*       Set an estimated r-dependent regular time-step (not quantized).
          STEPR(I) = SMAX/8.0D0*SQRT(1.0 + RI2)
          STEPR(I) = MIN(STEPR(I),SMAX)
*       Initialize F & FDOT.
          DO 5 K = 1,3
              F(K,I) = 0.0
              FDOT(K,I) = 0.0
    5     CONTINUE
   10 CONTINUE
*!$omp end parallel do
*
*       Send all particles (X & XDOT) and zero F & FDOT to GPU.
!$omp parallel do private(I)
      DO 15 I = IFIRST,NTOT
          CALL GPUIRR_SET_JP(I,X(1,I),XDOT(1,I),F(1,I),FDOT(1,I),
     &                                          BODY(I),T0(I))
   15 CONTINUE
!$omp end parallel do
*
      CALL GPUIRR_PRED_ALL(IFIRST,NTOT,TIME)  ! may not be needed.
*
*       Send all single particles and c.m. bodies (NFI) to GPU.
      CALL GPUNB_SEND(NFI,BODY(IFIRST),X(1,IFIRST),XDOT(1,IFIRST))
*
*       Define maximum GPU neighbour number and initialize counters.
      NBMAX = MIN(NNBMAX + 150,LMAX-5)
      JNEXT = 0
      NOFL2 = 0
*
*       Loop over all particles split into NIMAX blocks.
      DO 100 II = IFIRST,NTOT,NIMAX
   20     NI = MIN(NFI-JNEXT,NIMAX)
*       Copy neighbour radius, STEPR and state vector for each block.
!$omp parallel do private(LL, JNEXT, I, K)
          DO 30 LL = 1,NI
              I = IREG(JNEXT+LL)
              H2I(LL) = RS(I)**2/BODYM
              DTR(LL) = STEPR(I)
   30     CONTINUE
!$omp end parallel do
*
*       Evaluate forces, first derivatives and neighbour lists for new block.
          I = IREG(JNEXT + 1)
          CALL GPUNB_REGF(NI,H2I,DTR,X(1,I),XDOT(1,I),GPUACC,GPUJRK,
     &                                     GPUPHI,LMAX,NBMAX,LISTGP)
*       Check neighbour lists for overflow or zero membership (NNB = 1).
          DO 50 LL = 1,NI
              NNB = LISTGP(1,LL)
*       Repeat last block with reduced RS(I) on NNB < 2 (at end of loop).
              IF (NNB.LT.2) THEN
                  I = IREG(JNEXT + LL)
                  RI2 = (X(1,I)-RDENS(1))**2 + (X(2,I)-RDENS(2))**2 +
     &                                         (X(3,I)-RDENS(3))**2
                  WRITE (41,40)  NSTEPR, NAME(I), NNB,
     &                           RS(I), SQRT(RI2)
   40             FORMAT (' OVERFLOW!   #R NAME NB RS ri ',
     &                                  I11,I7,I5,2F8.3)
                  CALL FLUSH(41)
                  IF (NNB.LT.0) THEN
                      RS(I) = 0.9*RS(I)
                  ELSE
                      RS(I) = 1.5*RS(I)
                  END IF
                  H2I(LL) = RS(I)**2/BODYM
                  NOFL2 = NOFL2 + 1
              END IF
   50     CONTINUE
*
*       Repeat the last block for rare exceptions (NNB < 0 and = 1).
          IF (NOFL2.GT.0) THEN
              NOFL2 = 0
              GO TO 20
           END IF
*
*      Copy regular force (FR & D1R) and neighbour list from GPU.
!$omp parallel do private(LL, I, ITEMP, NNB, L1, L)
          DO 70 LL = 1,NI
              I = IREG(JNEXT + LL)
              DO 55 K = 1,3
                  FR(K,I) = GPUACC(K,LL)
                  D1R(K,I) = GPUJRK(K,LL)
   55         CONTINUE
              NNB = LISTGP(1,LL)
              L1 = 1
              DO 60 L = 2,NNB+1
*       Note GPU address starts from 0 (hence add IFIRST to neighbour list).
                  ITEMP = LISTGP(L,LL) + IFIRST
                  IF (ITEMP.NE.I) THEN
                      L1 = L1 + 1
                      LISTGP(L1,LL) = ITEMP
                  END IF
   60         CONTINUE
              LISTGP(1,LL) = L1 - 1
              CALL GPUIRR_SET_LIST(I,LISTGP(1,LL))
              DO 65 L = 1,L1
                  LIST(L,I) = LISTGP(L,LL)
   65         CONTINUE
   70     CONTINUE
!$omp end parallel do
*
*       Evaluate irregular forces by vector procedure (IREG(1) is first).
          JJ = JNEXT + 1
          CALL GPUIRR_FIRR_VEC(NI,IREG(JJ),GF(1,1),GFD(1,1))
*
*       Copy irregular force and first derivative.
!$omp parallel do private(I)
          DO 80 LL = 1,NI
              I = IREG(JNEXT + LL)
              DO 75 K = 1,3
                  FI(K,I) = GF(K,LL)
                  D1(K,I) = GFD(K,LL)
   75         CONTINUE
   80     CONTINUE
!$omp end parallel do
*
*       Advance to next block (may not be used).
          JNEXT = JNEXT + NI
  100 CONTINUE
*
*       Check option for external force.
      IF (KZ(14).GT.0) THEN
          CALL XTRNLD(IFIRST,NTOT,1)
      END IF
*
*       Form total force & force derivative and extra variables for XVPRED.
!$omp parallel do private(I)
      DO 110 I = IFIRST,NTOT
          DO 105 K = 1,3
              F(K,I) = FI(K,I) + FR(K,I)
              FDOT(K,I) = D1(K,I) + D1R(K,I)
              D0(K,I) = FI(K,I)
              D0R(K,I) = FR(K,I)
              FIDOT(K,I) = D1(K,I)
              FRDOT(K,I) = D1R(K,I)
  105     CONTINUE
  110 CONTINUE
!$omp end parallel do
*
*       Close the GPU libraries (limits change in INTGRT).
      CALL GPUNB_CLOSE
      CALL GPUIRR_CLOSE
*
      RETURN
*
      END
