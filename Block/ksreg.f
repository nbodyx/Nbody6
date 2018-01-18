      SUBROUTINE KSREG
*
*
*       New KS regularization.
*       ----------------------
*
      INCLUDE 'common6.h'
      COMMON/FSAVE/  SAVEIT(6)
      REAL*8  SAVE(15)
      EXTERNAL RENAME
*
*
*       Ensure ICOMP < JCOMP.
      IF (ICOMP.GT.JCOMP) THEN
          ISAVE = ICOMP
          ICOMP = JCOMP
          JCOMP = ISAVE
      END IF
*
*       Save updated regular force & derivative for new KSINIT procedure.
      DTR = TIME - T0R(ICOMP)
      DO K = 1,3
          SAVEIT(K) = FR(K,ICOMP) + FRDOT(K,ICOMP)*DTR
          SAVEIT(K+3) = D1R(K,ICOMP) + D2R(K,ICOMP)*DTR
      END DO
*
*       Replace #JCOMP by arbitrary body in case it is the only neighbour.
      NNB = LIST(1,ICOMP)
      DO L = 2,NNB+1
          J = LIST(L,ICOMP)
          IF (J.NE.JCOMP.AND.BODY(J).GT.0d0) THEN
              GO TO 101
          ENDIF
      END DO
*
      DO J = IFIRST+2,NTOT
          IF (J.NE.ICOMP.AND.J.NE.JCOMP.AND.BODY(J).GT.0d0) THEN
              LIST(1,ICOMP) = LIST(1,ICOMP) + 1
              LIST(2,ICOMP) = J
              GO TO 101
          END IF
      END DO
  101 CONTINUE
*
*       Copy neighbour list of #ICOMP without JCOMP.
      NNB1 = 1
      DO 1 L = 1,NNB
          IF (LIST(L+1,ICOMP).EQ.JCOMP) GO TO 1
          NNB1 = NNB1 + 1
          ILIST(NNB1) = LIST(L+1,ICOMP)
    1 CONTINUE
      ILIST(1) = NNB1 - 1
*
*       Save basic variables for components unless in correct location.
      DO 10 KCOMP = 1,2
*       Treat the first & second component in turn.
          IF (KCOMP.EQ.1) THEN
              I = ICOMP
          ELSE
              I = JCOMP
          END IF
          J = 2*NPAIRS + KCOMP
          IF (I.EQ.J) GO TO 10
*
          DO 2 K = 1,3
              SAVE(K) = X(K,I)
              SAVE(K+3) = X0DOT(K,I)
    2     CONTINUE
          SAVE(7) = BODY(I)
          SAVE(8) = RS(I)
          SAVE(9) = RADIUS(I)
          SAVE(10) = TEV(I)
          SAVE(11) = BODY0(I)
          SAVE(12) = TEV0(I)
          SAVE(13) = EPOCH(I)
          SAVE(14) = SPIN(I)
          SAVE(15) = ZLMSTY(I)
          NAMEI = NAME(I)
          KSI = KSTAR(I)
*
*       Exchange first & second single particle with ICOMP & JCOMP.
          DO 4 K = 1,3
              X(K,I) = X(K,J)
              X0(K,I) = X0(K,J)
              X0DOT(K,I) = X0DOT(K,J)
              XDOT(K,I) = XDOT(K,J)
              F(K,I) = F(K,J)
              FDOT(K,I) = FDOT(K,J)
              FI(K,I) = FI(K,J)
              FIDOT(K,I) = FIDOT(K,J)
              D0(K,I) = D0(K,J)
              D1(K,I) = D1(K,J)
              D2(K,I) = D2(K,J)
              D3(K,I) = D3(K,J)
              FR(K,I) = FR(K,J)
              FRDOT(K,I) = FRDOT(K,J)
              D0R(K,I) = D0R(K,J)
              D1R(K,I) = D1R(K,J)
              D2R(K,I) = D2R(K,J)
              D3R(K,I) = D3R(K,J)
              X(K,J) = SAVE(K)
              X0DOT(K,J) = SAVE(K+3)
    4     CONTINUE
*
          BODY(I) = BODY(J)
          RS(I) = RS(J)
          RADIUS(I) = RADIUS(J)
          TEV(I) = TEV(J)
          TEV0(I) = TEV0(J)
          BODY0(I) = BODY0(J)
          EPOCH(I) = EPOCH(J)
          SPIN(I) = SPIN(J)
          ZLMSTY(I) = ZLMSTY(J)
          NAME(I) = NAME(J)
          KSTAR(I) = KSTAR(J)
          STEP(I) = STEP(J)
          STEPR(I) = STEPR(J)
          T0(I) = T0(J)
          T0R(I) = T0R(J)
          K = LIST(1,J) + 1
          DO 5 L = 1,K
              LIST(L,I) = LIST(L,J)
    5     CONTINUE
          BODY(J) = SAVE(7)
          RS(J) = SAVE(8)
          RADIUS(J) = SAVE(9)
          TEV(J) = SAVE(10)
          BODY0(J) = SAVE(11)
          TEV0(J) = SAVE(12)
          EPOCH(J) = SAVE(13)
          SPIN(J) = SAVE(14)
          ZLMSTY(J) = SAVE(15)
          NAME(J) = NAMEI
          KSTAR(J) = KSI
   10 CONTINUE
*
*       Save neighbour list of first component for RENAME & FPOLY.
      DO 15 L = 1,NNB1
          LIST(L,IFIRST) = ILIST(L)
   15 CONTINUE
*
*       Increase pair index, total number & single particle index.
      NPAIRS = NPAIRS + 1
      NTOT = N + NPAIRS
      IFIRST = 2*NPAIRS + 1
*
*       Update all relevant COMMON list arrays.
      CALL RENAME
*
*       Check replacing of single KS component by corresponding c.m.
      DO 30 I = IFIRST-2,NTOT-1
   20     IF (LIST(1,I).GT.0.AND.LIST(2,I).LT.IFIRST) THEN
              J2 = LIST(2,I)
              J = KVEC(J2) + N
              IF (LIST(1,I).EQ.1) THEN
                  LIST(2,I) = J
              ELSE
                  L = 2
   22             JNEXT = LIST(L+1,I)
                  IF (JNEXT.LT.J) THEN
                      LIST(L,I) = JNEXT
                      L = L + 1
                      IF (L.LE.LIST(1,I)) GO TO 22
                      LIST(L,I) = J
                  ELSE IF (JNEXT.EQ.J) THEN
                      DO 25 LL = L,LIST(1,I)
                          LIST(LL,I) = LIST(LL+1,I)
   25                 CONTINUE
                      LIST(1,I) = LIST(1,I) - 1
                  ELSE
                      LIST(L,I) = J
                  END IF
*       Check again until first neighbour > ICOMP.
                  GO TO 20
              END IF
          END IF
   30 CONTINUE
*
*       Copy neighbour list for second component & c.m. (NNB1 = LIST(1,I)+1).
      I = 2*NPAIRS - 1
      DO 40 L = 1,NNB1
          LIST(L,I+1) = LIST(L,I)
          LIST(L,NTOT) = LIST(L,I)
   40 CONTINUE
*
*       Initialize the regularized solution.
      CALL KSINIT
*
*       Check optional binary analysis after merger or multiple collision.
*     IF (KZ(4).GT.0.AND.IPHASE.GT.3) THEN
*         CALL EVOLVE(NPAIRS,-1)
*     END IF
*
*       Check updating of global index for chain c.m.
      IF (NCH.GT.0) THEN
          CALL CHFIND
      END IF
*
      RETURN
*
      END
