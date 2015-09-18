      SUBROUTINE CMFIRR2(I,NNB,KLIST,XI,XID,FIRR,FD)
*
*
*       Irregular force correction from J > N (singles or unpert KS).
*       -------------------------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/PREDICT/ TPRED(NMAX)
      REAL*8  XI(3),XID(3),FIRR(3),FD(3),DX(3),DV(3),
     &        FP(6),FPD(6),FCM(3),FCMD(3),CMX(3),CMV(3),XK(6),VK(6)
      INTEGER KLIST(LMAX)
*
*
*       Set neighbour number & list index of the last single particle.
      NNB1 = NNB + 1
      NNB2 = NNB1
   25 IF (KLIST(NNB2).LE.N) GO TO 30
      NNB2 = NNB2 - 1
      IF (NNB2.GT.1) GO TO 25
*       Include special case of only c.m. neighbours.
      GO TO 40
*
*       Exit if no binaries on neighbour list.
   30 IF (NNB2.EQ.NNB1) GO TO 60

*       Perform differential correction for regularized c.m. neighbours.
   40 DO 50 LL = NNB2+1,NNB1
          J = KLIST(LL)
*       See whether to sum over binary components.
          JPAIR = J - N
          J1 = 2*JPAIR - 1
*       Skip unperturbed binary (treated as single particle).
          IF (LIST(1,J1).EQ.0) THEN
              GO TO 50
          END IF
*
*       Obtain coordinates & velocity by prediction or copying.
          IPRED = 1
          IF (TPRED(J).EQ.TIME) IPRED = 0
          CALL KSRES3(JPAIR,J1,J2,IPRED,CMX,CMV,XK,VK)
*
*       Obtain individual c.m. force with single particle approximation.
          A1 = CMX(1) - X(1,I)
          A2 = CMX(2) - X(2,I)
          A3 = CMX(3) - X(3,I)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
          DV(1) = CMV(1) - XDOT(1,I)
          DV(2) = CMV(2) - XDOT(2,I)
          DV(3) = CMV(3) - XDOT(3,I)
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FCM(1) = A1*DR3I
          FCM(2) = A2*DR3I
          FCM(3) = A3*DR3I
          FCMD(1) = (DV(1) - A1*DRDV)*DR3I
          FCMD(2) = (DV(2) - A2*DRDV)*DR3I
          FCMD(3) = (DV(3) - A3*DRDV)*DR3I
*
*       Evaluate perturbation on first component of body #J.
          dr2 = 0.0
          drdv = 0.0
          DO 42 L = 1,3
              dx(L) = XK(L) - X(L,I)
              dv(L) = VK(L) - XDOT(L,I)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
   42     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(J1)*dr2i*SQRT(dr2i)
          drdv = 3.0*drdv*dr2i
*
          DO 45 L = 1,3
              FP(L) = dx(L)*dr3i
              FPD(L) = (dv(L) - dx(L)*drdv)*dr3i
   45     CONTINUE
*
*       Evaluate perturbation on second component of body #J.
          dr2 = 0.0
          drdv = 0.0
          DO 46 L = 1,3
              dx(L) = XK(L+3) - X(L,I)
              dv(L) = VK(L+3) - XDOT(L,I)
              dr2 = dr2 + dx(L)**2
              drdv = drdv + dx(L)*dv(L)
   46     CONTINUE
*
          dr2i = 1.0/dr2
          dr3i = BODY(J2)*dr2i*SQRT(dr2i)
          drdv = 3.0*drdv*dr2i
*
          DO 48 L = 1,3
              FP(L+3) = dx(L)*dr3i
              FPD(L+3) = (dv(L) - dx(L)*drdv)*dr3i
   48     CONTINUE
*
*       Accumulate individual correction terms together.
          DO 49 L = 1,3
              FIRR(L) = FIRR(L) + (FP(L) + FP(L+3) - FCM(L))
              FD(L) = FD(L) + (FPD(L) + FPD(L+3) - FCMD(L))
   49     CONTINUE
   50 CONTINUE
*
   60 RETURN
*
      END
