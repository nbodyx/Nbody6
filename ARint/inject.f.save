      SUBROUTINE INJECT(ISUB)
*
*
*       Injection of chain member.
*       --------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX3=3*NMX,NMX4=4*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC
      COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &           XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &           XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
*
*
*       Determine closest approaching neighbour (JCLOSE for ABSORB).
      RX = 1.0
      DO 10 I = IFIRST,N
          IF (I.EQ.ICH) GO TO 10
          RIJ2 = 0.0
          DO 5 K = 1,3
              RIJ2 = RIJ2 + (X(K,I) - X(K,ICH))**2
    5     CONTINUE
          IF (RIJ2.LT.RX) THEN
              RD = 0.0
              DO 6 K = 1,3
                  RD = RD + (X(K,I)-X(K,ICH))*(XDOT(K,I)-XDOT(K,ICH))
    6         CONTINUE
              IF (RD.LT.0.0) THEN
                  JCLOSE = I
                  RX = RIJ2
              END IF
          END IF
   10 CONTINUE
*
*       Save KS component for termination.
      JCOMP = JCLOSE
      WRITE (6,20)  TIME+TOFF, NN, NAME(JCLOSE), SQRT(RX), 1.0/RINV(1)
   20 FORMAT (/,' CHAIN INJECT    T NN NM RX RBH ',F8.3,I4,I8,1P,2E10.2)
      CALL FLUSH(6)
*
*     NNB = LISTC(1) + 2
*     LISTC(NNB) = ICH
*     DO 30 L = 2,NNB
*     J = LISTC(L)
*     WRITE (6,25) J, NAME(J),BODY(J), STEP(J), (X(K,J),K=1,3)
*  25 FORMAT (' LISTC   J NM S M X  ',2I7,1P,2E10.2,2X,3E10.2)
*  30 CONTINUE
*
*       Increase chain membership in the usual way.
      MASS = MASS + BODY(JCLOSE)
      CALL ABSORB(ISUB)
*
*       Update perturber list.
      CALL CHLIST(ICH)
*
      RETURN
*
      END
