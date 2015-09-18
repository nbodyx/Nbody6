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
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
*
*
*       Determine closest approaching perturber (JCLOSE for ABSORB).
      RX = 1.0
      JCLOSE = 0
      NP = LISTC(1)
      DO 10 L = 2,NP+1
          J = LISTC(L)
          RIJ2 = 0.0
          DO 5 K = 1,3
              RIJ2 = RIJ2 + (X(K,J) - X(K,ICH))**2
    5     CONTINUE
          IF (RIJ2.LT.RX) THEN
              RD = 0.0
              DO 6 K = 1,3
                  RD = RD + (X(K,J)-X(K,ICH))*(XDOT(K,J)-XDOT(K,ICH))
    6         CONTINUE
              IF (RD.LT.0.0) THEN
                  JCLOSE = J
                  RX = RIJ2
              END IF
          END IF
   10 CONTINUE
*
*       Impose realistic distance limit (closest perturber may be receding).
      IF (JCLOSE.EQ.0.OR.(RX.GT.RMIN22.AND.GPERT.LT.0.2)) GO TO 50
      IF (RX.GT.4.0*RMIN22) GO TO 50
      WRITE (6,20)  TIME+TOFF, NSTEP1, NN, NPERT, NAME(JCLOSE),
     &              SQRT(RX), GPERT, 1.0/RINV(1)
   20 FORMAT (/,' CHAIN INJECT    T # NN NP NM RX GP RB ',
     &                            F9.3,I6,2I4,I8,1P,3E10.2)
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
   50 RETURN
*
      END
