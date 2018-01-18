      SUBROUTINE CHLIST(II)
*
*
*       Perturber list for ARC regularization.
*       --------------------------------------
*
      INCLUDE 'common6.h'
        REAL*8  M,MASS,MC
        PARAMETER (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &  NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
         COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &   XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      SAVE ITER
      DATA ITER /0/
*
*
*     IF (II.EQ.0) I = ICH
      RPERT = 0.5*RSUM
      RCRIT2 = 2.0*RPERT**2/BODY(ICH)
      RCRIT3 = RCRIT2*RPERT/(0.03*GMIN)
      RCRIT2 = CMSEP2*RPERT**2
      GPERT = 0.0
*
*       Select new perturbers from c.m. neighbour list (search to 100*RPERT).
    1 NNB1 = 1
      NNB2 = LIST(1,ICH) + 1
      DO 10 L = 2,NNB2
          J = LIST(L,ICH)
*       Note call from routine CHAIN with II = 0 (hence use ICH).
          W1 = X(1,J) - X(1,ICH)
          W2 = X(2,J) - X(2,ICH)
          W3 = X(3,J) - X(3,ICH)
          RSEP2 = W1*W1 + W2*W2 + W3*W3
*       Employ fast test followed by actual condition.
          IF (RSEP2.GT.RCRIT2) GO TO 10
          RIJ3 = RSEP2*SQRT(RSEP2)
          IF (RIJ3.GT.BODY(J)*RCRIT3) GO TO 10
          GAM = 0.0
*       Loop over chain members to get maximum for each neighbour.
          DO 5 LL = 1,NCH-1
*       Estimate perturbation for each chain distance using same RIJ3.
              GX = 2.0*BODY(J)/(M(LL) + M(LL+1))/(RIJ3*RINV(LL)**3)
              GAM = MAX(GX,GAM)
    5     CONTINUE
*       Form perturber list and set maximum relative perturbation > 10^{-7}.
          IF (GAM.GT.1.0D-07) THEN
              NNB1 = NNB1 + 1
              LISTC(NNB1) = J
              GPERT = MAX(GAM,GPERT)
          END IF
   10 CONTINUE
*
*       Limit the perturber number to 20 using volume factor.
      IF (NNB1.GT.20) THEN
          RCRIT3 = (20.0/FLOAT(NNB1))*RCRIT3
          GO TO 1
      END IF
          
*       Save perturber membership.
      LISTC(1) = NNB1 - 1
      NPERT = LISTC(1)
*
*       Ensure consistent perturbation indicator.
      IF (GPERT.GT.1.0D-07) THEN
          IPERT = 1
      ELSE
          IPERT = 0
      END IF
*
      IF (NPERT.GT.6.AND.ITER.LT.10) THEN
        ITER = ITER + 1
        WRITE (6,20)  NPERT, NNB2, GPERT, RPERT, (1.0/RINV(K),K=1,NCH-1)
   20   FORMAT (' CHLIST!   NP NB2 GP RPERT R  ',I4,I5,1P,8E10.2)
      END IF
*
      RETURN
*
      END
