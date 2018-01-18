      SUBROUTINE CHLIST(I)
*
*
*       Perturber list for chain regularization.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  M,MASS,MC,MIJ,MKK
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      SAVE ITER
      DATA ITER /0/
*
*
*       Replace zero argument from CHAIN by #ICH.
      IF (I.EQ.0) I = ICH
*       Use 5 x harmonic mean of RMIN & RGRAV for basic search distance.
*     RPERT = 5.0*RMIN*ABS(RGRAV)/(RMIN + ABS(RGRAV))
      RPERT = 0.5*RSUM      ! adopted 01/16.
      RPERT = MIN(RPERT,2.0*RMIN)
      RCRIT2 = 2.0*RPERT**2/BODY(ICH)
      RCRIT3 = RCRIT2*RPERT/(0.01*GMIN)
      RCRIT2 = CMSEP2*RPERT**2  ! maybe a bit much
      GPERT = 0.0
*
*       Select new perturbers from c.m. neighbour list (search to 100*RPERT).
    1 NNB1 = 1
      NNB2 = LIST(1,ICH) + 1
      DO 10 L = 2,NNB2
          J = LIST(L,ICH)
*       Note call from routine CHAIN with I = 1 (hence use ICH).
          W1 = X(1,J) - X(1,ICH)
          W2 = X(2,J) - X(2,ICH)
          W3 = X(3,J) - X(3,ICH)
          RSEP2 = W1*W1 + W2*W2 + W3*W3
*       Employ fast test followed by actual condition.
          IF (RSEP2.GT.RCRIT2) GO TO 10
          RIJ3 = RSEP2*SQRT(RSEP2)
          IF (RIJ3.GT.BODY(J)*RCRIT3) GO TO 10
          GAM = 0.0
*       Loop over chain members for each neighbour.
          DO 5 LL = 1,NCH-1
*       Estimate perturbation for each chain distance using same RIJ3.
              GX = 2.0*BODY(J)/(M(LL) + M(LL+1))/(RIJ3*RINV(LL)**3)
              GAM = MAX(GX,GAM)
    5     CONTINUE
*       Form perturber list and set maximum relative perturbation > 10^{-8}.
          IF (GAM.GT.1.0D-08) THEN
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
      IF (GPERT.GT.1.0D-08) THEN
          IPERT = 1
      ELSE
          IPERT = 0
      END IF
*
      IF (NPERT.GT.6.AND.ITER.LT.10) THEN
        ITER = ITER + 1
        WRITE (6,20)  NPERT, NNB2, GPERT, RPERT, (1.0/RINV(K),K=1,NCH-1)
   20   FORMAT (' CHLIST!   NP NB2 GP RPERT R  ',I4,I5,1P,6E10.2)
      END IF
*
      RETURN
*
      END
