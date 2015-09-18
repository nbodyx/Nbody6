      SUBROUTINE CHLIST(II)
*
*
*       Perturber list for chain regularization.
*       -----------------------------------------
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
      COMMON/POSTN2/  SEMIGR,ECCGR,DEGR,ISPIN
      SAVE ITER
      DATA ITER /0/
*
*
*    IF (I.EQ.0) I = ICH
*       Form perturber list based on RMIN, RGRAV and GMIN/50.
      RPERT = MIN(RMIN,SEMIGR*(1.0 + ECCGR))
      RPERT = 6.0*RPERT
      IF (SEMIGR.LT.0.2*RMIN) RPERT = 10.0*SEMIGR
*       Use 5 x harmonic mean of RMIN & RGRAV for basic search distance.
*     RPERT = 5.0*RMIN*ABS(RGRAV)/(RMIN + ABS(RGRAV))
      RCRIT2 = 2.0*RPERT**2/BODY(ICH)
      RCRIT3 = RCRIT2*RPERT/(0.001*GMIN)
      RCRIT2 = CMSEP2*RPERT**2  ! maybe a bit much
      PMAX = 0.0
*
*       Select new perturbers from c.m. neighbour list (search to 100*RMIN).
    5 NNB1 = 1
      NNB2 = LIST(1,ICH) + 1
      DO 10 L = 2,NNB2
          J = LIST(L,ICH)
*       Note call from routine CHAIN with I = 1 (hence use ICH).
          W1 = X(1,J) - X(1,ICH)
          W2 = X(2,J) - X(2,ICH)
          W3 = X(3,J) - X(3,ICH)
          RSEP2 = W1*W1 + W2*W2 + W3*W3
          IF (RSEP2.LT.RCRIT2) THEN
              RIJ3 = RSEP2*SQRT(RSEP2)
*       Estimate perturbed distance from tidal limit approximation.
              IF (RIJ3.LT.BODY(J)*RCRIT3) THEN
                  NNB1 = NNB1 + 1
                  LISTC(NNB1) = J
                  PMAX = MAX(BODY(J)/RIJ3,PMAX)
              END IF
          END IF
   10 CONTINUE
*
*       Limit the perturber number to 20 using volume factor.
      IF (NNB1.GT.20) THEN
          RCRIT3 = (20.0/FLOAT(NNB1))*RCRIT3
          GO TO 5
      END IF
          
*       Save perturber membership.
      LISTC(1) = NNB1 - 1
*
*       Form current perturbation from nearest body (effective value RPERT).
      RP = MIN(RGRAV,RMIN)
      IF (NCH.EQ.2) RP = RGRAV
      GPERT = 2.0*PMAX*RP**3/BODY(ICH)
*
*       Obtain actual perturbation (note GAMX zero first time).
*     CALL CHPERT(GAMX)
      GAMX = 0.0
*
      IF (GAMX.GT.0.0D0) THEN
          GPERT = GAMX
      END IF
      NPERT = LISTC(1)
*
*       Ensure consistent perturbation indicator.
      IF (GPERT.GT.1.0D-09) THEN
          IPERT = 1
      ELSE
          IPERT = 0
      END IF
*
      IF (NPERT.GT.6.AND.ITER.LT.10) THEN
      ITER = ITER + 1
      WRITE (6,20)  NPERT, NNB2, GPERT, RPERT, SEMIGR
   20 FORMAT (' CHLIST!   NP NB2 GP RPERT AGR  ',I4,I5,1P,4E10.2)
      END IF
*
      RETURN
*
      END
