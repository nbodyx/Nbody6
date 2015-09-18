      SUBROUTINE KSPINIT
*     
*     
*     Perturber lists for primordial binaries.
*     ----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  UI(4),VI(4)
*     
*
      TIME0 = TIME
!$omp parallel do 
!$omp& private(ipair,icm,i1,i2,eb,semi,tk,vi,ui,tp,imod,K)
      DO 50 IPAIR = 1,NPAIRS
         ICM = IPAIR + N
         I2 = 2*IPAIR
         I1 = I2 - 1
*
*       Form perturber list for each pair.
         CALL KSLIST(IPAIR)
*     
*       Transform unperturbed hard binary to apocentre and set time-step.
         EB = H(IPAIR)*BODY(I1)*BODY(I2)/BODY(ICM)
         SEMI = -0.5*BODY(ICM)/H(IPAIR)
         IF (LIST(1,I1).EQ.0.AND.EB.LT.EBH) THEN  ! EBH set in zero.f.
            TK = TWOPI*SEMI*SQRT(SEMI/BODY(ICM))
*       Note TIME is not commensurate after KSPERI (cf. CHTERM & STEPS).
            DO 30 K = 1,4
               UI(K) = U(K,IPAIR)
               VI(K) = UDOT(K,IPAIR)
   30       CONTINUE
*       Determine pericentre time (TP < 0 if TDOT2 < 0) and add TK/2.
            CALL TPERI(SEMI,UI,VI,BODY(ICM),TP)
*       Note: apocentre to apocentre gives almost zero step.
            STEP(I1) = 0.5*MIN(TK,STEP(ICM)) - TP
*       Transform KS variables to peri and by pi/2 to apocentre (skip apo).
            IF (ABS(TDOT2(IPAIR)).GT.1.0E-12.OR.R(IPAIR).LT.SEMI) THEN
               TIME = TIME0
               CALL KSPERI(IPAIR)
               CALL KSAPO(IPAIR)
               TIME = TIME0
               TIME = MAX(TIME,0.0D0)
*       Reset TIME to quantized value (small > 0 or < 0 possible initially).
            ELSE IF (TDOT2(IPAIR).GT.0.0) THEN
               TDOT2(IPAIR) = -1.0E-20
            END IF
         END IF
*     
*       Estimate an appropriate KS slow-down index for G < GMIN.
         IMOD = 1
         IF (LIST(1,I1).EQ.0.AND.SEMI.GT.0.0) THEN
            TK = TWOPI*SEMI*SQRT(SEMI/BODY(ICM))
            IF (KZ(26).GT.0.AND.STEP(ICM).GT.TK) THEN
               IMOD = 1 + INT(LOG(STEP(ICM)/TK)/0.69)
               IMOD = MIN(IMOD,5)
            END IF
         END IF
*
*       Obtain polynomials for perturbed KS motion (standard case & merger).
         CALL KSPOLY(IPAIR,IMOD)
*
*       Set primordial binary indicator in second component (needs updating).
         LIST(2,I2) = -1
*       Print diagnostics.
c            WRITE (6,40) IPAIR, NAME(ICM), NAME(I1),NAME(I2),STEP(ICM),
c     &            STEP(I1),LIST(1,ICM),LIST(1,I1),GAMMA(IPAIR),IMOD
c  40        FORMAT ('KS PAIR: IPAIR',I10,' NAME(ICM,I1,I2)',3I10,
c     &           ' STEP(ICM,I1)',1P,2E10.2,0P,' NNB',I4,' NP',I4,
c     &           ' GAMMA',1P,E10.2,0P,' IMOD',I3)
   50 CONTINUE    
!$omp end parallel do      
*
      RETURN
*
      END
