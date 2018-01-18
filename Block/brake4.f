      SUBROUTINE BRAKE4(I1,I2,KCASE,DW)
*
*
*       GR analytical orbit shrinkage.
*       ------------------------------
*
      INCLUDE 'common6.h'
      COMMON/POSTN/  CVEL,TAUGR,RZ1,GAMMAZ,TKOZ,EMAX,TSP,KZ24,IGR,IPN
      COMMON/KSPAR/  ISTAT(KMAX)
      REAL*8 M1,M2,UI(4),UIDOT(4)
      SAVE ITER,HIP,ADOT0,EDOT0
      DATA ITER,HIP,ADOT0,EDOT0 /0,0.0D0,0.0D0,0.0D0/
*
*
*       Check relativistic conditions (at least one >= NS).
      IF (MAX(KSTAR(I1),KSTAR(I2)).LT.13) GO TO 100
*
*       Specify the basic elements from BH or KS treatments.
      M1 = BODY(I1)
      M2 = BODY(I2)
      IPAIR = KVEC(I1)
      I = N + IPAIR
      SEMI = -0.5*BODY(I)/H(IPAIR)
      E2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(E2)
*
*       Replace the BH value of CVEL by CLIGHT for standard binaries.
      IF (CVEL.EQ.0.0D0) CVEL = CLIGHT
*
*       Employ Einstein shift above 1.0D-04 per orbit for decision-making.
      DW = 3.0*TWOPI*(BODY(I))/(SEMI*CVEL**2*(1.0-E2))
      IF (DW.LT.1.0D-04.OR.DW.GT.1.1D-03) GO TO 100   ! Note limit 0.0011.
*
*       Form da/dt & de/dt according to Peters 1964.
      ADOT =  64.0/5.0*M1*M2*(M1+M2)/(CVEL**5*SEMI**3*(1.0 - E2)**3.5)
      ADOT = ADOT*(1.0 + 73.0/24.0*E2 + 37.0/96.0*E2**2)
      EDOT = 304.0/15.0*ECC*M1*M2*(M1 + M2)/(CVEL**5*SEMI**4)
      EDOT = EDOT/(1.0 - E2)**2.5*(1.0 + 121.0/304.0*E2)
      TGR = SEMI/ADOT
*
*       Set local KS vectors.
      DO 2 K = 1,4
          UI(K) = U0(K,IPAIR)
          UIDOT(K) = UDOT(K,IPAIR)
    2 CONTINUE
*
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
      DT1 = STEP(I1)
      THETA = DW*DT1/TK
*       Impose limit if THETA > TWOPI but use full STEP(I1) for change.
      IF (THETA.GT.TWOPI) THEN
          THETA = DMOD(THETA,TWOPI)
      END IF
      DT = DT1
*
*       Rotate KS orbit by THETA/2 (period halving).
      CALL KSROT(UI,UIDOT,THETA)
*
*       Copy back to KS common variables.
      DO 15 K = 1,4
          U(K,IPAIR) = UI(K)
          U0(K,IPAIR) = UI(K)
          UDOT(K,IPAIR) = UIDOT(K)
   15 CONTINUE
*
*       Define RZ and increment the change in SEMI & ECC to second order.
      RZ = 8.0*(M1 + M2)/CVEL**2
      ADOT2 = (ADOT - ADOT0)/DT
      EDOT2 = (EDOT - EDOT0)/DT
      DH = H(IPAIR) - HIP
      IF (ABS(DH).LT.1.0D-04*ABS(HIP)) THEN
          SEMI1 = SEMI - (ADOT + 0.5*ADOT2*DT)*DT
          ECC1 = ECC - (EDOT + 0.5*EDOT2*DT)*DT
      ELSE
          SEMI1 = SEMI
          ECC1 = ECC
      END IF
*
*       Save the derivatives for constructing second order.
      ADOT0 = ADOT
      EDOT0 = EDOT
      IF (SEMI1.LT.0.9*SEMI) THEN
          WRITE (6,20)  NAME(I1), SEMI, SEMI1, ADOT*DT
   20     FORMAT (' PN WARNING    NM A A1 ADOT*DT ',I6,1P,3E10.2)
      END IF
*
*       Update binding energy and collision energy (zero on first call).
      HI = H(IPAIR)
      H(IPAIR) = -0.5*BODY(I)/SEMI1
      ZMU = BODY(I1)*BODY(I2)/BODY(I) 
      ECOLL = ECOLL + ZMU*(HI - H(IPAIR))
      HIP = H(IPAIR)
*
*       Change KS variables at original ECC and modify ECC at H = const.
      CALL EXPAND2(IPAIR,SEMI)
      CALL KSPERI(IPAIR)
*       Note simplified versions of standard routines (STEP(I1) is given).
      CALL DEFORM2(IPAIR,ECC,ECC1)
*
*       Re-initialize the KS binary (unperturbed case also needed).
*     IMOD = 1
      CALL RESOLV(IPAIR,1)
*     CALL KSPOLY(IPAIR,IMOD)
*
      ITER = ITER + 1
      IF (ITER.LT.1000.OR.MOD(ITER,1000).EQ.0) THEN
          WRITE (94,30)  TIME+TOFF, ECC, THETA, DT, TGR, SEMI, DW
   30     FORMAT (' GR SHRINK    T E TH DT TGR A DW ',
     &                           F11.4,F9.5,1P,3E9.1,E12.4,E10.2)
          CALL FLUSH(94)
          IF (ITER.GT.2000000000) ITER = 1000
      END IF
*
*       Check KS termination for weak PN.
      IF (DW.GT.1.0D-03) THEN
*       Note that first order Peters formulation is not valid for strong GR.
*         JP = LIST(2,I)
*         LIST(1,I1) = 1      ! Suppressed 06/16.
*         LIST(2,I1) = JP
*       Set PN indicator for ARCHAIN (DW limit means small TGR).
*         IPN = 2
*         WRITE (6,40)  JP, NAME(JP), STEP(I1), STEP(I), SEMI, TGR, DW
*  40     FORMAT (' ENFORCED PERTURB    JP NM S1 SI A TGR DW ',
*    &                                  2I6,1P,5E10.2)
          GO TO 100
      END IF
*
*       Activate coalescence condition using local index.
      JPHASE = 0
      IF ((SEMI1.LT.100.0*RZ.AND.TGR.LT.0.1).OR.TGR.LT.0.01) THEN
          WRITE (6,45)  KSTAR(I1), KSTAR(I2), RADIUS(I1), RADIUS(I2),
     &                  SEMI1
   45     FORMAT (' PN COAL    K* R1 R2 A  ',2I4,1P,3E10.2)
          CALL FLUSH(6)
          IQCOLL = -2
          JPHASE = -1
*         KSPAIR = IPAIR
*         CALL CMBODY(SEMI1,2)
*
*       Include optional kick velocity of 5*VRMS km/s for coalescence recoil.
          IF (KZ(43).GT.0.AND.MIN(KSTAR(I1),KSTAR(I2)).GT.13) THEN
*       Initialize basic variables at start of new step.
              VI20 = 0.0
              DO 48 K = 1,3
                  X0(K,I) = X(K,I)
                  X0DOT(K,I) = XDOT(K,I)
                  VI20 = VI20 + XDOT(K,I)**2
   48         CONTINUE
*
              VF = 5.0*(VRMS/VSTAR)/SQRT(VI20)
              DO 50 K = 1,3
                  XDOT(K,I) = VF*XDOT(K,I)
                  X0DOT(K,I) = XDOT(K,I)
   50         CONTINUE
              ECD0 = ECDOT
              ECDOT = ECDOT + 0.5*BODY(I)*VI20*(1.0 - VF**2)
              VESC = 5.0*VRMS
              WRITE (6,60)  VF, ECD0-ECDOT, VESC
   60         FORMAT (' COALESCENCE KICK    VF ECDOT VESC ',
     &                                      F7.3,F10.6,F6.1)
*       Form neighbour list and new polynomials.
              RS0 = RS(I)
              CALL NBLIST(I,RS0)
              CALL FPOLY1(I,I,0)
              CALL FPOLY2(I,I,0)
          END IF
      END IF
*
      IF (ISTAT(KCASE).EQ.0) ISTAT(KCASE) = JPHASE
*
  100 RETURN
*
      END
