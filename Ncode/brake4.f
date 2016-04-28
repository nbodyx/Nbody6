      SUBROUTINE BRAKE4(I1,I2,DT)
*
*
*       GR analytical orbit shrinkage.
*       ------------------------------
*
      INCLUDE 'common6.h'
      COMMON/POSTN/  CVEL,TAUGR,RZ1,GAMMAZ,TKOZ,EMAX,TSP,KZ24,IGR,IPN
      REAL*8 M1,M2,UI(4),UIDOT(4)
      SAVE ADOT0,EDOT0,ITER
      DATA ADOT0,EDOT0,ITER /0.0D0,0.0D0,0/
*
*
*       Check relativistic conditions (at least one >= NS).
      IF (MAX(KSTAR(I1),KSTAR(I2)).LT.13) GO TO 100
*
*       See whether CLIGHT has been initialized in ARchain.
      IF (ITER.EQ.0) THEN
          READ (5,*)  CLIGHT
          ITER = 1
      END IF
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
      DW = 3.0*TWOPI*(BODY(I1) + BODY(I2))/(SEMI*CVEL**2*(1.0-E2))
      IF (DW.LT.1.0D-04) GO TO 100
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
      TK = TWOPI*SEMI*SQRT(SEMI/(BODY(I1) + BODY(I2)))
      THETA = DW*DT/TK
*     DT1 = MIN(DT,1.0D-05*SEMI/ADOT)
      DT1 = 1.0D-05*SEMI/ADOT
*       Impose limit of time-step if THETA > TWOPI.
      IF (THETA.GT.TWOPI) THEN
          DT1 = TWOPI*TK/DW
          THETA = DMOD(THETA,TWOPI)
      END IF
      DT = DT1
      T00 = T0(I1)
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
*       Define RZ and obtain the incremental change in SEMI & ECC.
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
*       Update binding energy and collision energy.
      HI = H(IPAIR)
      H(IPAIR) = -0.5*BODY(I)/SEMI1
      ZMU = BODY(I1)*BODY(I2)/BODY(I) 
      ECOLL = ECOLL + ZMU*(HI - H(IPAIR))
      HIP = H(IPAIR)
      NSTEPQ = NSTEPQ + 1
*
*       Change KS variables at original ECC and modify ECC at H = const.
      CALL EXPAND2(IPAIR,SEMI)
      CALL KSPERI(IPAIR)
*       Note simplified versions of standard routines (new STEP(I1) given).
      CALL DEFORM2(IPAIR,ECC,ECC1)
*
*       Re-initialize the KS binary (unperturbed case also needed).
      IMOD = 1
      CALL RESOLV(IPAIR,1)
      CALL KSPOLY(IPAIR,IMOD)
*       Restore the original values of STEP & T0.
      STEP(I1) = DT1
      T0(I1) = T00
*
      ITER = ITER + 1
      IF (ITER.LT.1000.OR.MOD(ITER,1000).EQ.0) THEN
          WRITE (94,30)  TIME+TOFF, ECC, THETA, DT, TGR, DW, SEMI
   30     FORMAT (' GR SHRINK    T E TH DT TGR DW A ',
     &                           F11.4,F9.5,1P,4E9.1,E12.4)
          CALL FLUSH(94)
      END IF
*
*       Check KS termination with PN activation in KSINT.
      IF (DW.GT.2.0D-03) THEN
*       Note that first order Peters formulation is not valid for strong GR.
          IF (ITER.LT.1000.OR.MOD(ITER,1000).EQ.0) THEN
              WRITE (6,40)  IPN, STEP(I1), STEP(I), SEMI, TGR
   40         FORMAT (' ENFORCED CHAIN    IPN S1 SI A TZ',
     &                                    I4,1P,5E10.2)
              CALL FLUSH(6)
          END IF
*       Set PN indicator for ARCHAIN (DW limit means small TGR).
          IPN = 2
          GO TO 100
      END IF
*
*       Activate coalescence condition using COMMON index.
      IF ((SEMI1.LT.100.0*RZ.AND.TGR.LT.0.1).OR.TGR.LT.0.01) THEN
          WRITE (6,45)  KSTAR(I1), KSTAR(I2), RADIUS(I1), RADIUS(I2),
     &                  SEMI1
   45     FORMAT (' PN COAL    K* R1 R2 A  ',2I4,1P,3E10.2)
          CALL FLUSH(6)
          IQCOLL = -2
          IPHASE = -1
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
  100 RETURN
*
      END
