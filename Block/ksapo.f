      SUBROUTINE KSAPO(IPAIR)
*
*
*       Apocentre/pericentre/random KS variables.
*       -----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  RAN2
*
*
*       Specify half regularized period (PI/2) or random phase (unperturbed).
      IF (IPAIR.GT.0) THEN
          THETA = 0.5*TWOPI
          IKICK = 0
      ELSE IF (IPAIR.LT.0) THEN
          THETA = 0.5*TWOPI*RAN2(IDUM1)
          IPAIR = -IPAIR
*       Note new type not known here but WD case kick decided by option #25.
          IKICK = 1
          IF (LIST(1,2*IPAIR-1).GT.0) THETA = 0.0D0
*       Initialize time because HDOT & TDOT2 not updated for RESOLV.
          T0(2*IPAIR-1) = TIME
*       Skip hyperbolic orbit (i.e. kick for second binary component).
          IF (H(IPAIR).GT.0.0) GO TO 30
      ELSE
          IKICK = -1
          IPAIR = KSPAIR
          SEMI = -0.5D0*BODY(N+IPAIR)/H(IPAIR)
          ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
     &                          TDOT2(IPAIR)**2/(BODY(N+IPAIR)*SEMI)
          ECC = SQRT(ECC2)
*
*       Increase angle near small pericentre (be careful for SEMI < 0).
          IF (ECC.LT.0.99) THEN
              TH = -ECC
              THETA = ACOS(TH)
          ELSE IF (SEMI.GT.0.0) THEN
              THETA = 0.5*TWOPI
          ELSE
              THETA = 1.0
          END IF
*
*         WRITE (6,5)  ECC, THETA, SEMI, SEMI*(1.0 - ECC), R(IPAIR),
*    &                 TDOT2(IPAIR)
*   5     FORMAT (' KSAPO    E THETA A PM R TD2 ',F10.6,F7.3,1P,4E10.2)
      END IF
*
*       Use half the angle for KS to reach apocentre (bug fix 11/2015).
      THETA = 0.5*THETA
*       Form transformation coefficients (Stiefel & Scheifele p. 85).
      XC = COS(THETA)
      YS = SIN(THETA)
      FF = SQRT(0.5D0*ABS(H(IPAIR)))
      R(IPAIR) = 0.0D0
      TDOT2(IPAIR) = 0.0D0
*
*       Generate analytical solutions for U & UDOT using old U0 & UDOT.
      DO 10 K = 1,4
          U(K,IPAIR) = U0(K,IPAIR)*XC + UDOT(K,IPAIR)*YS/FF
          UDOT(K,IPAIR) = UDOT(K,IPAIR)*XC - U0(K,IPAIR)*YS*FF
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0D0*U(K,IPAIR)*UDOT(K,IPAIR)
   10 CONTINUE
*
*       Impose R' < 0 for apocentre procedures (IKICK = 0).
      SEMI = -0.5D0*BODY(N+IPAIR)/H(IPAIR)
      IF (TDOT2(IPAIR).GT.0.0D0.AND.R(IPAIR).GT.SEMI) THEN
          IF (IKICK.EQ.0) THEN
              TDOT2(IPAIR) = -1.0E-20
          END IF
      END IF
*
*       Include diagnostic check that correct apocentre has been set.
*     WRITE (6,20)  SEMI, R(IPAIR), H(IPAIR), SEMI*(1.0 + ECC)
*  20 FORMAT (' APOCENTRE:    A RA H APO ',1P,4E12.4,6E10.1)
*
*       Save KS parameters for WD or neutron star kick (routine FCORR).
   30 IF (IKICK.GT.0) THEN
          CALL KICK(IPAIR,0,0,0.0D0)
      END IF
*
      RETURN
*
      END
