      SUBROUTINE CSTAB4(ITERM)
*
*
*       Four-body ARC stability & termination test.
*       -------------------------------------------
*
      IMPLICIT REAL*8 (A-H,M,O-Z)
      PARAMETER (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2,NMX4=4*NMX)
      COMMON/ARCHAIN/X(NMX3),V(NMX3),WTTL,M(NMX), 
     &   XC(NMX3),WC(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),N 
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,JCOLL,NDISS1
      REAL*8  M,MB,MB1,MB2,R2(NMX,NMX),XCM(3),VCM(3),XX(3,3),VV(3,3),
     &        XCM2(3),VCM2(3)
      INTEGER  IJ(NMX)
      DATA NAME34,I34 /0,0/
      SAVE
*
*
*       Sort particle separations (Note I3 & I4 form closest pair).
      CALL R2SORT(IJ,R2)
*
*       Save indices with I1 & I2 as wide binary and I3 & I4 as outer body.
      I1 = IJ(3)
      I2 = IJ(4)
      I3 = IJ(1)
      I4 = IJ(2)
      MB = M(I1) + M(I2)
      MB2 = M(I3) + M(I4)
      MB1 = MB + MB2
*
*       Use perturbed triple stability instead for small middle distance.
      IF (1.0/RINV(2).LT.0.01*RSUM) THEN
          CALL CSTAB3(ITERM)
          IF (ITERM.LT.0) GO TO 40
      END IF
*
*       Form orbital parameters with c.m. of I3 & I4 as third body.
      VREL2 = 0.0D0
      VREL21 = 0.0D0
      VREL34 = 0.0D0
      RDOT = 0.0D0
      RDOT3 = 0.0D0
      RB0 = 0.0
      RI2 = 0.0D0
      D1 = 0.0
      D2 = 0.0
      VR21 = 0.0
      VR22 = 0.0
      RDOT1 = 0.0
      RDOT2 = 0.0
      DO 10 K = 1,3
          J1 = 3*(I1-1) + K
          J2 = 3*(I2-1) + K
          J3 = 3*(I3-1) + K
          J4 = 3*(I4-1) + K
          RB0 = RB0 + (X(J1) - X(J2))**2
          VREL2 = VREL2 + (V(J1) - V(J2))**2
          RDOT = RDOT + (X(J1) - X(J2))*(V(J1) - V(J2))
          XCM(K) = (M(I1)*X(J1) + M(I2)*X(J2))/MB
          VCM(K) = (M(I1)*V(J1) + M(I2)*V(J2))/MB
          XCM2(K) = (M(I3)*X(J3) + M(I4)*X(J4))/MB2
          VCM2(K) = (M(I3)*V(J3) + M(I4)*V(J4))/MB2
          RI2 = RI2 + (XCM2(K) - XCM(K))**2
          VREL21 = VREL21 + (VCM2(K) - VCM(K))**2
          VREL34 = VREL34 + (V(J3) - V(J4))**2
          RDOT3 = RDOT3 + (XCM2(K) - XCM(K))*(VCM2(K) - VCM(K))
          XX(K,1) = X(J1)
          XX(K,2) = X(J2)
          XX(K,3) = XCM2(K)
          VV(K,1) = V(J1)
          VV(K,2) = V(J2)
          VV(K,3) = VCM2(K)
          D1 = D1 + (XCM2(K) - X(J1))**2
          D2 = D2 + (XCM2(K) - X(J2))**2
          VR21 = VR21 + (VCM2(K) - V(J1))**2
          VR22 = VR22 + (VCM2(K) - V(J2))**2
          RDOT1 = RDOT1 + (XCM2(K) - X(J1))*(VCM2(K) - V(J1))
          RDOT2 = RDOT2 + (XCM2(K) - X(J2))*(VCM2(K) - V(J2))
   10 CONTINUE
*
*       Evaluate orbital elements for inner and outer motion.
      RB = SQRT(R2(I1,I2))
      R3 = SQRT(RI2)
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      ECC = SQRT((1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
      SEMI1 = 2.0/R3 - VREL21/MB1
      SEMI1 = 1.0/SEMI1
      ECC1 = SQRT((1.0D0 - R3/SEMI1)**2 + RDOT3**2/(SEMI1*MB1))
      RB2 = SQRT(R2(I3,I4))
      SEMI2 = 2.0/RB2 - VREL34/MB2
      SEMI2 = 1.0/SEMI2
*       Obtain the inclination.
      CALL INCLIN(XX,VV,XCM,VCM,ALPHA)
*
*       Skip if innermost triple is not stable (including ECC & ECC1 > 1).
      IF (SEMI.LT.0.0.OR.ECC.GT.1.0.OR.ECC1.GT.0.9999) THEN  ! Use E1 > 0.999.
          ITERM = 1
          GO TO 40
      END IF
*       Evaluate the Valtonen stability criterion.
      QST = QSTAB(ECC,ECC1,ALPHA,M(I1),M(I2),MB2)
      PMIN = SEMI1*(1.0D0 - ECC1)
*     WRITE (6,99)  ECC1, ALPHA, QST, PMIN, QST*SEMI, SEMI2
*  99 FORMAT (' ECC1 ALPHA QST PM Q*A A2 ',3F8.4,1P,3E10.2)
*     CALL FLUSH(6)
      IF (QST*SEMI.LT.PMIN) THEN
*       Allow for small theoretical correction factor.
          ZF = ABS(SEMI2)/SEMI
          ZF = MIN(ZF,0.2)
          PCRIT = QST*SEMI*(1.0 + 0.1*ZF)
      ELSE
          ITERM = 1
          GO TO 35
      END IF
*
*       Check hierarchical stability condition for two separated binaries.
      ITERM = 0
      IF (PMIN.GT.PCRIT.AND.RB.GT.SEMI.AND.RDOT3.GT.0.0.AND.
     &    R3.GT.4.0*SEMI) THEN
          ITERM = -1
          ALPHA = 180.0*ALPHA/3.1415
          WRITE (6,20)  ECC, ECC1, SEMI, SEMI1, PMIN, PCRIT, ALPHA
   20     FORMAT (' QUAD HIARCH    E =',F6.3,'  E1 =',F6.3,
     &                     '  A =',1P,E8.1,'  A1 =',E8.1,'  PM =',E9.2,
     &                     '  PC =',E9.2,'  IN =',0P,F6.1)
          RI = SQRT(CM(1)**2 + CM(2)**2 + CM(3)**2)
          WRITE (81,30)  TIMEC, RI, NAMEC(I3), ECC, ECC1, SEMI, SEMI1,
     &                   PCRIT/PMIN, ALPHA
   30     FORMAT (' CSTAB4  ',F9.5,F5.1,I6,2F6.3,1P,2E10.2,0P,F5.2,F6.1)
          CALL FLUSH(81)
      END IF
*
*       Include scheme for switching from ARC (no BH) to KS with perturbers.
   35 IF (ISTAR(I3).GT.13.OR.ISTAR(I4).GT.13) GO TO 40
*
*       Evaluate two-body parameters with respect to small binary I3 & I4.
      D1 = SQRT(D1)
      D2 = SQRT(D2)
      GAM1 = 2.0*M(I1)/(MB2 + M(I1))*(SEMI2/D1)**3
      GAM2 = 2.0*M(I2)/(MB2 + M(I2))*(SEMI2/D2)**3
      IF (MAX(GAM1,GAM2).GT.1.0D-04) GO TO 40
*
*       Form outer semi-major axis, eccentricity and pericentre.
      A1 = 2.0/D1 - VR21/(MB2 + M(I1))
      A2 = 2.0/D2 - VR22/(MB2 + M(I2))
      A1 = 1.0/A1
      A2 = 1.0/A2
      E1 = (1.0 - D1/A1)**2 + RDOT1**2/((MB2 + M(I1))*A1)
      E2 = (1.0 - D2/A2)**2 + RDOT2**2/((MB2 + M(I2))*A2)
      E1 = SQRT(E1)
      E2 = SQRT(E2)
      PM1 = A1*(1.0 - E1)
      PM2 = A2*(1.0 - E2)
*       Restrict pericentre to 20*SEMI (tighter than GAM1/2).
      IF (PM1.LT.20.0*SEMI2.OR.PM2.LT.20.0*SEMI2) GO TO 40
*
*       Consider the energy budget.
      NAM34 = NAMEC(I3) + NAMEC(I4)
      IF (I34.EQ.0) THEN
          NAME34 = NAM34
          I34 = I34 + 1
      ELSE IF (NAM34.EQ.NAME34) THEN
          I34 = I34 + 1
          IF (I34.GE.10) THEN
              EB = -0.5*M(I3)*M(I4)/SEMI2
*       Define energetic binary as 4 % of total initial binding energy.
              IF (EB.LT.-0.01) THEN
                  ITERM = -1
                  WRITE (6,38)  NAMEC(I3), NAMEC(I4), I34, EB, SEMI2
   38             FORMAT (' ENFORCED CSTAB4    NAME I34 EB A2 ',
     &                                         2I6,I4,F8.3,1P,E10.2)
                  CALL FLUSH(6)
                  I34 = 0
                  NAME34 = 0
              END IF
          END IF
      ELSE
*       Initialize new binary name or make new start.
          I34 = I34 + 1
          IF (I34.GT.50) THEN
              I34 = 0
              NAME34 = 0
          END IF
      END IF
*
   40 RETURN
*
      END
