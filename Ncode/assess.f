      SUBROUTINE ASSESS(IPAIR,IM,ECC,SEMI,ITERM)
*
*
*       Assessment of hierarchical stability.
*       -------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      REAL*8  XX(3,3),VV(3,3)
      SAVE ITIME
      DATA ITIME /0/
*
*
*       Define indices for components and c.m.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      ICM = N + IPAIR
*
*       Determine inclination between inner relative motion and outer orbit.
      RV = 0.0
      DO 10 K = 1,3
          XX(K,1) = XREL(K,IM)
          XX(K,2) = 0.0
          XX(K,3) = X(K,I2)
          VV(K,1) = VREL(K,IM)
          VV(K,2) = 0.0
          VV(K,3) = XDOT(K,I2)
          RV = RV + XREL(K,IM)*VREL(K,IM)
   10 CONTINUE
*
      CALL INCLIN(XX,VV,X(1,ICM),XDOT(1,ICM),ANGLE)
*
*       Form inner semi-major axis and eccentricity.
      RIN = SQRT(XREL(1,IM)**2 + XREL(2,IM)**2 + XREL(3,IM)**2)
      SEMI0 = -0.5*BODY(I1)/HM(IM)
      ECC2 = (1.0 - RIN/SEMI0)**2 + RV**2/(BODY(I1)*SEMI0)
      ECC0 = SQRT(ECC2)
      PMIN = SEMI*(1.0 - ECC)
*
*       Evaluate the general stability function after perturbation effect.
      IF (ECC.LT.1.0) THEN
*       Include de/dt corrections inside generous stability boundary.
          EOUT = ECC
          IF (PMIN.LT.6.0*SEMI0.OR.EOUT.GT.0.9) THEN
              NNB1 = LIST(1,I1) + 1
              DE = 0.0
              DO 20 L = 2,NNB1      ! Sum positive perturber contributions.
                  JP = LIST(L,I1)
                  CALL EDOT(I1,I2,JP,SEMI,ECC,ECCDOT)
                  IF (ECCDOT.GT.0.0) THEN
                      TK = TWOPI*SEMI*(SQRT(SEMI)/BODY(ICM))
                      DE = DE - MIN(ECCDOT*TK,0.001D0)
                  END IF
   20         CONTINUE
              EOUT = EOUT - DE
              EOUT = MAX(EOUT,0.0D0)
*       Obtain new stability boundary from Valtonen's criterion (MN 2017).
              QST = QSTAB(ECC0,EOUT,ANGLE,CM(1,IM),CM(2,IM),BODY(I2))
          ELSE
*       Employ basic criterion without EDOT outside safety boundary.
              QST = QSTAB(ECC0,EOUT,ANGLE,CM(1,IM),CM(2,IM),BODY(I2))
          END IF
*
          ITIME = ITIME + 1
          IF (ITIME.GT.2000000000) ITIME = 0
          IF (MOD(ITIME,1000).EQ.0) THEN
              ALPH = 360.0*ANGLE/TWOPI
              PCRIT = QST*SEMI0
              WRITE (6,30) ECC0, ECC, ALPH, SEMI, PCRIT, PMIN
   30         FORMAT (' ASSESS    E0 E1 INC A1 PCR PMIN ',
     &                            2F8.4,F7.1,1P,3E9.1)
          END IF
*       Define termination by QSTAB.
          IF (QST*SEMI0.LT.PMIN) THEN
              ITERM = 0
          ELSE
              ITERM = 1
          END IF
      END IF
*
      RETURN
*
      END
