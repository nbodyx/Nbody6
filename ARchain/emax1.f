      SUBROUTINE EMAX1(MB,XX,VV,XCM,VCM,ECC2,EX,EM)
*
*       Maximum & minimum eccentricity.
*       -------------------------------
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  XX(3,3),VV(3,3),XCM(3),VCM(3),XREL(3),VREL(3),
     &        A1(3),A2(3),EI(3),HI(3),HO(3),MB
*
*
*       Prepare evaluation of maximum eccentricity (see HIGROW).
      DO 10 K = 1,3
          XREL(K) = XX(K,1) - XX(K,2)
          VREL(K) = VV(K,1) - VV(K,2)
   10 CONTINUE
      A12 = 0.0
      A22 = 0.0
      A1A2 = 0.0
      RI2 = 0.0
      VI2 = 0.0
      RVI = 0.0
      DO 12 K = 1,3
          K1 = K + 1
          IF (K1.GT.3) K1 = 1
          K2 = K1 + 1
          IF (K2.GT.3) K2 = 1
          A1(K) = XREL(K1)*VREL(K2) - XREL(K2)*VREL(K1)
          A2(K) = (XX(K1,3) - XCM(K1))*(VV(K2,3) - VCM(K2))
     &          - (XX(K2,3) - XCM(K2))*(VV(K1,3) - VCM(K1))
          A12 = A12 + A1(K)**2
          A22 = A22 + A2(K)**2
          A1A2 = A1A2 + A1(K)*A2(K)
          RI2 = RI2 + XREL(K)**2
          VI2 = VI2 + VREL(K)**2
          RVI = RVI + XREL(K)*VREL(K)
   12 CONTINUE
*
*       Construct the Runge-Lenz vector (Heggie & Rasio 1995, Eq.(5)).
      EI2 = 0.0
      DO 15 K = 1,3
          EI(K) = (VI2*XREL(K) - RVI*VREL(K))/MB - XREL(K)/SQRT(RI2)
          EI2 = EI2 + EI(K)**2
   15 CONTINUE
      EI2 = MIN(EI2,0.99999D0)
*
*       Define unit vectors for inner eccentricity and angular momenta.
      COSJ = 0.0
      SJSG = 0.0
      DO 20 K = 1,3
          EI(K) = EI(K)/SQRT(EI2)
          HI(K) = A1(K)/SQRT(A12)
          HO(K) = A2(K)/SQRT(A22)
          COSJ = COSJ + HI(K)*HO(K)
          SJSG = SJSG + EI(K)*HO(K)
   20 CONTINUE
*
*       Evaluate the expressions A & Z (relativistic ECC2 introduced 3/7/11).
      A = COSJ*SQRT(1.0 - ECC2)
*     A = COSJ*SQRT(1.0 - EI2)
*     Z = (1.0 - EI2)*(2.0 - COSJ**2) + 5.0*EI2*SJSG**2
      Z = (1.0 - ECC2)*(2.0 - COSJ**2) + 5.0*ECC2*SJSG**2
*
*       Obtain maximum inner eccentricity (Douglas Heggie, Sept. 1995).
      Z2 = Z**2 + 25.0 + 16.0*A**4 - 10.0*Z - 20.0*A**2 - 8.0*A**2*Z
      Z2 = MAX(Z2,0.0D0)
      EX2 = (Z + 1.0 - 4.0*A**2 + SQRT(Z2))/6.0
      EX2 = MAX(EX2,0.0001D0)
      EX = SQRT(EX2)
*
*       Form minimum eccentricity (Douglas Heggie, Sept. 1996).
      AZ = A**2 + Z - 2.0
      IF (AZ.GE.0.0) THEN
          AZ1 = 1.0 + Z - 4.0*A**2
          EMIN2 = (AZ1 - SQRT(AZ1**2 - 12.0*AZ))/6.0
      ELSE
          EMIN2 = 1.0 - 0.5*(A**2 + Z)
      END IF
      EMIN2 = MAX(EMIN2,0.0D0)
      EM = SQRT(EMIN2)
*
      RETURN
*
      END
