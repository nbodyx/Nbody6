      SUBROUTINE DEFORM2(IPAIR,ECC0,ECC)
*
*
*       Deformation of elliptic orbit.
*       ------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Evaluate square regularized velocity.
      V20 = 0.0
      DO 10 K = 1,4
          V20 = V20 + UDOT(K,IPAIR)**2
   10 CONTINUE
*
*       Form KS coordinate & velocity scaling factors at peri or apocentre.
      I = N + IPAIR
      SEMI = -0.5*BODY(I)/H(IPAIR)
      IF (R(IPAIR).LT.SEMI) THEN
          EFAC = (1.0 - ECC)/(1.0 - ECC0)
          RNEW = SEMI*(1.0 - ECC)
      ELSE
          EFAC = (1.0 + ECC)/(1.0 + ECC0)
          RNEW = SEMI*(1.0 + ECC)
          RNEW = EFAC*R(IPAIR)
      END IF
      C1 = SQRT(EFAC)
      V2 = 0.5*(BODY(I) + H(IPAIR)*RNEW)
      IF(V2.LE.0.D0)THEN
         C2 = 1.0D-06
      ELSE
         C2 = SQRT(V2/V20)
      ENDIF
*
*       Re-scale KS variables at constant energy to new eccentricity (ECC).
      R(IPAIR) = 0.0D0
*       Retain sign of radial velocity for unperturbed KS (apo or peri).
      TDOT2(IPAIR) = 0.0D0
      DO 20 K = 1,4
          U(K,IPAIR) = C1*U(K,IPAIR)
          UDOT(K,IPAIR) = C2*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0*U(K,IPAIR)*UDOT(K,IPAIR)
   20 CONTINUE
*
      RETURN
*
      END
