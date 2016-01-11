      SUBROUTINE PNPERT2(M1,M2,X,V,ACC,DER,DT,IC,CLIGHT,DE)
*
*
*       Post-Newtonian perturbations.
*       -----------------------------
*
*       L. Blanchet & B. Iyer, Class. Quantum Grav. 20 (2003), 755.
*       T. Mora & Clifford Will, gr-qc/0312082.
*       ---------------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/POSTN3/  ESAVE(10)
      PARAMETER (ONE12=1.0/12.0D0)
      SAVE
      REAL*8  ACC(3),DER(3),VDOT(3),X(3),V(3)
      DATA ED0,ED20 /2*0.0D0/
      DATA IIC /0/
*
*
*       Save inverse powers of CLIGHT.
      IF (IIC.EQ.0) THEN
          CL2 = 1.0/CLIGHT**2
          CL3 = 1.0/CLIGHT**3
          ESAVE(1) = 0.0
          ESAVE(2) = 0.0
          ESAVE(3) = 0.0
          IIC = 1
      END IF
*
*       Set total mass and reduced mass parameter.
      M = M1 + M2
      ETA = M1*M2/M**2
      ZMU = M*ETA
*
*       Form useful scalars.
      R2 = X(1)**2 + X(2)**2 + X(3)**2
      R = SQRT(R2)
      MR = M/R
      V2 = V(1)**2 + V(2)**2 + V(3)**2
      RD = (X(1)*V(1) + X(2)*V(2) + X(3)*V(3))/R
*
*     A1 = 2*(2+ETA)*MR - (1+3*ETA)*V2 + 1.5*ETA*RD**2
*     B1 = 2*(2-ETA)*RD
*
      A1 = 0.0   ! including first order PN leads to large error (10/14).
      B1 = 0.0
      HT = 0.5*V2 - MR
      SEMI = -0.5*M/HT
*       Define unscaled PN2.5 (positive sign convention).
      A25 = 8.0/5.0*ETA*MR*RD*(17.0/3.0*MR + 3.0*V2)
      A = (A1 + A25*CL3)*CL2
*
      B25 = -8.0/5.0*ETA*MR*(3.0*MR + V2)
      B = (B1 + B25*CL3)*CL2
*
*       Set the standard PN2.5 perturbing acceleration for VDOT.
      DO 5 K = 1,3
          ACC(K) = M/R2*(A*X(K)/R + B*V(K))
    5 CONTINUE
*
*     WRITE (15,6) R, A, B, CLIGHT
*   6 FORMAT (' R A B C  ',1P,4E10.2)
*       Note ACC scaled by M/R**2/C**3 for actual perturbing force.
      RVD = 0.0
      VVDOT = 0.0
      DO 10 K = 1,3
          VDOT(K) = ACC(K) - MR*X(K)/R2
          RVD = RVD + X(K)*VDOT(K)
          VVDOT = VVDOT + V(K)*VDOT(K)
   10 CONTINUE
*       Form d2R/dt**2 from dR/dt & d(R*V)/dt and absorb 1/R in ADOT.
      RD2 = V2 + RVD - RD**2
*
*     AD1 = -2*(2+ETA)*MR*RD/R - 2*(1+3*ETA)*VVDOT + 3*ETA*RD*RD2/R
*     BD1 = 2*(2-ETA)*RD2/R
      AD1 = 0.0
      BD1 = 0.0
*
*       Adopt Peter Berczik's A & B derivatives (these have correct sign!).
      ADOT = (17.0/3.0*MR + 3.0*V2)*(RD2 - RD**2)*M/R2 +
     &       MR*RD*(-17.0/3.0*M*RD/R2 + 6.0*VVDOT)
      BDOT = MR*((3.0*MR + V2)*RD/R + 3.0*M*RD/R2 - 2.0*VVDOT)
      ADOT = 1.6*ETA*ADOT
      BDOT = 1.6*ETA*BDOT
*
      ADOT = (AD1 + ADOT*CL3)*CL2
      BDOT = (BD1 + BDOT*CL3)*CL2
*
*       Use equation of motion dV/dt = M/R2*((-1 + A/C**5)*X/R + B*V/C**5).
      DO 20 K = 1,3
*       Note that M/R2 is omitted in all terms (see final scaling).
          DER(K) = -2.0*RD/R*(A*X(K)/R + B*V(K)) + ADOT*X(K)/R -
     &              A*RD*X(K)/R2 + A*V(K)/R + BDOT*V(K) + B*VDOT(K)
   20 CONTINUE
*
*       Scale the acceleration and derivative by the leading term M/R**2.
      GMC = M/R2
      DO 30 K = 1,3
          DER(K) = GMC*DER(K)
   30 CONTINUE
*
*     WRITE (6,78)  ACC, DER
* 78  FORMAT (' ACC DER  ',1P,8E10.2)
*       Employ Hermite 4th-order integration for the energy equation.
      ED = 0.0
      ED2 = 0.0
*       Form the first and second energy derivative.
      DO 40 K = 1,3
          ED = ED + ACC(K)*V(K)
          ED2 = ED2 + ACC(K)*VDOT(K) + DER(K)*V(K)
   40 CONTINUE
*
*       Copy previous values from COMMON.
      ED0 = ESAVE(1)
      ED20 = ESAVE(2)
*
*       Construct standard Hermite third and 4th-order derivatives.
***   ED3 = 2.0*(-3.0*(ED0 - ED) - (2.0*ED20 + ED2)*DT)/DT**2
***   ED4 = 6.0*(2.0*(ED0 - ED) + (ED20 + ED2)*DT)/DT**3
*
*       Get final result from Keigo's Hermite corrector.
      DE = 0.5*(ED0 + ED)*DT + ONE12*(ED20 - ED2)*DT**2
      DE = ZMU*DE
*
*     WRITE (6,81)  R, ED, ED2, DE, ESAVE(3), SEMI
*  81 FORMAT (' FINAL    R ED ED2 DE SUM A ',1P,5E10.2,E12.4)
*       Save the two derivatives for next step (safer using COMMON).
      ED0 = ED
      ED20 = ED2
      ESAVE(1) = ED
      ESAVE(2) = ED2
      IF (IC.EQ.0) DE = 0.0
      ESAVE(3) = ESAVE(3) + DE
      IF (MOD(IC,1000).EQ.0) THEN
          HT = 0.5*V2 - MR
          SEMI = -0.5*M/HT
          WRITE (6,50)  IC, R, ESAVE(3), SEMI
   50     FORMAT (' PNPERT2 ENERGY    IC R EGR A ',I9,1P,2E10.2,E14.6)
      END IF
*
      RETURN
*
      END
