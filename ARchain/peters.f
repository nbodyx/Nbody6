      IMPLICIT REAL*8 (A-H,O-Z)
      REAL M1,M2
*
      M1 = 3.0D-04
      M2 = 3.0D-04
      ECCGR = 0.999
      SEMIGR = 2.0D-05
      CVEL = 18000.0
      CVEL = 10000.0
      M1 = 7.6
      M2 = 0.8
      ECCGR = 0.54
      SEMIGR = 3.0D-08
      E2 = ECCGR**2
      ADOT =  64.0/5.0*M1*M2*(M1+M2)/(CVEL**5*SEMIGR**3*(1.0-E2)**3.5)
      ADOT = ADOT*(1.0 + 73.0/24.0*E2 + 37.0/96.0*E2**2)
      EDOT = 304.0/15.0*ECCGR*M1*M2*(M1+M2)/(CVEL**5*SEMIGR**4)
      EDOT = EDOT/(1.0-E2)**2.5*(1.0 + 121.0/304.0*E2)
      TGR = SEMIGR/ADOT
      EGR = ECCGR/EDOT
*
      SMU = 6.0D+04
      ASUN = SEMIGR*2.0D+05*215.0
*       Include Banerjee formula in yrs (SEMI in solar radii & SMU in Msun).
      TM = 1.5D+08*ASUN**4/(M1*SMU)**3*(1.0-E2)**3.5
      THETA = 3.0*6.2830*(M1 + M2)/(SEMIGR*CVEL**2*(1.0-E2))
      WRITE (6,10)  ECCGR, SEMIGR, TGR, EGR, THETA, TM
   10 FORMAT (' E A TGR EGR THETA TM ',F10.5,1P,5E10.2)
      END
