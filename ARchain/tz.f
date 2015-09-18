      
      CLIGHT = 18000.0
      ECC = 0.308
      SEMI = 7.0D-07
      ECC = 0.99
      SEMI = 2.0D-05
      RAU = 1.0*2.0D+05
      SMU = 6.1D+04
      CLIGHT = 10000.0
      ECC = 0.54
      SEMI = 3.0D-08
      SMU = 4.4D+07
      ECC = 0.0
      SEMI = 4.3D-10
      SMU = 8.9D+07
      TAUGR = 1.3D+18*RAU**4/SMU**3
      WRITE (6,1) TAUGR
    1 FORMAT (' TAUGR   ',1P,E10.2)
          ECC2 = ECC**2
          FE = 1.0 + (73.0/24.0 + 37.0*ECC2/96.0)*ECC2
          GE = (1.0 - ECC2)**3.5/FE
          ZX = 3.0D-04
      ZX = 8.0D-04
          RATIO = 1.0
      ZX = 2.0D-05
      RATIO = 0.5
*         RATIO = 7.6/0.8
*       Replace physical time-scale by N-body units (cf. Lee 1993).
*         TZ = TAUGR*GE*SEMI**4/(RATIO*(1.0 + RATIO)*ZX**3)
          TZ = GE*SEMI**4/(RATIO*(1.0 + RATIO)*ZX**3)
      WRITE (6,3)  SEMI, ZX, TZ
    3 FORMAT (' SEMI ZX TZ  ',1P,3E10.2)
          TZ = 5.0/64.0*CLIGHT**5*TZ
      WRITE (6,6)  TZ
    6 FORMAT (' TZ  ',1P,E10.2)
      STOP
      END 
