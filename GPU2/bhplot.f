      SUBROUTINE BHPLOT
*
*
*       Black hole plotting data (also NCH & KS).
*       -----------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/POSTN2/  SEMIGR,ECCGR,DEGR,ISPIN
      REAL*8  XI(8),XCM(3),VCM(3)
*
*
*       Decide between standard BH, KS or ARC.
      IF (NCH.GT.0.OR.KZ(45).GT.1) GO TO 30
*
*       Search for one or two BHs (NAME = 1 or 2 but #24 holds number).
      IBH = 0
      JBH = 0
      IPAIR = 0
      DO 6 L = 1,2
      DO 5 I = 1,N
          IF (NAME(I).EQ.1.AND.IBH.EQ.0) IBH = I
          IF (NAME(I).EQ.2.AND.IBH.GT.0.AND.JBH.EQ.0) JBH = I
    5 CONTINUE
    6 CONTINUE
      IF (KZ(24).LT.2) JBH = 0
*
*       Copy coordinates and velocities (singles or c.m. for simplicity).
      I = IBH
      KK = 0
*       Switch to second component if IBH not identified.
      IF (I.EQ.0) THEN
          I = JBH
          KK = 0
      END IF
   10 IF (I.GE.IFIRST) THEN
          DO 12 K = 1,3
              XI(KK+K) = X(K,I)
   12     CONTINUE
      ELSE IF (I.GT.0) THEN
          IPAIR = KVEC(I)
          DO 15 K = 1,3
              XI(KK+K) = X(K,N+IPAIR)
   15     CONTINUE
      END IF
*       Include central distance as 4th entry.
      XI(KK+4) = SQRT(XI(KK+1)**2 + XI(KK+2)**2 + XI(KK+3)**2)
*       Use c.m. if both BHs in the same pair.
      IF (IPAIR.GT.0) THEN
          IF (NAME(2*IPAIR-1) + NAME(2*IPAIR).EQ.3) GO TO 20
      END IF
*
*       Repeat for possible second BH.
      IF (I.EQ.IBH.AND.JBH.GT.0) THEN
          I = JBH
          KK = KK + 4
          GO TO 10
      END IF
*
*       Produce plotting output in Myr on unit 45 (to match option #45).
   20 IF (IBH.GT.0.OR.JBH.GT.0) THEN
          K1 = 0
          K2 = KK
          WRITE (45,25)  (TIME+TOFF)*TSTAR, (XI(K),K=K1+1,K2+4)
   25     FORMAT (' ',F10.4,8F9.4)
          CALL FLUSH(45)
      END IF
      GO TO 200
*
*       Treat case of BH in KS binary (#45 > 1).
   30 IF (NCH.GT.0.OR.KZ(45).LE.1) GO TO 70
*
      LZ = 0
      DO 35 L = 1,2*NPAIRS
          IF (NAME(L).EQ.1) LZ = L
   35 CONTINUE
      IF (LZ.EQ.0) GO TO 70
      IP = KVEC(LZ)
      I1 = 2*IP - 1
      IF (LIST(1,I1).EQ.0) GO TO 70
      I2 = I1 + 1
      RIJ2 = 0.0
      VI2 = 0.0
      RD = 0.0
      DO 40 K = 1,3
          RIJ2 = RIJ2 + (X(K,I1) - X(K,I2))**2
          VI2 = VI2 + (XDOT(K,I1) - XDOT(K,I2))**2
          RD = RD + (X(K,I1) - X(K,I2))*(XDOT(K,I1) - XDOT(K,I2))
   40 CONTINUE
      RIJ = SQRT(RIJ2)
      ZMB = BODY(I1) + BODY(I2)
      SEMI = 2.0/RIJ - VI2/ZMB
      SEMI = 1.0/SEMI
      ECC2 = (1.0 - RIJ/SEMI)**2 + RD**2/(ZMB*SEMI)
      ECC = SQRT(ECC2)
*
*       Consider the perturbers (if any).
      ICM = N + IP
      NP = LIST(1,I1)
      AX = 1.0
      EX = 0.0
      NMJ = 0
      DO 60 L = 2,NP+1
          J = LIST(L,I1)
          RIJ2 = 0.0
          VI2 = 0.0
          RD = 0.0
          DO 50 K = 1,3
              RIJ2 = RIJ2 + (X(K,J) - X(K,ICM))**2
              VI2 = VI2 + (XDOT(K,J) - XDOT(K,ICM))**2
              RD = RD + (X(K,J) - X(K,ICM))*(XDOT(K,J) - XDOT(K,ICM))
   50     CONTINUE
          ZMB = ZMB + BODY(J)
          RIJ = SQRT(RIJ2)
          AJ = 2.0/RIJ - VI2/ZMB
          AJ = 1.0/AJ
          ECC2 = (1.0 - RIJ/AJ)**2 + RD**2/(ZMB*AJ)
          EJ = SQRT(ECC2)
*       Determine the smallest positive semi-major axis of perturbers.
          IF (AJ.GT.0.0.AND.AJ.LT.AX) THEN
              AX = AJ
              EX = EJ
              NMJ = NAME(J)
          END IF
   60 CONTINUE
*
      IF (ECC.GT.0.90.OR.EX.GT.0.90) THEN
          IF (BODY(I1).GT.BODY(I2)) THEN
              NM1 = NAME(I2)
          ELSE
              NM1 = NAME(I1)
          END IF
          RP = SEMI*(1.0 - ECC)
          RX = AX*(1.0 - EX)
          WRITE (45,65)  (TIME+TOFF)*TSTAR, NM1, ECC, SEMI, NMJ, EX,
     &                   AX, RX
   65     FORMAT (' BHPLOT    T NM E A NMJ EX AX RX ',
     &                        F8.1,I7,F9.5,1P,E9.1,0P,I7,F9.5,1P,2E9.1)
          CALL FLUSH(45)
      END IF
*
*       Check ARC configuration.
   70 IF (NCH.EQ.0.OR.KZ(45).LE.2) GO TO 200
*
      NP = LISTC(1)
      AX = 1.0
      EX = 0.0
      DO 80 L = 2,NP+1
          J = LISTC(L)
          RIJ2 = 0.0
          VI2 = 0.0
          RD = 0.0
          DO 75 K = 1,3
              RIJ2 = RIJ2 + (X(K,J) - X(K,ICH))**2
              VI2 = VI2 + (XDOT(K,J) - XDOT(K,ICH))**2
              RD = RD + (X(K,J) - X(K,ICH))*(XDOT(K,J) - XDOT(K,ICH))
   75     CONTINUE
          ZMB = BODY(ICH) + BODY(J)
          RIJ = SQRT(RIJ2)
          AJ = 2.0/RIJ - VI2/ZMB
          AJ = 1.0/AJ
          ECC2 = (1.0 - RIJ/AJ)**2 + RD**2/(ZMB*AJ)
          EJ = SQRT(ECC2)
*       Determine the smallest positive semi-major axis of perturbers.
          IF (AJ.GT.0.0.AND.AJ.LT.AX) THEN
              AX = AJ
              EX = EJ
              NMJ = NAME(J)
          END IF
   80 CONTINUE
*
      IF (EX.GT.0.90) THEN
          PX = AX*(1.0 - EX)
          WRITE (46,85)  (TIME+TOFF)*TSTAR, NMJ, EX, AX, PX
   85     FORMAT (' BHPLOT    T NMJ EX AX PX ',F8.1,I7,F9.5,1P,2E10.2)
          CALL FLUSH(46)
      END IF
*
*       Include special output for second innermost member (unit #49).
      IF (KZ(45).LE.3) GO TO 200
*       Determine the most distant body.
      RX2 = 0.0
      DO 100 L = 1,NCH
          RIJ2 = 0.0
          DO 95 K = 1,3
              RIJ2 = RIJ2 + XC(K,L)**2
   95     CONTINUE
          IF (RIJ2.GT.RX2) THEN
              RX2 = RIJ2
              I3 = L
          END IF
  100 CONTINUE
*
      DO 105 K = 1,3
          XCM(K) = 0.0
          VCM(K) = 0.0
  105 CONTINUE
*
*       Form c.m. of inner binary.
      ZMB = 0.0
      DO 110 L = 1,NCH
          IF (L.EQ.I3) GO TO 110
          DO 108 K = 1,3
              XCM(K) = XCM(K) + BODYC(L)*XC(K,L)
              VCM(K) = VCM(K) + BODYC(L)*UC(K,L)
  108     CONTINUE
          ZMB = ZMB + BODYC(L)
  110 CONTINUE
      DO 115 K = 1,3
          XCM(K) = XCM(K)/ZMB
          VCM(K) = VCM(K)/ZMB
  115 CONTINUE
*
*       Obtain SEMI & ECC for body #I3.
      RIJ2 = 0.0
      VIJ2 = 0.0
      RD = 0.0
      DO 120 K = 1,3
          RIJ2 = RIJ2 + (XC(K,I3) - XCM(K))**2
          VIJ2 = VIJ2 + (UC(K,I3) - VCM(K))**2
          RD = RD + (XC(K,I3) - XCM(K))*(UC(K,I3) - VCM(K))
  120 CONTINUE
      RIJ = SQRT(RIJ2)
      SEMI = 2.0/RIJ - VIJ2/ZMB
      SEMI = 1.0/SEMI
      ECC2 = (1.0 - RIJ/SEMI)**2 + RD**2/(ZMB*SEMI)
      ECC = SQRT(ECC2)
      STAB = SEMI*(1.0 - ECC)/(SEMIGR*(1.0 + ECCGR))
*
      WRITE (49,130)  (TIME+TOFF)*TSTAR, NCH, BODYC(I3)*SMU, STAB,
     &                                   ECC, ECCGR, SEMI, SEMIGR
  130 FORMAT (' BHPLOT    T NCH M3 APO/PM E EGR A AGR ',
     &                    F8.1,I4,F6.1,F7.1,2F9.5,1P,2E10.2)
      CALL FLUSH(49)
*
  200 CONTINUE
*
      IF (IBH.GT.0.AND.JBH.GT.0) THEN
*       Exclude free-floating component or hierarchy.
          IF (JBH.GT.2*NPAIRS.OR.IABS(IBH-JBH).NE.1) GO TO 300
          IPAIR = KVEC(IBH)
          SEMI = -0.5*BODY(N+IPAIR)/H(IPAIR)
          ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
     &                           TDOT2(IPAIR)**2/(BODY(N+IPAIR)*SEMI)
          ECC = SQRT(ECC2)
          EB = -0.5*BODY(IBH)*BODY(JBH)/SEMI
          IPN = -1
          NP = LIST(1,2*IPAIR-1)
          WRITE (57,210)  IPN, TIME+TOFF, ECC, SEMI, EB, NP,GAMMA(IPAIR)
  210     FORMAT (' BBH   IP T E A EB NP G ',
     &                    I3,F12.5,F9.5,1P,E12.4,0P,F12.6,I3,1P,E9.1)
          CALL FLUSH(57)
      END IF
*
  300 RETURN
*
      END
