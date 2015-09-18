      SUBROUTINE BHSTAB(EB,RSUB,ISTAB)
*
*
*       BH binary stability.
*       --------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/POSTN2/ SEMIGR,ECCGR,DEGR,ISPIN
*
*
*       Look for bound pericentres inside 3*SEMIGR.
      EB = 0.0
      ZMM = 0.0
      ISTAB = 10
      NP = LISTC(1)
      DO 20 L = 2,NP+1
          J = LISTC(L)
          RIJ2 = 0.0
          VIJ2 = 0.0
          RDOT = 0.0
          DO 10 K = 1,3
              RIJ2 = RIJ2 + (X(K,ICH) - X(K,J))**2
              VIJ2 = VIJ2 + (XDOT(K,ICH) - XDOT(K,J))**2
              RDOT = RDOT + (X(K,ICH) - X(K,J))*(XDOT(K,ICH)-XDOT(K,J))
   10     CONTINUE
          RIJ = SQRT(RIJ2)
          SEMI = 2.0/RIJ - VIJ2/(BODY(ICH) + BODY(J))
          SEMI = 1.0/SEMI
          ECC2 = (1.0-RIJ/SEMI)**2 + RDOT**2/((BODY(ICH)+BODY(J))*SEMI)
          ECC = SQRT(ECC2)
          PMIN = SEMI*(1.0 - ECC)
          IF (SEMI.GT.0.0) THEN
              ALPH = 0.0
              ZMM = ZMM + BODY(ICH)*BODY(J)
*       Employ the stability criterion for dominant binary.
              NST = NSTAB(SEMIGR,SEMI,ECCGR,ECC,ALPH,BODYC(1),BODYC(2),
     &                                                        BODY(J))
              EB = EB - 0.5*BODY(ICH)*BODY(J)/SEMI
*       Define stability index.
              ISTAB = PMIN/SEMIGR
          ELSE
              NST = 0
              ISTAB = 0
          END IF
*         IF (PMIN.LT.3.0*SEMIGR.AND.ECC.LT.1.0) THEN
          IF (ECC.LT.1.0) THEN
              WRITE (77,15)  TIME+TOFF, NAME(J), NP, NST, ECCGR, ECC,
     &                       PMIN/SEMIGR
   15         FORMAT (' BHSTAB    T NM NP NST E E1 PM/AGR ',
     &                            F9.3,I7,2I4,F8.4,F7.3,F7.1)
              CALL FLUSH(77)
          END IF
   20 CONTINUE
*
*       Form characteristic size of subsystem (cf. NEWREG diagnostics).
      IF (ZMM.GT.0.0) RSUB = -ZMM/EB
*
      RETURN
*
      END
