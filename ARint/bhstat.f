      SUBROUTINE BHSTAT
*
*
*       Black hole statistics.
*       ----------------------
*
      INCLUDE 'common6.h'
      PARAMETER (NMX=10)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CX(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/POSTN2/ SEMIGR,ECCGR,DEGR,ISPIN
      INTEGER NAMEB(10),IS(10)
      REAL*8  CM(3),EI(10)
*
*
*       Identify heavy single or binary in case of no chain.
      IF (NCH.EQ.0) THEN
          ICH = 0
          DO 3 I = 1,NTOT
          IF (BODY(I).GT.0.0001) ICH = I
    3     CONTINUE
          IF (ICH.EQ.0) GO TO 100
      END IF
*
*       Use c.m. as reference point unless case of strong recoil.
      RICH = SQRT((X(1,ICH)-RDENS(1))**2 + (X(2,ICH)-RDENS(2))**2 +
     &                                     (X(3,ICH)-RDENS(3))**2)
      IF (RICH.LT.0.2.AND.NCH.GT.0) THEN
          DO 1 K = 1,3
              CM(K) = X(K,ICH)
    1     CONTINUE
      ELSE
          DO 2 K = 1,3
              CM(K) = 0.0
    2     CONTINUE
      END IF
*
*       Define names of all BHs.
      Do 5 L = 1,NBH0
          NAMEB(L) = L
    5 CONTINUE
*
*       Determine global location but skip ghosts.
      LX = 0
      DO 10 I = IFIRST,N
          IF (BODY(I).EQ.0.0D0) GO TO 10
          DO 8 L = 1,NBH0
              IF (NAME(I).EQ.NAMEB(L)) THEN
                  LX = LX + 1
                  IS(LX) = I
                  IF (LX.EQ.NBH0) GO TO 12
                  GO TO 10
              END IF 
    8     CONTINUE
   10 CONTINUE
*
*       Evaluate binding energy of at most LX single BHs and chain c.m.
   12 IF (LX.EQ.0) GO TO 100
      LX = LX + 1
*       Ignore extra case if NCH = 0.
      IS(LX) = ICH
      DO 20 L = 1,LX
          I = IS(L)
          POTJ = 0.0D0
          DO 16 J = IFIRST,NTOT
              IF (J.EQ.I) GO TO 16
              RIJ2 = 0.0D0
              DO 15 K = 1,3
                  RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
   15         CONTINUE
              POTJ = POTJ + BODY(J)/SQRT(RIJ2)
   16     CONTINUE
          VI2 = 0.0
          DO 18 K = 1,3
              VI2 = VI2 + XDOT(K,I)**2
   18     CONTINUE
*       Note chain c.m. added for L = LX.
          EI(L) = 0.5*VI2 - POTJ
          IF (KZ(14).EQ.1) THEN
              CALL XTRNLV(I,I)
              EI(L) = EI(L) + HT/BODY(I)
          END IF
   20 CONTINUE
*
*       Save central distance and binding energy for single BHs.
      DO 30 L = 1,LX-1
          I = IS(L)
          RI = SQRT((X(1,I) - CM(1))**2 + (X(2,I) - CM(2))**2 +
     &                                    (X(3,I) - CM(3))**2)
          RD = (X(1,I) - CM(1))*XDOT(1,I) + (X(2,I) - CM(2))*XDOT(2,I)
     &                                    + (X(3,I) - CM(3))*XDOT(3,I)
          RD = RD/RI
          RIJ2 = 0.0
          VIJ2 = 0.0
          RDOT = 0.0
          DO 22 K = 1,3
              RIJ2 = RIJ2 + (X(K,I) - X(K,ICH))**2
              VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,ICH))**2
              RDOT = RDOT + (X(K,I) - X(K,ICH))*(XDOT(K,I)-XDOT(K,ICH))
   22     CONTINUE
          RIJ = SQRT(RIJ2)
          SEMI = 2.0/RIJ - VIJ2/(BODY(I) + BODY(ICH))
          SEMI = 1.0/SEMI
          ECC2 = (1.0-RIJ/SEMI)**2 + RDOT**2/(SEMI*(BODY(I)+BODY(ICH)))
          ECC = SQRT(ECC2)
          IF (RI.GT.100.0*RMIN) ECC = 0.0
          PMIN = SEMI*(1.0 - ECC)/SEMIGR
          IF (RI.GT.100.0*RMIN) PMIN = 0.0
          IF (NAME(I).EQ.1) THEN
              WRITE (91,25)  TTOT, ECC, PMIN, RI, EI(L)
   25         FORMAT (' ',F8.2,F7.3,F6.1,1P,2E10.2)
              CALL FLUSH(91)
          ELSE IF (NAME(I).EQ.2) THEN
              WRITE (92,25)  TTOT, ECC, PMIN, RI, EI(L)
              CALL FLUSH(92)
          ELSE IF (NAME(I).EQ.3) THEN
              WRITE (93,25)  TTOT, ECC, PMIN, RI, EI(L)
              CALL FLUSH(93)
          ELSE IF (NAME(I).EQ.4) THEN
              WRITE (94,25)  TTOT, ECC, PMIN, RI, EI(L)
              CALL FLUSH(94)
          END IF
   30 CONTINUE
*
*       Assign small values inside the chain for plotting purposes (log!).
      ECC = 0.0
      PMIN = 0.0
      DO 40 L = 1,NCH
          RI = 1.0D-04
          IF (NAMEC(L).EQ.1) THEN
              WRITE (91,25)  TTOT, ECC, PMIN, RI, EI(LX)
              CALL FLUSH(91)
          ELSE IF (NAMEC(L).EQ.2) THEN
              WRITE (92,25)  TTOT, ECC, PMIN, RI, EI(LX)
              CALL FLUSH(92)
          ELSE IF (NAMEC(L).EQ.3) THEN
              WRITE (93,25)  TTOT, ECC, PMIN, RI, EI(LX)
              CALL FLUSH(93)
          ELSE IF (NAMEC(L).EQ.4) THEN
              WRITE (94,25)  TTOT, ECC, PMIN, RI, EI(LX)
              CALL FLUSH(94)
          END IF
   40 CONTINUE
*
*       Include data for any chain c.m.
      IF (NCH.GT.0) THEN
          VI2 = XDOT(1,ICH)**2 + XDOT(2,ICH)**2 + XDOT(3,ICH)**2
          WRITE (95,50)  TTOT, ECC, PMIN, RICH, EI(LX), VI2
   50     FORMAT (' ',F8.2,F7.3,F6.1,1P,2E10.2,0P,F7.3)
          CALL FLUSH(95)
      END IF
*
  100 RETURN
*
      END
