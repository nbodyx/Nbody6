      SUBROUTINE INFALL(IBH,IESC,NBH2,ISUB)
*
*       Disruption of star by BH.
*       -------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX3=3*NMX,NMX4=4*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,M1,MASS,MC,MMIJ,XCM(3),VCM(3),CG(6)
      COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &           XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &           XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/ARCHAIN2/ MMIJ,CMX(3),CMV(3),ENERGY,EnerGR,CHTIME
      common/TIMECOMMON/Taika,timecomparison
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/ARZERO/  ISTAR0(NMX),SIZE0(NMX)
      COMMON/POSTN2/  SEMIGR,ECCGR,DEGR,ISPIN
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
      REAL*8  RVEC(3),VVEC(3)
*
*
*       Copy chain variables to standard form.
      LK = 0
      DO 4 L = 1,NCH
          DO 1 K = 1,3
              LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
    1     CONTINUE
    4 CONTINUE
*
*       Ensure c.m. values of X & XDOT are updated (just in case; 05/16).
      DO 2 K = 1,3
          X(K,ICH) = X0(K,ICH)
          XDOT(K,ICH) = X0DOT(K,ICH)
    2 CONTINUE
*
*       Obtain original c.m. kinetic energy.
      ZK1 = 0.0
      DO 5 K = 1,3
          ZK1 = ZK1 + XDOT(K,ICH)**2
    5 CONTINUE
      ZK1 = 0.5*BODY(ICH)*ZK1
*
*       Define new local c.m. for two-body system IBH & IESC.
      LX = IBH
      LK = 3*(IESC - 1)
      LN = 3*(LX - 1)
*       Choose 1/10 mass for ghost or whole star to be swallowed.
      IF (KZ(43).GE.2.AND.ISTAR(IESC).LT.10) THEN
          M1 = 0.1*M(IESC)
      ELSE
          M1 = M(IESC)
      END IF
*       Note use of XCM & VCM even if adopting direct escape.
      DO 10 K = 1,3
          LK = LK + 1
          LN = LN + 1
          XCM(K) = (M1*XCH(LK) + M(LX)*XCH(LN))/(M1 + M(LX))
          VCM(K) = (M1*VCH(LK) + M(LX)*VCH(LN))/(M1 + M(LX))
*       NB! XCH is not at actual pericentre (determined by A*(1-ECC)).
   10 CONTINUE
*
*       Ensure largest stellar type and reduce membership.
      ISTAR(LX) = MAX(ISTAR(IBH),ISTAR(IESC))
      ISTAR0(LX) = ISTAR(LX)
      NCH = NCH - 1
      NN = NCH
*
*       Check optional treatment of tidal disruption by BH.
      IF (KZ(43).GE.2.AND.ISTAR(IESC).LT.10) THEN
          NDISR = NDISR + 1
          ZMB =  BODYC(IESC) + BODYC(LX)
          M(IESC) = 0.1*M(IESC)
*       Accrete 1/10 the mass of #IESC and copy new BH variables.
          M(LX) = M(LX) + M(IESC)
          BODYC(LX) = BODYC(LX) + M(IESC)
          ZMASS = ZMASS - 9.0*M(IESC)
*       Update system masses.
          BODY(ICH) = BODY(ICH) - 9.0*M(IESC)
          MASS = MASS - 9.0*M(IESC)
          DO 12 K = 1,3
              X4(K,LX) = XCM(K)
              XDOT4(K,LX) = VCM(K)
   12     CONTINUE
*
*       Obtain external energy for ghost (internal pot in ECH).
          POT1 = 0.0
          DO 20 J = IFIRST,NTOT
              IF (J.EQ.ICH) GO TO 20
              RIJ2 = 0.0
              DO 15 K = 1,3
                  RIJ2 = RIJ2 + (X(K,ICH) - X(K,J))**2
   15         CONTINUE
              POT1 = POT1 + BODY(J)/SQRT(RIJ2)
   20     CONTINUE
*
*       Subtract ghost mass contribution to the potential energy.
          ECOLL = ECOLL - 9.0*M(IESC)*POT1
          WRITE (24,25)  TIME+TOFF, NDISR, NAMEC(IESC), ISTAR(IESC),
     &                   ECCGR, 10.0*M(IESC)*SMU, M(LX)*SMU, SEMIGR
   25     FORMAT (' DISRUPT2    T NDISR NM K* E M1 M2 SEMI ',
     &                          F8.1,I5,I7,I4,F10.6,2F6.1,1P,E10.2)
*       Note that modified M(LX) and another component will be new KS.
          CALL FLUSH(24)
*
          VR2 = 0.0
          RR2 = 0.0
          DO 30 K = 1,3
              RR2 = RR2 + (X4(K,IESC) - X4(K,LX))**2
              VI2 = VI2 + (XDOT4(K,IESC) - XDOT4(K,LX))**2
              RVEC(K) = X4(K,IESC) - X4(K,LX)
              VVEC(K) = XDOT4(K,IESC) - XDOT4(K,LX)
   30     CONTINUE
          RR = SQRT(RR2)
          SEMI = 2.0/RR - VI2/ZMB
          SEMI = 1.0/SEMI
          VA2 = ZMB/SEMI
          XFAC = SEMI/RR
          VFAC = SQRT(RR/SEMI)
          DO 35 K = 1,3
              X4(K,IESC) = XFAC*BODYC(LX)*RVEC(K)/ZMB
              X4(K,LX) = -XFAC*BODYC(IESC)*RVEC(K)/ZMB
              XDOT4(K,IESC) = VFAC*BODYC(LX)*VVEC(K)/ZMB
              XDOT4(K,LX) = -VFAC*BODYC(IESC)*VVEC(K)/ZMB
              XCM(K) = X4(K,LX)
              VCM(K) = XDOT4(K,LX)
   35     CONTINUE     
          RI2 = 0.0
          VI2 = 0.0
          DO 40 K = 1,3
              RI2 = RI2 + (X4(K,IESC) - X4(K,LX))**2
              VI2 = VI2 + (XDOT4(K,IESC) - XDOT4(K,LX))**2
   40     CONTINUE
          RIJ = SQRT(RI2)
          ANEW = 2.0/RIJ - VI2/BODYC(LX)
          ANEW = 1.0/ANEW
*         WRITE (6,45)  RVEC, RIJ, ANEW, VA2, VI2
*  45     FORMAT (' INFALL TRANSF    RVEC R A VA2 VI2 ',1P,7E10.2)
      ELSE
*       Implement complete swallowing of compact star or BH.
          M(LX) = M(LX) + M(IESC)
          BODYC(LX) = BODYC(LX) + M(IESC)
          NCOLL = NCOLL + 1
          WRITE (6,50)  NAMEC(IESC), ISTAR(IESC), BODYC(IESC)*SMU,
     &                  M(LX)*SMU
   50     FORMAT (' SWALLOWED STAR/BH    NAMC K* M1 M2 ',I6,I4,2F7.2)
          NDISR = NDISR + 1
          WRITE (24,25)  TIME+TOFF, NDISR, NAMEC(IESC), ISTAR(IESC),
     &                   ECCGR, BODYC(IESC)*SMU, M(LX)*SMU, SEMIGR
          CALL FLUSH(24)
      END IF
*
*       Check possible reduction of dominant body index.
      IF (LX.GT.IESC) LX = LX - 1
*
*       Identify global index of coalescence/escaper body.
      I = 0
      DO 55 J = IFIRST,NTOT
          IF (NAME(J).EQ.NAMEC(IESC).OR.NAME(J).EQ.0) THEN
              I = J
*       Exit on first identification (bug fix 9/16; next bit was also wrong).
              IF (NAME(I).EQ.NAMEC(IESC)) GO TO 58
          END IF
   55 CONTINUE
*
*       Switch to NAMEC(IBH) in case NAMEC(IESC) is current c.m. (rare case).
   58 IF (I.EQ.ICH) THEN
          ICH0 = ICH
*       Restore NAME(ICH) and search for NAME(IBH) as new c.m.
          NAME(I) = NAME0
          DO 60 J = IFIRST,N
              IF (NAME(J).EQ.NAMEC(IBH)) THEN
                  ICH = J
              END IF
*       Determine new global index of disrupted star.
              IF (NAME(J).EQ.NAMEC(IESC)) THEN
                  I = J
              END IF
   60     CONTINUE
*       Save new NAME(ICH) and copy BODYC(IBH) & c.m. variables to ICH.
          NAME0 = NAME(ICH)
          NAME(ICH) = 0
          BODY(ICH) = BODYC(IBH)
          DO 62 K = 1,3
              X(K,ICH) = X(K,ICH0)
              XDOT(K,ICH) = XDOT(K,ICH0)
   62     CONTINUE
      END IF
*
*       Define M(IESC) as ghost (partial or complete swallowing).
      CALL GHOST(I)
*       Ensure escape condition is satisfied (R**2 < 1D+10).
      X0(1,I) = 1.0D+04
      X(1,I) = 1.0D+04
*
*       Remove chain (and clump) mass & reference name of absorbed member.
      LI = 3*(IESC - 1)
      DO 65 L = IESC,NCH
          M(L) = M(L+1)
          BODYC(L) = BODYC(L+1)
          NAMEC(L) = NAMEC(L+1)
          SIZE(L) = SIZE(L+1)
          ISTAR(L) = ISTAR(L+1)
          SIZE0(L) = SIZE0(L+1)
          ISTAR0(L) = ISTAR0(L+1)
          BODYS(L,ISUB) = BODYS(L+1,ISUB)
          NAMES(L,ISUB) = NAMES(L+1,ISUB)
          DO 64 K = 1,3
              XCH(LI+K) = XCH(LI+K+3)
              VCH(LI+K) = VCH(LI+K+3)
   64     CONTINUE 
          LI = LI + 3
   65 CONTINUE
*
*       Set XCM & VCM in XCH & VCH of the current IBH.
      LK = 3*(LX - 1)
      DO 70 K = 1,3
          LK = LK + 1
          XCH(LK) = XCM(K)
          VCH(LK) = VCM(K)
   70 CONTINUE
*
*       Copy new chain coordinates & velocities to standard variables.
      LK = 0
      DO 80 L = 1,NCH
          DO 75 K = 1,3
              LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
   75     CONTINUE
   80 CONTINUE
*
*       Perform c.m. check (not perfect before correction).
      DO 85 K = 1,6
          CG(K) = 0.0
   85 CONTINUE
*
      LK = 0
      DO 95 L = 1,NCH
          DO 90 K = 1,3
              LK = LK + 1
              CG(K) = CG(K) + BODYC(L)*XCH(LK)
              CG(K+3) = CG(K+3) + BODYC(L)*VCH(LK)
   90     CONTINUE
   95 CONTINUE
*
      DO 96 K = 1,6
          CG(K) = CG(K)/BODY(ICH)
   96 CONTINUE
*
*       Adopt c.m. condition and copy to X4 & XDOT4.
      LK = 0
      DO 98 L = 1,NCH
          DO 97 K = 1,3
              LK = LK + 1
              XCH(LK) = XCH(LK) - CG(K)
              VCH(LK) = VCH(LK) - CG(K+3)
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
   97     CONTINUE
   98 CONTINUE
*
      IF (KZ(30).NE.0) THEN
          WRITE (6,99)  TIME+TOFF, (CG(K),K=1,6)
   99     FORMAT (' INFALL CHECK:   T CG ',F10.2,1P,6E9.1)
      END IF
*
*       Form new c.m. kinetic energy.
      ZK2 = 0.0
      DO 120 K = 1,3
          ZK2 = ZK2 + XDOT(K,ICH)**2
  120 CONTINUE
      ZK2 = 0.5*BODY(ICH)*ZK2
*
*       Add c.m. kinetic energy change for conservation.
      ECOLL = ECOLL + (ZK1 - ZK2)
*
*       Re-initialize force polynomials.
      TIME = TBLOCK
      CALL FPOLY1(ICH,ICH,0)
      CALL FPOLY2(ICH,ICH,0)
*
      RETURN
*
      END
