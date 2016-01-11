      SUBROUTINE KICK(I,ICASE,KW,DM)
*
*
*       Velocity kick for WD, neutron stars or black holes.
*       ---------------------------------------------------
*
* There are various choices that can be made for kicks. 
* Make sure you are comfortable with the choices set (see below) 
* as this will critically affect retention statistics. 
*
* For regular NSs the kick choice is controlled by the value of 
* DISP (sigma in BSE). Choices are: 
*    DISP < 0.0 - kick drawn randomly between 0 - ABS(DISP) km/s
*    DISP = 0.0 - no kick
*    DISP > 0.0 - kick drawn from a Maxwellian with dispersion DISP. 
* You may also choose to have the kick velocity distribution scaled 
* by VSTAR (i.e. scaled by the initial escape velocity). 
* To do this set VFAC to a non-zero value and VFAC*VSTAR will be 
* either the maximum of the flat distribution (DISP < 0) or 
* the dispersion of the Maxwellian (DISP > 0). 
*
* Then for an electron capture supernova or an accretion-induced 
* collapse the choice is determined by the value of ECSIG set 
* internally here. Choices are: 
*    ECSIG = 0.0 - no kick
*    ECSIG > 0.0 - kick drawn from a Maxwellian, dispersion ECSIG 
*    ECSIG < 0.0 - same as for the regular NSs but scaled by ABS(ECSIG).
* These supernova are identified by their mass of 1.26 Msun. 
*
* For BHs the kick choice is controlled by the value of BHFLAG.  
* Choices are: 
*    BHFLAG = 0 - no kick
*    BHFLAG = 1 - same as for the regular NSs
*    BHFLAG = 2 - same as for the regular NSs but scaled by fallback. 
*
* Small kicks for WDs can also be set if KZ(25) > 0 in the input file. 
* In this case you can distinguish: 
*   WDSIG1 - He and COWDs
*   WDSIG2 - ONeWDs 
* as the dispersion in a Maxwellian for the different WD types. 
* A limit of WDKMAX is set. 
* See Fellhauer et al. 2003, ApJ, 595, L53 for more on WD kicks. 
*
      INCLUDE 'common6.h'
      INTEGER BHFLAG
      REAL*8  RAN2, VK(4)
      SAVE  IPAIR, KC, VDIS, RI
      DATA  IPAIR,KC /0,0/
      REAL*8 DISP, VFAC, ECSIG, WDSIG1, WDSIG2, WDKMAX
      REAL*8 DISP0, VK2
      LOGICAL IFLAT
*
*       Choose the kick settings. 
*       Some suggested combinations are included: 
*
*       Adopt the Maxwellian of Hansen & Phinney (MN 291, 569, 1997) 
*       for regular NSs, with EC kicks from a Maxwellian with a lower 
*       peak and BH kicks scaled by fallback. WD kicks depend on #25. 
*     DISP = 190.D0
*     VFAC = 0.D0
*     ECSIG = 20.D0
*     WDSIG1 = 2.D0
*     WDSIG2 = 2.D0
*     WDKMAX = 6.D0
*     BHFLAG = 2
*
*       Take a flat distribution between 0-100 km/s for regular NSs, 
*       scale EC kicks down by a factor of 4 and do not give kicks 
*       to BHs. WD kicks depend on #25. 
      DISP = -100.D0
      VFAC = 0.D0
      ECSIG = -0.25D0
      WDSIG1 = 2.D0
      WDSIG2 = 2.D0
      WDKMAX = 6.D0
      BHFLAG = 0
      BHFLAG = 2
      DISP = 20.0
*
*       Save orbital parameters in case of KS binary (called from KSAPO).
      IF (ICASE.EQ.0) THEN
          IPAIR = I
          KC = KSTAR(N+IPAIR)
*       Identify the correct component (KSTAR reversed in MDOT or EXPEL).
          I1 = 2*IPAIR - 1
          IF (KSTAR(I1).LE.0) THEN
              IN = I1
          ELSE IF (KSTAR(I1+1).LE.0) THEN
              IN = I1 + 1
          END IF
          KSTAR(IN) = -KSTAR(IN)
*
*       Determine mass loss and actual disruption velocity.
*         DM = BODY(IN) - 1.4/ZMBAR
*         IF (KW.LT.13) DM = 0.0
*         VD2 = 2.0*(BODY(N+IPAIR) - DM)/R(IPAIR)
*       Determine disruption velocity after mass-loss.
          VD2 = 2.0*BODY(N+IPAIR)/R(IPAIR)
          VDIS = SQRT(VD2)*VSTAR
*       Set cluster escape velocity (add twice central potential).
          VESC = SQRT(VD2 + 4.0)*VSTAR
          SEMI = -0.5*BODY(N+IPAIR)/H(IPAIR)
          ZM1 = BODY(I1)*SMU
          ZM2 = BODY(I1+1)*SMU
*       Skip case of massless binary.
          IF (BODY(N+IPAIR).GT.0.0D0) THEN
              EB = BODY(I1)*BODY(I1+1)/BODY(N+IPAIR)*H(IPAIR)
          ELSE
              EB = 0.0
          END IF
          RI = R(IPAIR)
*       Skip on #25 = 0/1 for consistent net WD modification of EKICK.
          IF (KW.LT.13.AND.KZ(25).EQ.0) GO TO 30
*       Sum whole binding energy (used by BINOUT for net change).
          EKICK = EKICK + EB
          EGRAV = EGRAV + EB
          I2 = I1 + 1
          WRITE (6,1)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2), ZM1,
     &                 ZM2, VESC, VDIS, R(IPAIR)/SEMI, EB, R(IPAIR)
    1     FORMAT (' BINARY KICK:    NAM K* M1 M2 VESC VDIS R/A EB R ',
     &                              2I6,2I4,4F7.1,F6.2,F9.4,1P,E10.2)
          NBKICK = NBKICK + 1
*       Remove any circularized KS binary from the CHAOS/SYNCH table.
          IF (KSTAR(N+IPAIR).GE.10.AND.NCHAOS.GT.0) THEN
              II = -(N + IPAIR)
              CALL SPIRAL(II)
              KSTAR(N+IPAIR) = 0
          END IF
          GO TO 30
      END IF
*
*       Generate velocity kick for neutron star (Gordon Drukier Tokyo paper).
*     IT = 0
*     V0 = 330.0
*     VCUT = 1000.0
*   2 VT = VCUT/V0*RAN2(IDUM1)
*     VP = VT*(2.0*RAN2(IDUM1) - 1.0)
*     VN = SQRT(VT**2 - VP**2)
*     FAC = 1.0/0.847*VN**0.3/(1.0 + VN**3.3)
*     IF (FAC.LT.RAN2(IDUM1).AND.IT.LT.10) GO TO 2
*     VKICK = V0*VT
*
      ZM = BODY(I)*ZMBAR
      VKICK = 0.D0
*
      IFLAT = .FALSE.
      IF(DISP.LT.-0.01)THEN
         IF(KW.EQ.13.AND.ZM.GE.1.28) IFLAT = .TRUE.
         IF(KW.EQ.13.AND.ZM.LT.1.28.AND.ECSIG.LT.-0.01) IFLAT = .TRUE.
         IF(KW.EQ.14.AND.BHFLAG.GT.0) IFLAT = .TRUE.
      ENDIF
*
      IF(IFLAT)THEN
*
*       Generate the kick velocity from a flat distribution. 
         DISP0 = ABS(DISP)
         IF(VFAC.GT.0.001D0) DISP0 = VFAC*VSTAR
         IF(KW.EQ.13.AND.ZM.LT.1.28) DISP0 = DISP0*ABS(ECSIG)
*
         VKICK = RAN2(IDUM1)*DISP0
         THETA = RAN2(IDUM1)*TWOPI
         SPHI = RAN2(IDUM1)
         X1 = ASIN(SPHI)
         CPHI = COS(X1)
         VK(1) = COS(THETA)*CPHI*VKICK
         VK(2) = SIN(THETA)*CPHI*VKICK
         VK(3) = SPHI*VKICK
         VK2 = VKICK*VKICK
*
      ELSE
*
*       Generate the kick velocity using a Maxwellian distribution. 
         DISP0 = MAX(DISP,0.D0)
         IF(VFAC.GT.0.001D0) DISP0 = VFAC*VSTAR
         IF(KW.EQ.10.OR.KW.EQ.11) DISP0 = MAX(WDSIG1,0.D0)
         IF(KW.EQ.12) DISP0 = MAX(WDSIG2,0.D0)
         IF(KW.EQ.13.AND.ZM.LT.1.28)THEN
            IF(ECSIG.LT.-0.01)THEN
               DISP0 = DISP0*ABS(ECSIG)
            ELSE
               DISP0 = MAX(ECSIG,0.D0)
            ENDIF
         ENDIF
         IF(KW.EQ.14.AND.BHFLAG.EQ.0) DISP0 = 0.D0
*
*       Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
         DO 2 K = 1,2
             X1 = RAN2(IDUM1)
             X2 = RAN2(IDUM1)
*       Generate two velocities from polar coordinates S & THETA.
             S = DISP0*SQRT(-2.0*LOG(1.0 - X1))
             THETA = TWOPI*X2
             VK(2*K-1) = S*COS(THETA)
             VK(2*K) = S*SIN(THETA)
    2    CONTINUE
*
         IF(DISP0.GT.0.001D0)THEN
            VK2 = VK(1)**2 + VK(2)**2 + VK(3)**2
            VKICK = SQRT(VK2)
         ELSE
            VK2 = 0.D0
            VKICK = 0.D0
         ENDIF
*
      ENDIF
*
      VK(4) = VKICK
*
* Impose the maximum WD kick velocity. 
*
      IF(KW.GE.10.AND.KW.LE.12.AND.VKICK.GT.WDKMAX)THEN
         VKICK = WDKMAX
      ENDIF
*
* Restrict the BH kick velocity by fallback. 
* This could be done better but in the N-body code we only have 
* limited information. 
*
      IF(KW.EQ.14.AND.BHFLAG.GT.1)THEN
         VKICK = VK(4)*DM/(BODY(I)+DM)
      ENDIF
*
*       Limit kick velocity to VDIS+10*VSTAR/10*VST.
*       (disabled)
      IF (IPAIR.GT.0) THEN
          VBF = SQRT(VDIS**2 + 100.0*VSTAR**2)
*       Include large kick velocity to ensure escape of disrupted star.
          IF (KW.LT.10.OR.(KW.LT.13.AND.KZ(25).EQ.0)) VKICK = 1.0*VBF
*         VKICK = MIN(VKICK,VBF)
*         VKICK = MAX(VKICK,VDIS+3.0*VSTAR)
      ELSE
*         VKICK = MIN(VKICK,10.0D0*VSTAR)
*       Ensure escape of massless star.
          IF (BODY(I).EQ.0.0D0) VKICK = 10.0*VSTAR
      END IF
*
      IF (VKICK.NE.VK(4)) THEN
         DO K = 1,3
            VK(K) = VK(K)*VKICK/VK(4)
         ENDDO
         VK2 = VKICK*VKICK
         VK(4) = VKICK
      END IF
*
*       Skip case of zero kick velocity.
      IF (VKICK.EQ.0.0D0.OR.DISP0.EQ.0.0D0) GO TO 30
*
*       Add truncated/full kick velocity and initialize X0DOT.
      VKICK = VKICK/VSTAR
      VI2 = 0.0
      VF2 = 0.0
      DO 10 K = 1,3
          VI2 = VI2 + XDOT(K,I)**2
          XDOT(K,I) = XDOT(K,I) + VK(K)/VSTAR
          X0DOT(K,I) = XDOT(K,I)
          VF2 = VF2 + XDOT(K,I)**2
   10 CONTINUE
*
*       Modify energy loss due to increased velocity of single particle.
      ECDOT = ECDOT - 0.5*BODY(I)*(VF2 - VI2)
      NKICK = NKICK + 1
*
*       Evaluate binary kick energy from relative velocity (for
*       diagnostics).
      IF (IPAIR.GT.0) THEN
          JP = KVEC(I)
          J = I + 1
          IF (I.EQ.2*JP) J = I - 1
          VF2 = 0.0
          DO 15 K = 1,3
              VF2 = VF2 + (XDOT(K,I) - XDOT(K,J))**2
   15     CONTINUE
          HNEW = 0.5*VF2 - (BODY(I) + BODY(J))/RI
*       Exclude colliding WDs.
          IF (BODY(I) + BODY(J).EQ.0.0D0) THEN
              EB1 = 0.0
          ELSE
              EB1 = BODY(I)*BODY(J)/(BODY(I) + BODY(J))*HNEW
          END IF
          IF (EB1.LT.0.0) THEN
              EKICK = EKICK - EB1
              EGRAV = EGRAV - EB1
          END IF
          IPAIR = 0
      END IF
*
      IF (NKICK.LT.50.OR.NAME(I).LE.2*NBIN0.OR.
     &    (KW.GE.13.AND.TTOT*TSTAR.GT.100.0)) THEN
          WRITE (6,20)  I, NAME(I), KSTAR(I), KC, BODY0(I)*ZMBAR, ZM,
     &                  SQRT(VI2)*VSTAR, VKICK*VSTAR, SQRT(VF2)*VSTAR
   20     FORMAT (' VELOCITY KICK:    I NAM K* KC* M0 M VI VK VF ',
     &                                2I6,2I4,2F7.2,3F7.1)
          KC = 0
      END IF
*
*       Highlight BH/NS velocities below 4 times rms velocity.
      IF (VKICK.LT.4.0*SQRT(0.5).AND.KW.GE.13.AND.VKICK.GT.0.05) THEN
          WRITE (6,25)  I, NAME(I), KW, VKICK*VSTAR, SQRT(VF2)*VSTAR
   25     FORMAT (' LOW KICK:    I NAM K* VK VF ',2I7,I4,2F7.2)
      END IF
*
*       Include optional list of high-velocity particles (double counting!).
*     IF (KZ(37).GT.0) THEN
*         CALL HIVEL(I)
*     END IF
*
   30 RETURN
*
      END
