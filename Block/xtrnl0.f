      SUBROUTINE XTRNL0
*
*
*       External force initialization.
*       ------------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)
*
*
*       Check option for cluster in circular galactic orbit.
      IF (KZ(14).NE.1) GO TO 20
*
*       Specify Oort's constants (units of km/sec/kpc).
      A = 14.4
      B = -12.0
*       Adopt local density from Gilmore & Kuijken (solar mass/pc**3).
      RHO = 0.11
*       Convert rotation constants to units of cm/sec/pc.
      A = 100.0*A
      B = 100.0*B
*
*       Specify the tidal term in star cluster units (solar mass & pc).
      TIDAL(1) = 4.0*A*(A - B)*(PC/GM)
*
*       Initialize the Y-component to zero.
      TIDAL(2) = 0.0
*
*       Specify the vertical force gradient.
      TIDAL(3) = -(2.0*TWOPI*RHO + 2.0*(A - B)*(A + B)*(PC/GM))
*
      FAC = 1.0E-10/(PC/GM)
      WRITE (6,5)  ZMBAR*ZMASS, FAC*TIDAL(1), FAC*TIDAL(3), PC/GM
    5 FORMAT (/,12X,'TOTAL MASS =',F8.1,'  TIDAL(1&3) =',1P,2E10.2,
     &              '  PC/GM =',E10.2)
*
*       Adopt twice the angular velocity for Coriolis terms.
      TIDAL(4) = 2.0*(A - B)*SQRT(PC/GM)
*
*       Define time scale in seconds and velocity scale in km/sec.
      TSCALE = SQRT(PC/GM)*PC
*
*       Convert time scale from units of seconds to million years.
      TSCALE = TSCALE/(3.15D+07*1.0D+06)
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
      DO 10 K = 1,3
          TIDAL(K) = TIDAL(K)*RBAR**3/ZMBAR
   10 CONTINUE
      TIDAL(4) = TIDAL(4)*SQRT(RBAR**3/ZMBAR)
      TSCALE = TSCALE*SQRT(RBAR**3/(ZMASS*ZMBAR))
*
*       Perform an energy scaling iteration to include ETIDE (once is OK).
      E0 = -0.25
      CALL XTRNLV(1,N)
*       Form total energy based on current values and ETIDE (cf. SCALE).
      ETOT = ZKIN - POT + ETIDE
      SX = E0/ETOT
*       Exclude re-scaling for positive energy or binary Plummer cluster.
      IF (ETOT.GE.0.0) SX = 1.0
      IF (KZ(5).EQ.2) SX = 1.0
*
*       Scale coordinates & velocities to yield ETOT = -0.25.
      DO 15 I = 1,N
          DO 12 K = 1,3
              X(K,I) = X(K,I)/SX
              XDOT(K,I) = XDOT(K,I)*SQRT(SX)
   12     CONTINUE
   15 CONTINUE
*
*       Consider alternatives: circular point-mass orbit or 3D galaxy model.
   20 ZMTOT = ZMASS*ZMBAR
      IF (KZ(14).EQ.2) THEN
*
*       Read galaxy mass and central distance (solar units and kpc).
          READ (5,*)  GMG, RG0
*
*       Set circular velocity in km/sec and angular velocity in cgs units.
          VG0 = 1.0D-05*SQRT(GMG/(1000.0*RG0))*SQRT(GM/PC)
          OMEGA = 100.0*VG0/RG0
*
*       Obtain King tidal radius in pc (eq. (9) of Fukushige & Heggie, 1995).
          RT = (ZMTOT/(3.0*GMG))**0.3333*(1000.0*RG0)
*
          IF (RTIDE.GT.0.0) THEN
*       Determine RBAR (N-body units) from RT (pc) and King model (see SCALE).
              IF(KZ(22).EQ.2) RBAR = RT/RTIDE
          ELSE
              RTIDE = RT/RBAR
          END IF
*
*       Convert from cgs to N-body units.
          OMEGA = OMEGA*SQRT(PC/GM)*SQRT(RBAR**3/ZMBAR)
*
*       Specify the galactic parameters for equations of motion.
          TIDAL(1) = 3.0*OMEGA**2
          TIDAL(2) = 0.0D0
          TIDAL(3) = -OMEGA**2
          TIDAL(4) = 2.0*OMEGA
          GMG = GMG/ZMTOT
*
*       Check re-scaling units to current RBAR (i.e. TSTAR & VSTAR).
          IF (KZ(22).EQ.2) THEN
              CALL UNITS
          END IF
*
          WRITE (6,35)  GMG, RG0, OMEGA, RTIDE, RBAR
*
*       Perform an energy scaling iteration to include ETIDE (once is OK).
      E0 = -0.25
      CALL XTRNLV(1,N)
*       Form total energy based on current values and ETIDE (cf. SCALE).
      ETOT = ZKIN - POT + ETIDE
      SX = E0/ETOT
      IF (ETOT.GE.0.0) SX = 1.0
*
*       Scale coordinates & velocities to yield ETOT = -0.25.
      DO 24 I = 1,N
          DO 22 K = 1,3
              X(K,I) = X(K,I)/SX
              XDOT(K,I) = XDOT(K,I)*SQRT(SX)
   22     CONTINUE
   24 CONTINUE
*
*       Treat the case of 3D orbit for point-mass, bulge, disk and/or halo.
      ELSE IF (KZ(14).EQ.3) THEN
*
*       Read all parameters (NB! Do not confuse with Oort's constants A, B).
          READ (5,*)  GMG, DISK, A, B, VCIRC, RCIRC, GMB, AR, GAM
*       Gamma/eta model: GMB = mass, AR = scale radius, GAM = exponent.
          READ (5,*)  (RG(K),K=1,3), (VG(K),K=1,3)
*
*       Specify planar motion from SEMI & ECC for no disk & halo if VZ = 0.
          IF (DISK + VCIRC + GMB.EQ.0.0.AND.VG(3).EQ.0.0D0) THEN
              RAP = RG(1)
              ECC = RG(2)
              SEMI = RAP/(1.0 + ECC)
              VG2 = GMG/(1000.0*SEMI)*(1.0 - ECC)/(1.0 + ECC)
              DO 25 K = 1,3
                  RG(K) = 0.0
                  VG(K) = 0.0
   25         CONTINUE
*       Initialize 2D orbit with given eccentricity at apocentre.
              RG(1) = RAP
              VG(2) = 1.0D-05*SQRT(VG2)*SQRT(GM/PC)
          END IF
*
*       Convert from kpc and km/sec to N-body units.
          DO 30 K = 1,3
              RG(K) = 1000.0*RG(K)/RBAR
              VG(K) = VG(K)/VSTAR
   30     CONTINUE
*
*       Define the angular velocity (z-component) and mass in N-body units.
          R02 = RG(1)**2 + RG(2)**2
          OMEGA = (RG(1)*VG(2) - RG(2)*VG(1))/R02
          TIDAL(4) = 2.0*OMEGA
          GMG = GMG/ZMTOT
          GMB = GMB/ZMTOT
          AR = 1000.0*AR/RBAR
*
          IF (GMG.GT.0.0) THEN
              WRITE (6,35)  GMG, SQRT(R02), OMEGA, RTIDE, RBAR
   35         FORMAT (/,12X,'POINT-MASS MODEL:    GMG =',1P,E9.1,
     &                      '  RG =',E9.1,'  OMEGA =',E9.1,
     &                      '  RT =',0P,F6.2,'  RBAR =',F6.2)
          END IF
          IF (GMB.GT.0.0D0) THEN
              WRITE (6,36)  GMB, AR, GAM
   36         FORMAT (/,12X,'GAMMA/ETA MODEL:    GMB =',1P,E9.1,
     &                      '  AR =',E9.1,'  GAM =',E9.1)
          END IF
*
*       Define disk and/or logarithmic halo parameters in N-body units.
          IF (DISK.GT.0.0D0) THEN
              DISK = DISK/ZMTOT
              A = 1000.0*A/RBAR
              B = 1000.0*B/RBAR
              WRITE (6,40)  DISK, A, B
   40         FORMAT (/,12X,'DISK MODEL:    MD =',1P,E9.1,
     &                                   '  A =',E9.1,'  B =',E9.1)
          END IF
*
*       Determine local halo velocity from total circular velocity.
          IF (VCIRC.GT.0.0D0) THEN
              VCIRC = VCIRC/VSTAR
              RCIRC = 1000.0*RCIRC/RBAR
              A2 = RCIRC**2 + (A + B)**2
              V02 = VCIRC**2 - (GMG/RCIRC + DISK*RCIRC**2/A2**1.5)
*       Include any contribution from bulge potential (V2 = R*F).
              IF (GMB.GT.0.0D0) THEN
                  VB2 = GMB/RCIRC*(1.0 + AR/RCIRC)**(GAM-3.0)
                  V02 = V02 - VB2
              END IF
              IF (V02.LT.0.0D0) THEN
                  WRITE (6,45)  V02, 0.001*RCIRC*RBAR
   45             FORMAT (' ',' NEGATIVE HALO VELOCITY!    V02 RCIRC ',
     &                                                     1P,2E10.2)
                  STOP
              END IF
*       Specify the corresponding scale length of logarithmic halo.
              RL2 = RCIRC**2*(VCIRC**2 - V02)/V02
*       Define the asymptotic circular velocity due to halo.
              V02 = VCIRC**2
*
*       Include table of circular velocity on unit #52 (km/sec & kpc).
              RI = 1000.0/RBAR
              DR = 1000.0/RBAR
              DO 60 K = 1,30
                  RI2 = RI**2
                  A2 = RI2 + (A + B)**2
                  VB2 = GMB/RI*(1.0 + AR/RI)**(GAM-3.0)
                  VCIRC2 = GMG/SQRT(RI2) + DISK*RI2/A2**1.5 +
     &                     V02*RI2/(RL2 + RI2) + VB2
                  WRITE (52,50)  SQRT(VCIRC2)*VSTAR, RI*RBAR/1000.0
   50             FORMAT (' CIRCULAR VELOCITY:    VC R ',F7.1,F7.2)
                  RI = RI + DR
   60         CONTINUE
              CALL FLUSH(52)
*
              A2 = R02 + (A + B)**2
              VB2 = GMB/SQRT(R02)*(1.0 + AR/SQRT(R02))**(GAM-3.0)
              VCIRC2 = GMG/SQRT(R02) + DISK*R02/A2**1.5 +
     &                 V02*R02/(RL2 + R02) + VB2
              VCIRC = SQRT(VCIRC2)*VSTAR
              WRITE (6,62)  VCIRC, SQRT(R02)/1000.0, SQRT(RL2)/1000.0
   62         FORMAT (/,12X,'CIRCULAR VELOCITY:    VC RG RL',F7.1,2F7.2)
          ELSE
              V02 = 0.0
          END IF
*
*       Initialize F & FDOT of reference frame (point-mass galaxy is OK).
          CALL GCINIT
*
*       Form tidal radius from circular angular velocity (assumes apocentre).
          IF (RTIDE.EQ.0.0D0) RTIDE = (0.5/OMEGA**2)**0.3333
*
          WRITE (6,65)  (RG(K),K=1,3), (VG(K),K=1,3), SQRT(V02)
   65     FORMAT (/,12X,'SCALED ORBIT:    RG = ',3F7.2,
     &                                 '  VG = ',3F7.1,'  V0 =',F6.1)
      END IF
*
*       Include Plummer potential for 2D and 3D (use MP = 0 if not needed).
      IF (KZ(14).EQ.3.OR.KZ(14).EQ.4) THEN
*       Check input for Plummer potential.
          READ (5,*)  MP, AP2, MPDOT, TDELAY
          WRITE (6,70)  MP, AP2, MPDOT, TDELAY
   70     FORMAT (/,12X,'PLUMMER POTENTIAL:    MP =',F7.3,'  AP =',F6.2,
     &                        '  MPDOT =',F7.3,'  TDELAY =',F5.1)
          MP0 = MP
          AP2 = AP2**2
*       Rescale velocities by including the Plummer & galactic virial energy.
          IF (ZKIN.GT.0.0D0) THEN
*       Note that QVIR = Q is saved in routine SCALE and VIR < 0 with GPU.
              CALL ENERGY
              VIR = POT - VIR
              QV = SQRT(QVIR*VIR/ZKIN)
              DO 80 I = 1,N
                  DO 78 K = 1,3
                      XDOT(K,I) = XDOT(K,I)*QV
   78             CONTINUE
   80         CONTINUE
          END IF
          IF (RTIDE.EQ.0.0D0) RTIDE = 10.0*RSCALE
      ELSE
          MP = 0.0
          MPDOT = 0.0
      END IF
      RTIDE0 = RTIDE
      TSCALE = TSTAR
*
*       Define tidal radius in scaled units for linearized field.
      IF (KZ(14).LE.2) THEN
          RTIDE = (ZMASS/TIDAL(1))**0.3333
          WRITE (6,90)  (TIDAL(K),K=1,4), TSCALE, RTIDE
   90     FORMAT (/,12X,'TIDAL PARAMETERS:  ',1P,4E10.2,
     &                  '  TSCALE =',E9.2,'  RTIDE =',0P,F6.2,/)
      END IF
*
      RETURN
*
      END
