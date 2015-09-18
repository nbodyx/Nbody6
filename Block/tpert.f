      SUBROUTINE TPERT(IPAIR,GA,DT)
*
*
*       Perturbation time scale.
*       ------------------------
*
      INCLUDE 'common6.h'
*
*
*       Set c.m. index and initialize scalars.
      I = N + IPAIR
      FMAX = 0.0
      DTIN = 1.0E+20
      JCL = 0    ! thread-safe local variable.
      NNB1 = LIST(1,I) + 1
      SEMI = -0.5*BODY(I)/H(IPAIR)
*       Define a safe boundary for omitting intruder search using square.
      RS2 = 2.0*BODY1*CMSEP2*SEMI**2/BODY(I)
*
*       Check rare case of no neighbours.
      IF (NNB1.LE.1) THEN
          DT = 2.0*STEPR(I)
          GO TO 20
      END IF
*
*       Find the most likely perturbers (first approach & maximum force).
      DO 10 L = 2,NNB1
          J = LIST(L,I)
          RIJ2 = (X(1,J)-X(1,I))**2 + (X(2,J)-X(2,I))**2+
     &                                (X(3,J)-X(3,I))**2
          IF (RIJ2.GT.RS2) GO TO 10
          RDOT = 0.0
*
          DO 6 K = 1,3
              XREL = X(K,J) - X(K,I)
              VREL = XDOT(K,J) - XDOT(K,I)
*             RIJ2 = RIJ2 + XREL**2
              RDOT = RDOT + XREL*VREL
    6     CONTINUE
*
          VR = RDOT/RIJ2
          IF (VR.LT.DTIN) THEN
              DTIN = VR
*       Note DTIN is inverse travel time to include case of no RDOT < 0.
              RCRIT2 = RIJ2
              JCRIT = J
          END IF
          FIJ = (BODY(I) + BODY(J))/RIJ2
          IF (FIJ.GT.FMAX) THEN
              FMAX = FIJ
              RJMIN2 = RIJ2
              JCL = J
          END IF
   10 CONTINUE
*
      IF (JCL.EQ.0) THEN
          DT = 4.0D0*STEP(I)
          GO TO 15
      END IF
*
*       Form radial velocity of body with shortest approach time (if any).
      RCRIT = SQRT(RCRIT2)
      RDOT = RCRIT*ABS(DTIN)
      A1 = 2.0/(BODY(I)*GA)
*     SEMI = -0.5*BODY(I)/H(IPAIR)
*
*       Use the actual apocentre for unperturbed travel time.
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      RI = SEMI*(1.0 + SQRT(ECC2))
*
*       Estimate time interval to reach tidal perturbation of GA.
      DT = (RCRIT - RI*(BODY(JCRIT)*A1)**0.3333)/RDOT
*
*       Compare the travel time based on acceleration only.
      DTMAX = SQRT(2.0D0*ABS(DT)*RDOT*RCRIT2/(BODY(I) + BODY(JCRIT)))
      DT = MIN(DT,DTMAX)
*
*       Skip dominant force test if there is only one critical body.
      IF (JCRIT.NE.JCL) THEN
*       Form the return time of the dominant body and choose the minimum.
          DR = SQRT(RJMIN2) - RI*(BODY(JCL)*A1)**0.3333
          DTMAX = SQRT(2.0D0*ABS(DR)/FMAX)
          DT = MIN(DT,DTMAX)
      END IF
*
*       Apply safety test in case background force dominates c.m. motion.
   15 DT = MIN(DT,4.0D0*STEP(I))
      DT = MAX(DT,0.0D0)
*
*       Save JCL for use in UNPERT (thread-safe procedure).
   20 IPAIR = JCL
*
      RETURN
*
      END
