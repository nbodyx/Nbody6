      SUBROUTINE CIRC
*
*
*       Circularization of inner binary.
*       --------------------------------
*
      INCLUDE 'common6.h'
      COMMON/KOZTMP/  EMAX
*
*
*       Determine inner eccentricity after merger termination.
      IP = NPAIRS
      I1 = 2*IP - 1
      A0 = -0.5*BODY(N+IP)/H(IP)
      ECC2 = (1.0 - R(IP)/A0)**2 + TDOT2(IP)**2/(BODY(N+IP)*A0)
      ECC0 = SQRT(ECC2)
*
*       Adopt circularization based on constant angular momentum.
*     SEMI = A0*(1.0 - ECC2)
      SEMI = A0*(1.0 - EMAX**2)
      SEMI0 = SEMI
      RX = MAX(RADIUS(I1),RADIUS(2*IP))
*       Impose a limit to avoid Roche case (but note BRAKE2 reduction).
      SEMI = MAX(SEMI,3.0*RX)
      ZMU = BODY(I1)*BODY(2*IP)/BODY(N+IP)
      HI = H(IP)
*       Set new binding energy and update ECOLL for conservation.
      H(IP) = -0.5*BODY(N+IP)/SEMI
      ECOLL = ECOLL + ZMU*(HI - H(IP))
*
*       Transform to new semi-major axis and small eccentricity (KSTAR=10).
      CALL EXPAND(IP,A0)
      CALL KSPERI(IP)
      ECC1 = 0.001
      CALL DEFORM(IP,ECC0,ECC1)
      KSTAR(N+IP) = MAX(10,KSTAR(N+IP))
      TEV0(N+IP) = TIME
      TEV(N+IP) = TIME + 0.1  ! important for exceeding TEV0.
*
*       Estimate circularization time in Myr (C. Tout, Cambody lecture 2006).
      IF (RADIUS(2*IP).GT.RADIUS(I1)) THEN
          J1 = I1 + 1
          J2 = I1
      ELSE
          J1 = I1
          J2 = I1 + 1
      END IF
      Q1 = BODY(J1)/BODY(J2)
      TCIRC = 2.0*Q1**2/(1.0 + Q1)*(SEMI/RADIUS(J1))**8
      TCIRC = 1.0D-06*TCIRC
*
      WRITE (6,15)  TIME+TOFF, NAME(I1), ECC0, A0, SEMI, TCIRC
      WRITE (45,15)  TIME+TOFF, NAME(I1), ECC0, A0, SEMI, TCIRC
   15 FORMAT (' CIRCULARIZED    T NM E A0 A TC ',F9.2,I6,F8.4,1P,3E10.2)
      CALL FLUSH(45)
*
*       Check standard collision condition.
      IF (RADIUS(I1) + RADIUS(2*IP).GT.SEMI) THEN
          KSPAIR = NPAIRS
          WRITE (6,20) KSPAIR, ECC0, A0, SEMI, RX
   20     FORMAT (' KOZAI COAL/COLL    KS E0 A0 A RX ',
     &                                 I4,F8.4,1P,3E10.2)
          IQCOLL = -2
          CALL CMBODY(R(IP),2)
      END IF
*
      RETURN
*
      END
