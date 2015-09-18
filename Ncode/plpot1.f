      SUBROUTINE PLPOT1(PHI)
*
*
*       Plummer potential energy.
*       -------------------------
*
      INCLUDE 'common6.h'
*
*
*       Evaluate tidal energy for Plummer potential.
          PHI = 0.0
          DO 10 I = IFIRST,NTOT
              RI2 = AP2
              DO 5 K = 1,3
                  RI2 = RI2 + X(K,I)**2
    5         CONTINUE
              PHI = PHI - BODY(I)*MP/SQRT(RI2)
   10 CONTINUE
*
      RETURN
*
      END
