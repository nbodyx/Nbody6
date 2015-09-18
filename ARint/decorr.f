      SUBROUTINE DECORR(DE)
*
*
*       Collision/coalescence energy correction.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
*
*       Accumulate energy change for conservation purposes.
      ECOLL = ECOLL + DE
*
      RETURN
*
      END
