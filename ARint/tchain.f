      SUBROUTINE TCHAIN(ISUB,TSMIN)
*
*
*       Time interval for next chain c.m.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
*
*
*       Set the maximum chain integration step to c.m. value.
      TSMIN = STEP(ICH)
*       Include the block-step as safety check (exclude first time).
      IF (TBLOCK.GT.TPREV) THEN
          TSMIN = MIN(TSMIN,TBLOCK - TPREV)
      END IF
*
*       Specify next interval unless termination.
      IF (STEPS(ISUB).GT.0.0D0) THEN
          STEPS(ISUB) = TSMIN
      END IF
*
      RETURN
*
      END
