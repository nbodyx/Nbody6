      SUBROUTINE REPAIR(I1,I2,I3,I4,NQ,LISTQ,TMIN)
*
*
*      Repair of LISTQ after new or terminated KS.
*       ------------------------------------------
*
      INCLUDE 'common6.h'
      INTEGER LISTQ(NMAX)
*
*
*       Treat case of NEW KS first (I4 = IFIRST+1 for IPHASE = 2).
      IF (IPHASE.EQ.1) THEN
*       Remove original KS components ICOMP0 & JCOMP0 in turn if < IFIRST.
          II = I1
*       Skip ICOMP0 >= IFIRST and try second exchanged component or dummy.
          IF (I1.GE.IFIRST) THEN
              II = I2
              IF (I2.GE.IFIRST) II = -1
          END IF
    1     LI = 0
*       Identify particle II in LISTQ and remove unless last list member.
          DO 5 L = 1,NQ
              IF (LISTQ(L).EQ.II) LI = L
    5     CONTINUE
          IF (LI.GT.0.AND.LI.LT.NQ) THEN
              DO 10 LL = LI,NQ
                  LISTQ(LL) = LISTQ(LL+1)
   10         CONTINUE
          END IF
*       Reduce membership on LI > 0 and select second component < IFIRST.
          IF (LI.GT.0) NQ = NQ - 1
          IF (II.EQ.I1.AND.I2.LT.IFIRST) THEN
              II = I2
              GO TO 1
          END IF
*       Check removal of one or more KS components.
          DO 12 L = 1,NQ
              IF (LISTQ(L).LT.IFIRST) THEN
                  II = LISTQ(L)
                  GO TO 1
              END IF
   12     CONTINUE
*       Add current c.m. as extra member at NTOT (typically small step).
          NQ = NQ + 1
          LISTQ(NQ) = I3
      ELSE
*       Add two first single particles after KS termination (small steps).
          LISTQ(NQ+1) = I3 
          LISTQ(NQ+2) = I4
*       Note: terminated KS components usually have small steps.
          NQ = NQ + 2
*       Find and remove terminated c.m. of N + KSPAIR and reduce higher IDs.
          LI = 0
          DO 15 L = 1,NQ
              IF (LISTQ(L).EQ.I1) LI = L
*       Note LISTQ is not sequential near the end (I3 & I4 added).
              IF (LISTQ(L).GT.I1) LISTQ(L) = LISTQ(L) - 1
   15     CONTINUE
*       Compress LISTQ array by removing current c.m.
          IF (LI.GT.0.AND.LI.LT.NQ) THEN
              DO 20 LL = LI,NQ
                  LISTQ(LL) = LISTQ(LL+1)
   20         CONTINUE
          END IF
*       Reduce membership on identified c.m.
          IF (LI.GT.0) NQ = NQ - 1
      END IF
*
*       Determine new TMIN from final LISTQ.
      TMIN = 1.0D+06
      DO 30 L = 1,NQ
          J = LISTQ(L)
          TMIN = MIN(TNEW(J),TMIN)
   30 CONTINUE
*
      RETURN
*
      END
