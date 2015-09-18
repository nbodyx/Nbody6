      SUBROUTINE NBREM(ICM,NSYS,NP)
*
*
*       Removal of ghosts from neighbour lists.
*       ---------------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Examine all members inside neighbour radius of body #ICM.
      DO 100 LL = 1,NP
          I = JPERT(LL) 
   10     NNB1 = LIST(1,I) + 1
          IF (NNB1.EQ.1) GO TO 100
*
*       First see whether body #ICM is a neighbour.
          DO 20 L = 2,NNB1
              IF (LIST(L,I).EQ.ICM) THEN
                  GO TO 30
              ELSE
                  IF (LIST(L,I).GT.ICM) GO TO 100
              END IF
   20     CONTINUE
*
*       Remove any other members of the subsystem.
   30     DO 60 K = 1,NSYS
              J = JLIST(K)
              IF (J.EQ.ICM) GO TO 60
*
*       Determine location of body #J.
              DO 50 L = 2,NNB1
                  IF (LIST(L,I).EQ.J) THEN
*       Reduce membership and move all subsequent members up by one.
                      LIST(1,I) = LIST(1,I) - 1
                      NNB1 = NNB1 - 1
                      DO 40 LJ = L,NNB1
                          LIST(LJ,I) = LIST(LJ+1,I)
   40                 CONTINUE
                      GO TO 60
                  ELSE
                      IF (LIST(L,I).GT.J) GO TO 60
                  END IF
   50         CONTINUE
   60     CONTINUE
*
*       Also check any KS perturber list.
          IF (I.GT.N) THEN
              I = 2*(I - N) - 1
              IF (LIST(1,I).GT.0) GO TO 10
          END IF
  100 CONTINUE
*
      RETURN
*
      END
