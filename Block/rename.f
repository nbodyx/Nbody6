      SUBROUTINE RENAME
*
*
*       Renaming of list arrays.
*       ------------------------
*
      INCLUDE 'common6.h'
*
*
      I1 = 2*NPAIRS-1
      I2 = 2*NPAIRS
*
*
*       Form the sequential names of exchanged & regularized bodies.
*       NL = total number of renaming bodies
      IF (ICOMP.EQ.I1) THEN
*       Skip correctly placed
          IF (JCOMP.EQ.I2) THEN
              NL = 0
              GO TO 41
          END IF
          NL = 3
          ILIST(1) = I1
          ILIST(2) = I2
          ILIST(3) = JCOMP
*       I2 <=> JCOMP; I1 is included to correctly calculate NIJ
          ILIST(5) = ILIST(1)
          ILIST(6) = ILIST(3)
          ILIST(7) = ILIST(2)
      ELSE IF (ICOMP.EQ.I2) THEN
          NL = 3
          ILIST(1) = I1
          ILIST(2) = I2
          ILIST(3) = JCOMP
*       I1 => JCOMP, I2 => I1, JCOMP => I2
          ILIST(5) = ILIST(3)
          ILIST(6) = ILIST(1)
          ILIST(7) = ILIST(2)
      ELSE
          NL = 4
          ILIST(1) = I1
          ILIST(2) = I2
          ILIST(3) = ICOMP
          ILIST(4) = JCOMP
*       I1 <=> ICOMP, I2 <=> JCOMP
          ILIST(5) = ILIST(3)
          ILIST(6) = ILIST(4)
          ILIST(7) = ILIST(1)
          ILIST(8) = ILIST(2)
      END IF
*
*       Rename neighbour lists with exchanged or regularized components.
      DO 40 J = 1,NTOT-1
*       See whether any neighbours need to be renamed.
          NNB = LIST(1,J)
          IF (NNB.EQ.0) GO TO 40
*       Skip modification if first neighbour comes after second component.
          IF (LIST(2,J).GT.ILIST(NL)) GO TO 40
*       No special precaution needed for case of no neighbours.
          NNB1 = 0
          NIJ = 0
          LL = 2
          DO 25 K = 1,NL
              IK = ILIST(K)
              LR = NNB + 1
              IF (IK.LT.LIST(LL,J)) GO TO 25
              IF (IK.GT.LIST(LR,J)) GO TO 30
*       Binary search
   20         IF (LL.GT.LR) GO TO 25
              LC = (LL+LR)/2
              JC = LIST(LC,J)
              IF (JC.LT.IK) THEN
                  LL = LC + 1
                  GO TO 20
              ELSEIF (JC.GT.IK) THEN
                  LR = LC - 1
                  GO TO 20
              ELSE
                  IF (IK.NE.ILIST(K+4)) THEN
                      NNB1 = NNB1 + 1
*       Save corresponding list location for the modification loop.
                      JLIST(NNB1) = K
                      JLIST(NNB1+4) = LC
                  END IF
*       Count number of identified KS components.
                  IF (IK.EQ.ICOMP.OR.IK.EQ.JCOMP) NIJ = NIJ + 1
                  IF (LC.EQ.NNB+1) GO TO 30
                  LL = LC + 1
              ENDIF
   25     CONTINUE
*
*       Skip modification if no neighbours need to be renamed.
   30     IF (NNB1.EQ.0) GO TO 40
*
*       Include c.m. renaming procedure if both components are neighbours.
          DO 38 K = 1,NNB1
              KTIME = JLIST(K)
              LS = JLIST(K+4)
*       Indicator for identified members of search list.
              IOLD = ILIST(KTIME)
              INEW = ILIST(KTIME+4)
*       Replace identified components (or single perturber) by c.m.
              IF (INEW.LT.IFIRST) THEN
                  IF (NIJ.EQ.2.OR.J.LT.IFIRST - 3) INEW = NTOT
              END IF
*       Start with the saved list index.
              LL = LS
              IF (INEW.LT.IOLD) THEN
*       Move list members down by one until location for new name is vacated.
                  DO 32 L = LL,3,-1
                      IF (LIST(L-1,J).GT.INEW) THEN
                          LIST(L,J) = LIST(L-1,J)
                          LS = L - 1
                      ELSE
                          GO TO 36
                      END IF
   32             CONTINUE
              ELSE
*       Move list members up by one until sequential location is reached.
                  DO 34 L = LL,NNB
                      IF (LIST(L+1,J).LT.INEW) THEN
                          LIST(L,J) = LIST(L+1,J)
                          LS = L + 1
                      ELSE
                          GO TO 36
                      END IF
   34             CONTINUE
              END IF
*       Set renamed neighbour in sequential location.
   36         LIST(LS,J) = INEW
*       Update positions
              DO 37 L = K+1,NNB1
                  IF (JLIST(L+4).LE.LS)
     &                JLIST(L+4) = JLIST(L+4) - 1
   37         CONTINUE
   38     CONTINUE
*
*       Reduce membership by one if c.m. set in last two locations.
          IF (NIJ.EQ.2) LIST(1,J) = NNB - 1
   40 CONTINUE
*
*       Update the list of recently regularized particles.
   41 NNB = LISTR(1)
      DO 44 L = 2,NNB+1
*       First see whether either component has been regularized before.
          IF (LISTR(L).EQ.ICOMP.OR.LISTR(L).EQ.JCOMP) THEN 
*       Remove corresponding pair even if only one component is present.
              J = 2*KVEC(L-1)
*       Move up the subsequent pairs and reduce membership by two.
              DO 42 K = J,NNB-1
                  LISTR(K) = LISTR(K+2)
   42         CONTINUE
              LISTR(1) = LISTR(1) - 2
*       Make a new search otherwise LISTR -> -2 if NNB = 2.
              GO TO 41
          END IF
   44 CONTINUE
*
*       Rename exchanged components if present in the list.
      IF (NL.GT.0) THEN
          DO 46 L = 2,NNB+1
              DO 45 K = 1,2
                  IF (LISTR(L).EQ.ILIST(K)) THEN
                      LISTR(L) = ILIST(K+4)
                      GO TO 46
                  END IF
   45         CONTINUE
   46     CONTINUE
      END IF
*
*       Update the list of high velocity particles.
      IF (LISTV(1).EQ.0) GO TO 70
      NNB = LISTV(1)
      L = 1
   50 L = L + 1
*       Check for removal of regularized component.
      IF (LISTV(L).EQ.ICOMP.OR.LISTV(L).EQ.JCOMP) THEN
          DO 52 K = L,NNB
              LISTV(K) = LISTV(K+1)
   52     CONTINUE
          LISTV(1) = LISTV(1) - 1
          NNB = NNB - 1
*       Consider the same location again.
          L = L - 1
      END IF
      IF (L.LE.NNB) GO TO 50
*
*       Rename exchanged components.
      IF (NL.GT.0) THEN
          DO 56 L = 2,NNB+1
              DO 55 K = 1,2
                  IF (LISTV(L).EQ.ILIST(K)) THEN
                      LISTV(L) = ILIST(K+4)
                      GO TO 56
                  END IF
   55         CONTINUE
   56     CONTINUE
      END IF
*
*       Remove any fast particles which have slowed down or are outside 2<R>.
      L = 1
   60 L = L + 1
      IF (LISTV(1).EQ.0) GO TO 70
      J = LISTV(L)
      A0 = XDOT(1,J)**2 + XDOT(2,J)**2 + XDOT(3,J)**2
      A2 = (X(1,J) - RDENS(1))**2 + (X(2,J) - RDENS(2))**2 +
     &                              (X(3,J) - RDENS(3))**2
      IF (A0.LT.16.0*ECLOSE.OR.A2.GT.4.0*RSCALE**2) THEN 
          DO 65 K = L,NNB
              LISTV(K) = LISTV(K+1)
   65     CONTINUE
          LISTV(1) = LISTV(1) - 1
          NNB = NNB - 1
*       Consider the same location again.
          L = L - 1
      END IF
      IF (L.LE.NNB) GO TO 60
*
   70 RETURN
*
      END
