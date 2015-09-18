      SUBROUTINE BHPLOT
*
*
*       Black hole plotting data (1 or 2).
*       ----------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  XI(8)
*
*
*       Search for one or two BHs (NAME = 1 or 2).
      IBH = 0
      JBH = 0
      DO 5 I = 1,N
          IF (NAME(I).EQ.1.AND.IBH.EQ.0) IBH = I
          IF (NAME(I).EQ.2.AND.IBH.GT.0) JBH = I
    5 CONTINUE
*
*       Copy coordinates and velocities (singles or c.m. for simplicity).
      I = IBH
      KK = 0
*       Switch to second component if IBH not identified.
      IF (I.EQ.0) THEN
          I = JBH
          KK = 4
      END IF
   10 IF (I.GE.IFIRST) THEN
          DO 12 K = 1,3
              XI(KK+K) = X(K,I)
   12     CONTINUE
      ELSE IF (I.GT.0) THEN
          IPAIR = KVEC(I)
          DO 15 K = 1,3
              XI(KK+K) = X(K,N+IPAIR)
   15     CONTINUE
      END IF
*       Include central distance as 4th entry.
      XI(KK+4) = SQRT(XI(KK+1)**2 + XI(KK+2)**2 + XI(KK+3)**2)
*
*       Repeat for possible second BH.
      IF (I.EQ.IBH.AND.JBH.GT.0) THEN
          I = JBH
          KK = KK + 4
          GO TO 10
      END IF
*
*       Produce plotting output in Myr on unit 45 (to match opotion #45).
      IF (IBH.GT.0) THEN
          WRITE (45,20)  (TIME+TOFF)*TSTAR, (XI(K),K=KK+1,KK+4)
   20     FORMAT (' ',F10.4,4F9.4)
          CALL FLUSH(45)
      END IF
*
      RETURN
*
      END
