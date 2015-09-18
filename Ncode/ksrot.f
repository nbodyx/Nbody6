      SUBROUTINE KSROT(U,UDOT,THETA)
*
*       KS Orbit rotation.
*       -----------------
*       Authors: D.J. Urminsky & D.C. Heggie (2005).
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 U(4),UDOT(4),UHAT(4),V(4),VHAT(4),UHAT2(4),VHAT2(4)
*
      US = SQRT(U(1)**2 + U(2)**2 + U(3)**2 + U(4)**2)
      UDOTR = 0.0
      DO 5 K = 1,4
          UHAT(K) = U(K)/US
          UDOTR = UDOTR + UDOT(K)*UHAT(K)
    5 CONTINUE
*
      VS = 0.0
      DO 10 K = 1,4
          V(K) = UDOT(K) - UDOTR*UHAT(K)
          VS = VS + V(K)**2
   10 CONTINUE
      VS = SQRT(VS)
      DO 15 K = 1,4
          VHAT(K) = V(K)/VS
   15 CONTINUE
*
*       Perform the rotation of physical angle THETA.
      DO 20 K = 1,4
          UHAT2(K) = UHAT(K)*COS(0.5D0*THETA) + VHAT(K)*SIN(0.5D0*THETA)
          VHAT2(K) = VHAT(K)*COS(0.5D0*THETA) - UHAT(K)*SIN(0.5D0*THETA)
          U(K) = US*UHAT2(K)
          UDOT(K) = UDOTR*UHAT2(K) + VS*VHAT2(K)
   20 CONTINUE
*
      RETURN
      END
