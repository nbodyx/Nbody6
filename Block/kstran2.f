      SUBROUTINE KSTRAN2(UI,UIDOT,Q1,Q2,Q3,RDOT)
*
      INCLUDE 'common6.h'
      REAL *8 UI(4),UIDOT(4),RDOT(3),A1(3,4)
*
*       Form relative coordinates obtained from explicit KS transformation.
      Q1 = UI(1)**2 - UI(2)**2 - UI(3)**2 + UI(4)**2
      Q2 = UI(1)*UI(2) - UI(3)*UI(4)
      Q3 = UI(1)*UI(3) + UI(2)*UI(4)
      Q2 = Q2 + Q2
      Q3 = Q3 + Q3
      RI = UI(1)**2+UI(2)**2+UI(3)**2+UI(4)**2
*
*       Set current transformation matrix.
      CALL MATRIX(UI,A1)
*
*       Obtain relative velocity from KS transformation.
      RINV = 2.0D0/RI
      DO 30 L = 1,3
          RDOT(L) = 0.0D0
          DO 25 K = 1,4
              RDOT(L) = RDOT(L) + A1(L,K)*UIDOT(K)
   25     CONTINUE
          RDOT(L) = RDOT(L)*RINV
   30 CONTINUE
*
      RETURN
*
      END
