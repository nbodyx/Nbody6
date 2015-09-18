      SUBROUTINE FFDOT(I,ZMS,XS,VS,SS)
*
*
*       Force & first derivative corrections.
*       -------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  A(9),F1(3),F1DOT(3),XS(3),VS(3)
*
*
*       Obtain force and first derivative due to body #J.
      DO 15 K = 1,3
          A(K) = XS(K) - X(K,I)
          A(K+3) = VS(K) - XDOT(K,I)
   15 CONTINUE
*
      A(7) = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
      A(8) = ZMS*A(7)*SQRT(A(7))
      A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
*
      DO 20 K = 1,3
          F1(K) = A(K)*A(8)
          F1DOT(K) = (A(K+3) - A(K)*A(9))*A(8)
   20 CONTINUE
*
*       Subtract dominant contributions from F/2 and FDOT/6.
      DO 25 K = 1,3
          F(K,I) = F(K,I) - SS*0.5D0*F1(K)
          FDOT(K,I) = FDOT(K,I) - SS*ONE6*F1DOT(K)
          FI(K,I) = FI(K,I) - SS*F1(K)
          FIDOT(K,I) = FIDOT(K,I) - SS*F1DOT(K)
   25 CONTINUE
*
      DE = SS*BODY(I)*ZMS*A(7)
      ECOLL = ECOLL + DE
      BE(3) = BE(3) + DE
*
      RETURN
*
      END
