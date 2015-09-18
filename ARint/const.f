      SUBROUTINE CONST(X,XD,M,NB,ENERGY,G,ALAG)
*
*
*       Constants of motion.
*       --------------------
*
      IMPLICIT REAL*8 (A-H,M,O-Z)
*     COMMON/SOFT/  EPS2
      REAL*8 X(3*NB),XD(3*NB),G(3),M(NB)
*
      EPS2 = 0.0
*
      U=0.0
      T=0.0
      UG=0.0
      G(1)=0.
      G(2)=0.
      G(3)=0.
      DO 10 I=1,NB
      K1=(I-1)*3+1
      K2=K1+1
      K3=K2+1
      T=T+.5D0*M(I)*(XD(K1)**2+XD(K2)**2+XD(K3)**2)
      G(1)=G(1)+M(I)*(X(K2)*XD(K3)-X(K3)*XD(K2))
      G(2)=G(2)-M(I)*(X(K1)*XD(K3)-X(K3)*XD(K1))
      G(3)=G(3)+M(I)*(X(K1)*XD(K2)-X(K2)*XD(K1))
      J1=I+1
      DO 9 J=J1,NB
      KI=(I-1)*3
      KJ=(J-1)*3
      R2=0.
      DO 8 K=1,3
      KI=KI+1
      KJ=KJ+1
      R2=R2+(X(KI)-X(KJ))**2
8     CONTINUE
*       Include possible softening for singular interaction.
      IF(J.LT.NB)R2=R2+EPS2
      U=U+M(I)*M(J)/SQRT(R2)
9     CONTINUE
10    CONTINUE
*
*       Form the total energy and Lagrangean energy.
      ENERGY=T-U
      ALAG=T+U
*     G(1) = U
*
      RETURN
      END
