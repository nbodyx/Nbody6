      SUBROUTINE R2SORT(IJ,RIJ2)
*
*
*       Sorting of square chain distances.
*       ---------------------------------
*
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
      COMMON/ARCHAIN/X(NMX3),V(NMX3),WTTL,M(NMX),
     &   XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),N
      REAL*8  M,MC,MASS,RIJ2(NMX,NMX)
      INTEGER IJ(NMX)
*
*
      RIJMIN=1.E30
      FMAX = 0.0
      DO I=1,N-1
      DO J=I+1,N
      RIJ2(I,J) = (X(3*I-2)-X(3*J-2))**2+(X(3*I-1)-X(3*J-1))**2+
     &            (X(3*I)-X(3*J))**2
      RIJ2(J,I)=RIJ2(I,J)
*     IF(RIJ2(I,J).LT.RIJMIN)THEN
*     RIJMIN=RIJ2(I,J)
*       Save dominant pair instead of smallest distance (08/09).
      IF ((M(I)+M(J))/RIJ2(I,J).GT.FMAX)THEN
      FMAX=(M(I)+M(J))/RIJ2(I,J)
      IJ(1)=I
      IJ(2)=J
      END IF
      END DO
      END DO
*
      L=2
      DO I=1,N
      IF(I.NE.IJ(1).AND.I.NE.IJ(2))THEN
      L=L+1
      IJ(L)=I
      END IF
      END DO
*
      RETURN
      END
