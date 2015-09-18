c       COMARC.CH =filename
        IMPLICIT REAL*8 (A-h,M,O-Z)
        PARAMETER (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &  NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
         COMMON/ARCHAIN/X(NMX3),V(NMX3),WTTL,M(NMX),
     &   XC(NMX3),WC(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/ARCHAIN2/MMIJ,CMX(3),CMV(3),ENERGY,CHTIME
      common/softening/ee,cmethod(3),Clight,NofBH
      common/TIMECOMMON/Taika
