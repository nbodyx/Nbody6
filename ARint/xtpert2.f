       SUBROUTINE XTPERT2(JPAIR,ACC)
*
*       Chain perturbation for KS components.
*       -------------------------------------
*
        INCLUDE 'common6.h'
        REAL*8  M,MASS,MC
        PARAMETER (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &  NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
         COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &   XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
        REAL*8 ACC(NMX3),dx(3)
*       DATA IT /0/
        SAVE
*
*       Resolve KS pair (note possibility low/high order using GAMMA).
*       IF (GAMMA(JPAIR).GT.1.0D-03) THEN
            CALL KSRES(JPAIR,J1,J2,0.0D0)
*       ELSE
*           CALL KSRES(JPAIR,J1,J2,100.0D0)
*       END IF
*
        ic0=0
*       Loop over the chain members and components.
        do ic=1,NN
        I = J1
    1   rr=0.0
        do k=1,3
        dx(k)=XC(k,ic)-X(k,I)
        rr=rr+dx(k)**2
        end do
        rrr=sqrt(rr)*rr
        do k=1,3
        ACC(ic0+k)=ACC(ic0+k)-BODY(I)/rrr*dx(k)
        end do
        IF (I.EQ.J1) THEN
            I = I + 1
            GO TO 1
        END IF
        ic0=ic0+3
        end do
*
*     IT = IT + 1
*     IF (MOD(IT,50).EQ.0) THEN
*     WRITE (77,2)  NAME(N+JPAIR), TIME, R(JPAIR),H(JPAIR),GAMMA(JPAIR)
*   2 FORMAT (' XTPERT2   NM T R H G  ',I8,F10.4,1P,3E10.2)
*     CALL FLUSH(77)
*     END IF
        return
        end
