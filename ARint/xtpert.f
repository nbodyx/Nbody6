*     subroutine EXTERNAL ACCELERATIONS(ACC)
      SUBROUTINE XTPERT(ACC)
*
*       Chain perturbation.
*       -------------------
*
        INCLUDE 'common6.h'
        REAL*8  M,MASS,MC,MMIJ
        PARAMETER (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &  NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
         COMMON/ARCHAIN/XCH(NMX3),VCH(NMX3),WTTL,M(NMX),
     &   XCDUM(NMX3),WCDUM(NMX3),MC(NMX),
     &   XI(NMX3),VI(NMX3),MASS,RINV(NMXm),RSUM,INAME(NMX),NN
      COMMON/ARCHAIN2/ MMIJ,CMX(3),CMV(3),ENERGY,EnerGR,CHTIME
      common/TIMECOMMON/Taika,timecomparison
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
        REAL*8 ACC(NMX3),dx(3)
        SAVE
*
c       HERE ONE MUST EVALUATE THE ACCELERATIONS DUE TO THE PERTURBERS.
C       Physical positions and velocities (in the inertial coordinate)
C       system are in vectors X and V
C       (X(1)=X_1,X(2)=Y_1,X(3)=Z_1, X(4)=X_2, X(5)=Y_2,...)
C       After a call to this routine the EXAccerations
C       are assumed to be in the vector ACC.
*
        TIME0 = TIME
        TIME = Taika + CHTIME ! this is the time to be used for prediction.
c                         ! chtime is measured from beginning of CHAIN call.
        TIME = TIME + T0S(1)
*       Predict perturbers, c.m. and resolve chain members.
        CALL XCPRED(1)
        DO I=1,3*NN
        ACC(I)=0
        END DO
*       if(1.eq.1) return ! NO PERTURBATINS
*
        NPC = LISTC(1) + 1
*       IT = 0
*       Check whether to resolve any KS pair (c.m. perturbers are at end).
    1   IF (LISTC(NPC).GT.N) THEN
            JPAIR = LISTC(NPC) - N
            CALL XTPERT2(JPAIR,ACC)
            NPC = NPC - 1
*           IT = IT + 1
*           IF (IT.GT.1) THEN
*           WRITE (78,2)  TIME, NAME(N+JPAIR), R(JPAIR), GAMMA(JPAIR)
*   2       FORMAT (' CM PERT   T NAM R G  ',F10.4,I7,1P,3E10.2)
*           CALL FLUSH(78)
*           END IF
            IF (NPC.GT.1) GO TO 1
        END IF
*
        ic0=0
        do ic=1,NN ! ic= chainparticle index
        do L=2,NPC ! L=perturber list index
        I = LISTC(L)
        rr=0
        do k=1,3 ! coordinate index
        dx(k)=XC(k,ic)-X(k,I)
        rr=rr+dx(k)**2
        end do ! k
        rrr=sqrt(rr)*rr
        do k=1,3
        ACC(ic0+k)=ACC(ic0+k)-BODY(I)/rrr*dx(k)
        end do ! k
        end do ! i
        ic0=ic0+3
        end do ! ic         

        TIME = TIME0
        return
        end
