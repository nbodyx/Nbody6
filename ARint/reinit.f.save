      SUBROUTINE REINIT(ISUB)
*
*
*       Re-initialization of chain system.
*       ----------------------------------
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
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(NMX),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/ECHAIN/  ECH
      REAL*8  FIRR(3),FD(3)
*
*
*       Set phase indicator for step reduction (routine STEPS).
      IPHASE = 10
*
*       Re-initialize neighbour list of c.m. and form perturber list.
      RS0 = RS(ICH)
      CALL NBLIST(ICH,RS0)
      CALL CHLIST(ICH)
      NNB1 = LIST(1,ICH) + 1
*
*       Predict neighbours to FDOT and c.m. to FDOT3.
      DO 1 L = 2,NNB1
          J = LIST(L,ICH)
          CALL XVPRED(J,0)
    1 CONTINUE
      CALL XVPRED(ICH,-1)
*
*       Initialize force polynomials (include differential force correction).
      CALL FPOLY1(ICH,ICH,0)
      DO 5 K = 1,3
          FIRR(K) = 0.0
          FD(K) = 0.0
    5 CONTINUE
      CALL CHFIRR(ICH,2,X(1,ICH),XDOT(1,ICH),FIRR,FD)
      DO 10 K = 1,3
          F(K,ICH) = F(K,ICH) + 0.5*FIRR(K)       !factorial bugs 4/7/07
          FDOT(K,ICH) = FDOT(K,ICH) + ONE6*FD(K)
   10 CONTINUE
      CALL FPOLY2(ICH,ICH,0)
*
*       Obtain maximum unperturbed separation based on dominant neighbour.
      CALL EXTEND(ISUB)
*
*       Update chain variables (just in case).
      CALL XCPRED(0)
*
*       Update decision-making variables for chain regularization.
      TS(ISUB) = TIME
*       Note: probably done on next CALL CHAIN.
***   STEPS(ISUB) = 0.01*STEP(ICH)
*
*       Re-calculate new energy of chain system (not needed!).
*     CALL CONST(XCH,VCH,M,NN,ECH1,ANG,ALAG)
*
*       Set phase indicator < 0 to ensure new NLIST in routine INTGRT.
      IPHASE = -1
*
      RETURN
*
      END
