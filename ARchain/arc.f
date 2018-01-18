      SUBROUTINE ARC   !CALLED FROM CHAIN
     &    (NN,XX,VX,MX,TIME,DELTAT,TOL,NEWREG,KSMX,soft,cl,Ixc,NBH,
     &    spini,CMXX,CMVX)
c       BETTER TO USE CM-coords & vels for XX & VX and CMXX CMVX
c       FOR CM-position (needed in the Perturbations routine).
c-----------------------------------------------------------------
c       NOTE: some variables (eg. Energy and EnerGR) are only in the
c       COMMON. The internal NB-energy = ENERGY+EnerGR  (should be)
c       Energy =  integrated E-value (excluding grav.radiation)
c       EnerGr =  Energy radiated away (grav.radiation IF Clight.ne.0.0)
C       CHAIN INTEGRATION. Perturbations & CM-motion included (in principle).
c       NN = # of bodies; XX = (cm)coords, VX = (cm)vels, MX = masses,
c       CMXX = coords of CM, CMVX = vels of CM ! removed
c       TIME = time, deltaT = time interval
c       STEP = stepsize (set = 0 initially)
c       NEWREG = .true. iff chain membership has changed
c       KSMX = max # of steps without RETURN (use some large # )
c       soft  = optional softening( U = 1/sqrt(r**2+soft**2) )
c       cmethod = 3-d vector that determines the method:
c       (1,0,0)  = logH, (0,1,0) = TTL,(0,0,1) = DIFSY2 without t-tranformation
c       cl = speed of light
c       NOTE: cl = 0  = > no relativistic terms !!!
c       Ixc = 1  = > exact time,  = 0 no exact time but RETURN after CHTIME>DELTAT
*       Note: EnerGR = 0 in version of June/15.

      INCLUDE 'ARCCOM2e2.ch'
      COMMON/DerOfTime/         GTIME
      COMMON/DIAGNOSTICS/       GAMMA,H,IWR
      COMMON/omegacoefficients/ OMEC(NMX,NMX)
      COMMON/collision/         icollision,ione,itwo,iwarning
!       COMMON/itemaxcommon/      aitemax,itemax,itemax_used

      REAL*8 G0(3),XX(*),VX(*),MX(*),spini(3),CMXX(3),CMVX(3)
      REAL*8 Y(1500),SY(1500),Yold(1500)
      LOGICAL MUSTSWITCH,NEWREG
      SAVE

      CHTIME = 0.0
      icollision = 0
      Taika = TIME ! to common
      NofBH = NBH  ! - " -

      IF (NEWREG) THEN
          step = 0
          iwarning = 0
!           itemax = 12
!           itemax_used = 0
          ee = soft**2  ! to common
          DO k = 1,3
              spin(k) = spini(k) ! SPIN
          END DO
          clight = cl   ! -"-
          N = NN
          mass = 0
          DO I = 1,N
              M(I) = MX(I)
              mass = mass+m(i)
          END DO

          MMIJ = 0.0
          DO I = 1,N-1
              DO J = I+1,N
                  MMIJ = MMIJ+M(I)*M(J)
              END DO
          END DO
          MMIJ = 2.0*MMIJ/(N*(N-1))
          DO I = 1,3*N
              X(I) = XX(I)
              V(I) = VX(I)
          END DO
          CALL FindChainIndices
*          IF (IWR.GT.0) WRITE (6,1232)time,(INAME(KW),KW = 1,N)
          CALL InitializeXcAndWc
          CALL ConstantsOfMotion(ENERGY,G0,ALAG)
          EnerGr = 0 ! energy radiated away
          gtime = 1/ALAG
          DO K = 1,3
              CMX(K) = CMXX(K)
              CMV(K) = CMVX(K)
          END DO
          CALL omegacoef
          STIME = 0.0
          NEWREG = .FALSE.
          WTTL = Wfunction()
          mmss = 0
          DO i = 1,n-1
              DO j = i+1,n
                  mmss = mmss+m(i)*m(j)
              END DO
          END DO
          CALL TakeYfromXcWc(Y,Nvar)
          DO i = 1,Nvar
              SY(i) = 0
          END DO
          IF (step.EQ.0.0) CALL InitialStepsize(X,V,M,N,ee,step) ! New step
          stimex = step
          EPS = TOL
          NCALL = 0
      END IF ! End NEWREG

      KSTEPS = 0
      nzero = 0
      stw = stimex
      step = min(step,2*ABS(stimex))
      stimex = 0
  777 KSTEPS = KSTEPS+1
      CALL TakeYfromXcWc(Y,Nvar)
      CALL ObtainOrderOfY(SY)
      stime = 0
      f1 = chtime-deltaT ! for exact time
      d1 = gtime
      dltime = -f1
      CALL TakeYfromXcWc(Yold,Nvar)
      CALL DIFSYAB(Nvar,EPS,SY,step,stime,Y)
      I_switch = 1
      CALL PutYtoXCWC(Y,Nvar)
      CALL CheckSwitchingConditions(MUST SWITCH)
      IF (MUST SWITCH) THEN
          I_switch = 0
          CALL ChainTransformation !
          WTTL = Wfunction() ! this may not be necessary, but probably OK.
          CALL TakeYfromXcWc(Y,Nvar)
*          IF (IWR.GT.0) WRITE (6,1232)time+chtime,(INAME(KW),KW = 1,N)
* 1232     FORMAT (1X,g12.4,' I-CHAIN',20I3)
      END IF ! MUST SWITCH
      f2 = chtime-deltaT ! for exact time iteration
      d2 = gtime
      x1 = -stime
      x2 = 0.0
      DLT = DELTAT! for short
      IF (CHTIME.LT.DLT.AND.(KSTEPS.LT.KSMX)
     &    .AND.(icollision.EQ.0)) GOTO 777
      IF (KSTEPS.LT.KSMX.AND.Ixc.EQ.1.AND.icollision.EQ.0) THEN
        ! Integrate to approximate EXACT OUTPUT TIME
          IF (ABS(f1).LT.ABS(f2)*I_switch) THEN ! I_switch prevents use of
              CALL PutYtoXCWC(yold,nvar)   ! f1 IF just SWITCHed
              CALL ObtainOrderOfY(sy)
              CALL EstimateStepsize(-f1,step2)
              cht_0 = chtime
              CALL DIFSYAB(Nvar,EPS,SY,step2,stime,Yold)
              CALL PutYtoXCWC(Yold,Nvar)
c              WRITE (6,*)'T: f1,f2,dlt_cht ',f1,f2,deltat-chtime,deltat-cht0
          ELSE
              CALL EstimateStepsize(-f2,step2)
              CALL ObtainOrderOfY(sy)
              cht_0 = chtime
              CALL DIFSYAB(Nvar,EPS,SY,step2,stime,Y)
              CALL PutYtoXCWC(Y,Nvar)
c              WRITE (6,*)'T: f2,f1,dlt_cht ',f2,f1,deltat-chtime,deltat-cht0
          END IF
          stimex = stimex+stime! for estimating maximum next step
c          CALL Iterate2ExactTime(Y,Nvar,deltaT,f1,d1,f2,d2,x1,x2,eps)
      END IF
      IF (stimex.EQ.0.0) stimex = step
      CALL UpdateXandV
      DO I = 1,3*N
          XX(I) = X(I)
          VX(I) = V(I)
      END DO
      DO I = 1,3
          spini(I) = spin(I)
          CMXX(I) = CMX(I)
          CMVX(I) = CMV(I)
      END DO
      TIME = TIME+CHTIME

      RETURN
      END

      SUBROUTINE Iterate2ExactTime
     &           (Y0,Nvar,deltaT,f1,d1,f2,d2,x1,x2,eps)
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/DerOfTime/ GTIME
      COMMON/collision/ icollision,Ione,Itwo,iwarning
      REAL*8 Y(1500),SY(1500),Y0(*)
      DATA tiny/1.d-6/
      SAVE

      iskeleita = 0
      it = 0
      hs = ABS(x1-x2)
 1111 CONTINUE
      it = it+1
      DO i = 1,nvar
          y(i) = y0(i)
      END DO
      stime = 0
      dx1 = -f1/d1
      dx2 = -f2/d2
      IF (ABS(dx1).LT.ABS(dx2))THEN
          xnew = x1+dx1
      ELSE
      xnew = x2+dx2
      END IF
c
      test = (x1-xnew)*(xnew-x2)
      IF (test.LT.(-tiny*hs).OR.(it+1).EQ.(it+1)/5*5) THEN
          xnew = (x1+x2)/2 ! bisect IF out of interval
      END IF

      sfinal = xnew

      CALL PutYtoXCWC(Y,Nvar)
c----------------------------------------------------------------------
      CALL ObtainOrderOfY(SY)
      steppi = 0
      DO k = 1,5
          step = sfinal-stime
          IF (ABS(step).GT.1.e-3*ABS(hs).OR.k.EQ.1) THEN !!!!
              steppi = step
              CALL DIFSYAB(Nvar,EPS,SY,step,stime,Y)
              iskeleita = iskeleita+1
c              it = it+1
          ELSE
              GOTO 222
          END IF
      END DO
  222 CONTINUE
      CALL PutYtoXCWC(Y,Nvar)
      CALL UpdateXandV
      fnew = chtime-deltaT
      dfnew = gtime
c       keep it bracketed
      IF (f1*fnew.LE.0.0) THEN
          f2 = fnew
          d2 = dfnew
          x2 = xnew
      ELSE
          f1 = fnew
          d1 = dfnew
          x1 = xnew
      END IF
      IF ((ABS(deltaT-chtime).GT.1.e-3*deltat).AND.(it.LT.100))
     &    GOTO 1111
c       ONE FINAL STEP SHOULD BE HERE (IF above not-so-accurate test)
c--------------------------------------------------------------------
      DO i = 1,Nvar
          y0(i) = y(i)
      END DO
      CALL PutYtoXCWC(Y,Nvar)
      CALL UpdateXandV

      RETURN
      END

      SUBROUTINE LEAPFROG(STEP,Leaps,stime)
      IMPLICIT REAL*8 (a-h,M,o-z)
      SAVE

      CALL PUTV2W
      hs = step
      h2 = hs/2
      CALL XCmotion(h2)
      stime = stime+h2
      DO k = 1,Leaps-1
          CALL WCmotion(hs)
          CALL XCmotion(hs)
          stime = stime+hs
      END DO
      CALL WCmotion(hs)
      CALL XCmotion(h2)
      stime = stime+h2

      RETURN
      END

      function Wfunction()
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/omegacoefficients/ OMEC(NMX,NMX)
      SAVE
      OMEGA = 0.0d0
      DO I = 1,N-1
          DO J = I+1,N
              IF (omec(i,j).ne.0.0) THEN
                  RIJ = SQRT(SQUARE(X(3*I-2),X(3*J-2)))
                  OMEGA = OMEGA+omec(i,j)/RIJ
              END IF
          END DO
      END DO
      Wfunction = OMEGA

      RETURN
      END

      SUBROUTINE omegacoef
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/omegacoefficients/ OMEC(NMX,NMX)
      SAVE

      icount = 0
      DO i = 1,N-1
          DO j = i+1,N
              IF (m(i)+m(j).GT.0.0 .AND. cmethod(2).ne.0.0) THEN
                  OMEC(I,J) = mmij
                  OMEC(J,I) = mmij
                  icount = icount+1
              ELSE
                  OMEC(I,J) = 0
                  OMEC(J,I) = 0
              END IF
          END DO
      END DO
      IF (icount.EQ.0.0) cmethod(2) = 0 ! all terms zero anyway

      RETURN
      END

      SUBROUTINE XCMOTION(hs)
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/IncrementCommon/ WTTLinc,XCinc(NMX3),WCinc(NMX3),
     &                        CMXinc(3),CMVinc(3),ENERGYinc,
     &                        Energrinc,CHTIMEinc,spin inc(3)
      COMMON/DerOfTime/       G
      COMMON/DIAGNOSTICS/     GAMMA,H,IWR
      SAVE

**      Te = -ENERGY-EnerGR
      Te = -ENERGY
      IF (cmethod(1).ne.0.0d0) THEN
          CALL EVALUATEV(V,WC)
          DO I = 1,N
              I0 = 3*I-3
              Te = Te+M(I)*(V(I0+1)**2+V(I0+2)**2+V(I0+3)**2)/2
          END DO
      END IF ! cmethod(1).ne.0.0d0
      G = 1/(Te*cmethod(1)+WTTL*cmethod(2)+cmethod(3)) ! = t'
      IF (G.LT.0.0.AND.iwr.GT.0) THEN
          WRITE (6,*)1/G,' tdot <0 ! '
          RETURN ! seriously wrong, but may work (this step gets rejected)
      END IF
      dT =  hs*G
      DO I = 1,N-1
          L = 3*(I-1)
          DO K = 1,3
              XCinc(L+K) = XCinc(L+k)+WC(L+K)*dT
              XC(L+K) = XC(L+K)+WC(L+K)*dT
          END DO
      END DO
      CHTIMEinc = CHTIMEinc+dT
      CHTIME = CHTIME+dT
      DO k = 1,3
          CMXinc(k) = CMXinc(k)+dt*cmv(k)
          cmx(k) = cmx(k)+dt*cmv(k)
      END DO

      RETURN
      END

      SUBROUTINE PUTV2W
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/vwcommon/ Ww(nmx3),WTTLw,cmvw(3),spinw(3)
      SAVE

      DO i = 1,3*(N-1)
          Ww(i) = WC(I)
      END DO
      WTTLw = WTTL
      DO k = 1,3
          spinw(k) = spin(k)
          cmvw(k) = cmv(k)
      END DO

      RETURN
      END

      SUBROUTINE VelocityDependentPerturbations
     &           (dT,Va,spina,acc,dcmv,df,dfGR,dspin)
      INCLUDE 'ARCCOM2e2.ch'
      REAL*8 df(*),Va(*),dcmv(3),dfGR(*),dfR(nmx3),acc(nmx3)
      REAL*8 dspin(3),spina(3)
      SAVE

      IF (Clight.ne.0.0) THEN ! INCLUDE only IF Clight set >0
          CALL RelativisticAccelerations(dfr,dfGR,Va,spina,dspin)
      ELSE
          DO i = 1,3*n
              dfr(i) = 0
              dfgr(i) = 0
          END DO
          DO k = 1,3
              dspin(k) = 0
          END DO
      END IF
      DO i = 1,3*n
          df(i) = acc(i)+dfr(i)
      END DO
      CALL Reduce2cm(df,m,n,dcmv)

      RETURN
      END

      SUBROUTINE CheckSwitchingConditions(MUSTSWITCH)
      INCLUDE 'ARCCOM2e2.ch'
      LOGICAL MUSTSWITCH
      DATA NCALL,NSWITCH/0,200000000/
      SAVE

      MUSTSWITCH = .FALSE.
      NCALL = NCALL+1
C       Switch anyway after every NSWITCHth step.
      IF (NCALL.GE.NSWITCH) THEN
          NCALL = 0
          MUSTSWITCH = .TRUE.
          RETURN
      END IF
C       Inspect the structure of the chain.
C       NOTE: Inverse values 1/r are used instead of r itself.
      ADISTI = 0.5*(N-1)/RSUM
      LRI = N-1
      DO I = 1,N-2
          DO J = I+2,N
              LRI = LRI+1
C       do not inspect IF 1/r is small.
              IF (RINV(LRI).GT.ADISTI) THEN
                  IF (J-I.GT.2) THEN
C       Check for a dangerous long loop.
C                      RINVMX = MAX(RINV(I-1),RINV(I),RINV(J-1),RINV(J))
                      IF (I.GT.1) THEN
                          RINVMX = MAX(RINV(I-1),RINV(I))
                      ELSE
                          RINVMX = RINV(1)
                      END IF
                      RINVMX = MAX(RINVMX,RINV(J-1))
                      IF (J.LT.N) RINVMX = MAX(RINVMX,RINV(J))
                      IF (RINV(LRI).GT.RINVMX) THEN ! 0.7*RINVMX may be more careful
                          MUSTSWITCH = .TRUE.
                          NCALL = 0
                          RETURN
                      END IF
                  ELSE
C       Is this a triangle with smallest size not regularised?
                      IF (RINV(LRI).GT.MAX(RINV(I),RINV(I+1))) THEN
                          MUSTSWITCH = .TRUE.
                          NCALL = 0
                          RETURN
                      END IF
                  END IF
              END IF
          END DO
      END DO

      RETURN
      END

      SUBROUTINE FindChainIndices
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/EXTRA/ KNAME(NMX)
      REAL*8 RIJ2(NMXM)
      INTEGER IC(NMX2),IJ(NMXM,2),IND(NMXM)
      LOGICAL USED(NMXM),SUC,LOOP
      SAVE

      L = 0
      DO I = 1,N-1
          DO J = I+1,N
              L = L+1
              RIJ2(L) = SQUARE(X(3*I-2),X(3*J-2))
              IJ(L,1) = I
              IJ(L,2) = J
              USED(L) = .FALSE.
          END DO
      END DO
      CALL ARRANGE(L,RIJ2,IND)
      LMIN = 1+NMX
      LMAX = 2+NMX
      IC(LMIN) = IJ(IND(1),1)
      IC(LMAX) = IJ(IND(1),2)
      USED(IND(1)) = .TRUE.
1     DO I = 2,L
          LI = IND(I)
          IF (.NOT.USED(LI)) THEN
              CALL CheckConnection(IC,LMIN,LMAX,IJ,LI,SUC,LOOP)
              IF (SUC) THEN
                  USED(LI) = .TRUE.
                  GOTO 2
              ELSE
                  USED(LI) = LOOP
              END IF
          END IF
      END DO
2     IF (LMAX-LMIN+1.LT.N) GO TO 1
      L = 0
      DO I = LMIN,LMAX
          L = L+1
          INAME(L) = IC(I)
      END DO
*       Form the inverse table (Seppo 13/6/10).
      DO I = 1,N
          KNAME(INAME(I)) = I
      END DO

      RETURN
      END

      SUBROUTINE CheckConnection(IC,LMIN,LMAX,IJ,LI,SUC,LOOP)
      INCLUDE 'ARCCOM2e2.ch'
      INTEGER IC(*),ICC(2),IJ(NMXM,2)
      LOGICAL SUC,LOOP
      SAVE

      SUC = .FALSE.
      LOOP = .FALSE.
      ICC(1) = IC(LMIN)
      ICC(2) = IC(LMAX)
      DO I = 1,2
          DO J = 1,2
              IF (ICC(I).EQ.IJ(LI,J)) THEN
                  JC = 3-J
                  LOOP = .TRUE.
                  DO L = LMIN,LMAX
                      IF (IC(L).EQ.IJ(LI,JC)) RETURN
                  END DO
                  SUC = .TRUE.
                  LOOP = .FALSE.
                  IF (I.EQ.1) THEN
                      LMIN = LMIN-1
                      IC(LMIN) = IJ(LI,JC)
                      RETURN
                  ELSE
                      LMAX = LMAX+1
                      IC(LMAX) = IJ(LI,JC)
                      RETURN
                  END IF
              END IF
          END DO
      END DO

      RETURN
      END

      SUBROUTINE ARRANGE(N,Array,Indx)
      IMPLICIT REAL*8 (a-h,o-z)
      DIMENSION Array(*),Indx(*)
      SAVE

      DO j = 1,N
          Indx(j) = j
      END DO
      IF (N.LT.2) RETURN
      l = N/2+1
      ir = N
10    IF (l.GT.1) THEN
          l = l-1
          Indxt = Indx(l)
          q = Array(Indxt)
      ELSE
          Indxt = Indx(ir)
          q = Array(Indxt)
          Indx(ir) = Indx(1)
          ir = ir-1
          IF (ir.EQ.1) THEN
              Indx(1) = Indxt
              RETURN
          END IF
      END IF
      i = l
      j = l+l
20    IF (j.LE.ir) THEN
          IF (j.LT.ir) THEN
              IF (Array(Indx(j)).LT.Array(Indx(j+1))) j = j+1
          END IF
          IF (q.LT.Array(Indx(j))) THEN
              Indx(i) = Indx(j)
              i = j
              j = j+j
          ELSE
              j = ir+1
          END IF
          GOTO 20
      END IF
      Indx(i) = Indxt
      GO TO 10

      END

      SUBROUTINE InitializeXcAndWc  ! Find CHAIN coords & velocities.
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,JCOLL,NDISS1
      COMMON/ARZERO/  ISTAR0(NMX),SIZE0(NMX)    ! New procedure 10/16.
      SAVE

C       Initialize centre of mass.
      DO K = 1,3
          CMX(K) = 0.0
          CMV(K) = 0.0
      END DO
      MASS = 0.0
      DO I = 1,N
          L = 3*(I-1)
          MC(I) = M(INAME(I)) ! masses along the chain
          MASS = MASS+MC(I)
          DO K = 1,3
              CMX(K) = CMX(K)+M(I)*X(L+K)
              CMV(K) = CMV(K)+M(I)*V(L+K)
          END DO
      END DO
      DO K = 1,3
          CMX(K) = CMX(K)/MASS
          CMV(K) = CMV(K)/MASS
      END DO
c       Rearange according to chain indices.
      DO I = 1,N
          L = 3*(I-1)
          LF = 3*INAME(I)-3
          DO K = 1,3
              XI(L+K) = X(LF+K)
              VI(L+K) = V(LF+K)
          END DO
      END DO

!       Chain coordinates & vels and initialize WTTL.!
      WTTL = 0.0
      DO I = 1,N-1
      L = 3*(I-1)
      DO K = 1,3
      XC(L+K) = XI(L+K+3)-XI(L+K)
      WC(L+K) = VI(L+K+3)-VI(L+K)
      END DO
      L1 = L+1
      r2 = cdot(XC(L1),XC(L1))
      RI = sqrt(r2)
      RINV(I) = 1/RI
      END DO

*       Initialize stellar evolution indices and radii.
      DO I = 1,N
          ISTAR(I) = ISTAR0(INAME(I)) ! This is correct if ISTAR0 and SIZE0
          SIZE(I) = SIZE0(INAME(I))   ! are in common and related to the 
                                      ! true star indices, not chain indices.
      END DO

      RETURN
      END

      SUBROUTINE UpdateXandV
      INCLUDE 'ARCCOM2e2.ch'
      REAL*8 X0(3),V0(3)
      SAVE

C       Obtain physical variables from chain quantities.
      DO K = 1,3
          XI(K) = 0.0
          VI(k) = 0.0
          X0(K) = 0.0
          V0(k) = 0.0
      END DO
      DO I = 1,N-1
          L = 3*(I-1)
          DO K = 1,3
              VI(L+3+K) = VI(L+K)+WC(L+K)
              XI(L+3+K) = XI(L+K)+XC(L+K)
          END DO
      END DO
      DO I = 1,N
          L = 3*(I-1)
          DO K = 1,3
              V0(K) = V0(K)+VI(L+K)*MC(I)/MASS
              X0(K) = X0(K)+XI(L+K)*MC(I)/MASS
          END DO
      END DO
C       Rearrange according to INAME(i) and add CM.
      DO I = 1,N
          L = 3*(I-1)
          LF = 3*(INAME(I)-1)
          DO K = 1,3
              X(LF+K) = XI(L+K)-X0(K)!+CMX(K) ! CM-coords
              V(LF+K) = VI(L+K)-V0(K)!+CMV(K) ! CM-vels
          END DO
      END DO

      RETURN
      END

      SUBROUTINE ChainTransformation

      INCLUDE 'ARCCOM2e2.ch'

      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,JCOLL,NDISS1
      COMMON/ARZERO/  ISTAR0(NMX),SIZE0(NMX)    ! New procedure 10/16.
      REAL*8 XCNEW(NMX3),WCNEW(NMX3)
      INTEGER IOLD(NMX)
      SAVE

      L2 = 3*(INAME(1)-1)
      DO K = 1,3
          X(L2+K) = 0.0
      END DO
C       Xs are needed when determining new chain indices.
      DO I = 1,N-1
          L = 3*(I-1)
          L1 = L2
          L2 = 3*(INAME(I+1)-1)
          DO K = 1,3
              X(L2+K) = X(L1+K)+XC(L+K)
          END DO
      END DO
C       Save the old chain indices, stellar types and radii.
      DO I = 1,N
          IOLD(I) = INAME(I)
      END DO

C       Find new ones.
      CALL FindChainIndices

C       Construct new chain coordinates. Transformation matrix
C       (from old to new) has only coefficients -1, 0 or +1.
      DO I = 1,3*(N-1)
          XCNEW(I) = 0.0
          WCNEW(I) = 0.0
      END DO
      DO ICNEW = 1,N-1
C       Obtain K0 &  K1 such that iold(k0) = iname(icnew)
c                                 iold(k1) = iname(icnew+1)
          LNEW = 3*(ICNEW-1)
          DO I = 1,N
              IF (IOLD(I).EQ.INAME(ICNEW)) K0 = I
              IF (IOLD(I).EQ.INAME(ICNEW+1)) K1 = I
          END DO
          DO ICOLD = 1,N-1
              LOLD = 3*(ICOLD-1)
              IF ((K1.GT.ICOLD).AND.(K0.LE.ICOLD)) THEN
C       Add
                  DO K = 1,3
                      XCNEW(LNEW+K) = XCNEW(LNEW+K)+XC(LOLD+K)
                      WCNEW(LNEW+K) = WCNEW(LNEW+K)+WC(LOLD+K)
                  END DO
              ELSEIF ((K1.LE.ICOLD).AND.(K0.GT.ICOLD)) THEN
C        Subtract
                  DO K = 1,3
                      XCNEW(LNEW+K) = XCNEW(LNEW+K)-XC(LOLD+K)
                      WCNEW(LNEW+K) = WCNEW(LNEW+K)-WC(LOLD+K)
                  END DO
              END IF
          END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I = 1,3*(N-1)   !!!!!!!!!!!!!!!!!!
          xc(i) = xcnew(i)   !!!!!!!!!!!!!!!!!!!
          wc(i) = wcnew(i)
      END DO           !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       Auxiliary quantities, including stellar types and radii.
      MASS = 0.0
      DO I = 1,N
          MC(I) = M(INAME(I))
          MASS = MASS+MC(I)
          ISTAR(I) = ISTAR0(INAME(I)) ! This is correct if ISTAR0 and SIZE0
          SIZE(I) = SIZE0(INAME(I))   ! are in common and related to the 
                                      ! true star indices, not chain indices.
      END DO

      RETURN
      END

      FUNCTION SQUARE(X,Y)
      IMPLICIT REAL*8 (a-h,m,o-z)
      COMMON/softening/ ee,cmethod(3),clight,NofBH ! only ee needed here
      REAL*8 X(3),Y(3),SQUARE
      SAVE

      SQUARE = (X(1)-Y(1))**2+(X(2)-Y(2))**2+(X(3)-Y(3))**2+ee

      RETURN
      END

      SUBROUTINE DIFSYAB(N,EPS,S,h,t,Y)!,Jmax)
      IMPLICIT REAL*8 (a-h,o-z)
c       N = coordin. maara ( = 3*NB)
c       F = funktion nimi (FORCE)
      PARAMETER (NMX = 1500,NMX2 = 2*NMX,nmx27 = nmx2*7) ! NMX = MAX(N),N = 3*NB
      REAL*8 Y(N),YR(NMX2),YS(NMX2),y0(NMX),DT(NMX2,7),D(7),S(N),EP(4)
      LOGICAL KONV,BO,KL,GR
      DATA EP/.4D-1,.16D-2,.64D-4,.256D-5/
      DATA dt/nmx27*0.0d0/
      SAVE

      Jmax = 10 ! JMAX set here
      IF (EPS.LT.1.D-14) EPS = 1.D-14
      IF (N.GT.NMX) WRITE (6,*) ' too many variables!', char(7)
      IF (jmax.LT.4) WRITE (6,*)' too small Jmax ( = ',jmax,')'
      JTI = 0
      FY = 1
      redu = 0.8d0
      Odot7 = 0.7
      DO i = 1,N
          y0(i) = y(i)
          s(i) = MAX(ABS(y0(i)),s(i))
      END DO
10    tN = t+H
      BO = .FALSE.
C
      M = 1
      JR = 2
      JS = 3
      DO J = 1,Jmax  ! 10
          DO i = 1,N
              ys(i) = y(i)
              s(i) = MAX(ABS(ys(i)),s(i))
          END DO
C
          IF (BO) THEN
              D(2) = 1.777777777777778D0
              D(4) = 7.111111111111111D0
              D(6) = 2.844444444444444D1
          ELSE
              D(2) = 2.25D0
              D(4) = 9.D0
              D(6) = 36.0D0
          END IF

          IF (J.GT.7) THEN
              L = 7
              D(7) = 6.4D1
          ELSE
              L = J
              D(L) = M*M
          END IF

          KONV = L.GT.3
          subH = H/M
          CALL SubSteps(Y0,YS,subH,M) ! M substeps of size H/M.
          KL = L.LT.2
          GR = L.GT.5
          FS = 0.

          DO I = 1,N
              V = DT(I,1)
              C = YS(I)
              DT(I,1) = C
              TA = C
              IF (.NOT.KL) THEN
                  DO K = 2,L
                      B1 = D(K)*V
                      B = B1-C
                      W = C-V
                      U = V
                      IF (B.ne.0.0) THEN
                          B = W/B
                          U = C*B
                          C = B1*B
                      END IF
                      V = DT(I,K)
                      DT(I,K) = U
                      TA = U+TA
                  END DO ! K = 2,L
                  SI = MAX(S(I),ABS(TA))
                  IF (DABS(YR(I)-TA).GT.SI*EPS) KONV = .FALSE.
                  IF (.NOT.(GR.OR.SI.EQ.0.D0)) THEN
                      FV = DABS(W)/SI
                      IF (FS.LT.FV) FS = FV
                  END IF
              END IF ! .NOT.KL.
              YR(I) = TA
          END DO ! I = 1,N
c       end of I-loop
          IF (FS.NE.0.D0) THEN
              FA = FY
              K = L-1
              FY = (EP(K)/FS)**(1.d0/FLOAT(L+K))
              FY = min(FY,1.4) !1.4 ~ 1/0.7 ; where 0.7 = initial reduction factor
              IF (.NOT.((L.NE.2.AND.FY.LT.Odot7*FA)
     &            .OR.FY.GT.Odot7)) THEN
                  H = H*FY
                  JTI = JTI+1
                  IF (JTI.GT.25) THEN
                      H = 0.0
                      RETURN
                  END IF
                  GO TO 10 ! Try again with a smaller step.
              END IF
          END IF

          IF (KONV) THEN
              t = tN
              H = H*FY
              DO I = 1,N
                  Y(I) = YR(I)+y0(i) !!!!!!!
              END DO
              RETURN
          END IF

          D(3) = 4.D0
          D(5) = 1.6D1
          BO = .NOT.BO
          M = JR
          JR = JS
          JS = M+M
      END DO ! J = 1,Jmax
      redu = redu*redu+.001d0 ! square reduction factor (minimum near 0.001).
      H = H*redu
      GO TO 10 ! Try again with smaller step.

      END

      SUBROUTINE SubSteps(Y0,Y,H,Leaps)
      IMPLICIT REAL*8 (a-h,m,o-z)
      COMMON/softening/ ee,cmethod(3),Clight,NofBh
      COMMON/collision/ icollision,Ione,Itwo,iwarning
      REAL*8 Y(*),Y0(*)!,ytest(1000)
      SAVE

      icollision = 0
      CALL PutYtoXCWC(Y0,Nvar) ! Y -> XC, WTTL, WC
      CALL InitializeIncrements2Zero
      CALL LEAPFROG(H,Leaps,stime) ! advance
      CALL TakeIncrements2Y(y)

      RETURN
      END

      SUBROUTINE InitializeIncrements2Zero
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/IncrementCommon/ WTTLinc,XCinc(NMX3),WCinc(NMX3),
     &                        CMXinc(3),CMVinc(3),ENERGYinc,
     &                        Energrinc,CHTIMEinc,spin inc(3)

      DO i = 1,3*(N-1)
          XCinc(i) = 0
          WCinc(i) = 0
      END DO
      DO k = 1,3
          CMXinc(k) = 0
          CMVinc(k) = 0
          spin inc(k) = 0
      END DO
      WTTLinc = 0
      ENERGYinc = 0
      EnerGRinc = 0
      CHTIMEinc = 0

      RETURN
      END

      SUBROUTINE TakeIncrements2Y(Y)
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/IncrementCommon/ WTTLinc,XCinc(NMX3),WCinc(NMX3),
     &                        CMXinc(3),CMVinc(3),ENERGYinc,
     &                        Energrinc,CHTIMEinc,spin inc(3)
      REAL*8 Y(*)
      SAVE

      L = 1
      Y(L) = CHTIMEinc
      DO i = 1,3*(N-1)
          L = L+1
          Y(L) = XCinc(I)
      END DO
      L = L+1
      Y(L) = WTTLinc
      DO i = 1,3*(N-1)
          L = L+1
          Y(L) = WCinc(I)
      END DO
      DO i = 1,3
          L = L+1
          Y(L) = CMXinc(I)
      END DO
      DO i = 1,3
          L = L+1
          Y(L) = CMVinc(I)
      END DO
      L = L+1
      Y(L) = ENERGYinc
      L = L+1
      Y(L) = EnerGRinc
      DO k = 1,3
          L = L+1
          Y(L) = spin inc(k)
      END DO
c       Nvar = L

      RETURN
      END

      SUBROUTINE PutYtoXCWC(Y,Lmx)
      INCLUDE 'ARCCOM2e2.ch'
      REAL*8 Y(*)
      SAVE

      L = 1
      CHTIME = Y(L)
      DO i = 1,3*(N-1)
          L = L+1
          XC(I) = Y(L)
      END DO
      L = L+1
      WTTL = Y(L)
      DO i = 1,3*(N-1)
          L = L+1
          WC(I) = Y(L)
      END DO
      DO i = 1,3
          L = L+1
          CMX(I) = Y(L)
      END DO
      DO i = 1,3
          L = L+1
          CMV(I) = Y(L)
      END DO
      L = L+1
      ENERGY = Y(L)
      L = L+1
      EnerGR = Y(L)
      DO k = 1,3
          L = L+1
          spin(k) = Y(L)
      END DO
      Lmx = L

      RETURN
      END

      SUBROUTINE TakeYfromXcWc(Y,Nvar)
      INCLUDE 'ARCCOM2e2.ch'
      REAL*8 Y(*)
      SAVE

      L = 1
      Y(L) = CHTIME
      DO i = 1,3*(N-1)
          L = L+1
          Y(L) = XC(I)
      END DO
      L = L+1
      Y(L) = WTTL
      DO i = 1,3*(N-1)
          L = L+1
          Y(L) = WC(I)
      END DO
      DO i = 1,3
          L = L+1
          Y(L) = CMX(I)
      END DO
      DO i = 1,3
          L = L+1
          Y(L) = CMV(I)
      END DO
      L = L+1
      Y(L) = ENERGY
      L = L+1
      Y(L) = EnerGR
      DO k = 1,3
          L = L+1
          Y(L) = spin(k)
      END DO
      Nvar = L

      RETURN
      END

      SUBROUTINE ObtainOrderOfY(SY)
      INCLUDE 'ARCCOM2e2.ch'
      REAL*8 SY(*)
      SAVE

      w_old = 0.90
      w_new = 1-w_old
      L = 1
      SY(L) = ABS(CHTIME)*w_new+sy(L)*w_old
      SR = 0
      XCmin = 1.d99
      UPO = 0
      DO i = 1,N-1
          i0 = 3*i-3
          XCA = ABS(XC(I0+1))+ABS(XC(I0+2))+ABS(XC(I0+3))
          SR = SR+XCA
          UPO = UPO+MMIJ/XCA
          XCmin = min(XCA,XCmin)
          DO k = 1,3
              L = L+1
              SY(L) = XCA*w_new+sy(L)*w_old
          END DO ! k
      END DO  ! I
      L = L+1
      SY(L) = (ABS(WTTL*1.e2)+mass**2/XCmin)*w_new+sy(L)*w_old
      SW0 = SQRT(ABS(Energy/mass))
      SW = 0
      DO i = 1,N-1
          i0 = 3*i-3
          WCA = ABS(WC(I0+1))+ABS(WC(I0+2))+ABS(WC(I0+3))
          SW = SW+WCA
          DO k = 1,3
              L = L+1

              IF (WCA.ne.0.0) THEN
                  SY(L) = WCA*w_new+sy(L)*w_old
              ELSE
                  SY(L) = SW0*w_new+sy(L)*w_old
              END IF
          END DO ! k
      END DO ! i

      L = 1
      DO i = 1,N-1
          i0 = 3*i-3
          DO k = 1,3
              L = L+1
              IF (SY(L).EQ.0.0) SY(L) = SR/N*w_new+sy(L)*w_old
          END DO ! k
      END DO  ! I
      L = L+1 ! WTTL
      DO i = 1,N-1
          i0 = 3*i-3
          DO k = 1,3
              L = L+1
              IF (SY(L).EQ.0.0)
     &            SY(L) = (SW/N+SQRT(UPO/mass))*w_new+sy(L)*w_old
c              IF (SY(L).EQ.0.0)SY(L) = 1
          END DO ! k
      END DO ! i

      CMXA = ABS(cmx(1))+ABS(cmx(2))+ABS(cmx(3))+SR/N
      CMVA = ABS(cmv(1))+ABS(cmv(2))+ABS(cmv(3))+SW/N

      DO i = 1,3
          L = L+1
          SY(L) = CMXA*w_new+sy(L)*w_old ! cmx
      END DO

      DO i = 1,3
          L = L+1
          SY(L) = CMVA*w_new+sy(L)*w_old ! cmv
      END DO

      L = L+1
      SY(L) = (ABS(ENERGY)+0.1*UPO)*w_new+sy(L)*w_old ! E
      L = L+1
      SY(L) = SY(L-1)*w_new+sy(L)*w_old
      IF (SY(1).EQ.0.0)
     &    SY(1) = (SQRT(sr/mass)*sr*1.d-2)*w_new+sy(1)*w_old
*       Comment for the line above says 'time'.
      DO k = 1,3
          L = L+1
          SY(L) = 1. ! spin components.
      END DO

      RETURN
      END

      SUBROUTINE EVALUATEX
      INCLUDE 'ARCCOM2e2.ch'
      REAL*8 X0(3)
      SAVE

C       Obtain physical variables from chain quantities.
      DO K = 1,3
          XI(K) = 0.0
          X0(K) = 0.0
      END DO
      DO I = 1,N-1
          L = 3*(I-1)
          DO K = 1,3
              XI(L+3+K) = XI(L+K)+XC(L+K)
          END DO
      END DO
      DO I = 1,N
          L = 3*(I-1)
          DO K = 1,3
              X0(K) = X0(K)+XI(L+K)*MC(I)/MASS
          END DO
      END DO
C       Rearrange according to INAME(i) and add CM.
      DO I = 1,N
          L = 3*(I-1)
          LF = 3*(INAME(I)-1)
          DO K = 1,3
              X(LF+K) = XI(L+K)-X0(K)!+CMX(K) ! CM-coords
          END DO
      END DO

      RETURN
      END

      SUBROUTINE EVALUATEV(VN,WI)
      INCLUDE 'ARCCOM2e2.ch'
      REAL*8 V0(3),VN(*),WI(*)
      SAVE

C       Obtain physical V's from chain quantities.
      DO K = 1,3
          V0(k) = 0.0
          VI(k) = 0.0
      END DO
      DO I = 1,N-1
          L = 3*(I-1)
          DO K = 1,3
              VI(L+3+K) = VI(L+K)+WI(L+K)!WC(L+K)
          END DO
      END DO
      DO I = 1,N
          L = 3*(I-1)
          DO K = 1,3
              V0(K) = V0(K)+VI(L+K)*MC(I)/MASS
          END DO
      END DO
C       Rearrange according to INAME(i) and add CM.
      DO I = 1,N
          L = 3*(I-1)
          LF = 3*(INAME(I)-1)
          DO K = 1,3
              VN(LF+K) = VI(L+K)-V0(K)!+CMV(K)
              V(LF+K) = VN(LF+K) !
          END DO
      END DO

      RETURN
      END

      SUBROUTINE RelativisticAccelerations(ACC,ACCGR,Va,spina,dspin)
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/EXTRA/     KNAME(NMX)
      COMMON/collision/ icollision,IBH,JBH,iwarning
      COMMON/notneeded/ rijnotneeded
      COMMON/deeveet/   dv2(3),dv4(3),dv5(3)
      REAL*8 ACC(*),dX(3),dW(3),dF(3),Va(*),ACCGR(*),dfGR(3),dsp(3),
     &       spina(3),dspin(3)
      SAVE
*
      Cl = Clight! SPEED OF LIGHT
C       INITIALIZE the relativistic acceration(s) here.
      DO  I = 1,3*N
          ACC(I) = 0.0
          ACCGR(I) = 0.0
      END DO
      DO k = 1,3
          dspin(k) = 0
      END DO
      IF (IBH.LE.0) RETURN
*      DO IK = 1,N      ! Seppo's original scheme.
      I = IBH          ! Explicit scheme, 14/6/10.
      IK = KNAME(I)
      I3 = 3*I
      I2 = I3-1
      I1 = I3-2
*      DO JK = IK+1,N
      J = JBH
      JK = KNAME(J)
      J3 = J+J+J
      J2 = J3-1
      J1 = J3-2
      IF (JK.NE.IK+1) THEN
          dx(1) = X(J1)-X(I1)
          dx(2) = X(J2)-X(I2)
          dx(3) = X(J3)-X(I3)
          dw(1) = Va(J1)-Va(I1)
          dw(2) = Va(J2)-Va(I2)
          dw(3) = Va(J3)-Va(I3)
      ELSE
          K1 = 3*IK-2
          K2 = K1+1
          K3 = K2+1
          dx(1) = XC(K1)
          dx(2) = XC(K2)
          dx(3) = XC(K3)
          dw(1) = Va(J1)-Va(I1)
          dw(2) = Va(J2)-Va(I2)
          dw(3) = Va(J3)-Va(I3)
      END IF
      vij2 = dw(1)**2+dw(2)**2+dw(3)**2

c       This (cheating) avoids vij>cl and produces only O(1/c^6) 'errors'.
      DO k = 1,3
          dw(k) = dw(k)/(1+(vij2/cl**2)**4)**.125d0 !  avoid V_ij > c !!
c       dw(k) = dw(k)/(1+(vij2/cl**2)**2)**.25d0 ! not so good
      END DO

      vij2 = dw(1)**2+dw(2)**2+dw(3)**2
      RS = 2.d0*(m(i)+m(j))/CL**2

      RIJ2 = dx(1)**2+dx(2)**2+dx(3)**2
      rij = SQRT(rij2)
      rdotv = dx(1)*dw(1)+dx(2)*dw(2)+dx(3)*dw(3)
*       Ii = min(i,j)  ! activated 18/6 (check with Seppo).
*       Jx = MAX(i,j)
      IF (M(I).GT.M(J))THEN
          Ii = I
          Jx = J
      ELSE
          Ii = J
          Jx = I
      END IF
*       CALL RelativisticTerms  ! replaced by 2010 NBODY7 version.
*    &  (Ii,dX,dW,rij,rdotv,vij2,m(Ii),m(Jx),cl,DF,dfGR,spina,dsp)
      CALL PNPERT
     &     (dX,dW,rij,rdotv,vij2,m(Ii),m(Jx),DF,dfGR,spina,dsp)
*       Note m(Ii), m(Jx) changed from m(i), m(j) 18/6/10 (OK by Seppo).
      RS = 2.d0*(m(i)+m(j))/CL**2
      IF (rij.LT.4.*RS.AND.iwarning.EQ.1) THEN
          WRITE (6,10)  r/RS, rij/RS, I, J, rij, RS
   10     FORMAT (' CHAIN COLLISION    r/RS rij/RS I J rij RS ',
     &                               2F6.2,2I4,1P,2E10.2)
      END IF
      IF (rij.LT.4*RS) THEN!
          iwarning = iwarning+1
          icollision = 1   ! collision indicator
      END IF
      DO k = 1,3
          dspin(k) = dspin(k)+dsp(k)
      END DO
*       Note all signs changed 24/6/10 consistent with Rel Terms and PNPERT.
      ACC(I1) = ACC(I1)-m(j)*dF(1) ! here I assume action = reaction
      ACC(I2) = ACC(I2)-m(j)*dF(2) ! which is not really true for
      ACC(I3) = ACC(I3)-m(j)*dF(3) ! relativistic terms (but who cares)
      ACC(J1) = ACC(J1)+m(i)*dF(1)
      ACC(J2) = ACC(J2)+m(i)*dF(2)
      ACC(J3) = ACC(J3)+m(i)*dF(3)
c        Grav.Rad.-terms split according to masses (allows c.m. velocity).
      ACCgr(I1) = ACCgr(I1)-m(j)*dFgr(1) ! here I assume action = reaction
      ACCgr(I2) = ACCgr(I2)-m(j)*dFgr(2) ! which is not really true for
      ACCgr(I3) = ACCgr(I3)-m(j)*dFgr(3) ! relativistic terms (but who cares)
      ACCgr(J1) = ACCgr(J1)+m(i)*dFgr(1)
      ACCgr(J2) = ACCgr(J2)+m(i)*dFgr(2)
      ACCgr(J3) = ACCgr(J3)+m(i)*dFgr(3)
*
*       END DO ! J  Suppressed after Seppo's new suggestion.
*       END DO ! I

      RETURN
      END

      SUBROUTINE Reduce2cm(x,m,nb,cm)
      IMPLICIT REAL*8 (a-h,m,o-z)
      REAL*8 x(*),m(*),cm(3)
      SAVE

      cm(1) = 0
      cm(2) = 0
      cm(3) = 0
      sm = 0
      DO i = 1,nb
          sm = sm+m(i)
          DO k = 1,3
              cm(k) = cm(k)+m(i)*x(k+3*(i-1))
          END DO
      END DO
      DO k = 1,3
          cm(k) = cm(k)/sm
      END DO
      DO i = 1,nb
          DO k = 1,3
              x(k+3*(i-1)) = x(k+3*(i-1))-cm(k)
          END DO
      END DO

      RETURN
      END

      SUBROUTINE cross2(a,b,c)
      REAL*8 a(3),b(3),c(3)
      SAVE

      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)

      RETURN
      END

      SUBROUTINE gopu_SpinTerms(X,V,r,M1,m2,c,alpha,dv3,dalpha)
      IMPLICIT REAL*8 (a-h,m,n,o-z)
      REAL*8 x(3),v(3),dv3(3),n(3)
      REAL*8 dalpha(3),w(3),alpha(3)
      REAL*8 nxa(3),vxa(3),J(3)
      REAL*8 dv_q(3) ! TEST
      SAVE
       ! This routine assumes: BH mass M1>>m2. Spin of m2 is neglected

      DO k = 1,3
          n(k) = x(k)/r
      END DO
      m = m1+m2
      eta = m1*m2/m**2
      SQ = SQRT(1-4*eta)
      Aq = -12/(1+sq)
      Bq =  -6/(1+sq)-3
      Cq = 1+6/(1+sq)
      rdot = cdot(n,v)
      CALL cross2(n,v,w)
      anxv = cdot(alpha,w)
      CALL cross2(n,alpha,nxa)
      CALL cross2(v,alpha,vxa)
      DO k = 1,3
          dv3(k) = -m1**2/(c*r)**3*
     &             (Aq*anxv*n(k)+rdot*Bq*nxa(k)+Cq*vxa(k))
      END DO
      coeff = eta*m/(c*r)**2*(3/(1+sq)+.5d0)
      CALL cross2(w,alpha,dalpha)
      DO k = 1,3
          dalpha(k) = coeff*dalpha(k)
      END DO
c       C. Will Q2-terms
      sjj = 0
      DO k = 1,3
          j(k) = M1**2/c*alpha(k)
          sjj = sjj+j(k)**2
      END DO
      sj = SQRT(sjj)
      IF (sj.ne.0.0) THEN  ! IF sj = 0, THEN J(k) = 0 and Q-term  = 0 anyway
          DO k = 1,3
              j(k) = j(k)/sj
          END DO
      END IF
      Q2 = -sjj/M1/c**2!  X = X_j-X_i in this code

c      DO k = 1,3 ! OLD (and wrong!)
c      dv3(k) = dv3(k)  ! add Quadrupole terms
c     & +1.5*Q2/r**4*(n(k)*(5*cdot(n,j)**2-1)-2*cdot(n,j)*j(k))
c      END DO

      Q2 = -Q2 ! earlier we had Q2 grad Q-Potential,
           ! now grad Q-ForceFunction  = > different sign
      CALL Q2term(m1,r,x,v,c,Q2,j,dv_q) ! Obtain the Q2 terms
      DO k = 1,3
          dv3(k) = dv3(k)+dv_q(k) ! add quadrupole terms (these are more correct)
      END DO

      RETURN
      END

      SUBROUTINE Q2term(m,r,x,v,c,Q2,e,dvq)
      IMPLICIT REAL*8 (a-h,m,o-z)
      REAL*8 x(3),v(3),dvq(3),Rx(3),Ux(3),e(3)
      ! m = m1+m2 (?),vv = v**2
      ! e = spin direction;  Q2 = m**3/c**4*xi**2, xi = |spin| = Kerr PARAMETER

      vv = cdot(v,v)
      er = cdot(e,x)
      RQ2 = (-1+3*(er/r)**2)/(2*r**3) ! quadrupole potential (exept factor Q2)
      U2b = m/r
      oc = 1/c
      DO k = 1,3
          Ux(k) = -x(k)*m/r**3 ! two-body acceleration
          Rx(k) = (3*e(k)*er)/r**5+
     &     (x(k)*(-3*er**2/r**6-(3*(-1+(3*(er)**2)/r**2))/(2*r**4)))/r ! quadrupole pot gradient
      END DO
      vRx = cdot(v,Rx)
      DO k = 1,3 ! complete quadrupole term in \dot v
          dvq(k) = Q2*(Rx(k)*(1 + oc**2*(-4*(Q2*RQ2 + U2b) + vv))
     &             -4*oc**2*(RQ2*Ux(k)+vRx*v(k)))
      END DO

      RETURN
      END

      SUBROUTINE InitialStepsize(X,V,M,NB,ee,step)
      IMPLICIT REAL*8 (A-H,m,O-Z)
      DIMENSION X(*),V(*),M(*)
      SAVE

      T = 0.0
      U = 0.0
      RMIN = 1.D30
      mass = M(NB)
      time_step2 = 1.e30
      DO I = 1,NB-1
          mass = mass+M(I)
          DO J = I+1,Nb
              MIJ = M(I)*M(J)
              KI = (I-1)*3
              KJ = (J-1)*3
              xx = X(KI+1)-X(KJ+1)
              yy = X(KI+2)-X(KJ+2)
              zz = X(KI+3)-X(KJ+3)
              R2 = xx*xx+yy*yy+zz*zz+ee
              vx = V(KI+1)-V(KJ+1)
              vy = V(KI+2)-V(KJ+2)
              vz = V(KI+3)-V(KJ+3)
              vv = vx*vx+vy*vy+vz*vz
              R1 = Sqrt(R2)
              time_step2 = min(time_step2,R2/(vv+(M(I)+M(J))/R1)) ! ~2B radius of convergence^2
              U = U+MIJ/R1
              T = T+MIJ*(vx*vx+vy*vy+vz*vz)
          END DO
      END DO
      T = T/(2*mass)
      ENERGY = T-U
      Alag = T+U
      STEP = 0.1*U*SQRT(time_step2)

      RETURN
      END

      function oot(alfa,eta,zeta,q,e,sqaf) ! oot = pericentre time
c       alfa = 1/a; eta = sqrt(a) e sin(E); zeta = e Cos(E),
c       q = a(1-e), e = ecc, sqaf = sqrt(|a|)
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER(tiny = 1.d-18)
      SAVE

      IF (zeta.GT.0.0) THEN
c        ellipse (near peri), parabola or hyperbola.
          ecc = MAX(e,tiny)
          X = eta/ecc
          Z = alfa*X*X
          oot = X*(q+X*X*g3(Z))
      ELSE
c       upper half of an elliptic orbit.
          oot = (atan2(eta*sqaf,zeta)/sqaf-eta)/alfa
      END IF

      RETURN
      END

      function g3(z)
      IMPLICIT REAL*8 (a-h,o-z)
      COMMON/mita/ zero
      SAVE

      IF (z.GT.0.025d0) THEN ! elliptic
          x = SQRT(z)
          g3 = (asin(x)-x)/x**3
      ELSEIF(z.LT.-0.025d0) THEN ! hyperbolic
          x = SQRT(-z)
          g3 = (log(x+SQRT(1+x*x))-x )/x/z
      ELSE ! Pade approximation for small  |z|
c       g3 = (1/6.d0-19177*z/170280 + 939109*z*z/214552800)/
c     &      (1-7987*z/7095 + 54145*z*z/204336)
          g3 = (1+6*(-19177*z/170280 + 939109*z*z/214552800))/
     &         (6*(1-7987*z/7095 + 54145*z*z/204336))
          zero = 0
      END IF

      RETURN
      END

      SUBROUTINE ConstantsOfMotion(ENE_NB,G,Alag)
c       IMPLICIT REAL*8 (A-H,m,O-Z)
c       DIMENSION G(3)
      INCLUDE 'ARCCOM2e2.ch'
      REAL*8 g(3)
      COMMON/justforfun/ Tkin,Upot,dSkin,dSpot
      SAVE

c       Contants of motion in the centre-of-mass system.
      T = 0.0
      U = 0.0
      G(1) = 0.
      G(2) = 0.
      G(3) = 0.
      RMIN = 1.D30
c       mass = M(N)
      DO Ik = 1,N-1
          I = INAME(IK)      ! along the chain
c       mass = mass+M(I)
          DO Jk = Ik+1,N
              J = INAME(JK)      !  -"-
              MIJ = M(I)*M(J)
              KI = (I-1)*3
              KJ = (J-1)*3
              IF (JK.NE.IK+1) THEN
                  xx = X(KI+1)-X(KJ+1)
                  yy = X(KI+2)-X(KJ+2)
                  zz = X(KI+3)-X(KJ+3)
                  vx = V(KI+1)-V(KJ+1)
                  vy = V(KI+2)-V(KJ+2)
                  vz = V(KI+3)-V(KJ+3)
              ELSE
                  K1 = 3*IK-2
                  K2 = K1+1
                  K3 = K2+1
                  XX = XC(K1)   ! use chain vectors when possible,
                  YY = XC(K2)   ! this often reduces roundoff
                  ZZ = XC(K3)
                  VX = WC(K1)
                  VY = WC(K2)
                  VZ = WC(K3)
              END IF
              R2 = xx*xx+yy*yy+zz*zz+ee
              U = U+MIJ/SQRT(R2)
              T = T+MIJ*(vx*vx+vy*vy+vz*vz)
              G(1) = G(1)+MIJ*(yy*vz-zz*vy)
              G(2) = G(2)+MIJ*(zz*vx-xx*vz)
              G(3) = G(3)+MIJ*(xx*vy-yy*vx)
          END DO
      END DO
      T = T/(2*mass)
      G(1) = G(1)/mass
      G(2) = G(2)/mass
      G(3) = G(3)/mass
      ENE_NB = T-U
      Alag = T+U
      Tkin = T ! to justforfun
      Upot = U ! to justforfun
      OmegaB = Wfunction()
      dSkin = cmethod(1)*(T-ENERGY-ENERGR)+cmethod(2)*WTTL+cmethod(3)
      dSpot = cmethod(1)*U+cmethod(2)*OmegaB+cmethod(3)

      RETURN
      END

      SUBROUTINE  WCMOTION(hs)
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/IncrementCommon/   WTTLinc,XCinc(NMX3),WCinc(NMX3),
     &                          CMXinc(3),CMVinc(3),ENERGYinc,
     &                          Energrinc,CHTIMEinc,spin inc(3)
      COMMON/vwcommon/          Ww(nmx3),WTTLw,cmvw(3),spinw(3)
      COMMON/omegacoefficients/ OMEC(NMX,NMX)
      COMMON/apuindex/          ikir
      COMMON/DerOfTime/         G
      COMMON/DIAGNOSTICS/       GAMMA,H,IWR
      REAL*8 FC(NMX3),XAUX(3),acc(nmx3)
      REAL*8 F(NMX3),!df(nmx3),dfGR(nmx3),
     &  GOM(nmx3)!,dcmv(3),Va(nmx3),afc(nmx3),dfE(3),dspin(3)
      DATA IT /0/
      SAVE

      CALL EVALUATEX
      RSUM = 0.0
      OMEGA = 0.0d0
      U = 0
      DO i = 1,3*N
          f(i) = 0
          GOM(i) = 0
      END DO
      DO I = 1,N-1
          L = 3*(I-1)
          RIJL2 = xc(L+1)**2+xc(L+2)**2+xc(L+3)**2+ee
          RIJL = SQRT(RIJL2)
C        Evaluate RSUM for decision-making.
          RSUM = RSUM+RIJL
          RINV(I) = 1.d0/RIJL
          U = U+MC(I)*MC(I+1)*RINV(I)
          A = RINV(I)**3
          i0 = 3*i-3
          j = i+1
          j0 = 3*j-3
          omeker = omec(iname(i),iname(j))
          DO K = 1,3
              AF = A*XC(I0+K)
              f(I0+k) = f(i0+k)+MC(J)*AF
              f(j0+k) = f(j0+k)-MC(I)*AF
              IF (cmethod(2).ne.0.0d0.AND.omeker.ne.0.0) THEN
                  GOM(I0+k) = GOM(I0+k)+AF*omeker
                  GOM(J0+k) = GOM(J0+k)-AF*omeker
              END IF
          END DO
          IF (cmethod(2).ne.0.0.AND.omeker.ne.0.0) THEN
              OMEGA = OMEGA+omeker*RINV(I)
          END IF
      END DO

      LRI = N-1
C       Physical coordinates
      DO K = 1,3
          XI(K) = 0.0
      END DO
      DO I = 1,N-1
          L = 3*(I-1)
          DO K = 1,3
              XI(L+3+K) = XI(L+K)+XC(L+K)
          END DO
      END DO
C       Non-chained contribution
      DO I = 1,N-2
          LI = 3*(I-1)
          DO J = I+2,N
              LJ = 3*(J-1)
              RIJ2 = 0.0+ee
              IF (J.GT.I+2) THEN
                  DO K = 1,3
                      XAUX(K) = XI(LJ+K)-XI(LI+K)
                      RIJ2 = RIJ2+XAUX(K)**2
                  END DO
              ELSE
                  DO K = 1,3
                      XAUX(K) = XC(LI+K)+XC(LI+K+3)
                      RIJ2 = RIJ2+XAUX(K)**2
                  END DO
              END IF
              RIJ2INV = 1/RIJ2
              LRI = LRI+1
              RINV(LRI) = SQRT(RIJ2INV)
              U = U+MC(I)*MC(J)*RINV(LRI)
              omeker = omec(iname(i),iname(j))
              IF (omeker.ne.0.0.AND.cmethod(2).ne.0.0) THEN
                  OMEGA = OMEGA+omeker*RINV(LRI)
              END IF
              DO K = 1,3
                  A = RINV(LRI)**3*XAUX(K)
                  f(LI+K) = f(LI+K)+MC(J)*A
                  f(LJ+K) = f(LJ+K)-MC(I)*A
                  IF (cmethod(2).ne.0.0d0.AND.omeker.ne.0.0) THEN
                      GOM(LI+K) = GOM(LI+K)+A*omeker
                      GOM(LJ+K) = GOM(LJ+K)-A*omeker
                  END IF
              END DO
          END DO ! J = I+2,N
      END DO  ! I = 1,N-2
      dT = hs/(U*cmethod(1)+OMEGA*cmethod(2)+cmethod(3)) ! time interval
*       Obtain Newtonian perturbation and copy to new labelled common (09/16).
      CALL XTPERT(ACC)
      DO i = 1,n-1
          DO k = 1,3
              L = 3*(i-1)
              FC(L+k) = f(3*i+k)-f(3*i+k-3)
          END DO
      END DO
      IF (clight.GT.0.0) THEN ! V-dependent ACC
          CALL V_jump(Ww,spinw,cmvw,WTTLw,WC,spin,FC,acc,dt/2,
     &                gom,energyj,energrj,1) ! Auxiliary W ( = Ww) etc
          CALL V_jump(WC,spin,cmv,WTTL,Ww,spinw,FC,acc,dt,
     &                gom,energy,energr,2)   ! 'true' W  etc
          CALL V_jump(Ww,spinw,cmvw,WTTLw,WC,spin,FC,acc,dt/2,
     &                gom,energyj,energrj,3) ! Auxiliary W ( = Ww) etc
      ELSE ! c>0
          CALL V_jACConly(WC,cmv,WTTL,FC,acc,dt,
     &                    gom,energy,energrj)  ! here ACC depends ONLY on COORDINATES
      END IF

      RETURN
      END

      SUBROUTINE V_jump(WCj,spinj,cmvj,wttlj,WCi,spini,FCj,acc,dt,
     &                  gom,energyj,energrj,ind)
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/IncrementCommon/ WTTLinc,XCinc(NMX3),WCinc(NMX3),
     &                        CMXinc(3),CMVinc(3),ENERGYinc,
     &                        Energrinc,CHTIMEinc,spin inc(3)
      REAL*8 wcj(*),fcj(*),df(nmx3),dcmv(3),afc(nmx3),gom(*),
     &  dfe(nmx3),dfgr(nmx3),dspin(3),spinj(3),cmvj(3),wci(nmx3),
     &  spini(3),acc(*)
      SAVE

      CALL EVALUATEV(V,WCi)
c       add V-dependent perts.
      IF (clight.GT.0.0) THEN
          CALL VelocityDependentPerturbations
     &         (dT,V,spini,acc,dcmv,df,dfGR,dspin)
      ELSE
          DO i = 1,3*n
              df(i) = acc(i)
          END DO
      END IF
      DO i = 1,n-1
          L = 3*I-3
          I1 = 3*INAME(I)-3
          I2 = 3*INAME(I+1)-3
          DO k = 1,3
              afc(L+k) = df(I2+k)-df(I1+k)
          END DO
      END DO
      IF (IND.EQ.2) THEN
          dotE = 0
          dotEGR = 0
          DO I = 1,N
              I0 = 3*I-3
              DO k = 1,3
                  dfE(k) = df(i0+k)!  +dfGR(i0+k)!
              END DO
              dotE = dotE+! NB-Energy change.
     &          M(I)*(V(I0+1)*dfE(1)+V(I0+2)*dfE(2)+V(I0+3)*dfE(3)) !
                DO k = 1,3
                    dfE(k) = dfGR(I0+k)
                END DO
                dotEGR = dotEGR+ ! radiated energy
     &          M(I)*(V(I0+1)*dfE(1)+V(I0+2)*dfE(2)+V(I0+3)*dfE(3))
          END DO
**          dotE = dotE+dotEGR
          ENERGYj = ENERGYj+dotE*dT
            EnerGrj = EnerGRj+dotEGR*dT
          IF (ind.EQ.2) THEN
              ENERGYinc = ENERGYinc+dotE*dt
                EnerGRinc = EnerGRinc+dotEGR*dT
          END IF !ind.EQ.2
      END IF ! IND = 2
      IF (cmethod(2).ne.0.0d0) THEN
          dotW = 0
          DO I = 1,N
              k0 = 3*I-3
              i0 = 3*iname(i)-3
              dotW = dotW+
     &          (V(I0+1)*GOM(k0+1)+V(I0+2)*GOM(K0+2)+V(I0+3)*GOM(K0+3))
          END DO
          WTTLj = WTTLj+dotW*dT
          IF (ind.EQ.2) WTTLinc = WTTLinc+dotW*dT
      END IF ! cmethod(2).ne.0.0
      DO I = 1,N-1
          L = 3*(I-1)
          DO K = 1,3
              IF (ind.EQ.2)
     &          WCinc(L+K) = WCinc(L+K)+(FCj(L+K)+afc(L+K))*dT
              WCj(L+K) = WCj(L+K)+(FCj(L+K)+afc(L+K))*dT
          END DO
      END DO

      DO k = 1,3
          spinj(k) = spinj(k)+dT*dspin(k)
          cmvj(k) = cmvj(k)+dT*dcmv(k)
      END DO
      IF (ind.EQ.2) THEN
          DO k = 1,3
              spin inc(k) = spin inc(k)+dT*dspin(k)
              cmv inc(k) = cmv inc(k)+dT*dcmv(k)
          END DO
      END IF ! ind.EQ.2

      RETURN
      END

      SUBROUTINE V_jACConly(WCj,CMVj,WTTLj,FC,acc,dt,
     &                      gom,energyj,energrj)
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/IncrementCommon/ WTTLinc,XCinc(NMX3),WCinc(NMX3),
     &                        CMXinc(3),CMVinc(3),ENERGYinc,
     &                        Energrinc,CHTIMEinc,spin inc(3)
      REAL*8 wcj(*),fc(*),dcmv(3),afc(nmx3),gom(*),
     &  dfe(nmx3),cmvj(3),acc(*),WCi(NMX3)
      SAVE

      DO I = 1,N-1
          L = 3*(I-1)
          DO K = 1,3
              WCi(L+K) = WC(L+K)+FC(L+K)*dT/2  !( no inc here!)
          END DO
      END DO
      CALL EVALUATEV(V,WCi)
      CALL Reduce2cm(acc,m,n,dcmv)
      DO I = 1,3*N
          V(I) = V(I)+acc(I)*dT/2 ! average Velocity
      END DO
c       add V-dependent perts.
      DO i = 1,n-1
          L = 3*I-3
          I1 = 3*INAME(I)-3
          I2 = 3*INAME(I+1)-3
          DO k = 1,3
              afc(L+k) = acc(I2+k)-acc(I1+k) ! CHAIN vector accelerations
          END DO
      END DO
      dotE = 0
      dotEGR = 0
      DO I = 1,N
          I1 = 3*I-2
          DO k = 1,3
              dfE(k) = acc(i0+k)
          END DO
          dotE = dotE+M(I)*cdot(V(I1),acc(i1))
      END DO
      ENERGYj = ENERGYj+dotE*dT
        EnerGrj = EnerGRj+dotEGR*dT
      ENERGYinc = ENERGYinc+dotE*dT
        EnerGRinc = EnerGRinc+dotEGR*dT
      IF (cmethod(2).ne.0.0d0) THEN
          dotW = 0
          DO I = 1,N
              k1 = 3*I-2
              i1 = 3*iname(i)-2
              dotW = dotW+cdot(V(I1),GOM(K1))
          END DO
          WTTLinc = WTTLinc+dotW*dT
          WTTLj = WTTLj+dotW*dT
      END IF ! cmethod(2).ne.0.0
      DO I = 1,N-1
          L = 3*(I-1)
          DO K = 1,3
              WCinc(L+K) = WCinc(L+K)+(FC(L+K)+afc(L+K))*dT
              WCj(L+K) = WCj(L+K)+(FC(L+K)+afc(L+K))*dT
          END DO
      END DO

      DO k = 1,3
          cmv inc(k) = cmv inc(k)+dT*dcmv(k)
          cmvj(k) = cmvj(k)+dT*dcmv(k)
      END DO

      RETURN
      END

      SUBROUTINE EstimateStepsize(dtime,step)
      INCLUDE 'ARCCOM2e2.ch'
      COMMON/collision/         icollision,ione,itwo,iwarning
      COMMON/omegacoefficients/ OMEC(NMX,NMX) ! not part of ARCCOM2e2.h
      COMMON/eitoimi/           iei
      COMMON/toolarge/          beta,ma,mb,itoo,iw,jw,n_alku
      PARAMETER (twopi = 6.283185307179586d0)
      REAL*8 xij(3),vij(3),gx(5)
      SAVE

      nr = 0
      nx = 0
c      evaluate lenght of chain
      CALL UpdateXandV  ! we need x and v
      step = cmethod(3)*dtime   ! contribution from cmethod(3)
      DO IK = 1,N-1
          DO JK = IK+1,N
              I = INAME(IK)
              J = INAME(JK)
              iw = i
              jw = j
              KI = (I-1)*3
              KJ = (J-1)*3
              IF (JK.NE.IK+1) THEN
                  xij(1) = X(KI+1)-X(KJ+1)
                  xij(2) = X(KI+2)-X(KJ+2)
                  xij(3) = X(KI+3)-X(KJ+3)
                  vij(1) = V(KI+1)-V(KJ+1)
                  vij(2) = V(KI+2)-V(KJ+2)
                  vij(3) = V(KI+3)-V(KJ+3)
                  ind = 0
              ELSE
                  ind = 123
                  K1 = 3*IK-2
                  K2 = K1+1
                  K3 = K2+1
                  xij(1) = -XC(K1)   ! use chain vectors when possible,
                  xij(2) = -XC(K2)   ! this often reduces roundoff
                  xij(3) = -XC(K3)
                  vij(1) = -WC(K1)
                  vij(2) = -WC(K2)
                  vij(3) = -WC(K3)
              END IF

              i0 = 3*i-3
              j0 = 3*j-3
              DO k = 1,3
                  xijk = x(i0+k)-x(j0+k)
                  vijk = v(i0+k)-v(j0+k)
              END DO
              rr = cdot(xij,xij)
              r = SQRT(rr)
              alfa = cmethod(1)*m(i)*m(j)+cmethod(2)*OMEC(I,J) ! potential & 'TTL' terms

              mipj = m(i)+m(j)
              vv = cdot(vij,vij)
              oa = 2/r-vv/mipj
              dltrr = dtime**2*vv
              IF (dltrr.LT..001*rr) THEN
                  nr = nr+1
                  step = step+dtime*alfa/r ! add contributions from large distances
              ELSE                   ! in this case use Stumpff-Weiss method
                  nx = nx+1
                  eta = cdot(xij,vij)
                  beta = mipj*oa
                  zeta = mipj-beta*r
                  period = 0
                  IF (oa.GT.0.0) THEN
                      period = twopi/(oa*SQRT(oa*mipj))
                      kp = INT(dtime/period)
                      delta_t = dtime-kp*period ! take periods into account differently
                  ELSE
                      kp = 0
                      delta_t = dtime !!!
                  END IF
                  ma = m(i)
                  mb = m(j)
                  Xa = 0
                  CALL Xanom(mipj,r,eta,zeta,beta,delta_t,Xa,rx,gx) ! Solve KEPLER-eqs.
                  step = step+alfa*(Xa+oa*kp*period) ! The Stumpff-Weiss principle.
              END IF
          END DO
      END DO

      RETURN
      END

      SUBROUTINE gfunc(xb,al,g)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 c(5),g(5)

      z = al*xb*xb
      CALL cfun(z, c)
      s = xb
      DO 1 i = 1,5
          g(i) = c(i)*s
          s = s*xb
    1 CONTINUE

      RETURN
      END

      SUBROUTINE cfun(z,c)!Stumpff(Z,C)
      IMPLICIT REAL*8 (A-H,m,O-Z)
      COMMON/toolarge/ beta,ma,mb,itoo,iw,jw,n_alku
      COMMON/diagno/   ncfunc
      PARAMETER(o2 = 1.d0/2,o6 = 1.d0/6,o8 = 1.d0/8,o16 = 1.d0/16)
      REAL*8 C(5)
      SAVE

      ncfunc = ncfunc+1
      itoo = 0
      h = z
      DO  K = 0,7
          IF (ABS(h).LT.0.9d0) GOTO 2
          h = h/4 ! divide by 4 untill h<.9
      END DO
      akseli = (ma+mb)/beta
      WRITE (6,106) Z,iw,jw,ma,mb,beta,akseli,n_alku
      WRITE (41,106) Z,iw,jw,ma,mb,beta,akseli,n_alku
  106 FORMAT(' too large Z = ',1p,g12.4, '4 c-functions',
     &       0p,2i5,1p,4g12.4,i5,' ijmab_beta ax n_a')

      c(1) = 0!1.
      DO k = 2,5
          c(k) = 0!c(k-1)/k ! something
      END DO
      itoo = 1
      RETURN
 2    C(4) =     ! use Pade-approximations for c_4 & c_5
     &        (201859257600.d0+h*(-3741257520.d0
     &      + (40025040.d0-147173.d0*h)*h))/
     &        (240.d0*(20185925760.d0 + h*(298738440.d0
     &      + h*(1945020.d0 + 5801.d0*h))))
      C(5) =
     &       (3750361655040.d0 + h*(-40967886960.d0
     &     + (358614256.d0 - 1029037.d0*h)*h))/
     &       (55440.d0*(8117665920.d0 + h*(104602680.d0
     &     + h*(582348.d0 + 1451.d0*h))))

      DO  I = 1,K  ! 4-fold argument K times
          C3 = o6-h*C(5)
          C2 = o2-h*C(4)
          C(5) = (C(5)+C(4)+C2*C3)*o16
          C(4) = C3*(2.D0-h*C3)*o8
          h = 4.d0*h
      END DO

      C(3) = o6-Z*C(5)
      C(2) = o2-Z*C(4)
      C(1) = 1-Z*C(3)

      RETURN
      END

      SUBROUTINE Xanom(m,r,eta,zet,beta,t,x,rx,g)
c      Kepler solver.
      IMPLICIT REAL*8 (a-h,m,o-z)
      COMMON/diagno/    ncfunc
      COMMON/collision/ icollision,ione,itwo,iwarning
      COMMON/eitoimi/   iei
      COMMON/toolarge/  betaa,ma,mb,itoo,iw,jw,n_alku
      REAL*8 g(5)
c      Solution of the `universal' form of Kepler's equation.
c      input: m = mass, r  = r(0) = dist, eta = r.v, zet = m-r*beta, beta = m/a, t = time-incr
c      { note:  eta = sqrt[m a]*e Sin[E],  zeta = m e Cos[E] }
c      output: x = \int dt/r, rx = r(t), g(k) = x^k*c_k(beta*x^2); c_k = Stumpff-funcs
c      recommend: IF a fairly good initial estimate is not available, use X = 0.
      SAVE

      betaa = beta
      iei = 0
      IF (t.EQ.0.0) THEN ! IF CALLed with t = 0
          x = 0
          DO k = 1,5
              g(k) = 0
          END DO
          rx = r
          RETURN
      END IF

c        initial estimate (IF not given as input i.e. IF not x*t>0 )
      IF (x*t.LE.0.0) THEN ! no initial estimate
          IF (zet.GT.0.0) THEN ! near pericentre
c              x = t/(r**3+m*t**2/6)**.333333333d0
              X = t/SQRT(r*r+(m*t**2/6)**.666667d0)
              Xens = X
          ELSE ! far from pericentre
              x = t/r
          END IF
      END IF

c        first bracket the root by stepping forwards using difference eqs.
      n_alku = 0
   66 r0 = r
      n_alku = n_alku+1
      eta0 = eta
      zet0 = zet
      tau0 = -t
      CALL gfunc(x,beta,g) ! 1.
      IF (itoo.EQ.1) WRITE (6,16) Xens,X,n_alku,iw,jw,
     &                           beta*x*x,tau0,tau1,beta
   16 FORMAT (1x,1p,2g12.4,3i5,1p,4g12.4,' X0 X n i j Z tau0 tau1 beta')
      xg = x
      g0 = 1-beta*g(2)
      tau1 = r0*x+eta0*g(2)+zet0*g(3)-t
      r1 = r0+eta0*g(1)+zet0*g(2)
      eta1 = eta0*g0+zet0*g(1)
      zet1 = zet0*g0-beta*eta0*g(1)
      x0 = 0
      x1 = x
      hhc2 = 2*g(2)
      DO k = 1,8 !!!!!!!!!!!!!
          IF (tau0*tau1.GT.0.0) THEN
              ddtau = hhc2*eta1
              ddr = hhc2*zet1
              r2 = 2*r1-r0+ddr
              zet2 = 2*zet1-zet0-beta*ddr
              tau2 = 2*tau1-tau0+ddtau
              eta2 = 2*eta1-eta0-beta*ddtau
              eta0 = eta1
              eta1 = eta2
              zet0 = zet1
              zet1 = zet2
              r0 = r1
              r1 = r2
              tau0 = tau1
              tau1 = tau2
              x0 = x1
              x1 = x1+x
          ELSE
              GOTO 77
          END IF
      END DO
      x = 1.5d0*x1
      GOTO 66 ! initial estimate much too small!
   77 CONTINUE
c       iterate to final solution
      dx = x
      DO i = 1,300 ! usually i_max  = 2 or 3 only
          itera = i
          IF (ABS(tau0*r1).LT.ABS(tau1*r0))THEN
              dx = -tau0/r0
c              dx = -tau0/(r0+eta0*dx/2)
c              dx = -tau0/(r0+eta0*dx/2+zet0*dx*dx/6)
              x = x0+dx
              dzeit = dx*(r0+eta0*dx/2+zet0*dx*dx/6)+tau0
              x00 = x0
              icase = 0
              tau = tau0
          ELSE
              dx = -tau1/r1
c              dx = -tau1/(r1+eta1*dx/2)
c              dx = -tau1/(r1+eta1*dx/2+zet1*dx*dx/6)
              x = x1+dx
              dzeit = dx*(r1+eta1*dx/2+zet1*dx*dx/6)+tau1
              x00 = x1
              icase = 1
              tau = tau1
          END IF

          IF((x1-x)*(x-x0).LT.0.0.OR.i.EQ.i/5*5) THEN ! IF out_of_brackets or
              x = (x0+x1)/2                                ! slow use bisection
              icase = -1
              GOTO 11
          END IF

          IF (ABS(dzeit).LT.1.d-3*ABS(t)
     &         .AND.ABS(dx).LT.1.e-3*ABS(x)) GOTO 99
   11     CONTINUE
          CALL gfunc(x,beta,g) ! 2.,...
          xg = x
          g0 = 1-beta*g(2)
          rpr = eta*g0+zet*g(1)
          rpp = zet*g0-beta*eta*g(1)
          rlast = r+eta*g(1)+zet*g(2)
          f = r*x+eta*g(2)+zet*g(3)-t

          IF (f*tau0.GT.0.0) THEN ! keep it bracketed
              x0 = x
              tau0 = f
              eta0 = rpr
              zet0 = rpp
              r0 = rlast
          ELSE
              x1 = x
              tau1 = f
              eta1 = rpr
              zet1 = rpp
              r1 = rlast
          END IF
      END DO ! i
      aks = m/beta
      periodi = 6.28*aks*SQRT(ABS(aks)/m)
      WRITE (6,166)aks,r0,r1,t,periodi,x,f/(r0+r1)*2
  166 FORMAT (1x,'NO CONV',1p,7g12.4,' a r0 r1 t prd x dx')
      iei = 1
   99 CONTINUE
c       final correction of g's  & r-evaluation
      IF (X00.ne.xg) THEN
          CALL gfunc(x,beta,g)
          xg = x
      ELSE
          g(5) = g(5)+dx*(g(4)+dx*g(3)/2.d0)
          g(4) = g(4)+dx*(g(3)+dx*g(2)/2.d0)
          g(3) = x**3/6.d0-beta*g(5)
          g(2) = x**2/2.d0-beta*g(4)
          g(1) = x -beta*g(3)
      END IF
      rx = r+eta*g(1)+zet*g(2)

      RETURN
      END

      function cdot(a,b)
      REAL*8  a(3),b(3),cdot

      cdot = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

      RETURN
      END
