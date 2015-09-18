      SUBROUTINE EREL(IM,EB,SEMI)
*
*
*       Dominant two-body energy in chain system.
*       -----------------------------------------
*
*     INCLUDE 'COMARC.CH'
      INCLUDE 'ARCCOM2e2.CH'
*
*
*       Determine indices of dominant motion.
      RX = 0.0
      DO 1 K = 1,N-1
          IF (RINV(K).GT.RX) THEN
              RX = RINV(K)
              IM = K
          END IF
    1 CONTINUE
      K1 = INAME(IM)
      K2 = INAME(IM+1)
*
*       Form kinetic energy terms of dominant c.m. (K1 + K2) and the rest.
      ZK = 0.0D0
      P0 = 0.0D0
      DO 10 K = 1,3
          PS = 0.0
*       Exclude dominant bodies using names of chain members.
          DO 5 I = 1,N
              IF (INAME(I).NE.K1.AND.INAME(I).NE.K2) THEN
                  J = 3*(I - 1)
                  PS = PS + M(I)*V(J+K)
                  ZK = ZK + M(I)*V(J+K)**2
              END IF
    5     CONTINUE
          P0 = P0 + PS**2
   10 CONTINUE
*
*       Obtain potential energy due to non-dominant chain distances.
      POT = 0.0D0
      DO 20 I = 1,N-1
          LI = 3*(I - 1)
          DO 15 J = I+1,N
              IF (I.EQ.K1.AND.J.EQ.K2) GO TO 15
              LJ = 3*(J - 1)
              RIJ2 = 0.0
              DO 12 K = 1,3
                  RIJ2 = RIJ2 + (X(LI+K) - X(LJ+K))**2
   12         CONTINUE
              POT = POT + M(I)*M(J)/SQRT(RIJ2)
   15     CONTINUE
   20 CONTINUE
*
*       Derive binding energy from total energy and perturbing function.
      MB = M(K1) + M(K2)
      VP = 0.5D0*(P0/MB + ZK) - POT
      EB = ENERGY - VP
*
*       Set semi-major axis.
      SEMI = -M(K1)*M(K2)/(2.0D0*EB)
*
      RETURN
*
      END
