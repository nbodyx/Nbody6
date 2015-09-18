      REAL*8 FUNCTION PREMSF(M,T,RR)        
      REAL*8 RZAMS,M,T,RR
      REAL*8 RZAMSF
      EXTERNAL RZAMSF
*      
      RZAMS = RZAMSF(M)
*       Coefficients in PREMSR
      IF(M.LE.1.d0)THEN
         PRE1 = 0d0
         PRE2 = 0d0
         PRE3 = 7.432d-2 - 9.43d-2*M + 7.439d-2*M**2
      ENDIF
      IF(M.GT.1d0.AND.M.LT.2d0)THEN
         PRE1 = -4.00772d0 + 4.00772d0*M
         PRE2 = 8.5656d0 - 8.5656d0*M
         PRE3 = -4.50678d0 + 4.56118d0*M
      ENDIF
      IF(M.GE.2d0)THEN
         PRE1 = 1.60324d0 + 2.20401d0*M - 0.60433d0*M**2 
     &        + 5.172d-2*M**3
         PRE2 = -4.56878d0 - 4.05305d0*M + 1.24575*M**2 
     &        - 0.10922d0*M**3
         PRE3 = 3.01153 + 1.85745*M - 0.64290d0*M**2 
     &        + 5.759d-2*M**3
      END IF
      PREMSF = RZAMS*10**((PRE1*T**3 + PRE2*T**4 
     &       + PRE3*T**5)/(1.05 - T)) - RR
      RETURN
      END
