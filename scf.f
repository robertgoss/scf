      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /MATHS  / PI
      PI = 3.1415926535898d0
      
      IOP   = 2
      N     = 3
      R     = 1.4632D0
      Z1    = 2.0925D0
      Z2    = 1.24D0
      ZA    = 2.0D0
      ZB    = 1.0D0
      
      CALL HFCALC(IOP,N,R,Z1,Z2,ZA,ZB)
      END
     
CCCCCC************************************************************
      SUBROUTINE HFCALC(IOP,N,R,Z1,Z2,ZA,ZB)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /MATHS  / PI
      
      IF (IOP.NE.0) PRINT 10, N, ZA, ZB
      CALL INTEGRALS(IOP,N,R,Z1,Z2,ZA,ZB)
      CALL COLLECT() 
      CALL SCF(IOP,N,R,Z1,Z2,ZA,ZB) 
   
10    FORMAT(1H1,2X,4HSTO-,I1,21HG FOR ATOMIC NUMBERS , F5.2, 5H AND , 
     $F5.2//)
    
      RETURN 
      END
CCCCCC************************************************************
      SUBROUTINE INTEGRALS(IOP,N,R,Z1,Z2,ZA,ZB)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /INTMAT / S(2,2), T(2,2), V(2,2,2), Q(2,2,2,2) 
      COMMON /MATHS  / PI
    
      DIMENSION COEFF(3,3), EXPON(3,3), D(3,2), A(3,2)      
     
      DATA COEFF /1.0D0,      0.0D0,      0.0D0, 
     $            0.678914D0, 0.430129D0, 0.0D0, 
     $            0.444635D0, 0.535328D0, 0.154329D0 /
    
      DATA EXPON /0.270950D0,      0.0D0, 0.0D0, 
     $            0.151623D0, 0.851819D0, 0.0D0, 
     $            0.109818D0, 0.405771D0, 2.22766D0 /
         
      S(1,2)  = 0.0D0  
     
      T(1,1)  = 0.0D0
      T(1,2)  = 0.0D0
      T(2,2)  = 0.0D0
        
      V(1,1,1)= 0.0D0
      V(1,2,1)= 0.0D0
      V(2,2,1)= 0.0D0
	  
      V(1,1,2)= 0.0D0
      V(1,2,2)= 0.0D0
      V(2,2,2)= 0.0D0

      Q(1,1,1,1) = 0.0D0 
      Q(2,1,1,1) = 0.0D0
      Q(2,1,2,1) = 0.0D0
      Q(2,2,1,1) = 0.0D0
      Q(2,2,2,1) = 0.0D0
      Q(2,2,2,2) = 0.0D0
     
      R2 = R*R
     
      DO 10 I = 1, N
       A(I,1) = EXPON(I,N)*(Z1**2)
       D(I,1) = COEFF(I,N)*((2.0D0*A(I,1)/PI)**0.75D0)
       A(I,2) = EXPON(I,N)*(Z2**2)
       D(I,2) = COEFF(I,N)*((2.0D0*A(I,2)/PI)**0.75D0)
10    CONTINUE
   
c      WRITE (6,*) A(1,1), A(1,2), A(2,1), A(2,2), A(3,1), A(3,2), PI
c      WRITE (6,*) D(1,1), D(1,2), D(2,1), D(2,2), D(3,1), D(3,2)  

CCCCC Calculate one electron integrals
	  
      DO 20 I=1,N
      DO 20 J=1,N
     
      RAP  = A(J,2) * R / (A(I,1) + A(J,2) )
      RAP2 = RAP**2
      RBP2 = (R-RAP)**2
      
      D11 = D(I,1) * D(J,1)
      D12 = D(I,1) * D(J,2)
      D22 = D(I,2) * D(J,2)
     
      S(1,2)  = S(1,2)   + S_RAW(A(I,1),A(J,2),R2   ) * D12  
     
      T(1,1)  = T(1,1)   + T_RAW(A(I,1),A(J,1),0.0D0) * D11
      T(1,2)  = T(1,2)   + T_RAW(A(I,1),A(J,2),R2   ) * D12
      T(2,2)  = T(2,2)   + T_RAW(A(I,2),A(J,2),0.0D0) * D22
        
      V(1,1,1)= V(1,1,1) + V_RAW(A(I,1),A(J,1),0.0D0, 0.0D0, ZA) * D11
      V(1,2,1)= V(1,2,1) + V_RAW(A(I,1),A(J,2),R2   , RAP2 , ZA) * D12 	  
      V(2,2,1)= V(2,2,1) + V_RAW(A(I,2),A(J,2),0.0D0, R2   , ZA) * D22 
	  
      V(1,1,2)= V(1,1,2) + V_RAW(A(I,1),A(J,1),0.0D0, R2,    ZB) * D11
      V(1,2,2)= V(1,2,2) + V_RAW(A(I,1),A(J,2),R2   , RBP2 , ZB) * D12 	  
      V(2,2,2)= V(2,2,2) + V_RAW(A(I,2),A(J,2),0.0D0, 0.0D0, ZB) * D22 
	  
20    CONTINUE     
     
      S(2,1) = S(1,2)
      T(2,1) = T(1,2)
      
      V(2,1,1) = V(1,2,1) 
      V(2,1,2) = V(1,2,2) 
     
CCCCC Calculate two electron integrals
	  
      DO 30 I=1,N
      DO 30 J=1,N     
      DO 30 K=1,N     
      DO 30 L=1,N     
	  
      RAP  = A(I,2) * R / (A(I,2) + A(J,1) )
      RAQ  = A(K,2) * R / (A(K,2) + A(L,1) )
      RBP  = R-RAP
      RBQ  = R-RAQ
      RPQ  = RAP - RAQ
	  
      RAP2  = RAP**2
      RAQ2  = RAQ**2
      RBP2  = RBP**2
      RBQ2  = RBQ**2
      RPQ2  = RPQ**2
	  
      Q(1,1,1,1) = Q(1,1,1,1) 
     $           +       D(I,1)*D(J,1)*D(K,1)*D(L,1) 
     $           * Q_RAW(A(I,1),A(J,1),A(K,1),A(L,1),0.0D0,0.0D0,0.0D0)
      Q(2,1,1,1) = Q(2,1,1,1) 
     $           +       D(I,2)*D(J,1)*D(K,1)*D(L,1) 
     $           * Q_RAW(A(I,2),A(J,1),A(K,1),A(L,1),R2   ,0.0D0,RAP2 )
      Q(2,1,2,1) = Q(2,1,2,1) 
     $           +       D(I,2)*D(J,1)*D(K,2)*D(L,1) 
     $           * Q_RAW(A(I,2),A(J,1),A(K,2),A(L,1),R2   ,R2   ,RPQ2 )
      Q(2,2,1,1) = Q(2,2,1,1) 
     $           +       D(I,2)*D(J,2)*D(K,1)*D(L,1) 
     $           * Q_RAW(A(I,2),A(J,2),A(K,1),A(L,1),0.0D0,0.0D0,R2   )
      Q(2,2,2,1) = Q(2,2,2,1) 
     $           +       D(I,2)*D(J,2)*D(K,2)*D(L,1) 
     $           * Q_RAW(A(I,2),A(J,2),A(K,2),A(L,1),0.0D0,R2,RBQ2 )
      Q(2,2,2,2) = Q(2,2,2,2) 
     $           +       D(I,2)*D(J,2)*D(K,2)*D(L,2) 
     $           * Q_RAW(A(I,2),A(J,2),A(K,2),A(L,2),0.0D0,0.0d0,0.0d0 )
	  
30    CONTINUE     
     
CCCCC Annoying Printing
  
      IF (IOP.NE.0) THEN
        PRINT 40 
        PRINT 50, R, Z1, Z2, S(1,2), T(1,1) 
        PRINT 60
        PRINT 50, T(1,2), T(2,2), V(1,1,1), V(1,2,1), V(2,2,1)
        PRINT 70
        PRINT 50, V(1,1,2), V(1,2,2), V(2,2,2)
        WRITE (6,*)
        WRITE (6,*) "   ","J(1,1,1,1), J(2,1,1,1), J(2,1,2,1), ",
     $                   "J(2,2,1,1), J(2,2,2,1), J(2,2,2,2)" 
        PRINT 80, Q(1,1,1,1), Q(2,1,1,1), Q(2,1,2,1), 
     $            Q(2,2,1,1), Q(2,2,2,1), Q(2,2,2,2)
      ENDIF 
     
CCCCC Build Coulomb Matrix

      Q(1,2,1,1) = Q(2,1,1,1)
      Q(1,1,2,1) = Q(2,1,1,1)
      Q(1,1,1,2) = Q(2,1,1,1)
      
      Q(1,2,1,2) = Q(2,1,2,1)
      Q(1,2,2,1) = Q(2,1,2,1)
      Q(2,1,1,2) = Q(2,1,2,1)
      
      Q(1,1,2,2) = Q(2,2,1,1)

      Q(1,2,2,2) = Q(2,2,2,1)
      Q(2,1,2,2) = Q(2,2,2,1)
      Q(2,2,1,2) = Q(2,2,2,1)      
   
      RETURN
     
CCCCC Annoying Fortran Print Codes

40    FORMAT(/,3X,1HR,10X,5HZETA1,6X,5HZETA2,6X, 3HS12, 8X, 3HT11)
50    FORMAT(5F11.6)
60    FORMAT(/,3X,3HT11,8X,3HT11,8X,4HV11A,7X,4HV12A,7X,4HV22A)
70    FORMAT(/,3X,4HV11B,7X,4HV12B,7X,4HV22B)
80    FORMAT(6F12.6)

      END
CCCCCC************************************************************
      FUNCTION F0(t)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /MATHS  / PI
     
      IF (t.LT.1.0D-6) THEN
        F0 = 1.0D0 - t / 3.0D0
      ELSE 
        F0 = 0.5D0 * DSQRT(PI/t)* ZERF(DSQRT(t))
c        write (6,*) "F0(",X,"): ", F0
      ENDIF

      RETURN
      END
CCCCCC************************************************************
      FUNCTION ZERF(t)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(5)
      DATA W/ 0.3275911D0   /
      DATA Y/ 0.254829592D0,
     $       -0.284496736D0,
     $        1.421413741D0,
     $       -1.453152027D0,
     $        1.061405429D0 /
     
      s    = 1.0D0 / (1.0D0 + W * t)
      snew = s
      POLY = Y(1) * snew
      
      DO 10 I = 2,5
        snew = snew * s
        POLY = POLY + Y(I) * snew
10    CONTINUE      
      
      ZERF = 1.0D0 - POLY * DEXP(-t**2)
      
c      write (6,*) 'p ', poly, zerf, X
      
      RETURN
      END
      
CCCCCC************************************************************
      FUNCTION S_RAW(A,B,RAB2)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /MATHS  / PI
	   
      S_RAW = (PI/(A+B))**1.5D0 * DEXP(-A*B*RAB2/(A+B))
      
      RETURN
      END

CCCCCC************************************************************
      FUNCTION T_RAW(A,B,RAB2)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /MATHS  / PI
      
      X = A*B/(A+B)
      
C      WRITE(6,*) A,B,X, PI
      
      T_RAW = X * (3.0D0 - 2.0D0*(X*RAB2))
     $      * ( PI / (A+B) )**1.5D0
     $      * DEXP(-X*RAB2)  
     
      RETURN
      END

CCCCCC************************************************************
      FUNCTION V_RAW(A,B,RAB2,RCP2,Z)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /MATHS  / PI
  
      V_RAW = - Z *2.0D0*PI*F0((A+B)*RCP2)* DEXP(-A*B*RAB2/(A+B))/(A+B)
      
      
      RETURN
      END
      
CCCCCC************************************************************
      FUNCTION Q_RAW(A,B,C,D,RAB2,RCD2,RPQ2)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /MATHS  / PI
      
      Q_RAW = 2.0D0 * (PI**2.5D0) / ((A+B)*(C+D)*DSQRT(A+B+C+D))
     $  * F0((A+B)*(C+D)*RPQ2/(A+B+C+D))  
     $  * DEXP(-A*B*RAB2/(A+B) - C*D*RCD2/(C+D))
     
      RETURN
      END
      
CCCCCC************************************************************
      SUBROUTINE COLLECT()
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /INTMAT   / S(2,2), T(2,2), V(2,2,2), Q(2,2,2,2) 
      COMMON /SCFMAT   / H(2,2), G(2,2), F(2,2)
      COMMON /TRANSMAT / X(2,2), P(2,2)
            
C Core Hamiltonian    
      H(1,1) = T(1,1) + V(1,1,1) + V(1,1,2)
      H(1,2) = T(1,2) + V(1,2,1) + V(1,2,2)
      H(2,1) = H(1,2)
      H(2,2) = T(2,2) + V(2,2,1) + V(2,2,2)
      
      CALL MATOUT(H,2,2,"H")
      
C S^-1/2      
      X(1,1) = 1.0D0 / DSQRT(2.0d0 * (1.0d0 + S(1,2))) 
      X(2,1) = X(1,1)
      X(1,2) = 1.0d0 / DSQRT(2.0d0 * (1.0d0 - S(1,2))) 
      X(2,2) = -X(1,2)

      CALL MATOUT(X,2,2,"X")
       
      RETURN 
      END
            
CCCCCC************************************************************
      SUBROUTINE SCF(IOP,N,R,Z1,Z2,ZA,ZB)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /INTMAT / S(2,2), T(2,2), V(2,2,2), Q(2,2,2,2) 
      COMMON /SCFMAT   / H(2,2), G(2,2), F(2,2)
      COMMON /TRANSMAT / X(2,2), P(2,2)
      COMMON /MATHS  / PI
      
      DIMENSION XT(2,2), FPRIME(2,2), CPRIME(2,2), C(2,2), E(2,2)
      DIMENSION POLD(2,2)
      
      CALL ZEROMATRIX(P,2,2)
      CALL TRANSPOSEMATRIX(X,XT,2,2)
      
      DATA TOL/1.0d-4/
      DATA MAXIT/25/
      
      
      ITER = 0
1000  ITER = ITER + 1
      IF (IOP.GT.0) THEN 
        WRITE (6,*) "Iteration Number: ", ITER
      ENDIF
       
      CALL FORMG
       
      EN=0.0D0
       
      DO 10 I=1,2
      DO 10 J=1,2        
10     EN = EN + 0.5d0 *P(I,J)*(H(I,J)+F(I,J))  
     
      IF (IOP.GE.2) THEN 
         CALL MATOUT(P,2,2,"P")
         CALL MATOUT(G,2,2,"G")
         CALL MATOUT(F,2,2,"F")
         PRINT 40, EN
      ENDIF
      
      CALL MULT(F,X,G,2)
      CALL MULT(XT,G,FPRIME,2)
      
      CALL DIAG(FPRIME,CPRIME,E)
      CALL MULT(X,CPRIME,C,2)
      

      CALL COPYMATRIX(P,POLD,2,2)
      DELTA = 0.0d0
      
      DO 21 I=1,2
      DO 21 J=1,2        
       P(I,J) = 0.0d0
      DO 20 K=1,1        
20     P(I,J) = P(I,J) + 2.0d0 *C(I,K)*C(J,K)
21     DELTA = DELTA + ( P(I,J) - POLD(I,J))**2

       DELTA = 0.50d0 * DSQRT(DELTA)
 
       IF (IOP.GE.2) THEN 
         CALL MATOUT(FPRIME,2,2,"F")
         CALL MATOUT(CPRIME,2,2,"CP")
         CALL MATOUT(E,2,2,"E")
         CALL MATOUT(C,2,2,"C")
         CALL MATOUT(P,2,2,"P")
         PRINT 40, EN
      ENDIF
 
       IF (ITER.LT.MAXIT) THEN
        IF (DELTA.GT.TOL) GOTO 1000
       ENDIF
       
       WRITE (6,*) "Converged Density Matrix"
              
       CALL MATOUT(P,2,2,"P")
       PRINT 40, EN
 
       EN = EN + ZA*ZB/R 
       
       PRINT 50, EN
       
c Calculate delta

      
      
30    FORMAT(/,28HStart of iteration number = , I2)      
40    FORMAT(///,20HElectronic Energy = ,D20.12)
50    FORMAT(///,20HTotal Energy      = ,D20.12)
      
      RETURN
      END      
      
CCCCCC************************************************************
      SUBROUTINE FORMG
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /INTMAT / S(2,2), T(2,2), V(2,2,2), Q(2,2,2,2) 
      COMMON /SCFMAT   / H(2,2), G(2,2), F(2,2)
      COMMON /TRANSMAT / X(2,2), P(2,2)
       
      DO 10 I=1,2 
      DO 10 J=1,2
       G(I,J) = 0d0      
      DO 20 K=1,2 
      DO 20 L=1,2
20     G(I,J) = G(I,J) + P(K,L)*(Q(I,J,K,L)-0.5D0*Q(I,L,K,J))
10     F(I,J) = H(I,J) + G(I,J)

      RETURN
      END      

CCCCCC************************************************************
      SUBROUTINE DIAG(A,C,E)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(2,2),C(2,2), E(2,2)
      COMMON /MATHS  / PI
      
      A1122 = A(1,1) - A(2,2)
      
      IF (DABS(A1122).LT.1.0D-20) THEN
        PHI = 0.25D0 * PI 
      ELSE 
        PHI = 0.5D0 * DATAN(2.0d0 * A(1,2) / A1122)
      ENDIF      
      
      C(1,1) = DCOS(PHI)
      C(1,2) = DSIN(PHI)
      C(2,1) = DSIN(PHI)
      C(2,2) = - DCOS(PHI)
      
      
      E(1,1) = A(1,1)*DCOS(PHI)**2
     $       + A(2,2)*DSIN(PHI)**2
     $       + A(1,2)*DSIN(2.0d0*PHI)
      
      E(2,2) = A(1,1)*DSIN(PHI)**2
     $       + A(2,2)*DCOS(PHI)**2
     $       - A(1,2)*DSIN(2.0d0*PHI)
	      
      E(1,2) = 0.0d0
      E(2,1) = 0.0d0	
				
      IF (E(1,1).gt.E(2,2)) THEN
       TMP = E(2,2)
       E(2,2) = E(1,1)
       E(1,1) = TMP
		 
       TMP = C(1,2)
       C(1,2) = C(1,1)
       C(1,1) = TMP
       
    	 TMP = C(2,2)
       C(2,2) = C(2,1)
       C(2,1) = TMP
      ENDIF
      
      RETURN
      END

      
CCCCCC************************************************************
      SUBROUTINE MULT(A,B,C,M)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,M), B(M,M), C(M,M)
      
      DO 10 I=1,M
      DO 10 J=1,M 
       C(I,J) = 0.0D0
      DO 10 K=1,M
10     C(I,J) = C(I,J) + A(I,K) * B(K,J)
      
      RETURN
      END
      
CCCCCC************************************************************
      SUBROUTINE MATOUT(A,M,N,LABEL)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N)
      CHARACTER LABEL
      
      WRITE (6,*) LABEL
      
      DO 40 I=1,M 
       PRINT 20,(A(I,J),J=1,N)
40    CONTINUE
 
20    FORMAT (1X,10D18.10)
30    FORMAT (//)
    
      RETURN
      END
      
CCCCCC************************************************************
      SUBROUTINE COPYMATRIX(A,ANEW,N,M)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
     
      DIMENSION A(N,M), ANEW(N,M)
     
      DO 10 I=1,N
      DO 10 J=1,N
 10     ANEW(I,J) = A(I,J)
 
      RETURN
      END
      
CCCCCC************************************************************
      SUBROUTINE TRANSPOSEMATRIX(A,ANEW,N,M)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
     
      DIMENSION A(N,M), ANEW(M,N)
     
      DO 10 I=1,N
      DO 10 J=1,M
 10     ANEW(J,I) = A(I,J)
 
      RETURN
      END
      
      
CCCCCC************************************************************
      SUBROUTINE ZEROMATRIX(A,N,M)
CCCCCC************************************************************
     
CCCCCC************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
     
      DIMENSION A(N,M)
     
      DO 10 I=1,N
      DO 10 J=1,M
 10     A(I,J) = 0.0d0
 
      RETURN
      END
      
