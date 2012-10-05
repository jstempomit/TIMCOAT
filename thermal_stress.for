C***********************************************************************
C                                                                      *
C  Subroutine TSTRESS_3D(DT, SIGR, SIGT)                               *
C                                                                      *
C                                                                      *
C                                                                      *
C  Discription: Use IMSL math libraries to calculate thermal stress    *
C               distributions in TRISO partilces.                      *
C                                                                      *
C                                                                      *
C  Assumption:  Spherical symmetry applies, and the analysis is 3D.    *
C               The pyrocarbon layers are anisotropic and SiC layer    *
C               is isotropic.							                 *
C                                                                      *
C                                                                      *
C  Actual argument description                                         *
C                                                                      *
C      DT:  Temperature load (C)                                       *
C                                                                      *
C  Description of returned quantities                                  *
C                                                                      *
C  SIGR(0:NDIV):  Radial stress distribution over three layers         *
C  SIGT(0:NDIV):  Tangential stress distribution over three layers     *
C                                                                      *
C  Dummy argument description                                          *
C                                                                      *
C                                                                      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
	SUBROUTINE TSTRESS_3D(DT, SIGR, SIGT)
C  Number of divisions for stress distribution
	PARAMETER ( NDIV = 30)
	PARAMETER ( IPATH = 1, LDA = 6, NDIM = 6 )
	DOUBLE PRECISION DT, SIGR, SIGT
	DOUBLE PRECISION P(LDA,LDA), Q(NDIM), X(NDIM) !arrays for linear equations
	DOUBLE PRECISION R1, R2, R3, R4, R5
	DOUBLE PRECISION IPYCE1,IPYCE3,IPYCNU13,IPYCNU12,IPYCNU31,
     &                 OPYCE1,OPYCE3,OPYCNU13,OPYCNU12,OPYCNU31,
     &                 SICE, SICNU, SICALPHA,
     &                 IPYCALPHA1,IPYCALPHA3,OPYCALPHA1,OPYCALPHA3
	DOUBLE PRECISION IPYCE, IPYCNU, OPYCE, OPYCNU, 
     &				 IPYCALPHA, OPYCALPHA
      DOUBLE PRECISION M1I, M2I, M1O, M2O, KPI, KPO, HI, HO, 
     &                 GI1, GI2, GI3, GI4, GS1, GS2, GO1, GO2, GO3, GO4
	DOUBLE PRECISION RANGE, R
	DOUBLE PRECISION A, B, C, D, E, F             !coefficients for stresses
	INTEGER I
	DIMENSION SIGR(0:NDIV), SIGT(0:NDIV)
C
      COMMON /RADII/ R1, R2, R3, R4, R5   !particle geometry
      COMMON /LAYERS/ IPYCE1,IPYCE3,IPYCNU13,IPYCNU12,IPYCNU31,
     &                OPYCE1,OPYCE3,OPYCNU13,OPYCNU12,OPYCNU31,
     &                SICE, SICNU, SICALPHA,
     &                IPYCALPHA1,IPYCALPHA3,OPYCALPHA1,OPYCALPHA3
	DATA X0/0.0 D+00/, X1/1.0 D+00/, X2/2.0 D+00/, X3/3.0 D+00/
C
C  Prepare quantities based on material properties
C
      IF((IPYCE1.NE.IPYCE3).AND.(OPYCE1.NE.OPYCE3)) THEN
	  M1I = (-X1 + SQRT(X1+8.0D0*IPYCNU13*(X1-IPYCNU31)/IPYCNU31/
     &		(X1-IPYCNU12)))/X2
	  M2I = (-X1 - SQRT(X1+8.0D0*IPYCNU13*(X1-IPYCNU31)/IPYCNU31/
     &		(X1-IPYCNU12)))/X2
	  M1O = (-X1 + SQRT(X1+8.0D0*OPYCNU13*(X1-OPYCNU31)/OPYCNU31/
     &		(X1-OPYCNU12)))/X2
	  M2O = (-X1 - SQRT(X1+8.0D0*OPYCNU13*(X1-OPYCNU31)/OPYCNU31/
     &		(X1-OPYCNU12)))/X2
        KPI = (IPYCNU31*(X1-IPYCNU12-IPYCNU13)*IPYCALPHA3-IPYCNU13*
     &		(X1-X2*IPYCNU31)*IPYCALPHA1)*DT/(IPYCNU31*(X1-IPYCNU12)
     &		 -IPYCNU13*(X1-IPYCNU31))
        KPO = (OPYCNU31*(X1-OPYCNU12-OPYCNU13)*OPYCALPHA3-OPYCNU13*
     &		(X1-X2*OPYCNU31)*OPYCALPHA1)*DT/(OPYCNU31*(X1-OPYCNU12)
     &		 -OPYCNU13*(X1-OPYCNU31))
        HI = IPYCNU13*IPYCE3*(IPYCALPHA3-IPYCALPHA1)*DT/
     &       (IPYCNU31*(X1-IPYCNU12)-IPYCNU13*(X1-IPYCNU31))
        HO = OPYCNU13*OPYCE3*(OPYCALPHA3-OPYCALPHA1)*DT/
     &       (OPYCNU31*(X1-OPYCNU12)-OPYCNU13*(X1-OPYCNU31))
        GI1 = ((X1-IPYCNU12)*M1I+X2*IPYCNU13)*IPYCE3/
     &		(X1-IPYCNU12-X2*IPYCNU13*IPYCNU31)
	  GI2 = ((X1-IPYCNU12)*M2I+X2*IPYCNU13)*IPYCE3/
     &		(X1-IPYCNU12-X2*IPYCNU13*IPYCNU31)
        GI3 = (X1+M1I)*GI1/X2
        GI4 = (X1+M2I)*GI2/X2
        GS1 = SICE/(X1 - X2*SICNU)
	  GS2 = SICE/(X1 + SICNU)
        GO1 = ((X1-OPYCNU12)*M1O+X2*OPYCNU13)*OPYCE3/
     &		(X1-OPYCNU12-X2*OPYCNU13*OPYCNU31)
	  GO2 = ((X1-OPYCNU12)*M2O+X2*OPYCNU13)*OPYCE3/
     &		(X1-OPYCNU12-X2*OPYCNU13*OPYCNU31)
        GO3 = (X1+M1O)*GO1/X2
        GO4 = (X1+M2O)*GO2/X2
	  RANGE = R5 - R2
C
C    Step 1: Assign values to P(LDA,LDA)
        P(1,1) = GI1*(R2**(M1I-X1))
	  P(1,2) = GI2*(R2**(M2I-X1))
	  P(2,1) = GI1*(R3**(M1I-X1))
	  P(2,2) = GI2*(R3**(M2I-X1))
        P(2,3) = -GS1
	  P(2,4) = X2*GS2*(R3**(-X3))
	  P(3,1) = R3**(M1I-X1)
	  P(3,2) = R3**(M2I-X1)
	  P(3,3) = -X1
	  P(3,4) = -R3**(-X3)
	  P(4,3) = -X1
	  P(4,4) = -R4**(-X3)
	  P(4,5) = R4**(M1O-X1)
	  P(4,6) = R4**(M2O-X1)
	  P(5,3) = -GS1
	  P(5,4) = X2*GS2*(R4**(-X3))
	  P(5,5) = GO1*(R4**(M1O-X1))
	  P(5,6) = GO2*(R4**(M2O-X1))
	  P(6,5) = GO1*(R5**(M1O-X1))
	  P(6,6) = GO2*(R5**(M2O-X1))
C    Step 2: Assign values to Q(NDIM)
        Q(1) = -HI
        Q(2) = -HI - GS1*SICALPHA*DT
	  Q(3) = -KPI
        Q(4) = -KPO
	  Q(5) = -HO - GS1*SICALPHA*DT
	  Q(6) = -HO
C    Step 3: Use DLSARG of IMSL library to calculate X
        CALL DLSARG (NDIM, P, LDA, Q, IPATH, X)
        A = X(1)
	  B = X(2)
	  C = X(3)
	  D = X(4)
	  E = X(5)
	  F = X(6)
C    Step 4: Calculate thermal stresses based on A to F
C
	  DO 600 I=0, NDIV
	    R = I*RANGE/FLOAT(NDIV) + R2
	    IF (R .LE. R3) THEN
	      SIGR(I) = HI + GI1*A*(R**(M1I-X1)) + GI2*B*(R**(M2I-X1))
	      SIGT(I) = HI + GI3*A*(R**(M1I-X1)) + GI4*B*(R**(M2I-X1))
	    ELSE IF (R .LT. R4) THEN
	      SIGR(I) = GS1*(C-SICALPHA*DT) - X2*GS2*D*(R**(-X3))
	      SIGT(I) = GS1*(C-SICALPHA*DT) + GS2*D*(R**(-X3))
	    ELSE
	      SIGR(I) = HO + GO1*E*(R**(M1O-X1)) + GO2*F*(R**(M2O-X1))
	      SIGT(I) = HO + GO3*E*(R**(M1O-X1)) + GO4*F*(R**(M2O-X1))
	    END IF
600     CONTINUE
      ELSE IF((IPYCE1.EQ.IPYCE3).AND.(OPYCE1.EQ.OPYCE3)) THEN
C    When PyC layers are isotropic, the formulations are different from above
	  IPYCE = IPYCE1
	  IPYCNU = IPYCNU13
	  OPYCE = OPYCE1
	  OPYCNU = OPYCNU13
	  IPYCALPHA = IPYCALPHA1
	  OPYCALPHA = OPYCALPHA1
C
        GI1 = IPYCE/(X1 - X2*IPYCNU)
	  GI2 = IPYCE/(X1 + IPYCNU)
	  GS1 = SICE/(X1 - X2*SICNU)
	  GS2 = SICE/(X1 + SICNU)
        GO1 = OPYCE/(X1 - X2*OPYCNU)
	  GO2 = OPYCE/(X1 + OPYCNU)
	  RANGE = R5 - R2
C
C    Step 1: Assign values to P(LDA,LDA)
        P(1,1) = GI1
	  P(1,2) = -X2*GI2*(R2**(-X3))
	  P(2,1) = GI1
	  P(2,2) = -X2*GI2*(R3**(-X3))
        P(2,3) = -GS1
	  P(2,4) = X2*GS2*(R3**(-X3))
	  P(3,1) = X1
	  P(3,2) = R3**(-X3)
	  P(3,3) = -X1
	  P(3,4) = -R3**(-X3)
	  P(4,3) = -X1
	  P(4,4) = -R4**(-X3)
	  P(4,5) = X1
	  P(4,6) = R4**(-X3)
	  P(5,3) = -GS1
	  P(5,4) = X2*GS2*(R4**(-X3))
	  P(5,5) = GO1
	  P(5,6) = -X2*GO2*(R4**(-X3))
	  P(6,5) = GO1
	  P(6,6) = -X2*GO2*(R5**(-X3))
C    Step 2: Assign values to Q(NDIM)
        Q(1) = GI1*IPYCALPHA*DT
        Q(2) = (GI1*IPYCALPHA - GS1*SICALPHA)*DT
	  Q(3) = 0
        Q(4) = 0
	  Q(5) = (GO1*OPYCALPHA - GS1*SICALPHA)*DT
	  Q(6) = GO1*OPYCALPHA*DT
C    Step 3: Use DLSARG of IMSL library to calculate X
        CALL DLSARG (NDIM, P, LDA, Q, IPATH, X)
        A = X(1)
	  B = X(2)
	  C = X(3)
	  D = X(4)
	  E = X(5)
	  F = X(6)
C    Step 4: Calculate thermal stresses based on A to F
C
	  DO 610 I=0, NDIV
	    R = I*RANGE/FLOAT(NDIV) + R2
	    IF (R .LE. R3) THEN
	      SIGR(I) = GI1*(A-IPYCALPHA*DT) - X2*GI2*B*(R**(-X3))
	      SIGT(I) = GI1*(A-IPYCALPHA*DT) + GI2*B*(R**(-X3))
	    ELSE IF (R .LT. R4) THEN
	      SIGR(I) = GS1*(C-SICALPHA*DT) - X2*GS2*D*(R**(-X3))
	      SIGT(I) = GS1*(C-SICALPHA*DT) + GS2*D*(R**(-X3))
	    ELSE
	      SIGR(I) = GO1*(E-OPYCALPHA*DT) - X2*GO2*F*(R**(-X3))
	      SIGT(I) = GO1*(E-OPYCALPHA*DT) + GO2*F*(R**(-X3))
	    END IF
610     CONTINUE
	END IF
      RETURN
	END
C                                                                      *
C***********************************************************************
C
C