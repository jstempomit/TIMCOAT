
C***********************************************************************
C                                                                      *
C  Subroutine M_FUEL                                                   *
C                                                                      *
C***********************************************************************
C
C  Function: Calculate stress distribution at certain fluence based on 
C            FUEL code.
C
C
C***********************************************************************
C  Description of entrance parameters of subroutine
C
C  PRESS      :  Internal gas pressure
C  PAMB       :  External ambient pressure
C  PYCE       :  Young's modulus of PyC (PYCE=2.0D+04+(5.0D+03)*FLU(STEP))
C  PYCNU      :  Poisson's ratio (elastic) for the pyrocarbons
C  PYCREEP    :  Creep coefficient for the pyrocarbons
C  SICE       :  Young's modulus for the SiC
C  SICNU      :  Poisson's ratio for the SiC
C  R2         :  Actual buffer outer radius
C  R3         :  Actual IPyC outer radius
C  R4         :  Actual SiC outer radius
C  R5         :  Actual OPyC outer radius
C  FLU(STEP)  :  Fluence at time of calculation
C  IPYCD      :  Density of the IPyC
C  OPYCD      :  Density of the OPyC
C
C  Description of returned quantities
C
C  Description of local variables
C
C
C***********************************************************************
C
	SUBROUTINE M_FUEL(PRESS, PAMB, FLU, SIGR, SIGT)
C  Number of intervals for stress distribution
	PARAMETER ( NDIV = 9)
	DOUBLE PRECISION PRESS, PAMB, FLU
	DOUBLE PRECISION R1, R2, R3, R4, R5
	DOUBLE PRECISION IPYCALPHA, IPYCCNU, IPYCE, IPYCNU, IPYCREEP,
     &                 OPYCALPHA, OPYCCNU, OPYCE, OPYCNU, OPYCREEP,
     &                 SICALPHA, SICE, SICNU, IPYCD, OPYCD
      DOUBLE PRECISION ISWR, ISWT, OSWR, OSWT
	DOUBLE PRECISION SIGR, SIGT
	DOUBLE PRECISION SIGRI, SIGRO, SIGTO, FI, FO
	DOUBLE PRECISION R
      INTEGER NDIVI, NDIVS, NDIVO
	INTEGER I
	DIMENSION SIGR(1:NDIV), SIGT(1:NDIV)
      DIMENSION OSWR(0:3), OSWT(0:3), ISWR(0:3), ISWT(0:3)
C
      COMMON /PAR_R/ R1, R2, R3, R4, R5
	COMMON /PAR_DIV/ NDIVI, NDIVS, NDIVO
      COMMON /PAR_M/ IPYCALPHA, IPYCCNU, IPYCE, IPYCNU, IPYCREEP,
     &               OPYCALPHA, OPYCCNU, OPYCE, OPYCNU, OPYCREEP,
     &               SICALPHA, SICE, SICNU, IPYCD, OPYCD,
     &               ISWR, ISWT, OSWR, OSWT
C
	DATA X0/0.0 D+00/, X1/1.0 D+00/, X2/2.0 D+00/, X3/3.0 D+00/
C
C FLU cannot be zero in INTRFACE.
	IF(FLU .EQ. 0.0D0) RETURN
C  Calculate stresses at two inferfaces
C  Here the Subroutine INTRFACE treats IPyC and OPyC as the same material
      CALL INTRFACE(.FALSE.,.FALSE., PRESS, PAMB, FLU,
     B              SIGRI, SIGRO, SIGTO, FI, FO)
	DO 902 I = 1, NDIVI+2
	  R = (I-1)*(R3-R2)/FLOAT(NDIVI+1) + R2
        SIGR(I) = (R3**(-3)-R**(-3))*(-PRESS)/(R3**(-3)-R2**(-3)) +
     A             (R**(-3)-R2**(-3))*SIGRI/(R3**(-3)-R2**(-3)) +
     B             X2*((R2**(-3)-R**(-3))*DLOG(R3)/(R3**(-3)-R2**(-3))
     C             + (R**(-3)-R3**(-3))*DLOG(R2)/(R3**(-3)-R2**(-3))
     D             + DLOG(R))*FI/X3
	  SIGT(I)=(R**(-3)+X2*R3**(-3))*(-PRESS)/(X2*(R3**(-3)-R2**(-3)))
     A          -(R**(-3)+X2*R2**(-3))*SIGRI/(X2*(R3**(-3)-R2**(-3)))
     B      +((R**(-3)+X2*R2**(-3))*DLOG(R3)/(X3*(R3**(-3)-R2**(-3)))
     C      - (R**(-3)+X2*R3**(-3))*DLOG(R2)/(X3*(R3**(-3)-R2**(-3)))
     D      + X2*DLOG(R)/X3 + X1/X3)*FI
902	CONTINUE
	DO 903 I = 1, NDIVS+2
	  R = (I-1)*(R4-R3)/FLOAT(NDIVS+1) + R3
	  SIGR(NDIVI+2+I) = (R4**(-3)-R**(-3))*SIGRI/(R4**(-3)-R3**(-3))+
     A            (R**(-3)-R3**(-3))*SIGRO/(R4**(-3)-R3**(-3))
	  SIGT(NDIVI+2+I) = (R**(-3)+X2*R4**(-3))*SIGRI/
     A			(X2*(R4**(-3)-R3**(-3)))-(R**(-3)+X2*R3**(-3))
     B            *SIGRO/(X2*(R4**(-3)-R3**(-3)))
903	CONTINUE
	DO 904 I = 1, NDIVO+2
	  R = (I-1)*(R5-R4)/FLOAT(NDIVO+1) + R4
	  SIGR(NDIVI+NDIVS+4+I) = (R4**(-3)-R**(-3))*(-PAMB)/
     A			(R4**(-3)-R5**(-3)) + (R**(-3)-R5**(-3))*SIGRO
     B            /(R4**(-3)-R5**(-3)) + X2*((R5**(-3)-R**(-3))
     C            *DLOG(R4)/(R4**(-3)-R5**(-3)) + (R**(-3)-R4**(-3))
     D            *DLOG(R5)/(R4**(-3)-R5**(-3)) + DLOG(R))*FO/X3
	  SIGT(NDIVI+NDIVS+4+I) = (R**(-3)+X2*R4**(-3))*(-PAMB)/
     A			(X2*(R4**(-3)-R5**(-3)))-(R**(-3)+X2*R5**(-3))
     B            *SIGRO/(X2*(R4**(-3)-R5**(-3)))
     C       +((R**(-3)+X2*R5**(-3))*DLOG(R4)/(X3*(R4**(-3)-R5**(-3)))
     D       - (R**(-3)+X2*R4**(-3))*DLOG(R5)/(X3*(R4**(-3)-R5**(-3)))
     E       + X2*DLOG(R)/X3 + X1/X3)*FO
904   CONTINUE
      RETURN
	END
C
C
C***********************************************************************
C                                                                      *
C  Subroutine INTRFACE                                                 *
C                                                                      *
C***********************************************************************
C
C  Calculate the stresses at two interfaces, IPyC/SiC and OPyC/SiC
C	
C  Assumption: IPyC and OPyC have to be the same material
C
C***********************************************************************
C
C  Description of input variables
C 
C  DEBONDI    :  Logical variable with the following outcomes
C                .TRUE.  IPyC layer debonds
C                .FALSE. No debonding of IPyC
C  DEBONDO    :  Logical variable with the following outcomes
C             :  .TRUE.  OPyC layer fails or debonds
C             :  .FALSE. No debonding or failure of OPyC
C
C  Note: Subroutine does not currently allow for debonding both PyC
C        layers in the same calculation.
C
C  PRESS      :  Internal gas pressure
C  PAMB       :  External pressure
C  PYCE       :  Young's modulus for the pyrocarbons
C  PYCNU      :  Poisson's ratio (elastic) for the pyrocarbons
C  R2         :  Actual buffer outer radius
C  R3         :  Actual IPyC outer radius
C  R4         :  Actual SiC outer radius
C  R5         :  Actual OPyC outer radius
C  SICE       :  Young's modulus for the SiC
C  SICNU      :  Poisson's ratio for the SiC
C  PYCREEP    :  Creep coefficient for the pyrocarbons
C  FLU        :  Fluence at time of calculation
C  IPYCD      :  Density of the IPyC
C  OPYCD      :  Density of the OPyC
C 
C  Description of calculated quantities
C
C  SIGRI      :  Radial stress at the IPyC - SiC interface due to creep,
C                 swelling, and pressure
C  SIGRO      :  Radial stress at the SiC - OPyC interface due to creep,
C                 swelling, and pressure
C  SIGTO      :  Tangential stress at the inner surface of the OPyC due
C                 to creep, swelling, and pressure
C
C  Description of local variables
C
C  ISWR0      :  Zeroth order IPyC radial swelling coefficients
C  ISWT0      :  Zeroth order IPyC tangential swelling coefficients
C  OSWR0      :  Zeroth order OPyC radial swelling coefficients
C  OSWT0      :  Zeroth order OPyC tangential swelling coefficients
C  ISWR       :  Higher order IPyC radial swelling coefficients
C  ISWT       :  Higher order IPyC tangential swelling coefficients
C  OSWR       :  Higher order OPyC radial swelling coefficients
C  OSWT       :  Higher order OPyC tangential swelling coefficients
C  A          :
C  A0         :
C  A1         :
C  A2         :
C  A3         :
C  A23        :  Coefficients associated with SiC
C  A24        :  Coefficients associated with SiC
C  A54        :  Coefficients associated with SiC
C  B          :
C  B11        :
C  B12        :
C  B21        :
C  B22        :
C  C          :  Creep coefficient for the pyrocarbons
C  C1         :
C  C2         :
C  Q1         :
C  Q2         :
C  H0         :
C  H1         :
C  H2         :
C  H3         :
C  G0         :
C  G1         :
C  G2         :
C  G3         :
C  C1         :  ! Why are C1 and C2 listed twice?
C  C2         :
C  F          :
C  FI         :
C  FO         :
C  K          :
C  K0         :
C  K1         :
C  K2         :
C  K3         :
C  L          :
C  L0         :
C  L1         :
C  L2         :
C  L3         :
C  M          :
C  M1         :
C  M2         :
C  N1         :
C  N2         :
C  N3         :
C  N4         :
C  N5         :
C  N6         :
C  N7         :
C  N8         :
C  M1         :  ! Why are M1 and M2 listed twice?
C  M2         :
C  RND        :  1/(B11*B22 - B21*B12)
C  R45L       :  DLOG(R4/R5)
C  R32L       :  DLOG(R3/R2)
C  R23        :  R2**3
C  R33        :  R3**3
C  R43        :  R4**3
C  R53        :  R5**3
C  X0         :  Numerical constant 0
C  X1         :  Numerical constant 1
C  X2         :  Numerical constant 2
C  X3         :  Numerical constant 3
C  X4         :  Numerical constant 4
C  X6         :  Numerical constant 6
C  X8         :  Numerical constant 8
C  X9         :  Numerical constant 9
C  X16        :  Numerical constant 16
C  Z          :
C  R23        :  R2**3
C  R33        :  R3**3
C  R43        :  R4**3
C  R53        :  R5**3
C
C***********************************************************************
C
C  Subroutines and functions called:
C
C  None
C
C  Intrinsic functions:
C
C  DEXP
C  DLOG
C  DSQRT
C
C***********************************************************************
C
      SUBROUTINE INTRFACE (DEBONDI, DEBONDO, PRESS, PAMB, FLU,
     &                 SIGRI, SIGRO, SIGTO, FI, FO)
      DOUBLE PRECISION PRESS, PAMB, FLU, SIGRI, SIGRO, SIGTO, FI, FO
	DOUBLE PRECISION R1, R2, R3, R4, R5
	DOUBLE PRECISION IPYCALPHA, IPYCCNU, IPYCE, IPYCNU, IPYCREEP,
     &                 OPYCALPHA, OPYCCNU, OPYCE, OPYCNU, OPYCREEP,
     &                 SICALPHA, SICE, SICNU, IPYCD, OPYCD
	DOUBLE PRECISION PYCE, PYCNU, PYCREEP
      DOUBLE PRECISION R23, R33, R43, R53
      DOUBLE PRECISION ISWR, ISWT, OSWR, OSWT
      DOUBLE PRECISION C
      DOUBLE PRECISION A, A0, A1, A2, A3, A23, A24, A54
      DOUBLE PRECISION B, B11, B12, B21, B22
      DOUBLE PRECISION Q1, Q2, H0, H1, H2, H3, G0, G1, G2, G3, C1, C2
      DOUBLE PRECISION F, Z
      DOUBLE PRECISION R45L, R32L, RND
      DOUBLE PRECISION X0, X1, X2, X3, X4, X6, X8, X9, X16
      LOGICAL DEBONDI,DEBONDO
C
C  Unconventional declarations
C
      DOUBLE PRECISION K, K0, K1, K2, K3, L, L0, L1, L2, L3, H
      DOUBLE PRECISION N1, N2, N3, N4, N5, N6, N7, N8, M, M1, M2

      DIMENSION OSWR(0:3), OSWT(0:3), ISWR(0:3), ISWT(0:3)

C
      COMMON /PAR_R/ R1, R2, R3, R4, R5
      COMMON /PAR_M/ IPYCALPHA, IPYCCNU, IPYCE, IPYCNU, IPYCREEP,
     &               OPYCALPHA, OPYCCNU, OPYCE, OPYCNU, OPYCREEP,
     &               SICALPHA, SICE, SICNU, IPYCD, OPYCD,
     &               ISWR, ISWT, OSWR, OSWT
C
      DATA X0/0.0 D+00/, X1/1.0 D+00/, X2/2.0 D+00/, X3/3.0 D+00/
      DATA X4/4.0 D+00/, X6/6.0 D+00/, X8/8.0 D+00/, X9/9.0 D+00/
      DATA X16/16.0 D+00/
C
	PYCE = IPYCE
	PYCNU = IPYCNU
	PYCREEP = IPYCREEP    
C  Determine the following four quantities based on the material
C  properties of the pyrocarbon layers
C
      L = X1/PYCE
      H = X2*PYCNU/PYCE
      M = PYCNU/PYCE
      K = (X1 - PYCNU)/PYCE
C
C  Following are the creep and swelling coefficients.  Swelling
C  coefficients are given for third order polynomials.
C
      C = PYCREEP
C
C  Radii**3                                        
C
      R23 = R2**3
      R33 = R3**3
      R43 = R4**3
      R53 = R5**3
C
C  Following are matrix coefficients needed to satisfy boundary
C  conditions
C
      A23 = -SICE/(X1 - X2*SICNU)
      A24 = X2*SICE/(R33*(X1 + SICNU))
      A54 = X2*SICE/(R43*(X1 + SICNU))
C
      SIGRI = X0
      SIGRO = X0
      SIGTO = X0
C
      A = (X2*(K - M)*R33/R23 + L + M)/X2
      B = (X2*(K - M)*R43/R53 + L + M)/X2
      B11 = - (X1 - R33/R23)/A + A23
      B12 = - (X1/R33 - X1/R23)/A + A24
      B21 = - (X1 - R43/R53)/B + A23
      B22 = - (X1/R43 - X1/R53)/B + A54
      N5 = (A23*A54 - A24*A23)/(- B11*A54 + A23*B12)
      N6 = (A23*B12 - A24*B11)/(- B11*A54 + A23*B12)
      N7 = (A23*A24 - A54*A23)/(- B21*A24 + A23*B22)
      N8 = (A23*B22 - A54*B21)/(- B21*A24 + A23*B22)
C
C  The following intrinsic functions will reduce execution time
C
      R45L = DLOG(R4/R5)
      R32L = DLOG(R3/R2)
C
C  If IPyC layer debonds
C
      IF (DEBONDI .EQV. .TRUE.) THEN
        A3 = (X2/C)*(OSWR(3) - OSWT(3))
        A2 = (X2/C)*((OSWR(2) - OSWT(2)) - X3*K*A3)
        A1 = (X2/C)*((OSWR(1) - OSWT(1)) - X2*K*A2)
        A0 = (X2/C)*((OSWR(0) - OSWT(0)) - K*A1)
        F = ((A3*FLU + A2)*FLU + A1)*FLU + A0 -
     A        A0*DEXP(- C*FLU/(X2*K))
        G0 = - (N7*R45L/B)*(OSWR(0) - OSWT(0)) -
     A         (N7*(X1 - R43/R53)/(X3*B))*(OSWR(0) + X2*OSWT(0)) +
     B         (X3*C*N7*PAMB)/(X4*B) - N8*PRESS/FLU
        G1 = - (N7*R45L/B)*(OSWR(1) - OSWT(1)) -
     A         (N7*(X1 - R43/R53)/(X3*B))*(OSWR(1) + X2*OSWT(1))
        G2 = - (N7*R45L/B)*(OSWR(2) - OSWT(2)) -
     A         (N7*(X1 - R43/R53)/(X3*B))*(OSWR(2) + X2*OSWT(2))
        G3 = - (N7*R45L/B)*(OSWR(3) - OSWT(3)) -
     A         (N7*(X1 - R43/R53)/(X3*B))*(OSWR(3) + X2*OSWT(3))
        M1 = X3*C*N7/(X4*B)
        L3 = - X4*B*G3/(X3*C*N7)
        L2 = - X4*B*(G2 - X3*L3)/(X3*C*N7)
        L1 = - X4*B*(G1 - X2*L2)/(X3*C*N7)
        L0 = - X4*B*(G0 - L1)/(X3*C*N7)
        C1 = X3*PAMB*K*N7/(X2*B) - L0
C
C  Radial and tangential stress in SiC and OPyC
C
        SIGRI = - PRESS
        SIGRO = ((L3*FLU + L2)*FLU + L1)*FLU + L0 + C1*DEXP(M1*FLU)
        SIGTO = SIGRO + (X3/(X2*(R43/R53 - X1)))*
     A         (PAMB + SIGRO - X2*F*R45L/X3) + F/X3
C
C  If IPyC layer is intact
C
      ELSE IF (DEBONDI .EQV. .FALSE.) THEN
        H0 = (X3*K*N5/(X2*A))*PRESS/FLU -
     A       (N5*R32L/A)*(ISWR(0) - ISWT(0)) -
     B       (N5*(X1 - R33/R23)/(X3*A))*(ISWR(0) + X2*ISWT(0))
        H1 = - (N5*R32L/A)*(ISWR(1) - ISWT(1)) -
     A         (N5*(X1 - R33/R23)/(X3*A))*(ISWR(1) + X2*ISWT(1)) +
     B         (X3*C*N5*PRESS/FLU)/(X4*A)
        H2 = - (N5*R32L/A)*(ISWR(2) - ISWT(2)) -
     A         (N5*(X1 - R33/R23)/(X3*A))*(ISWR(2) + X2*ISWT(2))
        H3 = - (N5*R32L/A)*(ISWR(3) - ISWT(3)) -
     A         (N5*(X1 - R33/R23)/(X3*A))*(ISWR(3) + X2*ISWT(3))
C
C  If OPyC layer debonds
C
        IF (DEBONDO .EQV. .TRUE.) THEN
          M1 = X3*C*N5/(X4*A)
          K3 = - X4*A*H3/(X3*C*N5)
          K2 = - X4*A*(H2 - X3*K3)/(X3*C*N5)
          K1 = - X4*A*(H1 - X2*K2)/(X3*C*N5)
          K0 = - X4*A*(H0 - K1)/(X3*C*N5)
          C1 = - PAMB*N6 - K0
          SIGRI = ((K3*FLU + K2)*FLU + K1)*FLU + C1*DEXP(M1*FLU) + K0
          SIGRO = - PAMB
          SIGTO = X0
C
C  If OPyC layer is intact
C
        ELSE IF (DEBONDO .EQV. .FALSE.) THEN
          A3 = (X2/C)*(OSWR(3) - OSWT(3))
          A2 = (X2/C)*((OSWR(2) - OSWT(2)) - X3*K*A3)
          A1 = (X2/C)*((OSWR(1) - OSWT(1)) - X2*K*A2)
          A0 = (X2/C)*((OSWR(0) - OSWT(0)) - K*A1)
          F = ((A3*FLU + A2)*FLU + A1)*FLU + A0 -
     A          A0*DEXP(- C*FLU/(X2*K))
          G0 = - (N7*R45L/B)*(OSWR(0) - OSWT(0)) -
     A           (N7*(X1 - R43/R53)/(X3*B))*(OSWR(0) + X2*OSWT(0)) +
     B           (X3*C*N7*PAMB)/(X4*B) - N8*PRESS/FLU
          G1 = - (N7*R45L/B)*(OSWR(1) - OSWT(1)) -
     A           (N7*(X1 - R43/R53)/(X3*B))*(OSWR(1) + X2*OSWT(1))
          G2 = - (N7*R45L/B)*(OSWR(2) - OSWT(2)) -
     A           (N7*(X1 - R43/R53)/(X3*B))*(OSWR(2) + X2*OSWT(2))
          G3 = - (N7*R45L/B)*(OSWR(3) - OSWT(3)) -
     A           (N7*(X1 - R43/R53)/(X3*B))*(OSWR(3) + X2*OSWT(3))

          RND = X1/(B11*B22 - B21*B12)
          N1 = (- A23*B22 + A24*B21)*RND
          N2 = (  A23*B12 - A24*B11)*RND
          N3 = (- A23*B22 + A54*B21)*RND
          N4 = (  A23*B12 - A54*B11)*RND
          H0 = (N1/N5)*H0 - ((N2*R45L)/B)*(OSWR(0) - OSWT(0)) -
     A         (N2*(X1 - R43/R53)/(X3*B))*(OSWR(0) + X2*OSWT(0)) + 
     B         (X3*C*N2*PAMB)/(X4*B)
          H1 = (N1/N5)*H1 - ((N2*R45L)/B)*(OSWR(1) - OSWT(1)) -
     A         (N2*(X1 - R43/R53)/(X3*B))*(OSWR(1) + X2*OSWT(1))
          H2 = (N1/N5)*H2 - ((N2*R45L)/B)*(OSWR(2) - OSWT(2)) -
     A         (N2*(X1 -  R43/R53)/(X3*B))*(OSWR(2) + X2*OSWT(2))
          H3 = (N1/N5)*H3 - ((N2*R45L)/B)*(OSWR(3) - OSWT(3)) -
     A         (N2*(X1 - R43/R53)/(X3*B))*(OSWR(3) + X2*OSWT(3))
          G0 = (N4/N7)*G0 - ((N3*R32L)/A)*(ISWR(0) - ISWT(0)) -
     A         (N3*(X1 - R33/R23)/(X3*A))*(ISWR(0) + X2*ISWT(0)) +
     B         (X3*K*N3/(X2*A))*(PRESS/FLU) + N4*N8*PRESS/FLU/N7
          G1 = (N4/N7)*G1 - ((N3*R32L)/A)*(ISWR(1) - ISWT(1)) -
     A         (N3*(X1 - R33/R23)/(X3*A))*(ISWR(1) + X2*ISWT(1)) +
     B         (X3*C*N3*PRESS/FLU)/(X4*A)
          G2 = (N4/N7)*G2 - ((N3*R32L)/A)*(ISWR(2) - ISWT(2)) -
     A         (N3*(X1 - R33/R23)/(X3*A))*(ISWR(2) + X2*ISWT(2))
          G3 = (N4/N7)*G3 - ((N3*R32L)/A)*(ISWR(3) - ISWT(3)) -
     A         (N3*(X1 - R33/R23)/(X3*A))*(ISWR(3) + X2*ISWT(3))
          Q1 = (X3*C/X4)*(N4/B + N1/A)
          Q2 = (X9*(C**2)/(X16*A*B))*(N4*N1 - N2*N3)
          K3 = (- N4*H3 + N2*G3)*(X4*A)/(X3*C*(N4*N1 - N2*N3))
          K2 = (X3*H3 - X3*C*N4*H2/(X4*B) + X3*C*N2*G2/(X4*B) +
     A          X3*Q1*K3)/Q2
          K1 = (X2*H2 - X3*C*N4*H1/(X4*B) + X3*C*N2*G1/(X4*B) -
     A          X6*K3 + X2*Q1*K2)/Q2
          K0 = (H1 - X3*C*N4*H0/(X4*B) + X3*C*N2*G0/(X4*B) - 
     A          X2*K2 + Q1*K1)/Q2
C
C  Roots of quadratic:  Z^2 + Q1 Z + Q2 = 0, are real
C
          Z = DSQRT(Q1*Q1 - X4*Q2)
          M1 = (Q1 + Z)/X2
          M2 = (Q1 - Z)/X2
          C2 = ((X3*PAMB*K*N2/(X2*B) - K0)*(- M1 + X3*C*N1/(X4*A)) +
     A           X9*C*N2*K*N4*PAMB/(X8*B*B) + H0 - K1 +
     B           X3*C*N1*K0/(X4*A))/(M2 - M1)
          C1 = - C2 + X3*PAMB*K*N2/(X2*B) - K0
C
C  Radial and tangential stress in SiC and OPyC
C
          SIGRI = C1*DEXP(M1*FLU) + C2*DEXP(M2*FLU) + K0 +
     A            ((K3*FLU + K2)*FLU + K1)*FLU
          SIGRO = (- X4*B/(X3*C*N2))*(C1*( - M1 + X3*C*N1/(X4*A))*
     A             DEXP(M1*FLU) + C2*( - M2 + X3*C*N1/(X4*A))*
     B             DEXP(M2*FLU) + H0 - K1 + X3*C*N1*K0/(X4*A) +
     C            (H1 - X2*K2 + X3*C*N1*K1/(X4*A))*FLU +
     D            (H2 - X3*K3 + X3*C*N1*K2/(X4*A))*FLU**2 +
     E            (H3 + X3*C*N1*K3/(X4*A))*FLU**3)
          SIGTO = SIGRO + (X3/(X2*(R43/R53 - X1)))*
     A            (PAMB + SIGRO - X2*F*R45L/X3) + F/X3
        END IF
      END IF
      A3 = (X2/C)*(ISWR(3) - ISWT(3))
      A2 = (X2/C)*((ISWR(2) - ISWT(2)) - X3*K*A3)
      A1 = (X2/C)*((ISWR(1) - ISWT(1)) - X2*K*A2)
      A0 = (X2/C)*((ISWR(0) - ISWT(0)) - K*A1)
      FI = ((A3*FLU + A2)*FLU + A1)*FLU + A0 -
     A      A0*DEXP(- C*FLU/(X2*K))
      A3 = (X2/C)*(OSWR(3) - OSWT(3))
      A2 = (X2/C)*((OSWR(2) - OSWT(2)) - X3*K*A3)
      A1 = (X2/C)*((OSWR(1) - OSWT(1)) - X2*K*A2)
      A0 = (X2/C)*((OSWR(0) - OSWT(0)) - K*A1)
      FO = ((A3*FLU + A2)*FLU + A1)*FLU + A0 -
     A      A0*DEXP(- C*FLU/(X2*K))
      RETURN
      END
C
