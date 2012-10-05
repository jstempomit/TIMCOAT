C***********************************************************************
C                                                                      *
C  Subroutine M_ANALYSIS(MCODE,PRESS,PAMB,FLU,SIGR,SIGT,EPIR,EPIT,UR)  *
C                                                                      *
C  Author:      Jing Wang, NED MIT, April 21, 2003                     *
C                                                                      *
C  Description: Use IMSL math libraries to apply series solutions and  *
C               calculate stress, strain and displacement distributions*
C			  in TRISO particles.                                    *
C                                                                      *
C  Improvements: Poisson's ratio for creep could be changed instead of *
C               being fixed at 0.5 in earlier formulation.             *
C                                                                      *
C  Assumption:  The analysis holds spherical symmetry                  *
C               The external pressure is constant during irradiation.  *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual argument description                                         *
C                                                                      *
C    MCODE           C: Indicates the type of analysis to perform      *
C					 'ISO3': full three-layer analysis               *
C					 'IS2' : IPyC/SiC two-layer analysis             *
C					 'SO2' : SiC/OPyC two-layer analysis             *
C					 'S1'  : SiC single-layer analysis               *
C    PRESS (MPa)    D: Internal gas pressure                           *
C    PAMB (MPa)     D: External ambient pressure                       *
C    FLU (10^21nvt) D: Fluence at time of calculation                  *
C    SIGR(1:NDIV) (MPa) -- returned quantity                           *
C                   D: Radial stress distribution over three layers    *
C    SIGT(1:NDIV) (MPa) -- returned quantity                           *
C                   D: Tangential stress distribution over three layers*
C    EPIR(1:NDIV) -- returned quantity                                 *
C                   D: Radial strain distribution over three layers    *
C    EPIT(1:NDIV) -- returned quantity                                 *
C                   D: Tangential strain distribution over three layers*
C    UR(1:NDIV) (um) -- returned quantity                              *
C				  D: Radial displacement distribution over three     *
C					 layers                                          *
C                                                                      *
C  Variables in COMMON blocks                                          *
C  /PAR_R/ :  fuel particle current geometry                           *
C    R1  - R5 (um)  D: Current radii of layers                         *
C  /PAR_DIV/ --                                                        *
C    NDIVI          I: Number of divisions in IPyC (For stress distr.) *
C    NDIVO          I: Number of divisions in OPyC (For stress distr.) *
C    NDIVS          I: Number of divisions in SiC  (For stress distr.) *
C  /PAR_M/ :  properties of layers used in mechanical analysis         *
C    IPYCALPHA (1/C)D: Thermal expansion coefficient of IPyC           *
C    IPYCCNU        D: Creep Poisson's ratio of IPyC                   *
C    IPYCE (MPa)    D: Young's modulus of IPyC                         *
C    IPYCNU         D: Poisson's ratio of IPyC                         *
C    IPYCREEP (strain/MPa.10**21 nvt) --                               *
C                   D: Creep coefficient of IPyC                       *
C    OPYCD (g/cm^3) D: Density of OPyC                                 *
C    OPYCALPHA (1/C)D: Thermal expansion coefficient of OPyC           *
C    OPYCCNU        D: Creep Poisson's ratio of OPyC                   *
C    OPYCE (MPa)    D: Young's modulus of OPyC                         *
C    OPYCNU         D: Poisson's ratio of OPyC                         *
C    OPYCREEP (strain/MPa.10**21 nvt) --                               *
C                   D: Creep coefficient of OPyC                       *
C    OPYCD (g/cm^3) D: Density of OPyC                                 *
C    SICALPHA (1/C) D: Thermal expansion coefficient of SiC            *
C    SICE (MPa)     D: Young's modulus of SiC                          *
C    SICNU          D: Poisson's ratio of SiC                          *
C    ISWR(0:NDEG)   D: Array of length NDEG+1 storing coefficients of  *
C                      the polynomial of IPyC radial swelling rate     *
C    ISWT(0:NDEG)   D: Array of length NDEG+1 storing coefficients of  *
C                      the polynomial of IPyC tangential swelling rate *
C    OSWR(0:NDEG)   D: Array of length NDEG+1 storing coefficients of  *
C                      the polynomial of OPyC radial swelling rate     *
C    OSWT(0:NDEG)   D: Array of length NDEG+1 storing coefficients of  *
C                      the polynomial of OPyC tangential swelling rate *
C  /CRACKED_PYC/ :  Quantities related to cracked PyC layers           *
C    EPIRCP(1:NDIV) D: Array recording the elastic radial strains of   *
C					 fully relaxed PyC layers at point of cracking.  *
C					 The symbol means 'Epsilon R of C prime'         *
C    EPITCP(1:NDIV) D: Array recording the elastic tangential strains  *
C					 of fully relaxed PyC layers at point of cracking*
C					 The symbol means 'Epsilon T of C prime'         *
C    URCP(1:NDIV)   D: Array recording the elastic radial displacement *
C					 of fully relaxed PyC layers at point of cracking*
C					 The symbol means 'Ur of C prime'                *
C    KIIPYC(MPa.um^1/2) D:                                             *
C					 Stress intensity factor in IPyC layer           *
C    KIOPYC(MPa.um^1/2) D:                                             *
C					 Stress intensity factor in OPyC layer           *
C    KI1 (MPa.um^1/2)D: Stress intensity factor from IPyC crack        *
C    KI2 (MPa.um^1/2)D: Stress intensity factor from OPyC crack        *
C    SHEARIPYC(MPa.um)  D:                                             *
C					 The shear force per unit length on SiC surface  *
C                      induced by IPyC crack                           *
C    SHEAROPYC(MPa.um)  D:                                             *
C					 The shear force per unit length on SiC surface  *
C                      induced by OPyC crack                           *
C    DF(10^21nvt)   D: Incremental fluence since last step             *
C                                                                      *
C  Local argument description                                          *
C    P(LDA,LDA), Q(NDIM), X(NDIM) --                                   *
C                   D: Arrays for linear equations                     *
C    A(0:IORDER+1), B, D, F, M, N --                                   *
C                   D: Coefficients for stresses                       *
C    SIG* (MPa)     D: Intermediate stress variables for calculations  *
C    EPIIR, EPIIT   D: The i-th order term of radial and tangential    *
C					 strains                                         *
C    UIR            D: The i-th order term of radial displacement      *
C    EPITP          D: Tangential strain induced by shear force at     *
C                      PyC/SiC interfaces.                             *
C    URP (um)       D: Radial displacement induced by shear force at   *
C					 PyC/SiC interfaces.                             *
C    AIPYC (um)     D: The crack length in IPyC layer                  *
C    AOPYC (um)     D: The crack length in OPyC layer                  *
C    SIGTPIBAR      D: Average tangential stress due to shear force at *
C					 IPyC/SiC interface.                             *
C    SIGTPOBAR      D: Average tangential stress due to shear force at *
C					 OPyC/SiC interface.                             *
C    SWELLRIPYC (1/10^21nvt) D:                                        *
C                      Current radial swelling rate in IPyC.           *
C    SWELLROPYC (1/10^21nvt) D:                                        *
C                      Current radial swelling rate in OPyC.           *
C    SWELLTIPYC (1/10^21nvt) D:                                        *
C                      Current tangential swelling rate in IPyC.       *
C    SWELLTOPYC (1/10^21nvt) D:                                        *
C                      Current tangential swelling rate in OPyC.       *
C    WSHEARIPYC (MPa.um^3) D:                                          *
C					 The shear strain energy associated with a crack *
C					 in IPyC layer                                   *
C    WSHEAROPYC (MPa.um^3) D:                                          *
C					 The shear strain energy associated with a crack *
C					 in OPyC layer                                   *
C    L*, G*, V*, C*, P*, F*, S* --                                     *
C                   D: Intermediate variables                          *
C                                                                      *
C  Parameters, counters and others                                     *
C    NDIV           I: Number of divisions for stress distributions    *
C    IORDER         I: Power of the series solution for stresses       *
C    NDEG           I: Degree of polynomial of swelling rate           *
C    IPATH, LDA, NDIM --                                               *
C                   I: Variables for calling IMSL subroutine DLSARG    *
C    CRITERIA       D: Convergence criterium for series solution       *
C    R (um)         D: Radius                                          *
C    FLUI           D: Variable for series summation                   *
C    PROBER, PROBET D: Variables for statistics                        *
C    CONV           L: Flag of whether the series solution converge    *
C                                                                      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
	SUBROUTINE M_ANALYSIS(MCODE,PRESS,PAMB,FLU,SIGR,SIGT,EPIR,EPIT,UR)
C
      USE IMSL_LIBRARIES          ! Allow usage of IMSL libraries
	PARAMETER ( NDIV = 30, IORDER = 200 )
      PARAMETER ( NDEG = 3)
	PARAMETER ( IPATH = 1, LDA = 6, NDIM = 6, NDIM2 = 4)
      PARAMETER ( PIE = 3.1415926535897932385 D0)
	PARAMETER ( SCALE_R = 3.0D2, SCALE_F = 1.5D0, SCALE_E = 1.0D3)
	PARAMETER ( CRITERIA = 1.0D-4 )
C
	DOUBLE PRECISION PRESS, PAMB, FLU
	DOUBLE PRECISION R1, R2, R3, R4, R5
	DOUBLE PRECISION IPYCALPHA, IPYCCNU, IPYCE, IPYCNU, IPYCREEP,
     &                 OPYCALPHA, OPYCCNU, OPYCE, OPYCNU, OPYCREEP,
     &                 SICALPHA, SICE, SICNU, IPYCD, OPYCD
      DOUBLE PRECISION ISWR, ISWT, OSWR, OSWT
	DOUBLE PRECISION EPIRCP, EPITCP, URCP, KIIPYC, KIOPYC, KI1, KI2,
     &				 SHEARIPYC, SHEAROPYC, DF
	DOUBLE PRECISION EPITP, URP, SIGTPIBAR, SIGTPOBAR, AIPYC, AOPYC,
     &				 SWELLTIPYC, SWELLTOPYC, SWELLRIPYC, SWELLROPYC,
     &				 WSHEARIPYC, WSHEAROPYC
	DOUBLE PRECISION SIGR, SIGT, SIGIR, SIGIT, SIGINTER
	DOUBLE PRECISION EPIR, EPIT, EPIIR, EPIIT, UR, UIR
	DOUBLE PRECISION P(LDA,LDA), Q(NDIM), X(NDIM)
	DOUBLE PRECISION A, B, D, F, M, N
	DOUBLE PRECISION L1I, L2I, L3I, L4I, L1O, L2O, L3O, L4O
	DOUBLE PRECISION G1I, G2I, G3I, G1O, G2O, G3O, G1S, G2S
      DOUBLE PRECISION V1I, V2I, V3I, V1O, V2O, V3O
	DOUBLE PRECISION SIGR2, SIGT2, SIGR3, SIGT3, SIGR4, SIGT4
	DOUBLE PRECISION SIGR5, SIGT5, SIG2, SIG3, SIG4, SIG5
	DOUBLE PRECISION SI, SO, SI1, SO1, FI, FO, GI, GO
	DOUBLE PRECISION PROBER, PROBET
	DOUBLE PRECISION TEMP
C  Note: SIGR2, SIGT2, SIGR3, SIGT3, SIGR4, SIGT4, SIGR5, SIGT5, SIG2, SIG3, SIG4, SIG5, 
C        SI, SO, FI, FO, GI, GO change with the order i of stress components
	DOUBLE PRECISION PI, PO, CI, CO, R, DR, FLUI
      INTEGER NDIVI, NDIVS, NDIVO
	INTEGER I, J, K
	CHARACTER*4 MCODE
	LOGICAL CONV
C
	DIMENSION SIGR(1:NDIV), SIGT(1:NDIV)
	DIMENSION SIGIR(1:NDIV), SIGIT(1:NDIV)
	DIMENSION EPIR(1:NDIV), EPIT(1:NDIV), UR(1:NDIV)
	DIMENSION EPIIR(1:NDIV), EPIIT(1:NDIV), UIR(1:NDIV)
      DIMENSION ISWR(0:NDEG), ISWT(0:NDEG), OSWR(0:NDEG), OSWT(0:NDEG)
	DIMENSION EPIRCP(1:NDIV), EPITCP(1:NDIV), URCP(1:NDIV)
	DIMENSION EPITP(1:NDIV), URP(1:NDIV)
	DIMENSION A(0:IORDER+1), B(0:IORDER+1), D(0:IORDER+1),
     &          F(0:IORDER+1), M(0:IORDER+1), N(0:IORDER+1)
C
      COMMON /PAR_R/ R1, R2, R3, R4, R5
	COMMON /PAR_DIV/ NDIVI, NDIVS, NDIVO
      COMMON /PAR_M/ IPYCALPHA, IPYCCNU, IPYCE, IPYCNU, IPYCREEP,
     &               OPYCALPHA, OPYCCNU, OPYCE, OPYCNU, OPYCREEP,
     &               SICALPHA, SICE, SICNU, IPYCD, OPYCD,
     &               ISWR, ISWT, OSWR, OSWT
	COMMON /CRACKED_PYC/ EPIRCP, EPITCP, URCP, KIIPYC, KIOPYC,
     &					 KI1, KI2, SHEARIPYC, SHEAROPYC, DF
C
C
	DATA X0/0.0 D+00/, X1/1.0 D+00/, X2/2.0 D+00/, X3/3.0 D+00/
C  Choose the type of analysis to perform
	IF(MCODE .EQ. 'ISO3') THEN
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments here.
	  R1 = R1/SCALE_R
	  R2 = R2/SCALE_R
	  R3 = R3/SCALE_R
	  R4 = R4/SCALE_R
	  R5 = R5/SCALE_R
	  DO 2805 K = 1,NDIV
	    UR(K) = UR(K)/SCALE_R
2805	  CONTINUE
	  IPYCE = IPYCE/SCALE_E
	  OPYCE = OPYCE/SCALE_E
	  SICE = SICE/SCALE_E
	  PRESS = PRESS/SCALE_E
	  PAMB = PAMB/SCALE_E
	  FLU = FLU/SCALE_F
	  IPYCREEP = IPYCREEP*SCALE_F*SCALE_E
	  OPYCREEP = OPYCREEP*SCALE_F*SCALE_E
	  DO 2806 I=0,NDEG
	    ISWR(I) = ISWR(I)*SCALE_F**(I+1)
	    ISWT(I) = ISWT(I)*SCALE_F**(I+1)
	    OSWR(I) = OSWR(I)*SCALE_F**(I+1)
	    OSWT(I) = OSWT(I)*SCALE_F**(I+1)
2806    CONTINUE
	  DO 2807 K=1,NDIV
	    SIGR(K) = SIGR(K)/SCALE_E
	    SIGT(K) = SIGT(K)/SCALE_E
2807    CONTINUE
C
C  Prepare quantities based on material properties
C
        L1I = X1/IPYCE
        L2I = X2*IPYCNU/IPYCE
        L3I = IPYCNU/IPYCE
        L4I = (X1 - IPYCNU)/IPYCE
        L1O = X1/OPYCE
        L2O = X2*OPYCNU/OPYCE
        L3O = OPYCNU/OPYCE
        L4O = (X1 - OPYCNU)/OPYCE
C
	  G1I = X1/(L4I - L3I)
	  G2I = X1/(L1I + L3I)
	  G3I = L4I/(L1I + L3I)
	  G1O = X1/(L4O - L3O)
	  G2O = X1/(L1O + L3O)
	  G3O = L4O/(L1O + L3O)
	  G1S = SICE/(X1 - X2*SICNU)
	  G2S = X2*SICE/(X1 + SICNU)
C
	  V1I = X1 - X2*IPYCCNU
	  V2I = X2*(X1 + IPYCCNU)/X3
	  V3I = X2*(L1I*IPYCCNU-L3I)*G1I
	  V1O = X1 - X2*OPYCCNU
	  V2O = X2*(X1 + OPYCCNU)/X3
	  V3O = X2*(L1O*OPYCCNU-L3O)*G1O
C
C  Prepare fluence, pressure and creep quantities
	  FLUI = 1.0D+00
C    Stresses are additive; boundary conditions should exclude residual stresses
        PI = -PRESS-SIGR(1)
	  PO = -PAMB-SIGR(NDIV)
	  CI = IPYCREEP
	  CO = OPYCREEP
C
C  Initialize arrays P and Q
        DO 2808 I=1,NDIM
	    Q(I)=0.0D0
2808    CONTINUE
C  Calculate initial stresses, strains, and displacement, i.e., I = 0
C    Step 1: Assign values to P(LDA,LDA). It's fixed for all orders of I.
        P(1,1) = 1.0D0
	  P(1,2) = -X2*G2I*(R2**(-3.0))/G1I
	  P(1,3) = 0.0D0
	  P(1,4) = 0.0D0
	  P(1,5) = 0.0D0
	  P(1,6) = 0.0D0
	  P(2,1) = 1.0D0
	  P(2,2) = -X2*G2I*(R3**(-3.0))/G1I
        P(2,3) = -G1S/G1I
	  P(2,4) = G2S*(R3**(-3.0))/G1I
	  P(2,5) = 0.0D0
	  P(2,6) = 0.0D0
	  P(3,1) = 1.0D0
	  P(3,2) = R3**(-3.0)
	  P(3,3) = -1.0D0
	  P(3,4) = -R3**(-3.0)
        P(3,5) = 0.0D0
        P(3,6) = 0.0D0
        P(4,1) = 0.0D0
        P(4,2) = 0.0D0
	  P(4,3) = -1.0D0
	  P(4,4) = -R4**(-3.0)
	  P(4,5) = 1.0D0
	  P(4,6) = R4**(-3.0)
	  P(5,1) = 0.0D0
	  P(5,2) = 0.0D0
	  P(5,3) = -G1S/G1O
	  P(5,4) = G2S*(R4**(-3.0))/G1O
	  P(5,5) = 1.0D0
	  P(5,6) = -X2*G2O*(R4**(-3.0))/G1O
	  P(6,1) = 0.0D0
	  P(6,2) = 0.0D0
	  P(6,3) = 0.0D0
	  P(6,4) = 0.0D0
	  P(6,5) = 1.0D0
	  P(6,6) = -X2*G2O*(R5**(-3.0))/G1O
C
C    Step 2: Assign values to Q(NDIM) for i=0
        Q(1) = PI/G1I
        Q(6) = PO/G1O
C    Step 3: Use DLSARG of IMSL library to calculate X
        CALL DLSARG (NDIM, P, LDA, Q, IPATH, X)
        A(0) = X(1)
	  B(0) = X(2)
	  D(0) = X(3)
	  F(0) = X(4)
	  M(0) = X(5)
	  N(0) = X(6)
C
C    Register initial stresses and strains
        DO 2820 J = 1, NDIVI+2
	    R = (J-1)*(R3-R2)/FLOAT(NDIVI+1) + R2
	    SIGIR(J) = G1I*A(0) - X2*G2I*B(0)*(R**(-3.0))
	    SIGIT(J) = SIGIR(J) + X3*G2I*B(0)*(R**(-3.0))
	    EPIIR(J) = A(0) - X2*B(0)*(R**(-3.0))
	    EPIIT(J) = A(0) + B(0)*(R**(-3.0))
	    UIR(J) = A(0)*R + B(0)*(R**(-2.0))
2820    CONTINUE
        DO 2830 J = 1, NDIVS+2
	    R = (J-1)*(R4-R3)/FLOAT(NDIVS+1) + R3
	    SIGIR(NDIVI+2+J) = G1S*D(0) - G2S*F(0)*(R**(-3.0))
	    SIGIT(NDIVI+2+J) = SIGIR(NDIVI+2+J) + 
     &                       X3*G2S*F(0)*(R**(-3.0))/X2
	    EPIIR(NDIVI+2+J) = D(0) - X2*F(0)*(R**(-3.0))
	    EPIIT(NDIVI+2+J) = D(0) + F(0)*(R**(-3.0))
	    UIR(NDIVI+2+J) = D(0)*R + F(0)*(R**(-2.0))
2830    CONTINUE
        DO 2840 J = 1, NDIVO+2
	    R = (J-1)*(R5-R4)/FLOAT(NDIVO+1) + R4
	    SIGIR(NDIVI+NDIVS+4+J) = G1O*M(0) - X2*G2O*N(0)*(R**(-3.0))
	    SIGIT(NDIVI+NDIVS+4+J) = SIGIR(NDIVI+NDIVS+4+J) + 
     &                             X3*G2O*N(0)*(R**(-3.0))
	    EPIIR(NDIVI+NDIVS+4+J) = M(0) - X2*N(0)*(R**(-3.0))
	    EPIIT(NDIVI+NDIVS+4+J) = M(0) + N(0)*(R**(-3.0))
	    UIR(NDIVI+NDIVS+4+J) = M(0)*R + N(0)*(R**(-2.0))
2840    CONTINUE
C
C    Add to residual quantities from last fuel cycle
	  DO 2850 J = 1, NDIV
	    SIGR(J) = SIGIR(J) + SIGR(J)
          SIGT(J) = SIGIT(J) + SIGT(J)
	    EPIR(J) = EPIIR(J) + EPIR(J)
	    EPIT(J) = EPIIT(J) + EPIT(J)
	    UR(J)   = UIR(J)   + UR(J)
2850    CONTINUE
C
        IF(FLU.EQ.0.0D0) THEN
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments at the beginning of this subroutine and 
C  change them back here before returning to the main program.
	    R1 = R1*SCALE_R
	    R2 = R2*SCALE_R
	    R3 = R3*SCALE_R
	    R4 = R4*SCALE_R
	    R5 = R5*SCALE_R
	    DO 2809 K = 1,NDIV
	      UR(K) = UR(K)*SCALE_R
2809	    CONTINUE
	    IPYCE = IPYCE*SCALE_E
	    OPYCE = OPYCE*SCALE_E
	    SICE = SICE*SCALE_E
	    PRESS = PRESS*SCALE_E
	    PAMB = PAMB*SCALE_E
	    FLU = FLU*SCALE_F
	    IPYCREEP = IPYCREEP/SCALE_F/SCALE_E
	    OPYCREEP = OPYCREEP/SCALE_F/SCALE_E
	    DO 2810 K=0,NDEG
	      ISWR(K) = ISWR(K)/(SCALE_F**(K+1))
	      ISWT(K) = ISWT(K)/(SCALE_F**(K+1))
	      OSWR(K) = OSWR(K)/(SCALE_F**(K+1))
	      OSWT(K) = OSWT(K)/(SCALE_F**(K+1))
2810      CONTINUE
	    DO 2811 K=1,NDIV
	      SIGR(K) = SIGR(K)*SCALE_E
	      SIGT(K) = SIGT(K)*SCALE_E
2811      CONTINUE
	    RETURN  !no need for calculating higher order terms
	  END IF
C    Step 4: Prepare stress entries for higher order calculations
        DO 2860 J=1,NDIV
	    SIGIR(J) = SIGR(J)
	    SIGIT(J) = SIGT(J)
2860    CONTINUE
        SIGR2 = SIGIR(1)
	  SIGT2 = SIGIT(1)
        SIGR3 = SIGIR(NDIVI+2)
	  SIGT3 = SIGIT(NDIVI+2)
	  SIGR4 = SIGIR(NDIVI+NDIVS+5)
	  SIGT4 = SIGIT(NDIVI+NDIVS+5)
	  SIGR5 = SIGIR(NDIV)
	  SIGT5 = SIGIT(NDIV)
	  SIG2 = (X1+L3I*G1I*V1I)*SIGR2 + (L1I*G1I*V1I-X1)*SIGT2
     	  SIG3 = (X1+L3I*G1I*V1I)*SIGR3 + (L1I*G1I*V1I-X1)*SIGT3
        SIG4 = (X1+L3O*G1O*V1O)*SIGR4 + (L1O*G1O*V1O-X1)*SIGT4
        SIG5 = (X1+L3O*G1O*V1O)*SIGR5 + (L1O*G1O*V1O-X1)*SIGT5
C
C  Begin to calculate higher order terms
	  GI = 0.0
	  GO = 0.0
C  Order loop
        CONV = .FALSE.
        DO 2800 I=0, IORDER
C  Initialize probes for testing convergence
          PROBER=0.
	    PROBET=0.
C
	    FLUI = FLUI*FLU
C
          IF (I .LE. NDEG) THEN
	      SI = (L4I*ISWR(I)+L2I*ISWT(I))/((I+1.0)*(L1I*L4I-L2I*L3I))
	      SO = (L4O*OSWR(I)+L2O*OSWT(I))/((I+1.0)*(L1O*L4O-L2O*L3O))
	      SI1 = G2I*(ISWR(I)-ISWT(I))/(I+1.0)
	      SO1 = G2O*(OSWR(I)-OSWT(I))/(I+1.0)
	      FI = (-CI)*V3I*GI/(X2*L4I*(I+1.0))+(ISWR(I)-ISWT(I))/
     &           (L4I*(I+1.0))
	      FO = (-CO)*V3O*GO/(X2*L4O*(I+1.0))+(OSWR(I)-OSWT(I))/
     &           (L4O*(I+1.0))
	      GI = FI - CI*G1I*V1I*GI/(I+1.0)
	      GO = FO - CO*G1O*V1O*GO/(I+1.0)
	    ELSE
	      SI = 0.0
	      SO = 0.0
	      SI1 = 0.0
	      SO1 = 0.0
	      FI = (-CI)*V3I*GI/(X2*L4I*(I+1.0))
	      FO = (-CO)*V3O*GO/(X2*L4O*(I+1.0))
	      GI = FI - CI*G1I*V1I*GI/(I+1.0)
	      GO = FO - CO*G1O*V1O*GO/(I+1.0)
	    END IF
C
C  Calculate A(I+1), B(I+1), D(I+1), F(I+1), M(I+1), N(I+1)
C    Step 1: Assign new values to Q(NDIM)
	    Q(1) = (-X2*FI*(DLOG(R2)+G3I)/X3+CI*G2I*SIG2+SI)/G1I
	    Q(2) = (-X2*FI*(DLOG(R3)+G3I)/X3+CI*G2I*SIG3+SI)/G1I
	    Q(3) = -X2*FI*DLOG(R3)/(X3*G1I)
	    Q(4) = -X2*FO*DLOG(R4)/(X3*G1O)
	    Q(5) = (-X2*FO*(DLOG(R4)+G3O)/X3+CO*G2O*SIG4+SO)/G1O
	    Q(6) = (-X2*FO*(DLOG(R5)+G3O)/X3+CO*G2O*SIG5+SO)/G1O
C    Step 2: Use DLSARG of IMSL library to calculate X
          CALL DLSARG (NDIM, P, LDA, Q, IPATH, X)
          A(I+1) = X(1)
	    B(I+1) = X(2)
	    D(I+1) = X(3)
	    F(I+1) = X(4)
	    M(I+1) = X(5)
	    N(I+1) = X(6)
C
C  Calculate stresses, strains, and displacement of order i+1
          DO 2870 J = 1, NDIVI+2
	      R = (J-1)*(R3-R2)/FLOAT(NDIVI+1) + R2
	      SIGINTER = SIGIR(J)
	      SIGIR(J) = G1I*A(I+1) - X2*G2I*B(I+1)*(R**(-3.0)) + X2*FI*
     A               (DLOG(R)+G3I)/X3 - CI*G2I*((X1+L3I*G1I*V1I)*
     B               SIGIR(J)+(L1I*G1I*V1I-X1)*SIGIT(J))/(I+1.0)- SI
	      SIGIT(J) = SIGIR(J) + X3*G2I*B(I+1)*(R**(-3.0)) - X2*FI*
     A               G2I/(X3*G1I) + X3*CI*G2I*V2I*(SIGINTER-SIGIT(J))/
     B               (X2*(I+1.0)) + SI1
	      EPIIR(J) = A(I+1) - X2*B(I+1)*(R**(-3.0)) + 
     A	  			  X2*(L4I-L3I)*FI*(DLOG(R)+X1)/X3
	      EPIIT(J) = A(I+1) + B(I+1)*(R**(-3.0)) +
     A				  X2*(L4I-L3I)*FI*DLOG(R)/X3
	      UIR(J)   = A(I+1)*R + B(I+1)*(R**(-2.0)) +
     A				  X2*(L4I-L3I)*FI*R*DLOG(R)/X3
2870      CONTINUE
          DO 2880 J = 1, NDIVS+2
	      R = (J-1)*(R4-R3)/FLOAT(NDIVS+1) + R3
	      SIGINTER = SIGIR(NDIVI+2+J)
	      SIGIR(NDIVI+2+J) = G1S*D(I+1) - G2S*F(I+1)*(R**(-3.0))
	      SIGIT(NDIVI+2+J) = SIGIR(NDIVI+2+J) + X3*G2S*F(I+1)*
     &                         (R**(-3.0))/X2
	      EPIIR(NDIVI+2+J) = D(I+1) - X2*F(I+1)*(R**(-3.0))
	      EPIIT(NDIVI+2+J) = D(I+1) + F(I+1)*(R**(-3.0))
	      UIR(NDIVI+2+J)   = D(I+1)*R + F(I+1)*(R**(-2.0))
2880      CONTINUE
	    DO 2890 J = 1, NDIVO+2
	      R = (J-1)*(R5-R4)/FLOAT(NDIVO+1) + R4
	      SIGINTER = SIGIR(NDIVI+NDIVS+4+J)
	      SIGIR(NDIVI+NDIVS+4+J) = G1O*M(I+1) - X2*G2O*N(I+1)*
     A             (R**(-3.0)) + X2*FO*(DLOG(R)+G3O)/X3 - CO*G2O*
     B             ((X1+L3O*G1O*V1O)*SIGIR(NDIVI+NDIVS+4+J)+
     C             (L1O*G1O*V1O-X1)*SIGIT(NDIVI+NDIVS+4+J))/(I+1.0) - SO
	      SIGIT(NDIVI+NDIVS+4+J) = SIGIR(NDIVI+NDIVS+4+J) + 
     A             X3*G2O*N(I+1)*(R**(-3.0)) - X2*FO*G2O/(X3*G1O) + 
     B             X3*CO*G2O*V2O*(SIGINTER-SIGIT(NDIVI+NDIVS+4+J))/
     C             (X2*(I+1.0)) + SO1
	      EPIIR(NDIVI+NDIVS+4+J) = M(I+1) - X2*N(I+1)*(R**(-3.0)) + 
     A							  X2*(L4O-L3O)*FO*(DLOG(R)+X1)/X3
	      EPIIT(NDIVI+NDIVS+4+J) = M(I+1) + N(I+1)*(R**(-3.0)) +
     A							  X2*(L4O-L3O)*FO*DLOG(R)/X3
	      UIR(NDIVI+NDIVS+4+J)   = M(I+1)*R + N(I+1)*(R**(-2.0)) +
     A							  X2*(L4O-L3O)*FO*R*DLOG(R)/X3
2890      CONTINUE
C
	    DO 2900 J=1, NDIV
	      SIGR(J) = SIGR(J) + SIGIR(J)*FLUI
	      SIGT(J) = SIGT(J) + SIGIT(J)*FLUI
		  EPIR(J) = EPIR(J) + EPIIR(J)*FLUI
	      EPIT(J) = EPIT(J) + EPIIT(J)*FLUI
		  UR(J)   = UR(J)   + UIR(J)*FLUI
C
C  Calculate the average-square-root of the term of order i
            PROBER = PROBER + (SIGIR(J)*FLUI/SIGR(J))**2
            PROBET = PROBET + (SIGIT(J)*FLUI/SIGT(J))**2
2900      CONTINUE
          PROBER = DSQRT(PROBER/NDIV)
          PROBET = DSQRT(PROBET/NDIV)
	    IF((PROBER.LT.CRITERIA).AND.(PROBET.LT.CRITERIA)) THEN
	      IF(CONV.EQ..FALSE.) THEN
	        CONV = .TRUE.
	      ELSE
C	        WRITE(*,*) 'Convergence achieved at order i = ', I
C
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments at the beginning of this subroutine and 
C  change them back here before returning to the main program.
	        R1 = R1*SCALE_R
	        R2 = R2*SCALE_R
	        R3 = R3*SCALE_R
	        R4 = R4*SCALE_R
	        R5 = R5*SCALE_R
		    DO 2812 K = 1,NDIV
			  UR(K) = UR(K)*SCALE_R
2812		    CONTINUE
	        IPYCE = IPYCE*SCALE_E
	        OPYCE = OPYCE*SCALE_E
	        SICE = SICE*SCALE_E
	        PRESS = PRESS*SCALE_E
	        PAMB = PAMB*SCALE_E
	        FLU = FLU*SCALE_F
	        IPYCREEP = IPYCREEP/SCALE_F/SCALE_E
	        OPYCREEP = OPYCREEP/SCALE_F/SCALE_E
	        DO 2813 K=0,NDEG
	          ISWR(K) = ISWR(K)/(SCALE_F**(K+1))
	          ISWT(K) = ISWT(K)/(SCALE_F**(K+1))
	          OSWR(K) = OSWR(K)/(SCALE_F**(K+1)) 
	          OSWT(K) = OSWT(K)/(SCALE_F**(K+1))
2813          CONTINUE
	        DO 2814 K=1,NDIV
	          SIGR(K) = SIGR(K)*SCALE_E
	          SIGT(K) = SIGT(K)*SCALE_E
2814          CONTINUE
	        RETURN
	      ENDIF
	    ELSE IF(CONV.EQ..TRUE.) THEN
C           WRITE(*,*) 'At Order = ',I, 'abormal'
	      CONV = .FALSE.
	    ENDIF
C
C  Update stress entries for higher orders
          SIGINTER = SIGR2
          SIGR2 = G1I*A(I+1) - X2*G2I*B(I+1)*(R2**(-3.0)) + X2*FI*
     A           (DLOG(R2)+G3I)/X3 - CI*G2I*((X1+L3I*G1I*V1I)*SIGR2 +
     B           (L1I*G1I*V1I-X1)*SIGT2)/(I+1.0) - SI
          SIGT2 = SIGR2 + X3*G2I*B(I+1)*(R2**(-3.0)) - X2*FI*
     A           G2I/(X3*G1I) + X3*CI*G2I*V2I*(SIGINTER-SIGT2)/
     B           (X2*(I+1.0)) + SI1
          SIGINTER = SIGR3
          SIGR3 = G1I*A(I+1) - X2*G2I*B(I+1)*(R3**(-3.0)) + X2*FI*
     A           (DLOG(R3)+G3I)/X3 - CI*G2I*((X1+L3I*G1I*V1I)*SIGR3 +
     B           (L1I*G1I*V1I-X1)*SIGT3)/(I+1.0) - SI
          SIGT3 = SIGR3 + X3*G2I*B(I+1)*(R3**(-3.0)) - X2*FI*
     A           G2I/(X3*G1I) + X3*CI*G2I*V2I*(SIGINTER-SIGT3)/
     B           (X2*(I+1.0)) + SI1
	    SIGINTER = SIGR4
	    SIGR4 = G1O*M(I+1) - X2*G2O*N(I+1)*(R4**(-3.0)) + X2*FO*
     A           (DLOG(R4)+G3O)/X3 - CO*G2O*((X1+L3O*G1O*V1O)*SIGR4 +
     B           (L1O*G1O*V1O-X1)*SIGT4)/(I+1.0) - SO
	    SIGT4 = SIGR4 + X3*G2O*N(I+1)*(R4**(-3.0)) - X2*FO*
     A           G2O/(X3*G1O) + X3*CO*G2O*V2O*(SIGINTER-SIGT4)/
     B           (X2*(I+1.0)) + SO1
	    SIGINTER = SIGR5
	    SIGR5 = G1O*M(I+1) - X2*G2O*N(I+1)*(R5**(-3.0)) + X2*FO*
     A           (DLOG(R5)+G3O)/X3 - CO*G2O*((X1+L3O*G1O*V1O)*SIGR5 +
     B           (L1O*G1O*V1O-X1)*SIGT5)/(I+1.0) - SO
	    SIGT5 = SIGR5 + X3*G2O*N(I+1)*(R5**(-3.0)) - X2*FO*
     A           G2O/(X3*G1O) + X3*CO*G2O*V2O*(SIGINTER-SIGT5)/
     B           (X2*(I+1.0)) + SO1
	    SIG2 = ((X1+L3I*G1I*V1I)*SIGR2+(L1I*G1I*V1I-X1)*SIGT2)/(I+2.0)
	    SIG3 = ((X1+L3I*G1I*V1I)*SIGR3+(L1I*G1I*V1I-X1)*SIGT3)/(I+2.0)
	    SIG4 = ((X1+L3O*G1O*V1O)*SIGR4+(L1O*G1O*V1O-X1)*SIGT4)/(I+2.0)
	    SIG5 = ((X1+L3O*G1O*V1O)*SIGR5+(L1O*G1O*V1O-X1)*SIGT5)/(I+2.0)
2800    CONTINUE
        IF (I.GE.IORDER) WRITE(*,*) 'Stress has not converge at IORDER.'
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments at the beginning of this subroutine and 
C  change them back here before returning to the main program.
	  R1 = R1*SCALE_R
	  R2 = R2*SCALE_R
	  R3 = R3*SCALE_R
	  R4 = R4*SCALE_R
	  R5 = R5*SCALE_R
	  DO 2815 K = 1,NDIV
	    UR(K) = UR(K)*SCALE_R
2815	  CONTINUE
	  IPYCE = IPYCE*SCALE_E
	  OPYCE = OPYCE*SCALE_E
	  SICE = SICE*SCALE_E
	  PRESS = PRESS*SCALE_E
	  PAMB = PAMB*SCALE_E
	  FLU = FLU*SCALE_F
	  IPYCREEP = IPYCREEP/SCALE_F/SCALE_E
	  OPYCREEP = OPYCREEP/SCALE_F/SCALE_E
	  DO 2816 K=0,NDEG
	    ISWR(K) = ISWR(K)/(SCALE_F**(K+1))
	    ISWT(K) = ISWT(K)/(SCALE_F**(K+1))
	    OSWR(K) = OSWR(K)/(SCALE_F**(K+1))
	    OSWT(K) = OSWT(K)/(SCALE_F**(K+1))
2816    CONTINUE
	  DO 2817 K=1,NDIV
	    SIGR(K) = SIGR(K)*SCALE_E
	    SIGT(K) = SIGT(K)*SCALE_E
2817    CONTINUE
C
	ELSE IF(MCODE .EQ. 'IS2') THEN
C  1) First evaluate the stresses, strains, and displacements in the cracked OPyC layer,
C  and calculate the shear force it applies to the SiC layer.
	  AOPYC = 0.99D0*(R5-R4)
	  SWELLROPYC = 0.0D0
	  SWELLTOPYC = 0.0D0
	  DO 3116 K = 0, NDEG
		SWELLROPYC = SWELLROPYC + OSWR(K)*FLU**K
		SWELLTOPYC = SWELLTOPYC + OSWT(K)*FLU**K
3116	  CONTINUE
	  DO 3136 K = 1, NDIVO + 2
		URP(NDIVI+NDIVS+4+K) = 0.0D0
3136	  CONTINUE
	  DR = (R5-R4)/FLOAT(NDIVO+1)
	  CONV = .FALSE.
	  DO WHILE (.NOT.CONV)
		TEMP = KIOPYC
		SIGTPOBAR = 0.0D0
		DO 3117 J = 1, NDIVO + 2
		  R = (J-1)*DR + R4
		  EPITP(NDIVI+NDIVS+4+J)=((URCP(NDIVI+NDIVS+4+J)+
     &		URP(NDIVI+NDIVS+4+J))/R -
     &		EPITCP(NDIVI+NDIVS+4+J)-SWELLTOPYC*DF-8.0D0*
     &		(X1-OPYCNU**2.0)*KIOPYC*SQRT(ABS(R-R4)/(X2*PIE*R**2.0))
     &		/PIE/OPYCE)/(X1+OPYCE*OPYCREEP*DF)
		  SIGTPOBAR = SIGTPOBAR + OPYCE*EPITP(NDIVI+NDIVS+4+J)
3117		CONTINUE
		SIGTPOBAR = SIGTPOBAR/FLOAT(NDIVO+2)
	    KIOPYC = 0.413*(X1+R5/R4)*SIGTPOBAR*SQRT(PIE*AOPYC)/
     &			 SQRT(X1 - AOPYC/(R5-R4))
C  Approaches the right value of KIOPYC
		KIOPYC = (KIOPYC + TEMP)/X2
	    IF(ABS((TEMP-KIOPYC)/KIOPYC).LE.CRITERIA) CONV = .TRUE.
C  Recalculate URP
		DO 3137 J = 1, NDIVO + 2
		  URP(NDIVI+NDIVS+4+J) = 0.0D0
		  DO 3138 K = 1, J - 1
			URP(NDIVI+NDIVS+4+J) = URP(NDIVI+NDIVS+4+J) - 
     &				OPYCNU*EPITP(NDIVI+NDIVS+4+K)*DR
3138		  CONTINUE
3137		CONTINUE
	  END DO
C  Update EPITCP, and record resulting stresses, strains, and displacements in the cracked OPyC layer
	  DO 3118 J = 1, NDIVO + 2
		EPITCP(NDIVI+NDIVS+4+J) = EPITCP(NDIVI+NDIVS+4+J) + 
     &				OPYCREEP*OPYCE*EPITP(NDIVI+NDIVS+4+J)*DF +
     &				SWELLTOPYC*DF
		EPIRCP(NDIVI+NDIVS+4+J) = EPIRCP(NDIVI+NDIVS+4+J) +
     &				SWELLROPYC*DF
		URCP(NDIVI+NDIVS+4+J) = 0.0D0
		DO 3139 K = 1, J -1
		  URCP(NDIVI+NDIVS+4+J) = URCP(NDIVI+NDIVS+4+J) +
     &				EPIRCP(NDIVI+NDIVS+4+K)*DR
3139		CONTINUE
		SIGT(NDIVI+NDIVS+4+J) = OPYCE*EPITP(NDIVI+NDIVS+4+J)
		SIGR(NDIVI+NDIVS+4+J) = -PAMB
		EPIR(NDIVI+NDIVS+4+J) = EPIRCP(NDIVI+NDIVS+4+J) - 
     &				OPYCNU*EPITP(NDIVI+NDIVS+4+J)
		EPIT(NDIVI+NDIVS+4+J) = EPITCP(NDIVI+NDIVS+4+J) + 
     &				EPITP(NDIVI+NDIVS+4+J)
		UR(NDIVI+NDIVS+4+J) = URCP(NDIVI+NDIVS+4+J) +
     &				URP(NDIVI+NDIVS+4+J)
3118	  CONTINUE
C  Calculate the strain energy stored in cracked OPyC layer due to shear force at the interface.
	  TEMP = 0.0D0
	  DO 3120 J = 1, NDIVO + 1
		R = (J-1)*DR + R4
		TEMP = TEMP + (EPITP(NDIVI+NDIVS+4+J)*R)**2.0*DR
3120	  CONTINUE
	  TEMP = TEMP*OPYCE/X2
	  WSHEAROPYC = (X1+OPYCNU)*KIOPYC**2.0*AOPYC*
     &			   (R4*(0.25D0*PIE-X1/X3)-0.5D0*AOPYC*
     &				(0.0625D0*PIE-X1/X3))/(16.0D0*OPYCE)
	  SHEAROPYC = X2*(TEMP+WSHEAROPYC/X2/PIE)/
     &			  (-PIE*R4**2.0*EPITCP(NDIVI+NDIVS+5))
C
C  2) Calculate the stresses in the intact SiC and IPyC layers
C
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments here.
	  R1 = R1/SCALE_R
	  R2 = R2/SCALE_R
	  R3 = R3/SCALE_R
	  R4 = R4/SCALE_R
	  R5 = R5/SCALE_R
	  DO 2905 K = 1,NDIV
	    UR(K) = UR(K)/SCALE_R
2905	  CONTINUE
	  IPYCE = IPYCE/SCALE_E
	  OPYCE = OPYCE/SCALE_E
	  SICE = SICE/SCALE_E
	  PRESS = PRESS/SCALE_E
	  PAMB = PAMB/SCALE_E
	  FLU = FLU/SCALE_F
	  IPYCREEP = IPYCREEP*SCALE_F*SCALE_E
	  OPYCREEP = OPYCREEP*SCALE_F*SCALE_E
	  DO 2906 I=0,NDEG
	    ISWR(I) = ISWR(I)*SCALE_F**(I+1)
	    ISWT(I) = ISWT(I)*SCALE_F**(I+1)
	    OSWR(I) = OSWR(I)*SCALE_F**(I+1)
	    OSWT(I) = OSWT(I)*SCALE_F**(I+1)
2906    CONTINUE
	  DO 2907 K=1,NDIV
	    SIGR(K) = SIGR(K)/SCALE_E
	    SIGT(K) = SIGT(K)/SCALE_E
2907    CONTINUE
C
C  Prepare quantities based on material properties
C
        L1I = X1/IPYCE
        L2I = X2*IPYCNU/IPYCE
        L3I = IPYCNU/IPYCE
        L4I = (X1 - IPYCNU)/IPYCE
C
	  G1I = X1/(L4I - L3I)
	  G2I = X1/(L1I + L3I)
	  G3I = L4I/(L1I + L3I)
	  G1S = SICE/(X1 - X2*SICNU)
	  G2S = X2*SICE/(X1 + SICNU)
C
	  V1I = X1 - X2*IPYCCNU
	  V2I = X2*(X1 + IPYCCNU)/X3
	  V3I = X2*(L1I*IPYCCNU-L3I)*G1I
C
C  Prepare fluence, pressure and creep quantities
	  FLUI = 1.0D+00
C    Stresses are additive; boundary conditions should exclude residual stresses
        PI = -PRESS-SIGR(1)
	  PO = -PAMB-SIGR(NDIVI+NDIVS+4)
	  CI = IPYCREEP
C
C  Initialize arrays P and Q
        DO 2908 I=1,NDIM2
	    Q(I)=0.0D0
2908    CONTINUE
C  Calculate initial stresses, strains, and displacement, i.e., I = 0
C    Step 1: Assign values to P(LDA,LDA). It's fixed for all orders of I.
        P(1,1) = 1.0D0
	  P(1,2) = -X2*G2I*(R2**(-3.0))/G1I
	  P(1,3) = 0.0D0
	  P(1,4) = 0.0D0
	  P(2,1) = 1.0D0
	  P(2,2) = -X2*G2I*(R3**(-3.0))/G1I
        P(2,3) = -G1S/G1I
	  P(2,4) = G2S*(R3**(-3.0))/G1I
	  P(3,1) = 1.0D0
	  P(3,2) = R3**(-3.0)
	  P(3,3) = -1.0D0
	  P(3,4) = -R3**(-3.0)
        P(4,1) = 0.0D0
        P(4,2) = 0.0D0
	  P(4,3) = 1.0D0
	  P(4,4) = -G2S*(R4**(-3.0))/G1S
C
C    Step 2: Assign values to Q(NDIM) for i=0
        Q(1) = PI/G1I
        Q(4) = PO/G1S
C    Step 3: Use DLSARG of IMSL library to calculate X
        CALL DLSARG (NDIM2, P(1:NDIM2,1:NDIM2), NDIM2, Q(1:NDIM2), 
     &			   IPATH, X(1:NDIM2))
        A(0) = X(1)
	  B(0) = X(2)
	  D(0) = X(3)
	  F(0) = X(4)
C
C    Register initial stresses and strains
        DO 2920 J = 1, NDIVI+2
	    R = (J-1)*(R3-R2)/FLOAT(NDIVI+1) + R2
	    SIGIR(J) = G1I*A(0) - X2*G2I*B(0)*(R**(-3.0))
	    SIGIT(J) = SIGIR(J) + X3*G2I*B(0)*(R**(-3.0))
	    EPIIR(J) = A(0) - X2*B(0)*(R**(-3.0))
	    EPIIT(J) = A(0) + B(0)*(R**(-3.0))
	    UIR(J) = A(0)*R + B(0)*(R**(-2.0))
2920    CONTINUE
        DO 2930 J = 1, NDIVS+2
	    R = (J-1)*(R4-R3)/FLOAT(NDIVS+1) + R3
	    SIGIR(NDIVI+2+J) = G1S*D(0) - G2S*F(0)*(R**(-3.0))
	    SIGIT(NDIVI+2+J) = SIGIR(NDIVI+2+J) + 
     &                       X3*G2S*F(0)*(R**(-3.0))/X2
	    EPIIR(NDIVI+2+J) = D(0) - X2*F(0)*(R**(-3.0))
	    EPIIT(NDIVI+2+J) = D(0) + F(0)*(R**(-3.0))
	    UIR(NDIVI+2+J) = D(0)*R + F(0)*(R**(-2.0))
2930    CONTINUE
C
C    Add to residual quantities from last fuel cycle
	  DO 2950 J = 1, NDIVI + NDIVS + 4
	    SIGR(J) = SIGIR(J) + SIGR(J)
          SIGT(J) = SIGIT(J) + SIGT(J)
	    EPIR(J) = EPIIR(J) + EPIR(J)
	    EPIT(J) = EPIIT(J) + EPIT(J)
	    UR(J)   = UIR(J)   + UR(J)
2950    CONTINUE
C
        IF(FLU.EQ.0.0D0) THEN
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments at the beginning of this subroutine and 
C  change them back here before returning to the main program.
	    R1 = R1*SCALE_R
	    R2 = R2*SCALE_R
	    R3 = R3*SCALE_R
	    R4 = R4*SCALE_R
	    R5 = R5*SCALE_R
	    DO 2909 K = 1,NDIV
	      UR(K) = UR(K)*SCALE_R
2909	    CONTINUE
	    IPYCE = IPYCE*SCALE_E
	    OPYCE = OPYCE*SCALE_E
	    SICE = SICE*SCALE_E
	    PRESS = PRESS*SCALE_E
	    PAMB = PAMB*SCALE_E
	    FLU = FLU*SCALE_F
	    IPYCREEP = IPYCREEP/SCALE_F/SCALE_E
	    OPYCREEP = OPYCREEP/SCALE_F/SCALE_E
	    DO 2910 K=0,NDEG
	      ISWR(K) = ISWR(K)/(SCALE_F**(K+1))
	      ISWT(K) = ISWT(K)/(SCALE_F**(K+1))
	      OSWR(K) = OSWR(K)/(SCALE_F**(K+1))
	      OSWT(K) = OSWT(K)/(SCALE_F**(K+1))
2910      CONTINUE
	    DO 2911 K=1,NDIV
	      SIGR(K) = SIGR(K)*SCALE_E
	      SIGT(K) = SIGT(K)*SCALE_E
2911      CONTINUE
	    RETURN  !no need for calculating higher order terms
	  END IF
C    Step 4: Prepare stress entries for higher order calculations
        DO 2960 J=1,NDIVI+NDIVS+4
	    SIGIR(J) = SIGR(J)
	    SIGIT(J) = SIGT(J)
2960    CONTINUE
        SIGR2 = SIGIR(1)
	  SIGT2 = SIGIT(1)
        SIGR3 = SIGIR(NDIVI+2)
	  SIGT3 = SIGIT(NDIVI+2)
	  SIG2 = (X1+L3I*G1I*V1I)*SIGR2 + (L1I*G1I*V1I-X1)*SIGT2
     	  SIG3 = (X1+L3I*G1I*V1I)*SIGR3 + (L1I*G1I*V1I-X1)*SIGT3
C
C  Begin to calculate higher order terms
	  GI = 0.0
C  Order loop
        CONV = .FALSE.
        DO 2990 I=0, IORDER
C  Initialize probes for testing convergence
          PROBER=0.
	    PROBET=0.
C
	    FLUI = FLUI*FLU
C
          IF (I .LE. NDEG) THEN
	      SI = (L4I*ISWR(I)+L2I*ISWT(I))/((I+1.0)*(L1I*L4I-L2I*L3I))
	      SI1 = G2I*(ISWR(I)-ISWT(I))/(I+1.0)
	      FI = (-CI)*V3I*GI/(X2*L4I*(I+1.0))+(ISWR(I)-ISWT(I))/
     &           (L4I*(I+1.0))
	      GI = FI - CI*G1I*V1I*GI/(I+1.0)
	    ELSE
	      SI = 0.0
	      SI1 = 0.0
	      FI = (-CI)*V3I*GI/(X2*L4I*(I+1.0))
	      GI = FI - CI*G1I*V1I*GI/(I+1.0)
	    END IF
C
C  Calculate A(I+1), B(I+1), D(I+1), F(I+1)
C    Step 1: Assign new values to Q(NDIM2)
	    Q(1) = (-X2*FI*(DLOG(R2)+G3I)/X3+CI*G2I*SIG2+SI)/G1I
	    Q(2) = (-X2*FI*(DLOG(R3)+G3I)/X3+CI*G2I*SIG3+SI)/G1I
	    Q(3) = -X2*FI*DLOG(R3)/(X3*G1I)
	    Q(4) = 0.0D0
C    Step 2: Use DLSARG of IMSL library to calculate X
          CALL DLSARG (NDIM2, P(1:NDIM2,1:NDIM2), NDIM2, Q(1:NDIM2), 
     &				 IPATH, X(1:NDIM2))
          A(I+1) = X(1)
	    B(I+1) = X(2)
	    D(I+1) = X(3)
	    F(I+1) = X(4)
C
C  Calculate stresses, strains, and displacement of order i+1
          DO 2970 J = 1, NDIVI+2
	      R = (J-1)*(R3-R2)/FLOAT(NDIVI+1) + R2
	      SIGINTER = SIGIR(J)
	      SIGIR(J) = G1I*A(I+1) - X2*G2I*B(I+1)*(R**(-3.0)) + X2*FI*
     A               (DLOG(R)+G3I)/X3 - CI*G2I*((X1+L3I*G1I*V1I)*
     B               SIGIR(J)+(L1I*G1I*V1I-X1)*SIGIT(J))/(I+1.0)- SI
	      SIGIT(J) = SIGIR(J) + X3*G2I*B(I+1)*(R**(-3.0)) - X2*FI*
     A               G2I/(X3*G1I) + X3*CI*G2I*V2I*(SIGINTER-SIGIT(J))/
     B               (X2*(I+1.0)) + SI1
	      EPIIR(J) = A(I+1) - X2*B(I+1)*(R**(-3.0)) + 
     A	  			  X2*(L4I-L3I)*FI*(DLOG(R)+X1)/X3
	      EPIIT(J) = A(I+1) + B(I+1)*(R**(-3.0)) +
     A				  X2*(L4I-L3I)*FI*DLOG(R)/X3
	      UIR(J)   = A(I+1)*R + B(I+1)*(R**(-2.0)) +
     A				  X2*(L4I-L3I)*FI*R*DLOG(R)/X3
2970      CONTINUE
          DO 2980 J = 1, NDIVS+2
	      R = (J-1)*(R4-R3)/FLOAT(NDIVS+1) + R3
	      SIGINTER = SIGIR(NDIVI+2+J)
	      SIGIR(NDIVI+2+J) = G1S*D(I+1) - G2S*F(I+1)*(R**(-3.0))
	      SIGIT(NDIVI+2+J) = SIGIR(NDIVI+2+J) + X3*G2S*F(I+1)*
     &                         (R**(-3.0))/X2
	      EPIIR(NDIVI+2+J) = D(I+1) - X2*F(I+1)*(R**(-3.0))
	      EPIIT(NDIVI+2+J) = D(I+1) + F(I+1)*(R**(-3.0))
	      UIR(NDIVI+2+J)   = D(I+1)*R + F(I+1)*(R**(-2.0))
2980      CONTINUE
C
	    DO 3000 J=1, NDIVI+NDIVS+4
	      SIGR(J) = SIGR(J) + SIGIR(J)*FLUI
	      SIGT(J) = SIGT(J) + SIGIT(J)*FLUI
		  EPIR(J) = EPIR(J) + EPIIR(J)*FLUI
	      EPIT(J) = EPIT(J) + EPIIT(J)*FLUI
		  UR(J)   = UR(J)   + UIR(J)*FLUI
C
C  Calculate the average-square-root of the term of order i
            PROBER = PROBER + (SIGIR(J)*FLUI/SIGR(J))**2
            PROBET = PROBET + (SIGIT(J)*FLUI/SIGT(J))**2
3000      CONTINUE
          PROBER = DSQRT(PROBER/(NDIVI+NDIVS+4))
          PROBET = DSQRT(PROBET/(NDIVI+NDIVS+4))
	    IF((PROBER.LT.CRITERIA).AND.(PROBET.LT.CRITERIA)) THEN
	      IF(CONV.EQ..FALSE.) THEN
	        CONV = .TRUE.
	      ELSE
C	        WRITE(*,*) 'Convergence achieved at order i = ', I
C
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments at the beginning of this subroutine and 
C  change them back here before returning to the main program.
	        R1 = R1*SCALE_R
	        R2 = R2*SCALE_R
	        R3 = R3*SCALE_R
	        R4 = R4*SCALE_R
	        R5 = R5*SCALE_R
		    DO 2912 K = 1,NDIV
			  UR(K) = UR(K)*SCALE_R
2912		    CONTINUE
	        IPYCE = IPYCE*SCALE_E
	        OPYCE = OPYCE*SCALE_E
	        SICE = SICE*SCALE_E
	        PRESS = PRESS*SCALE_E
	        PAMB = PAMB*SCALE_E
	        FLU = FLU*SCALE_F
	        IPYCREEP = IPYCREEP/SCALE_F/SCALE_E
	        OPYCREEP = OPYCREEP/SCALE_F/SCALE_E
	        DO 2913 K=0,NDEG
	          ISWR(K) = ISWR(K)/(SCALE_F**(K+1))
	          ISWT(K) = ISWT(K)/(SCALE_F**(K+1))
	          OSWR(K) = OSWR(K)/(SCALE_F**(K+1)) 
	          OSWT(K) = OSWT(K)/(SCALE_F**(K+1))
2913          CONTINUE
	        DO 2914 K=1,NDIV
	          SIGR(K) = SIGR(K)*SCALE_E
	          SIGT(K) = SIGT(K)*SCALE_E
2914          CONTINUE
	        RETURN
	      ENDIF
	    ELSE IF(CONV.EQ..TRUE.) THEN
C           WRITE(*,*) 'At Order = ',I, 'abormal'
	      CONV = .FALSE.
	    ENDIF
C
C  Update stress entries for higher orders
          SIGINTER = SIGR2
          SIGR2 = G1I*A(I+1) - X2*G2I*B(I+1)*(R2**(-3.0)) + X2*FI*
     A           (DLOG(R2)+G3I)/X3 - CI*G2I*((X1+L3I*G1I*V1I)*SIGR2 +
     B           (L1I*G1I*V1I-X1)*SIGT2)/(I+1.0) - SI
          SIGT2 = SIGR2 + X3*G2I*B(I+1)*(R2**(-3.0)) - X2*FI*
     A           G2I/(X3*G1I) + X3*CI*G2I*V2I*(SIGINTER-SIGT2)/
     B           (X2*(I+1.0)) + SI1
          SIGINTER = SIGR3
          SIGR3 = G1I*A(I+1) - X2*G2I*B(I+1)*(R3**(-3.0)) + X2*FI*
     A           (DLOG(R3)+G3I)/X3 - CI*G2I*((X1+L3I*G1I*V1I)*SIGR3 +
     B           (L1I*G1I*V1I-X1)*SIGT3)/(I+1.0) - SI
          SIGT3 = SIGR3 + X3*G2I*B(I+1)*(R3**(-3.0)) - X2*FI*
     A           G2I/(X3*G1I) + X3*CI*G2I*V2I*(SIGINTER-SIGT3)/
     B           (X2*(I+1.0)) + SI1
	    SIG2 = ((X1+L3I*G1I*V1I)*SIGR2+(L1I*G1I*V1I-X1)*SIGT2)/(I+2.0)
	    SIG3 = ((X1+L3I*G1I*V1I)*SIGR3+(L1I*G1I*V1I-X1)*SIGT3)/(I+2.0)
2990    CONTINUE
        IF (I.GE.IORDER) WRITE(*,*) 'Stress has not converge at IORDER.'
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments at the beginning of this subroutine and 
C  change them back here before returning to the main program.
	  R1 = R1*SCALE_R
	  R2 = R2*SCALE_R
	  R3 = R3*SCALE_R
	  R4 = R4*SCALE_R
	  R5 = R5*SCALE_R
	  DO 2915 K = 1,NDIV
	    UR(K) = UR(K)*SCALE_R
2915	  CONTINUE
	  IPYCE = IPYCE*SCALE_E
	  OPYCE = OPYCE*SCALE_E
	  SICE = SICE*SCALE_E
	  PRESS = PRESS*SCALE_E
	  PAMB = PAMB*SCALE_E
	  FLU = FLU*SCALE_F
	  IPYCREEP = IPYCREEP/SCALE_F/SCALE_E
	  OPYCREEP = OPYCREEP/SCALE_F/SCALE_E
	  DO 2916 K=0,NDEG
	    ISWR(K) = ISWR(K)/(SCALE_F**(K+1))
	    ISWT(K) = ISWT(K)/(SCALE_F**(K+1))
	    OSWR(K) = OSWR(K)/(SCALE_F**(K+1))
	    OSWT(K) = OSWT(K)/(SCALE_F**(K+1))
2916    CONTINUE
	  DO 2917 K=1,NDIV
	    SIGR(K) = SIGR(K)*SCALE_E
	    SIGT(K) = SIGT(K)*SCALE_E
2917    CONTINUE
C
	ELSE IF(MCODE .EQ. 'SO2') THEN
C  1) First evaluate the stresses, strains, and displacements in the cracked IPyC layer,
C  and calculate the shear force it applies to the SiC layer.
	  AIPYC = 0.99D0*(R3-R2)
	  SWELLRIPYC = 0.0D0
	  SWELLTIPYC = 0.0D0
	  DO 3111 K=0,NDEG
		SWELLRIPYC = SWELLRIPYC + ISWR(K)*FLU**K
		SWELLTIPYC = SWELLTIPYC + ISWT(K)*FLU**K
3111	  CONTINUE
	  DO 3131 K = 1, NDIVI + 2
		URP(K) = 0.0D0
3131	  CONTINUE
	  DR = (R3-R2)/FLOAT(NDIVI+1)
	  CONV = .FALSE.
	  DO WHILE (.NOT.CONV)
		TEMP = KIIPYC
		SIGTPIBAR = 0.0D0
		DO 3112 J = 1, NDIVI + 2
		  R = (J-1)*DR + R2
		  EPITP(J)=((URCP(J)+URP(J))/R-EPITCP(J)-SWELLTIPYC*DF-8.0D0*
     &		(X1-IPYCNU**2.0)*KIIPYC*SQRT(ABS(R3-R)/(X2*PIE*R**2.0))
     &		/PIE/IPYCE)/(X1+IPYCE*IPYCREEP*DF)
		  SIGTPIBAR = SIGTPIBAR + IPYCE*EPITP(J)
3112		CONTINUE
		SIGTPIBAR = SIGTPIBAR/FLOAT(NDIVI+2)
	    KIIPYC = 0.413*(X1+R2/R3)*SIGTPIBAR*SQRT(PIE*AIPYC)/
     &			 SQRT(X1 - AIPYC/(R3-R2))
C  Approaches the right value of KIIPYC
		KIIPYC = (KIIPYC + TEMP)/X2
	    IF(ABS((TEMP-KIIPYC)/KIIPYC).LE.CRITERIA) CONV = .TRUE.
C  Recalculate URP
		DO 3132 J = 1, NDIVI + 2
		  URP(J) = 0.0D0
		  DO 3133 K = J, NDIVI + 1
			URP(J) = URP(J) + IPYCNU*EPITP(K)*DR
3133		  CONTINUE
3132		CONTINUE
	  END DO
C  Update EPITCP, and record resulting stresses, strains, and displacements in the cracked IPyC layer
	  DO 3113 J = 1, NDIVI + 2
		EPITCP(J) = EPITCP(J) + IPYCREEP*IPYCE*EPITP(J)*DF +
     &				SWELLTIPYC*DF
		EPIRCP(J) = EPIRCP(J) + SWELLRIPYC*DF
		URCP(J) = 0.0D0
		DO 3134 K = J, NDIVI + 1
		  URCP(J) = URCP(J) - EPIRCP(K)*DR
3134		CONTINUE
		SIGT(J) = IPYCE*EPITP(J)
		SIGR(J) = -PRESS
		EPIR(J) = EPIRCP(J) - IPYCNU*EPITP(J)
		EPIT(J) = EPITCP(J) + EPITP(J)
		UR(J) = URCP(J) + URP(J)
3113	  CONTINUE
C  Calculate the strain energy stored in cracked IPyC layer due to shear force at the interface.
	  TEMP = 0.0D0
	  DO 3115 J = 1, NDIVI + 1
		R = (J-1)*DR + R2
		TEMP = TEMP + (EPITP(J)*R)**2.0*DR
3115	  CONTINUE
	  TEMP = TEMP*IPYCE/X2
	  WSHEARIPYC = (X1+IPYCNU)*KIIPYC**2.0*AIPYC*
     &			   (R3*(0.25D0*PIE-X1/X3)+0.5D0*AIPYC*
     &				(0.0625D0*PIE-X1/X3))/(16.0D0*IPYCE)
	  SHEARIPYC = X2*(TEMP+WSHEARIPYC/X2/PIE)/
     &			  (-PIE*R3**2.0*EPITCP(NDIVI+2))
C
C  2) Calculate the stresses in the intact SiC and OPyC layers
C
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments here.
	  R1 = R1/SCALE_R
	  R2 = R2/SCALE_R
	  R3 = R3/SCALE_R
	  R4 = R4/SCALE_R
	  R5 = R5/SCALE_R
	  DO 3005 K = 1,NDIV
	    UR(K) = UR(K)/SCALE_R
3005	  CONTINUE
	  IPYCE = IPYCE/SCALE_E
	  OPYCE = OPYCE/SCALE_E
	  SICE = SICE/SCALE_E
	  PRESS = PRESS/SCALE_E
	  PAMB = PAMB/SCALE_E
	  FLU = FLU/SCALE_F
	  IPYCREEP = IPYCREEP*SCALE_F*SCALE_E
	  OPYCREEP = OPYCREEP*SCALE_F*SCALE_E
	  DO 3006 I=0,NDEG
	    ISWR(I) = ISWR(I)*SCALE_F**(I+1)
	    ISWT(I) = ISWT(I)*SCALE_F**(I+1)
	    OSWR(I) = OSWR(I)*SCALE_F**(I+1)
	    OSWT(I) = OSWT(I)*SCALE_F**(I+1)
3006    CONTINUE
	  DO 3007 K=1,NDIV
	    SIGR(K) = SIGR(K)/SCALE_E
	    SIGT(K) = SIGT(K)/SCALE_E
3007    CONTINUE
C
C  Prepare quantities based on material properties
C
        L1O = X1/OPYCE
        L2O = X2*OPYCNU/OPYCE
        L3O = OPYCNU/OPYCE
        L4O = (X1 - OPYCNU)/OPYCE
C
	  G1O = X1/(L4O - L3O)
	  G2O = X1/(L1O + L3O)
	  G3O = L4O/(L1O + L3O)
	  G1S = SICE/(X1 - X2*SICNU)
	  G2S = X2*SICE/(X1 + SICNU)
C
	  V1O = X1 - X2*OPYCCNU
	  V2O = X2*(X1 + OPYCCNU)/X3
	  V3O = X2*(L1O*OPYCCNU-L3O)*G1O
C
C  Prepare fluence, pressure and creep quantities
	  FLUI = 1.0D+00
C    Stresses are additive; boundary conditions should exclude residual stresses
        PI = -PRESS-SIGR(NDIVI+3)
	  PO = -PAMB-SIGR(NDIV)
	  CO = OPYCREEP
C
C  Initialize arrays P and Q
        DO 3008 I=1,NDIM2
	    Q(I)=0.0D0
3008    CONTINUE
C  Calculate initial stresses, strains, and displacement, i.e., I = 0
C    Step 1: Assign values to P(LDA,LDA). It's fixed for all orders of I.
        P(1,1) = 1.0D0
	  P(1,2) = -G2S*(R3**(-3.0))/G1S
	  P(1,3) = 0.0D0
	  P(1,4) = 0.0D0
	  P(2,1) = -G1S/G1O
	  P(2,2) = G2S*(R4**(-3.0))/G1O
        P(2,3) = 1.0D0
	  P(2,4) = -X2*G2O*(R4**(-3.0))/G1O
	  P(3,1) = -1.0D0
	  P(3,2) = -R4**(-3.0)
	  P(3,3) = 1.0D0
	  P(3,4) = R4**(-3.0)
        P(4,1) = 0.0D0
        P(4,2) = 0.0D0
	  P(4,3) = 1.0D0
	  P(4,4) = -X2*G2O*(R5**(-3.0))/G1O
C
C    Step 2: Assign values to Q(NDIM) for i=0
        Q(1) = PI/G1S
        Q(4) = PO/G1O
C    Step 3: Use DLSARG of IMSL library to calculate X
        CALL DLSARG (NDIM2, P(1:NDIM2,1:NDIM2), NDIM2, Q(1:NDIM2), 
     &			   IPATH, X(1:NDIM2))
        D(0) = X(1)
	  F(0) = X(2)
	  M(0) = X(3)
	  N(0) = X(4)
C
C    Register initial stresses and strains
        DO 3030 J = 1, NDIVS+2
	    R = (J-1)*(R4-R3)/FLOAT(NDIVS+1) + R3
	    SIGIR(NDIVI+2+J) = G1S*D(0) - G2S*F(0)*(R**(-3.0))
	    SIGIT(NDIVI+2+J) = SIGIR(NDIVI+2+J) + 
     &                       X3*G2S*F(0)*(R**(-3.0))/X2
	    EPIIR(NDIVI+2+J) = D(0) - X2*F(0)*(R**(-3.0))
	    EPIIT(NDIVI+2+J) = D(0) + F(0)*(R**(-3.0))
	    UIR(NDIVI+2+J) = D(0)*R + F(0)*(R**(-2.0))
3030    CONTINUE
        DO 3040 J = 1, NDIVO+2
	    R = (J-1)*(R5-R4)/FLOAT(NDIVO+1) + R4
	    SIGIR(NDIVI+NDIVS+4+J) = G1O*M(0) - X2*G2O*N(0)*(R**(-3.0))
	    SIGIT(NDIVI+NDIVS+4+J) = SIGIR(NDIVI+NDIVS+4+J) + 
     &                             X3*G2O*N(0)*(R**(-3.0))
	    EPIIR(NDIVI+NDIVS+4+J) = M(0) - X2*N(0)*(R**(-3.0))
	    EPIIT(NDIVI+NDIVS+4+J) = M(0) + N(0)*(R**(-3.0))
	    UIR(NDIVI+NDIVS+4+J) = M(0)*R + N(0)*(R**(-2.0))
3040    CONTINUE
C
C    Add to residual quantities from last fuel cycle
	  DO 3050 J = NDIVI+3, NDIV
	    SIGR(J) = SIGIR(J) + SIGR(J)
          SIGT(J) = SIGIT(J) + SIGT(J)
	    EPIR(J) = EPIIR(J) + EPIR(J)
	    EPIT(J) = EPIIT(J) + EPIT(J)
	    UR(J)   = UIR(J)   + UR(J)
3050    CONTINUE
C
        IF(FLU.EQ.0.0D0) THEN
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments at the beginning of this subroutine and 
C  change them back here before returning to the main program.
	    R1 = R1*SCALE_R
	    R2 = R2*SCALE_R
	    R3 = R3*SCALE_R
	    R4 = R4*SCALE_R
	    R5 = R5*SCALE_R
	    DO 3009 K = 1,NDIV
	      UR(K) = UR(K)*SCALE_R
3009	    CONTINUE
	    IPYCE = IPYCE*SCALE_E
	    OPYCE = OPYCE*SCALE_E
	    SICE = SICE*SCALE_E
	    PRESS = PRESS*SCALE_E
	    PAMB = PAMB*SCALE_E
	    FLU = FLU*SCALE_F
	    IPYCREEP = IPYCREEP/SCALE_F/SCALE_E
	    OPYCREEP = OPYCREEP/SCALE_F/SCALE_E
	    DO 3010 K=0,NDEG
	      ISWR(K) = ISWR(K)/(SCALE_F**(K+1))
	      ISWT(K) = ISWT(K)/(SCALE_F**(K+1))
	      OSWR(K) = OSWR(K)/(SCALE_F**(K+1))
	      OSWT(K) = OSWT(K)/(SCALE_F**(K+1))
3010      CONTINUE
	    DO 3011 K=1,NDIV
	      SIGR(K) = SIGR(K)*SCALE_E
	      SIGT(K) = SIGT(K)*SCALE_E
3011      CONTINUE
	    RETURN  !no need for calculating higher order terms
	  END IF
C    Step 4: Prepare stress entries for higher order calculations
        DO 3060 J = NDIVI+3, NDIV
	    SIGIR(J) = SIGR(J)
	    SIGIT(J) = SIGT(J)
3060    CONTINUE
	  SIGR4 = SIGIR(NDIVI+NDIVS+5)
	  SIGT4 = SIGIT(NDIVI+NDIVS+5)
	  SIGR5 = SIGIR(NDIV)
	  SIGT5 = SIGIT(NDIV)
        SIG4 = (X1+L3O*G1O*V1O)*SIGR4 + (L1O*G1O*V1O-X1)*SIGT4
        SIG5 = (X1+L3O*G1O*V1O)*SIGR5 + (L1O*G1O*V1O-X1)*SIGT5
C
C  Begin to calculate higher order terms
	  GO = 0.0
C  Order loop
        CONV = .FALSE.
        DO 3190 I=0, IORDER
C  Initialize probes for testing convergence
          PROBER=0.
	    PROBET=0.
C
	    FLUI = FLUI*FLU
C
          IF (I .LE. NDEG) THEN
	      SO = (L4O*OSWR(I)+L2O*OSWT(I))/((I+1.0)*(L1O*L4O-L2O*L3O))
	      SO1 = G2O*(OSWR(I)-OSWT(I))/(I+1.0)
	      FO = (-CO)*V3O*GO/(X2*L4O*(I+1.0))+(OSWR(I)-OSWT(I))/
     &           (L4O*(I+1.0))
	      GO = FO - CO*G1O*V1O*GO/(I+1.0)
	    ELSE
	      SO = 0.0
	      SO1 = 0.0
	      FO = (-CO)*V3O*GO/(X2*L4O*(I+1.0))
	      GO = FO - CO*G1O*V1O*GO/(I+1.0)
	    END IF
C
C  Calculate A(I+1), B(I+1), D(I+1), F(I+1)
C    Step 1: Assign new values to Q(NDIM2)
	    Q(1) = 0.0D0
	    Q(2) = (-X2*FO*(DLOG(R4)+G3O)/X3+CO*G2O*SIG4+SO)/G1O
	    Q(3) = -X2*FO*DLOG(R4)/(X3*G1O)
	    Q(4) = (-X2*FO*(DLOG(R5)+G3O)/X3+CO*G2O*SIG5+SO)/G1O
C    Step 2: Use DLSARG of IMSL library to calculate X
          CALL DLSARG (NDIM2, P(1:NDIM2,1:NDIM2), NDIM2, Q(1:NDIM2), 
     &				 IPATH, X(1:NDIM2))
          D(I+1) = X(1)
	    F(I+1) = X(2)
	    M(I+1) = X(3)
	    N(I+1) = X(4)
C
C  Calculate stresses, strains, and displacement of order i+1
          DO 3080 J = 1, NDIVS+2
	      R = (J-1)*(R4-R3)/FLOAT(NDIVS+1) + R3
	      SIGINTER = SIGIR(NDIVI+2+J)
	      SIGIR(NDIVI+2+J) = G1S*D(I+1) - G2S*F(I+1)*(R**(-3.0))
	      SIGIT(NDIVI+2+J) = SIGIR(NDIVI+2+J) + X3*G2S*F(I+1)*
     &                         (R**(-3.0))/X2
	      EPIIR(NDIVI+2+J) = D(I+1) - X2*F(I+1)*(R**(-3.0))
	      EPIIT(NDIVI+2+J) = D(I+1) + F(I+1)*(R**(-3.0))
	      UIR(NDIVI+2+J)   = D(I+1)*R + F(I+1)*(R**(-2.0))
3080      CONTINUE
	    DO 3090 J = 1, NDIVO+2
	      R = (J-1)*(R5-R4)/FLOAT(NDIVO+1) + R4
	      SIGINTER = SIGIR(NDIVI+NDIVS+4+J)
	      SIGIR(NDIVI+NDIVS+4+J) = G1O*M(I+1) - X2*G2O*N(I+1)*
     A             (R**(-3.0)) + X2*FO*(DLOG(R)+G3O)/X3 - CO*G2O*
     B             ((X1+L3O*G1O*V1O)*SIGIR(NDIVI+NDIVS+4+J)+
     C             (L1O*G1O*V1O-X1)*SIGIT(NDIVI+NDIVS+4+J))/(I+1.0) - SO
	      SIGIT(NDIVI+NDIVS+4+J) = SIGIR(NDIVI+NDIVS+4+J) + 
     A             X3*G2O*N(I+1)*(R**(-3.0)) - X2*FO*G2O/(X3*G1O) + 
     B             X3*CO*G2O*V2O*(SIGINTER-SIGIT(NDIVI+NDIVS+4+J))/
     C             (X2*(I+1.0)) + SO1
	      EPIIR(NDIVI+NDIVS+4+J) = M(I+1) - X2*N(I+1)*(R**(-3.0)) + 
     A							  X2*(L4O-L3O)*FO*(DLOG(R)+X1)/X3
	      EPIIT(NDIVI+NDIVS+4+J) = M(I+1) + N(I+1)*(R**(-3.0)) +
     A							  X2*(L4O-L3O)*FO*DLOG(R)/X3
	      UIR(NDIVI+NDIVS+4+J)   = M(I+1)*R + N(I+1)*(R**(-2.0)) +
     A							  X2*(L4O-L3O)*FO*R*DLOG(R)/X3
3090      CONTINUE
C
	    DO 3100 J = NDIVI+3, NDIV
	      SIGR(J) = SIGR(J) + SIGIR(J)*FLUI
	      SIGT(J) = SIGT(J) + SIGIT(J)*FLUI
		  EPIR(J) = EPIR(J) + EPIIR(J)*FLUI
	      EPIT(J) = EPIT(J) + EPIIT(J)*FLUI
		  UR(J)   = UR(J)   + UIR(J)*FLUI
C
C  Calculate the average-square-root of the term of order i
            PROBER = PROBER + (SIGIR(J)*FLUI/SIGR(J))**2
            PROBET = PROBET + (SIGIT(J)*FLUI/SIGT(J))**2
3100      CONTINUE
          PROBER = DSQRT(PROBER/(NDIV-NDIVI-2))
          PROBET = DSQRT(PROBET/(NDIV-NDIVI-2))
	    IF((PROBER.LT.CRITERIA).AND.(PROBET.LT.CRITERIA)) THEN
	      IF(CONV.EQ..FALSE.) THEN
	        CONV = .TRUE.
	      ELSE
C	        WRITE(*,*) 'Convergence achieved at order i = ', I
C
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments at the beginning of this subroutine and 
C  change them back here before returning to the main program.
	        R1 = R1*SCALE_R
	        R2 = R2*SCALE_R
	        R3 = R3*SCALE_R
	        R4 = R4*SCALE_R
	        R5 = R5*SCALE_R
		    DO 3012 K = 1,NDIV
			  UR(K) = UR(K)*SCALE_R
3012		    CONTINUE
	        IPYCE = IPYCE*SCALE_E
	        OPYCE = OPYCE*SCALE_E
	        SICE = SICE*SCALE_E
	        PRESS = PRESS*SCALE_E
	        PAMB = PAMB*SCALE_E
	        FLU = FLU*SCALE_F
	        IPYCREEP = IPYCREEP/SCALE_F/SCALE_E
	        OPYCREEP = OPYCREEP/SCALE_F/SCALE_E
	        DO 3013 K=0,NDEG
	          ISWR(K) = ISWR(K)/(SCALE_F**(K+1))
	          ISWT(K) = ISWT(K)/(SCALE_F**(K+1))
	          OSWR(K) = OSWR(K)/(SCALE_F**(K+1)) 
	          OSWT(K) = OSWT(K)/(SCALE_F**(K+1))
3013          CONTINUE
	        DO 3014 K=1,NDIV
	          SIGR(K) = SIGR(K)*SCALE_E
	          SIGT(K) = SIGT(K)*SCALE_E
3014          CONTINUE
	        RETURN
	      ENDIF
	    ELSE IF(CONV.EQ..TRUE.) THEN
C           WRITE(*,*) 'At Order = ',I, 'abormal'
	      CONV = .FALSE.
	    ENDIF
C
C  Update stress entries for higher orders
	    SIGINTER = SIGR4
	    SIGR4 = G1O*M(I+1) - X2*G2O*N(I+1)*(R4**(-3.0)) + X2*FO*
     A           (DLOG(R4)+G3O)/X3 - CO*G2O*((X1+L3O*G1O*V1O)*SIGR4 +
     B           (L1O*G1O*V1O-X1)*SIGT4)/(I+1.0) - SO
	    SIGT4 = SIGR4 + X3*G2O*N(I+1)*(R4**(-3.0)) - X2*FO*
     A           G2O/(X3*G1O) + X3*CO*G2O*V2O*(SIGINTER-SIGT4)/
     B           (X2*(I+1.0)) + SO1
	    SIGINTER = SIGR5
	    SIGR5 = G1O*M(I+1) - X2*G2O*N(I+1)*(R5**(-3.0)) + X2*FO*
     A           (DLOG(R5)+G3O)/X3 - CO*G2O*((X1+L3O*G1O*V1O)*SIGR5 +
     B           (L1O*G1O*V1O-X1)*SIGT5)/(I+1.0) - SO
	    SIGT5 = SIGR5 + X3*G2O*N(I+1)*(R5**(-3.0)) - X2*FO*
     A           G2O/(X3*G1O) + X3*CO*G2O*V2O*(SIGINTER-SIGT5)/
     B           (X2*(I+1.0)) + SO1
	    SIG4 = ((X1+L3O*G1O*V1O)*SIGR4+(L1O*G1O*V1O-X1)*SIGT4)/(I+2.0)
	    SIG5 = ((X1+L3O*G1O*V1O)*SIGR5+(L1O*G1O*V1O-X1)*SIGT5)/(I+2.0)
3190    CONTINUE
        IF (I.GE.IORDER) WRITE(*,*) 'Stress has not converge at IORDER.'
C  Because the calculation encouters error accumulation when some quantities are large,
C  we make some dimension adjustments at the beginning of this subroutine and 
C  change them back here before returning to the main program.
	  R1 = R1*SCALE_R
	  R2 = R2*SCALE_R
	  R3 = R3*SCALE_R
	  R4 = R4*SCALE_R
	  R5 = R5*SCALE_R
	  DO 3015 K = 1,NDIV
	    UR(K) = UR(K)*SCALE_R
3015	  CONTINUE
	  IPYCE = IPYCE*SCALE_E
	  OPYCE = OPYCE*SCALE_E
	  SICE = SICE*SCALE_E
	  PRESS = PRESS*SCALE_E
	  PAMB = PAMB*SCALE_E
	  FLU = FLU*SCALE_F
	  IPYCREEP = IPYCREEP/SCALE_F/SCALE_E
	  OPYCREEP = OPYCREEP/SCALE_F/SCALE_E
	  DO 3016 K=0,NDEG
	    ISWR(K) = ISWR(K)/(SCALE_F**(K+1))
	    ISWT(K) = ISWT(K)/(SCALE_F**(K+1))
	    OSWR(K) = OSWR(K)/(SCALE_F**(K+1))
	    OSWT(K) = OSWT(K)/(SCALE_F**(K+1))
3016    CONTINUE
	  DO 3017 K=1,NDIV
	    SIGR(K) = SIGR(K)*SCALE_E
	    SIGT(K) = SIGT(K)*SCALE_E
3017    CONTINUE
C
	ELSE IF(MCODE .EQ. 'S1') THEN
C  1) First evaluate the stresses, strains, and displacements in the cracked IPyC layer,
C  and calculate the shear force it applies to the SiC layer.
	  AIPYC = 0.99D0*(R3-R2)
	  SWELLRIPYC = 0.0D0
	  SWELLTIPYC = 0.0D0
	  DO 3121 K=0,NDEG
		SWELLRIPYC = SWELLRIPYC + ISWR(K)*FLU**K
		SWELLTIPYC = SWELLTIPYC + ISWT(K)*FLU**K
3121	  CONTINUE
	  DO 3141 K = 1, NDIVI + 2
		URP(K) = 0.0D0
3141	  CONTINUE
	  DR = (R3-R2)/FLOAT(NDIVI+1)
	  CONV = .FALSE.
	  DO WHILE (.NOT.CONV)
		TEMP = KIIPYC
		SIGTPIBAR = 0.0D0
		DO 3122 J = 1, NDIVI + 2
		  R = (J-1)*DR + R2
		  EPITP(J)=((URCP(J)+URP(J))/R-EPITCP(J)-SWELLTIPYC*DF-8.0D0*
     &		(X1-IPYCNU**2.0)*KIIPYC*SQRT(ABS(R3-R)/(X2*PIE*R**2.0))
     &		/PIE/IPYCE)/(X1+IPYCE*IPYCREEP*DF)
		  SIGTPIBAR = SIGTPIBAR + IPYCE*EPITP(J)
3122		CONTINUE
		SIGTPIBAR = SIGTPIBAR/FLOAT(NDIVI+2)
	    KIIPYC = 0.413*(X1+R2/R3)*SIGTPIBAR*SQRT(PIE*AIPYC)/
     &			 SQRT(X1 - AIPYC/(R3-R2))
C  Approaches the right value of KIIPYC
		KIIPYC = (KIIPYC + TEMP)/X2
	    IF(ABS((TEMP-KIIPYC)/KIIPYC).LE.CRITERIA) CONV = .TRUE.
C  Recalculate URP
		DO 3142 J = 1, NDIVI + 2
		  URP(J) = 0.0D0
		  DO 3143 K = J, NDIVI + 1
			URP(J) = URP(J) + IPYCNU*EPITP(K)*DR
3143		  CONTINUE
3142		CONTINUE
	  END DO
C  Update EPITCP, and record resulting stresses, strains, and displacements in the cracked IPyC layer
	  DO 3123 J = 1, NDIVI + 2
		EPITCP(J) = EPITCP(J) + IPYCREEP*IPYCE*EPITP(J)*DF +
     &				SWELLTIPYC*DF
		EPIRCP(J) = EPIRCP(J) + SWELLRIPYC*DF
		URCP(J) = 0.0D0
		DO 3144 K = J, NDIVI + 1
		  URCP(J) = URCP(J) - EPIRCP(K)*DR
3144		CONTINUE
		SIGT(J) = IPYCE*EPITP(J)
		SIGR(J) = -PRESS
		EPIR(J) = EPIRCP(J) - IPYCNU*EPITP(J)
		EPIT(J) = EPITCP(J) + EPITP(J)
		UR(J) = URCP(J) + URP(J)
3123	  CONTINUE
C  Calculate the strain energy stored in cracked IPyC layer due to shear force at the interface.
	  TEMP = 0.0D0
	  DO 3125 J = 1, NDIVI + 1
		R = (J-1)*DR + R2
		TEMP = TEMP + (EPITP(J)*R)**2.0*DR
3125	  CONTINUE
	  TEMP = TEMP*IPYCE/X2
	  WSHEARIPYC = (X1+IPYCNU)*KIIPYC**2.0*AIPYC*
     &			   (R3*(0.25D0*PIE-X1/X3)+0.5D0*AIPYC*
     &				(0.0625D0*PIE-X1/X3))/(16.0D0*IPYCE)
	  SHEARIPYC = X2*(TEMP+WSHEARIPYC/X2/PIE)/
     &			  (-PIE*R3**2.0*EPITCP(NDIVI+2))
C
C  2) Secondly evaluate the stresses, strains, and displacements in the cracked OPyC layer,
C  and calculate the shear force it applies to the SiC layer.
	  AOPYC = 0.99D0*(R5-R4)
	  SWELLROPYC = 0.0D0
	  SWELLTOPYC = 0.0D0
	  DO 3126 K=0,NDEG
		SWELLROPYC = SWELLROPYC + OSWR(K)*FLU**K
		SWELLTOPYC = SWELLTOPYC + OSWT(K)*FLU**K
3126	  CONTINUE
	  DO 3146 K = 1, NDIVO + 2
		URP(NDIVI+NDIVS+4+K) = 0.0D0
3146	  CONTINUE
	  DR = (R5-R4)/FLOAT(NDIVO+1)
	  CONV = .FALSE.
	  DO WHILE (.NOT.CONV)
		TEMP = KIOPYC
		SIGTPOBAR = 0.0D0
		DO 3127 J = 1, NDIVO + 2
		  R = (J-1)*DR + R4
		  EPITP(NDIVI+NDIVS+4+J)=((URCP(NDIVI+NDIVS+4+J)+
     &		URP(NDIVI+NDIVS+4+J))/R -
     &		EPITCP(NDIVI+NDIVS+4+J)-SWELLTOPYC*DF-8.0D0*
     &		(X1-OPYCNU**2.0)*KIOPYC*SQRT(ABS(R-R4)/(X2*PIE*R**2.0))
     &		/PIE/OPYCE)/(X1+OPYCE*OPYCREEP*DF)
		  SIGTPOBAR = SIGTPOBAR + OPYCE*EPITP(NDIVI+NDIVS+4+J)
3127		CONTINUE
		SIGTPOBAR = SIGTPOBAR/FLOAT(NDIVO+2)
	    KIOPYC = 0.413*(X1+R5/R4)*SIGTPOBAR*SQRT(PIE*AOPYC)/
     &			 SQRT(X1 - AOPYC/(R5-R4))
C  Approaches the right value of KIOPYC
		KIOPYC = (KIOPYC + TEMP)/X2
	    IF(ABS((TEMP-KIOPYC)/KIOPYC).LE.CRITERIA) CONV = .TRUE.
C  Recalculate URP
		DO 3147 J = 1, NDIVO + 2
		  URP(NDIVI+NDIVS+4+J) = 0.0D0
		  DO 3148 K = 1, J - 1
			URP(NDIVI+NDIVS+4+J) = URP(NDIVI+NDIVS+4+J) - 
     &				OPYCNU*EPITP(NDIVI+NDIVS+4+K)*DR
3148		  CONTINUE
3147		CONTINUE
	  END DO
C  Update EPITCP, and record resulting stresses, strains, and displacements in the cracked OPyC layer
	  DO 3128 J = 1, NDIVO + 2
		EPITCP(NDIVI+NDIVS+4+J) = EPITCP(NDIVI+NDIVS+4+J) + 
     &				OPYCREEP*OPYCE*EPITP(NDIVI+NDIVS+4+J)*DF +
     &				SWELLTOPYC*DF
		EPIRCP(NDIVI+NDIVS+4+J) = EPIRCP(NDIVI+NDIVS+4+J) +
     &				SWELLROPYC*DF
		URCP(NDIVI+NDIVS+4+J) = 0.0D0
		DO 3149 K = 1, J -1
		  URCP(NDIVI+NDIVS+4+J) = URCP(NDIVI+NDIVS+4+J) +
     &				EPIRCP(NDIVI+NDIVS+4+K)*DR
3149		CONTINUE
		SIGT(NDIVI+NDIVS+4+J) = OPYCE*EPITP(NDIVI+NDIVS+4+J)
		SIGR(NDIVI+NDIVS+4+J) = -PAMB
		EPIR(NDIVI+NDIVS+4+J) = EPIRCP(NDIVI+NDIVS+4+J) - 
     &				OPYCNU*EPITP(NDIVI+NDIVS+4+J)
		EPIT(NDIVI+NDIVS+4+J) = EPITCP(NDIVI+NDIVS+4+J) + 
     &				EPITP(NDIVI+NDIVS+4+J)
		UR(NDIVI+NDIVS+4+J) = URCP(NDIVI+NDIVS+4+J) +
     &				URP(NDIVI+NDIVS+4+J)
3128	  CONTINUE
C  Calculate the strain energy stored in cracked OPyC layer due to shear force at the interface.
	  TEMP = 0.0D0
	  DO 3130 J = 1, NDIVO + 1
		R = (J-1)*DR + R4
		TEMP = TEMP + (EPITP(NDIVI+NDIVS+4+J)*R)**2.0*DR
3130	  CONTINUE
	  TEMP = TEMP*OPYCE/X2
	  WSHEAROPYC = (X1+OPYCNU)*KIOPYC**2.0*AOPYC*
     &			   (R4*(0.25D0*PIE-X1/X3)-0.5D0*AOPYC*
     &				(0.0625D0*PIE-X1/X3))/(16.0D0*OPYCE)
	  SHEAROPYC = X2*(TEMP+WSHEAROPYC/X2/PIE)/
     &			  (-PIE*R4**2.0*EPITCP(NDIVI+NDIVS+5))
C
C  3) Finally calculate the stresses in SiC layer
C
C  Prepare quantities based on material properties
C
	  G1S = SICE/(X1 - X2*SICNU)
	  G2S = X2*SICE/(X1 + SICNU)
C  Prepare pressure
        PI = -PRESS
	  PO = -PAMB
C  Calculate stresses, strains, and displacement
	  DO 3018 K = 1, NDIVS+2
	    R = (K-1)*(R4-R3)/FLOAT(NDIVS+1) + R3
	    SIGR(NDIVI+2+K)=(PI*(R4**(-3.0)-R**(-3.0)) +
     &			PO*(R**(-3.0)-R3**(-3.0)))/(R4**(-3.0)-R3**(-3.0))
	    SIGT(NDIVI+2+K)=(PI*(X2*R4**(-3.0)+R**(-3.0))-
     &		PO*(R**(-3.0)+X2*R3**(-3.0)))/X2/(R4**(-3.0)-R3**(-3.0))
	    EPIR(NDIVI+2+K)=(PI*(R4**(-3.0))-PO*(R3**(-3.0)))
     &			/G1S/(R4**(-3.0)-R3**(-3.0)) -
     &			(X2*(PI-PO)*R**(-3.0))/G2S/(R4**(-3.0)-R3**(-3.0))
	    EPIT(NDIVI+2+K)=(PI*(R4**(-3.0))-PO*(R3**(-3.0)))
     &			/G1S/(R4**(-3.0)-R3**(-3.0)) +
     &			((PI-PO)*R**(-3.0))/G2S/(R4**(-3.0)-R3**(-3.0))
	    UR(NDIVI+2+K)  =(PI*(R4**(-3.0))-PO*(R3**(-3.0)))*R
     &			/G1S/(R4**(-3.0)-R3**(-3.0)) +
     &			((PI-PO)*R**(-2.0))/G2S/(R4**(-3.0)-R3**(-3.0)) 
3018	  CONTINUE
	END IF
C
      RETURN
	END
C                                                                      *
C***********************************************************************