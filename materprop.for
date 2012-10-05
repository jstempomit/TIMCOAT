C
C This library gives the dependence of coating ceramic material properties
C on various parameters. The reference is F. Ho, "Material Models of Pyro-
C carbon and Pyrolytic Silicon Carbide", CEGA-002820, July 1993
C
C 1. Dense pyrocarbon (density = 1.9 +/- 0.1g/cm^3)
C
C   1) Young's modulus (MPa)
C
C***********************************************************************
C                                                                      *
C  Function E_PYC                                                      *
C                                                                      *
C    Young's modulus of dense pyrocarbon as a function of density,     *
C  BAF0, crystallite size, neutron fluence and temperature.            *
C                                                                      *
C  Requirements:                                                       *
C    Density:        : 1.9 +/- 0.1 g/cm^3                              *
C    Crystallite size: 30 +/- 5 A                                      *
C    Neutron fluence : < 4.0*10^21 nvt                                 *
C    Temperature     : < 2000 C                                        *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual variable description                                         *
C    E_PYC(MPa)D: Young's modulus of dense pyrocarbon                  *
C    D (g/cm^3)     D: Density of dense pyrocarbon                     *
C    BAF0           D: Unirradiated BAF of pyrocarbon                  *
C    LC (A)         D: Crystallite size                                *
C    FLUENCE (10^21nvt) --                                             *
C                   D: Neutron fluence                                 *
C    T (C)          D: Temperature in dense PyC layer                  *
C  Local variable description                                          *
C                                                                      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
	FUNCTION E_PYC(D, BAF0, LC, FLUENCE, T)
C
	DOUBLE PRECISION E_PYC, D, BAF0, LC, FLUENCE, T
C
	IF(BAF0 .LT. 2.15D0) THEN
	  E_PYC = 2.55D4*(0.384D0 + 0.324D0*D)*(0.481D0 + 0.519D0*BAF0)
     &          *(2.985D0 - 0.0662D0*LC)*(1.0D0 + 0.23D0*FLUENCE)
     &          *(1.0D0 + 1.5D-4*(T-20.D0))
	ELSE
	  E_PYC = 2.55D4*(0.384D0 + 0.324D0*D)*1.605D0
     &          *(2.985D0 - 0.0662D0*LC)*(1.0D0 + 0.23D0*FLUENCE)
     &          *(1.0D0 + 1.5D-4*(T-20.0D0))
	END IF
C  Fix it temporarily
C	E_PYC = 2.55D4
	RETURN
	END
C                                                                      *
C***********************************************************************
C
C   2) Irradiated Bacon Anisotropy Factor (BAFi)
C
C***********************************************************************
C                                                                      *
C  Function BAFI_PYC                                                   *
C                                                                      *
C    Irradiated Bacon Anisotropy Factor of pyrocarbon as a function of *
C  neutron fluence under unrestrained condition.                       *
C                                                                      *
C  Actual variable description                                         *
C    BAFI_PYC       D: Irradiated BAF of pyrocarbon                    *
C    BAF0           D: Unirradiated BAF of pyrocarbon                  *
C    FLUENCE (10^21nvt) --                                             *
C                   D: Neutron fluence                                 *
C  Local variable description                                          *
C    BAFI_F         D: Array giving the dependence of BAFI on fluence  *
C                      The fluence range is 0 - 5D21 nvt.              *
C    BAFI           D: Irradiated BAF for density = 1.9g/cm^3          *
C    C              D: the ratio of BAFI/BAF0                          *
C    INDEX          I: Array index                                     *
C                                                                      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
	FUNCTION BAFI_PYC(BAF0, FLUENCE)
C
	DOUBLE PRECISION BAFI_PYC, BAF0, FLUENCE
	DOUBLE PRECISION BAFI_F(1:8, -1:0)
	DOUBLE PRECISION BAFI, C
	INTEGER  INDEX
	DATA BAFI_F / 0.0D0,   1.0D0,   2.0D0,   3.0D0,   3.5D0,   4.0D0,
     &              4.5D0,   5.0D0,  1.05D0, 1.052D0, 1.062D0, 1.073D0,
     &             1.08D0, 1.084D0, 1.088D0,  1.09D0/
C
	IF(FLUENCE .LE. BAFI_F(1,-1)) THEN
	  INDEX = 1
        CALL ERR_HANDLER('FUNCTION BAFI_PYC: Neutron fluence 
     &               below lower limit in BAFI_F', 61, 0, 0, IERR)  !lower limit is 0 here.
	ELSE IF(FLUENCE .GE. BAFI_F(8,-1)) THEN
	  INDEX = 8 - 1
        CALL ERR_HANDLER('FUNCTION BAFI_PYC: Neutron fluence 
     &               above upper limit in BAFI_F', 61, 0, 0, IERR)  !upper limit is 5 here.
	ELSE
	  CALL LOCATE_ARRAY(BAFI_F(:,-1),8,1,FLUENCE,INDEX)
	END IF
C      
	BAFI = ((FLUENCE-BAFI_F(INDEX,-1))*BAFI_F(INDEX+1,0)+
     &        (BAFI_F(INDEX+1,-1)-FLUENCE)*BAFI_F(INDEX,0))/
     &       (BAFI_F(INDEX+1,-1)-BAFI_F(INDEX,-1))
	C = BAFI/BAFI_F(1,0)
	BAFI_PYC = BAF0*C
	RETURN
	END
C                                                                      *
C***********************************************************************
C
C   3.a) Irradiation-induced dimensional changes (unrestrained)
C
C***********************************************************************
C  Subroutine SWELLU(TEMP,FLUENCE,DENSITY,BAF0,CRATE,SR_DOT,ST_DOT)    *
C                                                                      *
C    This subroutine provides the UNRESTRAINED radial and tangential   *
C  irradiation-induced dimensional changes in IPyC or OPyC layers of   *
C  TRISO fuel with temperature, fast fluence, density, BAF0 and coating*
C  rate given for that layer. The algorithm of picking out the data is *
C  from p2-30 Fig. 5 of 'Material model of pyrocarbon and pyrolytic    *
C  silicon carbide', CEGA-002820.                                      *
C                                                                      *
C  Requirements:                                                       *
C    1) The fast fluence MUST be lower than 4.0*10^21nvt because       *
C       the fitting curves come from data below this limit.            *
C    2) The BAF should be in the range of 1.00 to 1.2787.              *
C    3) The density should be in the range of 1.66 to 1.99 g/cm^3.     *
C    4) The coating rate should be in the range of 0 to 12 um/min.     *
C    5) The temperature that could be modeled is 600 - 1350 C.         *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual argument description                                         *
C    TEMP (C)       D: Current temperature of that layer (degree C)    *
C    FLUENCE (10^21n/cm^2) --                                          *
C                   D: Current fast fluence                            *
C    DENSITY (g/cm^3) --                                               *
C                   D: Density of the layer                            *
C    BAF0           D: Initial Bacon Anisotropy Factor (dimensionless) *
C    CRATE (um/min) D: Coating rate of the layer                       *
C    SR_DOT (/10^21nvt) --                                             *
C                   D: Returned value of radial strain rate            *
C    ST_DOT (/10^21nvt) --                                             *
C                   D: Returned value of tangential strain rate        *
C                                                                      *
C  Local variable description                                          *
C  (1) Arrays:                                                         *
C    BAF_D_CR       D: Array of data giving the dependency of BAF      *
C                      on density and average coating rate             *
C                      (from Figure A-8, CEGA-002820)                  *
C    BAF_CR         D: BAF0 as a function of coating rate at the       *
C                      inputed density                                 *
C    EISO_T_F       D: Array of data giving denpendency of Epi_ISO on  *
C                      temperature and fast fluence (from Table 1,CEGA)*
C    EISO_F         D: Epi_ISO as a function of fast fluence at the    *
C                      inputed temperature                             *
C    EISO_D         D: Epi_ISO as a function of density at T=1100C and *
C                      fluence = 3.7*10^21nvt (from Table on pp. 2.25  *
C                      of CEGA)                                        *
C    ER_T_F_BAF     D: The dependency of Er on temperature, fast       *
C                      fluence and BAF01 at density=1.96               *
C                      (Table 2 of CEGA-002820)                        *
C    ET_T_F_BAF     D: The dependency of Et on temperature, fast       *
C                      fluence and BAF01 at density=1.96               *
C                      (Table 2 of CEGA-002820)                        *
C    ER_T_F         D: The dependency of Er on temperature, fast       *
C                      fluence at derived BAF01 and density=1.96       *
C    ET_T_F         D: The dependency of Et on temperature, fast       *
C                      fluence at derived BAF01 and density=1.96       *
C    ER_F           D: Er as a function of fluence at given            *
C                      temperature, BAF01 and density = 1.96           *
C    ET_F           D: Et as a function of fluence at given            *
C                      temperature, BAF01 and density = 1.96           *
C    DE_D           D: Delta E as a function of density at T=1100C and *
C                      fluence = 3.7*10^21nvt (From Fig. 4 of CEGA)    *
C    COEFF          D: Array used in solving non-linear equations with *
C                      IMSL                                            *
C    ROOT           D: Array used in solving non-linear equations with *
C                      IMSL                                            *
C  (2) Simple variables:                                               *
C    BAF01          D: The corresponding BAF at density = 1.96g/cm^3   *
C                      given the coating rate remains unchanged        *
C    BAF0B          D: The copy of original BAF0 when BAF0 is out of   *
C                      the range of given database                     *
C    DE (%)         D: Delta E at all specified quantities             *
C    DE_DOT (%)     D: The rate of Delta E at all specified quantities *
C    EISO (%)       D: Epi_ISO at all specified quantities             *
C    EISO_DOT (%)   D: The rate of Epi_ISO at all specified quantities *
C    ER (%)         D: Er at density = 1.96                            *
C    ER_DOT (%)     D: The rate of Er at density = 1.96                *
C    ET (%)         D: Et at density = 1.96                            *
C    ET_DOT (%)     D: The rate of Et at density = 1.96                *
C    SCALING        D: Scalling factor for adjust EISO and DE due to   *
C                      density difference                              *
C    TEMP0 (C)      D: Reference temperature in case the given         *
C                      temperature is out of range                     *
C                                                                      *
C  Parameters and counters                                             *
C    BGROUPS        I: Number of BAF entries in related arrays         *
C    CRMAX          I: Maximum allowable coating rate                  *
C    DGROUPS        I: Number of density entries in related arrays     *
C    IEER           I: Lun for error message file  ('ERROR.msg')       *
C    INDEX          I: General-purpose indices for array entry locating*
C    INDEXR         I: Indices for array entry locating when dealing   *
C                      with data array of material radial direction    *
C    INDEXT         I: Indices for array entry locating when dealing   *
C                      with data array of material tangential direction*
C    NDEG           I: Degree of polynomial for swelling (not rate)    *
C    TGROUPS        I: Number of temperature entries in related arrays *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE SWELLU(TEMP,FLUENCE,DENSITY,BAF0,CRATE,SR_DOT,ST_DOT)
C
      DOUBLE PRECISION TEMP,FLUENCE,DENSITY,BAF0,CRATE,SR_DOT,ST_DOT
      DOUBLE PRECISION BAF_D_CR, BAF_CR, BAF01, BAF0B
      DOUBLE PRECISION EISO_T_F, EISO_F, EISO, TEMP0
      DOUBLE PRECISION EISO_D, DE_D, SCALING
      DOUBLE PRECISION ER_T_F_BAF, ET_T_F_BAF, ER_T_F, ET_T_F, 
     &                 ER_F, ET_F, ER, ET, DE, SR, ST
	DOUBLE PRECISION EISO_DOT, DE_DOT, ER_DOT, ET_DOT
      DOUBLE PRECISION COEFF(0:3)
      DOUBLE COMPLEX ROOT(3)
      INTEGER  CRMAX, DGROUPS, TGROUPS, BGROUPS
      INTEGER  INDEX, INDEXR, INDEXT, I
      PARAMETER (NDEG=4, CRMAX=13, DGROUPS=8, 
     &           TGROUPS=3, BGROUPS=6, IERR = 12)
      DIMENSION ER_T_F_BAF(1:TGROUPS,1:BGROUPS,-1:NDEG)
      DIMENSION ET_T_F_BAF(1:TGROUPS,1:BGROUPS,-1:NDEG)
      DIMENSION ER_T_F(1:TGROUPS,-1:NDEG), ET_T_F(1:TGROUPS,-1:NDEG)
      DIMENSION ER_F(0:NDEG), ET_F(0:NDEG)
      DIMENSION BAF_D_CR(1:DGROUPS,-1:CRMAX), BAF_CR(0:CRMAX)
      DIMENSION EISO_T_F(1:TGROUPS,-1:NDEG), EISO_F(0:NDEG)
      DIMENSION EISO_D(1:9,-1:0), DE_D(1:9,-1:0)
      EXTERNAL DZPLRC, WRCRN
      DATA ER_T_F_BAF /1.00D0,1.00D0,1.00D0,1.0212D0,1.0212D0,1.0212D0,
     &           1.0488D0,1.0488D0,1.0488D0,1.0769D0,1.0769D0,1.0769D0,
     &           1.1746D0,1.1746D0,1.1746D0,1.2787D0,1.2787D0,1.2787D0,
     &           0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,
     &           0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,
     &           0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,
     &     -1.2408D0,-1.5239D0,-1.4284D0,-1.1064D0,-2.0752D0,-1.5433D0,
     &  -0.94333D0,-2.0047D0,-1.4964D0,-0.78045D0,-1.8169D0,-0.89522D0,
     &     -0.15714D0,-1.1854D0,1.2093D0,0.40265D0,-0.45900D0,3.7162D0,
     &   0.00175D0,0.13048D0,-0.19563D0,-0.03128D0,1.37845D0,0.59804D0,
     &     -0.03589D0,1.3038D0,1.16621D0,-0.02975D0,1.1085D0,0.80331D0,
     &  -0.14889D0,0.64995D0,-0.53861D0,-0.16501D0,0.51172D0,-2.7042D0,
     &   0.08533D0,0.06299D0,0.18991D0,0.09184D0,-0.48993D0,-0.09997D0,
     &  0.08184D0,-0.3728D0,-0.30106D0,0.06655D0,-0.23868D0,-0.09009D0,
     &     0.07546D0,0.01380D0,0.43114D0,0.03676D0,-0.03245D0,1.1799D0,
     & -0.01253D0,-0.01072D0,-0.02591D0,-0.01220D0,0.06602D0,0.00978D0,
     &   -0.00958D0,0.04538D0,0.03475D0,-0.00626D0,0.02484D0,0.00467D0,
     & -0.00293D0,-0.01284D0,-0.05590D0,0.00706D0,-0.00142D0,-0.1391D0/
C
      DATA ET_T_F_BAF /1.00D0,1.00D0,1.00D0,1.0303D0,1.0303D0,1.0303D0,
     &           1.0769D0,1.0769D0,1.0769D0,1.1250D0,1.1250D0,1.1250D0,
     &           1.2258D0,1.2258D0,1.2258D0,1.3333D0,1.3333D0,1.3333D0,
     &           0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,
     &           0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,
     &           0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,0.0000D0,
     &     -1.2408D0,-1.5239D0,-1.4284D0,-1.3855D0,-1.5759D0,-2.2468D0,
     &     -1.4679D0,-1.3220D0,-2.8293D0,-1.6466D0,-1.1870D0,-3.2555D0,
     &   -1.8499D0,-0.96963D0,-4.4478D0,-2.1919D0,-0.81239D0,-5.6714D0,
     &    0.00175D0,0.13048D0,-0.19563D0,0.05307D0,0.09019D0,0.48243D0,
     &  -0.02836D0,-0.51928D0,0.76088D0,0.03928D0,-0.90635D0,0.90423D0,
     &      -0.09358D0,-1.5911D0,1.6032D0,0.02675D0,-2.2076D0,2.4192D0,
     &    0.08533D0,0.06299D0,0.18991D0,0.07620D0,0.05306D0,-0.07687D0,
     &   0.12139D0,0.27603D0,-0.22314D0,0.10067D0,0.41046D0,-0.33175D0,
     &   0.18119D0,0.64689D0,-0.58683D0,0.15352D0,0.88496D0,-0.86155D0,
     &-0.01253D0,-0.01072D0,-0.02591D0,-0.01245D0,-0.00815D0,0.00464D0,
     & -0.01948D0,-0.03465D0,0.02431D0,-0.01764D0,-0.05067D0,0.04329D0,
     & -0.03036D0,-0.07682D0,0.07458D0,-0.02972D0,-0.10457D0,0.10668D0/
C
      DATA ER_T_F(:,-1) /600.0D0,1032.0D0,1350.0D0/
      DATA ET_T_F(:,-1) /600.0D0,1032.0D0,1350.0D0/
C
C  The second column of array BAF_D_CR, which corresponds to BAF0 at coating rate of
C  0, was changed all to the value of 1.34, so that BAF0 range of 1.025 to 1.34 is
C  valid for density of range 1.66 - 1.99. Also, the column of coating rate 13um/min
C  was arbitrarily added and set to value 1.0, so that BAF0 is kept above 1.0.
C  (Modification occured on 03/03/03)
      DATA BAF_D_CR 
     &   /1.0D0,1.660D0,1.740D0,1.830D0,1.880D0,1.940D0,1.960D0,1.990D0,
     &    1.34D0,1.34D0,1.340D0,1.340D0,1.340D0,1.340D0,1.340D0,1.340D0,
     &    1.0D0,1.03000,1.03622,1.04000,1.05000,1.05944,1.07139,1.08333,
     &    1.0D0,1.01611,1.01978,1.02267,1.02867,1.03300,1.04089,1.04878,
     &    1.0D0,1.01067,1.01311,1.01656,1.02044,1.02411,1.03072,1.03733,
     &    1.0D0,1.00889,1.01156,1.01433,1.01767,1.02000,1.02605,1.03211,
     &    1.0D0,1.00833,1.01056,1.01333,1.01633,1.01889,1.02445,1.03000,
     &    1.0D0,1.00780,1.01000,1.01311,1.01556,1.01844,1.02356,1.02867,
     &    1.0D0,1.00779,1.00978,1.01278,1.01501,1.01778,1.02273,1.02767,
     &    1.0D0,1.00778,1.00977,1.01277,1.01500,1.01777,1.02250,1.02722,
     &    1.0D0,1.00777,1.00976,1.01276,1.01456,1.01776,1.02205,1.02644,
     &    1.0D0,1.00776,1.00975,1.01275,1.01422,1.01775,1.02184,1.02600,
     &    1.0D0,1.00775,1.00974,1.01222,1.01421,1.01774,1.02161,1.02544,
     &    1.0D0,1.00774,1.00973,1.01221,1.01420,1.01773,1.02150,1.02500,
     &	1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0/
C
      DATA EISO_T_F /600.0D0,     1032.0D0,     1350.0D0,
     &               0.0D0,       0.0D0,        0.0D0,
     &               -1.4764D0,   -1.7111D0,    -1.8886D0,
     &               0.29487D0,   0.38139D0,    0.38982D0,
     &               -0.02790D0,  -0.03674D0,   -0.03691D0,
     &               0.00105D0,   0.00135D0,    0.00130D0/
      DATA EISO_D  /1.000D0,1.200D0,1.400D0,1.500D0,1.600D0,1.800D0,
     &              1.900D0,1.960D0,2.000D0,-16.15D0,-13.11D0,-9.98D0,
     &              -8.50D0,-6.97D0,-4.42D0,-3.41D0,-2.75D0,-2.33D0/
      DATA DE_D    /1.00D0, 1.20D0, 1.40D0, 1.50D0, 1.60D0, 1.80D0,
     &              1.90D0, 1.96D0, 2.00D0, 0.00D0, 0.50D0, 1.10D0,
     &              1.65D0, 2.45D0, 6.305D0, 7.9D0, 9.60D0, 11.1D0/
      DATA COEFF   /0.0D0,0.0D0,0.0D0,1.0D0/
C
C  Step 1: Convert the given set of density, BAF0 and coating rate to BAF01 at density = 1.96g/cm^3
      IF(DENSITY.LT.BAF_D_CR(1,-1)) THEN
c        DENSITY = BAF_D_CR(1,-1)         !Currently BAF_D_CR(1,-1)=1.00
        INDEX = 1
        CALL ERR_HANDLER('SUBROUTINE SWELLU: Density of IPyC/OPyC layer
     & below lower limit in BAF_D_CR', 75, 0, 0, IERR)
      ELSE IF(DENSITY.GT.BAF_D_CR(DGROUPS,-1)) THEN
c        DENSITY = BAF_D_CR(DGROUPS,-1)   !Currently BAF_D_CR(DGROUPS,-1)=1.99
        INDEX = DGROUPS - 1
        CALL ERR_HANDLER('SUBROUTINE SWELLU: Density of IPyC/OPyC layer
     & above upper limit in BAF_D_CR', 75, 0, 0, IERR)
      ELSE
        CALL LOCATE_ARRAY(BAF_D_CR(:,-1), DGROUPS, 1, DENSITY, INDEX)
        IF((INDEX.LT.1).OR.(INDEX.GT.(DGROUPS-1))) THEN
          CALL ERR_HANDLER('SUBROUTINE SWELLU: Out of Array BAF_D_CR 
     &					  range', 45, 0, 2, IERR)
        END IF
      END IF
      DO 2910 I=0,CRMAX
        BAF_CR(I)=((DENSITY-BAF_D_CR(INDEX,-1))*BAF_D_CR(INDEX+1,I)+
     &            (BAF_D_CR(INDEX+1,-1)-DENSITY)*BAF_D_CR(INDEX,I))/
     &            (BAF_D_CR(INDEX+1,-1)-BAF_D_CR(INDEX,-1))
2910  CONTINUE
C   look for coating rate based on the given BAF0
	BAF0B = BAF0
      IF(BAF0.GT.BAF_CR(0)) THEN   !BAF0 descends as coating rate increases
c        BAF0B = BAF_CR(0)
        INDEX = 0
        CALL ERR_HANDLER('SUBROUTINE SWELLU: BAF0 of IPyC/OPyC layer
     & above upper limit in BAF_D_CR', 72, 0, 0, IERR)
      ELSE IF(BAF0.LT.BAF_CR(CRMAX)) THEN
c        BAF0B = BAF_CR(CRMAX)
        INDEX = CRMAX - 1
        CALL ERR_HANDLER('SUBROUTINE SWELLU: BAF0 of IPyC/OPyC layer
     & below lower limit in BAF_D_CR', 72, 0, 0, IERR)
      ELSE
C        BAF0B = BAF0
	  CALL LOCATE_ARRAY(BAF_CR,CRMAX+1,0,BAF0B,INDEX)
        IF((INDEX.LT.0).OR.(INDEX.GT.(CRMAX-1))) THEN
          CALL ERR_HANDLER('SUBROUTINE SWELLU: Out of Array BAF_D_CR 
     &					  range', 45, 0, 2, IERR)
        END IF
      END IF
      CRATE = ((BAF0B-BAF_CR(INDEX))*FLOAT(INDEX+1) + (BAF_CR(INDEX+1)
     &         -BAF0B)*FLOAT(INDEX))/(BAF_CR(INDEX+1)-BAF_CR(INDEX))
C   convert to BAF01 at density = 1.96g/cm^3, which is given by BAF_D_CR(7,*)
      BAF01 = (CRATE-INDEX)*BAF_D_CR(7,INDEX+1) + 
     &        (INDEX+1-CRATE)*BAF_D_CR(7,INDEX)
C  Step 2: Obtain Epi_ISO and the rate of Epi_ISO for the given T and Phi at the density of 1.96g/cm^3
	TEMP0 = TEMP
      IF(TEMP.LE.EISO_T_F(1,-1)) THEN
c        TEMP0 = EISO_T_F(1,-1)           !Currently EISO_T_F(1,-1)=600.0
        INDEX = 1
c	  CALL ERR_HANDLER('SUBROUTINE SWELLU: Temperature
c     & below lower limit in EISO_T_F', 60, 0, 0, IERR)
      ELSE IF(TEMP.GE.EISO_T_F(TGROUPS,-1)) THEN
c        TEMP0 = EISO_T_F(TGROUPS,-1)     !Currently EISO_T_F(TGROUPS,-1)=1350.0
        INDEX = TGROUPS - 1
c	  CALL ERR_HANDLER('SUBROUTINE SWELLU: Temperature
c     & above upper limit in EISO_T_F', 60, 0, 0, IERR)
      ELSE
c        TEMP0 = TEMP
        CALL LOCATE_ARRAY(EISO_T_F(:,-1), TGROUPS, 1, TEMP, INDEX)
        IF((INDEX.LT.1).OR.(INDEX.GT.(TGROUPS-1))) THEN
          CALL ERR_HANDLER('SUBROUTINE SWELLU: Out of Array EISO_T_F 
     &                range', 45, 0, 2, IERR)
        END IF
      END IF
      EISO = 0.0
	EISO_DOT = 0.0
      DO 2920 I=0,NDEG
        EISO_F(I)=((TEMP0-EISO_T_F(INDEX,-1))*EISO_T_F(INDEX+1,I)+
     &            (EISO_T_F(INDEX+1,-1)-TEMP0)*EISO_T_F(INDEX,I))/
     &            (EISO_T_F(INDEX+1,-1)-EISO_T_F(INDEX,-1))
        EISO = EISO + EISO_F(I)*FLUENCE**I
	  IF(I.NE.0) EISO_DOT = EISO_DOT + I*EISO_F(I)*FLUENCE**(I-1)
2920  CONTINUE
C  Step 3: Adjust Epi_ISO and the rate of Epi_ISO for density difference
      CALL LOCATE_ARRAY(EISO_D(:,-1),9,1,DENSITY,INDEX)
      SCALING = ((DENSITY-EISO_D(INDEX,-1))*EISO_D(INDEX+1,0)+
     &          (EISO_D(INDEX+1,-1)-DENSITY)*EISO_D(INDEX,0))/
     &          (EISO_D(INDEX+1,-1)-EISO_D(INDEX,-1))
C   Take the ratio of Epi_ISO at this density with respect to that at 1.96g/cm^3, given by EISO_D(8.0)
      SCALING = SCALING/EISO_D(8,0)
      EISO = EISO*SCALING
	EISO_DOT = EISO_DOT*SCALING
C
	IF(BAF0.EQ.1.0D0) THEN
	  SR_DOT = EISO_DOT
	  ST_DOT = EISO_DOT
C   Change from percentage to real dimensional change rate
        SR_DOT = SR_DOT/100.0D0
        ST_DOT = ST_DOT/100.0D0
        RETURN
      END IF
C  Step 4: Obtain Epi_r and Epi_t and Delta_E and their rates at T, Phi, BAF01 and density = 1.96
      DO 2940 I=1,TGROUPS
	  IF(BAF01 .LE. ER_T_F_BAF(I,1,-1)) THEN
	    INDEXR = 1
	    CALL ERR_HANDLER('SUBROUTINE SWELLU: BAF0 of IPyC/OPyC layer
     & below lower limit in ER_T_F_BAF', 74, 0, 0, IERR)
	  ELSE IF(BAF01 .GE. ER_T_F_BAF(I,BGROUPS,-1)) THEN
	    INDEXR = BGROUPS - 1
	    CALL ERR_HANDLER('SUBROUTINE SWELLU: BAF0 of IPyC/OPyC layer
     & above upper limit in ER_T_F_BAF', 74, 0, 0, IERR)
	  ELSE
	    CALL LOCATE_ARRAY(ER_T_F_BAF(I,:,-1),BGROUPS,1,BAF01,INDEXR)
	  END IF
        DO 2930 J=0,NDEG
          ER_T_F(I,J)=((BAF01-ER_T_F_BAF(I,INDEXR,-1))*
     &              ER_T_F_BAF(I,INDEXR+1,J)+(ER_T_F_BAF(I,INDEXR+1,-1)-
     &              BAF01)*ER_T_F_BAF(I,INDEXR,J))/
     &              (ER_T_F_BAF(I,INDEXR+1,-1)-ER_T_F_BAF(I,INDEXR,-1))
2930    CONTINUE
	  IF(BAF01 .LE. ET_T_F_BAF(I,1,-1)) THEN
	    INDEXT = 1
	    CALL ERR_HANDLER('SUBROUTINE SWELLU: BAF0 of IPyC/OPyC layer
     & below lower limit in ET_T_F_BAF', 74, 0, 0, IERR)
	  ELSE IF(BAF01 .GE. ET_T_F_BAF(I,BGROUPS,-1)) THEN
	    INDEXT = BGROUPS - 1
	    CALL ERR_HANDLER('SUBROUTINE SWELLU: BAF0 of IPyC/OPyC layer
     & above upper limit in ET_T_F_BAF', 74, 0, 0, IERR)
	  ELSE
	    CALL LOCATE_ARRAY(ET_T_F_BAF(I,:,-1),BGROUPS,1,BAF01,INDEXT)
	  END IF
        DO 2935 J=0,NDEG
          ET_T_F(I,J)=((BAF01-ET_T_F_BAF(I,INDEXT,-1))*
     &              ET_T_F_BAF(I,INDEXT+1,J)+(ET_T_F_BAF(I,INDEXT+1,-1)-
     &              BAF01)*ET_T_F_BAF(I,INDEXT,J))/
     &              (ET_T_F_BAF(I,INDEXT+1,-1)-ET_T_F_BAF(I,INDEXT,-1))
2935    CONTINUE
2940  CONTINUE
      IF(TEMP0.LE.ER_T_F(1,-1)) THEN
C        TEMP0 = ER_T_F(1,-1)          !Currently ER_T_F(1,-1)=ET_T_F(1,-1)=600.0
        INDEXR = 1
        INDEXT = 1
c	  CALL ERR_HANDLER('SUBROUTINE SWELLU: Temperature
c     & below lower limit in ER_T_F and ET_T_F', 69, 0, 0, IERR)
      ELSE IF(TEMP0.GE.ER_T_F(TGROUPS,-1)) THEN
C        TEMP0 = ER_T_F(TGROUPS,-1)    !Currently ER_T_F(TGROUPS,-1)=ET_T_F(TGROUPS,-1)=1350.0
        INDEXR = TGROUPS -1
        INDEXT = TGROUPS -1
c	  CALL ERR_HANDLER('SUBROUTINE SWELLU: Temperature
c     & above upper limit in ER_T_F and ET_T_F', 69, 0, 0, IERR)
      ELSE
        CALL LOCATE_ARRAY(ER_T_F(:,-1),TGROUPS,1,TEMP0,INDEXR)
        CALL LOCATE_ARRAY(ET_T_F(:,-1),TGROUPS,1,TEMP0,INDEXT)
        IF((INDEXR.LT.1).OR.(INDEXR.GT.(TGROUPS-1))) THEN
          CALL ERR_HANDLER('SUBROUTINE SWELLU: Out of Array ER_T_F 
     &                range', 45, 0, 2, IERR)
        END IF
        IF((INDEXT.LT.1).OR.(INDEXT.GT.(TGROUPS-1))) THEN
          CALL ERR_HANDLER('SUBROUTINE SWELLU: Out of Array ET_T_F 
     &                range', 45, 0, 2, IERR)
        END IF
      END IF
      ER = 0.0D0
      ET = 0.0D0
      DE = 0.0D0
	ER_DOT = 0.0D0
	ET_DOT = 0.0D0
	DE_DOT = 0.0D0
      DO 2950 I=0,NDEG
        ER_F(I)=((TEMP0-ER_T_F(INDEXR,-1))*ER_T_F(INDEXR+1,I) + 
     &          (ER_T_F(INDEXR+1,-1)-TEMP0)*ER_T_F(INDEXR,I))/
     &          (ER_T_F(INDEXR+1,-1)-ER_T_F(INDEXR,-1))
        ET_F(I)=((TEMP0-ET_T_F(INDEXT,-1))*ET_T_F(INDEXT+1,I) + 
     &          (ET_T_F(INDEXT+1,-1)-TEMP0)*ET_T_F(INDEXT,I))/
     &          (ET_T_F(INDEXT+1,-1)-ET_T_F(INDEXT,-1))
        ER = ER + ER_F(I)*FLUENCE**I
        ET = ET + ET_F(I)*FLUENCE**I
	  IF(I.NE.0) THEN
	    ER_DOT = ER_DOT + I*ER_F(I)*FLUENCE**(I-1)
	    ET_DOT = ET_DOT + I*ET_F(I)*FLUENCE**(I-1)
	  END IF
2950  CONTINUE
      DE = ER - ET
	DE_DOT = ER_DOT - ET_DOT
C  Step 5: Adjust Delta_E and the rate of it back to the input density and BAF0
      CALL LOCATE_ARRAY(DE_D(:,-1),9,1,DENSITY,INDEX)
      SCALING = ((DENSITY-DE_D(INDEX,-1))*DE_D(INDEX+1,0)+
     &          (DE_D(INDEX+1,-1)-DENSITY)*DE_D(INDEX,0))/
     &          (DE_D(INDEX+1,-1)-DE_D(INDEX,-1))
C   Take the ratio of Delta_E at this density with respect to that at 1.96g/cm^3, given by DE_D(8.0)
      SCALING = SCALING/DE_D(8,0)
      DE = DE*SCALING
	DE_DOT = DE_DOT*SCALING
C  Step 6: Use Delta_E and Epi_ISO to derive Er and Et
C    (100 - Et)**3 - Delta_E*(100- Et)**2 - (100 - Epi_ISO)**3 = 0
      IF(FLUENCE.NE.0.0D0) THEN
	  COEFF(0)=(-1.0D0)*(100.0D0-EISO)**3
        COEFF(2)=(-1.0D0)*DE
        CALL DZPLRC(3,COEFF,ROOT)
        ST = -100.0D0
        DO 2960 I=1,3
          IF(ABS(AIMAG(ROOT(I))).LT.1.0E-10)  ST=100.0D0-REAL(ROOT(I))
2960    CONTINUE
        IF(ST.NE.-100.0D0) THEN
          SR = DE + ST
        ELSE
          CALL ERR_HANDLER('SUBROUTINE SWELLU: Cannot determine 
     & dimensional changes', 54, 0, 2, IERR)
        END IF
	ELSE
	  SR = 0.0D0
	  ST = 0.0D0
	END IF
	ST_DOT = (3.0D0*EISO_DOT - (100.0D0-EISO)*DE_DOT/(100.0D0-SR))/
     &         ((100.0D0-EISO)/(100.0D0-SR)+2.0D0*(100.0D0-EISO)/
     &          (100.0D0-ST))
      SR_DOT = DE_DOT + ST_DOT
C  Change from percentage to real dimensional change rate
      SR_DOT = SR_DOT/100.0D0
      ST_DOT = ST_DOT/100.0D0
C  Done!
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C   3.b) Irradiation-induced dimensional changes (restrained)
C
C***********************************************************************
C  Subroutine SWELLR(BAFI, EC, DFLU, SR_DOT, ST_DOT)                   *
C                                                                      *
C    This subroutine calculates the restrained irradiation-induced     *
C  dimensional changes in PyC given the unrestrained strains due to    *
C  creep in stressed material. The algorithm of conversion is          *
C  from p2-33 to p2-34, CEGA-002820.                                   *
C                                                                      *
C  Requirements:                                                       *
C    1) The fast fluence MUST be lower than 4.0*10^21nvt because       *
C       the fitting curves come from data below this limit.            *
C    2) The BAF should be in the range of 1.00 to 1.33.                *
C                                                                      *
C  Actual argument description                                         *
C    BAFI           D: Irradiated BAF (dimensionless)                  *
C    EC             D: Array of apparent creep strain                  *
C                      EC(1): Radial component                         *
C                      EC(2): Tangential component                     *
C    DFLU (10^21nvt)D: Incremental fluence gained in one step          *
C    SR_DOT (/10^21nvt) --                                             *
C                   D: radial strain rate                              *
C    ST_DOT (/10^21nvt) --                                             *
C                   D: tangential strain rate                          *
C                                                                      *
C  Local variable description                                          *
C    BAFIA          D: Adjusted Irradiated BAF due to restrained       *
C                      condition (dimensionless)                       *
C    DSR            D: Incremental unrestrained radial irradiation     *
C                      strain over the last fluence step               *
C    DST            D: Incremental unrestrained tangential irradiation *
C                      strain over the last fluence step               *
C    SAR            D: reorientation radial strain                     *
C    SAT            D: reorientation tangential strain                 *
C    SBR            D: densification radial strain                     *
C    SBT            D: densification tangential strain                 *
C    RI             D: Preferred orientation parameter                 *
C    XA             D: Macroscopic dimensional strains of crystallites *
C                      in the a direction of PyC                       *
C    XC             D: Macroscopic dimensional strains of crystallites *
C                      in the c direction of PyC                       *
C                                                                      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE SWELLR(BAFI, EC, DFLU, SR_DOT, ST_DOT)
C
      DOUBLE PRECISION BAFI, EC, DFLU, SR_DOT, ST_DOT
	DOUBLE PRECISION BAFIA, RI, DSR, DST, SAR, SAT, SBR, SBT, XA, XC
	DIMENSION EC(1:2)
C  Preparation for calculation
	IF(((EC(1).EQ.0.0D0).AND.(EC(2).EQ.0.0D0)).OR.(DFLU.EQ.0.0D0)
     &	.OR.(BAFI.EQ.1.0D0)) RETURN
	RI = 2.0D0/(2.0D0 + BAFI)
C  Adjust for restrained condition
	DSR = SR_DOT*DFLU
	DST = ST_DOT*DFLU
	SAR = 2.0D0*(DSR - DST)/3.0D0
	SAT = (DST - DSR)/3.0D0
	SBR = (DSR + 2.0D0*DST)/3.0D0
	SBT = SBR
	XA = 2.0D0*(DST-DSR)/(3.0D0*(2.0D0-3.0D0*RI))
	XC = 4.0D0*(DSR-DST)/(3.0D0*(2.0D0-3.0D0*RI))
C    adjust restrained BAFI
C     The following formula is from p2-36 CEGA-002820 by assuming the principle
C     of independent action and the material remains transversely isotropic.
C
C	BAFIA = BAFI*(1.0D0+5.25D0*EC(2))*DSQRT(1.0D0-5.25D0*EC(1))
C
C     However, the above makes a huge increase in BAFI, which is unrealistic.
C     Considering the shear strain may be effective in changing pyrocarbon's
C     crystallite alignment, we use effective or octahedral-shear creep strain
C     instead of apparent strain.
      BAFIA = BAFI*DSQRT(1.0D0 + 5.25*(DSQRT(2.0D0/3.0D0)*
     &                   ABS(EC(1)+EC(2))))
	RI = 2.0D0/(2.0D0 + BAFIA)
	SAR = RI*XA + (1.0D0 - RI)*XC
	SAT = (1.0D0 - RI/2.0D0)*XA + (RI/2.0D0)*XC
	DSR = SAR + SBR
	DST = SAT + SBT
	SR_DOT = DSR/DFLU
	ST_DOT = DST/DFLU
C  Done!
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C   4) Steady state irradiation creep
C
C***********************************************************************
C                                                                      *
C  Subroutine SSICREEP_PYC                                             *
C                                                                      *
C    Steady state irradiation creep coefficient and Poisson's ratio    *
C  for PyC.                                                            *
C                                                                      *
C  Requirements:                                                       *
C    Density:        : <2.05 g/cm^3                                    *
C                                                                      *
C  Actual variable description                                         *
C    CREEP(10^-21/MPa/nvt) --                                          *
C                   D: Steady state irradiation creep coefficient      *
C    NU_C           D: Poisson's ratio for irradiation creep           *
C    D (g/cm^3)     D: Density of dense pyrocarbon                     *
C    FLUENCE (10^21nvt) --                                             *
C                   D: Neutron fluence                                 *
C    T (C)          D: Temperature in dense PyC layer                  *
C    EC             D: Apparent creep strains                          *
C  Local variable description                                          *
C    EEC            D: Effective creep strain                          *
C                                                                      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
	SUBROUTINE SSICREEP_PYC(FLUENCE, D, T, EC, CREEP, NU_C)
C
	DOUBLE PRECISION FLUENCE, D, T, EC, CREEP, NU_C
	DOUBLE PRECISION EEC
	DIMENSION EC(1:2)
C
	EEC = DSQRT(2.0D0/3.0D0)*ABS(EC(1)+EC(2))
C
	CREEP = (1.0D0 + (1.9D0 - D)*2.38D0)*
     &        (2.193D-4 - 4.85D-7*T + 4.0147D-10*T**2)
C
C	IF (EEC.LE.0.013D0) THEN
C	  NU_C = 0.5D0 - EEC/0.13D0
C	ELSE
C	  NU_C = 0.4D0
C	END IF
      IF (FLUENCE.LE.0.3D0) THEN
        NU_C = 0.5D0 - FLUENCE/3.0D0
      ELSE
        NU_C = 0.4D0
      ENDIF
C for comparison with INEEL's results, fix NU_C to 0.5
c	NU_C = 0.5D0
	RETURN
	END
C                                                                      *
C***********************************************************************
C
C
C   5) PyC strength
C
C***********************************************************************
C                                                                      *
C  Subroutine STRENGTH_PYC                                             *
C                                                                      *
C    Calculate mean fracture strength and Weibull modulus for PyC.     *
C                                                                      *
C  Requirements:                                                       *
C    Density:        : <2.00 g/cm^3                                    *
C                                                                      *
C  Effects:                                                            *
C    Calculate the initial characteristic strength, mean fracture      *
C  strength at any time and Weibull modulus, return the results via    *
C  parameter list.                                                     *
C                                                                      *
C  Actual variable description                                         *
C    LAYER          C: The layer to be examined (IPYC or OPYC)         *
C    FLAG           C: Indicates if the user wants initial strength    *
C                      or irradiated strength                          *
C    D (g/cm^3)     D: Density of dense pyrocarbon                     *
C    BAF0           D: Initial BAF                                     *
C    FLU (10^21nvt) --                                                 *
C                   D: Neutron fluence                                 *
C    T (C)          D: Temperature in dense PyC layer                  *
C    SIGT (MPa)     D: Tangential stress distribution in layers        *
C    SIGMA0 (MPa.meter^3/m) (m is Weibull modulus here)                *
C                   D: Initial characteristic strength of PyC          *                                                              *
C    SIGF (MPa)     D: Mean fracture strength of PyC                   *
C    M              D: Weibull modulus of PyC                          *
C                                                                      *
C  Local variable description                                          *
C    R1 - R5 (um)   D: Layer radii                                     *
C    SIGMA0I (MPa.meter^3/m) (m is Weibull modulus here)               *
C                   D: Irradiated characteristic strength              *
C    SIGTBAR (MPa)  D: Average tangential stress in certain PyC layer  *
C    SIGTIBAR (MPa) D: Average tangential stress of each integration   *
C                      shell in certain PyC layer                      *
C    V (m^3)        D: The volume of certain PyC layer                 *
C    VI (m^3)       D: The volume of each integration shell in         *
C                      certain PyC layer                               *
C    INTEGRAL (MPA^m.meter^3)                                          *
C                   D: An intergral in the calculation of mean fracture*
C                      strength                                        *
C    RI (um)        D: The inner radius of the integration shell       *
C    RO (um)        D: The outer radius of the integration shell       *
C    DR (um)        D: The thickness of the intergration shell         *
C    NDIV?          I: The number of sampling points in each layer     *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
	SUBROUTINE STRENGTH_PYC(LAYER, FLAG, D, BAF0, FLU, T, SIGT,
     &					    SIGMA0, SIGF, M)
C
	PARAMETER ( PIE = 3.1415926535897932385 D0 )
	PARAMETER ( SCALE = 1.0D-18)
	DOUBLE PRECISION D, BAF0, FLU, T, SIGT, SIGMA0, SIGF, M
	DOUBLE PRECISION R1, R2, R3, R4, R5
	DOUBLE PRECISION SIGMA0I, SIGTBAR, SIGTIBAR, V, VI, INTEGRAL
	DOUBLE PRECISION RI, RO, DR
	INTEGER NDIVI, NDIVS, NDIVO
	CHARACTER*4 LAYER
	CHARACTER*3 FLAG
C
	DIMENSION SIGT(1:)
C
      COMMON /PAR_R/ R1, R2, R3, R4, R5         !particle geometry
	COMMON /PAR_DIV/ NDIVI, NDIVS, NDIVO
C
	IF(FLAG.EQ.'INI') THEN
	  M = 9.5D0
C       The characteristic strength of PyC at room temperature and density
C       of 1.90 (Weibull modulus equals to 9.5).
	  SIGMA0 = 154.46D0*BAF0**2 - 141.1D0*BAF0
C       Adjust characteristic strength for specified density 
	  IF(D.GE.1.8D0) THEN
		CONTINUE
	  ELSE IF(D.GE.1.1D0) THEN
		SIGMA0 = SIGMA0*(1.241D0*D - 1.234D0)
	  ELSE IF(D.GT.1.0D0) THEN
		SIGMA0 = SIGMA0*(0.664D0*D - 0.577D0)
	  ELSE IF(D.EQ.1.0D0) THEN
		SIGMA0 = SIGMA0*0.0667D0
	  END IF
	  IF(LAYER .EQ. 'IPYC') THEN
	    V = 4.0D0*PIE*SCALE*(R3**3.0 - R2**3.0)/3.0D0
	  ELSE IF(LAYER .EQ. 'OPYC') THEN
	    V = 4.0D0*PIE*SCALE*(R5**3.0 - R4**3.0)/3.0D0
	  END IF
C       Calculate mean fracture strength based on the characteristic strength
	  SIGF = SIGMA0/((2.0D0*V)**(1.0D0/M))
	ELSE IF(FLAG.EQ.'IRR') THEN
C       Calculate characteristic strength at temperature T and irradiation FLU
	  SIGMA0I = SIGMA0*DSQRT((1.0D0 + 0.23D0*FLU)
     &                         *(1.0D0 + 1.5D-4*(T-20.D0)))
	  SIGTBAR = 0.0D0
	  INTEGRAL = 0.0D0
	  IF(LAYER .EQ. 'IPYC') THEN
	    DO 4001 I = 1, NDIVI+2
	      SIGTBAR = SIGTBAR + SIGT(I)
4001      CONTINUE
          SIGTBAR = SIGTBAR/FLOAT(NDIVI+2)
		V = 4.0D0*PIE*SCALE*(R3**3.0 - R2**3.0)/3.0D0
		DR = (R3 - R2)/(NDIVI + 1)
		RI = R2
		RO = RI + DR
	    DO 4002 I = 1, NDIVI + 1
	      SIGTIBAR = (SIGT(I)+SIGT(I+1))/2.0D0
C		  If the stress is negative, intergration can not be performed because
C		  SIGTIBAR**M can be a complex value.
		  IF(SIGTIBAR .LE. 0.0D0) THEN
			SIGF = SIGMA0I/((2.0D0*V)**(1.0D0/M))
			RETURN
		  END IF
	      VI = 4.0D0*PIE*SCALE*(RO**3.0 - RI**3.0)/3.0D0
		  INTEGRAL = INTEGRAL + (SIGTIBAR**M)*VI
		  RI = RO
	      RO = RI + DR
4002      CONTINUE
		SIGF = SIGMA0I*SIGTBAR/((2.0D0*INTEGRAL)**(1.0D0/M))
	  ELSE IF(LAYER .EQ. 'OPYC') THEN
	    DO 4003 I = 1, NDIVO+2
	      SIGTBAR = SIGTBAR + SIGT(NDIVI+NDIVS+4+I)
4003      CONTINUE
		SIGTBAR = SIGTBAR/FLOAT(NDIVO+2)
		V = 4.0D0*PIE*SCALE*(R5**3.0 - R4**3.0)/3.0D0
		DR = (R5 - R4)/(NDIVO + 1)
		RI = R4
		RO = RI + DR
		DO 4004 I = 1, NDIVO + 1
	      SIGTIBAR = (SIGT(NDIVI+NDIVS+4+I)+SIGT(NDIVI+NDIVS+4+I+1))
     &				 /2.0D0
C		  If the stress is negative, intergration can not be performed because
C		  SIGTIBAR**M can be a complex value.
		  IF(SIGTIBAR .LE. 0.0D0) THEN
			SIGF = SIGMA0I/((2.0D0*V)**(1.0D0/M))
			RETURN
		  END IF
		  VI = 4.0D0*PIE*SCALE*(RO**3.0 - RI**3.0)/3.0D0
	      INTEGRAL = INTEGRAL + (SIGTIBAR**M)*VI
		  RI = RO
	      RO = RI + DR
4004		CONTINUE
		SIGF = SIGMA0I*SIGTBAR/((2.0D0*INTEGRAL)**(1.0D0/M))
	  END IF
	END IF
	RETURN
	END
C                                                                      *
C***********************************************************************
C
C
C 2. Silicon Carbide
C
C
C   1) Young's modulus (MPa)
C
C
C***********************************************************************
C                                                                      *
C  Function E_SIC                                                      *
C                                                                      *
C    Young's modulus vs Temperature of SiC                             *
C                                                                      *
C  Actual variable description                                         *
C    E_SIC (MPa)    D: Young's modulus of SiC                          *
C    T (C)          D: Temperature of SiC layer                        *
C  Local variable description                                          *
C    E_T_SIC        D: Array of E_SIC dependence on temperature        *
C    INDEX          D: Array index                                     *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      FUNCTION E_SIC(T)
C
      DOUBLE PRECISION E_SIC, T
	DOUBLE PRECISION E_T_SIC(1:4,-1:0)
	INTEGER INDEX
C
	DATA E_T_SIC /  25.0D0,   940.0D0,  1215.0D0,  1400.0D0,
     &			  4.2718D5, 3.74816D5, 3.39677D5, 2.71466D5/
C
      CALL LOCATE_ARRAY(E_T_SIC(:,-1),4,1,T,INDEX)
	IF(INDEX.LT.4) THEN
	  E_SIC = ((T-E_T_SIC(INDEX,-1))*E_T_SIC(INDEX+1,0)+
     &           (E_T_SIC(INDEX+1,-1)-T)*E_T_SIC(INDEX,0))/
     &          (E_T_SIC(INDEX+1,-1)-E_T_SIC(INDEX,-1))
	ELSE
	  E_SIC = ((T-E_T_SIC(INDEX-1,-1))*E_T_SIC(INDEX,0)+
     &           (E_T_SIC(INDEX,-1)-T)*E_T_SIC(INDEX-1,0))/
     &          (E_T_SIC(INDEX,-1)-E_T_SIC(INDEX-1,-1))
	END IF
C  fix it temporarily
C	E_SIC = 3.0D5
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C   5) SiC strength
C
C***********************************************************************
C                                                                      *
C  Subroutine STRENGTH_SIC                                             *
C                                                                      *
C    Calculate mean fracture strength and Weibull modulus for SiC.     *
C                                                                      *
C  Notes: Due to inconsistency of data, the effects of irradiation     *
C         and temperature on strength are not included for now.        *
C         However, their entries in the parameter list are reserved    *
C         for future improvements.                                     *
C                                                                      *
C  Requirements:                                                       *
C                                                                      *
C  Effects:                                                            *
C    Calculate the initial characteristic strength, mean fracture      *
C  strength at any time and Weibull modulus, return the results via    *
C  parameter list.                                                     *
C                                                                      *
C  Actual variable description                                         *
C    FLAG           C: Indicates if the user wants initial strength    *
C                      or irradiated strength                          *
C    FLU (10^21nvt) --                                                 *
C                   D: Neutron fluence                                 *
C    T (C)          D: Temperature in SiC layer                        *
C    SIGT (MPa)     D: Tangential stress distribution in layers        *
C    SIGMA0 (MPa.meter^3/m) (m is Weibull modulus here)                *
C                   D: Initial characteristic strength of SiC          *                                                              *
C    SIGF (MPa)     D: Mean fracture strength of SiC                   *
C    M              D: Weibull modulus of SiC                          *
C                                                                      *
C  Local variable description                                          *
C    R1 - R5 (um)   D: Layer radii                                     *
C    SIGTBAR (MPa)  D: Average tangential stress in SiC layer          *
C    SIGTIBAR (MPa) D: Average tangential stress of each integration   *
C                      shell in SiC layer                              *
C    V (m^3)        D: The volume of SiC layer                         *
C    VI (m^3)       D: The volume of each integration shell in         *
C                      SiC layer                                       *
C    INTEGRAL (MPA^m.meter^3)                                          *
C                   D: An intergral in the calculation of mean fracture*
C                      strength                                        *
C    RI (um)        D: The inner radius of the integration shell       *
C    RO (um)        D: The outer radius of the integration shell       *
C    DR (um)        D: The thickness of the intergration shell         *
C    NDIV?          I: The number of sampling points in each layer     *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
	SUBROUTINE STRENGTH_SIC(FLAG, FLU, T, SIGT, SIGMA0, SIGF, M)
C
	PARAMETER ( PIE = 3.1415926535897932385 D0 )
	PARAMETER ( SCALE = 1.0D-18)
	DOUBLE PRECISION FLU, T, SIGT, SIGMA0, SIGF, M
	DOUBLE PRECISION R1, R2, R3, R4, R5
	DOUBLE PRECISION SIGTBAR, SIGTIBAR, V, VI, INTEGRAL
	DOUBLE PRECISION RI, RO, DR
	INTEGER NDIVI, NDIVS, NDIVO
	CHARACTER*3 FLAG
C
	DIMENSION SIGT(1:)
C
      COMMON /PAR_R/ R1, R2, R3, R4, R5         !particle geometry
	COMMON /PAR_DIV/ NDIVI, NDIVS, NDIVO
C
	IF(FLAG.EQ.'INI') THEN
	  M = 6.0D0
C       The characteristic strength of SiC
	  SIGMA0 = 9.64D0
	  V = 4.0D0*PIE*SCALE*(R4**3.0 - R3**3.0)/3.0D0
C       Calculate mean fracture strength based on the characteristic strength
	  SIGF = SIGMA0/((2.0D0*V)**(1.0D0/M))
	ELSE IF(FLAG.EQ.'IRR') THEN
	  SIGTBAR = 0.0D0
	  INTEGRAL = 0.0D0
	  DO 4011 I = 1, NDIVS+2
	    SIGTBAR = SIGTBAR + SIGT(NDIVI+2+I)
4011    CONTINUE
        SIGTBAR = SIGTBAR/FLOAT(NDIVS+2)
	  V = 4.0D0*PIE*SCALE*(R4**3.0 - R3**3.0)/3.0D0
	  DR = (R4 - R3)/(NDIVS+1)
	  RI = R3
	  RO = RI + DR
	  DO 4012 I = 1, NDIVS+1
	    SIGTIBAR = (SIGT(NDIVI+2+I)+SIGT(NDIVI+2+I+1))/2.0D0
C		If the stress is negative, intergration can not be performed because
C		SIGTIBAR**M can be a complex value.
		IF(SIGTIBAR .LE. 0.0D0) THEN
		  SIGF = SIGMA0/((2.0D0*V)**(1.0D0/M))
		  RETURN
		END IF
	    VI = 4.0D0*PIE*SCALE*(RO**3.0 - RI**3.0)/3.0D0
		INTEGRAL = INTEGRAL + (SIGTIBAR**M)*VI
		RI = RO
	    RO = RI + DR
4012    CONTINUE
	  SIGF = SIGMA0*SIGTBAR/((2.0D0*INTEGRAL)**(1.0D0/M))
	END IF
	RETURN
	END
C                                                                      *
C***********************************************************************
C
C
C 3. Apparent creep strain of PyC
C
C
C***********************************************************************
C                                                                      *
C  Subroutine EPI_C                                                    *
C                                                                      *
C    Calculate apparent creep strains in pyrocarbon                    *
C                                                                      *
C  Actual variable description                                         *
C    SIGR (MPa)     D: Radial stress distribution                      *
C    SIGT (MPa)     D: Tangential stress distribution                  *
C    CREEP ((strain/MPa.10**21 nvt) --                                 *
C                   D: Irradiation creep coefficient                   *
C    CNU            D: Poisson's ratio for creep                       *
C    DFLU (10^21nvt)D: Incremental fluence gained in last step         *
C    EC             D: Array of apparent creep strain                  *
C                      EC(1): Radial component                         *
C                      EC(2): Tangential component                     *
C  Local variable description                                          *
C    DEC            D: Incremental creep strain from last step         *
C    SIGR_BAR (MPa) D: Average radial stress in certain layer          *
C    SIGT_BAR (MPa) D: Average tangential stress in certain layer      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
	SUBROUTINE EPI_C(LAYER, SIGR, SIGT, CREEP, CNU, DFLU, EC)
C
	DOUBLE PRECISION SIGR, SIGT, CREEP, CNU, DFLU, EC
	DOUBLE PRECISION DEC, SIGR_BAR, SIGT_BAR
	INTEGER   NDIVI, NDIVS, NDIVO
	INTEGER   I
	CHARACTER*(*) LAYER
	DIMENSION SIGR(1:), SIGT(1:)
	DIMENSION EC(1:2)
	COMMON /PAR_DIV/ NDIVI, NDIVS, NDIVO
C
	SIGR_BAR = 0.0D0
	SIGT_BAR = 0.0D0
	IF(LAYER.EQ.'IPYC') THEN
	  DO 2711 I = 1, NDIVI+2
	    SIGR_BAR = SIGR_BAR + SIGR(I)
		SIGT_BAR = SIGT_BAR + SIGT(I)
2711	  CONTINUE
	  SIGR_BAR = SIGR_BAR/(NDIVI+2)
	  SIGT_BAR = SIGT_BAR/(NDIVI+2)
	  DEC = CREEP*(SIGR_BAR-2.0D0*CNU*SIGT_BAR)*DFLU
C	  EC(1) = DEC
	  EC(1) = EC(1) + DEC
	  DEC = CREEP*((1.0D0-CNU)*SIGT_BAR-CNU*SIGR_BAR)*DFLU
C	  EC(2) = DEC
	  EC(2) = EC(2) + DEC
	ELSE IF(LAYER.EQ.'OPYC') THEN
	  DO 2712 I = NDIVI+NDIVS+5, NDIVI+NDIVS+NDIVO+6
		SIGR_BAR = SIGR_BAR + SIGR(I)
		SIGT_BAR = SIGT_BAR + SIGT(I)
2712	  CONTINUE
	  SIGR_BAR = SIGR_BAR/(NDIVO+2)
	  SIGT_BAR = SIGT_BAR/(NDIVO+2)
	  DEC = CREEP*(SIGR_BAR-2.0D0*CNU*SIGT_BAR)*DFLU
C	  EC(1) = DEC
	  EC(1) = EC(1) + DEC
	  DEC = CREEP*((1.0D0-CNU)*SIGT_BAR-CNU*SIGR_BAR)*DFLU
C	  EC(2) = DEC
	  EC(2) = EC(2) + DEC	  
	ELSE  !No creep if the layer is not pyrocarbon, or used to clear EC.
	  EC(1) = 0.0D0
	  EC(2) = 0.0D0
	END IF
C    When doing some study on the effects on apparent creep strain, it is
C    turned off sometimes.
C	EC(1) = 0.0D0
C	EC(2) = 0.0D0
	RETURN
	END
C                                                                      *
C***********************************************************************
C