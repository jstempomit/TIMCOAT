C***********************************************************************
C                                                                      *
C                        PROGRAM TIMCOAT v.2                           *
C                                                                      *
C  Author:                    Jing Wang                                *
C                 Department of Nuclear Engineering                    *
C               Massachusetts Institute of Technology                  *
C                         Tel: 617-253-3247                            *
C                      Email: jingwang@mit.edu                         *
C                                                                      *
C  Date: February, 2004                                                *
C                                                                      *
C  General Description                                                 *
C     The program name is a acronym for The Integrated Model of        *
C  COATed particle fuel performance.                                   *
C     This FORTRAN program is the implementation of the fuel           *
C  performance model for High Temperature Gas-cooled Reactors (HTGRs). *
C  It samples TRISO fuel particles by Monte Carlo sampling method      *
C  and allows the particles to be cycled in a steady-state pebble bed  *
C  reactor core. By this method, the irradiation history of each       *
C  particle is created and then mechanical and chemical analyses       *
C  are performed to predict the behavior of fuel particles, then       *
C  particle failures are evulated based on fuel failure models.        *
C                                                                      *
C  Highlights                                                          *
C  1.  Fracture mechanics is used to develop a pyrocarbon crack        *
C      induced particle failure model to work together with the        *
C      classical pressure vessel model to predict fuel failure         *
C  2.  The model consists of correlated mechancial, chemical           *
C      and thermal analyses, which makes it a well-intergrated model.  *
C  3.  The mechanical model considers irradiation-induced dimensional  *
C      change and creep in PyC layers, and pressure load in all three  *
C      structural layers in a TRISO particle. Analytical solution is   *
C      used here because of speed-demanding MC sampling. Thermal       *
C      stresses due to thermal expansion is dealt with in a separate   *
C      subroutine. Most of the data for material properties come from  *
C      the GA report -- CEGA002820, 1993.                              *
C  4.  This fuel performance model can simulate the irradiation        *
C      conditions in real pebble bed reactor environments. Also, it    *
C      can evaluate the performance of fuel particles with arbitrary   *
C      power histories.                                                *
C                                                                      *  
C                          Modified by:Jane Diecker                    * 	
C                    Department of Nuclear Engineering                 *
C                  Massachusetts Institute of Technology               *
C                         Email: janed@mit.edu                         *
C                                                                      *
C                            Date: May, 2005                           *	  
C                                                                      *
C  Changes:                                                            *
C  1.  The model now includes fission product attack of the SiC        *	                                                    
C  2.  The model now accounts for failure due to the amoeba effect     *
C	         	                                                       *
C***********************************************************************
C
C  Description of main program variables
C
C    In the following table, the column labeled 'Declared' has values:
C           
C           R*8      : REAL*8 (Double Precision)
C           I*4      : INTEGER*4
C           LOG      : LOGICAL
C           C*n      : CHARACTER*n
C
C    In the following table, the column labeled 'Type' is coded:
C           
C             P      : FORTRAN PARAMETER; typically a physical constant
C                      or a standard array dimension
C             C      : Code version or control parameters
C             I      : Input from the user
C             V      : Important program variables
C             S      : Batch statistical or stress histrogram variables
C             M      : Miscellaneous variables, usually intermediate
C                      values in a calculation or tally, loop 
C                      counters, etc.
C             F      : Function names; values returned by a function
C           
C           
C  Variable  Declared  Type                  Description
C  --------  --------  ----  -------------------------------------------
C
C  AVGT        R*8      M  : Average temperature in particle (K)
C  AWC         R*8      P  : Molecular Weight of Carbon [gm/gm-mol]
C  AWO         R*8      P  : Molecular Weight of Oxygen [gm/gm-mol]
C  AWU         R*8      V  : Atomic weight of Uranium mix [gm/gm-mol]
C  AWU235      R*8      P  : Molecular Weight of U235 [gm/gm-mol]
C  AWU238      R*8      P  : Molecular Weight of U238 [gm/gm-mol]
C  BDEN        R*8      V  : Actual buffer density [g/cc]
C  BLOCKMAP    I*4      I  : Specifies number of layers in every channel
C  BLOCKS      R*8      I  : Defines blocks of every channel in MPBR core (from VSOP outputs)
C  BPRATE      R*8      V  : Burnup rate [MW/T/day]
C  BU_CONV     R*8      P  : The factor that convert fima to MWd/T
C  BUFFD       R*8      I  : Buffer mean density [g/cc]
C  BUFFDVAR    R*8      I  : Buffer density std deviation [g/cc]
C  BUFFT       R*8      P  : Theoretical density of PyC [g/cm^3]
C  BUFFTHK     R*8      I  : Nominal Buffer thickness [microns]
C  BUFFVAR     R*8      I  : Buffer std deviation [microns]
C  BURNUP      R*8      V  : Current value of burnup [% FIMA]
C  CDATE       C*8      M  : Date format MM/DD/YY
C  CHANNELS    R*8      I  : Defines channels of MPBR core (from VSOP output)
C  CL          R*8      M  : Critical crack length of SiC (used to calculate
C                            KICSIC if it's not a given input)
C  CONFIDIF    R*8      S  : +/- 1 sigma limit on PROBIF
C  CONFIDOF    R*8      S  : +/- 1 sigma limit on PROBOF
C  CONFIDPF    R*8      S  : +/- 1 sigma limit on PROBPF
C  CONFIDSF    R*8      S  : +/- 1 sigma limit on PROBSF
C  COMP        C*3      C  : Compiler type
C  CORE_HEIGHT R*8      I  : Height of reactor core [m]
C  CORE_RADIUS R*8      I  : Radius of reactor core [m]
C  COREMODEL   I*4      V  : Indicating whether new or old VSOP core model is employed
C  CTIME       C*8      M  : Time format HH/MM/SS
C  CURAT       R*8      I  : Carbon to Uranium atom ratio
C  CV          R*8      P  : Numerical constant (4/3) Pi X 10^-12
C  DCORR       R*8      M  : Distance of fission product corrosion of SiC (microm)
C  DEBUG       LOG      I  : Flag for debugging
C  D_FLUENCE   R*8      M  : Incremental neutron fluence gained in one step [10^21nvt]
C  DF          R*8      M  : Incremental fluence since last step [10^21nvt].
C						   It is used in common block /CRACKED_PYC/.
C  DIFFUSION   LOG      I  : Flag indicating whether to take into account
C                            diffusion in gas release model or not
C  DPD         R*8      M  : Distance of Pd crack penetration 
C  DT          R*8      I  : Time increment in the reactor core (sec)
C  DUMMY       R*8      M  : Useful dummy variable for subroutine calls
C  E_IPYCREEP  R*8      V  : Apparent creep strain in IPyC
C  E_IPYCREEP0 R*8      V  : Back up of E_IPYCREEP at certain reference point
C  E_OPYCREEP  R*8      V  : Apparent creep strain in OPyC
C  E_OPYCREEP0 R*8      V  : Back up of E_OPYCREEP at certain reference point
C  E_PYC       R*8      F  : Function that returns the PyC elastic modulus
C  E_SIC       R*8      F  : Function that returns the SiC elastic modulus
C  ELAPTIME    R*4      F  : Function that returns the number of seconds
C                            since midnight of the current day
C  EOLBUP      R*8      I  : End of Life Burnup [% FIMA]
C  EOLFLU      R*8      I  : End of Life Fast Fluence [10^21 nvt]
C  EPIR        R*8      V  : Array of radial strain distribution
C  EPIR0       R*8      V  : Array of radial residual strain distribution after each cycle
C  EPIRCP      R*8      V  : Array recording the elastic radial strains of
C					       fully relaxed PyC layers at point of cracking.
C					       The symbol means 'Epsilon R of C prime'
C  EPIT        R*8      V  : Array of tangential strain distribution
C  EPIT0       R*8      V  : Array of tangential residual strain distribution after each cycle
C  EPITCP      R*8      V  : Array recording the elastic tangential strains
C					       of fully relaxed PyC layers at point of cracking.
C					       The symbol means 'Epsilon T of C prime'
C  ERROR_CODE  I*4      C  : Code symbolizes the degree of error:
C                            2 means this is an unconditionally fatal error.
C                            1 means this is a recoverable error.
C                            0 means this is a warning message only.
C  FAIL        LOG      V  : Indicates whether a particle is failed or not
C  FAILTYPE    I*4      V  : array counting the number of each type of failure
C                            FAILTYPE(0): fracture induced by IPyC cracking
C                            FAILTYPE(1): overpressure rupture
C                            FAILTYPE(2): amoeba effect
C  FAILUREPATH C*42     C  : Records the path a failed particle undertake
C  FF          C*1      M  : ASCII character
C  FFLUX       R*8      V  : Current fast flux (/cm^2.sec)
C  FMODE       I*4      V  : Failure mode of the current failed particle
C  FOPEN1      LOG      M  : Flags whether the input file INPFILE1 is open (0/1)
C  FOPEN2      LOG      M  : Flags whether the input file INPFILE2 is open (0/1)
C  FOPEN3      LOG      M  : Flags whether the input file INPFILE3 is open (0/1)
C  FOPEN4      LOG      M  : Flags whether the input file INPFILE4 is open (0/1)
C  FLUENCE     R*8      V  : Fast fluence vs time [nvt]
C  FLUENCE_R   R*8      V  : Fast fluence gained in every fuel cycle [nvt]
C  FLUX        R*8      F  : Function that calculates the current neutron flux
C  FUELTYPE    C*3      V  : Type of fuel evaluated ('UCO','UO2', or '(Th,U)O2')
C  GAUSSIAN    R*8      F  : Function that samples a Gaussian variate
C                            from a given mean and standard deviation
C  HCARD       R*8      V  : History card for particles, the columns of which 
C                            represent the following variables: (HCARD(*,C*))
C                            C1: OPERTIME	(second)
C                            C2: PEBBLE_R (cm)
C                            C3: PEBBLE_Z (cm)
C                            C4: FLUENCE/1.0D21 (10^21nvt)
C                            C5: FLUENCE_R/1.0D21 (10^21nvt)
C                            C6: QPPP	    (W/m3)
C                            C7: BURNUP   (fima)
C                            C8: T_HE     (C)
C                            C9-C14: T_PARTICLE(0:5) (C)
C                            C15: PRESS   (MPa)
C                            C16: SRDOT_IPYC (dimensionless)
C                            C17: STDOT_IPYC (dimensionless)
C                            C18: SRDOT_OPYC (dimensionless)
C                            C19: STDOT_OPYC (dimensionless)
C  HCARDB      R*8      V  : History card for particles under constant power history;
C                            the structure is the same as HCARD
C  HISTGRMI    I*4      S  : Array recording the IPyC failure histogram
C                            HISTGRM*(:,1) : Time at which the failure happens
C                            HISTGRM*(:,2) : Tangential stress at which the failure happens
C                            HISTGRM*(:,3) : Fluence at which the failure happens
C                            HISTGRM*(:,4) : Burnup at which the failure happens
C  HISTGRMO    I*4      S  : Array recording the OPyC failure histogram
C  HISTGRMP    I*4      S  : Array recording the particle failure histogram
C  HISTGRMS    I*4      S  : Array recording the SiC failure histogram
C  HISTOGRAM   LOG      V  : Flag for ouptput histogram
C  HSIGR       R*8      V  : Array of radial stress evolution [MPa]
C  HSIGRB      R*8      V  : Array of radial stress evolution for constant power [MPa]
C  HSIGT       R*8      V  : Array of tangential stress evolution [MPa]
C  HSIGTB      R*8      V  : Array of tangential stress evolution for constant power [MPa]
C  I           I*4      M  : General purpose loop counter variable
C  IDAT1       I*4      C  : Lun of input data file 'INPFILE1.DAT'
C  IDAT2       I*4      C  : Lun of input data file 'INPFILE2.DAT'
C  IDAT3       I*4      C  : Lun of input data file 'INPFILE3.DAT'
C  IDAT4       I*4      C  : Lun of input data file 'INPFILE4.DAT'
C  IERR        I*4      C  : Lun of error file 'ERROR.MSG'
C  IEVOLR      I*4      C  : Lun of output file of radial stress evolution
C  IEVOLT      I*4      C  : Lun of output file of tangential stress evolution
C  IHIS        I*4      C  : Lun of histogram output file
C  IKEY        I*4      C  : Lun of data entered by keyboard
C  IMSG        I*4      C  : Error code returned by IOSTAT during READ
C  INDEX       I*4      M  : Array index
C  INDEX2      I*4      M  : Array index
C  INITSEED    I*4      F  : Function to init ISEED with system clock
C  INPFILE1    C*8      I  : User-input filename ('.dat' filetype omitted)
C  INPFILE2    C*8      I  : Input file of power distribution
C  INPFILE3    C*8      I  : Input file of power history
C  INPFILE4    C*8      I  : Input file of fast flux history
C  IOUT        I*4      C  : Lun of general output file 'OSPEC.OUT'
C  IOS         I*4      C  : Store IOSTATS of READ or WRITE statements
C  IPYCALPHA   R*8      V  : Thermal expansion coeff of IPyC [per K]
C  IPYCBAF0    R*8      V  : Sampled initial Bacon Anisotropy Factor of IPyC [dimensionless]
C  IPYCBAF0I   R*8      I  : Initial mean Bacon Anisotropy Factor of IPyC [dimensionless]
C  IPYCBAFI    R*8      V  : Irradiated Bacon Anisotropy Factor of IPyC [dimensionless]
C  IPYCBAFVAR  R*8      I  : Standard deviation of IPYCBAF0I
C  IPYCCNU     R*8      V  : Creep Poisson's ratio of IPyC
C  IPYCCRATE   R*8      I  : Coating rate of IPyC layer [um/min]
C  IPYCD       R*8      I  : IPyC density [g/cc]
C  IPYCE       R*8      V  : Elastic modulus of IPyC [MPa]
C  IPYCF       R*8      I  : Initial IPyC characteristic strength [MPa.meter^3/m]
C  IPYCFAIL    I*4      S  : Count of number of particles with IPyC layer failure
C  IPYCFAIL0   I*4      S  : Count of number of particles with IPyC layer failure
C						   during irradiation one group prior to IPYCFAIL. (For statistical purpose)
C  IPYCFAILED  I*2      V  : Indicates if IPyC of current particle failed
C  IPYCLC      R*8      V  : Crystallite size of IPyC (A)
C  IPYCM       R*8      I  : IPyC Weibull Modulus
C  IPYCNU      R*8      V  : Poisson's ratio of IPyC
C  IPYCREEP    R*8      V  : Creep coefficient of IPyC [strain per MPa-10**21 nvt]
C  IPYCTHK     R*8      I  : Nominal IPyC thickness [microns]
C  IPYCVAR     R*8      I  : IPyC std deviation [microns]
C  IRRHISTRY   R*8      I  : Array for irradiation history of particles, used in
C						   the branch PSWITCH = 2.
C						   IRRHISTRY(,1) : Ellapse of irradiation test [days]
C						   IRRHISTRY(,2) : Full power days [days]
C						   IRRHISTRY(,3) : Irradiation temperature [degree C]
C						   IRRHISTRY(,4) : Neutron fast fluence [10^21nvt]
C						   IRRHISTRY(,5) : Burnup [% fima]
C  IRRTIME     R*8      I  : Effective irradiation time [days]
C  ISEED       I*4      I  : Current random number seed
C                            (user inputs only the first seed)
C  ISWR(0:NDEG)R*8      V  : Array of length NDEG+1 storing coefficients of
C                            the polynomial of IPyC radial swelling rate
C  ISWT(0:NDEG)R*8      V  : Array of length NDEG+1 storing coefficients of
C                            the polynomial of IPyC tangential swelling rate
C  ITERM       I*4      C  : Logical unit number (lun) of terminal
C  ITEST       I*4      C  : Lun of test output file 'TEST.DAT'
C  J           I*4      M  : General purpose loop counter variable
C  K_PM        R*8      V  : Thermal conductivity of pebble matrix [W/m.K]
C  K_PFM       R*8      V  : Thermal conductivity of pebble fuel zone matrix [W/m.K]
C  K_PFZ       R*8      V  : Thermal conductivity of pebble fuel zone [W/m.K]
C  KDEN        R*8      V  : Actual kernel density [g/cc]
C  KERND       R*8      I  : Kernel mean density [g/cc]
C  KERNDIA     R*8      I  : Nominal kernel diameter [microns]
C  KERNDVAR    R*8      I  : Kernel density std deviation [g/cc]
C  KERNT       R*8      I  : Kernel theoretical density [g/cc]
C  KERNVAR     R*8      I  : Kernel diameter std deviation [microns]
C  KI1         R*8      V  : Stress intensity factor from IPyC crack [MPa.um^1/2]
C  KI2         R*8      V  : Stress intensity factor from OPyC crack [MPa.um^1/2]
C  KICSIC      R*8      M  : Actual critical stress intensity factor of SiC [MPa*um1/2]
C  KICSICI     R*8      M  : Unirradiated critical stress intensity factor of SiC [MPa*um1/2]
C  KIIPYC      R*8      M  : Stress intensity factor in IPyC layer
C  KIOPYC      R*8      M  : Stress intensity factor in OPyC layer
C  KMC         R*8      M  : Kernel migration coefficient (m^2K/s)
C  LENGTH_OF_FILE  R*8  M  : Length of input files
C  MACH        C*3      C  : Machine type
C  MACHTIME    R*8      V  : Program running time (real time)
C  MCODE       C*4      C  : Indicates the type of analysis to perform
C							'ISO3': full three-layer analysis
C							'IS2' : IPyC/SiC two-layer analysis
C							'SO2' : SiC/OPyC two-layer analysis
C							'S1'  : SiC single-layer analysis
C  MD          R*8      M  : Kernel migration distance 
C  MF_HE       R*8      I  : Mass flow rate of Helium [kg/s]
C  MOLES       R*8      V  : Moles of gas in particle [gm-mol]
C  MOLESU      R*8      V  : Moles of uranium in kernel [gm-mol]
C  N           I*4      M  : Current case (particle) number  (1..NCASES)
C  NAXIAL      I*4      I  : Number of axial divisions in the reactor core specified by VSOP model
C  NBLOCK      I*4      V  : Number of blocks/layers in the channel
C  NBURP       I*4      I  : Frequency of user status messages (burps)
C                            during sampling of particles
C  NCASES      I*4      I  : Number of Monte Carlo particle cases to run:
C                            Maximum Size is 2 billion (2,147,483,647)
C  NCHANNEL    I*4	  I  : Number of channels in reactor core specified by VSOP model
C  NCHAR       I*2      I  : Counts number of characters in FAILUREPATH
C  NDIV        I*4      P  : Number of divisions for stress distributions
C  NDIVI       I*2      V  : Number of divisions in IPyC
C  NDIVO       I*2      V  : Number of divisions in OPyC
C  NDIVS       I*2      V  : Number of divisions in SiC
C  NDUMPED     I*4      V  : Number of particles who exceed end-of-life burnup
C  NHIS        I*4      P  : Number of divisions for histograms
C  NLAYER      I*4      I  : Number of layers/blocks in the reactor core specified by VSOP model
C  NOMINAL     LOG      I  : Flag to override statistical sampling (0/1)
C  NPEBBLE     I*4      I  : Number of pebbles in the reactor core
C  NPARTICLE   I*4      I  : Average number of particles in each pebble
C  NSTEPINCYCLE  I*4    M  : Number of steps in each irradiation cycle
C						   (This variable is used particularly in PSWITCH = 2)
C  OPERTIME    R*8      V  : Reactor operation time (virtual time ) [seconds]
C  OPYCALPHA   R*8      V  : Thermal expansion coeff of OPyC [per K]
C  OPYCBAF0    R*8      V  : Sampled initial Bacon Anisotropy Factor of OPyC [dimensionless]
C  OPYCBAF0I   R*8      I  : Initial mean Bacon Anisotropy Factor of OPyC [dimensionless]
C  OPYCBAFI    R*8      V  : Irradiated Bacon Anisotropy Factor of OPyC [dimensionless]
C  OPYCBAFVAR  R*8      I  : Standard deviation of OPYCBAF0I
C  OPYCCNU     R*8      V  : Creep Poisson's ratio of OPyC
C  OPYCCRATE   R*8      I  : Coating rate of OPyC layer [um/min]
C  OPYCD       R*8      I  : OPyC density [g/cc]
C  OPYCE       R*8      V  : Elastic modulus of OPyC [MPa]
C  OPYCF       R*8      I  : Initial OPyC characteristic strength [MPa.meter^3/m]
C  OPYCFAIL    I*4      S  : Count of number of particles with OPyC layer failure
C  OPYCFAIL0   I*4      S  : Count of number of particles with OPyC layer failure
C						   during irradiation one group prior to OPYCFAIL. (For statistical purpose)
C  OPYCFAILED  I*2      V  : Indicates if OPyC of current particle failed
C  OPYCLC      R*8      V  : Crystallite size of OPyC (A)
C  OPYCM       R*8      I  : OPyC Weibull modulus
C  OPYCNU      R*8      V  : Poisson's ratio of OPyC
C  OPYCREEP    R*8      V  : Creep coefficient of OPyC [strain per MPa-10**21 nvt]
C  OPYCTHK     R*8      I  : Nominal OPyC thickness [microns]
C  OPYCVAR     R*8      I  : OPyC std deviation [microns]
C  OUTTIME     R*8      I  : Period of time pebbles outside of core (seconds)
C  OSPEC       C*8      I  : Output file specifier for the input set
C  OSWR(0:NDEG)R*8      V  : Array of length NDEG+1 storing coefficients of
C                            the polynomial of OPyC radial swelling rate
C  OSWT(0:NDEG)R*8      V  : Array of length NDEG+1 storing coefficients of
C                            the polynomial of OPyC tangential swelling rate
C  OURAT       R*8      I  : Oxygen to Uranium atom ratio
C  PARASET     R*8      S  : Set of parameters for locating maximum and minimum stresses in parametric study
C  P_BLOCK     R*8      V  : Power of certain block [W]
C  P_CORE      R*8      I  : Thermal power of reactor [MWth]
C  P_CHANNEL   R*8      V  : The power of a channel [W]
C  P_PAR       R*8      V  : The power of a particle [MW]
C  PACKING     R*8      V  : Packing fraction of pebbles in the reactor core
C  PAMB        R*8      I  : Ambient Pressure [MPa]
C  PARAMETRIC_STUDY     L  : Flag for parametric study
C  PARFAIL     I*4      S  : Count of number of particles that fail
C  PARFAIL0    I*4      S  : Count of number of particles that fail
C                            during irradiation one group prior to PARFAIL. (For statistical purpose)
C  PARFAILED   I*2      V  : Indicates if current particle failed
C  PATH        R*8      V  : The stream line of the pebble
C  PEBBLE_R    R*8      V  : Pebble radial position in the reactor core
C  PEBBLE_Z    R*8      V  : Pebble axial position in the reactor core
C  PEBRADIUS   R*8      I  : Pebble radius [m]
C  PERTURBATION_ANALYSIS L : Flag for one parameter perturbation analysis
C						   in parametric study
C  PFZRADIUS   R*8      I  : Pebble fuel zone radius [m]
C  PIE         R*8      P  : Value of pi
C  POWDISTR    R*8      I  : Array of power distribution in reactor core, used
C						   in the branch of PSWITCH = 1.
C  POWER       R*8      V  : Current power of pebbles in the reactor core
C  PRESS       R*8      V  : Pressure of gas in particle [MPa]
C  PROBIF      R*8      S  : Failure probability of IPyC during irradiation
C  PROBOF      R*8      S  : Failure probability of OPyC during irradiation
C  PROBPF      R*8      S  : Failure probability of particles during irradiation
C  PROBSF      R*8      S  : Failure probability of SiC during irradiation
C  PSTATE      C*3      C  : Indicates the current particle state
C  PSWITCH     I*4      C  : Flag for choosing how to get power history
C                            PSWITCH = 1: user provides power distribution
C                            PSWITCH = 2: user provides power history
C                            PSWITCH = 3: user doesn't provide power; assume uniform power density
C  QPPP        R*8      V  : Current power density (W/m^3)
C  QPPP_AVG    R*8      I  : Average power density (W/m^3)
C  R_IN_PEBBLE R*8      V  : The position of the sampled particle in pebble [m]
C  R           R*8      P  : Gas Constant [cc-MPa/gm-mol-K]
C  R1          R*8      V  : Actual Kernel radius [microns]
C  R13         R*8      M  : Cube of R1  [microns^3]
C  R2          R*8      V  : Actual Buffer outer radius [microns]
C  R23         R*8      M  : Cube of R2  [microns^3]
C  R3          R*8      V  : Actual IPyC outer radius [microns]
C  R33         R*8      M  : Cube of R3  [microns^3]
C  R4          R*8      V  : Actual SiC outer radius [microns]
C  R43         R*8      M  : Cube of R4  [microns^3]
C  R5          R*8      V  : Actual OPyC outer radius [microns]
C  R53         R*8      M  : Cube of R5  [microns^3]
C  RADIUS      R*8      M  : Position in the layers of particles
C  RAND        R*8      F  : Random number generator
C  RNCASES     R*8      M  : Number of cases (NCASES) as a REAL*8
C  RRNCASES    R*8      M  : Reciprocal of RNCASES
C  RUNIRR      C*8      I  : Treatments during irradiation
C							'STRESS': Stress calculation without failure evaluation
C							'FAILURE': Stress calculation with failure evaluation
C  SHEARIPYC   R*8      V  : The shear force per unit length on SiC surface
C                            induced by IPyC crack [MPa.um]
C  SHEAROPYC   R*8      V  : The shear force per unit length on SiC surface
C                            induced by OPyC crack [MPa.um]
C  SHUFFLE     I*4      I  : times that pebbles enter the reactor
C  SICALPHA    R*8      V  : Thermal expansion coeff of SiC [per K]
C  SICE        R*8      V  : Elastic modulus of SiC [MPa]
C  SICF        R*8      I  : Initial SiC characteristic strength [MPa.meter^3/m]
C  SICFAIL     I*4      S  : Count of number of particles with SiC failure
C  SICFAIL0    I*4      S  : Count of number of particles with SiC failure
C						   during irradiation one group prior to SICFAIL. (For statistical purpose)
C  SICFAILED   I*2      V  : Indicates if SiC of current particle failed
C  SICKIC0     R*8      I  : SiC Mean Fracture Toughness [MPa*um1/2]
C  SICKVAR     R*8      I  : SiC Fracture Toughness std. deviation [MPa*um1/2]
C  SICM        R*8      I  : SiC Weibull modulus
C  SICNU       R*8      V  : Poisson's ratio of SiC
C  SICTHK      R*8      I  : Nominal SiC thickness [microns]
C  SICVAR      R*8      I  : SiC std deviation [microns]
C  SIG_LOWER   R*8      P  : Stress lower limit for histogram
C  SIG_UPPER   R*8      P  : Stress upper limit for histogram
C  SIGFCIPYC   R*8      S  : Actual fracture strength of IPyC when maximum IPyC stress prevails [MPa]
C  SIGFCOPYC   R*8      S  : Actual fracture strength of OPyC when maximum OPyC stress prevails [MPa]
C  SIGFCSIC    R*8      S  : Actual fracture strength of SiC when maximum SiC stress prevails [MPa]
C  SIGFSIC     R*8      V  : Actual fracture strength of SiC [MPa]
C  SIGFSICI    R*8      V  : Unirradiated fracture strength of SiC [MPa]
C  SIGFIPYC    R*8      V  : Actual fracture strength of IPyC [MPa]
C  SIGFIPYCI   R*8      V  : Unirradiated fracture strength of IPyC [MPa]
C  SIGFOPYC    R*8      V  : Actual fracture strength of OPyC [MPa]
C  SIGFOPYCI   R*8      V  : Unirradiated fracture strength of OPyC [MPa]
C  SIGLBARI    R*8      S  : Population mean of the tangential stress
C                            in IPyC at the end of irradiation  [MPa]
C  SIGLBARO    R*8      S  : Population mean of the tangential stress
C                            in OPyC at the end of irradiation  [MPa]
C  SIGLBARS    R*8      S  : Population mean of the tangential stress
C                            in SiC at the end of irradiation  [MPa]
C  SIGLIPYC    R*8      V  : The tangential stress in IPyC at the end of irradiation [MPa]
C  SIGLOPYC    R*8      V  : The tangential stress in OPyC at the end of irradiation [MPa]
C  SIGLSIC     R*8      V  : The tangential stress in SiC at the end of irradiation [MPa]
C  SIGLVARI    R*8      S  : Standard deviation of the tangential stress
C                            in IPyC at the end of irradiation  [MPa]
C  SIGLVARO    R*8      S  : Standard deviation of the tangential stress
C                            in OPyC at the end of irradiation  [MPa]
C  SIGLVARS    R*8      S  : Standard deviation of the tangential stress
C                            in SiC at the end of irradiation  [MPa]
C  SIGMBARI    R*8      S  : Population mean of minimum tangential stress
C                            in IPyC during irradiation  [MPa]
C  SIGMBARO    R*8      S  : Population mean of minimum tangential stress
C                            in OPyC during irradiation  [MPa]
C  SIGMBARS    R*8      S  : Population mean of minimum tangential stress
C                            in SiC during irradiation  [MPa]
C  SIGMIPYC    R*8      V  : Minimum tangential stress in IPyC [MPa]
C  SIGMOPYC    R*8      V  : Minimum tangential stress in OPyC [MPa]
C  SIGMSIC     R*8      V  : Minimum tangential stress in SiC [MPa]
C  SIGMSICM    R*8      S  : The minimum value of SIGMSIC in parametric study
C  SIGMSICX    R*8      S  : The maximum value of SIGMSIC in parametric study
C  SIGMVARI    R*8      S  : Standard deviation of min tangential stress
C                            in IPyC during irradiation  [MPa]
C  SIGMVARO    R*8      S  : Standard deviation of min tangential stress
C                            in OPyC during irradiation  [MPa]
C  SIGMVARS    R*8      S  : Standard deviation of min tangential stress
C                            in SiC during irradiation  [MPa]
C  SIGR        R*8      V  : Array of radial stress distribution
C  SIGR0       R*8      V  : Array of radial residual stress distribution after each cycle
C  SIGT        R*8      V  : Array of tangential stress distribution
C  SIGT0       R*8      V  : Array of tangential residual stress distribution after each cycle
C  SIGTIPYC    R*8      V  : Positional mean tangential stress across IPyC [MPa]
C  SIGTOPYC    R*8      V  : Positional mean tangential stress across OPyC [MPa]
C  SIGTSIC     R*8      V  : Positional mean tangential stress in across SiC [MPa]
C  SIGXBARI    R*8      S  : Population mean of maximum tangential stress
C                            in IPyC during irradiation  [MPa]
C  SIGXBARO    R*8      S  : Population mean of maximum tangential stress
C                            in OPyC during irradiation  [MPa]
C  SIGXBARS    R*8      S  : Population mean of maximum tangential stress
C                            in SiC during irradiation  [MPa]
C  SIGXIPYC    R*8      V  : Maximum tangential stress in IPyC [MPa]
C  SIGXIPYCM   R*8      S  : The minimum value of SIGXIPYC in parametric study
C  SIGXIPYCX   R*8      S  : The maximum value of SIGXIPYC in parametric study
C  SIGXOPYC    R*8      V  : Maximum tangential stress in OPyC [MPa]
C  SIGXSIC     R*8      V  : Maximum tangential stress in SiC [MPa]
C  SIGXVARI    R*8      S  : Standard deviation of max tangential stress
C                            in IPyC during irradiation  [MPa]
C  SIGXVARO    R*8      S  : Standard deviation of max tangential stress
C                            in OPyC during irradiation  [MPa]
C  SIGXVARS    R*8      S  : Standard deviation of max tangential stress
C                            in SiC during irradiation  [MPa]
C  SRDOT_IPYC	 R*8      V  : Radial swelling rate of IPyC [/10^21nvt]
C  SRDOT_OPYC  R*8      V  : Radial swelling rate of OPyC [/10^21nvt]
C  SSPOLY      R*8      M  : Array for polynomial fitting by IMSL
C  STAT        R*8      M  : Array for polynomial fitting by IMSL
C  STDOT_IPYC  R*8      V  : Tangential swelling rate of IPyC [/10^21nvt]
C  STDOT_OPYC  R*8      V  : Tangential swelling rate of OPyC [/10^21nvt]
C  STATUS      C*11     C  : Program FUEL status identifier
C  SURFACE_ANALYSIS     L  : Flag for two parameter surface analysis
C						   in parametric study
C  T_GASIN     R*8      P  : Coolant (He) entry temperature [C]
C  T_GASOUT    R*8      P  : Coolant (He) exit temperature  [C]
C  T_IRR       R*8      I  : Irradiation temperature [C]
C  TAB         C*1      M  : ASCII tab character
C  T_PARTICLE  R*8      V  : Temperature distribution in the particles (C)
C                            T_PARTICLE(0) is the value at the center of fuel
C                            T_PARTICLE(1-5) correspond to the values at R1-R5, respectively.
C  TEMPERORY   R*8      M  : Variable for temperory usage
C  TGRAD       R*8      M  : Temperature gradient in paticle (K/m)
C  T_HE        R*8      V  : Bulk Helium temperature [C]
C  TIME        R*4      S  : Elapsed wall clock time [sec]
C  TIME0       R*4      S  : Starting wall clock time [sec]
C  TIMELIMIT   R*8      P  : Upper time limit for histogram
C  TIMESTEP    I*4      M  : Time step in the channel
C  TIMESTEP_A  I*4      M  : Accumulative time step used for registering
C                            variables in the history card -- HCARD
C  TITLE       C*40     I  : Input set descriptive title
C  TRIANGLE    R*8      F  : Function which generates triangular distribution
C  U235E       R*8      V  : Actual U235 enrichment [%]
C  U235ENR     R*8      I  : U235 Enrichment [%]
C  U235VAR     R*8      I  : Enrichment std deviation [%]
C  UR          R*8      V  : Array of radial displacement distribution
C  UR0         R*8      V  : Array of radial residual displacement distribution after each cycle
C  URCP        R*8      V  : Array recording the elastic radial displacement
C					       of fully relaxed PyC layers at point of cracking.
C					       The symbol means 'Ur of C prime'.
C  USERSEED    LOG      I  : Flags user input of random number seed
C  VERSION     C*12     C  : Program FUEL version identifier
C  VOIDVOL     R*8      V  : Void volume [cc]
C  WEIBULL     R*8      F  : Function that samples a Weibull variate
C                             from a given mean and modulus
C  WFU         R*8      V  : Weight fraction of uranium in kernel
C  WHICH_BTH   I*4      V  : The batch which the pebble is in
C  WHICH_BLK   I*4      V  : The block/layer which the pebble is in
C  WHICH_CHN   I*4      V  : The channel which the pebble is in
C  X0          R*8      M  : Numerical constant 0
C  X1          R*8      M  : Numerical constant 1
C  X2          R*8      M  : Numerical constant 2
C
C    
C  Functions and subroutines called
C
C  Intrinsic functions:
C     CHAR
C     DBLE
C     DSQRT
C     INT
C
C  External functions:
C     BAFI_PYC
C     E_PYC
C     E_SIC
C     ELAPTIME
C     FLUX
C     GAUSSIAN
C	INITSEED
C     LEAP
C     NDAYS
C	RAND
C     TIMEON
C     TRIANGLE
C     WEIBULL
C
C
C  Subroutines:
C     CORE
C     CURRENTDATE
C     CURRENTTIME
C     EPI_C
C     FAILURE
C     FEEDPEBBLE
C     GASRLS
C     LOCATE
C     LOCATE_ARRAY
C     M_ANALYSIS
C     MACHINE
C     OPENFILE
C     SSICREEP_PYC
C     STRENGTH_PYC
C     STRENGTH_SIC
C     STRESS
C     SWELLR
C     SWELLU
C     TEMPERATURE
C     TPARTICLE
C     TPEBBLE
C	ERR_HANDLER
C
C************************************************************************
C
C************************************************************************
C                                                                       *
C  Main program starts
C
      PROGRAM TIMCOAT
      INCLUDE 'link_fnl_static.h'
      !DEC$ OBJCOMMENT LIB:'libiomp5md.lib'

      USE IMSL_LIBRARIES ! Use IMSL math and statistics libraries
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C  1. Reactor specifications
      DOUBLE PRECISION EOLBUP, EOLFLU, POWER, QPPP_AVG, QPPP, T_HE,
     &                 DT, OUTTIME, IRRTIME,
     &                 PEBBLE_R, PEBBLE_Z, P_CHANNEL, P_BLOCK, P_PAR
	DOUBLE PRECISION CORE_HEIGHT, CORE_RADIUS, P_CORE, T_IRR, T_GASIN,
     &                 T_GASOUT, MF_HE, PACKING
	DOUBLE PRECISION PEBRADIUS, PFZRADIUS, K_PM, K_PFM, K_PFZ,
     &                 R_IN_PEBBLE
C  2. Particle geometry
	DOUBLE PRECISION KERNDIA, KERNVAR, BUFFTHK, BUFFVAR, IPYCTHK,
     &                 IPYCVAR, SICTHK, SICVAR, OPYCTHK, OPYCVAR,
     &                 R10, R20, R30, R40, R50, R1, R2, R3, R4, R5,
     &                 R13, R23, R33, R43, R53
C  3. Layer properties
	DOUBLE PRECISION KDEN, KERND, KERNDVAR, KERNT, BDEN, BUFFD,
     &                 BUFFDVAR, BUFFT, IPYCD, OPYCD, IPYCBAF0, 
     &                 IPYCBAF0I, IPYCBAFVAR, IPYCBAFI, IPYCCRATE, 
     &                 IPYCLC, OPYCBAF0, OPYCBAF0I, OPYCBAFVAR, 
     &                 OPYCBAFI, OPYCCRATE, OPYCLC
C  4. Thermal expansion coefficients, creep coefficients and elastic moduli
	DOUBLE PRECISION IPYCALPHA, SICALPHA, OPYCALPHA, 
     &                 IPYCREEP, IPYCCNU, OPYCREEP, OPYCCNU,
     &                 IPYCE, IPYCNU, SICE, SICNU, OPYCE, OPYCNU
C     Swelling variables of pyrocarbon layers
      DOUBLE PRECISION SRDOT_IPYC, STDOT_IPYC, SRDOT_OPYC, STDOT_OPYC,
     &                 ISWR, ISWT, OSWR, OSWT
C  5. Fracture strengths and toughness of layers
	DOUBLE PRECISION IPYCF, IPYCM, SIGFIPYCI, SIGFIPYC, 
     &                 SICF, SICM, SIGFSICI, SIGFSIC,
     &                 OPYCF, OPYCM, SIGFOPYCI, SIGFOPYC, 
     &				 SIGFCIPYC, SIGFCSIC, SIGFCOPYC,
     &                 SICKIC0, SICKVAR, KICSICI, KICSIC, 
     &				 KIIPYC, KIOPYC, KI1, KI2, CL
C  6. Fuel specifications
      DOUBLE PRECISION CURAT, OURAT, U235ENR, U235VAR, U235E, AWU, WFU,
     &                 PAMB, T_PARTICLE
C  7. Stresses
	DOUBLE PRECISION SIGTIPYC, SIGTSIC, SIGTOPYC, SIGXIPYC, SIGXSIC,
     &                 SIGXOPYC, SIGMIPYC,SIGMSIC,  SIGMOPYC, SIGLIPYC,
     &				 SIGLOPYC, SIGLSIC, SIGR, SIGT, SIGR0, SIGT0
C  8. Strains & Displacement
	DOUBLE PRECISION EPIR, EPIT, EPIR0, EPIT0, UR, UR0
     &				 E_IPYCREEP, E_OPYCREEP, E_IPYCREEP0, E_OPYCREEP0
C  9. Quantities of cracked PyC layers
	DOUBLE PRECISION EPIRCP, EPITCP, URCP, SHEARIPYC, SHEAROPYC
C  10. Statistical variables
      DOUBLE PRECISION PROBIF, PROBSF, PROBOF, PROBPF,
     &                 CONFIDIF, CONFIDSF, CONFIDOF, CONFIDPF,
     &                 SIGXBARI, SIGXBARS, SIGXBARO,
     &                 SIGXVARI, SIGXVARS, SIGXVARO,
     &				 SIGMBARI, SIGMBARS, SIGMBARO,
     &                 SIGMVARI, SIGMVARS, SIGMVARO,
     &				 SIGLBARI, SIGLBARS, SIGLBARO,
     &                 SIGLVARI, SIGLVARS, SIGLVARO
C  11. Intimediate variables, seeds and counters
      DOUBLE PRECISION MOLES, MOLESU, VOIDVOL, PRESS, BURNUP,
     &                 BPRATE, FFLUX, FLUENCE, FLUENCE_R, D_FLUENCE, 
     &                 DF, RADIUS, DUMMY, RNCASES, RRNCASES, MACHTIME,
     &				 OPERTIME, TEMPERORY
C  12. Array for polynomial fitting by IMSL
      DOUBLE PRECISION SSPOLY, STAT(10)
C  13. Fission Product Attack Variables
      DOUBLE PRECISION DPD, DCORR
C  14. Amoeba Effect Variables
      DOUBLE PRECISION MD, AVGT, TGRAD, KMC
C  Other variables
      DOUBLE PRECISION NA        !Avogadro Number
      DOUBLE PRECISION TIME, TIME0, ELAPTIME
	DOUBLE PRECISION TIMELIMIT, SIG_UPPER, SIG_LOWER
	DOUBLE PRECISION PARASET(1:2, 1:4), SIGXIPYCX, SIGXIPYCM,
     &				 SIGMSICX, SIGMSICM
	INTEGER*4 NPEBBLE, NPARTICLE
      INTEGER*4 PARFAIL, IPYCFAIL, SICFAIL, OPYCFAIL, NCHAR, FMODE,
     &		  PARFAIL0, IPYCFAIL0, SICFAIL0, OPYCFAIL0, FAILTYPE,
     &          NBURP, NCASES, NDUMPED, LENGTH_OF_FILE, I, J, K, N
	INTEGER*4 HISTGRMP, HISTGRMI, HISTGRMO, HISTGRMS
	INTEGER   IKISIC, CORR, MSWITCH  
C  Number of radial and axial divisions for power distribution, and number of total layers/blocks
	INTEGER   NCHANNEL, NAXIAL, NLAYER
	INTEGER   ISEED
	INTEGER   NDIVI, NDIVS, NDIVO
      INTEGER   PSWITCH, COREMODEL
	INTEGER   IPYCFAILED, SICFAILED, OPYCFAILED, PARFAILED
      INTEGER   TIMESTEP, TIMESTEP_A, NBLOCK, NSTEPINCYCLE
      INTEGER   SHUFFLE, WHICH_CHN, WHICH_BLK, WHICH_BTH
	INTEGER   INDEX, INDEX2, ERROR_CODE
	CHARACTER*3  FUELTYPE
	CHARACTER*3  PSTATE
	CHARACTER*4  MCODE
	CHARACTER*42 FAILUREPATH
      CHARACTER*12 VERSION
      CHARACTER*11 STATUS
      CHARACTER*8  CTIME,CDATE
      CHARACTER*70 TITLE
	CHARACTER*256 INPFILE1
      CHARACTER*256  INPFILE2, INPFILE3, INPFILE4
      CHARACTER*256  OSPEC
      CHARACTER*8  RUNIRR
      CHARACTER*3  MACH,COMP,FILESTAT
      INTEGER INPSTAT
      CHARACTER*1  CHOICE1
      CHARACTER*1  TAB     
      CHARACTER*1  FF
      LOGICAL NOMINAL       ! Overrides statistical sampling
	LOGICAL DIFFUSION
	LOGICAL HISTOGRAM
	LOGICAL PARAMETRIC_STUDY, PERTURBATION_ANALYSIS, SURFACE_ANALYSIS
      LOGICAL FAIL   ! Flags particle failure
      LOGICAL FOPEN1, FOPEN2, FOPEN3, FOPEN4  ! Flags whether input file is open
      LOGICAL USERSEED      ! Flags user input of random number seed
	LOGICAL DEBUG, EX
	DOUBLE PRECISION RN
C
      PARAMETER (  NA = 6.02D23,
     &             R = 8.3144 D0,
     &             AWC = 12.01115 D0,
     &             AWO = 15.9994 D0,
     &             AWU235 = 235.0439 D0,
     &             AWU238 = 238.0508 D0,
C     &             BUFFT = 2.25 D0,
     &             BU_CONV = 9.02298D5,    !for U235
     &             PIE = 3.1415926535897932385 D0,
     &             CV = 4.1887 90204 78639 09846 D-12,
     &             IKEY = 5,
     &             ITERM = 6,
     &             IOUT = 7,
     &             IDAT1 = 8,
     &             IDAT2 = 9,
     &             IDAT3 = 10,
     &             IDAT4 = 11,
     &             IERR = 12,
     &             IEVOLR = 13,
     &             IEVOLT = 14,
     &             IHIS = 15,
     &             ITEST = 16,
     &             IOUTR = 17,
     &             IOUTT = 18,
     &             IOUTSW = 19,
     &             IOUTSR = 20,
     &             IOUTST = 21,
     &             IOUTER = 22,
     &             IOUTET = 23,
     &             IOUTUR = 24,
     &             IDBG = 25,
     &			 IPRMI = 26,
     &			 IPRMS = 27,
     &			 IOUTFL = 28)
C  Number of divisions for stress distribution         |Modification|
	PARAMETER ( NDIV = 30)    !Must be larger than 9
C  Boundaries for fuel failure histograms
      PARAMETER ( NHIS = 200,
     &            TIMELIMIT = 8.640D7, !1000 days
     &            SIG_UPPER = 1.0D3,
     &            SIG_LOWER = -1.0D3)
C  Degree of polynomial for swelling rate
      PARAMETER ( NDEG = 3)
      INTEGER, PARAMETER :: END_OF_FILE = -1
C
C  Arrays to DIMENSION
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: POWDISTR
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: IRRHISTRY
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: HCARDB
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: HSIGRB, HSIGTB
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: CHANNELS
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: BLOCKS
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: PATH
C  History card for the particle -- a big array
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: HCARD
C  History stresses for the particle
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: HSIGR, HSIGT
	INTEGER, DIMENSION (:,:), ALLOCATABLE :: BLOCKMAP
      DIMENSION ISWR(0:NDEG), ISWT(0:NDEG), OSWR(0:NDEG), OSWT(0:NDEG) 
	DIMENSION T_PARTICLE(0:5)
	DIMENSION SIGR(1:NDIV), SIGT(1:NDIV), SIGR0(1:NDIV), SIGT0(1:NDIV)
	DIMENSION EPIR(1:NDIV), EPIT(1:NDIV), EPIR0(1:NDIV), EPIT0(1:NDIV)
	DIMENSION UR(1:NDIV), UR0(1:NDIV)
	DIMENSION EPIRCP(1:NDIV), EPITCP(1:NDIV), URCP(1:NDIV)
	DIMENSION E_IPYCREEP(1:2), E_OPYCREEP(1:2)
	DIMENSION E_IPYCREEP0(1:2), E_OPYCREEP0(1:2)
	DIMENSION HISTGRMP(0:NHIS,4), HISTGRMI(0:NHIS,4), 
     &          HISTGRMO(0:NHIS,4), HISTGRMS(0:NHIS,4)
      DIMENSION SSPOLY(NDEG+1)
      DIMENSION FAILTYPE(0:2)
C
C  Define interfaces of external procedures
	INTERFACE
        SUBROUTINE CORE(COREMODEL, CHANNELS, BLOCKS, BLOCKMAP)
          DOUBLE PRECISION CHANNELS, BLOCKS
	    INTEGER COREMODEL, BLOCKMAP
	    DIMENSION CHANNELS(1:,-1:), BLOCKS(1:,1:)
          DIMENSION BLOCKMAP(1:,1:)
	  END SUBROUTINE
	END INTERFACE
C
	INTERFACE
	  SUBROUTINE FEEDPEBBLE(COREMODEL, CHANNELS, ENTRANCE, 
     &						WHICH_CHN, PATH)
          DOUBLE PRECISION CHANNELS, ENTRANCE, PATH
          INTEGER   COREMODEL, WHICH_CHN
	    DIMENSION CHANNELS(1:,-1:), PATH(1:,1:)
	  END SUBROUTINE
	END INTERFACE
C
	INTERFACE
	  SUBROUTINE STRENGTH_PYC(LAYER, FLAG, D, BAF0, FLU, T, SIGT,
     &					      SIGMA0, SIGF, M)
		DOUBLE PRECISION D, BAF0, FLU, T, SIGT, SIGMA0, SIGF, M
		CHARACTER*4 LAYER
		CHARACTER*3 FLAG		
		DIMENSION SIGT(1:)
	  END SUBROUTINE
	END INTERFACE
C
	INTERFACE
	  SUBROUTINE STRENGTH_SIC(FLAG, FLU, T, SIGT, SIGMA0, SIGF, M)
		DOUBLE PRECISION FLU, T, SIGT, SIGMA0, SIGF, M
		CHARACTER*3 FLAG		
		DIMENSION SIGT(1:)
	  END SUBROUTINE
	END INTERFACE
C
	INTERFACE
	  SUBROUTINE EPI_C(LAYER, SIGR, SIGT, CREEP, CNU, DFLU, EC)
		DOUBLE PRECISION SIGR, SIGT, CREEP, CNU, DFLU, EC
		CHARACTER*(*) LAYER
		DIMENSION SIGR(1:), SIGT(1:)
		DIMENSION EC(1:2)
	  END SUBROUTINE
	END INTERFACE
C  Common blocks
      COMMON /MTYPE/ MACH, COMP   !Common block used by date and timing routines
      COMMON /PBED/  CORE_HEIGHT, CORE_RADIUS, P_CORE, T_GASIN, 
     &               T_GASOUT, MF_HE, PACKING, NPEBBLE
      COMMON /PEBBLE/ PEBRADIUS, PFZRADIUS, K_PM, K_PFM, K_PFZ,
     &                R_IN_PEBBLE, NPARTICLE
      COMMON /PAR_R0/ R10, R20, R30, R40, R50   !initial particle geometry
      COMMON /PAR_R/ R1, R2, R3, R4, R5         !particle geometry
	COMMON /PAR_DIV/ NDIVI, NDIVS, NDIVO
      COMMON /PAR_K/ KDEN, BDEN, KERNT, BUFFT, U235E, CURAT, OURAT,
     &               FUELTYPE
      COMMON /PAR_M/ IPYCALPHA, IPYCCNU, IPYCE, IPYCNU, IPYCREEP,
     &               OPYCALPHA, OPYCCNU, OPYCE, OPYCNU, OPYCREEP,
     &               SICALPHA, SICE, SICNU, IPYCD, OPYCD,
     &               ISWR, ISWT, OSWR, OSWT
      COMMON /PAR_F/ SIGFIPYC, SIGFOPYC, SIGFSIC, KICSIC
	COMMON /CRACKED_PYC/ EPIRCP, EPITCP, URCP, KIIPYC, KIOPYC,
     &					 KI1, KI2, SHEARIPYC, SHEAROPYC, DF
      COMMON /NFAIL/ IPYCFAIL, SICFAIL, OPYCFAIL, PARFAIL,
     &               IPYCFAILED, SICFAILED, OPYCFAILED, PARFAILED,
     &			   NCHAR, FMODE, PSTATE, MCODE, FAILUREPATH
	COMMON /ERRHANDLE/ ERROR_CODE
	COMMON /OUTPUT/ IKISIC
C      
C  Functions declared EXTERNAL
	EXTERNAL BAFI_PYC
	EXTERNAL E_PYC
      EXTERNAL E_SIC
      !EXTERNAL ELAPTIME
	EXTERNAL FLUX
      EXTERNAL GAUSSIAN
      EXTERNAL INITSEED
C      EXTERNAL RAND
      EXTERNAL TRIANGLE
      EXTERNAL WEIBULL
C
C  User input gathered into NAMELIST for easy handling
C
      NAMELIST /INPUT/  BUFFD,		BUFFDVAR,		BUFFT,
     &                  BUFFTHK,		BUFFVAR,		CORE_HEIGHT,
     &                  CORE_RADIUS,	CURAT,			DEBUG,
     &                  DT,			EOLBUP,			EOLFLU,
     &                  IPYCBAF0I,	IPYCBAFVAR, 	IPYCCRATE,
     &                  IPYCD,		IPYCF,			IPYCLC,
     &                  IPYCM,		IPYCTHK,		IPYCVAR,
     &				  IRRTIME,		ISEED,			KERND,
     &                  KERNDIA,		KERNDVAR,		KERNT,
     &                  KERNVAR,		MF_HE,			NBURP,
     &                  NCASES,		NOMINAL,		NPEBBLE,
     &                  NPARTICLE,	QPPP_AVG,		OPYCBAF0I,
     &                  OPYCBAFVAR,	OPYCCRATE,		OPYCD,
     &                  OPYCF,		OPYCLC,			OPYCM,	
     &                  OPYCTHK,		OPYCVAR,		OSPEC,
     &                  OURAT,		OUTTIME,		P_CORE,
     &                  PAMB,			PEBRADIUS,		PFZRADIUS,
     &                  RUNIRR,		SHUFFLE,		SICF,    
     &                  SICKIC0,		SICKVAR,		SICM,
     &				  SICTHK,		SICVAR,			T_GASIN,
     &				  T_GASOUT,		T_IRR,			TITLE,
     &                  U235ENR,		U235VAR,		USERSEED,
     &                  FUELTYPE,		DIFFUSION,		HISTOGRAM
C  Important reactor core parameters, output to DEBUG file
C
	NAMELIST /REACTOR/ CORE_HEIGHT, CORE_RADIUS,	P_CORE, 
     &				   T_IRR,		T_GASIN,		T_GASOUT,
     &                   MF_HE,		PACKING,		NPEBBLE,
     &				   EOLBUP,		EOLFLU,			QPPP_AVG,
     &				   DT,			OUTTIME,		OPERTIME,
     &				   IRRTIME
C  Important pebble parameters, output to DEBUG file
C
	NAMELIST /PEBBLES/ PEBRADIUS,	PFZRADIUS,		K_PM,
     &				   K_PFM,		K_PFZ,			NPARTICLE
C  Important particle parameters, output to DEBUG file
C
      NAMELIST /PART0/  FUELTYPE,		U235E,			KERNDIA,
     &				  KERNVAR,			
     &				  BUFFTHK,		BUFFVAR,		
     &				  IPYCTHK,		IPYCVAR,		SICTHK,
     &                  SICVAR,		OPYCTHK,		OPYCVAR,
     &                  R10,			R20,			R30,
     &				  R40,			R50, 
     &				  R1,			R2,				R3, 
     &				  R4,			R5,
     &				  KDEN,			BDEN,
     &                  IPYCALPHA,	IPYCBAF0I,		IPYCBAF0,
     &                  IPYCBAFI,		IPYCCNU,		IPYCD, 
     &				  IPYCE,		IPYCNU,			IPYCREEP, 
     &				  IPYCCRATE,	IPYCLC,
     &                  OPYCALPHA,	OPYCBAF0I,		OPYCBAF0,
     &                  OPYCBAFI,		OPYCCNU,		OPYCD, 
     &				  OPYCE,		OPYCNU,			OPYCREEP, 
     &				  OPYCCRATE,	OPYCLC,
     &       			  SICALPHA,		SICE,			SICNU,
     &				  IPYCF,		IPYCM,			SIGFIPYCI, 
     &				  SIGFIPYC,		SICF,			SICM,
     &                  SIGFSICI,		SIGFSIC,		OPYCF,
     &                  OPYCM,		SIGFOPYCI,		SIGFOPYC, 
     &				  SICKIC0,		SICKVAR,		KICSICI, 
     &				  KICSIC,
     &				  NDIVI,		NDIVS,			NDIVO
C
	NAMELIST /PARTI/  R1,			R2,				R3, 
     &				  R4,			R5,
     &                  PRESS,		BURNUP, 		FLUENCE,
     &                  IPYCALPHA,	IPYCBAFI,		IPYCD,
     &                  IPYCE,		IPYCNU,			IPYCCNU,		 
     &				  IPYCREEP, 
     &                  OPYCALPHA,	OPYCBAFI,		OPYCD,
     &                  OPYCE,		OPYCNU,			OPYCCNU,		 
     &				  OPYCREEP, 
     &       			  SICALPHA,		SICE,			SICNU,
     &				  SIGFIPYC,		SIGFOPYC,		SIGFSIC,
     &				  KICSIC,		KIIPYC,			KIOPYC,		CL,
     &                  SIGTIPYC,		SIGTSIC,		SIGTOPYC,
     &                  SIGXIPYC,		SIGXSIC,		SIGXOPYC,
     &                  SIGMIPYC,		SIGMSIC,		SIGMOPYC,
     &                  SIGLIPYC,		SIGLSIC,		SIGLOPYC,
     &				  SIGFCIPYC,	SIGFCSIC,		SIGFCOPYC,
     &                  IPYCFAIL,		OPYCFAIL,		SICFAIL,
     &				  PARFAIL
C
C  Version data
      DATA VERSION /'v 2.0   '/     ! May 2005
      DATA STATUS  /'DEBUG'/
      DATA CTIME   /'        '/, CDATE/'        '/
C  Numerical constants
      DATA X0/0.0 D+00/, X1/1.0 D+00/, X2/2.0 D+00/
	DATA CL/15.0/
C  Miscellaneous variables
      DATA FOPEN1 /.FALSE./, FOPEN2 /.FALSE./
      DATA FOPEN3 /.FALSE./, FOPEN4 /.FALSE./
C
C
C  Start of Executables
C
C  Initialize machine and compiler dependent routines
      CALL MACHINE
C
C  Special characthers
      TAB = CHAR(9)      ! Tab character for histograms
      FF = CHAR(12)
C 
C  ************************* Initialization ****************************
C
C  Display greeting window
C	Call Greeting()
C
	IPYCFAIL = 0
	SICFAIL = 0
	OPYCFAIL = 0
      PARFAIL = 0
	IPYCFAIL0 = 0
	SICFAIL0 = 0
	OPYCFAIL0 = 0
      PARFAIL0 = 0
	NDUMPED = 0
	IKISIC = 29
	CORR = 30
	DO 110  I = 0, 2
	  FAILTYPE(I) = 0
110   CONTINUE
C
      PROBIF = 0.0D0
	PROBSF = 0.0D0
	PROBOF = 0.0D0
	PROBPF = 0.0D0
	CONFIDIF = 0.0D0   ! Initialization of +/- 1 sigma limit on PROBIF
	CONFIDSF = 0.0D0   ! Initialization of +/- 1 sigma limit on PROBSF
      CONFIDOF = 0.0D0   ! Initialization of +/- 1 sigma limit on PROBOF
	CONFIDPF = 0.0D0   
	SIGXBARI = 0.0D0   ! Population mean of maximum tangential stress in IPyC
	SIGXBARS = 0.0D0   ! Population mean of maximum tangential stress in SiC
	SIGXBARO = 0.0D0   ! Population mean of maximum tangential stress in OPyC
      SIGXVARI = 0.0D0   ! Standard deviation of max tangential stress in IPyC
	SIGXVARS = 0.0D0   ! Standard deviation of max tangential stress in SiC
	SIGXVARO = 0.0D0   ! Standard deviation of max tangential stress in OPyC
	SIGMBARI = 0.0D0   ! Population mean of minimum tangential stress in IPyC
	SIGMBARS = 0.0D0   ! Population mean of minimum tangential stress in SiC
	SIGMBARO = 0.0D0   ! Population mean of minimum tangential stress in OPyC
      SIGMVARI = 0.0D0   ! Standard deviation of min tangential stress in IPyC
	SIGMVARS = 0.0D0   ! Standard deviation of min tangential stress in SiC
	SIGMVARO = 0.0D0   ! Standard deviation of min tangential stress in OPyC
	SIGLBARI = 0.0D0   ! Population mean of EOL tangential stress in IPyC
	SIGLBARS = 0.0D0   ! Population mean of EOL tangential stress in SiC
	SIGLBARO = 0.0D0   ! Population mean of EOL tangential stress in OPyC
      SIGLVARI = 0.0D0   ! Standard deviation of EOL tangential stress in IPyC
	SIGLVARS = 0.0D0   ! Standard deviation of EOL tangential stress in SiC
	SIGLVARO = 0.0D0   ! Standard deviation of EOL tangential stress in OPyC
C
C  Reactor core and TRISO fuel input data dialog boxes
C 115	Call InputDialog(INPFILE1)
      INPFILE1='DS_MPBR1.dat'
	FILESTAT = 'OLD'
      OPEN(FILE = INPFILE1,STATUS = FILESTAT,UNIT = IDAT1, 
     &     IOSTAT=INPSTATUS)
	IF(INPSTATUS>0) THEN
	  WRITE(*,*) 'FILENAME: '//TRIM(INPFILE1)//' CANNOT BE FOUND!'
	  STOP
	END IF
	FOPEN1 = .TRUE.
C
C  Read input data file:
C
C  In general, an IOSTAT < 0 indicates a fatal error condition.
C  In particular, IOSTAT = -1 is an end of file condition.
C  An IOSTAT = 0 indicates no errors or warnings.  On some systems
C  such as the Dec 5000 class workstations, IOSTAT returns a
C  positive integer value.  This positive IOSTAT is a message.
C     
      READ(UNIT = IDAT1, NML = INPUT, IOSTAT=IMSG)
	IF (IMSG .NE. 0) THEN
        CALL ERR_HANDLER(TRIM('Main: IOSTAT error from input file'),
     &					34,IMSG,2,IERR)
	END IF
C
      IF (NCASES .EQ. 0) STOP  ! Detect null set input, END program
C
C
C
	IF (.NOT. USERSEED) ISEED = INITSEED(0)
C
C  Initialize random number generator with ISEED
C
      IF (ISEED .EQ. 0) THEN
	  CALL ERR_HANDLER (
     A 	TRIM('MAIN: RNG seed, ISEED = 0, using default seed 305'),
     B     50, 1, 1, IERR)
        ISEED = 305
	END IF
	CALL RANDOM_SEED (PUT = (/ISEED/))
	CALL RANDOM_NUMBER(RN)
	Z = RN    !Initialize random number generator
C
C  Run TIMCOAT as Version 1 or as Version 2
C  Version 2 adds Pd migration, corrosion (thinning) of the SiC layer, and the Amoeba effect.
C      WRITE(ITERM,*) '(v1 = 1, v2 = 2) why is this thing stupid? dumbo',
C     &                    ' me no likey' 
C      READ(IKEY,*) VERSIONSWITCH
C      WRITE(ITERM,*) '(v1 = 1, v2 = 2) why is this thing stupid? dumbo'
C      READ(IKEY,*) VERSIONSWITCH
      WRITE(ITERM,*) 'Run TIMCOAT as Version 1 or 2 (v1 = 1, v2 = 2)'
      READ(IKEY,*) MSWITCH
C
C     MSWITCH = 2
C
C  Select the type of simulation to run (pebble bed reactor core simulation,
C  irradiation experiment simulation, or constant irradiation simulation)
C
C	Call PowerTypeDialog(PSWITCH)
C
      PSWITCH = 1          !Select reactor power history
                               ! 1 - user provides power distribution
                               ! 2 - user provides power history
                               ! 3 - user doesn't provide power; assume uniform power density
      IF(PSWITCH .EQ. 1) THEN
	  WRITE(ITERM,*) ' Is it from new or old core model?',
     &				 ' Please choose a number (new = 1, old = 2):'
	  READ(IKEY,*) COREMODEL
        WRITE(ITERM,*) ' Enter filename of power distribution ',
     &                 '(Omit .dat extension): '
        READ(IKEY,602) INPFILE2
        CALL OPENFILE(INPFILE2,'.dat','OLD',IDAT2)
        FOPEN2 = .TRUE.
        LENGTH_OF_FILE = 0
C     Scan the length of the file
        DO
          READ(UNIT=IDAT2,FMT=*,IOSTAT=IOS) TEMPERORY
          IF(IOS.EQ.END_OF_FILE) THEN
            EXIT
          ELSE
            LENGTH_OF_FILE = LENGTH_OF_FILE + 1
          END IF
        END DO
        REWIND (UNIT = IDAT2)
C     Determine the dimension of array of power distribution and read in data
        ALLOCATE (POWDISTR(LENGTH_OF_FILE,6))
        DO 120 I=1, LENGTH_OF_FILE
          READ(UNIT=IDAT2,FMT=*)  POWDISTR(I,1), POWDISTR(I,2),
     &      POWDISTR(I,3),POWDISTR(I,4),POWDISTR(I,5),POWDISTR(I,6)
120      CONTINUE
        CLOSE (UNIT = IDAT2, STATUS = 'KEEP')
        FOPEN2 = .FALSE.
C     Need geometric information on the reactor core specified by VSOP model
C     These numbers must be consistent with the input files above.
	  WRITE(ITERM,*) ' Please tell me the number of channels,',
     &				 ' axial divisions, and blocks in the core: '
	  READ(IKEY,*) NCHANNEL, NAXIAL, NLAYER
	  WRITE(ITERM,*) ' Thanks.'
C
	  ALLOCATE (CHANNELS(NAXIAL,-1:NCHANNEL))
	  ALLOCATE (BLOCKS(NLAYER,5))
	  ALLOCATE (PATH(NAXIAL,2))
	  ALLOCATE (BLOCKMAP(2,NCHANNEL))
C     Read in geometric information of channels and blocks ASSOCIATED with INPFILE2.dat
        CALL CORE(COREMODEL, CHANNELS, BLOCKS, BLOCKMAP)
C     Assign arrays for particle history card and history stresses
	  ALLOCATE (HCARD(NAXIAL,19))
	  ALLOCATE (HSIGR(0:NAXIAL,1:NDIV))
	  ALLOCATE (HSIGT(0:NAXIAL,1:NDIV))
      ELSE IF(PSWITCH .EQ. 2) THEN
        WRITE(ITERM,*) ' Please enter the filename ',
     &		       ' for the irradiation history ',
     &                 '(include .dat extension): '
C     (Omit .dat)
       READ(IKEY,602) INPFILE3
       FILESTAT = 'OLD'
C       WRITE(*,*) 'HELLO'
C            STOP
       OPEN(FILE = INPFILE3,STATUS = FILESTAT,UNIT = IDAT3, 
     &     IOSTAT=INPSTATUS)
        FOPEN3 = .TRUE.
        WRITE(ITERM,*) ' Thanks.'
        LENGTH_OF_FILE = 0
C     Scan the length of the file
        DO
          READ(UNIT=IDAT3,FMT=*,IOSTAT=IOS) TEMPERORY
          IF(IOS.EQ.END_OF_FILE) THEN
            EXIT
          ELSE
            LENGTH_OF_FILE = LENGTH_OF_FILE + 1
          END IF
        END DO
        REWIND (UNIT = IDAT3)
C     Determine the dimension of array of irradiaton history and read in data
        ALLOCATE (IRRHISTRY(0:LENGTH_OF_FILE-1,5))
	  DO 121 I = 0, LENGTH_OF_FILE - 1
		READ(UNIT=IDAT3,FMT=*) IRRHISTRY(I,1), IRRHISTRY(I,2),
     &	    IRRHISTRY(I,3), IRRHISTRY(I,4), IRRHISTRY(I,5)
121	  CONTINUE
	  CLOSE (UNIT = IDAT3, STATUS = 'KEEP')
	  FOPEN3 = .FALSE.
      ELSE IF(PSWITCH .EQ. 3) THEN
        WRITE(ITERM,*) ' Simulations will be run assuming uniform ',
     &               'power density and fast flux for the fuel.'
	ENDIF
      WRITE(ITERM,*)
C
C  In order to perform parametric study, place the flag PARAMETRIC_STUDY as true
	PARAMETRIC_STUDY = .FALSE.
C  If perturbation analysis is to be performed, set PERTURBATION_ANALYSIS to .TRUE.
C  If surface analysis is to be performed, set SURFACE_ANALYSIS to .TRUE., instead
	PERTURBATION_ANALYSIS = .TRUE.
c	SURFACE_ANALYSIS = .TRUE.
C
C  Open output files
	OPEN (FILE = TRIM(OSPEC)//'.out',STATUS = "REPLACE",UNIT = IOUT)	
C
	IF(.NOT. PARAMETRIC_STUDY) THEN
C  Open an output file for debugging imformation
	  IF(DEBUG) THEN
	    OPEN(FILE = TRIM(OSPEC)//'.dbg',STATUS="REPLACE",UNIT = IDBG)
C  Open intermediate variable output files for chemistry model
          OPEN (FILE ='KI_SIC'//'.out',STATUS ="REPLACE", UNIT = IKISIC)
	      WRITE (IKISIC,*) 'N','OPERTIME','DT','KI1'
	  OPEN (FILE = 'CORR'//'.out',STATUS = "REPLACE", UNIT = CORR)
	      WRITE (CORR,*) 'N','OPERTIME','R3'
	  END IF         
C  Open an output file for particle histogram
        IF(HISTOGRAM)  THEN
	    OPEN (FILE=TRIM(OSPEC)//'.his',STATUS="REPLACE",UNIT=IHIS)
C    Clear arrays for histograms
	    DO 125 I = 0, NHIS
	      DO 126 J = 1, 4
	        HISTGRMP(I,J) = 0
	        HISTGRMI(I,J) = 0
	        HISTGRMO(I,J) = 0
	        HISTGRMS(I,J) = 0
126         CONTINUE
125       CONTINUE
	  END IF
C  Open an output file to record information about failed particles
C  if Monte Carlo sampling is performed.
	  IF(.NOT.NOMINAL) THEN
	    OPEN(FILE='failures'//'.dat',STATUS = "NEW",UNIT=IOUTFL)
	    IF(PSWITCH.EQ.1) THEN
	      WRITE(IOUTFL, 643)
	    ELSE
	      WRITE(IOUTFL, 644)
	    END IF
	  END IF
C  Open an output file for parametric study
	ELSE
	  OPEN(FILE='_IPyCStress'//'.prm',STATUS="NEW",UNIT=IPRMI)
	  IF(SURFACE_ANALYSIS) THEN
	    OPEN(FILE='_SiCStress'//'.prm',STATUS="NEW",UNIT=IPRMS)
	  END IF
	  DO 131 I = 1, 2
		DO 132 J = 1, 4
		  PARASET(I,J)= 0.0D0
132		CONTINUE
131	  CONTINUE
	  SIGXIPYCX = 0.0D0
	  SIGXIPYCM = 3000.0D0
	  SIGMSICX = -5000.0D0
	  SIGMSICM = 0.0D0
	END IF
C
!C  Initialize date and time routines
!      CALL CURRENTDATE (CDATE)
!      CALL CURRENTTIME (CTIME)
!      TIME0 = ELAPTIME(0.0D0)      ! Starting wall clock time for this run
!	MACHTIME = 0.0D0
C
C  Write headers to terminal
      WRITE(ITERM,601)
	WRITE(ITERM,*)
      WRITE(ITERM,*) '   THE MIT FUEL PERFORMANCE MODEL FOR HTGR -- ',
     A               VERSION,'  ',STATUS
      WRITE(ITERM,*)
      WRITE(ITERM,625) 'The calculation starts at: '
	WRITE(ITERM,617) CDATE,CTIME
C
C  Write headers to output files
      WRITE(IOUT,617) CDATE,CTIME
	WRITE(IOUT,601)
      WRITE(IOUT,*)
      WRITE(IOUT,*) '   THE MIT FUEL PERFORMANCE MODEL FOR HTGR -- ',
     A              VERSION
      WRITE(IOUT,*)
      WRITE(UNIT = IOUT, NML = INPUT)
	WRITE(IOUT,*) ' PSWITCH = ', PSWITCH
      WRITE(IOUT,*)
      WRITE(IOUT,*) ' Initial Random Number SEED = ',ISEED
      WRITE(IOUT,*)
	WRITE(IOUT,*) ' The number of divisions for stress distributions:'
     &              , NDIV
	WRITE(IOUT,*)
      WRITE(IOUT,625) 'The calculation starts at: '
	WRITE(IOUT,617) CDATE,CTIME
C
      WRITE(ITERM,*) 'If you run hundreds of thousands of particles, ',
     &               'it will take hours. Please stand by....'
      IF (NCASES .GT. NBURP) THEN
        WRITE(ITERM,604) NCASES
        WRITE(IOUT,604) NCASES
      END IF
C  The loop below is for parametric study, tuning one or several parameters
C  from a low value to a high one. If parametric study is not to be performed,
C  just execute what is inside once.
	IF(.NOT. PARAMETRIC_STUDY) GO TO 511
C
	IF(PERTURBATION_ANALYSIS) THEN
	  WRITE(IPRMI, *) 'OPyC Density perturbation on nominal ',
     &				  'particle'
	  WRITE(IPRMI, *)
	  WRITE(IPRMI, 639)
	ELSE IF(SURFACE_ANALYSIS) THEN
	  WRITE(IPRMI, *) 'Kernel Diameter - IPyC BAF0  Surface'
	  WRITE(IPRMI, 640) ' '
	  WRITE(IPRMS, *) 'Kernel Diameter - IPyC BAF0  Surface'
	  WRITE(IPRMS, 640) ' '
	  DO 399 IPYCBAF0I = 1.0D0, 1.331D0, 0.01D0
	    WRITE(IPRMI, 629) IPYCBAF0I
	    WRITE(IPRMS, 629) IPYCBAF0I
399	  CONTINUE
	  WRITE(IPRMI, *)
	  WRITE(IPRMS, *)
	END IF
C
	DO 400	OPYCD = 1.8D0, 1.991D0, 0.01D0
	WRITE(IPRMI, 629) OPYCD
	IF(SURFACE_ANALYSIS) WRITE(IPRMS, 629) OPYCD	
	DO 401 IPYCBAF0I = 1.08D0, 1.08D0, 0.01D0
C
C  Prepare several irradiation dependent parameters based on inputs for constant power history
511   IF(PSWITCH.EQ.3) THEN
	  BPRATE = EOLBUP*BU_CONV/IRRTIME
	  FFLUX = EOLFLU*1.0D21/(IRRTIME*86400.0D0)
	  LENGTH_OF_FILE = INT(EOLFLU*1.0D21/(FFLUX*DT)+1.5D0)
        ALLOCATE (HCARDB(0:LENGTH_OF_FILE,19))
	  ALLOCATE (HSIGRB(0:LENGTH_OF_FILE,NDIV))
	  ALLOCATE (HSIGTB(0:LENGTH_OF_FILE,NDIV))
	END IF
C
C  Prepare a few materials properties not dependent on particles and fluence
C
      IPYCNU = 0.235D0                 ! Elastic Poisson's ratio of IPyC
	IPYCCNU = 0.5D0
      OPYCNU = 0.235D0                 ! Elastic Poisson's ratio of OPyC
	OPYCCNU = 0.5D0
      OPYCALPHA = (4.167D-6 * OPYCD) - 2.033D-6  ! OPyC expansion
      IPYCALPHA = (4.167D-6 * IPYCD) - 2.033D-6  ! IPyC expansion
      SICNU = 0.13D0
      SICALPHA = 4.9D-6            ! Thermal expansion of SiC [per C] 
C     
C  Determine the number of divisions in IPyC, SiC and OPyC for outputing
C  stress distributions across the structural layers
      NDIVI = (NDIV-6)/3        ! Reserve 6 points for layer surfaces (IPyC/SiC/OPyC)
	NDIVS = (NDIV-6-NDIVI)/2
	NDIVO = NDIV-6-NDIVI-NDIVS
C
C#######################################################################
C  Main Monte Carlo Loop - sample fuel particles
C       
C
C  Open output files which will be appended after each cycle in the loop below
      IF(.NOT. PARAMETRIC_STUDY) THEN
	    IF(PSWITCH.EQ.1) THEN
	      OPEN(FILE='test    '//'.out',STATUS="REPLACE",UNIT=ITEST)
	      OPEN(FILE='out_core'//'.dat',STATUS="REPLACE",UNIT=IOUTR)
	       WRITE(IOUTR,637) 
	      OPEN(FILE='out_temp'//'.dat',STATUS="REPLACE",UNIT=IOUTT)
	      OPEN(FILE='out_swel'//'.dat',STATUS="REPLACE",UNIT=IOUTSW)
	      OPEN(FILE='out_sigr'//'.dat',STATUS="REPLACE",UNIT=IOUTSR)
	      OPEN(FILE='out_sigt'//'.dat',STATUS="REPLACE",UNIT=IOUTST)
	      OPEN(FILE='out_epir'//'.dat',STATUS="REPLACE",UNIT=IOUTER)
	      OPEN(FILE='out_epit'//'.dat',STATUS="REPLACE",UNIT=IOUTET)
	      OPEN(FILE='out_ur'//'.dat',STATUS="REPLACE",UNIT=IOUTUR)
	    ELSE IF(PSWITCH.EQ.2) THEN
	      OPEN(FILE='irr_histry'//'.out',STATUS="REPLACE",UNIT=ITEST)
	      OPEN(FILE='out_swel'//'.dat',STATUS="REPLACE",UNIT=IOUTSW)
	      OPEN(FILE='out_sigr'//'.dat',STATUS="REPLACE",UNIT=IOUTSR)
	      OPEN(FILE='out_sigt'//'.dat',STATUS="REPLACE",UNIT=IOUTST)
	      OPEN(FILE='out_epir'//'.dat',STATUS="REPLACE",UNIT=IOUTER)
	      OPEN(FILE='out_epit'//'.dat',STATUS="REPLACE",UNIT=IOUTET)
	      OPEN(FILE='out_ur'//'.dat',STATUS="REPLACE",UNIT=IOUTUR)
	    ELSE IF(PSWITCH.EQ.3) THEN
	      OPEN(FILE='cap_test'//'.out',STATUS="REPLACE",UNIT=ITEST)
	      OPEN(FILE='out_swel'//'.dat',STATUS="REPLACE",UNIT=IOUTSW)
	      OPEN(FILE='out_sigr'//'.dat',STATUS="REPLACE",UNIT=IOUTSR)
	      OPEN(FILE='out_sigt'//'.dat',STATUS="REPLACE",UNIT=IOUTST)
	      OPEN(FILE='out_epir'//'.dat',STATUS="REPLACE",UNIT=IOUTER)
	      OPEN(FILE='out_epit'//'.dat',STATUS="REPLACE",UNIT=IOUTET)
	      OPEN(FILE='out_ur'//'.dat',STATUS="REPLACE",UNIT=IOUTUR)
	    END IF
      END IF
C
      DO 1000 N = 1, NCASES
C  Sample fuel particles
C    NOMINAL = .TRUE.  --> The same particle runs NCASES power histories
C    NOMINAL = .FALSE. --> Sample NCASES particles (NCASES about 1,000,000)
      IF (NOMINAL) THEN
        R10 = KERNDIA/X2
        R20 = R10 + BUFFTHK
        R30 = R20 + IPYCTHK
        R40 = R30 + SICTHK
        R50 = R40 + OPYCTHK
	  R1 = R10
	  R2 = R20
	  R3 = R30
	  R4 = R40
	  R5 = R50
        KDEN = KERND
        BDEN = BUFFD
        U235E = U235ENR
	  IPYCBAF0 = IPYCBAF0I
	  OPYCBAF0 = OPYCBAF0I
	  CALL STRENGTH_PYC('IPYC','INI',IPYCD,IPYCBAF0,FLUENCE,T_IRR,
     &                    SIGT,IPYCF,SIGFIPYCI,IPYCM)
	  CALL STRENGTH_PYC('OPYC','INI',OPYCD,OPYCBAF0,FLUENCE,T_IRR,
     &                    SIGT,OPYCF,SIGFOPYCI,OPYCM)
        CALL STRENGTH_SIC('INI',FLUENCE,T_IRR,SIGT,SICF,
     &                    SIGFSICI,SICM)
        KICSICI = SICKIC0
C
      IF(.NOT. PARAMETRIC_STUDY) THEN
C         Write headings to stress, strains, and displacement output files
          WRITE(IOUTSR,629) RADIUS
          WRITE(IOUTST,629) RADIUS
          WRITE(IOUTER,642) RADIUS
          WRITE(IOUTET,642) RADIUS
          WRITE(IOUTUR,629) RADIUS
          DO 127 K = 0, NDIVI+1
	      RADIUS = K*(R3-R2)/FLOAT(NDIVI+1) + R2
	      WRITE(IOUTSR,629) RADIUS
            WRITE(IOUTST,629) RADIUS
            WRITE(IOUTER,642) RADIUS
            WRITE(IOUTET,642) RADIUS
            WRITE(IOUTUR,629) RADIUS
127       CONTINUE
          DO 128 K = 0, NDIVS+1
	      RADIUS = K*(R4-R3)/FLOAT(NDIVS+1) + R3
	      WRITE(IOUTSR,629) RADIUS
            WRITE(IOUTST,629) RADIUS
            WRITE(IOUTER,642) RADIUS
            WRITE(IOUTET,642) RADIUS
            WRITE(IOUTUR,629) RADIUS
128       CONTINUE
          DO 129 K = 0, NDIVO+1
	      RADIUS = K*(R5-R4)/FLOAT(NDIVO+1) + R4
	      WRITE(IOUTSR,629) RADIUS
            WRITE(IOUTST,629) RADIUS
            WRITE(IOUTER,642) RADIUS
            WRITE(IOUTET,642) RADIUS
            WRITE(IOUTUR,629) RADIUS
129       CONTINUE
          WRITE(IOUTSR,*)
          WRITE(IOUTST,*)
          WRITE(IOUTER,*)
          WRITE(IOUTET,*)
          WRITE(IOUTUR,*)
	  END IF
C
      ELSE           ! Perform sampling and check limits
C130     R10 = GAUSSIAN(KERNDIA,  KERNVAR)/X2
C        IF (R10 .LE. X0) GOTO 130
C        R1 = R10
C         
C140     R20 = R10 + GAUSSIAN(BUFFTHK,  BUFFVAR)
C        IF (R20 .LE. R10) GOTO 140
C	  R2 = R20
C         
C150     R30 = R20 + GAUSSIAN(IPYCTHK,  IPYCVAR)
C        IF (R30 .LE. R20) GOTO 150
C	  R3 = R30
C         
C160     R40 = R30 + GAUSSIAN(SICTHK,   SICVAR)
C        IF (R40 .LE. R30) GOTO 160
C	  R4 = R40
C         
C170     R50 = R40 + GAUSSIAN(OPYCTHK,  OPYCVAR)
C        IF (R50 .LE. R40) GOTO 170
C	  R5 = R50
C         
C180     KDEN = GAUSSIAN(KERND, KERNDVAR)
C        IF ((KDEN .GT. KERNT).OR.(KDEN .LE. 0.0D0)) GOTO 180
C         
C190     BDEN = GAUSSIAN(BUFFD, BUFFDVAR)
C        IF ((BDEN .GT. BUFFT).OR.(BDEN .LE. 0.0D0)) GOTO 190
C 
C200     U235E = GAUSSIAN(U235ENR,  U235VAR)
C        IF ((U235E .GT. 100.0D0).OR.(U235E .LE. 0.0D0)) GOTO 200
C         
C201	  IPYCBAF0 = GAUSSIAN(IPYCBAF0I, IPYCBAFVAR)
C	  IF (IPYCBAF0 .LE. 1.0D0) GOTO 201
C202	  OPYCBAF0 = GAUSSIAN(OPYCBAF0I, OPYCBAFVAR)
C	  IF (OPYCBAF0 .LE. 1.0D0) GOTO 202
	  R10 = TRIANGLE(KERNDIA,  KERNVAR)/X2
	  R20 = R10 + TRIANGLE(BUFFTHK,  BUFFVAR)
	  R30 = R20 + TRIANGLE(IPYCTHK,  IPYCVAR)
	  R40 = R30 + TRIANGLE(SICTHK,   SICVAR)
	  R50 = R40 + TRIANGLE(OPYCTHK,  OPYCVAR)
	  R1 = R10
	  R2 = R20
	  R3 = R30
	  R4 = R40
	  R5 = R50
	  KDEN = TRIANGLE(KERND, KERNDVAR)
	  BDEN = TRIANGLE(BUFFD, BUFFDVAR)
	  U235E = TRIANGLE(U235ENR,  U235VAR)
	  IPYCBAF0 = TRIANGLE(IPYCBAF0I, IPYCBAFVAR)
	  OPYCBAF0 = TRIANGLE(OPYCBAF0I, OPYCBAFVAR)
C
	  CALL STRENGTH_PYC('IPYC','INI',IPYCD,IPYCBAF0,FLUENCE,T_IRR,
     &                    SIGT,IPYCF,SIGFIPYCI,IPYCM)
	  CALL STRENGTH_PYC('OPYC','INI',OPYCD,OPYCBAF0,FLUENCE,T_IRR,
     &                    SIGT,OPYCF,SIGFOPYCI,OPYCM)
	  TEMPERORY = IPYCF
	  IPYCF = WEIBULL(IPYCF,  IPYCM)
	  SIGFIPYCI = SIGFIPYCI*IPYCF/TEMPERORY
	  TEMPERORY = OPYCF
        OPYCF = WEIBULL(OPYCF,  OPYCM)
	  SIGFOPYCI = SIGFOPYCI*OPYCF/TEMPERORY
C
        CALL STRENGTH_SIC('INI',FLUENCE,T_IRR,SIGT,SICF,SIGFSICI,SICM)
	  TEMPERORY = SICF
	  SICF = WEIBULL(SICF,  SICM)
	  SIGFSICI = SIGFSICI*SICF/TEMPERORY
C
        KICSICI = TRIANGLE(SICKIC0,  SICKVAR)
      END IF
C
      SIGFSIC = SIGFSICI
	SIGFIPYC = SIGFIPYCI
	SIGFOPYC = SIGFOPYCI
	KICSIC = KICSICI
C
	FAIL = .FALSE.        !the new sampled particle is fine :)
	MCODE = 'ISO3'		  !the mechanical analysis for a fine particle
	PSTATE = 'ISO'		  !the particle state that indicates intact particle
	FMODE = -1			  !particle failure mode is reset
	FAILUREPATH = ''      !clear Failure Path
	NCHAR = 0			  !reset character count for Failure Path
	PARFAILED = 0  
	IPYCFAILED = 0
	OPYCFAILED = 0
	SICFAILED = 0
	ERROR_CODE = -1  !capture errors when analyze this particle
	OPERTIME = 0.0D0
C    For constant power irradiation, namely, PSWITCH=3, IRRTIME is known beforehand
      IF((PSWITCH.EQ.1).OR.(PSWITCH.EQ.2)) THEN
	  IRRTIME = 0.0D0
	END IF
      BURNUP = 0.0D0    !the particle/pebble is fresh.
      FLUENCE = 0.0D0
	D_FLUENCE = 0.0D0
	DF = 0.0D0
	QPPP = 0.0D0
	PRESS = 0.0D0
	CALL RANDOM_NUMBER(RN)
       R_IN_PEBBLE = RN
C	CALL RANDOM_NUMBER(RN)  !WHY IS THIS HERE A SECOND TIME?????
	R_IN_PEBBLE = R_IN_PEBBLE*PFZRADIUS
C    Clear kernel migration distance
      MD = 0.0D0
C    Reset statistical variables
      SIGXIPYC = -1000.0D0    !Reset maximum stress in IPyC to lower limit
	SIGXSIC  = -1000.0D0    !Reset maximum stress in SiC to lower limit
	SIGXOPYC = -1000.0D0    !Reset maximum stress in OPyC to lower limit
	SIGMIPYC =  1000.0D0    !Reset minimum stress in IPyC to upper limit
	SIGMSIC  =  1000.0D0    !Reset minimum stress in SiC to upper limit
	SIGMOPYC =  1000.0D0    !Reset minimum stress in OPyC to upper limit
	SIGLIPYC =  0.0D0       !Reset end-of-life stress in IPyC
	SIGLSIC  =  0.0D0       !Reset end-of-life stress in SiC
	SIGLOPYC =  0.0D0       !Reset end-of-life stress in OPyC
	SIGFCIPYC = 0.0D0		!Reset IPyC fracture strength at maximum IPyC stress
	SIGFCSIC =  0.0D0		!Reset SiC fracture strength at maximum SiC stress
	SIGFCOPYC = 0.0D0		!Reset OPyC fracture strength at maximum OPyC stress
C    Clear stresses, strains, and displacement
	DO 205 I = 1, NDIV
	  SIGR0(I) = 0.0D0
	  SIGT0(I) = 0.0D0
	  SIGR(I) = 0.0D0
	  SIGT(I) = 0.0D0
	  EPIR0(I) = 0.0D0
	  EPIT0(I) = 0.0D0
	  EPIR(I) = 0.0D0
	  EPIT(I) = 0.0D0
	  UR0(I) = 0.0D0
	  UR(I) = 0.0D0
	  EPIRCP(I) = 0.0D0
	  EPITCP(I) = 0.0D0
	  URCP(I) = 0.0D0
205   CONTINUE
C    Clear stress intensity factors and shear forces from cracked layers
	KIIPYC = 0.0D0
	KIOPYC = 0.0D0
	SHEARIPYC = 0.0D0
	SHEAROPYC = 0.0D0
C    Clear creep strains
	DO 206 I = 1, 2
	  E_IPYCREEP(I) = 0.0D0
	  E_OPYCREEP(I) = 0.0D0
	  E_IPYCREEP0(I) = 0.0D0
	  E_OPYCREEP0(I) = 0.0D0
206	CONTINUE
C
C  Write to debug file various info. at the initial state
	IF((.NOT.PARAMETRIC_STUDY) .AND. DEBUG .AND. NOMINAL) THEN
	  WRITE(IDBG, *) 'Start of debug file'
	  WRITE(IDBG, *)
	  WRITE(IDBG, *) 'INITIAL PARTICLE PROPERTIES:'
        WRITE(UNIT = IDBG, NML = PART0)
	  WRITE(IDBG, *)
	  WRITE(IDBG, *) '**** Strength data of layers ****'
	  WRITE(IDBG, *)
	  WRITE(IDBG, 635)
	  WRITE(IDBG, *)
	  WRITE(IDBG, 623) FLUENCE, IPYCF, SIGFIPYC, SICF, SIGFSIC,
     &				   OPYCF, SIGFOPYC, KICSIC, KIIPYC, KIOPYC,
     &				   KI1, KI2
	END IF
C
C  Choose the type of power history
C
      IF(PSWITCH.EQ.1) THEN
C  PSWITCH = 1:  run the refueling system
C  Clear history card of particles
        DO 220 I = 1, NAXIAL
	    DO 210 J = 1, 19
	      HCARD(I,J) = 0.0D0
210       CONTINUE
220     CONTINUE
        TIMESTEP_A = 0
C  Inner loop -- shuffle the pebbles into the reactor core and generate the history card
        DO 500 I = 1, SHUFFLE
C    Sample one channel into which the pebble goes and the path it flows
          CALL FEEDPEBBLE(COREMODEL, CHANNELS, ENTRANCE, WHICH_CHN,PATH)
	    NBLOCK = BLOCKMAP(1,WHICH_CHN) !Number of blocks in the selected channel
	    IF ((.NOT.PARAMETRIC_STUDY) .AND. NOMINAL) THEN
 	      WRITE(ITEST, *) OPERTIME, WHICH_CHN, NBLOCK
            END IF
          PEBBLE_Z = PATH(1,1)        !z-position at the entrance
          PEBBLE_R = PATH(1,2)        !r-position at the entrance
	    WHICH_BTH = I
	    TIMESTEP = 0
	    FLUENCE_R = 0.0D0
	    T_HE = T_GASIN
          CALL TEMPERATURE(0.0D0, T_HE, BURNUP, T_PARTICLE)  !Reset T_PARTICLE at the entrance
C    Calculate the total power of this channel for the purpose of scaling He temp.
          P_CHANNEL = 0.0D0
          DO 230 J = 1, NBLOCK
	      WHICH_BLK = BLOCKMAP(2,WHICH_CHN)-BLOCKMAP(1,WHICH_CHN)+J
	      P_BLOCK = 0.0D0
            DO 235 K = 1, SHUFFLE
	        P_BLOCK = P_BLOCK + QPPP_AVG*POWDISTR((WHICH_BLK-1)*
     &				 (SHUFFLE+1)+K, 4)*POWDISTR((WHICH_BLK-1)*
     &				 (SHUFFLE+1)+K, 2)*1.0D-6
235         CONTINUE
            HCARD(J,8) = P_BLOCK
	      P_CHANNEL = P_CHANNEL + P_BLOCK
230       CONTINUE
C
          DO WHILE (TIMESTEP.LT.NBLOCK)  !while in the core
	      TIMESTEP = TIMESTEP + 1
	      TIMESTEP_A = TIMESTEP_A + 1
	      OPERTIME = OPERTIME + DT   !update operation time
		  IRRTIME = IRRTIME + DT
C    Identify which block the pebble is in
	      WHICH_BLK = BLOCKMAP(2,WHICH_CHN)-BLOCKMAP(1,WHICH_CHN)
     &				  + TIMESTEP
C    Locate where on the streamline the pebble is
            CALL LOCATE_ARRAY(PATH(1:NAXIAL,1),NAXIAL,1,
     &                        BLOCKS(WHICH_BLK,1),INDEX)
C    Update new position of pebble by linear interpolation along the streamline
		  PEBBLE_Z=BLOCKS(WHICH_BLK,1)
		  PEBBLE_R=((PEBBLE_Z-PATH(INDEX,1))*PATH(INDEX+1,2)+
     &                (PATH(INDEX+1,1)-PEBBLE_Z)*PATH(INDEX,2))/
     &               (PATH(INDEX+1,1)-PATH(INDEX,1))
            HCARD(TIMESTEP,1) = OPERTIME
	      HCARD(TIMESTEP,2) = PEBBLE_R
	      HCARD(TIMESTEP,3) = PEBBLE_Z
C    Pick out flux from this block/layer
C    Temporarily use FLUX 1 (Fast flux) to obtain flux
            FFLUX = BLOCKS(WHICH_BLK, 2)
            FLUENCE = FLUENCE + FFLUX*DT          !nvt
	      FLUENCE_R = FLUENCE_R + FFLUX*DT      !nvt
	      HCARD(TIMESTEP,4) = FLUENCE/1.0D21  !convert nvt to 10^21nvt
	      HCARD(TIMESTEP,5) = FLUENCE_R/1.0D21
C    Determine current burnup
            IF(TIMESTEP.EQ.1) THEN
	        IF(WHICH_BTH.EQ.1) THEN
	          BURNUP = POWDISTR((WHICH_BLK-1)*(SHUFFLE+1)+
     &					WHICH_BTH,5)/BU_CONV
	        ELSE
	          TEMPERORY = (POWDISTR((WHICH_BLK-1)*(SHUFFLE+1)
     &					   +WHICH_BTH,5)-POWDISTR((WHICH_BLK-1
     &					   +NBLOCK-1)*(SHUFFLE+1)+WHICH_BTH-1,5))
     &					   /BU_CONV
	          IF(TEMPERORY.GT.0.0D0) BURNUP = BURNUP + TEMPERORY
	        END IF
	      ELSE
              BURNUP = BURNUP+(POWDISTR((WHICH_BLK-1)*(SHUFFLE+1)+
     &				 WHICH_BTH,5)-POWDISTR((WHICH_BLK-1-1)*(SHUFFLE+1)
     &				 +WHICH_BTH,5))/BU_CONV
	      END IF
            QPPP = QPPP_AVG*POWDISTR((WHICH_BLK-1)*(SHUFFLE+1)+
     &			 WHICH_BTH, 4)   !pick out power at that position
	      HCARD(TIMESTEP,6) = QPPP
	      HCARD(TIMESTEP,7) = BURNUP
            HCARD(TIMESTEP,8) = T_HE + HCARD(TIMESTEP,8)*
     &                          (T_GASOUT-T_GASIN)/P_CHANNEL
	      T_HE = HCARD(TIMESTEP,8)
C    Calculate temperature distribution in particles
            CALL TEMPERATURE(QPPP, T_HE, BURNUP, T_PARTICLE)  !calculate T distribution
	      HCARD(TIMESTEP,9)  = T_PARTICLE(0)
	      HCARD(TIMESTEP,10)  = T_PARTICLE(1)
	      HCARD(TIMESTEP,11) = T_PARTICLE(2)
	      HCARD(TIMESTEP,12) = T_PARTICLE(3)
	      HCARD(TIMESTEP,13) = T_PARTICLE(4)
	      HCARD(TIMESTEP,14) = T_PARTICLE(5)
		  CALL GASRLS(T_PARTICLE, BURNUP, OPERTIME, DIFFUSION, PRESS)  !calculate internal pressure
	      HCARD(TIMESTEP,15) = PRESS
C	      
C  Switch to run TIMCOAT as mode 1 or mode 2
C  Calculate the kernel migration distance if mode 2 is ON (MSWITCH = 2):				
      IF (MSWITCH .EQ. 2)  THEN
           IF (FUELTYPE .EQ. 'UO2') THEN
            KMC =1.7E-7*exp(-9.21E4/(8.314*(T_PARTICLE(3)+273.15)))
           ELSE
	        KMC = 0.62*exp(-3.11E5/(8.314*(T_PARTICLE(3)+273.15)))
           END IF
           AVGT = 273.15+((T_PARTICLE(0) + T_PARTICLE(5))/2)
	       TGRAD = (T_PARTICLE(0) - T_PARTICLE(5))/(R5*1E-6) 
           MD = MD + (KMC*DT*(1/AVGT**2)*TGRAD)/1E-6
      ELSE
      END IF
C    Account for FP corrosion of SiC, from equation 3.18 in Diecker 2005
      IF (MWSWITCH .EQ. 2)  THEN
          DCORR = (255.2*DT/3600)*exp(-159.9/(0.008314*
     &              (T_PARTICLE(3)+273.15)))	       
	      R3=R3+DCORR
	      R2=R2+DCORR
	      WRITE(CORR,*) N, OPERTIME, R3
      ELSE
      END IF
C
C    Calculate unrestrained swelling rates in PyC
            CALL SWELLU(T_PARTICLE(3), FLUENCE/1.0D21, IPYCD, IPYCBAF0,
     &                  IPYCCRATE, SRDOT_IPYC, STDOT_IPYC)
            CALL SWELLU(T_PARTICLE(5), FLUENCE/1.0D21, OPYCD, OPYCBAF0,
     &                  OPYCCRATE, SRDOT_OPYC, STDOT_OPYC)
	      HCARD(TIMESTEP,16) = SRDOT_IPYC
	      HCARD(TIMESTEP,17) = STDOT_IPYC
	      HCARD(TIMESTEP,18) = SRDOT_OPYC
	      HCARD(TIMESTEP,19) = STDOT_OPYC
	      IF(BURNUP.GT.EOLBUP) EXIT
          END DO
C
C  Do analyses after the pebble goes through the reactor core once
C  and the history card is created.
C
C  First time stress analysis: unrestrained condition
C    Fit polynomials for swelling rate to be used in stress analysis
	    IF(TIMESTEP.EQ.1) THEN
	      ISWR(0) = HCARD(TIMESTEP,16)
	      ISWT(0) = HCARD(TIMESTEP,17)
	      OSWR(0) = HCARD(TIMESTEP,18)
	      OSWT(0) = HCARD(TIMESTEP,19)
	      DO 241 K = 1, NDEG
	        ISWR(K) = 0.0D0
	        ISWT(K) = 0.0D0
	        OSWR(K) = 0.0D0
	        OSWT(K) = 0.0D0
241         CONTINUE
          ELSE IF(TIMESTEP.EQ.2) THEN  !linear interpolation
	      ISWR(0) = (HCARD(TIMESTEP,5)*HCARD(1,16) - 
     &                 HCARD(1,5)*HCARD(TIMESTEP,16))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      ISWT(0) = (HCARD(TIMESTEP,5)*HCARD(1,17) - 
     &                 HCARD(1,5)*HCARD(TIMESTEP,17))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      OSWR(0) = (HCARD(TIMESTEP,5)*HCARD(1,18) - 
     &                 HCARD(1,5)*HCARD(TIMESTEP,18))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      OSWT(0) = (HCARD(TIMESTEP,5)*HCARD(1,19) - 
     &                 HCARD(1,5)*HCARD(TIMESTEP,19))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      ISWR(1) = (HCARD(TIMESTEP,16) - HCARD(1,16))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      ISWT(1) = (HCARD(TIMESTEP,17) - HCARD(1,17))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      OSWR(1) = (HCARD(TIMESTEP,18) - HCARD(1,18))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      OSWT(1) = (HCARD(TIMESTEP,19) - HCARD(1,19))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      IF(NDEG.GE.2) THEN
              DO 242 K = 2, NDEG
	          ISWR(K) = 0.0D0
	          ISWT(K) = 0.0D0
	          OSWR(K) = 0.0D0
	          OSWT(K) = 0.0D0
242           CONTINUE
            END IF
          ELSE IF((TIMESTEP-2).LT.NDEG) THEN
	      CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,16),TIMESTEP-2,
     &                  ISWR(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,17),TIMESTEP-2,
     &                  ISWT(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,18),TIMESTEP-2,
     &                  OSWR(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,19),TIMESTEP-2,
     &                  OSWT(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
            DO 243 K = TIMESTEP-1, NDEG
	        ISWR(K) = 0.0D0
	        ISWT(K) = 0.0D0
	        OSWR(K) = 0.0D0
	        OSWT(K) = 0.0D0
243         CONTINUE
          ELSE
C    IMSL library DRCURV is to do curve fitting
            CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,16),NDEG,ISWR,SSPOLY,STAT)
C	      IF (STAT(5).LE.80.0D0) THEN  !R-squared less than 80%
C	        CALL ERR_HANDLER('Main: Curve fitting for ISWR not good enough',
C     &                    44,0,0,IERR)
C	      ELSE IF(STAT(10).NE.0) THEN  !Data contain NaN
C	        CALL ERR_HANDLER('Main: Fitting data for ISWR contain NaN',
C     &                    39,1,1,IERR)
C	      END IF
C
            CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,17),NDEG,ISWT,SSPOLY,STAT)
C	      IF (STAT(5).LE.80.0D0) THEN  !R-squared less than 80%
C	        CALL ERR_HANDLER('Main: Curve fitting for ISWT not good enough',
C     &                    44,0,0,IERR)
C	      ELSE IF(STAT(10).NE.0) THEN  !Data contain NaN
C	        CALL ERR_HANDLER('Main: Fitting data for ISWT contain NaN',
C     &                    39,1,1,IERR)
C	      END IF
C
            CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,18),NDEG,OSWR,SSPOLY,STAT)
C	      IF (STAT(5).LE.80.0D0) THEN  !R-squared less than 80%
C	        CALL ERR_HANDLER('Main: Curve fitting for OSWR not good enough',
C     &                    44,0,0,IERR)
C	      ELSE IF(STAT(10).NE.0) THEN  !Data contain NaN
C	        CALL ERR_HANDLER('Main: Fitting data for OSWR contain NaN',
C     &                    39,1,1,IERR)
C	      END IF
C
            CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,19),NDEG,OSWT,SSPOLY,STAT)
C	      IF (STAT(5).LE.80.0D0) THEN  !R-squared less than 80%
C	        CALL ERR_HANDLER('Main: Curve fitting for OSWT not good enough',
C     &                    44,0,0,IERR)
C	      ELSE IF(STAT(10).NE.0) THEN  !Data contain NaN
C	        CALL ERR_HANDLER('Main: Fitting data for OSWT contain NaN',
C     &                    39,1,1,IERR)
C	      END IF
C    End of polynomial fitting
          END IF
C    Examine the fitting results
C          DO 248 J = 1, TIMESTEP
C	      WRITE(ITEST,623) HCARD(J,4), HCARD(J,5), HCARD(J,16),
C     &                       HCARD(J,17),HCARD(J,18),HCARD(J,19)
C248       CONTINUE
C          WRITE(ITEST,*) 'Fitting coefficients for ISWR, ISWT, OSWR',
C     &                   ' and OSWT --'
C	    WRITE(ITEST,624) ISWR(0), ISWT(0), OSWR(0), OSWT(0)
C	    WRITE(ITEST,624) ISWR(1), ISWT(1), OSWR(1), OSWT(1)
C	    WRITE(ITEST,624) ISWR(2), ISWT(2), OSWR(2), OSWT(2)
C	    WRITE(ITEST,624) ISWR(3), ISWT(3), OSWR(3), OSWT(3)
C	    WRITE(ITEST,*)
C
C
C    Restore apparent creep strains from last cycle
		E_IPYCREEP(1) = E_IPYCREEP0(1)
		E_IPYCREEP(2) = E_IPYCREEP0(2)
		E_OPYCREEP(1) = E_OPYCREEP0(1)
		E_OPYCREEP(2) = E_OPYCREEP0(2)
C
          DO 240 J = 1, TIMESTEP
            IPYCE = E_PYC(IPYCD,IPYCBAF0,IPYCLC,HCARD(J,4),HCARD(J,12))
            OPYCE = E_PYC(OPYCD,OPYCBAF0,OPYCLC,HCARD(J,4),HCARD(J,13))
C
	      CALL SSICREEP_PYC(HCARD(J,4),IPYCD,HCARD(J,12),E_IPYCREEP,
     &						IPYCREEP,IPYCCNU)
            CALL SSICREEP_PYC(HCARD(J,4),OPYCD,HCARD(J,13),E_OPYCREEP,
     &						OPYCREEP,OPYCCNU)
C
            SICE = E_SIC(HCARD(J,13))
C
		  IF(J.GT.1) THEN
			DF = HCARD(J,5) - HCARD(J-1,5)
		  ELSE
			DF = HCARD(J,5)
		  END IF
C    Restore residual stresses, strains, and displacement from last cycle
	      DO 251 K = 1, NDIV
              SIGR(K) = SIGR0(K)
	        SIGT(K) = SIGT0(K)
			EPIR(K) = EPIR0(K)
			EPIT(K) = EPIT0(K)
			UR(K)   = UR0(K)
251         CONTINUE
C    Do mechanical analyses
	      CALL M_ANALYSIS(MCODE, HCARD(J,15), PAMB, HCARD(J,5), SIGR,
     &					  SIGT, EPIR, EPIT, UR)
		  HSIGR(J,:) = SIGR
		  HSIGT(J,:) = SIGT
C    Accumulate apparent creep strains
		  D_FLUENCE = DF
		  CALL EPI_C('IPYC',SIGR,SIGT,IPYCREEP,
     &                 IPYCCNU,D_FLUENCE,E_IPYCREEP)
		  CALL EPI_C('OPYC',SIGR,SIGT,OPYCREEP,
     &                 OPYCCNU,D_FLUENCE,E_OPYCREEP)
240       CONTINUE
C  End of first time stress analysis
C
C  Second time stress analysis: restrained condition
C    Restore apparent creep strains from last cycle
		E_IPYCREEP(1) = E_IPYCREEP0(1)
		E_IPYCREEP(2) = E_IPYCREEP0(2)
		E_OPYCREEP(1) = E_OPYCREEP0(1)
		E_OPYCREEP(2) = E_OPYCREEP0(2)
C
C    convert unrestrained dimensional changes to restrained ones
		DO 255 J = 1, TIMESTEP
C    Calculate irradiated BAF, apparent creep strains and restrained swelling rates in PyC
		  IPYCBAFI = BAFI_PYC(IPYCBAF0, HCARD(J, 4))
		  OPYCBAFI = BAFI_PYC(OPYCBAF0, HCARD(J, 4))
C
		  CALL SSICREEP_PYC(HCARD(J,4),IPYCD,HCARD(J,12),E_IPYCREEP,
     &						IPYCREEP,IPYCCNU)
            CALL SSICREEP_PYC(HCARD(J,4),OPYCD,HCARD(J,13),E_OPYCREEP,
     &						OPYCREEP,OPYCCNU)
C
		  IF(J.GT.1) THEN
			D_FLUENCE = HCARD(J,5) - HCARD(J-1,5)
		  ELSE
			D_FLUENCE = HCARD(J,5)
		  END IF
C
		  CALL EPI_C('IPYC',HSIGR(J,:),HSIGT(J,:),IPYCREEP,
     &                 IPYCCNU,D_FLUENCE,E_IPYCREEP)
		  CALL EPI_C('OPYC',HSIGR(J,:),HSIGT(J,:),OPYCREEP,
     &                 OPYCCNU,D_FLUENCE,E_OPYCREEP)
C
		  CALL SWELLR(IPYCBAFI,E_IPYCREEP,D_FLUENCE,
     &				  HCARD(J,16),HCARD(J,17))
		  CALL SWELLR(OPYCBAFI,E_OPYCREEP,D_FLUENCE,
     &				  HCARD(J,18),HCARD(J,19))
255		CONTINUE
C    Fit polynomials for swelling rate to be used in stress analysis
	    IF(TIMESTEP.EQ.1) THEN
	      ISWR(0) = HCARD(TIMESTEP,16)
	      ISWT(0) = HCARD(TIMESTEP,17)
	      OSWR(0) = HCARD(TIMESTEP,18)
	      OSWT(0) = HCARD(TIMESTEP,19)
	      DO 245 K = 1, NDEG
	        ISWR(K) = 0.0D0
	        ISWT(K) = 0.0D0
	        OSWR(K) = 0.0D0
	        OSWT(K) = 0.0D0
245         CONTINUE
          ELSE IF(TIMESTEP.EQ.2) THEN  !linear interpolation
	      ISWR(0) = (HCARD(TIMESTEP,5)*HCARD(1,16) - 
     &                 HCARD(1,5)*HCARD(TIMESTEP,16))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      ISWT(0) = (HCARD(TIMESTEP,5)*HCARD(1,17) - 
     &                 HCARD(1,5)*HCARD(TIMESTEP,17))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      OSWR(0) = (HCARD(TIMESTEP,5)*HCARD(1,18) - 
     &                 HCARD(1,5)*HCARD(TIMESTEP,18))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      OSWT(0) = (HCARD(TIMESTEP,5)*HCARD(1,19) - 
     &                 HCARD(1,5)*HCARD(TIMESTEP,19))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      ISWR(1) = (HCARD(TIMESTEP,16) - HCARD(1,16))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      ISWT(1) = (HCARD(TIMESTEP,17) - HCARD(1,17))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      OSWR(1) = (HCARD(TIMESTEP,18) - HCARD(1,18))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      OSWT(1) = (HCARD(TIMESTEP,19) - HCARD(1,19))/
     &                (HCARD(TIMESTEP,5)-HCARD(1,5))
	      IF(NDEG.GE.2) THEN
              DO 246 K = 2, NDEG
	          ISWR(K) = 0.0D0
	          ISWT(K) = 0.0D0
	          OSWR(K) = 0.0D0
	          OSWT(K) = 0.0D0
246           CONTINUE
            END IF
          ELSE IF((TIMESTEP-2).LT.NDEG) THEN
	      CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,16),TIMESTEP-2,
     &                  ISWR(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,17),TIMESTEP-2,
     &                  ISWT(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,18),TIMESTEP-2,
     &                  OSWR(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,19),TIMESTEP-2,
     &                  OSWT(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
            DO 247 K = TIMESTEP-1, NDEG
	        ISWR(K) = 0.0D0
	        ISWT(K) = 0.0D0
	        OSWR(K) = 0.0D0
	        OSWT(K) = 0.0D0
247         CONTINUE
          ELSE
C    IMSL library DRCURV is to do curve fitting
            CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,16),NDEG,ISWR,SSPOLY,STAT)
C	      IF (STAT(5).LE.80.0D0) THEN  !R-squared less than 80%
C	        CALL ERR_HANDLER('Main: Curve fitting for ISWR not good enough',
C     &                    44,0,0,IERR)
C	      ELSE IF(STAT(10).NE.0) THEN  !Data contain NaN
C	        CALL ERR_HANDLER('Main: Fitting data for ISWR contain NaN',
C     &                    39,1,1,IERR)
C	      END IF
C
            CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,17),NDEG,ISWT,SSPOLY,STAT)
C	      IF (STAT(5).LE.80.0D0) THEN  !R-squared less than 80%
C	        CALL ERR_HANDLER('Main: Curve fitting for ISWT not good enough',
C     &                    44,0,0,IERR)
C	      ELSE IF(STAT(10).NE.0) THEN  !Data contain NaN
C	        CALL ERR_HANDLER('Main: Fitting data for ISWT contain NaN',
C     &                    39,1,1,IERR)
C	      END IF
C
            CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,18),NDEG,OSWR,SSPOLY,STAT)
C	      IF (STAT(5).LE.80.0D0) THEN  !R-squared less than 80%
C	        CALL ERR_HANDLER('Main: Curve fitting for OSWR not good enough',
C     &                    44,0,0,IERR)
C	      ELSE IF(STAT(10).NE.0) THEN  !Data contain NaN
C	        CALL ERR_HANDLER('Main: Fitting data for OSWR contain NaN',
C     &                    39,1,1,IERR)
C	      END IF
C
            CALL DRCURV(TIMESTEP,HCARD(1:TIMESTEP,5),
     &                  HCARD(1:TIMESTEP,19),NDEG,OSWT,SSPOLY,STAT)
C	      IF (STAT(5).LE.80.0D0) THEN  !R-squared less than 80%
C	        CALL ERR_HANDLER('Main: Curve fitting for OSWT not good enough',
C     &                    44,0,0,IERR)
C	      ELSE IF(STAT(10).NE.0) THEN  !Data contain NaN
C	        CALL ERR_HANDLER('Main: Fitting data for OSWT contain NaN',
C     &                    39,1,1,IERR)
C	      END IF
C    End of polynomial fitting
          END IF
C    Examine the fitting results
C          DO 248 J = 1, TIMESTEP
C	      WRITE(ITEST,623) HCARD(J,4), HCARD(J,5), HCARD(J,16),
C     &                       HCARD(J,17),HCARD(J,18),HCARD(J,19)
C248       CONTINUE
C          WRITE(ITEST,*) 'Fitting coefficients for ISWR, ISWT, OSWR',
C     &                   ' and OSWT --'
C	    WRITE(ITEST,624) ISWR(0), ISWT(0), OSWR(0), OSWT(0)
C	    WRITE(ITEST,624) ISWR(1), ISWT(1), OSWR(1), OSWT(1)
C	    WRITE(ITEST,624) ISWR(2), ISWT(2), OSWR(2), OSWT(2)
C	    WRITE(ITEST,624) ISWR(3), ISWT(3), OSWR(3), OSWT(3)
C	    WRITE(ITEST,*)
C
C
C    Restore apparent creep strains from last cycle
		E_IPYCREEP(1) = E_IPYCREEP0(1)
		E_IPYCREEP(2) = E_IPYCREEP0(2)
		E_OPYCREEP(1) = E_OPYCREEP0(1)
		E_OPYCREEP(2) = E_OPYCREEP0(2)
C
          DO 250 J = 1, TIMESTEP
            IPYCE = E_PYC(IPYCD,IPYCBAF0,IPYCLC,HCARD(J,4),HCARD(J,12))
            OPYCE = E_PYC(OPYCD,OPYCBAF0,OPYCLC,HCARD(J,4),HCARD(J,13))
            SICE = E_SIC(HCARD(J,13))
C
	      CALL SSICREEP_PYC(HCARD(J,4),IPYCD,HCARD(J,12),E_IPYCREEP,
     &						IPYCREEP,IPYCCNU)
            CALL SSICREEP_PYC(HCARD(J,4),OPYCD,HCARD(J,13),E_OPYCREEP,
     &						OPYCREEP,OPYCCNU)
C
		  IF(J.GT.1) THEN
			DF = HCARD(J,5) - HCARD(J-1,5)
		  ELSE
			DF = HCARD(J,5)
		  END IF
C    Restore residual stresses, strains, and displacement from last cycle
	      DO 260 K = 1, NDIV
              SIGR(K) = SIGR0(K)
	        SIGT(K) = SIGT0(K)
			EPIR(K) = EPIR0(K)
			EPIT(K) = EPIT0(K)
			UR(K)   = UR0(K)
260         CONTINUE
C    Do mechanical analyses
	      CALL M_ANALYSIS(MCODE, HCARD(J,15), PAMB, HCARD(J,5), SIGR,
     &					  SIGT, EPIR, EPIT, UR)
C
C    Accumulate apparent creep strains
		  D_FLUENCE = DF
		  CALL EPI_C('IPYC',SIGR,SIGT,IPYCREEP,
     &                 IPYCCNU,D_FLUENCE,E_IPYCREEP)
		  CALL EPI_C('OPYC',SIGR,SIGT,OPYCREEP,
     &                 OPYCCNU,D_FLUENCE,E_OPYCREEP)
C
C    Update the mean fracture strength of PyC layers due to temp. variation and irradiation
		  CALL STRENGTH_PYC('IPYC','IRR',IPYCD,IPYCBAF0,HCARD(J,4),
     &				HCARD(J,12),SIGT,IPYCF,SIGFIPYC,IPYCM)
		  CALL STRENGTH_PYC('OPYC','IRR',OPYCD,OPYCBAF0,HCARD(J,4),
     &				HCARD(J,13),SIGT,OPYCF,SIGFOPYC,OPYCM)
C    Calculate the mean fracture strength of SiC layer
		  CALL STRENGTH_SIC('IRR',HCARD(J,4),HCARD(J,12),SIGT,
     &						SICF,SIGFSIC,SICM)
C    Register maximum and minimum tangential stresses in layers
            DO 261 K = 1, NDIVI+2
              IF (SIGT(K).GT.SIGXIPYC) THEN
			  SIGXIPYC = SIGT(K)
			  SIGFCIPYC = SIGFIPYC
			END IF
			IF (SIGT(K).LT.SIGMIPYC) SIGMIPYC = SIGT(K)
261         CONTINUE
            DO 262 K = 1, NDIVS+2
	        IF (SIGT(NDIVI+2+K).GT.SIGXSIC) THEN
			  SIGXSIC = SIGT(NDIVI+2+K)
			  SIGFCSIC = SIGFSIC
			END IF
	        IF (SIGT(NDIVI+2+K).LT.SIGMSIC)  SIGMSIC = SIGT(NDIVI+2+K)
262         CONTINUE
            DO 263 K = 1, NDIVO+2
	        IF (SIGT(NDIVI+NDIVS+4+K).GT.SIGXOPYC) THEN
                SIGXOPYC = SIGT(NDIVI+NDIVS+4+K)
			  SIGFCOPYC = SIGFOPYC
			END IF
	        IF (SIGT(NDIVI+NDIVS+4+K).LT.SIGMOPYC) 
     &           SIGMOPYC = SIGT(NDIVI+NDIVS+4+K)
263         CONTINUE
C    Register mean tangential stresses in layers
	      SIGTIPYC = 0.0D0
	      DO 264 K = 1, NDIVI+2
	        SIGTIPYC = SIGTIPYC + SIGT(K)
264         CONTINUE
            SIGTIPYC = SIGTIPYC/FLOAT(NDIVI+2)
	      SIGTSIC = 0.0D0
	      DO 265 K = 1, NDIVS+2
	        SIGTSIC = SIGTSIC + SIGT(NDIVI+2+K)
265         CONTINUE
            SIGTSIC = SIGTSIC/FLOAT(NDIVS+2)
            SIGTOPYC = 0.0D0
	      DO 266 K = 1, NDIVO+2
	        SIGTOPYC = SIGTOPYC + SIGT(NDIVI+NDIVS+4+K)
266         CONTINUE
            SIGTOPYC = SIGTOPYC/FLOAT(NDIVO+2)
C    Register end-of-life mean tangential stresses in layers
		  SIGLIPYC = SIGTIPYC
		  SIGLOPYC = SIGTOPYC
		  SIGLSIC  = SIGTSIC
C
C    Evaluate fuel failure
C    DPD is given by eqn 3.14 from Diecker 2005
      IF(RUNIRR .EQ. 'FAILURE') THEN
       DPD = 255.2*(OPERTIME/3600)*exp(-159.9/(0.008314*
     &        (T_PARTICLE(3)+273.15)))        
       CALL FAILURE(SIGR, SIGT, FAIL, FAILTYPE, DPD, N, 
     &                     OPERTIME, DT, MD)
      END IF
C
C    Write to debug file of material strength data
		  IF((.NOT.PARAMETRIC_STUDY).AND. DEBUG .AND. NOMINAL) THEN
			WRITE(IDBG, 623) HCARD(J,4), IPYCF, SIGFIPYC, SICF, 
     &			SIGFSIC, OPYCF, SIGFOPYC, KICSIC, KIIPYC, KIOPYC,
     &			KI1, KI2
		  END IF
C
C    Create histograms here
            IF(HISTOGRAM) THEN
C             Particle histogram
	        IF(PARFAILED.EQ.1) THEN
	          INDEX = INT(NHIS*HCARD(J,1)/TIMELIMIT+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMP(INDEX,1) = HISTGRMP(INDEX,1) + 1
C
	          INDEX = INT(NHIS*(SIGTSIC-SIG_LOWER)/
     &                     (SIG_UPPER-SIG_LOWER)+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMP(INDEX,2) = HISTGRMP(INDEX,2) + 1
C
	          INDEX = INT(NHIS*HCARD(J,4)/EOLFLU+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMP(INDEX,3) = HISTGRMP(INDEX,3) + 1
C
	          INDEX = INT(NHIS*HCARD(J,7)/EOLBUP+0.5D0)
                IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMP(INDEX,4) = HISTGRMP(INDEX,4) + 1
	        END IF
C             IPyC layer histogram
	        IF(IPYCFAILED.EQ.1) THEN
	          INDEX = INT(NHIS*HCARD(J,1)/TIMELIMIT+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMI(INDEX,1) = HISTGRMI(INDEX,1) + 1
C
	          INDEX = INT(NHIS*(SIGTIPYC-SIG_LOWER)/
     &                      (SIG_UPPER-SIG_LOWER)+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMI(INDEX,2) = HISTGRMI(INDEX,2) + 1
C
	          INDEX = INT(NHIS*HCARD(J,4)/EOLFLU+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMI(INDEX,3) = HISTGRMI(INDEX,3) + 1
C
	          INDEX = INT(NHIS*HCARD(J,7)/EOLBUP+0.5D0)
                IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMI(INDEX,4) = HISTGRMI(INDEX,4) + 1
	        END IF
C             OPyC layer histogram
	        IF(OPYCFAILED.EQ.1) THEN
	          INDEX = INT(NHIS*HCARD(J,1)/TIMELIMIT+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMO(INDEX,1) = HISTGRMO(INDEX,1) + 1
C
	          INDEX = INT(NHIS*(SIGTOPYC-SIG_LOWER)/
     &                      (SIG_UPPER-SIG_LOWER)+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMO(INDEX,2) = HISTGRMO(INDEX,2) + 1
C
	          INDEX = INT(NHIS*HCARD(J,4)/EOLFLU+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMO(INDEX,3) = HISTGRMO(INDEX,3) + 1
C
	          INDEX = INT(NHIS*HCARD(J,7)/EOLBUP+0.5D0)
                IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMO(INDEX,4) = HISTGRMO(INDEX,4) + 1
	        END IF
C             SiC layer histogram
	        IF(SICFAILED.EQ.1) THEN
	          INDEX = INT(NHIS*HCARD(J,1)/TIMELIMIT+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMS(INDEX,1) = HISTGRMS(INDEX,1) + 1
C
	          INDEX = INT(NHIS*(SIGTSIC-SIG_LOWER)/
     &                      (SIG_UPPER-SIG_LOWER)+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMS(INDEX,2) = HISTGRMS(INDEX,2) + 1
C
	          INDEX = INT(NHIS*HCARD(J,4)/EOLFLU+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMS(INDEX,3) = HISTGRMS(INDEX,3) + 1
C
	          INDEX = INT(NHIS*HCARD(J,7)/EOLBUP+0.5D0)
                IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMS(INDEX,4) = HISTGRMS(INDEX,4) + 1
	        END IF
	      END IF
C    End of creating histograms
C
C    Output a bunch of results to various output files
            IF ((.NOT.PARAMETRIC_STUDY) .AND. NOMINAL) THEN
              WRITE(IOUTR,626) HCARD(J,1)/86400.0D0, HCARD(J,4), 
     &                         HCARD(J,6), HCARD(J,7), HCARD(J,8)
	        WRITE(IOUTT,627) HCARD(J,1)/86400.0D0, HCARD(J,9),
     &                         HCARD(J,10), HCARD(J,11), HCARD(J,12),
     &                         HCARD(J,13), HCARD(J,14)
	        WRITE(IOUTSW,628) HCARD(J,4), HCARD(J,16), HCARD(J,17), 
     &                         HCARD(J,18), HCARD(J,19)
	        WRITE(IOUTSR,629) HCARD(J,4)
	 	    WRITE(IOUTST,629) HCARD(J,4)
	        WRITE(IOUTER,642) HCARD(J,4)
	 	    WRITE(IOUTET,642) HCARD(J,4)
	        WRITE(IOUTUR,629) HCARD(J,4)
		    DO 280 K = 1, NDIV
	          WRITE(IOUTSR, 629) SIGR(K)
	          WRITE(IOUTST, 629) SIGT(K)
	          WRITE(IOUTER, 641) EPIR(K)
	          WRITE(IOUTET, 641) EPIT(K)
	          WRITE(IOUTUR, 629) UR(K)
280           CONTINUE
              WRITE(IOUTSR,*)
              WRITE(IOUTST,*)
              WRITE(IOUTER,*)
              WRITE(IOUTET,*)
              WRITE(IOUTUR,*)
		  END IF
C
            IF(FAIL) THEN
C             write information about the failed particle
			IF(.NOT.NOMINAL) THEN
			  WRITE(IOUTFL, 645) PARFAIL, FMODE, R1, R2, R3, R4, R5,
     &						   BDEN, IPYCD, OPYCD, IPYCBAF0, OPYCBAF0,
     &						   SIGFIPYC, SIGFOPYC, SIGFSIC, KICSIC,
     &						   SIGLIPYC, SIGLOPYC, SIGLSIC,
     &						   HCARD(J,1)/86400.0D0, HCARD(J,4),
     &						   HCARD(J,7), HCARD(J,11), HCARD(J,15), 
     &						   WHICH_CHN, J, WHICH_BTH, FAILUREPATH
			END IF
		    GO TO 990       
	      END IF
C
250       CONTINUE
          OPERTIME = OPERTIME + OUTTIME  !account for the time outside
	    TIMESTEP_A = TIMESTEP_A + 1
          IF(BURNUP.GT.EOLBUP) THEN
		  NDUMPED = NDUMPED + 1
		  EXIT
	    END IF
C    Save residual stresses, strains, and displacement
	    DO 270 J = 1, NDIV
	      SIGR0(J) = SIGR(J)
	      SIGT0(J) = SIGT(J)
		  EPIR0(J) = EPIR(J)
		  EPIT0(J) = EPIT(J)
		  UR0(J)   = UR(J)
270       CONTINUE
C    Save apparent creep strains gained in this cycle
		E_IPYCREEP0(1) = E_IPYCREEP(1)
		E_IPYCREEP0(2) = E_IPYCREEP(2)
		E_OPYCREEP0(1) = E_OPYCREEP(1)
		E_OPYCREEP0(2) = E_OPYCREEP(2)		
C         maybe do thermal stress analysis and gas release calculation here.
500     CONTINUE
C
C***********************************************************************
C***********************************************************************
C                                                                      *
C                                                                      *
C
      ELSE IF(PSWITCH.EQ.2) THEN
C  PSWITCH = 2: use the irradiation history that user provides.
	  TIMESTEP_A = 0
C     Outer loop: represents the entire irradiation process
	  DO WHILE (TIMESTEP_A .LT. LENGTH_OF_FILE - 1)
		I = 0
C		Determine the number of steps in next irradiation cycle
		DO WHILE((IRRHISTRY(TIMESTEP_A+I+1,2)-
     &             IRRHISTRY(TIMESTEP_A+I,2) .GT. 1.0D-3))
		  I = I + 1
	      IF(TIMESTEP_A+I+1 .GE. LENGTH_OF_FILE) EXIT
		END DO
		NSTEPINCYCLE = I
C		Prepare arrays for storing data in this cycle        
		ALLOCATE (HCARDB(0:NSTEPINCYCLE,19))
		ALLOCATE (HSIGRB(0:NSTEPINCYCLE,NDIV))
		ALLOCATE (HSIGTB(0:NSTEPINCYCLE,NDIV))
		DO 710 I = 0, NSTEPINCYCLE
		  DO 720 J = 1, 19
			HCARDB(I,J) = 0.0D0
720		  CONTINUE
710		CONTINUE
C    Register initial non-zero quantities
		DO 730 I = 0, 5
		  T_PARTICLE(I) = IRRHISTRY(TIMESTEP_A,3)
730		CONTINUE
		TIMESTEP = 0
		OPERTIME = IRRHISTRY(TIMESTEP_A,1)*86400.0D0  !seconds
		IRRTIME = IRRHISTRY(TIMESTEP_A,2)		!days
		BURNUP = IRRHISTRY(TIMESTEP_A,5)/100.0D0  !convert from percent to absolute value
C
		HCARDB(TIMESTEP,1)  = OPERTIME
		HCARDB(TIMESTEP,4)  = IRRHISTRY(TIMESTEP_A,4) !Fluence from last cycle
		HCARDB(TIMESTEP,7)  = BURNUP
C
		HCARDB(TIMESTEP,9)  = T_PARTICLE(0)
		HCARDB(TIMESTEP,10) = T_PARTICLE(1)
		HCARDB(TIMESTEP,11) = T_PARTICLE(2)
		HCARDB(TIMESTEP,12) = T_PARTICLE(3)
		HCARDB(TIMESTEP,13) = T_PARTICLE(4)
		HCARDB(TIMESTEP,14) = T_PARTICLE(5)
		CALL GASRLS(T_PARTICLE, BURNUP, OPERTIME, DIFFUSION, PRESS)
		CALL SWELLU(T_PARTICLE(3), 0.0D0, IPYCD, IPYCBAF0,
     &                IPYCCRATE, SRDOT_IPYC, STDOT_IPYC)
		CALL SWELLU(T_PARTICLE(5), 0.0D0, OPYCD, OPYCBAF0,
     &                OPYCCRATE, SRDOT_OPYC, STDOT_OPYC)
		HCARDB(TIMESTEP,15) = PRESS
		HCARDB(TIMESTEP,16) = SRDOT_IPYC
		HCARDB(TIMESTEP,17) = STDOT_IPYC
		HCARDB(TIMESTEP,18) = SRDOT_OPYC
		HCARDB(TIMESTEP,19) = STDOT_OPYC
C
C		Register quantities in this irradiation cycle
		DO WHILE (TIMESTEP .LT. NSTEPINCYCLE)
		  TIMESTEP = TIMESTEP + 1
		  TIMESTEP_A = TIMESTEP_A + 1
		  OPERTIME = IRRHISTRY(TIMESTEP_A,1)*86400.0D0  !seconds
		  IRRTIME = IRRHISTRY(TIMESTEP_A,2)		!days
		  FLUENCE = IRRHISTRY(TIMESTEP_A,4)*1.0D21
		  FLUENCE_R = (IRRHISTRY(TIMESTEP_A,4)-HCARDB(0,4))*1.0D21
		  BURNUP = IRRHISTRY(TIMESTEP_A,5)/100.0D0  !convert from percent to absolute value
		  DO 740 I = 0, 5
			T_PARTICLE(I) = IRRHISTRY(TIMESTEP_A,3)
740		  CONTINUE
C
		  HCARDB(TIMESTEP,1) = OPERTIME
		  HCARDB(TIMESTEP,4) = FLUENCE/1.0D21
	      HCARDB(TIMESTEP,5) = FLUENCE_R/1.0D21
	      HCARDB(TIMESTEP,7) = BURNUP
	      HCARDB(TIMESTEP,9) = T_PARTICLE(0)
	      HCARDB(TIMESTEP,10) = T_PARTICLE(1)
	      HCARDB(TIMESTEP,11) = T_PARTICLE(2)
	      HCARDB(TIMESTEP,12) = T_PARTICLE(3)
	      HCARDB(TIMESTEP,13) = T_PARTICLE(4)
	      HCARDB(TIMESTEP,14) = T_PARTICLE(5)
		  CALL GASRLS(T_PARTICLE, BURNUP, OPERTIME, DIFFUSION, PRESS)
	      HCARDB(TIMESTEP,15) = PRESS
		  CALL SWELLU(T_PARTICLE(3), FLUENCE/1.0D21, IPYCD, IPYCBAF0,
     &                  IPYCCRATE, SRDOT_IPYC, STDOT_IPYC)
            CALL SWELLU(T_PARTICLE(5), FLUENCE/1.0D21, OPYCD, OPYCBAF0,
     &                  OPYCCRATE, SRDOT_OPYC, STDOT_OPYC)
	      HCARDB(TIMESTEP,16) = SRDOT_IPYC
		  HCARDB(TIMESTEP,17) = STDOT_IPYC
	      HCARDB(TIMESTEP,18) = SRDOT_OPYC
	      HCARDB(TIMESTEP,19) = STDOT_OPYC
		END DO
C		The history card for the current cycle is created.
C
C		First time stress analysis: unrestrained condition
C		Fit polynomials for swelling rate to be used in stress analysis
	    IF(TIMESTEP.EQ.1) THEN
	      ISWR(0) = HCARDB(TIMESTEP,16)
	      ISWT(0) = HCARDB(TIMESTEP,17)
	      OSWR(0) = HCARDB(TIMESTEP,18)
	      OSWT(0) = HCARDB(TIMESTEP,19)
	      DO 745 K = 1, NDEG
	        ISWR(K) = 0.0D0
	        ISWT(K) = 0.0D0
	        OSWR(K) = 0.0D0
	        OSWT(K) = 0.0D0
745         CONTINUE
          ELSE IF(TIMESTEP.EQ.2) THEN  !linear interpolation
	      ISWR(0) = (HCARDB(TIMESTEP,5)*HCARDB(1,16) - 
     &                 HCARDB(1,5)*HCARDB(TIMESTEP,16))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      ISWT(0) = (HCARDB(TIMESTEP,5)*HCARDB(1,17) - 
     &                 HCARDB(1,5)*HCARDB(TIMESTEP,17))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      OSWR(0) = (HCARDB(TIMESTEP,5)*HCARDB(1,18) - 
     &                 HCARDB(1,5)*HCARDB(TIMESTEP,18))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      OSWT(0) = (HCARDB(TIMESTEP,5)*HCARDB(1,19) - 
     &                 HCARDB(1,5)*HCARDB(TIMESTEP,19))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      ISWR(1) = (HCARDB(TIMESTEP,16) - HCARDB(1,16))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      ISWT(1) = (HCARDB(TIMESTEP,17) - HCARDB(1,17))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      OSWR(1) = (HCARDB(TIMESTEP,18) - HCARDB(1,18))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      OSWT(1) = (HCARDB(TIMESTEP,19) - HCARDB(1,19))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      IF(NDEG.GE.2) THEN
              DO 746 K = 2, NDEG
	          ISWR(K) = 0.0D0
	          ISWT(K) = 0.0D0
	          OSWR(K) = 0.0D0
	          OSWT(K) = 0.0D0
746           CONTINUE
            END IF
          ELSE IF((TIMESTEP-2).LT.NDEG) THEN
	      CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                  HCARDB(1:TIMESTEP,16),TIMESTEP-2,
     &                  ISWR(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                  HCARDB(1:TIMESTEP,17),TIMESTEP-2,
     &                  ISWT(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                  HCARDB(1:TIMESTEP,18),TIMESTEP-2,
     &                  OSWR(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                  HCARDB(1:TIMESTEP,19),TIMESTEP-2,
     &                  OSWT(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
            DO 747 K = TIMESTEP-1, NDEG
	        ISWR(K) = 0.0D0
	        ISWT(K) = 0.0D0
	        OSWR(K) = 0.0D0
	        OSWT(K) = 0.0D0
747         CONTINUE
          ELSE
		  CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                HCARDB(1:TIMESTEP,16),NDEG,ISWR,SSPOLY,STAT)
		  CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                HCARDB(1:TIMESTEP,17),NDEG,ISWT,SSPOLY,STAT)
		  CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                HCARDB(1:TIMESTEP,18),NDEG,OSWR,SSPOLY,STAT)
		  CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                HCARDB(1:TIMESTEP,19),NDEG,OSWT,SSPOLY,STAT)
		END IF
C
C		IF((.NOT.PARAMETRIC_STUDY) .AND. NOMINAL) THEN
C		Output swelling rate fitting of IPyC for inspection
C		  WRITE(IOUTSW, *) 'Radial and Tang. Swelling Coeffs: '
C	      WRITE(IOUTSW,628) ISWR(0), ISWR(1), ISWR(2), ISWR(3)
C	      WRITE(IOUTSW,628) ISWT(0), ISWT(1), ISWT(2), ISWT(3)
C	      WRITE(IOUTSW, *)
C		END IF
C		Restore apparent creep strains from last cycle
		E_IPYCREEP(1) = E_IPYCREEP0(1)
		E_IPYCREEP(2) = E_IPYCREEP0(2)
		E_OPYCREEP(1) = E_OPYCREEP0(1)
		E_OPYCREEP(2) = E_OPYCREEP0(2)		
C
C    Begin first time stress analysis with unrestrained condition
		DO 750 J = 0, TIMESTEP
		  IPYCE = E_PYC(IPYCD,IPYCBAF0,IPYCLC,HCARDB(J,4),
     &					HCARDB(J,12))
            OPYCE = E_PYC(OPYCD,OPYCBAF0,OPYCLC,HCARDB(J,4),
     &					HCARDB(J,13))
		  SICE = E_SIC(HCARDB(J,13))
C
	      CALL SSICREEP_PYC(HCARDB(J,4),IPYCD,HCARDB(J,12),E_IPYCREEP,
     &						IPYCREEP,IPYCCNU)
            CALL SSICREEP_PYC(HCARDB(J,4),OPYCD,HCARDB(J,13),E_OPYCREEP,
     &						OPYCREEP,OPYCCNU)
C
		  IF(J.GE.1) THEN
			DF = HCARDB(J,5) - HCARDB(J-1,5)
		  ELSE
			DF = HCARDB(J,5)
		  END IF
C		  Restore residual stresses, strains, and displacement if any
		  DO 751 K = 1, NDIV
			SIGR(K) = SIGR0(K)
			SIGT(K) = SIGT0(K)
			EPIR(K) = EPIR0(K)
			EPIT(K) = EPIT0(K)
			UR(K)   = UR0(K)
751		  CONTINUE
C		  Do mechanical analysis
		  CALL M_ANALYSIS(MCODE, HCARDB(J,15), PAMB, HCARDB(J,5),SIGR,
     &					  SIGT, EPIR, EPIT, UR)
		  HSIGRB(J,:) = SIGR
		  HSIGTB(J,:) = SIGT
C
C		  Accumulate apparent creep strains
		  D_FLUENCE = DF
		  CALL EPI_C('IPYC',SIGR,SIGT,IPYCREEP,
     &                 IPYCCNU,D_FLUENCE,E_IPYCREEP)
		  CALL EPI_C('OPYC',SIGR,SIGT,OPYCREEP,
     &                 OPYCCNU,D_FLUENCE,E_OPYCREEP)
C
C		  IF((.NOT.PARAMETRIC_STUDY) .AND. NOMINAL) THEN
C			Send stresses, strains, and displacement to output files
C			WRITE(IOUTSR,629) HCARDB(J,4)
C			WRITE(IOUTST,629) HCARDB(J,4)
C			WRITE(IOUTER,642) HCARDB(J,4)
C			WRITE(IOUTET,642) HCARDB(J,4)
C			WRITE(IOUTUR,629) HCARDB(J,4)
C			DO 752 K = 1, NDIV
C			  WRITE(IOUTSR, 629) SIGR(K)
C	          WRITE(IOUTST, 629) SIGT(K)
C			  WRITE(IOUTER, 641) EPIR(K)
C			  WRITE(IOUTET, 641) EPIT(K)
C			  WRITE(IOUTUR, 629) UR(K)
C752			CONTINUE
C			WRITE(IOUTSR,*)
C			WRITE(IOUTST,*)
C			WRITE(IOUTER,*)
C			WRITE(IOUTET,*)
C			WRITE(IOUTUR,*)
C			Output swelling rate of IPyC for inspection
C			WRITE(IOUTSW,628) HCARDB(J,4), HCARDB(J,16), HCARDB(J,17)
C		  END IF
C
750		CONTINUE
C    End of first time stress analysis
C
C		Restore apparent creep strains from last cycle
		E_IPYCREEP(1) = E_IPYCREEP0(1)
		E_IPYCREEP(2) = E_IPYCREEP0(2)
		E_OPYCREEP(1) = E_OPYCREEP0(1)
		E_OPYCREEP(2) = E_OPYCREEP0(2)
C		Calculate irradiated BAF, apparent creep strains and restrained swelling rates in PyC 
		DO 760 J = 1, TIMESTEP
		  IPYCBAFI = BAFI_PYC(IPYCBAF0, HCARDB(J, 4))
		  OPYCBAFI = BAFI_PYC(OPYCBAF0, HCARDB(J, 4))
C
		  CALL SSICREEP_PYC(HCARDB(J,4),IPYCD,HCARDB(J,12),E_IPYCREEP,
     &						IPYCREEP,IPYCCNU)
            CALL SSICREEP_PYC(HCARDB(J,4),OPYCD,HCARDB(J,13),E_OPYCREEP,
     &						OPYCREEP,OPYCCNU)
C
		  IF(J.GT.1) THEN
			D_FLUENCE = HCARDB(J,5) - HCARDB(J-1,5)
		  ELSE IF(J.EQ.1) THEN
		    D_FLUENCE = HCARDB(J,5)
		  ELSE
		    D_FLUENCE = 0.0D0
		  END IF
C
		  CALL EPI_C('IPYC',HSIGRB(J,1:NDIV),HSIGTB(J,1:NDIV),
     &               IPYCREEP,IPYCCNU,D_FLUENCE,E_IPYCREEP)
		  CALL EPI_C('OPYC',HSIGRB(J,1:NDIV),HSIGTB(J,1:NDIV),
     &               OPYCREEP,OPYCCNU,D_FLUENCE,E_OPYCREEP)
C
		  CALL SWELLR(IPYCBAFI,E_IPYCREEP,D_FLUENCE,
     &				  HCARDB(J,16),HCARDB(J,17))
		  CALL SWELLR(OPYCBAFI,E_OPYCREEP,D_FLUENCE,
     &				  HCARDB(J,18),HCARDB(J,19))
760		CONTINUE
C
C		Fit polynomials for swelling rate to be used in stress analysis
	    IF(TIMESTEP.EQ.1) THEN
	      ISWR(0) = HCARDB(TIMESTEP,16)
	      ISWT(0) = HCARDB(TIMESTEP,17)
	      OSWR(0) = HCARDB(TIMESTEP,18)
	      OSWT(0) = HCARDB(TIMESTEP,19)
	      DO 765 K = 1, NDEG
	        ISWR(K) = 0.0D0
	        ISWT(K) = 0.0D0
	        OSWR(K) = 0.0D0
	        OSWT(K) = 0.0D0
765         CONTINUE
          ELSE IF(TIMESTEP.EQ.2) THEN  !linear interpolation
	      ISWR(0) = (HCARDB(TIMESTEP,5)*HCARDB(1,16) - 
     &                 HCARDB(1,5)*HCARDB(TIMESTEP,16))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      ISWT(0) = (HCARDB(TIMESTEP,5)*HCARDB(1,17) - 
     &                 HCARDB(1,5)*HCARDB(TIMESTEP,17))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      OSWR(0) = (HCARDB(TIMESTEP,5)*HCARDB(1,18) - 
     &                 HCARDB(1,5)*HCARDB(TIMESTEP,18))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      OSWT(0) = (HCARDB(TIMESTEP,5)*HCARDB(1,19) - 
     &                 HCARDB(1,5)*HCARDB(TIMESTEP,19))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      ISWR(1) = (HCARDB(TIMESTEP,16) - HCARDB(1,16))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      ISWT(1) = (HCARDB(TIMESTEP,17) - HCARDB(1,17))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      OSWR(1) = (HCARDB(TIMESTEP,18) - HCARDB(1,18))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      OSWT(1) = (HCARDB(TIMESTEP,19) - HCARDB(1,19))/
     &                (HCARDB(TIMESTEP,5)-HCARDB(1,5))
	      IF(NDEG.GE.2) THEN
              DO 766 K = 2, NDEG
	          ISWR(K) = 0.0D0
	          ISWT(K) = 0.0D0
	          OSWR(K) = 0.0D0
	          OSWT(K) = 0.0D0
766           CONTINUE
            END IF
          ELSE IF((TIMESTEP-2).LT.NDEG) THEN
	      CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                  HCARDB(1:TIMESTEP,16),TIMESTEP-2,
     &                  ISWR(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                  HCARDB(1:TIMESTEP,17),TIMESTEP-2,
     &                  ISWT(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                  HCARDB(1:TIMESTEP,18),TIMESTEP-2,
     &                  OSWR(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
	      CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                  HCARDB(1:TIMESTEP,19),TIMESTEP-2,
     &                  OSWT(0:TIMESTEP-2),SSPOLY(1:TIMESTEP-1),STAT)
            DO 767 K = TIMESTEP-1, NDEG
	        ISWR(K) = 0.0D0
	        ISWT(K) = 0.0D0
	        OSWR(K) = 0.0D0
	        OSWT(K) = 0.0D0
767         CONTINUE
          ELSE
		  CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                HCARDB(1:TIMESTEP,16),NDEG,ISWR,SSPOLY,STAT)
		  CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                HCARDB(1:TIMESTEP,17),NDEG,ISWT,SSPOLY,STAT)
		  CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                HCARDB(1:TIMESTEP,18),NDEG,OSWR,SSPOLY,STAT)
		  CALL DRCURV(TIMESTEP,HCARDB(1:TIMESTEP,5),
     &                HCARDB(1:TIMESTEP,19),NDEG,OSWT,SSPOLY,STAT)
		END IF
C
		IF((.NOT.PARAMETRIC_STUDY) .AND. NOMINAL) THEN
		  WRITE(IOUTSW, *)
		  WRITE(IOUTSW, *) 'Radial and Tang. Swelling Coeffs for IPyC'
		  WRITE(IOUTSW,628) ISWR(0), ISWR(1), ISWR(2), ISWR(3)
	      WRITE(IOUTSW,628) ISWT(0), ISWT(1), ISWT(2), ISWT(3)
	      WRITE(IOUTSW, *)
		END IF
C
C    Begin second time stress analysis with restrained condition
C		Restore apparent creep strains from last cycle
		E_IPYCREEP(1) = E_IPYCREEP0(1)
		E_IPYCREEP(2) = E_IPYCREEP0(2)
		E_OPYCREEP(1) = E_OPYCREEP0(1)
		E_OPYCREEP(2) = E_OPYCREEP0(2)
C
		DO 770 J = 0, TIMESTEP
	      IPYCE = E_PYC(IPYCD,IPYCBAF0,IPYCLC,HCARDB(J,4),
     &					HCARDB(J,12))
            OPYCE = E_PYC(OPYCD,OPYCBAF0,OPYCLC,HCARDB(J,4),
     &					HCARDB(J,13))
            SICE = E_SIC(HCARDB(J,13))
C
	      CALL SSICREEP_PYC(HCARDB(J,4),IPYCD,HCARDB(J,12),E_IPYCREEP,
     &						IPYCREEP,IPYCCNU)
            CALL SSICREEP_PYC(HCARDB(J,4),OPYCD,HCARDB(J,13),E_OPYCREEP,
     &					    OPYCREEP,OPYCCNU)
C
		  IF(J.GE.1) THEN
		    DF = HCARDB(J,5) - HCARDB(J-1,5)
		  ELSE
		    DF = HCARDB(J,5)
		  END IF
C		  Restore residual stresses, strains, and displacement if any
		  DO 771 K = 1, NDIV
			SIGR(K) = SIGR0(K)
			SIGT(K) = SIGT0(K)
			EPIR(K) = EPIR0(K)
			EPIT(K) = EPIT0(K)
			UR(K)   = UR0(K)
771		  CONTINUE
C		  Do mechanical analysis
		  CALL M_ANALYSIS(MCODE, HCARDB(J,15), PAMB, HCARDB(J,5),SIGR,
     &					  SIGT, EPIR, EPIT, UR)
C
C		  Accumulate apparent creep strains
		  D_FLUENCE = DF
		  CALL EPI_C('IPYC',SIGR,SIGT,IPYCREEP,
     &                 IPYCCNU,D_FLUENCE,E_IPYCREEP)
		  CALL EPI_C('OPYC',SIGR,SIGT,OPYCREEP,
     &                 OPYCCNU,D_FLUENCE,E_OPYCREEP)
C
C           Update the fracture strength of PyC layers due to 
C           temp. variation and irradiation
		  CALL STRENGTH_PYC('IPYC','IRR',IPYCD,IPYCBAF0,HCARDB(J,4),
     &				HCARDB(J,12),SIGT,IPYCF,SIGFIPYC,IPYCM)
		  CALL STRENGTH_PYC('OPYC','IRR',OPYCD,OPYCBAF0,HCARDB(J,4),
     &				HCARDB(J,13),SIGT,OPYCF,SIGFOPYC,OPYCM)
C           Calculate the mean fracture strength of SiC layer
		  CALL STRENGTH_SIC('IRR',HCARDB(J,4),HCARDB(J,12),SIGT,
     &				SICF,SIGFSIC,SICM)
C           Register maximum and minimum tangential stresses in layers
		  DO 772 K = 1, NDIVI+2
			IF (SIGT(K).GT.SIGXIPYC) THEN 
			  SIGXIPYC = SIGT(K)
			  SIGFCIPYC = SIGFIPYC
			END IF
			IF (SIGT(K).LT.SIGMIPYC) SIGMIPYC = SIGT(K)
772         CONTINUE
            DO 773 K = 1, NDIVS+2
			IF (SIGT(NDIVI+2+K).GT.SIGXSIC) THEN
			  SIGXSIC = SIGT(NDIVI+2+K)
			  SIGFCSIC = SIGFSIC
			END IF
	        IF (SIGT(NDIVI+2+K).LT.SIGMSIC)  SIGMSIC = SIGT(NDIVI+2+K)
773		  CONTINUE
            DO 774 K = 1, NDIVO+2
	        IF (SIGT(NDIVI+NDIVS+4+K).GT.SIGXOPYC) THEN
                SIGXOPYC = SIGT(NDIVI+NDIVS+4+K)
			  SIGFCOPYC = SIGFOPYC
			END IF
	        IF (SIGT(NDIVI+NDIVS+4+K).LT.SIGMOPYC) 
     &           SIGMOPYC = SIGT(NDIVI+NDIVS+4+K)
774         CONTINUE
C
C		  Register mean tangential stresses in layers
		  SIGTIPYC = 0.0D0
		  DO 775 K = 1, NDIVI+2
	        SIGTIPYC = SIGTIPYC + SIGT(K)
775         CONTINUE
            SIGTIPYC = SIGTIPYC/FLOAT(NDIVI+2)
	      SIGTSIC = 0.0D0
	      DO 776 K = 1, NDIVS+2
	        SIGTSIC = SIGTSIC + SIGT(NDIVI+2+K)
776         CONTINUE
            SIGTSIC = SIGTSIC/FLOAT(NDIVS+2)
            SIGTOPYC = 0.0D0
	      DO 777 K = 1, NDIVO+2
	        SIGTOPYC = SIGTOPYC + SIGT(NDIVI+NDIVS+4+K)
777         CONTINUE
            SIGTOPYC = SIGTOPYC/FLOAT(NDIVO+2)
C    Register end-of-life mean tangential stresses in layers
		  SIGLIPYC = SIGTIPYC
		  SIGLOPYC = SIGTOPYC
		  SIGLSIC  = SIGTSIC
C
C    Evaluate fuel failure, do NOT do so if performing parametric study
		  IF((.NOT.PARAMETRIC_STUDY).AND.(RUNIRR .EQ. 'FAILURE')) THEN
	        DPD = 255.2*(OPERTIME/3600)*exp(-159.9/(0.008314*
     &        (T_PARTICLE(3)+273.15)))  
			CALL FAILURE(SIGR, SIGT, FAIL, FAILTYPE, DPD, N, 
     &                     OPERTIME, DT, MD)
		  END IF
C
C    Write to debug file of material strength data
		  IF((.NOT.PARAMETRIC_STUDY) .AND. DEBUG .AND. NOMINAL) THEN
		    WRITE(IDBG, 623) HCARDB(J,4), IPYCF, SIGFIPYC, SICF, 
     &			SIGFSIC, OPYCF, SIGFOPYC, KICSIC, KIIPYC, KIOPYC,
     &			KI1, KI2
		  END IF
C
C    Create histograms here
            IF(HISTOGRAM) THEN
C             Particle histogram
	        IF(PARFAILED.EQ.1) THEN
	          INDEX = INT(NHIS*HCARDB(J,1)/TIMELIMIT+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMP(INDEX,1) = HISTGRMP(INDEX,1) + 1
C
	          INDEX = INT(NHIS*(SIGTSIC-SIG_LOWER)/
     &                     (SIG_UPPER-SIG_LOWER)+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMP(INDEX,2) = HISTGRMP(INDEX,2) + 1
C
	          INDEX = INT(NHIS*HCARDB(J,4)/EOLFLU+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMP(INDEX,3) = HISTGRMP(INDEX,3) + 1
C
	          INDEX = INT(NHIS*HCARDB(J,7)/EOLBUP+0.5D0)
                IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMP(INDEX,4) = HISTGRMP(INDEX,4) + 1
	        END IF
C             IPyC layer histogram
	        IF(IPYCFAILED.EQ.1) THEN
	          INDEX = INT(NHIS*HCARDB(J,1)/TIMELIMIT+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMI(INDEX,1) = HISTGRMI(INDEX,1) + 1
C
	          INDEX = INT(NHIS*(SIGTIPYC-SIG_LOWER)/
     &                      (SIG_UPPER-SIG_LOWER)+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMI(INDEX,2) = HISTGRMI(INDEX,2) + 1
C
	          INDEX = INT(NHIS*HCARDB(J,4)/EOLFLU+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMI(INDEX,3) = HISTGRMI(INDEX,3) + 1
C
	          INDEX = INT(NHIS*HCARDB(J,7)/EOLBUP+0.5D0)
                IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMI(INDEX,4) = HISTGRMI(INDEX,4) + 1
	        END IF
C             OPyC layer histogram
	        IF(OPYCFAILED.EQ.1) THEN
	          INDEX = INT(NHIS*HCARDB(J,1)/TIMELIMIT+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMO(INDEX,1) = HISTGRMO(INDEX,1) + 1
C
	          INDEX = INT(NHIS*(SIGTOPYC-SIG_LOWER)/
     &                      (SIG_UPPER-SIG_LOWER)+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMO(INDEX,2) = HISTGRMO(INDEX,2) + 1
C
	          INDEX = INT(NHIS*HCARDB(J,4)/EOLFLU+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMO(INDEX,3) = HISTGRMO(INDEX,3) + 1
C
	          INDEX = INT(NHIS*HCARDB(J,7)/EOLBUP+0.5D0)
                IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMO(INDEX,4) = HISTGRMO(INDEX,4) + 1
	        END IF
C             SiC layer histogram
	        IF(SICFAILED.EQ.1) THEN
	          INDEX = INT(NHIS*HCARDB(J,1)/TIMELIMIT+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMS(INDEX,1) = HISTGRMS(INDEX,1) + 1
C
	          INDEX = INT(NHIS*(SIGTSIC-SIG_LOWER)/
     &                      (SIG_UPPER-SIG_LOWER)+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMS(INDEX,2) = HISTGRMS(INDEX,2) + 1
C
	          INDEX = INT(NHIS*HCARDB(J,4)/EOLFLU+0.5D0)
	          IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMS(INDEX,3) = HISTGRMS(INDEX,3) + 1
C
	          INDEX = INT(NHIS*HCARDB(J,7)/EOLBUP+0.5D0)
                IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &          HISTGRMS(INDEX,4) = HISTGRMS(INDEX,4) + 1
	        END IF
	      END IF
C    End of creating histograms
C
	      IF((.NOT.PARAMETRIC_STUDY) .AND. NOMINAL) THEN
C             Send stresses, strains, and displacement to output files
              WRITE(ITEST,*)  HCARDB(J,1), HCARDB(J,4), HCARDB(J,7)
	        WRITE(IOUTSR,629) HCARDB(J,4)
	        WRITE(IOUTST,629) HCARDB(J,4)
	        WRITE(IOUTER,642) HCARDB(J,4)
	        WRITE(IOUTET,642) HCARDB(J,4)
	        WRITE(IOUTUR,629) HCARDB(J,4)
              DO 778 K = 1, NDIV
	          WRITE(IOUTSR, 629) SIGR(K)
	          WRITE(IOUTST, 629) SIGT(K)
	          WRITE(IOUTER, 641) EPIR(K)
	          WRITE(IOUTET, 641) EPIT(K)
	          WRITE(IOUTUR, 629) UR(K)
778           CONTINUE
              WRITE(IOUTSR,*)
              WRITE(IOUTST,*)
              WRITE(IOUTER,*)
              WRITE(IOUTET,*)
              WRITE(IOUTUR,*)
C             Output swelling rate of IPyC for inspection
	        WRITE(IOUTSW,628) HCARDB(J,4), HCARDB(J,16), HCARDB(J,17)
	      END IF
C
		  IF(FAIL) THEN
C             write information about the failed particle
			IF(.NOT.NOMINAL) THEN
			  WRITE(IOUTFL, 646) PARFAIL, FMODE, R1, R2, R3, R4, R5,
     &						   BDEN, IPYCD, OPYCD, IPYCBAF0, OPYCBAF0,
     &						   SIGFIPYC, SIGFOPYC, SIGFSIC, KICSIC,
     &						   SIGLIPYC, SIGLOPYC, SIGLSIC,
     &						   HCARDB(J,1)/86400.0D0, HCARDB(J,4),
     &						   HCARDB(J,7),HCARDB(J,11),HCARDB(J,15),
     &						   FAILUREPATH
			END IF
C			Release the allocated arrays
			DEALLOCATE (HCARDB)
			DEALLOCATE (HSIGRB)
			DEALLOCATE (HSIGTB)
		    GO TO 990       
	      END IF
C
770		CONTINUE
C
C         Save residual stresses, strains, and displacement
		DO 780 J = 1, NDIV
	      SIGR0(J) = SIGR(J)
	      SIGT0(J) = SIGT(J)
		  EPIR0(J) = EPIR(J)
		  EPIT0(J) = EPIT(J)
		  UR0(J)   = UR(J)
780       CONTINUE
C		Save apparent creep strains gained in this cycle
		E_IPYCREEP0(1) = E_IPYCREEP(1)
		E_IPYCREEP0(2) = E_IPYCREEP(2)
		E_OPYCREEP0(1) = E_OPYCREEP(1)
		E_OPYCREEP0(2) = E_OPYCREEP(2)
C
C         maybe do thermal stress analysis and gas release calculation here.
		DEALLOCATE (HCARDB)
		DEALLOCATE (HSIGRB)
		DEALLOCATE (HSIGTB)
C    Shift TIMESTEP_A from the end of one cycle to the start of next cycle
		TIMESTEP_A = TIMESTEP_A + 1
	  END DO
C                                                                      *
C                                                                      *
C***********************************************************************
C***********************************************************************
C                                                                      *
C                                                                      *
	ELSE    !PSWITCH choice is important here
C
	  DO 1110 K = 0, 5
	    T_PARTICLE(K) = T_IRR
1110    CONTINUE
        DO 1120 I = 0, LENGTH_OF_FILE
	    DO 1130 J = 1, 19
	      HCARDB(I,J) = 0.0D0
1130      CONTINUE
1120    CONTINUE
C    Register initial non-zero quantities
        TIMESTEP = 0
	  HCARDB(TIMESTEP,9)  = T_PARTICLE(0)
	  HCARDB(TIMESTEP,10) = T_PARTICLE(1)
	  HCARDB(TIMESTEP,11) = T_PARTICLE(2)
	  HCARDB(TIMESTEP,12) = T_PARTICLE(3)
	  HCARDB(TIMESTEP,13) = T_PARTICLE(4)
	  HCARDB(TIMESTEP,14) = T_PARTICLE(5)
        CALL SWELLU(T_PARTICLE(3), FLUENCE/1.0D21, IPYCD, IPYCBAF0,
     &             IPYCCRATE, SRDOT_IPYC, STDOT_IPYC)
        CALL SWELLU(T_PARTICLE(5), FLUENCE/1.0D21, OPYCD, OPYCBAF0,
     &             OPYCCRATE, SRDOT_OPYC, STDOT_OPYC)
	  HCARDB(TIMESTEP,16) = SRDOT_IPYC
	  HCARDB(TIMESTEP,17) = STDOT_IPYC
	  HCARDB(TIMESTEP,18) = SRDOT_OPYC
	  HCARDB(TIMESTEP,19) = STDOT_OPYC
C
C    Register quantities during irradiation
        DO WHILE (((FLUENCE+FFLUX*DT-FLUENCE/1.0D6)/1.0D21).LE.EOLFLU)
	    TIMESTEP = TIMESTEP + 1
	    OPERTIME = OPERTIME + DT
          FLUENCE = FLUENCE + FFLUX*DT
	    BURNUP = BURNUP + BPRATE*(DT/86400.0D0)/BU_CONV   !convert sec to day
C         BU_CONV = (190.0D0*1.602D-13*NA)/(86400.0D0*AWU235)
	    HCARDB(TIMESTEP,1) = OPERTIME
          HCARDB(TIMESTEP,4) = FLUENCE/1.0D21
	    HCARDB(TIMESTEP,5) = FLUENCE/1.0D21
	    HCARDB(TIMESTEP,7) = BURNUP
	    HCARDB(TIMESTEP,9) = T_PARTICLE(0)
	    HCARDB(TIMESTEP,10) = T_PARTICLE(1)
	    HCARDB(TIMESTEP,11) = T_PARTICLE(2)
	    HCARDB(TIMESTEP,12) = T_PARTICLE(3)
	    HCARDB(TIMESTEP,13) = T_PARTICLE(4)
	    HCARDB(TIMESTEP,14) = T_PARTICLE(5)
	    CALL GASRLS(T_PARTICLE, BURNUP, OPERTIME, DIFFUSION, PRESS)
	    HCARDB(TIMESTEP,15) = PRESS
          CALL SWELLU(T_PARTICLE(3), FLUENCE/1.0D21, IPYCD, IPYCBAF0,
     &               IPYCCRATE, SRDOT_IPYC, STDOT_IPYC)
          CALL SWELLU(T_PARTICLE(5), FLUENCE/1.0D21, OPYCD, OPYCBAF0,
     &               OPYCCRATE, SRDOT_OPYC, STDOT_OPYC)
	    HCARDB(TIMESTEP,16) = SRDOT_IPYC
	    HCARDB(TIMESTEP,17) = STDOT_IPYC
	    HCARDB(TIMESTEP,18) = SRDOT_OPYC
	    HCARDB(TIMESTEP,19) = STDOT_OPYC
	    IF(TIMESTEP.EQ.LENGTH_OF_FILE) EXIT
        END DO
C
        CALL DRCURV(TIMESTEP,HCARDB(0:TIMESTEP,5),
     &              HCARDB(0:TIMESTEP,16),NDEG,ISWR,SSPOLY,STAT)
        CALL DRCURV(TIMESTEP,HCARDB(0:TIMESTEP,5),
     &              HCARDB(0:TIMESTEP,17),NDEG,ISWT,SSPOLY,STAT)
        CALL DRCURV(TIMESTEP,HCARDB(0:TIMESTEP,5),
     &              HCARDB(0:TIMESTEP,18),NDEG,OSWR,SSPOLY,STAT)
        CALL DRCURV(TIMESTEP,HCARDB(0:TIMESTEP,5),
     &              HCARDB(0:TIMESTEP,19),NDEG,OSWT,SSPOLY,STAT)
C
	  IF((.NOT.PARAMETRIC_STUDY) .AND. NOMINAL) THEN
C    Output swelling rate fitting of IPyC for inspection
		WRITE(IOUTSW, *) 'Radial and Tang. Swelling Coeffs: '
	    WRITE(IOUTSW,628) ISWR(0), ISWR(1), ISWR(2), ISWR(3)
	    WRITE(IOUTSW,628) ISWT(0), ISWT(1), ISWT(2), ISWT(3)
	    WRITE(IOUTSW, *)
	  END IF
C
C    Restore apparent creep strains from last cycle
		E_IPYCREEP(1) = E_IPYCREEP0(1)
		E_IPYCREEP(2) = E_IPYCREEP0(2)
		E_OPYCREEP(1) = E_OPYCREEP0(1)
		E_OPYCREEP(2) = E_OPYCREEP0(2)
C
C    Begin first time stress analysis with unrestrained condition
        DO 1145 J = 0, TIMESTEP
          IPYCE = E_PYC(IPYCD,IPYCBAF0,IPYCLC,HCARDB(J,4),HCARDB(J,12))
          OPYCE = E_PYC(OPYCD,OPYCBAF0,OPYCLC,HCARDB(J,4),HCARDB(J,13))
          SICE = E_SIC(HCARDB(J,13))
C
	    CALL SSICREEP_PYC(HCARDB(J,4),IPYCD,HCARDB(J,12),E_IPYCREEP,
     &					  IPYCREEP,IPYCCNU)
          CALL SSICREEP_PYC(HCARDB(J,4),OPYCD,HCARDB(J,13),E_OPYCREEP,
     &					  OPYCREEP,OPYCCNU)
C
		IF(J.GE.1) THEN
		  DF = HCARDB(J,5) - HCARDB(J-1,5)
		ELSE
		  DF = HCARDB(J,5)
		END IF
C    Restore residual stresses, strains, and displacement if any
	    DO 1141 K = 1, NDIV
            SIGR(K) = SIGR0(K)
	      SIGT(K) = SIGT0(K)
		  EPIR(K) = EPIR0(K)
		  EPIT(K) = EPIT0(K)
		  UR(K)   = UR0(K)
1141      CONTINUE
C    Do mechanical analyses
	    CALL M_ANALYSIS(MCODE, HCARDB(J,15), PAMB, HCARDB(J,5), SIGR,
     &					SIGT, EPIR, EPIT, UR)
		HSIGRB(J,:) = SIGR
		HSIGTB(J,:) = SIGT
C
C    Accumulate apparent creep strains
		D_FLUENCE = DF
		CALL EPI_C('IPYC',SIGR,SIGT,IPYCREEP,
     &               IPYCCNU,D_FLUENCE,E_IPYCREEP)
		CALL EPI_C('OPYC',SIGR,SIGT,OPYCREEP,
     &               OPYCCNU,D_FLUENCE,E_OPYCREEP)
C
	    IF((.NOT.PARAMETRIC_STUDY) .AND. NOMINAL) THEN
C    Send stresses, strains, and displacement to output files
	      WRITE(IOUTSR,629) HCARDB(J,4)
	      WRITE(IOUTST,629) HCARDB(J,4)
	      WRITE(IOUTER,642) HCARDB(J,4)
	      WRITE(IOUTET,642) HCARDB(J,4)
	      WRITE(IOUTUR,629) HCARDB(J,4)
            DO 1142 K = 1, NDIV
	        WRITE(IOUTSR, 629) SIGR(K)
	        WRITE(IOUTST, 629) SIGT(K)
	        WRITE(IOUTER, 641) EPIR(K)
	        WRITE(IOUTET, 641) EPIT(K)
	        WRITE(IOUTUR, 629) UR(K)
1142        CONTINUE
            WRITE(IOUTSR,*)
            WRITE(IOUTST,*)
            WRITE(IOUTER,*)
            WRITE(IOUTET,*)
            WRITE(IOUTUR,*)
C    Output swelling rate of IPyC for inspection
	      WRITE(IOUTSW,628) HCARDB(J,4), HCARDB(J,16), HCARDB(J,17)
	    END IF
C
1145    CONTINUE
C    End of first time stress analysis
C
C    Restore apparent creep strains from last cycle
	  E_IPYCREEP(1) = E_IPYCREEP0(1)
	  E_IPYCREEP(2) = E_IPYCREEP0(2)
	  E_OPYCREEP(1) = E_OPYCREEP0(1)
	  E_OPYCREEP(2) = E_OPYCREEP0(2)
C
C    Calculate irradiated BAF, apparent creep strains and restrained swelling rates in PyC
	  DO 1143 J = 1, TIMESTEP
		IPYCBAFI = BAFI_PYC(IPYCBAF0, HCARDB(J, 4))
		OPYCBAFI = BAFI_PYC(OPYCBAF0, HCARDB(J, 4))
C
	    CALL SSICREEP_PYC(HCARDB(J,4),IPYCD,HCARDB(J,12),E_IPYCREEP,
     &					  IPYCREEP,IPYCCNU)
          CALL SSICREEP_PYC(HCARDB(J,4),OPYCD,HCARDB(J,13),E_OPYCREEP,
     &					  OPYCREEP,OPYCCNU)
C
		IF(J.GT.1) THEN
		  D_FLUENCE = HCARDB(J,5) - HCARDB(J-1,5)
		ELSE IF(J.EQ.1) THEN
		  D_FLUENCE = HCARDB(J,5)
		ELSE
		  D_FLUENCE = 0.0D0
		END IF
C
		CALL EPI_C('IPYC',HSIGRB(J,1:NDIV),HSIGTB(J,1:NDIV),
     &               IPYCREEP,IPYCCNU,D_FLUENCE,E_IPYCREEP)
		CALL EPI_C('OPYC',HSIGRB(J,1:NDIV),HSIGTB(J,1:NDIV),
     &               OPYCREEP,OPYCCNU,D_FLUENCE,E_OPYCREEP)
C
		CALL SWELLR(IPYCBAFI,E_IPYCREEP,D_FLUENCE,
     &				HCARDB(J,16),HCARDB(J,17))
		CALL SWELLR(OPYCBAFI,E_OPYCREEP,D_FLUENCE,
     &				HCARDB(J,18),HCARDB(J,19))
1143	  CONTINUE
C
        CALL DRCURV(TIMESTEP,HCARDB(0:TIMESTEP,5),
     &              HCARDB(0:TIMESTEP,16),NDEG,ISWR,SSPOLY,STAT)
        CALL DRCURV(TIMESTEP,HCARDB(0:TIMESTEP,5),
     &              HCARDB(0:TIMESTEP,17),NDEG,ISWT,SSPOLY,STAT)
        CALL DRCURV(TIMESTEP,HCARDB(0:TIMESTEP,5),
     &              HCARDB(0:TIMESTEP,18),NDEG,OSWR,SSPOLY,STAT)
        CALL DRCURV(TIMESTEP,HCARDB(0:TIMESTEP,5),
     &              HCARDB(0:TIMESTEP,19),NDEG,OSWT,SSPOLY,STAT)
	  IF((.NOT.PARAMETRIC_STUDY) .AND. NOMINAL) THEN
C    Output swelling rate fitting of IPyC for inspection
C		
	    WRITE(IOUTSW, *)
		WRITE(IOUTSW, *) 'Radial and Tang. Swelling Coeffs: for IPyC'
	    WRITE(IOUTSW,628) ISWR(0), ISWR(1), ISWR(2), ISWR(3)
	    WRITE(IOUTSW,628) ISWT(0), ISWT(1), ISWT(2), ISWT(3)
	    WRITE(IOUTSW, *)
	  END IF
C
C    Begin second time stress analysis with restrained condition
C    Restore apparent creep strains from last cycle
	  E_IPYCREEP(1) = E_IPYCREEP0(1)
	  E_IPYCREEP(2) = E_IPYCREEP0(2)
	  E_OPYCREEP(1) = E_OPYCREEP0(1)
	  E_OPYCREEP(2) = E_OPYCREEP0(2)
C
        DO 1140 J = 0, TIMESTEP
          IPYCE = E_PYC(IPYCD,IPYCBAF0,IPYCLC,HCARDB(J,4),HCARDB(J,12))
          OPYCE = E_PYC(OPYCD,OPYCBAF0,OPYCLC,HCARDB(J,4),HCARDB(J,13))
          SICE = E_SIC(HCARDB(J,13))
C
	    CALL SSICREEP_PYC(HCARDB(J,4),IPYCD,HCARDB(J,12),E_IPYCREEP,
     &					  IPYCREEP,IPYCCNU)
          CALL SSICREEP_PYC(HCARDB(J,4),OPYCD,HCARDB(J,13),E_OPYCREEP,
     &					  OPYCREEP,OPYCCNU)
C
		IF(J.GE.1) THEN
		  DF = HCARDB(J,5) - HCARDB(J-1,5)
		ELSE
		  DF = HCARDB(J,5)
		END IF
C    Restore residual stresses, strains, and displacement if any
	    DO 1150 K = 1, NDIV
            SIGR(K) = SIGR0(K)
	      SIGT(K) = SIGT0(K)
		  EPIR(K) = EPIR0(K)
		  EPIT(K) = EPIT0(K)
		  UR(K)   = UR0(K)
1150      CONTINUE
C    Do mechanical analyses
	    CALL M_ANALYSIS(MCODE, HCARDB(J,15), PAMB, HCARDB(J,5), SIGR,
     &					SIGT, EPIR, EPIT, UR)
C
C    Accumulate apparent creep strains
		D_FLUENCE = DF
		CALL EPI_C('IPYC',SIGR,SIGT,IPYCREEP,
     &               IPYCCNU,D_FLUENCE,E_IPYCREEP)
		CALL EPI_C('OPYC',SIGR,SIGT,OPYCREEP,
     &               OPYCCNU,D_FLUENCE,E_OPYCREEP)
C
C    Update the fracture strength of PyC layers due to temp. variation and irradiation
		CALL STRENGTH_PYC('IPYC','IRR',IPYCD,IPYCBAF0,HCARDB(J,4),
     &				HCARDB(J,12),SIGT,IPYCF,SIGFIPYC,IPYCM)
		CALL STRENGTH_PYC('OPYC','IRR',OPYCD,OPYCBAF0,HCARDB(J,4),
     &				HCARDB(J,13),SIGT,OPYCF,SIGFOPYC,OPYCM)
C    Calculate the mean fracture strength of SiC layer
		CALL STRENGTH_SIC('IRR',HCARDB(J,4),HCARDB(J,12),SIGT,
     &					  SICF,SIGFSIC,SICM)
C
C    Register maximum tangential stresses in layers
          DO 1160 K = 1, NDIVI+2
            IF (SIGT(K).GT.SIGXIPYC) THEN
		    SIGXIPYC = SIGT(K)
			SIGFCIPYC = SIGFIPYC
		  END IF
            IF (SIGT(K).LT.SIGMIPYC) SIGMIPYC = SIGT(K)
1160       CONTINUE
          DO 1170 K = 1, NDIVS+2
	      IF (SIGT(NDIVI+2+K).GT.SIGXSIC) THEN
		    SIGXSIC = SIGT(NDIVI+2+K)
			SIGFCSIC = SIGFSIC
		  END IF
	      IF (SIGT(NDIVI+2+K).LT.SIGMSIC)  SIGMSIC = SIGT(NDIVI+2+K)
1170      CONTINUE
          DO 1180 K = 1, NDIVO+2
	      IF (SIGT(NDIVI+NDIVS+4+K).GT.SIGXOPYC) THEN
              SIGXOPYC = SIGT(NDIVI+NDIVS+4+K)
			SIGFCOPYC = SIGFOPYC
		  END IF
	      IF (SIGT(NDIVI+NDIVS+4+K).LT.SIGMOPYC) 
     &         SIGMOPYC = SIGT(NDIVI+NDIVS+4+K)
1180      CONTINUE
C
C    Register mean tangential stresses in layers
	    SIGTIPYC = 0.0D0
	    DO 284 K = 1, NDIVI+2
	      SIGTIPYC = SIGTIPYC + SIGT(K)
284       CONTINUE
          SIGTIPYC = SIGTIPYC/FLOAT(NDIVI+2)
	    SIGTSIC = 0.0D0
	    DO 285 K = 1, NDIVS+2
	      SIGTSIC = SIGTSIC + SIGT(NDIVI+2+K)
285       CONTINUE
          SIGTSIC = SIGTSIC/FLOAT(NDIVS+2)
          SIGTOPYC = 0.0D0
	    DO 286 K = 1, NDIVO+2
	      SIGTOPYC = SIGTOPYC + SIGT(NDIVI+NDIVS+4+K)
286       CONTINUE
          SIGTOPYC = SIGTOPYC/FLOAT(NDIVO+2)
C    Register end-of-life tangential stresses in layers
		SIGLIPYC = SIGTIPYC
		SIGLOPYC = SIGTOPYC
		SIGLSIC  = SIGTSIC
C
C    Evaluate fuel failure, do NOT do so if performing parametric study
		IF((.NOT.PARAMETRIC_STUDY).AND.(RUNIRR .EQ. 'FAILURE')) THEN
              DPD = 255.2*(OPERTIME/3600)*exp(-159.9/(0.008314*
     &        (T_PARTICLE(3)+273.15)))     
			CALL FAILURE(SIGR, SIGT, FAIL, FAILTYPE, DPD, N, 
     &                     OPERTIME, DT, MD)
		END IF
C
C    Write to debug file of material strength data
		IF((.NOT.PARAMETRIC_STUDY) .AND. DEBUG .AND. NOMINAL) THEN
		  WRITE(IDBG, 623) HCARDB(J,4), IPYCF, SIGFIPYC, SICF, 
     &			SIGFSIC, OPYCF, SIGFOPYC, KICSIC, KIIPYC, KIOPYC,
     &			KI1, KI2
		END IF
C
C    Create histograms here
          IF(HISTOGRAM) THEN
C           Particle histogram
	      IF(PARFAILED.EQ.1) THEN
	        INDEX = INT(NHIS*HCARDB(J,1)/TIMELIMIT+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMP(INDEX,1) = HISTGRMP(INDEX,1) + 1
C
	        INDEX = INT(NHIS*(SIGTSIC-SIG_LOWER)/
     &                   (SIG_UPPER-SIG_LOWER)+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMP(INDEX,2) = HISTGRMP(INDEX,2) + 1
C
	        INDEX = INT(NHIS*HCARDB(J,4)/EOLFLU+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMP(INDEX,3) = HISTGRMP(INDEX,3) + 1
C
	        INDEX = INT(NHIS*HCARDB(J,7)/EOLBUP+0.5D0)
              IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMP(INDEX,4) = HISTGRMP(INDEX,4) + 1
	      END IF
C           IPyC layer histogram
	      IF(IPYCFAILED.EQ.1) THEN
	        INDEX = INT(NHIS*HCARDB(J,1)/TIMELIMIT+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMI(INDEX,1) = HISTGRMI(INDEX,1) + 1
C
	        INDEX = INT(NHIS*(SIGTIPYC-SIG_LOWER)/
     &                    (SIG_UPPER-SIG_LOWER)+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMI(INDEX,2) = HISTGRMI(INDEX,2) + 1
C
	        INDEX = INT(NHIS*HCARDB(J,4)/EOLFLU+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMI(INDEX,3) = HISTGRMI(INDEX,3) + 1
C
	        INDEX = INT(NHIS*HCARDB(J,7)/EOLBUP+0.5D0)
              IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMI(INDEX,4) = HISTGRMI(INDEX,4) + 1
	      END IF
C             OPyC layer histogram
	      IF(OPYCFAILED.EQ.1) THEN
	        INDEX = INT(NHIS*HCARDB(J,1)/TIMELIMIT+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMO(INDEX,1) = HISTGRMO(INDEX,1) + 1
C
	        INDEX = INT(NHIS*(SIGTOPYC-SIG_LOWER)/
     &                    (SIG_UPPER-SIG_LOWER)+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMO(INDEX,2) = HISTGRMO(INDEX,2) + 1
C
	        INDEX = INT(NHIS*HCARDB(J,4)/EOLFLU+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMO(INDEX,3) = HISTGRMO(INDEX,3) + 1
C
	        INDEX = INT(NHIS*HCARDB(J,7)/EOLBUP+0.5D0)
              IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMO(INDEX,4) = HISTGRMO(INDEX,4) + 1
	      END IF
C           SiC layer histogram
	      IF(SICFAILED.EQ.1) THEN
	        INDEX = INT(NHIS*HCARDB(J,1)/TIMELIMIT+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMS(INDEX,1) = HISTGRMS(INDEX,1) + 1
C
	        INDEX = INT(NHIS*(SIGTSIC-SIG_LOWER)/
     &                    (SIG_UPPER-SIG_LOWER)+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMS(INDEX,2) = HISTGRMS(INDEX,2) + 1
C
	        INDEX = INT(NHIS*HCARDB(J,4)/EOLFLU+0.5D0)
	        IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMS(INDEX,3) = HISTGRMS(INDEX,3) + 1
C
	        INDEX = INT(NHIS*HCARDB(J,7)/EOLBUP+0.5D0)
              IF((INDEX.GE.0).AND.(INDEX.LE.NHIS))
     &        HISTGRMS(INDEX,4) = HISTGRMS(INDEX,4) + 1
	      END IF
	    END IF
C    End of creating histograms
C
	    IF((.NOT.PARAMETRIC_STUDY) .AND. NOMINAL) THEN
C    Send stresses, strains, and displacement to output files
            WRITE(ITEST,*)  HCARDB(J,1), HCARDB(J,4), HCARDB(J,7)
	      WRITE(IOUTSR,629) HCARDB(J,4)
	      WRITE(IOUTST,629) HCARDB(J,4)
	      WRITE(IOUTER,642) HCARDB(J,4)
	      WRITE(IOUTET,642) HCARDB(J,4)
	      WRITE(IOUTUR,629) HCARDB(J,4)
            DO 1190 K = 1, NDIV
	        WRITE(IOUTSR, 629) SIGR(K)
	        WRITE(IOUTST, 629) SIGT(K)
	        WRITE(IOUTER, 641) EPIR(K)
	        WRITE(IOUTET, 641) EPIT(K)
	        WRITE(IOUTUR, 629) UR(K)
1190        CONTINUE
            WRITE(IOUTSR,*)
            WRITE(IOUTST,*)
            WRITE(IOUTER,*)
            WRITE(IOUTET,*)
            WRITE(IOUTUR,*)
C    Output swelling rate of IPyC for inspection
	      WRITE(IOUTSW,628) HCARDB(J,4), HCARDB(J,16), HCARDB(J,17)
	    END IF
C
          IF(FAIL) THEN
C         write information about the failed particle
		  IF(.NOT.NOMINAL) THEN
			WRITE(IOUTFL, 646) PARFAIL, FMODE, R1, R2, R3, R4, R5,
     &						 BDEN, IPYCD, OPYCD, IPYCBAF0, OPYCBAF0,
     &						 SIGFIPYC, SIGFOPYC, SIGFSIC, KICSIC,
     &						 SIGLIPYC, SIGLOPYC, SIGLSIC,
     &						 HCARDB(J,1)/86400.0D0, HCARDB(J,4),
     &						 HCARDB(J,7),HCARDB(J,11),HCARDB(J,15),
     &						 FAILUREPATH
		  END IF
		  GO TO 990       
	    END IF
C
1140    CONTINUE
C                                                                      *
C                                                                      *
C***********************************************************************
C***********************************************************************
      ENDIF
C  Send to the terminal intermediate results per NBURP of particles. 
C  NCASES of particles are divided into groups NBURP particles for statistical purpose.
990   IF (MOD(N,NBURP) .EQ. 0) THEN
        TIME = 0.0
        WRITE(ITERM,605) N, TIME, PARFAIL, SICFAIL, IPYCFAIL, OPYCFAIL
        WRITE(IOUT,605) N, TIME, PARFAIL, SICFAIL, IPYCFAIL, OPYCFAIL
C  Gain statistical information on layer/particle failures every group of NBURP particles
	  PROBIF = PROBIF + DBLE(IPYCFAIL - IPYCFAIL0)
	  PROBSF = PROBSF + DBLE(SICFAIL - SICFAIL0)
	  PROBOF = PROBOF + DBLE(OPYCFAIL - OPYCFAIL0)
	  PROBPF = PROBPF + DBLE(PARFAIL - PARFAIL0)
	  CONFIDIF = CONFIDIF + DBLE((IPYCFAIL - IPYCFAIL0)**2)
	  CONFIDSF = CONFIDSF + DBLE((SICFAIL - SICFAIL0)**2)
	  CONFIDOF = CONFIDOF + DBLE((OPYCFAIL - OPYCFAIL0)**2)
	  CONFIDPF = CONFIDPF + DBLE((PARFAIL - PARFAIL0)**2)
	  IPYCFAIL0 = IPYCFAIL
	  SICFAIL0 = SICFAIL
	  OPYCFAIL0 = OPYCFAIL
	  PARFAIL0 = PARFAIL
      END IF
C  Gain statistical information on stresses
      SIGXBARI = SIGXBARI + SIGXIPYC
	SIGXBARS = SIGXBARS + SIGXSIC
	SIGXBARO = SIGXBARO + SIGXOPYC
	SIGXVARI = SIGXVARI + SIGXIPYC**2
	SIGXVARS = SIGXVARS + SIGXSIC**2
	SIGXVARO = SIGXVARO + SIGXOPYC**2
      SIGMBARI = SIGMBARI + SIGMIPYC
	SIGMBARS = SIGMBARS + SIGMSIC
	SIGMBARO = SIGMBARO + SIGMOPYC
	SIGMVARI = SIGMVARI + SIGMIPYC**2
	SIGMVARS = SIGMVARS + SIGMSIC**2
	SIGMVARO = SIGMVARO + SIGMOPYC**2
      SIGLBARI = SIGLBARI + SIGLIPYC
	SIGLBARS = SIGLBARS + SIGLSIC
	SIGLBARO = SIGLBARO + SIGLOPYC
	SIGLVARI = SIGLVARI + SIGLIPYC**2
	SIGLVARS = SIGLVARS + SIGLSIC**2
	SIGLVARO = SIGLVARO + SIGLOPYC**2
C
1000  CONTINUE  ! End of Monte Carlo Loop: Go get another particle
C  Release memory for data stored arrays
      IF(PSWITCH.EQ.1) THEN
        DEALLOCATE (POWDISTR)
	  DEALLOCATE (CHANNELS)
	  DEALLOCATE (BLOCKS)
	  DEALLOCATE (PATH)
	  DEALLOCATE (BLOCKMAP)
	  DEALLOCATE (HCARD)
	  DEALLOCATE (HSIGR)
	  DEALLOCATE (HSIGT)
	ELSE IF(PSWITCH.EQ.2) THEN
	  DEALLOCATE (IRRHISTRY)
	ELSE IF(PSWITCH.EQ.3) THEN
	  DEALLOCATE (HCARDB)
	  DEALLOCATE (HSIGRB)
	  DEALLOCATE (HSIGTB)
	END IF
C
C  Do statistics on stresses
      SIGXBARI = SIGXBARI/NCASES
      SIGXBARS = SIGXBARS/NCASES
      SIGXBARO = SIGXBARO/NCASES
      SIGMBARI = SIGMBARI/NCASES
      SIGMBARS = SIGMBARS/NCASES
      SIGMBARO = SIGMBARO/NCASES
      SIGLBARI = SIGLBARI/NCASES
      SIGLBARS = SIGLBARS/NCASES
      SIGLBARO = SIGLBARO/NCASES
	IF(NCASES .GT. 1) THEN
	  SIGXVARI = DSQRT(DABS(SIGXVARI-NCASES*SIGXBARI**2)/(NCASES-1))
	  SIGXVARS = DSQRT(DABS(SIGXVARS-NCASES*SIGXBARS**2)/(NCASES-1))
	  SIGXVARO = DSQRT(DABS(SIGXVARO-NCASES*SIGXBARO**2)/(NCASES-1))
	  SIGMVARI = DSQRT(DABS(SIGMVARI-NCASES*SIGMBARI**2)/(NCASES-1))
	  SIGMVARS = DSQRT(DABS(SIGMVARS-NCASES*SIGMBARS**2)/(NCASES-1))
	  SIGMVARO = DSQRT(DABS(SIGMVARO-NCASES*SIGMBARO**2)/(NCASES-1))
	  SIGLVARI = DSQRT(DABS(SIGLVARI-NCASES*SIGLBARI**2)/(NCASES-1))
	  SIGLVARS = DSQRT(DABS(SIGLVARS-NCASES*SIGLBARS**2)/(NCASES-1))
	  SIGLVARO = DSQRT(DABS(SIGLVARO-NCASES*SIGLBARO**2)/(NCASES-1))
	ELSE
	  SIGXVARI = 0.0D0
	  SIGXVARS = 0.0D0
	  SIGXVARO = 0.0D0
	  SIGMVARI = 0.0D0
	  SIGMVARS = 0.0D0
	  SIGMVARO = 0.0D0
	  SIGLVARI = 0.0D0
	  SIGLVARS = 0.0D0
	  SIGLVARO = 0.0D0
	END IF
C  Do statistics on layer/particle failures
	TEMPERORY = DBLE(NCASES)/DBLE(NBURP)
	IF(NCASES .GT. 1) THEN
	  CONFIDIF = DSQRT(DABS((TEMPERORY*CONFIDIF-PROBIF**2)
     &			  /TEMPERORY/(TEMPERORY-1)))
	  CONFIDSF = DSQRT(DABS((TEMPERORY*CONFIDSF-PROBSF**2)
     &			  /TEMPERORY/(TEMPERORY-1)))
	  CONFIDOF = DSQRT(DABS((TEMPERORY*CONFIDOF-PROBOF**2)
     &			  /TEMPERORY/(TEMPERORY-1)))
	  CONFIDPF = DSQRT(DABS((TEMPERORY*CONFIDPF-PROBPF**2)
     &			  /TEMPERORY/(TEMPERORY-1)))
	ELSE
	  CONFIDIF = 0.0D0
	  CONFIDSF = 0.0D0
	  CONFIDOF = 0.0D0
	  CONFIDPF = 0.0D0
	END IF
	PROBIF = PROBIF/NCASES
	PROBSF = PROBSF/NCASES
	PROBOF = PROBOF/NCASES
	PROBPF = PROBPF/NCASES
	CONFIDIF = CONFIDIF/NBURP
	CONFIDSF = CONFIDSF/NBURP
	CONFIDOF = CONFIDOF/NBURP
	CONFIDPF = CONFIDPF/NBURP
C
	IF(.NOT.PARAMETRIC_STUDY) GO TO 512
C  End of the parametric study loop
	IF(PERTURBATION_ANALYSIS) THEN
	  WRITE(IPRMI, 636) SIGXIPYC, SIGXSIC, SIGXOPYC, 
     &                    SIGMIPYC, SIGMSIC, SIGMOPYC,
     &					SIGLIPYC, SIGLSIC, SIGLOPYC,
     &					SIGFCIPYC,SIGFCSIC,SIGFCOPYC
	ELSE IF(SURFACE_ANALYSIS) THEN
	  WRITE(IPRMI, 629) SIGXIPYC
	  WRITE(IPRMS, 629) SIGMSIC
C
	  IF(SIGXIPYC .GT. SIGXIPYCX) THEN
	    SIGXIPYCX = SIGXIPYC
	    PARASET(1,1) = KERNDIA
	    PARASET(2,1) = IPYCBAF0I
	  END IF
	  IF(SIGXIPYC .LT. SIGXIPYCM) THEN
	    SIGXIPYCM = SIGXIPYC
	    PARASET(1,2) = KERNDIA
	    PARASET(2,2) = IPYCBAF0I
	  END IF
	  IF(SIGMSIC .GT. SIGMSICX) THEN
	    SIGMSICX = SIGMSIC
	    PARASET(1,3) = KERNDIA
	    PARASET(2,3) = IPYCBAF0I
	  END IF
	  IF(SIGMSIC .LT. SIGMSICM) THEN
	    SIGMSICM = SIGMSIC
	    PARASET(1,4) = KERNDIA
	    PARASET(2,4) = IPYCBAF0I
	  END IF
	END IF
C
401	CONTINUE
	WRITE(IPRMI, *)
	IF(SURFACE_ANALYSIS) WRITE(IPRMS, *)
400	CONTINUE
C
	IF(PARAMETRIC_STUDY .AND. SURFACE_ANALYSIS) THEN
	  WRITE(IPRMI, *)
	  WRITE(IPRMI, *) 'Lowest value of maximum IPyC stress: ', 
     &                 SIGXIPYCM
	  WRITE(IPRMI, *) '@  Kernel Diameter: ', PARASET(1,2)
	  WRITE(IPRMI, *) '&  IPyC BAF0: ', PARASET(2,2)
	  WRITE(IPRMI, *)
	  WRITE(IPRMI, *) 'Highest value of maximum IPyC stress: ', 
     &                 SIGXIPYCX
	  WRITE(IPRMI, *) '@  Kernel Diameter: ', PARASET(1,1)
	  WRITE(IPRMI, *) '&  IPyC BAF0: ', PARASET(2,1)
C
	  WRITE(IPRMS, *)
	  WRITE(IPRMS, *) 'Lowest value of minimum SiC stress: ', 
     &                 SIGMSICM
	  WRITE(IPRMS, *) '@  Kernel Diameter: ', PARASET(1,4)
	  WRITE(IPRMS, *) '&  IPyC BAF0: ', PARASET(2,4)
	  WRITE(IPRMS, *)
	  WRITE(IPRMS, *) 'Highest value of minimum SiC stress: ', 
     &                 SIGMSICX
	  WRITE(IPRMS, *) '@  Kernel Diameter: ', PARASET(1,3)
	  WRITE(IPRMS, *) '&  IPyC BAF0: ', PARASET(2,3)
	END IF
C
512	IF((.NOT.PARAMETRIC_STUDY) .AND. HISTOGRAM) THEN
	  WRITE(IHIS,632) 
	  DO 1010 I = 0, NHIS
	    WRITE(IHIS,633) I
	    WRITE(IHIS,634) HISTGRMP(I,1), HISTGRMP(I,2), 
     &                    HISTGRMP(I,3), HISTGRMP(I,4)
	    WRITE(IHIS,634) HISTGRMS(I,1), HISTGRMS(I,2), 
     &                    HISTGRMS(I,3), HISTGRMS(I,4)
	    WRITE(IHIS,634) HISTGRMI(I,1), HISTGRMI(I,2), 
     &                    HISTGRMI(I,3), HISTGRMI(I,4)
	    WRITE(IHIS,634) HISTGRMO(I,1), HISTGRMO(I,2), 
     &                    HISTGRMO(I,3), HISTGRMO(I,4)
	    WRITE(IHIS,*)
1010    CONTINUE
	END IF
C  Simulation all done, send results to user
C 
C  Write the debug file various info. after the calculation
	IF((.NOT.PARAMETRIC_STUDY) .AND. DEBUG) THEN
		WRITE(IDBG, *)
		WRITE(IDBG, *) 'REACTOR PARAMETERS:'
		WRITE(UNIT = IDBG, NML = REACTOR)
		WRITE(IDBG, *)
		WRITE(IDBG, *) 'FUEL PEBBLE PARAMETERS:'
		WRITE(UNIT = IDBG, NML = PEBBLES)
		WRITE(IDBG, *)
		WRITE(IDBG, *) 'FUEL PARTICLE PROPERTIES AND STATISTICS:'
		WRITE(UNIT = IDBG, NML = PARTI)
		WRITE(IDBG, *)
	    WRITE(IDBG, *) 'CONTROL VARIABLES FOR HISTOGRAMS:'
	    WRITE(IDBG, *) 'TIMELIMIT = ', TIMELIMIT
	    WRITE(IDBG, *) 'EOL FLUENCE = ', EOLFLU
	    WRITE(IDBG, *) 'EOL BURNUP = ', EOLBUP
	    WRITE(IDBG, *) 'STRESS UPPER LIMIT = ', SIG_UPPER
	    WRITE(IDBG, *) 'STRESS LOWER LIMIT = ', SIG_LOWER
	    WRITE(IDBG, *) 'HISTOGRAM INTEVALS = ', NHIS
		WRITE(IDBG, *)
		WRITE(IDBG, *) 'End of debug file.'
	END IF
C
C  Comment out the original time statements because they do not work.
C      TIME = 0
C	MACHTIME = TIME
C	WRITE(ITERM,*)
C	IF (TIME .LT. 60.0) THEN
C	  WRITE(ITERM,*) 'Elapsed calculation time: ',TIME,' seconds'
C	ELSE IF ((TIME .GT. 60.0) .AND. (TIME .LT. 3600.0)) THEN
C	  WRITE(ITERM,*) 'Elapsed calculation time: ',INT(TIME/60.0),
C     A                 ' minutes and ',MOD(INT(TIME),60),' seconds'
C	ELSE IF (TIME .GT. 3600.0) THEN
C	  WRITE(ITERM,*) 'Elapsed calculation time: ',INT(TIME/3600.0),
C     A                 ' hours ',INT(MOD(INT(TIME),3600)/60.0),
C     A                 'minutes and ',MOD(MOD(INT(TIME),3600),60),
C     A                 'seconds'
C	END IF
C
C  Print the elapsed CPU time to the terminal
      CPUTIME = MCLOCK()
      WRITE(ITERM,*)
      WRITE(ITERM,*)'Elapsed calculation time (sec):', CPUTIME/1000
      WRITE(ITERM,*)
C  Print the number of particles per second to the terminal
	WRITE(ITERM,*)
	WRITE(ITERM,*) 'Average number of particles per second = ',
     A                FLOAT(NCASES)/(CPUTIME/1000)
	WRITE(ITERM,*)
C
C	WRITE(IOUT,*)
C	IF (TIME .LT. 60.0) THEN
C	  WRITE(IOUT,*) 'Elapsed calculation time: ',TIME,' seconds'
C	ELSE IF ((TIME .GT. 60.0) .AND. (TIME .LT. 3600.0)) THEN
C	  WRITE(IOUT,*) 'Elapsed calculation time: ',INT(TIME/60.0),
C     A                 ' minutes and ',MOD(INT(TIME),60),' seconds'
C	ELSE IF (TIME .GT. 3600.0) THEN
C	  WRITE(IOUT,*) 'Elapsed calcuation time: ',INT(TIME/3600.0),
C     A                 ' hours ',INT(MOD(INT(TIME),3600)/60.0),
C     A                 'minutes and', MOD(MOD(INT(TIME),3600),60),
C     A                 'seconds'
C	END IF
	WRITE(IOUT,*)
	WRITE(IOUT,*) '======================= Statistical Report ',
     &              '======================='
	WRITE(IOUT,*)
	WRITE(IOUT,*) 'Average number of particles per second = ',
     A                FLOAT(NCASES)/MAX(TIME, 1.0 D0)
	WRITE(IOUT,*)
C
      !CALL CURRENTDATE (CDATE)
      !CALL CURRENTTIME (CTIME)
	WRITE(IOUT,638) 'Sample maximum tangential stress in IPyC',
     &				SIGXBARI,'+/-',SIGXVARI
	WRITE(IOUT,638) 'Sample maximum tangential stress in SiC ',
     &				SIGXBARS,'+/-',SIGXVARS
	WRITE(IOUT,638) 'Sample maximum tangential stress in OPyC',
     &				SIGXBARO,'+/-',SIGXVARO
	WRITE(IOUT,638) 'Sample minimum tangential stress in IPyC',
     &				SIGMBARI,'+/-',SIGMVARI
	WRITE(IOUT,638) 'Sample minimum tangential stress in SiC ',
     &				SIGMBARS,'+/-',SIGMVARS
	WRITE(IOUT,638) 'Sample minimum tangential stress in OPyC',
     &				SIGMBARO,'+/-',SIGMVARO
	WRITE(IOUT,638) 'Sample EOL tangential stress in IPyC    ',
     &				SIGLBARI,'+/-',SIGLVARI
	WRITE(IOUT,638) 'Sample EOL tangential stress in SiC     ',
     &				SIGLBARS,'+/-',SIGLVARS
	WRITE(IOUT,638) 'Sample EOL tangential stress in OPyC    ',
     &				SIGLBARO,'+/-',SIGLVARO
	WRITE(IOUT,*)
      WRITE(IOUT,630) 'Number of failed particles:           ',PARFAIL
	WRITE(IOUT,630) 'Number of particle with SiC crack:    ',SICFAIL
	WRITE(IOUT,630) 'Number of particle with IPyC crack:   ',IPYCFAIL
	WRITE(IOUT,630) 'Number of particle with OPyC crack:   ',OPYCFAIL
	WRITE(IOUT,638) 'Particle failure probability:           ',PROBPF
     &				,'+/-',CONFIDPF
	WRITE(IOUT,638) 'SiC layer failure probability:          ',PROBSF
     &				,'+/-',CONFIDSF
	WRITE(IOUT,638) 'IPyC layer failure probability:         ',PROBIF
     &				,'+/-',CONFIDIF
	WRITE(IOUT,638) 'OPyC layer failure probability:         ',PROBOF
     &				,'+/-',CONFIDOF
      WRITE(IOUT,630) 'Fracture induced by IPyC cracking:    ',
     &                FAILTYPE(0)
	WRITE(IOUT,630) 'Overpressure rupture:                 ',
     &                FAILTYPE(1)
	WRITE(IOUT,630) 'Failure due to the Amoeba Effect:     ',
     &                FAILTYPE(2)
	WRITE(IOUT,630) 'Dumped particles exceeding EOL burnup ',NDUMPED
	WRITE(IOUT,*)
	WRITE(IOUT,*) '=================== End of Statistical Report ',
     &              '===================='
	WRITE(IOUT,*)
      WRITE(IOUT,622) 'The calculation ends at: '
	WRITE(IOUT,617) CDATE,CTIME
C
      CLOSE(UNIT = IOUT, STATUS = 'KEEP')
	CLOSE(UNIT = IERR, STATUS = 'KEEP')
	IF (NOMINAL) THEN
	   CLOSE(UNIT = IEVOLR, STATUS = 'KEEP')
	   CLOSE(UNIT = IEVOLT, STATUS = 'KEEP')
	ENDIF
      WRITE(ITERM,*)
C
      CLOSE(UNIT = IDAT1, STATUS = 'KEEP')  ! Normal end of program here
     	FOPEN1 = .FALSE.
      IF(FOPEN2.EQ..TRUE.) THEN
        CLOSE(UNIT = IDAT2, STATUS = 'KEEP')
        FOPEN2 = .FALSE.
      ENDIF
      IF(FOPEN3.EQ..TRUE.) THEN
        CLOSE(UNIT = IDAT3, STATUS = 'KEEP')
        FOPEN3 = .FALSE.
      ENDIF
      IF(FOPEN4.EQ..TRUE.) THEN
        CLOSE(UNIT = IDAT4, STATUS = 'KEEP')
        FOPEN4 = .FALSE.
      ENDIF
	WRITE(ITERM,622) 'The calculation ends at: '
	WRITE(ITERM,617) CDATE, CTIME
	WRITE(ITERM,*) ' ** Please view the output files for details. **'
	WRITE(ITERM,*)
      WRITE(ITERM,*) ' One input set detected:  ',
     A               'Normal end of code'
      WRITE(ITERM,*)
      STOP
C
C  Formatting
C
601   FORMAT('  ##### ### #   #   ###   ##     ##   #####'/,
     A	   '    #    #  ## ##  #     #  #   #  #    #  '/,
     A	   '    #    #  # # # #     #    # #    #   #  '/,
     A	   '    #    #  #   # #     #    # ######   #  '/,
     A	   '    #    #  #   # #     #    # #    #   #  '/,
     A	   '    #    #  #   #  #     #  #  #    #   #  '/,
     A	   '    #   ### #   #   ###   ##   #    #   #  '/,/,
     A	   'Monte Carlo calculation of mechanical failures ',
     A       'in a statistical '/,
     B       'sample of TRISO-coated fuel particles for the High ',
     C       'Temperature Gas-cooled '/,
     D       'Reactors (HTGRs).'//)
602   FORMAT(A256)
603   FORMAT(/' Start of Sampling:',I9,' Cases'/
     A/'            Elapsed      SiC        OPyC        SiC       ',
     B   'OPyC',
     C/'   Cases     Time      Failures   Failures   Failures   ',
     D  'Failures',
     E/' Completed    sec       Irrad      Irrad      Anneal     ',
     F  'Anneal'/)
604   FORMAT(/' Start of Sampling:',I9,' Cases'/
     &/'             Elapsed    Particle     SiC        IPyC       ',
     &    'OPyC',
     &/'   Cases      Time      Failures   Failures    Failures   ',
     &  'Failures',
     &/' Completed    (sec)      (Irrad)   (Irrad)     (Irrad)    ',
     &   '(Irrad)'/)
605   FORMAT(I10,F9.0,I10,I11,I12,I12)
606   FORMAT(' ',' Debug output for case ',I7,' at step ',I3)
607   FORMAT(/I10,' Cases Sampled:'//
     A  '   SiC Failure probability in irradiation:',
     B                       1PE10.2,' +/- ',1PE10.2,//
     C  '  OPyC Failure probability in irradiation:',
     D                       1PE10.2,' +/- ',1PE10.2,//
     E  '     SiC Failure probability in annealing:',
     F                       1PE10.2,' +/- ',1PE10.2,//
     G  '    OPyC Failure probability in annealing:',
     H                       1PE10.2,' +/- ',1PE10.2,//)
608   FORMAT(/I10,' Cases Sampled:'//
     A  '   SiC Failure probability in irradiation:',
     B                       1PE10.2,' +/- ',1PE10.2,//
     C  '  OPyC Failure probability in irradiation:',
     D                       1PE10.2,' +/- ',1PE10.2,//)
609   FORMAT(' Summary of Statistics:',
     A   /'  Layer    Conditions        Mean Stress     Std. Dev.',
     B   /'                                [MPa]         [MPa]'   )
611   FORMAT(/A6,4X,A11,4X,F12.2,2X,F12.2)
612   FORMAT('   Step ',I2,'   SiC Failures = ',I9)
C 212  FORMAT(A6, A1, A3, A1, A3, A1, A4, A1, A4)
C 213  FORMAT(A3, A1, A5, A1, A6, A1, A5, A1, A6)
613   FORMAT(A6, A1, A1, A3, A1, A1, A3, A1, A1, A4, A1, A1, A4)
614   FORMAT(A3, A1, A1, A5, A1, A1, A6, A1, A1, A5, A1, A1, A6)
615   FORMAT(A1,1PE10.4,A1,1PE9.3,A1,1PE9.3,A1,1PE9.3,A1,1PE9.3)
616   FORMAT(   1PE10.4,A1,1PE9.3,A1,1PE9.3,A1,1PE9.3,A1,1PE9.3)
617   FORMAT(A8, A14,/)
618   FORMAT(74('_'),/)
619   FORMAT(I3,4X,F7.2,4X,F10.4,4X,F10.4)
621   FORMAT(2X,F8.3,\)
622   FORMAT(A25,\)
623   FORMAT(11X,F5.3,6(4X,E12.5),5(5X,F9.2))
624   FORMAT(32X,4(2X,E14.7))
625   FORMAT(A27,\)
626   FORMAT(3X,F9.3,11X,F9.3,9X,E12.5,5X,F9.4,8X,F9.3)
627   FORMAT(7(3X,F8.3))
628   FORMAT(5(4X,E12.5))
629   FORMAT(2X,F9.3,\)
630   FORMAT(A39,I11)
631   FORMAT(A39,E11.4)
632   FORMAT('  I',' | ',9('*'),' Particle ',9('*'),' |',
     &            ' | ',8('*'),' SiC Layer ',8('*'),' |',
     &           ' | ',7('*'),' IPyC Layer ',8('*'),' |',
     &           ' | ',7('*'),' OPyC Layer ',8('*'),' |',/,
     &       8X,'Time',3X,'Stress',4X,'Flu',4X,'Burp',
     &       3(4X,'Time',3X,'Stress',4X,'Flu',4X,'Burp'))
633   FORMAT(I3,\)
634   FORMAT(4(I8),\)
635	FORMAT('  Fluence (10^21nvt) IPyC Sigma0    IPyC strength',\,
     &	   '     SiC Sigma0     SiC strength    OPyC Sigma0',\,
     &	   '    OPyC strength      SiC KIC        IPyC KI',\,
     &	   '      OPyC KI    SiC KI from IPyC  SiC KI from OPyC')
636	FORMAT(12(4X,F8.2),\)
637   FORMAT('Operation Time (d)   Fluence (10E21nvt)   ',\,
     &	   'QPPP (W/m^3)    Burnup (FIMA)    Coolant T (C)',/)
638   FORMAT(A41,2X,E11.4,1X,A3,1X,E11.4)
639	FORMAT('                Max IPyC    Max SiC     ',\,
     &				  'Max OPyC    Min IPyC    Min SiC     Min OPyC'
     &	  ,\,'    EOL IPyC    EOL SiC     EOL OPyC',\
     &	  ,'   IPyC Strgth  SiC Strgth  OPyC Strgth',/)
640	FORMAT(A11,\)
641   FORMAT(2X,E10.3,\)
642   FORMAT(2X,F10.3,\)
643	FORMAT('Case No.  Failure Type    R1    R2    R3    R4    R5',\,
     &	   '    BuffD   IPyCD   OPyCD   IPyCBAF0   OPyCBAF0   ',\,
     &	   'IPyC SigF   OPyC SigF   SiC SigF   SiC KIc   ',\,
     &	   'IPyC Stress   OPyC Stress   SiC Stress   Irr. Time   ',\,
     &	   'Fluence     Burnup   IPyC T   Inner P   Channel   ',\,
     &	   'Block   Cycle     FailurePath',//
     &	   '                         (um)  (um)  (um)  (um)  (um)',\,
     &	   '   (g/cc)  (g/cc)  (g/cc)                        ',\,
     &	   '  (MPa)       (MPa)       (MPa)  (MPa.um^1/2)',\,
     &	   '   (MPa)         (MPa)        (MPa)       (days)     ',\,
     &	   '(10^21nvt)  (FIMA)     (C)     (MPa)              ',\,
     &	   '                ',/)
644	FORMAT('Case No.  Failure Type    R1    R2    R3    R4    R5',\,
     &	   '    BuffD   IPyCD   OPyCD   IPyCBAF0   OPyCBAF0   ',\,
     &	   'IPyC SigF   OPyC SigF   SiC SigF   SiC KIc   ',\,
     &	   'IPyC Stress   OPyC Stress   SiC Stress   Irr. Time   ',\,
     &	   'Fluence     Burnup   IPyC T   Inner P     FailurePath',//
     &	   '                         (um)  (um)  (um)  (um)  (um)',\,
     &	   '   (g/cc)  (g/cc)  (g/cc)                        ',\,
     &	   '  (MPa)       (MPa)       (MPa)  (MPa.um^1/2)',\,
     &	   '   (MPa)         (MPa)        (MPa)       (days)     ',\,
     &	   '(10^21nvt)  (FIMA)     (C)     (MPa)    ',/)
645   FORMAT(I7,7X,I2,8X,5(F5.1,1X),2X,3(F5.3,3X),2(1X,F6.4,4X),\,
     &	   2(2X,F6.1,4X),1X,F6.1,5X,F6.1,3X,2(3X,F6.1,5X),2X,F6.1,\,
     &	   7X,F6.1,4X,1X,F5.3,6X,F6.4,3X,F6.1,3X,F6.2,7X,I1,8X,\,
     &	   I2,6X,I2,2X,A42/)
646   FORMAT(I7,7X,I2,8X,5(F5.1,1X),2X,3(F5.3,3X),2(1X,F6.4,4X),\,
     &	   2(2X,F6.1,4X),1X,F6.1,5X,F6.1,3X,2(3X,F6.1,5X),2X,F6.1,\,
     &	   7X,F6.1,4X,1X,F5.3,6X,F6.4,3X,F6.1,3X,F6.2,2X,A42/)
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C  Subroutine CORE(COREMODEL, CHANNELS, BLOCKS, BLOCKMAP)              *
C                                                                      *
C    This subroutine reads in channel and block information for input  *
C  files channels.dat, blocks.dat.                                     *
C    Caution: The dimension of CHANNELS and BLOCKS are specified in the*
C             main program, but their stored data are read in this     *
C             subroutine from channels.dat and blocks.dat, therefore   *
C             they must be consistent. If the input files are changed, *
C             dimension should be changed accordingly.                 *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual argument description                                         *
C    COREMODEL      I: indicates whether new or old core model is      *
C					 employed. (currently 1 for new and 2 for old)   *
C    CHANNELS(NAXIAL,-1:NCHANNEL) --                                   *
C                   D: defines boundary lines of every channel         *
C                      format: (z-position,r-position)                 *
C                      (coordinates in cm)                             *
C    BLOCKS(NLAYER,5) --                                               *
C                   D: defines axial positions (cm) and neutron fluxes *
C                      (1/cm^2.sec) of blocks of every channel         *
C    BLOCKMAP(2,NCHANNEL) --                                           *
C                   D: specifies number of layers in every channel     *
C                               First row: num. of layers in each chn  *
C                               Second row: cumulative layers          *
C                                                                      *
C  Local variable description                                          *
C    I              I: loop counter                                    *
C                                                                      *
C  Parameters                                                          *
C    NCHANNEL1      I: number of channels according to new VSOP model  *
C    NAXIAL1        I: number of axial divisions from new VSOP model   *
C    NLAYER1        I: total number of layers from new VSOP model      *
C    NCHANNEL2      I: number of channels according to old VSOP model  *
C    NAXIAL2        I: number of axial divisions from old VSOP model   *
C    NLAYER2        I: total number of layers from old VSOP model      *
C    IDAT           I: lun of input file                               *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE CORE(COREMODEL, CHANNELS, BLOCKS, BLOCKMAP)
C
      DOUBLE PRECISION CHANNELS, BLOCKS
	INTEGER COREMODEL, BLOCKMAP, I
C
      PARAMETER ( NCHANNEL1 = 5, NAXIAL1 = 21, NLAYER1 = 93)
      PARAMETER ( NCHANNEL2 = 5, NAXIAL2 = 15, NLAYER2 = 57)
      PARAMETER ( IDAT = 16, IERR = 17)
C
	DIMENSION CHANNELS(1:,-1:), BLOCKS(1:,1:)
      DIMENSION BLOCKMAP(1:,1:)
C
	IF (COREMODEL .EQ. 1) THEN
C  Initialize BLOCKMAP based on new VSOP core model
        BLOCKMAP(1,1) = 20
        BLOCKMAP(1,2) = 18
        BLOCKMAP(1,3) = 18
        BLOCKMAP(1,4) = 18
        BLOCKMAP(1,5) = 19
        BLOCKMAP(2,1) = 20
        BLOCKMAP(2,2) = 38
        BLOCKMAP(2,3) = 56
        BLOCKMAP(2,4) = 74
        BLOCKMAP(2,5) = 93
C  Open channels.dat
        CALL OPENFILE('channels','.dat','OLD',IDAT)
C  First omit two-line headings in channels.dat, dependent on input file.
        READ(UNIT=IDAT,FMT=*)
        READ(UNIT=IDAT,FMT=*)
C  Read in data
        DO 2205 I = 1, NAXIAL1
	    READ(UNIT=IDAT,FMT=*) CHANNELS(I,-1),CHANNELS(I,0),
     &             CHANNELS(I,1),CHANNELS(I,2),CHANNELS(I,3),
     &			 CHANNELS(I,4),CHANNELS(I,5)
 2205   CONTINUE
 	  CLOSE (UNIT = IDAT, STATUS = 'KEEP')
C  Open blocks.dat
        CALL OPENFILE('blocks  ','.dat','OLD',IDAT)
C  Omit two-line heading of blocks.dat, dependent on input file.
        READ(UNIT=IDAT,FMT=*)
        READ(UNIT=IDAT,FMT=*)
C  Read in data
	  DO 2210 I = 1, NLAYER1
	    READ(UNIT=IDAT,FMT=*) BLOCKS(I,1), BLOCKS(I,2), BLOCKS(I,3),
     &                          BLOCKS(I,4), BLOCKS(I,5)
 2210   CONTINUE
        CLOSE (UNIT = IDAT, STATUS = 'KEEP')
	ELSE IF ( COREMODEL .EQ. 2) THEN
C  Initialize BLOCKMAP based on old VSOP core model
        BLOCKMAP(1,1) = 9
        BLOCKMAP(1,2) = 10
        BLOCKMAP(1,3) = 11
        BLOCKMAP(1,4) = 12
        BLOCKMAP(1,5) = 15
        BLOCKMAP(2,1) = 9
        BLOCKMAP(2,2) = 19
        BLOCKMAP(2,3) = 30
        BLOCKMAP(2,4) = 42
        BLOCKMAP(2,5) = 57
C  Open channels.dat
        CALL OPENFILE('channels','.dat','OLD',IDAT)
C  First omit two-line headings in channels.dat, dependent on input file.
        READ(UNIT=IDAT,FMT=*)
        READ(UNIT=IDAT,FMT=*)
C  Read in data
        DO 2215 I = 1, NAXIAL2
	    READ(UNIT=IDAT,FMT=*) CHANNELS(I,-1),CHANNELS(I,0),
     &             CHANNELS(I,1),CHANNELS(I,2),CHANNELS(I,3),
     &			 CHANNELS(I,4),CHANNELS(I,5)
 2215   CONTINUE
 	  CLOSE (UNIT = IDAT, STATUS = 'KEEP')
C  Open blocks.dat
        CALL OPENFILE('blocks  ','.dat','OLD',IDAT)
C  Omit two-line heading of blocks.dat, dependent on input file.
        READ(UNIT=IDAT,FMT=*)
        READ(UNIT=IDAT,FMT=*)
C  Read in data
	  DO 2220 I = 1, NLAYER2
	    READ(UNIT=IDAT,FMT=*) BLOCKS(I,1), BLOCKS(I,2), BLOCKS(I,3),
     &                          BLOCKS(I,4), BLOCKS(I,5)
 2220   CONTINUE
        CLOSE (UNIT = IDAT, STATUS = 'KEEP')
	ELSE
	  CALL ERR_HANDLER(TRIM('COREMODEL not specified when calling CORE'),
     &			  45,2,2,IERR)
	END IF
	RETURN
	END
C                                                                      *
C***********************************************************************
C
C
C
!C***********************************************************************
!C                                                                      *
!C  Subroutine CURRENTTIME                                              *
!C                                                                      *
!C     Get the current time hour:minutes:seconds                        *
!C                                                                      *
!C  It is assumed that the common block TYPE has been initialized by    *
!C  calling subroutine MACHINE.                                         *
!C                                                                      *
!C  Notation                                                            *
!C    D               : REAL*8 (Double Precision)                       *
!C    I               : INTEGER                                         *
!C    L               : LOGICAL                                         *
!C    C               : CHARACTER*n                                     *
!C                                                                      *
!C  Actual argument description                                         *
!C    CTIME          C: Time format HH:MM:SS                            *
!C                                                                      *
!C  Common blocks:                                                      *
!C  /MTYPE/                                                             *
!C                                                                      *
!C  Description of variables                                            *
!C    CTIME2         C: Time format HH:MM:SS.SS                         *
!C    COMP           C: Three character description of compiler type    *
!C    C1             I: Numerical constant, 2 (2^31 - 1) = 4 294 967 294*
!C    C2             I: Numerical constant, 86400, number of seconds    *
!C                      in a day                                        *
!C    IERR           I: Logical Unit Number of error file 'ERROR.MSG'   *
!C    MACH           C: Three character description of computer type    *
!C    NHOUR          I: Current time in hours                           *
!C    NMIN           I: Current time in minutes                         *
!C    NSEC           I: Current time in seconds                         *
!C    NTIME          I: Reference time in seconds                       *
!C    PERIOD         D: Processed reference time                        *
!C                                                                      *
!C  None-standard features:                                             *
!C    LONG            : 32 bit absolute memory addressing function for  *
!C                      Macintosh computers.  The memory address 524    *
!C                      (decimal) has the time in seconds since midnight*
!C                       January 1, 1904. LONG is an INTEGER function.  *
!C                                                                      *
!C    TIME()          : On the DEC workstations with DEC Ultrix Fortran,*
!C                      the call to TIME returns the time in seconds    *
!C                      since midnight January 1, 1970.  TIME is an     *
!C                      INTEGER function.                               *
!C                                                                      *
!C    TIME()          : On the IBM, TIME returns the time in the format *
!C                      Hour:Minute:Second:Hundredth of a Second. Time  *
!C                      is a CHARACTER function.                        *
!C                                                                      *
!C  Functions and subroutines called:                                   *
!C    ERR_HANDLER     : Processes an error (diagnostic) message.        *
!C                                                                      *
!C  Intrinsic functions called:                                         *
!C    CHAR                                                              *
!C    DBLE                                                              *
!C    ICHAR                                                             *
!C    IDINT                                                             *
!C    MOD                                                               *
!C                                                                      *
!C***********************************************************************
!C
!C***********************************************************************
!C                                                                      *
!      SUBROUTINE CURRENTTIME (CTIME)
!C
!      DOUBLE PRECISION PERIOD, C1, C2
!      INTEGER NTIME, NHOUR, NMIN, NSEC, IERR
!C
!C  On the DEC workstations with DEC Ultrix Fortran, the call to TIME
!C  returns the time in seconds since midnight January 1, 1970
!C  Greenwich Meridian Time (GMT).  The function TIME must be declared
!C  INTEGER.
!C
!C#######################################################################
!CDEC	INTEGER TIME
!C#######################################################################
!      CHARACTER*8 CTIME
!      CHARACTER*11 CTIME2
!      CHARACTER*3 MACH, COMP
!C
!      COMMON /MTYPE/MACH, COMP
!C
!      PARAMETER (IERR = 12)
!C
!C  Set CTIME2 to blanks to avoid compiler warnings
!      DATA CTIME2/'           '/
!C
!C  Set to zero to avoid compiler warnings
!      DATA NTIME/0/
!C
!C  Define numerical constants
!C  2 (2^31 - 1)
!      DATA C1/4 294 967 294.0 D+00/
!      DATA C2/86400.0 D+00/     ! Seconds in day
!C
!C  Machine dependent calculations
!C
!      IF ((MACH .EQ. 'MAC') .OR. (MACH .EQ. 'DEC')) THEN
!C
!C  Generic to all Macintosh compilers
!C
!C  The memory address 524 (decimal) has the time in seconds since
!C  midnight January 1, 1904.
!C
!C#######################################################################
!CMAC      IF (MACH .EQ. 'MAC') NTIME = LONG(524)
!C#######################################################################
!C
!C  On the DEC workstations with DEC Ultrix Fortran, the call to TIME
!C  returns the time in seconds since midnight January 1, 1970.
!C
!C#######################################################################
!CDEC      IF (MACH .EQ. 'DEC') NTIME = TIME()
!C#######################################################################
!C
!C  The number of seconds returned are in the range - 2^31 -> 2^31 - 1
!C  Negative values are used to extend the range of the system clocks.
!C  If a negative value for time is returned 2*(2^31 - 1) must be added
!C  to get the correct time.  This would cause INTEGER*4 overflow.
!C  To prevent this, the time is converted to double precision, then
!C  scaled to the number of days.
!C  
!        IF (NTIME .GE. 0) THEN
!          PERIOD = DBLE(NTIME)
!        ELSE
!          PERIOD = C1 + DBLE(NTIME)
!        END IF
!C
!C  Remaining fraction of day in seconds
!C
!        NTIME = IDINT(C2*(PERIOD/C2 - DBLE(IDINT(PERIOD/C2))))
!        NHOUR = NTIME/3600
!        NMIN = MOD(NTIME,3600)/60
!        NSEC = MOD(NTIME,3600) - 60*NMIN
!      ELSE IF (MACH .EQ. 'IBM') THEN
!C
!C  Only for IBM compilers
!C
!C  IBM routine TIME returns the time in the format
!C  'Hour:Minute:Second '
!        CALL TIME(CTIME2)
!        NHOUR = 10*(ICHAR(CTIME2(1:1)) - 48) + ICHAR(CTIME2(2:2)) - 48
!        NMIN = 10*(ICHAR(CTIME2(4:4)) - 48) + ICHAR(CTIME2(5:5)) -48
!        NSEC = 10*(ICHAR(CTIME2(7:7)) - 48) + ICHAR(CTIME2(8:8)) - 48
!      ELSE
!        CALL ERR_HANDLER('SUBROUTINE CURRENTTIME: unsupported 
!     &					machine type', 42, 2, 2, IERR)
!      END IF
!C
!C  Express time in HH:MM:SS form
!C
!C  ASCII number for 0 is 48, ... ASCII number for 9 is 57
!C
!C  Tens of hours
!      CTIME(1:1) = CHAR(48 + NHOUR/10) 
!C  Single hours
!      CTIME(2:2) = CHAR(48 + MOD(NHOUR,10))
!C  Break
!      CTIME(3:3) = ':' 
!C  Tens of minutes
!      CTIME(4:4) = CHAR(48 + NMIN/10)
!C  Single minutes
!      CTIME(5:5) = CHAR(48 + MOD(NMIN,10)) 
!C  Break
!      CTIME(6:6) = ':'  
!C  Tens of seconds
!      CTIME(7:7) = CHAR(48 + NSEC/10)
!C  Single seconds
!      CTIME(8:8) = CHAR(48 + MOD(NSEC,10)) 
!      RETURN
!      END
!C                                                                      *
!C***********************************************************************
C
C
C
C***********************************************************************
C  Subroutine FAILURE(SIGR, SIGT, FAIL, FAILTYPE, DPD, N, OPERTIME,    *
C                     DT, MD)                                          *
C                                                                      *
C    This subroutine evaluates fuel particle failures due to both      *
C  mechanical and chemical effects.                                    *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual argument description                                         *
C    SIGR(1:NDIV) (MPa) --                                             *
C                   D: Radial stress distribution in particles         *
C    SIGT(1:NDIV) (MPa) --                                             *
C                   D: Tangential stress distribution in particles     *
C    FAIL           L: Flag indicating if particle fails or not        *
C    FAILTYPE(0:2)  I: array counting the number of each type of       *
C                      failure                                         *
C                      FAILTYPE(0): fracture induced by IPyC cracking  *
C                      FAILTYPE(1): overpressure rupture               *
C                      FAILTYPE(2): amoeba effect                      *
C  DPD (micrometers)D: Distance of Pd crack penetration                *  
C  N                I: Current Case (Particle) Number                  *
C  OPERTIME (sec)   D: Total virtual reactor operation time            *
C  DT (sec)         D: Time increment in the reactor core              *
C  MD (micrometers) D: kernel migration distance                       *
C  Variables in COMMON blocks                                          *
C  /PAR_R/ --                                                          *
C    R1 - R5        D: Fuel particle geometry                          *
C  /PAR_DIV/ --                                                        *
C    NDIVI          I: Number of divisions in IPyC (For stress distr.) *
C    NDIVO          I: Number of divisions in OPyC (For stress distr.) *
C    NDIVS          I: Number of divisions in SiC  (For stress distr.) *
C  /PAR_F/ --                                                          *
C    SIGFIPYC (MPa) D: Fracture strength of IPyC                       *
C    SIGFOPYC (MPa) D: Fracture strength of OPyC                       *
C    SIGFSIC (MPa)  D: Fracture strength of SiC                        *
C    KICSIC (MPa.um^1/2) --                                            *
C                   D: Critical stress intensity factor of SiC         *
C  /NFAIL/ --                                                          *
C    IPYCFAIL       I: Number of particles with IPyC failure           *
C    OPYCFAIL       I: Number of particles with OPyC failure           *
C    SICFAIL        I: Number of particles with SiC failure            *
C    PARFAIL        I: Number of failed particles                      *
C    IPYCFAILED     I: Indicates if IPyC failed or not                 *
C                      0: means it's fine                              *
C                      1: means it just failed                         *
C                      2: means it already failed                      *
C    OPYCFAILED     I: Indicates if OPyC failed or not                 *
C    SICFAILED      I: Indicates if SiC failed or not                  *
C    PARFAILED      I: Indicates if particle failed or not             *
C    FMODE          I: The failure mode of current failed particle     *
C                      0:  fracture induced by IPyC cracking           *
C                      1:  overpressure rupture                        *
C                      2:  amoeba effect                               *
C    MCODE          C: Indicates the type of analysis to perform       *
C					 'ISO3': full three-layer analysis               *
C					 'IS2' : IPyC/SiC two-layer analysis             *
C					 'SO2' : SiC/OPyC two-layer analysis             *
C					 'S1'  : SiC single-layer analysis               *
C    PSTATE         C: Indicates the current particle state            *
C    FAILUREPATH    C: Records the path a failed particle undertake    *
C    NCHAR          I: Counts number of characters in FAILUREPATH      *
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
C                                                                      *
C  Local argument description                                          *
C    SIGTIPYC (MPa) D: Average tangential stress in IPyC               *
C    SIGTSIC (MPa)  D: Average tangential stress in SiC                *
C    SIGTOPYC (MPa) D: Average tangential stress in OPyC               *
C                                                                      *
C  Parameters                                                          *
C    NDIV           I: Number of sampled points in structural layers   *
C                      for stress analyses                             *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE FAILURE(SIGR, SIGT, FAIL, FAILTYPE, DPD, N, 
     &         	OPERTIME, DT, MD)
C
      DOUBLE PRECISION SIGR, SIGT
	DOUBLE PRECISION SIGTIPYC, SIGTSIC, SIGTOPYC
	DOUBLE PRECISION R1, R2, R3, R4, R5
	DOUBLE PRECISION SIGFIPYC, SIGFOPYC, SIGFSIC, KICSIC
	DOUBLE PRECISION IPYCALPHA, IPYCCNU, IPYCE, IPYCNU, IPYCREEP,
     &                 OPYCALPHA, OPYCCNU, OPYCE, OPYCNU, OPYCREEP,
     &                 SICALPHA, SICE, SICNU, IPYCD, OPYCD
      DOUBLE PRECISION ISWR, ISWT, OSWR, OSWT
	DOUBLE PRECISION EPIRCP, EPITCP, URCP, KIIPYC, KIOPYC, KI1, KI2,
     &				 SHEARIPYC, SHEAROPYC, DF
	DOUBLE PRECISION R, OPERTIME, DT
	DOUBLE PRECISION DPD, MD
      LOGICAL   FAIL
      INTEGER*4 FAILTYPE(0:2)
	INTEGER*4 IPYCFAIL, SICFAIL, OPYCFAIL, PARFAIL, NCHAR, FMODE
	INTEGER   NDIVI, NDIVS, NDIVO
	INTEGER   I, J, K
	INTEGER   IPYCFAILED, SICFAILED, OPYCFAILED, PARFAILED
	INTEGER   IKISIC
	CHARACTER*3 PSTATE
	CHARACTER*4 MCODE
	CHARACTER*42 FAILUREPATH
C
      PARAMETER (  NDIV = 30)
      PARAMETER (  NDEG = 3)
      PARAMETER (  PIE = 3.1415926535897932385 D0)
C
      DIMENSION SIGR(1:NDIV), SIGT(1:NDIV)
      DIMENSION ISWR(0:NDEG), ISWT(0:NDEG), OSWR(0:NDEG), OSWT(0:NDEG)
	DIMENSION EPIRCP(1:NDIV), EPITCP(1:NDIV), URCP(1:NDIV)
C
      COMMON /PAR_R/ R1, R2, R3, R4, R5   !particle geometry
	COMMON /PAR_DIV/ NDIVI, NDIVS, NDIVO
      COMMON /PAR_F/ SIGFIPYC, SIGFOPYC, SIGFSIC, KICSIC
      COMMON /NFAIL/ IPYCFAIL, SICFAIL, OPYCFAIL, PARFAIL,
     &               IPYCFAILED, SICFAILED, OPYCFAILED, PARFAILED,
     &			   NCHAR, FMODE, PSTATE, MCODE, FAILUREPATH
      COMMON /PAR_M/ IPYCALPHA, IPYCCNU, IPYCE, IPYCNU, IPYCREEP,
     &               OPYCALPHA, OPYCCNU, OPYCE, OPYCNU, OPYCREEP,
     &               SICALPHA, SICE, SICNU, IPYCD, OPYCD,
     &               ISWR, ISWT, OSWR, OSWT
	COMMON /CRACKED_PYC/ EPIRCP, EPITCP, URCP, KIIPYC, KIOPYC,
     &					 KI1, KI2, SHEARIPYC, SHEAROPYC, DF
	COMMON /OUTPUT/ IKISIC
C
C  Calculate average tangential stresses in layers
      SIGTIPYC = 0.0D0
	SIGTSIC = 0.0D0
	SIGTOPYC = 0.0D0
      DO 2310 I = 1, NDIVI+2
	  SIGTIPYC = SIGTIPYC + SIGT(I)
2310  CONTINUE
      SIGTIPYC = SIGTIPYC/FLOAT(NDIVI+2)
	DO 2320 I = 1, NDIVS+2
	  SIGTSIC = SIGTSIC + SIGT(NDIVI+2+I)
2320  CONTINUE
      SIGTSIC = SIGTSIC/FLOAT(NDIVS+2)
	DO 2330 I = 1, NDIVO+2
	  SIGTOPYC = SIGTOPYC + SIGT(NDIVI+NDIVS+4+I)
2330  CONTINUE
      SIGTOPYC = SIGTOPYC/FLOAT(NDIVO+2)
C
	IF (MD .GT. R3) THEN
	  SICFAIL = SICFAIL + 1
	  SICFAILED = 1
	  IPYCFAIL = IPYCFAIL + 1
	  IPYCFAILED = 1
	  OPYCFAIL = OPYCFAIL + 1
	  OPYCFAILED = 1
	  PARFAIL = PARFAIL + 1
	  PARFAILED = 1
	  PSTATE = 'FFF'
	  FAIL = .TRUE.
	  FAILTYPE(2) = FAILTYPE(2) +1
	  FMODE = 2
	ELSE IF(PSTATE .EQ. 'ISO') THEN
	  IF(SIGTSIC .GT. SIGFSIC) THEN
	    SICFAIL = SICFAIL + 1
	    SICFAILED = 1
	    IPYCFAIL = IPYCFAIL + 1
	    IPYCFAILED = 1
	    OPYCFAIL = OPYCFAIL + 1
	    OPYCFAILED = 1
	    PARFAIL = PARFAIL + 1
	    PARFAILED = 1
	    PSTATE = 'FFF'
	    FAILUREPATH(NCHAR+1:NCHAR+10) = 'ISO -> FFF'
	    NCHAR = NCHAR + 10
	    FAIL = .TRUE.
	    FAILTYPE(1) = FAILTYPE(1) + 1
	    FMODE = 1
	  ELSE IF(SIGTIPYC .GT. SIGFIPYC) THEN
	    IPYCFAIL = IPYCFAIL + 1
	    IPYCFAILED = 1
	    PSTATE = 'CSO'
	    FAILUREPATH(NCHAR+1:NCHAR+10) = 'ISO -> CSO'
	    NCHAR = NCHAR + 10
	    KI1 = 0.0D0
	    IF(KI1. GT. KICSIC) THEN
	      SICFAIL = SICFAIL + 1
	      SICFAILED = 1
	      IF(SIGTSIC .GT. 0.0D0) THEN
	        OPYCFAIL = OPYCFAIL + 1
	        OPYCFAILED = 1
	        PARFAIL = PARFAIL + 1
	        PARFAILED = 1
	        PSTATE = 'FFF'
	        FAILUREPATH(NCHAR+1:NCHAR+14) = ' -> FFO -> FFF'
	        NCHAR = NCHAR + 14
	        FAIL = .TRUE.
	        FAILTYPE(0) = FAILTYPE(0) + 1
	        FMODE = 0
	      ELSE
	        PSTATE = 'FCO'
	        FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> FCO'
	        NCHAR = NCHAR + 7
	        MCODE = 'SO2'
	      END IF
	    ELSE
	      PSTATE = 'FSO'
	      FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> FSO'
	      NCHAR = NCHAR + 7
	      MCODE = 'SO2'
C         Record the elastic response if cracked layer is fully relaxed
		  DO 2350 J = 1, NDIVI + 2
			EPIRCP(J) = (-(SIGR(J)-SIGR(1)) + IPYCNU*SIGT(J))/IPYCE
			EPITCP(J) = (-SIGT(J) + IPYCNU*(SIGR(J)-SIGR(1)))/IPYCE
			URCP(J) = 0.0D0
			DO 2351 K = J, NDIVI + 1
			  URCP(J) = URCP(J) - EPIRCP(K)*(R3-R2)/FLOAT(NDIVI+1)
2351			CONTINUE
2350		  CONTINUE
	    END IF
	  ELSE IF(SIGTOPYC .GT. SIGFOPYC) THEN
	    OPYCFAIL = OPYCFAIL + 1
	    OPYCFAILED = 1
	    PSTATE = 'ISC'
	    FAILUREPATH(NCHAR+1:NCHAR+10) = 'ISO -> ISC'
	    NCHAR = NCHAR + 10
		KI2 = 0
	    IF(KI2. GT. KICSIC) THEN
	      SICFAIL = SICFAIL + 1
	      SICFAILED = 1
	      IF(SIGTSIC .GT. 0.0D0) THEN
	        IPYCFAIL = IPYCFAIL + 1
	        IPYCFAILED = 1
	        PARFAIL = PARFAIL + 1
	        PARFAILED = 1
	        PSTATE = 'FFF'
	        FAILUREPATH(NCHAR+1:NCHAR+14) = ' -> IFF -> FFF'
	        NCHAR = NCHAR + 14
	        FAIL = .TRUE.
	        FAILTYPE(0) = FAILTYPE(0) + 1
	        FMODE = 0
	      ELSE
	        PSTATE = 'ICF'
	        FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> ICF'
	        NCHAR = NCHAR + 7
	        MCODE = 'IS2'
	      END IF
	    ELSE
	      PSTATE = 'ISF'
	      FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> ISF'
	      NCHAR = NCHAR + 7
	      MCODE = 'IS2'
C         Record the elastic response if cracked layer is fully relaxed
		  DO 2352 J = 1, NDIVO + 2
			EPIRCP(NDIVI+NDIVS+4+J) = (-(SIGR(NDIVI+NDIVS+4+J)-
     &				SIGR(NDIV)) + OPYCNU*SIGT(NDIVI+NDIVS+4+J))/OPYCE
			EPITCP(NDIVI+NDIVS+4+J) = (-SIGT(NDIVI+NDIVS+4+J) + 
     &				OPYCNU*(SIGR(NDIVI+NDIVS+4+J)-SIGR(NDIV)))/OPYCE
			URCP(NDIVI+NDIVS+4+J) = 0.0D0
			DO 2353 K = 1, J-1
			  URCP(NDIVI+NDIVS+4+J) = URCP(NDIVI+NDIVS+4+J) + 
     &				EPIRCP(NDIVI+NDIVS+4+K)*(R5-R4)/FLOAT(NDIVO+1)
2353			CONTINUE
2352		  CONTINUE
	    END IF
	  END IF
	ELSE IF(PSTATE .EQ. 'FSO') THEN
C   First update the layer flags
        IPYCFAILED = 2
C
C      KI1 is stress intensity contribution to the SiC layer from KIIPYC
C      The inherent flaw size of the SiC layer on the IPyC side is approximated as 2 micron
C      The inherent flaw size of the SiC layer on the OPyC side is approximated as 1 micron
C
	  KI1 = KIIPYC*SQRT((2+DPD)/(R3-R2)) + 0.413*(1.0D0+R2/R4)*SIGTSIC*
     &		SQRT(PIE*(R3-R2))/SQRT(1.0D0 - (R3-R2)/(R4-R2))
	  WRITE(IKISIC,*) N, OPERTIME, DT, KI1
C
	  IF(SIGTSIC .GT. SIGFSIC) THEN
	    SICFAIL = SICFAIL + 1
	    SICFAILED = 1
	    OPYCFAIL = OPYCFAIL + 1
	    OPYCFAILED = 1
	    PARFAIL = PARFAIL + 1
	    PARFAILED = 1
	    PSTATE = 'FFF'
	    FAILUREPATH(NCHAR+1:NCHAR+14) = ' -> FFO -> FFF'
	    NCHAR = NCHAR + 14
	    FAIL = .TRUE.
	    FAILTYPE(1) = FAILTYPE(1) + 1
	    FMODE = 1
	  ELSE IF(KI1 .GT. KICSIC) THEN
	    SICFAIL = SICFAIL + 1
	    SICFAILED = 1
	    OPYCFAIL = OPYCFAIL + 1
	    OPYCFAILED = 1
	    PARFAIL = PARFAIL + 1
	    PARFAILED = 1
	    PSTATE = 'FFF'
	    FAILUREPATH(NCHAR+1:NCHAR+14) = ' -> FFO -> FFF'
	    NCHAR = NCHAR + 14
	    FAIL = .TRUE.
	    FAILTYPE(0) = FAILTYPE(0) + 1
	    FMODE = 0	
        ELSE IF(SIGTOPYC .GT. SIGFOPYC) THEN
	    OPYCFAIL = OPYCFAIL + 1
	    OPYCFAILED = 1
	    PSTATE = 'FSC'
	    FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> FSC'
	    NCHAR = NCHAR + 7
		KI2 = 0.0D0
	    IF(KI2. GT. KICSIC) THEN
	      SICFAIL = SICFAIL + 1
	      SICFAILED = 1
	      PARFAIL = PARFAIL + 1
	      PARFAILED = 1
	      PSTATE = 'FFF'
	      FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> FFF'
	      NCHAR = NCHAR + 7
	      FAIL = .TRUE.
	      FAILTYPE(0) = FAILTYPE(0) + 1
	      FMODE = 0
	    ELSE
	      PSTATE = 'FSF'
	      FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> FSF'
	      NCHAR = NCHAR + 7
	      MCODE = 'S1'
C         Record the elastic response if cracked layer is fully relaxed
		  DO 2354 J = 1, NDIVO + 2
			EPIRCP(NDIVI+NDIVS+4+J) = (-(SIGR(NDIVI+NDIVS+4+J)-
     &				SIGR(NDIV)) + OPYCNU*SIGT(NDIVI+NDIVS+4+J))/OPYCE
			EPITCP(NDIVI+NDIVS+4+J) = (-SIGT(NDIVI+NDIVS+4+J) + 
     &				OPYCNU*(SIGR(NDIVI+NDIVS+4+J)-SIGR(NDIV)))/OPYCE
			URCP(NDIVI+NDIVS+4+J) = 0.0D0
			DO 2355 K = 1, J-1
			  URCP(NDIVI+NDIVS+4+J) = URCP(NDIVI+NDIVS+4+J) + 
     &				EPIRCP(NDIVI+NDIVS+4+K)*(R5-R4)/FLOAT(NDIVO+1)
2355			CONTINUE
2354		  CONTINUE
	    END IF
	  END IF
	ELSE IF(PSTATE .EQ. 'FCO') THEN
	  IPYCFAILED = 2
	  SICFAILED = 2
	  IF(SIGTSIC .GT. 0.0D0) THEN
	    OPYCFAIL = OPYCFAIL + 1
	    OPYCFAILED = 1
	    PARFAIL = PARFAIL + 1
	    PARFAILED = 1
	    PSTATE = 'FFF'
	    FAILUREPATH(NCHAR+1:NCHAR+14) = ' -> FFO -> FFF'
	    NCHAR = NCHAR + 14
	    FAIL = .TRUE.
	    FAILTYPE(0) = FAILTYPE(0) + 1
          FMODE = 0
	  ELSE IF(SIGTOPYC .GT. SIGFOPYC) THEN
	    OPYCFAIL = OPYCFAIL + 1
	    OPYCFAILED = 1
	    PARFAIL = PARFAIL + 1
	    PARFAILED = 1
	    PSTATE = 'FFF'
	    FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> FFF'
	    NCHAR = NCHAR + 7
	    FAIL = .TRUE.
	    FAILTYPE(1) = FAILTYPE(1) + 1
          FMODE = 1
	  END IF
	ELSE IF(PSTATE .EQ. 'ISF') THEN
C   First update the layer flags
        OPYCFAILED = 2
C
	  KI2 = KIOPYC*SQRT(1.0D0/(R5-R4)) + 0.413*(1.0D0+R5/R3)*SIGTSIC*
     &		SQRT(PIE*(R5-R4))/SQRT(1.0D0 - (R5-R4)/(R5-R3))
C
	  IF(SIGTSIC .GT. SIGFSIC) THEN
	    SICFAIL = SICFAIL + 1
	    SICFAILED = 1
	    OPYCFAIL = OPYCFAIL + 1
	    OPYCFAILED = 1
	    PARFAIL = PARFAIL + 1
	    PARFAILED = 1
	    PSTATE = 'FFF'
	    FAILUREPATH(NCHAR+1:NCHAR+14) = ' -> IFF -> FFF'
	    NCHAR = NCHAR + 14
	    FAIL = .TRUE.
	    FAILTYPE(1) = FAILTYPE(1) + 1
	    FMODE = 1
	  ELSE IF(KI2 .GT. KICSIC) THEN
	    SICFAIL = SICFAIL + 1
	    SICFAILED = 1
	    OPYCFAIL = OPYCFAIL + 1
	    OPYCFAILED = 1
	    PARFAIL = PARFAIL + 1
	    PARFAILED = 1
	    PSTATE = 'FFF'
	    FAILUREPATH(NCHAR+1:NCHAR+14) = ' -> IFF -> FFF'
	    NCHAR = NCHAR + 14
	    FAIL = .TRUE.
	    FAILTYPE(0) = FAILTYPE(0) + 1
	    FMODE = 0
        ELSE IF(SIGTIPYC .GT. SIGFIPYC) THEN
	    IPYCFAIL = IPYCFAIL + 1
	    IPYCFAILED = 1
	    PSTATE = 'CSF'
	    FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> CSF'
	    NCHAR = NCHAR + 7
		KI1 = 0.0D0
	    IF(KI1. GT. KICSIC) THEN
	      SICFAIL = SICFAIL + 1
	      SICFAILED = 1
	      PARFAIL = PARFAIL + 1
	      PARFAILED = 1
	      PSTATE = 'FFF'
	      FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> FFF'
	      NCHAR = NCHAR + 7
	      FAIL = .TRUE.
	      FAILTYPE(0) = FAILTYPE(0) + 1
	      FMODE = 0
	    ELSE
	      PSTATE = 'FSF'
	      FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> FSF'
	      NCHAR = NCHAR + 7
	      MCODE = 'S1'
C         Record to elastic response if cracked layer is fully relaxed
		  DO 2356 J = 1, NDIVI + 2
			EPIRCP(J) = (-(SIGR(J)-SIGR(1)) + IPYCNU*SIGT(J))/IPYCE
			EPITCP(J) = (-SIGT(J) + IPYCNU*(SIGR(J)-SIGR(1)))/IPYCE
			URCP(J) = 0.0D0
			DO 2357 K = J, NDIVI + 1
			  URCP(J) = URCP(J) - EPIRCP(K)*(R3-R2)/FLOAT(NDIVI+1)
2357			CONTINUE
2356		  CONTINUE
	    END IF
	  END IF
	ELSE IF(PSTATE .EQ. 'ICF') THEN
	  OPYCFAILED = 2
	  SICFAILED = 2
	  IF(SIGTSIC .GT. 0.0D0) THEN
	    IPYCFAIL = IPYCFAIL + 1
	    IPYCFAILED = 1
	    PARFAIL = PARFAIL + 1
	    PARFAILED = 1
	    PSTATE = 'FFF'
	    FAILUREPATH(NCHAR+1:NCHAR+14) = ' -> IFF -> FFF'
	    NCHAR = NCHAR + 14
	    FAIL = .TRUE.
	    FAILTYPE(0) = FAILTYPE(0) + 1
          FMODE = 0
	  ELSE IF(SIGTIPYC .GT. SIGFIPYC) THEN
	    IPYCFAIL = IPYCFAIL + 1
	    IPYCFAILED = 1
	    PARFAIL = PARFAIL + 1
	    PARFAILED = 1
	    PSTATE = 'FFF'
	    FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> FFF'
	    NCHAR = NCHAR + 7
	    FAIL = .TRUE.
	    FAILTYPE(1) = FAILTYPE(1) + 1
          FMODE = 1
	  END IF
	ELSE IF(PSTATE .EQ. 'FSF') THEN
	  IPYCFAILED = 2
	  OPYCFAILED = 2
C
	  KI1 = KIIPYC*SQRT((2+DPD)/(R3-R2)) + 0.413*(1.0D0+R2/R4)*
     &		SIGTSIC*SQRT(PIE*(R3-R2))/SQRT(1.0D0 - (R3-R2)/(R4-R2))
	  KI2 = KIOPYC*SQRT(1.0D0/(R5-R4)) + 0.413*(1.0D0+R5/R3)*
     &		SIGTSIC*SQRT(PIE*(R5-R4))/SQRT(1.0D0 - (R5-R4)/(R5-R3))
	  WRITE(IKISIC,*) N, OPERTIME, DT, KI1
C
	  IF(SIGTSIC .GT. SIGFSIC) THEN
	    SICFAIL = SICFAIL + 1
	    SICFAILED = 1
	    PARFAIL = PARFAIL + 1
	    PARFAILED = 1
	    PSTATE = 'FFF'
	    FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> FFF'
	    NCHAR = NCHAR + 7
	    FAIL = .TRUE.
	    FAILTYPE(1) = FAILTYPE(1) + 1
          FMODE = 1
	  ELSE IF((KI1 .GT. KICSIC).OR.(KI2 .GT. KICSIC)) THEN
	    SICFAIL = SICFAIL + 1
	    SICFAILED = 1
	    PARFAIL = PARFAIL + 1
	    PARFAILED = 1
	    PSTATE = 'FFF'
	    FAILUREPATH(NCHAR+1:NCHAR+7) = ' -> FFF'
	    NCHAR = NCHAR + 7
	    FAIL = .TRUE.
	    FAILTYPE(0) = FAILTYPE(0) + 1
          FMODE = 0
	  END IF
	ELSE IF(PSTATE .EQ. 'FFF') THEN
	  FAIL = .TRUE.
	END IF
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C  Subroutine FEEDPEBBLE(COREMODEL, CHANNELS, ENTRANCE, WHICH_CHN,PATH)*
C                                                                      *
C    This subroutine randomly samples the radial position where a      *
C  pebble is fed in from and gives the streamline the pebble will      *
C  flow along.                                                         *
C    Caution: The dimension of CHANNELS and PATH are specified in the  *
C             main program, but their stored data are read in this     *
C             subroutine from channels.dat and blocks.dat, therefore   *
C             they must be consistent. If the input files are changed, *
C             dimension should be changed accordingly.                 * 
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual argument description                                         *
C    COREMODEL      I: indicates whether new or old core model is      *
C					 employed. (currently 1 for new and 2 for old)   *
C    CHANNELS(NAXIAL,-1:NCHANNEL) --                                   *
C                   D: defines boundary lines of every channel         *
C                      format: (z-position,r-position)                 *
C                              (coordinates in cm)                     *
C    ENTRANCE       D: randomly determined radial entrance position    *
C                      (cm)                                            *
C    PATH(NAXIAL,2) --                                                 *
C                   D: defines the streamline of pebble in the format  *
C                      (z, r)                                          *
C    WHICH_CHN      I: indicates which channel the pebble goes in NOW  *
C                                                                      *
C  Local variable description                                          *
C    MODERATOR      D: thickness of center moderator (cm)              *
C    R              D: randomly sampled r position (cm)                *
C    WIDTH          D: width of reactor core (cm)                      *
C                                                                      *
C  Parameters                                                          *
C    NCHANNEL1      D: number of channels from new VSOP model          *
C    NAXIAL1        D: number of axial divisions from new VSOP model   *
C    NCHANNEL2      D: number of channels from old VSOP model          *
C    NAXIAL2        D: number of axial divisions from old VSOP model   *
C                                                                      *
C  Function called                                                     *
C    RAND            D: random number generator (0,1]                   *
C    DSQRT          D: calculate square root                           *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE FEEDPEBBLE(COREMODEL,CHANNELS,ENTRANCE,WHICH_CHN,PATH)
C
      DOUBLE PRECISION CHANNELS, ENTRANCE, PATH
      DOUBLE PRECISION MODERATOR, R, WIDTH
      INTEGER    COREMODEL, WHICH_CHN
C
	PARAMETER ( NCHANNEL1 = 5, NAXIAL1 = 21)
	PARAMETER ( NCHANNEL2 = 5, NAXIAL2 = 15)
C
	DIMENSION CHANNELS(1:,-1:), PATH(1:,1:)
C
C      EXTERNAL  RAND
C
	IF (COREMODEL .EQ. 1) THEN
        WIDTH = CHANNELS(1,NCHANNEL1)
        MODERATOR = CHANNELS(1,0)
        CALL RANDOM_NUMBER(RN)
2405    R = WIDTH*RN
        IF(R.LE.MODERATOR) THEN !Pebbles don't go into moderator, so sample again
          GO TO 2405
        END IF
        CALL LOCATE_ARRAY(CHANNELS(1,0:5),6,1,R,WHICH_CHN)
C  Calculate the streamline of pebble
        PATH(1,1) = CHANNELS(1,-1)
	  PATH(1,2) = R
        DO 2410 I = 2, NAXIAL1
          PATH(I,1) = CHANNELS(I,-1)     !z-position
          PATH(I,2) = ((PATH(1,2)-CHANNELS(1,WHICH_CHN-1))*
     &              CHANNELS(I,WHICH_CHN)+(CHANNELS(1,WHICH_CHN)-
     &              PATH(1,2))*CHANNELS(I,WHICH_CHN-1))/
     &              (CHANNELS(1,WHICH_CHN)-CHANNELS(1,WHICH_CHN-1))
2410    CONTINUE
        ENTRANCE = R
	ELSE IF (COREMODEL .EQ. 2) THEN
        WIDTH = CHANNELS(1,NCHANNEL2) - CHANNELS(1,0)
        MODERATOR = CHANNELS(1,1)
2415  CALL RANDOM_NUMBER(RN)  
      R = WIDTH*RN
        IF(R.LE.MODERATOR) THEN !Pebbles don't go into moderator, so sample again
          GO TO 2415
	  ELSE IF(R.LE.CHANNELS(1,2)) THEN
	  CALL RANDOM_NUMBER(RN)
          IF(RN.LE.0.5D0) GO TO 2415 !Channel two is graphite/pebble 50/50 mixed
        END IF
        CALL LOCATE_ARRAY(CHANNELS(1,1:5),5,2,R,WHICH_CHN)
C  Calculate the streamline of pebble
        PATH(1,1) = CHANNELS(1,-1)
	  PATH(1,2) = R
        DO 2420 I = 2, NAXIAL2
          PATH(I,1) = CHANNELS(I,-1)     !z-position
          PATH(I,2) = ((PATH(1,2)-CHANNELS(1,WHICH_CHN-1))*
     &              CHANNELS(I,WHICH_CHN)+(CHANNELS(1,WHICH_CHN)-
     &              PATH(1,2))*CHANNELS(I,WHICH_CHN-1))/
     &              (CHANNELS(1,WHICH_CHN)-CHANNELS(1,WHICH_CHN-1))
2420    CONTINUE
        ENTRANCE = R
	ELSE
	  CALL ERR_HANDLER(TRIM('COREMODEL not specified when calling 
     &					FEEDPEBBLE'), 45,2,2,IERR)
	END IF
	RETURN
	END
C                                                                      *
C***********************************************************************
C
C
C      
C***********************************************************************
C  Subroutine GASRLS(T_PARTICLE, BURNUP, OPERTIME, DIFFUSION, PRESSURE)*
C                                                                      *
C    This subroutine calculates the fraction of fission gas released,  *
C  the Oxygen produced, and the internal pressure by ideal gas law as  *
C  a function of temperature, burnup and reactor operating time.       *
C    The reference to this subroutine is R. Gontard, H. Nabielek, "    *
C  Performance Evaluation of Modern HTR TRISO Fuel," HTA-IB-05/90,     *
C  KFA, Julich, Germany, Appendix C.                                   *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual argument description                                         *
C    T_PARTICLE(0:5)D: Temperature distribution in the particles in C  *
C    BURNUP (FIMA)  D: Current burnup                                  *
C    OPERTIME (sec) D: Reactor operating time                          *
C    DIFFUSION      L: Flag indicating whether to take into account    *
C                      diffusion in gas release model or not           *
C  Returned value                                                      *
C    PRESSURE (MPa) D: Internal pressure                               *
C                                                                      *
C  Variables in COMMON blocks                                          *
C  /PAR_R/ --                                                          *
C    R1-R5 (um)     D: TRISO fuel dimension                            *
C  /PAR_K/ --                                                          *
C    KDEN (g/cm^3)  D: Fuel kernel density                             *
C    BDEN (g/cm^3)  D: Buffer density                                  *
C    KERNT (g/cm^3) D: Fuel kernel theoretical density                 *
C    BUFFT (g/cm^3) D: Buffer theoretical density                      *
C    FUELTYPE       C: Type of fuel evaluated ('UCO','UO2', or         *
C                      '(Th,U)O2')                                     *
C    U235E          D: U235 enrichment                                 *
C    CURAT          D: Carbon to Uranium weight ratio in kernel        *
C    OURAT          D: Oxygen to Uranium weight ratio in kernel        *
C  Local variable description                                          *
C    TKERNEL (C)    D: Fuel kernel temperature                         *
C    TBUFFER (C)    D: Buffer temperature                              *
C    DCOEFFI (sec^-1) --                                               *
C                   D: Reduced diffussion coefficient for noble gases  *
C                      during irradiation                              *
C    DCOEFFA (sec^-1) --                                               *
C                   D: Reduced diffussion coefficient for noble gases  *
C                      during annealing                                *
C    TAUI           D: Dimensionless varible for irradiation,          *
C                      TAUI=DCOEFFI*OPERTIME                           *
C    TAUA           D: Dimensionless varible for annealing,            *
C                      TAUA=DCOEFFA*ANNEALTIME                         *
C    OPF            D: Oxygen atoms per fission                        *
C    VKERNEL (m^3)  D: Kernel volume                                   *
C    AWU (g/mol)    D: Atomic weight of Uranium mix                    *
C    WFU            D: Weight fraction of Uranium in fuel kernel       *
C    MOLESU (mole)  D: Moles of Uranium in fuel kernel                 *
C    VMOLAR (m^3/mol) --                                               *
C                   D: Molar volume of heavy metal in the kernel       *
C    VOIDVOL (m^3)  D: Available void volume for fission gases         *
C    FD             D: Fractional release of fission gases             *
C  Parameters and counters:                                            *
C    ANNEALTIME (sec) --                                               *
C                   D: Annealing time required by the formulation.     *
C                      It cannot be zero, but we can set a trivial     *
C                      value like 1 second.                            *
C    TANNEAL (C)    D: Annealing temperature                           *
C    AWC (g/mol)    D: Atomic weight of C12                            *
C    AWO (g/mol)    D: Atomic weight of O16                            *
C    AWU235 (g/mol) D: Atomic weight of U235                           *
C    AWU238 (g/mol) D: Atomic weight of U238                           *
C    VOIDF          D: Fraction of voids in the buffer layer           *
C    YIELD          D: Yield of stable fission gases                   *
C    R              D: Gas constant (J/(mol*C))                        *
C    PIE            D: 3.14159...                                      *
C    CV             D: Constant for volume calculation of spherical    *
C                      geometry (4/3)*PIE*10^(-18) m^3/um^3            *
C    ORDERI,ORDERA  D: Summands in the calculation of FD               *
C    FI,FA          D: The summations in the calculation of FD         *
C    N              I: Order of summations                             *
C    LIMIT          I: Highest order taken into account                *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE GASRLS(T_PARTICLE,BURNUP,OPERTIME,DIFFUSION,PRESSURE)
C
      DOUBLE PRECISION T_PARTICLE, BURNUP, OPERTIME, PRESSURE
      DOUBLE PRECISION TKERNEL, TBUFFER, DCOEFFI, TAUI, OPF, VKERNEL, 
     &                 VOIDVOL, FD
      DOUBLE PRECISION R1, R2, R3, R4, R5
	DOUBLE PRECISION KDEN, BDEN, KERNT, BUFFT, U235E, CURAT, OURAT
      DOUBLE PRECISION R, AWC, AWO, AWU235, AWU238, YIELD, VOIDF,
     &                 OPFMAX, PIE, CV, AWU, WFU, MOLESU, VMOLAR
	DOUBLE PRECISION ORDERI, ORDERA, FI, FA
      DOUBLE PRECISION ANNEALTIME, TANNEAL, DCOEFFA, TAUA
      INTEGER*4 N, LIMIT
	LOGICAL DIFFUSION
      CHARACTER*3  FUELTYPE
C
      PARAMETER (  R = 8.3144 D0,         !(J/(mole*K))
     &             AWC = 12.01115 D0,
     &             AWO = 15.9994 D0,
     &             AWU235 = 235.0439 D0,
     &             AWU238 = 238.0508 D0,
     &             YIELD = 0.31 D0,
     &             VOIDF = 0.5 D0,
     &             OPFMAX = 0.625D0,
     &             PIE = 3.1415926535897932385 D0,
     &             CV = 4.1887 90204 78639 09846 D-18, !Including the conversion from um to m
     &             ANNEALTIME = 1.0D0,    !1.0 sec
     &             TANNEAL = 1250.0D0,    !degree C
     &             LIMIT = 100)
C
      DIMENSION T_PARTICLE(0:5)
C
      COMMON /PAR_R/ R1, R2, R3, R4, R5   !particle geometry
      COMMON /PAR_K/ KDEN, BDEN, KERNT, BUFFT, U235E, CURAT, OURAT,
     &               FUELTYPE
C
      IF (DIFFUSION) THEN
C  Pick out kernel temperature from temperature distribution T_PARTICLE(0:5)
        TKERNEL = T_PARTICLE(0)
	  TBUFFER = T_PARTICLE(1)
C  Calculate reduced diffusion coefficient for noble gases in the kernel
        IF((FUELTYPE.EQ.'UCO').OR.(FUELTYPE.EQ.'UO2')) THEN
          DCOEFFI = 10.0D0**(-2.30D0-8.116D3/(TKERNEL+273.0D0)) !Valid for UO2 and UCO
          DCOEFFA = 10.0D0**(-2.30D0-8.116D3/(TANNEAL+273.0D0))
        ELSE  !For (Th,U)O2 fuel
          DCOEFFI = 10.0D0**(-5.94D0+3.24D0/(1.0D0+0.11D0/BURNUP)-
     &                        5.460D3/(TKERNEL+273.0D0))
          DCOEFFA = 10.0D0**(-5.94D0+3.24D0/(1.0D0+0.11D0/BURNUP)-
     &                        5.460D3/(TANNEAL+273.0D0))
        END IF
        TAUI = DCOEFFI*OPERTIME                     !dimensionless argument
        TAUA = DCOEFFA*ANNEALTIME
C  Calculate fractional release of fission gases Xe and Kr (Allelein 1983)
        FI = 1.0D0
        FA = 1.0D0
        N=1
        DO WHILE(N.LE.LIMIT)
          ORDERI = (1.0D0-DEXP(-(TAUI+TAUA)*N**2*PIE**2)/((N*PIE)**4))
          FI = FI - 6.0D0*ORDERI/(TAUI+TAUA)
          ORDERA = (1.0D0-DEXP(-TAUA*N**2*PIE**2)/((N*PIE)**4))
          FA = FA - 6.0D0*ORDERA/TAUA
	    N = N + 1
        END DO
C       TAUI might be zero when OPERTIME is zero
	  IF(OPERTIME .EQ. 0.0D0) THEN
	    FD = 0.0D0
	  ELSE
          FD = ((TAUI+TAUA)*FI - TAUA*FA)/TAUI
	  END IF
C  Number of Oxygen atoms per fission OPF and molar volume of heavy metal
        IF(FUELTYPE.EQ.'UCO') THEN
          OPF = 0.0D0
          VMOLAR = 2.50654 D-5
        ELSE IF(FUELTYPE.EQ.'UO2') THEN
C          OPF = 10.0D0**(-10.08D0-8.5D3/(TKERNEL+273.0D0))*(OPERTIME)**2
C  The following OPF is from Kazuhiro Sawa for comparison with HTTR
          OPF = 10.0D0**(-0.21D0-8.5D3/(TKERNEL+273.0D0))*
     &          (OPERTIME/86400.0D0)**2
          VMOLAR = 2.43796 D-5
        ELSE  !For (Th,U)O2
          OPF = 10.0D0**(0.96D0-4.42D3/(TKERNEL+273.0D0) +
     &          0.4*DLOG10(10.0D0) + 0.3*DLOG10(BURNUP))
          VMOLAR = 2.51905 D-5
        END IF
	  IF(OPF.GT.OPFMAX) OPF = OPFMAX
C  Kernel volume and void volume
        VKERNEL = CV*R1**3
C        VOIDVOL = CV*(R2**3 - R1**3)*VOIDF
        VOIDVOL = CV*((BUFFT - BDEN)/BUFFT)*(R2**3 - R1**3)
C  Calculate internal pressure
        PRESSURE=(FD*YIELD+OPF)*BURNUP*(VKERNEL/VMOLAR)*R*
     &           (TBUFFER+273.0D0)/VOIDVOL
        PRESSURE = PRESSURE/1.0D6     !convert to MPa
	ELSE
	  TBUFFER = T_PARTICLE(1)
	  VOIDVOL = CV*(((BUFFT - BDEN)/BUFFT)*(R2**3 - R1**3) +
     &            (KERNT - KDEN)/KERNT*R1**3)
        AWU = AWU235*(U235E/100.D0) + AWU238*(100.D0-U235E)/100.D0
        WFU = AWU/(AWU+CURAT*AWC+OURAT*AWO)
        MOLESU = WFU*KDEN*(4.D0*PIE/3.D0)*R1**3*(1.D-12)/AWU  ! convert um to cm
        PRESSURE = YIELD*BURNUP*MOLESU*R*(TBUFFER+273.0D0)/VOIDVOL
	  PRESSURE = PRESSURE/1.0D6
	END IF
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
C  Subroutine LOCATE                                                   *
C                                                                      *
C  Indices for look-up table                                           *
C  Subroutine LOCATE comes from W. Press, Numerical Recipes,           *
C    Cambridge, Cambridge Univ. Press, 1986.                           *
C                                                                      *
C  Subroutine LOCATE is used to bracket the x-position for which f(x)  *
C  is sought by two adjecent tabulated positions.  That is, given a    *
C  monotonic array of Vi, and given a value x, it finds the two        *
C  values Vi, Vi+1 that bracket x by a bisection method.  If the       *
C  array V does not bracket x the returned index is zero.              *
C                                                                      *
C  Description of variables                                            *
C                                                                      *
C  IV        :  Integer input array (look-up table)                    *
C  N         :  Dimension of array IV                                  *
C  IX        :  Element of look-up table                               *
C  J         :  Returned lower index                                   *
C                                                                      *
C  Description of local variables                                      *
C                                                                      *
C  JL        :  Lower index                                            *
C  JU        :  Upper index                                            *
C  JM        :  Intermediate index                                     *
C                                                                      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE LOCATE (IV, N, IX, J)
C
      INTEGER IV, IX, N, J, JL, JU, JM
C
      DIMENSION IV(N)
C
C  Taken from Press, et.al.
C
C
C  Initialize upper and lower limits
C
      JL = 0
      JU = N + 1
C
C  If we are not yet done, compute a midpoint and replace either the
C  lower or the upper limit, as appropriate
C
2510    IF (JU - JL .GT. 1) THEN
        JM = (JU + JL)/2
        IF ((IV(N) .GT. IV(1)) .EQV. (IX .GE. IV(JM))) THEN
          JL = JM
        ELSE
          JU = JM
        ENDIF
C
C  Repeat until the test condition in line number 10 is satisfied.
C
        GO TO 2510
      ENDIF
C
C  Then set the output and return
C
      J = JL
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
C  Subroutine LOCATE_ARRAY                                             *
C                                                                      *
C  Indices for look-up table                                           *
C  Subroutine LOCATE_ARRAY comes from W. Press, Numerical Recipes,     *
C    Cambridge, Cambridge Univ. Press, 1986.                           *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  This subroutine is used to bracket the x-position for which f(x)    *
C  is sought by two adjecent tabulated positions.  That is, given a    *
C  monotonic array of Vi, and given a value x, it finds the two        *
C  values Vi, Vi+1 that bracket x by a bisection method.  If the       *
C  array V does not bracket x the returned index is zero.              *
C                                                                      *
C  Actual variable description                                         *
C    AI              D: Input array (look-up table)                    *
C    N               I: Dimension of array AI                          *
C    LI              I: Lower index of input array                     *
C    X               D: Element of look-up table                       *
C    J               I: Returned lower index                           *
C                                                                      *
C  Local variable description                                          *
C    JL              I: Lower limit                                    *
C    JU              I: Upper limit                                    *
C    JM              I: Intermediate index                             *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE LOCATE_ARRAY (AI, N, LI, X, J)
C
      DOUBLE PRECISION AI, X
      INTEGER N, LI, J, JL, JU, JM
C
      DIMENSION AI(1:N)
C
C  Initialize upper and lower limits
C
      JL = 0
      JU = N + 1
C
C  If we are not yet done, compute a midpoint and replace either the
C  lower or the upper limit, as appropriate
C
      DO WHILE(JU - JL .GT. 1)
        JM = (JU + JL)/2
        IF ((AI(N) .GT. AI(1)) .EQV. (X .GE. AI(JM))) THEN
          JL = JM
        ELSE
          JU = JM
        END IF
      END DO
C
C  If JL is below the lower bound or above the upper bound, adjust it.
	IF(JL .LT. 1) THEN
	  JL = JL + 1
	ELSE IF(JL .GE. N) THEN
	  JL = JL - 1
	END IF
C  Set the output and return
C  Notice: Array AI starts from indice 1, but input array starts from LI,
C          hence we need to adjust output to the right position in the 
C          original array.
      J = JL - (1 - LI)
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
C  Subroutine MACHINE                                                  *
C                                                                      *
C    Initializes computer and compiler type for use in Date and Time   *
C  routines                                                            *
C                                                                      *
C  To activate a computer and compiler type, remove the corresponding  *
C  comment cards.                                                      *
C                                                                      *
C  Variables in COMMON blocks                                          *
C  /MTYPE/ --                                                          *
C    MACH           C: Three character description of computer type    *
C    COMP           C: Three character description of compiler type    *
C                                                                      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE MACHINE
C
      CHARACTER*3 MACH, COMP
C
      COMMON /MTYPE/ MACH, COMP
C#######################################################################
C  Macintosh computer
C      MACH = 'MAC'  
C  IBM/IBM compatible computer
      MACH = 'IBM' 
C  DEC work station
C      MACH = 'DEC'
C
C  No compiler specified
C      COMP = '   ' 
C  Language system compiler for Macintosh
C      COMP = 'LSC'
C  Absoft 0/20 compiler for Macintosh
C      COMP = '020'
C  Absoft MPW compiler for Macintosh
C      COMP = 'MPW' 
C  Lahey compiler for IBM
C      COMP = 'LAH'
C  Visual Fortran 5.0 for WindowsNT 4.0
      COMP = 'VF5'
C  DEC Ultrix Fortran compiler, Version 3.2.1 or greater
CDEC  COMP = 'ULT'
C#######################################################################
C
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
C  SUBROUTINE OPENFILE                                                 *
C                                                                      *
C    Opens file without leading or trailing blanks in filename         *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual variable description                                         *
C    FILENAME       C: File name (up to 256 characters)                *
C    FILETYPE       C: File type (always 4 non-blank characters: '.xxx'*
C    STAT           C: Status of file: 'OLD' or 'NEW'                  *
C    LUN            L: Logical unit number                             *
C                                                                      *
C  Local variable description                                          *
C    I              I: Loop counter                                    *
C    J              I: Tally                                           *
C    NAME           C: Catenated file name:                            *
C                      NAME = (FILENAME - all blank spaces)//FILETYPE  *
C    NB             I: Number of blank spaces in FILENAME              *
C    IERR           I: Logical Unit Number of error file 'ERROR.MSG'   *
C    EX             L: Indicates if the file exists                    *
C    KEYSTROKE      L: Keep the character inputed from keyboard        *
C                                                                      *
C  Functions and subroutines called                                    *
C                                                                      *
C    ERR_HANDLER     : Processes an error (diagnostic) message.        *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE OPENFILE(FILENAME,FILETYPE,STAT,LUN)
C
      INTEGER*4 LUN, I, J, IERR
      CHARACTER*(*) FILENAME
      CHARACTER*4 FILETYPE
      CHARACTER*3 STAT
      CHARACTER*260 NAME
	CHARACTER*1 KEYSTROKE
	LOGICAL EX
C
	PARAMETER (IERR = 12)
C
C  Reconstruct the file name without blank spaces
C
      I = LEN(FILENAME)
	J=0
      DO WHILE(FILENAME(J+1:J+1).NE.' ')
	  J = J + 1
	  NAME(J:J) = FILENAME(J:J)
	  IF(J.EQ.I) EXIT
	END DO
      NAME(J+1:J+4) = FILETYPE
C
C  Generation of file name complete, open file NAME
C
	INQUIRE (FILE = NAME(1:J+4), EXIST = EX)
	IF (EX .AND. (STAT .EQ. 'NEW')) THEN 
	  WRITE(*,*) 'File ',NAME(1:J+4),
     &             ' already exists, overwrite it? (Y/N)'
	  READ(*,*) KEYSTROKE
	  IF((KEYSTROKE.EQ.'Y').OR.(KEYSTROKE.EQ.'y')) THEN
	    OPEN(UNIT = LUN, FILE = NAME(1:J+4), STATUS = 'OLD')
c		STAT = 'OLD'
	  ELSE
	    CALL ERR_HANDLER(TRIM('OPENFILE: File already exists'),30,2,2,IERR)
	  END IF
	ELSE IF((.NOT.EX) .AND. (STAT.EQ.'OLD')) THEN
	  WRITE(*,*) 'File ',NAME(1:J+4),' does not exist.'
	  CALL ERR_HANDLER(TRIM('OPENFILE: File does not exist'),29,1,1,IERR)
C	  STAT = 'ERR'
	ELSE
        OPEN(UNIT = LUN, FILE = NAME(1:J+4), STATUS = STAT)
	END IF
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C  Subroutine TEMPERATURE(QPPP, T_HE, BURNUP, T_PARTICLE)              *
C                                                                      *
C    This subroutine calculates the temperature distribution through   *
C  the fuel particles.                                                 *
C    Notes:                                                            *
C  1.  Call TEMPERATURE(0.0D0, T_HE, BURNUP, T_PARTICLE) for           *
C    initialization.                                                   * 
C  2.  Currently there is no gap considered in the model               *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual argument description                                         *
C    QPPP (W/m^3)   D: Current power density                           *
C    T_HE (C)       D: Current Helium temperature                      *
C    BURNUP (FIMA)  D: Current burnup of the pebble                    *
C    T_PARTICLE(0:5) D: (returned arguments)                           *
C                      Temperature distribution in the particles       *
C                      T_PARTICLE(0) is the value at the center of fuel*
C                      T_PARTICLE(1-5) correspond to the values at     *
C                      R1-R5, respectively.                            *
C                                                                      *
C  Variables in COMMON blocks                                          *
C  /PBED/ : reactor core parameters                                    *
C    CORE_HEIGHT(m) D: Height of reactor core                          *
C    CORE_RADIUS(m) D: Radius of reactor core                          *
C    P_CORE (MWth)  D: Thermal power of reactor                        *
C    T_GASIN (C)    D: Coolant (He) entry temperature                  *
C    T_GASOUT (C)   D: Coolant (He) exit temperature                   *
C    MF_HE (kg/s)   D: Mass flow rate of Helium                        *
C    PACKING        D: Packing fraction of pebbles in the reactor core *
C    NPEBBLE        I: Number of pebbles in the reactor core           *
C  /PEBBLE/ : pebble parameters                                        *
C    PEBRADIUS (m)  D: Pebble radius                                   *
C    PFZRADIUS (m)  D: Pebble fuel zone radius                         *
C    K_PM (W/m.K)   D: Thermal conductivity of pebble matrix           *
C    K_PFM (W/m.K)  D: Thermal conductivity of pebble fuel zone matrix *
C    K_PFZ (W/m.K)  D: Thermal conductivity of pebble fuel zone        *
C    R_IN_PEBBLE (m)D: The position of the sampled particle in pebble  *
C    NPARTICLE      I: Average number of particles in each pebble      *
C  /PAR_R0/ : fuel particle initial geometry                           *
C    R10 - R50 (um) D: Initial radii of layers                         *
C  /PAR_R/ :  fuel particle current geometry                           *
C    R1  - R5 (um)  D: Current radii of layers                         *
C  /PAR_K/ --                                                          *
C    KDEN (g/cm^3)  D: Fuel kernel density                             *
C    BDEN (g/cm^3)  D: Buffer density                                  *
C    KERNT (g/cm^3) D: Fuel kernel theoretical density                 *
C    BUFFT (g/cm^3) D: Buffer theoretical density                      *
C    FUELTYPE       C: Type of fuel evaluated ('UCO','UO2', or         *
C                      '(Th,U)O2')                                     *
C    U235E          D: U235 enrichment                                 *
C    CURAT          D: Carbon to Uranium weight ratio in kernel        *
C    OURAT          D: Oxygen to Uranium weight ratio in kernel        *
C  /PAR_T/ :  thermal properties of particle materials                 *
C    K_FUEL (W/m.K) D: Thermal conductivity of fuel kernel             *
C    K_BUFFER (W/m.K) D:                                               *
C                      Thermal conductivity of buffer                  *
C    K_IPYC (W/m.K) D: Thermal conductivity of IPyC layer              *
C    K_SIC (W/m.K)  D: Thermal conductivity of SiC layer               *
C    K_OPYC (W/m.K) D: Thermal conductivity of OPyC layer              *
C  /ERRHANDLE/ : handle errors                                         *
C    ERROR_CODE     I: Code symbolizes the degree of error:            *
C                      2 means this is an unconditionally fatal error. *
C                      1 means this is a recoverable error.            *
C                      0 means this is a warning message only.         *
C                                                                      *
C  Local variable description                                          *
C    A_CORE (m^2)   D: Cross-section area of pebble bed                *
C    V_CORE (m^3)   D: Volume of pebble bed                            *
C    PEBDIAMETER(m) D: Pebble diameter                                 *
C    V_PEBBLE (m^3) D: Volume of a pebble                              *
C    V_PFZ (m^3)    D: Volume of pebble fuel zone                      *
C    V_PFM (m^3)    D: Volume of pebble fuel zone matrix               *
C    V_PM (m^3)     D: Volume of pebble matrix                         *
C    Q_PFZ (W/m^3)   : Volumetric heat generation rate in pebble       *
C                      fuel zone                                       *
C    T_PEBBLE (C)   D: Array storing temperature profile in a pebble   *
C                      T_PEBBLE(0) is the value at the center of pebble*
C                      T_PEBBLE(1) is the value at fuel zone/matrix    *
C                      boundary.                                       *
C                      T_PEBBLE(2) is the value at pebble surface.     *
C                      R1-R5, respectively.                            *
C    V_PARTICLE (m^3) D:                                               *
C                      Volume of a fuel particle                       *
C    V_FUEL (m^3)   D: Volume of the fuel kernel                       *
C    SWELL_FUEL     D: Swelling of fuel kernel as a function of burnup *
C    V_BUFFER (m^3) D: Volume of buffer                                *
C    V_IPYC (m^3)   D: Volume of IPyC layer                            *
C    V_SIC (m^3)    D: Volume of SiC layer                             *
C    V_OPYC (m^3)   D: Volume of OPyC layer                            *
C    UO2_THERMAL (W/m.K) --                                            *
C                   D: Array of UO2 thermal conductivity as a function *
C                      of temperature                                  *
C    UCO_THERMAL (W/m.K) --                                            *
C                   D: Array of UCO thermal conductivity as a function *
C                      of temperature                                  *
C    PYC_THERMAL (W/m.K) --                                            *
C                   D: Array of PyC thermal conductivity as a function *
C                      of temperature                                  *
C    SIC_THERMAL (W/m.K) --                                            *
C                   D: Array of SiC thermal conductivity as a function *
C                      of temperature                                  *
C    Q_FUEL (W/m^3) D: Volumetric heat generation rate of fuel kernel  *
C    T_PS (C)       D: Calculated particle surface temperature         *
C    HE_THERMAL     D: Array storing thermal properties of Helium      *
C                      HE_THERMAL(:,0) (C): temperature                *
C                      HE_THERMAL(:,1) (W/m.K): thermal conductivity   *
C                      HE_THERMAL(:,2) (kg/m^3): density               *
C                      HE_THERMAL(:,3) (kg/m.s): viscosity             *
C    CP_HE (J/kg.K) D: Specific heat of Helium                         *
C    D_HE (kg/m^3)  D: Density of Helium                               *
C    K_HE (W/m.K)   D: Thermal conductivity of Helium                  *
C    MU_HE (kg/m.s) D: Viscosity of Helium                             *
C    VC_HE (m/s)    D: Characteristic velocity of Helium flow          *
C    RE_HE          D: Reynold's number of He                          *
C    PR_HE          D: Prandtl number of He                            *
C    H_HE (W/m^2.K) D: Heat transfer coefficient of He                 *
C                                                                      *
C  Others                                                              *
C    T (C)          D: temporary temperature                           *
C    R (m or um)    D: temporary radius                                *
C    INDEX          I: Index for position locating in arrays           *
C                                                                      *
C  Parameters, counters and constants                                  *
C    BU_CONV (MWd/T/FIMA) : Convert burnup from FIMA to MWd/T          *
C    PIE            D: 3.14159.....                                    *
C    IERR           I: Unit number of error message file 'ERROR.MSG'   *
C    I              I: Lun for loops                                   *
C                                                                      *
C  Functions and subroutines called                                    *
C    LOCATE_ARRAY   S:                                                 *
C    ERR_HANDLER    S:                                                 *
C    TPEBBLE        F:                                                 *
C    TPARTICLE      F:                                                 *
C    RAND           F:                                                 *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE TEMPERATURE(QPPP, T_HE, BURNUP, T_PARTICLE)
C
	DOUBLE PRECISION QPPP, T_HE, BURNUP, T_PARTICLE
C  Core parameters
	DOUBLE PRECISION CORE_HEIGHT, CORE_RADIUS, A_CORE, V_CORE,
     &                 P_CORE, T_GASIN, T_GASOUT, MF_HE, PACKING
C  Pebble related variables
	DOUBLE PRECISION PEBRADIUS, PEBDIAMETER, PFZRADIUS, V_PEBBLE,
     &                 V_PFZ, V_PFM, V_PM, Q_PFZ, K_PM, K_PFM, K_PFZ,
     &                 T_PEBBLE, R_IN_PEBBLE
C  Particle geometric variables
	DOUBLE PRECISION R10, R20, R30, R40, R50, R1, R2, R3, R4, R5
      DOUBLE PRECISION V_PARTICLE, V_FUEL, SWELL_FUEL, V_BUFFER,
     &                 V_IPYC, V_SIC, V_OPYC
C  Particle kernel specifications
	DOUBLE PRECISION KDEN, BDEN, KERNT, BUFFT, U235E, CURAT, OURAT
C  Particle thermal properties
	DOUBLE PRECISION UO2_THERMAL, UCO_THERMAL, 
     &                 PYC_THERMAL, SIC_THERMAL
	DOUBLE PRECISION Q_FUEL, K_FUEL, K_BUFFER, K_IPYC, K_SIC, K_OPYC,
     &                 T_PS
C  Helium thermal properties
	DOUBLE PRECISION HE_THERMAL, CP_HE, D_HE, MU_HE, K_HE, VC_HE,
     &                 RE_HE, PR_HE, H_HE
C  Others
	DOUBLE PRECISION T, R
	INTEGER*4 NPEBBLE, NPARTICLE
      INTEGER   ERROR_CODE, INDEX, I
      CHARACTER*3  FUELTYPE
C
      PARAMETER (BU_CONV = 9.02298D5, PIE = 3.1415926535897932385 D0)
	PARAMETER (IERR = 12)
C
      DIMENSION T_PARTICLE(0:5), T_PEBBLE(0:2)
	DIMENSION UO2_THERMAL(1:16,0:1), UCO_THERMAL(1:14,0:1),
     &          PYC_THERMAL(1:6,0:1), SIC_THERMAL(1:4,0:1),
     &          HE_THERMAL(1:7,0:3)
C  Common blocks
      COMMON /PBED/ CORE_HEIGHT, CORE_RADIUS, P_CORE, T_GASIN, 
     &              T_GASOUT, MF_HE, PACKING, NPEBBLE
      COMMON /PEBBLE/ PEBRADIUS, PFZRADIUS, K_PM, K_PFM, K_PFZ,
     &                R_IN_PEBBLE, NPARTICLE
      COMMON /PAR_K/ KDEN, BDEN, KERNT, BUFFT, U235E, CURAT, OURAT,
     &               FUELTYPE
      COMMON /PAR_R0/ R10, R20, R30, R40, R50   !initial particle geometry
      COMMON /PAR_R/ R1, R2, R3, R4, R5    !particle geometry
	COMMON /PAR_T/ K_FUEL, K_BUFFER, K_IPYC, K_SIC, K_OPYC
	COMMON /ERRHANDLE/ ERROR_CODE
C
C	EXTERNAL RAND
      EXTERNAL TPARTICLE
	EXTERNAL TPEBBLE
C  Save pebble temperature profile for next step usage.
	SAVE T_PEBBLE
C  Fuel UO2 or UCO conductivity data are from El-Wakil; temperatures are in C.
      DATA UO2_THERMAL /92.96,204.07,315.18,426.29,537.41,648.52,
     &                  759.63,870.74,981.85,1092.96,1204.07,1315.18,
     &                  1426.29,1537.41,1648.52,1759.63,
     &                  7.79,6.06,4.85,4.33,3.81,3.46,2.77,2.60,
     &                  2.42,2.25,2.08,1.90,1.90,1.90,1.90,1.90/
      DATA UCO_THERMAL /92.96,148.52,204.07,259.63,315.18,
     &                  370.74,426.29,481.85,537.41,648.52,
     &                  759.68,870.74,981.85,1092.96,
     &                  25.56,24.35,23.33,22.53,21.93,21.44,21.10,
     &                  20.80,20.61,20.46,20.35,20.25,20.20,20.02/
C  Helium thermal data are from A. F. Mills; temperatures are in C.
	DATA HE_THERMAL /-23.15,26.85,126.85,226.85,326.85,526.85,726.85,
     &                 0.133, 0.149, 0.178, 0.205, 0.229, 0.273, 0.313,
     &                 0.195,0.1624,0.1218,0.0974,0.0812,0.0609,0.0487,
     &         17.9D-6,20.1D-6,24.4D-6,28.2D-6,31.7D-6,37.8D-6,43.3D-6/
C  Dense PyC thermal conductivity data are from D. Dobranich; temperatures are in C.
	DATA PYC_THERMAL /76.85,116.85,216.85,326.85,456.85,556.85,
     &                  9.2D0, 7.5D0, 6.7D0, 5.9D0, 5.4D0, 5.2D0/
C  SiC thermal conductivity data are from Vollman & Hanson; temperatures are in C.
      DATA SIC_THERMAL /500.0D0,900.0D0,1100.0D0,1300.0D0,
     &                   39.8D0, 36.5D0,  35.7D0,  35.5D0/
C
C  First time run: initialize arrays storing temperature profiles
      IF (QPPP .EQ. 0.0D0) THEN
	  DO 3110 I = 0, 2
	    T_PEBBLE(I) = T_HE
3110    CONTINUE
        DO 3120 I = 0, 5
	    T_PARTICLE(I) = T_HE
3120    CONTINUE
        RETURN
	END IF
C  Modeling begins here ...
C
C  Thermal properties for coolant He (from A. F. Mills)
      CP_HE = 5200.0D0                                                 ! J/kg.K
      CALL LOCATE_ARRAY(HE_THERMAL(:,0), 7, 1, T_HE, INDEX)
	IF (INDEX .EQ. 0) INDEX = INDEX + 1   !exceeds lowerbound; make adjustion.
      IF (INDEX .EQ. 7) INDEX = INDEX - 1   !exceeds upperbound; make adjustion.
      K_HE = ((T_HE-HE_THERMAL(INDEX,0))*HE_THERMAL(INDEX+1,1)         ! W/m.K
     &        +(HE_THERMAL(INDEX+1,0)-T_HE)*HE_THERMAL(INDEX,1))
     &       /(HE_THERMAL(INDEX+1,0)-HE_THERMAL(INDEX,0))
      D_HE = ((T_HE-HE_THERMAL(INDEX,0))*HE_THERMAL(INDEX+1,2)         ! kg/m^3
     &        +(HE_THERMAL(INDEX+1,0)-T_HE)*HE_THERMAL(INDEX,2))
     &       /(HE_THERMAL(INDEX+1,0)-HE_THERMAL(INDEX,0))
      MU_HE = ((T_HE-HE_THERMAL(INDEX,0))*HE_THERMAL(INDEX+1,3)        ! kg/m.s
     &         +(HE_THERMAL(INDEX+1,0)-T_HE)*HE_THERMAL(INDEX,3))
     &        /(HE_THERMAL(INDEX+1,0)-HE_THERMAL(INDEX,0))
C
C  TRISO fuel geometries, including calculation of fuel swelling
      V_PARTICLE = 4.0D0*PIE*R5**3*1.0D-18/3.0D0
      V_FUEL = 4.0D0*PIE*R10**3/3.0D0
      SWELL_FUEL = BURNUP*BU_CONV*1.0D-6
      V_FUEL = V_FUEL*(1.0D0 + SWELL_FUEL)
      R1 = (V_FUEL*3.0D0/4.0D0/PIE)**(1.0D0/3.0D0)     !Update R1
	V_FUEL = V_FUEL*1.0D-18    !Convert V_FUEL from um^3 to m^3
      IF (R1 .GE. R2) THEN
	  ERROR_CODE = 1
	  CALL ERR_HANDLER(TRIM('SUBROUTINE TEMPERATURE: R1 > R2, bad sample or
     &               wrong calculation'), 64, 1, 1, IERR)
	  RETURN
      END IF
      V_BUFFER = 4.0D0*PIE*(R2**3-R1**3)*1.0D-18/3.0D0
      V_IPYC = 4.0D0*PIE*(R3**3-R2**3)*1.0D-18/3.0D0
      V_SIC = 4.0D0*PIE*(R4**3-R3**3)*1.0D-18/3.0D0
      V_OPYC = 4.0D0*PIE*(R5**3-R4**3)*1.0D-18/3.0D0
C
C  Coating material thermal properties
C   Use temperatures in fuel, buffer, IPyC and OPyC from last step for approximation
C   Find K of fuel: UO2, UCO or else.
      T = (T_PARTICLE(0)+T_PARTICLE(1))/2.0D0
	IF (FUELTYPE .EQ. 'UO2') THEN
        CALL LOCATE_ARRAY(UO2_THERMAL(:,0), 16, 1, T, INDEX)
	  IF (INDEX .EQ. 0) INDEX = INDEX + 1   !exceeds lowerbound; make adjustion.
	  IF (INDEX .EQ. 16) INDEX = INDEX - 1  !exceeds upperbound; make adjustion.
	  K_FUEL = ((T-UO2_THERMAL(INDEX,0))*UO2_THERMAL(INDEX+1,1)      ! W/m.K
     &            +(UO2_THERMAL(INDEX+1,0)-T)*UO2_THERMAL(INDEX,1))
     &           /(UO2_THERMAL(INDEX+1,0)-UO2_THERMAL(INDEX,0))
	ELSE IF(FUELTYPE .EQ. 'UCO') THEN
        CALL LOCATE_ARRAY(UCO_THERMAL(:,0), 14, 1, T, INDEX)
	  IF (INDEX .EQ. 0) INDEX = INDEX + 1   !exceeds lowerbound; make adjustion.
	  IF (INDEX .EQ. 14) INDEX = INDEX - 1  !exceeds upperbound; make adjustion.
	  K_FUEL = ((T-UCO_THERMAL(INDEX,0))*UCO_THERMAL(INDEX+1,1)      ! W/m.K
     &            +(UCO_THERMAL(INDEX+1,0)-T)*UCO_THERMAL(INDEX,1))
     &           /(UCO_THERMAL(INDEX+1,0)-UCO_THERMAL(INDEX,0))        
      ELSE   !Temperorily use this value. Don't know K of other fuels yet.
	  K_FUEL = 2.30D0                                                ! W/m.K
	END IF
C   Find K of buffer
      T = (T_PARTICLE(1)+T_PARTICLE(2))/2.0D0
      CALL LOCATE_ARRAY(PYC_THERMAL(:,0), 6, 1, T, INDEX)
	IF (INDEX .EQ. 0) INDEX = INDEX + 1     !exceeds lowerbound; make adjustion.
      IF (INDEX .EQ. 6) INDEX = INDEX - 1     !exceeds upperbound; make adjustion.
C    The value 0.25 below comes from the assumption that K of buffer is
C    one quarter of that of dense PyC.
	K_BUFFER = 0.25D0*((T-PYC_THERMAL(INDEX,0))*PYC_THERMAL(INDEX+1,1) ! W/m.K
     &            +(PYC_THERMAL(INDEX+1,0)-T)*PYC_THERMAL(INDEX,1))
     &            /(PYC_THERMAL(INDEX+1,0)-PYC_THERMAL(INDEX,0))
C   Find K of IPyC
      T = (T_PARTICLE(2)+T_PARTICLE(3))/2.0D0
      CALL LOCATE_ARRAY(PYC_THERMAL(:,0), 6, 1, T, INDEX)
	IF (INDEX .EQ. 0) INDEX = INDEX + 1     !exceeds lowerbound; make adjustion.
      IF (INDEX .EQ. 6) INDEX = INDEX - 1     !exceeds upperbound; make adjustion.
      K_IPYC = ((T-PYC_THERMAL(INDEX,0))*PYC_THERMAL(INDEX+1,1)        ! W/m.K
     &          +(PYC_THERMAL(INDEX+1,0)-T)*PYC_THERMAL(INDEX,1))
     &         /(PYC_THERMAL(INDEX+1,0)-PYC_THERMAL(INDEX,0))
C   Find K of OPyC
      T = (T_PARTICLE(4)+T_PARTICLE(5))/2.0D0
      CALL LOCATE_ARRAY(PYC_THERMAL(:,0), 6, 1, T, INDEX)
      IF (INDEX .EQ. 0) INDEX = INDEX + 1     !exceeds lowerbound; make adjustion.
      IF (INDEX .EQ. 6) INDEX = INDEX - 1     !exceeds upperbound; make adjustion.
      K_OPYC = ((T-PYC_THERMAL(INDEX,0))*PYC_THERMAL(INDEX+1,1)        ! W/m.K
     &          +(PYC_THERMAL(INDEX+1,0)-T)*PYC_THERMAL(INDEX,1))
     &         /(PYC_THERMAL(INDEX+1,0)-PYC_THERMAL(INDEX,0))
C   Find K of SiC
      T = (T_PARTICLE(3)+T_PARTICLE(4))/2.0D0
      CALL LOCATE_ARRAY(SIC_THERMAL(:,0), 4, 1, T, INDEX)
	IF (INDEX .EQ. 0) INDEX = INDEX + 1     !exceeds lowerbound; make adjustion.
      IF (INDEX .EQ. 4) INDEX = INDEX - 1     !exceeds upperbound; make adjustion.
C    The value of 0.5 below comes from the assumption that K of irradiated
C    SiC is half of the original value.
      K_SIC = 0.5D0*((T-SIC_THERMAL(INDEX,0))*SIC_THERMAL(INDEX+1,1)   ! W/m.K
     &          +(SIC_THERMAL(INDEX+1,0)-T)*SIC_THERMAL(INDEX,1))
     &         /(SIC_THERMAL(INDEX+1,0)-SIC_THERMAL(INDEX,0))
C
C  Pebble model is in this block
C   Volume of pebble zones (dimension: m^3)
      V_PEBBLE = 4.0D0*PIE*PEBRADIUS**3/3.0D0
      V_PFZ = 4.0D0*PIE*PFZRADIUS**3/3.0D0
      V_PFM = V_PFZ - NPARTICLE*V_PARTICLE
	V_PM = 4.0D0*PIE*(PEBRADIUS**3 - PFZRADIUS**3)/3.0D0
      IF (V_PFM .LT. 0.0D0) THEN
	  ERROR_CODE = 2
	  CALL ERR_HANDLER(TRIM('SUBROUTINE TEMPERATURE: Input error on fuel 
     &               particle quantities'), 63, 2, 2, IERR)
	END IF
	A_CORE = PIE*CORE_RADIUS**2
      V_CORE = A_CORE*CORE_HEIGHT
	PACKING = NPEBBLE*V_PEBBLE/V_CORE
      Q_PFZ = (QPPP/PACKING)*V_PEBBLE/V_PFZ
C   Thermal conductivity of matrix graphite is from Kania and Nickel; temperatures are in C.
      T = (T_PEBBLE(0)+T_PEBBLE(1))/2.0D0
	K_PFM = 47.4D0*(1.0D0-9.7556D-4*(T-100.0D0)*DEXP(-6.036D-4*T))
      T = (T_PEBBLE(1)+T_PEBBLE(2))/2.0D0
	K_PM = 47.4D0*(1.0D0-9.7556D-4*(T-100.0D0)*DEXP(-6.036D-4*T))
C   Heat transfer coefficient of He is calculated using Achenbach Correlation.
      PEBDIAMETER = 2.0D0*PEBRADIUS
      VC_HE = MF_HE/(D_HE*A_CORE*(1.0D0-PACKING)) !Characteristic velocity of He
      RE_HE = D_HE*PEBDIAMETER*VC_HE/MU_HE        !Reynold's number of He
      PR_HE = MU_HE*CP_HE/K_HE                    !Prandtl number of He
	H_HE = (K_HE/PEBDIAMETER)*(PR_HE**(1.0D0/3.0D0))*               ! W/m^2.K
     &       (((1.18D0*RE_HE**0.58D0)**4.0D0+
     &         (0.23D0*RE_HE**0.75D0)**4.0D0)**0.25D0)
C   Average K of pebble fuel zone
      K_PFZ = ((V_FUEL*K_FUEL+V_BUFFER*K_BUFFER+V_IPYC*K_IPYC+
     &         V_SIC*K_SIC+V_OPYC*K_OPYC)*NPARTICLE+V_PFM*K_PFM)/V_PFZ
C   Update pebble temperature profile
      T_PEBBLE(0) = TPEBBLE(T_HE, H_HE, Q_PFZ, 0.0D0)
      T_PEBBLE(1) = TPEBBLE(T_HE, H_HE, Q_PFZ, PFZRADIUS)
      T_PEBBLE(2) = TPEBBLE(T_HE, H_HE, Q_PFZ, PEBRADIUS)
C  Calculate the temperature profile of a particle randomly chosen in the pebble
      R = R_IN_PEBBLE
	T_PS = TPEBBLE(T_HE, H_HE, Q_PFZ, R)
      Q_FUEL = Q_PFZ*V_PFZ/(NPARTICLE*V_FUEL)
      T_PARTICLE(0) = TPARTICLE(T_PS, Q_FUEL, 0.0D0)
	T_PARTICLE(1) = TPARTICLE(T_PS, Q_FUEL, R1)
	T_PARTICLE(2) = TPARTICLE(T_PS, Q_FUEL, R2)
	T_PARTICLE(3) = TPARTICLE(T_PS, Q_FUEL, R3)
	T_PARTICLE(4) = TPARTICLE(T_PS, Q_FUEL, R4)
	T_PARTICLE(5) = TPARTICLE(T_PS, Q_FUEL, R5)
C
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
C  Subroutine CURRENTDATE                                              *
C                                                                      *
C    The function CURRENTDATE returns the current date in the form     *
C  Month/Day/Year.                                                     *
C                                                                      *
C    It is assumed that the common block TYPE has been initialized by  *
C  calling subroutine MACHINE.                                         *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual variable description                                         *
C    CDATE          C: Date in format MM/DD/YY                         *
C                                                                      *
C  Local variable description                                          *
C    CDATE2         C: Date in format DD/MMM/YYYY                      *
C    COMP           C: Three character description of compiler type    *
C    C1             D: Numerical constant, 2 (2^31 - 1) = 4294967294   *
C    C2             D: Numerical constant, 86400, number of seconds    *
C                      in a day                                        *
C    IERR           I: Logical Unit Number of error file 'ERROR.MSG'   *
C    IMONTH         I: Current month                                   *
C    IDAY           I: Current day                                     *
C    IYEAR          I: Current two digit year                          *
C    IY             I: Accumulate sum of days per month for normal     *
C                      length year.                                    *
C    IYL            I: Accumulate sum of days per month for leap year. *
C    J              I: Intermediate variable                           *
C    LPYR           L: Logical, TRUE if a leap year, FALSE otherwise.  *
C    LY             I: = 1 if NYEAR is a leap year, = 0 otherwise      *
C    MACH           C: Three character description of computer type    *
C    N1             C: String index for the first occurence of         *
C                      separator '/'                                   *
C    N2             C: String index for the second occurence of        *
C                      separator '/'                                   *
C    NYEAR          I: Running year variable                           *
C    NTIME          I: Estimate of the number of days elapsed          *
C    PERIOD         D: Processed reference time                        *
C                                                                      *
C  None-standard features                                              *
C    LONG            : 32 bit absolute memory addressing function for  *
C                      Macintosh computers.  The memory address 524    *
C                      (decimal) has the time in seconds since midnight*
C                      January 1, 1904. LONG is an INTEGER function.   *
C                                                                      *
C    TIME()          : On the DEC workstations with DEC Ultrix Fortran,*
C                      the call to TIME returns the time in seconds    *
C                      since midnight January 1, 1970.  TIME is an     *
C                      INTEGER function.                               *
C                                                                      *
C    DATE            : On the IBM, DATE returns the date in the        *
C                      CHARACTER*11 form MM:DD:YY, where MM (month),   *
C                      DD (day), YY (year), are two digits each.       *
C                      If the result is not exactly 8 characters in    *
C                      length, the date is left-justified and truncated*
C                      or space-filled as necessary.                   *
C                                                                      *
C  Functions and subroutines called                                    *
C    LEAP            : Logical function that returns TRUE if a year    *
C                      is a leap year, returns FALSE otherwise.        *
C    LOCATE          : Indices for integer look-up table.              *
C    ERR_HANDLER     : Processes an error (diagnostic) message.        *
C                                                                      *
C  Intrinsic functions                                                 *
C    CHAR                                                              *
C    DBLE                                                              *
C    ICHAR                                                             *
C    IDINT                                                             *
C    INDEX                                                             *
C    MOD                                                               *
C                                                                      *
C  Common blocks                                                       *
C    MTYPE                                                             *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE CURRENTDATE (CDATE)
C
      DOUBLE PRECISION PERIOD, C1, C2
      INTEGER IMONTH, IDAY, IYEAR, NTIME, IY, IYL, NYEAR, J, LY, IERR
C
C#######################################################################
CDEC     INTEGER TIME
C#######################################################################
      CHARACTER*8 CDATE
      CHARACTER*11 CDATE2      ! original
      CHARACTER*3 MACH, COMP
      LOGICAL LPYR, LEAP
C
      PARAMETER (IERR = 12)
C
      DIMENSION IY(13), IYL(13)
C
      COMMON /MTYPE/ MACH, COMP
C
      EXTERNAL LEAP
C
      SAVE C1, C2, IY, IYL
C
	DATA CDATE2/'        '/
C
C  Set to zero to avoid compiler warnings
      DATA NTIME/0/
C
C  Initialize logical variables
      DATA LPYR/.FALSE./
C  2 (2^31 - 1)
C
C  Define numerical constants
      DATA C1/4 294 967 294.0 D+00/  
      DATA C2/86400.0 D+00/
C
C  Define date tables
      DATA IY/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304,
     A        334, 365/
      DATA IYL/0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305,
     A         335, 366/
C
C  Machine dependent calculations
      IF ((MACH .EQ. 'MAC') .OR. (MACH .EQ. 'DEC')) THEN
C
C  Generic to all Macintosh compilers
C
C  The memory address 524 (decimal) has the time in seconds since
C  midnight January 1, 1904.
C
C#######################################################################
CMAC      IF (MACH .EQ. 'MAC') NTIME = LONG(524) 
C#######################################################################
C
C  On the DEC workstations with DEC Ultrix Fortran, the call to TIME
C  returns the time in seconds since midnight January 1, 1970.
C
C#######################################################################
CDEC      IF (MACH .EQ. 'DEC') NTIME = TIME()
C#######################################################################
C
C  The number of seconds returned are in the range - 2^31 -> 2^31 - 1
C  Negative values are used to extend the range of the system clocks.
C  If a negative value for time is returned 2*(2^31 - 1) must be added
C  to get the correct time.  This would cause INTEGER*4 overflow.
C  To prevent this, the time is converted to double precision, then
C  scaled to the number of days.
C
        IF (NTIME .GE. 0) THEN
          PERIOD = DBLE(NTIME)
        ELSE
          PERIOD = C1 + DBLE(NTIME)
        END IF
C
C  Convert to the number of days (86400 seconds) since midnight
C  January 1, 1904 for Macintosh or midnight January 1, 1970
C
C  Number of days
        NTIME = IDINT(PERIOD/C2) 
C
C  Determine Year
C
C  The current year is determined by subtracting the number of days
C  corresponding to each year from the day count, NTIME, until the
C  day count is negative.  Then adjust the day count.
C
C  Starting year
        IF (MACH .EQ. 'MAC') NYEAR = 1904
        IF (MACH .EQ. 'DEC') NYEAR = 1970
3210    LPYR = LEAP(NYEAR)
        IF (LPYR) THEN
          LY = 1
        ELSE
          LY = 0
        END IF
        NTIME = NTIME - 365 - LY
        IF (NTIME .GE. 0) THEN
C  Go to next year
          NYEAR = NYEAR + 1 
          GO TO 3210
C  Adjust time for over count
        ELSE  
          NTIME = NTIME + 365 + LY
        END IF
C
C  Convert year to two digit form
        IYEAR = MOD(NYEAR,100)
C
C  NTIME represents the number of days completed.
C
C  The current day adds 1 to the number of days completed.
        NTIME = NTIME + 1 
C
C  Determine Month
        IF (LPYR) THEN
          CALL LOCATE (IYL, 13, NTIME, J)
        ELSE
          CALL LOCATE (IY, 13, NTIME, J)
        END IF
        IMONTH = J
C
C  Determine Day
C
        IF (LPYR) THEN
          NTIME = NTIME - IYL(J)
        ELSE
          NTIME = NTIME - IY(J)
        END IF
        IDAY = NTIME
      ELSE IF (MACH .EQ. 'IBM') THEN
C  Only for IBM compilers
C
C  IBM routine DATE_AND_TIME returns the date in the format 'YYYYMMDD '
C
C#######################################################################
        CALL DATE_AND_TIME (CDATE2)
C#######################################################################
C	
        IYEAR = 10*(ICHAR(CDATE2(3:3)) - 48) + ICHAR(CDATE2(4:4)) - 48
	  IMONTH = 10*(ICHAR(CDATE2(5:5)) - 48) + ICHAR(CDATE2(6:6)) - 48
	  IDAY = 10*(ICHAR(CDATE2(7:7)) - 48) + ICHAR(CDATE2(8:8)) - 48
	ELSE
	  CALL ERR_HANDLER(TRIM('SUBROUTINE CURRENTDATE: unsupported 
     &					machine type'), 42, 2, 2, IERR)
	END IF

C  This was the original month scheme
C      IDAY = 10*ICHAR(CDATE2(1:1)) + ICHAR(CDATE2(2:2))
C      IYEAR = 1000*ICHAR(CDATE2(8:8)) + 100*ICHAR(CDATE2(9:9)) +
C     A         10*ICHAR(CDATE2(10:10)) + ICHAR(CDATE2(11:11))
C        IF (CDATE(4:6) .EQ. 'JAN') THEN
C          IMONTH = 1
C        ELSE IF (CDATE(4:6) .EQ. 'FEB') THEN
C          IMONTH = 2
C        ELSE IF (CDATE(4:6) .EQ. 'MAR') THEN
C          IMONTH = 3
C        ELSE IF (CDATE(4:6) .EQ. 'APR') THEN
C          IMONTH = 4
C       ELSE IF (CDATE(4:6) .EQ. 'MAY') THEN
C          IMONTH = 5
C        ELSE IF (CDATE(4:6) .EQ. 'JUN') THEN
C          IMONTH = 6
C        ELSE IF (CDATE(4:6) .EQ. 'JUL') THEN
C          IMONTH = 7
C        ELSE IF (CDATE(4:6) .EQ. 'AUG') THEN
C          IMONTH = 8
C        ELSE IF (CDATE(4:6) .EQ. 'SEP') THEN
C          IMONTH = 9
C        ELSE IF (CDATE(4:6) .EQ. 'OCT') THEN
C         IMONTH = 10
C       ELSE IF (CDATE(4:6) .EQ. 'NOV') THEN
C         IMONTH = 11
C       ELSE IF (CDATE(4:6) .EQ. 'DEC') THEN
C         IMONTH = 12
C       END IF
C        CALL ERR_HANDLER('CURRENTDATE:  unsupported machine type',33,2,2,IERR)
C      END IF
C
C  Express date in MM/DD/YY form
C
C  ASCII number for 0 is 48, ... ASCII number for 9 is 57
C
C  Tens of months
      CDATE(1:1) = CHAR(48 + IMONTH/10) 
C  Single months
      CDATE(2:2) = CHAR(48 + MOD(IMONTH,10))
C  Break
      CDATE(3:3) = '/' 
C  Tens of days
      CDATE(4:4) = CHAR(48 + IDAY/10)
C  Single days
      CDATE(5:5) = CHAR(48 + MOD(IDAY,10))
C  Break
      CDATE(6:6) = '/' 
C  Tens of years
      CDATE(7:7) = CHAR(48 + IYEAR/10) 
C  Single years
      CDATE(8:8) = CHAR(48 + MOD(IYEAR,10))
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
C  Subroutine ERR_HANDLER                                              *
C                                                                      *
C    Processes a diagnostic message, in a manner determined            *
C  by the value of LEVEL.                                              *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual argument description                                         *
C    MESSG          C: The Hollerith message to be processed,          *
C                      containing no more than 72 characters           *
C    LENGTH         I: The actual number of characters in MESSG.       *
C    NERR           I: The error number associated with this message   *
C    LEVEL          I: Error category                                  *
C                      2 means this is an unconditionally fatal error. *
C                      1 means this is a recoverable error.            *
C                      0 means this is a warning message only.         *
C    LUN            I: Unit number of error message file 'ERROR.MSG'   *
C                                                                      *
C***********************************************************************
C                                                                      *
C      From a user's point of view there are only a few important      *
C  things to know about ERR_HANDLER. When a library routine detects an *
C  error situation ERR_HANDLER will be called.  This call includes a   *
C  message which describes the difficulty, an error number and a       *
C  seriousness level.  The ERR_HANDLER package prints the message,     *
C  records the error number in a table for later summary.              *
C  By convention, the error message contains in its first few          *
C  characters the name of the routine causing the error.               *
C                                                                      *
C      The seriousness level can be FATAL, RECOVERABLE, or WARNING.    *
C  A FATAL error causes a program abort and in this situation          *
C  ERR_HANDLER never returns to the subroutine which called it.        *
C  Thus the math library subroutine writer who calls ERR_HANDLER with  *
C  a FATAL error does not need to deal with what to do next.  For      *
C  example, not dimensioning an array large enough is a typical        *
C  FATAL error. A RECOVERABLE error is something that the user might   *
C  be able to correct. For example, a data value out of range. By      *
C  declaring an error RECOVERABLE, the library subroutine writer must  *
C  anticipate that ERR_HANDLER will return and must code for the       *
C  possibility that the computation should be able to proceed.  A      *
C  WARNING error is just like a RECOVERABLE error except that only     *
C  unusual action by the user will prevent ERR_HANDLER from returning  *
C  to the library routine. In practice most errors are declared either *
C  FATAL or WARNING.                                                   *
C                                                                      *
C  Note:  IOSTAT error code returned on READ error.  Other errors      *
C           are self-explanatory (and just return an error 0).         *
C                                                                      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE ERR_HANDLER(MESSG,LENGTH,NERR,LEVEL,LUN)
      INTEGER*4 LENGTH, NERR, LEVEL, LUN
      CHARACTER MESSG*256
      LOGICAL EX, FOPEN
C
C  The error messages are written to file 'ERROR.MSG'.
C
      INQUIRE (FILE = 'error.msg', EXIST = EX)
	IF (EX) THEN
	  INQUIRE (FILE = 'error.msg', OPENED = FOPEN)
	  IF (.NOT. FOPEN) THEN
	    OPEN(UNIT = LUN,FILE = 'error.msg',STATUS = 'OLD',
     A         ACCESS = 'SEQUENTIAL')
	  END IF
      ELSE
	  OPEN(UNIT = LUN,FILE = 'error.msg',STATUS = 'NEW',
     A     ACCESS = 'SEQUENTIAL')	  
      END IF
      WRITE(LUN,3301) MESSG(1:LENGTH)
      IF (LEVEL .EQ. 2) THEN
	  WRITE(*,3301) MESSG(1:LENGTH)
        WRITE(LUN,3302)
	  WRITE(*,3302)
        WRITE(LUN,3305) NERR
        CLOSE(UNIT = LUN)
        STOP                          ! Terminate the program
      ELSE IF (LEVEL .EQ. 1) THEN
        WRITE(LUN,3303)
        WRITE(LUN,3305) NERR
	  WRITE(LUN,*)
      ELSE IF (LEVEL .EQ. 0) THEN
        WRITE(LUN,3304)
        WRITE(LUN,3305) NERR
	  WRITE(LUN,*)
      END IF
      RETURN
C
C  Format block:
3301  FORMAT(' ',A)
3302  FORMAT(' ','FATAL error, program terminating')
3303  FORMAT(' ','RECOVERABLE error, program continuing')
3304  FORMAT(' ','WARNING error, program continuing')
3305  FORMAT(' ','Error type:  ',I5)
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
!C  Function ELAPTIME                                                   *
!C                                                                      *
!C    Elapsed time since T0 in seconds                                  *
!C                                                                      *
!C  Notation                                                            *
!C    D               : REAL*8 (Double Precision)                       *
!C    I               : INTEGER                                         *
!C    L               : LOGICAL                                         *
!C    C               : CHARACTER*n                                     *
!C                                                                      *
!C  Actual variable description                                         *
!C    T0 (sec)       D: Reference time                                  *
!C    ELAPTIME (sec) D: Elapsed time since T0                           *
!C                                                                      *
!C  Functions and subroutines called                                    *
!C    TIMEON          : Time in seconds that the machine has been on    *
!C                      or since midnight depending on machine type.    *
!C***********************************************************************
!C
!C***********************************************************************
!C                                                                      *
!      FUNCTION ELAPTIME (T0)
!      DOUBLE PRECISION ELAPTIME, TIMEON, T0
!
!      EXTERNAL TIMEON
!
!      ELAPTIME = TIMEON(T0) - T0
!      RETURN
!      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C  Function  FLUX(U235E, QPPP, T_PARTICLE)                             *
C                                                                      *
C    This function calculates the fast neutron flux based on power     *
C  density and temperature.                                            *
C                                                                      *
C  Actual argument description                                         *
C    U235E (%)      D: U235 enrickment                                 *
C    QPPP (W/m^3)   D: Current power density                           *
C    T_PARTICLE(0:5) (C) --                                            *
C                   D: Temperature distribution in the particles       *
C    FLUX (n/cm^2)  D: Fast neutron flux                               *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      FUNCTION FLUX(U235E, QPPP, T_PARTICLE)
C
      DOUBLE PRECISION FLUX, U235E, QPPP, T_PARTICLE
      DIMENSION T_PARTICLE(0:5)
C
      FLUX = 1.17 D+14
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
C  Function GAUSSIAN                                                   *
C                                                                      *
C    Gaussian distributed kernel and coating thickness                 *
C                                                                      *
C    Adaptation of GASDEV modified to accept mean and standard         *
C  deviation.                                                          *
C    Function GASDEV comes from W. Press, Numerical Recipes, Cambridge,*
C  Cambridge Univ. Press, 1986.                                        *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual variable description                                         *
C    MEAN           D: Mean of distribution                            *
C    STDDEV         D: Standard deviation of distribution              *
C    GAUSSIAN       D: Returned value                                  *
C                                                                      *
C  Local variables description                                         *
C    X1             D: Numerical constant 1                            *
C    X2             D: Numerical constant 2                            *
C    V1, V2         D: Intermediate variables                          *
C    R              D: Intermediate variable                           *
C    FACTOR         D: Intermediate variable                           *
C                                                                      *
C  Functions and subroutines called                                    *
C    RAND            : Random number generator on the interval (0,1].  *
C                      Assumed initialized by calling routine.         *
C                                                                      *
C  Intrinsic functions called                                          *
C    DLOG                                                              *
C    DSQRT                                                             *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      FUNCTION GAUSSIAN(MEAN, STDDEV)
C
      DOUBLE PRECISION GAUSSIAN, MEAN, STDDEV             ! Global variables
      DOUBLE PRECISION V1, V2, R, FACTOR, X1, X2     ! Local variables
      INTEGER ISET
C
C      EXTERNAL RAND
C
      SAVE ISET, V1, FACTOR, X1, X2
C
C  Local flag to use second deviate
      DATA ISET /0/
	DATA X1/1.0 D+00/, X2/2.0 D+00/
C
      IF (ISET .EQ. 0) THEN
C  Cosine and Sin distributions
4110    CALL RANDOM_NUMBER(RN)    
        V1 = X2*RN - X1
        CALL RANDOM_NUMBER(RN)
        V2 = X2*RN - X1
        R = V1**2 + V2**2
        IF (R .GE. X1) GO TO 4110
        FACTOR = DSQRT(- X2*DLOG(R)/R)
        GAUSSIAN = MEAN + STDDEV*(V2*FACTOR)
        ISET = 1
      ELSE
C  Uses saved variate
        GAUSSIAN = MEAN + STDDEV*(V1*FACTOR)
        ISET = 0
      ENDIF
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
C  Function INITSEED                                                   *
C                                                                      *
C    Initial random number seed based on system clock                  *
C                                                                      *
C  Actual variable description                                         *
C    INITSEED       I: Initial random number seed                      *
C    IDUMMY         I: Dummy variable                                  *
C                                                                      *
C  Subroutines and functions called                                    *
C    TIMEON          : Time in seconds that the machine has been on,   *
C                      or since midnight depending on machine type.    *
C  Intrinsic functions                                                 *
C    INT                                                               *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      FUNCTION INITSEED (IDUMMY)
C
!      DOUBLE PRECISION TIMEON
      INTEGER INITSEED, IDUMMY
C
!      EXTERNAL TIMEON
C
      INITSEED = 5
C
C  The random number generator seed must be non-zero.  If INITSEED = 0,
C  set to documented default seed 305.
      IF (INITSEED .EQ. 0) INITSEED = 305
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
C  Function LEAP                                                       *
C                                                                      *
C    Logical function returns TRUE if a year is a leap year, returns   *
C  FALSE otherwise.                                                    *
C                                                                      *
C  Actual variable description                                         *
C    LEAP           L: Logical result of leap year test                *
C    NYEAR          I: Year to be tested                               *
C                                                                      *
C  Intrinsic functions                                                 *
C    MOD                                                               *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      FUNCTION LEAP (NYEAR)
C
      INTEGER NYEAR
      LOGICAL LEAP
C
C  A year is a leap year if it is evenly divisible by 4.  However, if
C  the year is evenly divisible by 100 it is NOT a leap year, unless
C  it is also evenly divisible by 400, in which case it IS.
C
      IF ((MOD(NYEAR, 4) .EQ. 0) .AND. (((MOD(NYEAR, 100) .NE. 0) .AND.
     A    (MOD(NYEAR, 400) .NE. 0)) .OR.
     B    (MOD(NYEAR, 400). EQ. 0))) THEN
        LEAP = .TRUE.
      ELSE
        LEAP = .FALSE.
      END IF
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
C  Function NDAYS (CDATE1,CDATE2)                                      *
C                                                                      *
C    Returns the number of days between CDATE1 and CDATE2              *
C                                                                      *
C  Actual variable description                                         *
C    CDATE1          C: Initial date                                   *
C    CDATE2          C: Finial date                                    *
C                                                                      *
C  Local variable description                                          *
C    I               I: Tally variable                                 *
C    IY              I: Accumulate sum of days per month for normal    *
C                       length year.                                   *
C    IYL             I: Accumulate sum of days per month for leap year.*
C    IMON1           I: Month corresponding to CDATE1                  *
C    IMON2           I: Month corresponding to CDATE2                  *
C    IDAY1           I: Day corresponding to CDATE1                    *
C    IDAY2           I: Day corresponding to CDATE2                    *
C    IYR1            I: Two digit year corresponding to CDATE1         *
C    IYR2            I: Two digit year corresponding to CDATE2         *
C    J               I: Tally and intermediate variable                *
C    LY              L: Logical result of leap year test               *
C    NDAY1           I: Number of days since 1900 for the initial date *
C    NDAY2           I: Number of days since 1900 for the final date   *
C                                                                      *
C  Subroutines and functions called                                    *
C    LEAP             : Logical function returns TRUE if a year is a   *
C                       leap year, returns FALSE otherwise             *
C    LOCATE           : Indicies for look-up table                     *
C                                                                      *
C  Intrinsic functions called                                          *
C    ICHAR                                                             *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      FUNCTION NDAYS (CDATE1,CDATE2)
C
      INTEGER NDAYS, NDAY1, NDAY2, IDAY1, IMON1, IYR1
      INTEGER IDAY2, IMON2, IYR2, IY, IYL, I, J
      CHARACTER*8 CDATE1, CDATE2
      LOGICAL LY, LEAP
C
      DIMENSION IY(13), IYL(13)
C
      EXTERNAL LEAP
C
      SAVE IY, IYL
C
C  Define date tables
C
      DATA IY/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304,
     A        334, 365/
      DATA IYL/0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305,
     A         335, 366/
C
C  Convert each date to the number of days since the last century,
C  this is a fudicial mark.
C
C  Date CDATE1
C
C  Convert character form of date into numerical form useful for
C  calculations.
C
      IDAY1 = 10*ICHAR(CDATE1(4:4)) + ICHAR(CDATE1(5:5))
      IMON1 = 10*ICHAR(CDATE1(1:1)) + ICHAR(CDATE1(2:2))
      IYR1 = 1900 + 10*ICHAR(CDATE1(7:7)) + ICHAR(CDATE1(8:8))
C
C  Accumulate the number of days corresponding to elapse years
C  since 1900.
C
      NDAY1 = 0
      DO 4210 I = 1900, IYR1 - 1
      LY = LEAP(I)
      IF (LY) THEN
        NDAY1 = NDAY1 + 366
      ELSE
        NDAY1 = NDAY1 + 365
      END IF
4210  CONTINUE
C
C  Leap year adjustment and the number of days of the input year
C
      LY = LEAP(IYR1)
      IF (LY) THEN
        CALL LOCATE(IYL,13,IMON1,J)
        NDAY1 = NDAY1 + IYL(J)
      ELSE
        CALL LOCATE(IY,13,IMON1,J)
        NDAY1 = NDAY1 + IY(J)
      END IF
C
C  Total number of days
C
      NDAY1 = NDAY1 + IDAY1
C
C  Convert each date to the number of days since the last century,
C  this is a fudicial mark.
C
C  Date CDATE2
C
C  Convert character form of date into numerical form useful for
C  calculations.
C
      IDAY2 = 10*ICHAR(CDATE2(4:4)) + ICHAR(CDATE2(5:5))
      IMON2 = 10*ICHAR(CDATE2(1:1)) + ICHAR(CDATE2(2:2))
      IYR2 = 1900 + 10*ICHAR(CDATE2(7:7)) + ICHAR(CDATE2(8:8))
C
C  Accumulate the number of days corresponding to elapse years
C  since 1900.
C
      NDAY2 = 0
      DO 4220 I = 1900, IYR2 - 1
      LY = LEAP(I)
      IF (LY) THEN
        NDAY2 = NDAY2 + 366
      ELSE
        NDAY2 = NDAY2 + 365
      END IF
4220  CONTINUE
C
C  Leap year adjustment and the number of days of the input year
C
      LY = LEAP(IYR2)
      IF (LY) THEN
        CALL LOCATE(IYL,13,IMON2,J)
        NDAY2 = NDAY2 + IYL(J)
      ELSE
        CALL LOCATE(IY,13,IMON2,J)
        NDAY2 = NDAY2 + IY(J)
      END IF
C
C  Total number of days
C
      NDAY2 = NDAY2 + IDAY2
C
C  Difference in days
C
      NDAYS = NDAY2 - NDAY1
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
!C  Function RAN                                                        *
!C                                                                      *
!C  This routine generates quasi uniform random numbers on the          *
!C  interval (0,1] with density (2^31 - 1) using the lagged Fibonacci   *
!C  sequence.                                                           *
!C                                                                      *
!C      Xn = [X(n-17) - X(n - 5)] Mod[m]                                *
!C                                                                      *
!C  where                                                               *
!C                                                                      *
!C      m = 2**(MDIG - 2) + (2**(MDIG - 2) - 1)                         *
!C      MDIG = number of binary digit available for representing        *
!C             integers, including the sign bit.                        *
!C                                                                      *
!C***********************************************************************
!C                                                                      *
!C  Author    Blue, James, Scientific Computing Division, NBS           *
!C            Kahaner, Davis, Scientific Computing Division, NBS        *
!C            Marsaglia, George, Computer Science Department,           *
!C            Washington State university.                              *
!C                                                                      *
!C  Modified by Lorenz H. Menke, Jr., (208) 526-8918 EG&G Idaho, Inc.,  *
!C  PO Box 1625, Idaho Falls, ID 83415-1575.                            *
!C                                                                      *
!C  Purpose   This routine generates two quasi uniform random numbers   *
!C            on the interval (0,1].  It can be used with any computer  *
!C            which allows integers at least as large as 32767.         *
!C                                                                      *
!C  Description                                                         *
!C                                                                      *
!C  This routine generates quasi uniform random numbers on the interval *
!C  (0,1].  It can be used with any computer which allows integers at   *
!C  least as large as 32767.                                            *
!C                                                                      *
!C   Use:                                                               *
!C                                                                      *
!C       First time ...                                                 *
!C                                                                      *
!C       Z = RAN (JD)                                                   *
!C                                                                      *
!C       Here JD is any  N O N - Z E R O  integer.  This causes         *
!C       initialization of the program and the first two random         *
!C       numbers to be returned as Z1 and Z2.                           *
!C                                                                      *
!C       Subsequent times...                                            *
!C                                                                      *
!C       Z = RAN (0)                                                    *
!C                                                                      *
!C       Causes the next random number to be returned as Z.             *
!C                                                                      *
!C***********************************************************************
!C                                                                      *
!C  Note:   Users who wish to transport this program from one computer  *
!C          to another should read the following information...         *
!C                                                                      *
!C  Machine dependencies...                                             *
!C                                                                      *
!C  MDIG    =  A lower bound on the number of binary digit available    *
!C             for representing integers, including the sign bit.       *
!C             This value must be at least 16, but may be increased     *
!C             in line with remark A below.                             *
!C                                                                      *
!C  Remarks...                                                          *
!C                                                                      *
!C      A. This program can be used in two ways:                        *
!C                                                                      *
!C          (1) To obtain repeatable results on different computers     *
!C              set 'MDIG' to the smallest of its values on each, or,   *
!C          (2) to allow the longest sequence of random numbers to be   *
!C              generated without cycling (repeating) set 'MDIG' to the *
!C              largest possible value.                                 *
!C                                                                      *
!C      B. The sequence of numbers depends on the initial input 'JD'    *
!C         as well as the value of 'MDIG'.                              *
!C                                                                      *
!C             If MDIG = 16 one should find that                        *
!C             The first evaluation                                     *
!C                 Z = UNI(305,Z1,Z2) gives Z1 = .027832881...          *
!C                 and Z2 = .56102176...                                *
!C             The second evaluation                                    *
!C                 Z = UNI(0,Z1,Z2) gives   Z1 = .41456343...           *
!C                 and Z2 = .19797357...                                *
!C                                                                      *
!C***********************************************************************
!C                                                                      *
!C  Actual variable description                                         *
!C                                                                      *
!C  Z1        :  First random number.                                   *
!C  Z2        :  Second random number.                                  *
!C  JD        :  Random number generator seed.  JD is used to           *
!C               initialize the random number generator.                *
!C  MDIG      :  A lower bound on the number of binary digit available  *
!C               for representing integers, including the sign bit.     *
!C  M         :  Array containing starting sequence.                    *
!C                                                                      *
!C  Local variables                                                     *
!C                                                                      *
!C  JSEED     :  Modified seed.                                         *
!C  I,J       :  Loop index.                                            *
!C  K         :  Differences between M array elements.                  *
!C  M1        :  Modulus.                                               *
!C  M2,K0,K1, :  intermediate values.                                   *
!C  J0,J1                                                               *
!C  R0        :  Numerical constant 0                                   *
!C  R1        :  Numerical constant 1                                   *
!C  Z0        :  1/M1, scales the random umber to (0,1]                 *
!C                                                                      *
!C***********************************************************************
!C                                                                      *
!C  Functions and subroutines called:                                   *
!C                                                                      *
!C  ERR_HANDLER :  Processes an error (diagnostic) message.             *
!C                                                                      *
!C  Intrinsic functions:                                                *
!C                                                                      *
!C  DBLE                                                                *
!C  IABS                                                                *
!C  MIN0                                                                *
!C  MOD                                                                 *
!C                                                                      *
!C***********************************************************************
!C
!C***********************************************************************
!C                                                                      *
!      FUNCTION RAN (JD)
!C
!      DOUBLE PRECISION RAN, Z0, Z1, Z2, R0, R1
!      INTEGER I, J, K, M, M1, M2, JD, MDIG, K0, K1, J0, J1, JSEED
!C
!      PARAMETER (IERR = 12)
!      DIMENSION M(17)
!C
!      SAVE I, J, M, M1, M2
!C
!	DATA M(1), M(2), M(3), M(4), M(5), M(6), M(7), M(8), M(9), M(10),
!     1     M(11), M(12), M(13), M(14), M(15), M(16), M(17)
!     2  / 30788, 23052, 2053, 19346, 10646, 19427, 23975, 19049, 10949,
!     3    19693, 29746, 26748, 2796, 23890, 29168, 31924, 16499 /
!      DATA M1, M2, I, J / 32767, 256, 5, 17 /
!C
!C  32 bit constant; 1/(2^31 - 1) = 4.6566 12875 24579 6923 X 10**-10
!C
!C      DATA Z0/4.6566 12875 24579 6923 D-10
!	DATA R0/0.0 D+00/, R1/1.0 D+00/
!C  
!C  First executable statement  RAN.
!C
!C  Check if random number generator has been initialized, if not
!C  initialize it.
!C
!      IF (JD .NE. 0) THEN
!C
!C***********************************************************************
!C
!C  Machine dependent variable, MDIG
!C
!C  Largest integer supported by 32 bit machines such as Macintosh,
!C  IBM/IBM compatible PCs and DEC 5000 class workstations is:
!C
!C      2**31 - 1 = 2,147,483,647
!C
!C      (2^31 - 1)  -> 31 + 1 (sign) = 32 = MDIG
!C
!        MDIG = 32
!C
!C***********************************************************************
!C
!C  Be sure that MDIG is at least 16.  If not, a fatal error condition
!C
!        IF (MDIG .LT. 16) CALL
!     AERR_HANDLER(TRIM('FUNCTION RAN --MDIG less than 16'),32,1,2,IERR)
!        M1 = 2**(MDIG - 2) + (2**(MDIG - 2) - 1)
!        M2 = 2**(MDIG/2)
!        Z0 = (1.0D+00)/DBLE(M1)
!	  JSEED = MIN0(IABS(JD),M1)
!C
!C  Check if JSEED is odd.
!C
!        IF (MOD(JSEED,2) .EQ. 0) JSEED = JSEED - 1
!        K0 = MOD(9069,M2)
!        K1 = 9069/M2
!        J0 = MOD(JSEED,M2)
!        J1 = JSEED/M2
!        DO 4310 I = 1, 17
!          JSEED = J0*K0
!          J1 = MOD(JSEED/M2 + J0*K1 + J1*K0, M2/2)
!          J0 = MOD(JSEED,M2)
!          M(I) = J0 + M2*J1
!4310      CONTINUE
!        I = 5
!        J = 17
!      END IF
!C
!C  Begin main loop here.
!C
!C  First random number, NOT normalized
!C
!4320  K = M(I) - M(J)
!      IF (K .LT. 0) K = K + M1
!      M(J) = K
!      I = I - 1
!      IF (I .EQ. 0) I = 17
!      J = J - 1
!      IF (J .EQ. 0) J = 17
!C     Z1 = DBLE(K)*Z3   ! original
!	Z1 = DBLE(K)
!C
!C  Second random number, normalized
!C      
!      K = M(I) - M(J)
!      IF (K .LT. 0) K = K + M1
!      M(J) = K
!      I = I - 1
!      IF (I .EQ. 0) I = 17
!      J = J - 1
!      IF (J .EQ. 0) J = 17
!C      Z2 = DBLE(K)*Z3    ! original
!	Z2 = DBLE(K)*Z0
!C
!C  Combination random number generator
!C
!C  Extend lower range to ~ Z0^2 = .....
!C
!	RAN = (Z1 - Z2)*Z0
!	IF ((RAN .GE. R1) .OR. (RAN .LE. R0)) GOTO 4320
!      RETURN
!      END
C                                                                      *
C***********************************************************************
C
C
C***********************************************************************
C                                                                      *
!C  Function TIMEON (TIME)                                              *
!C                                                                      *
!C    The function TIMEON returns the number of seconds elapse since    *
!C  the machine were turned on.  The starting reference point is        *
!C  arbitrary.  On Macintosh machines the starting point is since the   *
!C  machine was turned on.  For IBM and DEC machines the starting point *
!C  is the time since midnight.  Date changes are tracked for run times *
!C  that run through midnight.  The TIMEON function is initialized by   *
!C  setting the argument TIME to 0.0                                    *
!C                                                                      *
!C    It is assumed that the common block TYPE has been initialized by  *
!C  calling subroutine MACHINE.                                         *
!C                                                                      *
!C  Actual variable description                                         *
!C    TIME0 (sec)    D: Used to initialize the calendar for TIMEON      *
!C    TIMEON (sec)   D: Number of seconds elapse since the machine was  *
!C                      turned on, or since midnight.                   *
!C                                                                      *
!C  Local variable description                                          *
!C    CDATE          C: Initialization date                             *
!C    CDATE2         C: Date of most recent call                        *
!C    COMP           C: Three character description of compiler type    *
!C    CTIME          C: Initialization time                             *
!C    DTIME (sec)    D: Number of seconds since midnight of the         *
!C                      initialization date                             *
!C    IERR           I: Logical Unit Number of error file 'ERROR.MSG'  *
!C    MACH           C: Three character description of computer type    *
!C    NDAY           I: Number of elapsed days since initialization call*
!C    TADJ (sec)     D: Correction to TIMEON that accounts for date     *
!C                      changes                                         *
!C                                                                      *
!C  None-standard features                                              *
!C    LONG            : 32 bit absolute memory addressing function for  *
!C                      Macintosh computers.  The memory address 362    *
!C                      (decimal) has the number of 'tics',             *
!C                      (60 tics = 1 second), elapse since the machine  *
!C                      was turned on.  LONG is an INTEGER function.    *
!C                                                                      *
!C    TIMER           : Returns the number of 'tics', (100 tics = 1     *
!C                      second), elapse since midnight on an IBM.       *
!C                      TIMER is an INTEGER function.                   *
!C                                                                      *
!C    SECNDS          : Returns the number of seconds since midnight    *
!C                      on a DEC. SECNDS is a REAL function.            *
!C                                                                      *
!C  Subroutines and functions called                                    *
!C    CURRENTTIME     : Returns the current time in the format HH:MM:SS.*
!C    NDAYS           : Returns the number of elapse days between two   *
!C                      dates.                                          *
!C    CURRENTDATE     : Returns the current date in the format MM/DD/YY.*
!C    ERR_HANDLER     : Processes an error (diagnostic) message.        *
!C                                                                      *
!C  Intrinsic functions                                                 *
!C    FLOAT                                                             *
!C                                                                      *
!C  Common blocks                                                       *
!C    MTYPE                                                             *
!C***********************************************************************
!C
!C***********************************************************************
!C                                                                      *
!      FUNCTION TIMEON (TIME0)
!C
!      DOUBLE PRECISION TIMEON, TIME0, DTIME, DTIME2
!      INTEGER NTICKS, NDAY, NDAYS, IERR
!      CHARACTER*3 MACH, COMP
!      CHARACTER*8 CDATE, CDATE2, CTIME, CTIME2
!C
!      PARAMETER (IERR = 12)
!C
!      COMMON /MTYPE/ MACH, COMP
!C
!      EXTERNAL NDAYS
!C
!      SAVE NDAY, CDATE
!C
!	DATA CDATE2/'        '/
!C  Set to zero to avoid compiler warnings
!C
!      DATA NTICKS/0/
!C
!C  The Dec and IBM timing routines are based on the number of seconds
!C  elapsed since midnight.  Therefore, tracking the number of day changes
!C  is necessary.
!C
!C  Accumulation of elapse days since initialization call
!C
!      IF (TIME0 .EQ. 0.0D0) THEN
!        CALL CURRENTDATE(CDATE)
!	  CALL CURRENTTIME(CTIME)
!        NDAY = 0
!C
!C  Convert to the number of seconds since midnight of the
!C  initialization call.
!C
!	  DTIME = FLOAT(
!     A          3600*(10*ICHAR(CTIME(1:1)) + ICHAR(CTIME(2:2))) +
!     B            60*(10*ICHAR(CTIME(4:4)) + ICHAR(CTIME(5:5))) +
!     C               (10*ICHAR(CTIME(7:7)) + ICHAR(CTIME(8:8))) )
!      TIMEON = DTIME
!      RETURN
!	END IF
!C
!C  Find current date and calulate the number of date changes since
!C  the last call.  Then accumulate the total number of date changes
!C  since the initialization call.
!C
!      CALL CURRENTDATE(CDATE2)
!      CALL CURRENTTIME(CTIME2)
!      NDAY = NDAY + NDAYS(CDATE,CDATE2)
!      CDATE = CDATE2
!	DTIME2 = FLOAT(
!     A         3600*(10*ICHAR(CTIME2(1:1)) + ICHAR(CTIME2(2:2))) +
!     B           60*(10*ICHAR(CTIME2(4:4)) + ICHAR(CTIME2(5:5))) +
!     C              (10*ICHAR(CTIME2(7:7)) + ICHAR(CTIME2(8:8))) )
!C
!C  Machine specific elapse time calculations
!C
!      IF (MACH .EQ. 'MAC') THEN
!C
!C  Generic to all Macintosh compilers
!C
!C  The memory address 362 (decimal) has the time in seconds elapse
!C  since the machine was turned on based on the 60 Hz counter.
!C
!C#######################################################################
!CMAC        NTICKS = LONG(362)
!C#######################################################################
!        TIMEON = FLOAT(NTICKS)/60.0
!      ELSE IF (MACH .EQ. 'DEC') THEN
!C
!C  The function SECONDS returns the number of seconds since midnight.
!C  Account of the number of days for timing periods that are longer
!C  than one day or run through midnight are made.
!C
!        TIMEON = FLOAT(NTICKS)/60.0
!C#######################################################################
!CDEC        TIMEON = SECNDS(0.0) + TADJ
!C#######################################################################
!      ELSE IF (MACH .EQ. 'IBM') THEN
!C
!        TIMEON = FLOAT(86400*NDAY) + DTIME2
!      ELSE
!        CALL ERR_HANDLER(TRIM('FUNCTION TIMEON: unsupported machine type
!     &              '), 41, 2, 2, IERR)
!      END IF
!      RETURN
!      END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C  Function TPARTICLE(T_PS, Q_FUEL, R)                                 *
C                                                                      *
C    Attention!!! R in this function has the dimension of "um".        *
C                                                                      *
C    Calcuate the temperature (C) at radius R in a particle, given the *
C  particle surface temperature, average power density in fuel kernel. *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual argument description                                         *
C    T_PS (C)       D: Current particle surface temperature            *
C    Q_FUEL (W/m^3) D: Volumetric heat generation rate in fuel kernel  *
C    R (um)         D: Radius                                          *
C                                                                      *
C  Variables in COMMON blocks                                          *
C  /PAR_R/           : Layer radii of TRISO particles                  *
C    R1 - R5 (um)   D: Various radii of layers in TRISO particles      *
C  /PAR_T/           : Thermal conductivity of particle materials      *
C    K_FUEL (W/m.K) D: Thermal conductivity of fuel kernel             *
C    K_BUFFER (W/m.K) D:                                               *
C                      Thermal conductivity of pyrocarbon buffer       *
C    K_IPYC (W/m.K) D: Thermal conductivity of IPyC layer              *
C    K_SIC (W/m.K)  D: Thermal conductivity of SiC layer               *
C    K_OPYC (W/m.K) D: Thermal conductivity of OPyC layer              *
C                                                                      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      FUNCTION TPARTICLE(T_PS, Q_FUEL, R)
C
      DOUBLE PRECISION TPARTICLE, T_PS, Q_FUEL, R
	DOUBLE PRECISION R1, R2, R3, R4, R5
      DOUBLE PRECISION K_FUEL, K_BUFFER, K_IPYC, K_SIC, K_OPYC
C
      COMMON /PAR_R/ R1, R2, R3, R4, R5               !particle geometry
	COMMON /PAR_T/ K_FUEL, K_BUFFER, K_IPYC, K_SIC, K_OPYC
C
      IF ((R.GE.0.0D0).AND.(R.LT.R1)) THEN
	  TPARTICLE = T_PS + 1.0D-12*(Q_FUEL*R1**3.0/3.0D0)*
     &                     ((1.0D0/R4-1.0D0/R5)/K_OPYC + 
     &                      (1.0D0/R3-1.0D0/R4)/K_SIC + 
     &                      (1.0D0/R2-1.0D0/R3)/K_IPYC + 
     &                      (1.0D0/R1-1.0D0/R2)/K_BUFFER + 
     &                      (1.0D0/R1-R**2.0/R1**3.0)/(2.0D0*K_FUEL))
	ELSE IF ((R.GE.R1).AND.(R.LT.R2)) THEN
	  TPARTICLE = T_PS + 1.0D-12*(Q_FUEL*R1**3.0/3.0D0)*
     &                     ((1.0D0/R4-1.0D0/R5)/K_OPYC + 
     &                      (1.0D0/R3-1.0D0/R4)/K_SIC + 
     &                      (1.0D0/R2-1.0D0/R3)/K_IPYC + 
     &                      (1.0D0/R-1.0D0/R2)/K_BUFFER)
	ELSE IF ((R.GE.R2).AND.(R.LT.R3)) THEN
	  TPARTICLE = T_PS + 1.0D-12*(Q_FUEL*R1**3.0/3.0D0)*
     &                     ((1.0D0/R4-1.0D0/R5)/K_OPYC + 
     &                      (1.0D0/R3-1.0D0/R4)/K_SIC + 
     &                      (1.0D0/R-1.0D0/R3)/K_IPYC)
	ELSE IF ((R.GE.R3).AND.(R.LT.R4)) THEN
	  TPARTICLE = T_PS + 1.0D-12*(Q_FUEL*R1**3.0/3.0D0)*
     &                     ((1.0D0/R4-1.0D0/R5)/K_OPYC + 
     &                      (1.0D0/R-1.0D0/R4)/K_SIC)
      ELSE IF ((R.GE.R4).AND.(R.LT.R5)) THEN
	  TPARTICLE = T_PS + 1.0D-12*(Q_FUEL*R1**3.0/3.0D0)*
     &                     (1.0D0/R-1.0D0/R5)/K_OPYC
      ELSE
	  TPARTICLE = T_PS
	END IF
	RETURN
	END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C  Function TPEBBLE(T_HE, H_HE, Q_PFZ, R)                              *
C                                                                      *
C    Calcuate the temperature (C) at radius R in a pebble, given the   *
C  He temperature, He heat transfer coefficient and power density in   *
C  the pebble fuel zone.                                               *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual argument description                                         *
C    T_HE (C)       D: Current Helium temperature                      *
C    H_HE (W/m^2.K) D: Heat transfer coefficient of He                 *
C    Q_PFZ (W/m^3)  D: Volumetric heat generation rate in the pebble   *
C                      fuel zone                                       *
C    R (m)          D: Radius                                          *
C                                                                      *
C  Variables in COMMON blocks                                          *
C  /PEBBLE/          : Pebble geometry and material thermal properties *
C    PEBRADIUS (m)  D: Pebble radius                                   *
C    PFZRADIUS (m)  D: Pebble fuel zone radius                         *
C    K_PM (W/m.K)   D: Thermal conductivity of pebble graphite matrix  *
C    K_PFM (W/m.K)  D: Thermal conductivity of matrix in pebble fuel   *
C                      zone                                            *
C    K_PFZ (W/m.K)  D: Thermal conductivity of pebble fuel zone        *
C    NPARTICLE      I: Number of TRISO particles in each pebble        *
C                                                                      *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      FUNCTION TPEBBLE(T_HE, H_HE, Q_PFZ, R)
C
      DOUBLE PRECISION TPEBBLE, T_HE, H_HE, Q_PFZ, R
	DOUBLE PRECISION PEBRADIUS, PFZRADIUS, R_IN_PEBBLE,
     &                 K_PM, K_PFM, K_PFZ
	INTEGER*4 NPARTICLE
C
      COMMON /PEBBLE/ PEBRADIUS, PFZRADIUS, K_PM, K_PFM, K_PFZ,
     &                R_IN_PEBBLE, NPARTICLE
C
      IF ((R.GE.0.0D0).AND.(R.LT.PFZRADIUS)) THEN
	  TPEBBLE = T_HE + (Q_PFZ*PFZRADIUS**3.0/3.0D0)*
     &                   (1.0D0/H_HE/PEBRADIUS**2.0 + 
     &                    (1.0D0/PFZRADIUS-1.0D0/PEBRADIUS)/K_PM + 
     &                    (1.0D0/PFZRADIUS-R**2.0/PFZRADIUS**3.0)/
     &                    (2.0D0*K_PFZ))
	ELSE IF ((R.GE.PFZRADIUS).AND.(R.LE.PEBRADIUS)) THEN
	  TPEBBLE = T_HE + (Q_PFZ*PFZRADIUS**3.0/3.0D0)*
     &                   (1.0D0/H_HE/PEBRADIUS**2.0 + 
     &                    (1.0D0/R-1.0D0/PEBRADIUS)/K_PM)
      ELSE
	  TPEBBLE = T_HE
	END IF
      RETURN
	END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C  Function TRIANGLE                                                   *
C                                                                      *
C    Triangular distributed critical stress intensity factor of SiC    *
C                                                                      *
C  Actual argument description                                         *
C    MEAN           D: Mean of distribution.                           *
C    STDDEV         D: Standard deviation of distribution.             *
C    TRIANGLE       D: Returned value.                                 *
C                                                                      *
C  Local variable description                                          *
C    S              D: parameter related to standard deviation         *
C    F              D: number (0:1) sampled by uniform distribution    *
C                                                                      *
C  Functions and subroutines called                                    *
C    RANDOM_NUMBER     Random number generator on the interval (0,1]   *
C                      Assumed initialized by calling routine.         *
C                                                                      *
C  Intrinsic functions called                                          *
C    SQRT                                                              *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      FUNCTION TRIANGLE(MEAN, STDDEV)
C
	DOUBLE PRECISION TRIANGLE, MEAN, STDDEV, S, F
C
C	EXTERNAL RAN
C
      S = STDDEV*SQRT(6.0)
        CALL RANDOM_NUMBER(RN)
	F = RN
	IF(F.LE.0.5) THEN
	  TRIANGLE = MEAN - S + S*SQRT(2.0*F)
	ELSE
	  TRIANGLE = MEAN + S - S*SQRT(2.0*(1.0-F))
	ENDIF
	RETURN
	END
C                                                                      *
C***********************************************************************
C
C
C
C***********************************************************************
C                                                                      *
C  Function WEIBULL                                                    *
C                                                                      *
C    Weibull distributed SiC layer strength                            *
C                                                                      *
C    Samples Weibull distribution, based on cumulative distribution    *
C  function.                                                           *
C    Cumulative Weibull comes from Bourne and Green, Reliability       *
C  Technology, 1974.                                                   *
C                                                                      *
C  Actual argument description                                         *
C    MEAN           D: Mean of distribution                            *
C    MODULUS        D: Modulus of distribution                         *  
C    WEIBULL        D: Weibull distribution                            *
C                                                                      *
C  Local variable description                                          *
C    F              D: Random variant                                  *
C                                                                      *
C  Functions and subroutines called                                    *
C    RAN            D: Random number generator on the interval (0,1].  *
C                      Assumed initialized by calling routine.         *
C                                                                      *
C  Intrinsic functions called                                          *
C    DLOG                                                              *
C    DEXP                                                              *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      FUNCTION WEIBULL(MEAN, MODULUS)
C
      DOUBLE PRECISION WEIBULL, MEAN, MODULUS, F, RAN
C
C      EXTERNAL RAN
C  F is the Weibull cumulative distribution function  
        CALL RANDOM_NUMBER(RN)
        F = RN
C
C  The distributions 1 - F and F are equivalent
C      WEIBULL = MEAN*DEXP(DLOG(- DLOG(1.0D0 - F))/MODULUS)    ! Original
	WEIBULL = MEAN*DEXP(DLOG(- DLOG(F))/MODULUS)
      RETURN
      END
C                                                                      *
C***********************************************************************
C
C