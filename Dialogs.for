C***********************************************************************
C                                                                      *
C  Subroutine Greeting                                                 *
C                                                                      *
C  Displays the welcome screen.                                        *
C                                                                      *
C  Author: Jing Wang, 01/27/2003                                       *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual variable description                                         *
C    None                                                              *
C  Local variable description                                          *
C    None                                                              *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE Greeting()
C
	USE DFLOGM
	IMPLICIT NONE
	INCLUDE 'resource.fd'
C
	INTEGER retint
	LOGICAL retlog
	TYPE (dialog) dlg_mpbr
C  Display greeting window
	retlog = DlgInit( IDD_DIALOG_GREETING, dlg_mpbr )
	retint = DlgModal( dlg_mpbr )
	Call DlgUninit( dlg_mpbr )
C
	END
C                                                                      *
C***********************************************************************
C
C
C***********************************************************************
C                                                                      *
C  Subroutine InputDialog                                              *
C                                                                      *
C  Displays the input file request dialog and get input file name.     *
C                                                                      *
C  Author: Jing Wang, 01/27/2003                                       *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual variable description                                         *
C    InputFilename  C: The returned input file name.                   *
C  Local variable description                                          *
C    None                                                              *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE InputDialog(InputFilename)
C
	USE DFLOGM
	IMPLICIT NONE
	INCLUDE 'resource.fd'
C
	CHARACTER*256 InputFilename
	INTEGER retint
	LOGICAL retlog
	TYPE (dialog) dlg_mpbr
C  Reactor core and TRISO fuel input data dialog boxes
C    Initialize the dialog box
	retlog = DlgInit(IDD_DIALOG_INPUT, dlg_mpbr)
	retlog = DlgSet( dlg_mpbr, IDC_EDIT_MPBRINPUT, "PBMR" )
!    Activate the modal dialog.
	retint = DlgModal( dlg_mpbr )
	retlog = DlgGet( dlg_mpbr, IDC_EDIT_MPBRINPUT, InputFilename)
C Release dialog resources	
	CALL DlgUninit( dlg_mpbr )
C
	END
C                                                                      *
C***********************************************************************
C
C
C***********************************************************************
C                                                                      *
C  Subroutine PowerTypeDialog                                          *
C                                                                      *
C  Displays the power type selection dialog and get power type.        *
C                                                                      *
C  Author: Jing Wang, 01/27/2003                                       *
C                                                                      *
C  Notation                                                            *
C    D               : REAL*8 (Double Precision)                       *
C    I               : INTEGER                                         *
C    L               : LOGICAL                                         *
C    C               : CHARACTER*n                                     *
C                                                                      *
C  Actual variable description                                         *
C    PowerType      I: The returned power type.                        *
C  Local variable description                                          *
C    None                                                              *
C***********************************************************************
C
C***********************************************************************
C                                                                      *
      SUBROUTINE PowerTypeDialog(PowerType)
C
	USE DFLOGM
	IMPLICIT NONE
	INCLUDE 'resource.fd'
C
	INTEGER PowerType
	INTEGER retint
	LOGICAL retlog
	LOGICAL PushedState
	TYPE (dialog) dlg_ptype
C    Initialize the dialog box
	retlog = DlgInit(IDD_DIALOG_POWERTYPE, dlg_ptype)
	retlog = DlgSet( dlg_ptype, IDC_RADIO_POWER1, .TRUE. )
!    Activate the modal dialog.
	retint = DlgModal( dlg_ptype )
	retlog = DlgGet( dlg_ptype, IDC_RADIO_POWER1, PushedState )
	IF (.NOT. PushedState) THEN
	  retlog = DlgGet( dlg_ptype, IDC_RADIO_POWER2, PushedState )
	  IF (.NOT. PushedState) THEN
	    PowerType = 3
	  ELSE
	    PowerType = 2
	  END IF
	ELSE
	  PowerType = 1
	END IF
C Release dialog resources	
	CALL DlgUninit( dlg_ptype )
C
	END
C                                                                      *
C***********************************************************************