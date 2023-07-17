! F90 program HITBIN
!
! VERSION (update VERSID)
!   23NOV18 AD v3.02 Fix Bug using unassigned WNO value.
!                    Reorder HITTYP data type into 8-byte boundaries for ifort.
!   03JUL18 AD v3.01 Fix Bug with assigning CO2 isotope#10,11,12
!   01JAN18 AD v3.00 F90 conversion.
!
! DESCRIPTION
!   Convert HITRAN 160 character .par file to RFM direct access binary file
!   Apart from .par file, requires additional file vibh16.txt containing
!   assignments of indices to vibrational levels.
!
!   The structure of this file has the HITBIN module before the main program
!   in order to simplify the compilation (eg gfortran hitbin.f90 -o hitbin)
!   MODULE HITBIN
!     VARIABLE KINDS
!     GLOBAL DATA
!     CONTAINS
!       SUBROUTINE REAVIB    this reads vibh16.txt
!       FUNCTION IDXVIB      this assigns vib level indices
!   END MODULE
!   PROGRAM HITBIN
!
!   Variables for each transition record
!              Input  Output Description
!      LSTAT     -      I 4   Priority of transition information.
!      MOL      I2      I 8   Molecule number.
!      ISO      I1      I 12  Isotope number (=1 most abundant, 2=second etc)
!      WNO     F12.6    D 20  Line frequency [cm-1].
!      STR     D10.3    R 24  Line strength  I=[cm-1./(molec.cm-2)] @ 296K.
!                                            O=[cm-1./(kg.moles.cm-2)] @296K
!                           (Scale input by Avogadro No.to avoid underflows.)
!      TPR     E10.3    R 28  Transition probab. [Debyes2] (old) or Einstein A
!      ABR      F5.4    R 32  Air-broadened halfwidth  (HWHM) [cm-1/atm] @ 296K.
!      SBR      F5.3    R 36  Self-broadened halfwidth (HWHM) [cm-1/atm] @ 296K.
!      ELS     F10.4    R 40  Lower-state energy [cm-1].
!      ABC      F4.2    R 44  Coefficient of temperature dependance of ABROAD
!      TSP      F8.6    R 48  Transition shift due to pressure
!      IUV       A15    I 52  Upper state global (=Vib) quanta index / ID
!      ILV       A15    I 56  Lower state global quanta index / ID
!      IUR       A15   A9 65  Upper state local (=Rot) quanta.
!      ILR       A15   A9 74  Lower state local quanta.
!      SPARE9          A9 83  Spare (used to be quality and reference info)
!      IFP       -      I 87  Forward pointer on data line.
!      SPARE1          A1 88  Spare (used to be empty byte at end of record)
!
MODULE HITBIN_MOD
!
! VARIABLE KINDS
    INTEGER, PARAMETER :: I4 = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: R4 = SELECTED_REAL_KIND(6)
    INTEGER, PARAMETER :: R8 = SELECTED_REAL_KIND(15,200)
!
! GLOBAL CONSTANTS
    INTEGER(I4),  PARAMETER :: IFORM  = 1  ! format version# for HITRAN output
    INTEGER(I4),  PARAMETER :: LSTEND = -2 ! LSTAT value for file termination
    INTEGER(I4),  PARAMETER :: LSTFWD = -7 ! LSTAT value for Fwd Ptr record
    INTEGER(I4),  PARAMETER :: LSTHDR =  4 ! LSTAT value for header records
    INTEGER(I4),  PARAMETER :: LSTINI =  0 ! LSTAT value for file header record
    INTEGER(I4),  PARAMETER :: LSTLIN = 10 ! LSTAT value for line param. record

    INTEGER(I4),  PARAMETER :: LUNBIN = 2  ! LUN for new binary (output) file
    INTEGER(I4),  PARAMETER :: LUNHIT = 1  ! LUN for ASCII (input) file
    INTEGER(I4),  PARAMETER :: LUNVIB = 1  ! LUN for reading in Vib.levels
    INTEGER(I4),  PARAMETER :: MAXMOL = 50 ! Max number of HITRAN molecules
! Don't change MAXNFP since value is assumed in RFM as well
    INTEGER(I4),  PARAMETER :: MAXNFP = 14 ! Max no.of fwd pointers per record
    INTEGER(I4),  PARAMETER :: NFPREC = 1 + MAXMOL/MAXNFP ! No.recs in f.p.block
    INTEGER(I4),  PARAMETER :: NRECFP = 200           ! No recs between fp blks
    REAL(R4),     PARAMETER :: AVOG   = 6.02214199    ! Avogadro's No. *1e-26
    CHARACTER(*), PARAMETER :: VERSID = '3.02'        ! Program Version
    CHARACTER(*), PARAMETER :: VIBFIL = 'vibh16.txt'  ! File with Vib IDs
!
! GLOBAL TYPES
  TYPE :: CLSTYP           ! List of vib.level IDs for each class of molecules
    INTEGER(I4)            :: NVB     ! No. vibrational levels in class
    INTEGER(I4),   POINTER :: IDX(:)  ! Index for each Vib level
    CHARACTER(15), POINTER :: VIB(:)  ! HITRAN ID for each Vib level
  END TYPE CLSTYP
!
  TYPE :: HITTYP           ! Contents of output binary records
    INTEGER(I4)  :: LST    ! Type of binary file record
    INTEGER(I4)  :: MOL    ! Molecule ID
    INTEGER(I4)  :: ISO    ! Isotope ID
    INTEGER(I4)  :: IFP    ! Forward pointer read from binary file
    REAL(R8)     :: WNO    ! Wavenumber
    REAL(R4)     :: STR    ! Line strength
    REAL(R4)     :: TPR    ! Transition Probability
    REAL(R4)     :: ABR    ! Air Broadened halfwidth
    REAL(R4)     :: SBR    ! Self broadened halfwidth
    REAL(R4)     :: ELS    ! Energy of the lower state
    REAL(R4)     :: ABC    ! Air Broad. Temp. Coeff
    REAL(R4)     :: TSP    ! Pressure shift
    INTEGER(I4)  :: IUV    ! Upper level vibrational ID (old format HITRAN)
    INTEGER(I4)  :: ILV    ! Index of Vib.Temp profile affecting lower level
    CHARACTER(9) :: IUR    ! Upper level rotational ID (old format HITRAN)
    CHARACTER(9) :: ILR    ! Lower level rotational ID (old format HITRAN)
    CHARACTER(9) :: SPARE9 = '         ' ! Spare bytes in output record
    CHARACTER(1) :: SPARE1 = ' '         ! Spare byte in output record
  END TYPE HITTYP
!
! GLOBAL VARIABLES
    INTEGER(I4) :: ICLMOL(MAXMOL) = 0    ! Class# for each molecule
    TYPE(CLSTYP), ALLOCATABLE :: CLS(:)  ! Different classes of molecules
!
CONTAINS
SUBROUTINE REAVIB
!
! DESCRIPTION
!   Read in Vib level data from file into CLS structure
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    INTEGER(I4)   :: ICLS ! Counter for different molecule classes
    INTEGER(I4)   :: IDUM ! Dummy variable
    INTEGER(I4)   :: IMOL ! Counter for molecules within each class
    INTEGER(I4)   :: IVIB ! Counter for different vib levels
    INTEGER(I4)   :: MOL  ! Molecule index
    INTEGER(I4)   :: NCLS ! No. different molecule classes in file
    INTEGER(I4)   :: NMOL ! No. different molecules in each class
    INTEGER(I4)   :: NVIB ! No. different vibrational levels in class
    CHARACTER(80) :: REC  ! Record read from vib.data file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  OPEN ( UNIT=LUNVIB, FILE=VIBFIL, STATUS='OLD', ACTION='READ' )
  DO
    READ ( LUNVIB, '(A)' ) REC
    IF ( REC(1:1) .NE. '!' ) EXIT
  END DO
  READ ( REC, * ) NCLS  ! No. different molecule classes
  ALLOCATE ( CLS(NCLS) )
!
  DO ICLS = 1, NCLS
    DO
      READ ( LUNVIB, '(A)' ) REC
      IF ( REC(1:1) .NE. '!' ) EXIT
    END DO
    READ ( REC, * ) IDUM, NMOL, NVIB
    CLS(ICLS)%NVB = NVIB
    ALLOCATE ( CLS(ICLS)%VIB(NVIB), CLS(ICLS)%IDX(NVIB) )
    DO IMOL = 1, NMOL
      READ ( LUNVIB, * ) MOL
      ICLMOL(MOL) = ICLS
    END DO
    DO IVIB = 1, NVIB
      READ ( LUNVIB, '(A)' ) REC
      DO WHILE ( REC(1:1) .EQ. '!' )
        READ ( LUNVIB, '(A)' ) REC
      END DO
      READ ( REC, '(X,A15,X,I8)') CLS(ICLS)%VIB(IVIB), CLS(ICLS)%IDX(IVIB)
    END DO
  END DO
  CLOSE ( LUNVIB )
!
END SUBROUTINE REAVIB

INTEGER(I4) FUNCTION IDXVIB ( IDXMOL, STRVIB )
!
! DESCRIPTION
!   Translate HITRAN C*15 Global Quantum ID to integer index
!
  IMPLICIT NONE
!
! ARGUMENTS
    INTEGER(I4),   INTENT(IN) :: IDXMOL ! HITRAN index of molecule
    CHARACTER(15), INTENT(IN) :: STRVIB ! HITRAN vib.level ID
!
! LOCAL VARIABLES
    INTEGER(I4)   :: ICLS   ! Index of molecular class
    INTEGER(I4)   :: IVIB   ! Index of vib level within class
    INTEGER(I4)   :: NVIB   ! No. different vib levels for class
    CHARACTER(15) :: STRLOC ! Local version of STRVIB
    INTEGER(I4),   POINTER :: IDXSAV(:)  ! List of indices of vib levels
    CHARACTER(15), POINTER :: VIBSAV(:)  ! HITRAN codes for vib levels
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  ICLS = ICLMOL(IDXMOL)
!
! HITRAN2016 C2N2 (#48) has some vib level IDs misaligned
  IF ( IDXMOL .EQ. 48 ) THEN
    STRLOC = ADJUSTR ( STRVIB )       ! Right justify all C2N2 IDs
  ELSE
    STRLOC = STRVIB
  END IF
!
  IF ( STRVIB .EQ. '' ) THEN  ! Empty vib lev string, assign as IDXVIB=0
!  WRITE (*,*) 'W-IDXVIB: Setting IDXVIB=0 for missing vib lev, IDXMOL=',IDXMOL
   CONTINUE
  ELSE IF ( ICLS .LE. 0 ) THEN     ! 0 or -1
    IF ( ICLS .EQ. 0 ) THEN
      WRITE (*,*) 'W-IDXVIB: Setting IDXVIB=0 for unidentified IDXMOL=', IDXMOL
      ICLMOL(IDXMOL) = -1     ! -1 to flag that warning has been issued
    END IF
  ELSE
    DO IVIB = 1, CLS(ICLS)%NVB
      IF ( STRLOC .EQ. CLS(ICLS)%VIB(IVIB) ) THEN
        IDXVIB = CLS(ICLS)%IDX(IVIB)
        RETURN                ! Normal exit with identified vib level
      END IF
    END DO
    WRITE ( *, * ) 'W-IDXVIB: Unrecognised Vibrational Quanta=''' // STRVIB &
      // ''' for IDXMOL=', IDXMOL
    IDXSAV => CLS(ICLS)%IDX
    VIBSAV => CLS(ICLS)%VIB
    NVIB = CLS(ICLS)%NVB + 1
    ALLOCATE ( CLS(ICLS)%VIB(NVIB), CLS(ICLS)%IDX(NVIB) )
    CLS(ICLS)%NVB = NVIB
    CLS(ICLS)%VIB(1:NVIB-1) = VIBSAV
    CLS(ICLS)%VIB(NVIB) = STRLOC
    CLS(ICLS)%IDX(1:NVIB-1) = IDXSAV
    CLS(ICLS)%IDX(NVIB) = 0
    DEALLOCATE ( VIBSAV, IDXSAV )
  END IF
!
  IDXVIB = 0
!
END FUNCTION IDXVIB

END MODULE HITBIN_MOD

PROGRAM HITBIN
!
! DESCRIPTION
!   Main program.
!
  USE HITBIN_MOD
!
  IMPLICIT NONE
!
! LOCAL VARIABLES
    LOGICAL     :: USEVIB          ! T=Use vib level indices, F=set to 0
    INTEGER(I4) :: IFPREC          ! Record# within forward pointer block
    INTEGER(I4) :: IMOL            ! Molecule counter
    INTEGER(I4) :: IOFF            ! Offset for molecule# in f.p. records
    INTEGER(I4) :: IOSVAL          ! Saved value of IOSTAT
    INTEGER(I4) :: IREC            ! Record# in new binary file
    INTEGER(I4) :: IREC1           ! First record# in new file after headers
    INTEGER(I4) :: IREC2           ! Last record# in new file (termination rec)
    INTEGER(I4) :: IRECFP          ! First record# in new file after last fp.blk
    INTEGER(I4) :: JREC            ! Secondary record counter
    INTEGER(I4) :: LSTAT           ! Output record type
    INTEGER(I4) :: MOL             ! Molecule ID#
    INTEGER(I4) :: RECLEN          ! Record length of binary files
    INTEGER(I4) :: IRCMOL(NFPREC*MAXNFP)  ! Next record# for each molecule
    REAL(R8)    :: DSTR            ! Line strength allowing for < 1.0E-38
    REAL(R8)    :: WNO             ! Wavenumber
    REAL(R8)    :: WNOUPP          ! Last wavenumber written/read
    REAL(R8)    :: WNORQ1 = -1.0D0 ! Lowest wavenumber to select from .par file
    REAL(R8)    :: WNORQ2 = 1.0D6  ! Highest wavenumber to select from .par file
    CHARACTER(200) :: FILNAM       ! User-input name of input and output files
    CHARACTER(48)  :: HEAD48       ! User header for binary file
    CHARACTER(84)  :: HEADR2       ! HITBIN header record for binary file
    CHARACTER(15)  :: JLR          ! Lower level rotational ID
    CHARACTER(15)  :: JLV          ! Lower level vibrational ID
    CHARACTER(15)  :: JUR          ! Upper level rotational ID
    CHARACTER(15)  :: JUV          ! Upper level vibrational ID
    CHARACTER(80)  :: RECORD       ! User input for wavenumber range
!
  TYPE(HITTYP) :: REC  ! Output record for binary file
!
! EXECUTABLE CODE -------------------------------------------------------------
!
  WRITE ( *, '(A,I2)' ) 'R-HITBIN: Running HITBIN v' // VERSID // &
    ' generating output Format#', IFORM
!
! Open file containing vibrational level indices
  INQUIRE ( FILE=VIBFIL, EXIST=USEVIB )
  IF ( USEVIB ) THEN
    WRITE (*,*) 'I-HITBIN: Reading in Vib level indices from ' // VIBFIL
    CALL REAVIB
  ELSE
    WRITE (*,*) 'W-HITBIN: Vib index file not found: ' // VIBFIL
  END IF

! Get HITRAN (ASCII) input file and open
  WRITE ( *, '(A)', ADVANCE='NO' ) 'Input HITRAN .par file: '
  READ ( *, '(A)' ) FILNAM
  OPEN ( UNIT=LUNHIT, FILE=FILNAM, STATUS='OLD', ACTION='READ' )
!
! Select wavenumber range for ASCII file (default = use all)
  WRITE ( *, '(A)', ADVANCE='NO' ) 'Wavenumber range (cm-1) [<CR>=all]: '
  READ ( *, '(A)' ) RECORD
  IF ( RECORD .NE. '' ) THEN
    READ ( RECORD, * ) WNORQ1, WNORQ2
    WRITE ( *, '(A)' ) 'I-HITBIN: Finding first record of ASCII file...'
    DO
      READ ( LUNHIT, '(3X,F12.6)', END=900 ) WNO
      IF ( WNO .GE. WNORQ1 ) EXIT
    END DO
    IF ( WNO .GT. WNORQ2 ) GOTO 900           ! no overlap with reqd range
    BACKSPACE ( LUNHIT )
  END IF
!
! Get name of the new binary file to be created
  WRITE ( *, '(A)', ADVANCE='NO' ) 'New binary file: '
  READ ( *, '(A)' ) FILNAM
!
! Get record length for new binary file
  INQUIRE ( IOLENGTH=RECLEN ) REC
  WRITE ( *, * ) 'Opening file with RECL=', RECLEN
  OPEN ( UNIT=LUNBIN, FILE=FILNAM, STATUS='UNKNOWN', ACCESS='DIRECT', &
         ACTION='READWRITE', RECL=RECLEN )
!
! Get header for new file
  WRITE ( *, '(A)', ADVANCE='NO' ) 'Header for new file (up to 48 chars): '
  READ ( *, '(A)' ) HEAD48
!
! Write HITBIN header (record length and version ID) to rec#2
  WRITE ( HEADR2, '(A,I3,A,A)' ) 'RECL=', RECLEN, &
    ' Converted to binary format by HITBIN v.', VERSID
  WRITE ( LUNBIN, REC=2 ) LSTHDR, HEADR2
!
  IREC = 3
  IREC1 = IREC
!
  WRITE ( *, '(A)' ) 'I-HITBIN: Writing new binary file...'
! Begin output data with an empty forward pointer block
  DO JREC = 1, NFPREC
    WRITE ( LUNBIN, REC=IREC ) LSTFWD
    IREC = IREC + 1
  END DO
  IRECFP = IREC
!
! Repeat for each record, taking lowest wavenumber
  DO
!
    IF ( IREC - IRECFP .EQ. NRECFP ) THEN  ! insert new forward pointer block
      DO JREC = 1, NFPREC
        WRITE ( LUNBIN, REC=IREC ) LSTFWD
        IREC = IREC + 1
      END DO
      IRECFP = IREC
    END IF
!
! Copy record from ASCII file to new file
! Note that format descriptors Fw.d (eg F5.2) do not generally represent the
! actual HITRAN format (eg F5.3) but simply set to avoid warning messages from
! some compilers (eg ifort) which recommend w ge d+3
    READ ( LUNHIT, 1001, ERR=902, IOSTAT=IOSVAL, END=100 ) &
      REC%MOL, REC%ISO, REC%WNO, DSTR, REC%TPR, REC%ABR, REC%SBR, REC%ELS, &
      REC%ABC, REC%TSP, JUV, JLV, JUR, JLR
1001 FORMAT ( I2, Z1, F12.6, F10.3, E10.3, F5.2, F5.2, F10.4, F4.1, F8.5, 4A15 )
    IF ( REC%WNO .GT. WNORQ2 ) EXIT
    WNO = REC%WNO
! HITRAN uses Isotope Id '0','A','B' to represent #10,11,12, so adjust
    SELECT CASE ( REC%ISO )
    CASE ( 0 )     ; REC%ISO = 10
    CASE ( 10:15 ) ; REC%ISO = REC%ISO + 1
    CASE DEFAULT   ;
    END SELECT
!
! Convert new format C*15 codes for upper,lower vibrational level into old
! format integer(I4) indices required by RFM if particular Vib levels have to be
! identified (user-selected or non-LTE).
    IF ( USEVIB ) THEN
      REC%IUV = IDXVIB ( REC%MOL, JUV )
      REC%ILV = IDXVIB ( REC%MOL, JLV )
    ELSE
      REC%IUV = 0
      REC%ILV = 0
    END IF
!
! Reformat local quantum numbers for CO2 lines (required by RFM for line-mixing)
! The first 9 characters of the HITRAN04 C*15 strings seems to match the
! old C*9 fields except that in the old format the characters are shifted one
! space left and the new format has an extra character 'e' or 'f' which needs
! to be set blank (this was fixed in v2.21 of this code)
! No idea what this will do for local quantum numbers for other molecules (!)
! but hopefully this information won't be required
    REC%IUR = JUR(2:9)//' '
    REC%ILR = JLR(2:9)//' '
! To avoid underflow problems, STR read as double precision and scaled by
! Avogadro's number before converting to single precision
    REC%STR = SNGL ( DSTR * AVOG * 1.0E26 )
!
    REC%LST = LSTLIN
    REC%IFP = 0
    WRITE ( LUNBIN, REC=IREC ) REC%LST, REC%MOL, REC%ISO, REC%WNO, REC%STR, &
      REC%TPR, REC%ABR, REC%SBR, REC%ELS, REC%ABC, REC%TSP, REC%IUV, REC%ILV, &
      REC%IUR, REC%ILR, REC%SPARE9, REC%IFP, REC%SPARE1
    IREC = IREC + 1
    WNOUPP = WNO
!
    IF ( MOD ( IREC, 100000 ) .EQ. 0 ) WRITE ( *, '(A,I8,A,F12.6)' ) &
      'I-HITBIN: Record#', IREC, ' Wavenumber=', WNOUPP
!
  END DO
100 CONTINUE
!
! Save record after last data record
  IREC2 = IREC
!
! Write header record
  WRITE ( LUNBIN, REC=1 ) LSTINI, IFORM, IREC1, IREC2, HEAD48
!
! Write termination record
  WRITE ( LUNBIN, REC=IREC2 ) LSTEND, 0, 0, WNOUPP, (0,IMOL=1,14)
  WRITE ( *, '(A,I8,A,F12.6)' ) &
    'I-HITBIN: Last Record#', IREC2, ' Wavenumber=', WNOUPP
!
! Calculate forward pointers
  WRITE ( *, '(A)' ) 'I-HITBIN: Writing forward pointers...'
  IRCMOL(:) = IREC2
  IFPREC = NFPREC
  DO IREC = IREC2-1, IREC1, -1
    READ ( LUNBIN, REC=IREC ) LSTAT
    IF ( LSTAT .EQ. LSTFWD ) THEN
      IOFF = ( IFPREC - 1 ) * MAXNFP
      WRITE ( LUNBIN, REC=IREC ) LSTFWD, IOFF+1, 0, WNOUPP, &
        ( IRCMOL(IMOL)-IREC, IMOL = IOFF+1, IOFF+MAXNFP )
      IF ( IFPREC .EQ. 1 ) THEN
        IFPREC = NFPREC
      ELSE
        IFPREC = IFPREC - 1
      END IF
    ELSE
      READ ( LUNBIN, REC=IREC ) REC%LST, REC%MOL, REC%ISO, REC%WNO, REC%STR, &
      REC%TPR, REC%ABR, REC%SBR, REC%ELS, REC%ABC, REC%TSP, REC%IUV, REC%ILV, &
      REC%IUR, REC%ILR, REC%SPARE9, REC%IFP, REC%SPARE1
      MOL = REC%MOL
      WNOUPP = REC%WNO
      REC%IFP = IRCMOL(MOL) - IREC
      IRCMOL(MOL) = IREC
      WRITE ( LUNBIN, REC=IREC ) REC%LST, REC%MOL, REC%ISO, REC%WNO, REC%STR, &
      REC%TPR, REC%ABR, REC%SBR, REC%ELS, REC%ABC, REC%TSP, REC%IUV, REC%ILV, &
      REC%IUR, REC%ILR, REC%SPARE9, REC%IFP, REC%SPARE1
    END IF
    IF ( MOD ( IREC, 100000 ) .EQ. 0 ) WRITE ( *, '(A,I8,A,F12.6)' ) &
      'I-HITBIN: Record#', IREC, ' Wavenumber=', WNOUPP
  END DO
!
  STOP  'R-HITBIN: Successful completion'
!
900  CONTINUE
  STOP 'F-HITBIN: HITRAN file does not overlap required wno.range'
!
902  CONTINUE
  WRITE ( *, * ) &
    'F-HITBIN: Error decoding record from HITRAN file. IOSTAT=',IOSVAL
  WRITE ( *, * ) 'F-HITBIN: Last Wno read was ', WNO
  STOP
!
END PROGRAM HITBIN
