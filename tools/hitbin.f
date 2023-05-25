      PROGRAM HITBIN
C
C VERSION (update VERSID)
C     03FEB14 AD 2.21 Fix bug with copying CO2 rotational ID
C     10JAN14 AD 2.20 Allow for additional non-HITRAN molecules >#47
C     24OCT13 AD 2.10 Correction: Change ISO=0 from input to 10 for output
C                     (currently only applies to CO2 isotopologue 838)
C     22AUG13 AD 2.01 Limit number of warning messages within IDXVIB
C     26JUL13 AD 2.00 Simplified code - remove some options from old version
C                     Directly convert new or old HITRAN format files
C     30JUN09 AD 1.20 HITRAN2008: allow for SBROAD in F5.4 or F5.3 (2008)
C     09MAR02 AD 1.10 DUPCHK: if two lines have equal status change from
C                       termination condition to warning.
C                       Remove redundant MAXNEW parameter.
C     09SEP01 AD 1.00 Original. Based on HITLIN.
C
C DESCRIPTION
C     Convert HITRAN sequential ASCII file to direct access binary file
C     Note: the record length of the unformatted output file depends
C           on the no of bytes per word of the machine being used.
C           System specific RECL is used in the output binary file
C           open statement
C
C     Variables for each transition record (100 characters in input records)
C              Input  Output Description
C             Old/New
C      LSTAT     -      I 4   Priority of transition information.
C      MOL      I2      I 8   Molecule number.
C      ISO      I1      I 12  Isotope number (=1 most abundant, 2=second etc)
C      WNO     F12.6    D 20  Line frequency [cm-1].
C      STR     D10.3    R 24  Line strength  I=[cm-1./(molec.cm-2)] @ 296K.
C                                          O=[cm-1./(kg.moles.cm-2)] @296K
C                           (Scale input by Avogadro No.to avoid underflows.)
C      TPR     E10.3    R 28  Transition probab. [Debyes2] (old) or Einstein A
C      ABR      F5.4    R 32  Air-broadened halfwidth  (HWHM) [cm-1/atm] @ 296K.
C      SBR    F5.4/F5.3 R 36  Self-broadened halfwidth (HWHM) [cm-1/atm] @ 296K.
C      ELS     F10.4    R 40   Lower-state energy [cm-1].
C      ABC      F4.2    R 44   Coefficient of temperature dependance of ABROAD
C      TSP      F8.6    R 48   Transition shift due to pressure
C      IUV    I3/A15    I 52    Upper state global (=Vib) quanta index / ID
C      ILV    I3/A15    I 56    Lower state global quanta index / ID
C      IUR    A9/A15   A9 65   Upper state local (=Rot) quanta.
C      ILR    A9/A15   A9 74   Lower state local quanta.
C      SPARE9          A9 83   Spare (used to be quality and reference info)
C      IFP       -      I 87    Forward pointer on data line.
C      SPARE1          A1 88   Spare (used to be empty byte at end of record)
C
C     LSTAT values for different types of binary file record
C        -7 forward pointer block
C        -2 file termination record
C         0 file header record
C         4 other header records
C        10 line transition
C
      IMPLICIT NONE
C
      EXTERNAL
     &  IDXVIB ! Translate new HITRAN C*15 Vib. ID to old integer index
     &, IRECL  ! Return no.bytes assumed for RECL parameter in OPEN statement
      INTEGER IDXVIB
      INTEGER IRECL
C
C LOCAL CONSTANTS
      INTEGER IFORM               ! format version# for HITRAN output
        PARAMETER ( IFORM = 1 )
      INTEGER LUNHIT              ! LUN for ASCII (input) file
        PARAMETER ( LUNHIT = 1 )
      INTEGER LUNBIN              ! LUN for new binary (output) file
        PARAMETER ( LUNBIN = 2 )
      INTEGER MAXNFP              ! Max no.of forward pointers per record
        PARAMETER ( MAXNFP = 14 ) ! Don't change this
      INTEGER NFPREC              ! No. forward pointer records in a f.p.block
        PARAMETER ( NFPREC = 4 )
      INTEGER MAXMOL              ! Max HITRAN molecule index allowed
        PARAMETER ( MAXMOL = NFPREC*MAXNFP ) ! Change NFPREC to increase MAXMOL
      INTEGER NRECFP              ! No of records between f.p. blocks
        PARAMETER ( NRECFP = 200 )
      REAL AVOG                   ! Avogadro's number (*1e-26)
        PARAMETER ( AVOG = 6.02214199 )
      DOUBLE PRECISION WNOBIG     ! Larger Wno than likely to occur in HITRAN
        PARAMETER ( WNOBIG = 1.0D6 )
C
C LOCAL VARIABLES
      INTEGER IEX             ! STREN exponent read from ASCII file
      INTEGER IFP             ! Forward pointer read from binary file
      INTEGER IFPREC          ! Record# within forward pointer block
      INTEGER ILV             ! Lower level vibrational ID (old format HITRAN)
      INTEGER IMOL            ! Molecule counter
      INTEGER IOFF            ! Offset for molecule# in f.p. records
      INTEGER IOSVAL          ! Saved value of IOSTAT
      INTEGER IRCMOL(MAXMOL)  ! Next record# containing transition for each mol
      INTEGER IREC            ! Record# in new binary file
      INTEGER IREC1           ! First record# in new file after headers
      INTEGER IREC2           ! Last record# in new file (termination rec)
      INTEGER IRECFP          ! First record# in new file after last f.p.block
      INTEGER ISO             ! Isotope ID
      INTEGER IUV             ! Upper level vibrational ID (old format HITRAN)
      INTEGER JREC            ! Secondary record counter
      INTEGER LSTAT           ! Type of binary file record
      INTEGER MOL             ! Molecule ID
      INTEGER RECLEN          ! Record length of binary files
      LOGICAL OLDFMT          ! T=old HITRAN format (<2004), F=new format
      REAL    ABC             ! Air Broad. Temp. Coeff
      REAL    ABR             ! Air Broadened halfwidth
      REAL    ELS             ! Energy of the lower state
      REAL    SBR             ! Self broadened halfwidth
      REAL    STR             ! Line strength
      REAL    TPR             ! Transition Probability
      REAL    TSP             ! Pressure shift
      DOUBLE PRECISION DSTR   ! Line strength allowing for < 1.0E-38
      DOUBLE PRECISION WNO    ! Wavenumber
      DOUBLE PRECISION WNOUPP ! Last wavenumber written/read
      DOUBLE PRECISION WNORQ1 ! Lowest wavenumber to select from ASCII file
      DOUBLE PRECISION WNORQ2 ! Highest wavenumber to select from ASCII file
      CHARACTER*80 FILNAM     ! User-input name of input and output files
      CHARACTER*48 HEAD48     ! User header for binary file
      CHARACTER*84 HEADR2     ! HITBIN header record for binary file
      CHARACTER*9  ILR        ! Lower level rotational ID (old format HITRAN)
      CHARACTER*9  IUR        ! Upper level rotational ID (old format HITRAN)
      CHARACTER*15 JLR        ! Lower level rotational ID (new format HITRAN)
      CHARACTER*15 JLV        ! Lower level vibrational ID (new format HITRAN)
      CHARACTER*15 JUR        ! Upper level rotational ID (new format HITRAN)
      CHARACTER*15 JUV        ! Upper level vibrational ID (new format HITRAN)
      CHARACTER*80  RECORD    ! User input for wavenumber range
      CHARACTER*100 REC100    ! 100-character record for old HITRAN format
      CHARACTER*160 REC160    ! 160-character record for new HITRAN format
      CHARACTER*1  SPARE1     ! Spare byte in output record
      CHARACTER*9  SPARE9     ! Spare bytes in output record (was C*9 acc.index
C                               and data references)
      CHARACTER*4  VERSID     ! Program Version identifier
C
      DATA SPARE9 / '         ' /    ! 9 spaces
      DATA SPARE1 / ' ' /            ! 1 space
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      VERSID = '2.21'
      WRITE ( *, '(A,A,A,I2)' ) 'R-HITBIN: Running HITBIN v', VERSID,
     &                   ' generating output Format#', IFORM
C
C Get HITRAN (ASCII) input file and open
      WRITE ( *, '(A$)' ) 'Input HITRAN file: '
      READ ( *, '(A)' ) FILNAM
      OPEN ( UNIT=LUNHIT, FILE=FILNAM, STATUS='OLD' )
C Check format
      OLDFMT = .FALSE.
      READ ( LUNHIT, '(A)', IOSTAT=IOSVAL ) REC160
      IF ( IOSVAL .NE. 0 ) THEN
        REWIND ( LUNHIT )
        READ ( LUNHIT, '(A)', IOSTAT=IOSVAL ) REC100
        OLDFMT = IOSVAL .EQ. 0
        IF ( OLDFMT ) THEN
          WRITE (*,*) 'File identifed as old (pre 2004) format'
        ELSE
          STOP 'Unable to identify HITRAN input file format'
        ENDIF
      END IF
      REWIND ( LUNHIT )
C
C Get wavenumber range for ASCII file (default = use all)
      WRITE ( *, '(A$)' ) 'Wavenumber range (cm-1) [<CR>=all]: '
      READ ( *, '(A)' ) RECORD
      IF ( RECORD .EQ. ' ' ) THEN
        WNORQ1 = -1.0D0
        WNORQ2 = WNOBIG
      ELSE
        READ ( RECORD, * ) WNORQ1, WNORQ2
        WRITE ( *, '(A)' )
     &    'I-HITBIN: Finding first record of ASCII file...'
        WNO = -1.0D0
        DO WHILE ( WNO .LT. WNORQ1 )
          READ ( LUNHIT, '(3X,F12.6)', END=900 ) WNO
        END DO
        IF ( WNO .GT. WNORQ2 ) GOTO 900           ! no overlap with reqd range
        BACKSPACE ( LUNHIT )
      END IF
C
C Get name of the new binary file to be created
      WRITE ( *, '(A$)' ) 'New binary file: '
      READ ( *, '(A)' ) FILNAM
C Get record length for new binary file
      RECLEN = 22 * IRECL ( LUNBIN )
      WRITE ( *, * ) 'Assuming appropriate RECL value is:', RECLEN
      OPEN ( UNIT=LUNBIN, FILE=FILNAM, STATUS='NEW', ACCESS='DIRECT',
     &       RECL=RECLEN )
C
C Get header for new file
      WRITE ( *, '(A$)' ) 'Header for new file (up to 48 chars): '
      READ ( *, '(A)' ) HEAD48
C
C Write HITBIN header (record length and version ID) to rec#2
      WRITE ( HEADR2, '(A,I3,A,A)' ) 'RECL=', RECLEN,
     &  ' Converted to binary format by HITBIN v.', VERSID
      WRITE ( LUNBIN, REC=2 ) 4, HEADR2
C
      IREC = 3
      IREC1 = IREC
C
      WRITE ( *, '(A)' ) 'I-HITBIN: Writing new binary file...'
C Begin output data with an empty forward pointer block
      DO JREC = 1, NFPREC
        WRITE ( LUNBIN, REC=IREC ) -7
        IREC = IREC + 1
      END DO
      IRECFP = IREC
C
C Repeat for each record, taking lowest wavenumber
      DO WHILE ( WNO .LE. WNORQ2 )
C Check if forward pointer block required
        IF ( IREC - IRECFP .EQ. NRECFP ) THEN          ! insert f.p. block
          DO JREC = 1, NFPREC
            WRITE ( LUNBIN, REC=IREC ) -7
            IREC = IREC + 1
          END DO
          IRECFP = IREC
        END IF

C Copy record from ASCII file to new file
C To avoid underflow problems, STR read as double precision and scaled by
C Avogadro's number before converting to single precision
C Note that format descriptors Fw.d (eg F5.2) do not generally represent the
C actual HITRAN format (eg F5.3) but simply set to avoid warning messages from
C some compilers (eg ifort) which recommend w ge d+3
        IF ( OLDFMT ) THEN
          READ ( LUNHIT, 1000, END=800 ) MOL, ISO, WNO, DSTR, TPR,
     &        ABR, SBR, ELS, ABC, TSP, IUV, ILV, IUR, ILR
 1000     FORMAT( I2, I1, F12.6, F10.3, E10.3,
     &            f5.2, F5.2, F10.4, F4.1, F8.5, 2I3, 2A9 )
          STR = AVOG * STR * ( 10.0**( IEX + 26 ) )
        ELSE
          READ ( LUNHIT, 1001, END=800 ) MOL, ISO, WNO, DSTR, TPR,
     &        ABR, SBR, ELS, ABC, TSP, JUV, JLV, JUR, JLR
 1001       FORMAT( I2, I1, F12.6, F10.3, E10.3,
     &            F5.2, F5.2, F10.4, F4.1, F8.5, 4A15 )
          STR = SNGL ( DSTR * AVOG * 1.0E26 )

C Convert new format C*15 codes for upper,lower vibrational level into old
C format integer indices required by RFM if particular Vib levels have to be
C identified (user-selected or non-LTE).
          IUV = IDXVIB ( MOL, JUV )
          ILV = IDXVIB ( MOL, JLV )
C Reformat local quantum numbers for CO2 lines (required by RFM for line-mixing)
C The first 9 characters of the HITRAN04 C*15 strings seems to match the
C old C*9 fields except that in the old format the characters are shifted one
C space left and the new format has an extra character 'e' or 'f' which needs
C to be set blank (this was fixed in v2.21 of this code)
C No idea what this will do for local quantum numbers for other molecules (!)
C but hopefully this information won't be required
          IUR = JUR(2:9)//' '
          ILR = JLR(2:9)//' '
        ENDIF
C
C HITRAN 2012 has 10 CO2 isotopologues with '0' used for #10. Use 10 in bin file
        IF ( ISO .EQ. 0 ) ISO = 10
C
        IF ( WNO .LE. WNORQ2 ) THEN
          WRITE ( LUNBIN, REC=IREC ) 10, MOL, ISO, WNO, STR, TPR, ! 10=LSTAT
     &       ABR, SBR, ELS, ABC, TSP, IUV, ILV, IUR, ILR,
     &       SPARE9, 0, SPARE1   ! 0 = forward ptr
          IREC = IREC + 1
          WNOUPP = WNO
C
          IF ( MOD ( IREC, 100000 ) .EQ. 0 )
     &      WRITE ( *, '(A,I8,A,F12.6)' )
     &      'I-HITBIN: Record#', IREC, ' Wavenumber=', WNOUPP
        END IF
      END DO
  800 CONTINUE
C
C Save record after last data record
      IREC2 = IREC
C
C Write header record
      WRITE ( LUNBIN, REC=1 ) 0, IFORM, IREC1, IREC2, HEAD48
C
C Write termination record
      WRITE ( LUNBIN, REC=IREC2 ) -2, 0, 0, WNOUPP, (0,IMOL=1,14)
      WRITE ( *, '(A,I8,A,F12.6)' )
     &  'I-HITBIN: Last Record#', IREC2, ' Wavenumber=', WNOUPP
C
C Calculate forward pointers
      WRITE ( *, '(A)' ) 'I-HITBIN: Writing forward pointers...'
      DO IMOL = 1, MAXMOL
        IRCMOL(IMOL) = IREC2
      END DO
      IFPREC = NFPREC
      DO IREC = IREC2-1, IREC1, -1
        READ ( LUNBIN, REC=IREC ) LSTAT
        IF ( LSTAT .EQ. -7 ) THEN
          IOFF = ( IFPREC - 1 ) * MAXNFP
          WRITE ( LUNBIN, REC=IREC ) -7, IOFF+1, 0, WNOUPP,
     &      ( IRCMOL(IMOL)-IREC, IMOL = IOFF+1, IOFF+MAXNFP )
          IF ( IFPREC .EQ. 1 ) THEN
            IFPREC = NFPREC
          ELSE
            IFPREC = IFPREC - 1
          END IF
        ELSE
          READ ( LUNBIN, REC=IREC ) LSTAT, MOL, ISO, WNO, STR, TPR,
     &      ABR, SBR, ELS, ABC, TSP, IUV, ILV, IUR, ILR, SPARE9, IFP,
     &      SPARE1
          IFP = IRCMOL(MOL) - IREC
          IRCMOL(MOL) = IREC
          WNOUPP = WNO
          WRITE ( LUNBIN, REC=IREC ) LSTAT, MOL, ISO, WNO, STR, TPR,
     &      ABR, SBR, ELS, ABC, TSP, IUV, ILV, IUR, ILR, SPARE9, IFP,
     &      SPARE1
        END IF
        IF ( MOD ( IREC, 100000 ) .EQ. 0 )
     &    WRITE ( *, '(A,I8,A,F12.6)' )
     &    'I-HITBIN: Record#', IREC, ' Wavenumber=', WNOUPP
      END DO
C
      STOP  'R-HITBIN: Successful completion'
C
 900  CONTINUE
      STOP 'F-HITBIN: HITRAN file does not overlap required wno.range'

      END
C
      INTEGER FUNCTION IRECL ( LUN )
C
C VERSION
C     25JUL13  AD  Original.
C
C DESCRIPTION
C     Return no.bytes assumed for RECL parameter in OPEN statement
C     RECL specifies record lengths either in bytes or units of 4 bytes,
C     depending on the compiler and compilation options.
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER  LUN  !  I  Spare LUN used for testing (freed after use)
C
C LOCAL VARIABLES
      INTEGER  I, J ! dummy variables written to file
      INTEGER  IOS  ! saved value of IOSTAT
C
C EXECUTABLE CODE --------------------------------------------------------------
C
      I = 0
      J = 0
      IRECL = 1
      OPEN ( UNIT=LUN, RECL=2*IRECL, ACCESS='DIRECT', STATUS='SCRATCH' )
      WRITE ( LUN, REC=1, IOSTAT=IOS ) I, J
      CLOSE ( LUN )
      IF ( IOS .EQ. 0 ) RETURN

      IRECL = 4
      OPEN ( UNIT=LUN, RECL=2*IRECL, ACCESS='DIRECT', STATUS='SCRATCH' )
      WRITE ( LUN, REC=1, IOSTAT=IOS ) I, J
      CLOSE ( LUN )
      IF ( IOS .EQ. 0 ) RETURN
C
      STOP 'F-IRECL: Unable to establish RECL definition'
C
      END
      INTEGER FUNCTION IDXVIB ( IDXMOL, STRVIB )
C
C VERSION
C     22AUG13  AD Limit number of warning messages
C     26JUL13  AD Add new molec 43-47 & vib states found in HITRAN2012
C     30JUN09  AD Add new HITRAN molec#40-42
C     19DEC07  AD Set IVIB locally, and set unidentified levels to 999
C     24MAR06  AD Original.
C
C DESCRIPTION
C     Translate HITRAN C*15 Global Quantum ID to old integer index
C     Based on a program supplied by Javier Martin-Torres (NASA LaRC).
C
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      IVIB    !  Local value of IDXVIB
      INTEGER      IDXMOL  !  HITRAN code for molecule
      CHARACTER*15 STRVIB  !  HITRAN04 field defining vibration level
C
C LOCAL CONSTANTS
      INTEGER MAXMOL       !  .GE. Max number of different molecules
        PARAMETER ( MAXMOL = 47 )
      INTEGER MAXWRN       !  Max number of warning messages per molecules
        PARAMETER ( MAXWRN = 10 )
      INTEGER MAXNEW       !  Max number of new molecules beyond MAXMOL
        PARAMETER ( MAXNEW = 10 )
C
C LOCAL VARIABLES
      INTEGER INEW           ! Counter for new molecules
      INTEGER NNEW           ! No. of new molecules found
      INTEGER NWRN(MAXMOL) ! No. warnings per molecule.
      INTEGER NEWLST(MAXNEW) ! List of IDXMOL values for new molecules
C
C DATA STATEMENTS
      DATA NWRN / MAXMOL * 0 /
      DATA NNEW / 0 /
      DATA NEWLST / MAXNEW * 0 /
      SAVE NWRN, NNEW, NEWLST
C
C EXECUTABLE CODE -------------------------------------------------------------
C
      IVIB = 0
C
      IF ( IDXMOL .GT. MAXMOL ) THEN
        IVIB = 0
        DO INEW = 1, NNEW
          IF ( IDXMOL .EQ. NEWLST(INEW) ) RETURN
        END DO
        IF ( NNEW .LT. MAXNEW ) THEN
          NNEW = NNEW + 1
          NEWLST(NNEW) = IDXMOL
          WRITE (*,*)
     &      'W-IDXVIB: Setting IVIB=0 for unidentified IDXMOL=', IDXMOL
        END IF
        RETURN
      END IF
C
C Class 1: Diatomic Molecules
      IF ( IDXMOL .EQ.  5 .OR.                                    ! CO
     &     IDXMOL .EQ. 14 .OR.                                    ! HF
     &     IDXMOL .EQ. 15 .OR.                                    ! HCl
     &     IDXMOL .EQ. 16 .OR.                                    ! HBr
     &     IDXMOL .EQ. 17 .OR.                                    ! HI
     &     IDXMOL .EQ. 22 .OR.                                    ! N2
     &     IDXMOL .EQ. 36 .OR.                                    ! NO+
     &     IDXMOL .EQ. 45 .OR.                                    ! H2
     &     IDXMOL .EQ. 46      ) THEN                             ! CS
        IF ( STRVIB .EQ. '              0' ) IVIB = 1
        IF ( STRVIB .EQ. '              1' ) IVIB = 2
        IF ( STRVIB .EQ. '              2' ) IVIB = 3
        IF ( STRVIB .EQ. '              3' ) IVIB = 4
        IF ( STRVIB .EQ. '              4' ) IVIB = 5
        IF ( STRVIB .EQ. '              5' ) IVIB = 6
        IF ( STRVIB .EQ. '              6' ) IVIB = 7
        IF ( STRVIB .EQ. '              7' ) IVIB = 8
        IF ( STRVIB .EQ. '              8' ) IVIB = 9
        IF ( STRVIB .EQ. '              9' ) IVIB = 10
        IF ( STRVIB .EQ. '             10' ) IVIB = 11
        IF ( STRVIB .EQ. '             11' ) IVIB = 12
        IF ( STRVIB .EQ. '             12' ) IVIB = 13
        IF ( STRVIB .EQ. '             13' ) IVIB = 14
        IF ( STRVIB .EQ. '             14' ) IVIB = 15
        IF ( STRVIB .EQ. '             15' ) IVIB = 16
        IF ( STRVIB .EQ. '             16' ) IVIB = 17
        IF ( STRVIB .EQ. '             17' ) IVIB = 18
        IF ( STRVIB .EQ. '             18' ) IVIB = 19
        IF ( STRVIB .EQ. '             19' ) IVIB = 20
        IF ( STRVIB .EQ. '             20' ) IVIB = 21
C new levels for HF, HCl found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '             21' ) IVIB = 22
        IF ( STRVIB .EQ. '             22' ) IVIB = 23
        IF ( STRVIB .EQ. '             23' ) IVIB = 24
        IF ( STRVIB .EQ. '             24' ) IVIB = 25
        IF ( STRVIB .EQ. '             25' ) IVIB = 26
        IF ( STRVIB .EQ. '             26' ) IVIB = 27
C
C Class 2: Diatomic molecules with different electonic levels
C HITRAN 2012 uses different spacing to previous versions so allow for both
      ELSE IF ( IDXMOL .EQ.  7 ) THEN                              ! O2
	IF ( STRVIB .EQ. '            X 0' ) IVIB = 1
        IF ( STRVIB .EQ. '       X      0' ) IVIB = 1    ! HITRAN 2012
	IF ( STRVIB .EQ. '            X 1' ) IVIB = 2
        IF ( STRVIB .EQ. '       X      1' ) IVIB = 2    ! HITRAN 2012
	IF ( STRVIB .EQ. '            a 0' ) IVIB = 3
        IF ( STRVIB .EQ. '       a      0' ) IVIB = 3    ! HITRAN 2012
	IF ( STRVIB .EQ. '            a 1' ) IVIB = 4
        IF ( STRVIB .EQ. '       a      1' ) IVIB = 4    ! HITRAN 2012
	IF ( STRVIB .EQ. '            b 0' ) IVIB = 5
        IF ( STRVIB .EQ. '       b      0' ) IVIB = 5    ! HITRAN 2012
	IF ( STRVIB .EQ. '            b 1' ) IVIB = 6
        IF ( STRVIB .EQ. '       b      1' ) IVIB = 6    ! HITRAN 2012
	IF ( STRVIB .EQ. '            b 2' ) IVIB = 7
        IF ( STRVIB .EQ. '       b      2' ) IVIB = 7    ! HITRAN 2012
	IF ( STRVIB .EQ. '               ' ) IVIB = 8
	IF ( STRVIB .EQ. '            X 2' ) IVIB = 9
	IF ( STRVIB .EQ. '            B 0' ) IVIB = 10
	IF ( STRVIB .EQ. '            B 1' ) IVIB = 11
	IF ( STRVIB .EQ. '            B 2' ) IVIB = 12
	IF ( STRVIB .EQ. '            B 3' ) IVIB = 13
	IF ( STRVIB .EQ. '            B 4' ) IVIB = 14
	IF ( STRVIB .EQ. '            B 5' ) IVIB = 15
	IF ( STRVIB .EQ. '            B 6' ) IVIB = 16
	IF ( STRVIB .EQ. '            B 7' ) IVIB = 17
	IF ( STRVIB .EQ. '            B 8' ) IVIB = 18
	IF ( STRVIB .EQ. '            B 9' ) IVIB = 19
	IF ( STRVIB .EQ. '           B 10' ) IVIB = 20
	IF ( STRVIB .EQ. '           B 11' ) IVIB = 21
	IF ( STRVIB .EQ. '           B 12' ) IVIB = 22
	IF ( STRVIB .EQ. '           B 13' ) IVIB = 23
	IF ( STRVIB .EQ. '           B 14' ) IVIB = 24
	IF ( STRVIB .EQ. '           B 15' ) IVIB = 25
	IF ( STRVIB .EQ. '           B 16' ) IVIB = 26
	IF ( STRVIB .EQ. '           B 17' ) IVIB = 27
	IF ( STRVIB .EQ. '           B 18' ) IVIB = 28
	IF ( STRVIB .EQ. '           B 19' ) IVIB = 29
C
C Class 3: Diatomic molecules with Pi-doublet electronic state
      ELSE IF ( IDXMOL .EQ.  8 .OR.                                ! NO
     &          IDXMOL .EQ. 13 .OR.                                ! OH
     &          IDXMOL .EQ. 18      ) THEN                         ! ClO
	IF ( STRVIB .EQ. '       X3/2   0' ) IVIB = 1
	IF ( STRVIB .EQ. '       X3/2   1' ) IVIB = 2
        IF ( STRVIB .EQ. '       X3/2   2' ) IVIB = 3
	IF ( STRVIB .EQ. '       X3/2   3' ) IVIB = 4
	IF ( STRVIB .EQ. '       X3/2   4' ) IVIB = 5
	IF ( STRVIB .EQ. '       X3/2   5' ) IVIB = 6
	IF ( STRVIB .EQ. '       X3/2   6' ) IVIB = 7
	IF ( STRVIB .EQ. '       X3/2   7' ) IVIB = 8
	IF ( STRVIB .EQ. '       X3/2   8' ) IVIB = 9
	IF ( STRVIB .EQ. '       X3/2   9' ) IVIB = 10
	IF ( STRVIB .EQ. '       X1/2   0' ) IVIB = 11
	IF ( STRVIB .EQ. '       X1/2   1' ) IVIB = 12
        IF ( STRVIB .EQ. '       X1/2   2' ) IVIB = 13
	IF ( STRVIB .EQ. '       X1/2   3' ) IVIB = 14
	IF ( STRVIB .EQ. '       X1/2   4' ) IVIB = 15
	IF ( STRVIB .EQ. '       X1/2   5' ) IVIB = 16
	IF ( STRVIB .EQ. '       X1/2   6' ) IVIB = 17
	IF ( STRVIB .EQ. '       X1/2   7' ) IVIB = 18
	IF ( STRVIB .EQ. '       X1/2   8' ) IVIB = 19
	IF ( STRVIB .EQ. '       X1/2   9' ) IVIB = 20
	IF ( STRVIB .EQ. '       X1/2  10' ) IVIB = 21
	IF ( STRVIB .EQ. '       X1/2  11' ) IVIB = 22
	IF ( STRVIB .EQ. '       X1/2  12' ) IVIB = 23
	IF ( STRVIB .EQ. '       X3/2  10' ) IVIB = 24
	IF ( STRVIB .EQ. '       X3/2  11' ) IVIB = 25
	IF ( STRVIB .EQ. '       X3/2  12' ) IVIB = 26
	IF ( STRVIB .EQ. '       A1     0' ) IVIB = 27
	IF ( STRVIB .EQ. '       A1     1' ) IVIB = 28
	IF ( STRVIB .EQ. '       A1     2' ) IVIB = 29
	IF ( STRVIB .EQ. '       A1     3' ) IVIB = 30
	IF ( STRVIB .EQ. '       A2     0' ) IVIB = 31
	IF ( STRVIB .EQ. '       A2     1' ) IVIB = 32
	IF ( STRVIB .EQ. '       A2     2' ) IVIB = 33
	IF ( STRVIB .EQ. '       A2     3' ) IVIB = 34
	IF ( STRVIB .EQ. '       X3/2  13' ) IVIB = 35
	IF ( STRVIB .EQ. '       X3/2  14' ) IVIB = 36
	IF ( STRVIB .EQ. '       X1/2  13' ) IVIB = 37
	IF ( STRVIB .EQ. '       X1/2  14' ) IVIB = 38
C
C Class 4: Linear triatomic molecules
      ELSE IF ( IDXMOL .EQ.  4 .OR.                                 ! N2O
     &          IDXMOL .EQ. 19 .OR.                                 ! OCS
     &          IDXMOL .EQ. 23      ) THEN                          ! HCN
	IF ( STRVIB .EQ. '        0 0 0 0' ) IVIB = 1
 	IF ( STRVIB .EQ. '        0 1 1 0' ) IVIB = 2
 	IF ( STRVIB .EQ. '        0 2 0 0' ) IVIB = 3
	IF ( STRVIB .EQ. '        0 2 2 0' ) IVIB = 4
	IF ( STRVIB .EQ. '        1 0 0 0' ) IVIB = 5
 	IF ( STRVIB .EQ. '        0 3 1 0' ) IVIB = 6
 	IF ( STRVIB .EQ. '        0 3 3 0' ) IVIB = 7
	IF ( STRVIB .EQ. '        1 1 1 0' ) IVIB = 8
 	IF ( STRVIB .EQ. '        0 4 0 0' ) IVIB = 9
 	IF ( STRVIB .EQ. '        0 4 2 0' ) IVIB = 10
	IF ( STRVIB .EQ. '        1 2 0 0' ) IVIB = 11
	IF ( STRVIB .EQ. '        1 2 2 0' ) IVIB = 12
	IF ( STRVIB .EQ. '        2 0 0 0' ) IVIB = 13
	IF ( STRVIB .EQ. '        0 0 0 1' ) IVIB = 14
 	IF ( STRVIB .EQ. '        0 5 1 0' ) IVIB = 15
 	IF ( STRVIB .EQ. '        1 3 1 0' ) IVIB = 16
	IF ( STRVIB .EQ. '        1 3 3 0' ) IVIB = 17
 	IF ( STRVIB .EQ. '        2 1 1 0' ) IVIB = 18
 	IF ( STRVIB .EQ. '        0 1 1 1' ) IVIB = 19
 	IF ( STRVIB .EQ. '        1 4 0 0' ) IVIB = 20
 	IF ( STRVIB .EQ. '        1 4 2 0' ) IVIB = 21
 	IF ( STRVIB .EQ. '        2 2 0 0' ) IVIB = 22
 	IF ( STRVIB .EQ. '        2 2 2 0' ) IVIB = 23
 	IF ( STRVIB .EQ. '        3 0 0 0' ) IVIB = 24
 	IF ( STRVIB .EQ. '        0 2 0 1' ) IVIB = 25
 	IF ( STRVIB .EQ. '        0 2 2 1' ) IVIB = 26
 	IF ( STRVIB .EQ. '        1 0 0 1' ) IVIB = 27
	IF ( STRVIB .EQ. '        2 3 1 0' ) IVIB = 28
	IF ( STRVIB .EQ. '        3 1 1 0' ) IVIB = 29
	IF ( STRVIB .EQ. '        0 3 1 1' ) IVIB = 30
 	IF ( STRVIB .EQ. '        0 3 3 1' ) IVIB = 31
 	IF ( STRVIB .EQ. '        1 1 1 1' ) IVIB = 32
 	IF ( STRVIB .EQ. '        4 0 0 0' ) IVIB = 33
 	IF ( STRVIB .EQ. '        3 2 0 0' ) IVIB = 34
 	IF ( STRVIB .EQ. '        2 0 0 1' ) IVIB = 35
 	IF ( STRVIB .EQ. '        1 2 0 1' ) IVIB = 36
 	IF ( STRVIB .EQ. '        1 2 2 1' ) IVIB = 37
 	IF ( STRVIB .EQ. '        0 0 0 2' ) IVIB = 38
	IF ( STRVIB .EQ. '        2 1 1 1' ) IVIB = 39
	IF ( STRVIB .EQ. '        0 1 1 2' ) IVIB = 40
 	IF ( STRVIB .EQ. '               ' ) IVIB = 41
	IF ( STRVIB .EQ. '        0 6 0 0' ) IVIB = 42
	IF ( STRVIB .EQ. '        0 6 2 0' ) IVIB = 43
	IF ( STRVIB .EQ. '        0 4 4 0' ) IVIB = 44
	IF ( STRVIB .EQ. '        0 5 3 0' ) IVIB = 45
	IF ( STRVIB .EQ. '        0 4 4 1' ) IVIB = 46
 	IF ( STRVIB .EQ. '        0 4 2 1' ) IVIB = 47
 	IF ( STRVIB .EQ. '        0 4 0 1' ) IVIB = 48
 	IF ( STRVIB .EQ. '        1 5 1 0' ) IVIB = 49
	IF ( STRVIB .EQ. '        1 5 3 0' ) IVIB = 50
	IF ( STRVIB .EQ. '        2 3 3 0' ) IVIB = 51
	IF ( STRVIB .EQ. '        0 5 3 1' ) IVIB = 52
	IF ( STRVIB .EQ. '        0 5 1 1' ) IVIB = 53
 	IF ( STRVIB .EQ. '        1 0 0 2' ) IVIB = 54
 	IF ( STRVIB .EQ. '        1 3 1 1' ) IVIB = 55
 	IF ( STRVIB .EQ. '        1 3 3 1' ) IVIB = 56
 	IF ( STRVIB .EQ. '        0 7 3 0' ) IVIB = 57
	IF ( STRVIB .EQ. '        0 7 1 0' ) IVIB = 58
 	IF ( STRVIB .EQ. '        1 6 0 0' ) IVIB = 59
 	IF ( STRVIB .EQ. '        2 4 0 0' ) IVIB = 60
	IF ( STRVIB .EQ. '        2 4 2 0' ) IVIB = 61
 	IF ( STRVIB .EQ. '        4 1 1 0' ) IVIB = 62
	IF ( STRVIB .EQ. '        3 2 2 0' ) IVIB = 63
	IF ( STRVIB .EQ. '        0 2 2 2' ) IVIB = 64
	IF ( STRVIB .EQ. '        0 2 0 2' ) IVIB = 65
 	IF ( STRVIB .EQ. '        1 4 0 1' ) IVIB = 66
	IF ( STRVIB .EQ. '        1 4 2 1' ) IVIB = 67
	IF ( STRVIB .EQ. '        2 2 0 1' ) IVIB = 68
	IF ( STRVIB .EQ. '        3 0 0 1' ) IVIB = 69
 	IF ( STRVIB .EQ. '        2 5 1 0' ) IVIB = 70
	IF ( STRVIB .EQ. '        4 2 0 0' ) IVIB = 71
	IF ( STRVIB .EQ. '        3 3 1 0' ) IVIB = 72
	IF ( STRVIB .EQ. '        0 6 2 1' ) IVIB = 73
 	IF ( STRVIB .EQ. '        1 1 1 2' ) IVIB = 74
 	IF ( STRVIB .EQ. '        2 3 1 1' ) IVIB = 75
 	IF ( STRVIB .EQ. '        3 1 1 1' ) IVIB = 76
	IF ( STRVIB .EQ. '        3 4 0 0' ) IVIB = 77
	IF ( STRVIB .EQ. '        5 0 0 0' ) IVIB = 78
	IF ( STRVIB .EQ. '        0 1 1 3' ) IVIB = 79
 	IF ( STRVIB .EQ. '        0 0 0 3' ) IVIB = 80
 	IF ( STRVIB .EQ. '        2 0 0 2' ) IVIB = 81
 	IF ( STRVIB .EQ. '        3 2 0 1' ) IVIB = 82
	IF ( STRVIB .EQ. '        4 0 0 1' ) IVIB = 83
 	IF ( STRVIB .EQ. '        1 0 0 3' ) IVIB = 84
 	IF ( STRVIB .EQ. '        2 2 2 1' ) IVIB = 85
	IF ( STRVIB .EQ. '        0 9 1 0' ) IVIB = 86
	IF ( STRVIB .EQ. '        0 7 1 1' ) IVIB = 87
	IF ( STRVIB .EQ. '        0 2 0 2' ) IVIB = 88
	IF ( STRVIB .EQ. '        0 5 1 1' ) IVIB = 89
	IF ( STRVIB .EQ. '        1 1 1 2' ) IVIB = 90
	IF ( STRVIB .EQ. '        3 4 2 0' ) IVIB = 91
C new levels for OCS found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '        0 3 1 2' ) IVIB = 92
        IF ( STRVIB .EQ. '        0 4 0 2' ) IVIB = 93
        IF ( STRVIB .EQ. '        0 4 2 2' ) IVIB = 94
        IF ( STRVIB .EQ. '        1 2 0 2' ) IVIB = 95
        IF ( STRVIB .EQ. '        1 2 2 2' ) IVIB = 96
C
C Class 5: Linear triatomic molecules with large Fermi resonance
      ELSE IF ( IDXMOL .EQ.  2 ) THEN                                ! CO2
	IF ( STRVIB .EQ. '       0 0 0 01' ) IVIB = 1
	IF ( STRVIB .EQ. '       0 1 1 01' ) IVIB = 2
	IF ( STRVIB .EQ. '       1 0 0 02' ) IVIB = 3
	IF ( STRVIB .EQ. '       0 2 2 01' ) IVIB = 4
 	IF ( STRVIB .EQ. '       1 0 0 01' ) IVIB = 5
 	IF ( STRVIB .EQ. '       1 1 1 02' ) IVIB = 6
 	IF ( STRVIB .EQ. '       0 3 3 01' ) IVIB = 7
 	IF ( STRVIB .EQ. '       1 1 1 01' ) IVIB = 8
	IF ( STRVIB .EQ. '       0 0 0 11' ) IVIB = 9
	IF ( STRVIB .EQ. '       2 0 0 03' ) IVIB = 10
	IF ( STRVIB .EQ. '       1 2 2 02' ) IVIB = 11
	IF ( STRVIB .EQ. '       2 0 0 02' ) IVIB = 12
	IF ( STRVIB .EQ. '       0 4 4 01' ) IVIB = 13
	IF ( STRVIB .EQ. '       1 2 2 01' ) IVIB = 14
	IF ( STRVIB .EQ. '       2 0 0 01' ) IVIB = 15
	IF ( STRVIB .EQ. '       0 1 1 11' ) IVIB = 16
	IF ( STRVIB .EQ. '       2 1 1 03' ) IVIB = 17
	IF ( STRVIB .EQ. '       1 3 3 02' ) IVIB = 18
	IF ( STRVIB .EQ. '       2 1 1 02' ) IVIB = 19
	IF ( STRVIB .EQ. '       0 5 5 01' ) IVIB = 20
	IF ( STRVIB .EQ. '       1 3 3 01' ) IVIB = 21
	IF ( STRVIB .EQ. '       2 1 1 01' ) IVIB = 22
	IF ( STRVIB .EQ. '       1 0 0 12' ) IVIB = 23
	IF ( STRVIB .EQ. '       0 2 2 11' ) IVIB = 24
	IF ( STRVIB .EQ. '       1 0 0 11' ) IVIB = 25
	IF ( STRVIB .EQ. '       3 0 0 04' ) IVIB = 26
	IF ( STRVIB .EQ. '       2 2 2 03' ) IVIB = 27
	IF ( STRVIB .EQ. '       1 4 4 02' ) IVIB = 28
	IF ( STRVIB .EQ. '       3 0 0 03' ) IVIB = 29
	IF ( STRVIB .EQ. '       2 2 2 02' ) IVIB = 30
	IF ( STRVIB .EQ. '       0 6 6 01' ) IVIB = 31
	IF ( STRVIB .EQ. '       3 0 0 02' ) IVIB = 32
	IF ( STRVIB .EQ. '       1 4 4 01' ) IVIB = 33
	IF ( STRVIB .EQ. '       2 2 2 01' ) IVIB = 34
	IF ( STRVIB .EQ. '       3 0 0 01' ) IVIB = 35
	IF ( STRVIB .EQ. '       1 1 1 12' ) IVIB = 36
	IF ( STRVIB .EQ. '       0 3 3 11' ) IVIB = 37
	IF ( STRVIB .EQ. '       1 1 1 11' ) IVIB = 38
	IF ( STRVIB .EQ. '       0 0 0 21' ) IVIB = 39
	IF ( STRVIB .EQ. '       3 1 1 04' ) IVIB = 40
	IF ( STRVIB .EQ. '       3 1 1 03' ) IVIB = 41
	IF ( STRVIB .EQ. '       3 1 1 02' ) IVIB = 42
	IF ( STRVIB .EQ. '       2 0 0 13' ) IVIB = 43
	IF ( STRVIB .EQ. '       1 2 2 12' ) IVIB = 44
	IF ( STRVIB .EQ. '       2 3 3 01' ) IVIB = 45
	IF ( STRVIB .EQ. '       3 1 1 01' ) IVIB = 46
	IF ( STRVIB .EQ. '       0 4 4 11' ) IVIB = 47
	IF ( STRVIB .EQ. '       2 0 0 12' ) IVIB = 48
	IF ( STRVIB .EQ. '       1 2 2 11' ) IVIB = 49
	IF ( STRVIB .EQ. '       2 0 0 11' ) IVIB = 50
	IF ( STRVIB .EQ. '       0 1 1 21' ) IVIB = 51
	IF ( STRVIB .EQ. '       4 0 0 04' ) IVIB = 52
	IF ( STRVIB .EQ. '       3 2 2 03' ) IVIB = 53
	IF ( STRVIB .EQ. '       2 1 1 13' ) IVIB = 54
	IF ( STRVIB .EQ. '       4 0 0 02' ) IVIB = 55
	IF ( STRVIB .EQ. '       1 3 3 12' ) IVIB = 56
	IF ( STRVIB .EQ. '       0 5 5 11' ) IVIB = 57
	IF ( STRVIB .EQ. '       2 1 1 12' ) IVIB = 58
	IF ( STRVIB .EQ. '       1 3 3 11' ) IVIB = 59
	IF ( STRVIB .EQ. '       2 1 1 11' ) IVIB = 60
	IF ( STRVIB .EQ. '       1 0 0 22' ) IVIB = 61
	IF ( STRVIB .EQ. '       0 2 2 21' ) IVIB = 62
	IF ( STRVIB .EQ. '       1 0 0 21' ) IVIB = 63
	IF ( STRVIB .EQ. '       3 0 0 14' ) IVIB = 64
	IF ( STRVIB .EQ. '       2 2 2 13' ) IVIB = 65
	IF ( STRVIB .EQ. '       1 4 4 12' ) IVIB = 66
	IF ( STRVIB .EQ. '       4 1 1 02' ) IVIB = 67
	IF ( STRVIB .EQ. '       3 0 0 13' ) IVIB = 68
	IF ( STRVIB .EQ. '       0 6 6 11' ) IVIB = 69
	IF ( STRVIB .EQ. '       2 2 2 12' ) IVIB = 70
	IF ( STRVIB .EQ. '       3 0 0 12' ) IVIB = 71
	IF ( STRVIB .EQ. '       4 1 1 01' ) IVIB = 72
	IF ( STRVIB .EQ. '       1 4 4 11' ) IVIB = 73
	IF ( STRVIB .EQ. '       2 2 2 11' ) IVIB = 74
	IF ( STRVIB .EQ. '       3 0 0 11' ) IVIB = 75
	IF ( STRVIB .EQ. '       1 1 1 22' ) IVIB = 76
	IF ( STRVIB .EQ. '       0 3 3 21' ) IVIB = 77
	IF ( STRVIB .EQ. '       1 1 1 21' ) IVIB = 78
	IF ( STRVIB .EQ. '       0 0 0 31' ) IVIB = 79
	IF ( STRVIB .EQ. '       3 1 1 14' ) IVIB = 80
	IF ( STRVIB .EQ. '       2 3 3 13' ) IVIB = 81
	IF ( STRVIB .EQ. '       3 1 1 13' ) IVIB = 82
	IF ( STRVIB .EQ. '       2 3 3 12' ) IVIB = 83
	IF ( STRVIB .EQ. '       3 1 1 12' ) IVIB = 84
	IF ( STRVIB .EQ. '       1 5 5 11' ) IVIB = 85
	IF ( STRVIB .EQ. '       2 0 0 23' ) IVIB = 86
	IF ( STRVIB .EQ. '       2 3 3 11' ) IVIB = 87
	IF ( STRVIB .EQ. '       1 2 2 22' ) IVIB = 88
	IF ( STRVIB .EQ. '       3 1 1 11' ) IVIB = 89
	IF ( STRVIB .EQ. '       2 0 0 22' ) IVIB = 90
	IF ( STRVIB .EQ. '       1 2 2 21' ) IVIB = 91
	IF ( STRVIB .EQ. '       2 0 0 21' ) IVIB = 92
	IF ( STRVIB .EQ. '       0 1 1 31' ) IVIB = 93
	IF ( STRVIB .EQ. '       4 0 0 15' ) IVIB = 94
	IF ( STRVIB .EQ. '       3 2 2 14' ) IVIB = 95
	IF ( STRVIB .EQ. '       4 0 0 14' ) IVIB = 96
	IF ( STRVIB .EQ. '       3 2 2 13' ) IVIB = 97
	IF ( STRVIB .EQ. '       4 0 0 13' ) IVIB = 98
	IF ( STRVIB .EQ. '       5 1 1 02' ) IVIB = 99
	IF ( STRVIB .EQ. '       3 2 2 12' ) IVIB = 100
	IF ( STRVIB .EQ. '       4 0 0 12' ) IVIB = 101
	IF ( STRVIB .EQ. '       2 1 1 23' ) IVIB = 102
	IF ( STRVIB .EQ. '       3 2 2 11' ) IVIB = 103
	IF ( STRVIB .EQ. '       2 1 1 22' ) IVIB = 104
	IF ( STRVIB .EQ. '       4 0 0 11' ) IVIB = 105
	IF ( STRVIB .EQ. '       2 1 1 21' ) IVIB = 106
	IF ( STRVIB .EQ. '       1 0 0 32' ) IVIB = 107
	IF ( STRVIB .EQ. '       0 2 2 31' ) IVIB = 108
	IF ( STRVIB .EQ. '       1 0 0 31' ) IVIB = 109
	IF ( STRVIB .EQ. '       4 1 1 14' ) IVIB = 110
	IF ( STRVIB .EQ. '       4 1 1 13' ) IVIB = 111
	IF ( STRVIB .EQ. '       4 1 1 12' ) IVIB = 112
	IF ( STRVIB .EQ. '       1 1 1 32' ) IVIB = 113
	IF ( STRVIB .EQ. '       0 3 3 31' ) IVIB = 114
	IF ( STRVIB .EQ. '       1 1 1 31' ) IVIB = 115
	IF ( STRVIB .EQ. '       2 0 0 33' ) IVIB = 116
	IF ( STRVIB .EQ. '       1 2 2 32' ) IVIB = 117
	IF ( STRVIB .EQ. '       2 0 0 32' ) IVIB = 118
	IF ( STRVIB .EQ. '       1 2 2 31' ) IVIB = 119
	IF ( STRVIB .EQ. '       2 0 0 31' ) IVIB = 120
	IF ( STRVIB .EQ. '       2 1 1 33' ) IVIB = 121
	IF ( STRVIB .EQ. '       2 1 1 32' ) IVIB = 122
	IF ( STRVIB .EQ. '       2 1 1 31' ) IVIB = 123
	IF ( STRVIB .EQ. '       2 3 3 03' ) IVIB = 124
	IF ( STRVIB .EQ. '       1 5 5 02' ) IVIB = 125
	IF ( STRVIB .EQ. '       2 3 3 02' ) IVIB = 126
	IF ( STRVIB .EQ. '       0 7 7 01' ) IVIB = 127
	IF ( STRVIB .EQ. '               ' ) IVIB = 128
	IF ( STRVIB .EQ. '       1 0 0 41' ) IVIB = 129
	IF ( STRVIB .EQ. '       1 0 0 51' ) IVIB = 130
	IF ( STRVIB .EQ. '       1 0 0 52' ) IVIB = 131
	IF ( STRVIB .EQ. '       0 0 0 51' ) IVIB = 132
C new levels found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '       0 0 0 00' ) IVIB = 133
        IF ( STRVIB .EQ. '       0 0 0 41' ) IVIB = 134
        IF ( STRVIB .EQ. '       0 1 1 41' ) IVIB = 135
        IF ( STRVIB .EQ. '       0 1 1 51' ) IVIB = 136
        IF ( STRVIB .EQ. '       0 2 2 41' ) IVIB = 137
        IF ( STRVIB .EQ. '       0 2 2 51' ) IVIB = 138
        IF ( STRVIB .EQ. '       0 4 4 21' ) IVIB = 139
        IF ( STRVIB .EQ. '       0 4 4 31' ) IVIB = 140
        IF ( STRVIB .EQ. '       0 5 5 21' ) IVIB = 141
        IF ( STRVIB .EQ. '       0 5 5 31' ) IVIB = 142
        IF ( STRVIB .EQ. '       0 7 7 11' ) IVIB = 143
        IF ( STRVIB .EQ. '       0 8 8 01' ) IVIB = 144
        IF ( STRVIB .EQ. '       0 8 8 11' ) IVIB = 145
        IF ( STRVIB .EQ. '       0 9 9 01' ) IVIB = 146
        IF ( STRVIB .EQ. '       1 0 0 42' ) IVIB = 147
        IF ( STRVIB .EQ. '       1 1 1 41' ) IVIB = 148
        IF ( STRVIB .EQ. '       1 1 1 42' ) IVIB = 149
        IF ( STRVIB .EQ. '       1 1 1 51' ) IVIB = 150
        IF ( STRVIB .EQ. '       1 1 1 52' ) IVIB = 151
        IF ( STRVIB .EQ. '       1 2 2 51' ) IVIB = 152
        IF ( STRVIB .EQ. '       1 2 2 52' ) IVIB = 153
        IF ( STRVIB .EQ. '       1 3 3 21' ) IVIB = 154
        IF ( STRVIB .EQ. '       1 3 3 22' ) IVIB = 155
        IF ( STRVIB .EQ. '       1 3 3 31' ) IVIB = 156
        IF ( STRVIB .EQ. '       1 3 3 32' ) IVIB = 157
        IF ( STRVIB .EQ. '       1 4 4 21' ) IVIB = 158
        IF ( STRVIB .EQ. '       1 4 4 22' ) IVIB = 159
        IF ( STRVIB .EQ. '       1 4 4 31' ) IVIB = 160
        IF ( STRVIB .EQ. '       1 4 4 32' ) IVIB = 161
        IF ( STRVIB .EQ. '       1 5 5 01' ) IVIB = 162
        IF ( STRVIB .EQ. '       1 5 5 12' ) IVIB = 163
        IF ( STRVIB .EQ. '       1 6 6 01' ) IVIB = 164
        IF ( STRVIB .EQ. '       1 6 6 02' ) IVIB = 165
        IF ( STRVIB .EQ. '       1 6 6 11' ) IVIB = 166
        IF ( STRVIB .EQ. '       1 6 6 12' ) IVIB = 167
        IF ( STRVIB .EQ. '       1 7 7 01' ) IVIB = 168
        IF ( STRVIB .EQ. '       1 7 7 02' ) IVIB = 169
        IF ( STRVIB .EQ. '       1 7 7 11' ) IVIB = 170
        IF ( STRVIB .EQ. '       1 7 7 12' ) IVIB = 171
        IF ( STRVIB .EQ. '       2 0 0 41' ) IVIB = 172
        IF ( STRVIB .EQ. '       2 0 0 42' ) IVIB = 173
        IF ( STRVIB .EQ. '       2 0 0 51' ) IVIB = 174
        IF ( STRVIB .EQ. '       2 0 0 52' ) IVIB = 175
        IF ( STRVIB .EQ. '       2 0 0 53' ) IVIB = 176
        IF ( STRVIB .EQ. '       2 2 2 21' ) IVIB = 177
        IF ( STRVIB .EQ. '       2 2 2 22' ) IVIB = 178
        IF ( STRVIB .EQ. '       2 2 2 23' ) IVIB = 179
        IF ( STRVIB .EQ. '       2 2 2 31' ) IVIB = 180
        IF ( STRVIB .EQ. '       2 2 2 32' ) IVIB = 181
        IF ( STRVIB .EQ. '       2 2 2 33' ) IVIB = 182
        IF ( STRVIB .EQ. '       2 3 3 21' ) IVIB = 183
        IF ( STRVIB .EQ. '       2 3 3 22' ) IVIB = 184
        IF ( STRVIB .EQ. '       2 3 3 23' ) IVIB = 185
        IF ( STRVIB .EQ. '       2 3 3 31' ) IVIB = 186
        IF ( STRVIB .EQ. '       2 3 3 32' ) IVIB = 187
        IF ( STRVIB .EQ. '       2 3 3 33' ) IVIB = 188
        IF ( STRVIB .EQ. '       2 4 4 01' ) IVIB = 189
        IF ( STRVIB .EQ. '       2 4 4 02' ) IVIB = 190
        IF ( STRVIB .EQ. '       2 4 4 03' ) IVIB = 191
        IF ( STRVIB .EQ. '       2 4 4 11' ) IVIB = 192
        IF ( STRVIB .EQ. '       2 4 4 12' ) IVIB = 193
        IF ( STRVIB .EQ. '       2 4 4 13' ) IVIB = 194
        IF ( STRVIB .EQ. '       2 5 5 01' ) IVIB = 195
        IF ( STRVIB .EQ. '       2 5 5 02' ) IVIB = 196
        IF ( STRVIB .EQ. '       2 5 5 03' ) IVIB = 197
        IF ( STRVIB .EQ. '       2 5 5 11' ) IVIB = 198
        IF ( STRVIB .EQ. '       2 5 5 12' ) IVIB = 199
        IF ( STRVIB .EQ. '       2 5 5 13' ) IVIB = 200
        IF ( STRVIB .EQ. '       2 6 6 01' ) IVIB = 201
        IF ( STRVIB .EQ. '       2 6 6 02' ) IVIB = 202
        IF ( STRVIB .EQ. '       2 6 6 11' ) IVIB = 203
        IF ( STRVIB .EQ. '       2 6 6 12' ) IVIB = 204
        IF ( STRVIB .EQ. '       2 6 6 13' ) IVIB = 205
        IF ( STRVIB .EQ. '       3 0 0 21' ) IVIB = 206
        IF ( STRVIB .EQ. '       3 0 0 22' ) IVIB = 207
        IF ( STRVIB .EQ. '       3 0 0 23' ) IVIB = 208
        IF ( STRVIB .EQ. '       3 0 0 24' ) IVIB = 209
        IF ( STRVIB .EQ. '       3 0 0 31' ) IVIB = 210
        IF ( STRVIB .EQ. '       3 0 0 32' ) IVIB = 211
        IF ( STRVIB .EQ. '       3 0 0 33' ) IVIB = 212
        IF ( STRVIB .EQ. '       3 0 0 34' ) IVIB = 213
        IF ( STRVIB .EQ. '       3 1 1 21' ) IVIB = 214
        IF ( STRVIB .EQ. '       3 1 1 22' ) IVIB = 215
        IF ( STRVIB .EQ. '       3 1 1 23' ) IVIB = 216
        IF ( STRVIB .EQ. '       3 1 1 24' ) IVIB = 217
        IF ( STRVIB .EQ. '       3 1 1 31' ) IVIB = 218
        IF ( STRVIB .EQ. '       3 1 1 32' ) IVIB = 219
        IF ( STRVIB .EQ. '       3 1 1 33' ) IVIB = 220
        IF ( STRVIB .EQ. '       3 1 1 34' ) IVIB = 221
        IF ( STRVIB .EQ. '       3 2 2 01' ) IVIB = 222
        IF ( STRVIB .EQ. '       3 2 2 02' ) IVIB = 223
        IF ( STRVIB .EQ. '       3 2 2 04' ) IVIB = 224
        IF ( STRVIB .EQ. '       3 3 3 01' ) IVIB = 225
        IF ( STRVIB .EQ. '       3 3 3 02' ) IVIB = 226
        IF ( STRVIB .EQ. '       3 3 3 03' ) IVIB = 227
        IF ( STRVIB .EQ. '       3 3 3 04' ) IVIB = 228
        IF ( STRVIB .EQ. '       3 3 3 11' ) IVIB = 229
        IF ( STRVIB .EQ. '       3 3 3 12' ) IVIB = 230
        IF ( STRVIB .EQ. '       3 3 3 13' ) IVIB = 231
        IF ( STRVIB .EQ. '       3 3 3 14' ) IVIB = 232
        IF ( STRVIB .EQ. '       3 4 4 01' ) IVIB = 233
        IF ( STRVIB .EQ. '       3 4 4 02' ) IVIB = 234
        IF ( STRVIB .EQ. '       3 4 4 03' ) IVIB = 235
        IF ( STRVIB .EQ. '       3 4 4 04' ) IVIB = 236
        IF ( STRVIB .EQ. '       3 4 4 11' ) IVIB = 237
        IF ( STRVIB .EQ. '       3 4 4 12' ) IVIB = 238
        IF ( STRVIB .EQ. '       3 4 4 13' ) IVIB = 239
        IF ( STRVIB .EQ. '       3 4 4 14' ) IVIB = 240
        IF ( STRVIB .EQ. '       3 5 5 01' ) IVIB = 241
        IF ( STRVIB .EQ. '       3 5 5 02' ) IVIB = 242
        IF ( STRVIB .EQ. '       3 5 5 12' ) IVIB = 243
        IF ( STRVIB .EQ. '       3 5 5 13' ) IVIB = 244
        IF ( STRVIB .EQ. '       4 0 0 01' ) IVIB = 245
        IF ( STRVIB .EQ. '       4 0 0 03' ) IVIB = 246
        IF ( STRVIB .EQ. '       4 0 0 05' ) IVIB = 247
        IF ( STRVIB .EQ. '       4 0 0 21' ) IVIB = 248
        IF ( STRVIB .EQ. '       4 0 0 22' ) IVIB = 249
        IF ( STRVIB .EQ. '       4 0 0 23' ) IVIB = 250
        IF ( STRVIB .EQ. '       4 0 0 24' ) IVIB = 251
        IF ( STRVIB .EQ. '       4 0 0 33' ) IVIB = 252
        IF ( STRVIB .EQ. '       4 0 0 34' ) IVIB = 253
        IF ( STRVIB .EQ. '       4 1 1 03' ) IVIB = 254
        IF ( STRVIB .EQ. '       4 1 1 04' ) IVIB = 255
        IF ( STRVIB .EQ. '       4 1 1 05' ) IVIB = 256
        IF ( STRVIB .EQ. '       4 1 1 11' ) IVIB = 257
        IF ( STRVIB .EQ. '       4 1 1 15' ) IVIB = 258
        IF ( STRVIB .EQ. '       4 1 1 25' ) IVIB = 259
        IF ( STRVIB .EQ. '       4 2 2 01' ) IVIB = 260
        IF ( STRVIB .EQ. '       4 2 2 02' ) IVIB = 261
        IF ( STRVIB .EQ. '       4 2 2 03' ) IVIB = 262
        IF ( STRVIB .EQ. '       4 2 2 04' ) IVIB = 263
        IF ( STRVIB .EQ. '       4 2 2 05' ) IVIB = 264
        IF ( STRVIB .EQ. '       4 2 2 11' ) IVIB = 265
        IF ( STRVIB .EQ. '       4 2 2 12' ) IVIB = 266
        IF ( STRVIB .EQ. '       4 2 2 13' ) IVIB = 267
        IF ( STRVIB .EQ. '       4 2 2 14' ) IVIB = 268
        IF ( STRVIB .EQ. '       4 2 2 15' ) IVIB = 269
        IF ( STRVIB .EQ. '       4 3 3 01' ) IVIB = 270
        IF ( STRVIB .EQ. '       4 3 3 02' ) IVIB = 271
        IF ( STRVIB .EQ. '       4 3 3 03' ) IVIB = 272
        IF ( STRVIB .EQ. '       4 3 3 04' ) IVIB = 273
        IF ( STRVIB .EQ. '       4 3 3 11' ) IVIB = 274
        IF ( STRVIB .EQ. '       4 3 3 12' ) IVIB = 275
        IF ( STRVIB .EQ. '       4 3 3 13' ) IVIB = 276
        IF ( STRVIB .EQ. '       4 3 3 14' ) IVIB = 277
        IF ( STRVIB .EQ. '       4 3 3 15' ) IVIB = 278
        IF ( STRVIB .EQ. '       4 4 4 01' ) IVIB = 279
        IF ( STRVIB .EQ. '       4 4 4 02' ) IVIB = 280
        IF ( STRVIB .EQ. '       4 4 4 03' ) IVIB = 281
        IF ( STRVIB .EQ. '       5 0 0 01' ) IVIB = 282
        IF ( STRVIB .EQ. '       5 0 0 02' ) IVIB = 283
        IF ( STRVIB .EQ. '       5 0 0 03' ) IVIB = 284
        IF ( STRVIB .EQ. '       5 0 0 04' ) IVIB = 285
        IF ( STRVIB .EQ. '       5 0 0 05' ) IVIB = 286
        IF ( STRVIB .EQ. '       5 0 0 06' ) IVIB = 287
        IF ( STRVIB .EQ. '       5 0 0 11' ) IVIB = 288
        IF ( STRVIB .EQ. '       5 0 0 12' ) IVIB = 289
        IF ( STRVIB .EQ. '       5 0 0 13' ) IVIB = 290
        IF ( STRVIB .EQ. '       5 0 0 14' ) IVIB = 291
        IF ( STRVIB .EQ. '       5 0 0 15' ) IVIB = 292
        IF ( STRVIB .EQ. '       5 0 0 16' ) IVIB = 293
        IF ( STRVIB .EQ. '       5 1 1 01' ) IVIB = 294
        IF ( STRVIB .EQ. '       5 1 1 03' ) IVIB = 295
        IF ( STRVIB .EQ. '       5 1 1 04' ) IVIB = 296
        IF ( STRVIB .EQ. '       5 1 1 05' ) IVIB = 297
        IF ( STRVIB .EQ. '       5 1 1 06' ) IVIB = 298
        IF ( STRVIB .EQ. '       5 1 1 11' ) IVIB = 299
        IF ( STRVIB .EQ. '       5 1 1 12' ) IVIB = 300
        IF ( STRVIB .EQ. '       5 1 1 13' ) IVIB = 301
        IF ( STRVIB .EQ. '       5 1 1 14' ) IVIB = 302
        IF ( STRVIB .EQ. '       5 1 1 15' ) IVIB = 303
        IF ( STRVIB .EQ. '       5 1 1 16' ) IVIB = 304
        IF ( STRVIB .EQ. '       5 2 2 01' ) IVIB = 305
        IF ( STRVIB .EQ. '       5 2 2 02' ) IVIB = 306
        IF ( STRVIB .EQ. '       5 2 2 03' ) IVIB = 307
        IF ( STRVIB .EQ. '       5 2 2 04' ) IVIB = 308
        IF ( STRVIB .EQ. '       5 2 2 05' ) IVIB = 309
        IF ( STRVIB .EQ. '       5 2 2 11' ) IVIB = 310
        IF ( STRVIB .EQ. '       5 2 2 12' ) IVIB = 311
        IF ( STRVIB .EQ. '       5 2 2 13' ) IVIB = 312
        IF ( STRVIB .EQ. '       5 2 2 14' ) IVIB = 313
        IF ( STRVIB .EQ. '       5 2 2 15' ) IVIB = 314
        IF ( STRVIB .EQ. '       5 3 3 02' ) IVIB = 315
        IF ( STRVIB .EQ. '       5 4 4 03' ) IVIB = 316
        IF ( STRVIB .EQ. '       6 0 0 01' ) IVIB = 317
        IF ( STRVIB .EQ. '       6 0 0 02' ) IVIB = 318
        IF ( STRVIB .EQ. '       6 0 0 03' ) IVIB = 319
        IF ( STRVIB .EQ. '       6 0 0 04' ) IVIB = 320
        IF ( STRVIB .EQ. '       6 0 0 06' ) IVIB = 321
        IF ( STRVIB .EQ. '       6 0 0 12' ) IVIB = 322
        IF ( STRVIB .EQ. '       6 0 0 13' ) IVIB = 323
        IF ( STRVIB .EQ. '       6 0 0 14' ) IVIB = 324
        IF ( STRVIB .EQ. '       6 0 0 15' ) IVIB = 325
        IF ( STRVIB .EQ. '       6 0 0 16' ) IVIB = 326
        IF ( STRVIB .EQ. '       6 0 0 17' ) IVIB = 327
        IF ( STRVIB .EQ. '       6 1 1 01' ) IVIB = 328
        IF ( STRVIB .EQ. '       6 1 1 02' ) IVIB = 329
        IF ( STRVIB .EQ. '       6 1 1 03' ) IVIB = 330
        IF ( STRVIB .EQ. '       6 1 1 04' ) IVIB = 331
        IF ( STRVIB .EQ. '       6 1 1 14' ) IVIB = 332
        IF ( STRVIB .EQ. '       6 1 1 15' ) IVIB = 333
C
C Class 6: Non-linear triatomic molecules
      ELSE IF ( IDXMOL .EQ.  1 .OR.                                ! H2O
     &          IDXMOL .EQ.  3 .OR. 	                           ! O3
     &          IDXMOL .EQ.  9 .OR. 	                           ! SO2
     &          IDXMOL .EQ. 10 .OR. 	                           ! NO2
     &          IDXMOL .EQ. 21 .OR. 	                           ! HOCl
     &          IDXMOL .EQ. 31 .OR. 	                           ! H2S
     &          IDXMOL .EQ. 33 .OR. 	                           ! HO2
     &          IDXMOL .EQ. 37      ) THEN                         ! HOBr
	IF ( STRVIB .EQ. '          0 0 0' ) IVIB = 1
	IF ( STRVIB .EQ. '          0 1 0' ) IVIB = 2
 	IF ( STRVIB .EQ. '          0 2 0' ) IVIB = 3
	IF ( STRVIB .EQ. '          1 0 0' ) IVIB = 4
	IF ( STRVIB .EQ. '          0 0 1' ) IVIB = 5
 	IF ( STRVIB .EQ. '          0 3 0' ) IVIB = 6
 	IF ( STRVIB .EQ. '          1 1 0' ) IVIB = 7
	IF ( STRVIB .EQ. '          0 1 1' ) IVIB = 8
 	IF ( STRVIB .EQ. '          0 4 0' ) IVIB = 9
	IF ( STRVIB .EQ. '          1 2 0' ) IVIB = 10
	IF ( STRVIB .EQ. '          0 2 1' ) IVIB = 11
	IF ( STRVIB .EQ. '          2 0 0' ) IVIB = 12
	IF ( STRVIB .EQ. '          1 0 1' ) IVIB = 13
 	IF ( STRVIB .EQ. '          0 0 2' ) IVIB = 14
 	IF ( STRVIB .EQ. '          1 3 0' ) IVIB = 15
	IF ( STRVIB .EQ. '          0 3 1' ) IVIB = 16
	IF ( STRVIB .EQ. '          2 1 0' ) IVIB = 17
	IF ( STRVIB .EQ. '          1 1 1' ) IVIB = 18
 	IF ( STRVIB .EQ. '          0 1 2' ) IVIB = 19
	IF ( STRVIB .EQ. '          0 4 1' ) IVIB = 20
	IF ( STRVIB .EQ. '          2 2 0' ) IVIB = 21
	IF ( STRVIB .EQ. '          1 2 1' ) IVIB = 22
	IF ( STRVIB .EQ. '          0 2 2' ) IVIB = 23
 	IF ( STRVIB .EQ. '          3 0 0' ) IVIB = 24
	IF ( STRVIB .EQ. '          2 0 1' ) IVIB = 25
	IF ( STRVIB .EQ. '          1 0 2' ) IVIB = 26
 	IF ( STRVIB .EQ. '          0 0 3' ) IVIB = 27
	IF ( STRVIB .EQ. '          1 3 1' ) IVIB = 28
	IF ( STRVIB .EQ. '          3 1 0' ) IVIB = 29
	IF ( STRVIB .EQ. '          2 1 1' ) IVIB = 30
	IF ( STRVIB .EQ. '          1 1 2' ) IVIB = 31
 	IF ( STRVIB .EQ. '          0 1 3' ) IVIB = 32
	IF ( STRVIB .EQ. '          1 4 1' ) IVIB = 33
	IF ( STRVIB .EQ. '          0 4 2' ) IVIB = 34
	IF ( STRVIB .EQ. '          3 2 0' ) IVIB = 35
 	IF ( STRVIB .EQ. '          2 2 1' ) IVIB = 36
	IF ( STRVIB .EQ. '          3 0 1' ) IVIB = 37
 	IF ( STRVIB .EQ. '          2 0 2' ) IVIB = 38
	IF ( STRVIB .EQ. '          1 2 2' ) IVIB = 39
	IF ( STRVIB .EQ. '          0 2 3' ) IVIB = 40
	IF ( STRVIB .EQ. '          4 0 0' ) IVIB = 41
	IF ( STRVIB .EQ. '          1 0 3' ) IVIB = 42
	IF ( STRVIB .EQ. '          0 0 4' ) IVIB = 43
	IF ( STRVIB .EQ. '          1 5 1' ) IVIB = 44
	IF ( STRVIB .EQ. '          3 3 0' ) IVIB = 45
	IF ( STRVIB .EQ. '          2 3 1' ) IVIB = 46
	IF ( STRVIB .EQ. '          2 1 2' ) IVIB = 47
	IF ( STRVIB .EQ. '          3 1 1' ) IVIB = 48
	IF ( STRVIB .EQ. '          4 1 0' ) IVIB = 49
	IF ( STRVIB .EQ. '          1 1 3' ) IVIB = 50
	IF ( STRVIB .EQ. '          3 2 1' ) IVIB = 51
	IF ( STRVIB .EQ. '          2 2 2' ) IVIB = 52
	IF ( STRVIB .EQ. '          3 0 2' ) IVIB = 53
	IF ( STRVIB .EQ. '          4 0 1' ) IVIB = 54
	IF ( STRVIB .EQ. '          4 2 0' ) IVIB = 55
	IF ( STRVIB .EQ. '          1 2 3' ) IVIB = 56
	IF ( STRVIB .EQ. '          5 0 0' ) IVIB = 57
	IF ( STRVIB .EQ. '          2 0 3' ) IVIB = 58
	IF ( STRVIB .EQ. '          1 0 4' ) IVIB = 59
	IF ( STRVIB .EQ. '               ' ) IVIB = 60
	IF ( STRVIB .EQ. '          3 3 1' ) IVIB = 61
	IF ( STRVIB .EQ. '          2 1 3' ) IVIB = 62
	IF ( STRVIB .EQ. '          3 1 2' ) IVIB = 63
	IF ( STRVIB .EQ. '          4 1 1' ) IVIB = 64
	IF ( STRVIB .EQ. '          3 0 3' ) IVIB = 65
	IF ( STRVIB .EQ. '          4 0 2' ) IVIB = 66
	IF ( STRVIB .EQ. '          4 0 3' ) IVIB = 67
	IF ( STRVIB .EQ. '          4 2 1' ) IVIB = 68
	IF ( STRVIB .EQ. '          5 0 1' ) IVIB = 69
	IF ( STRVIB .EQ. '          3 1 3' ) IVIB = 70
	IF ( STRVIB .EQ. '          4 1 2' ) IVIB = 71
	IF ( STRVIB .EQ. '          2 3 2' ) IVIB = 72
	IF ( STRVIB .EQ. '          0 5 0' ) IVIB = 73
	IF ( STRVIB .EQ. '          0 6 0' ) IVIB = 74
	IF ( STRVIB .EQ. '          0 7 0' ) IVIB = 75
	IF ( STRVIB .EQ. '          0 3 2' ) IVIB = 76
	IF ( STRVIB .EQ. '          0 5 1' ) IVIB = 77
	IF ( STRVIB .EQ. '          0 6 1' ) IVIB = 78
	IF ( STRVIB .EQ. '          0 8 0' ) IVIB = 79
	IF ( STRVIB .EQ. '          1 4 0' ) IVIB = 80
	IF ( STRVIB .EQ. '          1 5 0' ) IVIB = 81
	IF ( STRVIB .EQ. '          0 3 3' ) IVIB = 82
	IF ( STRVIB .EQ. '          0 3 4' ) IVIB = 83
	IF ( STRVIB .EQ. '          0 4 3' ) IVIB = 84
	IF ( STRVIB .EQ. '          0 5 3' ) IVIB = 85
	IF ( STRVIB .EQ. '          0 6 1' ) IVIB = 86
	IF ( STRVIB .EQ. '          0 6 3' ) IVIB = 87
	IF ( STRVIB .EQ. '          0 7 1' ) IVIB = 88
	IF ( STRVIB .EQ. '          1 1 5' ) IVIB = 89
	IF ( STRVIB .EQ. '          1 3 2' ) IVIB = 90
	IF ( STRVIB .EQ. '          1 3 3' ) IVIB = 91
	IF ( STRVIB .EQ. '          1 4 2' ) IVIB = 92
	IF ( STRVIB .EQ. '          1 6 0' ) IVIB = 93
	IF ( STRVIB .EQ. '          1 7 0' ) IVIB = 94
	IF ( STRVIB .EQ. '          2 2 3' ) IVIB = 95
	IF ( STRVIB .EQ. '          2 4 0' ) IVIB = 96
	IF ( STRVIB .EQ. '          2 4 1' ) IVIB = 97
	IF ( STRVIB .EQ. '          3 2 2' ) IVIB = 98
	IF ( STRVIB .EQ. '          3 4 0' ) IVIB = 99
	IF ( STRVIB .EQ. '          3 4 1' ) IVIB = 100
	IF ( STRVIB .EQ. '          4 3 0' ) IVIB = 101
	IF ( STRVIB .EQ. '          4 3 1' ) IVIB = 102
	IF ( STRVIB .EQ. '          5 1 0' ) IVIB = 103
	IF ( STRVIB .EQ. '          5 1 1' ) IVIB = 104
	IF ( STRVIB .EQ. '          5 2 0' ) IVIB = 105
	IF ( STRVIB .EQ. '          6 0 0' ) IVIB = 106
	IF ( STRVIB .EQ. '          6 0 1' ) IVIB = 107
	IF ( STRVIB .EQ. '          6 1 0' ) IVIB = 108
	IF ( STRVIB .EQ. '          6 1 1' ) IVIB = 109
	IF ( STRVIB .EQ. '          6 2 0' ) IVIB = 110
	IF ( STRVIB .EQ. '          7 0 0' ) IVIB = 111
	IF ( STRVIB .EQ. '          7 0 1' ) IVIB = 112
	IF ( STRVIB .EQ. '          8 0 0' ) IVIB = 113
	IF ( STRVIB .EQ. '          0 5 0' ) IVIB = 114
 	IF ( STRVIB .EQ. '          0 6 0' ) IVIB = 115
C new levels for H2O found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '          0 0 5' ) IVIB = 116
        IF ( STRVIB .EQ. '          0 0 6' ) IVIB = 117
        IF ( STRVIB .EQ. '          0 0 7' ) IVIB = 118
        IF ( STRVIB .EQ. '          0 1 4' ) IVIB = 119
        IF ( STRVIB .EQ. '          0 1 5' ) IVIB = 120
        IF ( STRVIB .EQ. '          0 2 4' ) IVIB = 121
        IF ( STRVIB .EQ. '          0 2 5' ) IVIB = 122
        IF ( STRVIB .EQ. '          0 3 5' ) IVIB = 123
        IF ( STRVIB .EQ. '          0 4 4' ) IVIB = 124
        IF ( STRVIB .EQ. '          0 4 5' ) IVIB = 125
        IF ( STRVIB .EQ. '          0 5 2' ) IVIB = 126
        IF ( STRVIB .EQ. '          0 6 2' ) IVIB = 127
        IF ( STRVIB .EQ. '          0 7 2' ) IVIB = 128
        IF ( STRVIB .EQ. '          0 7 3' ) IVIB = 129
        IF ( STRVIB .EQ. '          0 8 1' ) IVIB = 130
        IF ( STRVIB .EQ. '          0 8 2' ) IVIB = 131
        IF ( STRVIB .EQ. '          0 9 0' ) IVIB = 132
        IF ( STRVIB .EQ. '          0 9 1' ) IVIB = 133
        IF ( STRVIB .EQ. '          0 9 2' ) IVIB = 134
        IF ( STRVIB .EQ. '          0 9 3' ) IVIB = 135
        IF ( STRVIB .EQ. '          010 0' ) IVIB = 136
        IF ( STRVIB .EQ. '          010 1' ) IVIB = 137
        IF ( STRVIB .EQ. '          011 0' ) IVIB = 138
        IF ( STRVIB .EQ. '          011 1' ) IVIB = 139
        IF ( STRVIB .EQ. '          012 0' ) IVIB = 140
        IF ( STRVIB .EQ. '          012 1' ) IVIB = 141
        IF ( STRVIB .EQ. '          013 0' ) IVIB = 142
        IF ( STRVIB .EQ. '          014 0' ) IVIB = 143
        IF ( STRVIB .EQ. '          015 0' ) IVIB = 144
        IF ( STRVIB .EQ. '          1 0 5' ) IVIB = 145
        IF ( STRVIB .EQ. '          1 0 6' ) IVIB = 146
        IF ( STRVIB .EQ. '          1 1 4' ) IVIB = 147
        IF ( STRVIB .EQ. '          1 2 4' ) IVIB = 148
        IF ( STRVIB .EQ. '          1 2 5' ) IVIB = 149
        IF ( STRVIB .EQ. '          1 3 4' ) IVIB = 150
        IF ( STRVIB .EQ. '          1 4 3' ) IVIB = 151
        IF ( STRVIB .EQ. '          1 5 2' ) IVIB = 152
        IF ( STRVIB .EQ. '          1 5 3' ) IVIB = 153
        IF ( STRVIB .EQ. '          1 5 4' ) IVIB = 154
        IF ( STRVIB .EQ. '          1 6 1' ) IVIB = 155
        IF ( STRVIB .EQ. '          1 6 2' ) IVIB = 156
        IF ( STRVIB .EQ. '          1 6 3' ) IVIB = 157
        IF ( STRVIB .EQ. '          1 7 1' ) IVIB = 158
        IF ( STRVIB .EQ. '          1 8 0' ) IVIB = 159
        IF ( STRVIB .EQ. '          1 8 1' ) IVIB = 160
        IF ( STRVIB .EQ. '          1 8 2' ) IVIB = 161
        IF ( STRVIB .EQ. '          1 9 0' ) IVIB = 162
        IF ( STRVIB .EQ. '          1 9 1' ) IVIB = 163
        IF ( STRVIB .EQ. '          110 0' ) IVIB = 164
        IF ( STRVIB .EQ. '          110 1' ) IVIB = 165
        IF ( STRVIB .EQ. '          111 0' ) IVIB = 166
        IF ( STRVIB .EQ. '          111 1' ) IVIB = 167
        IF ( STRVIB .EQ. '          114 0' ) IVIB = 168
        IF ( STRVIB .EQ. '          2 0 4' ) IVIB = 169
        IF ( STRVIB .EQ. '          2 1 4' ) IVIB = 170
        IF ( STRVIB .EQ. '          2 2 4' ) IVIB = 171
        IF ( STRVIB .EQ. '          2 3 0' ) IVIB = 172
        IF ( STRVIB .EQ. '          2 3 3' ) IVIB = 173
        IF ( STRVIB .EQ. '          2 4 2' ) IVIB = 174
        IF ( STRVIB .EQ. '          2 4 3' ) IVIB = 175
        IF ( STRVIB .EQ. '          2 5 0' ) IVIB = 176
        IF ( STRVIB .EQ. '          2 5 1' ) IVIB = 177
        IF ( STRVIB .EQ. '          2 5 3' ) IVIB = 178
        IF ( STRVIB .EQ. '          2 6 0' ) IVIB = 179
        IF ( STRVIB .EQ. '          2 6 1' ) IVIB = 180
        IF ( STRVIB .EQ. '          2 7 0' ) IVIB = 181
        IF ( STRVIB .EQ. '          2 7 1' ) IVIB = 182
        IF ( STRVIB .EQ. '          2 7 2' ) IVIB = 183
        IF ( STRVIB .EQ. '          2 8 0' ) IVIB = 184
        IF ( STRVIB .EQ. '          2 8 1' ) IVIB = 185
        IF ( STRVIB .EQ. '          2 9 0' ) IVIB = 186
        IF ( STRVIB .EQ. '          210 0' ) IVIB = 187
        IF ( STRVIB .EQ. '          210 1' ) IVIB = 188
        IF ( STRVIB .EQ. '          211 0' ) IVIB = 189
        IF ( STRVIB .EQ. '          3 2 3' ) IVIB = 190
        IF ( STRVIB .EQ. '          3 3 2' ) IVIB = 191
        IF ( STRVIB .EQ. '          3 3 3' ) IVIB = 192
        IF ( STRVIB .EQ. '          3 4 2' ) IVIB = 193
        IF ( STRVIB .EQ. '          3 5 0' ) IVIB = 194
        IF ( STRVIB .EQ. '          3 5 1' ) IVIB = 195
        IF ( STRVIB .EQ. '          3 5 2' ) IVIB = 196
        IF ( STRVIB .EQ. '          3 6 0' ) IVIB = 197
        IF ( STRVIB .EQ. '          3 6 1' ) IVIB = 198
        IF ( STRVIB .EQ. '          3 7 0' ) IVIB = 199
        IF ( STRVIB .EQ. '          3 7 1' ) IVIB = 200
        IF ( STRVIB .EQ. '          3 8 0' ) IVIB = 201
        IF ( STRVIB .EQ. '          3 8 1' ) IVIB = 202
        IF ( STRVIB .EQ. '          4 2 2' ) IVIB = 203
        IF ( STRVIB .EQ. '          4 3 2' ) IVIB = 204
        IF ( STRVIB .EQ. '          4 4 0' ) IVIB = 205
        IF ( STRVIB .EQ. '          4 4 1' ) IVIB = 206
        IF ( STRVIB .EQ. '          4 5 0' ) IVIB = 207
        IF ( STRVIB .EQ. '          4 5 1' ) IVIB = 208
        IF ( STRVIB .EQ. '          4 6 0' ) IVIB = 209
        IF ( STRVIB .EQ. '          4 7 0' ) IVIB = 210
        IF ( STRVIB .EQ. '          5 0 2' ) IVIB = 211
        IF ( STRVIB .EQ. '          5 1 2' ) IVIB = 212
        IF ( STRVIB .EQ. '          5 2 1' ) IVIB = 213
        IF ( STRVIB .EQ. '          5 3 0' ) IVIB = 214
        IF ( STRVIB .EQ. '          5 3 1' ) IVIB = 215
        IF ( STRVIB .EQ. '          5 5 0' ) IVIB = 216
        IF ( STRVIB .EQ. '          6 2 1' ) IVIB = 217
        IF ( STRVIB .EQ. '          6 3 0' ) IVIB = 218
        IF ( STRVIB .EQ. '          6 4 0' ) IVIB = 219
        IF ( STRVIB .EQ. '          7 1 0' ) IVIB = 220
        IF ( STRVIB .EQ. '         -2-2-2' ) IVIB = 221
C new levels for O3 found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '       I  1 0 5' ) IVIB = 222
        IF ( STRVIB .EQ. '       I  1 2 4' ) IVIB = 223
        IF ( STRVIB .EQ. '       I  2 0 5' ) IVIB = 224
        IF ( STRVIB .EQ. '       I  2 2 3' ) IVIB = 225
        IF ( STRVIB .EQ. '       I  2 3 3' ) IVIB = 226
        IF ( STRVIB .EQ. '      II  1 0 5' ) IVIB = 227
        IF ( STRVIB .EQ. '      II  1 2 4' ) IVIB = 228
        IF ( STRVIB .EQ. '      II  2 2 3' ) IVIB = 229
        IF ( STRVIB .EQ. '      II  2 3 3' ) IVIB = 230
C
C Class 7: Linear tetratomic molecules
      ELSE IF ( IDXMOL .EQ. 26 ) THEN                                ! C2H2
 	IF ( STRVIB .EQ. ' 0 0 0 0 1 1   ' ) IVIB = 1
 	IF ( STRVIB .EQ. ' 0 0 0 0 0 0+  ' ) IVIB = 2
 	IF ( STRVIB .EQ. ' 0 0 1 0 0 0+  ' ) IVIB = 3
 	IF ( STRVIB .EQ. ' 1 0 1 0 0 0+  ' ) IVIB = 4
 	IF ( STRVIB .EQ. ' 0 0 0 0 1 1  u' ) IVIB = 5
 	IF ( STRVIB .EQ. ' 0 0 0 0 0 0+ g' ) IVIB = 6
 	IF ( STRVIB .EQ. ' 0 0 0 0 3 1  u' ) IVIB = 7
 	IF ( STRVIB .EQ. ' 0 0 0 1 1 0+ u' ) IVIB = 8
 	IF ( STRVIB .EQ. ' 0 0 0 2 1 1 1u' ) IVIB = 9
 	IF ( STRVIB .EQ. ' 0 0 0 2 1 1 2u' ) IVIB = 10
 	IF ( STRVIB .EQ. ' 0 0 1 0 0 0+ u' ) IVIB = 11
 	IF ( STRVIB .EQ. ' 0 1 0 1 1 0+ u' ) IVIB = 12
 	IF ( STRVIB .EQ. ' 1 0 1 0 0 0+ u' ) IVIB = 13
 	IF ( STRVIB .EQ. ' 1 1 0 1 1 0+ u' ) IVIB = 14
 	IF ( STRVIB .EQ. ' 0 0 0 0 2 0+ g' ) IVIB = 15
 	IF ( STRVIB .EQ. ' 0 0 0 0 2 2  g' ) IVIB = 16
 	IF ( STRVIB .EQ. ' 0 0 0 0 4 0+ g' ) IVIB = 17
 	IF ( STRVIB .EQ. ' 0 0 0 0 4 2  g' ) IVIB = 18
 	IF ( STRVIB .EQ. ' 0 0 0 2 2 0+2g' ) IVIB = 19
 	IF ( STRVIB .EQ. ' 0 0 0 2 2 0- g' ) IVIB = 20
 	IF ( STRVIB .EQ. ' 0 0 0 2 2 2 2g' ) IVIB = 21
 	IF ( STRVIB .EQ. ' 0 1 0 1 0 1  g' ) IVIB = 22
 	IF ( STRVIB .EQ. ' 1 0 1 0 1 1  g' ) IVIB = 23
 	IF ( STRVIB .EQ. ' 0 0 0 1 0 1  g' ) IVIB = 24
 	IF ( STRVIB .EQ. ' 0 0 0 1 1 0- u' ) IVIB = 25
 	IF ( STRVIB .EQ. ' 0 0 0 1 1 2  u' ) IVIB = 26
 	IF ( STRVIB .EQ. ' 0 0 0 1 3 0+ u' ) IVIB = 27
 	IF ( STRVIB .EQ. ' 0 0 0 1 3 0- u' ) IVIB = 28
 	IF ( STRVIB .EQ. ' 0 0 0 1 3 2 1u' ) IVIB = 29
 	IF ( STRVIB .EQ. ' 0 0 0 1 3 2 2u' ) IVIB = 30
 	IF ( STRVIB .EQ. ' 0 0 0 3 1 0+ u' ) IVIB = 31
 	IF ( STRVIB .EQ. ' 0 0 0 3 1 0- u' ) IVIB = 32
 	IF ( STRVIB .EQ. ' 0 0 0 3 1 2 1u' ) IVIB = 33
 	IF ( STRVIB .EQ. ' 0 0 0 3 1 2 2u' ) IVIB = 34
 	IF ( STRVIB .EQ. ' 0 1 0 0 1 1  u' ) IVIB = 35
 	IF ( STRVIB .EQ. ' 1 0 1 1 0 1  u' ) IVIB = 36
C new levels for C2H2 found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. ' 0 0 0 0 2 0+  ' ) IVIB = 37
        IF ( STRVIB .EQ. ' 0 0 0 0 2 2   ' ) IVIB = 38
        IF ( STRVIB .EQ. ' 0 0 0 0 3 1   ' ) IVIB = 39
        IF ( STRVIB .EQ. ' 0 0 0 0 3 3   ' ) IVIB = 40
        IF ( STRVIB .EQ. ' 0 0 0 1 0 1   ' ) IVIB = 41
        IF ( STRVIB .EQ. ' 0 0 0 1 1 0+  ' ) IVIB = 42
        IF ( STRVIB .EQ. ' 0 0 0 1 1 0-  ' ) IVIB = 43
        IF ( STRVIB .EQ. ' 0 0 0 1 1 2   ' ) IVIB = 44
        IF ( STRVIB .EQ. ' 0 0 0 1 2 1   ' ) IVIB = 45
        IF ( STRVIB .EQ. ' 0 0 0 1 2 1 1g' ) IVIB = 46
        IF ( STRVIB .EQ. ' 0 0 0 1 2 1 2g' ) IVIB = 47
        IF ( STRVIB .EQ. ' 0 0 0 1 2 3   ' ) IVIB = 48
        IF ( STRVIB .EQ. ' 0 0 0 1 3 2  u' ) IVIB = 49
        IF ( STRVIB .EQ. ' 0 0 0 2 0 0+  ' ) IVIB = 50
        IF ( STRVIB .EQ. ' 0 0 0 2 0 0+ g' ) IVIB = 51
        IF ( STRVIB .EQ. ' 0 0 0 2 0 2   ' ) IVIB = 52
        IF ( STRVIB .EQ. ' 0 0 0 2 0 2  g' ) IVIB = 53
        IF ( STRVIB .EQ. ' 0 0 0 2 1 1   ' ) IVIB = 54
        IF ( STRVIB .EQ. ' 0 0 0 2 1 3   ' ) IVIB = 55
        IF ( STRVIB .EQ. ' 0 0 0 2 1 3  u' ) IVIB = 56
        IF ( STRVIB .EQ. ' 0 0 0 2 2 0+ g' ) IVIB = 57
        IF ( STRVIB .EQ. ' 0 0 0 3 0 1   ' ) IVIB = 58
        IF ( STRVIB .EQ. ' 0 0 0 3 0 3   ' ) IVIB = 59
        IF ( STRVIB .EQ. ' 0 0 0 3 1 2  u' ) IVIB = 60
        IF ( STRVIB .EQ. ' 0 0 1 0 1 1  g' ) IVIB = 61
        IF ( STRVIB .EQ. ' 0 0 1 0 2 0+ u' ) IVIB = 62
        IF ( STRVIB .EQ. ' 0 0 1 0 2 2  u' ) IVIB = 63
        IF ( STRVIB .EQ. ' 0 0 1 1 0 1  u' ) IVIB = 64
        IF ( STRVIB .EQ. ' 0 0 1 1 1 0+ g' ) IVIB = 65
        IF ( STRVIB .EQ. ' 0 0 1 1 1 0- g' ) IVIB = 66
        IF ( STRVIB .EQ. ' 0 0 1 1 1 2  g' ) IVIB = 67
        IF ( STRVIB .EQ. ' 0 0 1 2 0 0+ u' ) IVIB = 68
        IF ( STRVIB .EQ. ' 0 0 1 2 0 2  u' ) IVIB = 69
        IF ( STRVIB .EQ. ' 0 0 1 3 0 1  u' ) IVIB = 70
        IF ( STRVIB .EQ. ' 0 0 2 0 0 0+ g' ) IVIB = 71
        IF ( STRVIB .EQ. ' 0 0 2 0 1 1  u' ) IVIB = 72
        IF ( STRVIB .EQ. ' 0 0 2 1 0 1  g' ) IVIB = 73
        IF ( STRVIB .EQ. ' 0 0 2 1 1 0+ u' ) IVIB = 74
        IF ( STRVIB .EQ. ' 0 0 2 2 0 2  g' ) IVIB = 75
        IF ( STRVIB .EQ. ' 0 0 3 0 0 0+ u' ) IVIB = 76
        IF ( STRVIB .EQ. ' 0 0 3 1 0 1  u' ) IVIB = 77
        IF ( STRVIB .EQ. ' 0 1 0 0 0 0+ g' ) IVIB = 78
        IF ( STRVIB .EQ. ' 0 1 0 0 3 1  u' ) IVIB = 79
        IF ( STRVIB .EQ. ' 0 1 0 1 2 1 2g' ) IVIB = 80
        IF ( STRVIB .EQ. ' 0 1 0 1 3 0+2u' ) IVIB = 81
        IF ( STRVIB .EQ. ' 0 1 0 1 3 2 2u' ) IVIB = 82
        IF ( STRVIB .EQ. ' 0 1 0 2 1 1 2u' ) IVIB = 83
        IF ( STRVIB .EQ. ' 0 1 0 2 2 0- g' ) IVIB = 84
        IF ( STRVIB .EQ. ' 0 1 0 2 2 2 2g' ) IVIB = 85
        IF ( STRVIB .EQ. ' 0 1 0 2 3 1 3u' ) IVIB = 86
        IF ( STRVIB .EQ. ' 0 1 0 3 1 0+ u' ) IVIB = 87
        IF ( STRVIB .EQ. ' 0 1 0 3 1 0+2u' ) IVIB = 88
        IF ( STRVIB .EQ. ' 0 1 0 3 1 2 2u' ) IVIB = 89
        IF ( STRVIB .EQ. ' 0 1 0 4 1 1  u' ) IVIB = 90
        IF ( STRVIB .EQ. ' 0 1 0 4 1 1 2u' ) IVIB = 91
        IF ( STRVIB .EQ. ' 0 1 1 0 0 0+ u' ) IVIB = 92
        IF ( STRVIB .EQ. ' 0 1 1 0 2 0+ u' ) IVIB = 93
        IF ( STRVIB .EQ. ' 0 1 1 1 0 1  u' ) IVIB = 94
        IF ( STRVIB .EQ. ' 0 1 1 2 0 0+ u' ) IVIB = 95
        IF ( STRVIB .EQ. ' 0 2 0 1 1 0+ u' ) IVIB = 96
        IF ( STRVIB .EQ. ' 0 2 0 1 3 0+ u' ) IVIB = 97
        IF ( STRVIB .EQ. ' 0 2 0 2 1 1 2u' ) IVIB = 98
        IF ( STRVIB .EQ. ' 0 2 0 3 1 0+ u' ) IVIB = 99
        IF ( STRVIB .EQ. ' 1 0 0 0 0 0+ g' ) IVIB = 100
        IF ( STRVIB .EQ. ' 1 0 0 0 1 1  u' ) IVIB = 101
        IF ( STRVIB .EQ. ' 1 0 0 0 2 0+ g' ) IVIB = 102
        IF ( STRVIB .EQ. ' 1 0 0 0 2 2  g' ) IVIB = 103
        IF ( STRVIB .EQ. ' 1 0 0 0 3 1  u' ) IVIB = 104
        IF ( STRVIB .EQ. ' 1 0 0 1 0 1  g' ) IVIB = 105
        IF ( STRVIB .EQ. ' 1 0 0 1 1 0+ u' ) IVIB = 106
        IF ( STRVIB .EQ. ' 1 0 0 1 1 0- u' ) IVIB = 107
        IF ( STRVIB .EQ. ' 1 0 0 1 1 2  u' ) IVIB = 108
        IF ( STRVIB .EQ. ' 1 0 0 1 2 1  g' ) IVIB = 109
        IF ( STRVIB .EQ. ' 1 0 0 2 1 1  u' ) IVIB = 110
        IF ( STRVIB .EQ. ' 1 0 0 2 1 1 1u' ) IVIB = 111
        IF ( STRVIB .EQ. ' 1 0 1 0 2 0+ u' ) IVIB = 112
        IF ( STRVIB .EQ. ' 1 0 1 0 2 2  u' ) IVIB = 113
        IF ( STRVIB .EQ. ' 1 0 1 1 1 0+ g' ) IVIB = 114
        IF ( STRVIB .EQ. ' 1 0 1 1 1 0- g' ) IVIB = 115
        IF ( STRVIB .EQ. ' 1 0 1 1 1 2  g' ) IVIB = 116
        IF ( STRVIB .EQ. ' 1 0 1 2 0 0+ u' ) IVIB = 117
        IF ( STRVIB .EQ. ' 1 0 1 2 0 2  u' ) IVIB = 118
        IF ( STRVIB .EQ. ' 1 1 0 1 2 1 2g' ) IVIB = 119
        IF ( STRVIB .EQ. ' 1 1 0 2 0 0+ g' ) IVIB = 120
        IF ( STRVIB .EQ. ' 1 1 0 2 1 1  u' ) IVIB = 121
        IF ( STRVIB .EQ. ' 1 1 0 2 1 1 1u' ) IVIB = 122
        IF ( STRVIB .EQ. ' 1 1 0 2 1 1 2u' ) IVIB = 123
        IF ( STRVIB .EQ. ' 1 1 1 0 0 0+ u' ) IVIB = 124
        IF ( STRVIB .EQ. ' 1 1 1 2 0 0+ u' ) IVIB = 125
        IF ( STRVIB .EQ. ' 1 2 0 1 1 0+ u' ) IVIB = 126
        IF ( STRVIB .EQ. ' 2 0 0 0 0 0+ g' ) IVIB = 127
        IF ( STRVIB .EQ. ' 2 0 0 0 1 1  u' ) IVIB = 128
        IF ( STRVIB .EQ. ' 2 0 0 1 0 1  g' ) IVIB = 129
        IF ( STRVIB .EQ. ' 2 0 1 0 0 0+ u' ) IVIB = 130

C
C Class 8: Pyramidal tetratomic molecules
      ELSE IF ( IDXMOL .EQ. 11 .OR.                                ! NH3
     &          IDXMOL .EQ. 28      ) THEN                         ! PH3
 	IF ( STRVIB .EQ. '      0 0 0 0  ' ) IVIB = 1
 	IF ( STRVIB .EQ. '      0 1 0 0  ' ) IVIB = 2
 	IF ( STRVIB .EQ. '      0 2 0 0  ' ) IVIB = 3
 	IF ( STRVIB .EQ. '      0 0 0 1  ' ) IVIB = 4
 	IF ( STRVIB .EQ. '      0 1 0 1  ' ) IVIB = 5
 	IF ( STRVIB .EQ. '      0 0 0 2  ' ) IVIB = 6
 	IF ( STRVIB .EQ. '      0 0 1 0  ' ) IVIB = 7
 	IF ( STRVIB .EQ. '      1 0 0 0  ' ) IVIB = 8
 	IF ( STRVIB .EQ. '      0 0 0 0 a' ) IVIB = 9
 	IF ( STRVIB .EQ. '      0 1 0 0 a' ) IVIB = 10
 	IF ( STRVIB .EQ. '      0 2 0 0 a' ) IVIB = 11
 	IF ( STRVIB .EQ. '      0 0 0 1 a' ) IVIB = 12
 	IF ( STRVIB .EQ. '      0 0 0 0 s' ) IVIB = 13
 	IF ( STRVIB .EQ. '      0 1 0 0 s' ) IVIB = 14
 	IF ( STRVIB .EQ. '      0 2 0 0 s' ) IVIB = 15
 	IF ( STRVIB .EQ. '      0 0 0 1 s' ) IVIB = 16
 	IF ( STRVIB .EQ. '               ' ) IVIB = 17
 	IF ( STRVIB .EQ. '      0 3 0 0 s' ) IVIB = 18
 	IF ( STRVIB .EQ. '      0 1 0 1 s' ) IVIB = 19
 	IF ( STRVIB .EQ. '      0 1 0 1 a' ) IVIB = 20
 	IF ( STRVIB .EQ. '      0 3 0 0 a' ) IVIB = 21
 	IF ( STRVIB .EQ. '      0 2 0 1 s' ) IVIB = 22
 	IF ( STRVIB .EQ. '      0 0 0 2As' ) IVIB = 23
 	IF ( STRVIB .EQ. '      0 0 0 2Aa' ) IVIB = 24
 	IF ( STRVIB .EQ. '      0 0 0 2Es' ) IVIB = 25
 	IF ( STRVIB .EQ. '      0 0 0 2Ea' ) IVIB = 26
 	IF ( STRVIB .EQ. '      1 0 0 0 s' ) IVIB = 27
 	IF ( STRVIB .EQ. '      1 0 0 0 a' ) IVIB = 28
 	IF ( STRVIB .EQ. '      0 0 1 0 s' ) IVIB = 29
 	IF ( STRVIB .EQ. '      0 0 1 0 a' ) IVIB = 30
 	IF ( STRVIB .EQ. '      0 4 0 0 s' ) IVIB = 31
 	IF ( STRVIB .EQ. '      0 2 0 1 a' ) IVIB = 32
 	IF ( STRVIB .EQ. '      0 3 0 1 s' ) IVIB = 33
 	IF ( STRVIB .EQ. '      0 4 0 0 a' ) IVIB = 34
 	IF ( STRVIB .EQ. '      0 1 0 2As' ) IVIB = 35
 	IF ( STRVIB .EQ. '      0 1 0 2Es' ) IVIB = 36
 	IF ( STRVIB .EQ. '      0 1 0 2Aa' ) IVIB = 37
 	IF ( STRVIB .EQ. '      0 1 0 2Ea' ) IVIB = 38
 	IF ( STRVIB .EQ. '      1 1 0 0 s' ) IVIB = 39
 	IF ( STRVIB .EQ. '      1 1 0 0 a' ) IVIB = 40
 	IF ( STRVIB .EQ. '      0 1 1 0 s' ) IVIB = 41
 	IF ( STRVIB .EQ. '      0 1 1 0 a' ) IVIB = 42
 	IF ( STRVIB .EQ. '      0 3 0 1 a' ) IVIB = 43
 	IF ( STRVIB .EQ. '      1 0 0 1 s' ) IVIB = 44
 	IF ( STRVIB .EQ. '      1 0 0 1 a' ) IVIB = 45
 	IF ( STRVIB .EQ. '      0 0 1 1 s' ) IVIB = 46
 	IF ( STRVIB .EQ. '      0 0 1 1 a' ) IVIB = 47
C new levels for NH3 found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '    0000 00 0  ' ) IVIB = 48
        IF ( STRVIB .EQ. ' 0000 00 0     ' ) IVIB = 49
        IF ( STRVIB .EQ. ' 0000 00 0 A1'' ' ) IVIB = 50
        IF ( STRVIB .EQ. ' 0000 00 0 A2" ' ) IVIB = 51
        IF ( STRVIB .EQ. ' 0000 00 0 A2'' ' ) IVIB = 52
        IF ( STRVIB .EQ. ' 0001 01 1     ' ) IVIB = 53
        IF ( STRVIB .EQ. ' 0001 01 1 E"  ' ) IVIB = 54
        IF ( STRVIB .EQ. ' 0001 01 1 E''  ' ) IVIB = 55
        IF ( STRVIB .EQ. ' 0002 00 0     ' ) IVIB = 56
        IF ( STRVIB .EQ. ' 0002 00 0 A1'' ' ) IVIB = 57
        IF ( STRVIB .EQ. ' 0002 00 0 A2" ' ) IVIB = 58
        IF ( STRVIB .EQ. ' 0002 02 2     ' ) IVIB = 59
        IF ( STRVIB .EQ. ' 0002 02 2 E"  ' ) IVIB = 60
        IF ( STRVIB .EQ. ' 0002 02 2 E''  ' ) IVIB = 61
        IF ( STRVIB .EQ. ' 0003 01 1     ' ) IVIB = 62
        IF ( STRVIB .EQ. ' 0003 01 1 E"  ' ) IVIB = 63
        IF ( STRVIB .EQ. ' 0003 01 1 E''  ' ) IVIB = 64
        IF ( STRVIB .EQ. ' 0003 03 3     ' ) IVIB = 65
        IF ( STRVIB .EQ. ' 0003 03 3 A2'' ' ) IVIB = 66
        IF ( STRVIB .EQ. ' 0010 10 1     ' ) IVIB = 67
        IF ( STRVIB .EQ. ' 0010 10 1 E"  ' ) IVIB = 68
        IF ( STRVIB .EQ. ' 0010 10 1 E''  ' ) IVIB = 69
        IF ( STRVIB .EQ. ' 0011 11 0 A1" ' ) IVIB = 70
        IF ( STRVIB .EQ. ' 0011 11 2     ' ) IVIB = 71
        IF ( STRVIB .EQ. ' 0011 11 2 A1" ' ) IVIB = 72
        IF ( STRVIB .EQ. ' 0011 11 2 A1'' ' ) IVIB = 73
        IF ( STRVIB .EQ. ' 0011 11 2 A2" ' ) IVIB = 74
        IF ( STRVIB .EQ. ' 0011 11 2 A2'' ' ) IVIB = 75
        IF ( STRVIB .EQ. ' 0011 11 2 E"  ' ) IVIB = 76
        IF ( STRVIB .EQ. ' 0011 11 2 E''  ' ) IVIB = 77
        IF ( STRVIB .EQ. ' 0012 10 1 E"  ' ) IVIB = 78
        IF ( STRVIB .EQ. ' 0012 10 1 E''  ' ) IVIB = 79
        IF ( STRVIB .EQ. ' 0012 11 1 E"  ' ) IVIB = 80
        IF ( STRVIB .EQ. ' 0012 12 1 E"  ' ) IVIB = 81
        IF ( STRVIB .EQ. ' 0012 12 1 E''  ' ) IVIB = 82
        IF ( STRVIB .EQ. ' 0012 12 3     ' ) IVIB = 83
        IF ( STRVIB .EQ. ' 0020 00 0 A1'' ' ) IVIB = 84
        IF ( STRVIB .EQ. ' 0020 00 0 A2" ' ) IVIB = 85
        IF ( STRVIB .EQ. ' 0020 20 2 E"  ' ) IVIB = 86
        IF ( STRVIB .EQ. ' 0020 20 2 E''  ' ) IVIB = 87
        IF ( STRVIB .EQ. ' 0100 00 0     ' ) IVIB = 88
        IF ( STRVIB .EQ. ' 0100 00 0 A1'' ' ) IVIB = 89
        IF ( STRVIB .EQ. ' 0100 00 0 A2" ' ) IVIB = 90
        IF ( STRVIB .EQ. ' 0101 01 1     ' ) IVIB = 91
        IF ( STRVIB .EQ. ' 0101 01 1 E"  ' ) IVIB = 92
        IF ( STRVIB .EQ. ' 0101 01 1 E''  ' ) IVIB = 93
        IF ( STRVIB .EQ. ' 0102 00 0 A1'' ' ) IVIB = 94
        IF ( STRVIB .EQ. ' 0102 02 2 E''  ' ) IVIB = 95
        IF ( STRVIB .EQ. ' 0110 10 1     ' ) IVIB = 96
        IF ( STRVIB .EQ. ' 0110 10 1 E"  ' ) IVIB = 97
        IF ( STRVIB .EQ. ' 0110 10 1 E''  ' ) IVIB = 98
        IF ( STRVIB .EQ. ' 0200 00 0     ' ) IVIB = 99
        IF ( STRVIB .EQ. ' 0200 00 0 A1'' ' ) IVIB = 100
        IF ( STRVIB .EQ. ' 0200 00 0 A2" ' ) IVIB = 101
        IF ( STRVIB .EQ. ' 0201 01 1 E"  ' ) IVIB = 102
        IF ( STRVIB .EQ. ' 0201 01 1 E''  ' ) IVIB = 103
        IF ( STRVIB .EQ. ' 0202 00 0     ' ) IVIB = 104
        IF ( STRVIB .EQ. ' 0202 00 0 A1'' ' ) IVIB = 105
        IF ( STRVIB .EQ. ' 0202 00 0 A2" ' ) IVIB = 106
        IF ( STRVIB .EQ. ' 0202 02 2     ' ) IVIB = 107
        IF ( STRVIB .EQ. ' 0202 02 2 E"  ' ) IVIB = 108
        IF ( STRVIB .EQ. ' 0203 01 1 E"  ' ) IVIB = 109
        IF ( STRVIB .EQ. ' 0203 01 1 E''  ' ) IVIB = 110
        IF ( STRVIB .EQ. ' 0210 10 1     ' ) IVIB = 111
        IF ( STRVIB .EQ. ' 0210 10 1 E''  ' ) IVIB = 112
        IF ( STRVIB .EQ. ' 0300 00 0     ' ) IVIB = 113
        IF ( STRVIB .EQ. ' 0300 00 0 A1'' ' ) IVIB = 114
        IF ( STRVIB .EQ. ' 0300 00 0 A2" ' ) IVIB = 115
        IF ( STRVIB .EQ. ' 0301 01 1 E"  ' ) IVIB = 116
        IF ( STRVIB .EQ. ' 0400 00 0 A1'' ' ) IVIB = 117
        IF ( STRVIB .EQ. ' 0401 01 1     ' ) IVIB = 118
        IF ( STRVIB .EQ. ' 0401 01 1 E''  ' ) IVIB = 119
        IF ( STRVIB .EQ. ' 1000 00 0     ' ) IVIB = 120
        IF ( STRVIB .EQ. ' 1000 00 0 A1'' ' ) IVIB = 121
        IF ( STRVIB .EQ. ' 1000 00 0 A2" ' ) IVIB = 122
        IF ( STRVIB .EQ. ' 1001 01 1     ' ) IVIB = 123
        IF ( STRVIB .EQ. ' 1001 01 1 E"  ' ) IVIB = 124
        IF ( STRVIB .EQ. ' 1001 01 1 E''  ' ) IVIB = 125
        IF ( STRVIB .EQ. ' 1002 02 2 E"  ' ) IVIB = 126
        IF ( STRVIB .EQ. ' 1002 02 2 E''  ' ) IVIB = 127
        IF ( STRVIB .EQ. ' 1010 10 1 E"  ' ) IVIB = 128
        IF ( STRVIB .EQ. ' 1010 10 1 E''  ' ) IVIB = 129
        IF ( STRVIB .EQ. ' 1100 00 0     ' ) IVIB = 130
        IF ( STRVIB .EQ. ' 1100 00 0 A1'' ' ) IVIB = 131
        IF ( STRVIB .EQ. ' 1100 00 0 A2" ' ) IVIB = 132
        IF ( STRVIB .EQ. ' 1101 01 1 E''  ' ) IVIB = 133
        IF ( STRVIB .EQ. ' 1200 00 0     ' ) IVIB = 134
        IF ( STRVIB .EQ. ' 1200 00 0 A1'' ' ) IVIB = 135
C new levels for PH3 found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '      0 0 1 1  ' ) IVIB = 136
        IF ( STRVIB .EQ. '      0 1 0 2  ' ) IVIB = 137
        IF ( STRVIB .EQ. '      0 1 1 0  ' ) IVIB = 138
        IF ( STRVIB .EQ. '      0 2 0 1  ' ) IVIB = 139
        IF ( STRVIB .EQ. '      0 3 0 0  ' ) IVIB = 140
        IF ( STRVIB .EQ. '      0 4 0 0  ' ) IVIB = 141
        IF ( STRVIB .EQ. '      1 0 0 1  ' ) IVIB = 142
        IF ( STRVIB .EQ. '      1 1 0 0  ' ) IVIB = 143
C
C Class 9: Non-linear tetratomic molecules
      ELSE IF ( IDXMOL .EQ. 20 .OR.                                  ! H2CO
     &          IDXMOL .EQ. 25 .OR.                                  ! H2O2
     &          IDXMOL .EQ. 29      ) THEN                           ! COF2
 	IF ( STRVIB .EQ. '    0 0 0 0 0 0' ) IVIB = 1
 	IF ( STRVIB .EQ. '    0 0 0 0 0 2' ) IVIB = 2
 	IF ( STRVIB .EQ. '    0 0 1 1 0 0' ) IVIB = 3
 	IF ( STRVIB .EQ. '    0 0 1 0 0 1' ) IVIB = 4
 	IF ( STRVIB .EQ. '    1 0 0 0 0 0' ) IVIB = 5
 	IF ( STRVIB .EQ. '    0 0 0 0 1 0' ) IVIB = 6
 	IF ( STRVIB .EQ. '    0 1 0 1 0 0' ) IVIB = 7
 	IF ( STRVIB .EQ. '    0 1 0 0 0 1' ) IVIB = 8
 	IF ( STRVIB .EQ. '    0 0 0 0 0 1' ) IVIB = 9
 	IF ( STRVIB .EQ. '    0 1 0 0 0 0' ) IVIB = 10
 	IF ( STRVIB .EQ. '    0 0 0 1 0 0' ) IVIB = 11
 	IF ( STRVIB .EQ. '    0 2 0 0 0 0' ) IVIB = 12
 	IF ( STRVIB .EQ. '               ' ) IVIB = 13
 	IF ( STRVIB .EQ. '    0 0 2 0 0 1' ) IVIB = 14
 	IF ( STRVIB .EQ. '    0 0 0 0 2 0' ) IVIB = 15
 	IF ( STRVIB .EQ. '    0 0 000 0 0' ) IVIB = 16
 	IF ( STRVIB .EQ. '    0 0 000 0 1' ) IVIB = 17
 	IF ( STRVIB .EQ. '    0 0 001 0 0' ) IVIB = 18
 	IF ( STRVIB .EQ. '    0 0 002 0 0' ) IVIB = 19
 	IF ( STRVIB .EQ. '    0 0 011 0 0' ) IVIB = 20
 	IF ( STRVIB .EQ. '    0 0 012 0 0' ) IVIB = 21
 	IF ( STRVIB .EQ. '    0 0 021 0 0' ) IVIB = 22
 	IF ( STRVIB .EQ. '    0 0 022 0 0' ) IVIB = 23
 	IF ( STRVIB .EQ. '    0 0 031 0 0' ) IVIB = 24
 	IF ( STRVIB .EQ. '    0 0 032 0 0' ) IVIB = 25
 	IF ( STRVIB .EQ. '    0 0 101 0 0' ) IVIB = 26
 	IF ( STRVIB .EQ. '    0 0 102 0 0' ) IVIB = 27
 	IF ( STRVIB .EQ. '    0 0 111 0 0' ) IVIB = 28
 	IF ( STRVIB .EQ. '    0 0 112 0 0' ) IVIB = 29
 	IF ( STRVIB .EQ. '    0 0 003 0 0' ) IVIB = 30
 	IF ( STRVIB .EQ. '    0 0 004 0 0' ) IVIB = 31
 	IF ( STRVIB .EQ. '    0 0 013 0 0' ) IVIB = 32
 	IF ( STRVIB .EQ. '    0 0 014 0 0' ) IVIB = 33
 	IF ( STRVIB .EQ. '    0 0 023 0 0' ) IVIB = 34
 	IF ( STRVIB .EQ. '    0 0 024 0 0' ) IVIB = 35
 	IF ( STRVIB .EQ. '    0 0 033 0 0' ) IVIB = 36
 	IF ( STRVIB .EQ. '    0 0 034 0 0' ) IVIB = 37
 	IF ( STRVIB .EQ. '    0 0 103 0 0' ) IVIB = 38
 	IF ( STRVIB .EQ. '    0 0 104 0 0' ) IVIB = 39
C new levels for H2O2 found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '    0 0 001 0 1' ) IVIB = 40
        IF ( STRVIB .EQ. '    0 0 002 0 1' ) IVIB = 41
        IF ( STRVIB .EQ. '    0 0 003 0 1' ) IVIB = 42
        IF ( STRVIB .EQ. '    0 0 004 0 1' ) IVIB = 43
        IF ( STRVIB .EQ. '    0 0 011 0 1' ) IVIB = 44
        IF ( STRVIB .EQ. '    0 0 012 0 1' ) IVIB = 45
        IF ( STRVIB .EQ. '    0 0 013 0 1' ) IVIB = 46
        IF ( STRVIB .EQ. '    0 0 014 0 1' ) IVIB = 47
        IF ( STRVIB .EQ. '    0 0 021 0 1' ) IVIB = 48
        IF ( STRVIB .EQ. '    0 0 022 0 1' ) IVIB = 49
        IF ( STRVIB .EQ. '    0 0 023 0 1' ) IVIB = 50
        IF ( STRVIB .EQ. '    0 0 024 0 1' ) IVIB = 51
        IF ( STRVIB .EQ. '    0 0 113 0 0' ) IVIB = 52
        IF ( STRVIB .EQ. '    0 0 114 0 0' ) IVIB = 53
        IF ( STRVIB .EQ. '    0 0 131 0 0' ) IVIB = 54
        IF ( STRVIB .EQ. '    0 0 132 0 0' ) IVIB = 55
        IF ( STRVIB .EQ. '    0 0 133 0 0' ) IVIB = 56
        IF ( STRVIB .EQ. '    0 0 134 0 0' ) IVIB = 57
        IF ( STRVIB .EQ. '    0 1 001 0 0' ) IVIB = 58
        IF ( STRVIB .EQ. '    0 1 002 0 0' ) IVIB = 59
        IF ( STRVIB .EQ. '    0 1 003 0 0' ) IVIB = 60
        IF ( STRVIB .EQ. '    0 1 004 0 0' ) IVIB = 61
        IF ( STRVIB .EQ. '    0 1 011 0 0' ) IVIB = 62
        IF ( STRVIB .EQ. '    0 1 012 0 0' ) IVIB = 63
        IF ( STRVIB .EQ. '    0 1 013 0 0' ) IVIB = 64
        IF ( STRVIB .EQ. '    0 1 014 0 0' ) IVIB = 65
C new levels for COF2 found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '    0 0 0 0 1 1' ) IVIB = 66
        IF ( STRVIB .EQ. '    0 0 1 0 0 0' ) IVIB = 67
C new levels for H2CO found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '    0 0 0 1 0 1' ) IVIB = 68
        IF ( STRVIB .EQ. '    0 0 0 2 0 0' ) IVIB = 69
        IF ( STRVIB .EQ. '    0 0 2 0 0 0' ) IVIB = 70
        IF ( STRVIB .EQ. '    0 1 1 0 0 0' ) IVIB = 71
C
C Class 10: Pentatomic or greater polyatomic molecules
      ELSE IF ( IDXMOL .EQ.  6 .OR.                                ! CH4
     &          IDXMOL .EQ. 12 .OR.                                ! HNO3
     &          IDXMOL .EQ. 24 .OR.                                ! CH3Cl
     &          IDXMOL .EQ. 27 .OR.                                ! C2H6
     &          IDXMOL .EQ. 30 .OR.                                ! SF6
     &          IDXMOL .EQ. 32 .OR.                                ! HCOOH
     &          IDXMOL .EQ. 35 .OR.                                ! ClONO2
     &          IDXMOL .EQ. 38 .OR.                                ! C2H4
     &          IDXMOL .EQ. 39 .OR.                                ! CH3OH
     &          IDXMOL .EQ. 40 .OR.                                ! CH3Br
     &          IDXMOL .EQ. 41 .OR.                                ! CH3CN
     &          IDXMOL .EQ. 42      ) THEN                         ! CF4
 	IF ( STRVIB .EQ. '         GROUND' ) IVIB = 1
 	IF ( STRVIB .EQ. '             V1' ) IVIB = 2
 	IF ( STRVIB .EQ. '             V2' ) IVIB = 3
 	IF ( STRVIB .EQ. '             V4' ) IVIB = 4
 	IF ( STRVIB .EQ. '             V5' ) IVIB = 5
 	IF ( STRVIB .EQ. '             V9' ) IVIB = 6
 	IF ( STRVIB .EQ. '            2V5' ) IVIB = 7
 	IF ( STRVIB .EQ. '            2V9' ) IVIB = 8
 	IF ( STRVIB .EQ. '            3V6' ) IVIB = 9
 	IF ( STRVIB .EQ. '            3V9' ) IVIB = 10
 	IF ( STRVIB .EQ. '          V5+V9' ) IVIB = 11
 	IF ( STRVIB .EQ. '               ' ) IVIB = 12
 	IF ( STRVIB .EQ. '             V6' ) IVIB = 13
 	IF ( STRVIB .EQ. '             V3' ) IVIB = 14
 	IF ( STRVIB .EQ. '            2V6' ) IVIB = 15
 	IF ( STRVIB .EQ. '             V7' ) IVIB = 16
 	IF ( STRVIB .EQ. '             V8' ) IVIB = 17
 	IF ( STRVIB .EQ. '          V8+V9' ) IVIB = 18
 	IF ( STRVIB .EQ. '          V3+V6' ) IVIB = 19
 	IF ( STRVIB .EQ. '            2V3' ) IVIB = 20
 	IF ( STRVIB .EQ. '          V5+V6' ) IVIB = 21
 	IF ( STRVIB .EQ. '          V3+V5' ) IVIB = 22
 	IF ( STRVIB .EQ. '          V4+V9' ) IVIB = 23
 	IF ( STRVIB .EQ. '            V10' ) IVIB = 24
 	IF ( STRVIB .EQ. '            V11' ) IVIB = 25
 	IF ( STRVIB .EQ. '         V2+V12' ) IVIB = 26
 	IF ( STRVIB .EQ. '       2V10+V12' ) IVIB = 27
 	IF ( STRVIB .EQ. '         V9+V10' ) IVIB = 28
 	IF ( STRVIB .EQ. '            V12' ) IVIB = 29
 	IF ( STRVIB .EQ. '           2V12' ) IVIB = 30
 	IF ( STRVIB .EQ. '         V6+V12' ) IVIB = 31
 	IF ( STRVIB .EQ. '         V7+V12' ) IVIB = 32
 	IF ( STRVIB .EQ. '         V8+V12' ) IVIB = 33
 	IF ( STRVIB .EQ. '        V8+2V12' ) IVIB = 34
 	IF ( STRVIB .EQ. '           3V12' ) IVIB = 35
 	IF ( STRVIB .EQ. '           4V12' ) IVIB = 36
 	IF ( STRVIB .EQ. '    0 0 0 0    ' ) IVIB = 37
 	IF ( STRVIB .EQ. '    0 0 0 0 1A1' ) IVIB = 38
 	IF ( STRVIB .EQ. '    0 0 0 1 1F2' ) IVIB = 39
 	IF ( STRVIB .EQ. '    0 0 0 2 1 E' ) IVIB = 40
 	IF ( STRVIB .EQ. '    0 0 0 2 1A1' ) IVIB = 41
 	IF ( STRVIB .EQ. '    0 0 0 2 1F2' ) IVIB = 42
 	IF ( STRVIB .EQ. '    0 0 0 3 1A1' ) IVIB = 43
 	IF ( STRVIB .EQ. '    0 0 0 3 1F1' ) IVIB = 44
 	IF ( STRVIB .EQ. '    0 0 0 3 1F2' ) IVIB = 45
 	IF ( STRVIB .EQ. '    0 0 0 3 2F2' ) IVIB = 46
 	IF ( STRVIB .EQ. '    0 0 0 4    ' ) IVIB = 47
 	IF ( STRVIB .EQ. '    0 0 1 0 1F2' ) IVIB = 48
 	IF ( STRVIB .EQ. '    0 0 1 1 1 E' ) IVIB = 49
 	IF ( STRVIB .EQ. '    0 0 1 1 1A1' ) IVIB = 50
 	IF ( STRVIB .EQ. '    0 0 1 1 1F1' ) IVIB = 51
 	IF ( STRVIB .EQ. '    0 0 1 1 1F2' ) IVIB = 52
 	IF ( STRVIB .EQ. '    0 0 1 2    ' ) IVIB = 53
 	IF ( STRVIB .EQ. '    0 0 1 2 1F2' ) IVIB = 54
 	IF ( STRVIB .EQ. '    0 0 2 0 1F2' ) IVIB = 55
 	IF ( STRVIB .EQ. '    0 0 3 0    ' ) IVIB = 56
 	IF ( STRVIB .EQ. '    0 1 0 0 1 E' ) IVIB = 57
 	IF ( STRVIB .EQ. '    0 1 0 1 1F1' ) IVIB = 58
 	IF ( STRVIB .EQ. '    0 1 0 1 1F2' ) IVIB = 59
 	IF ( STRVIB .EQ. '    0 1 0 2 1 E' ) IVIB = 60
 	IF ( STRVIB .EQ. '    0 1 0 2 1A1' ) IVIB = 61
 	IF ( STRVIB .EQ. '    0 1 0 2 1A2' ) IVIB = 62
 	IF ( STRVIB .EQ. '    0 1 0 2 1F1' ) IVIB = 63
 	IF ( STRVIB .EQ. '    0 1 0 2 1F2' ) IVIB = 64
 	IF ( STRVIB .EQ. '    0 1 0 2 2 E' ) IVIB = 65
 	IF ( STRVIB .EQ. '    0 1 1 0 1F1' ) IVIB = 66
 	IF ( STRVIB .EQ. '    0 1 1 0 1F2' ) IVIB = 67
 	IF ( STRVIB .EQ. '    0 1 2 0    ' ) IVIB = 68
 	IF ( STRVIB .EQ. '    0 2 0 0 1 E' ) IVIB = 69
 	IF ( STRVIB .EQ. '    0 2 0 0 1A1' ) IVIB = 70
 	IF ( STRVIB .EQ. '    0 2 0 1 1F1' ) IVIB = 71
 	IF ( STRVIB .EQ. '    0 2 0 1 1F2' ) IVIB = 72
 	IF ( STRVIB .EQ. '    0 2 0 1 2F2' ) IVIB = 73
 	IF ( STRVIB .EQ. '    0 3 0 0    ' ) IVIB = 74
 	IF ( STRVIB .EQ. '    0 3 0 0 1 E' ) IVIB = 75
 	IF ( STRVIB .EQ. '    0 3 0 0 1A1' ) IVIB = 76
 	IF ( STRVIB .EQ. '    0 3 0 0 1A2' ) IVIB = 77
 	IF ( STRVIB .EQ. '    1 0 0 0 1A1' ) IVIB = 78
 	IF ( STRVIB .EQ. '    1 0 0 1 1F2' ) IVIB = 79
 	IF ( STRVIB .EQ. '    1 1 0 0 1 E' ) IVIB = 80
 	IF ( STRVIB .EQ. '    0 0 0 2 1 E' ) IVIB = 81
 	IF ( STRVIB .EQ. '    0 0 0 2 1A1' ) IVIB = 82
 	IF ( STRVIB .EQ. '    0 0 0 2 1F2' ) IVIB = 83
 	IF ( STRVIB .EQ. '    0 0 0 3 1A1' ) IVIB = 84
 	IF ( STRVIB .EQ. '    0 0 0 3 1F1' ) IVIB = 85
 	IF ( STRVIB .EQ. '    0 0 0 3 1F2' ) IVIB = 86
 	IF ( STRVIB .EQ. '    0 0 0 3 2F2' ) IVIB = 87
 	IF ( STRVIB .EQ. '    0 0 1 1 1 E' ) IVIB = 88
 	IF ( STRVIB .EQ. '    0 0 1 1 1A1' ) IVIB = 89
 	IF ( STRVIB .EQ. '    0 0 1 1 1F1' ) IVIB = 90
 	IF ( STRVIB .EQ. '    0 0 1 1 1F2' ) IVIB = 91
C new levels for CH4 found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '            3V2' ) IVIB = 92
        IF ( STRVIB .EQ. '            3V3' ) IVIB = 93
        IF ( STRVIB .EQ. '          V2+V3' ) IVIB = 94
        IF ( STRVIB .EQ. '          V2+V5' ) IVIB = 95
        IF ( STRVIB .EQ. '          V2+V6' ) IVIB = 96
        IF ( STRVIB .EQ. '         2V2 E ' ) IVIB = 97
        IF ( STRVIB .EQ. '         2V3+V6' ) IVIB = 98
        IF ( STRVIB .EQ. '         3V5 A1' ) IVIB = 99
        IF ( STRVIB .EQ. '         3V5 E ' ) IVIB = 100
        IF ( STRVIB .EQ. '         V3+2V6' ) IVIB = 101
        IF ( STRVIB .EQ. '         V5+2V6' ) IVIB = 102
        IF ( STRVIB .EQ. '       V1+V2+V6' ) IVIB = 103
        IF ( STRVIB .EQ. '       V1+V5 E ' ) IVIB = 104
        IF ( STRVIB .EQ. '       V1+V6 E ' ) IVIB = 105
        IF ( STRVIB .EQ. '       V2+V4+V6' ) IVIB = 106
        IF ( STRVIB .EQ. '       V3+V4 E ' ) IVIB = 107
        IF ( STRVIB .EQ. '       V3+V5+V6' ) IVIB = 108
        IF ( STRVIB .EQ. '       V4+V5 A1' ) IVIB = 109
        IF ( STRVIB .EQ. '       V4+V5 A2' ) IVIB = 110
        IF ( STRVIB .EQ. '       V4+V5 E ' ) IVIB = 111
        IF ( STRVIB .EQ. '       V4+V6 A1' ) IVIB = 112
        IF ( STRVIB .EQ. '       V4+V6 A2' ) IVIB = 113
        IF ( STRVIB .EQ. '       V4+V6 E ' ) IVIB = 114
        IF ( STRVIB .EQ. '      2V3+V5 E ' ) IVIB = 115
        IF ( STRVIB .EQ. '      2V5+V6 A1' ) IVIB = 116
        IF ( STRVIB .EQ. '      2V5+V6 A2' ) IVIB = 117
        IF ( STRVIB .EQ. '      2V5+V6 E ' ) IVIB = 118
        IF ( STRVIB .EQ. '      V2+2V5+V6' ) IVIB = 119
        IF ( STRVIB .EQ. '      V2+2V6 A1' ) IVIB = 120
        IF ( STRVIB .EQ. '      V2+2V6 E ' ) IVIB = 121
        IF ( STRVIB .EQ. '      V3+2V5 A1' ) IVIB = 122
        IF ( STRVIB .EQ. '      V3+2V5 E ' ) IVIB = 123
        IF ( STRVIB .EQ. '    0 0 0 2 1E ' ) IVIB = 124
        IF ( STRVIB .EQ. '    0 0 0 5  A1' ) IVIB = 125
        IF ( STRVIB .EQ. '    0 0 0 5  E ' ) IVIB = 126
        IF ( STRVIB .EQ. '    0 0 0 5  F1' ) IVIB = 127
        IF ( STRVIB .EQ. '    0 0 0 5  F2' ) IVIB = 128
        IF ( STRVIB .EQ. '    0 0 1 1 1E ' ) IVIB = 129
        IF ( STRVIB .EQ. '    0 0 1 2 1A1' ) IVIB = 130
        IF ( STRVIB .EQ. '    0 0 1 2 1E ' ) IVIB = 131
        IF ( STRVIB .EQ. '    0 0 1 2 1F1' ) IVIB = 132
        IF ( STRVIB .EQ. '    0 0 2 0  A1' ) IVIB = 133
        IF ( STRVIB .EQ. '    0 0 2 0  E ' ) IVIB = 134
        IF ( STRVIB .EQ. '    0 0 2 0  F2' ) IVIB = 135
        IF ( STRVIB .EQ. '    0 1 0 0 1E ' ) IVIB = 136
        IF ( STRVIB .EQ. '    0 1 0 2 1E ' ) IVIB = 137
        IF ( STRVIB .EQ. '    0 1 0 2 2E ' ) IVIB = 138
        IF ( STRVIB .EQ. '    0 1 0 3 1F1' ) IVIB = 139
        IF ( STRVIB .EQ. '    0 1 0 4  A1' ) IVIB = 140
        IF ( STRVIB .EQ. '    0 1 0 4  A2' ) IVIB = 141
        IF ( STRVIB .EQ. '    0 1 0 4  E ' ) IVIB = 142
        IF ( STRVIB .EQ. '    0 1 0 4  F1' ) IVIB = 143
        IF ( STRVIB .EQ. '    0 1 0 4  F2' ) IVIB = 144
        IF ( STRVIB .EQ. '    0 1 1 1    ' ) IVIB = 145
        IF ( STRVIB .EQ. '    0 1 1 1  A1' ) IVIB = 146
        IF ( STRVIB .EQ. '    0 1 1 1  A2' ) IVIB = 147
        IF ( STRVIB .EQ. '    0 1 1 1  E ' ) IVIB = 148
        IF ( STRVIB .EQ. '    0 1 1 1  F1' ) IVIB = 149
        IF ( STRVIB .EQ. '    0 1 1 1  F2' ) IVIB = 150
        IF ( STRVIB .EQ. '    0 1 1 1 1A1' ) IVIB = 151
        IF ( STRVIB .EQ. '    0 1 1 1 1A2' ) IVIB = 152
        IF ( STRVIB .EQ. '    0 1 1 1 1E ' ) IVIB = 153
        IF ( STRVIB .EQ. '    0 1 1 1 1F1' ) IVIB = 154
        IF ( STRVIB .EQ. '    0 1 1 1 1F2' ) IVIB = 155
        IF ( STRVIB .EQ. '    0 2 0 0 1E ' ) IVIB = 156
        IF ( STRVIB .EQ. '    0 2 0 2 1A1' ) IVIB = 157
        IF ( STRVIB .EQ. '    0 2 0 2 1E ' ) IVIB = 158
        IF ( STRVIB .EQ. '    0 2 0 2 1F1' ) IVIB = 159
        IF ( STRVIB .EQ. '    0 2 0 2 1F2' ) IVIB = 160
        IF ( STRVIB .EQ. '    0 2 1 0    ' ) IVIB = 161
        IF ( STRVIB .EQ. '    0 2 1 0  F1' ) IVIB = 162
        IF ( STRVIB .EQ. '    0 2 1 0  F2' ) IVIB = 163
        IF ( STRVIB .EQ. '    0 3 0 0 1E ' ) IVIB = 164
        IF ( STRVIB .EQ. '    0 3 0 1  F1' ) IVIB = 165
        IF ( STRVIB .EQ. '    0 3 0 1  F2' ) IVIB = 166
        IF ( STRVIB .EQ. '    0 3 0 1 1F1' ) IVIB = 167
        IF ( STRVIB .EQ. '    0 3 0 1 1F2' ) IVIB = 168
        IF ( STRVIB .EQ. '    1 0 0 2 1F2' ) IVIB = 169
        IF ( STRVIB .EQ. '    1 0 0 3  A1' ) IVIB = 170
        IF ( STRVIB .EQ. '    1 0 0 3  F2' ) IVIB = 171
        IF ( STRVIB .EQ. '    1 0 1 0  F2' ) IVIB = 172
        IF ( STRVIB .EQ. '    1 0 1 0 1F2' ) IVIB = 173
        IF ( STRVIB .EQ. '    1 1 0 0 1E ' ) IVIB = 174
        IF ( STRVIB .EQ. '    1 1 0 1 1F2' ) IVIB = 175
        IF ( STRVIB .EQ. '    1 2 0 0  E ' ) IVIB = 176
        IF ( STRVIB .EQ. '    1 2 0 0  F1' ) IVIB = 177
        IF ( STRVIB .EQ. '    2 0 0 0  A1' ) IVIB = 178
        IF ( STRVIB .EQ. '    2 0 0 0 1A1' ) IVIB = 179
C new levels for HNO3 found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '          V5+V7' ) IVIB = 180
        IF ( STRVIB .EQ. '          V6+V7' ) IVIB = 181
C new levels for C2H6 found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '            2V4' ) IVIB = 182
        IF ( STRVIB .EQ. '            3V4' ) IVIB = 183
        IF ( STRVIB .EQ. '          V4+V6' ) IVIB = 184
        IF ( STRVIB .EQ. '          V4+V8' ) IVIB = 185
        IF ( STRVIB .EQ. '         2V4+V9' ) IVIB = 186
        IF ( STRVIB .EQ. '         V4+V12' ) IVIB = 187
        IF ( STRVIB .EQ. '        2V4+V12' ) IVIB = 188
C new level for HCOOH found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '          V6+V9' ) IVIB = 189
C new level for CH3Cl found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '         2V3+V5' ) IVIB = 190
      ELSE IF ( IDXMOL .EQ. 43 ) THEN                               ! C4H2
C new levels found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. ' 000000000 e g ' ) IVIB = 1
        IF ( STRVIB .EQ. ' 000000000 e u ' ) IVIB = 2
        IF ( STRVIB .EQ. ' 000000001 e g ' ) IVIB = 3
        IF ( STRVIB .EQ. ' 000000001 e u ' ) IVIB = 4
        IF ( STRVIB .EQ. ' 000000001 f g ' ) IVIB = 5
        IF ( STRVIB .EQ. ' 000000001 f u ' ) IVIB = 6
        IF ( STRVIB .EQ. ' 000000002 e g ' ) IVIB = 7
        IF ( STRVIB .EQ. ' 000000002 e u ' ) IVIB = 8
        IF ( STRVIB .EQ. ' 000000002 f g ' ) IVIB = 9
        IF ( STRVIB .EQ. ' 000000002 f u ' ) IVIB = 10
        IF ( STRVIB .EQ. ' 000000003 e g ' ) IVIB = 11
        IF ( STRVIB .EQ. ' 000000003 e u ' ) IVIB = 12
        IF ( STRVIB .EQ. ' 000000003 f g ' ) IVIB = 13
        IF ( STRVIB .EQ. ' 000000003 f u ' ) IVIB = 14
        IF ( STRVIB .EQ. ' 000000004 e g ' ) IVIB = 15
        IF ( STRVIB .EQ. ' 000000004 e u ' ) IVIB = 16
        IF ( STRVIB .EQ. ' 000000004 f g ' ) IVIB = 17
        IF ( STRVIB .EQ. ' 000000004 f u ' ) IVIB = 18
        IF ( STRVIB .EQ. ' 000000005 e g ' ) IVIB = 19
        IF ( STRVIB .EQ. ' 000000005 e u ' ) IVIB = 20
        IF ( STRVIB .EQ. ' 000000005 f g ' ) IVIB = 21
        IF ( STRVIB .EQ. ' 000000005 f u ' ) IVIB = 22
        IF ( STRVIB .EQ. ' 000000006 e g ' ) IVIB = 23
        IF ( STRVIB .EQ. ' 000000006 e u ' ) IVIB = 24
        IF ( STRVIB .EQ. ' 000000006 f g ' ) IVIB = 25
        IF ( STRVIB .EQ. ' 000000006 f u ' ) IVIB = 26
        IF ( STRVIB .EQ. ' 000000007 e g ' ) IVIB = 27
        IF ( STRVIB .EQ. ' 000000007 e u ' ) IVIB = 28
        IF ( STRVIB .EQ. ' 000000007 f g ' ) IVIB = 29
        IF ( STRVIB .EQ. ' 000000007 f u ' ) IVIB = 30
        IF ( STRVIB .EQ. ' 000000008 e g ' ) IVIB = 31
        IF ( STRVIB .EQ. ' 000000008 e u ' ) IVIB = 32
        IF ( STRVIB .EQ. ' 000000008 f u ' ) IVIB = 33
        IF ( STRVIB .EQ. ' 000000010 e g ' ) IVIB = 34
        IF ( STRVIB .EQ. ' 000000010 e u ' ) IVIB = 35
        IF ( STRVIB .EQ. ' 000000010 f g ' ) IVIB = 36
        IF ( STRVIB .EQ. ' 000000010 f u ' ) IVIB = 37
        IF ( STRVIB .EQ. ' 000000011 e g ' ) IVIB = 38
        IF ( STRVIB .EQ. ' 000000011 e u ' ) IVIB = 39
        IF ( STRVIB .EQ. ' 000000011 f g ' ) IVIB = 40
        IF ( STRVIB .EQ. ' 000000011 f u ' ) IVIB = 41
        IF ( STRVIB .EQ. ' 000000012 e g ' ) IVIB = 42
        IF ( STRVIB .EQ. ' 000000012 e u ' ) IVIB = 43
        IF ( STRVIB .EQ. ' 000000012 f g ' ) IVIB = 44
        IF ( STRVIB .EQ. ' 000000012 f u ' ) IVIB = 45
        IF ( STRVIB .EQ. ' 000000013 e g ' ) IVIB = 46
        IF ( STRVIB .EQ. ' 000000013 e u ' ) IVIB = 47
        IF ( STRVIB .EQ. ' 000000013 f g ' ) IVIB = 48
        IF ( STRVIB .EQ. ' 000000013 f u ' ) IVIB = 49
        IF ( STRVIB .EQ. ' 000000014 e g ' ) IVIB = 50
        IF ( STRVIB .EQ. ' 000000014 e u ' ) IVIB = 51
        IF ( STRVIB .EQ. ' 000000014 f g ' ) IVIB = 52
        IF ( STRVIB .EQ. ' 000000014 f u ' ) IVIB = 53
        IF ( STRVIB .EQ. ' 000000015 e g ' ) IVIB = 54
        IF ( STRVIB .EQ. ' 000000015 e u ' ) IVIB = 55
        IF ( STRVIB .EQ. ' 000000015 f g ' ) IVIB = 56
        IF ( STRVIB .EQ. ' 000000015 f u ' ) IVIB = 57
        IF ( STRVIB .EQ. ' 000000016 e g ' ) IVIB = 58
        IF ( STRVIB .EQ. ' 000000016 e u ' ) IVIB = 59
        IF ( STRVIB .EQ. ' 000000016 f g ' ) IVIB = 60
        IF ( STRVIB .EQ. ' 000000016 f u ' ) IVIB = 61
        IF ( STRVIB .EQ. ' 000000020 e g ' ) IVIB = 62
        IF ( STRVIB .EQ. ' 000000020 e u ' ) IVIB = 63
        IF ( STRVIB .EQ. ' 000000020 f g ' ) IVIB = 64
        IF ( STRVIB .EQ. ' 000000020 f u ' ) IVIB = 65
        IF ( STRVIB .EQ. ' 000000021 e g ' ) IVIB = 66
        IF ( STRVIB .EQ. ' 000000021 e u ' ) IVIB = 67
        IF ( STRVIB .EQ. ' 000000021 f g ' ) IVIB = 68
        IF ( STRVIB .EQ. ' 000000021 f u ' ) IVIB = 69
        IF ( STRVIB .EQ. ' 000000022 e g ' ) IVIB = 70
        IF ( STRVIB .EQ. ' 000000022 e u ' ) IVIB = 71
        IF ( STRVIB .EQ. ' 000000022 f g ' ) IVIB = 72
        IF ( STRVIB .EQ. ' 000000022 f u ' ) IVIB = 73
        IF ( STRVIB .EQ. ' 000000023 e g ' ) IVIB = 74
        IF ( STRVIB .EQ. ' 000000023 e u ' ) IVIB = 75
        IF ( STRVIB .EQ. ' 000000023 f g ' ) IVIB = 76
        IF ( STRVIB .EQ. ' 000000023 f u ' ) IVIB = 77
        IF ( STRVIB .EQ. ' 000000030 e g ' ) IVIB = 78
        IF ( STRVIB .EQ. ' 000000030 e u ' ) IVIB = 79
        IF ( STRVIB .EQ. ' 000000030 f g ' ) IVIB = 80
        IF ( STRVIB .EQ. ' 000000030 f u ' ) IVIB = 81
        IF ( STRVIB .EQ. ' 000000100 e g ' ) IVIB = 82
        IF ( STRVIB .EQ. ' 000000100 e u ' ) IVIB = 83
        IF ( STRVIB .EQ. ' 000000100 f g ' ) IVIB = 84
        IF ( STRVIB .EQ. ' 000000100 f u ' ) IVIB = 85
        IF ( STRVIB .EQ. ' 000000101 e g ' ) IVIB = 86
        IF ( STRVIB .EQ. ' 000000101 e u ' ) IVIB = 87
        IF ( STRVIB .EQ. ' 000000101 f g ' ) IVIB = 88
        IF ( STRVIB .EQ. ' 000000101 f u ' ) IVIB = 89
        IF ( STRVIB .EQ. ' 000000102 e g ' ) IVIB = 90
        IF ( STRVIB .EQ. ' 000000102 e u ' ) IVIB = 91
        IF ( STRVIB .EQ. ' 000000102 f g ' ) IVIB = 92
        IF ( STRVIB .EQ. ' 000000102 f u ' ) IVIB = 93
        IF ( STRVIB .EQ. ' 000000103 e g ' ) IVIB = 94
        IF ( STRVIB .EQ. ' 000000103 e u ' ) IVIB = 95
        IF ( STRVIB .EQ. ' 000000103 f g ' ) IVIB = 96
        IF ( STRVIB .EQ. ' 000000103 f u ' ) IVIB = 97
        IF ( STRVIB .EQ. ' 000000105 e g ' ) IVIB = 98
        IF ( STRVIB .EQ. ' 000000105 e u ' ) IVIB = 99
        IF ( STRVIB .EQ. ' 000000105 f g ' ) IVIB = 100
        IF ( STRVIB .EQ. ' 000000105 f u ' ) IVIB = 101
        IF ( STRVIB .EQ. ' 000000106 e g ' ) IVIB = 102
        IF ( STRVIB .EQ. ' 000000106 e u ' ) IVIB = 103
        IF ( STRVIB .EQ. ' 000000106 f g ' ) IVIB = 104
        IF ( STRVIB .EQ. ' 000000106 f u ' ) IVIB = 105
        IF ( STRVIB .EQ. ' 000000110 e g ' ) IVIB = 106
        IF ( STRVIB .EQ. ' 000000110 e u ' ) IVIB = 107
        IF ( STRVIB .EQ. ' 000000110 f g ' ) IVIB = 108
        IF ( STRVIB .EQ. ' 000000110 f u ' ) IVIB = 109
        IF ( STRVIB .EQ. ' 000000111 e g ' ) IVIB = 110
        IF ( STRVIB .EQ. ' 000000111 e u ' ) IVIB = 111
        IF ( STRVIB .EQ. ' 000000111 f g ' ) IVIB = 112
        IF ( STRVIB .EQ. ' 000000111 f u ' ) IVIB = 113
        IF ( STRVIB .EQ. ' 000000112 e g ' ) IVIB = 114
        IF ( STRVIB .EQ. ' 000000112 e u ' ) IVIB = 115
        IF ( STRVIB .EQ. ' 000000112 f g ' ) IVIB = 116
        IF ( STRVIB .EQ. ' 000000112 f u ' ) IVIB = 117
        IF ( STRVIB .EQ. ' 000000113 e g ' ) IVIB = 118
        IF ( STRVIB .EQ. ' 000000113 e u ' ) IVIB = 119
        IF ( STRVIB .EQ. ' 000000113 f g ' ) IVIB = 120
        IF ( STRVIB .EQ. ' 000000113 f u ' ) IVIB = 121
        IF ( STRVIB .EQ. ' 000000120 e g ' ) IVIB = 122
        IF ( STRVIB .EQ. ' 000000120 e u ' ) IVIB = 123
        IF ( STRVIB .EQ. ' 000000120 f g ' ) IVIB = 124
        IF ( STRVIB .EQ. ' 000000120 f u ' ) IVIB = 125
        IF ( STRVIB .EQ. ' 000000200 e g ' ) IVIB = 126
        IF ( STRVIB .EQ. ' 000000200 e u ' ) IVIB = 127
        IF ( STRVIB .EQ. ' 000000200 f g ' ) IVIB = 128
        IF ( STRVIB .EQ. ' 000000200 f u ' ) IVIB = 129
        IF ( STRVIB .EQ. ' 000000201 e g ' ) IVIB = 130
        IF ( STRVIB .EQ. ' 000000201 e u ' ) IVIB = 131
        IF ( STRVIB .EQ. ' 000000201 f g ' ) IVIB = 132
        IF ( STRVIB .EQ. ' 000000201 f u ' ) IVIB = 133
        IF ( STRVIB .EQ. ' 000000202 e g ' ) IVIB = 134
        IF ( STRVIB .EQ. ' 000000202 e u ' ) IVIB = 135
        IF ( STRVIB .EQ. ' 000000202 f g ' ) IVIB = 136
        IF ( STRVIB .EQ. ' 000000202 f u ' ) IVIB = 137
        IF ( STRVIB .EQ. ' 000000210 e g ' ) IVIB = 138
        IF ( STRVIB .EQ. ' 000000210 e u ' ) IVIB = 139
        IF ( STRVIB .EQ. ' 000000210 f g ' ) IVIB = 140
        IF ( STRVIB .EQ. ' 000000210 f u ' ) IVIB = 141
        IF ( STRVIB .EQ. ' 000000211 e g ' ) IVIB = 142
        IF ( STRVIB .EQ. ' 000000211 e u ' ) IVIB = 143
        IF ( STRVIB .EQ. ' 000000211 f g ' ) IVIB = 144
        IF ( STRVIB .EQ. ' 000000211 f u ' ) IVIB = 145
        IF ( STRVIB .EQ. ' 000000212 e g ' ) IVIB = 146
        IF ( STRVIB .EQ. ' 000000212 e u ' ) IVIB = 147
        IF ( STRVIB .EQ. ' 000000212 f g ' ) IVIB = 148
        IF ( STRVIB .EQ. ' 000000212 f u ' ) IVIB = 149
        IF ( STRVIB .EQ. ' 000001000 e g ' ) IVIB = 150
        IF ( STRVIB .EQ. ' 000001000 e u ' ) IVIB = 151
        IF ( STRVIB .EQ. ' 000001000 f g ' ) IVIB = 152
        IF ( STRVIB .EQ. ' 000001000 f u ' ) IVIB = 153
        IF ( STRVIB .EQ. ' 000001001 e g ' ) IVIB = 154
        IF ( STRVIB .EQ. ' 000001001 e u ' ) IVIB = 155
        IF ( STRVIB .EQ. ' 000001001 f g ' ) IVIB = 156
        IF ( STRVIB .EQ. ' 000001001 f u ' ) IVIB = 157
        IF ( STRVIB .EQ. ' 000001002 e g ' ) IVIB = 158
        IF ( STRVIB .EQ. ' 000001002 e u ' ) IVIB = 159
        IF ( STRVIB .EQ. ' 000001002 f g ' ) IVIB = 160
        IF ( STRVIB .EQ. ' 000001002 f u ' ) IVIB = 161
        IF ( STRVIB .EQ. ' 000001003 e g ' ) IVIB = 162
        IF ( STRVIB .EQ. ' 000001003 e u ' ) IVIB = 163
        IF ( STRVIB .EQ. ' 000001003 f g ' ) IVIB = 164
        IF ( STRVIB .EQ. ' 000001003 f u ' ) IVIB = 165
        IF ( STRVIB .EQ. ' 000001005 f u ' ) IVIB = 166
        IF ( STRVIB .EQ. ' 000001006 e g ' ) IVIB = 167
        IF ( STRVIB .EQ. ' 000001006 e u ' ) IVIB = 168
        IF ( STRVIB .EQ. ' 000001006 f g ' ) IVIB = 169
        IF ( STRVIB .EQ. ' 000001010 e g ' ) IVIB = 170
        IF ( STRVIB .EQ. ' 000001010 e u ' ) IVIB = 171
        IF ( STRVIB .EQ. ' 000001010 f g ' ) IVIB = 172
        IF ( STRVIB .EQ. ' 000001010 f u ' ) IVIB = 173
        IF ( STRVIB .EQ. ' 000001011 e g ' ) IVIB = 174
        IF ( STRVIB .EQ. ' 000001011 e u ' ) IVIB = 175
        IF ( STRVIB .EQ. ' 000001011 f g ' ) IVIB = 176
        IF ( STRVIB .EQ. ' 000001011 f u ' ) IVIB = 177
        IF ( STRVIB .EQ. ' 000001012 e g ' ) IVIB = 178
        IF ( STRVIB .EQ. ' 000001012 e u ' ) IVIB = 179
        IF ( STRVIB .EQ. ' 000001012 f g ' ) IVIB = 180
        IF ( STRVIB .EQ. ' 000001012 f u ' ) IVIB = 181
        IF ( STRVIB .EQ. ' 000001013 e g ' ) IVIB = 182
        IF ( STRVIB .EQ. ' 000001013 e u ' ) IVIB = 183
        IF ( STRVIB .EQ. ' 000001013 f g ' ) IVIB = 184
        IF ( STRVIB .EQ. ' 000001013 f u ' ) IVIB = 185
        IF ( STRVIB .EQ. ' 000001020 e g ' ) IVIB = 186
        IF ( STRVIB .EQ. ' 000001020 e u ' ) IVIB = 187
        IF ( STRVIB .EQ. ' 000001020 f g ' ) IVIB = 188
        IF ( STRVIB .EQ. ' 000001020 f u ' ) IVIB = 189
        IF ( STRVIB .EQ. ' 000001100 e g ' ) IVIB = 190
        IF ( STRVIB .EQ. ' 000001100 e u ' ) IVIB = 191
        IF ( STRVIB .EQ. ' 000001100 f g ' ) IVIB = 192
        IF ( STRVIB .EQ. ' 000001100 f u ' ) IVIB = 193
        IF ( STRVIB .EQ. ' 000001110 e g ' ) IVIB = 194
        IF ( STRVIB .EQ. ' 000001110 e u ' ) IVIB = 195
        IF ( STRVIB .EQ. ' 000001110 f g ' ) IVIB = 196
        IF ( STRVIB .EQ. ' 000001110 f u ' ) IVIB = 197
        IF ( STRVIB .EQ. ' 000002000 e g ' ) IVIB = 198
        IF ( STRVIB .EQ. ' 000002000 e u ' ) IVIB = 199
        IF ( STRVIB .EQ. ' 000002000 f g ' ) IVIB = 200
        IF ( STRVIB .EQ. ' 000002000 f u ' ) IVIB = 201
        IF ( STRVIB .EQ. ' 000002001 e g ' ) IVIB = 202
        IF ( STRVIB .EQ. ' 000002001 e u ' ) IVIB = 203
        IF ( STRVIB .EQ. ' 000002001 f g ' ) IVIB = 204
        IF ( STRVIB .EQ. ' 000002001 f u ' ) IVIB = 205
        IF ( STRVIB .EQ. ' 000002002 e g ' ) IVIB = 206
        IF ( STRVIB .EQ. ' 000002002 e u ' ) IVIB = 207
        IF ( STRVIB .EQ. ' 000002002 f g ' ) IVIB = 208
        IF ( STRVIB .EQ. ' 000002002 f u ' ) IVIB = 209
        IF ( STRVIB .EQ. ' 000002003 e g ' ) IVIB = 210
        IF ( STRVIB .EQ. ' 000002003 e u ' ) IVIB = 211
        IF ( STRVIB .EQ. ' 000002003 f g ' ) IVIB = 212
        IF ( STRVIB .EQ. ' 000002003 f u ' ) IVIB = 213
        IF ( STRVIB .EQ. ' 000002010 e g ' ) IVIB = 214
        IF ( STRVIB .EQ. ' 000002010 e u ' ) IVIB = 215
        IF ( STRVIB .EQ. ' 000002010 f g ' ) IVIB = 216
        IF ( STRVIB .EQ. ' 000002010 f u ' ) IVIB = 217
        IF ( STRVIB .EQ. ' 000002100 e g ' ) IVIB = 218
        IF ( STRVIB .EQ. ' 000002100 e u ' ) IVIB = 219
        IF ( STRVIB .EQ. ' 000002100 f g ' ) IVIB = 220
        IF ( STRVIB .EQ. ' 000002100 f u ' ) IVIB = 221
        IF ( STRVIB .EQ. ' 000003000 e g ' ) IVIB = 222
        IF ( STRVIB .EQ. ' 000003000 e u ' ) IVIB = 223
        IF ( STRVIB .EQ. ' 000003000 f g ' ) IVIB = 224
        IF ( STRVIB .EQ. ' 000003000 f u ' ) IVIB = 225
        IF ( STRVIB .EQ. ' 001000000 e g ' ) IVIB = 226
        IF ( STRVIB .EQ. ' 001000000 e u ' ) IVIB = 227
        IF ( STRVIB .EQ. ' 001000001 e g ' ) IVIB = 228
        IF ( STRVIB .EQ. ' 001000001 e u ' ) IVIB = 229
        IF ( STRVIB .EQ. ' 001000001 f g ' ) IVIB = 230
        IF ( STRVIB .EQ. ' 001000001 f u ' ) IVIB = 231
        IF ( STRVIB .EQ. ' 001000002 e g ' ) IVIB = 232
        IF ( STRVIB .EQ. ' 001000002 e u ' ) IVIB = 233
        IF ( STRVIB .EQ. ' 001000002 f g ' ) IVIB = 234
        IF ( STRVIB .EQ. ' 001000002 f u ' ) IVIB = 235
        IF ( STRVIB .EQ. ' 001000003 e g ' ) IVIB = 236
        IF ( STRVIB .EQ. ' 001000003 e u ' ) IVIB = 237
        IF ( STRVIB .EQ. ' 001000003 f g ' ) IVIB = 238
        IF ( STRVIB .EQ. ' 001000003 f u ' ) IVIB = 239
        IF ( STRVIB .EQ. ' 001000004 e g ' ) IVIB = 240
        IF ( STRVIB .EQ. ' 001000004 e u ' ) IVIB = 241
        IF ( STRVIB .EQ. ' 001000004 f g ' ) IVIB = 242
        IF ( STRVIB .EQ. ' 001000004 f u ' ) IVIB = 243
        IF ( STRVIB .EQ. ' 001000010 e g ' ) IVIB = 244
        IF ( STRVIB .EQ. ' 001000010 e u ' ) IVIB = 245
        IF ( STRVIB .EQ. ' 001000010 f g ' ) IVIB = 246
        IF ( STRVIB .EQ. ' 001000010 f u ' ) IVIB = 247
        IF ( STRVIB .EQ. ' 001000011 e g ' ) IVIB = 248
        IF ( STRVIB .EQ. ' 001000011 e u ' ) IVIB = 249
        IF ( STRVIB .EQ. ' 001000011 f g ' ) IVIB = 250
        IF ( STRVIB .EQ. ' 001000011 f u ' ) IVIB = 251
        IF ( STRVIB .EQ. ' 001000012 e g ' ) IVIB = 252
        IF ( STRVIB .EQ. ' 001000012 e u ' ) IVIB = 253
        IF ( STRVIB .EQ. ' 001000012 f g ' ) IVIB = 254
        IF ( STRVIB .EQ. ' 001000012 f u ' ) IVIB = 255
        IF ( STRVIB .EQ. ' 001000100 e g ' ) IVIB = 256
        IF ( STRVIB .EQ. ' 001000100 e u ' ) IVIB = 257
        IF ( STRVIB .EQ. ' 001000100 f g ' ) IVIB = 258
        IF ( STRVIB .EQ. ' 001000100 f u ' ) IVIB = 259
        IF ( STRVIB .EQ. ' 001000101 e g ' ) IVIB = 260
        IF ( STRVIB .EQ. ' 001000101 e u ' ) IVIB = 261
        IF ( STRVIB .EQ. ' 001000101 f g ' ) IVIB = 262
        IF ( STRVIB .EQ. ' 001000101 f u ' ) IVIB = 263
        IF ( STRVIB .EQ. ' 001000102 e g ' ) IVIB = 264
        IF ( STRVIB .EQ. ' 001000102 e u ' ) IVIB = 265
        IF ( STRVIB .EQ. ' 001000102 f g ' ) IVIB = 266
        IF ( STRVIB .EQ. ' 001000102 f u ' ) IVIB = 267
        IF ( STRVIB .EQ. ' 001000200 e g ' ) IVIB = 268
        IF ( STRVIB .EQ. ' 001000200 e u ' ) IVIB = 269
        IF ( STRVIB .EQ. ' 001000200 f g ' ) IVIB = 270
        IF ( STRVIB .EQ. ' 001000200 f u ' ) IVIB = 271
        IF ( STRVIB .EQ. ' 001001000 e g ' ) IVIB = 272
        IF ( STRVIB .EQ. ' 001001000 e u ' ) IVIB = 273
        IF ( STRVIB .EQ. ' 001001000 f g ' ) IVIB = 274
        IF ( STRVIB .EQ. ' 001001000 f u ' ) IVIB = 275
        IF ( STRVIB .EQ. ' 001001001 e g ' ) IVIB = 276
        IF ( STRVIB .EQ. ' 001001001 e u ' ) IVIB = 277
        IF ( STRVIB .EQ. ' 001001001 f g ' ) IVIB = 278
        IF ( STRVIB .EQ. ' 001001001 f u ' ) IVIB = 279
        IF ( STRVIB .EQ. ' 001001002 e g ' ) IVIB = 280
        IF ( STRVIB .EQ. ' 001001002 e u ' ) IVIB = 281
        IF ( STRVIB .EQ. ' 001001002 f g ' ) IVIB = 282
        IF ( STRVIB .EQ. ' 001001002 f u ' ) IVIB = 283
        IF ( STRVIB .EQ. ' 002000000 e g ' ) IVIB = 284
        IF ( STRVIB .EQ. ' 002000000 e u ' ) IVIB = 285
        IF ( STRVIB .EQ. ' 002000001 e g ' ) IVIB = 286
        IF ( STRVIB .EQ. ' 002000001 e u ' ) IVIB = 287
        IF ( STRVIB .EQ. ' 002000001 f g ' ) IVIB = 288
        IF ( STRVIB .EQ. ' 002000001 f u ' ) IVIB = 289
      ELSE IF ( IDXMOL .EQ. 44 ) THEN                               ! HC3N
C new levels found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. '  0000000 0 0 0' ) IVIB = 1
        IF ( STRVIB .EQ. '  0000001 0 0 1' ) IVIB = 2
        IF ( STRVIB .EQ. '  0000002 0 0 0' ) IVIB = 3
        IF ( STRVIB .EQ. '  0000002 0 0 2' ) IVIB = 4
        IF ( STRVIB .EQ. '  0000003 0 0 1' ) IVIB = 5
        IF ( STRVIB .EQ. '  0000003 0 0 3' ) IVIB = 6
        IF ( STRVIB .EQ. '  0000004 0 0 0' ) IVIB = 7
        IF ( STRVIB .EQ. '  0000004 0 0 2' ) IVIB = 8
        IF ( STRVIB .EQ. '  0000004 0 0 4' ) IVIB = 9
        IF ( STRVIB .EQ. '  0000005 0 0 1' ) IVIB = 10
        IF ( STRVIB .EQ. '  0000005 0 0 3' ) IVIB = 11
        IF ( STRVIB .EQ. '  0000005 0 0 5' ) IVIB = 12
        IF ( STRVIB .EQ. '  0000006 0 0 0' ) IVIB = 13
        IF ( STRVIB .EQ. '  0000006 0 0 2' ) IVIB = 14
        IF ( STRVIB .EQ. '  0000006 0 0 4' ) IVIB = 15
        IF ( STRVIB .EQ. '  0000006 0 0 6' ) IVIB = 16
        IF ( STRVIB .EQ. '  0000007 0 0 1' ) IVIB = 17
        IF ( STRVIB .EQ. '  0000007 0 0 3' ) IVIB = 18
        IF ( STRVIB .EQ. '  0000007 0 0 5' ) IVIB = 19
        IF ( STRVIB .EQ. '  0000007 0 0 7' ) IVIB = 20
        IF ( STRVIB .EQ. '  0000008 0 0 0' ) IVIB = 21
        IF ( STRVIB .EQ. '  0000008 0 0 2' ) IVIB = 22
        IF ( STRVIB .EQ. '  0000008 0 0 4' ) IVIB = 23
        IF ( STRVIB .EQ. '  0000008 0 0 6' ) IVIB = 24
        IF ( STRVIB .EQ. '  0000009 0 0 1' ) IVIB = 25
        IF ( STRVIB .EQ. '  0000009 0 0 3' ) IVIB = 26
        IF ( STRVIB .EQ. '  0000009 0 0 5' ) IVIB = 27
        IF ( STRVIB .EQ. '  0000009 0 0 7' ) IVIB = 28
        IF ( STRVIB .EQ. '  0000010 0 1 0' ) IVIB = 29
        IF ( STRVIB .EQ. '  0000011 0 1 1' ) IVIB = 30
        IF ( STRVIB .EQ. '  0000011 0 1-1' ) IVIB = 31
        IF ( STRVIB .EQ. '  0000012 0 1 0' ) IVIB = 32
        IF ( STRVIB .EQ. '  0000012 0 1 2' ) IVIB = 33
        IF ( STRVIB .EQ. '  0000012 0-1 2' ) IVIB = 34
        IF ( STRVIB .EQ. '  0000013 0 1 1' ) IVIB = 35
        IF ( STRVIB .EQ. '  0000013 0 1 3' ) IVIB = 36
        IF ( STRVIB .EQ. '  0000013 0 1-1' ) IVIB = 37
        IF ( STRVIB .EQ. '  0000013 0-1 3' ) IVIB = 38
        IF ( STRVIB .EQ. '  0000014 0 1 0' ) IVIB = 39
        IF ( STRVIB .EQ. '  0000014 0 1 2' ) IVIB = 40
        IF ( STRVIB .EQ. '  0000014 0 1 4' ) IVIB = 41
        IF ( STRVIB .EQ. '  0000014 0-1 2' ) IVIB = 42
        IF ( STRVIB .EQ. '  0000014 0-1 4' ) IVIB = 43
        IF ( STRVIB .EQ. '  0000015 0 1 1' ) IVIB = 44
        IF ( STRVIB .EQ. '  0000015 0 1 3' ) IVIB = 45
        IF ( STRVIB .EQ. '  0000015 0 1 5' ) IVIB = 46
        IF ( STRVIB .EQ. '  0000015 0 1-1' ) IVIB = 47
        IF ( STRVIB .EQ. '  0000015 0-1 3' ) IVIB = 48
        IF ( STRVIB .EQ. '  0000015 0-1 5' ) IVIB = 49
        IF ( STRVIB .EQ. '  0000016 0 1 0' ) IVIB = 50
        IF ( STRVIB .EQ. '  0000016 0 1 2' ) IVIB = 51
        IF ( STRVIB .EQ. '  0000016 0 1 4' ) IVIB = 52
        IF ( STRVIB .EQ. '  0000016 0 1 6' ) IVIB = 53
        IF ( STRVIB .EQ. '  0000017 0 1 1' ) IVIB = 54
        IF ( STRVIB .EQ. '  0000017 0 1 3' ) IVIB = 55
        IF ( STRVIB .EQ. '  0000017 0 1 5' ) IVIB = 56
        IF ( STRVIB .EQ. '  0000017 0 1 7' ) IVIB = 57
        IF ( STRVIB .EQ. '  0000017 0 1-1' ) IVIB = 58
        IF ( STRVIB .EQ. '  0000020 0 0 0' ) IVIB = 59
        IF ( STRVIB .EQ. '  0000020 0 2 0' ) IVIB = 60
        IF ( STRVIB .EQ. '  0000021 0 0 1' ) IVIB = 61
        IF ( STRVIB .EQ. '  0000021 0 2 1' ) IVIB = 62
        IF ( STRVIB .EQ. '  0000021 0 2-1' ) IVIB = 63
        IF ( STRVIB .EQ. '  0000022 0 0 0' ) IVIB = 64
        IF ( STRVIB .EQ. '  0000022 0 0 2' ) IVIB = 65
        IF ( STRVIB .EQ. '  0000022 0 2 0' ) IVIB = 66
        IF ( STRVIB .EQ. '  0000022 0 2 2' ) IVIB = 67
        IF ( STRVIB .EQ. '  0000022 0 2-2' ) IVIB = 68
        IF ( STRVIB .EQ. '  0000023 0 0 1' ) IVIB = 69
        IF ( STRVIB .EQ. '  0000023 0 0 3' ) IVIB = 70
        IF ( STRVIB .EQ. '  0000023 0 2 1' ) IVIB = 71
        IF ( STRVIB .EQ. '  0000023 0 2 3' ) IVIB = 72
        IF ( STRVIB .EQ. '  0000023 0 2-1' ) IVIB = 73
        IF ( STRVIB .EQ. '  0000023 0-2 3' ) IVIB = 74
        IF ( STRVIB .EQ. '  0000024 0 0 0' ) IVIB = 75
        IF ( STRVIB .EQ. '  0000024 0 0 2' ) IVIB = 76
        IF ( STRVIB .EQ. '  0000024 0 0 4' ) IVIB = 77
        IF ( STRVIB .EQ. '  0000024 0 2 0' ) IVIB = 78
        IF ( STRVIB .EQ. '  0000024 0 2 2' ) IVIB = 79
        IF ( STRVIB .EQ. '  0000024 0 2 4' ) IVIB = 80
        IF ( STRVIB .EQ. '  0000024 0 2-2' ) IVIB = 81
        IF ( STRVIB .EQ. '  0000025 0 0 1' ) IVIB = 82
        IF ( STRVIB .EQ. '  0000025 0 0 3' ) IVIB = 83
        IF ( STRVIB .EQ. '  0000025 0 0 5' ) IVIB = 84
        IF ( STRVIB .EQ. '  0000025 0 2 1' ) IVIB = 85
        IF ( STRVIB .EQ. '  0000025 0 2 3' ) IVIB = 86
        IF ( STRVIB .EQ. '  0000025 0 2 5' ) IVIB = 87
        IF ( STRVIB .EQ. '  0000025 0 2-1' ) IVIB = 88
        IF ( STRVIB .EQ. '  0000030 0 1 0' ) IVIB = 89
        IF ( STRVIB .EQ. '  0000030 0 3 0' ) IVIB = 90
        IF ( STRVIB .EQ. '  0000031 0 1 1' ) IVIB = 91
        IF ( STRVIB .EQ. '  0000031 0 1-1' ) IVIB = 92
        IF ( STRVIB .EQ. '  0000031 0 3 1' ) IVIB = 93
        IF ( STRVIB .EQ. '  0000031 0 3-1' ) IVIB = 94
        IF ( STRVIB .EQ. '  0000032 0 1 0' ) IVIB = 95
        IF ( STRVIB .EQ. '  0000032 0 1 2' ) IVIB = 96
        IF ( STRVIB .EQ. '  0000032 0 3 0' ) IVIB = 97
        IF ( STRVIB .EQ. '  0000032 0 3 2' ) IVIB = 98
        IF ( STRVIB .EQ. '  0000032 0 3-2' ) IVIB = 99
        IF ( STRVIB .EQ. '  0000033 0 1 1' ) IVIB = 100
        IF ( STRVIB .EQ. '  0000033 0 1 3' ) IVIB = 101
        IF ( STRVIB .EQ. '  0000033 0 1-1' ) IVIB = 102
        IF ( STRVIB .EQ. '  0000033 0 3 1' ) IVIB = 103
        IF ( STRVIB .EQ. '  0000033 0 3 3' ) IVIB = 104
        IF ( STRVIB .EQ. '  0000033 0 3-1' ) IVIB = 105
        IF ( STRVIB .EQ. '  0000033 0 3-3' ) IVIB = 106
        IF ( STRVIB .EQ. '  0000040 0 0 0' ) IVIB = 107
        IF ( STRVIB .EQ. '  0000040 0 2 0' ) IVIB = 108
        IF ( STRVIB .EQ. '  0000040 0 4 0' ) IVIB = 109
        IF ( STRVIB .EQ. '  0000041 0 0 1' ) IVIB = 110
        IF ( STRVIB .EQ. '  0000041 0 2 1' ) IVIB = 111
        IF ( STRVIB .EQ. '  0000041 0 2-1' ) IVIB = 112
        IF ( STRVIB .EQ. '  0000041 0 4 1' ) IVIB = 113
        IF ( STRVIB .EQ. '  0000041 0 4-1' ) IVIB = 114
        IF ( STRVIB .EQ. '  0000100 1 0 0' ) IVIB = 115
        IF ( STRVIB .EQ. '  0000101 1 0 1' ) IVIB = 116
        IF ( STRVIB .EQ. '  0000101 1 0-1' ) IVIB = 117
        IF ( STRVIB .EQ. '  0000102 1 0 0' ) IVIB = 118
        IF ( STRVIB .EQ. '  0000102 1 0 2' ) IVIB = 119
        IF ( STRVIB .EQ. '  0000102-1 0 2' ) IVIB = 120
        IF ( STRVIB .EQ. '  0000103 1 0 1' ) IVIB = 121
        IF ( STRVIB .EQ. '  0000103 1 0 3' ) IVIB = 122
        IF ( STRVIB .EQ. '  0000103 1 0-1' ) IVIB = 123
        IF ( STRVIB .EQ. '  0000103-1 0 3' ) IVIB = 124
        IF ( STRVIB .EQ. '  0000104 1 0 0' ) IVIB = 125
        IF ( STRVIB .EQ. '  0000104 1 0 2' ) IVIB = 126
        IF ( STRVIB .EQ. '  0000104 1 0 4' ) IVIB = 127
        IF ( STRVIB .EQ. '  0000104-1 0 2' ) IVIB = 128
        IF ( STRVIB .EQ. '  0000104-1 0 4' ) IVIB = 129
        IF ( STRVIB .EQ. '  0000105 1 0 1' ) IVIB = 130
        IF ( STRVIB .EQ. '  0000105 1 0 3' ) IVIB = 131
        IF ( STRVIB .EQ. '  0000105 1 0 5' ) IVIB = 132
        IF ( STRVIB .EQ. '  0000105 1 0-1' ) IVIB = 133
        IF ( STRVIB .EQ. '  0000105-1 0 3' ) IVIB = 134
        IF ( STRVIB .EQ. '  0000105-1 0 5' ) IVIB = 135
        IF ( STRVIB .EQ. '  0000106 1 0 0' ) IVIB = 136
        IF ( STRVIB .EQ. '  0000106 1 0 2' ) IVIB = 137
        IF ( STRVIB .EQ. '  0000106 1 0 4' ) IVIB = 138
        IF ( STRVIB .EQ. '  0000106 1 0 6' ) IVIB = 139
        IF ( STRVIB .EQ. '  0000106-1 0 2' ) IVIB = 140
        IF ( STRVIB .EQ. '  0000106-1 0 4' ) IVIB = 141
        IF ( STRVIB .EQ. '  0000106-1 0 6' ) IVIB = 142
        IF ( STRVIB .EQ. '  0000110 1 1 0' ) IVIB = 143
        IF ( STRVIB .EQ. '  0000110 1-1 0' ) IVIB = 144
        IF ( STRVIB .EQ. '  0000111 1 1 1' ) IVIB = 145
        IF ( STRVIB .EQ. '  0000111 1 1-1' ) IVIB = 146
        IF ( STRVIB .EQ. '  0000111 1-1 1' ) IVIB = 147
        IF ( STRVIB .EQ. '  0000111-1 1 1' ) IVIB = 148
        IF ( STRVIB .EQ. '  0000112 1 1 0' ) IVIB = 149
        IF ( STRVIB .EQ. '  0000112 1 1 2' ) IVIB = 150
        IF ( STRVIB .EQ. '  0000112 1 1-2' ) IVIB = 151
        IF ( STRVIB .EQ. '  0000112 1-1 0' ) IVIB = 152
        IF ( STRVIB .EQ. '  0000112 1-1 2' ) IVIB = 153
        IF ( STRVIB .EQ. '  0000112-1 1 2' ) IVIB = 154
        IF ( STRVIB .EQ. '  0000113 1 1 1' ) IVIB = 155
        IF ( STRVIB .EQ. '  0000113 1 1 3' ) IVIB = 156
        IF ( STRVIB .EQ. '  0000113 1 1-1' ) IVIB = 157
        IF ( STRVIB .EQ. '  0000113-1 1 1' ) IVIB = 158
        IF ( STRVIB .EQ. '  0000113-1 1 3' ) IVIB = 159
        IF ( STRVIB .EQ. '  0000114 1 1 0' ) IVIB = 160
        IF ( STRVIB .EQ. '  0000114 1 1 2' ) IVIB = 161
        IF ( STRVIB .EQ. '  0000114 1 1 4' ) IVIB = 162
        IF ( STRVIB .EQ. '  0000114 1 1-2' ) IVIB = 163
        IF ( STRVIB .EQ. '  0000114-1 1 2' ) IVIB = 164
        IF ( STRVIB .EQ. '  0000114-1 1 4' ) IVIB = 165
        IF ( STRVIB .EQ. '  0000120 1 0 0' ) IVIB = 166
        IF ( STRVIB .EQ. '  0000120 1 2 0' ) IVIB = 167
        IF ( STRVIB .EQ. '  0000120-1 2 0' ) IVIB = 168
        IF ( STRVIB .EQ. '  0000121 1 0 1' ) IVIB = 169
        IF ( STRVIB .EQ. '  0000121 1 0-1' ) IVIB = 170
        IF ( STRVIB .EQ. '  0000121 1 2 1' ) IVIB = 171
        IF ( STRVIB .EQ. '  0000121 1 2-1' ) IVIB = 172
        IF ( STRVIB .EQ. '  0000121-1 2 1' ) IVIB = 173
        IF ( STRVIB .EQ. '  0000122 1 0 0' ) IVIB = 174
        IF ( STRVIB .EQ. '  0000122 1 0 2' ) IVIB = 175
        IF ( STRVIB .EQ. '  0000122 1 2 0' ) IVIB = 176
        IF ( STRVIB .EQ. '  0000122 1 2 2' ) IVIB = 177
        IF ( STRVIB .EQ. '  0000122 1 2-2' ) IVIB = 178
        IF ( STRVIB .EQ. '  0000122-1 0 2' ) IVIB = 179
        IF ( STRVIB .EQ. '  0000122-1 2 0' ) IVIB = 180
        IF ( STRVIB .EQ. '  0000122-1 2 2' ) IVIB = 181
        IF ( STRVIB .EQ. '  0000130 1 1 0' ) IVIB = 182
        IF ( STRVIB .EQ. '  0000130 1 3 0' ) IVIB = 183
        IF ( STRVIB .EQ. '  0000130-1 3 0' ) IVIB = 184
        IF ( STRVIB .EQ. '  0000200 0 0 0' ) IVIB = 185
        IF ( STRVIB .EQ. '  0000200 2 0 0' ) IVIB = 186
        IF ( STRVIB .EQ. '  0000201 0 0 1' ) IVIB = 187
        IF ( STRVIB .EQ. '  0000201 2 0 1' ) IVIB = 188
        IF ( STRVIB .EQ. '  0000201 2 0-1' ) IVIB = 189
        IF ( STRVIB .EQ. '  0000202 0 0 0' ) IVIB = 190
        IF ( STRVIB .EQ. '  0000202 0 0 2' ) IVIB = 191
        IF ( STRVIB .EQ. '  0000202 2 0 0' ) IVIB = 192
        IF ( STRVIB .EQ. '  0000202 2 0 2' ) IVIB = 193
        IF ( STRVIB .EQ. '  0000202 2 0-2' ) IVIB = 194
        IF ( STRVIB .EQ. '  0000203 0 0 1' ) IVIB = 195
        IF ( STRVIB .EQ. '  0000203 0 0 3' ) IVIB = 196
        IF ( STRVIB .EQ. '  0000203 2 0 1' ) IVIB = 197
        IF ( STRVIB .EQ. '  0000203 2 0 3' ) IVIB = 198
        IF ( STRVIB .EQ. '  0000203 2 0-1' ) IVIB = 199
        IF ( STRVIB .EQ. '  0000203-2 0 3' ) IVIB = 200
        IF ( STRVIB .EQ. '  0000210 0 1 0' ) IVIB = 201
        IF ( STRVIB .EQ. '  0000210 2 1 0' ) IVIB = 202
        IF ( STRVIB .EQ. '  0000211 0 1 1' ) IVIB = 203
        IF ( STRVIB .EQ. '  0000211 0 1-1' ) IVIB = 204
        IF ( STRVIB .EQ. '  0000211 2 1 1' ) IVIB = 205
        IF ( STRVIB .EQ. '  0000211 2 1-1' ) IVIB = 206
        IF ( STRVIB .EQ. '  0000300 1 0 0' ) IVIB = 207
        IF ( STRVIB .EQ. '  0000300 3 0 0' ) IVIB = 208
        IF ( STRVIB .EQ. '  0001000 0 0 0' ) IVIB = 209
        IF ( STRVIB .EQ. '  0001001 0 0 1' ) IVIB = 210
        IF ( STRVIB .EQ. '  0001002 0 0 0' ) IVIB = 211
        IF ( STRVIB .EQ. '  0001002 0 0 2' ) IVIB = 212
        IF ( STRVIB .EQ. '  0001003 0 0 1' ) IVIB = 213
        IF ( STRVIB .EQ. '  0001003 0 0 3' ) IVIB = 214
        IF ( STRVIB .EQ. '  0001004 0 0 0' ) IVIB = 215
        IF ( STRVIB .EQ. '  0001004 0 0 2' ) IVIB = 216
        IF ( STRVIB .EQ. '  0001004 0 0 4' ) IVIB = 217
        IF ( STRVIB .EQ. '  0001005 0 0 1' ) IVIB = 218
        IF ( STRVIB .EQ. '  0001005 0 0 3' ) IVIB = 219
        IF ( STRVIB .EQ. '  0001005 0 0 5' ) IVIB = 220
        IF ( STRVIB .EQ. '  0001010 0 1 0' ) IVIB = 221
        IF ( STRVIB .EQ. '  0001011 0 1 1' ) IVIB = 222
        IF ( STRVIB .EQ. '  0001011 0 1-1' ) IVIB = 223
        IF ( STRVIB .EQ. '  0001012 0 1 0' ) IVIB = 224
        IF ( STRVIB .EQ. '  0001012 0 1 2' ) IVIB = 225
        IF ( STRVIB .EQ. '  0001013 0 1 1' ) IVIB = 226
        IF ( STRVIB .EQ. '  0001013 0 1 3' ) IVIB = 227
        IF ( STRVIB .EQ. '  0001013 0 1-1' ) IVIB = 228
        IF ( STRVIB .EQ. '  0001020 0 0 0' ) IVIB = 229
        IF ( STRVIB .EQ. '  0001020 0 2 0' ) IVIB = 230
        IF ( STRVIB .EQ. '  0001021 0 0 1' ) IVIB = 231
        IF ( STRVIB .EQ. '  0001021 0 2 1' ) IVIB = 232
        IF ( STRVIB .EQ. '  0001021 0 2-1' ) IVIB = 233
        IF ( STRVIB .EQ. '  0001100 1 0 0' ) IVIB = 234
        IF ( STRVIB .EQ. '  0001101 1 0 1' ) IVIB = 235
        IF ( STRVIB .EQ. '  0001101 1 0-1' ) IVIB = 236
        IF ( STRVIB .EQ. '  0001102 1 0 0' ) IVIB = 237
        IF ( STRVIB .EQ. '  0001102 1 0 2' ) IVIB = 238
        IF ( STRVIB .EQ. '  0001102-1 0 2' ) IVIB = 239
        IF ( STRVIB .EQ. '  0001110 1 1 0' ) IVIB = 240
        IF ( STRVIB .EQ. '  0002000 0 0 0' ) IVIB = 241
        IF ( STRVIB .EQ. '  0002001 0 0 1' ) IVIB = 242
      ELSE IF ( IDXMOL .EQ. 47 ) THEN                               ! SO3
C new levels found in HITRAN2012 - arbitrarily assigned IVIB values
        IF ( STRVIB .EQ. ' 0 0 0 0 0 0A1''' ) IVIB = 1
        IF ( STRVIB .EQ. ' 0 0 0 0 1 1E'' ' ) IVIB = 2
        IF ( STRVIB .EQ. ' 0 0 0 0 2 0A1''' ) IVIB = 3
        IF ( STRVIB .EQ. ' 0 0 0 0 2 2E'' ' ) IVIB = 4
        IF ( STRVIB .EQ. ' 0 0 1 1 0 0E'' ' ) IVIB = 5
        IF ( STRVIB .EQ. ' 0 0 2 2 0 0E'' ' ) IVIB = 6
        IF ( STRVIB .EQ. ' 0 1 0 0 0 0A2"' ) IVIB = 7
        IF ( STRVIB .EQ. ' 0 1 0 0 1 1E" ' ) IVIB = 8
        IF ( STRVIB .EQ. ' 0 2 0 0 0 0A1''' ) IVIB = 9
        IF ( STRVIB .EQ. ' 1 0 0 0 0 0A1''' ) IVIB = 10
C
C Monatomic molecules
      ELSE IF ( IDXMOL .EQ. 34 ) THEN                               ! O
 	IF ( STRVIB .EQ. '              0' ) IVIB = 1
      ELSE
        WRITE ( *, * ) 'F-IVIB: Unrecognised value IDXMOL=', IDXMOL
        STOP
      END IF
C
C Check that level has been identified
      IF ( IVIB .EQ. 0 ) THEN
        IF ( NWRN(IDXMOL) .LT. MAXWRN ) THEN
          WRITE ( *, * ) 'W-IDXVIB: Unrecognised Vibrational Quanta=''',
     &      STRVIB, ''' for IDXMOL=', IDXMOL
          NWRN(IDXMOL) = NWRN(IDXMOL) + 1
          IF ( NWRN(IDXMOL) .EQ. MAXWRN ) WRITE ( *, * )
     &      'Further such warnings for this molecule suppressed'
        END IF
        IVIB = 999
      END IF
C
      IDXVIB = IVIB
C
      END
