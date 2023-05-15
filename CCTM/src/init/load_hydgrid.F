!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!


C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/work/rep/arc/CCTM/src/init/yamo/load_cgrid.F,v 1.9 2012/01/19 14:47:23 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE LOAD_DX1GRID ( FNAME, SPC_CAT, CGRID )

C-----------------------------------------------------------------------
C Function:
C   Initialize the model CGRID array from file data

C Revision history:
C   29 Nov 12 S.L.Napelenok: adapted from load_cgrid for cmaq 5.0
C   12 Dec 19 S.L.Napelenok: ddm-3d implementation for version 5.3.1
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID mechanism species
      USE UTILIO_DEFN
C     USE DDM3D_DEFN, ONLY: NP, NPMAX, SENNUM, SPCNAME, SEN_PAR, SENGRID
	  USE HDMod

      IMPLICIT NONE

      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:

      CHARACTER( 16 ) :: FNAME
      INTEGER      JDATE
      INTEGER      JTIME
      CHARACTER( 2 ) :: SPC_CAT
      TYPE(hyperdual), POINTER :: CGRID( :,:,:,: )  ! for initial CONC
C Parameters:

C minimum aerosol sulfate concentration [ ug/m**3 ]
c     REAL, PARAMETER :: AEROCONCMIN = 0.001

C The following two factors assume that sulfate density is 1.8e3 [ kg/m**3 ]
C and that the geometric mean diameter and geometric standard deviations
C for the Aitken mode are 0.01e-6 [ m ] and 1.7 respectively
C and are 0.07e-6 and 2.0 respectively for the accumulation mode.

C factor to calculate aerosol number concentration from aerosol sulfate mass
C concentration in the Aitken mode [ ug ].
c     REAL, PARAMETER :: NUMFACT_I = 2.988524 E11

C factor to calculate aerosol number concentration from aerosol sulfate mass
C concentration in the Accumulation mode [ ug ].
c     REAL, PARAMETER :: NUMFACT_J = 3.560191 E08

C fraction of sulfuric acid vapor taken as aerosol for first time step
c     REAL, PARAMETER :: SO4VAPTOAER = 0.999
C initial fraction of total aerosol sulfate in the Aitken mode
c     REAL, PARAMETER :: IFRACATKN = 0.04

      INTEGER, SAVE :: MXSPC
      INTEGER ASTAT

C File variables:

c     REAL      :: DENS( NCOLS,NROWS,NLAYS )       ! air density (kg/m^3)
      REAL      :: RHOJ( NCOLS,NROWS,NLAYS ) ! air density X Jacobian (kg/m^2)

C External Functions:

c     INTEGER, EXTERNAL :: FINDEX       !  looks up number in table.

C Local Variables

c     REAL         MWH2SO4                           ! H2SO4 molec. wt.
c     REAL         H2SO4CONV                         ! ppm -> ug/m**3
c     INTEGER      LSULF                             ! Gas chem CGRID index
c     INTEGER      ISO4AJ, ISO4AI, INUMATKN, INUMACC ! CGRID aerosol indices

      INTEGER      GXOFF, GYOFF               ! global origin offset from file

C for XTRACT3
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3
      INTEGER       :: STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI
      REAL      :: DBUFF( NCOLS,NROWS,NLAYS )

      INTEGER      SPC_STRT
      INTEGER      N_SPCS                     ! no. of species for this call
      INTEGER      NDX                        ! loop copy of INDX
c     INTEGER      ISUR                       ! surrogate index
      INTEGER, ALLOCATABLE, SAVE :: INDX( : ) ! Variable indices for all IC species
c     REAL,    ALLOCATABLE, SAVE :: ICBC_FAC( : ) ! Factor to be applied to ICs
      INTEGER      C, R, L, SPC, V            ! loop counters
      INTEGER      ASPC                       ! CGRID RHOJ pointer
c     INTEGER      NCOLSDENS, NROWSDENS       ! local for DENS

      CHARACTER( 16 ) :: PNAME = 'LOAD_DX1GRID'
      CHARACTER( 16 ) :: VNAME
c     CHARACTER( 16 ) :: CONCMIN
      CHARACTER( 96 ) :: XMSG = ' '
      CHARACTER( 40 ) :: CHWARN = 'Domain extents dfrnt from model for '
      CHARACTER( 24 ) :: ESTR1 = 'No IC found for species '
      CHARACTER( 34 ) :: ESTR2 = ' '
      CHARACTER( 34 ) :: ESTR3 = ' '

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

c     INTEGER, SAVE :: LOGDEV

c     INTEGER SENNUM
c     CHARACTER( 8 ) :: SPCNAME
c     CHARACTER( 16) :: SENNAME


C-----------------------------------------------------------------------
      JDATE = STDATE
      JTIME = STTIME
      
      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         CALL LOG_HEADING( LOGDEV, 'Loading Initial Hyd DX1 Sensitivities' ) 

         MXSPC = N_GC_SPC + N_AE_SPC + N_NR_SPC + N_TR_SPC
         ALLOCATE ( INDX( MXSPC ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'ERROR allocating INDX'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF
      END IF

	  ! This fname should be hyd dx1 information
      IF ( .NOT. OPEN3( FNAME, FSREAD3, PNAME ) ) THEN
         XMSG = 'Could not open ' // FNAME // ' file'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF
 
      IF ( .NOT. DESC3( FNAME ) ) THEN
         XMSG = 'Could not get ' // FNAME // ' file description'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF
 
      IF ( GL_NCOLS .NE. NCOLS3D ) THEN
         WRITE( LOGDEV,* ) ' '
         WRITE( LOGDEV,* ) '    WARNING: ' // CHWARN // FNAME
         WRITE( LOGDEV,* ) '>>  GL_NCOLS: ', GL_NCOLS, '  NCOLS3D: ', NCOLS3D
      END IF
 
      IF ( GL_NROWS .NE. NROWS3D ) THEN
         WRITE( LOGDEV,* ) ' '
         WRITE( LOGDEV,* ) '    WARNING: ' // CHWARN // FNAME
         WRITE( LOGDEV,* ) '>>  GL_NROWS: ', GL_NROWS, '  NROWS3D: ', NROWS3D
      END IF
 
      IF ( NLAYS .NE. NLAYS3D ) THEN
         XMSG = 'Wrong number of layers in ' // FNAME // ' file'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF

C     ESTR2 = ' in ' // TRIM( FNAME ) // '; Look for '
C     ESTR3 = ' in ' // TRIM( FNAME ) // '; set to ' // TRIM( CONCMIN )
      

C Get INDX
      DO SPC = 1, MXSPC
         INDX( SPC ) = 0
      END DO


      IF ( SPC_CAT .EQ. 'GC' ) THEN
         
         WRITE( XMSG,1009 ) 'transported gas-phase (reactive) species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         
         SPC_STRT = GC_STRT
         N_SPCS = N_GC_SPC
         
         ! Basically the code for reading a restart file in load_cgrid
         DO SPC = 1, N_SPCS
            NDX = INDEX1( GC_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( GC_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), GC_SPC( SPC )         

         END DO

      ELSE IF ( SPC_CAT .EQ. 'AE' ) THEN

         WRITE( XMSG,1009 ) 'transported aerosol species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         
         SPC_STRT = AE_STRT
         N_SPCS = N_AE_SPC
         DO SPC = 1, N_SPCS
            NDX = INDEX1( AE_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( AE_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), AE_SPC( SPC )

         END DO

      ELSE IF ( SPC_CAT .EQ. 'NR' ) THEN

         WRITE( XMSG,1009 ) 'transported non-reactive gas species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         SPC_STRT = NR_STRT
         N_SPCS = N_NR_SPC
         DO SPC = 1, N_SPCS

            NDX = INDEX1( NR_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( NR_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), NR_SPC( SPC )

         END DO

      ELSE IF ( SPC_CAT .EQ. 'TR' ) THEN

         WRITE( XMSG, 1009 ) 'transported inert tracer gas species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         SPC_STRT = TR_STRT
         N_SPCS = N_TR_SPC
         DO SPC = 1, N_SPCS

            NDX = INDEX1( TR_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( TR_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), TR_SPC( SPC )

         END DO

      ELSE IF ( SPC_CAT .EQ. 'RJ' ) THEN
         N_SPCS = 0
      ELSE
         XMSG = 'Species categories incorrect for CGRID '
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
      END IF
        
C Read into DX1GRID

      CALL SUBHFILE ( FNAME, GXOFF, GYOFF,
     &                STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI )
C IOFDESC common now loaded with FNAME header

      DO SPC = 1, N_SPCS
         V = SPC_STRT - 1 + SPC
         NDX = INDX( SPC )

            IF ( NDX .GT. 0 ) THEN
               IF ( .NOT. XTRACT3( FNAME, VNAME3D ( NDX ),
     &              1,NLAYS, STRTROWINI,ENDROWINI, STRTCOLINI,ENDCOLINI,
     &              JDATE, JTIME, DBUFF ) ) THEN
                  XMSG = 'Could not read ' // TRIM( VNAME3D( NDX ) )
     &                 // ' from ' // FNAME
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
               
               
               DO L = 1, NLAYS
                  DO R = 1, NROWS
                     DO C = 1, NCOLS
                        CGRID( C,R,L,V )%dx1 = REAL ( DBUFF( C,R,L ) , 8 )
                     END DO
                  END DO
               END DO
            ELSE
               DO L = 1, NLAYS
                  DO R = 1, NROWS
                     DO C = 1, NCOLS
                        CGRID( C,R,L,V )%dx1 = 0.0d0
                     END DO
                  END DO
               END DO
            END IF   ! INDX .GT. 0
      END DO

      IF ( N_SPCS .NE. 0 ) WRITE( LOGDEV,'(/ 5X, A)' )
     &                            SPC_CAT // ' loaded into DX1GRID'

      IF ( SPC_CAT .EQ. 'RJ' ) THEN
 
C Load RHOJ for transport and mixing ratio advection adjustment
 
C        VNAME = 'DENSA_J'
C        IF ( .NOT. XTRACT3( MET_CRO_3D, VNAME,
C    &              1,NLAYS, STRTROWMC3,ENDROWMC3, STRTCOLMC3,ENDCOLMC3,
C    &              JDATE, JTIME, RHOJ ) ) THEN
C            XMSG = 'Could not read DENSA_J from ' // MET_CRO_3D
C           CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
C        END IF
C
C        ASPC = GC_STRT - 1 + N_GC_SPCD
C        DO NP = 1, NPMAX
C           DO L = 1, NLAYS
C              DO R = 1, NROWS
C                 DO C = 1, NCOLS
C                    SENGRID( C,R,L,NP,ASPC ) = RHOJ( C,R,L )
C                 END DO
C              END DO
C           END DO
C        END DO
C
C        WRITE( LOGDEV,'(/ 5X, A)' ) 'Dens*Jacobian loaded into SENGRID'
	     DO L = 1, NLAYS
	     	DO R = 1, NROWS
	     		DO C = 1, NCOLS
	     			CGRID( C,R,L,V )%dx1 = 0.0d0
	     		ENDDO
	     	ENDDO
	     ENDDO
 
      END IF

C Close the file

!     IF ( .NOT. CLOSE3( FNAME ) ) THEN
!        XMSG = 'Could not close ' // FNAME // ' file'
!        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!     END IF

      RETURN

1009  FORMAT( 'Initial Condition Factors used for ', A )
1013  FORMAT( 5X, I3, 2X, A, 1PG13.5 )
      END
      
      
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE LOAD_DX2GRID ( FNAME, SPC_CAT, CGRID )

C-----------------------------------------------------------------------
C Function:
C   Initialize the model CGRID array from file data

C Revision history:
C   29 Nov 12 S.L.Napelenok: adapted from load_cgrid for cmaq 5.0
C   12 Dec 19 S.L.Napelenok: ddm-3d implementation for version 5.3.1
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID mechanism species
      USE UTILIO_DEFN
C     USE DDM3D_DEFN, ONLY: NP, NPMAX, SENNUM, SPCNAME, SEN_PAR, SENGRID
	  USE HDMod
      IMPLICIT NONE

      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:

      CHARACTER( 16 ) :: FNAME
      INTEGER      JDATE
      INTEGER      JTIME
      CHARACTER( 2 ) :: SPC_CAT
      TYPE(hyperdual), POINTER :: CGRID( :,:,:,: )  ! for initial CONC
C Parameters:

C minimum aerosol sulfate concentration [ ug/m**3 ]
c     REAL, PARAMETER :: AEROCONCMIN = 0.001

C The following two factors assume that sulfate density is 1.8e3 [ kg/m**3 ]
C and that the geometric mean diameter and geometric standard deviations
C for the Aitken mode are 0.01e-6 [ m ] and 1.7 respectively
C and are 0.07e-6 and 2.0 respectively for the accumulation mode.

C factor to calculate aerosol number concentration from aerosol sulfate mass
C concentration in the Aitken mode [ ug ].
c     REAL, PARAMETER :: NUMFACT_I = 2.988524 E11

C factor to calculate aerosol number concentration from aerosol sulfate mass
C concentration in the Accumulation mode [ ug ].
c     REAL, PARAMETER :: NUMFACT_J = 3.560191 E08

C fraction of sulfuric acid vapor taken as aerosol for first time step
c     REAL, PARAMETER :: SO4VAPTOAER = 0.999
C initial fraction of total aerosol sulfate in the Aitken mode
c     REAL, PARAMETER :: IFRACATKN = 0.04

      INTEGER, SAVE :: MXSPC
      INTEGER ASTAT

C File variables:

c     REAL      :: DENS( NCOLS,NROWS,NLAYS )       ! air density (kg/m^3)
      REAL      :: RHOJ( NCOLS,NROWS,NLAYS ) ! air density X Jacobian (kg/m^2)

C External Functions:

c     INTEGER, EXTERNAL :: FINDEX       !  looks up number in table.

C Local Variables

c     REAL         MWH2SO4                           ! H2SO4 molec. wt.
c     REAL         H2SO4CONV                         ! ppm -> ug/m**3
c     INTEGER      LSULF                             ! Gas chem CGRID index
c     INTEGER      ISO4AJ, ISO4AI, INUMATKN, INUMACC ! CGRID aerosol indices

      INTEGER      GXOFF, GYOFF               ! global origin offset from file

C for XTRACT3
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3
      INTEGER       :: STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI
      REAL      :: DBUFF( NCOLS,NROWS,NLAYS )

      INTEGER      SPC_STRT
      INTEGER      N_SPCS                     ! no. of species for this call
      INTEGER      NDX                        ! loop copy of INDX
c     INTEGER      ISUR                       ! surrogate index
      INTEGER, ALLOCATABLE, SAVE :: INDX( : ) ! Variable indices for all IC species
c     REAL,    ALLOCATABLE, SAVE :: ICBC_FAC( : ) ! Factor to be applied to ICs
      INTEGER      C, R, L, SPC, V            ! loop counters
      INTEGER      ASPC                       ! CGRID RHOJ pointer
c     INTEGER      NCOLSDENS, NROWSDENS       ! local for DENS

      CHARACTER( 16 ) :: PNAME = 'LOAD_DX2GRID'
      CHARACTER( 16 ) :: VNAME
c     CHARACTER( 16 ) :: CONCMIN
      CHARACTER( 96 ) :: XMSG = ' '
      CHARACTER( 40 ) :: CHWARN = 'Domain extents dfrnt from model for '
      CHARACTER( 24 ) :: ESTR1 = 'No IC found for species '
      CHARACTER( 34 ) :: ESTR2 = ' '
      CHARACTER( 34 ) :: ESTR3 = ' '
      
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

c     INTEGER, SAVE :: LOGDEV

c     INTEGER SENNUM
c     CHARACTER( 8 ) :: SPCNAME
c     CHARACTER( 16) :: SENNAME


C-----------------------------------------------------------------------
      JDATE = STDATE
      JTIME = STTIME
      
      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         CALL LOG_HEADING( LOGDEV, 'Loading Initial Hyd DX2 Sensitivities' ) 

         MXSPC = N_GC_SPC + N_AE_SPC + N_NR_SPC + N_TR_SPC
         ALLOCATE ( INDX( MXSPC ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'ERROR allocating INDX'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF
      END IF

	  ! This fname should be hyd sensitivities dx2 information
      IF ( .NOT. OPEN3( FNAME, FSREAD3, PNAME ) ) THEN
         XMSG = 'Could not open ' // FNAME // ' file'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF
 
      IF ( .NOT. DESC3( FNAME ) ) THEN
         XMSG = 'Could not get ' // FNAME // ' file description'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF
 
      IF ( GL_NCOLS .NE. NCOLS3D ) THEN
         WRITE( LOGDEV,* ) ' '
         WRITE( LOGDEV,* ) '    WARNING: ' // CHWARN // FNAME
         WRITE( LOGDEV,* ) '>>  GL_NCOLS: ', GL_NCOLS, '  NCOLS3D: ', NCOLS3D
      END IF
 
      IF ( GL_NROWS .NE. NROWS3D ) THEN
         WRITE( LOGDEV,* ) ' '
         WRITE( LOGDEV,* ) '    WARNING: ' // CHWARN // FNAME
         WRITE( LOGDEV,* ) '>>  GL_NROWS: ', GL_NROWS, '  NROWS3D: ', NROWS3D
      END IF
 
      IF ( NLAYS .NE. NLAYS3D ) THEN
         XMSG = 'Wrong number of layers in ' // FNAME // ' file'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF

C     ESTR2 = ' in ' // TRIM( FNAME ) // '; Look for '
C     ESTR3 = ' in ' // TRIM( FNAME ) // '; set to ' // TRIM( CONCMIN )
      

C Get INDX
      DO SPC = 1, MXSPC
         INDX( SPC ) = 0
      END DO


C      write(logdev,*) "load_dx2grid", SPC_CAT

      IF ( SPC_CAT .EQ. 'GC' ) THEN
         
         WRITE( XMSG,1009 ) 'transported gas-phase (reactive) species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         
         SPC_STRT = GC_STRT
         N_SPCS = N_GC_SPC
         
         ! Basically the code for reading a restart file in load_cgrid
         DO SPC = 1, N_SPCS
            NDX = INDEX1( GC_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( GC_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), GC_SPC( SPC )         

         END DO

      ELSE IF ( SPC_CAT .EQ. 'AE' ) THEN

         WRITE( XMSG,1009 ) 'transported aerosol species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         
         SPC_STRT = AE_STRT
         N_SPCS = N_AE_SPC
         DO SPC = 1, N_SPCS
            NDX = INDEX1( AE_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( AE_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), AE_SPC( SPC )

         END DO

      ELSE IF ( SPC_CAT .EQ. 'NR' ) THEN

         WRITE( XMSG,1009 ) 'transported non-reactive gas species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         SPC_STRT = NR_STRT
         N_SPCS = N_NR_SPC
         DO SPC = 1, N_SPCS

            NDX = INDEX1( NR_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( NR_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), NR_SPC( SPC )

         END DO

      ELSE IF ( SPC_CAT .EQ. 'TR' ) THEN

         WRITE( XMSG, 1009 ) 'transported inert tracer gas species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         SPC_STRT = TR_STRT
         N_SPCS = N_TR_SPC
         DO SPC = 1, N_SPCS

            NDX = INDEX1( TR_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( TR_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), TR_SPC( SPC )

         END DO

      ELSE IF ( SPC_CAT .EQ. 'RJ' ) THEN
         N_SPCS = 0
      ELSE
         XMSG = 'Species categories incorrect for CGRID '
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
      END IF
        
C Read into DX2GRID

      CALL SUBHFILE ( FNAME, GXOFF, GYOFF,
     &                STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI )
C IOFDESC common now loaded with FNAME header

      DO SPC = 1, N_SPCS
         V = SPC_STRT - 1 + SPC
         NDX = INDX( SPC )

            IF ( NDX .GT. 0 ) THEN
               IF ( .NOT. XTRACT3( FNAME, VNAME3D ( NDX ),
     &              1,NLAYS, STRTROWINI,ENDROWINI, STRTCOLINI,ENDCOLINI,
     &              JDATE, JTIME, DBUFF ) ) THEN
                  XMSG = 'Could not read ' // TRIM( VNAME3D( NDX ) )
     &                 // ' from ' // FNAME
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
               DO L = 1, NLAYS
                  DO R = 1, NROWS
                     DO C = 1, NCOLS
                        CGRID( C,R,L,V )%dx2 = REAL ( DBUFF( C,R,L ) , 8 ) 
                     END DO
                  END DO
               END DO
            ELSE
               DO L = 1, NLAYS
                  DO R = 1, NROWS
                     DO C = 1, NCOLS
                        CGRID( C,R,L,V )%dx2 = 0.0d0
                     END DO
                  END DO
               END DO
            END IF   ! INDX .GT. 0
      END DO

      IF ( N_SPCS .NE. 0 ) WRITE( LOGDEV,'(/ 5X, A)' )
     &                            SPC_CAT // ' loaded into DX2GRID'

      IF ( SPC_CAT .EQ. 'RJ' ) THEN
 
C Load RHOJ for transport and mixing ratio advection adjustment
 
C        VNAME = 'DENSA_J'
C        IF ( .NOT. XTRACT3( MET_CRO_3D, VNAME,
C    &              1,NLAYS, STRTROWMC3,ENDROWMC3, STRTCOLMC3,ENDCOLMC3,
C    &              JDATE, JTIME, RHOJ ) ) THEN
C            XMSG = 'Could not read DENSA_J from ' // MET_CRO_3D
C           CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
C        END IF
C
C        ASPC = GC_STRT - 1 + N_GC_SPCD
C        DO NP = 1, NPMAX
C           DO L = 1, NLAYS
C              DO R = 1, NROWS
C                 DO C = 1, NCOLS
C                    SENGRID( C,R,L,NP,ASPC ) = RHOJ( C,R,L )
C                 END DO
C              END DO
C           END DO
C        END DO
C
C        WRITE( LOGDEV,'(/ 5X, A)' ) 'Dens*Jacobian loaded into SENGRID'
	     DO L = 1, NLAYS
	     	DO R = 1, NROWS
	     		DO C = 1, NCOLS
	     			CGRID( C,R,L,V )%dx2 = 0.0d0
	     		ENDDO
	     	ENDDO
	     ENDDO
 
      END IF

C Close the file

!     IF ( .NOT. CLOSE3( FNAME ) ) THEN
!        XMSG = 'Could not close ' // FNAME // ' file'
!        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!     END IF

      RETURN

1009  FORMAT( 'Initial Condition Factors used for ', A )
1013  FORMAT( 5X, I3, 2X, A, 1PG13.5 )
      END
            
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE LOAD_DX1X2GRID ( FNAME, SPC_CAT, CGRID )

C-----------------------------------------------------------------------
C Function:
C   Initialize the model CGRID array from file data

C Revision history:
C   29 Nov 12 S.L.Napelenok: adapted from load_cgrid for cmaq 5.0
C   12 Dec 19 S.L.Napelenok: ddm-3d implementation for version 5.3.1
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID mechanism species
      USE UTILIO_DEFN
C     USE DDM3D_DEFN, ONLY: NP, NPMAX, SENNUM, SPCNAME, SEN_PAR, SENGRID
	  USE HDMod
      IMPLICIT NONE

      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:

      CHARACTER( 16 ) :: FNAME
      INTEGER      JDATE
      INTEGER      JTIME
      CHARACTER( 2 ) :: SPC_CAT
      TYPE(hyperdual), POINTER :: CGRID( :,:,:,: )  ! for initial CONC
C Parameters:

C minimum aerosol sulfate concentration [ ug/m**3 ]
c     REAL, PARAMETER :: AEROCONCMIN = 0.001

C The following two factors assume that sulfate density is 1.8e3 [ kg/m**3 ]
C and that the geometric mean diameter and geometric standard deviations
C for the Aitken mode are 0.01e-6 [ m ] and 1.7 respectively
C and are 0.07e-6 and 2.0 respectively for the accumulation mode.

C factor to calculate aerosol number concentration from aerosol sulfate mass
C concentration in the Aitken mode [ ug ].
c     REAL, PARAMETER :: NUMFACT_I = 2.988524 E11

C factor to calculate aerosol number concentration from aerosol sulfate mass
C concentration in the Accumulation mode [ ug ].
c     REAL, PARAMETER :: NUMFACT_J = 3.560191 E08

C fraction of sulfuric acid vapor taken as aerosol for first time step
c     REAL, PARAMETER :: SO4VAPTOAER = 0.999
C initial fraction of total aerosol sulfate in the Aitken mode
c     REAL, PARAMETER :: IFRACATKN = 0.04

      INTEGER, SAVE :: MXSPC
      INTEGER ASTAT

C File variables:

c     REAL      :: DENS( NCOLS,NROWS,NLAYS )       ! air density (kg/m^3)
      REAL      :: RHOJ( NCOLS,NROWS,NLAYS ) ! air density X Jacobian (kg/m^2)

C External Functions:

c     INTEGER, EXTERNAL :: FINDEX       !  looks up number in table.

C Local Variables

c     REAL         MWH2SO4                           ! H2SO4 molec. wt.
c     REAL         H2SO4CONV                         ! ppm -> ug/m**3
c     INTEGER      LSULF                             ! Gas chem CGRID index
c     INTEGER      ISO4AJ, ISO4AI, INUMATKN, INUMACC ! CGRID aerosol indices

      INTEGER      GXOFF, GYOFF               ! global origin offset from file

C for XTRACT3
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3
      INTEGER       :: STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI
      REAL      :: DBUFF( NCOLS,NROWS,NLAYS )

      INTEGER      SPC_STRT
      INTEGER      N_SPCS                     ! no. of species for this call
      INTEGER      NDX                        ! loop copy of INDX
c     INTEGER      ISUR                       ! surrogate index
      INTEGER, ALLOCATABLE, SAVE :: INDX( : ) ! Variable indices for all IC species
c     REAL,    ALLOCATABLE, SAVE :: ICBC_FAC( : ) ! Factor to be applied to ICs
      INTEGER      C, R, L, SPC, V            ! loop counters
      INTEGER      ASPC                       ! CGRID RHOJ pointer
c     INTEGER      NCOLSDENS, NROWSDENS       ! local for DENS

      CHARACTER( 16 ) :: PNAME = 'LOAD_DX1X2GRID'
      CHARACTER( 16 ) :: VNAME
c     CHARACTER( 16 ) :: CONCMIN
      CHARACTER( 96 ) :: XMSG = ' '
      CHARACTER( 40 ) :: CHWARN = 'Domain extents dfrnt from model for '
      CHARACTER( 24 ) :: ESTR1 = 'No IC found for species '
      CHARACTER( 34 ) :: ESTR2 = ' '
      CHARACTER( 34 ) :: ESTR3 = ' '

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

c     INTEGER, SAVE :: LOGDEV

c     INTEGER SENNUM
c     CHARACTER( 8 ) :: SPCNAME
c     CHARACTER( 16) :: SENNAME


C-----------------------------------------------------------------------
      
      JDATE = STDATE
      JTIME = STTIME
      
      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         CALL LOG_HEADING( LOGDEV, 'Loading Initial Hyd DX1X2 Sensitivities' ) 

         MXSPC = N_GC_SPC + N_AE_SPC + N_NR_SPC + N_TR_SPC
         ALLOCATE ( INDX( MXSPC ), STAT = ASTAT )
         IF ( ASTAT .NE. 0 ) THEN
            XMSG = 'ERROR allocating INDX'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF
      END IF

	  ! This fname should be hyd sensitivities dx1x2 information
      IF ( .NOT. OPEN3( FNAME, FSREAD3, PNAME ) ) THEN
         XMSG = 'Could not open ' // FNAME // ' file'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF
 
      IF ( .NOT. DESC3( FNAME ) ) THEN
         XMSG = 'Could not get ' // FNAME // ' file description'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF
 
      IF ( GL_NCOLS .NE. NCOLS3D ) THEN
         WRITE( LOGDEV,* ) ' '
         WRITE( LOGDEV,* ) '    WARNING: ' // CHWARN // FNAME
         WRITE( LOGDEV,* ) '>>  GL_NCOLS: ', GL_NCOLS, '  NCOLS3D: ', NCOLS3D
      END IF
 
      IF ( GL_NROWS .NE. NROWS3D ) THEN
         WRITE( LOGDEV,* ) ' '
         WRITE( LOGDEV,* ) '    WARNING: ' // CHWARN // FNAME
         WRITE( LOGDEV,* ) '>>  GL_NROWS: ', GL_NROWS, '  NROWS3D: ', NROWS3D
      END IF
 
      IF ( NLAYS .NE. NLAYS3D ) THEN
         XMSG = 'Wrong number of layers in ' // FNAME // ' file'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF

C     ESTR2 = ' in ' // TRIM( FNAME ) // '; Look for '
C     ESTR3 = ' in ' // TRIM( FNAME ) // '; set to ' // TRIM( CONCMIN )
      

C Get INDX
      DO SPC = 1, MXSPC
         INDX( SPC ) = 0
      END DO


C      write(logdev,*) "load_dx1x2grid", SPC_CAT

      IF ( SPC_CAT .EQ. 'GC' ) THEN
         
         WRITE( XMSG,1009 ) 'transported gas-phase (reactive) species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         
         SPC_STRT = GC_STRT
         N_SPCS = N_GC_SPC
         
         ! Basically the code for reading a restart file in load_cgrid
         DO SPC = 1, N_SPCS
            NDX = INDEX1( GC_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( GC_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), GC_SPC( SPC )         

         END DO

      ELSE IF ( SPC_CAT .EQ. 'AE' ) THEN

         WRITE( XMSG,1009 ) 'transported aerosol species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         
         SPC_STRT = AE_STRT
         N_SPCS = N_AE_SPC
         DO SPC = 1, N_SPCS
            NDX = INDEX1( AE_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( AE_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), AE_SPC( SPC )

         END DO

      ELSE IF ( SPC_CAT .EQ. 'NR' ) THEN

         WRITE( XMSG,1009 ) 'transported non-reactive gas species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         SPC_STRT = NR_STRT
         N_SPCS = N_NR_SPC
         DO SPC = 1, N_SPCS

            NDX = INDEX1( NR_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( NR_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), NR_SPC( SPC )

         END DO

      ELSE IF ( SPC_CAT .EQ. 'TR' ) THEN

         WRITE( XMSG, 1009 ) 'transported inert tracer gas species'
         WRITE( LOGDEV, * )
         CALL LOG_MESSAGE( LOGDEV, XMSG )
         SPC_STRT = TR_STRT
         N_SPCS = N_TR_SPC
         DO SPC = 1, N_SPCS

            NDX = INDEX1( TR_SPC( SPC ), NVARS3D, VNAME3D )
            IF ( NDX .NE. 0 ) THEN
            	INDX( SPC ) = NDX   ! index in the IC file
            ELSE
                XMSG = ESTR1 // TRIM( TR_SPC( SPC ) ) // ESTR3
                CALL M3MESG( XMSG )
            END IF
            
            IF ( INDX( SPC ) .GT. 0 )
     &         WRITE( LOGDEV,1013 ) INDX( SPC ), TR_SPC( SPC )

         END DO

      ELSE IF ( SPC_CAT .EQ. 'RJ' ) THEN
         N_SPCS = 0
      ELSE
         XMSG = 'Species categories incorrect for CGRID '
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
      END IF
        
C Read into DX1X2GRID

      CALL SUBHFILE ( FNAME, GXOFF, GYOFF,
     &                STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI )
C IOFDESC common now loaded with FNAME header

      DO SPC = 1, N_SPCS
         V = SPC_STRT - 1 + SPC
         NDX = INDX( SPC )

            IF ( NDX .GT. 0 ) THEN
               IF ( .NOT. XTRACT3( FNAME, VNAME3D ( NDX ),
     &              1,NLAYS, STRTROWINI,ENDROWINI, STRTCOLINI,ENDCOLINI,
     &              JDATE, JTIME, DBUFF ) ) THEN
                  XMSG = 'Could not read ' // TRIM( VNAME3D( NDX ) )
     &                 // ' from ' // FNAME
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
               DO L = 1, NLAYS
                  DO R = 1, NROWS
                     DO C = 1, NCOLS
                        CGRID( C,R,L,V )%dx1x2 = REAL ( DBUFF( C,R,L ) , 8 ) 
                     END DO
                  END DO
               END DO
            ELSE
               DO L = 1, NLAYS
                  DO R = 1, NROWS
                     DO C = 1, NCOLS
                        CGRID( C,R,L,V )%dx1x2 = 0.0d0
                     END DO
                  END DO
               END DO
            END IF   ! INDX .GT. 0
      END DO

      IF ( N_SPCS .NE. 0 ) WRITE( LOGDEV,'(/ 5X, A)' )
     &                            SPC_CAT // ' loaded into DX1X2GRID'

      IF ( SPC_CAT .EQ. 'RJ' ) THEN
 
C Load RHOJ for transport and mixing ratio advection adjustment
 
C        VNAME = 'DENSA_J'
C        IF ( .NOT. XTRACT3( MET_CRO_3D, VNAME,
C    &              1,NLAYS, STRTROWMC3,ENDROWMC3, STRTCOLMC3,ENDCOLMC3,
C    &              JDATE, JTIME, RHOJ ) ) THEN
C            XMSG = 'Could not read DENSA_J from ' // MET_CRO_3D
C           CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
C        END IF
C
C        ASPC = GC_STRT - 1 + N_GC_SPCD
C        DO NP = 1, NPMAX
C           DO L = 1, NLAYS
C              DO R = 1, NROWS
C                 DO C = 1, NCOLS
C                    SENGRID( C,R,L,NP,ASPC ) = RHOJ( C,R,L )
C                 END DO
C              END DO
C           END DO
C        END DO
C
C        WRITE( LOGDEV,'(/ 5X, A)' ) 'Dens*Jacobian loaded into SENGRID'
	     DO L = 1, NLAYS
	     	DO R = 1, NROWS
	     		DO C = 1, NCOLS
	     			CGRID( C,R,L,V )%dx1x2 = 0.0d0
	     		ENDDO
	     	ENDDO
	     ENDDO
 
      END IF

C Close the file

!     IF ( .NOT. CLOSE3( FNAME ) ) THEN
!        XMSG = 'Could not close ' // FNAME // ' file'
!        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!     END IF

      RETURN

1009  FORMAT( 'Initial Condition Factors used for ', A )
1013  FORMAT( 5X, I3, 2X, A, 1PG13.5 )
      END
           
