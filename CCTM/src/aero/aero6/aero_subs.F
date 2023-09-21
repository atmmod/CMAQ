
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

C routines for aerosol formation and transformation processes

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE AEROPROC( DT, COL, ROW, LAYER )

C-----------------------------------------------------------------------
C  SUBROUTINE AEROPROC advances the number, second moment, and mass
C   concentrations for each mode over the time interval DT.

C  KEY SUBROUTINES CALLED:
C     GETPAR, ORGAER5, GETCOAGS, INTERCOAG_GH, INTRACOAG_GH,
C     VOLINORG

C  KEY FUNCTIONS CALLED:  ERF, ERFC, GETAF

C  REVISION HISTORY:
C     Coded in December 1999 by Dr. Francis S. Binkowski
C      Modified from older versions used in CMAQ.

C FSB 05/17/00  new version of RPMARES included

C FSB 05/30/00  Fixed minor bug in awater.f

C FSB Correction to extinction coefficient in getbext.

C FSB 07/21/00 corrected units on ORGRATES and ORGBRATE_IN from
C      ppm/min to ppm/sec in AEROPROC and AEROSTEP.

C FSB 11/30/00 following changes from 07/28/00
C     Fixed problem for new emissions file version
C     Combined emissions for M3
C     Used a fixed value of Dpmin for plotting
C     Added variables OMEGA_AT & OMEGA_AC for partitioning
C     Eliminated the restriction on relative humidity for nucleation.
C     Added a branch in EQL for very low relative humidity (<1%).

C FSB Following changes to RPMARES:
C     Number of iterations reduced from 150 to 50.
C     Iterations are used only if 0.5 < = RATIO.
C     In calculating the molality of the bisulfate ion, a
C     MAX(1.0e-10, MSO4 ) is used.

C FSB 08/08/01 Changes to NEWPART to correct mass rate and to AEROSTEP
C     to correct sulfate, include emissions, and trap problem in
C     accumulation mode number calculation.

C FSB 09/19/01 Version AE3, major changes
C     Organics are done with Dr. Benedikt Schell's approach. (ORGAER3)
C     Particle production uses Kulmala approach (NEWPART3)
C     Condensational factors for sulfate and organics are calculated
C      separately.
C     Emissions are assumed to be input in vertical diffusion.
C     Include files replaced by Fortran 90 Modules.
C     The modules also contain subroutines.

C FSB 10/24/01 Added the treatment of gaseous N2O5 -> aerosol HNO3

C FSB 10/25/01 Changed the mass transfer calculation for Aitken to
C     accumulation mode by coagulation as recommended by Dr. Benedikt
C     Schell.

C SJR 04/24/02 Replaced thermodynamic code RPMARES with ISORROPIA
C     --Shawn J. Roselle

C GLG 04/04/03 Modifications to allow for evaporation of semi-volatile
C     organics from aerosol phase.  --Gerald L. Gipson

C FSB 11/18/03 Corrections to sulfate condensation.  Previously,
C     SCONDRATE was undefined when SO4RATE=0 and SCONDRATE was
C     negative when SO4RATE < DMDT_so4.

C PVB 01/08/04 Several changes in preparation for simulation of 2001
C      calendar year --Dr. Prakash V. Bhave
C    -Added interface to new subroutine, GETCOAGS, for calculating
C      coagulation rates.  GETCOAGS is used instead of Gauss-Hermite
C      quadrature for computational efficiency, by setting
C      FASTCOAG_FLAG = .TRUE.  --PVB
C    -Removed SOA from the definition of "DRY" aerosol.  Aerosol surface
C      area is now transported without SOA.  See notes in GETPAR, AERO,
C      and AERO_DEPV subroutines.  --SJR
C    -Moved EQL3 call from AEROPROC to AEROSTEP, immediately following
C      the ORGAER3 call.  This is a side effect of transporting aerosol
C      surface area without SOA.  --SJR
C    -New subroutine, HCOND3, to calculate condensation rates for
C      2nd and 3rd moments.  Results are unchanged.  --FSB
C    -Revised method of calculating SOA.  Partition SOA to the modes
C      in proportion to the amounts of total organic mass (SOA plus
C      primary) in each mode.  Modal geometric standard deviations are
C      now preserved during SOA condensation and evaporation. --FSB
C    -Combined the former subroutines AEROPROC and AEROSTEP into this
C      subroutine; retained the name AEROPROC.  --PVB

C PVB 09/21/04 added in-line documentation with input from FSB.
C     Changed MWH2SO4 from 98.07354 to 98.0 g/mol for consistency with
C     the mechanism files.

C PVB 09/27/04 removed the IF(XM3.GT.0.0) mode-merging precondition
C     because it caused significant erroneous mode crossover.  Fix
C     suggested by Dr. Chris Nolte.

C PVB 01/19/05 Added SO4RATE to the EQL3 call vector.  This is necessary
C     to ensure that the gas and inorganic fine PM (i+j) concentrations
C     are in thermodynamic equilibrium at the end of each time step.

C PVB 05/02/05 Modified ERF statement function for negative arguments such
C     that erf(-x) = -erf(x).  Previous version had erf(-x) = erf(x).

C PVB 04/06/06 Added GAMMA_N2O5 to the AEROPROC and EQL3 call vectors,
C     so it can be written to the aerosol diagnostic file.

C SLN 09/07/07 Several changes in preparation for new SOA module
C    -Removed Aitken-mode SOA species (CBLK(VORGAI),CBLK(VORGBI)) and
C     all of the related variables (OMASS_I, OMASS_J, FRACI, FRACJ,
C     OLD_M2_I, OLD_M3_I, NEW_M2_I, NEW_M3_I)
C    -Replaced the other SOA species (CBLK(VORGAJ),CBLK(VORGBJ)) with
C     precursor-specific SOA species (e.g., CBLK(VTOL1J)).
C    -Replaced SOA_A and SOA_B with an array, SOA_ALL.  Replaced
C     OLDSOA_A and OLDSOA_B with an array, OLDSOA.  Updated ORGAER3
C     call vector accordingly.
C    -Removed dependence of mode-merging criteria on SOA condensation
C     because that logic was flawed.  Deleted all variables needed
C     for that calculation (CHEMRATE_ORG, ORGRATE, ORGBRATE).

C PVB 11/02/07 Moved heterogenous N2O5 chemistry from EQL3 to a new
C     subroutine, HETCHEM.

C PVB 11/29/07 Implementation of new SOA formation mechanisms
C    -Replaced call to ORGAER3 with a call to ORGAER5
C    -Moved CBLK updates of SV species, SOA species, M2, and M3 into
C     ORGAER5 instead of updating those values in AEROPROC

C JOY 04/08/08 White space, alignment, readability

C JTK 04/17/08 Implemented coarse chemistry updates
C     -replaced condensation calculations with call to VOLINORG
C     -added code for coarse surface area and variable coarse std. dev.
C     -modified modal dynamics equations because new particle formation
C      and growth are now calculated in VOLINORG

C JOY 12/07/08 Removed unused EQL3 code, add column and row arguments to
C     VOLINORG for diagnostic logging

C SH  12/17/09 Major restructuring of all aerosol codes
C     -use new Fortran modules (e.g., AERO_DATA, SOA_DEFN); shortened call
C      vector; replaced CBLK with aerospc structure; eliminated AERO_INFO

C SH  03/10/11 Renamed met_data to aeromet_data

C HS  03/10/11 Added call to POAAGE between ORGAER and HETCHEM

C SR  03/25/11 Replaced I/O API include files with UTILIO_DEFN

C BH  09/30/13 Removed Call to HETCHEM and GAMMA_N2O5 argument
C              because of merging of gas and heterogeneous 
C              chemistry

C JB  02/07/14 Added Jai Xing's mass balance fix when there is excessive
C              condensation or evaporation under cold conditions

C BH  07/21/14 Removed call to POAAGE because reaction represented in
C              chemical mechanism

C GS/KF 09/23/14 Updated nucleation scheme with Vehkamaki et al. (2002)

C HP/BM  4/16 Removed m3_wet_flag since wet/dry status is now
C             tracked by AERO_DATA variable. Moments are always wet in AEROPROC
C             Updated treatment of aerosol moments

C SLN 12/30/19 Update for ddm-3d implementation in version 5.3.1

C  REFERENCES:
C   1. Binkowski, F.S. and U. Shankar, The regional particulate matter
C      model 1. Model description and preliminary results, J. Geophys.
C      Res., Vol 100, No D12, 26101-26209, 1995.

C   2. Binkowski, F.S. and S.J. Roselle, Models-3 Community
C      Multiscale Air Quality (CMAQ) model aerosol component 1:
C      Model Description.  J. Geophys. Res., Vol 108, No D6, 4183
C      doi:10.1029/2001JD001409, 2003.

C   3. Bhave, P.V., S.J. Roselle, F.S. Binkowski, C.G. Nolte, S. Yu,
C      G.L. Gipson, and K.L. Schere, CMAQ aerosol module development:
C      recent enhancements and future plans, Paper No. 6.8, CMAS Annual
C      Conference, Chapel Hill, NC, 2004.
C-----------------------------------------------------------------------
      
      USE AERO_DATA
      USE PRECURSOR_DATA      ! gas phase aero precursor data
      USE SOA_DEFN
      USE AEROMET_DATA
      USE UTILIO_DEFN
      USE AERO_BUDGET

#ifdef sens
      USE DDM3D_DEFN, ONLY : NP, NPMAX
      USE AERO_DDM3D, ONLY : s_aerospc_conc, s_precursor_conc, s_so4rate, s_vapor_conc, s_orgprod_conc
#endif


	  use HDMod 
      IMPLICIT NONE

C *** arguments:
      REAL,    INTENT( IN ) :: DT          ! synchronization time step, sec
      INTEGER, INTENT( IN ) :: COL         ! Column of cell
      INTEGER, INTENT( IN ) :: ROW         ! Row of cell
      INTEGER, INTENT( IN ) :: LAYER       ! Layer of cell

C *** Parameters
      REAL, PARAMETER :: DIFFSULF = 9.36E-06  ! molecular diffusiviity for sulfuric acid
      REAL, PARAMETER :: SQRT2 = 1.4142135623731
      REAL, PARAMETER :: T0 = 288.15   ! [ K ]
      REAL, PARAMETER :: TWOTHIRDS   =  2.0 / 3.0
      REAL, PARAMETER :: ONE_OVER_TICE =  1.0 / 273.16

C *** local variables
      REAL         DIFFCORR   ! Correction to DIFFSULF & DIFFORG for pressure
      REAL         DV_SO4     ! molecular diffusivity of H2SO4 vapor after correction for ambient conditions
      REAL         SQRT_TEMP  ! square root of ambient temperature
      REAL         XLM        ! atmospheric mean free path [m]
      REAL         AMU        ! atmospheric dynamic viscosity [kg/m s]

      REAL( 8 ) :: CGR( N_MODE-1 ) ! Aitken & Accum. modes

      REAL( 8 ) :: LAMDA      ! mean free path [ m ]
      REAL( 8 ) :: KNC        ! KNC = TWO3 * BOLTZMANN *  AIRTEMP / AMU

C *** Free Molecular regime (depends upon modal density)
      TYPE(hyperdual) :: KFMAT      ! = SQRT( 3.0*BOLTZMANN * AIRTEMP / PDENSAT )
      TYPE(hyperdual) :: KFMAC      ! = SQRT( 3.0*BOLTZMANN * AIRTEMP / PDENSAC )
      TYPE(hyperdual) :: KFMATAC    ! = SQRT( 6.0*BOLTZMANN * AIRTEMP / ( PDENSAT + PDENSAC ) )

C *** Intermodal coagulation rates [ m**3/s ] ( 0th & 2nd moments )
      TYPE(hyperdual) :: BATAC( 2 ) ! Aitken to accumulation
      TYPE(hyperdual) :: BACAT( 2 ) ! accumulation from Aitken

C *** Intramodal coagulation rates [ m**3/s ] ( 0th & 2nd moments )
      TYPE(hyperdual) :: BATAT( 2 ) ! Aitken mode
      TYPE(hyperdual) :: BACAC( 2 ) ! accumulation mode

C *** Intermodal coagulation rate [ m**3/s ] ( 3rd moment )
      TYPE(hyperdual) :: C3IJ       ! Aitken to accumulation
      TYPE(hyperdual) :: C30ATAC    ! Aitken to accumulation
      TYPE(hyperdual) :: DG_D ( N_MODE )
      TYPE(hyperdual) :: SG_D ( N_MODE )
      TYPE(hyperdual) :: XXL_D( N_MODE )

! *** Variables for computing the budget
      TYPE(hyperdual) :: CBUDGET0_NUM ( N_MODE )
      TYPE(hyperdual) :: CBUDGET0_SRF ( N_MODE )
      TYPE(hyperdual) :: CBUDGET0_MASS( N_AEROSPC,N_MODE )

C *** variables for advancing concentrations one time step
      TYPE(hyperdual) :: A, B
      TYPE(hyperdual) :: Y0, Y
      TYPE(hyperdual) :: EXPDT
      TYPE(hyperdual) :: LOSS, PROD, POL
      TYPE(hyperdual) :: TMASS
      TYPE(hyperdual) :: FACTRANS ! special factor to compute mass transfer
      TYPE(hyperdual) :: M20      ! for initial condidtions in time stepping

C *** Variables for mode merging
      TYPE(hyperdual) ::  GETAF
      TYPE(hyperdual) ::  AAA, XNUM, XM2, XM3, XXM2, XXM3
      TYPE(hyperdual) ::  FNUM, FM2, FM3, PHNUM, PHM2, PHM3

C *** Comment this out, wrote overloading function for hyperdual
C *** ERF and ERFC in util/util/HDMod.f90
c  *** REAL       ERF, ERFC 


C *** local variables
      INTEGER      SPC        ! loop counter
      INTEGER      N          ! loop counter
      INTEGER      MODE       ! loop counter (DO NOT REMOVE)


      logical, save :: firstime = .true.
C-----------------------------------------------------------------------

      if ( firstime ) then
         firstime = .false.
#ifdef nomm
         write( logdev,* ) 'aero_subs: nomm'
#endif
      end if

C *** square root of the ambient temperature for later use
      SQRT_TEMP = SQRT( AIRTEMP )
 
C *** Calculate mean free path [ m ]:
C     6.6328E-8 is the sea level value given in Table I.2.8
C     on page 10 of U.S. Standard Atmosphere 1962

      XLM = 6.6328E-8 * STDATMPA * AIRTEMP  / ( T0 * AIRPRES )

C *** Calculate dynamic viscosity [ kg m**-1 s**-1 ]:
C     U.S. Standard Atmosphere 1962 page 14 expression
C     for dynamic viscosity is:
C     dynamic viscosity =  beta * T * sqrt(T) / ( T + S)
C     where beta = 1.458e-6 [ kg sec^-1 K**-0.5 ], s = 110.4 [ K ].
      AMU = 1.458E-6 * AIRTEMP * SQRT_TEMP / ( AIRTEMP + 110.4 )

C *** Set minimums for coarse mode
      MOMENT0_CONC( N_MODE ) = MAX( AEROMODE_MINNUM( N_MODE ),
     &                              MOMENT0_CONC( N_MODE ) )

C *** Secondary Organics
C     Update the secondary organic aerosol (SOA) mass concentrations
C     and the SVOC mass concentrations by equilibrium absorptive
C     partitioning between the particle and vapor phases.  Assume all
C     SOA resides in the accumulation mode.
C
C     Aerosol is wet when it enters orgaer
      CALL ORGAER( DT, LAYER )

C *** Secondary Inorganics
C     The VOLINORG subroutine includes the treatment of new particle
C     production and a fully dynamic treatment of inorganic gas-to-
C     particle mass transfer.

C *** Compute H2SO4 diffusivity, correct for temperature and pressure
      DIFFCORR = ( STDATMPA / AIRPRES ) * ( ONE_OVER_TICE * AIRTEMP ) ** 1.75
      DV_SO4 = DIFFSULF * DIFFCORR

C *** Update size parameters (distribution is wet from ORGAER)
      CALL GETPAR( FIXED_sg )

! *** Mass transfer of inorganic constituents      
! *** Process Budget is calculated inside this routine because NPF and
! *** condensation are calculated simultaneously
      CALL VOLINORG( DT, COL, ROW, LAYER, DV_SO4, CGR )

! *** Coagulation
!     Calculate coagulation coefficients using a method dictated by
!     the value of FASTCOAG_FLAG.  If TRUE, the computationally-
!     efficient GETCOAGS routine is used.  If FALSE, the more intensive
!     Gauss-Hermite numerical quadrature method is used.  See Section
!     2.1 of Bhave et al. (2004) for further discussion.
! *** Initialize Budget Variables
      CBUDGET0_NUM  = MOMENT0_CONC
      CBUDGET0_SRF  = MOMENT2_CONC
      CBUDGET0_MASS = AEROSPC_CONC
 
! *** set atmospheric mean free path in double precision
      LAMDA    = XLM

! *** calculate term used in Equation A6 of Binkowski & Shankar (1995)
      KNC      = TWOTHIRDS * BOLTZMANN *  AIRTEMP / AMU

! *** calculate terms used in Equation A5 of Binkowski & Shankar (1995)
      KFMAT    = SQRT( 3.0 * BOLTZMANN * AIRTEMP / AEROMODE_DENS( 1 ) )
      KFMAC    = SQRT( 3.0 * BOLTZMANN * AIRTEMP / AEROMODE_DENS( 2 ) )
      KFMATAC  = SQRT( 6.0 * BOLTZMANN * AIRTEMP
     &         / ( AEROMODE_DENS( 1 ) + AEROMODE_DENS( 2 ) ) )

! *** transfer of number to accumulation mode from Aitken mode is zero
       BACAT( 1 ) = 0.0d0

       IF ( FASTCOAG_FLAG ) THEN ! Solve coagulation analytically
 
C *** set geometric mean diameters, geometric standard deviations, and
C     ln(GSD) in double precision
         DO N = 1, N_MODE
            DG_D( N ) = AEROMODE_DIAM( N )
            SG_D( N ) = EXP( AEROMODE_LNSG( N ) )
            XXL_D( N ) = AEROMODE_LNSG( N )
         END DO

C *** calculate intermodal and intramodal coagulation coefficients
C     for zeroth and second moments, and intermodal coagulation
C     coefficient for third moment
         CALL GETCOAGS( LAMDA, KFMATAC, KFMAT, KFMAC, KNC,
     &                  DG_D(1), DG_D(2), SG_D(1), SG_D(2),
     &                  XXL_D(1),XXL_D(2),
     &                  BATAT( 2 ), BATAT( 1 ), BACAC( 2 ), BACAC( 1 ),
     &                  BATAC( 2 ), BACAT( 2 ), BATAC( 1 ), C3IJ )

      ELSE                 ! Use Gauss-Hermite numerical quadrature

C *** calculate Aitken-mode intramodal coagulation coefficients
C     for zeroth and second moments
         CALL INTRACOAG_GH( LAMDA, KFMAT, KNC, AEROMODE_DIAM( 1 ),
     &                      AEROMODE_LNSG( 1 ), BATAT( 2 ), BATAT( 1 ) )

C *** calculate accumulation-mode intramodal coagulation coefficients
C     for zeroth and second moments
         CALL INTRACOAG_GH( LAMDA, KFMAC, KNC, AEROMODE_DIAM( 2 ),
     &                      AEROMODE_LNSG( 2 ), BACAC( 2 ), BACAC( 1 ) )

C *** calculate intermodal coagulation coefficients for zeroth, second,
C     and third moments
         CALL INTERCOAG_GH( LAMDA, KFMATAC, KNC,
     &                      AEROMODE_DIAM( 1 ), AEROMODE_DIAM( 2 ),
     &                      AEROMODE_LNSG( 1 ), AEROMODE_LNSG( 2 ),
     &                      BATAC( 2 ), BACAT( 2 ), BATAC( 1 ), C3IJ )

      END IF   ! FASTCOAG_FLAG

C *** calculate 3rd moment intermodal transfer rate by coagulation
      C30ATAC = C3IJ * MOMENT0_CONC( 1 ) * MOMENT0_CONC( 2 )

C *** TAKE ONE FORWARD TIME STEP - Solve Modal Dynamics Equations
C     This code implements Section 1.4 of Binkowski and Roselle (2003)
C     with two notable exceptions.  1) emissions are treated in
C     CMAQ`s vertical diffusion routine, so they do not appear in the
C     following equations. 2) new particle formation and condensational
C     growth are now treated in the VOLINORG subroutine, so they do not
C     appear in the following equations.
 
C     M2 is updated before M0 because the intermodal transfer rate of
C     M2 is a function of the number concentrations.  In contrast,
C     production and loss rates of M0 are independent of M2.  Advancing
C     M2 before M0 avoids operator splitting within the modal-dynamic-
C     equation solution.  A similar rearrangement would be necessary
C     for the M3 update, but the dependence of M3 on number
C     concentrations already is accounted for in the C30ATAC term.

C *** UPDATE SECOND MOMENT
C     For each lognormal mode, solve equations of form:
C        dM2/dt = P2 - L2*M2   ! if L2 > 0
C     with solution
C        M2(t) = P2/L2 + ( M2(t0) - P2/L2 ) * exp( -L2*dt )
C     or
C        dM2/dt = P2           ! if L2 = 0
C     with solution
C        M2(t) = M2(t0) + P2*dt

C *** Aitken mode: initial value of M2
      M20 = MOMENT2_CONC( 1 )

C *** Loss of 2nd moment from Aitken mode is due to intermodal
C     coagulation with accumulation mode and intramodal coagulation.
C     Production term is removed, because new particle formation
C     and condensational growth are accounted for in VOLINORG.
      LOSS = (
     &        ( BATAT( 2 ) * MOMENT0_CONC( 1 )
     &        + BATAC( 2 ) * MOMENT0_CONC( 2 ) ) * MOMENT0_CONC( 1 )
     &       ) / M20

C *** Solve for M2_Aitken based on LOSS during this time step
C     Note: LOSS is assumed to be non-negative.
      IF ( LOSS .GT. 0.0 ) THEN
         Y = M20 * EXP( -LOSS * DT )
      ELSE
         Y = M20
      END IF ! test on loss

C *** Transfer new value of M2_Aitken to the array
      MOMENT2_CONC( 1 ) = MAX( AEROMODE_MINM2( 1 ), Y )

C *** Accumulation mode: initial value of M2
      M20 = MOMENT2_CONC( 2 )

C *** Production of 2nd moment in accumulation mode is due to
C     intermodal coagulation Aitken mode
      PROD = BACAT( 2 ) * MOMENT0_CONC( 1 ) * MOMENT0_CONC( 2 )

C *** Loss of 2nd moment from accumulation mode is due only to
C     intramodal coagulation
      LOSS = ( BACAC( 2 ) * MOMENT0_CONC( 2 ) * MOMENT0_CONC( 2 ) ) / M20

C *** Solve for M2_accum based on PROD and LOSS during this time step
C     Note: LOSS is assumed to be non-negative.
      IF ( LOSS .GT. 0.0 ) THEN
         POL = PROD / LOSS
         Y = POL + ( M20 - POL ) * EXP( -LOSS * DT )
      ELSE
         Y = M20 + PROD * DT
      END IF ! test on loss

C *** Transfer new value of M2_accum to moment array
      MOMENT2_CONC( 2 ) = MAX( AEROMODE_MINM2( 2 ) , Y  )

C *** Coarse mode: no change because coagulation of coarse particles
C     is neglected in current model version.

! ** Update Budget Variable for Second Moment
      COAG_BUDGET( AEROSRF_MAP( 1 ),1 ) = MOMENT2_CONC( 1 ) - CBUDGET0_SRF( 1 )
      COAG_BUDGET( AEROSRF_MAP( 2 ),2 ) = MOMENT2_CONC( 2 ) - CBUDGET0_SRF( 2 )

C *** end of update for second moment

C *** Update Zeroth Moment (i.e. number concentration)

C *** Aitken mode: initial value of M0

      Y0 = MOMENT0_CONC( 1 )

C *** The rate of change for M0_Aitken is described in Equation 8a of
C     Binkowski & Roselle (2003), with the c_i term equal to 0.

      A = BATAT( 1 )                      ! intramodal coagulation
      B = BATAC( 1 ) * moment0_conc( 2 )  ! intermodal coagulation

      EXPDT = EXP( - B * DT )
      IF ( EXPDT .LT. 1.0D0 ) THEN
         Y = B * Y0 * EXPDT / ( B + A * Y0 * ( 1.0D0 - EXPDT ) )
      ELSE
         Y = Y0                 ! solution in the limit that B approaches zero
      END IF

C *** Transfer new value of M0_Aitken to the moment array
      MOMENT0_CONC( 1 ) = MAX( AEROMODE_MINNUM( 1 ), Y )

C *** Accumulation mode: initial value of M0
      Y0 = MOMENT0_CONC( 2 )

C *** The rate of change for M0_accum is described in Equation 8b of
C     Binkowski & Roselle (2003), except the coefficient C is zero
C     because emissions are treated outside the CMAQ aerosol module.
C     The equation reduces to the form: dY/dt = -A * Y**2 , where
      A = BACAC( 1 )                 ! intramodal coagulation

C *** Solve for M0_accum using Smoluchowski`s solution
      Y = Y0 / ( 1.0D0 + A * Y0 * DT )

C *** Transfer new value of M0_accum to the moment array
       MOMENT0_CONC( 2 ) = MAX( AEROMODE_MINNUM( 2 ), Y )

! ** Update Budget Variable for Second Moment
      COAG_BUDGET( AERONUM_MAP( 1 ),1 ) = MOMENT0_CONC( 1 ) - CBUDGET0_NUM( 1 )
      COAG_BUDGET( AERONUM_MAP( 2 ),2 ) = MOMENT0_CONC( 2 ) - CBUDGET0_NUM( 2 )

C *** end of update for zeroth moment - note that the coarse mode number does
C     not change because coarse-mode coagulation is neglected in the model

C *** UPDATE MASS CONCENTRATIONS (for each species)
C     The following procedure is described in Paragraphs 21-23
C     of Binkowski & Roselle (2003), except the Ei,n and Ej,n terms
C     are excluded here because emissions are treated outside the
C     CMAQ aerosol module.

C     Aitken mode mass concentration rates of change are of the form:
C        dc/dt = P - L*c    ! Equation 9a of Binkowski & Roselle (2003)
C     with solution
C        c(t0 + dt) = P/L + ( c(t0) - P/L ) * exp(-L*dt)

C     For all species, loss of Aitken mode mass is due to intermodal
C     coagulation.
C        LOSSn = PI/6 * RHOn * C30ATAC / MASSn
C        RHOn  = MASSn / (M3 * PI/6)
C     When above equations are combined, the PI/6 terms cancel yielding
C        LOSSn = C30ATAC / M3
C     where LOSSn is the loss rate of species n, RHOn is the mass of
C     species n per unit of particle volume, C30ATAC is the 3rd moment
C     loss rate due to intermodal coagulation, MASSn is the mass
C     concentration of species n, and M3 is the 3rd moment
C     concentration.

      LOSS = C30ATAC / MOMENT3_CONC( 1 )

C *** Set up extra variables to solve for Aitken mode mass concentrations
      FACTRANS = LOSS * DT
      EXPDT = EXP( -FACTRANS )

C *** Transfer mass from Aitken to accumulation mode, resulting from 
C     intermodal coagulation.

      DO SPC = 1, N_AEROSPC
         IF ( AERO_MISSING(SPC,1) ) CYCLE
 
         TMASS = AEROSPC_CONC( SPC,1 ) + AEROSPC_CONC( SPC,2 )
         IF( SPC .EQ. APHGJ_IDX  )THEN 
! assumes all production adsorb onto accumulation mode
             TMASS = TMASS + PHG_RATE * DT
         END IF

         AEROSPC_CONC( SPC,1 ) = MAX( AEROSPC( SPC )%MIN_CONC( 1 ),
     &                                AEROSPC_CONC( SPC,1 ) * EXPDT )
         AEROSPC_CONC( SPC,2 ) = MAX( AEROSPC( SPC )%MIN_CONC( 2 ),
     &                                TMASS - AEROSPC_CONC( SPC,1 ) )


         ! Update Budget Variable for Second Moment
         IF ( AEROSPC_MAP( SPC,1 ) .NE. 0 ) 
     &      COAG_BUDGET( AEROSPC_MAP( SPC,1 ),1 ) = 
     &           AEROSPC_CONC( SPC,1 ) - CBUDGET0_MASS( SPC,1 )
         IF ( AEROSPC_MAP( SPC,2 ) .NE. 0 ) 
     &      COAG_BUDGET( AEROSPC_MAP( SPC,2 ),1 ) = 
     &           AEROSPC_CONC( SPC,2 ) - CBUDGET0_MASS( SPC,2 )

      END DO

C *** end of update for species mass concentrations

C *** Mode Merging
C     This code implements Section 1.5 of Binkowski and Roselle (2003).
C     If the Aitken mode mass is growing faster than accumulation mode
C     mass and the Aitken mode number concentration exceeds the
C     accumulation mode number concentration, then modes are merged by
C     renaming.

! *** Initialize Budget Variables
      CBUDGET0_NUM  = MOMENT0_CONC
      CBUDGET0_SRF  = MOMENT2_CONC
      CBUDGET0_MASS = AEROSPC_CONC

#ifdef nomm
      if ( .false. ) then
#else
      IF ( CGR( 1 ) .GT. CGR( 2 ) .AND.
     &     MOMENT0_CONC( 1 ) .GT. MOMENT0_CONC( 2 ) ) THEN
#endif

C *** Before mode merging, update the third moments, geometric mean
C     diameters, geometric standard deviations, modal mass totals, and
C     particle densities, based on the new concentrations of M2, M0, and
C     speciated masses calculated above. This is still the wet
C     distribution.
         CALL GETPAR( FIXED_sg )

C *** Calculate AAA = ln( Dij / DGATK ) / ( SQRT2 * XXLSGAT ), where Dij
C     is the diameter at which the Aitken-mode and accumulation-mode
C     number distributions intersect (i.e., overlap).  AAA is equivalent
C     to the "Xnum" term described below Equation 10a by Binkowski and
C     Roselle (2003).
         AAA = GETAF( MOMENT0_CONC( 1 ), MOMENT0_CONC( 2 ),
     &                AEROMODE_DIAM( 1 ), AEROMODE_DIAM( 2 ), 
     &                AEROMODE_LNSG( 1 ), AEROMODE_LNSG( 2 ),
     &                SQRT2 ) 

C *** Ensure that Xnum is large enough so that no more than half of
C     the Aitken mode mass is merged into the accumulation mode during
C     any given time step.  This criterion is described in Paragraph 26
C     of Binkowski and Roselle (2003).
         XXM3 = 3.0 * AEROMODE_LNSG( 1 ) / SQRT2
         XNUM = MAX( AAA, XXM3 )

C *** Factors used in error function calls for M2 and M3 mode merging
         XXM2 = TWOTHIRDS * XXM3
         XM2  = XNUM - XXM2 ! set up for 2nd moment transfer
         XM3  = XNUM - XXM3 ! set up for 3rd moment and mass transfers

C *** Calculate the fractions of the number, 2nd, and 3rd moment
C     distributions with diameter greater than the intersection diameter
         FNUM  = 0.5 * ERFC( XNUM )            ! Eq 10a of B&R 2003
         FM2   = 0.5 * ERFC( XM2 )             ! Eq 10b of B&R 2003
         FM3   = 0.5 * ERFC( XM3 )             ! Eq 10b of B&R 2003

C *** Calculate the fractions of the number, 2nd, and 3rd moment
C     distributions with diameters less than the intersection diameter.
         PHNUM = 0.5 * ( 1.0 + ERF( XNUM ) )  ! Eq 10c of B&R 2003
         PHM2  = 0.5 * ( 1.0 + ERF( XM2 ) )   ! Eq 10d of B&R 2003
         PHM3  = 0.5 * ( 1.0 + ERF( XM3 ) )   ! Eq 10d of B&R 2003

C *** Update accumulation-mode moment concentrations using
C     Equations 11a - 11c of Binkowski and Roselle (2003).
         MOMENT0_CONC( 2 ) = MOMENT0_CONC( 2 ) + MOMENT0_CONC( 1 ) * FNUM
         MOMENT2_CONC( 2 ) = MOMENT2_CONC( 2 ) + MOMENT2_CONC( 1 ) * FM2
         MOMENT3_CONC( 2 ) = MOMENT3_CONC( 2 ) + MOMENT3_CONC( 1 ) * FM3

C *** Update Aitken-mode moment concentrations using
C     Equations 11d - 11f of Binkowski and Roselle (2003).
         MOMENT0_CONC( 1 ) = MOMENT0_CONC( 1 ) * PHNUM
         MOMENT2_CONC( 1 ) = MOMENT2_CONC( 1 ) * PHM2
         MOMENT3_CONC( 1 ) = MOMENT3_CONC( 1 ) * PHM3
 
C *** Rename masses of each species from Aitken mode to acumulation mode
C     using Equation 11b of Binkowski and Roselle (2003). Do this for
C     all species, even the aerosol water.
         DO SPC = 1, N_AEROSPC
            IF ( AERO_MISSING(SPC,1) ) CYCLE

            AEROSPC_CONC( SPC,2 ) = AEROSPC_CONC( SPC,2 ) + AEROSPC_CONC( SPC,1 ) * FM3
            AEROSPC_CONC( SPC,1 ) = AEROSPC_CONC( SPC,1 ) * PHM3
         END DO

      END IF ! end check on necessity for merging

C *** end of update for mode merging

C *** Update the third moments, geometric mean diameters, geometric
C     standard deviations, modal mass totals, and particle densities,
C     based on the final concentrations of M2, M0, and speciated masses
C     after mode merging is complete. This should be done for the wet
C     distribution.
      CALL GETPAR( FIXED_sg )

C *** Set minimum value for all concentrations in the CBLK array

      DO N = 1, N_MODE
         DO SPC = 1, N_AEROSPC
            AEROSPC_CONC( SPC,N ) = MAX( AEROSPC_CONC( SPC,N ),
     &                                   AEROSPC( SPC )%MIN_CONC( N ) )
         END DO
      END DO

! *** Propagate Mode Merging and Minimum value impacts to budget variables
	  
C     GROWTH_BUDGET( AERONUM_MAP(:) ) = MOMENT0_CONC(:) - CBUDGET0_NUM(:)
C     GROWTH_BUDGET( AEROSRF_MAP(:) ) = MOMENT0_CONC(:) - CBUDGET0_NUM(:)
C     DO SPC = 1,N_AEROSPC
C        WHERE( AEROSPC_MAP( SPC,: ) .NE. 0 )
C    &          GROWTH_BUDGET( AEROSPC_MAP(SPC,:) )=
C    &              AEROSPC_CONC(SPC,:) - CBUDGET0_MASS(SPC,:)
C     END DO
	  
      DO MODE = 1, N_MODE
 	      GROWTH_BUDGET(AERONUM_MAP(MODE)) = MOMENT0_CONC(MODE) - CBUDGET0_NUM(MODE)
 	      GROWTH_BUDGET(AEROSRF_MAP(MODE)) = MOMENT0_CONC(MODE) - CBUDGET0_NUM(MODE)
      ENDDO 
   
      DO SPC = 1,N_AEROSPC
            DO MODE = 1,N_MODE
                  IF (AEROSPC_MAP(SPC, MODE) .NE. 0) THEN
                  GROWTH_BUDGET( AEROSPC_MAP(SPC, MODE) ) = AEROSPC_CONC(SPC, MODE) - CBUDGET0_MASS(SPC, MODE)
 	            ENDIF 
            ENDDO
      ENDDO 
	  
		  
      RETURN
      END SUBROUTINE AEROPROC


C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE VOLINORG( DT, COL, ROW, LAYER, DV_SO4, CGR )

C *** Calculates the partitioning of inorganic components (CL,NO3,NH4,SO4)
C     between the aerosol and gas phase over the operator synchronization
C     timestep (DT). Partitioning is calculated using the Hybrid approach,
C     where dynamic mass transfer of species to/from the coarse mode is
C     calculated using multiple sub-operator time steps (TSTEP) and the
C     fine modes are equilibrated with the gas phase. The mass transfer
C     calculations are made using the H+ flux-limiting approach of Pilinis
C     et al. (2000). If 'OPTIONFLAG' is not set to 'Hybrid', the mass
C     transfer calculations for the coarse mode are skipped, and the fine
C     modes are equilibrated with the gas phase.

C     Returns updated volatile inorganic species concentrations in the gas
C     and particulate phase, and the aerosol modal parameters

C *** Revision history: 4/07 - Moved HCOND3 and NEWPART3 calls from 
C                              AEROPROC to this subroutine for 
C                              mass transfer calculation  
C     15 Jul 08 J.Young, P.Bhave: increased cutoff to hybrid from .01 to .05 ug/m**3
C               J.Young: change 'OPTIONFLAG' to just a logical variable, 'Hybrid'

C     10 Mar 11 S.Howard Renamed met_data to aeromet_data
C     25 Mar 11 S.Roselle Replaced I/O API include files with UTILIO_DEFN
C     26 Apr 11 G.Sarwar Replaced existing ISORROPIA with ISORROPIA 2.1
C               Updated coarse-mode aerosol speciation and H+ calculation  
C     12 Apr 16 H. Pye and B. Murphy: Update with consistent treatment for wet particles.
C     15 Apr 16 J.Young: Use aerosol factors from AERO_DATA module named constants

C *** References
C 1. Pilinis C, Capaldo KP, Nenes A, Pandis SN (2000) MADM - A new
C    multicomponent aerosol dynamics model. AEROSOL SCIENCE AND TECHNOLOGY.
C    32(5):482-502
C
C 2. Capaldo KP, Pilinis C, Pandis SN (2000) A computationally efficient hybrid
C    approach for dynamic gas/aerosol transfer in air quality models. ATMOSPHERIC
C    ENVIRONMENT. 34(21):3617-3627
C
C 3. Fountoukis C, Nenes, A (2007) ISORROPIA II: a comnputationally efficient 
C    thermodynamic equilibrium model for K+-Ca2+-Mg2+-NH4+-SO42--NO3--Cl-H2O 
C    aerosols. ATMOSPHERIC CHEMISTRY AND PHYSICS. 7, 4639-4659

      USE AERO_DATA
      USE PRECURSOR_DATA
      USE SOA_DEFN
      USE AEROMET_DATA
      USE UTILIO_DEFN
	  USE HDMod 
   	  USE ISORROPIAII_MAIN_MOD, ONLY : ISOROPIA
	  
      IMPLICIT NONE

C *** Arguments:
      REAL    DT              ! time step [sec]
      INTEGER COL             ! grid column index
      INTEGER ROW             ! grid row index
      INTEGER LAYER           ! model layer index
      REAL    DV_SO4          ! molecular diffusivity of H2SO4 vapor 
                              ! after correction for ambient conditions
      REAL( 8 ) :: CGR( N_MODE-1 ) ! 3rd moment SO4 growth rate [m^3/m^3-s]

C *** Parameters: 
      INTEGER, PARAMETER :: NINORG = 9      ! number of inorganic species
      INTEGER, PARAMETER :: NVOLINORG = 3   ! number of volatile inorganic species

      ! indices for inorganic species
!     INTEGER, PARAMETER :: KNH4 = 1, KNO3 = 2, KCL = 3, KSO4 = 4, KNA = 5, KHP = 6, KMG = 7, KK = 8, KCA = 9
      INTEGER, PARAMETER :: KNH4 = 1, KNO3 = 2, KCL = 3, KSO4 = 4, KNA = 5, KMG = 6, KK = 7, KCA = 8, KHP = 9

      REAL( 8 ), PARAMETER :: D_TWOTHIRDS = 2.0D0 / 3.0D0

      REAL, PARAMETER :: CUTOFF = 0.05  ! [ug/m**3]
      REAL, PARAMETER :: ALPHSULF = 0.1 ! Accommodation coefficient for sulfuric acid
                                        ! see Capaldo et al. (2000)

C *** Local Variables:

C *** Inputs to subroutine HCOND3

      REAL, SAVE :: COFCBAR_SO4  ! Temperature-independent coefficients
                                 ! for caculating molecular vel [m/s]
                                 ! = sqrt((8*Rgas)/(pi*MW)) 
      REAL         CBAR_SO4      ! molecular velocity of H2SO4                      

      TYPE(hyperdual) :: AM0( N_MODE ) ! zeroth moments
      TYPE(hyperdual) :: AM1( N_MODE ) ! first moments
      TYPE(hyperdual) :: AM2( N_MODE ) ! second moments

      TYPE(hyperdual) :: M2DRY_INIT( N_MODE )  ! Dry Second Moment Initial Array
      TYPE(hyperdual) :: M3DRY_INIT( N_MODE )  ! Dry Third Moment Initial Array
      TYPE(hyperdual) :: M2WET_INIT( N_MODE )  ! Wet Second Moment Initial Array
      TYPE(hyperdual) :: M3WET_INIT( N_MODE )  ! Wet Third Moment Initial Array
      TYPE(hyperdual) :: M2WET_FINAL( N_MODE ) ! Wet Second Moment Final Array
      TYPE(hyperdual) :: M3WET_FINAL( N_MODE ) ! Wet Third Moment Final Array

C *** Outputs from HCOND3: size-dependent term in the condensational-growth 
C     expressions defined in Equations A13-A14 of [Binkowski & Shankar,1995]
      TYPE(hyperdual) :: FCONC_SO4( N_MODE,2 )  ! All sizes 2nd and 3rd moments
      TYPE(hyperdual) :: FCONC_OUT( 2 )         ! Output for HCOND3, 2nd and 3rd moments
      TYPE(hyperdual) :: FCONCM1_SO4       ! reciprocals of total SO4 cond rates

C *** Modal partition factors [ dimensionless ]
C     defined in Equations A17-A18 of [Binkowski & Shankar,1995]
      TYPE(hyperdual) :: OMEGA_AT_SO4  ! Aitken mode 2nd and 3rd moments
      TYPE(hyperdual) :: OMEGA_AC_SO4  ! Accumulation mode 2nd and 3rd moments
      TYPE(hyperdual) :: OMEGA( 2 )    ! partitioning coefficient for equilibrium PM mass
      TYPE(hyperdual) :: PHI( NINORG,2 ) ! Mass Fraction of each component in each aerosol mode
      TYPE(hyperdual) :: TOTAER( NINORG )! Total aerosol component across multiple modes 

C *** Variables for new particle formation:
      TYPE(hyperdual) :: XH2SO4    ! steady state H2SO4 concentration
      TYPE(hyperdual) :: DMDT_SO4  ! particle mass production rate [ ug/m**3 s ]
      TYPE(hyperdual) :: DNDT      ! particle number production rate [ # / m**3 s ]
      TYPE(hyperdual) :: DM2DT     ! second moment production rate [ m**2 / m**3 s]
      TYPE(hyperdual) :: SCONDRATE ! SO4 condensation rate [ ug/m**3 s ]

C *** Mode-specific sulfate production rate [ ug/m**3 s ]
      TYPE(hyperdual) :: CONDSO4( N_MODE )    ! sulfate condensation rate [ ug/m**3 s ]
      TYPE(hyperdual) :: RATE                 ! CONDSO4 or cond+nucl rate

C *** Size-dependent portion of mass-transfer rate equation
      TYPE(hyperdual) :: GRFAC1( N_MODE )     ! 2nd moment [ m**2/m**3-s ] 
      TYPE(hyperdual) :: GRFAC2( N_MODE )     ! 3rd moment [ m**3/m**3-s ] 
      
C *** ISORROPIA input variables
      TYPE(hyperdual) :: WI( NINORG - 1 )     ! species array [ mol/m**3 ]
      REAL( 8 ) :: RHI                  ! relative humidity [ fraction ]
      REAL( 8 ) :: TEMPI                ! temperature   [ deg K]
      REAL( 8 ) :: CNTRL( 2 )           ! ISOROPIA control parameters 

C *** ISORROPIA output variables
      TYPE(hyperdual) :: WT( NINORG - 1 )     ! species output array (unused)
      TYPE(hyperdual) :: GAS( 3 )             ! gas-phase   "     " 
      TYPE(hyperdual) :: AERLIQ( 15 )         ! liq aerosol "     " 
      TYPE(hyperdual) :: AERSLD( 19 )         ! solid "     "     "  (unused)
      TYPE(hyperdual) :: OTHER( 9 )           ! supplmentary output array (unused)
      CHARACTER( 15 ) :: SCASI          ! subcase number output (unused)

C *** Variables to account for mass conservation violations in ISRP3F
!     LOGICAL TRUSTNH4                  ! false if ISOROPIA's partitioning
                                        !  of NH4/NH3 is to be ignored
      LOGICAL TRUSTCL                   ! false if ISOROPIA's partitioning       
                                        !  of Cl/HCl is to be ignored

C *** Initial (double-precision) concentrations [ug/m3]
      TYPE(hyperdual) :: GNH3R8               ! gas-phase ammonia
      TYPE(hyperdual) :: GNO3R8               ! gas-phase nitric acid
      TYPE(hyperdual) :: GCLR8                ! gas-phase hydrochloric acid

C *** Variables for volatile species mass transfer between gas and aerosol and
C     mass partitioning between the modes 
      LOGICAL         HYBRID ! mass transfer option flag (mass transfer if .TRUE.)
      REAL( 8 )       :: DELT                 ! time step DT [s]
      TYPE(hyperdual) :: HPLUS( N_MODE )      ! scratch var for H+ [umol/m**3]

      TYPE(hyperdual), SAVE :: H2SO4RATM1     ! Mol. wt. ratio of SO4/H2SO4

      TYPE(hyperdual) :: DVOLINORG( NVOLINORG ) ! vol inorg spcs mass to be xferred [mol/m3]
      TYPE(hyperdual) :: DVOLMAX                ! max value for DVOLINORG 
      TYPE(hyperdual) :: CINORG( NINORG,N_MODE ) ! scratch array for inorg spcs [ug/m**3]
      TYPE(hyperdual) :: SEACAT                 ! coarse sea-salt cations [ug/m**3]
      TYPE(hyperdual) :: SOILwVOLS              ! windblown dust before removal of
                                          ! SO4,NO3,CL,H2O [ug/m**3]
      TYPE(hyperdual) :: PMCwVOLS               ! anthrop coarse material before removal
                                          ! of SO4,NO3,CL,H2O [ug/m**3]
      REAL( 8 ) :: INT_TIME               ! internal mass transfer time (s)
      REAL( 8 ) :: TSTEP                  ! mass transfer time step [s]
      TYPE(hyperdual) :: DRYM20, Y              ! scratch vars for 2nd moment [m**2/m**3]
      TYPE(hyperdual) :: M3OTHR                 ! vars for 3rd moment calculation [m**3/m**3]
      REAL( 8 ), SAVE :: DF( NVOLINORG )  ! scratch array for mole -> ug conversion
      REAL( 8 ), SAVE :: DFH2OR8          ! mole -> ug conversion for H2O
      TYPE(hyperdual) :: J( NVOLINORG )         ! condensation/evaporation flux [mol/m**3-s]
      TYPE(hyperdual) :: CFINAL( NVOLINORG,N_MODE ) ! conc after mass xfer step [ug/m**3]
      TYPE(hyperdual) :: H2O                    ! Scratch LWC variable for output
      TYPE(hyperdual) :: H2O_NEW                ! Update of LWC after condensation 
      TYPE(hyperdual) :: SO4                    ! modal SO4 after condensation or cond + nucl
      TYPE(hyperdual) :: DDRYM3DT               ! rate of 3rd moment transfer - dry inorg spcs
      TYPE(hyperdual) :: DDRYM2DT               ! rate of 2nd moment transfer -  "     "    "
      TYPE(hyperdual) :: DRYM3, WETM3           ! scratch vars for 3rd moment calc [m**3/m**3]
      TYPE(hyperdual) :: DRYM2, WETM2           ! scratch vars for 2nd moment calc [m**2/m**3]
      TYPE(hyperdual) :: LOSS                   ! rate of loss of second moment [1/s] 
      TYPE(hyperdual) :: EQLBHIJ                ! H+ concentration from isoropia for I plus J-modes [ug/m**3]

      TYPE(hyperdual) :: DELNUM( N_MODE ) ! Change in Number due to some process
      TYPE(hyperdual) :: DELSRF( N_MODE ) ! Change in Surface Area due to some process
      TYPE(hyperdual) :: DELSO4( N_MODE ) ! Change in Sulfate due to some process
      TYPE(hyperdual) :: DELNO3( N_MODE ) ! Change in Nitrate due to some process
      TYPE(hyperdual) :: DELNH4( N_MODE ) ! Change in Ammonium due to some process
      TYPE(hyperdual) :: DELCL( N_MODE )  ! Change in Chloride due to some process
      TYPE(hyperdual) :: DELH2O( N_MODE ) ! Change in Water due to some process
      TYPE(hyperdual) :: DELH3OP( N_MODE )! Change in Proton due to some process
      TYPE(hyperdual) :: CBUDGET0_NUM( N_MODE ) ! Initial Number before some process
      TYPE(hyperdual) :: CBUDGET0_SRF( N_MODE ) ! Initial Surface Area before some process
      TYPE(hyperdual) :: CBUDGET0_SO4( N_MODE ) ! Initial Sulfate before some process
      TYPE(hyperdual) :: CBUDGET0_NO3( N_MODE ) ! Initial Nitrate before some process
      TYPE(hyperdual) :: CBUDGET0_NH4( N_MODE ) ! Initial Ammonium before some process
      TYPE(hyperdual) :: CBUDGET0_CL( N_MODE )  ! Initial Chloride before some process
      TYPE(hyperdual) :: CBUDGET0_H2O( N_MODE ) ! Initial Particle Water before some process
      TYPE(hyperdual) :: CBUDGET0_H3OP( N_MODE )! Initial Proton before some process

      TYPE(hyperdual) :: TMP
      REAL( 8 )       :: TEMP_CUTOFF      ! Cutoff temperature for whether to use FwdIso
      REAL( 8 )       :: PRES_CUTOFF      ! Cutoff pressure for whether to use FwdIso
      INTEGER I, N                        ! loop and array indices
      INTEGER IMODE                       ! mode loop index  
      INTEGER MODE                        ! mode loop index_el
      INTEGER ISTEP                       ! loop index, mass transfer time step loop
      INTEGER ISP                         ! loop index, species loop
      LOGICAL TrustIso                    ! For negative vap. press., TrustIso = F
      LOGICAL FwdIso                      ! J.E.L. for not using ISORROPIA  
      LOGICAL HybridSens                  ! J.E.L. for not using rev ISORROPIA hyd parts
      LOGICAL FwdSens                     ! 
      
      ! Additional variables for resetting sensitivities for asymptotic steps
      TYPE(hyperdual)       :: HPLUS_PREV
      TYPE(hyperdual)       :: moment0_conc3_orig
      TYPE(hyperdual)       :: moment2_conc3_orig
      TYPE(hyperdual)       :: moment3_conc3_orig
      TYPE(hyperdual)       :: aeromode_dens3_orig
      TYPE(hyperdual)       :: aeromode_lnsg3_orig
      TYPE(hyperdual)       :: aeromode_diam3_orig
      TYPE(hyperdual)       :: aeromode_mass3_orig
      TYPE(hyperdual)       :: anh4k_orig
      TYPE(hyperdual)       :: ano3k_orig
      TYPE(hyperdual)       :: aclk_orig
      TYPE(hyperdual)       :: aso4k_orig
      TYPE(hyperdual)       :: ah2ok_orig
      TYPE(hyperdual)       :: nh3_orig
      TYPE(hyperdual)       :: hno3_orig
      TYPE(hyperdual)       :: hcl_orig
      TYPE(hyperdual)       :: hplus_orig
      TYPE(hyperdual)       :: seacat_orig
      TYPE(hyperdual)       :: soil_orig
      TYPE(hyperdual)       :: cors_orig
      
      
      
C *** Local Saved Variables 
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      LOGICAL, SAVE :: FIRSTWRITE = .TRUE.
      REAL( 8 ), SAVE :: SO4FAC         !  F6DPIM9 / RHOSO4
      REAL( 8 ), SAVE :: SOILFAC        !  F6DPIM9 / RHOSOIL
      REAL( 8 ), SAVE :: ANTHFAC        !  F6DPIM9 / RHOANTH
      REAL( 8 ), SAVE :: H2SO4RAT       ! Mol. wt ratios H2SO4/SO4
      REAL( 8 ), SAVE :: NH3RAT         ! Mol. wt ratios NH3/NH4
      REAL( 8 ), SAVE :: HNO3RAT        ! Mol. wt ratios HNO3/NO3
      REAL( 8 ), SAVE :: HCLRAT         ! Mol. wt ratios HCL/CL 
      REAL( 8 ), SAVE :: MWH2SO4        ! molecular weight for H2SO4
      REAL( 8 ), SAVE :: FAERH2SO4      ! 1e-6 / mw(h2so4)

      REAL( 8 ), SAVE :: CFAC_ANA
      REAL( 8 ), SAVE :: CFAC_ASO4
      REAL( 8 ), SAVE :: CFAC_ANH4
      REAL( 8 ), SAVE :: CFAC_ANO3
      REAL( 8 ), SAVE :: CFAC_ACL
      REAL( 8 ), SAVE :: CFAC_ACA
      REAL( 8 ), SAVE :: CFAC_AK
      REAL( 8 ), SAVE :: CFAC_AMG
      REAL( 8 ), SAVE :: CFAC_ASEACAT
      REAL( 8 ), SAVE :: CFAC_GNH3
      REAL( 8 ), SAVE :: CFAC_GHNO3
      REAL( 8 ), SAVE :: CFAC_GHCL

      REAL( 8 ), SAVE :: M3FAC_ANA
      REAL( 8 ), SAVE :: M3FAC_ASO4
      REAL( 8 ), SAVE :: M3FAC_ANH4
      REAL( 8 ), SAVE :: M3FAC_ANO3
      REAL( 8 ), SAVE :: M3FAC_ACL
      REAL( 8 ), SAVE :: M3FAC_H2O
	  
      INTEGER         :: C_GAS      ! loop counter (E.L.)
      LOGICAL         :: SAME_ALL_AFTER = .false.
	  
	  
      logical, save :: write1 = .true.

C-----------------------------------------------------------------------
 
      IF ( FIRSTIME ) THEN
      
         FIRSTIME = .FALSE.
#ifdef noiso
         write( logdev,* ) 'aero_subs, volinorg: noiso'
#endif

         CFAC_ANA  = 1.0D-6 / REAL( AEROSPC_MW( ANA_IDX ), 8 )
         CFAC_ASO4 = 1.0D-6 / REAL( AEROSPC_MW( ASO4_IDX ), 8 )
         CFAC_ANH4 = 1.0D-6 / REAL( AEROSPC_MW( ANH4_IDX ), 8 )
         CFAC_ANO3 = 1.0D-6 / REAL( AEROSPC_MW( ANO3_IDX ), 8 )
         CFAC_ACL  = 1.0D-6 / REAL( AEROSPC_MW( ACL_IDX ), 8 )
         CFAC_ACA  = 1.0D-6 / REAL( AEROSPC_MW( ACA_IDX ), 8 )
         CFAC_AK   = 1.0D-6 / REAL( AEROSPC_MW( AK_IDX ), 8 )
         CFAC_AMG  = 1.0D-6 / REAL( AEROSPC_MW( AMG_IDX ), 8 )
         CFAC_ASEACAT = 1.0D-6 / REAL( AEROSPC_MW( ASEACAT_IDX ), 8 )
         CFAC_GNH3 = 1.0D-6 / PRECURSOR_MW( NH3_IDX )
         CFAC_GHNO3= 1.0D-6 / PRECURSOR_MW( HNO3_IDX )
         CFAC_GHCL = 1.0D-6 / PRECURSOR_MW( HCL_IDX )

         M3FAC_ANA  = 1.0D-9 * F6DPI / REAL( AEROSPC( ANA_IDX )%DENSITY, 8 )
         M3FAC_ASO4 = 1.0D-9 * F6DPI / REAL( AEROSPC( ASO4_IDX )%DENSITY, 8 )
         M3FAC_ANH4 = 1.0D-9 * F6DPI / REAL( AEROSPC( ANH4_IDX )%DENSITY, 8 )
         M3FAC_ANO3 = 1.0D-9 * F6DPI / REAL( AEROSPC( ANO3_IDX )%DENSITY, 8 )
         M3FAC_ACL  = 1.0D-9 * F6DPI / REAL( AEROSPC( ACL_IDX )%DENSITY, 8 )
         M3FAC_H2O  = 1.0D-9 * F6DPI / REAL( AEROSPC( AH2O_IDX )%DENSITY, 8 )

         MWH2SO4 = PRECURSOR_MW( SULF_IDX ) ! molecular weight for H2SO4
         FAERH2SO4 = 1.0E-6 / MWH2SO4

         COFCBAR_SO4 = SQRT( 8.0 * RGASUNIV / ( PI * REAL( MWH2SO4,4 ) * 1.0E-3 ) )
         H2SO4RATM1 = AEROSPC_MW( ASO4_IDX ) / MWH2SO4
         SO4FAC  = 1.0D-9 * F6DPI / REAL( AEROSPC( ASO4_IDX )%DENSITY, 8 )
         SOILFAC = 1.0D-9 * F6DPI / REAL( AEROSPC( ASOIL_IDX )%DENSITY, 8 )
         ANTHFAC = 1.0D-9 * F6DPI / REAL( AEROSPC( ACORS_IDX )%DENSITY, 8 )

         DF( KNH4 ) = 1.0D0 / CFAC_ANH4
         DF( KNO3 ) = 1.0D0 / CFAC_ANO3
         DF( KCL )  = 1.0D0 / CFAC_ACL
         DFH2OR8    = 1.0D6 * MWWAT      ! aerospc_mw(AH2O_IDX) (18.0 != 18.0153)

         ! Mol. wt ratios H2SO4/SO4, NH3/NH4, HNO3/NO3, HCL/CL
         H2SO4RAT= PRECURSOR_MW( SULF_IDX ) / REAL( AEROSPC_MW( ASO4_IDX ), 8 )
         NH3RAT  = PRECURSOR_MW( NH3_IDX )  / REAL( AEROSPC_MW( ANH4_IDX ), 8 )
         HNO3RAT = PRECURSOR_MW( HNO3_IDX ) / REAL( AEROSPC_MW( ANO3_IDX ), 8 )
         HCLRAT  = PRECURSOR_MW( HCL_IDX )  / REAL( AEROSPC_MW( ACL_IDX ), 8 )
      END IF

C *** Determine if Hybrid
      TMP = 0.0D0
      
#ifdef noiso
	hybrid = .false.
#else
	  TEMP_CUTOFF = 260.d0
      PRES_CUTOFF = 2.D4
      FWDISO = ( AIRTEMP .GE. TEMP_CUTOFF ) .AND. ( AIRPRES .GE. PRES_CUTOFF ) .AND. (AIRRH .GE. 0.18)
C     FWDISO = .TRUE.
      
      DO I = 1, N_AEROSPC
         IF(  AEROSPC( I )%tracer )CYCLE
         IF ( AEROSPC( I )%CHARGE .NE. 0 ) TMP = TMP + AEROSPC_CONC( I,N_MODE )
      END DO
      HYBRID = ( TMP .GE. CUTOFF ) .AND. ( AIRRH .GE. 0.18 ) .AND. FWDISO
#endif

      DELT  = REAL( DT, 8 )
      TEMPI = REAL( AIRTEMP , 8 )             ! assume const within synch timestep
      RHI   = MIN( 0.95D0, REAL ( AIRRH, 8 ) )  ! "        "     "      "     "

C *** Calculate molecular velocities (temperature dependent) and
C     H+ concentration

      CBAR_SO4 = COFCBAR_SO4 * SQRT( AIRTEMP )

      HPLUS = 0.0D0
      DO I = 1, N_MODE
         DO N = 1, N_AEROSPC
            IF(  AEROSPC( N )%tracer )CYCLE
            HPLUS( I ) = HPLUS( I )
     &                 - AEROSPC( N )%CHARGE * AEROSPC_CONC( N,I ) / AEROSPC_MW( N )
         END DO
      END DO

C *** Condensational Growth (Size-dependent terms)
C     Calculate intermediate variables needed to determine the 2nd and
C     3rd moment condensational-growth rates.  3rd moment terms are 
C     needed for the calculation of new particle production.  See 
C     Section 3.3 of Jiang & Roth (2003) for a detailed discussion.
C    
C *** Set moments using Equation 4 of Binkowski & Shankar
C     (1995) or Equation 3 of Binkowski and Roselle (2003).
C     N.B: these are for a "wet" size distribution

      DO I = 1, N_MODE
         AM0( I ) = MOMENT0_CONC( I ) 
         AM1( I ) = MOMENT0_CONC( I ) * AEROMODE_DIAM( I )
     &            * EXP( 0.5 * AEROMODE_LNSG( I ) * AEROMODE_LNSG( I ) )
         AM2( I ) = MOMENT2_CONC( I )
      END DO
      
C *** Calculate the size-dependent terms in the condensational-
C     growth factor expressions for sulfate using 
C     Equations A13-A14 of Binkowski & Shankar (1995). 
       
      DO I = 1, N_MODE
         CALL HCOND3( AM0( I ), AM1( I ), AM2( I ),
     &                DV_SO4, ALPHSULF, CBAR_SO4, FCONC_OUT )
         FCONC_SO4( I, : ) = FCONC_OUT( : )
      END DO 

      IF ( .NOT. HYBRID ) THEN
         FCONC_SO4( N_MODE,1 ) = 0.0D0
         FCONC_SO4( N_MODE,2 ) = 0.0D0
      END IF

      DO I = 1, N_MODE
         GRFAC1( I ) = FCONC_SO4( I,1 )
         GRFAC2( I ) = FCONC_SO4( I,2 )
      END DO

C *** New Particle Production
C     Calculate the new particle production rate due to binary
C     nucleation of H2O and H2SO4.  These calculations are performed 
C     only when the gas-phase production rate of H2SO4 (i.e., SO4RATE) 
C     is non-zero.  The condensation rate of H2SO4 is calculated as the
C     gas-phase production rate minus the new particle production rate.

C *** Initialize Variables
      DMDT_SO4  = 0.0D0
      DNDT      = 0.0D0
      DM2DT     = 0.0D0
      SCONDRATE = 0.0D0

C *** Produce new particles only during time steps when the gas-phase 
C     production rate of H2SO4 is non-zero

      IF ( SO4RATE .NE. 0.0D0 ) THEN

C *** Adjust sulfuric acid vapor concentration to a value in
C     equilibrium with the production of new particles and the
C     condensation of sulfuric acid vapor on existing particles, based 
C     on Equations A21 and A23 of Binkowski & Shankar (1995).
         TMP = 0.0d0
         DO I = 1, N_MODE
            TMP = TMP + FCONC_SO4( I,2 )
         END DO

         XH2SO4 = SO4RATE / TMP
         XH2SO4 = MAX( XH2SO4, CONMIN )
         PRECURSOR_CONC( SULF_IDX ) = XH2SO4 

C *** Calculate new particle production rate for 0th, 2nd, & 3rd moments
         CALL NEWPART3 ( AIRRH, AIRTEMP, XH2SO4, SO4RATE,
     &                   DNDT, DMDT_SO4, DM2DT )
         
C *** Calculate sulfate condensation rate as the gas-phase production 
C     rate minus the new particle production rate, following Equation
C     3.23 of Jiang & Roth (2003).
         SCONDRATE = MAX( SO4RATE - DMDT_SO4, 0.0D0 )

! *** Propagate NPF Change to budget variables
         NPF_BUDGET( PRECURSOR_MAP( SULF_IDX ) ) = -DMDT_SO4 
     &             * H2SO4RAT  * DELT ! ug m-3 of Sulfuric Acid
                                      ! (converted to ppmv in update_precursor)
         NPF_BUDGET( AEROSPC_MAP( ASO4_IDX,1 ) ) =  DMDT_SO4 * DELT ! ug m-3
         NPF_BUDGET( AERONUM_MAP( 1 ) ) =  DNDT * DELT  ! N m-3
         NPF_BUDGET( AEROSRF_MAP( 1 ) ) =  DM2DT * DELT ! m2 m-3

      END IF   ! SO4RATE .NE. 0

C *** Sulfate Condensation (Size-resolved)
C     Calculate rate at which condensing sulfate should be added to each
C     mode.  The "omega" factors are defined in Equations 7a and 7b of
C     Binkowski & Shankar (1995). The i-mode and j-mode factors are 
C     calculated using Equation A17 of Binkowski & Shankar (1995). The 
C     condensation rate for accumulation mode (fine-equilibrium scheme) or 
C     coarse mode (hybrid and dynamic schemes) is computed by difference, 
C     to avoid mass conservation violations arising from numerical error.
      
      TMP = 0.d0

      DO I = 1, N_MODE
         TMP = TMP + FCONC_SO4( I,2 )
      END DO

      FCONCM1_SO4  = 1.0D0 / TMP
      OMEGA_AT_SO4 = FCONCM1_SO4 * FCONC_SO4( 1,2 )
      OMEGA_AC_SO4 = FCONCM1_SO4 * FCONC_SO4( 2,2 )

C *** Growth values for mode merge condition
      CGR( 1 ) = SO4FAC * SCONDRATE%x * OMEGA_AT_SO4%x
      CGR( 2 ) = SO4FAC * SCONDRATE%x * OMEGA_AC_SO4%x

! *** Initialize Budget Counter Variables         
      CBUDGET0_NUM( : ) = MOMENT0_CONC( : )
      CBUDGET0_NUM( 1 ) = CBUDGET0_NUM( 1 ) + DNDT * DELT
      CBUDGET0_SRF( : ) = MOMENT2_CONC( : )
      CBUDGET0_SRF( 1 ) = CBUDGET0_SRF( 1 ) + DM2DT * DELT
      CBUDGET0_SO4( : ) = AEROSPC_CONC( ASO4_IDX, : )
      CBUDGET0_SO4( 1 ) = CBUDGET0_SO4( 1 ) + DMDT_SO4 * DELT
      CBUDGET0_NO3( : ) = AEROSPC_CONC( ANO3_IDX, : )
      CBUDGET0_NH4( : ) = AEROSPC_CONC( ANH4_IDX, : )
      CBUDGET0_CL ( : ) = AEROSPC_CONC( ACL_IDX , : )
      CBUDGET0_H2O( : ) = AEROSPC_CONC( AH2O_IDX, : )
      CBUDGET0_H3OP( : )= AEROSPC_CONC( AH3OP_IDX, : )

! *** Diagnose Condensaiton Sink to Each Mode
      CONDSO4( 1 ) = OMEGA_AT_SO4 * SCONDRATE      
      
      IF ( HYBRID ) THEN 
         CONDSO4( 2 ) = OMEGA_AC_SO4 * SCONDRATE      
         CONDSO4( 3 ) = SCONDRATE - ( CONDSO4( 1 ) + CONDSO4( 2 ) ) 
      ELSE                                  ! fine equilibrium
         CONDSO4( 2 ) = SCONDRATE - CONDSO4( 1 )
         CONDSO4( 3 ) = 0.0D0               ! no coarse mode chemistry
      END IF

C *** For Hybrid approach, calculate dynamic mass trasfer for
C     semi-volatile components of coarse mode (NO3, CL, NH4)   
      HybridSens = .TRUE. 
      IF ( HYBRID ) THEN 

         CNTRL( 1 ) = 1.0D0 ! reverse problem
         CNTRL( 2 ) = 1.0D0 ! aerosol in metastable state
		 
         INT_TIME = 0.0D0
         TSTEP    = 90.0D0
         
         ISTEP    = 1
         IMODE    = 3
         TrustIso = .TRUE.
         
 		 
         DO WHILE ( INT_TIME .LT. DELT ) 

            IF ( INT_TIME + TSTEP .GT. DELT ) TSTEP = DELT - INT_TIME 
            INT_TIME = INT_TIME + TSTEP 
            ISTEP = ISTEP + 1   
            		
 			IF ( INT_TIME .GT. 90.0D0 ) THEN
 				HPLUS_PREV = AERLIQ( 1 )
 			ENDIF
			
C *** Calculate first moments using Equation 4 of Binkowski & Shankar
C     (1995) or Equation 3 of Binkowski and Roselle (2003).
C     N.B: these are for a "wet" size distribution
            AM0( IMODE ) = MOMENT0_CONC( IMODE )
            AM1( IMODE ) = MOMENT0_CONC( IMODE ) * AEROMODE_DIAM( IMODE )
     &                  * EXP( 0.5 * AEROMODE_LNSG( IMODE ) * AEROMODE_LNSG( IMODE ) )
            AM2( IMODE ) = MOMENT2_CONC( IMODE )

C *** Calculate the size-dependent terms in the condensational-
C     growth factor expressions for sulfate using 
C     Equations A13-A14 of Binkowski & Shankar (1995). 
            CALL HCOND3( AM0( IMODE ), AM1( IMODE ), AM2( IMODE ), DV_SO4, ALPHSULF, 
     &                      CBAR_SO4, FCONC_OUT )  ! adapted from Eq A14
            FCONC_SO4( IMODE, : ) = FCONC_OUT( : )

            GRFAC1( IMODE ) = FCONC_SO4( IMODE,1 ) 
            GRFAC2( IMODE ) = FCONC_SO4( IMODE,2 ) 

C *** Set conc array to aerosol concentrations prior to mass transfer
C     (The _RENORM and cation_FAC constants are set in the AERO_DATA module)

            SEACAT    = AEROSPC_CONC( ASEACAT_IDX,IMODE )              
            SOILwVOLS = AEROSPC_CONC( ASOIL_IDX,IMODE ) / ASOIL_RENORM 
            PMCwVOLS  = AEROSPC_CONC( ACORS_IDX,IMODE ) / ACORSEM_RENORM 

            CINORG( KNH4,IMODE ) = AEROSPC_CONC( ANH4_IDX,IMODE )  ! KNH4 = 1, ANH4_IDX = 4
            CINORG( KNO3,IMODE ) = AEROSPC_CONC( ANO3_IDX,IMODE )  ! KNO3 = 2, ANO3_IDX = 2
            CINORG( KCL, IMODE ) = AEROSPC_CONC( ACL_IDX,IMODE )   ! KCL = 3, ACL_IDX = 3
            CINORG( KSO4,IMODE ) = AEROSPC_CONC( ASO4_IDX,IMODE ) ! KSO4 = 4, ASO4_IDX = 1
            CINORG( KNA, IMODE ) = ASCAT_NA_FAC * SEACAT                      ! KNA = 5
     &                           + ASOIL_NA_FAC * SOILwVOLS
     &                           + ACORS_NA_FAC * PMCwVOLS
            CINORG( KMG, IMODE ) = ASCAT_MG_FAC * SEACAT                      ! KMG = 6
     &                           + ASOIL_MG_FAC * SOILwVOLS
     &                           + ACORS_MG_FAC * PMCwVOLS
            CINORG( KK,  IMODE ) = ASCAT_K_FAC  * SEACAT                      ! KK = 7
     &                           + ASOIL_K_FAC  * SOILwVOLS
     &                           + ACORS_K_FAC  * PMCwVOLS
            CINORG( KCA, IMODE ) = ASCAT_CA_FAC * SEACAT                      ! KCA = 8
     &                           + ASOIL_CA_FAC * SOILwVOLS
     &                           + ACORS_CA_FAC * PMCwVOLS
            CINORG( KHP, IMODE ) = HPLUS( IMODE )                  ! KHP = 9

            M3OTHR = SOILFAC * AEROSPC_CONC( ASOIL_IDX,IMODE )
     &             + ANTHFAC * AEROSPC_CONC( ACORS_IDX,IMODE )
            WETM3  = MOMENT3_CONC( IMODE )
            WETM2  = MOMENT2_CONC( IMODE )
            DRYM3  = WETM3 - H2OFAC * AEROSPC_CONC( AH2O_IDX, IMODE )  ! Assume no SOA in coarse mode
            DRYM20 = WETM2 * ( DRYM3 / WETM3 ) ** D_TWOTHIRDS

C *** Initial vapor-phase concentrations [ug/m3]
            GNO3R8 = PRECURSOR_CONC( HNO3_IDX )
            GNH3R8 = PRECURSOR_CONC( NH3_IDX )
            GCLR8  = PRECURSOR_CONC( HCL_IDX )
            
C *** Save the sensitivities and real values to new variable, before the first run of ISORROPIA            
            IF ( INT_TIME .EQ. 90.0D0 ) THEN
            	moment0_conc3_orig  = MOMENT0_CONC( IMODE )
            	moment2_conc3_orig  = MOMENT2_CONC( IMODE )
            	moment3_conc3_orig  = MOMENT3_CONC( IMODE )
            	aeromode_dens3_orig = aeromode_dens( IMODE )
            	aeromode_lnsg3_orig = aeromode_lnsg( IMODE )
            	aeromode_diam3_orig = aeromode_diam( IMODE )
            	aeromode_mass3_orig = aeromode_mass( IMODE )
            	anh4k_orig          = AEROSPC_CONC( ANH4_IDX,IMODE )
            	ano3k_orig          = AEROSPC_CONC( ANO3_IDX,IMODE )
            	aclk_orig           = AEROSPC_CONC( ACL_IDX,IMODE )
            	aso4k_orig          = AEROSPC_CONC( ASO4_IDX,IMODE )
            	ah2ok_orig          = AEROSPC_CONC( AH2O_IDX,IMODE )
            	nh3_orig            = PRECURSOR_CONC( NH3_IDX )
            	hno3_orig           = PRECURSOR_CONC( HNO3_IDX )
            	hcl_orig            = PRECURSOR_CONC( HCL_IDX )
            	hplus_orig          = HPLUS( IMODE )
            	seacat_orig         = AEROSPC_CONC( ASEACAT_IDX,IMODE ) 
            	soil_orig           = AEROSPC_CONC( ASOIL_IDX,IMODE )
            	cors_orig           = AEROSPC_CONC( ACORS_IDX,IMODE )
            ENDIF
            
C *** Before the last run of ISORROPIA, after the value "converges", 
C *** replace the calculated sensitivities of reverse ISORROPIA     

C *** Compute sulfate production rate [ug/m3 s] for coarse mode

            RATE = CONDSO4( IMODE )
            SO4  = CINORG( KSO4,IMODE ) + RATE * TSTEP * H2SO4RATM1
			
            IF ( TrustIso ) THEN

C *** Double Precision vars for ISORROPIA [mole/m3]
C              Compute equilibrium vapor pressures [mol/m3] of NH3, HNO3, and HCL
C              at the gas/particle interface of coarse mode aerosol.
C                 GAS(1) = NH3, GAS(2) = HNO3, GAS(3) = HCl
               WI( 1 ) = CINORG( KNA, IMODE ) * CFAC_ANA
               WI( 2 ) =                  SO4 * CFAC_ASO4
               WI( 3 ) = CINORG( KNH4,IMODE ) * CFAC_ANH4               
               WI( 4 ) = CINORG( KNO3,IMODE ) * CFAC_ANO3
               WI( 5 ) = CINORG( KCL, IMODE ) * CFAC_ACL
               WI( 6 ) = CINORG( KCA, IMODE ) * CFAC_ACA
               WI( 7 ) = CINORG( KK,  IMODE ) * CFAC_AK
               WI( 8 ) = CINORG( KMG, IMODE ) * CFAC_AMG
C              Also obtain the aqueous H+ concentration, AERLIQ(1) [mol/m3]

#ifdef verbose_aero
         if ( write1 ) then
         write( logdev,'(a, 8e13.3)' ) "VOLINORG,WI's C:",
     &                                  wi( 1 ), wi( 2 ), wi( 3 ), wi( 4 ),
     &                                  wi( 5 ), wi( 6 ), wi( 7 ), wi( 8 )
         write( logdev,* ) "VOLINORG,RH,T:", rhi, tempi
         end if
#endif

               CALL ISOROPIA( WI, RHI, TEMPI, CNTRL, WT, GAS, AERLIQ,  
     &                        AERSLD, SCASI, OTHER, COL, ROW, LAYER, FwdSens )

			   	 
               IF ( GAS( 1 ) .LT. 0.0D0 .OR. GAS( 2 ) .LT. 0.0D0 .OR.
     &              GAS( 3 ) .LT. 0.0D0 ) THEN
                  IF ( FIRSTWRITE ) THEN
                     FIRSTWRITE = .FALSE.
                     WRITE( LOGDEV,2023 )
                  END IF
                  WRITE( LOGDEV,2029 ) COL, ROW, LAYER, GAS( 1 ), GAS( 2 ), GAS( 3 )
                  TrustIso = .FALSE.
               END IF

            END IF   ! TrustIso
		 
	 ! If there is a change in H+ concentration through rev iso, ignore the sensitivities		 		 
 	 IF ( INT_TIME .GT. 90.0D0 ) THEN
 	 	
 	 	! Criteria 1: Significant Change of pH during the 4 steps of reverse ISORROPIA
 	 	! Used to be 1.25, AERLIQ(1) is the H+ concentration from the new reverse ISORROPIA run, and 
 	 	! HPLUS_PREV is the saved H+ concentration
 	 	IF  ( ( ( AERLIQ(1) / HPLUS_PREV ) .GE. 1.25d0 ) .OR. ( HPLUS_PREV / AERLIQ(1) ) .GE. 1.25d0 ) THEN
 	 		HybridSens = .false.
 	 	ENDIF
 	 	
 	 	! Criteria 2: Close to Neutral but not corrected by the SQRT(AKW)
 	 	IF ( AERLIQ( 1 ) /= 0.0D0 .AND. AERLIQ( 12 ) /= 0.0D0 ) THEN
 	 		IF (AERLIQ( 1 ) /= AERLIQ( 12 ) .AND. ( (AERLIQ( 1 ) / AERLIQ( 12 ) .LE. 1.D3 ) .OR. (AERLIQ( 12 ) / AERLIQ( 1 ) ) .LE. 1.D3 ) ) THEN
 	 			HybridSens = .false. 
 	 		ENDIF
 	 	ENDIF
 	 
 	 ELSE    ! INT_TIME .EQ. 90.0D0 for this case, only record down HPLUS here from first rev iso
 	 	HPLUS_PREV = AERLIQ( 1 )
 	 ENDIF
		 
		 
		 
C *** Change in volatile inorganic PM concentration to achieve
C     equilibrium, calculated as initial-gas-phase concentration minus 
C     equilibrium gas-phase concentration.  DVOLINORG is positive for
C     condensation and negative for evaporation.

#ifdef verbose_aero
         if ( write1 ) then
         write( logdev,'(a, 3e13.3)' ) "GASes:", gas( 1 ), gas( 2 ), gas( 3 )

         dvolinorg( knh4 ) = gnh3r8 * CFAC_GNH3  - gas( 1 )  ! [mol/m**3]
         dvolinorg( kno3 ) = gno3r8 * CFAC_GHNO3 - gas( 2 )  ! [mol/m**3]
         dvolinorg( kcl )  = gclr8  * CFAC_GHCL  - gas( 3 )  ! [mol/m**3]
           
         write( logdev,'(a, 3e13.3)' ) "DVOLINORG_coarse:",
     &      dvolinorg( knh4 ) / ( gnh3r8 * CFAC_GNH3 ),
     &      dvolinorg( kno3 ) / ( gno3r8 * CFAC_GHNO3),
     &      dvolinorg( kcl )  / ( gclr8  * CFAC_GHCL )
         end if
#endif

C *** Calculate condensation/evaporation flux for this time step and update 
C     volatile species concentrations.  Final aerosol conc set to be no less
C     than the minimum aerosol conc.
            IF ( TrustIso ) THEN
               CALL COMPUTE_FLUX( NVOLINORG, GNH3R8, GNO3R8, GCLR8, KNH4,
     &                            KNO3, KCL, GAS( 1:3 ), GRFAC2( IMODE ),
     &                            AERLIQ( 1 ), RATE, J, COL, ROW, LAYER )
            ELSE
               J( : ) = 0.0D0
            END IF 

            IF ( J( KNH4 ) * TSTEP * DF( KNH4 ) * NH3RAT .GT. GNH3R8 ) THEN
               WRITE( LOGDEV,* ) 'Condensed amt. exceeds NH3 conc: aero_subs.f'
               J( KNH4 ) = GNH3R8 / ( TSTEP * DF( KNH4 ) * NH3RAT )
            END IF
            IF ( J( KNO3 ) * TSTEP * DF( KNO3 ) * HNO3RAT .GT. GNO3R8 ) THEN
               WRITE( LOGDEV,* ) 'Condensed amt. exceeds HNO3 conc: aero_subs.f'
               J( KNO3 ) = GNO3R8 / ( TSTEP * DF( KNO3 ) * HNO3RAT )
            END IF
            IF ( J( KCL ) * TSTEP * DF(KCL) * HCLRAT .GT. GCLR8 ) THEN
               WRITE( LOGDEV,* ) 'Condensed amt. exceeds HCl conc: aero_subs.f'
               J( KCL ) = GCLR8 / ( TSTEP * DF( KCL ) * HCLRAT )
            END IF

C *** Integrate mass transfer equation, convert flux from molar to mass

            DO ISP = 1, NVOLINORG
               CFINAL( ISP,IMODE ) = MAX( 0.0D0,
     &                                    CINORG( ISP,IMODE )
     &                                    + J( ISP ) * TSTEP * DF( ISP ) )
            END DO 
            
C *** Calculate updated H+ concentration 

            HPLUS( IMODE ) = 0.0
     &                     - AEROSPC( ASO4_IDX )%CHARGE * SO4                  / AEROSPC_MW( ASO4_IDX )
     &                     - AEROSPC( ANO3_IDX )%CHARGE * CFINAL( KNO3,IMODE ) / AEROSPC_MW( ANO3_IDX )
     &                     - AEROSPC( ACL_IDX )%CHARGE  * CFINAL( KCL, IMODE ) / AEROSPC_MW( ACL_IDX )
     &                     - AEROSPC( ANH4_IDX )%CHARGE * CFINAL( KNH4,IMODE ) / AEROSPC_MW( ANH4_IDX )
!    &                     - AEROSPC( ANA_IDX )%CHARGE  * CINORG( KNA, IMODE ) / AEROSPC_MW( ANA_IDX )
     &                     - AEROSPC( ASEACAT_IDX )%CHARGE  * CINORG( KNA, IMODE ) / AEROSPC_MW( ASEACAT_IDX )
     &                     - AEROSPC( AMG_IDX )%CHARGE  * CINORG( KMG, IMODE ) / AEROSPC_MW( AMG_IDX )
     &                     - AEROSPC( AK_IDX )%CHARGE   * CINORG( KK,  IMODE ) / AEROSPC_MW( AK_IDX )
     &                     - AEROSPC( ACA_IDX )%CHARGE  * CINORG( KCA, IMODE ) / AEROSPC_MW( ACA_IDX )

C *** Equilibrate aerosol LWC with CFINAL by calling CALC_H2O
            WI( 1 ) = CINORG( KNA, IMODE ) * CFAC_ASEACAT
            WI( 2 ) =                  SO4 * CFAC_ASO4
            WI( 3 ) = CFINAL( KNH4,IMODE ) * CFAC_ANH4
            WI( 4 ) = CFINAL( KNO3,IMODE ) * CFAC_ANO3
            WI( 5 ) = CFINAL( KCL, IMODE ) * CFAC_ACL
            WI( 6 ) = CINORG( KCA, IMODE ) * CFAC_ACA
            WI( 7 ) = CINORG( KK,  IMODE ) * CFAC_AK
            WI( 8 ) = CINORG( KMG, IMODE ) * CFAC_AMG

            CALL CALC_H2O( WI, RHI, TEMPI, H2O_NEW ) 

            H2O = H2O_NEW * DFH2OR8 

C *** Update all Local Aerosol Mass and Vapor Concentrations 
            !Aerosol
            AEROSPC_CONC( ANH4_IDX,IMODE ) = CFINAL( KNH4,IMODE )
            AEROSPC_CONC( ANO3_IDX,IMODE ) = CFINAL( KNO3,IMODE )
            AEROSPC_CONC( ACL_IDX,IMODE )  = CFINAL( KCL, IMODE )
            AEROSPC_CONC( ASO4_IDX,IMODE ) = SO4
            AEROSPC_CONC( AH2O_IDX,IMODE ) = H2O

            !Gas
            PRECURSOR_CONC( NH3_IDX ) = GNH3R8 + ( CINORG( KNH4,IMODE )
     &                                 -CFINAL( KNH4,IMODE ) ) * NH3RAT 
            PRECURSOR_CONC( HNO3_IDX )= GNO3R8 + ( CINORG( KNO3,IMODE )
     &                                 -CFINAL( KNO3,IMODE) ) * HNO3RAT  
            PRECURSOR_CONC( HCL_IDX ) = GCLR8 + ( CINORG( KCL,IMODE )
     &                                 -CFINAL( KCL,IMODE) ) * HCLRAT 
     		 			
C *** Compute net change in 3rd moment due to dry inorganic mass transfer

            DDRYM3DT = ( ( CFINAL( KNH4,IMODE ) - CINORG( KNH4,IMODE ) ) * M3FAC_ANH4
     &                 + ( CFINAL( KNO3,IMODE ) - CINORG( KNO3,IMODE ) ) * M3FAC_ANO3
     &                 + ( CFINAL( KCL, IMODE ) - CINORG( KCL,IMODE ) )  * M3FAC_ACL
     &                 + ( SO4                  - CINORG( KSO4,IMODE ) ) * M3FAC_ASO4 ) / TSTEP

C *** Compute net change in 2nd moment due to dry inorganic mass transfer
C     (including nucleation) using equation A7 of Binkowski & Shankar (1995)
            DDRYM2DT = D_TWOTHIRDS * GRFAC1( IMODE ) / GRFAC2( IMODE ) * DDRYM3DT   

C *** Update dry 2nd moment for condensation/evaporation based on whether
C     net change in dry 2nd moment is production or loss
            IF ( DDRYM2DT .LT. 0.0D0 ) THEN
               LOSS = DDRYM2DT / DRYM20
               Y = DRYM20 * EXP( LOSS * TSTEP )
            ELSE
               Y = DRYM20 + DDRYM2DT * TSTEP
            END IF

C *** Add water (no SOA) 2nd moment while preserving standard deviation

            !Calculate 3rd Moment
            DRYM3 = ( M3FAC_ASO4 ) * SO4
     &            + M3FAC_ANH4 * CFINAL( KNH4,IMODE )
     &            + M3FAC_ANO3 * CFINAL( KNO3,IMODE )
     &            + M3FAC_ACL  * CFINAL( KCL,IMODE )
     &            + M3FAC_ANA  * SEACAT
     &            + M3OTHR                   
            WETM3 = DRYM3 + H2O * M3FAC_H2O

            !Calculate 2nd moment
            DRYM2 = MAX( REAL( AEROMODE_MINM2( IMODE ), 8 ), Y )
            WETM2 = DRYM2 * ( WETM3 / DRYM3 ) ** D_TWOTHIRDS

            MOMENT2_CONC( IMODE ) = WETM2

C *** Update the third moments, geometric mean diameters, geometric 
C     standard deviations, modal mass totals, and modal particle 
C     densities. It is a waste of time updating the aitken and
C     accumulation modes but the coarse mode does need to be updated
C     each sub-time step. This should be for the wet distribution
               
            CALL GETPAR( .TRUE. )
         END DO   ! end mass transfer time step loop
         
      END IF   ! for 'Hybrid' method
	  
	  ! If there are drastic changes in pH during rev iso, the sensitivities have asymptotic behaviors
	  ! See Henningan et. al. 2015 ACP paper for details
	
C  moment0_conc3_orig = 0.0d0
C  moment2_conc3_orig = 0.0d0
C  moment3_conc3_orig = 0.0d0
C  aeromode_dens3_orig = 0.0d0
C  aeromode_lnsg3_orig = 0.0d0
C  aeromode_diam3_orig = 0.0d0
C  aeromode_mass3_orig = 0.0d0
C  anh4k_orig = 0.0d0
C  ano3k_orig = 0.0d0
C  aso4k_orig = 0.0d0
C  nh3_orig = 0.0d0
C  hno3_orig = 0.0d0
C  hcl_orig = 0.0d0
C  hplus_orig = 0.0d0
C  seacat_orig = 0.0d0
C  soil_orig = 0.0d0
C  cors_orig = 0.0d0
	

      IF ( HYBRID .AND. ( HybridSens .eq. .false.) ) THEN
   	   ! Reset all moment related stuff
   		call reset_hyd_sens( moment0_conc( 3 ), moment0_conc3_orig)
   		call reset_hyd_sens( moment2_conc( 3 ), moment2_conc3_orig)
   		call reset_hyd_sens( moment3_conc( 3 ), moment3_conc3_orig)
   		call reset_hyd_sens( aeromode_dens( 3 ), aeromode_dens3_orig)
   		call reset_hyd_sens( aeromode_lnsg( 3 ), aeromode_lnsg3_orig)
   		call reset_hyd_sens( aeromode_diam( 3 ), aeromode_diam3_orig)
   		call reset_hyd_sens( aeromode_mass( 3 ), aeromode_mass3_orig)
   		call reset_hyd_sens( AEROSPC_CONC( ANH4_IDX, 3 ), anh4k_orig)
   		call reset_hyd_sens( AEROSPC_CONC( ANO3_IDX, 3 ), ano3k_orig)
   		call reset_hyd_sens( AEROSPC_CONC( ACL_IDX, 3 ), aclk_orig)
   		call reset_hyd_sens( AEROSPC_CONC( ASO4_IDX, 3 ), aso4k_orig)
   		call reset_hyd_sens( AEROSPC_CONC( AH2O_IDX, 3), ah2ok_orig)
   		call reset_hyd_sens( PRECURSOR_CONC( NH3_IDX ), nh3_orig)	
   		call reset_hyd_sens( PRECURSOR_CONC( HNO3_IDX ), hno3_orig)
   		call reset_hyd_sens( PRECURSOR_CONC( HCL_IDX ), hcl_orig)  	
      ENDIF
      write1 = .false.
C *** End of Coarse Mode dynamic mass transfer calculations

C *** Fine Aerosol Modes: Call ISORROPIA in forward mode to calculate gas-particle equilibrium
      !Get Precursors Vapor concentrations [ug m-3]
      GNH3R8 = PRECURSOR_CONC( NH3_IDX )   
      GNO3R8 = PRECURSOR_CONC( HNO3_IDX )  
      GCLR8  = PRECURSOR_CONC( HCL_IDX )

C *** Diagnose all total gas+particle concentrations to passed to
C     ISORROPIA. Convert everything to [mol m-3].
      WI( 1 ) = SUM( AEROSPC_CONC( ANA_IDX,1:2 ))   * CFAC_ANA
      !Compute sulfate from total sulfate production rate [ug/m3-s] for fine 
      !modes; add in H2SO4 nucleated in model timestep
      WI( 2 ) = ( SUM( AEROSPC_CONC( ASO4_IDX,1:2 ) ) 
     &           +( DMDT_SO4 + SUM( CONDSO4( 1:2 ) ) ) * DELT * H2SO4RATM1 ) * CFAC_ASO4
      WI( 3 ) = PRECURSOR_CONC( NH3_IDX )  * CFAC_GNH3
     &         +SUM( AEROSPC_CONC( ANH4_IDX,1:2 ) ) * CFAC_ANH4 
      WI( 4 ) = PRECURSOR_CONC( HNO3_IDX ) * CFAC_GHNO3
     &         +SUM( AEROSPC_CONC( ANO3_IDX,1:2 ) ) * CFAC_ANO3 
      WI( 5 ) = PRECURSOR_CONC( HCL_IDX )  * CFAC_GHCL 
     &         +SUM( AEROSPC_CONC( ACL_IDX,1:2 ) )  * CFAC_ACL 
      WI( 6 ) = SUM( AEROSPC_CONC( ACA_IDX,1:2 ) )  * CFAC_ACA 
      WI( 7 ) = SUM( AEROSPC_CONC( AK_IDX,1:2 ) )   * CFAC_AK 
      WI( 8 ) = SUM( AEROSPC_CONC( AMG_IDX,1:2 ) )  * CFAC_AMG 

      CNTRL( 1 ) = 0.0D0   ! Forward Problem
      CNTRL( 2 ) = 1.0D0   ! Aerosol in Metastable State

C *** Set flags to account for mass conservation violations in ISRP3F
      TRUSTCL  = .TRUE.
      IF ( (WI( 1 ) + WI( 5 )) .LT. 1.0D-20 .or. WI( 5 ) .LT. 1.0D-10 ) THEN
         TRUSTCL = .FALSE.
      END IF
         
#ifndef noiso
      
      IF ( FWDISO ) THEN 
           
      CALL ISOROPIA( WI, RHI, TEMPI, CNTRL, WT, GAS, AERLIQ,
     &               AERSLD, SCASI, OTHER, COL, ROW, LAYER, FwdSens)
     

C *** Save H+ concentration information in microgram/m3 for consistency
      EQLBHIJ = AERLIQ(1) * 1.0D6 * AEROSPC_MW( ah3op_idx)
      
      ELSE 
      
      gas( 1 ) = precursor_conc( nh3_idx ) * CFAC_GNH3
      gas( 2 ) = precursor_conc( hno3_idx) * CFAC_GHNO3
      gas( 3 ) = precursor_conc( hcl_idx ) * CFAC_GHCL

      EQLBHIJ =  HPLUS( 1 ) + HPLUS( 2 )  ! use charge balance if lacking isoropia info 
      
      ENDIF

#else
      gas( 1 ) = precursor_conc( nh3_idx ) * CFAC_GNH3
      gas( 2 ) = precursor_conc( hno3_idx) * CFAC_GHNO3
      gas( 3 ) = precursor_conc( hcl_idx ) * CFAC_GHCL

      EQLBHIJ =  HPLUS( 1 ) + HPLUS( 2 )  ! use charge balance if lacking isoropia info 

#endif

C *** Change in volatile inorganic PM concentration to achieve
C     equilibrium, calculated as initial-gas-phase concentration minus 
C     equilibrium gas-phase concentration.  DVOLINORG is positive for
C     condensation and negative for evaporation.
	
	  ! If the real concentration of NH3 fluctuates too much in forward ISORROPIA 
C  FwdSens = .TRUE. 
C  IF ( FwdIso .AND. (FwdSens .eq. .FALSE.) ) THEN
C  	call reset_hyd_sens( gas( 1 ), precursor_conc( nh3_idx ) * CFAC_GNH3)
C  	call reset_hyd_sens( gas( 2 ), precursor_conc( hno3_idx ) * CFAC_GHNO3)
C  	call reset_hyd_sens( gas( 3 ), precursor_conc( hcl_idx ) * CFAC_GHCL)
C  ENDIF	  	
	  
	  
      DVOLINORG( KNH4 ) = GNH3R8 * CFAC_GNH3  - GAS( 1 )   ! mol m-3
      DVOLINORG( KNO3 ) = GNO3R8 * CFAC_GHNO3 - GAS( 2 )   ! mol m-3
      DVOLINORG( KCL )  = GCLR8  * CFAC_GHCL  - GAS( 3 )   ! mol m-3
      

      IF ( DVOLINORG( KNH4 ) .LT. 0.0D0 ) THEN
         DVOLMAX = -SUM(AEROSPC_CONC( ANH4_IDX,1:2 ) ) * CFAC_ANH4 + EVAPMIND
         DVOLINORG( KNH4 ) = MAX( DVOLINORG( KNH4 ), DVOLMAX )
      END IF

      IF ( DVOLINORG( KNO3 ) .LT. 0.0D0 ) THEN
         DVOLMAX = -SUM(AEROSPC_CONC( ANO3_IDX,1:2 ) ) * CFAC_ANO3 + EVAPMIND
         DVOLINORG( KNO3 ) = MAX( DVOLINORG( KNO3 ), DVOLMAX)
      END IF

      IF ( .not.TRUSTCL ) THEN  
         DVOLINORG( KCL ) = 0.0D0
      ELSEIF ( DVOLINORG( KCL ) .LT. 0.0D0 ) THEN
         DVOLMAX = -SUM( AEROSPC_CONC( ACL_IDX,1:2 ) ) * CFAC_ACL + EVAPMIND
         DVOLINORG( KCL ) = MAX( DVOLINORG( KCL ), DVOLMAX )
      END IF

C *** Apply modal partitioning of equilibrium aerosol mass
      ! Calculate Distribution of Mass Transfer Among Modes
      OMEGA( 1 ) = GRFAC2( 1 ) / ( GRFAC2( 1 ) + GRFAC2( 2 ) )
      OMEGA( 2 ) = 1.0D0 - OMEGA( 1 )

      ! Save Initial Concentrations
      DO IMODE = 1, 2  
         CINORG( KSO4,IMODE ) = AEROSPC_CONC( ASO4_IDX, IMODE )
         CINORG( KNH4,IMODE ) = AEROSPC_CONC( ANH4_IDX, IMODE )
         CINORG( KNO3,IMODE ) = AEROSPC_CONC( ANO3_IDX, IMODE )
         CINORG( KNA, IMODE ) = AEROSPC_CONC( ANA_IDX,  IMODE )
         CINORG( KCL, IMODE ) = AEROSPC_CONC( ACL_IDX,  IMODE )
         CINORG( KCA, IMODE ) = AEROSPC_CONC( ACA_IDX,  IMODE )
         CINORG( KK,  IMODE ) = AEROSPC_CONC( AK_IDX,   IMODE )
         CINORG( KMG, IMODE ) = AEROSPC_CONC( AMG_IDX,  IMODE )
      END DO

      ! Calculate Initial Distribution of Mass Composition Among Modes
      DO ISP = 1, NVOLINORG
         TOTAER( ISP ) = MAX( SUM( CINORG( ISP,1:2 ) ), CONMIND )
         DO IMODE = 1, 2
            PHI( ISP, IMODE ) = CINORG( ISP,IMODE ) / TOTAER( ISP )
         ENDDO
      ENDDO

      ! Initialize Final Concentrations
      CFINAL = 0.0D0

      ! Calculate Final Concentrations
      DO ISP = 1, NVOLINORG
         IF ( DVOLINORG( ISP ) .LT. 0.0 ) THEN
            ! Evaporate Mass Using Condensed-Phase Fraction in each Mode
            CFINAL( ISP,1:2 ) = CFINAL( ISP,1:2 ) + CINORG( ISP,1:2 )
     &                        + PHI( ISP,1:2 ) * DVOLINORG( ISP ) * DF( ISP )
         ELSE
            ! Condense Mass Using Condensation Sink Factors
            CFINAL( ISP,1:2 ) = CFINAL( ISP,1:2 ) + CINORG( ISP,1:2 )
     &                        + OMEGA( 1:2 ) * DVOLINORG( ISP ) * DF( ISP )
         END IF
      END DO

C *** Apply Final Concentrations to Moment Variables in Aerosol Scheme
      ! Store Initial Wet Moments
      M3WET_INIT( : ) = MOMENT3_CONC( : )
      M2WET_INIT( : ) = MOMENT2_CONC( : )

      ! Calculate and Store Initial Dry Moments and set wet_moments_flag to false/dry
      ! There has not yet been an update to AEROSPC_CONC so this will
      ! calculate the old 3rd moment.
      call calcmoments( .false. )
      M3DRY_INIT( : ) = MOMENT3_CONC( : )
      M2DRY_INIT( : ) = MOMENT2_CONC( : )

      DO IMODE = 1, 2  ! modal partitioning of equilibrium aerosol mass

         IF ( IMODE .EQ. 1 ) THEN
            ! Update Number Concentration with NPF
            MOMENT0_CONC( IMODE ) =  AM0(IMODE) + DNDT * DELT 
            ! Add NPF to total condensation rate
            RATE = DMDT_SO4 + CONDSO4( IMODE )
            SO4 = CINORG( KSO4,IMODE ) + RATE * DELT * H2SO4RATM1 
         ELSE
            ! Ignore NPF. The small particles are not in these modes
            SO4 = CINORG( KSO4,IMODE ) + CONDSO4( IMODE ) * DELT * H2SO4RATM1
         END IF

C *** Double precision vars for CALC_H2O
         WI( 1 ) = CINORG( KNA, IMODE ) * CFAC_ANA
         WI( 2 ) =                  SO4 * CFAC_ASO4
         WI( 3 ) = CFINAL( KNH4,IMODE ) * CFAC_ANH4
         WI( 4 ) = CFINAL( KNO3,IMODE ) * CFAC_ANO3
         WI( 5 ) = CFINAL( KCL, IMODE ) * CFAC_ACL
         WI( 6 ) = CINORG( KCA, IMODE ) * CFAC_ACA
         WI( 7 ) = CINORG( KK,  IMODE ) * CFAC_AK
         WI( 8 ) = CINORG( KMG, IMODE ) * CFAC_AMG

         CALL CALC_H2O( WI, RHI, TEMPI, H2O_NEW ) 
         
         ! Update All Aerosol Concentrations
         AEROSPC_CONC( AH2O_IDX, IMODE ) = H2O_NEW * DFH2OR8 ! 4
         AEROSPC_CONC( ANH4_IDX, IMODE ) = CFINAL( KNH4,IMODE ) ! 4
         AEROSPC_CONC( ANO3_IDX, IMODE ) = CFINAL( KNO3,IMODE ) ! 4
         AEROSPC_CONC( ACL_IDX, IMODE )  = CFINAL( KCL ,IMODE ) ! 4
         AEROSPC_CONC( ASO4_IDX, IMODE ) = SO4                  ! 4

C *** Compute net change in 3rd moment due to dry inorganic mass transfer
C     (includes nucleated sulfate mass). This is for projecting the
C     change to the second moment due to dry inorganic condensation
         
         DDRYM3DT = ( ( CFINAL( KNH4,IMODE ) - CINORG( KNH4,IMODE ) ) * M3FAC_ANH4
     &              + ( CFINAL( KNO3,IMODE ) - CINORG( KNO3,IMODE ) ) * M3FAC_ANO3
     &              + ( CFINAL( KCL, IMODE ) - CINORG( KCL,IMODE ) )  * M3FAC_ACL
     &              + ( SO4                  - CINORG( KSO4,IMODE ) ) * M3FAC_ASO4 ) 
     &              / DELT 

C *** Compute net change in 2nd moment due to dry inorganic mass transfer
C     (including nucleation) using equation A7 of Binkowski & Shankar (1995)
         DDRYM2DT = D_TWOTHIRDS * GRFAC1( IMODE ) / GRFAC2( IMODE ) * DDRYM3DT

C *** Update dry 2nd moment for condensation/evaporation based on whether
C     net change in dry 2nd moment is production or loss
         IF ( DDRYM2DT .LT. 0.0D0 ) THEN
            LOSS = DDRYM2DT / M2DRY_INIT( IMODE )
            Y = M2DRY_INIT( IMODE ) * EXP( LOSS * DELT )
         ELSE
            Y = M2DRY_INIT( IMODE ) + DDRYM2DT * DELT
         END IF
         moment2_conc( IMODE ) = MAX( REAL(AEROMODE_MINM2( IMODE ), 8), Y )

      END DO

C *** Add water and SOA to 2nd moment while preserving standard deviation
      call calcmoments( .true. )
      M3WET_FINAL = moment3_conc
      M2WET_FINAL = moment2_conc

C *** Assign H+ Concentration to each Mode
      HPLUS( 1:2 ) = 0.0D0
      DO I = 1, N_AEROSPC
         IF(  AEROSPC( I )%tracer )CYCLE
         HPLUS( 1:2 ) = HPLUS( 1:2 )
     &                  - AEROSPC( I )%CHARGE * AEROSPC_CONC( I,1:2 ) / AEROSPC_MW( I )
      END DO

      H2O = AEROSPC_CONC( AH2O_IDX, 1 ) + AEROSPC_CONC( AH2O_IDX, 2 )
      IF(  H2O .GT. CONMIN )THEN
          H2O = 1.0 / H2O
          AEROSPC_CONC( AH3OP_IDX, 1:2 ) = EQLBHIJ * H2O * AEROSPC_CONC( AH2O_IDX, 1:2 )  ! 4
      ELSE
          AEROSPC_CONC( AH3OP_IDX, 1:2 ) = CONMIND
      END IF
      AEROSPC_CONC( AH3OP_IDX, 3 ) = HPLUS( 3 ) * AEROSPC_MW( ah3op_idx ) ! Coarse mode H+ concentration in ug/m3, 4

C *** Update the third moments, geometric mean diameters, geometric 
C     standard deviations, modal mass totals, and modal particle 
C     densities. Note that moment2_conc needs to be up to date when this
C     routine is called. Moment3_conc does not need to be up to date
C     because it will be recalculated inside GETPAR as the sum of 
C     aerospc_conc variables. This should be for the wet distribution.
      CALL GETPAR( .TRUE. )
       
! *** Propagate Concentration Changes from Coarse-Mode Mass
! *** Transfer to Budget Vectors
      DELNUM( : )  = MOMENT0_CONC( : ) - CBUDGET0_NUM( : )
      DELSRF( : )  = MOMENT2_CONC( : ) - CBUDGET0_SRF( : )
      DELSO4( : )  = AEROSPC_CONC( ASO4_IDX, : ) - CBUDGET0_SO4( : )
      DELNH4( : )  = AEROSPC_CONC( ANH4_IDX, : ) - CBUDGET0_NH4( : )
      DELNO3( : )  = AEROSPC_CONC( ANO3_IDX, : ) - CBUDGET0_NO3( : )
      DELCL( : )   = AEROSPC_CONC( ACL_IDX , : ) - CBUDGET0_CL( : )
      DELH2O( : )  = AEROSPC_CONC( AH2O_IDX, : ) - CBUDGET0_H2O( : )
      DELH3OP( : ) = AEROSPC_CONC( AH3OP_IDX,: ) - CBUDGET0_H3OP( : )

      ! precision = 4
      DO MODE = 1,N_MODE
      	COND_BUDGET(AERONUM_MAP(MODE)) = DELNUM(MODE)
      	COND_BUDGET(AEROSRF_MAP(MODE)) = DELSRF(MODE)
      	COND_BUDGET(AEROSPC_MAP(ASO4_IDX,MODE)) = DELSO4(MODE)
      	COND_BUDGET( AEROSPC_MAP( ANH4_IDX,MODE ) ) = DELNH4( MODE )
      	COND_BUDGET( AEROSPC_MAP( ANO3_IDX,MODE ) ) = DELNO3( MODE )
     	COND_BUDGET( AEROSPC_MAP( ACL_IDX ,MODE ) ) = DELCL( MODE )
      	COND_BUDGET( AEROSPC_MAP( AH2O_IDX,MODE ) ) = DELH2O( MODE )
      	COND_BUDGET( AEROSPC_MAP( AH3OP_IDX,MODE) ) = DELH3OP( MODE )
      ENDDO
       
     
C     COND_BUDGET( AERONUM_MAP( : ) ) =  DELNUM( : )
C     COND_BUDGET( AEROSRF_MAP( : ) ) =  DELSRF( : )
C     COND_BUDGET( AEROSPC_MAP( ASO4_IDX,: ) ) = DELSO4( : )
C     COND_BUDGET( AEROSPC_MAP( ANH4_IDX,: ) ) = DELNH4( : )
C     COND_BUDGET( AEROSPC_MAP( ANO3_IDX,: ) ) = DELNO3( : )
C     COND_BUDGET( AEROSPC_MAP( ACL_IDX ,: ) ) = DELCL( : )
C     COND_BUDGET( AEROSPC_MAP( AH2O_IDX,: ) ) = DELH2O( : )
C     COND_BUDGET( AEROSPC_MAP( AH3OP_IDX,:) ) = DELH3OP( : )

      COND_BUDGET( PRECURSOR_MAP( SULF_IDX ) ) = SUM( -DELSO4(:) * H2SO4RAT)
      COND_BUDGET( PRECURSOR_MAP( NH3_IDX  ) ) = SUM( -DELNH4(:) * NH3RAT )
      COND_BUDGET( PRECURSOR_MAP( HNO3_IDX ) ) = SUM( -DELNO3(:) * HNO3RAT )
      COND_BUDGET( PRECURSOR_MAP( HCL_IDX  ) ) = SUM( -DELCL(:)  * HCLRAT  )

 
C *** Update gas-phase concentrations from Aitken and Accumulation Mode Partitioning
      PRECURSOR_CONC( NH3_IDX )  = GNH3R8 + sum( CINORG( KNH4,1:2 )
     &                            -CFINAL( KNH4,1:2 ) ) * NH3RAT 
      PRECURSOR_CONC( HNO3_IDX ) = GNO3R8 + sum( CINORG( KNO3,1:2 )
     &                            -CFINAL( KNO3,1:2) ) * HNO3RAT  
      PRECURSOR_CONC( HCL_IDX )  = GCLR8 + sum( CINORG( KCL,1:2 )
     &                            -CFINAL( KCL,1:2) ) * HCLRAT 
     
         
2023  FORMAT( 1X, 'VOLINORG returning negative gas concentrations from ISOROPIA:'
     &       /10X, 'GAS(1) = NH3, GAS(2) = HNO3, GAS(3) = HCl' )
2029  FORMAT( 1X, '[see VOLINORG msg]'
     &        1X, 'at (C,R,L): ', 3I4, 1X, 'GAS Conc:', 3( 1PE11.3 ) )

      RETURN
      END SUBROUTINE VOLINORG

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE HCOND3( AM0, AM1, AM2, DV, ALPHA, CBAR, F )

C  calculates the size-dependent term in the condensational-growth rate
C  expression for the 2nd and 3rd moments of a lognormal aerosol mode
C  using the harmonic mean method.  This code follows Section A2 of
C  Binkowski & Shankar (1995).
 
C  Key Subroutines/Functions called:  none
 
C  Revision History:
C     coded November 7, 2003 by Dr. Francis S. Binkowski
C     Revised November 20, 2003 by F. Binkowski to have am1 and
C     am2 as inputs
 
C  Reference:
C   1. Binkowski, F.S. and U. Shankar, The regional particulate matter
C      model 1. Model description and preliminary results, J. Geophys.
C      Res., Vol 100, No D12, 26101-26209, 1995.

	  use HDMod 
      IMPLICIT NONE

C *** Includes:

      INCLUDE SUBST_CONST     ! physical and mathematical constants

C *** Arguments:

      TYPE(hyperdual), INTENT( IN ) :: AM0   ! zeroth moment of mode  [ #/m**3 ]
      TYPE(hyperdual), INTENT( IN ) :: AM1   ! first moment of mode   [ m/m**3 ]
      TYPE(hyperdual), INTENT( IN ) :: AM2   ! second moment of mode  [ m**2/m**3 ]
      REAL,            INTENT( IN ) :: Dv    ! molecular diffusivity of the
                                       ! condensing vapor  [ m**2/s ]
      REAL,      INTENT( IN ) :: ALPHA ! accommodation coefficient
      REAL,      INTENT( IN ) :: CBAR  ! kinetic velocity of condensing vapor [ m/s ]

      TYPE(hyperdual), INTENT( OUT ) :: F( 2 ) ! size-dependent term in condensational-growth
                                         ! rate: F(1) = 2nd moment [ m**2/m**3 s ]
                                         !       F(2) = 3rd moment [ m**3/m**3 s ]

C *** Local Variables:

      TYPE(hyperdual) :: GNC2 ! integrals used to calculate F(1) [m^2 / m^3 s]
      TYPE(hyperdual) :: GFM2 !

      TYPE(hyperdual) :: GNC3 ! integrals used to calculate F(2) [m^3 / m^3 s]
      TYPE(hyperdual) :: GFM3 !

      REAL( 8 ), PARAMETER :: TWOPI = 2.0D0 * PI
      REAL( 8 ), PARAMETER :: PI4 = 0.25D0 * PI

C-----------------------------------------------------------------------

C *** Implement equation A15 of Binkowski & Shankar (1995) for the
C     2nd and 3rd moments of a lognormal mode of arbitrary size.

      GNC2 = TWOPI * DV * AM0          ! 2nd moment, near-continuum
      GNC3 = TWOPI * DV * AM1          ! 3rd moment, near-continuum
      GFM2 = PI4 * ALPHA * CBAR * AM1  ! 2nd moment, free-molecular
      GFM3 = PI4 * ALPHA * CBAR * AM2  ! 3rd moment, free-molecular

C *** Implement equation A13 of Binkowski & Shankar (1995) for a
C     lognormal mode of arbitrary size.  These are the size-dependent
C     terms in the condensational-growth rate expression, given in
C     equation 7a of B&S (1995).

      F( 1 ) = GNC2 * GFM2 / ( GNC2 + GFM2 )  ! 2nd moment
      F( 2 ) = GNC3 * GFM3 / ( GNC3 + GFM3 )  ! 3rd moment

      RETURN
      END SUBROUTINE HCOND3

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE NEWPART3 ( RH, T, XH2SO4, SO4RATE, DNDT, DMDT_so4, DM2DT )

      USE AERO_DATA
      USE AEROMET_DATA   ! Includes CONST.EXT
      USE PRECURSOR_DATA, ONLY: PRECURSOR_MW, SULF_IDX
      USE UTILIO_DEFN
      
	  use HDMod 
      IMPLICIT NONE

C   REVISION HISTORY:
C     Replacement of Kulmala et al., 1998 nucleation scheme with scheme 
C     of Vehkamaki et al. (2002) by G. Sarwar and K. Fahey - 03/2014

C.. References:
C     Vehkamaki, H., Kulmala, M, Napari, I., Lehtinen, K.E.J., Timmreck, C.,
C     Noppel, M., and A. Laaksonen. (2002) An improved parameterization for 
C     sulfuric acid-water nucleation rates for tropospheric and stratospheric
C     conditions.  JGR, v107(D22).    
   
C..Inputs: 
C     The RH and T here are AIRRH and AIRTEMP from AEROMET_DATA
C     They are real numbers and not the same one in ISORROPIA
c     These should be single precision reals

      REAL, INTENT(IN ) :: RH              ! fractional relative humidity      
      REAL, INTENT(IN)  :: T               ! Ambient temperature [ K ]
      TYPE(hyperdual), INTENT(IN)  :: XH2SO4          ! sulfuric acid concentration [ ug/m**3 ]
      TYPE(hyperdual), INTENT(IN)  :: SO4RATE         ! gas-phase H2SO4 production rate [ ug/m**3 s ]

C.. Outputs:
      TYPE(hyperdual), INTENT( OUT ) :: DNDT     ! particle number production rate [ m^-3/s ]
      TYPE(hyperdual), INTENT( OUT ) :: DMDT_so4 ! SO4 mass production rate [ ug/m**3 s ]
      TYPE(hyperdual), INTENT( OUT ) :: DM2DT    ! second moment production rate [ m**2/m**3 s ]

C.. Parameters
      CHARACTER( 16 ), PARAMETER :: PNAME = 'NEWPART'                 

C.. Particle size parameters:
      REAL, PARAMETER :: d20 = 2.0E-07                ! diameter of a new particle [cm]
      REAL, PARAMETER :: d20sq = d20 * d20            ! new-particle diameter squared [cm**2]
      REAL, PARAMETER :: m2_20 = 1.0E-4 * d20sq       ! new-particle diameter squared [m**2]
      REAL, PARAMETER :: v20 = PI * d20 * d20sq /6.0  ! volume of a new particle [cm**3]
      
      TYPE(hyperdual)       :: sulfmass                     ! mass of a new particle [ug]
      TYPE(hyperdual)       :: sulfmass1                    ! inverse of sulfmass [ug**-1]

C.. Set constants and local variables for Vehkamaki et al. (2002) scheme
      REAL, PARAMETER :: C1 = 0.740997
      REAL, PARAMETER :: C2 = -0.00266379
      REAL, PARAMETER :: C3 = -0.00349998
      REAL, PARAMETER :: C4 = 0.0000504022
      REAL, PARAMETER :: C5 = 0.00201048
      REAL, PARAMETER :: C6 = -0.000183289
      REAL, PARAMETER :: C7 = 0.00157407
      REAL, PARAMETER :: C8 = -0.0000179059
      REAL, PARAMETER :: C9 = 0.000184403
      REAL, PARAMETER :: C10 = -1.50345E-6

      TYPE(hyperdual) :: XSTAR           ! mole fraction of sulfuric acid in the critical cluster
      TYPE(hyperdual) :: NA              ! total gas phase concentration of H2SO4  [ #/cm**3 ]
      REAL :: TEMP            ! ambient temperature
      REAL            :: LNRH, LNRH2, LNRH3
      TYPE(hyperdual) :: LNNA, LNNA2, LNNA3
      REAL :: TEMP2, TEMP3    ! TEMP**2, TEMP**3
      TYPE(hyperdual) :: XFAC       ! exponential term for the nucleation rate  
      TYPE(hyperdual) :: Jnuc       ! nucleation rate [ #/cm**3 s ]
      
      TYPE(hyperdual) :: A, B, C, D, E, F, G, H, I, J

      REAL :: MW_H2SO4        ! MW of H2SO4 in [ g / mole ]        
      REAL :: DENSITY_H2SO4   ! DENSITY of H2SO4 in [ kg / m**3 ]

      REAL, PARAMETER :: SCALEFAC = 1.0E-06           ! for [ 1 / m**3 ] ->  [ 1 / cm**3 ]
      REAL, PARAMETER :: MUG2G = 1.0E-6               ! [ ug  ] -> [ g ]                 

C.. Initialize variables
      DNDT     = 0.0D0
      DMDT_so4 = 0.0D0
      DM2DT    = 0.0D0

C.. Calculate molecular weight of H2SO4 [ g / mole ] 
      MW_H2SO4 = PRECURSOR_MW( SULF_IDX )
         
C.. Calculate density of sulfuric acid [ kg / m**3 ] 
      DENSITY_H2SO4 = AEROSPC( ASO4_IDX )%DENSITY           

C.. Calculate mass of a new particle [ug]
      sulfmass = 1.0D+3 * DENSITY_H2SO4 * v20

C.. Calculate inverse of sulfmass [ug**-1]      
      sulfmass1 = 1.0D0 / sulfmass   
     
C.. Calculate sulfuric acid concentration in molecules/cm3 
      NA = XH2SO4 * MUG2G * AVO * SCALEFAC / MW_H2SO4     

C.. The parameterization is valid at sulfuric acid concentrations of 1.0E4 - 1.0E11 molecules cm-3
      NA  = MAX (NA, 1.0D4)
      NA  = MIN (1.0D11, NA)

C.. The parameterization is valid at temperatures of 190.00-305.15 K
      TEMP = MAX (T, 190.00)
      TEMP = MIN (305.15, TEMP)

C.. The parameterization is valid at RH of 0.0001-1.0
C   aero_driver.f limits RH to 0.005-0.99; thus no additional constraint is needed

C.. Define convenient constants
      TEMP2 = TEMP * TEMP
      TEMP3 = TEMP * TEMP2
      
      LNRH = LOG( RH )  
      LNNA = LOG( NA )
      LNRH2 = LNRH * LNRH
      LNRH3 = LNRH * LNRH2
      LNNA2 = LNNA * LNNA
      LNNA3 = LNNA * LNNA2

C.. Calculate mole fraction of sulfuric acid in the critical cluster
      XSTAR = C1 + C2 * TEMP + C3 * LNNA + C4 * TEMP * LNNA + 
     &      C5 * LNRH + C6 * TEMP * LNRH +
     &      C7 * LNRH2 + C8 * TEMP * LNRH2 + C9 * LNRH3 + 
     &      C10 * TEMP * LNRH3

C.. Calculate coefficients needed for the nucleation rate [Eq-12 - Vehkamaki et al., 2002]           
      A = 0.14309 + 2.21956 * TEMP - 0.0273911 * TEMP2 + 
     &      0.0000722811 * TEMP3 + 5.91822 / XSTAR
     
      B = 0.117489 + 0.462532 * TEMP - 0.0118059 * TEMP2 + 
     &      0.0000404196 * TEMP3 + 15.7963 / XSTAR
     
      C = -0.215554 - 0.0810269 * TEMP + 0.00143581 * TEMP2 - 
     &      4.7758E-6 * TEMP3 - 2.91297 / XSTAR
     
      D = -3.58856 + 0.049508 * TEMP - 0.00021382 * TEMP2 + 
     &      3.10801E-7 * TEMP3 - 0.0293333 / XSTAR
     
      E = 1.14598 - 0.600796 * TEMP + 0.00864245 * TEMP2 - 
     &      0.0000228947 * TEMP3 - 8.44985 / XSTAR
     
      F = 2.15855 + 0.0808121 * TEMP - 0.000407382 * TEMP2 - 
     &      4.01957E-7 * TEMP3 + 0.721326 / XSTAR

      G = 1.6241 - 0.0160106 * TEMP + 0.0000377124 * TEMP2 + 
     &      3.21794E-8 * TEMP3 - 0.0113255 / XSTAR
      
      H = 9.71682 - 0.115048 * TEMP + 0.000157098 * TEMP2 + 
     &      4.00914E-7 * TEMP3 + 0.71186 / XSTAR
     
      I = -1.05611 + 0.00903378 * TEMP - 0.0000198417 * TEMP2 +
     &      2.46048E-8 * TEMP3 - 0.0579087 / XSTAR

      J = -0.148712 + 0.00283508 * TEMP - 9.24619E-6 * TEMP2 + 
     &      5.00427E-9 * TEMP3 - 0.0127081 / XSTAR

C.. Calculate the exponential term for the nucleation rate [Eq-12 - Vehkamaki et al., 2002] 
      XFAC = A + B * LNRH + C * LNRH2 + D * LNRH3 + 
     &      E * LNNA + F * LNRH * LNNA + G * LNRH2 * LNNA +     
     &      H * LNNA2 + I * LNRH * LNNA2 + J * LNNA3

C.. Calculate particle nucleation rate: unit [ 1 / cm**3 s] [Eq-12 - Vehkamaki et al., 2002] 
      Jnuc = EXP(XFAC)  
      
C.. The parameterization is valid for nucleation rates of 1.0E-7-1.0E10 [ 1 / cm**3 s]
      Jnuc = MAX (Jnuc, 1.0D-7)
      Jnuc = MIN (1.0D10, Jnuc)

C.. Convert the unit of particle nucleation rate into [ 1 / m**3 s] by multiplying it by 1.0E6  
       DNDT = Jnuc * 1.0E06  ! (1/(m**3 s))

C.. Calculate mass production rate [ ug / (m**3 s) ] analogous to
C   Equation 6a of Binkowski & Roselle (2003). Set the upper limit
C   of the mass production rate as the gas-phase production rate of
C   H2SO4, and adjust the number production rate accordingly.
      DMDT_so4 = sulfmass * DNDT 

      IF ( DMDT_so4 .GT. SO4RATE ) THEN
         DMDT_so4 = SO4RATE
         DNDT = DMDT_SO4 * sulfmass1
      END IF

C.. Calculate the production rate of 2nd moment [ m**2 / (m**3 s) ]
C   This is similar to Equation 6c of Binkowski & Roselle (2003),
C   except the factor of PI is removed and the assumed particle
C   diameter is different.
      DM2DT = DNDT * m2_20

      RETURN    
      END SUBROUTINE NEWPART3 

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE Compute_Flux ( nvolinorg, GNH3R8, GNO3R8, GCLR8, KNH4,
     &                          KNO3, KCL, Ceq, CondRate, Hplus, rate, J, 
     &							COL, ROW, LAYER )
 
C Description
C   Determines the evaporative/condensational flux of volatile
C   inorganic species to aerosol modes. In cases where the resulting H+
C   flux is greater than a specified limit, the Pilinis et al. (2000)
C   AS&T approach is used to modify species vapor pressures such that
C   cond./evap. produces an H+ flux equal to the limit (which is
C   proportional to the current mode concentration of H+).
C   Routine called by VOLINORG.
 
C Arguments
C   Inputs
C     nvolinorg: Number of Volatile inorganic species
C     GNH3R8   : NH3(g) concentration (ug/m3)
C     GNO3R8   : HNO3(g) concentration (ug/m3)
C     GCLR8    : HCl(g) concentration (ug/m3)
C     KNH4     : Index to NH4 species
C     KNO3     : Index to NO3 species
C     KCL      : Index to NO3 species 
C     Ceq      : vapor concentration (mol/m3)
C     CondRate : effective condensation rate (I) of 3rd moment to mode
C              : [treat units as (1/s)]
C     Hplus    : aerosol hydrogen ion concentration (mol/m3) for mode
C     rate     : H2SO4(g) condensation rate (ug/m3/s) for mode
 
C   Output
C     Ceq      : modified vapor concentration (mol/m3)
C     J        : molar cond./evap. flux of volatile inorganics (mol/m3-s)
 
C-----------------------------------------------------------------------

      USE AERO_DATA
      USE PRECURSOR_DATA
      USE AEROMET_DATA

	  use HDMod 
      IMPLICIT NONE

C     Arguments:
      INTEGER      nvolinorg
      TYPE(hyperdual) :: GNH3R8, GNO3R8, GCLR8 ! gas concentrations [ug/m3]
      INTEGER      KNH4, KNO3, KCL       ! Indices to species
      TYPE(hyperdual) :: Ceq( nvolinorg )      ! vapor concentrations [mol/m3]
      TYPE(hyperdual) :: CondRate              ! effective condensation rate (I) for 3rd moment
      TYPE(hyperdual) :: Hplus                 ! hydrogen ion concentration for mode [mol/m3]
      TYPE(hyperdual) :: rate
      TYPE(hyperdual) :: J( nvolinorg )        ! molar cond./evap. flux [mol/m3-s]

C     Local Variables:
      REAL( 8 ),  PARAMETER :: Afact = 1.0D-01  ! factor for H+ limiter
      REAL( 8 ),  PARAMETER :: small = 1.0D-25
      REAL( 8 ),  PARAMETER :: small_prod = 2.0D-22
      TYPE(hyperdual) :: Cinf( nvolinorg ) ! gas concentration in mol/m3
      TYPE(hyperdual) :: Qk              ! factor for modifying vapor press. based on H+ limit
      TYPE(hyperdual) :: Hflux           ! flux of H+ to mode from cond/evap
      TYPE(hyperdual) :: Hlim            ! maximum allowable H+ flux to mode
      TYPE(hyperdual) :: aa, bb, cc      ! terms in quadratic equation
      TYPE(hyperdual) :: JH2SO4          ! molar flux of H2SO4(g) [mol/m3/s]
      TYPE(hyperdual) :: CH2SO4          ! effective H2SO4(g) concentration [mol/m3]
      INTEGER      isp             ! inorganic species index
      INTEGER         :: COL, ROW, LAYER

C-----------------------------------------------------------------------

C     Convert gas concentration from ug/m3 to mol/m3
      Cinf( KNH4 ) = GNH3R8 * 1.0D-6 / PRECURSOR_MW( NH3_IDX )
      Cinf( KNO3 ) = GNO3R8 * 1.0D-6 / PRECURSOR_MW( HNO3_IDX )
      Cinf( KCL )  = GCLR8  * 1.0D-6 / PRECURSOR_MW( HCL_IDX )

C     Calculate cond/evap fluxes (no H+ limiting)
      DO isp = 1, nvolinorg
         J( isp ) = CondRate * ( Cinf( isp ) - Ceq( isp ) )
      END DO

C     Convert rate to mol/m3/s and get effective Cinf for H2SO4(g)
      JH2SO4  = rate * 1.0D-6 / PRECURSOR_MW( SULPRD_IDX )
      CH2SO4  = JH2SO4 / CondRate

C     Limit H+ flux (Pilinis et al., 2000, AS&T). Note: J is flux
C     to entire mode, not one particle
      Hlim  = Afact * Hplus
      Hflux = 2.0D0 * JH2SO4 + J( KNO3 ) + J( KCL ) - J( KNH4 )

C     If Hflux is too large, limit the flux by modifying species
C     vapor pressures with Qk factor (Pilinis et al., 2000, AS&T).
      IF ( ABS( Hflux ) .GT. Hlim ) THEN
         Hlim = SIGN( Hlim, Hflux )

C        Solve quadratic for Qk: aa*Qk^2 + bb*Qk + cc = 0
         aa = Ceq( KCL ) + Ceq( KNO3 )
C
         bb = Hlim / CondRate
     &      + Cinf( KNH4) - Cinf( KNO3 ) - Cinf( KCL ) - 2.0D0 * CH2SO4
         cc = -Ceq( KNH4 )
C
         Qk = 0.0D0 ! initialize Qk

         IF ( aa .LT. small .AND. 0.0D0 .LT. bb ) THEN ! bb*Qk + cc = 0
            Qk = -cc / bb
         ELSE IF (aa .LT. small .AND. bb .LE. 0.0D0 ) THEN
            Qk = 0.0D0
         ELSE IF (-cc .LT. small .AND. bb .LT. 0.0D0 ) THEN  ! aa*Qk^2 + bb*Qk = 0
            Qk = -bb / aa
         ELSE IF (-cc .LT. small .AND. 0.0D0 .LE. bb ) THEN
            Qk = 0.0D0
         ELSE
            Qk = ( -bb + SQRT ( bb**2 - 4.0D0 * aa * cc ) ) / ( 2.0D0 * aa )
            IF ( bb ** 2 - 4.0D0 * aa * cc .LT. 0.0D0 ) THEN
               PRINT *, 'Compute_Flux, sqrt<0'
               Qk = 0.0D0
            END IF
         END IF
		 
C     Modify vapor pressures and get new fluxes
         IF ( Qk .GT. small ) THEN
            Ceq( KNH4 ) = Ceq( KNH4 ) / Qk
            Ceq( KNO3 ) = Ceq( KNO3 ) * Qk
            Ceq( KCl )  = Ceq( KCl )  * Qk
            DO isp = 1, nvolinorg
               J( isp ) = CondRate * ( Cinf( isp ) - Ceq( isp ) )
            END DO
         END IF

      END IF   ! |Hflux| > Hlim

      END SUBROUTINE Compute_Flux

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CALC_H2O ( WI, RH, T, H2O_NEW )

C Description
C   Calculate the water content of aerosol at the new time step.  Water
C   calculations use the ZSR mixing rule with salts determined by the
C   ISORROPIA approach.
C   Routine called by VOLINORG.
 
C Arguments
C   Input
C     WI      : Concentration of components [mol/m^3] at new step
C     RH      : Relative humidity [0-1]
C     T       : Temperature [K]
 
C   Output
C     H2O_NEW : Water [mol/m^3] content at new time step
 
C-----------------------------------------------------------------------

	  use HDMod 
      IMPLICIT NONE

C Parameters:
      INTEGER, PARAMETER :: NCMP = 8, NPAIR = 23
      REAL( 8 ),  PARAMETER :: SMALL = 1.0D-20
      REAL( 8 ),  PARAMETER :: Mw = 0.018D0   ! molar mass H2O (kg/mol)

C Arguments:
      TYPE(hyperdual), INTENT( IN )  :: WI( NCMP )
      REAL( 8 ), INTENT( IN )  :: RH, T
      TYPE(hyperdual), INTENT( OUT ) :: H2O_NEW

C Local Variables:
      CHARACTER( 3 ) :: SC ! subcase for composition
      TYPE(hyperdual) :: FSO4, FNH4, FNA, FNO3, FCL ! "free" ion amounts
      TYPE(hyperdual) :: FCA, FK, FMG            
      TYPE(hyperdual) :: CASO4         ! amount of calcium sulfate, does not participate in ZSR calc
      TYPE(hyperdual) :: WATER         ! kg of water for new time step 
      TYPE(hyperdual) :: X, Y
      TYPE(hyperdual) :: CONC( NCMP )    ! concentration (mol/m^3)
      TYPE(hyperdual) :: CONCR( NPAIR )  ! concentration (mol/m^3) ion "pairs" 
      REAL( 8 ) :: M0I( NPAIR )    ! single-solute molalities
      INTEGER :: J

C-----------------------------------------------------------------------

C     Return if small concentration
      IF ( WI( 1 ) + WI( 2 ) + WI( 3 ) + WI( 4 )
     &   + WI( 5 ) + WI( 6 ) + WI( 7 ) + WI( 8 ) .LE. SMALL) THEN
         H2O_NEW = SMALL
         RETURN
      END IF

C     Set component array (mol/m^3) for determining salts
      CONC = WI   ! array assignment

C     Get the sub-case to use in determining salts
      CALL GETSC ( CONC, RH, T, SC )

#ifdef verbose_aero
!     write( logdev,* ) 'CALC_H2O -SC: ', sc
#endif

C     Initialize ion "pairs" (i.e., salts) used in ZSR
      CONCR( : ) = 0.0D0

C     Depending on case, determine moles of salts in solution (i.e., CONCR)
C     for ZSR calculation below

      IF ( SC .EQ. 'S2' ) THEN    ! sulfate poor (NH4-SO4 system), old K2
         CONCR( 4 )= MIN ( CONC( 2 ), 0.5D0 * CONC( 3 ) )  ! (NH4)2SO4

      ELSE IF ( SC .EQ. 'B4' ) THEN  ! sulfate rich (no acid), old L4, O4
         X = 2.0D0 * CONC( 2 ) - CONC( 3 )     ! 2SO4 - NH4
         Y = CONC( 3 ) - CONC( 2 )            ! NH4 - SO4
         IF ( X .LE. Y ) THEN
            CONCR( 13 ) = X      ! (NH4)3H(SO4)2 is MIN (X,Y)
            CONCR(  4 ) = Y - X  ! (NH4)2SO4
         ELSE
            CONCR( 13 ) = Y      ! (NH4)3H(SO4)2 is MIN (X,Y)
            CONCR(  9 ) = X - Y  ! NH4HSO4
         END IF

      ELSE IF ( SC .EQ. 'C2' ) THEN  ! sulfate rich (free acid), old M2, P2
         CONCR( 9 ) = CONC( 3 )                      ! NH4HSO4
         CONCR( 7 ) = MAX( CONC( 2 ) - CONC( 3 ), 0.0D0 )   ! H2SO4

      ELSE IF ( SC .EQ. 'N3' ) THEN    ! sulfate poor (NH4-SO4-NO3 system)
         CONCR( 4 ) = MIN ( CONC( 2 ), 0.5D0 * CONC( 3 ) )           ! (NH4)2SO4
         FNH4       = MAX ( CONC( 3 ) - 2.0D0 * CONCR( 4 ), 0.0D0 )  ! available NH4
         CONCR( 5 ) = MAX ( MIN ( FNH4, CONC( 4 ) ), 0.0D0 )         ! NH4NO3=MIN(NH4,NO3)

      ELSE IF ( SC .EQ. 'Q5' ) THEN    ! sulfate poor, sodium poor (NH4-SO4-NO3-Cl-Na)
         CONCR( 2 ) = 0.5D0 * CONC( 1 )                              ! Na2SO4
         FSO4       = MAX ( CONC( 2 ) - CONCR( 2 ), 0.0D0 )          ! available SO4
         CONCR( 4 ) = MAX ( MIN ( FSO4, 0.5D0 * CONC( 3 ) ), SMALL ) ! NH42S4=MIN(NH4,S4)
         FNH4       = MAX ( CONC( 3 ) - 2.0D0 * CONCR( 4 ), 0.0D0 )  ! available NH4
         CONCR( 5 ) = MIN ( FNH4, CONC( 4 ) )                        ! NH4NO3=MIN(NH4,NO3)
         FNH4       = MAX ( FNH4 - CONCR( 5 ), 0.0D0 )               ! avaialable NH4
         CONCR( 6 ) = MIN ( FNH4, CONC( 5 ) )                        ! NH4Cl=MIN(NH4,Cl)

      ELSE IF ( SC .EQ. 'R6' ) THEN   ! sulfate poor, sodium rich (NH4-SO4-NO3-Cl-Na)
         CONCR( 2 ) = CONC( 2 )                            ! Na2SO4
         FNA        = MAX ( CONC( 1 ) - 2.0D0 * CONCR( 2 ), 0.0D0 )

         CONCR( 3 ) = MIN ( FNA, CONC( 4 ) )               ! NaNO3
         FNO3       = MAX ( CONC( 4 ) - CONCR( 3 ), 0.0D0 )
         FNA        = MAX ( FNA - CONCR( 3 ), 0.0D0 )

         CONCR( 1 ) = MIN ( FNA, CONC( 5 ) )               ! NaCl
         FCL        = MAX ( CONC( 5 ) - CONCR( 1 ), 0.0D0 )
         FNA        = MAX ( FNA - CONCR( 1 ), 0.0D0 )

         CONCR( 5 ) = MIN ( FNO3, CONC( 3 ) )              ! NH4NO3
         FNO3       = MAX ( FNO3 - CONCR( 5 ), 0.0D0 )
         FNH4       = MAX ( CONC( 3 ) - CONCR( 5 ), 0.0D0 )

         CONCR( 6 ) = MIN (FCL, FNH4 )                     ! NH4Cl

      ELSE IF ( SC .EQ. 'I6' ) THEN   ! sulfate rich (no acid) (NH4-SO4-NO3-Cl-Na)
         CONCR(  2 ) = 0.5D0 * CONC( 1 )                          ! Na2SO4
         FSO4        = MAX ( CONC( 2 ) - CONCR( 2 ), 0.0D0 )
         CONCR( 13 ) = MIN ( CONC( 3 ) / 3.0D0, FSO4 / 2.0D0 )    ! (NH4)3H(SO4)2
         FSO4        = MAX ( FSO4 - 2.0D0 * CONCR( 13 ), 0.0D0 )
         FNH4        = MAX ( CONC( 3 ) - 3.0D0 * CONCR( 13 ), 0.0D0 )

         IF ( FSO4 .LE. SMALL ) THEN    ! reduce (NH4)3H(SO4)2, add (NH4)2SO4
            CONCR( 13 ) = MAX ( CONCR( 13 ) - FNH4, 0.0D0 )   ! (NH4)3H(SO4)2
            CONCR(  4 ) = 2.0D0 * FNH4                  ! (NH4)2SO4
         ELSE IF ( FNH4 .LE. SMALL ) THEN ! reduce (NH4)3H(SO4)2, add NH4HSO4
            CONCR(  9 ) = 3.0D0 * MIN ( FSO4, CONCR( 13 ) ) ! NH4HSO4
            CONCR( 13 ) = MAX ( CONCR( 13 ) - FSO4, 0.0D0 )
            IF ( CONCR( 2 ) .GT. SMALL ) THEN ! reduce Na2SO4, add NaHSO4
               FSO4        = MAX ( FSO4 - CONCR( 9 ) / 3.0D0, 0.0D0 )
               CONCR( 12 ) = 2.0D0 * FSO4                ! NaHSO4
               CONCR(  2 ) = MAX ( CONCR( 2 ) - FSO4, 0.0D0 )  ! Na2SO4
             END IF
         END IF

      ELSE IF ( SC .EQ. 'J3' ) THEN   ! sulfate rich (free acid) (NH4-SO4-NO3-Cl-Na)
         CONCR(  9 ) = CONC( 3 )                             ! NH4HSO4
         CONCR( 12 ) = CONC( 1 )                             ! NAHSO4
         CONCR(  7 ) = MAX ( CONC( 2 ) - CONC( 3 ) - CONC( 1 ), 0.0D0 ) ! H2SO4

      ! Crustal cases
      ELSE IF ( SC .EQ. 'V7' ) THEN  ! sulfate poor, sodium+crustal poor
         CASO4     = MIN ( CONC( 6 ), CONC( 2 ) )            ! CCASO4
         FSO4      = MAX ( CONC( 2 ) - CASO4, 0.0D0 )
         FCA       = MAX ( CONC( 6 ) - CASO4, 0.0D0 )

         CONCR( 17 ) = MIN ( 0.5D0 * CONC( 7 ), FSO4 )       ! CK2SO4
         FK          = MAX ( CONC( 7 ) - 2.D0 * CONCR( 17 ), 0.0D0 )
         FSO4        = MAX ( FSO4 - CONCR( 17 ), 0.0D0 )

         CONCR( 2 )  = MIN ( 0.5D0 * CONC( 1 ), FSO4 )       ! CNA2SO4
         FNA         = MAX ( CONC( 1 ) - 2.0D0 * CONCR( 2 ), 0.0D0 )
         FSO4        = MAX ( FSO4 - CONCR( 2 ), 0.0D0 )

         CONCR( 21 ) = MIN ( CONC( 8 ), FSO4 )               ! CMGSO4
         FMG         = MAX ( CONC( 8 ) - CONCR( 21 ), 0.0D0 )
         FSO4        = MAX ( FSO4 - CONCR( 21 ), 0.0D0 )

         CONCR( 4 )  = MAX ( MIN ( FSO4 , 0.5D0 * CONC( 3 ) ) , SMALL ) ! CNH42S4
         FNH4        = MAX ( CONC( 3 ) - 2.0D0 * CONCR( 4 ), 0.0D0 )

         CONCR( 5 )  = MIN ( FNH4, CONC( 4 ) )               ! CNH4NO3
         FNH4        = MAX ( FNH4 - CONCR( 5 ), 0.0D0 )

         CONCR( 6 )  = MIN ( FNH4, CONC( 5 ) )               ! CNH4CL

      ELSE IF ( SC .EQ. 'U8' ) THEN  ! sulfate poor, crustal+sodium rich, crustal poor
         CASO4       = MIN ( CONC( 6 ), CONC( 2 ) )          ! CCASO4
         FSO4        = MAX ( CONC( 2 ) - CASO4, 0.0D0 )
         FCA         = MAX ( CONC( 6 ) - CASO4, 0.0D0 )

         CONCR( 17 ) = MIN ( 0.5D0 * CONC( 7 ), FSO4 )       ! CK2SO4
         FK          = MAX ( CONC( 7 ) - 2.0D0 * CONCR( 17 ), 0.0D0 )
         FSO4        = MAX ( FSO4 - CONCR( 17 ), 0.0D0 )

         CONCR( 21 ) = MIN ( CONC( 8 ), FSO4 )               ! CMGSO4
         FMG         = MAX ( CONC( 8 ) - CONCR( 21 ), 0.0D0 )
         FSO4        = MAX ( FSO4 - CONCR( 21 ), 0.0D0 )

         CONCR( 2 )  = MAX ( FSO4, 0.0D0 )                   ! CNA2SO4
         FNA         = MAX ( CONC( 1 ) - 2.0D0 * CONCR( 2 ), 0.0D0 )

         CONCR( 3 )  = MIN ( FNA, CONC( 4 ) )                ! NaNO3
         FNO3        = MAX ( CONC( 4 ) - CONCR( 3 ), 0.0D0 )
         FNA         = MAX ( FNA - CONCR( 3 ), 0.0D0 )

         CONCR( 1 )  = MIN ( FNA, CONC( 5 ) )                ! NaCl
         FCL         = MAX ( CONC( 5 ) - CONCR( 1 ), 0.0D0 )
         FNA         = MAX ( FNA - CONCR( 1 ), 0.0D0 )

         CONCR( 5 )  = MIN ( FNO3, CONC( 3 ) )               ! NH4NO3
         FNO3        = MAX ( FNO3 - CONCR( 5 ), 0.0D0 )
         FNH4        = MAX ( CONC( 3 ) - CONCR( 5 ), 0.0D0 )

         CONCR( 6 )  = MIN ( FCL, FNH4 )                     ! NH4Cl
         FCL         = MAX ( FCL - CONCR( 6 ), 0.0D0 )
         FNH4        = MAX ( FNH4 - CONCR( 6 ), 0.0D0 )

      ELSE IF ( SC .EQ. 'W13' ) THEN  ! sulfate poor, crustal+sodium rich
         CASO4       = MIN ( CONC( 2 ), CONC( 6 ) )          ! CASO4
         FCA         = MAX ( CONC( 6 ) - CASO4, 0.0D0 )
         FSO4        = MAX ( CONC( 2 ) - CASO4, 0.0D0 )

         CONCR( 17 ) = MIN ( FSO4, 0.5D0 * CONC( 7 ) )       ! K2SO4
         FK          = MAX ( CONC( 7 ) - 2.0D0 * CONCR( 17 ), 0.0D0 )
         FSO4        = MAX ( FSO4 - CONCR( 17 ), 0.0D0 )

         CONCR( 21 ) = FSO4                                  ! MGSO4
         FMG         = MAX ( CONC( 8 ) - CONCR( 21 ), 0.0D0 )

         CONCR( 1 )  = MIN ( CONC( 1 ), CONC( 5 ) )          ! NACL
         FNA         = MAX ( CONC( 1 ) - CONCR( 1 ), 0.0D0 )
         FCL         = MAX ( CONC( 5 ) - CONCR( 1 ), 0.0D0 )

         CONCR( 16 ) = MIN ( FCA, 0.5D0 * FCL )              ! CACL2
         FCA         = MAX ( FCA - CONCR( 16 ), 0.0D0 )
         FCL         = MAX ( CONC( 5 ) - 2.0D0 * CONCR( 16 ), 0.0D0 )

         CONCR( 15 ) = MIN ( FCA, 0.5D0 * CONC( 4 ) )        ! CA(NO3)2
         FCA         = MAX ( FCA - CONCR( 15 ), 0.0D0 )
         FNO3        = MAX ( CONC( 4 ) - 2.0D0 * CONCR( 15 ), 0.0D0 )

         CONCR( 23 ) = MIN ( FMG, 0.5D0 * FCL )              ! MGCL2
         FMG         = MAX ( FMG - CONCR( 23 ), 0.0D0 )
         FCL         = MAX ( FCL - 2.0D0 * CONCR( 23 ), 0.0D0 )

         CONCR( 22 ) = MIN ( FMG, 0.5D0 * FNO3 )             ! MG(NO3)2
         FMG         = MAX ( FMG - CONCR( 22 ), 0.0D0 )
         FNO3        = MAX ( FNO3 - 2.0D0 * CONCR( 22 ), 0.0D0 )

         CONCR( 3 )  = MIN ( FNA, FNO3 )                     ! NANO3
         FNA         = MAX ( FNA - CONCR( 3 ), 0.0D0 )
         FNO3        = MAX ( FNO3 - CONCR( 3 ), 0.0D0 )

         CONCR( 20 ) = MIN ( FK, FCL )                       ! KCL
         FK          = MAX ( FK - CONCR( 20 ), 0.0D0 )
         FCL         = MAX ( FCL - CONCR( 20 ), 0.0D0 )

         CONCR( 19 ) = MIN ( FK, FNO3 )                      ! KNO3
         FK          = MAX ( FK - CONCR( 19 ), 0.0D0 )
         FNO3        = MAX ( FNO3 - CONCR( 19 ), 0.0D0 )

      ELSE IF ( SC .EQ. 'L9' ) THEN  ! sulfate rich, no free acid
         CASO4       = MIN ( CONC( 6 ), CONC( 2 ) )          ! CCASO4
         FSO4        = MAX ( CONC( 2 ) - CASO4, 0.0D0 )
         FCA         = MAX ( CONC( 6 ) - CASO4, 0.0D0 )

         CONCR( 17 ) = MIN ( 0.5D0 * CONC( 7 ), FSO4 )       ! CK2SO4
         FK          = MAX ( CONC( 7 ) - 2.0D0 * CONCR( 17 ), 0.0D0 )
         FSO4        = MAX ( FSO4 - CONCR( 17 ), 0.0D0 )

         CONCR( 2 )  = MIN ( 0.5D0 * CONC( 1 ), FSO4 )       ! CNA2SO4
         FNA         = MAX ( CONC( 1 ) - 2.0D0 * CONCR( 2 ), 0.0D0 )
         FSO4        = MAX ( FSO4 - CONCR( 2 ), 0.0D0 )

         CONCR( 21 ) = MIN ( CONC( 8 ), FSO4 )               ! CMGSO4
         FMG         = MAX ( CONC( 8 ) - CONCR( 21 ), 0.0D0 )
         FSO4        = MAX ( FSO4 - CONCR( 21 ), 0.0D0 )

         CONCR( 13 ) = MIN ( CONC( 3 ) / 3.0D0, FSO4 / 2.0D0 ) ! CLC
         FSO4        = MAX ( FSO4 - 2.0D0 * CONCR( 13 ), 0.0D0 )
         FNH4        = MAX ( CONC( 3 )- 3.0D0 * CONCR( 13 ),  0.0D0 )

         IF ( FSO4 .LE. SMALL ) THEN                           ! convert (NH4)3H(SO4)2 to (NH4)2SO4
            CONCR( 13 ) = MAX( CONCR( 13 ) - FNH4, 0.0D0 )
            CONCR(  4 ) = 2.0D0 * FNH4                         ! CNH42S4 

         ELSE IF ( FNH4 .LE. SMALL ) THEN                      ! convert (NH4)3H(SO4)2 to NH4HSO4
            CONCR(  9 ) = 3.0D0 * MIN( FSO4, CONCR( 13 ) )     ! CNH4HS4
            CONCR( 13 ) = MAX( CONCR( 13 ) - FSO4, 0.0D0 )     ! CLC, (NH4)3H(SO4)2
            IF ( CONCR( 2 ) .GT. SMALL ) THEN                  ! convert Na2SO4 to NaHSO4
               FSO4        = MAX( FSO4 - CONCR( 9 ) / 3.0D0, 0.0D0 )
               CONCR( 12 ) = 2.0D0 * FSO4                      ! CNAHSO4
               CONCR(  2 ) = MAX( CONCR( 2 ) - FSO4, 0.0D0 )   ! CNA2SO4
            END IF
            IF ( CONCR( 17 ) .GT. SMALL ) THEN                 ! convert K2SO4 to KHSO4
               FSO4        = MAX( FSO4 - CONCR( 9 ) / 3.0D0, 0.0D0 )
               CONCR( 18 ) = 2.0D0 * FSO4                      ! CKHSO4
               CONCR( 17 ) = MAX( CONCR( 17 ) - FSO4, 0.0D0 )  ! CK2SO4
            END IF
         END IF
         
      ELSE IF ( SC .EQ. 'K4' ) THEN ! sulfate super rich, free acid
         CONCR(  9 ) = CONC( 3 )                               ! NH4HSO4 = NH3
         CONCR( 12 ) = CONC( 1 )                               ! NaHSO4  = Na
         CONCR( 18 ) = CONC( 7 )                               ! KHSO4   = K
         CONCR( 21 ) = CONC( 8 )                               ! MgSO4   = Mg
         CONCR(  7 ) = MAX( CONC( 2 ) - CONC( 3 ) - CONC( 1 )
     &                    - CONC( 6 ) - CONC( 7 ) - CONC( 8 ), 0.0D0 ) ! H2SO4 = SO4 - NH4 - Na - Ca - K - Mg

      ELSE
         PRINT*, 'aero_subs.f: case not supported ',
     &          '(metastable reverse only)'
      END IF

C     Get single-solute molalities for ZSR calculation
      CALL GETM0I ( RH, M0I )

C     Calculate H2O with ZSR and determine delta water
      WATER = 0.0D0
      DO J = 1, NPAIR
         WATER = WATER + CONCR( J ) / M0I( J )
      END DO

      WATER = MAX ( WATER, SMALL )
      H2O_NEW = WATER / Mw

      END SUBROUTINE CALC_H2O  

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE GETM0I ( RHIN, M0I )

!!!!!!!!!!!!!!! We want to get rid of this dependency !!!!!!!!!!!!!!!!

C Description
C   Determines single-solute molalities for the 13 possible salts at
C   the ambient RH.  These molalities are used in the ZSR calculation
C   in CALC_H2O. Note that the molalities were determined at the beginning
C   of the time step, and so they are available in the IONS common block
C   of isrpia.inc.
C   Routine called by CALC_H2O.

C Revision History
C   ??? ???? Prakash Bhave(?), Jim Kelly(?): initial revision
C   Apr 2011 Havala Pye: Removed use of IONS common block since it requires
C                        ISORROPIA to be called to setup the initial array values.
C                        Now uses single-solute molalities from ZSR common block
C                        in isrpia.inc that is defined in BLKISO in isocom.f
C                        BE VERY CAREFUL ABOUT IMPLICIT VARIABLES HERE!
 
C Arguments
C   Output
C     M0I : Single-solute molalities (mol/kg-H2O) for 13 salts
 
C-----------------------------------------------------------------------

     

	  use HDMod 
      use ISORROPIAII_MAIN_MOD
	  implicit none
    ! INCLUDE 'h_isrpia.inc'

C Arguments
      INTEGER :: IZ
      REAL( 8 ), INTENT( IN )  :: RHIN

    
      REAL( 8 ),       INTENT( OUT ) :: M0I( NPAIR )
      
CC Location in pure molality array (function of RH)
CC force convert for IZ here     
      IZ = MIN( INT( RHIN * REAL( NZSR, 8 ) + 0.5D0 ), NZSR )
      IZ = MAX( IZ, 1 )
C
C Default value
      M0I = 1.0D+5   ! array assignment

C Actual values (10,11 not provided)

      AWAS(1:10)    = 187.72
      AWAS(11:NZSR) = (/ 
     & 158.13,134.41,115.37,100.10, 87.86, 78.00, 70.00, 63.45, 58.02,
     &  53.46,
     &  49.59, 46.26, 43.37, 40.84, 38.59, 36.59, 34.79, 33.16, 31.67,
     &  30.31,
     &  29.07, 27.91, 26.84, 25.84, 24.91, 24.03, 23.21, 22.44, 21.70,
     &  21.01,
     &  20.34, 19.71, 19.11, 18.54, 17.99, 17.46, 16.95, 16.46, 15.99,
     &  15.54,
     &  15.10, 14.67, 14.26, 13.86, 13.47, 13.09, 12.72, 12.36, 12.01,
     &  11.67,
     &  11.33, 11.00, 10.68, 10.37, 10.06,  9.75,  9.45,  9.15,  8.86,
     &   8.57,
     &   8.29,  8.01,  7.73,  7.45,  7.18,  6.91,  6.64,  6.37,  6.10,
     &   5.83,
     &   5.56,  5.29,  5.02,  4.74,  4.47,  4.19,  3.91,  3.63,  3.34,
     &   3.05,
     &   2.75,  2.45,  2.14,  1.83,  1.51,  1.19,  0.87,  0.56,  0.26,
     &  0.1/)
!
! awsn= sodium nitrate
!
      AWSN(1:10)    = 394.54
      AWSN(11:NZSR) = (/ 
     & 338.91,293.01,254.73,222.61,195.56,172.76,153.53,137.32,123.65,
     & 112.08,
     & 102.26, 93.88, 86.68, 80.45, 75.02, 70.24, 66.02, 62.26, 58.89,
     &  55.85,
     &  53.09, 50.57, 48.26, 46.14, 44.17, 42.35, 40.65, 39.06, 37.57,
     &  36.17,
     &  34.85, 33.60, 32.42, 31.29, 30.22, 29.20, 28.22, 27.28, 26.39,
     &  25.52,
     &  24.69, 23.89, 23.12, 22.37, 21.65, 20.94, 20.26, 19.60, 18.96,
     &  18.33,
     &  17.72, 17.12, 16.53, 15.96, 15.40, 14.85, 14.31, 13.78, 13.26,
     &  12.75,
     &  12.25, 11.75, 11.26, 10.77, 10.29,  9.82,  9.35,  8.88,  8.42,
     &   7.97,
     &   7.52,  7.07,  6.62,  6.18,  5.75,  5.32,  4.89,  4.47,  4.05,
     &   3.64,
     &   3.24,  2.84,  2.45,  2.07,  1.70,  1.34,  0.99,  0.65,  0.31,
     &  0.1/)
!
! awsc= sodium chloride
!
      AWSC(1:10) = 28.16
      AWSC(11:NZSR) = (/
     &  27.17, 26.27, 25.45, 24.69, 23.98, 23.33, 22.72, 22.14, 21.59,
     &  21.08,
     &  20.58, 20.12, 19.67, 19.24, 18.82, 18.43, 18.04, 17.67, 17.32,
     &  16.97,
     &  16.63, 16.31, 15.99, 15.68, 15.38, 15.08, 14.79, 14.51, 14.24,
     &  13.97,
     &  13.70, 13.44, 13.18, 12.93, 12.68, 12.44, 12.20, 11.96, 11.73,
     &  11.50,
     &  11.27, 11.05, 10.82, 10.60, 10.38, 10.16,  9.95,  9.74,  9.52,
     &   9.31,
     &   9.10,  8.89,  8.69,  8.48,  8.27,  8.07,  7.86,  7.65,  7.45,
     &   7.24,
     &   7.04,  6.83,  6.62,  6.42,  6.21,  6.00,  5.79,  5.58,  5.36,
     &   5.15,
     &   4.93,  4.71,  4.48,  4.26,  4.03,  3.80,  3.56,  3.32,  3.07,
     &   2.82,
     &   2.57,  2.30,  2.04,  1.76,  1.48,  1.20,  0.91,  0.61,  0.30,
     &  0.1/)
!
! awac= ammonium chloride
!
      AWAC(1:10) = 1209.00
      AWAC(11:NZSR) = (/
     & 1067.60,949.27,848.62,761.82,686.04,619.16,559.55,505.92,457.25,
     & 412.69,
     & 371.55,333.21,297.13,262.81,229.78,197.59,165.98,135.49,108.57,
     &  88.29,
     &  74.40, 64.75, 57.69, 52.25, 47.90, 44.30, 41.27, 38.65, 36.36,
     &  34.34,
     &  32.52, 30.88, 29.39, 28.02, 26.76, 25.60, 24.51, 23.50, 22.55,
     &  21.65,
     &  20.80, 20.00, 19.24, 18.52, 17.83, 17.17, 16.54, 15.93, 15.35,
     &  14.79,
     &  14.25, 13.73, 13.22, 12.73, 12.26, 11.80, 11.35, 10.92, 10.49,
     &  10.08,
     &   9.67,  9.28,  8.89,  8.51,  8.14,  7.77,  7.42,  7.06,  6.72,
     &   6.37,
     &   6.03,  5.70,  5.37,  5.05,  4.72,  4.40,  4.08,  3.77,  3.45,
     &   3.14,
     &   2.82,  2.51,  2.20,  1.89,  1.57,  1.26,  0.94,  0.62,  0.31,
     &  0.1/)
!
! awss= sodium sulfate
!
      AWSS(1:10) = 24.10
      AWSS(11:NZSR) = (/
     &  23.17, 22.34, 21.58, 20.90, 20.27, 19.69, 19.15, 18.64, 18.17,
     &  17.72,
     &  17.30, 16.90, 16.52, 16.16, 15.81, 15.48, 15.16, 14.85, 14.55,
     &  14.27,
     &  13.99, 13.73, 13.47, 13.21, 12.97, 12.73, 12.50, 12.27, 12.05,
     &  11.84,
     &  11.62, 11.42, 11.21, 11.01, 10.82, 10.63, 10.44, 10.25, 10.07,
     &   9.89,
     &   9.71,  9.53,  9.36,  9.19,  9.02,  8.85,  8.68,  8.51,  8.35,
     &   8.19,
     &   8.02,  7.86,  7.70,  7.54,  7.38,  7.22,  7.06,  6.90,  6.74,
     &   6.58,
     &   6.42,  6.26,  6.10,  5.94,  5.78,  5.61,  5.45,  5.28,  5.11,
     &   4.93,
     &   4.76,  4.58,  4.39,  4.20,  4.01,  3.81,  3.60,  3.39,  3.16,
     &   2.93,
     &   2.68,  2.41,  2.13,  1.83,  1.52,  1.19,  0.86,  0.54,  0.25,
     &  0.1/)
!
! awab= ammonium bisulfate
!
      AWAB(1:10) = 312.84
      AWAB(11:NZSR) = (/
     & 271.43,237.19,208.52,184.28,163.64,145.97,130.79,117.72,106.42,
     &  96.64,
     &  88.16, 80.77, 74.33, 68.67, 63.70, 59.30, 55.39, 51.89, 48.76,
     &  45.93,
     &  43.38, 41.05, 38.92, 36.97, 35.18, 33.52, 31.98, 30.55, 29.22,
     &  27.98,
     &  26.81, 25.71, 24.67, 23.70, 22.77, 21.90, 21.06, 20.27, 19.52,
     &  18.80,
     &  18.11, 17.45, 16.82, 16.21, 15.63, 15.07, 14.53, 14.01, 13.51,
     &  13.02,
     &  12.56, 12.10, 11.66, 11.24, 10.82, 10.42, 10.04,  9.66,  9.29,
     &   8.93,
     &   8.58,  8.24,  7.91,  7.58,  7.26,  6.95,  6.65,  6.35,  6.05,
     &   5.76,
     &   5.48,  5.20,  4.92,  4.64,  4.37,  4.09,  3.82,  3.54,  3.27,
     &   2.99,
     &   2.70,  2.42,  2.12,  1.83,  1.52,  1.22,  0.90,  0.59,  0.28,
     &  0.1/)
!
! awsa= sulfuric acid
!
      AWSA = (/34.00, 33.56, 29.22, 26.55, 24.61, 23.11, 21.89, 20.87,
     &  19.99, 18.45,
     &  17.83, 17.26, 16.73, 16.25, 15.80, 15.38, 14.98, 14.61, 14.26,
     &  13.93,
     &  13.61, 13.30, 13.01, 12.73, 12.47, 12.21, 11.96, 11.72, 11.49,
     &  11.26,
     &  11.04, 10.83, 10.62, 10.42, 10.23, 10.03,  9.85,  9.67,  9.49,
     &   9.31,
     &   9.14,  8.97,  8.81,  8.65,  8.49,  8.33,  8.18,  8.02,  7.87,
     &   7.73,
     &   7.58,  7.44,  7.29,  7.15,  7.01,  6.88,  6.74,  6.61,  6.47,
     &   6.34,
     &   6.21,  6.07,  5.94,  5.81,  5.68,  5.55,  5.43,  5.30,  5.17,
     &   5.04,
     &   4.91,  4.78,  4.65,  4.52,  4.39,  4.26,  4.13,  4.00,  3.86,
     &   3.73,
     &   3.59,  3.45,  3.31,  3.17,  3.02,  2.87,  2.71,  2.56,  2.39,
     &   2.22,
     &   2.05,  1.87,  1.68,  1.48,  1.27,  1.04,  0.80,  0.55,  0.28,
     &  0.1/)
!
! awlc= (NH4)3H(SO4)2
!
      AWLC(1:10) = 125.37
      AWLC(11:NZSR) = (/
     & 110.10, 97.50, 86.98, 78.08, 70.49, 63.97, 58.33, 53.43, 49.14,
     &  45.36,
     &  42.03, 39.07, 36.44, 34.08, 31.97, 30.06, 28.33, 26.76, 25.32,
     &  24.01,
     &  22.81, 21.70, 20.67, 19.71, 18.83, 18.00, 17.23, 16.50, 15.82,
     &  15.18,
     &  14.58, 14.01, 13.46, 12.95, 12.46, 11.99, 11.55, 11.13, 10.72,
     &  10.33,
     &   9.96,  9.60,  9.26,  8.93,  8.61,  8.30,  8.00,  7.72,  7.44,
     &   7.17,
     &   6.91,  6.66,  6.42,  6.19,  5.96,  5.74,  5.52,  5.31,  5.11,
     &   4.91,
     &   4.71,  4.53,  4.34,  4.16,  3.99,  3.81,  3.64,  3.48,  3.31,
     &   3.15,
     &   2.99,  2.84,  2.68,  2.53,  2.37,  2.22,  2.06,  1.91,  1.75,
     &   1.60,
     &   1.44,  1.28,  1.12,  0.95,  0.79,  0.62,  0.45,  0.29,  0.14,
     &  0.1/)
!
! awan= ammonium nitrate
!
      AWAN(1:10) = 960.19
      AWAN(11:NZSR) = (/
     & 853.15,763.85,688.20,623.27,566.92,517.54,473.91,435.06,400.26,
     & 368.89,
     & 340.48,314.63,291.01,269.36,249.46,231.11,214.17,198.50,184.00,
     & 170.58,
     & 158.15,146.66,136.04,126.25,117.24,108.97,101.39, 94.45, 88.11,
     &  82.33,
     &  77.06, 72.25, 67.85, 63.84, 60.16, 56.78, 53.68, 50.81, 48.17,
     &  45.71,
     &  43.43, 41.31, 39.32, 37.46, 35.71, 34.06, 32.50, 31.03, 29.63,
     &  28.30,
     &  27.03, 25.82, 24.67, 23.56, 22.49, 21.47, 20.48, 19.53, 18.61,
     &  17.72,
     &  16.86, 16.02, 15.20, 14.41, 13.64, 12.89, 12.15, 11.43, 10.73,
     &  10.05,
     &   9.38,  8.73,  8.09,  7.47,  6.86,  6.27,  5.70,  5.15,  4.61,
     &   4.09,
     &   3.60,  3.12,  2.66,  2.23,  1.81,  1.41,  1.03,  0.67,  0.32,
     &  0.1/)
!
! awsb= sodium bisulfate
!
      AWSB(1:10) = 55.99
      AWSB(11:NZSR) = (/
     &  53.79, 51.81, 49.99, 48.31, 46.75, 45.28, 43.91, 42.62, 41.39,
     &  40.22,
     &  39.10, 38.02, 36.99, 36.00, 35.04, 34.11, 33.21, 32.34, 31.49,
     &  30.65,
     &  29.84, 29.04, 28.27, 27.50, 26.75, 26.01, 25.29, 24.57, 23.87,
     &  23.17,
     &  22.49, 21.81, 21.15, 20.49, 19.84, 19.21, 18.58, 17.97, 17.37,
     &  16.77,
     &  16.19, 15.63, 15.08, 14.54, 14.01, 13.51, 13.01, 12.53, 12.07,
     &  11.62,
     &  11.19, 10.77, 10.36,  9.97,  9.59,  9.23,  8.87,  8.53,  8.20,
     &   7.88,
     &   7.57,  7.27,  6.97,  6.69,  6.41,  6.14,  5.88,  5.62,  5.36,
     &   5.11,
     &   4.87,  4.63,  4.39,  4.15,  3.92,  3.68,  3.45,  3.21,  2.98,
     &   2.74,
     &   2.49,  2.24,  1.98,  1.72,  1.44,  1.16,  0.87,  0.57,  0.28,
     &  0.1/)
!
! awpc= potassium chloride
!
      AWPC = (/172.62, 165.75, 159.10, 152.67, 146.46, 140.45, 134.64,
     &          129.03, 123.61, 118.38, 113.34, 108.48, 103.79, 99.27,
     &          94.93, 90.74, 86.71, 82.84, 79.11, 75.53, 72.09, 68.79,
     &          65.63, 62.59, 59.68, 56.90, 54.23, 51.68, 49.24, 46.91,
     &          44.68, 42.56, 40.53, 38.60, 36.76, 35.00, 33.33, 31.75,
     &          30.24, 28.81, 27.45, 26.16, 24.94, 23.78, 22.68, 21.64,
     &          20.66, 19.74, 18.86, 18.03, 17.25, 16.51, 15.82, 15.16,
     &          14.54, 13.96, 13.41, 12.89, 12.40, 11.94, 11.50, 11.08,
     &          10.69, 10.32, 9.96, 9.62, 9.30, 8.99, 8.69, 8.40, 8.12,
     &          7.85, 7.59, 7.33, 7.08, 6.83, 6.58, 6.33, 6.08, 5.84,
     &          5.59, 5.34, 5.09, 4.83, 4.57, 4.31, 4.04, 3.76, 3.48,
     &          3.19, 2.90, 2.60, 2.29, 1.98, 1.66, 1.33, 0.99, 0.65,
     &          0.30, 0.1/)
!
! awps= potassium sulfate
!
      AWPS = (/1014.82, 969.72, 926.16, 884.11, 843.54, 804.41, 766.68,
     &          730.32, 695.30, 661.58, 629.14, 597.93, 567.92, 539.09,
     &          511.41, 484.83, 459.34, 434.89, 411.47, 389.04, 367.58,
     &          347.05, 327.43, 308.69, 290.80, 273.73, 257.47, 241.98,
     &          227.24, 213.22, 199.90, 187.26, 175.27, 163.91, 153.15,
     &          142.97, 133.36, 124.28, 115.73, 107.66, 100.08, 92.95,
     &          86.26, 79.99, 74.12, 68.63, 63.50, 58.73, 54.27, 50.14,
     &          46.30, 42.74, 39.44, 36.40, 33.59, 31.00, 28.63, 26.45,
     &          24.45, 22.62, 20.95, 19.43, 18.05, 16.79, 15.64, 14.61,
     &          13.66,  12.81, 12.03, 11.33, 10.68, 10.09, 9.55, 9.06,
     &          8.60, 8.17, 7.76, 7.38, 7.02, 6.66, 6.32, 5.98, 5.65,
     &          5.31, 4.98, 4.64, 4.31, 3.96, 3.62, 3.27, 2.92, 2.57,
     &          2.22, 1.87, 1.53, 1.20, 0.87, 0.57, 0.28, 0.1/)
!
! awpn= potassium nitrate
!
      AWPN(1:44) = 1000.00
      AWPN(45:NZSR) = (/ 953.05, 881.09, 813.39,
     &          749.78, 690.09, 634.14, 581.77, 532.83, 487.16, 444.61,
     &          405.02, 368.26, 334.18, 302.64, 273.51, 246.67, 221.97,
     &          199.31, 178.56, 159.60, 142.33, 126.63, 112.40, 99.54,
     &          87.96, 77.55, 68.24, 59.92, 52.53, 45.98, 40.2, 35.11,
     &          30.65, 26.75, 23.35, 20.40, 17.85, 15.63, 13.72, 12.06,
     &          10.61, 9.35, 8.24, 7.25, 6.37, 5.56, 4.82, 4.12, 3.47,
     &          2.86, 2.28, 1.74, 1.24, 0.79, 0.40, 0.1/)
!
! awpb= potassium bisulfate
!
      AWPB(1:10) = 55.99
      AWPB(11:NZSR) = (/
     &  53.79, 51.81, 49.99, 48.31, 46.75, 45.28, 43.91, 42.62, 41.39,
     &  40.22,
     &  39.10, 38.02, 36.99, 36.00, 35.04, 34.11, 33.21, 32.34, 31.49,
     &  30.65,
     &  29.84, 29.04, 28.27, 27.50, 26.75, 26.01, 25.29, 24.57, 23.87,
     &  23.17,
     &  22.49, 21.81, 21.15, 20.49, 19.84, 19.21, 18.58, 17.97, 17.37,
     &  16.77,
     &  16.19, 15.63, 15.08, 14.54, 14.01, 13.51, 13.01, 12.53, 12.07,
     &  11.62,
     &  11.19, 10.77, 10.36,  9.97,  9.59,  9.23,  8.87,  8.53,  8.20,
     &   7.88,
     &   7.57,  7.27,  6.97,  6.69,  6.41,  6.14,  5.88,  5.62,  5.36,
     &   5.11,
     &   4.87,  4.63,  4.39,  4.15,  3.92,  3.68,  3.45,  3.21,  2.98,
     &   2.74,
     &   2.49,  2.24,  1.98,  1.72,  1.44,  1.16,  0.87,  0.57,  0.28,
     &  0.1/)
!
! awcc= calcium chloride
!
      AWCC = (/19.9, 19.0, 18.15, 17.35, 16.6, 15.89, 15.22, 14.58,
     &          13.99, 13.43, 12.90, 12.41, 11.94, 11.50, 11.09, 10.7,
     &          10.34, 9.99, 9.67, 9.37, 9.09, 8.83, 8.57, 8.34, 8.12,
     &          7.91, 7.71, 7.53, 7.35, 7.19, 7.03, 6.88, 6.74, 6.6,
     &          6.47, 6.35, 6.23, 6.12, 6.01, 5.90, 5.80, 5.70, 5.61,
     &          5.51, 5.42, 5.33, 5.24, 5.16, 5.07, 4.99, 4.91, 4.82,
     &          4.74, 4.66, 4.58, 4.50, 4.42, 4.34, 4.26, 4.19, 4.11,
     &          4.03, 3.95, 3.87, 3.79, 3.72, 3.64, 3.56, 3.48, 3.41,
     &          3.33, 3.25, 3.17, 3.09, 3.01, 2.93, 2.85, 2.76, 2.68,
     &          2.59, 2.50, 2.41, 2.32, 2.23, 2.13, 2.03, 1.93, 1.82,
     &          1.71, 1.59, 1.47, 1.35, 1.22, 1.07, 0.93, 0.77, 0.61,
     &          0.44, 0.25, 0.1/)
!
! awcn= calcium nitrate
!
      AWCN = (/32.89, 31.46, 30.12, 28.84, 27.64, 26.51, 25.44, 24.44,
     &          23.49, 22.59, 21.75, 20.96, 20.22, 19.51, 18.85, 18.23,
     &          17.64, 17.09, 16.56, 16.07, 15.61, 15.17, 14.75, 14.36,
     &          13.99, 13.63, 13.3, 12.98, 12.68, 12.39, 12.11, 11.84,
     &          11.59, 11.35, 11.11, 10.88, 10.66, 10.45, 10.24, 10.04,
     &          9.84, 9.65, 9.46, 9.28, 9.1, 8.92, 8.74, 8.57, 8.4,
     &          8.23, 8.06, 7.9, 7.73, 7.57, 7.41, 7.25, 7.1,6.94, 6.79,
     &          6.63, 6.48, 6.33, 6.18, 6.03, 5.89, 5.74, 5.60, 5.46,
     &          5.32, 5.17, 5.04, 4.9, 4.76, 4.62, 4.49, 4.35, 4.22,
     &          4.08, 3.94, 3.80, 3.66, 3.52, 3.38, 3.23, 3.08, 2.93,
     &          2.77, 2.60, 2.43, 2.25, 2.07, 1.87, 1.67, 1.45, 1.22,
     &          0.97, 0.72, 0.44, 0.14, 0.1/)
!
! awmc= magnesium chloride
!
      AWMC = (/11.24, 10.99, 10.74, 10.5, 10.26, 10.03, 9.81, 9.59,
     &          9.38, 9.18, 8.98, 8.79, 8.60, 8.42, 8.25, 8.07, 7.91,
     &          7.75, 7.59, 7.44, 7.29, 7.15, 7.01, 6.88, 6.75, 6.62,
     &          6.5, 6.38, 6.27, 6.16, 6.05, 5.94, 5.85, 5.75, 5.65,
     &          5.56, 5.47, 5.38, 5.30, 5.22, 5.14, 5.06, 4.98, 4.91,
     &          4.84, 4.77, 4.7, 4.63, 4.57, 4.5, 4.44, 4.37, 4.31,
     &          4.25, 4.19, 4.13, 4.07, 4.01, 3.95, 3.89, 3.83, 3.77,
     &          3.71, 3.65, 3.58, 3.52, 3.46, 3.39, 3.33, 3.26, 3.19,
     &          3.12, 3.05, 2.98, 2.9, 2.82, 2.75, 2.67, 2.58, 2.49,
     &          2.41, 2.32, 2.22, 2.13, 2.03, 1.92, 1.82, 1.71, 1.60,
     &          1.48, 1.36, 1.24, 1.11, 0.98, 0.84, 0.70, 0.56, 0.41,
     &          0.25, 0.1/)
!
! awmn= magnesium nitrate
!
      AWMN = (/12.00, 11.84, 11.68, 11.52, 11.36, 11.2, 11.04, 10.88,
     &          10.72, 10.56, 10.40, 10.25, 10.09, 9.93, 9.78, 9.63,
     &          9.47, 9.32, 9.17, 9.02, 8.87, 8.72, 8.58, 8.43, 8.29,
     &          8.15, 8.01, 7.87, 7.73, 7.59, 7.46, 7.33, 7.2, 7.07,
     &          6.94, 6.82, 6.69, 6.57, 6.45, 6.33, 6.21, 6.01, 5.98,
     &          5.87, 5.76, 5.65, 5.55, 5.44, 5.34, 5.24, 5.14, 5.04,
     &          4.94, 4.84, 4.75, 4.66, 4.56, 4.47, 4.38, 4.29, 4.21,
     &          4.12, 4.03, 3.95, 3.86, 3.78, 3.69, 3.61, 3.53, 3.45,
     &          3.36, 3.28, 3.19, 3.11, 3.03, 2.94, 2.85, 2.76, 2.67,
     &          2.58, 2.49, 2.39, 2.3, 2.2, 2.1, 1.99, 1.88, 1.77, 1.66,
     &          1.54, 1.42, 1.29, 1.16, 1.02, 0.88, 0.73, 0.58, 0.42,
     &          0.25, 0.1/)
!
! awmn= magnesium sulfate
!
      AWMS = (/0.93, 2.5, 3.94, 5.25, 6.45, 7.54, 8.52, 9.40, 10.19,
     &          10.89, 11.50, 12.04, 12.51, 12.90, 13.23, 13.50, 13.72,
     &          13.88, 13.99, 14.07, 14.1, 14.09, 14.05, 13.98, 13.88,
     &          13.75, 13.6, 13.43, 13.25, 13.05, 12.83, 12.61, 12.37,
     &          12.13, 11.88, 11.63, 11.37, 11.12, 10.86, 10.60, 10.35,
     &          10.09, 9.85, 9.6, 9.36, 9.13, 8.9, 8.68, 8.47, 8.26,
     &          8.07, 7.87, 7.69, 7.52, 7.35, 7.19, 7.03, 6.89, 6.75,
     &          6.62, 6.49, 6.37, 6.26, 6.15, 6.04, 5.94, 5.84, 5.75,
     &          5.65, 5.56, 5.47, 5.38, 5.29, 5.20, 5.11, 5.01, 4.92,
     &          4.82, 4.71, 4.60, 4.49, 4.36, 4.24, 4.10, 3.96, 3.81,
     &          3.65, 3.48, 3.30, 3.11, 2.92, 2.71, 2.49, 2.26, 2.02,
     &          1.76, 1.50, 1.22, 0.94, 0.64/)


      M0I( 01 ) = AWSC( IZ )   ! NACl
      M0I( 02 ) = AWSS( IZ )   ! (NA)2SO4
      M0I( 03 ) = AWSN( IZ )   ! NANO3
      M0I( 04 ) = AWAS( IZ )   ! (NH4)2SO4
      M0I( 05 ) = AWAN( IZ )   ! NH4NO3
      M0I( 06 ) = AWAC( IZ )   ! NH4CL
      M0I( 07 ) = AWSA( IZ )   ! 2H-SO4
      M0I( 08 ) = AWSA( IZ )   ! H-HSO4
      M0I( 09 ) = AWAB( IZ )   ! NH4HSO4


      M0I( 12 ) = AWSB( IZ )   ! NAHSO4
      M0I( 13 ) = AWLC( IZ )   ! (NH4)3H(SO4)2

      M0I( 15 ) = AWCN( IZ )   ! CA(NO3)2
      M0I( 16 ) = AWCC( IZ )   ! CACl2
      M0I( 17 ) = AWPS( IZ )   ! K2SO4
      M0I( 18 ) = AWPB( IZ )   ! KHSO4
      M0I( 19 ) = AWPN( IZ )   ! KNO3
      M0I( 20 ) = AWPC( IZ )   ! KCl
      M0I( 21 ) = AWMS( IZ )   ! MGSO4
      M0I( 22 ) = AWMN( IZ )   ! MG(NO3)2
      M0I( 23 ) = AWMC( IZ )   ! MGCL2
    
      END SUBROUTINE GETM0I               

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE GETSC ( CONC, RH, T, SC )
 
C Description
C   Determines the sub-case to use for water uptake calculations.
C   Follows the procedure of ISORROPIA.
C   Routine called by CALC_H2O.
 
C ArgumentS
C   Inputs
C     CONC : Concentration [mol/m^3] of aerosol components. This routine
C            sets minimum CONC to 1.0D-20
C     RH   : Relative humidity
C     T    : Temperature (K)
     
C   Output
C     SC   : Sub-case for aerosol composition
 
C-----------------------------------------------------------------------

	  use HDMod 
      use ISORROPIAII_MAIN_MOD, ONLY : GETASR 
	  
      IMPLICIT NONE

      INTEGER, PARAMETER :: NCMP = 8    ! was NCMP  = 5    ! number of aerosol components
      REAL( 8 ), PARAMETER :: SMALL = 1.0D-20

C Arguments:
!     REAL( 8 ), INTENT( IN )    :: CONC(  NCMP )
      TYPE(hyperdual), INTENT( INOUT ) :: CONC(  NCMP )
      REAL( 8 ), INTENT( IN )    :: RH, T
      CHARACTER( 3 ), INTENT( OUT ) :: SC
            
C Local Variables:
      REAL( 8 ) :: T0, TCF                     ! DRH(T) factor
      TYPE(hyperdual) :: S4RAT, S4RATW, NaRAT, SRI   ! sulfate & sodium ratios
      TYPE(hyperdual) :: CRAT                        ! crustals ratio
      TYPE(hyperdual) :: FSO4                        ! "free" sulfate
      REAL( 8 ) :: DNACL, DNH4CL, DNANO3, DNH4NO3, DNH42S4 ! DRH values

C     TYPE(hyperdual) :: GETASR    ! ISORROPIA function for sulfate ratio

      LOGICAL :: SCP1R, SCP2R, SCP3R, SCP4R ! concentration regime

C-----------------------------------------------------------------------

      SCP1R = .FALSE.
      SCP2R = .FALSE.
      SCP3R = .FALSE.
      SCP4R = .FALSE.

C     See if any components are negligible (see isocom.for)
      IF ( CONC( 1 ) + CONC( 4 ) + CONC( 5 ) + 
     &     CONC( 6 ) + CONC( 7 ) + CONC( 8 ) .LE. SMALL ) THEN       ! Ca,K,Mg,Na,Cl,NO3=0
         SCP1R = .TRUE.                                    
      ELSE IF ( CONC( 1 ) +        CONC( 5 ) +
     &          CONC( 6 ) + CONC( 7 ) + CONC( 8 ) .LE. SMALL ) THEN  ! Ca,K,Mg,Na,Cl=0
         SCP2R = .TRUE.                                     
      ELSE IF ( CONC( 6 ) + CONC( 7 ) + CONC( 8 ) .LE. SMALL ) THEN  ! Ca,K,Mg=0
         SCP3R = .TRUE.                                     
      ELSE                                                           ! all species
         SCP4R = .TRUE.
      END IF

      CONC( : ) = MAX ( CONC( : ), SMALL )

C     Deliquescence RH calculations
      DNH42S4 = 0.7997D0
      DNH4NO3 = 0.6183D0
      IF ( INT( T ) .NE. 298 ) THEN
         T0      = 298.15D0
         TCF     = 1.0D0 / T - 1.0D0 / T0
         DNH4NO3 = DNH4NO3 * EXP( 852.0D0 * TCF )
         DNH42S4 = DNH42S4 * EXP(  80.0D0 * TCF )
         DNH4NO3 = MIN ( DNH4NO3, DNH42S4 ) ! adjust for curves crossing T<271K
      END IF

C     Find sub-case "SC"
      IF ( SCP1R ) THEN ! NH4-S04 system

         IF ( RH .GE. DNH42S4 ) THEN
            S4RATW = GETASR( CONC( 2 ), RH ) ! aerosol sulfate ratio
         ELSE
            S4RATW = 2.0D0                ! dry aerosol sulfate ratio
         END IF
         S4RAT  = CONC( 3 ) / CONC( 2 )     ! sulfate ratio (NH4/SO4)

         IF ( S4RATW .LE. S4RAT ) THEN      ! sulfate poor
            SC = 'S2'
         ELSE IF ( 1.0D0 .LE. S4RAT .AND. S4RAT .LT. S4RATW ) THEN ! sulfate rich (no acid)
            SC = 'B4'
         ELSE IF ( S4RAT .LT. 1.0D0 ) THEN   ! sulfate rich (free acid)
            SC = 'C2'
         END IF

      ELSE IF ( SCP2R ) THEN ! NH4-SO4-NO3 system

         IF ( RH .GE. DNH4NO3 ) THEN
            S4RATW = GETASR( CONC( 2 ), RH )
         ELSE
            S4RATW = 2.0D0               ! dry aerosol ratio
         END IF
         S4RAT = CONC( 3 ) / CONC( 2 )

         IF ( S4RATW .LE. S4RAT ) THEN     ! sulfate poor
            SC = 'N3'
         ELSE IF ( 1.0D0 .LE. S4RAT .AND. S4RAT .LT. S4RATW ) THEN  ! sulfate rich (no acid)
            SC = 'B4'
         ELSE IF ( S4RAT .LT. 1.0D0 ) THEN    ! sulfate rich (free acid)
            SC = 'C2'
         END IF

      ELSE IF ( SCP3R )  THEN ! NH4-SO4-NO3-Na-Cl system

C        Adjust DRH of NH4NO3 for low temperature
         DNACL  = 0.7528D0
         DNANO3 = 0.7379D0
         DNH4CL = 0.7710D0
         IF ( INT( T ) .NE. 298 ) THEN
            DNACL   = DNACL  * EXP(  25.0D0 * TCF )
            DNANO3  = DNANO3 * EXP( 304.0D0 * TCF )
            DNH4CL  = DNH4Cl * EXP( 239.0D0 * TCF )
            DNH4NO3 = MIN ( DNH4NO3, DNH4CL, DNANO3, DNACL )
         END IF

         IF ( RH .GE. DNH4NO3 ) THEN
            FSO4   = CONC( 2 ) - CONC( 1 ) / 2.0D0   ! sulfate unbound by Na+
            FSO4   = MAX ( FSO4, SMALL )
            SRI    = GETASR ( FSO4, RH )
            S4RATW = ( CONC( 1 ) + FSO4 * SRI ) / CONC( 2 )
            S4RATW = MIN ( S4RATW, 2.0D0 )
         ELSE
            S4RATW = 2.0D0                       ! ratio for dry aerosol
         END IF
         S4RAT = ( CONC( 1 ) + CONC( 3 ) ) / CONC( 2 )
         NaRAT = CONC( 1 ) / CONC( 2 )

         IF ( S4RATW .LE. S4RAT .AND. NaRAT .LT. 2.0D0 ) THEN ! sulfate poor, sodium poor
            SC = 'Q5'
         ELSE IF ( S4RAT .GE. S4RATW .AND. NaRAT .GE. 2.0D0 ) THEN ! SO4 poor, Na rich
            SC = 'R6'
         ELSE IF ( 1.0D0 .LE. S4RAT .AND. S4RAT .LT. S4RATW ) THEN ! SO4 rich, no acid
            SC = 'I6'
         ELSE IF ( S4RAT .LT. 1.0D0 ) THEN ! sulfate rich, free acid
            SC = 'J3'
         END IF

      ELSE IF ( SCP4R ) THEN ! NH4-SO4-Na-Cl-Ca-K-Mg system

         ! Do I need an RH if check here????
         FSO4   = CONC( 2 ) - CONC( 1 ) / 2.0D0
     &          - CONC( 6 ) - CONC( 7 ) / 2.0D0 - CONC( 8 )  ! sulfate unbound by sodium,calcium,pottasium,magnesium
         FSO4   = MAX ( FSO4, SMALL )
         SRI    = GETASR( FSO4, RH )                         ! sulfate ratio for NH4+
         S4RATW = ( CONC( 1 ) + FSO4 * SRI + CONC( 6 )
     &            + CONC( 7 ) + CONC( 8 ) ) / CONC( 2 )      ! limiting sulfate ratio
         S4RATW = MIN ( S4RATW, 2.0D0 )
         S4RAT = ( CONC( 1 ) + CONC( 3 ) + CONC( 6 ) + CONC( 7 ) + CONC( 8 ) ) / CONC( 2 ) ! sulfate ratio
         NaRAT = ( CONC( 1 ) + CONC( 6 ) + CONC( 7 ) + CONC( 8 ) ) / CONC( 2 ) ! crustals+sodium ratio
         CRAT  = ( CONC( 6 ) + CONC( 7 ) + CONC( 8 ) ) / CONC( 2 )             ! crustals ratio

         IF ( S4RATW .LE. S4RAT .AND. NaRAT .LT. 2.0D0 ) THEN ! sulfate, sodium, crustal poor
            SC = 'V7'
         ELSE IF ( S4RAT .GE. S4RATW .AND. NaRAT .GE. 2.0D0 ) THEN
            IF ( CRAT .LE. 2.0D0 ) THEN       ! sulfate poor, dust+sodium rich, dust poor
               SC = 'U8'
            ELSE                              ! sulfate poor, dust+sodium rich, dust rich
               SC = 'W13'
            END IF
         ELSE IF ( 1.0D0 .LE. S4RAT .AND. S4RAT .LT. S4RATW ) THEN ! sulfate rich, no acid
            SC = 'L9'
         ELSE IF ( S4RAT .LT. 1.0D0 ) THEN     ! sulfate rich, free acid
            SC = 'K4'
         END IF
      END IF

      !print*,'SUBCASE identified in calc_h2o', SC

      END SUBROUTINE GETSC

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      TYPE(hyperdual) FUNCTION GETAF( NI, NJ, DGNI, DGNJ, XLSGI, XLSGJ, SQRT2)

C  Returns the value of "Xnum" in Equations 10a and 10c
C  of Binkowski and Roselle (2003), given the number concentrations,
C  median diameters, and natural logs of the geometric standard
C  deviations, in two lognormal modes.  The value returned by GETAF
C  is used subsequently in the mode merging calculations:
C       GETAF = ln( Dij / Dgi ) / ( SQRT2 * ln(Sgi) )
C  where Dij is the diameter of intersection,
C        Dgi is the median diameter of the smaller size mode, and
C        Sgi is the geometric standard deviation of smaller mode.
C  A quadratic equation is solved to obtain GETAF, following the
C  method of Press et al.
C 
C  References:
C   1. Binkowski, F.S. and S.J. Roselle, Models-3 Community
C      Multiscale Air Quality (CMAQ) model aerosol component 1:
C      Model Description.  J. Geophys. Res., Vol 108, No D6, 4183
C      doi:10.1029/2001JD001409, 2003.
C   2. Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P.
C      Flannery, Numerical Recipes in Fortran 77 - 2nd Edition.
C      Cambridge University Press, 1992.

	  use HDMod 
      IMPLICIT NONE

      TYPE(hyperdual) :: NI, NJ, DGNI, DGNJ, XLSGI, XLSGJ
      REAL :: SQRT2
      TYPE(hyperdual) :: AA, BB, CC, DISC, QQ, ALFA, L, YJI

C-----------------------------------------------------------------------

C *** Store intermediate values used for the quadratic solution
C     to reduce computational burden
      ALFA = XLSGI / XLSGJ
      YJI = LOG( DGNJ / DGNI ) / ( SQRT2 * XLSGI )
      L = LOG( ALFA * NJ / NI)

C *** Calculate quadratic equation coefficients & discriminant
      AA = 1.0 - ALFA * ALFA
      BB = 2.0 * YJI * ALFA * ALFA
      CC = L - YJI * YJI * ALFA * ALFA
      DISC = BB * BB - 4.0 * AA * CC

C *** If roots are imaginary, return a negative GETAF value so that no
C     mode merging takes place.
      IF ( DISC .LT. 0.0 ) THEN
         GETAF = - 5.0D0       ! error in intersection
         RETURN
      END IF

C *** Equation 5.6.4 of Press et al.
      QQ = -0.5 * ( BB + SIGN( 1.0D0, BB ) * SQRT( DISC ) )

C *** Return solution of the quadratic equation that corresponds to a
C     diameter of intersection lying between the median diameters of
C     the 2 modes.
      GETAF = CC / QQ       ! See Equation 5.6.5 of Press et al.

      RETURN
      END FUNCTION GETAF

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE AERO_INLET ( DGN, XXLSG, RHOP, fPM1, fPM25, fPM10 )

C  Calculates the volume fraction of a given aerosol mode that would enter
C  a sharp-cut inlet, using equations from Jiang et al (2006). The
C  subroutine calculates transmission factors for PM1, PM2.5 and PM10,
C  after first converting those cutoffs from Aerodynamic Diameter to
C  Stokes Diameter (used internally within CMAQ).
 
C  CMAQ aerosols are represented with Stokes Diameter. There are many
C  alternate forms of diameter including the following:
C       Volume Equivalent Diameter (D_ve) - the diameter of a sphere with the
C          same volume as the particle.
C       Stokes Diameter (D_st) - the diameter of a sphere with the same
C          terminal velocity as the real particle. If the shape factor
C          is 1, then the Stokes and Volume Equivalent Diameters are
C          equal.
C       Aerodynamic Diameter (D_ad) - the diameter of a sphere with unit
C          density and the same terminal velocity as the real particle.
C          If the density is 1, then the Aerodynamic Diameter and Stokes
C          Diameter are equal.
C       Vacuum Aerodynamic Diameter (D_va) - the diameter of a sphere with very
C          high Knudsen number and thus a simplified relationship with
C          Stokes diameter that is only dependent on density and shape
C          factor, not Slip Correction. The very high Knudsen number
C          could be because the particle is very small or because the
C          system pressure is very low.
C       Electrical Mobility Diameter (D_em) - the diameter of a sphere with the
C          same migration velocity as the real particle. Requires a
C          slip-correction and a shape factor correction when converting
C          from volume equivalent diameter.
 
C  Conversions:
C       Parameters
C         rho  - Density
C         rho0 - Unit Density (1 g cm-3)
C         Cc( ) - Cunningham Slip Correction of a particular diameter
C         X   - Shape factor = (NonSpherical Drag Force)/(Spherical Drag Force)
C       Aerodynamic -> Stokes
C         D_st = D_ad * sqrt(rho0/rho) * sqrt( Cc(D_ad) / Cc(D_st) )
C       Volume Equivalent -> Stokes
C         D_st = D_ve * sqrt(1/X) * sqrt( Cc(D_ve) / Cc(D_st) )
C       Electrical Mobility -> Volume Equivalent
C         D_ve = D_em * [Cc(D_ve) / Cc(D_em)] * 1/X
C             (If X = 1 (i.e. spherical), then D_st = D_ve = D_em )
C       Vacuum Aerodynamic -> Volume Equivalent (DeCarlo et al., 2004)
C         D_ve = D_va * rho0/rho * X
C             (If X = 1 (i.e. spherical), then D_st = D_ve = D_va * rho0/rho )
 
C  Key Subroutines called: none
 
C  Key Functions called:  ERF
 
C  Revision History:
C    Coded Jul 2005 by Prakash Bhave
C          Apr 2008 J.Kelly: corrected equation for Dst25 calculation
C          Feb 2016 B.Murphy: modified to output PM1 and PM10 mode
C                             fractions in addition to PM2.5
 
C  References:
C   1. Jiang, W., Smyth, S., Giroux, E., Roth, H., Yin, D., Differences
C   between CMAQ fine mode particle and PM2.5 concentrations and their
C   impact on model performance evaluation in the Lower Fraser Valley,
C   Atmos. Environ., 40:4973-4985, 2006.
C   2. Meng, Z., Seinfeld, J.H., On the source of the submicrometer
C   droplet mode of urban and regional aerosols, Aerosol Sci. and
C   Technology, 20:253-265, 1994.
 
C-----------------------------------------------------------------------

	use HDMod 
      IMPLICIT NONE

      INCLUDE SUBST_CONST    ! for PI

C    Input variables
      TYPE(hyperdual)   DGN     ! geometric mean Stokes diameter by number [ m ]
      TYPE(hyperdual)   XXLSG   ! natural log of geometric standard deviation
      TYPE(hyperdual)   RHOP    ! average particle density [ kg/m**3 ]

C    Output variable
      TYPE(hyperdual), INTENT( OUT ) :: fPM1   ! fraction of particulate volume transmitted through 1
      TYPE(hyperdual), INTENT( OUT ) :: fPM25  ! fraction of particulate volume transmitted through 2
      TYPE(hyperdual), INTENT( OUT ) :: fPM10  ! fraction of particulate volume transmitted through 10

C    Internal variables
      REAL, PARAMETER :: SQRT2 = 1.4142136  !  SQRT( 2 )
      REAL, PARAMETER :: D_AD1  = 1.0  ! aerodynamic diameter cut point [ um ]
      REAL, PARAMETER :: D_AD25 = 2.5  ! aerodynamic diameter cut point [ um ]
      REAL, PARAMETER :: D_AD10 = 10.0 ! aerodynamic diameter cut point [ um ]
      REAL, PARAMETER :: B = 0.21470 ! Cunningham slip-correction approx. param [ um ]
                                     ! This factor works well applied to the entire particle size-range
                                     ! The approximation is: Cc(Dp) = 1 + B/Dp
      TYPE(hyperdual) D_ST1, D_ST25, D_ST10     ! Stokes diameter equivalent of DCA..
      TYPE(hyperdual) DG                        ! DGN converted to [ um ]
      TYPE(hyperdual) ERFARG                    ! argument of ERF, from Step#6 of Jiang et al. (2006)

C *** Error function approximation, from Meng & Seinfeld (1994)
      ! Comment these out, wrote ERF_APPROX to replace the statement function
      
      !REAL ERF        ! Error function
      !REAL XX         ! dummy argument for ERF
      !ERF( XX )  = SIGN( 1.0, XX ) * SQRT( 1.0 - EXP( -4.0 * XX * XX / PI ) ) ! THIS IS A STATEMENT FUNCTION, WON'T WORK WITH 

C ----------------------- Begin solution -------------------------------

C *** Calculate Transmission Fractions for Inlets with Aerodynamic
C     Cutoffs
      DG = DGN * 1.0E+06 ! [um] The units need to be equivalent with the B parameter

      ! Convert size cut to equivalent Stokes diameter using
      ! equation 2 of Jiang et al. (2006). Note: the equation in Step 5
      ! of this paper has a typo (i.e., eq. 2 is correct).
      D_ST1  = 0.5D0 * ( SQRT( B ** 2 + 4.0 * D_AD1 *
     &                     ( D_AD1 + B ) * 1.0E+03 / RHOP ) - B )
      D_ST25 = 0.5D0 * ( SQRT( B ** 2 + 4.0 * D_AD25 *
     &                     ( D_AD25 + B ) * 1.0E+03 / RHOP ) - B )
      D_ST10 = 0.5D0 * ( SQRT( B ** 2 + 4.0 * D_AD10 *
     &                     ( D_AD10 + B ) * 1.0E+03 / RHOP ) - B )

      ! Calculate mass fraction with Dca < SizeCut, using ERF approximation
      ! from Meng & Seinfeld (1994) and modified form of Fk(X) equation
      ! in Step#6 of Jiang et al. (2006).
      ERFARG = ( LOG( D_ST1 ) - LOG( DG ) ) / ( SQRT2 * XXLSG ) - 3.0 * XXLSG / SQRT2
      fPM1   = 0.5 * ( 1.0 + ERF( ERFARG ) )
      ERFARG = ( LOG( D_ST25 ) - LOG( DG ) ) / ( SQRT2 * XXLSG ) - 3.0 * XXLSG / SQRT2
      fPM25  = 0.5 * ( 1.0 + ERF( ERFARG ) )
      ERFARG = ( LOG( D_ST10 ) - LOG( DG ) ) / ( SQRT2 * XXLSG ) - 3.0 * XXLSG / SQRT2
      fPM10  = 0.5 * ( 1.0 + ERF( ERFARG ) )

      END SUBROUTINE AERO_INLET

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE AERO_AMS ( M3_WET, M2_WET, M0, M_H2O, RHOP, RHO_H2O, fAMS )

C  The subroutine calculates transmission factors applicable for
C  comparison with AMS measurements using the curve suggested by Jimenez
C  and coworkers  (http://cires1.colorado.edu/jimenez-group/wiki/index.php/
C  FAQ_for_AMS_Data_Users#What_is_the_size_cut_of_AMS_measurements.3F).
C  The AMS measures Vacuum Aerodynamic Diameter (D_va) because the pressure is
C  so low in the instrument. This curve (also described in Ensberg et
C  al., 2013)  includes the following (Note: before applying the curve,
C  the bounds must be converted to Stokes Diameter):
C     D_va < 0.04 um: 0% transmission (T)
C     0.04 < D_va < 0.1 um:  T% = 1666.6 * (D_va - 0.04)   (1)
C     0.1  < D_va < 0.55 um: T% = 100%                     (2)
C     0.55 < D_va < 2.0 um:  T% = 100-68.965*(D_va-0.55)   (3)
C     2.0 um < D_va: T = 0%
 
C  CMAQ aerosols are represented with Stokes Diameter. There are many
C  alternate froms of diameter including the following:
C       Volume Equivalent Diameter (D_ve) - the diameter of a sphere with the
C          same volume as the particle.
C       Stokes Diameter (D_st) - the diameter of a sphere with the same
C          terminal velocity as the real particle. If the shape factor
C          is 1, then the Stokes and Volume Equivalent Diameters are
C          equal.
C       Vacuum Aerodynamic Diameter (D_va) - the diameter of a sphere with very
C          high Knudsen number and thus a simplified relationship with
C          Stokes diameter that is only dependent on density and shape
C          factor, not Slip Correction. The very high Knudsen number
C          could be because the particle is very small or because the
C          system pressure is very low.
 
C  Conversions:
C       Parameters
C         ro  - Density
C         ro0 - Unit Density (1 g cm-3)
C         Cc( ) - Cunningham Slip Correction of a particular diameter
C         X   - Shape factor = (NonSpherical Drag Force)/(Spherical Drag Force)
C       Vacuum Aerodynamic -> Volume Equivalent (DeCarlo et al., 2004)
C         D_ve = D_va * ro0/ro * X
C             (If X = 1 (i.e. spherical), then D_st = D_ve = D_va * ro0/ro )
 
C  Key Functions called:  ERF
 
C  Revision History:
C    Coded Feb 2016 B.Murphy: Created

C  References:
C   1. DeCarlo et al., Particle Morphology and Density Characterization
C   by Combined Mobility and Aerodynamic Diameter Measurements. Part 1:
C   Theory, Aerosol Sci. and Technology, 38:1185-1205, 2004
C   2. Ensberg et al., Inorganic and black carbon aerosols in the Los
C   Angeles Basin during CalNex, Journ. Geophys. Res., 2013.

C-----------------------------------------------------------------------

      USE AERO_DATA, ONLY : MIN_SIGMA_G, MAX_SIGMA_G
      USE AEROMET_DATA, ONLY : F6PI, PI6

      use HDMod 
      IMPLICIT NONE

      INCLUDE SUBST_CONST    ! for PI

C    Input variables
      TYPE(hyperdual), INTENT( IN ) :: M3_WET  ! Third Moment of Wet Distribution (m3/m3)
      TYPE(hyperdual), INTENT( IN ) :: M2_WET  ! Second Moment of Wet Distribution (m2/m3)
      TYPE(hyperdual), INTENT( IN ) :: M0      ! Number of Particles in Distribution (N/m3)
      TYPE(hyperdual), INTENT( IN ) :: M_H2O   ! Mass Conc. of Water in Particles (ug/m3)
      TYPE(hyperdual), INTENT( IN ) :: RHOP    ! average particle density [ kg/m**3 ]
      REAL, INTENT( IN ) :: RHO_H2O ! Water Density [ kg/m3 ]

C    Output variable
      TYPE(hyperdual), INTENT(OUT) :: fAMS   ! fraction of particulate volume transmitted through AMS Inlet

C    Internal variables
      REAL, PARAMETER :: SQRT2 = 1.4142136    !  SQRT( 2 )
      REAL( 8 ), PARAMETER :: DGMIN = 1.0D-9       !  min(Dp) in [m]
      REAL, PARAMETER :: ONETHIRD = 1.0 / 3.0
      REAL, PARAMETER :: TWOTHIRDS = 2.0 * ONETHIRD

      TYPE(hyperdual) DG                        ! DGN converted to [ um ]
      TYPE(hyperdual) DGv                       ! Volume Median Diameter
      TYPE(hyperdual) XXLSG                     ! ln(StndDev) for current mode
      TYPE(hyperdual) M3_DRY, M2_DRY, M3SUBT, M_WET, M_DRY, DRY_DENS, DENSFAC
      TYPE(hyperdual) XFSUM, LXFM2, L2SG, ES36
      TYPE(hyperdual) DBlo_st, DBhi_st, LOG_HI, LOG_LO
      REAL            LOG_HILO, LOG_LOHI
      TYPE(hyperdual) SQRT2LSG, TERM1, TERM2, TERM3, TERM4, TERM5

C *** Error function approximation, from Meng & Seinfeld (1994)
C *** Commented out, rewrite the statement function as ERF_APPROX
C     TYPE(hyperdual) ERF        ! Error function
C     TYPE(hyperdual) XX         ! dummy argument for ERF
C     ERF( XX )  = SIGN( 1.0, XX ) * SQRT( 1.0 - EXP( -4.0 * XX * XX / PI ) )

C ----------------------- Begin solution -------------------------------

C *** First Calculate Parameters of Dry Distribution since this is more
C     applicable in general to AMS measurements.
      M3SUBT = ( 1.0E-9 * F6PI / RHO_H2O ) * M_H2O        ! m3 m-3cd
      M3_DRY = Max(M3_WET - M3SUBT, 0.0) + TINY(0.0)   ! m3 m-3 
      M2_DRY = M2_WET * ( M3_DRY / M3_WET ) ** TWOTHIRDS  ! m2 m-3

      M_WET = M3_WET * 1.0E+9 * PI6 * RHOP ! ug m-3
      M_DRY = Max(M_WET - M_H2O, 0.0) + TINY(0.0) ! m3 m-3 

      DRY_DENS = M_DRY / M3_DRY * F6PI * 1.0E-9   ! kg m-3

      XFSUM = ONETHIRD * Log( M0 ) + TWOTHIRDS * Log( M3_DRY )
      LXFM2 = Log( M2_DRY )
      L2SG = XFSUM - LXFM2   ! ( ln(sigma) )^2

      L2SG = Min( Max( L2SG, LOG( MIN_SIGMA_G ) ** 2 ), LOG( MAX_SIGMA_G ) ** 2 )
      LXFM2 = XFSUM - L2SG
      ES36 = Exp( 4.5 * L2SG )

      DG = Max( DGMIN, ( M3_DRY / ( M0 * ES36 ) ) ** ONETHIRD ) * 1.0E+06 ![um] Units should correspond to D
      XXLSG = Sqrt( L2SG )  ! ln(sigma)

      DGv = EXP( LOG( DG ) + 3.0 * L2SG )

C *** Calculate Transmission Fractions for AMS with Vacuum Aerodynamic
C     Cutoffs. This approximation is split into a piecewise function (see
C     Appendix B in Ensberg et al., 2013).
      fAMS = 0.0D0
      DENSFAC = 1.0E+03 / DRY_DENS  ! Density Correction Factor
                                    !   = rho0 / rho
                                    !   = 1000.0 / rho
      SQRT2LSG = SQRT2 * XXLSG

      ! First Piece of Function [T(%) = 1666.6 * (Dva - 0.04) ]
      DBlo_st = 0.040 * DENSFAC  ! Stokes Lower Bound of Piece [um]
      DBhi_st = 0.100 * DENSFAC  ! Stokes Upper Bound of Piece [um]

      LOG_HILO = LOG( 100.0/40.0 )
      LOG_HI   = LOG( DBhi_st/DGv )
      LOG_LO   = LOG( DBlo_st/DGv )

      TERM1 = LOG( DGv/DBlo_st ) / LOG_HILO
      TERM2 =  ERF( LOG_HI / SQRT2LSG )
     &        -ERF( LOG_LO / SQRT2LSG )
      TERM3 = XXLSG / ( LOG_HILO * (pi/2.0) ** 0.5 )
      TERM4 = EXP( -1.0 * ( LOG_LO / SQRT2LSG ) ** 2 )
      TERM5 = EXP( -1.0 * ( LOG_HI / SQRT2LSG ) ** 2 )

      fAMS = ( TERM1*TERM2  +  TERM3*(TERM4-TERM5) )

      ! Second Piece of Function [T(%) = 100]
      DBlo_st = DBhi_st         ! Stokes Lower Bound [um]
      DBhi_st = 0.55 * DENSFAC  ! Stokes Upper Bound [um]

      TERM1 = ERF( LOG( DBhi_st/DGv ) / SQRT2LSG )
      TERM2 = ERF( LOG( DBlo_st/DGv ) / SQRT2LSG )
      fAMS = fAMS +  (TERM1 - TERM2)

      ! Third Piece of Function [T(%) = 1.0 - 0.6805 * (Dva - 0.55) ]
      DBlo_st = DBhi_st        ! Stokes Lower Bound [um]
      DBhi_st = 2.0 * DENSFAC  ! Stokes Upper Bound [um]

      LOG_LOHI = LOG( 550.0/2000.0 )
      LOG_HI   = LOG( DBhi_st/DGv )
      LOG_LO   = LOG( DBlo_st/DGv )

      TERM1 = LOG( DGv/DBhi_st ) / LOG_LOHI
      TERM2 =  ERF( LOG_HI / SQRT2LSG )
     &        -ERF( LOG_LO / SQRT2LSG )
      TERM3 = XXLSG / ( LOG_LOHI * (pi/2.0) ** 0.5 )
      TERM4 = EXP( -1.0 * ( LOG_LO / SQRT2LSG ) ** 2 )
      TERM5 = EXP( -1.0 * ( LOG_HI / SQRT2LSG ) ** 2 )

      fAMS = fAMS + TERM1 * TERM2  +  TERM3 * (TERM4-TERM5)

      ! Apply the Factor of 0.5 Consistent with Appendix B in Ensberg et
      ! al., 2013. Omit the total mass quantity since these are fractions
      ! we want.
      fAMS = 0.5 * fAMS

      END SUBROUTINE AERO_AMS

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE HETCHEM_EXTRACT_AERO ( CONCVEC, C, R, L, WET_M2, WET_M3, 
     &                                  DE_WET )

c  Calculates the change in aerosol surface area due to heterogeneous
C  reactions and updates the CGRID array with this information.
  
C  History
C  May 2016: BM Created

C  Key Subroutines Called: none
C
C  Called by: Gas-Phase Chemistry Driver
      USE AERO_DATA, ONLY: EXTRACT_AERO, CALCMOMENTS, N_MODE, MOMENT2_CONC,
     &                     MOMENT3_CONC, AEROMODE_DIAM, AEROMODE_LNSG,
     &                     CHEM_M2DRY_INIT, CHEM_M3DRY_INIT, FIXED_SG
      USE AEROMET_DATA, ONLY: pi
      USE GRID_CONF, ONLY: NLAYS, NCOLS, NROWS
      USE PRECURSOR_DATA, ONLY: EXTRACT_PRECURSOR

#ifdef sens
      USE DDM3D_DEFN, ONLY: SENGRID

#endif
	  use HDMod 
      IMPLICIT NONE
      
      TYPE(hyperdual), INTENT( IN ) :: CONCVEC ( : )  ! pointer to model concentrations
      REAL, PARAMETER    :: TWOTHIRDS = 2.0 / 3.0
      TYPE(hyperdual)               :: M2DRY_FINAL( N_MODE ) ! Dry 2nd Aerosol Moment
      TYPE(hyperdual)               :: M3DRY_FINAL( N_MODE ) ! Dry 3rd Aerosol Moment
      INTEGER            :: C, R, L, N
      TYPE(hyperdual), INTENT(OUT) :: WET_M2( : ), WET_M3( : ), DE_WET( : )
      LOGICAL, SAVE      :: FIRSTIME = .TRUE.

      IF ( FIRSTIME ) THEN
          FIRSTIME = .FALSE.
          ALLOCATE( CHEM_M2DRY_INIT( NCOLS,NROWS,NLAYS,N_MODE ),
     &              CHEM_M3DRY_INIT( NCOLS,NROWS,NLAYS,N_MODE ) )

      END IF

C *** extract grid cell concentrations of aero species from CGRID
C     into aerospc_conc in aero_data module
C     also converts dry surface area to wet second moment
#ifdef sens
               CALL EXTRACT_AERO ( CONCVEC(:), .TRUE., SENGRID(NCOLS,NROWS,NLAYS,:,:), .FALSE. )
#else
               CALL EXTRACT_AERO ( CONCVEC(:), .TRUE. )
#endif

C *** extract in inorganic aerosol processors

               CALL EXTRACT_PRECURSOR( CONCVEC(:) ) 
 
C *** Update geometric mean diameters, geometric
C     standard deviations, modal mass totals, and modal particle
C     densities, based on the concentrations of M2, M0, and speciated
C     masses.
               CALL getpar( FIXED_sg )
     
C *** set up variables needed for calculating KN2O5 and YIELD_CLNO2

               DO N = 1, N_MODE
C *** estimate the "wet third moments" from moment3_conc
C     Note: this is the H2O concentration from previous time step
                  WET_M3( N ) = MOMENT3_CONC( N )
C *** calculate "wet second moment" assuming that H2O does not
C     affect the geometric standard deviation
                  WET_M2( N ) = MOMENT2_CONC( N )
C *** The "wet" geometric mean (same as median) diameter was updated in
C     getpar. It is stored in aeromode_diam
C *** calculate effective diameter (this is actually the mean) using Eq 3 of Pleim et al (1995)
                  DE_WET( N ) = AEROMODE_DIAM( N ) * EXP( 1.5 * AEROMODE_LNSG( N ) ** 2 )
               END DO

C *** Retrieve and Save the Dry 3rd and 2nd Moments
               CALL calcmoments( .False. )
               CHEM_M2DRY_INIT( C,R,L,: ) = MOMENT2_CONC( : )
               CHEM_M3DRY_INIT( C,R,L,: ) = MOMENT3_CONC( : )

      END SUBROUTINE HETCHEM_EXTRACT_AERO

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE HETCHEM_UPDATE_AERO ( CGRID )

c  Calculates the change in aerosol surface area due to heterogeneous
C  reactions and updates the CGRID array with this information.
C  
C  History
C  May 2016: BM Created

C  Key Subroutines Called: none
C
C  Called by: Gas-Phase Chemistry Driver
      USE AERO_DATA
      USE AEROMET_DATA, ONLY: pi
      USE GRID_CONF, ONLY: NLAYS, NROWS, NCOLS
	  use HDMod 
	  
      IMPLICIT NONE
      
      TYPE(hyperdual), POINTER   :: CGRID ( :,:,:,: )  ! pointer to model concentrations
      REAL, PARAMETER :: TWOTHIRDS = 2.0 / 3.0
      TYPE(hyperdual)            :: M2DRY_FINAL( N_MODE ) ! Dry 2nd Aerosol Moment
      TYPE(hyperdual)            :: M3DRY_FINAL( N_MODE ) ! Dry 3rd Aerosol Moment
      INTEGER L, R, C, M 

      LOOP_LAY: DO L = 1, NLAYS
         LOOP_ROW: DO R = 1, NROWS
            LOOP_COL: DO C = 1, NCOLS   
      
               ! Extract grid cell concentrations of aero species from CGRID
               ! into aerospc_conc in aero_data module
               ! also converts dry surface area to wet second moment
#ifdef sens
               CALL EXTRACT_AERO ( CGRID( C,R,L,: ), .TRUE., CGRID(C,R,:,: ), .FALSE. )
#else
               CALL extract_aero ( CGRID( C,R,L,: ), .True. )
#endif

               ! Recalculate Dry 3rd Moment After Chemistry Processes
               ! 2nd moment hasnt been update yet.
               CALL calcmoments( .False. )
               M2DRY_FINAL( : ) = MOMENT2_CONC( : )
               M3DRY_FINAL( : ) = MOMENT3_CONC( : )

               ! Calculate new Second Moment Manually
               ! Assume standard deviation is fixed.
               M2DRY_FINAL( : ) = CHEM_M2DRY_INIT( C,R,L,: )  
     &             * ( M3DRY_FINAL( : ) 
     &             / CHEM_M3DRY_INIT( C,R,L,: ) ) ** TWOTHIRDS

               ! Update the Aerosol Surface Area making sure to multiply
               ! the 2nd moment by Pi to convert properly
               DO M = 1, N_MODE
                 CGRID( C,R,L, aerosrf_map( M ) ) = 
     &                        PI * M2DRY_FINAL( M )
               END DO

            END DO LOOP_COL
         END DO LOOP_ROW
      END DO LOOP_LAY

      END SUBROUTINE HETCHEM_UPDATE_AERO
C
 
