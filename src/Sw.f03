!************************************************************************
!************************************************************************
! *     Subroutine SW
!-----------------------------------------------------------------------
! *     This subroutine calculates the soil water availability for the plant,
! *     considering the rain, runoff, deep percolation (drainage) and water
! *     use by the plant (evapotranspiration). It is divided in subroutines
! *     that calculate those parameters separately. Daily data from climate
! *     comes from WEATHER and daily LAI from PLANT subroutines. SW supplies
! *     PLANT with daily soil water factor of availability (SWFAC)
! *
!************************************************************************
! *
! *          LIST OF VARIABLES
! *
! *     CN      = runoff curve number
! *     DATE    = date of irrigation applications (YYDDD)
! *     DOY     = day of year
! *     DP      = depth of the profile considered in the simulation (cm)
! *     DRN     = daily subsurface drainage (mm)
! *     DRNp    = daily drainage percentage (fraction of void space)
! *     DYN     = dynamic control variable
! *     EPa     = actual daily plant transpiration (mm)
! *     EPp     = potential plant transpiration (mm)
! *     ESa     = daily soil evaporation (mm)
! *     ESp     = potential soil evaporation (mm)
! *     ETp     = daily potential evapotranspiration (mm)
! *     FC      = soil water storage at field capacity (mm)
! *     FCp     = water content at field capacity (fraction of void space)
! *     INF     = daily infiltration (mm)
! *     IRR     = daily irrigation (mm)
! *     LAI     = leaf area index (m2/m2)
! *     POTINF  = potential infiltration (mm)
! *     RAIN    = daily rainfall (mm)
! *     ROF     = daily runoff (mm)
! *     SRAD    = solar radiation (mj/m2/day)
! *     ST      = soil water storage at saturation (mm)
! *     STp     = water content saturation (fraction of void space)
! *     SWC     = actual soil water storage in the profile (mm)
! *     SWC_ADJ = cumulative adjustment factor for soil water content (mm)
! *     SWC_INIT= initial soil water content (mm)
! *     SWFAC1  = soil water deficit stress factor
! *     SWFAC2  = soil water excess stress factor
! *     TDRN    = cumulative vertical drainage (mm)
! *     TEPA    = cumulative plant transpiration (mm)
! *     TESA    = cumulative soil evaporation (mm)
! *     TINF    = cumulative infiltration (mm)
! *     TIRR    = cumulative irrigation applied (mm)
! *     TRAIN   = cumulative precipitation (mm)
! *     TROF    = cumulative runoff (mm)
! *     TMAX    = daily maximum temperature (c)
! *     TMIN    = daily minimum temperature (c)
! *     WP      = soil water storage at wilting point (mm)
! *     WPp     = water content at wilting point (fraction of void space)

!************************************************************************
SUBROUTINE SW(&
        DOY, LAI, RAIN, SRAD, TMAX, TMIN, &             !Input
        SWFAC1, SWFAC2, &                               !Output
        DYN)                                            !Control

    !-----------------------------------------------------------------------
    use SoilComponent
    IMPLICIT NONE
    SAVE

    integer :: doy
    real :: lai, rain, srad, tmax, tmin, swfac1, swfac2
    character(len=10) :: dyn
    type(SoilModel) :: m
    type(IrrigationInput) :: irrigation

    m%doy = doy
    m%lai = lai
    m%rain = rain
    m%srad = srad
    m%tmax = tmax
    m%tmin = tmin
    m%swfac1 = swfac1
    m%swfac2 = swfac2

    !************************************************************************
    !************************************************************************
    !     INITIALIZATION
    !************************************************************************
    if (index(dyn, 'INITIAL') .ne. 0) then
        !************************************************************************
        call initialize_from_files(m)

        !************************************************************************
        !************************************************************************
        !     RATE CALCULATIONS
        !************************************************************************
    ELSEIF (INDEX(DYN, 'RATE') .NE. 0) THEN
        !************************************************************************
        call read_irrig(irrigation)
        m%date = irrigation%date
        m%irr  = irrigation%irr
        call rate(m)

        !************************************************************************
        !************************************************************************
        !     INTEGRATION
        !************************************************************************
    ELSEIF (INDEX(DYN, 'INTEG') .NE. 0) THEN
        !************************************************************************
        call integ(m)

        !************************************************************************
        !************************************************************************
        !     OUTPUT
        !************************************************************************
    ELSEIF (INDEX(DYN, 'OUTPUT    ') .NE. 0) THEN
        !************************************************************************
        call output_to_files(m)
        !************************************************************************
        !************************************************************************
        !     CLOSE
        !************************************************************************
    ELSEIF (INDEX(DYN, 'CLOSE') .NE. 0) THEN
        !************************************************************************
        call close_files(m)

        !************************************************************************
        !************************************************************************
        !     End of dynamic 'IF' construct
        !************************************************************************
    ENDIF
    !************************************************************************
    doy = m%doy
    lai = m%lai
    rain = m%rain
    srad = m%srad
    tmax = m%tmax
    tmin = m%tmin
    swfac1 = m%swfac1
    swfac2 = m%swfac2
END SUBROUTINE SW
!************************************************************************


!************************************************************************
!     Subroutine DRAINE
!     Calculates vertical drainage.
!-----------------------------------------------------------------------
!     Input:  SWC, FC, DRNp
!     Output: DRN
!************************************************************************

SUBROUTINE DRAINE(SWC, FC, DRNp, DRN)

    !-----------------------------------------------------------------------
    IMPLICIT NONE
    SAVE
    REAL SWC, FC, DRN, DRNp
    !-----------------------------------------------------------------------

    IF (SWC .GT. FC) THEN
        DRN = (SWC - FC) * DRNp
    ELSE
        DRN = 0
    ENDIF

    !-----------------------------------------------------------------------
    RETURN
END SUBROUTINE DRAINE
!************************************************************************



!************************************************************************
! *     Subroutine ESaS
! *     Calculates the actual daily soil evaporation.
!-----------------------------------------------------------------------
! *     Input:  SWC, WP, FC, ESp
! *     Output: ESa
!************************************************************************

SUBROUTINE ESaS(ESp, SWC, FC, WP, ESa)

    !-----------------------------------------------------------------------
    IMPLICIT NONE
    SAVE
    REAL a, SWC, WP, FC, ESa, ESp
    !-----------------------------------------------------------------------
    IF (SWC .LT. WP) THEN
        a = 0
    ELSEIF (SWC .GT. FC) THEN
        a = 1
    ELSE
        a = (SWC - WP) / (FC - WP)
    ENDIF

    ESa = ESp * a

    !-----------------------------------------------------------------------
    RETURN
END SUBROUTINE ESAS
!************************************************************************



!************************************************************************
! *     Subroutine ETpS
! *     Calculates the daily potential evapotranspiration.
!-----------------------------------------------------------------------
! *     Input:  LAI, TMAX, TMIN, SRAD
! *     Output: ETp
!************************************************************************
! C
! *     Local Variables
! *     ALB  =  ALBEDO OF CROP-SOIL SURFACE
! *     EEQ  =  EQUILIBRIUM EVAPOTRANSPIRATION (mm)
! *     Tmed =  ESTIMATED AVERAGE DAILY TEMPERATURE (C)
! *     f    =

!-----------------------------------------------------------------------
SUBROUTINE ETpS(SRAD, TMAX, TMIN, LAI, ETp)

    !-----------------------------------------------------------------------
    IMPLICIT NONE
    SAVE
    REAL    ALB, EEQ, f, Tmed, LAI
    REAL TMAX, TMIN, SRAD, ETP

    !-----------------------------------------------------------------------
    ALB = 0.1 * EXP(-0.7 * LAI) + 0.2 * (1 - EXP(-0.7 * LAI))
    Tmed = 0.6 * TMAX + 0.4 * TMIN
    EEQ = SRAD * (4.88E-03 - 4.37E-03 * ALB) * (Tmed + 29)

    IF (TMAX .LT. 5) THEN
        f = 0.01 * EXP(0.18 * (TMAX + 20))
    ELSEIF (TMAX .GT. 35) THEN
        f = 1.1 + 0.05 * (TMAX - 35)
    ELSE
        f = 1.1
    ENDIF

    ETp = f * EEQ
    !-----------------------------------------------------------------------
    RETURN
END SUBROUTINE ETPS
!************************************************************************



!************************************************************************
! *     Subroutine RUNOFF
! *     Calculates the daily runoff.
!************************************************************************
! *     Input:  POTINF, CN
! *     Output: ROF
! !-----------------------------------------------------------------------
! *     Local Variables
! *     CN = CURVE NUMBER SCS EQUATION
! *     S  = WATERSHED STORAGE SCS EQUATION (MM)

!-----------------------------------------------------------------------

subroutine runoff_rate(potinf, rof, s)
    implicit none
    real :: potinf, rof, s

    IF (POTINF .GT. 0.2 * S)  THEN
        ROF = ((POTINF - 0.2 * S)**2) / (POTINF + 0.8 * S)
    ELSE
        ROF = 0
    ENDIF
end subroutine runoff_rate

!************************************************************************
! *     Sub-subroutine STRESS calculates soil water stresses.
! *     Today's stresses will be applied to tomorrow's rate calcs.
!-----------------------------------------------------------------------
! *     Input:  SWC, DP, FC, ST, WP
! *     Output: SWFAC1, SWFAC2
!************************************************************************

subroutine stress_integ(SWC, DP, FC, ST, WP, SWFAC1, SWFAC2, the)
    implicit none
    real :: swc, dp, fc, st, wp, swfac1, swfac2, the
    real, parameter :: STRESS_DEPTH = 250
    real :: wtable, dwt
    IF (SWC .LT. WP) THEN
        SWFAC1 = 0.0
    ELSEIF (SWC .GT. THE) THEN
        SWFAC1 = 1.0
    ELSE
        SWFAC1 = (SWC - WP) / (THE - WP)
        SWFAC1 = MAX(MIN(SWFAC1, 1.0), 0.0)
    ENDIF

    !-----------------------------------------------------------------------
    !     Excess water stress factor - SWFAC2
    !-----------------------------------------------------------------------
    IF (SWC .LE. FC) THEN
        WTABLE = 0.0
        DWT = DP * 10.              !DP in cm, DWT in mm
        SWFAC2 = 1.0
    ELSE
        !FC water is distributed evenly throughout soil profile.  Any
        !  water in excess of FC creates a free water surface
        !WTABLE - thickness of water table (mm)
        !DWT - depth to water table from surface (mm)
        WTABLE = (SWC - FC) / (ST - FC) * DP * 10.
        DWT = DP * 10. - WTABLE

        IF (DWT .GE. STRESS_DEPTH) THEN
            SWFAC2 = 1.0
        ELSE
            SWFAC2 = DWT / STRESS_DEPTH
        ENDIF
        SWFAC2 = MAX(MIN(SWFAC2, 1.0), 0.0)
    ENDIF
end subroutine stress_integ

!************************************************************************
! *     Subroutine WBAL
! *     Seasonal water balance
!-----------------------------------------------------------------------
!     Input:  SWC, SWC_INIT, TDRN, TEPA,
!                 TESA, TIRR, TRAIN, TROF
!     Output: None
!************************************************************************

SUBROUTINE WBAL(SWC_INIT, SWC, TDRN, TEPA, &
        TESA, TIRR, TRAIN, TROF, SWC_ADJ, TINF)

    !-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, PARAMETER :: LSWC = 21
    REAL SWC, SWC_INIT
    REAL TDRN, TEPA, TESA, TIRR, TRAIN, TROF
    REAL WATBAL, SWC_ADJ, TINF
    REAL CHECK
    !-----------------------------------------------------------------------
    OPEN (LSWC, FILE = 'output/wbal.out', STATUS = 'REPLACE')

    WATBAL = (SWC_INIT - SWC) + (TRAIN + TIRR) - &
            !      0.0   =(Change in storage)+  (Inflows)    -
            (TESA + TEPA + TROF + TDRN)
    !                         (Outflows)

    WRITE(*, 100)   SWC_INIT, SWC, TRAIN, TIRR, TESA, TEPA, TROF, TDRN
    WRITE(LSWC, 100)SWC_INIT, SWC, TRAIN, TIRR, TESA, TEPA, TROF, TDRN
    100 FORMAT(//, 'SEASONAL SOIL WATER BALANCE', //, &
            'Initial soil water content (mm):', F10.3, /, &
            'Final soil water content (mm):  ', F10.3, /, &
            'Total rainfall depth (mm):      ', F10.3, /, &
            'Total irrigation depth (mm):    ', F10.3, /, &
            'Total soil evaporation (mm):    ', F10.3, /, &
            'Total plant transpiration (mm): ', F10.3, /, &
            'Total surface runoff (mm):      ', F10.3, /, &
            'Total vertical drainage (mm):   ', F10.3, /)

    IF (SWC_ADJ .NE. 0.0) THEN
        WRITE(*, 110) SWC_ADJ
        WRITE(LSWC, 110) SWC_ADJ
        110   FORMAT('Added water for SWC<0 (mm):     ', E10.3, /)
    ENDIF

    WRITE(*, 200) WATBAL
    WRITE(LSWC, 200) WATBAL
    200 FORMAT('Water Balance (mm):             ', F10.3, //)

    CHECK = TRAIN + TIRR - TROF
    IF ((CHECK - TINF) .GT. 0.0005) THEN
        WRITE(*, 300) CHECK, TINF, (CHECK - TINF)
        WRITE(LSWC, 300) CHECK, TINF, (CHECK - TINF)
        300   FORMAT(/, 'Error: TRAIN + TIRR - TROF = ', F10.4, /, &
                'Total infiltration =         ', F10.4, /, &
                'Difference =                 ', F10.4)
    ENDIF

    CLOSE (LSWC)

    !-----------------------------------------------------------------------
    RETURN
END SUBROUTINE WBAL


!************************************************************************
!************************************************************************
