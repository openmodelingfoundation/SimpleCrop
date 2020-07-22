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
