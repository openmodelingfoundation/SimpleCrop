module SoilComponent
    use, intrinsic :: iso_c_binding
    implicit none

    type, public, bind(c) :: SoilInput
        integer(c_int) :: doy
        real(c_float) :: lai, rain, srad, tmax, tmin, swfac1, swfac2
    end type SoilInput

    type, public, bind(c) :: SoilModel
        integer(c_int) :: date, doy
        real(c_float) ::  srad, tmax, tmin, rain, swc, inf, irr, rof, esa, epa, drnp
        real(c_float) ::  drn, dp, wpp, fcp, stp, wp, fc, st, esp, epp, etp, lai
        real(c_float) ::  cn, swfac1, swfac2, potinf
        real(c_float) ::  swc_init, train, tirr, tesa, tepa, trof, tdrn
        real(c_float) ::  tinf, swc_adj
        real(c_float) ::  s, the
    end type SoilModel

    type, public, bind(c) :: IrrigationInput
        integer(c_int) :: date
        real(c_float) :: irr
    end type IrrigationInput
contains
    subroutine initialize_from_files(m)
        implicit none
        type(SoilModel) :: m

        OPEN(3, FILE = 'data/soil.inp', action = 'read', STATUS = 'OLD')
        OPEN(10, FILE = 'output/soil.out', action = 'write', STATUS = 'REPLACE')
        OPEN(11, FILE = 'data/irrig.inp', action = 'read', STATUS = 'OLD')

        READ(3, '(5X,F5.2,5X,F5.2,5X,F5.2,5X,F7.2,5X,F5.2,5X,F5.2,5X,F5.2)') &
                m%WPp, m%FCp, m%STp, m%DP, m%DRNp, m%CN, m%SWC
        CLOSE(3)

        WRITE(10, 15)
        15  FORMAT('Results of soil water balance simulation:', &
                /, 105X, 'Soil', /, 73X, 'Pot.  Actual  Actual    Soil   Water', &
                10X, 'Excess', /, '  Day   Solar     Max     Min', 42X, &
                'Evapo-    Soil   Plant   Water Content Drought   Water', &
                /, '   of    Rad.    Temp    Temp    Rain   Irrig  Runoff', &
                '   Infil   Drain   Trans   Evap.  Trans. content   (mm3/', &
                '  Stress  Stress', /, ' Year (MJ/m2)    (oC)    (oC)    (mm)', &
                '    (mm)    (mm)    (mm)    (mm)    (mm)    (mm)    (mm)', &
                '    (mm)    mm3)  Factor  Factor')

        m%WP = m%DP * m%WPp * 10.0
        m%FC = m%DP * m%FCp * 10.0
        m%ST = m%DP * m%STp * 10.0
        m%SWC_INIT = m%SWC

        m%s = 254 * (100 / m%CN - 1)
        m%the = m%WP + 0.75 * (m%FC - m%WP)
        call stress_integ(m%SWC, m%DP, m%FC, m%ST, m%WP, m%SWFAC1, m%SWFAC2, m%the)

        !     Keep totals for water balance
        m%TRAIN = 0.0
        m%TIRR = 0.0
        m%TESA = 0.0
        m%TEPA = 0.0
        m%TROF = 0.0
        m%TDRN = 0.0
        m%TINF = 0.0
        m%SWC_ADJ = 0.0
    end subroutine initialize_from_files

    subroutine read_irrig(irrigation)
        implicit none
        type(IrrigationInput) :: irrigation
        read(11, 25) irrigation%DATE, irrigation%IRR
        25 format(i5, 2x, f4.1)
    end subroutine read_irrig

    subroutine output_to_files(m)
        implicit none
        type(SoilModel) :: m
        WRITE(10, '(I5,3F8.1,9F8.2,3F8.3)') m%DOY, m%SRAD, m%TMAX, m%TMIN, m%RAIN, m%IRR, m%ROF, m%INF, m%DRN, &
                m%ETP, m%ESa, m%EPa, m%SWC, m%SWC / m%DP, m%SWFAC1, m%SWFAC2
    end subroutine output_to_files

    subroutine close_files(m)
        implicit none
        type(SoilModel) :: m
        CALL WBAL(m%SWC_INIT, m%SWC, m%TDRN, m%TEPA, &
                m%TESA, m%TIRR, m%TRAIN, m%TROF, m%SWC_ADJ, m%TINF)
        CLOSE(10)
        CLOSE(11)
    end subroutine close_files

    subroutine rate(m)
        implicit none
        type(SoilModel) :: m

        m%TIRR = m%TIRR + m%IRR
        m%POTINF = m%RAIN + m%IRR
        m%TRAIN = m%TRAIN + m%RAIN
        CALL DRAINE(m%SWC, m%FC, m%DRNp, m%DRN)

        IF (m%POTINF .GT. 0.0) THEN
            CALL runoff_rate(m%POTINF, m%ROF, m%s)
            m%INF = m%POTINF - m%ROF
        ELSE
            m%ROF = 0.0
            m%INF = 0.0
        ENDIF

        !     Potential evapotranspiration (ETp), soil evaporation (ESp) and
        !       plant transpiration (EPp)
        CALL ETpS(m%SRAD, m%TMAX, m%TMIN, m%LAI, m%ETp)
        m%ESp = m%ETp * EXP(-0.7 * m%LAI)
        m%EPp = m%ETp * (1 - EXP(-0.7 * m%LAI))

        !     Actual soil evaporation (ESa), plant transpiration (EPa)
        CALL ESaS(m%ESp, m%SWC, m%FC, m%WP, m%ESa)
        m%EPa = m%EPp * MIN(m%SWFAC1, m%SWFAC2)
    end subroutine rate

    subroutine c_rate(m) bind(c, name = 'soil_rate')
        implicit none
        type(SoilModel) :: m
        call rate(m)
    end subroutine c_rate

    subroutine integ(m)
        type(SoilModel) :: m
        m%SWC = m%SWC + (m%INF - m%ESa - m%EPa - m%DRN)

        IF (m%SWC .GT. m%ST) THEN
            m%ROF = m%ROF + (m%SWC - m%ST)
            m%SWC = m%ST
        ENDIF

        IF (m%SWC .LT. 0.0) THEN
            m%SWC_ADJ = m%SWC_ADJ - m%SWC
            m%SWC = 0.0
        ENDIF

        m%TINF = m%TINF + m%INF
        m%TESA = m%TESA + m%ESA
        m%TEPA = m%TEPA + m%EPA
        m%TDRN = m%TDRN + m%DRN
        m%TROF = m%TROF + m%ROF

        call stress_integ(m%SWC, m%DP, m%FC, m%ST, m%WP, m%SWFAC1, m%SWFAC2, m%the)
    end subroutine integ

    subroutine c_integ(m) bind(c, name = 'soil_integ')
        implicit none
        type(SoilModel) :: m
        call integ(m)
    end subroutine c_integ
end module SoilComponent