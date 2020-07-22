!***********************************************************************
!***********************************************************************
!     Subroutine PLANT
!     This subroutine simulates the growth of the plant using pre-determined
!     conditions.Hourly values of temperature and photosyntetically active
!     radiation come from WEATHER subroutine and daily values of availability
!     of water in the soil come from SW subroutine. This subroutine supplies
!     the SW subroutine with daily values of leaf area index (LAI).
!****************************************************************************

!                  LIST OF VARIABLES
!     di    = daily accumulated temperature above tb (degree days)
!     dLAI  = daily increase in leaf area index (m2/m2/d)
!     dN    = incremental leaf number
!     DOY   = day of the year
!     DYN   = dynamic control variable
!     dw    = incremental total plant dry matter weight (g m-2)
!     dwc   = incremental canopy dry matter weight (g m-2)
!     dwf   = incremental fruit dry matter weight (g m-2)
!     dwr   = incremental root dry matter weight (g m-2)
!     E     = conversion efficiency of CH2O to plant tissue (g g-1)
!     EMP1  = empirical coef. for expoilinear eq.
!     EMP2  = empirical coef. for expoilinear eq.
!     endsim= code signifying physiological maturity (end of simulation)
!     Fc    = fraction of total crop growth partitioned to canopy
!     FL    = code for development phase (1=vegetative phase,
!                 2=reproductive phase)
!     int   = accumulated temperature after reproductive phase starts (c)
!     INTOT = duration of reproductive stage (degree days)
!     LAI   = canopy leaf area index (m2 m-2)
!     Lfmax = maximum number of leaves
!     N     = leaf number
!     nb    = empirical coef. for expoilinear eq.
!     p1    = dry matter of leaves removed per plant per unit development after
!              maximum number of leaves is reached (g)
!     PD    = plant density m-2
!     Pg    = canopy gross photosynthesis rate (g plant-1 day-1)
!     PT    = photosynthesis reduction factor for temp.
!     rm    = maximum rate of leaf appearearance (day-1)
!     sla   = specific leaf area (m2 g-1)
!     SRAD  = Daily solar radiation (MJ m-2)
!     SWFAC1= soil water deficit stress factor
!     SWFAC2= soil water excess stress factor
!     tb    = base temperature above which reproductive growth occurs (c)
!     TMAX  = Daily maximum temperature (c)
!     TMIN  = Daily manimum temperature (c)
!     TMN   = Daily mean temperature (c)
!     W     = total plant dry matter weight (g m-2)
!     Wc    = canopy dry matter weight (g m-2)
!     Wf    = fruit dry matter weight (g m-2)
!     Wr    = root dry matter weight (g m-2)


!***********************************************************************
subroutine plant(&
        doy, endsim, tmax, tmin, par, swfac1, swfac2, & !input
        lai, & !output
        dyn)                                           !control

    !-----------------------------------------------------------------------
    use PlantComponent, only: initialize_from_file, rate, integ, output_to_file, close_file, PlantModel, PlantInput
    use iso_c_binding
    implicit none
    save

    integer(c_int) :: doy, endsim
    real(c_float) :: tmax, tmin, par, swfac1, swfac2, lai
    character(len=10) dyn
    type(PlantModel) :: model

    model%doy = doy
    model%endsim = endsim
    model%tmax = tmax
    model%tmin = tmin
    model%par = par
    model%swfac1 = swfac1
    model%swfac2 = swfac2
    model%lai = lai

    !************************************************************************
    !************************************************************************
    !     initialization
    !************************************************************************
    if (index(dyn, 'INITIAL') /= 0) then
        !************************************************************************
        call initialize_from_file(model)

        !************************************************************************
        !************************************************************************
        !     rate calculations
        !************************************************************************
    elseif (index(dyn, 'RATE') /= 0) then
        !************************************************************************
        call rate(model)

        !************************************************************************
        !************************************************************************
        !     integration
        !************************************************************************
    elseif (index(dyn, 'INTEG') /= 0) then
        !************************************************************************
        call integ(model)

        !************************************************************************
        !************************************************************************
        !     output
        !************************************************************************
    elseif (index(dyn, 'OUTPUT') /= 0) then
        !************************************************************************
        call output_to_file(model)

        !************************************************************************
        !************************************************************************
        !     close
        !************************************************************************
    elseif (index(dyn, 'CLOSE') /= 0) then
        !************************************************************************
        call close_file(model)

        !************************************************************************
        !************************************************************************
        !     end of dynamic 'if' construct
        !************************************************************************
    endif
    !************************************************************************
    doy = model%doy
    endsim = model%endsim
    tmax = model%tmax
    tmin = model%tmin
    par = model%par
    swfac1 = model%swfac1
    swfac2 = model%swfac2
    lai = model%lai
end subroutine plant
!***********************************************************************



!***********************************************************************
!     subroutine lais
!     calculates the canopy leaf area index (lai)
!-----------------------------------------------------------------------
!     input:  fl, di, pd, emp1, emp2, n, nb, swfac1, swfac2, pt, dn
!     output: dlai
!************************************************************************
! pure
subroutine lais(fl, di, pd, emp1, emp2, n, nb, swfac1, swfac2, pt, &
        dn, p1, sla, dlai)

    !-----------------------------------------------------------------------
    implicit none
    save
    real :: pd, emp1, emp2, n, nb, dlai, swfac, a, dn, p1, sla
    real :: swfac1, swfac2, pt, di, fl
    !-----------------------------------------------------------------------

    swfac = min(swfac1, swfac2)
    if (fl == 1.0) then
        a = exp(emp2 * (n - nb))
        dlai = swfac * pd * emp1 * pt * (a / (1 + a)) * dn
    elseif (fl == 2.0) then

        dlai = - pd * di * p1 * sla

    endif
    !-----------------------------------------------------------------------
    return
end subroutine lais
!***********************************************************************



!****************************************************************************
!     subroutine pgs
!     calculates the canopy gross photosysntesis rate (pg)
!*****************************************************************************
! pure
subroutine pgs(swfac1, swfac2, par, pd, pt, lai, pg)

    !-----------------------------------------------------------------------
    implicit none
    save
    real :: par, lai, pg, pt, y1
    real :: swfac1, swfac2, swfac, rowspc, pd

    !-----------------------------------------------------------------------
    !     rowsp = row spacing
    !     y1 = canopy light extinction coefficient

    swfac = min(swfac1, swfac2)
    rowspc = 60.0
    y1 = 1.5 - 0.768 * ((rowspc * 0.01)**2 * pd)**0.1
    pg = pt * swfac * 2.1 * par / pd * (1.0 - exp(-y1 * lai))

    !-----------------------------------------------------------------------
    return
end subroutine pgs
!***********************************************************************



!***********************************************************************
!     subroutine pts
!     calculates the factor that incorporates the effect of temperature
!     on photosynthesis
!************************************************************************
subroutine pts(tmax, tmin, pt)
    !-----------------------------------------------------------------------
    implicit none
    save
    real :: pt, tmax, tmin

    !-----------------------------------------------------------------------
    pt = 1.0 - 0.0025 * ((0.25 * tmin + 0.75 * tmax) - 26.0)**2

    !-----------------------------------------------------------------------
    return
end subroutine pts
!***********************************************************************
!***********************************************************************