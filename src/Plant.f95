module PlantComponent
    use, intrinsic :: iso_c_binding

    type, public, bind(C) :: PlantInput
        integer(c_int) :: doy, endsim
        real(c_float) :: tmax, tmin, par, swfac1, swfac2
    end type PlantInput

    type, public, bind(C) :: PlantModel
        real(c_float) :: e, fc, lai, nb, n, pt, pg, di, par
        real(c_float) :: rm, dwf, intc, tmax, tmin, p1, sla
        real(c_float) :: pd, emp1, emp2, lfmax, dwc, tmn
        real(c_float) :: dwr, dw, dn, w, wc, wr, wf, tb, intot, dlai, fl
        real(c_float) :: swfac1, swfac2
        integer(c_int) :: doy, endsim, count
    end type PlantModel
contains
    subroutine initialize_from_file(this)
        use iso_c_binding
        implicit none
        type(PlantModel) :: this
        this%endsim = 0

        open (2, file = 'data/plant.inp', status = 'UNKNOWN')
        open (1, file = 'output/plant.out', status = 'REPLACE')

        read(2, 10) this%lfmax, this%emp2, this%emp1, this%pd, this%nb, this%rm, this%fc, this%tb &
                , this%intot, this%n, this%lai, this%w, this%wr, this%wc &
                , this%p1, this%sla
        10 format(17(1x, f7.4))
        close(2)

        write(1, 11)
        write(1, 12)
        11 format('results of plant growth simulation: ')
        12 format(/ &
                /, '                accum', &
                /, '       number    temp                                    leaf', &
                /, '  day      of  during   plant  canopy    root   fruit    area', &
                /, '   of    leaf  reprod  weight  weight  weight  weight   index', &
                /, ' year   nodes    (oc)  (g/m2)  (g/m2)  (g/m2)  (g/m2) (m2/m2)', &
                /, ' ----  ------  ------  ------  ------  ------  ------  ------')

        write(*, 11)
        write(*, 12)

        this%count = 0
    end subroutine initialize_from_file

    subroutine c_inititialize_from_file(this) bind(c, name='pm_initialize_from_file')
        use iso_c_binding
        implicit none
        type(PlantModel) :: this
        call initialize_from_file(this)
    end subroutine c_inititialize_from_file

    subroutine output_to_file(this)
        use iso_c_binding
        implicit none
        type(PlantModel) :: this

        write(1, 20) this%doy, this%n, this%intc, this%w, this%wc, this%wr, this%wf, this%lai
        20 format(i5, 7f8.2)

        if (this%count == 23) then
            this%count = 0
            write(*, 30)
            30 format(2/)
            31 format(/ &
                /, '                accum', &
                /, '       number    temp                                    leaf', &
                /, '  day      of  during   plant  canopy    root   fruit    area', &
                /, '   of    leaf  reprod  weight  weight  weight  weight   index', &
                /, ' year   nodes    (oc)  (g/m2)  (g/m2)  (g/m2)  (g/m2) (m2/m2)', &
                /, ' ----  ------  ------  ------  ------  ------  ------  ------')
            write(*, 31)
        endif

        this%count = this%count + 1
        write(*, 20) this%doy, this%n, this%intc, this%w, this%wc, this%wr, this%wf, this%lai
    end subroutine output_to_file

    subroutine c_output_to_file(this) bind(c, name='pm_output_to_file')
        use iso_c_binding
        implicit none
        type(PlantModel) :: this
        call output_to_file(this)
    end subroutine c_output_to_file

    subroutine close_file(this)
        implicit none
        type(PlantModel) :: this
        close(1)
    end subroutine close_file

    subroutine c_close(this) bind(c, name='pm_close')
        implicit none
        type(PlantModel) :: this
        call close_file(this)
    end subroutine c_close

    subroutine rate(this)
        use iso_c_binding
        implicit none
        type(PlantModel) :: this

        this%tmn = 0.5 * (this%tmax + this%tmin)
        call pts(this%tmax, this%tmin, this%pt)
        call pgs(this%swfac1, this%swfac2, this%par, this%pd, this%pt, this%lai, this%pg)

        if (this%n < this%lfmax) then
            !         vegetative phase
            this%fl = 1.0
            this%e = 1.0
            this%dn = this%rm * this%pt

            call lais(this%fl, this%di, this%pd, this%emp1, this%emp2, this%n, this%nb, this%swfac1, this%swfac2, this%pt, &
                    this%dn, this%p1, this%sla, this%dlai)
            this%dw = this%e * (this%pg) * this%pd
            this%dwc = this%fc * this%dw
            this%dwr = (1 - this%fc) * this%dw
            this%dwf = 0.0

        else
            !         reproductive phase
            this%fl = 2.0

            if (this%tmn >= this%tb .and. this%tmn <= 25) then
                this%di = (this%tmn - this%tb)
            else
                this%di = 0.0
            endif

            this%intc = this%intc + this%di
            this%e = 1.0
            call lais(this%fl, this%di, this%pd, this%emp1, this%emp2, this%n, this%nb, this%swfac1, this%swfac2, this%pt, &
                    this%dn, this%p1, this%sla, this%dlai)
            this%dw = this%e * (this%pg) * this%pd
            this%dwf = this%dw
            this%dwc = 0.0
            this%dwr = 0.0
            this%dn = 0.0
        endif
    end subroutine rate

    subroutine c_rate(this) bind(c, name='pm_rate')
        implicit none
        type(PlantModel) :: this
        call rate(this)
    end

    subroutine update(this, input)
        use iso_c_binding
        implicit none
        type(PlantModel) :: this
        type(PlantInput) :: input
        this%tmax = input%tmax
        this%tmin = input%tmin
        this%par  = input%par
        this%swfac1 = input%swfac1
        this%swfac2 = input%swfac2
    end subroutine update

    subroutine c_update(this, input) bind(c, name='pm_update')
        implicit none
        type(PlantModel) :: this
        type(PlantInput) :: input
        call update(this, input)
    end subroutine c_update

    subroutine integ(this)
        use iso_c_binding
        implicit none
        type(PlantModel) :: this

        this%lai = this%lai + this%dlai
        this%w = this%w + this%dw
        this%wc = this%wc + this%dwc
        this%wr = this%wr + this%dwr
        this%wf = this%wf + this%dwf

        this%lai = max(this%lai, 0.0)
        this%w = max(this%w, 0.0)
        this%wc = max(this%wc, 0.0)
        this%wr = max(this%wr, 0.0)
        this%wf = max(this%wf, 0.0)

        this%n = this%n + this%dn
        if (this%intc > this%intot) then
            this%endsim = 1
            return
        endif
    end subroutine integ

    subroutine c_integ(this) bind(c, name='pm_integ')
        implicit none
        type(PlantModel) :: this
        call integ(this)
    end subroutine c_integ
end module PlantComponent

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