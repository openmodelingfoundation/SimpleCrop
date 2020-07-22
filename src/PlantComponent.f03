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
    subroutine initialize_from_file(m)
        use iso_c_binding
        implicit none
        type(PlantModel) :: m
        m%endsim = 0

        open (2, file = 'data/plant.inp', status = 'OLD')
        open (1, file = 'output/plant.out', status = 'REPLACE')

        read(2, 10) m%lfmax, m%emp2, m%emp1, m%pd, m%nb, m%rm, m%fc, m%tb &
                , m%intot, m%n, m%lai, m%w, m%wr, m%wc &
                , m%p1, m%sla
        10 format(17(1x, f7.4))
        close(2)

        write(1, 11)
        write(1, 12)
        11 format('results of plant growth simulation: ')
        12 format(/ &
                /, '                Accum', &
                /, '       Number    Temp                                    Leaf', &
                /, '  Day      of  during   Plant  Canopy    Root   Fruit    Area', &
                /, '   of    Leaf  Reprod  Weight  Weight  Weight  weight   Index', &
                /, ' Year   Nodes    (oC)  (g/m2)  (g/m2)  (g/m2)  (g/m2) (m2/m2)', &
                /, ' ----  ------  ------  ------  ------  ------  ------  ------')

        write(*, 11)
        write(*, 12)

        m%count = 0
    end subroutine initialize_from_file

    subroutine c_inititialize_from_file(m) bind(c, name='pm_initialize_from_file')
        use iso_c_binding
        implicit none
        type(PlantModel) :: m
        call initialize_from_file(m)
    end subroutine c_inititialize_from_file

    subroutine output_to_file(m)
        use iso_c_binding
        implicit none
        type(PlantModel) :: m

        write(1, 20) m%doy, m%n, m%intc, m%w, m%wc, m%wr, m%wf, m%lai
        20 format(i5, 7f8.2)

        if (m%count == 23) then
            m%count = 0
            write(*, 30)
            30 format(2/)
            31 format(/ &
                /, '                Accum', &
                /, '       Number    Temp                                    Leaf', &
                /, '  Day      of  During   Plant  Canopy    Root   Fruit    Area', &
                /, '   of    Leaf  Reprod  Weight  Weight  Weight  Weight   Index', &
                /, ' Year   Nodes    (oC)  (g/m2)  (g/m2)  (g/m2)  (g/m2) (m2/m2)', &
                /, ' ----  ------  ------  ------  ------  ------  ------  ------')
            write(*, 31)
        endif

        m%count = m%count + 1
        write(*, 20) m%doy, m%n, m%intc, m%w, m%wc, m%wr, m%wf, m%lai
    end subroutine output_to_file

    subroutine c_output_to_file(m) bind(c, name='pm_output_to_file')
        use iso_c_binding
        implicit none
        type(PlantModel) :: m
        call output_to_file(m)
    end subroutine c_output_to_file

    subroutine close_file(m)
        implicit none
        type(PlantModel) :: m
        close(1)
    end subroutine close_file

    subroutine c_close(m) bind(c, name='pm_close')
        implicit none
        type(PlantModel) :: m
        call close_file(m)
    end subroutine c_close

    subroutine rate(m)
        use iso_c_binding
        implicit none
        type(PlantModel) :: m

        m%tmn = 0.5 * (m%tmax + m%tmin)
        call pts(m%tmax, m%tmin, m%pt)
        call pgs(m%swfac1, m%swfac2, m%par, m%pd, m%pt, m%lai, m%pg)

        if (m%n < m%lfmax) then
            !         vegetative phase
            m%fl = 1.0
            m%e = 1.0
            m%dn = m%rm * m%pt

            call lais(m%fl, m%di, m%pd, m%emp1, m%emp2, m%n, m%nb, m%swfac1, m%swfac2, m%pt, &
                    m%dn, m%p1, m%sla, m%dlai)
            m%dw = m%e * (m%pg) * m%pd
            m%dwc = m%fc * m%dw
            m%dwr = (1 - m%fc) * m%dw
            m%dwf = 0.0

        else
            !         reproductive phase
            m%fl = 2.0

            if (m%tmn >= m%tb .and. m%tmn <= 25) then
                m%di = (m%tmn - m%tb)
            else
                m%di = 0.0
            endif

            m%intc = m%intc + m%di
            m%e = 1.0
            call lais(m%fl, m%di, m%pd, m%emp1, m%emp2, m%n, m%nb, m%swfac1, m%swfac2, m%pt, &
                    m%dn, m%p1, m%sla, m%dlai)
            m%dw = m%e * (m%pg) * m%pd
            m%dwf = m%dw
            m%dwc = 0.0
            m%dwr = 0.0
            m%dn = 0.0
        endif
    end subroutine rate

    subroutine c_rate(m) bind(c, name='pm_rate')
        implicit none
        type(PlantModel) :: m
        call rate(m)
    end subroutine c_rate

    subroutine update(m, input)
        type(PlantModel) :: m
        type(PlantInput) :: input
        m%tmax = input%tmax
        m%tmin = input%tmin
        m%par  = input%par
        m%swfac1 = input%swfac1
        m%swfac2 = input%swfac2
    end subroutine update

    subroutine c_update(m, input) bind(c, name='pm_update')
        implicit none
        type(PlantModel) :: m
        type(PlantInput) :: input
        call update(m, input)
    end subroutine c_update

    subroutine integ(m)
        use iso_c_binding
        implicit none
        type(PlantModel) :: m

        m%lai = m%lai + m%dlai
        m%w = m%w + m%dw
        m%wc = m%wc + m%dwc
        m%wr = m%wr + m%dwr
        m%wf = m%wf + m%dwf

        m%lai = max(m%lai, 0.0)
        m%w = max(m%w, 0.0)
        m%wc = max(m%wc, 0.0)
        m%wr = max(m%wr, 0.0)
        m%wf = max(m%wf, 0.0)

        m%n = m%n + m%dn
        if (m%intc > m%intot) then
            m%endsim = 1
            return
        endif
    end subroutine integ

    subroutine c_integ(m) bind(c, name='pm_integ')
        implicit none
        type(PlantModel) :: m
        call integ(m)
    end subroutine c_integ
end module PlantComponent