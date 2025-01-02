module m_daytim

    !! Date and Time routines
    !!
    !! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use m_std
    use iso_fortran_env, only: error_unit

    implicit none
    private
    save

    public :: daytim__isLeapYear         ! Check Leap Year
    public :: daytim__ymd2jul            ! Convert Year/Month/Day to Day of Year (Julian Day)
    public :: daytim__jul2md             ! Convert Day of Year (Julian Day) to Year/Month/Day
    public :: daytim__dayweek            ! Day of the week
    public :: daytim__timelocal          ! Convert Date and Time to Seconds form the reference time
    public :: daytim__localtime          ! Convert Seconds form the reference time to Date and Time
    public :: daytim__getDate            ! Obtain formatted date and time

    interface daytim__getdate

        !! Obtain current date and time, either formatted or integer form

        module procedure getdate_c, getdate_i1, getdate_i2

    end interface daytim__getdate

    interface daytim__dayweek

        !! Obtain day of the week

        module procedure dayweek_i, dayweek_a

    end interface daytim__dayweek

contains

    subroutine daytim__isLeapYear(year, isLeap)

        !! Return logical value whether if the given year is leap year

        integer, intent(in)  :: year
        logical, intent(out) :: isLeap

        if (mod(year, 100) == 0 .and. mod(year, 400) == 0) then
            isLeap = .true.
        else if (mod(year, 100) /= 0 .and. mod(year, 4) == 0) then
            isLeap = .true.
        else
            isLeap = .false.
        end if

    end subroutine daytim__isLeapYear


    subroutine daytim__ymd2jul(year, month, day, julday)

        !! Convert Year/Month/Day to Day of Year (Julian Day)

        integer, intent(in)  :: year
        integer, intent(in)  :: month
        integer, intent(in)  :: day
        integer, intent(out) :: julday

        integer :: dom(12)  ! day of month
        logical :: is_leap  ! leap year ?
        integer :: i

        call daytim__isLeapYear(year, is_leap)

        if (is_leap) then
            dom = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
        else
            dom = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
        end if

        !! initialize
        julday = 0

        !! error handling
        if (day > dom(month)) then
            write (error_unit, *) 'ymd2jul: input day excees the last day of the month'
            julday = 0
            return
        end if
        if (month > 12) then
            write (error_unit, *) 'ymd2jul: month excees 12 !'
            julday = 0
            return
        end if

        !! month
        do i = 1, month - 1
            julday = julday + dom(i)
        end do

        !! day
        julday = julday + day

    end subroutine daytim__ymd2jul


    subroutine daytim__jul2md(julday, year, month, day)

        !! Convert Day of Year (Julian Day) to Year/Month/Day

        integer, intent(in)  :: julday
        integer, intent(in)  :: year
        integer, intent(out) :: month
        integer, intent(out) :: day

        integer :: dom(12)   !! day of month
        logical :: is_leap   !! check for uruu doshi

        call daytim__isLeapYear(year, is_leap)

        if (is_leap) then
            dom = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
        else
            dom = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
        end if

        if (is_leap .and. julday > 366) then
            write (error_unit, *) 'jul2md: day of year excees Dec. 31'
            month = 0
            day = 0
            return
        else if (.not. is_leap) then
            if (julday > 365) then
                write (error_unit, *) 'jul2md: day of year excees Dec. 31'
                month = 0
                day = 0
                return
            end if

        end if

        month = 1

        do while (sum(dom(1:month)) < julday)
            month = month + 1
        end do

        day = julday - sum(dom(1:month - 1))

    end subroutine daytim__jul2md


    subroutine dayweek_i(year, month, day, dw)

        !! Returns day of the week from date information
        !! Return dw (integer): 0 (Sunday) to 6 (Satureday)

        integer, intent(in)  :: year
        integer, intent(in)  :: month
        integer, intent(in)  :: day
        integer, intent(out) :: dw

        integer :: y, m, d

        if (month <= 0 .or. month >= 14) then
            write (error_unit, '(A)') 'subroutine dayweek: invalied argument'
            dw = -1
            return
        end if
        if (month == 1 .or. month == 2) then
            y = year - 1
            m = month + 12
            d = day
        else
            y = year
            m = month
            d = day
        end if

        dw = y + int(y / 4.) - int(y / 100.) + int(y / 400.) + &
             int((26 * m + 16) / 10.) + d
        dw = mod(dw, 7)

    end subroutine dayweek_i


    subroutine dayweek_a(year, month, day, dw)

        !! Returns character formatted day of the week from the given date

        integer, intent(in)  :: year
        integer, intent(in)  :: month
        integer, intent(in)  :: day
        character(*), intent(out) :: dw

        integer :: idw
        character(9) :: dwname(0:6) = (/'Sunday   ', &
                                        'Monday   ', &
                                        'Tuesday  ', &
                                        'Wednesday', &
                                        'Thursday ', &
                                        'Friday   ', &
                                        'Saturday '/)

        call dayweek_i(year, month, day, idw)
        dw = dwname(idw)

    end subroutine dayweek_a


    subroutine daytim__timelocal(year, month, day, hour, min, sec, tim)

        !! Return the elapsed seconds from 1970-01-01 00:00:00 (UNIX time/POSIX time) from given date/time.
        !! Almost compatible with timelocal function in perl, but the month takes 1-12, not 0-11 as perl.

        integer, intent(in)  :: year
        integer, intent(in)  :: month
        integer, intent(in)  :: day
        integer, intent(in)  :: hour
        integer, intent(in)  :: min
        integer, intent(in)  :: sec
        integer, intent(out) :: tim

        integer :: doy
        integer :: jday
        integer :: i

        tim = 0
        if (year > 1970) then
            do i = 1970, year - 1
                call daytim__ymd2jul(i, 12, 31, doy)
                tim = tim + doy * 24 * 60 * 60
            end do
        else
            do i = year, 1970 - 1
                call daytim__ymd2jul(i, 12, 31, doy)
                tim = tim - doy * 24 * 60 * 60
            end do
        end if
        call daytim__ymd2jul(year, month, day, jday)
        tim = tim + ((jday - 1) * 24 * 60 + hour * 60 + min) * 60 + sec

    end subroutine daytim__timelocal


    subroutine daytim__localtime(tim, year, month, day, hour, min, sec)

        !! Inversely convert UNIX/POSIX time (seconds from 1970-01-01 00:00:00) to date and time
        !! This is an inverse routine of daytim__timelocal

        integer, intent(in)  :: tim
        integer, intent(out) :: year
        integer, intent(out) :: month
        integer, intent(out) :: day
        integer, intent(out) :: hour
        integer, intent(out) :: min
        integer, intent(out) :: sec

        integer :: doy, soy, ttim, jday

        ttim = tim

        if (ttim >= 0) then
            year = 1970
            do
                call daytim__ymd2jul(year, 12, 31, doy)
                soy = doy * 24 * 60 * 60
                if (ttim < soy) then
                    exit
                else
                    year = year + 1
                    ttim = ttim - soy
                end if

            end do

            jday = ttim / 86400 + 1
            call daytim__jul2md(jday, year, month, day)

            ttim = ttim - (jday - 1) * 86400
            hour = ttim / 3600
            ttim = ttim - hour * 3600
            min = ttim / 60
            sec = ttim - min * 60

        else
            year = 1969
            do
                call daytim__ymd2jul(year, 12, 31, doy)
                soy = doy * 24 * 60 * 60
                ttim = ttim + soy
                if (ttim >= 0) then
                    exit
                else
                    year = year - 1
                end if

            end do

            jday = ttim / 86400 + 1
            call daytim__jul2md(jday, year, month, day)

            ttim = ttim - (jday - 1) * 86400
            hour = ttim / 3600
            ttim = ttim - hour * 3600
            min = ttim / 60
            sec = ttim - min * 60

        end if

    end subroutine daytim__localtime


    subroutine getDate_c(date)

        !! Return date and time in formatted characters
        !! Example: "2005/07/31 15:32:53"

        character(20), intent(out) :: date

        character(8)  :: ymd
        character(10) :: hms

        call date_and_time(ymd, hms)
        date = ymd(1:4)//'/'//ymd(5:6)//'/'//ymd(7:8)//' '//hms(1:2)//':'//hms(3:4)//':'//hms(5:6)

    end subroutine getDate_c


    subroutine getDate_i1(yr, mo, dy, hr, mi, sc)

        !! Return date and time as intergers in the following order:

        integer, intent(out) :: yr, mo, dy, hr, mi, sc

        character(8)  :: ymd
        character(10) :: hms

        call date_and_time(ymd, hms)

        read (ymd(1:4), *) yr
        read (ymd(5:6), *) mo
        read (ymd(7:8), *) dy
        read (hms(1:2), *) hr
        read (hms(3:4), *) mi
        read (hms(5:6), *) sc

    end subroutine getDate_i1


    subroutine getDate_i2(itim)

        !! Return date and time as intergers in the localtime integer

        integer, intent(out) :: itim

        integer :: yr, mo, dy, hr, mi, sc
        character(8)  :: ymd
        character(10) :: hms

        call date_and_time(ymd, hms)
        read (ymd(1:4), *) yr
        read (ymd(5:6), *) mo
        read (ymd(7:8), *) dy
        read (hms(1:2), *) hr
        read (hms(3:4), *) mi
        read (hms(5:6), *) sc

        call daytim__timelocal(yr, mo, dy, hr, mi, sc, itim)

    end subroutine getDate_i2

    
end module m_daytim
