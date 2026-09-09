#pragma once

#include <iostream>
#include <stdexcept>
#include <cstdint>

/**
 * @file datetime_utils.hpp
 * @brief Small date and date-range utilities (civil calendar, Gregorian rules).
 *
 * This header provides:
 * - Date      : a simple civil date (day, month, year) with basic arithmetic.
 * - DateRange : an inclusive range [start..end] with an iterator and length().
 *
 * Notes / assumptions
 * -------------------
 * - Uses proleptic Gregorian calendar rules (leap years: divisible by 4, except centuries
 *   not divisible by 400).
 * - No timezone logic. This is pure calendar arithmetic.
 * - daysSince/hoursSince/secondsSince are computed via an O(1) day-number conversion
 *   (not by looping day-by-day).
 */

/**
 * @brief Represents a calendar date (day/month/year) in the Gregorian calendar.
 *
 * This is intentionally minimal:
 * - No time-of-day.
 * - No timezone or locale.
 *
 * Arithmetic:
 * - increment()/decrement(): +/- 1 day with proper month/year rollover.
 * - addDays(delta): returns a new Date shifted by delta days (uses repeated +/-1 steps).
 *   This is O(delta). A closed-form version would be faster for large deltas, but
 *   the hot paths are secondsSince() / daysSince(), which are O(1).
 * - addMonths(delta): shifts month, clamps day to last valid day of target month.
 *
 * Comparison:
 * - Lexicographic by (year, month, day).
 *
 * Safety:
 * - Constructor validates month/day ranges and throws if invalid.
 */
class Date
{
    private:
    int day,month,year;

    // Helper: leap year test (Gregorian).
    bool isLeapYear(int y) const
    {
        return (y%4==0) && ((y%100!=0) || (y%400==0));
    }

    // Helper: days in month (1..12) for given year.
    int daysInMonth(int month,int year) const
    {
        switch(month)
        {
            case 1: case 3: case 5: case 7: case 8: case 10: case 12:
                return 31;
            case 4: case 6: case 9: case 11:
                return 30;
            case 2:
                return isLeapYear(year) ? 29 : 28;
            default:
                throw std::runtime_error("Invalid month");
        }
    }

    // Helper: convert civil date to a monotonically increasing day number.
    //
    // This is an O(1) algorithm (Howard Hinnant / "civil_from_days" family).
    // It works for a proleptic Gregorian calendar and avoids looping by days.
    //
    // Returns days since 1970-01-01 (can be negative for earlier dates).
    static std::int64_t daysFromCivil(int y,int m,int d)
    {
        // Shift March to the start of the year to make leap handling simple.
        y -= (m <= 2);
        const int era = (y >= 0 ? y : y-399) / 400;
        const unsigned yoe = static_cast<unsigned>(y - era * 400);              // [0, 399]
        const unsigned doy = (153u*(static_cast<unsigned>(m + (m > 2 ? -3 : 9))) + 2u)/5u
                           + static_cast<unsigned>(d) - 1u;                     // [0, 365]
        const unsigned doe = yoe*365u + yoe/4u - yoe/100u + doy;                 // [0, 146096]
        return static_cast<std::int64_t>(era)*146097 + static_cast<std::int64_t>(doe) - 719468;
    }

    public:
    /**
     * @brief Construct a Date.
     * @throws std::runtime_error if month/day are invalid.
     */
    Date(int day,int month,int year) : day(day),month(month),year(year)
    {
        if (month<1 || month>12)
            throw std::runtime_error("Date: month out of range");

        const int dim=daysInMonth(month,year);
        if (day<1 || day>dim)
            throw std::runtime_error("Date: day out of range for month/year");
    }

    // Increment the date by one day (with rollover).
    void increment()
    {
        day++;
        if (day>daysInMonth(month,year))
        {
            day=1;
            month++;
            if (month>12)
            {
                month=1;
                year++;
            }
        }
    }

    // Decrement the date by one day (with rollover).
    void decrement()
    {
        if (day>1)
        {
            --day;
            return;
        }

        if (month==1)
        {
            month=12;
            --year;
        }
        else
        {
            --month;
        }

        day=daysInMonth(month,year);
    }

    // Comparisons (lexicographic by year, month, day).
    bool operator<=(const Date& other) const
    {
        if (year<other.year) return true;
        if (year==other.year && month<other.month) return true;
        if (year==other.year && month==other.month && day<=other.day) return true;
        return false;
    }

    bool operator<(const Date& other) const
    {
        if (year<other.year) return true;
        if (year==other.year && month<other.month) return true;
        if (year==other.year && month==other.month && day<other.day) return true;
        return false;
    }

    bool operator==(const Date& other) const
    {
        return (year==other.year && month==other.month && day==other.day);
    }

    bool operator>(const Date& other) const {return other<*this;}

    bool operator>=(const Date& other) const {return !(*this<other);}

    bool operator!=(const Date& other) const {return !(*this==other);}

    // Stream output (debug-friendly).
    friend std::ostream& operator<<(std::ostream& os,const Date& date)
    {
        os<<"Month: "<<date.month<<", Day: "<<date.day<<", Year: "<<date.year;
        return os;
    }

    int getYear() const {return year;}

    int getMonth() const {return month;}

    int getDay() const {return day;}

    /**
     * @brief Days difference: (*this) - other in days.
     *
     * Returns:
     * - positive if this date is later than other
     * - negative if earlier
     * - zero if equal
     *
     * Complexity: O(1).
     */
    int daysSince(const Date& other) const
    {
        const std::int64_t a=daysFromCivil(year,month,day);
        const std::int64_t b=daysFromCivil(other.year,other.month,other.day);

        // Cast to int: safe for typical simulation ranges (years ~ 1900..2100).
        // If you plan to handle centuries-wide spans, switch return type to int64_t.
        return static_cast<int>(a-b);
    }

    int hoursSince(const Date& other) const
    {
        int days=daysSince(other);
        return days*24;
    }

    int secondsSince(const Date& other) const
    {
        int days=daysSince(other);
        return days*24*3600;
    }

    /**
     * @brief Return a new Date shifted by delta days.
     *
     * Note: This uses repeated increment()/decrement(). It is simple and safe.
     */
    Date addDays(int delta) const
    {
        if (delta==0) return *this;

        Date tmp=*this;
        if (delta>0)
        {
            for (int i=0;i<delta;i++) tmp.increment();
        }
        else
        {
            for (int i=0;i<-delta;i++) tmp.decrement();
        }
        return tmp;
    }

    /**
     * @brief Return a new Date shifted by delta months.
     *
     * The day is clamped to the last valid day in the resulting month.
     * Example: 31 Jan + 1 month -> 28 Feb (or 29 in leap years).
     */
    Date addMonths(int delta) const
    {
        if (delta==0) return *this;

        int y=year;
        int m=month;
        int d=day;

        int index=y*12+(m-1)+delta;
        int new_y=index/12;
        int new_m=index%12+1;

        int max_day=daysInMonth(new_m,new_y);
        if (d>max_day) d=max_day;

        return Date(d,new_m,new_y);
    }

    // First day of current month.
    Date beginOfMonth() const
    {
        return Date(1,month,year);
    }
};

/**
 * @brief Inclusive date range [startDate .. endDate] with iterator support.
 *
 * Example:
 *   for (const Date& d : DateRange(Date(1,1,2000), Date(10,1,2000))) { ... }
 */
class DateRange
{
    Date startDate,endDate;

public:
    // inclusive range [startDate ... endDate]
    DateRange(const Date& s,const Date& e) : startDate(s),endDate(e)
    {
        if (e<s)
            throw std::runtime_error("DateRange: endDate is before startDate");
    }

    struct iterator
    {
        Date current,last;                       // 'last' is inclusive

        iterator(const Date& cur,const Date& last) : current(cur),last(last) {}

        const Date& operator*()  const {return current;}

        iterator& operator++()
        {
            current.increment();
            return *this;
        }

        bool operator!=(const iterator& other) const
        {
            return current!=other.current;
        }
    };

    iterator begin() const {return iterator(startDate,endDate);}

    iterator end() const
    {
        Date onePastEnd=endDate;
        onePastEnd.increment();
        return iterator(onePastEnd,endDate);
    }

    // number of days in the range (inclusive)
    int length() const
    {
        return endDate.daysSince(startDate)+1;
    }
};
