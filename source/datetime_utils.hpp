#include <iostream>
#include <stdexcept>

class Date
{
	private:
	int day,month,year;

	bool isLeapYear(int year) const
	{
		return (year%4==0);
	}

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

	public:
	Date(int day,int month,int year) : day(day),month(month),year(year) {}

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
		if (year==other.year && month==other.month && day==other.day) return true;
		return false;
	}

	friend std::ostream& operator<<(std::ostream& os,const Date& date)
	{
		os<<"Month: "<<date.month<<", Day: "<<date.day<<", Year: "<<date.year;
		return os;
	}

	int getYear() const
	{
		return year;
	}

	int getMonth() const
	{
		return month;
	}

	int getDay() const
	{
		return day;
	}

	int daysSince(const Date& other) const
	{
		if (*this==other)
			return 0;
		else if (*this<other)
			return -other.calculateDaysSince(*this);
		else
			return calculateDaysSince(other);
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

	private:
	int calculateDaysSince(const Date& other) const
	{
		int days=0;
		Date temp=other;

		while (temp<*this)
		{
			days++;
			temp.increment();
		}

		return days;
	}
};

class DateRange
{
	private:
	Date startDate,endDate,currentDate;

	public:
	DateRange(const Date& startDate,const Date& endDate) : startDate(startDate),endDate(endDate),currentDate(startDate) {}

	struct DateIterator
	{
		Date currentDate,endDate;
		DateIterator(const Date &currentDate,const Date &endDate) : currentDate(currentDate),endDate(endDate) {}

		bool operator!=(const DateIterator& other) const
		{
			return currentDate<=other.endDate;
		}

		const Date& operator*() const
		{
			return currentDate;
		}

		DateIterator& operator++()
		{
			currentDate.increment();
			return *this;
		}
	};

	DateIterator begin() const
	{
		return DateIterator(startDate,endDate);
	}

	DateIterator end() const
	{
		return DateIterator(endDate,endDate);
	}
};
