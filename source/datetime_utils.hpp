#include <iostream>
#include <stdexcept>

// Class representing a calendar date.
class Date
{
	private:
	// Data members to store the day, month, and year.
	int day,month,year;

	// Helper function to check if a given year is a leap year.
	bool isLeapYear(int y) const
	{
		return (y%4==0) && ((y%100!=0) || (y%400==0));
	}

	// Returns the number of days in a given month for a specific year.
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
	// Constructor to initialize the Date object with given day, month, and year.
	Date(int day,int month,int year) : day(day),month(month),year(year) {}

	// Increment the date by one day.
	// Rolls over the day, month, and year as necessary.
	void increment()
	{
		day++;
		// If the day exceeds the number of days in the current month, reset day and increment month.
		if (day>daysInMonth(month,year))
		{
			day=1;
			month++;
			// If the month exceeds 12, reset month and increment year.
			if (month>12)
			{
				month=1;
				year++;
			}
		}
	}

	void decrement()
	{
		if (day>1)
		{
			--day;
			return;
		}

		/* need to borrow from previous month */
		if (month==1)
		{
			month=12;
			--year;
		}
		else
		{
			--month;
		}
		day=daysInMonth(month,year);   // set to last day of new month
	}

	// Overloaded operator to compare if this Date is less than or equal to another Date.
	bool operator<=(const Date& other) const
	{
		if (year<other.year) return true;
		if (year==other.year && month<other.month) return true;
		if (year==other.year && month==other.month && day<=other.day) return true;
		return false;
	}

	// Overloaded operator to compare if this Date is strictly less than another Date.
	bool operator<(const Date& other) const
	{
		if (year<other.year) return true;
		if (year==other.year && month<other.month) return true;
		if (year==other.year && month==other.month && day<other.day) return true;
		return false;
	}

	// Overloaded equality operator to check if two Dates are exactly equal.
	bool operator==(const Date& other) const
	{
		if (year==other.year && month==other.month && day==other.day) return true;
		return false;
	}

	bool operator>(const Date& other) const {return other<*this;}

	bool operator>=(const Date& other) const {return !(*this<other);}

	bool operator!=(const Date& other) const {return !(*this==other);}

        // Overloaded stream insertion operator to print a Date in a formatted manner.
	friend std::ostream& operator<<(std::ostream& os,const Date& date)
	{
		os<<"Month: "<<date.month<<", Day: "<<date.day<<", Year: "<<date.year;
		return os;
	}

	// Getter for the year.
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
	// Compute the number of days between this Date and another Date.
	// Returns a positive number if this Date is later than the other,
	// negative if earlier, and zero if they are the same.
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

	// Shift the date by +-delta days and return a new Date.
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

	// Shift the date by +-delta months and return a new Date.
	// The day is clamped to the last day of the target month.
	Date addMonths(int delta) const
	{
		if (delta==0) return *this;         // fast path

		int y=year;
		int m=month;                       // 1 ... 12
		int d=day;

		int index=y*12+(m-1)+delta; // months since year 0
		int new_y=index/12;
		int new_m=index%12+1;           // back to 1 ... 12

		int max_day=daysInMonth(new_m,new_y); // uses fixed leap logic above
		if (d>max_day) d=max_day;

		return Date(d,new_m,new_y);
	}

	// Return the first day of the current month.
	Date beginOfMonth() const
	{
		return Date(1,month,year);          // month & year are members
	}

	private:
	// Helper function to calculate the number of days between two Dates.
	// It increments a temporary date from the 'other' date until it reaches this Date.
	int calculateDaysSince(const Date& other) const
	{
		int days=0;
		Date temp=other;	// Start from the other date.

		// Increment temp until it reaches the current date.
		while (temp<*this)
		{
			days++;
			temp.increment();
		}

		return days;
	}
};

// Class representing a range of dates with an iterator for traversal.
/*class DateRange
{
	private:
	// Start date, end date, and a current date for iteration.
	Date startDate,endDate,currentDate;

	public:
	// Constructor that initializes the date range with a start and end date.
	// The current date is initially set to the start date.
	DateRange(const Date& startDate,const Date& endDate) : startDate(startDate),endDate(endDate),currentDate(startDate) {}

        // Nested iterator class to allow iteration over the date range.
	struct DateIterator
	{
		// Overloaded inequality operator to compare iterators.
		// This implementation returns true as long as the current date is less than or equal to the end date.
		Date currentDate,endDate;
		DateIterator(const Date &currentDate,const Date &endDate) : currentDate(currentDate),endDate(endDate) {}

		bool operator!=(const DateIterator& other) const
		{
			return currentDate<=other.endDate;
		}

		// Dereference operator returns the current date.
		const Date& operator*() const
		{
			return currentDate;
		}

		// Pre-increment operator to move to the next date in the range.
		DateIterator& operator++()
		{
			currentDate.increment();
			return *this;
		}
	};

	// Returns an iterator pointing to the beginning of the date range.
	DateIterator begin() const
	{
		return DateIterator(startDate,endDate);
	}

	// Returns an iterator representing the end of the date range.
	// Note: This iterator uses endDate as both current and end.
	DateIterator end() const
	{
		return DateIterator(endDate,endDate);
	}
};*/

// ---------------------------------------------------------------------------
// Class representing a range of dates with an iterator for traversal.
// ---------------------------------------------------------------------------
class DateRange
{
	Date startDate,endDate;

public:
	// inclusive range [startDate … endDate]
	DateRange(const Date& s,const Date& e) : startDate(s),endDate(e) {}

	// -------------------------------------------------------------
	// iterator
	// -------------------------------------------------------------
	struct iterator
	{
		Date current,last;                       // 'last' is inclusive

		iterator(const Date& cur,const Date& last) : current(cur),last(last) {}

		const Date& operator*()  const {return current;}

		iterator& operator++()
		{                  // pre-increment
			current.increment();
			return *this;
		}

		bool operator!=(const iterator& other) const
		{
			return current!=other.current;      // compare current positions
		}
	};

	iterator begin() const {return iterator(startDate,endDate);}

	iterator end() const
	{
		Date onePastEnd=endDate;
		onePastEnd.increment();                   // sentinel iterator
		return iterator(onePastEnd,endDate);
	}

	// -------------------------------------------------------------
	// NEW: number of days in the range (inclusive)
	// -------------------------------------------------------------
	int length() const
	{
		return endDate.daysSince(startDate)+1;  // inclusive count
	}
};
