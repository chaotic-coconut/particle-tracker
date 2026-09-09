// Round-trip and boundary tests for the date and argument helpers.
//
// These functions are on the hot path of every run: parseDate turns a
// command-line argument into a Date, secondsSince2000 turns that Date into
// simulation time, and beginOfMonth is what the month cache keys on. The dates
// exercised here span the range for which velocity fields actually exist
// (roughly 1992 onwards) rather than exotic calendar cases that cannot occur.

#include "test_helpers.hpp"

#include <particle_tracker/data_prep_utils.hpp>

#include <string>
#include <vector>

namespace {

using particle_tracker::dateFromSecondsSince2000;
using particle_tracker::formatDate;
using particle_tracker::parseDate;
using particle_tracker::secondsSince2000;
using particle_tracker::splitArgs;

// First and last day of every month over the data era, plus every leap day.
std::vector<Date> interestingDates() {
  std::vector<Date> dates;
  for (int year = 1992; year <= 2030; ++year) {
    for (int month = 1; month <= 12; ++month) {
      const Date first(1, month, year);
      dates.push_back(first);
      dates.push_back(first.addMonths(1).addDays(-1)); // last day of the month
    }
  }
  return dates;
}

std::string describe(const Date &d) { return formatDate(d); }

void check_format_shape() {
  require(formatDate(Date(1, 1, 2000)) == "2000_01_01",
          "formatDate zero-pads day and month");
  require(formatDate(Date(31, 12, 1999)) == "1999_12_31",
          "formatDate renders a two-digit day and month");
}

void check_parse_accepts_both_separators() {
  const Date dashed = parseDate("2003-07-09");
  const Date underscored = parseDate("2003_07_09");
  require(dashed == Date(9, 7, 2003), "parseDate reads YYYY-MM-DD");
  require(dashed == underscored,
          "parseDate treats '-' and '_' as the same separator");
}

void check_parse_format_round_trip() {
  for (const Date &d : interestingDates())
    require(parseDate(formatDate(d)) == d,
            "parseDate(formatDate(d)) == d for " + describe(d));
}

void check_seconds_round_trip() {
  for (const Date &d : interestingDates())
    require(dateFromSecondsSince2000(secondsSince2000(d)) == d,
            "seconds round trip for " + describe(d));
}

void check_seconds_are_ordered() {
  // Includes dates before the epoch: the data era starts in the 1990s, so the
  // negative branch of the conversion is reached in every real run.
  const auto dates = interestingDates();
  for (std::size_t i = 1; i < dates.size(); ++i)
    require(secondsSince2000(dates[i - 1]) < secondsSince2000(dates[i]),
            "secondsSince2000 is strictly increasing at " + describe(dates[i]));

  require(secondsSince2000(Date(31, 12, 1999)) < 0,
          "secondsSince2000 is negative before the epoch");
  require(secondsSince2000(Date(1, 1, 2000)) == 0,
          "secondsSince2000 is zero at the epoch");
}

void check_begin_of_month() {
  for (const Date &d : interestingDates()) {
    const Date begin = d.beginOfMonth();
    require(begin.getDay() == 1, "beginOfMonth returns day 1 for " + describe(d));
    require(begin.getMonth() == d.getMonth() && begin.getYear() == d.getYear(),
            "beginOfMonth preserves month and year for " + describe(d));
    require(begin.beginOfMonth() == begin,
            "beginOfMonth is idempotent for " + describe(d));
  }
  require(Date(29, 2, 2000).beginOfMonth() == Date(1, 2, 2000),
          "beginOfMonth of a leap day");
}

void check_split_args() {
  const auto three = splitArgs("lon,lat,water_u");
  require(three.size() == 3 && three[0] == "lon" && three[2] == "water_u",
          "splitArgs splits on commas");
  require(splitArgs("lon").size() == 1, "splitArgs with a single element");
  require(splitArgs("").empty(), "splitArgs of an empty string is empty");
  require(splitArgs(",,").empty(), "splitArgs drops empty fields");
  require(splitArgs("lon,,lat").size() == 2,
          "splitArgs skips a repeated separator");
  require(splitArgs("lon,lat,").size() == 2,
          "splitArgs ignores a trailing separator");
}

} // namespace

int main() {
  check_format_shape();
  check_parse_accepts_both_separators();
  check_parse_format_round_trip();
  check_seconds_round_trip();
  check_seconds_are_ordered();
  check_begin_of_month();
  check_split_args();
}
