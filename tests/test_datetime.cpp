#include "test_helpers.hpp"

#include <particle_tracker/detail/datetime_utils.hpp>

int main() {
  const Date epoch(1, 1, 2000);

  require(epoch.addDays(0) == epoch, "Date addDays(0)");
  require(epoch.addDays(3653) == Date(1, 1, 2010),
          "Date conversion across a decade");
  require(epoch.addDays(-1095) == Date(1, 1, 1997),
          "Negative date conversion");
  require(Date(28, 2, 2000).addDays(1) == Date(29, 2, 2000),
          "Leap-day addition");
  require(Date(1, 3, 1900).addDays(-1) == Date(28, 2, 1900),
          "Non-leap century subtraction");

  const DateRange range(Date(28, 2, 2000), Date(1, 3, 2000));
  require(range.length() == 3, "Inclusive leap-day range length");
}
