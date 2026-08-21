#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

inline void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

template <typename Actual, typename Expected, typename Tolerance>
inline void require_near(Actual actual, Expected expected, Tolerance tolerance,
                         const std::string &message) {
  if (std::abs(actual - expected) <= tolerance)
    return;

  std::ostringstream out;
  out.precision(17);
  out << message << ": actual=" << actual << ", expected=" << expected
      << ", tolerance=" << tolerance;
  throw std::runtime_error(out.str());
}
