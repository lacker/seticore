#include "catch/catch.hpp"

#include "fil_reader.h"

TEST_CASE("converting from sigproc ra", "[fil]") {
  double seconds_per_hour = 60 * 60;
  double sigproc = 123000.5;
  double hours = convertFromSigprocRaOrDec(sigproc);
  REQUIRE(hours == Approx(12.5 + 0.5 / seconds_per_hour));

  hours = convertFromSigprocRaOrDec(-sigproc);
  REQUIRE(hours == Approx(-12.5 - 0.5 / seconds_per_hour));
}
