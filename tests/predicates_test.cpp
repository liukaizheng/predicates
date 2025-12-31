#include "predicates/predicates.hpp"

#include <simde-common.h>
#include <x86/avx2.h>

#include <iostream>
#include "predicates/interval_number.hpp"

int main() {
    IntervalNumber a(2.0);
    IntervalNumber b(-3.0);
    auto c = a * b;
    std::cout << c.data_[0] << " " << c.data_[1] << std::endl;
    int failures = 0;
  const auto check = [&](bool ok, const char* expr) {
    if (!ok) {
      std::cerr << "FAILED: " << expr << '\n';
      ++failures;
    }
  };

  check(predicates::is_even(0), "predicates::is_even(0)");
  check(predicates::is_even(2), "predicates::is_even(2)");
  check(!predicates::is_even(3), "!predicates::is_even(3)");

  check(predicates::is_odd(1), "predicates::is_odd(1)");
  check(predicates::is_odd(3), "predicates::is_odd(3)");
  check(!predicates::is_odd(2), "!predicates::is_odd(2)");

  if (failures != 0) {
    std::cerr << failures << " test(s) failed.\n";
    return 1;
  }
  return 0;
}
