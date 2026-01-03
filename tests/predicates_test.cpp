#include "predicates/predicates.hpp"

#include <iostream>

using namespace predicates;

int main() {

    int failures = 0;
    const auto check = [&](bool ok, const char* expr) {
        if (!ok) {
            std::cerr << "FAILED: " << expr << '\n';
            ++failures;
        }
    };

  // check(predicates::is_even(0), "predicates::is_even(0)");
  // check(predicates::is_even(2), "predicates::is_even(2)");
  // check(!predicates::is_even(3), "!predicates::is_even(3)");

  // check(predicates::is_odd(1), "predicates::is_odd(1)");
  // check(predicates::is_odd(3), "predicates::is_odd(3)");
  // check(!predicates::is_odd(2), "!predicates::is_odd(2)");

    if (failures != 0) {
        std::cerr << failures << " test(s) failed.\n";
        return 1;
    }
    return 0;
}
