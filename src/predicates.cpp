#include "predicates/predicates.hpp"
#include <predicates/natural_number.hpp>
#include <predicates/float_number.hpp>

namespace predicates {

bool is_even(int value) noexcept {
  return value % 2 == 0;
}

bool is_odd(int value) noexcept {
  return value % 2 != 0;
}

}  // namespace predicates
