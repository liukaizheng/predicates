#include "predicates/predicates.hpp"

#include <iostream>

int main() {
  std::cout << "2 is even: " << predicates::is_even(2) << '\n';
  std::cout << "3 is odd: " << predicates::is_odd(3) << '\n';
  return 0;
}

