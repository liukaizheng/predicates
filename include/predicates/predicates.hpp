#pragma once

#include <array>

namespace predicates {

double orient2d(const double* pa, const double* pb, const double* pc) noexcept;
double orient3d(const double* pa, const double* pb, const double* pc, const double* pd) noexcept;
double orient4d(const double* pa, const double* pb, const double* pc, const double* pd, const double* pe) noexcept;

std::array<double, 2> mi_orient0d(const double a, const double b) noexcept;
std::array<double, 2> mi_orient1d(const double* pa, const double* pb, const double* pc) noexcept;
std::array<double, 2> mi_orient2d(const double* pa, const double* pb, const double* pc, const double* pd) noexcept;
std::array<double, 2> mi_orient4d(const double* pa, const double* pb, const double* pc, const double* pd, const double* pe) noexcept;

}  // namespace predicates
