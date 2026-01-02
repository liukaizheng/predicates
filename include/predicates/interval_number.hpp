#pragma once

#include <cfenv>
#include <x86/avx2.h>

namespace predicates {

inline void set_upward() noexcept {
    std::fesetround(FE_UPWARD);
}

inline void set_nearest() noexcept {
    std::fesetround(FE_TONEAREST);
}


class IntervalNumber {
public:
    explicit IntervalNumber(double value) noexcept
        : data_(simde_mm_set_pd(-value, value)) {}

    explicit IntervalNumber(simde__m128d data) noexcept
        : data_(data) {}

    IntervalNumber(const IntervalNumber&) noexcept = default;
    IntervalNumber(IntervalNumber&&) noexcept = default;
    IntervalNumber& operator=(const IntervalNumber&) noexcept = default;
    IntervalNumber& operator=(IntervalNumber&&) noexcept = default;
    ~IntervalNumber() = default;

    [[nodiscard]] IntervalNumber operator+(const IntervalNumber& other) const noexcept {
        return IntervalNumber(simde_mm_add_pd(data_, other.data_));
    }

    [[nodiscard]] IntervalNumber operator-(const IntervalNumber& other) const noexcept {
        return IntervalNumber(simde_mm_add_pd(data_, simde_mm_shuffle_pd(other.data_, other.data_, 1)));
    }

    [[nodiscard]] IntervalNumber operator*(const IntervalNumber& other) const noexcept {
        simde__m256d this_expanded = simde_mm256_castpd128_pd256(data_);
        this_expanded = simde_mm256_insertf128_pd(this_expanded, data_, 1);
        simde__m256d other_expanded = simde_mm256_castpd128_pd256(other.data_);
        other_expanded = simde_mm256_insertf128_pd(other_expanded, other.data_, 1);

        simde__m256d this_swizzled = simde_mm256_shuffle_pd(this_expanded, this_expanded, 5);
        simde__m256d other_negated_high = simde_mm256_xor_pd(other_expanded, simde_mm256_set_pd(-0.0, -0.0, 0.0, 0.0));
        simde__m256d other_negated_low = simde_mm256_xor_pd(other_expanded, simde_mm256_set_pd(0.0, 0.0, -0.0, -0.0));

        simde__m256d products_a = simde_mm256_mul_pd(this_expanded, other_negated_high);
        simde__m256d products_b = simde_mm256_mul_pd(this_swizzled, other_negated_low);

        simde__m256d max_products = simde_mm256_max_pd(products_a, products_b);
        simde__m256d max_shuffled = simde_mm256_shuffle_pd(max_products, max_products, 5);
        max_products = simde_mm256_max_pd(max_products, max_shuffled);

        simde__m256d result = simde_mm256_permute4x64_pd(max_products, 72);
        return IntervalNumber(simde_mm256_castpd256_pd128(result));
    }

private:
    simde__m128d data_;
};

}  // namespace predicates
