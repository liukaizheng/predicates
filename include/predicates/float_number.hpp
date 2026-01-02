#pragma once

#include "natural_number.hpp"
#include <cstdint>
#include <cstring>
#include <limits>
#include <bit>
#include <string>

namespace predicates {

// A FloatNumber is a floating point number with arbitrarily large mantissa.
// A FloatNumber f evaluates to: f = sign * mantissa * 2^exponent
//
// The mantissa is a NaturalNumber whose least significant bit is 1 (normalized).
// Number is zero if mantissa is empty (sign == 0).

class FloatNumber {
private:
    NaturalNumber mantissa_;
    int32_t exponent_;
    int32_t sign_;  // -1, 0, or 1

    // Right-shift mantissa as long as the least significant bit is zero
    void normalize() {
        if (mantissa_.is_zero()) {
            sign_ = 0;
            exponent_ = 0;
            return;
        }

        // Count and remove trailing zeros
        uint32_t trailing_zeros = mantissa_.count_trailing_zeros();
        if (trailing_zeros > 0) {
            mantissa_ >>= trailing_zeros;
            exponent_ += static_cast<int32_t>(trailing_zeros);
        }

        // Final check
        if (mantissa_.is_zero()) {
            sign_ = 0;
            exponent_ = 0;
        }
    }

    // Left-shift the mantissa by n bits and reduce the exponent accordingly
    void left_shift(uint32_t n) {
        if (n == 0 || mantissa_.is_zero()) {
            return;
        }
        mantissa_ <<= n;
        exponent_ -= static_cast<int32_t>(n);
    }

public:
    // Default constructor creates zero
    FloatNumber() : exponent_(0), sign_(0) {}

    // Construct from double (lossless)
    explicit FloatNumber(double d) {
        sign_ = (d > 0) - (d < 0);

        if (sign_ == 0) {
            exponent_ = 0;
            return;
        }

        // Extract bits from double
        const auto bits = std::bit_cast<uint64_t>(d);

        // Extract mantissa (52 bits) and add implicit bit
        uint64_t mantissa_bits = (bits & 0x000FFFFFFFFFFFFFull) | 0x0010000000000000ull;

        // Extract exponent
        int32_t exp = static_cast<int32_t>((bits >> 52) & 0x7FF);
        exponent_ = exp - 1023 - 52;  // Adjust for bias and mantissa position

        // Create mantissa from bits
        mantissa_ = NaturalNumber(mantissa_bits);
    }

    // Construct from components
    FloatNumber(const NaturalNumber& m, int32_t e, int32_t s)
        : mantissa_(m), exponent_(e), sign_(s) {
        if (mantissa_.is_zero()) {
            sign_ = 0;
            exponent_ = 0;
        } else {
            normalize();
        }
    }

    // Copy and move constructors/assignments
    FloatNumber(const FloatNumber&) = default;
    FloatNumber(FloatNumber&&) noexcept = default;
    FloatNumber& operator=(const FloatNumber&) = default;
    FloatNumber& operator=(FloatNumber&&) noexcept = default;
    ~FloatNumber() = default;

    // Get sign (-1, 0, or 1)
    [[nodiscard]] int32_t sgn() const noexcept { return sign_; }

    // Check if zero
    [[nodiscard]] bool is_zero() const noexcept { return sign_ == 0; }

    // Get mantissa
    [[nodiscard]] const NaturalNumber& get_mantissa() const noexcept { return mantissa_; }

    // Get exponent
    [[nodiscard]] int32_t get_exponent() const noexcept { return exponent_; }

    // Negate in place
    void negate() { sign_ = -sign_; }

    // Return negated copy
    [[nodiscard]] FloatNumber operator-() const {
        FloatNumber result = *this;
        result.negate();
        return result;
    }

    // Truncated base-2 logarithm
    [[nodiscard]] int32_t log2() const {
        return exponent_ + static_cast<int32_t>(mantissa_.num_bits()) - 1;
    }

    // Convert to double (truncated)
    [[nodiscard]] double to_double() const {
        if (sign_ == 0) {
            return 0.0;
        }
        if (mantissa_.is_zero()) {
            return 0.0;
        }

        // Get number of significant bits
        auto num_bits = static_cast<int32_t>(mantissa_.num_bits());
        int32_t exp = exponent_ + num_bits - 1;

        // Check for overflow/underflow
        if (exp > 1023) {
            return sign_ < 0 ? -std::numeric_limits<double>::infinity()
                             : std::numeric_limits<double>::infinity();
        }
        if (exp < -1022) {
            return 0.0;
        }

        // Extract top 53 bits for mantissa
        uint64_t mantissa_bits = 0;
        NaturalNumber temp = mantissa_;

        // Shift to get exactly 53 bits (52 + implicit)
        if (num_bits > 53) {
            temp >>= static_cast<uint32_t>(num_bits - 53);
        } else if (num_bits < 53) {
            temp <<= static_cast<uint32_t>(53 - num_bits);
        }

        // Convert to uint64
        temp.to_uint64(mantissa_bits);

        // Remove implicit bit
        mantissa_bits &= 0x000FFFFFFFFFFFFFull;

        // Build double
        uint64_t result_bits = 0;
        if (sign_ < 0) result_bits |= 0x8000000000000000ull;
        result_bits |= static_cast<uint64_t>(exp + 1023) << 52;
        result_bits |= mantissa_bits;
        return std::bit_cast<double>(result_bits);
    }

    // Addition
    [[nodiscard]] FloatNumber operator+(const FloatNumber& other) const {
        if (is_zero()) return other;
        if (other.is_zero()) return *this;

        if (exponent_ == other.exponent_) {
            FloatNumber result;

            if (sign_ == other.sign_) {
                result.mantissa_ = mantissa_ + other.mantissa_;
                result.sign_ = sign_;
            } else if (other.mantissa_ >= mantissa_) {
                result.mantissa_ = other.mantissa_ - mantissa_;
                result.sign_ = other.sign_;
            } else {
                result.mantissa_ = mantissa_ - other.mantissa_;
                result.sign_ = sign_;
            }

            result.exponent_ = exponent_;
            result.normalize();
            return result;
        } else if (exponent_ > other.exponent_) {
            FloatNumber adjusted = *this;
            adjusted.left_shift(static_cast<uint32_t>(exponent_ - other.exponent_));
            return adjusted + other;
        } else {
            FloatNumber adjusted = other;
            adjusted.left_shift(static_cast<uint32_t>(other.exponent_ - exponent_));
            return adjusted + *this;
        }
    }

    // Subtraction
    [[nodiscard]] FloatNumber operator-(const FloatNumber& other) const {
        if (is_zero()) return -other;
        if (other.is_zero()) return *this;

        if (exponent_ == other.exponent_) {
            FloatNumber result;

            if (sign_ != other.sign_) {
                result.mantissa_ = mantissa_ + other.mantissa_;
                result.sign_ = sign_;
            } else if (other.mantissa_ >= mantissa_) {
                result.mantissa_ = other.mantissa_ - mantissa_;
                result.sign_ = -sign_;
            } else {
                result.mantissa_ = mantissa_ - other.mantissa_;
                result.sign_ = sign_;
            }

            result.exponent_ = exponent_;
            result.normalize();
            return result;
        } else if (exponent_ > other.exponent_) {
            FloatNumber adjusted = *this;
            adjusted.left_shift(static_cast<uint32_t>(exponent_ - other.exponent_));
            return adjusted - other;
        } else {
            FloatNumber adjusted = other;
            adjusted.left_shift(static_cast<uint32_t>(other.exponent_ - exponent_));
            return *this - adjusted;
        }
    }

    // Multiplication
    [[nodiscard]] FloatNumber operator*(const FloatNumber& other) const {
        if (is_zero() || other.is_zero()) return FloatNumber{};

        FloatNumber result;
        result.mantissa_ = mantissa_ * other.mantissa_;
        result.sign_ = sign_ * other.sign_;
        result.exponent_ = exponent_ + other.exponent_;
        result.normalize();
        return result;
    }

    // Compound assignment operators
    FloatNumber& operator+=(const FloatNumber& other) {
        *this = *this + other;
        return *this;
    }

    FloatNumber& operator-=(const FloatNumber& other) {
        *this = *this - other;
        return *this;
    }

    FloatNumber& operator*=(const FloatNumber& other) {
        *this = *this * other;
        return *this;
    }

    // Comparison operators
    [[nodiscard]] bool operator==(const FloatNumber& other) const {
        return (*this - other).sign_ == 0;
    }

    [[nodiscard]] bool operator!=(const FloatNumber& other) const {
        return (*this - other).sign_ != 0;
    }

    // String representation
    [[nodiscard]] std::string to_string() const {
        if (sign_ == 0) return "0";

        std::string s;
        if (sign_ < 0) s += "-";
        s += mantissa_.to_decimal_string();
        s += " * 2^";
        s += std::to_string(exponent_);
        return s;
    }
};

// Free function for sign
[[nodiscard]] inline int32_t sgn(const FloatNumber& f) {
    return f.sgn();
}

// Operators with left-hand double
[[nodiscard]] inline FloatNumber operator+(double lhs, const FloatNumber& rhs) {
    return FloatNumber(lhs) + rhs;
}

[[nodiscard]] inline FloatNumber operator-(double lhs, const FloatNumber& rhs) {
    return FloatNumber(lhs) - rhs;
}

[[nodiscard]] inline FloatNumber operator*(double lhs, const FloatNumber& rhs) {
    return FloatNumber(lhs) * rhs;
}

} // namespace predicates
