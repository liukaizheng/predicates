#pragma once

#include <algorithm>
#include <compare>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <compare>
#include <vector>
#include <bit>

namespace predicates {

// A NaturalNumber is an arbitrarily large non-negative integer.
// It is made of a sequence of digits in base 2^32.
// Leading zero-digits are not allowed.
// The value 'zero' is represented by an empty digit sequence.
// digits_[0] is most significant, digits_.back() is least significant.

class NaturalNumber {
public:
    using digit_t = uint32_t;
    using double_digit_t = uint64_t;
    static constexpr int kDigitBits = 32;
    static constexpr double_digit_t kDigitBase = static_cast<double_digit_t>(1) << kDigitBits;
    static constexpr digit_t kDigitMax = static_cast<digit_t>(~0);

private:
    std::vector<digit_t> digits_;

    // Remove leading zeros
    void normalize() {
        auto it = std::ranges::find_if_not(digits_, [](digit_t d) { return d == 0; });
        digits_.erase(digits_.begin(), it);
    }

public:
    // Default constructor creates zero
    NaturalNumber() = default;

    // Construct from uint32_t
    explicit NaturalNumber(digit_t value) {
        if (value != 0) {
            digits_.push_back(value);
        }
    }

    // Construct from uint64_t
    explicit NaturalNumber(double_digit_t value) {
        if (value != 0) {
            auto high = static_cast<digit_t>(value >> kDigitBits);
            auto low = static_cast<digit_t>(value);
            if (high != 0) {
                digits_.push_back(high);
            }
            digits_.push_back(low);
        }
    }

    // Copy constructor
    NaturalNumber(const NaturalNumber&) = default;

    // Move constructor
    NaturalNumber(NaturalNumber&&) noexcept = default;

    // Copy assignment
    NaturalNumber& operator=(const NaturalNumber&) = default;

    // Move assignment
    NaturalNumber& operator=(NaturalNumber&&) noexcept = default;

    // Destructor
    ~NaturalNumber() = default;

    // Check if zero
    [[nodiscard]] bool is_zero() const noexcept {
        return digits_.empty();
    }

    // Number of digits
    [[nodiscard]] std::size_t size() const noexcept {
        return digits_.size();
    }

    // Access digit (0 = most significant)
    [[nodiscard]] digit_t operator[](std::size_t index) const {
        return digits_[index];
    }

    // Get the least significant digit
    [[nodiscard]] digit_t least_significant() const noexcept {
        return digits_.empty() ? 0 : digits_.back();
    }

    // Get number of significant bits
    [[nodiscard]] std::size_t num_bits() const noexcept {
        if (digits_.empty()) {
            return 0;
        }
        return digits_.size() * kDigitBits - std::countl_zero(digits_.front());
    }

    // Convert to uint64_t if possible
    bool to_uint64(double_digit_t& result) const noexcept {
        if (digits_.empty()) {
            result = 0;
            return true;
        }
        if (digits_.size() == 1) {
            result = digits_[0];
            return true;
        }
        if (digits_.size() == 2) {
            result = (static_cast<double_digit_t>(digits_[0]) << kDigitBits) | digits_[1];
            return true;
        }
        return false;
    }

    // Comparison operators
    [[nodiscard]] auto operator<=>(const NaturalNumber& other) const noexcept {
        if (digits_.size() != other.digits_.size()) {
            return digits_.size() <=> other.digits_.size();
        }
        for (std::size_t i = 0; i < digits_.size(); ++i) {
            if (digits_[i] != other.digits_[i]) {
                return digits_[i] <=> other.digits_[i];
            }
        }
        return std::strong_ordering::equal;
    }

    [[nodiscard]] bool operator==(const NaturalNumber& other) const noexcept {
        return digits_ == other.digits_;
    }

    // Addition: this += other
    NaturalNumber& operator+=(const NaturalNumber& other) {
        if (other.is_zero()) return *this;
        if (is_zero()) {
            *this = other;
            return *this;
        }

        // Handle self-addition
        if (this == &other) {
            NaturalNumber copy = other;
            return *this += copy;
        }

        const std::size_t max_size = std::max(digits_.size(), other.digits_.size());
        std::vector<digit_t> result(max_size + 1, 0);

        double_digit_t carry = 0;

        // Add from least significant to most significant
        for (std::size_t i = 0; i < max_size; ++i) {
            std::size_t this_idx = digits_.size() - 1 - i;
            std::size_t other_idx = other.digits_.size() - 1 - i;
            std::size_t result_idx = max_size - i - 1;

            double_digit_t a = (i < digits_.size()) ? digits_[this_idx] : 0;
            double_digit_t b = (i < other.digits_.size()) ? other.digits_[other_idx] : 0;

            double_digit_t sum = a + b + carry;
            result[result_idx] = static_cast<digit_t>(sum);
            carry = sum >> kDigitBits;
        }
        if (carry != 0) {
            result.emplace(result.begin(), static_cast<digit_t>(carry));
        }

        digits_ = std::move(result);
        return *this;
    }

    // Addition with digit_t
    NaturalNumber& operator+=(digit_t value) {
        if (value == 0) return *this;
        if (is_zero()) {
            digits_.push_back(value);
            return *this;
        }

        double_digit_t carry = value;
        for (auto it = digits_.rbegin(); it != digits_.rend() && carry != 0; ++it) {
            double_digit_t sum = static_cast<double_digit_t>(*it) + carry;
            *it = static_cast<digit_t>(sum);
            carry = sum >> kDigitBits;
        }

        if (carry != 0) {
            digits_.insert(digits_.begin(), static_cast<digit_t>(carry));
        }

        return *this;
    }

    // Addition operators
    [[nodiscard]] NaturalNumber operator+(const NaturalNumber& other) const {
        NaturalNumber result = *this;
        result += other;
        return result;
    }

    [[nodiscard]] NaturalNumber operator+(digit_t value) const {
        NaturalNumber result = *this;
        result += value;
        return result;
    }

    // Subtraction: this -= other
    // Precondition: this >= other
    NaturalNumber& operator-=(const NaturalNumber& other) noexcept {
        if (other.is_zero()) return *this;

        // Handle self-subtraction
        if (this == &other) {
            digits_.clear();
            return *this;
        }

        double_digit_t borrow = 0;

        for (std::size_t i = 0; i < digits_.size(); ++i) {
            std::size_t this_idx = digits_.size() - 1 - i;
            std::size_t other_idx = other.digits_.size() - 1 - i;

            double_digit_t a = digits_[this_idx];
            double_digit_t b = (i < other.digits_.size()) ? other.digits_[other_idx] : 0;
            b += borrow;

            if (a < b) {
                a += kDigitBase;
                borrow = 1;
            } else {
                borrow = 0;
            }

            digits_[this_idx] = static_cast<digit_t>(a - b);
        }

        normalize();
        return *this;
    }

    // Subtraction with digit_t
    NaturalNumber& operator-=(digit_t value) noexcept {
        if (value == 0) {
            return *this;
        }
        return *this -= NaturalNumber(value);
    }

    // Subtraction operators
    [[nodiscard]] NaturalNumber operator-(const NaturalNumber& other) const {
        NaturalNumber result = *this;
        result -= other;
        return result;
    }

    [[nodiscard]] NaturalNumber operator-(digit_t value) const {
        NaturalNumber result = *this;
        result -= value;
        return result;
    }

    // Multiplication: this *= other
    NaturalNumber& operator*=(const NaturalNumber& other) {
        *this = *this * other;
        return *this;
    }

    // Multiplication with digit_t
    NaturalNumber& operator*=(digit_t value) {
        if (value == 0 || is_zero()) {
            digits_.clear();
            return *this;
        }
        if (value == 1) {
            return *this;
        }

        double_digit_t carry = 0;
        for (auto it = digits_.rbegin(); it != digits_.rend(); ++it) {
            double_digit_t product = static_cast<double_digit_t>(*it) * value + carry;
            *it = static_cast<digit_t>(product);
            carry = product >> kDigitBits;
        }

        if (carry != 0) {
            digits_.emplace(digits_.begin(), static_cast<digit_t>(carry));
        }

        return *this;
    }

    // Multiplication operators
    [[nodiscard]] NaturalNumber operator*(const NaturalNumber& other) const {
        if (is_zero() || other.is_zero()) {
            return NaturalNumber{};
        }

        // Use standard long multiplication
        std::vector<digit_t> result(digits_.size() + other.digits_.size(), 0);

        for (std::size_t i = 0; i < digits_.size(); ++i) {
            std::size_t this_idx = digits_.size() - 1 - i;
            double_digit_t carry = 0;

            for (std::size_t j = 0; j < other.digits_.size(); ++j) {
                std::size_t other_idx = other.digits_.size() - 1 - j;
                std::size_t result_idx = result.size() - 1 - i - j;

                double_digit_t product = static_cast<double_digit_t>(digits_[this_idx]) *
                                         static_cast<double_digit_t>(other.digits_[other_idx]) +
                                         result[result_idx] + carry;

                result[result_idx] = static_cast<digit_t>(product);
                carry = product >> kDigitBits;
            }

            // Propagate remaining carry
            std::size_t carry_idx = result.size() - 1 - i - other.digits_.size();
            while (carry != 0) {
                double_digit_t sum = static_cast<double_digit_t>(result[carry_idx]) + carry;
                result[carry_idx] = static_cast<digit_t>(sum);
                carry = sum >> kDigitBits;
                if (carry_idx == 0) break;
                --carry_idx;
            }
        }

        NaturalNumber res;
        res.digits_ = std::move(result);
        res.normalize();
        return res;
    }

    [[nodiscard]] NaturalNumber operator*(digit_t value) const {
        NaturalNumber result = *this;
        result *= value;
        return result;
    }

    // Left shift by n bits
    NaturalNumber& operator<<=(uint32_t n) {
        if (is_zero() || n == 0) return *this;

        const uint32_t digit_shift = n / kDigitBits;
        const uint32_t bit_shift = n % kDigitBits;

        // First handle bit shift within digits
        if (bit_shift > 0) {
            const uint32_t reverse_shift = kDigitBits - bit_shift;
            digit_t carry = 0;

            for (auto it = digits_.rbegin(); it != digits_.rend(); ++it) {
                digit_t new_carry = *it >> reverse_shift;
                *it = (*it << bit_shift) | carry;
                carry = new_carry;
            }

            if (carry != 0) {
                digits_.insert(digits_.begin(), carry);
            }
        }

        // Then add zero digits for digit shift
        if (digit_shift > 0) {
            digits_.resize(digits_.size() + digit_shift, 0);
        }

        return *this;
    }

    [[nodiscard]] NaturalNumber operator<<(uint32_t n) const {
        NaturalNumber result = *this;
        result <<= n;
        return result;
    }

    // Right shift by n bits
    NaturalNumber& operator>>=(uint32_t n) {
        if (is_zero() || n == 0) {
            return *this;
        }

        const uint32_t total_bits = static_cast<uint32_t>(num_bits());
        if (n >= total_bits) {
            digits_.clear();
            return *this;
        }

        const uint32_t digit_shift = n / kDigitBits;
        const uint32_t bit_shift = n % kDigitBits;

        // First remove whole digits
        if (digit_shift > 0) {
            digits_.resize(digits_.size() - digit_shift);
        }

        // Then handle bit shift within digits
        if (bit_shift > 0 && !digits_.empty()) {
            const uint32_t reverse_shift = kDigitBits - bit_shift;
            digit_t carry = 0;

            for (auto& digit : digits_) {
                digit_t new_carry = digit << reverse_shift;
                digit = (digit >> bit_shift) | carry;
                carry = new_carry;
            }
        }

        normalize();
        return *this;
    }

    [[nodiscard]] NaturalNumber operator>>(uint32_t n) const {
        NaturalNumber result = *this;
        result >>= n;
        return result;
    }

    // Get bit at position i (0 = least significant)
    [[nodiscard]] bool get_bit(uint32_t i) const {
        if (is_zero()) return false;
        const uint32_t digit_idx = static_cast<uint32_t>(digits_.size()) - 1 - (i / kDigitBits);
        if (digit_idx >= digits_.size()) return false;
        const uint32_t bit_idx = i % kDigitBits;
        return (digits_[digit_idx] & (static_cast<digit_t>(1) << bit_idx)) != 0;
    }

    // Count trailing zeros (zeros in least significant positions)
    [[nodiscard]] uint32_t count_trailing_zeros() const {
        if (is_zero()) {
            return 0;
        }

        uint32_t count = 0;
        for (auto it = digits_.rbegin(); it != digits_.rend(); ++it) {
            if (*it == 0) {
                count += kDigitBits;
            } else {
                // Count trailing zeros in this digit
                return count + std::countr_zero(*it);
            }
        }
        return count;
    }

    // String representation in decimal
    [[nodiscard]] std::string to_decimal_string() const {
        if (is_zero()) return "0";

        std::string result;
        NaturalNumber temp = *this;

        while (!temp.is_zero()) {
            digit_t remainder = temp.divide_by_digit(10);
            result.push_back(static_cast<char>('0' + remainder));
        }

        std::reverse(result.begin(), result.end());
        return result;
    }

    // String representation in binary
    [[nodiscard]] std::string to_binary_string() const {
        if (is_zero()) return "0";

        std::string result;
        for (const auto& digit : digits_) {
            for (int i = kDigitBits - 1; i >= 0; --i) {
                result.push_back((digit & (static_cast<digit_t>(1) << i)) ? '1' : '0');
            }
        }

        // Remove leading zeros
        auto first_one = result.find('1');
        if (first_one != std::string::npos) {
            result = result.substr(first_one);
        }

        return result;
    }

private:
    // Helper: divide by single digit, returns remainder
    digit_t divide_by_digit(digit_t divisor) {
        if (divisor == 0) {
            throw std::domain_error("Division by zero");
        }

        double_digit_t remainder = 0;
        for (auto& digit : digits_) {
            double_digit_t dividend = (remainder << kDigitBits) | digit;
            digit = static_cast<digit_t>(dividend / divisor);
            remainder = dividend % divisor;
        }

        normalize();
        return static_cast<digit_t>(remainder);
    }
};

// Free function operators for symmetry
[[nodiscard]] inline NaturalNumber operator+(NaturalNumber::digit_t lhs, const NaturalNumber& rhs) {
    return rhs + lhs;
}

[[nodiscard]] inline NaturalNumber operator*(NaturalNumber::digit_t lhs, const NaturalNumber& rhs) {
    return rhs * lhs;
}

// Stream output
inline std::ostream& operator<<(std::ostream& os, const NaturalNumber& num) {
    return os << num.to_decimal_string();
}

} // namespace predicates
