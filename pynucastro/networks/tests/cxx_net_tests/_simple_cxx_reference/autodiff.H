#ifndef AUTODIFF_H
#define AUTODIFF_H

#include <cassert>
#include <cmath>
#include <compare>
#include <type_traits>

namespace dual {

    /// Derivatives with respect to two independent variables, without arrays.
    template <typename T>
    struct DerivativePair {

        T d1{};
        T d2{};

        friend bool operator== (const DerivativePair&, const DerivativePair&) = default;

        friend DerivativePair operator+ (const DerivativePair& a, const DerivativePair& b) {
            return {a.d1 + b.d1, a.d2 + b.d2};
        }

        friend DerivativePair operator- (const DerivativePair& a, const DerivativePair& b) {
            return {a.d1 - b.d1, a.d2 - b.d2};
        }

        DerivativePair operator- () const {return {-d1, -d2};}

        friend DerivativePair operator* (const DerivativePair& a, T scale) {
            return {a.d1 * scale, a.d2 * scale};
        }

        friend DerivativePair operator* (T scale, const DerivativePair& a) {
            return a * scale;
        }

        friend DerivativePair operator/ (const DerivativePair& a, T scale) {
            return {a.d1 / scale, a.d2 / scale};
        }

    };

    ///
    /// A scalar value and a derivative seed (a scalar or DerivativePair).
    ///
    template <typename T, typename D = T>
    class Dual {

    public:

        T value{};
        D deriv{};

        Dual (T _value, D _deriv=D{})
          : value{_value}, deriv{_deriv} {}

        // For a pair, the default seed selects d1; pass {0, 1} for d2.
        void seed(D derivative=D{1}) {
            deriv = derivative;
        }

        // note, since our constructor is not explicit, these operators
        // work for Dual + Dual, Dual + T, or T + Dual.  The ideal is that
        // if we do Dual + double, then the double will first be cast to
        // a Dual.

        friend auto operator<=> (const Dual& lhs, const Dual& rhs) {
            return lhs.value <=> rhs.value;
        }

        friend bool operator== (const Dual& lhs, const Dual& rhs) {
            return lhs.value == rhs.value;
        }

        friend Dual operator+ (const Dual& lhs, const Dual& rhs) {
            return {lhs.value + rhs.value,
                    lhs.deriv + rhs.deriv};
        }

        friend Dual operator- (const Dual& lhs, const Dual& rhs) {
            return {lhs.value - rhs.value,
                    lhs.deriv - rhs.deriv};
        }

        Dual operator- () const {return {-value, -deriv};}

        friend Dual operator* (const Dual& lhs, const Dual& rhs) {
            return {lhs.value * rhs.value,
                    lhs.value * rhs.deriv + lhs.deriv * rhs.value};
        }

        friend Dual operator/ (const Dual& lhs, const Dual& rhs) {
            return {lhs.value / rhs.value,
                    -lhs.value / (rhs.value * rhs.value) * rhs.deriv +
                    lhs.deriv / rhs.value};
        }

    };

    template <typename T, typename D>
    Dual<T, D> constant_like (const Dual<T, D>&, std::type_identity_t<T> value) {
        return {value, D{}};
    }

    // Ties select the first operand. Scalar bounds convert to zero-seeded Duals.
    template <typename T, typename D>
    Dual<T, D> min (const Dual<T, D>& a, const std::type_identity_t<Dual<T, D>>& b) {
        return b.value < a.value ? b : a;
    }

    template <typename U, typename T, typename D>
        requires std::is_arithmetic_v<U>
    Dual<T, D> min (U a, const Dual<T, D>& b) {
        return min(Dual<T, D>(static_cast<T>(a)), b);
    }

    template <typename T, typename D>
    Dual<T, D> max (const Dual<T, D>& a, const std::type_identity_t<Dual<T, D>>& b) {
        return b.value > a.value ? b : a;
    }

    template <typename U, typename T, typename D>
        requires std::is_arithmetic_v<U>
    Dual<T, D> max (U a, const Dual<T, D>& b) {
        return max(Dual<T, D>(static_cast<T>(a)), b);
    }

    // A constant exponent does not need the logarithm used for a Dual exponent.
    template <typename T, typename D, typename U>
        requires (std::is_floating_point_v<U> || std::is_integral_v<U>)
    Dual<T, D> pow (const Dual<T, D>& x, U y) {

        if (y == U{0}) {
            return {T{1}, D{}};
        }
        if (y == U{1}) {
            return x;
        }

        // Subtract in floating point to avoid overflowing a minimum integer.
        using Exponent = std::conditional_t<std::is_integral_v<U>, long double, U>;
        auto exponent = static_cast<Exponent>(y);
        T value = std::pow(x.value, exponent);
        D deriv = exponent * std::pow(x.value, exponent - Exponent{1}) * x.deriv;
        return {value, deriv};

    }

    // Deduce T from x, allowing convertible exponents as well as Dual<T, D>.
    template <typename T, typename D>
    Dual<T, D> pow (const Dual<T, D>& x, const std::type_identity_t<Dual<T, D>>& y) {

        if (y.deriv == D{}) {
            return pow(x, y.value);
        }

        // At zero with a positive exponent, the exponent partial is zero
        // (the limit from positive bases). Use the constant-exponent partial.
        if (x.value == T{} && y.value > T{}) {
            return pow(x, y.value);
        }

        // A varying exponent otherwise requires a strictly positive base.
        assert (x.value > T{});

        if (y.value == T{}) {
            return {T{1}, y.deriv * std::log(x.value)};
        }

        T value = std::pow(x.value, y.value);
        D deriv = value * (y.deriv * std::log(x.value) +
                           y.value * x.deriv / x.value);
        return {value, deriv};

    }

    template <typename T, typename D>
    Dual<T, D> abs (const Dual<T, D>& x) {
        if (std::isnan(x.value)) {
            return x;
        }
        if (x.value < T{}) {
            return {-x.value, -x.deriv};
        } else if (x.value > T{}) {
            return x;
        }
        return {T{}, D{}};
    }

    template <typename T, typename D>
    Dual<T, D> exp (const Dual<T, D>& x) {
        auto value = std::exp(x.value);
        return {value, value * x.deriv};
    }

    template <typename T, typename D>
    Dual<T, D> log (const Dual<T, D>& x) {
        auto value = std::log(x.value);
        return {value, x.deriv / x.value};
    }

    template <typename T, typename D>
    Dual<T, D> log1p (const Dual<T, D>& x) {
        auto value = std::log1p(x.value);
        return {value, x.deriv / (T{1} + x.value)};
    }

    template <typename T, typename D>
    Dual<T, D> sqrt (const Dual<T, D>& x) {
        auto value = std::sqrt(x.value);
        return {value, T{0.5} * x.deriv / value};
    }

    template <typename T, typename D>
    Dual<T, D> cbrt (const Dual<T, D>& x) {
        auto value = std::cbrt(x.value);
        return {value, (T{1} / T{3}) * x.deriv / (value * value)};
    }

    template <typename T, typename D>
    Dual<T, D> cos (const Dual<T, D>& x) {
        return {std::cos(x.value), -std::sin(x.value) * x.deriv};
    }

    template <typename T, typename D>
    Dual<T, D> sin (const Dual<T, D>& x) {
        return {std::sin(x.value), std::cos(x.value) * x.deriv};
    }

    template <typename T, typename D>
    Dual<T, D> atan (const Dual<T, D>& x) {
        return {std::atan(x.value), x.deriv / (T{1} + x.value * x.value)};
    }

}

#endif
