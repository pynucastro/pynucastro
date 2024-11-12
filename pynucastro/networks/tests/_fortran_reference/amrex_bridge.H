#ifndef AMREX_BRIDGE_H
#define AMREX_BRIDGE_H

#include <cassert>

// this is a set of classes and data that allow us to use some
// AMReX concepts in the base C++ network.  Note: many of these
// are stripped down from their AMReX equivalents


// adapted from AMReX_REAL.H

using Real = double;

inline namespace literals {

    constexpr Real
    operator"" _rt( long double x )
    {
        return Real( x );
    }

    constexpr Real
    operator"" _rt( unsigned long long int x )
    {
        return Real( x );
    }

} // namespace literals


// adapted from AMReX_Array.H

template <class T, int XLO, int XHI>
struct Array1D
{

    // Returns the number of elements in the Array1D object as an unsigned
    // integer.
    [[nodiscard]] inline
    static constexpr unsigned int size () noexcept { return (XHI-XLO+1); }

    // Returns the number of elements in the Array1D object as an unsigned
    // integer.
    [[nodiscard]] inline
    static constexpr unsigned int len () noexcept { return (XHI-XLO+1); }


    // The elements of an Array1D object are accessed using parentheses,
    // array(i), instead of using square brackets.
    [[nodiscard]] inline
    const T& operator() (int i) const noexcept {
        assert (i >= XLO && i <= XHI);
        return arr[i-XLO];
    }

    [[nodiscard]] inline
    T& operator() (int i) noexcept {
        assert(i >= XLO && i <= XHI);
        return arr[i-XLO];
    }

    // Array1D is implemented as a fixed-size array.
    T arr[(XHI-XLO+1)];
};

// a 2D array with column-major ordering
// adapted from AMReX-Microphysics
template <int XLO, int XHI, int YLO, int YHI>
struct MathArray2D
{
    inline
    void zero() {
        for (int i = 0; i < (YHI-YLO+1)*(XHI-XLO+1); ++i) {
            arr[i] = 0.0_rt;
        }
    }

    inline
    void set (const int i, const int j, const Real x) noexcept {
        assert(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
        arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)] = x;
    }

    [[nodiscard]] inline
    Real get (const int i, const int j) const noexcept {
        assert(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
        return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
    }

    inline
    const Real& operator() (int i, int j) const noexcept {
        assert(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
        return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
    }

    inline
    Real& operator() (int i, int j) noexcept {
        assert(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
        return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
    }

    Real arr[(XHI-XLO+1)*(YHI-YLO+1)];
};


// copied from AMReX.H

namespace amrex {

template <class... Ts>
inline
void ignore_unused (const Ts&...) {}

} // namespace amrex

#endif
