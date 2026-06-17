#ifndef INTERP_TOOLS_H
#define INTERP_TOOLS_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

namespace interp_net {

    // index locator
    // return the index i such that x_array[i] <= x0 < x_array[i+1]

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    int
    find_index(const amrex::Real x0, const T& x_array) {

        int left = x_array.lo();
        int right = x_array.hi();

        int idx = -1;

        if (x0 >= x_array(left) && x0 < x_array(right)) {

            // find the largest x element <= x0 using a
            // binary search

            while (left < right) {
                int mid = (left + right) / 2;
                if (x_array(mid) > x0) {
                    right = mid;
                } else {
                    left = mid + 1;
                }
            }

            idx = right - 1;
        }

        return idx;
    }

    // an index locator that extrapolates at the boundaries instead
    // of returning -1

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    int
    find_index_extrap(const amrex::Real x0, const T& x_array) {

        int left = x_array.lo();
        int right = x_array.hi();

        int idx;

        if (x0 < x_array(left)) {
            idx = left;
        } else if (x0 >= x_array(right)) {
            idx = right - 1;
        } else {
            idx = find_index(x0, x_array);
        }

        return idx;
    }

    template <int do_derivative, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    std::pair<amrex::Real, amrex::Real>
    linear_extrap(const amrex::Real xv, const int idx0, const int idx1,
                  const T& x_array, const T& f_array) {

        const auto x0 = x_array(idx0);
        const auto x1 = x_array(idx1);

        const auto f0 = f_array(idx0);
        const auto f1 = f_array(idx1);

        const amrex::Real slope = (f1 - f0) / (x1 - x0);
        const amrex::Real f = f0 + slope * (xv - x0);

        amrex::Real dfdx {0.0_rt};
        if constexpr (do_derivative) {
            dfdx = slope;
        }

        return {f, dfdx};
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    amrex::Real
    limit_slope(const amrex::Real hm, const amrex::Real hp,
                const amrex::Real dm, const amrex::Real dp) {
        // Slope Limiter to preserve monotonicity
        // See C. Moler, Numerical Computing with Matlab, 2004.
        // :doi:`10.1137/1.9780898717952`, Chapter 3.4

        if (dm * dp <= 0.0_rt) {
            return 0.0_rt;
        }

        const amrex::Real w1 = 2.0_rt * hp + hm;
        const amrex::Real w2 = hp + 2.0_rt * hm;
        return (w1 + w2) / (w1 / dm + w2 / dp);
    }

    // Monotone cubic Hermite interpolation
    // Uses Hermite polynomials as basis to construct the cubic interpolant
    //
    // H(x) = f0 H_{1,0}(x) + f1 H_{1,1}(x) + f'0 Ĥ_{1,0}(x) + f'1 Ĥ_{1,0}(x)
    //
    // The slope f' is limited to preserve monotonicity.
    // Also use simple one-sided slope for edge cases
    // this version works for unevenly spaced data points
    // it will return f and f' at x = xv
    //
    // See C. Moler, Numerical Computing with Matlab, 2004.
    // :doi:`10.1137/1.9780898717952`, Chapter 3.4 and 3.6

    template <int do_derivative, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    std::pair<amrex::Real, amrex::Real>
    monotone_cubic_interp(const amrex::Real xv, const int idx,
                          const T& x_array, const T& f_array) {

        const int ilo = x_array.lo();
        const int ihi = x_array.hi();

        // Clamp index within table boundary
        // Note idx should already be clamped within [ilo, ihi-1]
        AMREX_ASSERT(idx >= ilo && idx < ihi);

        const int im1 = amrex::max(idx - 1, ilo);
        const int ip1 = idx + 1;
        const int ip2 = amrex::min(idx + 2, ihi);

        const auto x_im1 = x_array(im1);
        const auto x_idx = x_array(idx);
        const auto x_ip1 = x_array(ip1);
        const auto x_ip2 = x_array(ip2);

        const auto f_im1 = f_array(im1);
        const auto f_idx = f_array(idx);
        const auto f_ip1 = f_array(ip1);
        const auto f_ip2 = f_array(ip2);

        // Compute the slope at i0 and i0 + 1
        amrex::Real m0, m1;
        m1 = m0 = (f_ip1 - f_idx) / (x_ip1 - x_idx);

        // If NOT on the edge, compute the slope from the other side
        // and limit the slope to preserve monotonicity
        if (idx != ilo) {
            amrex::Real d_im1 = (f_idx - f_im1) / (x_idx - x_im1);
            m0 = limit_slope(x_idx - x_im1, x_ip1 - x_idx, d_im1, m0);
        }

        if (ip1 != ihi) {
            amrex::Real d_ip1 = (f_ip2 - f_ip1) / (x_ip2 - x_ip1);
            m1 = limit_slope(x_ip1 - x_idx, x_ip2 - x_ip1, m1, d_ip1);
        }

        // Compute cubic Hermite polynomial basis
        const amrex::Real h = x_ip1 - x_idx;
        const amrex::Real t = (xv - x_idx) / h;
        const amrex::Real t2 = t * t;
        const amrex::Real t3 = t2 * t;

        const amrex::Real H1 = -2.0_rt * t3 + 3.0_rt * t2;
        const amrex::Real H0 = -H1 + 1.0_rt;
        const amrex::Real Hhat1 = h * (t3 - t2);
        const amrex::Real Hhat0 = h * (t3 - 2.0_rt * t2 + t);

        // Construct the interpolant
        // H_3(x) = f0 H_{1,0}(x) + f1 H_{1,1}(x) + f'0 Ĥ_{1,0}(x) + f'1 Ĥ_{1,0}(x)

        const amrex::Real fv = f_idx * H0 + f_ip1 * H1 + m0 * Hhat0 + m1 * Hhat1;

        // Compute the derivative if necessary
        amrex::Real dfdx {0.0_rt};
        if constexpr (do_derivative) {
            const amrex::Real dH1dx = -6.0_rt * (t2 - t) / h;
            const amrex::Real dH0dx = -dH1dx;
            const amrex::Real dHhat1dx = 3.0_rt * t2 - 2.0_rt * t;
            const amrex::Real dHhat0dx = 3.0_rt * t2 - 4.0_rt * t + 1.0_rt;
            dfdx = f_idx * dH0dx + f_ip1 * dH1dx + m0 * dHhat0dx + m1 * dHhat1dx;
        }

        return {fv, dfdx};
    }

    template <int do_derivative, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    std::pair<amrex::Real, amrex::Real>
    monotone_1d_interp(const amrex::Real xv, const T& x_array, const T& f_array) {
        // Main driver for interpolating and extrapolating 1D table data.
        // Uses monotone Hermite cubic interpolation and linear extrapolation
        // to preserve monotonicity.

        // find the index such that x_array[idx] <= xv < x_array[idx+1]
        // We already ensure that ilo <= i <= ihi - 1
        const int idx = interp_net::find_index_extrap(xv, x_array);

        const int ilo = x_array.lo();
        const int ihi = x_array.hi();

        if (xv < x_array(ilo)) {
            return linear_extrap<do_derivative>(xv, ilo, ilo+1, x_array, f_array);
        }

        if (xv > x_array(ihi)){
            return linear_extrap<do_derivative>(xv, ihi, ihi-1, x_array, f_array);
        }

        return monotone_cubic_interp<do_derivative>(xv, idx, x_array, f_array);
    }

}
#endif
