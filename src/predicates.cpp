#include "predicates/predicates.hpp"
#include "predicates/interval_number.hpp"
#include "predicates/float_number.hpp"

#include <algorithm>
#include <type_traits>

namespace predicates {
template<typename T>
auto orient1d_one_impl(const T& p0, const T& p1, const T& q0, const T& q1) noexcept {
    T a_ = q0 - p0;
    T b_ = q1 - p1;
    T result = a_ - b_;
    if constexpr (std::is_floating_point_v<T>) {
        return std::array<T, 2>{{result, std::max(std::abs(a_), std::abs(b_))}};
    } else {
        return result;
    }
}

template<typename T>
auto orient2d_impl(const T& p1x, const T& p1y, const T& p2x, const T& p2y, const T& p3x, const T& p3y) noexcept {
    T dl  = (p2x - p1x) * (p3y - p1y);
    T dr  = (p2y - p1y) * (p3x - p1x);
    T det = dl - dr;
    if constexpr (std::is_floating_point_v<T>) {
        return std::array<T, 2>{{det, std::abs(dl) + std::abs(dr)}};
    } else {
        return det;
    }
}

template<typename T>
auto orient2d_one_impl(const T& p0, const T& p1, const T& p2, const T& q0, const T& q1, const T& q2, const T& r0, const T& r1, const T& r2) noexcept {
    T a_ = q0 - p0;
    T b_ = q1 - p1;
    T c_ = q2 - p2;
    T d_ = r0 - p0;
    T e_ = r1 - p1;
    T f_ = r2 - p2;
    T bf = b_ * f_;
    T ce = c_ * e_;
    T af = a_ * f_;
    T cd = c_ * d_;
    T ae = a_ * e_;
    T bd = b_ * d_;
    T bf_ce = bf - ce;
    T af_cd = af - cd;
    T ae_bd = ae - bd;
    T result0 = bf_ce - af_cd;
    T result = result0 + ae_bd;
    if constexpr (std::is_floating_point_v<T>) {
        return std::array<T, 2>{{result, std::max({std::abs(a_), std::abs(b_), std::abs(c_), std::abs(d_), std::abs(e_), std::abs(f_)})}};
    } else {
        return result;
    }
}

template<typename T>
auto orient3d_impl(
    const T& px, const T& py, const T& pz,
    const T& qx, const T& qy, const T& qz,
    const T& rx, const T& ry, const T& rz,
    const T& sx, const T& sy, const T& sz
) noexcept {
    T fadx = qx - px;
    T fbdx = rx - px;
    T fcdx = sx - px;
    T fady = qy - py;
    T fbdy = ry - py;
    T fcdy = sy - py;
    T fadz = qz - pz;
    T fbdz = rz - pz;
    T fcdz = sz - pz;

    T fbdxcdy = fbdx * fcdy * fadz;
    T fcdxbdy = fcdx * fbdy * fadz;
    T fcdxady = fcdx * fady * fbdz;
    T fadxcdy = fadx * fcdy * fbdz;
    T fadxbdy = fadx * fbdy * fcdz;
    T fbdxady = fbdx * fady * fcdz;

    T det = (fbdxcdy - fcdxbdy) + (fcdxady - fadxcdy) + (fadxbdy - fbdxady);
    if constexpr (std::is_floating_point_v<T>) {
        return std::array<T, 2>{{det, std::abs(fbdxcdy) + std::abs(fcdxbdy) + std::abs(fcdxady) +
                                std::abs(fadxcdy) + std::abs(fadxbdy) + std::abs(fbdxady)}};
    } else {
        return det;
    }
}

template<typename T>
auto orient3d_one_impl(
    const T& p0, const T& p1, const T& p2, const T& p3,
    const T& q0, const T& q1, const T& q2, const T& q3,
    const T& r0, const T& r1, const T& r2, const T& r3,
    const T& s0, const T& s1, const T& s2, const T& s3
) noexcept {
    T a_ = q0 - p0;
    T b_ = q1 - p1;
    T c_ = q2 - p2;
    T d_ = q3 - p3;
    T e_ = r0 - p0;
    T f_ = r1 - p1;
    T g_ = r2 - p2;
    T h_ = r3 - p3;
    T i_ = s0 - p0;
    T j_ = s1 - p1;
    T k_ = s2 - p2;
    T l_ = s3 - p3;
    T af = a_ * f_;
    T be = b_ * e_;
    T ce = c_ * e_;
    T ag = a_ * g_;
    T ah = a_ * h_;
    T de = d_ * e_;
    T bg = b_ * g_;
    T cf = c_ * f_;
    T df = d_ * f_;
    T bh = b_ * h_;
    T ch = c_ * h_;
    T dg = d_ * g_;
    T d1 = af - be;
    T d2 = k_ - l_;
    T d3 = ce - ag;
    T d4 = j_ - l_;
    T d5 = ah - de;
    T d6 = j_ - k_;
    T d7 = bg - cf;
    T d8 = i_ - l_;
    T d9 = df - bh;
    T d10 = i_ - k_;
    T d11 = ch - dg;
    T d12 = i_ - j_;
    T term1 = d1 * d2;
    T term2 = d3 * d4;
    T term3 = d5 * d6;
    T term4 = d7 * d8;
    T term5 = d9 * d10;
    T term6 = d11 * d12;
    T r12 = term1 + term2;
    T r34 = term3 + term4;
    T r56 = term5 + term6;
    T r1234 = r12 + r34;
    T r123456 = r1234 + r56;
    if constexpr (std::is_floating_point_v<T>) {
        return std::array<T, 2>{{r123456, std::max({std::abs(a_), std::abs(b_), std::abs(c_), std::abs(d_), std::abs(e_), std::abs(f_), std::abs(g_), std::abs(h_), std::abs(i_), std::abs(j_), std::abs(k_), std::abs(l_)})}};
    } else {
        return r123456;
    }
}
template<typename T>
auto orient4d_impl(
    const T& p0, const T& p1, const T& p2, const T& p3,
    const T& q0, const T& q1, const T& q2, const T& q3,
    const T& r0, const T& r1, const T& r2, const T& r3,
    const T& s0, const T& s1, const T& s2, const T& s3,
    const T& t0, const T& t1, const T& t2, const T& t3
) noexcept {
    T a_ = q0 - p0;
    T b_ = q1 - p1;
    T c_ = q2 - p2;
    T d_ = q3 - p3;
    T e_ = r0 - p0;
    T f_ = r1 - p1;
    T g_ = r2 - p2;
    T h_ = r3 - p3;
    T i_ = s0 - p0;
    T j_ = s1 - p1;
    T k_ = s2 - p2;
    T l_ = s3 - p3;
    T m_ = t0 - p0;
    T n_ = t1 - p1;
    T o_ = t2 - p2;
    T p_ = t3 - p3;
    T af = a_ * f_;
    T be = b_ * e_;
    T kp = k_ * p_;
    T lo = l_ * o_;
    T ce = c_ * e_;
    T ag = a_ * g_;
    T jp = j_ * p_;
    T ln = l_ * n_;
    T ah = a_ * h_;
    T de = d_ * e_;
    T jo = j_ * o_;
    T kn = k_ * n_;
    T bg = b_ * g_;
    T cf = c_ * f_;
    T ip = i_ * p_;
    T lm = l_ * m_;
    T df = d_ * f_;
    T bh = b_ * h_;
    T io = i_ * o_;
    T km = k_ * m_;
    T ch = c_ * h_;
    T dg = d_ * g_;
    T in = i_ * n_;
    T jm = j_ * m_;
    T d1 = af - be;
    T d2 = kp - lo;
    T d3 = ce - ag;
    T d4 = jp - ln;
    T d5 = ah - de;
    T d6 = jo - kn;
    T d7 = bg - cf;
    T d8 = ip - lm;
    T d9 = df - bh;
    T d10 = io - km;
    T d11 = ch - dg;
    T d12 = in - jm;
    T term1 = d1 * d2;
    T term2 = d3 * d4;
    T term3 = d5 * d6;
    T term4 = d7 * d8;
    T term5 = d9 * d10;
    T term6 = d11 * d12;
    T r12 = term1 + term2;
    T r34 = term3 + term4;
    T r56 = term5 + term6;
    T r1234 = r12 + r34;
    T r123456 = r1234 + r56;
    if constexpr (std::is_floating_point_v<T>) {
        return std::array<T, 2>{{r123456, std::max({std::abs(a_), std::abs(b_), std::abs(c_), std::abs(d_), std::abs(e_), std::abs(f_), std::abs(g_), std::abs(h_), std::abs(i_), std::abs(j_), std::abs(k_), std::abs(l_), std::abs(m_), std::abs(n_), std::abs(o_), std::abs(p_)})}};
    } else {
        return r123456;
    }
}

double orient1d_one(const double* pa, const double* pb) noexcept {
    const auto [det, max] = orient1d_one_impl<double>(pa[0], pa[1], pb[0], pb[1]);
    if (std::abs(det) > max * 4.440892098500627e-16) {
        return det;
    }
    return orient1d_one_impl<FloatNumber>(pa[0], pa[1], pb[0], pb[1]).to_double();
}

double orient2d(const double* pa, const double* pb, const double* pc) noexcept {
    const auto [det, max] = orient2d_impl<double>(pa[0], pa[1], pb[0], pb[1], pc[0], pc[1]);
    if (std::abs(det) > max * 3.3306690738754706e-016) {
        return det;
    }
    return orient2d_impl<FloatNumber>(pa[0], pa[1], pb[0], pb[1], pc[0], pc[1]).to_double();
}

double orient2d_one(const double* pa, const double* pb, const double* pc) noexcept {
    const auto [det, max] = orient2d_one_impl<double>(pa[0], pa[1], pa[2], pb[0], pb[1], pb[2], pc[0], pc[1], pc[2]);
    if (std::abs(det) > max * max * 3.552713678800501e-15 ) {
        return det;
    }
    return orient2d_one_impl<FloatNumber>(pa[0], pa[1], pa[2], pb[0], pb[1], pb[2], pc[0], pc[1], pc[2]).to_double();
}

double orient3d(const double* pa, const double* pb, const double* pc, const double* pd) noexcept {
    {
        const auto [det, max] = orient3d_impl<double>(pa[0], pa[1], pa[2], pb[0], pb[1], pb[2], pc[0], pc[1], pc[2], pd[0], pd[1], pd[2]);
        if (std::abs(det) > max * 7.7715611723761027e-016) {
            return det;
        }
    }
    {
        const auto det = orient3d_impl<IntervalNumber>(pa[0], pa[1], pa[2], pb[0], pb[1], pb[2], pc[0], pc[1], pc[2], pd[0], pd[1], pd[2]);
        if (det.defined()) {
            return det.to_double();
        }
    }
    {
        return orient3d_impl<FloatNumber>(pa[0], pa[1], pa[2], pb[0], pb[1], pb[2], pc[0], pc[1], pc[2], pd[0], pd[1], pd[2]).to_double();
    }
}

double orient3d_one(const double* pa, const double* pb, const double* pc, const double* pd) noexcept {
    {
        const auto [det, max] = orient3d_one_impl<double>(pa[0], pa[1], pa[2], pa[3], pb[0], pb[1], pb[2], pb[3], pc[0], pc[1], pc[2], pc[3], pd[0], pd[1], pd[2], pd[3]);
        if (std::abs(det) > max * max * max * 2.486899575160352e-14) {
            return det;
        }
    }
    {
        const auto det = orient3d_one_impl<IntervalNumber>(pa[0], pa[1], pa[2], pa[3], pb[0], pb[1], pb[2], pb[3], pc[0], pc[1], pc[2], pc[3], pd[0], pd[1], pd[2], pd[3]);
        if (det.defined()) {
            return det.to_double();
        }
    }
    {
        return orient3d_one_impl<FloatNumber>(pa[0], pa[1], pa[2], pa[3], pb[0], pb[1], pb[2], pb[3], pc[0], pc[1], pc[2], pc[3], pd[0], pd[1], pd[2], pd[3]).to_double();
    }
}

double orient4d(const double* pa, const double* pb, const double* pc, const double* pd, const double* pe) noexcept {
    {
        auto [det, max] = orient4d_impl<double>(pa[0], pa[1], pa[2], pa[3], pb[0], pb[1], pb[2], pb[3], pc[0], pc[1], pc[2], pc[3], pd[0], pd[1], pd[2], pd[3], pe[0], pe[1], pe[2], pe[3]);
        max *= max;
        max *= max;
        if (std::abs(det) > max * 3.019806626980427e-14) {
            return det;
        }
    }
    {
        const auto det = orient4d_impl<IntervalNumber>(pa[0], pa[1], pa[2], pa[3], pb[0], pb[1], pb[2], pb[3], pc[0], pc[1], pc[2], pc[3], pd[0], pd[1], pd[2], pd[3], pe[0], pe[1], pe[2], pe[3]);
        if (det.defined()) {
            return det.to_double();
        }
    }
    {
        return orient4d_impl<FloatNumber>(pa[0], pa[1], pa[2], pa[3], pb[0], pb[1], pb[2], pb[3], pc[0], pc[1], pc[2], pc[3], pd[0], pd[1], pd[2], pd[3], pe[0], pe[1], pe[2], pe[3]).to_double();
    }
}

std::array<double, 2> mi_orient0d(const double a, const double b) noexcept {
    return {{a - b, 1.0}};
}

std::array<double, 2> mi_orient1d(const double* pa, const double* pb, const double* pc) noexcept {
    const auto numerator = orient2d(pa, pb, pc);
    const auto denominator = orient1d_one(pa, pb);
    if (denominator > 0.0) {
        return {{numerator, denominator}};
    } else {
        return {{ -numerator, -denominator }};
    }
}
std::array<double, 2> mi_orient2d(const double* pa, const double* pb, const double* pc, const double* pd) noexcept {
    const auto numerator = orient3d(pa, pb, pc, pd);
    const auto denominator = orient2d_one(pa, pb, pc);
    if (denominator > 0.0) {
        return {{numerator, denominator}};
    } else {
        return {{ -numerator, -denominator }};
    }
}

std::array<double, 2> mi_orient4d(const double* pa, const double* pb, const double* pc, const double* pd, const double* pe) noexcept {
    const auto numerator = orient4d(pa, pb, pc, pd, pe);
    const auto denominator = orient3d_one(pa, pb, pc, pd);
    if (denominator > 0.0) {
        return {{numerator, denominator}};
    } else {
        return {{ -numerator, -denominator }};
    }
}
}  // namespace predicates
