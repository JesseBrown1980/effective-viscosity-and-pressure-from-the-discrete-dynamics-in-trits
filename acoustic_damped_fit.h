#ifndef ASOLARIA_ACOUSTIC_DAMPED_FIT_H
#define ASOLARIA_ACOUSTIC_DAMPED_FIT_H

// Measurement-only floating-point least squares. This does not modify the
// integer/trit collision or streaming dynamics.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace acoustic {

struct FitResult {
    bool valid = false;
    double omega = std::numeric_limits<double>::quiet_NaN();
    double alpha = std::numeric_limits<double>::quiet_NaN();
    double offset = std::numeric_limits<double>::quiet_NaN();
    double cosine = std::numeric_limits<double>::quiet_NaN();
    double sine = std::numeric_limits<double>::quiet_NaN();
    double r_squared = std::numeric_limits<double>::quiet_NaN();
    double rmse = std::numeric_limits<double>::quiet_NaN();
    double cycles = std::numeric_limits<double>::quiet_NaN();
    bool frequency_at_boundary = false;
    bool damping_at_boundary = false;
};

// Coefficients are referenced to t.front(), so callers supply relative time.
inline double predict(const FitResult& fit_result, double relative_time) {
    return fit_result.offset + std::exp(-fit_result.alpha * relative_time) *
        (fit_result.cosine * std::cos(fit_result.omega * relative_time) +
         fit_result.sine * std::sin(fit_result.omega * relative_time));
}

namespace detail {
constexpr double pi = 3.141592653589793238462643383279502884;

struct Candidate {
    double w = 0.0; // omega * observed duration
    double a = 0.0; // alpha * observed duration
    double sse = std::numeric_limits<double>::infinity();
    double offset = 0.0;
    double cosine = 0.0;
    double sine = 0.0;
};

struct Projector {
    const std::vector<double>& times; // relative times divided by duration
    const std::vector<double>& values; // centred and scaled observations

    Candidate operator()(double w, double a) const {
        Candidate result;
        result.w = w;
        result.a = a;
        const std::size_t n = times.size();
        const double nd = static_cast<double>(n);
        std::vector<double> qc(n), qs(n);
        double mean_c = 0.0, mean_s = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double envelope = std::exp(-a * times[i]);
            qc[i] = envelope * std::cos(w * times[i]);
            qs[i] = envelope * std::sin(w * times[i]);
            mean_c += qc[i] / nd;
            mean_s += qs[i] / nd;
        }
        // QR variable projection: remove the constant column, normalize the
        // cosine column, then orthogonalize the sine column twice. Do not form
        // X'X: normal equations square the condition number of this design.
        double norm_c = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            qc[i] -= mean_c;
            qs[i] -= mean_s;
            norm_c += qc[i] * qc[i];
        }
        norm_c = std::sqrt(norm_c);
        const double rank_tolerance = 1e-10 * std::sqrt(nd);
        if (!(norm_c > rank_tolerance)) return result;
        for (double& value : qc) value /= norm_c;
        double cross = 0.0;
        for (int pass = 0; pass < 2; ++pass) {
            double correction = 0.0;
            for (std::size_t i = 0; i < n; ++i) correction += qc[i] * qs[i];
            cross += correction;
            for (std::size_t i = 0; i < n; ++i) qs[i] -= correction * qc[i];
        }
        double norm_s = 0.0;
        for (double value : qs) norm_s += value * value;
        norm_s = std::sqrt(norm_s);
        if (!(norm_s > rank_tolerance)) return result;
        for (double& value : qs) value /= norm_s;
        double rhs_c = 0.0, rhs_s = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            rhs_c += qc[i] * values[i];
            rhs_s += qs[i] * values[i];
        }
        result.sine = rhs_s / norm_s;
        result.cosine = (rhs_c - cross * result.sine) / norm_c;
        result.offset = -mean_c * result.cosine - mean_s * result.sine;
        result.sse = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double residual = values[i] - rhs_c * qc[i] - rhs_s * qs[i];
            result.sse += residual * residual;
        }
        if (!std::isfinite(result.sse) || !std::isfinite(result.offset) ||
            !std::isfinite(result.cosine) || !std::isfinite(result.sine)) {
            result.sse = std::numeric_limits<double>::infinity();
        }
        return result;
    }
};

inline Candidate refine(const Projector& project, Candidate start,
                        double w_min, double w_max, double a_max) {
    auto evaluate = [&](double w, double a) {
        return project(std::clamp(w, w_min, w_max), std::clamp(a, 0.0, a_max));
    };
    const double w_step = std::min(0.4, (w_max - w_min) / 4.0);
    std::array<Candidate, 3> simplex = {
        start,
        evaluate(start.w + ((start.w + w_step <= w_max) ? w_step : -w_step), start.a),
        evaluate(start.w, start.a + ((start.a + 0.25 <= a_max) ? 0.25 : -0.25))
    };
    // Bounded, multi-start Nelder-Mead on the two nonlinear parameters only.
    // Linear offset and quadratures are re-solved at every evaluation.
    for (int iteration = 0; iteration < 220; ++iteration) {
        std::sort(simplex.begin(), simplex.end(), [](const Candidate& x, const Candidate& y) {
            return x.sse < y.sse;
        });
        double spread = 0.0;
        for (std::size_t i = 1; i < simplex.size(); ++i) {
            spread = std::max(spread, std::abs(simplex[i].w - simplex[0].w));
            spread = std::max(spread, std::abs(simplex[i].a - simplex[0].a));
        }
        if (spread < 1e-9) break;
        const double mw = (simplex[0].w + simplex[1].w) / 2.0;
        const double ma = (simplex[0].a + simplex[1].a) / 2.0;
        Candidate reflected = evaluate(2.0 * mw - simplex[2].w, 2.0 * ma - simplex[2].a);
        if (reflected.sse < simplex[0].sse) {
            Candidate expanded = evaluate(3.0 * mw - 2.0 * simplex[2].w,
                                          3.0 * ma - 2.0 * simplex[2].a);
            simplex[2] = (expanded.sse < reflected.sse) ? expanded : reflected;
        } else if (reflected.sse < simplex[1].sse) {
            simplex[2] = reflected;
        } else {
            const bool outside = reflected.sse < simplex[2].sse;
            const Candidate& target = outside ? reflected : simplex[2];
            Candidate contracted = evaluate((mw + target.w) / 2.0, (ma + target.a) / 2.0);
            if (contracted.sse < target.sse) {
                simplex[2] = contracted;
            } else {
                for (std::size_t i = 1; i < simplex.size(); ++i) {
                    simplex[i] = evaluate((simplex[0].w + simplex[i].w) / 2.0,
                                          (simplex[0].a + simplex[i].a) / 2.0);
                }
            }
        }
    }
    Candidate best = *std::min_element(simplex.begin(), simplex.end(),
        [](const Candidate& x, const Candidate& y) { return x.sse < y.sse; });
    // A clipped simplex can collapse onto alpha=0. An independent pattern
    // polish reopens that direction instead of treating collapse as evidence
    // for zero damping.
    double dw = 0.08, da = 0.08;
    for (int iteration = 0; iteration < 180 && std::max(dw, da) > 1e-9; ++iteration) {
        Candidate origin = best;
        for (int iw = -1; iw <= 1; ++iw) {
            for (int ia = -1; ia <= 1; ++ia) {
                if (iw == 0 && ia == 0) continue;
                Candidate candidate = evaluate(origin.w + iw * dw, origin.a + ia * da);
                if (candidate.sse < best.sse) best = candidate;
            }
        }
        if (!(best.sse < origin.sse)) { dw *= 0.5; da *= 0.5; }
    }
    return best;
}
} // namespace detail

// Fits C + exp(-alpha*t) [B*cos(omega*t) + D*sin(omega*t)]. Input must
// contain >=8 finite, strictly increasing, uniformly spaced samples. The
// search covers 0.25 cycle/window <= omega*T/(2*pi) <= 0.495*T/dt and
// 0 <= alpha*T <= 20. The upper frequency is 99% of temporal Nyquist,
// NOT a grid-velocity or thermodynamic sound-speed constraint. The lower
// frequency and damping bounds are explicit finite numerical search limits.
// The coarse grid has <=1/6 cycle/window spacing, followed by independent
// refinements of up to 12 separated candidates. This is a deterministic
// numerical search, not a proof of global optimality or mode identification.
// valid means a finite numerical fit; callers must separately gate R^2,
// >=3 observed fitted cycles, boundary solutions and cross-window stability.
inline FitResult fit(const std::vector<double>& t, const std::vector<double>& y) {
    FitResult output;
    const std::size_t n = t.size();
    if (n < 8 || y.size() != n) return output;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(t[i]) || !std::isfinite(y[i])) return output;
    }
    const double duration = t.back() - t.front();
    const double dt = duration / static_cast<double>(n - 1);
    if (!(duration > 0.0) || !std::isfinite(duration) || !(dt > 0.0)) return output;
    for (std::size_t i = 1; i < n; ++i) {
        const double step = t[i] - t[i - 1];
        if (!(step > 0.0) || std::abs(step - dt) > 1e-8 * dt) return output;
    }
    double mean = 0.0;
    const double nd = static_cast<double>(n);
    for (double value : y) mean += value / nd;
    if (!std::isfinite(mean)) return output;
    double scale = 0.0, max_abs = 0.0;
    for (double value : y) {
        scale = std::max(scale, std::abs(value - mean));
        max_abs = std::max(max_abs, std::abs(value));
    }
    if (!std::isfinite(scale) || !(scale > 8.0 * std::numeric_limits<double>::epsilon() * max_abs)) {
        return output; // A constant (or numerically unresolved) signal has no frequency.
    }
    std::vector<double> times(n), values(n);
    double sst = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        times[i] = (t[i] - t.front()) / duration;
        values[i] = (y[i] - mean) / scale;
        sst += values[i] * values[i];
    }
    const detail::Projector project{times, values};
    const double w_min = 0.5 * detail::pi;
    const double w_max = 0.99 * detail::pi * static_cast<double>(n - 1);
    constexpr double a_max = 20.0;
    const int intervals = static_cast<int>(std::ceil((w_max - w_min) / (detail::pi / 3.0)));
    const std::array<double, 16> damping_grid = {
        0.0, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 20.0
    };
    std::vector<detail::Candidate> candidates;
    candidates.reserve(static_cast<std::size_t>(intervals + 1));
    for (int iw = 0; iw <= intervals; ++iw) {
        const double w = w_min + (w_max - w_min) * static_cast<double>(iw) / intervals;
        detail::Candidate best;
        for (double a : damping_grid) {
            const detail::Candidate candidate = project(w, a);
            if (candidate.sse < best.sse) best = candidate;
        }
        if (std::isfinite(best.sse)) candidates.push_back(best);
    }
    if (candidates.empty()) return output;
    std::sort(candidates.begin(), candidates.end(), [](const detail::Candidate& x, const detail::Candidate& z) {
        return x.sse < z.sse;
    });
    detail::Candidate best = candidates.front();
    std::vector<double> started;
    for (const detail::Candidate& candidate : candidates) {
        bool separated = true;
        for (double prior : started) if (std::abs(prior - candidate.w) < 1.5) separated = false;
        if (!separated) continue;
        started.push_back(candidate.w);
        const detail::Candidate refined = detail::refine(project, candidate, w_min, w_max, a_max);
        if (refined.sse < best.sse) best = refined;
        if (started.size() == 12) break;
    }
    output.omega = best.w / duration;
    output.alpha = best.a / duration;
    output.offset = mean + scale * best.offset;
    output.cosine = scale * best.cosine;
    output.sine = scale * best.sine;
    output.r_squared = 1.0 - best.sse / sst;
    output.rmse = scale * std::sqrt(best.sse / nd);
    output.cycles = best.w / (2.0 * detail::pi);
    output.frequency_at_boundary = best.w - w_min < 1e-5 || w_max - best.w < 1e-5;
    output.damping_at_boundary = best.a < 1e-5 || a_max - best.a < 1e-5;
    output.valid = std::isfinite(output.omega) && std::isfinite(output.alpha) &&
        std::isfinite(output.offset) && std::isfinite(output.cosine) &&
        std::isfinite(output.sine) && std::isfinite(output.r_squared) && std::isfinite(output.rmse);
    return output;
}
} // namespace acoustic

#endif
