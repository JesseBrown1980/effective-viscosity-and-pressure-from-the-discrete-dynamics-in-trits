#include "acoustic_damped_fit.h"
#include <cstdlib>
#include <iostream>
#include <string>

namespace {
int failures = 0;
int checks = 0;

void check(bool condition, const std::string& name) {
    ++checks;
    if (!condition) { ++failures; std::cerr << "FAIL: " << name << '\n'; }
}

struct Signal { std::vector<double> t, y; };

Signal make_signal(int count, double dt, double origin, double omega, double alpha,
                   double offset, double cosine, double sine, double noise = 0.0) {
    Signal signal;
    for (int i = 0; i < count; ++i) {
        const double relative = static_cast<double>(i) * dt;
        signal.t.push_back(origin + relative);
        signal.y.push_back(offset + std::exp(-alpha * relative) *
            (cosine * std::cos(omega * relative) + sine * std::sin(omega * relative)) +
            noise * (std::sin(1.731 * relative + 0.3) + 0.4 * std::cos(2.313 * relative)));
    }
    return signal;
}

void exact_recovery(const std::string& name, int count, double dt, double origin,
                    double omega, double alpha, double offset, double cosine, double sine) {
    const auto signal = make_signal(count, dt, origin, omega, alpha, offset, cosine, sine);
    const auto result = acoustic::fit(signal.t, signal.y);
    check(result.valid, name + " finite fit");
    check(std::abs(result.omega - omega) < 1e-6 * omega, name + " free frequency");
    check(std::abs(result.alpha - alpha) < 1e-7, name + " damping");
    check(std::abs(result.offset - offset) < 1e-5, name + " offset");
    check(std::abs(result.cosine - cosine) < 1e-5, name + " cosine quadrature");
    check(std::abs(result.sine - sine) < 1e-5, name + " sine quadrature");
    check(result.r_squared > 1.0 - 1e-10, name + " residual quality");
    double max_error = 0.0;
    for (std::size_t i = 0; i < signal.t.size(); ++i) {
        max_error = std::max(max_error, std::abs(acoustic::predict(result, signal.t[i] - origin) - signal.y[i]));
    }
    check(max_error < 1e-5, name + " prediction and time origin");
}
} // namespace

int main() {
    const double pi = 3.14159265358979323846;
    exact_recovery("off-grid period 81.37", 370, 1.0, 137.0,
                   2.0 * pi / 81.37, 0.0043, 12.7, 4.2, -2.3);
    exact_recovery("high temporal frequency", 150, 0.5, -8.0,
                   2.2173, 0.0183, -1.4, -3.2, 2.8);
    exact_recovery("undamped interior frequency", 150, 1.0, 0.0,
                   0.35127, 0.0, 8.2, 4.5, 1.7);

    auto noise = make_signal(340, 1.0, 0.0, 2.0 * pi / 79.71, 0.0037, 2.4, 3.0, -1.2, 0.025);
    auto noisy_fit = acoustic::fit(noise.t, noise.y);
    check(noisy_fit.valid && std::abs(noisy_fit.omega - 2.0 * pi / 79.71) < 2e-4,
          "deterministic residual disturbance frequency recovery");
    check(std::abs(noisy_fit.alpha - 0.0037) < 2e-4, "disturbed signal damping recovery");
    check(noisy_fit.r_squared > 0.995 && noisy_fit.r_squared < 1.0, "disturbance remains in residual");
    check(noisy_fit.cycles >= 3.0 && !noisy_fit.frequency_at_boundary, "multi-period coverage reported");

    auto short_signal = make_signal(55, 1.0, 0.0, 2.0 * pi / 80.3, 0.004, 2.0, 3.0, -0.8);
    auto short_fit = acoustic::fit(short_signal.t, short_signal.y);
    check(short_fit.valid, "short trajectory may yield a numerical fit");
    check(short_fit.cycles < 3.0, "short trajectory cannot pass caller's three-cycle gate");
    check(short_fit.r_squared > 0.999, "excellent R-squared does not imply adequate coverage");

    auto undamped = make_signal(160, 1.0, 0.0, 0.2319, 0.0, 4.0, 2.0, 1.0);
    auto undamped_fit = acoustic::fit(undamped.t, undamped.y);
    check(undamped_fit.valid && undamped_fit.damping_at_boundary,
          "zero-damping boundary is explicitly flagged");

    auto lower_bound = make_signal(130, 1.0, 0.0, 0.5 * pi / 129.0, 0.01, 4.0, 2.0, 1.0);
    auto lower_fit = acoustic::fit(lower_bound.t, lower_bound.y);
    check(lower_fit.valid && lower_fit.frequency_at_boundary, "lower-frequency search boundary flagged");

    auto flat = make_signal(30, 1.0, 0.0, 0.2, 0.01, 9.0, 0.0, 0.0);
    check(!acoustic::fit(flat.t, flat.y).valid, "constant signal has no fitted frequency");
    flat.y.assign(flat.y.size(), 0.0);
    check(!acoustic::fit(flat.t, flat.y).valid, "zero signal rejected");
    flat.y[4] = std::numeric_limits<double>::quiet_NaN();
    check(!acoustic::fit(flat.t, flat.y).valid, "NaN observations rejected");
    auto invalid = make_signal(30, 1.0, 0.0, 0.2, 0.01, 1.0, 2.0, 1.0);
    invalid.t[4] = invalid.t[3];
    check(!acoustic::fit(invalid.t, invalid.y).valid, "duplicate times rejected");
    invalid.t[4] = 4.1;
    check(!acoustic::fit(invalid.t, invalid.y).valid, "irregular times rejected by documented contract");
    invalid.t[4] = std::numeric_limits<double>::infinity();
    check(!acoustic::fit(invalid.t, invalid.y).valid, "infinite time rejected");
    check(!acoustic::fit({0, 1, 2, 3}, {3, 2, 1, 2}).valid, "too few samples rejected");
    check(!acoustic::fit({0, 1, 2}, {3, 2}).valid, "mismatched input lengths rejected");

    std::cout << "Acoustic fitter numerical checks: " << (checks - failures) << '/' << checks
              << " passed\n";
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
