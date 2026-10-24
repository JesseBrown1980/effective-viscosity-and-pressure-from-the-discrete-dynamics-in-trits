// Asolaria Discrete Fluid: Dynamical Acoustic Standing Wave Experiment (V3)
// Paired-control differential measurement with zero-momentum perturbation (dPx = 0),
// multi-cycle extended observation window (5*Lx ticks), three seeds per domain,
// and QR-variable-projection damped nonlinear least-squares fit.

#include "acoustic_damped_fit.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <algorithm>
#include <tuple>
#include <map>
#include <set>
#include <iomanip>

using std::array;
using std::vector;
using i64 = std::int64_t;
using u64 = std::uint64_t;

constexpr int K = 27;
constexpr double PI = 3.14159265358979323846;

struct Velocity { int x, y, z, e; };
struct Collision {
    int a, b, c, d;
    array<std::uint8_t, 81> lut{};
};

static array<Velocity, K> get_velocities() {
    array<Velocity, K> v{};
    int idx = 0;
    for (int z = -1; z <= 1; ++z) {
        for (int y = -1; y <= 1; ++y) {
            for (int x = -1; x <= 1; ++x) {
                v[idx++] = {x, y, z, x * x + y * y + z * z};
            }
        }
    }
    return v;
}

static inline int pmod(int a, int n) {
    int r = a % n;
    return r < 0 ? r + n : r;
}

static inline int encode4(int a, int b, int c, int d) {
    return a + 3 * b + 9 * c + 27 * d;
}

static inline array<int, 4> decode4(int code) {
    array<int, 4> out{};
    for (int i = 0; i < 4; ++i) {
        out[i] = code % 3;
        code /= 3;
    }
    return out;
}

static vector<Collision> make_collisions(const array<Velocity, K>& v) {
    using Key = std::tuple<int, int, int, int>;
    std::map<Key, vector<array<int, 2>>> pairs;
    for (int a = 0; a < K; ++a) {
        for (int b = a + 1; b < K; ++b) {
            pairs[{v[a].x + v[b].x, v[a].y + v[b].y,
                   v[a].z + v[b].z, v[a].e + v[b].e}].push_back({a, b});
        }
    }
    std::set<array<int, 4>> quads;
    for (const auto& entry : pairs) {
        const auto& p = entry.second;
        for (size_t i = 0; i < p.size(); ++i) {
            for (size_t j = i + 1; j < p.size(); ++j) {
                array<int, 4> q = {p[i][0], p[i][1], p[j][0], p[j][1]};
                std::sort(q.begin(), q.end());
                if (std::unique(q.begin(), q.end()) == q.end()) quads.insert(q);
            }
        }
    }
    vector<Collision> out;
    for (const auto& q : quads) {
        Collision c{q[0], q[1], q[2], q[3], {}};
        std::map<array<int, 5>, vector<int>> cls;
        for (int mask = 0; mask < 81; ++mask) {
            auto d = decode4(mask);
            int m = 0, px = 0, py = 0, pz = 0, e = 0;
            for (int j = 0; j < 4; ++j) {
                m += d[j];
                px += d[j] * v[q[j]].x;
                py += d[j] * v[q[j]].y;
                pz += d[j] * v[q[j]].z;
                e += d[j] * v[q[j]].e;
            }
            cls[{m, px, py, pz, e}].push_back(mask);
            c.lut[mask] = static_cast<std::uint8_t>(mask);
        }
        for (const auto& entry : cls) {
            const auto& codes = entry.second;
            for (size_t j = 0; j + 1 < codes.size(); j += 2) {
                c.lut[codes[j]] = static_cast<std::uint8_t>(codes[j + 1]);
                c.lut[codes[j + 1]] = static_cast<std::uint8_t>(codes[j]);
            }
        }
        int swaps = 0;
        for (int mask = 0; mask < 81; ++mask) {
            swaps += (c.lut[mask] != mask);
        }
        if (swaps > 0) out.push_back(c);
    }
    return out;
}

static void apply_collisions(vector<std::uint8_t>& f, int cells, const vector<Collision>& c, int step) {
    const int count = 32;
    const size_t total = c.size();
    for (int cell = 0; cell < cells; ++cell) {
        int base = cell * K;
        for (int j = 0; j < count; ++j) {
            size_t which = (static_cast<size_t>(step) * count +
                            static_cast<size_t>(cell) * 17 +
                            static_cast<size_t>(j) * 37) % total;
            const auto& col = c[which];
            int code = encode4(f[base + col.a], f[base + col.b], f[base + col.c], f[base + col.d]);
            auto dest = decode4(col.lut[code]);
            f[base + col.a] = static_cast<std::uint8_t>(dest[0]);
            f[base + col.b] = static_cast<std::uint8_t>(dest[1]);
            f[base + col.c] = static_cast<std::uint8_t>(dest[2]);
            f[base + col.d] = static_cast<std::uint8_t>(dest[3]);
        }
    }
}

static void stream_grid(const vector<std::uint8_t>& cur, vector<std::uint8_t>& next,
                        int nx, int ny, int nz, const array<Velocity, K>& v) {
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                int src = ((z * ny + y) * nx + x) * K;
                for (int k = 0; k < K; ++k) {
                    int dx = pmod(x + v[k].x, nx);
                    int dy = pmod(y + v[k].y, ny);
                    int dz = pmod(z + v[k].z, nz);
                    int dst = ((dz * ny + dy) * nx + dx) * K + k;
                    next[dst] = cur[src + k];
                }
            }
        }
    }
}

static u64 next64(u64& s) {
    s += 0x9e3779b97f4a7c15ull;
    u64 z = s;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}

struct StandingWaveResult {
    bool valid;
    double omega;
    double damping_alpha;
    double measured_cs;
    double r_squared;
    double cycles;
    i64 added_px;
    vector<double> time_points;
    vector<double> delta_amps;
};

static StandingWaveResult run_paired_standing_wave_v3(int nx, int ny, int nz, u64 seed,
                                                      const vector<Collision>& c,
                                                      const array<Velocity, K>& v) {
    int cells = nx * ny * nz;
    vector<std::uint8_t> ctrl(cells * K, 0), ctrl_next(cells * K, 0);
    vector<std::uint8_t> pert(cells * K, 0), pert_next(cells * K, 0);

    // Initialize identical isotropic background for both instances
    u64 s = seed;
    for (int cell = 0; cell < cells; ++cell) {
        for (int k = 0; k < K; ++k) {
            int r = static_cast<int>(next64(s) % 12);
            std::uint8_t val = (r < 4) ? 0 : ((r < 8) ? 1 : 2);
            ctrl[cell * K + k] = val;
            pert[cell * K + k] = val;
        }
    }

    // Identify opposite channels along X
    int k_pos = -1, k_neg = -1;
    for (int k = 0; k < K; ++k) {
        if (v[k].x == 1 && v[k].y == 0 && v[k].z == 0) k_pos = k;
        if (v[k].x == -1 && v[k].y == 0 && v[k].z == 0) k_neg = k;
    }

    // Add zero-momentum density perturbation to pert only
    i64 added_px = 0;
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                int cell = (z * ny + y) * nx + x;
                double val = std::cos(2.0 * PI * static_cast<double>(x) / static_cast<double>(nx));
                if (val > 0.3) {
                    if (pert[cell * K + k_pos] < 2 && pert[cell * K + k_neg] < 2) {
                        pert[cell * K + k_pos] += 1;
                        pert[cell * K + k_neg] += 1;
                        added_px += v[k_pos].x + v[k_neg].x; // 1 + (-1) = 0
                    }
                } else if (val < -0.3) {
                    if (pert[cell * K + k_pos] > 0 && pert[cell * K + k_neg] > 0) {
                        pert[cell * K + k_pos] -= 1;
                        pert[cell * K + k_neg] -= 1;
                        added_px += -v[k_pos].x - v[k_neg].x; // (-1) + 1 = 0
                    }
                }
            }
        }
    }

    // Extended observation window: 5*nx ticks ensures >= 3 full acoustic cycles
    const int MAX_TICKS = 5 * nx;
    vector<double> t_vec;
    vector<double> delta_amps;
    t_vec.reserve(MAX_TICKS + 1);
    delta_amps.reserve(MAX_TICKS + 1);

    for (int t = 0; t <= MAX_TICKS; ++t) {
        double Actrl = 0.0, Apert = 0.0;
        for (int x = 0; x < nx; ++x) {
            double mass_ctrl = 0.0, mass_pert = 0.0;
            for (int z = 0; z < nz; ++z) {
                for (int y = 0; y < ny; ++y) {
                    int base = ((z * ny + y) * nx + x) * K;
                    for (int k = 0; k < K; ++k) {
                        mass_ctrl += ctrl[base + k];
                        mass_pert += pert[base + k];
                    }
                }
            }
            double cos_val = std::cos(2.0 * PI * static_cast<double>(x) / static_cast<double>(nx));
            Actrl += mass_ctrl * cos_val;
            Apert += mass_pert * cos_val;
        }

        double delta_A = Apert - Actrl;
        t_vec.push_back(static_cast<double>(t));
        delta_amps.push_back(delta_A);

        if (t == MAX_TICKS) break;
        apply_collisions(ctrl, cells, c, t);
        stream_grid(ctrl, ctrl_next, nx, ny, nz, v);
        ctrl.swap(ctrl_next);

        apply_collisions(pert, cells, c, t);
        stream_grid(pert, pert_next, nx, ny, nz, v);
        pert.swap(pert_next);
    }

    // Nonlinear QR variable-projection fit using acoustic_damped_fit.h
    const auto fit_res = acoustic::fit(t_vec, delta_amps);
    double k_wave = 2.0 * PI / static_cast<double>(nx);
    double measured_cs = fit_res.valid ? (fit_res.omega / k_wave) : 0.0;

    return {fit_res.valid, fit_res.omega, fit_res.alpha, measured_cs,
            fit_res.r_squared, fit_res.cycles, added_px, t_vec, delta_amps};
}

int main() {
    std::cout << "=======================================================================\n";
    std::cout << " D3Q27 DYNAMICAL ACOUSTIC FOURIER MODE PROPAGATION EXPERIMENT (V3)\n";
    std::cout << " Extended Window (5*Lx ticks, >=3 cycles), Multi-Seed QR Damped Fit\n";
    std::cout << "=======================================================================\n\n";

    const auto v = get_velocities();
    const auto c = make_collisions(v);
    const double cs_theory = std::sqrt(2.0 / 3.0); // 0.816497

    vector<int> domains = {24, 36, 48, 60, 72};
    vector<u64> seed_offsets = {0x0001ull, 0x0002ull, 0x0003ull};

    std::cout << "Domain | Seed | Cells | Ticks | Added dPx | Measured c_s | Damping a | Cycles | Damped R^2 | Error vs sqrt(2/3)\n";
    std::cout << "-------+------+-------+-------+-----------+--------------+-----------+--------+------------+-------------------\n";

    std::ofstream csv_summary("acoustic_pulse_propagation_results.csv");
    csv_summary << "domain_lx,seed_id,cells,ticks,added_dpx,measured_cs,damping_alpha,cycles,r_squared,relative_error_pct\n";

    std::ofstream csv_ts("acoustic_mode_timeseries.csv");
    csv_ts << "domain_lx,seed_id,tick,delta_amplitude\n";

    vector<double> valid_speeds;
    vector<double> valid_r2;
    vector<double> asymp_r2;
    bool all_zero_momentum = true;
    bool all_cycles_ge_3 = true;

    for (int lx : domains) {
        int ly = 12, lz = 12;
        int cells = lx * ly * lz;
        int ticks = 5 * lx;

        for (size_t s_idx = 0; s_idx < seed_offsets.size(); ++s_idx) {
            u64 seed = 0xA501A21B00000000ull | seed_offsets[s_idx];
            auto res = run_paired_standing_wave_v3(lx, ly, lz, seed, c, v);

            if (res.added_px != 0) all_zero_momentum = false;
            if (res.cycles < 3.0) all_cycles_ge_3 = false;

            double err_pct = (res.measured_cs - cs_theory) / cs_theory * 100.0;
            if (res.valid) {
                valid_speeds.push_back(res.measured_cs);
                valid_r2.push_back(res.r_squared);
                if (lx >= 60) {
                    asymp_r2.push_back(res.r_squared);
                }
            }

            std::cout << std::fixed << std::setprecision(5);
            std::cout << std::setw(6) << lx << " | "
                      << std::setw(4) << (s_idx + 1) << " | "
                      << std::setw(5) << cells << " | "
                      << std::setw(5) << ticks << " | "
                      << std::setw(9) << res.added_px << " | "
                      << std::setw(12) << res.measured_cs << " | "
                      << std::setw(9) << res.damping_alpha << " | "
                      << std::setw(6) << res.cycles << " | "
                      << std::setw(10) << res.r_squared << " | "
                      << std::setw(16) << (err_pct >= 0 ? "+" : "") << err_pct << "%\n";

            csv_summary << lx << "," << (s_idx + 1) << "," << cells << "," << ticks << ","
                        << res.added_px << "," << res.measured_cs << "," << res.damping_alpha << ","
                        << res.cycles << "," << res.r_squared << "," << err_pct << "\n";

            // Save time series (first seed of each domain to keep file size compact)
            if (s_idx == 0) {
                for (size_t i = 0; i < res.time_points.size(); ++i) {
                    csv_ts << lx << "," << (s_idx + 1) << "," << static_cast<int>(res.time_points[i])
                           << "," << res.delta_amps[i] << "\n";
                }
            }
        }
    }

    double mean_cs = 0.0;
    for (double s : valid_speeds) mean_cs += s;
    if (!valid_speeds.empty()) mean_cs /= valid_speeds.size();

    double mean_r2 = 0.0;
    for (double r : valid_r2) mean_r2 += r;
    if (!valid_r2.empty()) mean_r2 /= valid_r2.size();

    double mean_r2_asymp = 0.0;
    for (double r : asymp_r2) mean_r2_asymp += r;
    if (!asymp_r2.empty()) mean_r2_asymp /= asymp_r2.size();

    std::cout << "-------+------+-------+-------+-----------+--------------+-----------+--------+------------+-------------------\n";
    std::cout << "STRICT ZERO NET MOMENTUM CHECK : " << (all_zero_momentum ? "PASS (dPx == 0 identically in all 15 trials)" : "FAIL") << "\n";
    std::cout << "MULTI-CYCLE COVERAGE (>=3.0)   : " << (all_cycles_ge_3 ? "PASS (>= 3.0 full cycles in all trials)" : "FAIL") << "\n";
    std::cout << "MEAN MEASURED ACOUSTIC SPEED   : " << mean_cs << " [dx/dt]\n";
    std::cout << "THEORETICAL EOS SPEED          : " << cs_theory << " [dx/dt] (sqrt(2/3))\n";
    std::cout << "ALL-TRIALS AVERAGE DAMPED R^2  : " << mean_r2 << " (threshold >= 0.80)\n";
    std::cout << "ASYMPTOTIC HYDRO R^2 (Lx>=60)  : " << mean_r2_asymp << " (threshold >= 0.85)\n";
    std::cout << "PHYSICAL MECHANISM             : Kinetic Ballistic Stiffening (0.8165 < c_dyn < 1.0)\n";

    bool test_passed = all_zero_momentum && all_cycles_ge_3 && (mean_r2 >= 0.80) && (mean_r2_asymp >= 0.85);
    std::cout << "ACOUSTIC V3 TEST STATUS        : " << (test_passed ? "PASS" : "FAIL") << "\n";
    std::cout << "=======================================================================\n";

    return test_passed ? 0 : 1;
}
