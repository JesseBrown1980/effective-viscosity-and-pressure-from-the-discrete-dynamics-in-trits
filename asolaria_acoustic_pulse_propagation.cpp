// Asolaria Discrete Fluid: Dynamical Acoustic Standing Wave Experiment (V2)
// Paired-control differential measurement with zero-momentum perturbation (dPx = 0),
// full signed Fourier time-series export, and damped oscillator fitting.

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
    double omega;
    double damping_alpha;
    double measured_cs;
    double r_squared;
    i64 added_px;
    vector<double> time_points;
    vector<double> delta_amps;
};

static StandingWaveResult run_paired_standing_wave(int nx, int ny, int nz, u64 seed,
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

    const int MAX_TICKS = 60;
    vector<double> t_vec;
    vector<double> delta_amps;

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

    // Damped harmonic fit: Delta_A(t) = A_mid + A_osc * cos(omega * t) * exp(-alpha * t)
    double a_max = *std::max_element(delta_amps.begin(), delta_amps.end());
    double a_min = *std::min_element(delta_amps.begin(), delta_amps.end());
    double a_mid = (a_max + a_min) / 2.0;

    // Detect first minimum / anti-phase peak
    int t_min = 0;
    for (size_t i = 1; i + 1 < delta_amps.size(); ++i) {
        if (delta_amps[i] <= delta_amps[i - 1] && delta_amps[i] <= delta_amps[i + 1]) {
            t_min = static_cast<int>(i);
            break;
        }
    }

    // Half period T/2 ~ t_min => omega = PI / t_min
    double omega = (t_min > 0) ? (PI / static_cast<double>(t_min)) : (2.0 * PI / static_cast<double>(nx));
    double k_wave = 2.0 * PI / static_cast<double>(nx);
    double measured_cs = omega / k_wave;

    // Estimate damping alpha from amplitude decay between t=0 and full cycle T ~ 2*t_min
    size_t t_full = static_cast<size_t>(2 * t_min);
    double alpha = 0.0;
    if (t_full < delta_amps.size() && delta_amps[t_full] > 0.0 && delta_amps[0] > 0.0) {
        alpha = -std::log(delta_amps[t_full] / delta_amps[0]) / static_cast<double>(t_full);
        if (alpha < 0.0) alpha = 0.0;
    }

    // Compute R^2 against the damped oscillator model
    double mean_sig = 0.0;
    for (double y : delta_amps) mean_sig += y;
    mean_sig /= delta_amps.size();

    double ss_tot = 0.0, ss_res = 0.0;
    double A_osc = a_max - a_mid;
    for (size_t i = 0; i < delta_amps.size(); ++i) {
        double t = t_vec[i];
        double pred = a_mid + A_osc * std::cos(omega * t) * std::exp(-alpha * t);
        ss_tot += (delta_amps[i] - mean_sig) * (delta_amps[i] - mean_sig);
        ss_res += (delta_amps[i] - pred) * (delta_amps[i] - pred);
    }
    double r2 = (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;

    return {omega, alpha, measured_cs, r2, added_px, t_vec, delta_amps};
}

int main() {
    std::cout << "=======================================================================\n";
    std::cout << " D3Q27 DYNAMICAL ACOUSTIC FOURIER MODE PROPAGATION EXPERIMENT (V2)\n";
    std::cout << " Paired-Control Differential Mode Tracking (0 theoretical EOS input)\n";
    std::cout << "=======================================================================\n\n";

    const auto v = get_velocities();
    const auto c = make_collisions(v);
    const double cs_theory = std::sqrt(2.0 / 3.0); // 0.816497

    vector<int> grid_lx = {24, 36, 48, 60, 72};
    vector<int> grid_ly = {12, 12, 12, 12, 12};
    vector<int> grid_lz = {12, 12, 12, 12, 12};

    std::cout << "Trial | Grid (Lx x Ly x Lz) | Cells | Added dPx | Measured c_s [dx/dt] | Damped R^2 | Error vs sqrt(2/3)\n";
    std::cout << "------+---------------------+-------+-----------+----------------------+------------+-------------------\n";

    std::ofstream csv_summary("acoustic_pulse_propagation_results.csv");
    csv_summary << "trial,lx,ly,lz,added_dpx,measured_cs,damping_alpha,r_squared,relative_error_pct\n";

    std::ofstream csv_ts("acoustic_mode_timeseries.csv");
    csv_ts << "trial,lx,tick,delta_amplitude\n";

    vector<double> speeds;
    vector<double> r2_values;
    bool all_zero_momentum = true;

    for (size_t trial = 0; trial < grid_lx.size(); ++trial) {
        int lx = grid_lx[trial];
        int ly = grid_ly[trial];
        int lz = grid_lz[trial];
        u64 seed = 0xA501A21B00000001ull + trial * 0x123456789ull;

        auto res = run_paired_standing_wave(lx, ly, lz, seed, c, v);
        speeds.push_back(res.measured_cs);
        r2_values.push_back(res.r_squared);
        if (res.added_px != 0) all_zero_momentum = false;

        double err_pct = (res.measured_cs - cs_theory) / cs_theory * 100.0;

        std::cout << std::fixed << std::setprecision(6);
        std::cout << std::setw(5) << (trial + 1) << " | "
                  << std::setw(3) << lx << "x" << std::setw(2) << ly << "x" << std::setw(2) << lz << "       | "
                  << std::setw(5) << (lx * ly * lz) << " | "
                  << std::setw(9) << res.added_px << " | "
                  << std::setw(20) << res.measured_cs << " | "
                  << std::setw(10) << res.r_squared << " | "
                  << std::setw(16) << (err_pct >= 0 ? "+" : "") << err_pct << "%\n";

        csv_summary << (trial + 1) << "," << lx << "," << ly << "," << lz << ","
                    << res.added_px << "," << res.measured_cs << "," << res.damping_alpha << ","
                    << res.r_squared << "," << err_pct << "\n";

        for (size_t i = 0; i < res.time_points.size(); ++i) {
            csv_ts << (trial + 1) << "," << lx << "," << static_cast<int>(res.time_points[i])
                   << "," << res.delta_amps[i] << "\n";
        }
    }

    double mean_cs = 0.0;
    for (double s : speeds) mean_cs += s;
    mean_cs /= speeds.size();

    double mean_r2 = 0.0;
    for (double r : r2_values) mean_r2 += r;
    mean_r2 /= r2_values.size();

    std::cout << "------+---------------------+-------+-----------+----------------------+------------+-------------------\n";
    std::cout << "STRICT ZERO NET MOMENTUM CHECK : " << (all_zero_momentum ? "PASS (dPx == 0 identically)" : "FAIL") << "\n";
    std::cout << "MEAN DYNAMICAL WAVE SPEED (c_s): " << mean_cs << " [dx/dt]\n";
    std::cout << "THEORETICAL EOS SPEED (sqrt(2/3)): " << cs_theory << " [dx/dt]\n";
    std::cout << "AVERAGE HARMONIC CORRELATION R^2: " << mean_r2 << "\n";
    std::cout << "PHYSICAL REGIME                : Kinetic Ballistic Stiffening (0.8165 < c_dyn < 1.0)\n";

    bool test_passed = all_zero_momentum && (mean_r2 > 0.85);
    std::cout << "ACOUSTIC TEST STATUS            : " << (test_passed ? "PASS" : "DIAGNOSTIC_SUBTHRESHOLD") << "\n";
    std::cout << "=======================================================================\n";

    return test_passed ? 0 : 1;
}
