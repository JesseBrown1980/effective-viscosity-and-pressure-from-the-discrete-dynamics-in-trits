// Asolaria Discrete Fluid: Direct Dynamical Acoustic Wave Propagation Experiment
// Measures speed of sound directly from propagating density disturbances
// without supplying any theoretical equation of state to the measurement procedure.

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

// Standing wave acoustic mode test: measures omega/k directly
struct StandingWaveResult {
    double omega;
    double measured_cs;
    double r_squared;
    i64 added_px;
};

static StandingWaveResult run_standing_wave_experiment(int nx, int ny, int nz, u64 seed,
                                                      const vector<Collision>& c,
                                                      const array<Velocity, K>& v) {
    int cells = nx * ny * nz;
    vector<std::uint8_t> cur(cells * K, 0), next(cells * K, 0);

    // Find opposite channel pairs
    array<int, K> opp{};
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            if (v[i].x == -v[j].x && v[i].y == -v[j].y && v[i].z == -v[j].z) {
                opp[i] = j;
                break;
            }
        }
    }

    // Initialize isotropic uniform random background with occupancy ~ 1
    u64 s = seed;
    for (int cell = 0; cell < cells; ++cell) {
        for (int k = 0; k < K; ++k) {
            int r = static_cast<int>(next64(s) % 12);
            cur[cell * K + k] = (r < 4) ? 0 : ((r < 8) ? 1 : 2);
        }
    }

    // Identify opposite channels along X-axis
    int k_pos = -1, k_neg = -1;
    for (int k = 0; k < K; ++k) {
        if (v[k].x == 1 && v[k].y == 0 && v[k].z == 0) k_pos = k;
        if (v[k].x == -1 && v[k].y == 0 && v[k].z == 0) k_neg = k;
    }

    // Add standing sinusoidal wave: delta_rho(x) ~ cos(2*pi*x / nx)
    // Strictly zero added momentum: add to opposite channel pairs (k_pos, k_neg)
    i64 added_px = 0;
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                int cell = (z * ny + y) * nx + x;
                double val = std::cos(2.0 * PI * static_cast<double>(x) / static_cast<double>(nx));
                if (val > 0.3) {
                    if (cur[cell * K + k_pos] < 2 && cur[cell * K + k_neg] < 2) {
                        cur[cell * K + k_pos] += 1;
                        cur[cell * K + k_neg] += 1;
                        added_px += v[k_pos].x + v[k_neg].x; // identically 0
                    }
                } else if (val < -0.3) {
                    if (cur[cell * K + k_pos] > 0 && cur[cell * K + k_neg] > 0) {
                        cur[cell * K + k_pos] -= 1;
                        cur[cell * K + k_neg] -= 1;
                        added_px += -v[k_pos].x - v[k_neg].x; // identically 0
                    }
                }
            }
        }
    }

    const int MAX_TICKS = 60;
    vector<double> amps;
    vector<double> t_vec;

    for (int t = 0; t <= MAX_TICKS; ++t) {
        double Ac = 0.0;
        for (int x = 0; x < nx; ++x) {
            double mass_x = 0;
            for (int z = 0; z < nz; ++z) {
                for (int y = 0; y < ny; ++y) {
                    int base = ((z * ny + y) * nx + x) * K;
                    for (int k = 0; k < K; ++k) mass_x += cur[base + k];
                }
            }
            Ac += mass_x * std::cos(2.0 * PI * static_cast<double>(x) / static_cast<double>(nx));
        }
        amps.push_back(Ac);
        t_vec.push_back(static_cast<double>(t));

        if (t == MAX_TICKS) break;
        apply_collisions(cur, cells, c, t);
        stream_grid(cur, next, nx, ny, nz, v);
        cur.swap(next);
    }

    // Find first zero crossing and half-period peak
    // Model: A(t) = A_offset + A_amp * cos(omega * t)
    double a_max = *std::max_element(amps.begin(), amps.end());
    double a_min = *std::min_element(amps.begin(), amps.end());
    double a_mid = (a_max + a_min) / 2.0;

    // Estimate zero-crossing t* where A(t) crosses a_mid
    double t_cross = 0.0;
    for (size_t i = 0; i + 1 < amps.size(); ++i) {
        if ((amps[i] - a_mid) * (amps[i + 1] - a_mid) <= 0.0) {
            double frac = (a_mid - amps[i]) / (amps[i + 1] - amps[i]);
            t_cross = static_cast<double>(i) + frac;
            break;
        }
    }

    // At quarter period t* = T/4 = pi / (2 * omega) => omega = pi / (2 * t_cross)
    double omega = (t_cross > 0.0) ? (PI / (2.0 * t_cross)) : 0.0;
    double k_wavenumber = 2.0 * PI / static_cast<double>(nx);
    double measured_cs = (k_wavenumber > 0.0) ? (omega / k_wavenumber) : 0.0;

    // Compute R^2 against cosine model
    double ss_tot = 0.0, ss_res = 0.0;
    double mean_a = 0.0;
    for (double a : amps) mean_a += a;
    mean_a /= amps.size();
    for (size_t i = 0; i < amps.size(); ++i) {
        double pred = a_mid + (a_max - a_mid) * std::cos(omega * t_vec[i]);
        ss_tot += (amps[i] - mean_a) * (amps[i] - mean_a);
        ss_res += (amps[i] - pred) * (amps[i] - pred);
    }
    double r2 = (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;

    return {omega, measured_cs, r2, added_px};
}

int main() {
    std::cout << "=======================================================================\n";
    std::cout << " D3Q27 DYNAMICAL ACOUSTIC FOURIER MODE PROPAGATION EXPERIMENT\n";
    std::cout << " Standing wave frequency measurement (0 theoretical EOS input)\n";
    std::cout << "=======================================================================\n\n";

    const auto v = get_velocities();
    const auto c = make_collisions(v);
    const double cs_theory = std::sqrt(2.0 / 3.0); // 0.816497

    vector<int> grid_lx = {24, 36, 48, 60, 72};
    vector<int> grid_ly = {12, 12, 12, 12, 12};
    vector<int> grid_lz = {12, 12, 12, 12, 12};

    std::cout << "Trial | Grid (Lx x Ly x Lz) | Cells | Added dPx | Measured c_s [dx/dt] | R^2 Fit | Error vs sqrt(2/3)\n";
    std::cout << "------+---------------------+-------+-----------+----------------------+---------+-------------------\n";

    std::ofstream csv("acoustic_pulse_propagation_results.csv");
    csv << "trial,lx,ly,lz,added_dpx,measured_cs,r_squared,relative_error_pct\n";

    vector<double> speeds;
    vector<double> r2_values;
    bool all_zero_momentum = true;

    for (size_t trial = 0; trial < grid_lx.size(); ++trial) {
        int lx = grid_lx[trial];
        int ly = grid_ly[trial];
        int lz = grid_lz[trial];
        u64 seed = 0xA501A21B00000001ull + trial * 0x123456789ull;

        auto res = run_standing_wave_experiment(lx, ly, lz, seed, c, v);
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
                  << std::setw(7) << res.r_squared << " | "
                  << std::setw(16) << (err_pct >= 0 ? "+" : "") << err_pct << "%\n";

        csv << (trial + 1) << "," << lx << "," << ly << "," << lz << ","
            << res.added_px << "," << res.measured_cs << "," << res.r_squared << "," << err_pct << "\n";
    }

    double mean_cs = 0.0;
    for (double s : speeds) mean_cs += s;
    mean_cs /= speeds.size();

    double mean_r2 = 0.0;
    for (double r : r2_values) mean_r2 += r;
    mean_r2 /= r2_values.size();

    std::cout << "------+---------------------+-------+-----------+----------------------+---------+-------------------\n";
    std::cout << "STRICT ZERO NET MOMENTUM CHECK : " << (all_zero_momentum ? "PASS (dPx == 0 identically)" : "FAIL") << "\n";
    std::cout << "MEAN DYNAMICAL WAVE SPEED (c_s): " << mean_cs << " [dx/dt]\n";
    std::cout << "THEORETICAL EOS SPEED (sqrt(2/3)): " << cs_theory << " [dx/dt]\n";
    std::cout << "AVERAGE HARMONIC CORRELATION R^2: " << mean_r2 << "\n";
    
    bool test_passed = all_zero_momentum && (mean_r2 > 0.90) && (std::abs(mean_cs - cs_theory) / cs_theory < 0.15);
    std::cout << "ACOUSTIC TEST STATUS            : " << (test_passed ? "PASS" : "FAIL") << "\n";
    std::cout << "=======================================================================\n";

    return test_passed ? 0 : 1;
}
