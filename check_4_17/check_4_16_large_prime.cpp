#include <algorithm>
#include <chrono>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Experimental direct-modulus checker for Finite arithmetic verification 4.16.
//
// The production computation is carried out in
//
//     F_p,  p = 2013265921 = 15 * 2^27 + 1,
//
// using primitive root 31.  Since p > 50000, there is no modular translation
// n = q p + s: throughout the finite rectangle q = 0 and s = n.  The program
// therefore reconstructs only
//
//     D_10(g,n) = [X^10] Q_{g,n}(X) / n!  (mod p).
//
// The formulas below are a direct-modulus specialization of the archived
// modular_low_truncation_table_large.cpp.  The large n-range is processed in
// overlapping blocks, so no 2000 x 50000 x 11 output table is stored or
// written to disk.

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;

class WallClock {
public:
    WallClock() : start_(std::chrono::steady_clock::now()) {}

    double seconds() const {
        return std::chrono::duration<double>(
                   std::chrono::steady_clock::now() - start_)
            .count();
    }

private:
    std::chrono::steady_clock::time_point start_;
};

std::string duration_string(double seconds) {
    const auto total = static_cast<u64>(seconds + 0.5);
    const u64 hours = total / 3600;
    const u64 minutes = (total % 3600) / 60;
    const u64 secs = total % 60;
    std::string result;
    if (hours) {
        result += std::to_string(hours) + "h ";
    }
    if (hours || minutes) {
        result += std::to_string(minutes) + "m ";
    }
    result += std::to_string(secs) + "s";
    return result;
}

std::size_t next_power_of_two(std::size_t required_length) {
    std::size_t result = 1;
    while (result < required_length) {
        if (result > std::numeric_limits<std::size_t>::max() / 2) {
            throw std::overflow_error("power-of-two length overflow");
        }
        result <<= 1;
    }
    return result;
}

std::string gibibytes(std::size_t bytes) {
    const double gib = static_cast<double>(bytes) / (1024.0 * 1024.0 * 1024.0);
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << gib << " GiB";
    return out.str();
}

template <u32 Modulus>
struct Field {
    static u32 normalize(std::int64_t value) {
        value %= static_cast<std::int64_t>(Modulus);
        if (value < 0) {
            value += Modulus;
        }
        return static_cast<u32>(value);
    }

    static u32 add(u32 a, u32 b) {
        const u64 sum = static_cast<u64>(a) + b;
        return static_cast<u32>(sum >= Modulus ? sum - Modulus : sum);
    }

    static u32 subtract(u32 a, u32 b) {
        return a >= b ? a - b : static_cast<u32>(static_cast<u64>(a) + Modulus - b);
    }

    static u32 multiply(u32 a, u32 b) {
        return static_cast<u32>(static_cast<u64>(a) * b % Modulus);
    }

    static u32 power(u32 base, u64 exponent) {
        u32 result = 1;
        while (exponent) {
            if (exponent & 1U) {
                result = multiply(result, base);
            }
            base = multiply(base, base);
            exponent >>= 1U;
        }
        return result;
    }

    static u32 inverse(u32 value) {
        if (value == 0) {
            throw std::runtime_error("attempted division by zero in finite field");
        }
        return power(value, Modulus - 2ULL);
    }
};

template <u32 Modulus, u32 PrimitiveRoot>
struct Ntt {
    using F = Field<Modulus>;

    static void transform(std::vector<u32>& values, bool inverse) {
        const std::size_t n = values.size();
        if (n == 0 || (n & (n - 1)) != 0) {
            throw std::runtime_error("NTT length is not a positive power of two");
        }
        if ((Modulus - 1ULL) % n != 0) {
            throw std::runtime_error("NTT length does not divide p-1");
        }

        for (std::size_t i = 1, j = 0; i < n; ++i) {
            std::size_t bit = n >> 1U;
            while (j & bit) {
                j ^= bit;
                bit >>= 1U;
            }
            j ^= bit;
            if (i < j) {
                std::swap(values[i], values[j]);
            }
        }

        for (std::size_t length = 2; length <= n; length <<= 1U) {
            u32 root = F::power(PrimitiveRoot, (Modulus - 1ULL) / length);
            if (inverse) {
                root = F::inverse(root);
            }
            const std::size_t half = length >> 1U;
            for (std::size_t start = 0; start < n; start += length) {
                u32 root_power = 1;
                for (std::size_t offset = 0; offset < half; ++offset) {
                    const u32 even = values[start + offset];
                    const u32 odd = F::multiply(values[start + offset + half], root_power);
                    values[start + offset] = F::add(even, odd);
                    values[start + offset + half] = F::subtract(even, odd);
                    root_power = F::multiply(root_power, root);
                }
            }
        }

        if (inverse) {
            const u32 inverse_n = F::inverse(static_cast<u32>(n % Modulus));
            for (u32& value : values) {
                value = F::multiply(value, inverse_n);
            }
        }
    }
};

std::vector<int> mobius_table(int maximum) {
    std::vector<int> mobius(maximum + 1, 0);
    std::vector<int> least_prime(maximum + 1, 0);
    std::vector<int> primes;
    mobius[1] = 1;
    for (int n = 2; n <= maximum; ++n) {
        if (least_prime[n] == 0) {
            least_prime[n] = n;
            primes.push_back(n);
            mobius[n] = -1;
        }
        for (int prime : primes) {
            if (prime > least_prime[n] || 1LL * n * prime > maximum) {
                break;
            }
            least_prime[n * prime] = prime;
            mobius[n * prime] = (n % prime == 0) ? 0 : -mobius[n];
        }
    }
    return mobius;
}

std::vector<std::vector<int>> divisor_table(int maximum) {
    std::vector<std::vector<int>> divisors(maximum + 1);
    for (int divisor = 1; divisor <= maximum; ++divisor) {
        for (int multiple = divisor; multiple <= maximum; multiple += divisor) {
            divisors[multiple].push_back(divisor);
        }
    }
    return divisors;
}

bool is_prime_trial_division(u32 number) {
    if (number < 2) {
        return false;
    }
    if (number % 2 == 0) {
        return number == 2;
    }
    for (u32 divisor = 3; static_cast<u64>(divisor) * divisor <= number; divisor += 2) {
        if (number % divisor == 0) {
            return false;
        }
    }
    return true;
}

template <u32 Modulus, u32 PrimitiveRoot>
void validate_field_and_ntt() {
    using F = Field<Modulus>;
    std::cout << "Validating modulus p=" << Modulus << " and primitive root "
              << PrimitiveRoot << " ... " << std::flush;
    if (!is_prime_trial_division(Modulus)) {
        throw std::runtime_error("the verification modulus is not prime");
    }

    // The production modulus has p-1 = 2^27 * 3 * 5; the smoke-test modulus
    // has p-1 = 2^13 * 5.  Testing the distinct prime factors certifies that
    // PrimitiveRoot generates the full multiplicative group.
    u32 remaining = Modulus - 1;
    std::vector<u32> prime_factors;
    for (u32 divisor = 2; static_cast<u64>(divisor) * divisor <= remaining; ++divisor) {
        if (remaining % divisor == 0) {
            prime_factors.push_back(divisor);
            while (remaining % divisor == 0) {
                remaining /= divisor;
            }
        }
    }
    if (remaining > 1) {
        prime_factors.push_back(remaining);
    }
    for (u32 factor : prime_factors) {
        if (F::power(PrimitiveRoot, (Modulus - 1ULL) / factor) == 1) {
            throw std::runtime_error("the stated NTT root is not primitive");
        }
    }

    for (u32 value : {1U, 2U, 3U, 17U, 1000U}) {
        if (F::multiply(value, F::inverse(value)) != 1) {
            throw std::runtime_error("finite-field inverse self-test failed");
        }
    }

    std::vector<u32> left = {1, 2, 3, 4, 5};
    std::vector<u32> right = {7, 11, 13, 17};
    std::vector<u32> expected(left.size() + right.size() - 1, 0);
    for (std::size_t i = 0; i < left.size(); ++i) {
        for (std::size_t j = 0; j < right.size(); ++j) {
            expected[i + j] = F::add(expected[i + j], F::multiply(left[i], right[j]));
        }
    }
    const std::size_t transform_length = next_power_of_two(expected.size());
    left.resize(transform_length, 0);
    right.resize(transform_length, 0);
    Ntt<Modulus, PrimitiveRoot>::transform(left, false);
    Ntt<Modulus, PrimitiveRoot>::transform(right, false);
    for (std::size_t i = 0; i < transform_length; ++i) {
        left[i] = F::multiply(left[i], right[i]);
    }
    Ntt<Modulus, PrimitiveRoot>::transform(left, true);
    if (!std::equal(expected.begin(), expected.end(), left.begin())) {
        throw std::runtime_error("NTT convolution self-test failed");
    }
    std::cout << "PASS\n";
}

struct RunParameters {
    int maximum_genus;
    int weight_degree;
    int n_count;
    int block_width;
    bool manuscript_domain;
    std::string unresolved_path;
};

template <u32 Modulus, u32 PrimitiveRoot>
class DirectCoefficientChecker {
public:
    using F = Field<Modulus>;
    using Transform = Ntt<Modulus, PrimitiveRoot>;

    explicit DirectCoefficientChecker(RunParameters parameters)
        : parameters_(std::move(parameters)),
          genus_bound_(parameters_.maximum_genus),
          degree_(parameters_.weight_degree),
          n_count_(parameters_.n_count),
          block_width_(parameters_.block_width),
          inverses_(std::max(2 * genus_bound_ + 2, n_count_) + 2, 0),
          mobius_(mobius_table(2 * genus_bound_)),
          divisors_(divisor_table(2 * genus_bound_)) {
        if (genus_bound_ < 1 || degree_ < 0 || n_count_ < 1 || block_width_ < 1) {
            throw std::invalid_argument("all run parameters must be positive");
        }
        if (block_width_ <= genus_bound_) {
            throw std::invalid_argument("block width must exceed the genus bound");
        }
        if (Modulus <= static_cast<u32>(std::max(2 * genus_bound_ + degree_ + 2,
                                                 n_count_ - 1))) {
            throw std::invalid_argument("modulus is not large enough for direct evaluation");
        }

        build_inverse_table();
        build_factorials();
        build_bernoulli_numbers();

        packed_base_ = next_power_of_two(
            static_cast<std::size_t>(block_width_) + 2ULL * genus_bound_);
        packed_transform_length_ = next_power_of_two(
            static_cast<std::size_t>(2 * genus_bound_ - 1) * packed_base_);
        if ((Modulus - 1ULL) % packed_transform_length_ != 0) {
            throw std::runtime_error("chosen block size needs an unsupported NTT length");
        }
    }

    u64 run() {
        const WallClock total_clock;
        print_run_header();

        build_r_series();
        build_h_series();
        build_cumulative_r_transforms();

        std::ofstream unresolved_output;
        if (!parameters_.unresolved_path.empty()) {
            unresolved_output.open(parameters_.unresolved_path);
            if (!unresolved_output) {
                throw std::runtime_error("cannot open unresolved-pair output file");
            }
            unresolved_output << "g\tn\tresidue_mod_" << Modulus << '\n';
        }

        std::vector<std::vector<u32>> polynomial_checkpoint(
            genus_bound_, std::vector<u32>(degree_ + 1, 0));
        for (int h = 0; h < genus_bound_; ++h) {
            polynomial_checkpoint[h][0] = 1;
        }
        int checkpoint_n = 0;

        u64 tested = 0;
        u64 zeros = 0;
        u64 checksum = 0;
        const int block_count = (n_count_ + block_width_ - 1) / block_width_;
        const WallClock block_phase_clock;

        for (int block_index = 0; block_index < block_count; ++block_index) {
            const int output_begin = block_index * block_width_;
            const int output_end = std::min(n_count_, output_begin + block_width_);
            const int input_begin = std::max(0, output_begin - genus_bound_);
            if (input_begin != checkpoint_n) {
                throw std::runtime_error("internal polynomial checkpoint mismatch");
            }

            const WallClock block_clock;
            std::cout << "\n[block " << (block_index + 1) << '/' << block_count << "] "
                      << "output n=" << output_begin << ".." << (output_end - 1)
                      << ", A-input n=" << input_begin << ".." << (output_end - 1)
                      << '\n' << std::flush;

            const int local_width = output_end - input_begin;
            auto compact_a = build_a_block(input_begin, output_end,
                                           polynomial_checkpoint, checkpoint_n);
            std::cout << "  constructed A coefficients in "
                      << duration_string(block_clock.seconds()) << '\n' << std::flush;

            std::vector<u32> accumulator(packed_transform_length_, 0);
            std::vector<u32> transform_buffer(packed_transform_length_, 0);

            for (int a_degree = 0; a_degree <= degree_; ++a_degree) {
                std::fill(transform_buffer.begin(), transform_buffer.end(), 0);
                const auto& compact = compact_a[a_degree];
                for (int h = 0; h < genus_bound_; ++h) {
                    const std::size_t compact_offset =
                        static_cast<std::size_t>(h) * local_width;
                    const std::size_t packed_offset =
                        static_cast<std::size_t>(h) * packed_base_;
                    std::copy_n(compact.data() + compact_offset,
                                local_width,
                                transform_buffer.data() + packed_offset);
                }

                Transform::transform(transform_buffer, false);
                const auto& r_transform = cumulative_r_transforms_[a_degree];
                for (std::size_t index = 0; index < packed_transform_length_; ++index) {
                    accumulator[index] = F::add(
                        accumulator[index],
                        F::multiply(transform_buffer[index], r_transform[index]));
                }
                std::cout << "  transformed A degree " << a_degree << '/' << degree_
                          << " (block elapsed " << duration_string(block_clock.seconds())
                          << ")\n" << std::flush;
            }

            compact_a.clear();
            compact_a.shrink_to_fit();
            transform_buffer.clear();
            transform_buffer.shrink_to_fit();

            std::cout << "  inverse transform ... " << std::flush;
            Transform::transform(accumulator, true);
            std::cout << "done\n" << std::flush;

            u64 block_tested = 0;
            u64 block_zeros = 0;
            for (int n = output_begin; n < output_end; ++n) {
                const std::size_t local_n = static_cast<std::size_t>(n - input_begin);
                for (int genus = 1; genus <= genus_bound_; ++genus) {
                    if (!belongs_to_checked_domain(genus, n)) {
                        continue;
                    }
                    u32 residue = accumulator[
                        static_cast<std::size_t>(genus - 1) * packed_base_ + local_n];
                    if (genus == 1 && n == 0) {
                        residue = F::subtract(residue, 1);
                    }
                    ++tested;
                    ++block_tested;
                    checksum = (checksum * 1000003ULL + residue) % 1000000007ULL;
                    if (residue == 0) {
                        ++zeros;
                        ++block_zeros;
                        if (unresolved_output) {
                            unresolved_output << genus << '\t' << n << '\t' << residue << '\n';
                        }
                    }
                }
            }

            const double elapsed_blocks = block_phase_clock.seconds();
            const double estimated_total_blocks =
                elapsed_blocks * block_count / static_cast<double>(block_index + 1);
            const double eta = std::max(0.0, estimated_total_blocks - elapsed_blocks);
            std::cout << "  block result: tested=" << block_tested
                      << ", zero residues=" << block_zeros
                      << ", block time=" << duration_string(block_clock.seconds()) << '\n'
                      << "  cumulative: tested=" << tested
                      << ", zero residues=" << zeros
                      << ", estimated remaining block time=" << duration_string(eta)
                      << '\n' << std::flush;
        }

        if (parameters_.manuscript_domain && tested != 99999950ULL) {
            throw std::runtime_error("domain count mismatch: expected 99,999,950 pairs");
        }

        std::cout << "\nSUMMARY\n"
                  << "  modulus       = " << Modulus << '\n'
                  << "  tested pairs  = " << tested << '\n'
                  << "  zero residues = " << zeros << '\n'
                  << "  checksum      = " << checksum << '\n'
                  << "  total time    = " << duration_string(total_clock.seconds()) << '\n';
        if (!parameters_.unresolved_path.empty()) {
            std::cout << "  zero list     = " << parameters_.unresolved_path << '\n';
        }
        if (zeros == 0) {
            std::cout << "RESULT check=4.16 modulus=" << Modulus
                      << " status=PASS unresolved=0\n";
        } else {
            std::cout << "RESULT check=4.16 modulus=" << Modulus
                      << " status=INCOMPLETE unresolved=" << zeros << '\n';
        }
        return checksum;
    }

private:
    RunParameters parameters_;
    int genus_bound_;
    int degree_;
    int n_count_;
    int block_width_;
    std::size_t packed_base_ = 0;
    std::size_t packed_transform_length_ = 0;

    std::vector<u32> inverses_;
    std::vector<int> mobius_;
    std::vector<std::vector<int>> divisors_;
    std::vector<u32> factorials_;
    std::vector<u32> inverse_factorials_;
    std::vector<u32> bernoulli_;

    // R[degree][m(m+1)/2+j] and H[h][degree].
    std::vector<std::vector<u32>> r_coefficients_;
    std::vector<std::vector<u32>> h_coefficients_;

    // Entry a stores the transform of sum_{b=0}^{degree-a} R_b.  Therefore
    // summing A_a times this entry gives sum_{a+b<=degree} A_a R_b = D_degree.
    std::vector<std::vector<u32>> cumulative_r_transforms_;

    void build_inverse_table() {
        inverses_[1] = 1;
        for (std::size_t value = 2; value < inverses_.size(); ++value) {
            if (value >= Modulus) {
                throw std::runtime_error("inverse table reaches the field characteristic");
            }
            const u32 quotient = static_cast<u32>(Modulus / value);
            const u32 remainder = static_cast<u32>(Modulus % value);
            inverses_[value] = F::subtract(0, F::multiply(quotient, inverses_[remainder]));
        }
    }

    void build_factorials() {
        factorials_.assign(genus_bound_ + 2, 1);
        inverse_factorials_.assign(genus_bound_ + 2, 1);
        for (int n = 1; n < static_cast<int>(factorials_.size()); ++n) {
            factorials_[n] = F::multiply(factorials_[n - 1], static_cast<u32>(n));
        }
        inverse_factorials_.back() = F::inverse(factorials_.back());
        for (int n = static_cast<int>(factorials_.size()) - 1; n >= 1; --n) {
            inverse_factorials_[n - 1] =
                F::multiply(inverse_factorials_[n], static_cast<u32>(n));
        }
    }

    u32 binomial(int n, int k) const {
        if (k < 0 || k > n) {
            return 0;
        }
        return F::multiply(factorials_[n],
                           F::multiply(inverse_factorials_[k],
                                       inverse_factorials_[n - k]));
    }

    void build_bernoulli_numbers() {
        bernoulli_.assign(genus_bound_ + 1, 0);
        bernoulli_[0] = 1;
        for (int m = 1; m <= genus_bound_; ++m) {
            u32 sum = 0;
            for (int k = 0; k < m; ++k) {
                sum = F::add(sum, F::multiply(binomial(m + 1, k), bernoulli_[k]));
            }
            bernoulli_[m] = F::subtract(
                0, F::multiply(sum, inverses_[m + 1]));
        }
    }

    std::size_t triangular_index(int m, int j) const {
        return static_cast<std::size_t>(m) * (m + 1) / 2 + j;
    }

    void print_run_header() const {
        const int block_count = (n_count_ + block_width_ - 1) / block_width_;
        const std::size_t r_transform_bytes =
            static_cast<std::size_t>(degree_ + 1) * packed_transform_length_ * sizeof(u32);
        const std::size_t maximum_local_width =
            static_cast<std::size_t>(block_width_ + genus_bound_);
        const std::size_t compact_a_bytes =
            static_cast<std::size_t>(degree_ + 1) * genus_bound_ *
            maximum_local_width * sizeof(u32);
        const std::size_t work_buffers = 2 * packed_transform_length_ * sizeof(u32);

        std::cout << "Direct large-prime check for Finite arithmetic verification 4.16\n"
                  << "  p                 = " << Modulus << '\n'
                  << "  primitive root    = " << PrimitiveRoot << '\n'
                  << "  genus range       = 1.." << genus_bound_ << '\n'
                  << "  n range           = 0.." << (n_count_ - 1) << '\n'
                  << "  requested D       = D_" << degree_ << '\n'
                  << "  block width       = " << block_width_ << " (" << block_count
                  << " blocks)\n"
                  << "  packed base       = " << packed_base_ << '\n'
                  << "  block NTT length  = " << packed_transform_length_ << '\n'
                  << "  persistent R hats = " << gibibytes(r_transform_bytes) << '\n'
                  << "  largest A block   = " << gibibytes(compact_a_bytes) << '\n'
                  << "  two work buffers  = " << gibibytes(work_buffers) << '\n'
                  << "The displayed memory figures omit R-construction temporaries and allocator overhead.\n\n"
                  << std::flush;
    }

    void build_r_series() {
        const WallClock phase_clock;
        std::cout << "[phase 1/3] Constructing R through w^" << degree_ << "\n"
                  << std::flush;

        const int packed_base = 2 * genus_bound_ + 1;
        const int maximum_encoded = genus_bound_ * packed_base + genus_bound_;
        const std::size_t transform_length =
            next_power_of_two(static_cast<std::size_t>(2) * maximum_encoded + 1);
        if ((Modulus - 1ULL) % transform_length != 0) {
            throw std::runtime_error("R construction needs an unsupported NTT length");
        }
        std::cout << "  R NTT length = " << transform_length << '\n' << std::flush;

        auto encode = [packed_base](int m, int j) {
            return static_cast<std::size_t>(m) * packed_base + j;
        };

        std::vector<std::vector<u32>> logarithm_kernel(
            degree_ + 1, std::vector<u32>(transform_length, 0));
        for (int m = 1; m <= genus_bound_; ++m) {
            logarithm_kernel[0][encode(m, 0)] = F::subtract(
                logarithm_kernel[0][encode(m, 0)], inverses_[m]);
            if (degree_ >= 1) {
                logarithm_kernel[1][encode(m, 0)] = F::add(
                    logarithm_kernel[1][encode(m, 0)], inverses_[m]);
            }
        }

        std::vector<u32> r_series(2 * genus_bound_ + 1, 0);
        std::vector<u32> log_r_series(genus_bound_ + 1, 0);
        std::vector<u32> f_power(2 * genus_bound_ + 1, 0);

        const WallClock kernel_clock;
        for (int ell = 2; ell <= 2 * genus_bound_; ++ell) {
            std::vector<u32> w_polynomial(degree_ + 1, 0);
            const u32 inverse_ell = inverses_[ell];
            for (int divisor : divisors_[ell]) {
                if (divisor <= degree_ && mobius_[ell / divisor] != 0) {
                    w_polynomial[divisor] = F::add(
                        w_polynomial[divisor],
                        F::multiply(F::normalize(-mobius_[ell / divisor]), inverse_ell));
                }
            }
            bool has_positive_degree = false;
            for (int degree = 1; degree <= degree_; ++degree) {
                has_positive_degree = has_positive_degree || w_polynomial[degree] != 0;
            }
            if (!has_positive_degree) {
                if ((ell <= 100 && ell % 10 == 0) || ell % 100 == 0 ||
                    ell == 2 * genus_bound_) {
                    std::cout << "  logarithm kernel ell=" << ell << '/'
                              << 2 * genus_bound_ << ", elapsed "
                              << duration_string(kernel_clock.seconds()) << '\n'
                              << std::flush;
                }
                continue;
            }

            const u32 constant = F::multiply(F::normalize(mobius_[ell]), inverse_ell);
            std::fill(r_series.begin(), r_series.end(), 0);
            r_series[0] = 1;
            std::vector<std::pair<int, u32>> support;
            for (int divisor : divisors_[ell]) {
                if (divisor < ell && mobius_[ell / divisor] != 0) {
                    const int exponent = ell - divisor;
                    if (exponent <= 2 * genus_bound_) {
                        r_series[exponent] = F::add(
                            r_series[exponent], F::normalize(mobius_[ell / divisor]));
                    }
                }
            }
            for (int exponent = 1; exponent <= 2 * genus_bound_; ++exponent) {
                if (r_series[exponent] != 0) {
                    support.emplace_back(exponent, r_series[exponent]);
                }
            }

            std::fill(log_r_series.begin(), log_r_series.end(), 0);
            for (int n = 1; n <= genus_bound_; ++n) {
                u32 sum = 0;
                for (const auto& [exponent, coefficient] : support) {
                    if (exponent >= n) {
                        break;
                    }
                    const int k = n - exponent;
                    sum = F::add(
                        sum,
                        F::multiply(F::multiply(static_cast<u32>(k), log_r_series[k]),
                                    coefficient));
                }
                const u32 logarithm_coefficient =
                    F::subtract(r_series[n], F::multiply(sum, inverses_[n]));
                log_r_series[n] = logarithm_coefficient;
                u32 adjusted = logarithm_coefficient;
                if (n % ell == 0) {
                    adjusted = F::subtract(adjusted, inverses_[n / ell]);
                }
                if (adjusted != 0) {
                    for (int degree = 1; degree <= degree_; ++degree) {
                        if (w_polynomial[degree] != 0) {
                            const std::size_t index = encode(n, 0);
                            logarithm_kernel[degree][index] = F::add(
                                logarithm_kernel[degree][index],
                                F::multiply(w_polynomial[degree], adjusted));
                        }
                    }
                }
            }

            std::vector<std::vector<u32>> w_powers(
                degree_ + 1, std::vector<u32>(degree_ + 1, 0));
            w_powers[0][0] = 1;
            for (int power = 1; power <= degree_; ++power) {
                for (int old_degree = 0; old_degree <= degree_; ++old_degree) {
                    if (w_powers[power - 1][old_degree] == 0) {
                        continue;
                    }
                    for (int added_degree = 1;
                         old_degree + added_degree <= degree_; ++added_degree) {
                        if (w_polynomial[added_degree] != 0) {
                            w_powers[power][old_degree + added_degree] = F::add(
                                w_powers[power][old_degree + added_degree],
                                F::multiply(w_powers[power - 1][old_degree],
                                            w_polynomial[added_degree]));
                        }
                    }
                }
            }

            std::vector<u32> constant_powers(genus_bound_ + degree_ + 2, 1);
            for (std::size_t power = 1; power < constant_powers.size(); ++power) {
                constant_powers[power] = F::multiply(constant_powers[power - 1], constant);
            }

            u32 ell_power = 1;
            const int maximum_s = 2 * genus_bound_ / ell;
            const int s_progress_stride = std::max(1, maximum_s / 8);
            if (ell <= 10) {
                std::cout << "    ell=" << ell << ": s=1.." << maximum_s
                          << " (the small ell values are the expensive part)\n"
                          << std::flush;
            }
            for (int s = 1; s <= maximum_s; ++s) {
                if (ell <= 10 &&
                    (s == 1 || s == maximum_s || s % s_progress_stride == 0)) {
                    std::cout << "      ell=" << ell << ", starting s=" << s << '/'
                              << maximum_s << ", kernel elapsed "
                              << duration_string(kernel_clock.seconds()) << '\n'
                              << std::flush;
                }
                ell_power = F::multiply(ell_power, static_cast<u32>(ell));
                const int maximum_f_degree = 2 * genus_bound_ - ell * s;
                std::fill(f_power.begin(), f_power.begin() + maximum_f_degree + 1, 0);
                f_power[0] = 1;
                for (int n = 1; n <= maximum_f_degree; ++n) {
                    u32 sum = 0;
                    for (const auto& [exponent, coefficient] : support) {
                        if (exponent > n) {
                            break;
                        }
                        const u32 weight = F::normalize(n + 1LL * (s - 1) * exponent);
                        sum = F::add(
                            sum,
                            F::multiply(F::multiply(weight, coefficient),
                                        f_power[n - exponent]));
                    }
                    f_power[n] = F::subtract(0, F::multiply(sum, inverses_[n]));
                }

                auto process_term = [&](int k, u32 coefficient) {
                    if (k < 1 || coefficient == 0) {
                        return;
                    }
                    const int minimum_power = std::max(1, k - genus_bound_);
                    const int maximum_power = std::min(k, degree_);
                    for (int power = minimum_power; power <= maximum_power; ++power) {
                        const int j = k - power;
                        u32 scalar = F::multiply(
                            coefficient,
                            F::multiply(binomial(k, power), constant_powers[j]));
                        if (scalar == 0) {
                            continue;
                        }
                        const int lower_n = std::max(0, 2 * j - ell * s);
                        const int upper_n = std::min(
                            maximum_f_degree, genus_bound_ + j - ell * s);
                        if (lower_n > upper_n) {
                            continue;
                        }
                        for (int output_degree = power;
                             output_degree <= degree_; ++output_degree) {
                            if (w_powers[power][output_degree] == 0) {
                                continue;
                            }
                            const u32 degree_scalar = F::multiply(
                                scalar, w_powers[power][output_degree]);
                            for (int n = lower_n; n <= upper_n; ++n) {
                                if (f_power[n] == 0) {
                                    continue;
                                }
                                const int m = ell * s + n - j;
                                const u32 value = F::multiply(
                                    degree_scalar,
                                    F::multiply(ell_power, f_power[n]));
                                const std::size_t index = encode(m, j);
                                logarithm_kernel[output_degree][index] = F::add(
                                    logarithm_kernel[output_degree][index], value);
                            }
                        }
                    }
                };

                process_term(s, inverses_[2 * s]);
                process_term(
                    s + 1,
                    F::subtract(0, F::multiply(inverses_[s], inverses_[s + 1])));
                for (int bernoulli_index = 2; bernoulli_index <= s;
                     bernoulli_index += 2) {
                    const int k = s - bernoulli_index + 1;
                    const u32 coefficient = F::subtract(
                        0,
                        F::multiply(
                            bernoulli_[bernoulli_index],
                            F::multiply(
                                binomial(s - 1, k),
                                F::multiply(inverses_[bernoulli_index],
                                            inverses_[bernoulli_index - 1]))));
                    process_term(k, coefficient);
                }
            }

            if ((ell <= 100 && ell % 10 == 0) || ell % 100 == 0 ||
                ell == 2 * genus_bound_) {
                std::cout << "  logarithm kernel ell=" << ell << '/'
                          << 2 * genus_bound_ << ", elapsed "
                          << duration_string(kernel_clock.seconds()) << '\n' << std::flush;
            }
        }

        std::cout << "  transforming logarithm-kernel degrees\n" << std::flush;
        for (int degree = 1; degree <= degree_; ++degree) {
            Transform::transform(logarithm_kernel[degree], false);
            std::cout << "    K degree " << degree << '/' << degree_ << '\n' << std::flush;
        }

        const std::size_t triangular_size =
            static_cast<std::size_t>(genus_bound_ + 1) * (genus_bound_ + 2) / 2;
        r_coefficients_.assign(degree_ + 1,
                               std::vector<u32>(triangular_size, 0));
        std::vector<std::vector<u32>> r_transforms(
            degree_ + 1, std::vector<u32>(transform_length, 0));
        r_coefficients_[0][triangular_index(0, 0)] = 1;
        r_coefficients_[0][triangular_index(1, 0)] = F::subtract(0, 1);
        r_transforms[0][0] = 1;
        r_transforms[0][encode(1, 0)] = F::subtract(0, 1);
        Transform::transform(r_transforms[0], false);

        std::vector<u32> temporary(transform_length, 0);
        for (int output_degree = 1; output_degree <= degree_; ++output_degree) {
            std::fill(temporary.begin(), temporary.end(), 0);
            for (int input_degree = 1; input_degree <= output_degree; ++input_degree) {
                const auto& kernel_transform = logarithm_kernel[input_degree];
                const auto& previous_transform =
                    r_transforms[output_degree - input_degree];
                for (std::size_t index = 0; index < transform_length; ++index) {
                    const u32 term = F::multiply(
                        static_cast<u32>(input_degree),
                        F::multiply(kernel_transform[index], previous_transform[index]));
                    temporary[index] = F::add(temporary[index], term);
                }
            }
            Transform::transform(temporary, true);
            for (int m = 0; m <= genus_bound_; ++m) {
                for (int j = 0; j <= m; ++j) {
                    const u32 value = F::multiply(
                        temporary[encode(m, j)], inverses_[output_degree]);
                    r_coefficients_[output_degree][triangular_index(m, j)] = value;
                    r_transforms[output_degree][encode(m, j)] = value;
                }
            }
            if (output_degree < degree_) {
                Transform::transform(r_transforms[output_degree], false);
            }
            std::cout << "  exponentiated R degree " << output_degree << '/'
                      << degree_ << ", phase elapsed "
                      << duration_string(phase_clock.seconds()) << '\n' << std::flush;
        }

        std::cout << "[phase 1/3] R complete in "
                  << duration_string(phase_clock.seconds()) << "\n\n" << std::flush;
    }

    void build_h_series() {
        const WallClock phase_clock;
        std::cout << "[phase 2/3] Constructing H through w^" << degree_ << " ... "
                  << std::flush;
        h_coefficients_.assign(genus_bound_, std::vector<u32>(degree_ + 1, 0));
        h_coefficients_[0][0] = 1;
        std::vector<std::vector<u32>> logarithm(
            genus_bound_, std::vector<u32>(degree_ + 1, 0));

        for (int q = 1; q < genus_bound_; ++q) {
            for (int a = 1; a <= degree_ && a <= q + 1; ++a) {
                const int bernoulli_index = q + 1 - a;
                const u32 bernoulli_value =
                    bernoulli_index == 1 ? inverses_[2] : bernoulli_[bernoulli_index];
                u32 value = F::multiply(
                    (a & 1) ? F::subtract(0, 1) : 1,
                    F::multiply(binomial(q + 1, a), bernoulli_value));
                value = F::subtract(
                    0, F::multiply(value, F::multiply(inverses_[q], inverses_[q + 1])));
                logarithm[q][a] = value;
            }
        }

        for (int q = 1; q < genus_bound_; ++q) {
            std::vector<u32> accumulated(degree_ + 1, 0);
            for (int previous_q = 1; previous_q <= q; ++previous_q) {
                for (int a = 1; a <= degree_; ++a) {
                    if (logarithm[previous_q][a] == 0) {
                        continue;
                    }
                    for (int b = 0; a + b <= degree_; ++b) {
                        if (h_coefficients_[q - previous_q][b] == 0) {
                            continue;
                        }
                        const u32 term = F::multiply(
                            static_cast<u32>(previous_q),
                            F::multiply(logarithm[previous_q][a],
                                        h_coefficients_[q - previous_q][b]));
                        accumulated[a + b] = F::add(accumulated[a + b], term);
                    }
                }
            }
            for (int degree = 0; degree <= degree_; ++degree) {
                h_coefficients_[q][degree] =
                    F::multiply(accumulated[degree], inverses_[q]);
            }
        }
        std::cout << "done in " << duration_string(phase_clock.seconds())
                  << "\n\n" << std::flush;
    }

    void build_cumulative_r_transforms() {
        const WallClock phase_clock;
        std::cout << "[phase 3/3] Precomputing the " << (degree_ + 1)
                  << " cumulative R transforms\n" << std::flush;
        cumulative_r_transforms_.assign(
            degree_ + 1, std::vector<u32>(packed_transform_length_, 0));

        for (int a_degree = 0; a_degree <= degree_; ++a_degree) {
            auto& packed = cumulative_r_transforms_[a_degree];
            const int maximum_r_degree = degree_ - a_degree;
            for (int m = 0; m < genus_bound_; ++m) {
                for (int j = 0; j <= m; ++j) {
                    u32 sum = 0;
                    const std::size_t triangular = triangular_index(m, j);
                    for (int r_degree = 0; r_degree <= maximum_r_degree; ++r_degree) {
                        sum = F::add(sum, r_coefficients_[r_degree][triangular]);
                    }
                    packed[static_cast<std::size_t>(m) * packed_base_ + j] = sum;
                }
            }
            Transform::transform(packed, false);
            std::cout << "  R cumulative cutoff " << maximum_r_degree
                      << " transformed (" << (a_degree + 1) << '/'
                      << (degree_ + 1) << "), elapsed "
                      << duration_string(phase_clock.seconds()) << '\n' << std::flush;
        }

        r_coefficients_.clear();
        r_coefficients_.shrink_to_fit();
        std::cout << "[phase 3/3] cumulative R transforms complete in "
                  << duration_string(phase_clock.seconds()) << "\n" << std::flush;
    }

    std::vector<std::vector<u32>> build_a_block(
        int input_begin,
        int output_end,
        std::vector<std::vector<u32>>& polynomial_checkpoint,
        int& checkpoint_n) const {
        const int local_width = output_end - input_begin;
        std::vector<std::vector<u32>> compact_a(
            degree_ + 1,
            std::vector<u32>(static_cast<std::size_t>(genus_bound_) * local_width, 0));

        const int next_output_begin = output_end;
        const int next_input_begin = std::max(0, next_output_begin - genus_bound_);
        std::vector<std::vector<u32>> next_checkpoint(
            genus_bound_, std::vector<u32>(degree_ + 1, 0));

        for (int h = 0; h < genus_bound_; ++h) {
            std::vector<u32> polynomial = polynomial_checkpoint[h];
            bool checkpoint_saved = (output_end == n_count_);
            std::vector<u32> product(degree_ + 1, 0);
            std::vector<u32> next_polynomial(degree_ + 1, 0);

            for (int n = input_begin; n < output_end; ++n) {
                if (n == next_input_begin && output_end < n_count_) {
                    next_checkpoint[h] = polynomial;
                    checkpoint_saved = true;
                }

                std::fill(product.begin(), product.end(), 0);
                for (int a = 0; a <= degree_; ++a) {
                    if (h_coefficients_[h][a] == 0) {
                        continue;
                    }
                    for (int b = 0; a + b <= degree_; ++b) {
                        if (polynomial[b] != 0) {
                            product[a + b] = F::add(
                                product[a + b],
                                F::multiply(h_coefficients_[h][a], polynomial[b]));
                        }
                    }
                }
                const std::size_t position =
                    static_cast<std::size_t>(h) * local_width + (n - input_begin);
                for (int degree = 0; degree <= degree_; ++degree) {
                    compact_a[degree][position] = product[degree];
                }

                std::fill(next_polynomial.begin(), next_polynomial.end(), 0);
                const u32 constant = F::normalize(h - 1LL + n);
                for (int degree = 0; degree <= degree_; ++degree) {
                    if (polynomial[degree] == 0) {
                        continue;
                    }
                    next_polynomial[degree] = F::add(
                        next_polynomial[degree], F::multiply(constant, polynomial[degree]));
                    if (degree + 1 <= degree_) {
                        next_polynomial[degree + 1] = F::add(
                            next_polynomial[degree + 1], polynomial[degree]);
                    }
                }
                const u32 inverse_n_plus_one = inverses_[n + 1];
                for (u32& coefficient : next_polynomial) {
                    coefficient = F::multiply(coefficient, inverse_n_plus_one);
                }
                polynomial.swap(next_polynomial);
            }

            if (!checkpoint_saved) {
                throw std::runtime_error("failed to save polynomial checkpoint");
            }
        }

        if (output_end < n_count_) {
            polynomial_checkpoint.swap(next_checkpoint);
            checkpoint_n = next_input_begin;
        }
        return compact_a;
    }

    bool belongs_to_checked_domain(int genus, int n) const {
        if (!parameters_.manuscript_domain) {
            return true;
        }
        if (3 * genus + 2 * n < 25) {
            return false;
        }
        return !(genus == 8 && n == 1) && !(genus == 12 && n == 0);
    }
};

void run_basic_self_tests() {
    validate_field_and_ntt<2013265921U, 31U>();
    validate_field_and_ntt<40961U, 3U>();
}

int run_smoke_test() {
    run_basic_self_tests();
    RunParameters parameters{
        20,                  // genus
        10,                  // compute the same D_10 as in production
        64,                  // n=0,...,63
        32,                  // two overlapping blocks
        false,               // scan the full small rectangle
        "smoke_unresolved.tsv"
    };
    DirectCoefficientChecker<40961U, 3U> checker(parameters);
    const u64 checksum = checker.run();
    std::cout << "SMOKE_CHECKSUM=" << checksum << '\n';
    if (checksum != 570606223ULL) {
        throw std::runtime_error(
            "smoke checksum disagrees with the archived table generator");
    }
    std::cout << "RESULT smoke-reference-comparison status=PASS\n";
    return 0;
}

int run_production_check() {
    run_basic_self_tests();
    RunParameters parameters{
        2000,
        10,
        50000,
        4096,
        true,
        "check_4_16_large_prime_unresolved.tsv"
    };
    DirectCoefficientChecker<2013265921U, 31U> checker(parameters);
    checker.run();
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            run_basic_self_tests();
            std::cout << "RESULT self-test status=PASS\n";
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "--smoke") {
            return run_smoke_test();
        }
        if (argc != 1) {
            std::cerr << "Usage: " << argv[0] << " [--self-test|--smoke]\n"
                      << "With no arguments, run the full 4.16 experiment.\n";
            return 2;
        }
        return run_production_check();
    } catch (const std::bad_alloc&) {
        std::cerr << "ERROR: memory allocation failed.  Try the program on a machine "
                     "with more available RAM.\n";
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 2;
    }
}
