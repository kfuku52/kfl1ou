#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

namespace {

struct GroupSpec {
    std::vector<int> cols;
    double hessian;
    double sqrt_size;
};

inline void check_user_interrupt_periodically(std::size_t& counter) {
    ++counter;
    if ((counter & 0x3fffU) == 0U) {
        Rcpp::checkUserInterrupt();
    }
}

double checked_double(long double value, const char* context) {
    const long double max_double =
        static_cast<long double>(std::numeric_limits<double>::max());
    if (!std::isfinite(value) || std::fabs(value) > max_double) {
        Rcpp::stop(
            std::string("numerical overflow in group-lasso ") + context +
            "; rescale the inputs"
        );
    }
    return static_cast<double>(value);
}

class ScaledSumSquares {
public:
    ScaledSumSquares() : scale_(0.0L), sum_squares_(1.0L) {}

    void add(double value) {
        if (!std::isfinite(value)) {
            Rcpp::stop("non-finite value encountered in group-lasso norm");
        }
        const long double magnitude =
            std::fabs(static_cast<long double>(value));
        if (magnitude == 0.0L) {
            return;
        }
        if (scale_ < magnitude) {
            const long double ratio = scale_ / magnitude;
            sum_squares_ = 1.0L + sum_squares_ * ratio * ratio;
            scale_ = magnitude;
        } else {
            const long double ratio = magnitude / scale_;
            sum_squares_ += ratio * ratio;
        }
    }

    double norm(const char* context) const {
        if (scale_ == 0.0L) {
            return 0.0;
        }
        return checked_double(
            scale_ * std::sqrt(sum_squares_), context
        );
    }

    double squared_norm(const char* context) const {
        if (scale_ == 0.0L) {
            return 0.0;
        }
        return checked_double(
            scale_ * scale_ * sum_squares_, context
        );
    }

private:
    long double scale_;
    long double sum_squares_;
};

double stable_dot_ptr(const double* x, const std::vector<double>& y, int n,
                      std::size_t& interrupt_counter) {
    long double sum = 0.0L;
    long double correction = 0.0L;
    for (int i = 0; i < n; ++i) {
        check_user_interrupt_periodically(interrupt_counter);
        const long double product =
            static_cast<long double>(x[i]) *
            static_cast<long double>(y[static_cast<std::size_t>(i)]);
        if (!std::isfinite(product)) {
            Rcpp::stop(
                "numerical overflow in group-lasso dot product; "
                "rescale the inputs"
            );
        }
        const long double updated = sum + product;
        if (!std::isfinite(updated)) {
            Rcpp::stop(
                "numerical overflow in group-lasso dot product; "
                "rescale the inputs"
            );
        }
        if (std::fabs(sum) >= std::fabs(product)) {
            correction += (sum - updated) + product;
        } else {
            correction += (product - updated) + sum;
        }
        sum = updated;
    }
    return checked_double(sum + correction, "dot product");
}

double stable_vector_dot(const std::vector<double>& left,
                         const std::vector<double>& right, int size,
                         std::size_t& interrupt_counter) {
    long double sum = 0.0L;
    long double correction = 0.0L;
    for (int i = 0; i < size; ++i) {
        check_user_interrupt_periodically(interrupt_counter);
        const std::size_t index = static_cast<std::size_t>(i);
        const long double product =
            static_cast<long double>(left[index]) *
            static_cast<long double>(right[index]);
        if (!std::isfinite(product)) {
            Rcpp::stop(
                "numerical overflow in group-lasso dot product; "
                "rescale the inputs"
            );
        }
        const long double updated = sum + product;
        if (!std::isfinite(updated)) {
            Rcpp::stop(
                "numerical overflow in group-lasso dot product; "
                "rescale the inputs"
            );
        }
        if (std::fabs(sum) >= std::fabs(product)) {
            correction += (sum - updated) + product;
        } else {
            correction += (product - updated) + sum;
        }
        sum = updated;
    }
    return checked_double(sum + correction, "dot product");
}

double l2_norm_n(const std::vector<double>& x, int size,
                 std::size_t& interrupt_counter) {
    ScaledSumSquares accumulator;
    for (int i = 0; i < size; ++i) {
        check_user_interrupt_periodically(interrupt_counter);
        accumulator.add(x[static_cast<std::size_t>(i)]);
    }
    return accumulator.norm("norm");
}

double squared_norm(const std::vector<double>& x,
                    std::size_t& interrupt_counter) {
    ScaledSumSquares accumulator;
    for (double value : x) {
        check_user_interrupt_periodically(interrupt_counter);
        accumulator.add(value);
    }
    return accumulator.squared_norm("residual sum of squares");
}

double group_norm(const std::vector<double>& coef, const GroupSpec& group,
                  std::size_t& interrupt_counter) {
    ScaledSumSquares accumulator;
    for (int col : group.cols) {
        check_user_interrupt_periodically(interrupt_counter);
        accumulator.add(coef[static_cast<std::size_t>(col)]);
    }
    return accumulator.norm("group norm");
}

double penalty_value(const std::vector<GroupSpec>& groups,
                     const std::vector<double>& norms, double lambda,
                     std::size_t& interrupt_counter) {
    long double out = 0.0L;
    for (std::size_t g = 0; g < groups.size(); ++g) {
        check_user_interrupt_periodically(interrupt_counter);
        out += static_cast<long double>(lambda) *
            static_cast<long double>(groups[g].sqrt_size) *
            static_cast<long double>(norms[g]);
        if (!std::isfinite(out)) {
            Rcpp::stop(
                "numerical overflow in group-lasso penalty; "
                "rescale the inputs"
            );
        }
    }
    return checked_double(out, "penalty");
}

double checked_multiply_add(double current, double left, double right,
                            const char* context) {
    return checked_double(
        static_cast<long double>(current) +
        static_cast<long double>(left) * static_cast<long double>(right),
        context
    );
}

double candidate_rss(double rss, double step_scale, double dot_res_xjd,
                     double dot_xjd_xjd) {
    long double value = static_cast<long double>(rss) -
        2.0L * static_cast<long double>(step_scale) *
            static_cast<long double>(dot_res_xjd) +
        static_cast<long double>(step_scale) *
            static_cast<long double>(step_scale) *
            static_cast<long double>(dot_xjd_xjd);
    if (!std::isfinite(value)) {
        Rcpp::stop(
            "numerical overflow in group-lasso residual sum of squares; "
            "rescale the inputs"
        );
    }
    if (value < 0.0L) {
        const long double scale = std::max(
            1.0L,
            static_cast<long double>(rss) +
                2.0L * std::fabs(
                    static_cast<long double>(step_scale) *
                    static_cast<long double>(dot_res_xjd)
                ) +
                static_cast<long double>(step_scale) *
                    static_cast<long double>(step_scale) *
                    static_cast<long double>(dot_xjd_xjd)
        );
        const long double rounding_tolerance =
            64.0L * static_cast<long double>(
                std::numeric_limits<double>::epsilon()
            ) * scale;
        if (value >= -rounding_tolerance) {
            value = 0.0L;
        } else {
            Rcpp::stop(
                "negative residual sum of squares encountered in "
                "group-lasso optimization"
            );
        }
    }
    return checked_double(value, "residual sum of squares");
}

std::vector<GroupSpec> build_groups(const Rcpp::NumericMatrix& x,
                                    const Rcpp::IntegerVector& group) {
    const int n = x.nrow();
    const int p = x.ncol();

    if (n < 1 || p < 1) {
        Rcpp::stop("x must have at least one row and one column.");
    }
    if (group.size() != p) {
        Rcpp::stop("length(group) must equal ncol(x).");
    }
    std::size_t interrupt_counter = 0;
    for (R_xlen_t i = 0; i < x.size(); ++i) {
        check_user_interrupt_periodically(interrupt_counter);
        if (!std::isfinite(x[i])) {
            Rcpp::stop("x must contain only finite values.");
        }
    }

    std::unordered_map<int, std::size_t> positions;
    positions.reserve(static_cast<std::size_t>(p));

    std::vector<GroupSpec> groups;
    groups.reserve(static_cast<std::size_t>(p));
    const double* x_ptr = x.begin();

    for (int col = 0; col < p; ++col) {
        check_user_interrupt_periodically(interrupt_counter);
        if (Rcpp::IntegerVector::is_na(group[col]) || group[col] < 1) {
            Rcpp::stop("group indices must be positive non-missing integers.");
        }

        const int key = group[col];
        auto it = positions.find(key);
        if (it == positions.end()) {
            positions.emplace(key, groups.size());
            groups.push_back(GroupSpec());
            groups.back().hessian = 1e-2;
            groups.back().sqrt_size = 0.0;
            it = positions.find(key);
        }
        groups[it->second].cols.push_back(col);
    }

    for (GroupSpec& spec : groups) {
        double max_diag = 1e-2;
        for (int col : spec.cols) {
            check_user_interrupt_periodically(interrupt_counter);
            const double* x_col = x_ptr + static_cast<std::size_t>(n) * col;
            ScaledSumSquares column_norm;
            for (int i = 0; i < n; ++i) {
                check_user_interrupt_periodically(interrupt_counter);
                column_norm.add(x_col[i]);
            }
            const double sumsq =
                column_norm.squared_norm("design-column sum of squares");
            const double diagonal = checked_double(
                2.0L * static_cast<long double>(sumsq),
                "group Hessian"
            );
            max_diag = std::max(max_diag, diagonal);
        }
        spec.hessian = max_diag;
        spec.sqrt_size = std::sqrt(static_cast<double>(spec.cols.size()));
    }

    return groups;
}

}  // namespace

// [[Rcpp::export]]
double linreg_group_lasso_lambda_max_cpp(Rcpp::NumericMatrix x,
                                         Rcpp::NumericVector y,
                                         Rcpp::IntegerVector group) {
    const int n = x.nrow();
    if (y.size() != n) {
        Rcpp::stop("length(y) must equal nrow(x).");
    }
    std::size_t interrupt_counter = 0;
    for (R_xlen_t i = 0; i < y.size(); ++i) {
        check_user_interrupt_periodically(interrupt_counter);
        if (!std::isfinite(y[i])) {
            Rcpp::stop("y must contain only finite values.");
        }
    }

    std::vector<GroupSpec> groups = build_groups(x, group);
    const double* x_ptr = x.begin();
    std::vector<double> residual(y.begin(), y.end());
    double max_value = 0.0;

    for (const GroupSpec& spec : groups) {
        check_user_interrupt_periodically(interrupt_counter);
        ScaledSumSquares gradient_norm;
        for (int col : spec.cols) {
            const double* x_col = x_ptr + static_cast<std::size_t>(n) * col;
            const double grad = checked_double(
                -2.0L * static_cast<long double>(
                    stable_dot_ptr(
                        x_col, residual, n, interrupt_counter
                    )
                ),
                "gradient"
            );
            gradient_norm.add(grad);
        }
        const double group_value = checked_double(
            static_cast<long double>(
                gradient_norm.norm("gradient norm")
            ) / static_cast<long double>(spec.sqrt_size),
            "maximum penalty"
        );
        max_value = std::max(max_value, group_value);
    }

    return max_value;
}

// [[Rcpp::export]]
Rcpp::List linreg_group_lasso_path_cpp(Rcpp::NumericMatrix x,
                                       Rcpp::NumericVector y,
                                       Rcpp::IntegerVector group,
                                       Rcpp::NumericVector lambda,
                                       double tol = 5e-8,
                                       int max_iter = 500,
                                       int inner_loops = 10,
                                       double beta_ls = 0.5,
                                       double sigma_ls = 0.1,
                                       bool line_search = true) {
    const int n = x.nrow();
    const int p = x.ncol();
    const R_xlen_t lambda_count = lambda.size();
    if (lambda_count > static_cast<R_xlen_t>(
            std::numeric_limits<int>::max())) {
        Rcpp::stop("lambda is too long for an R matrix solution path.");
    }
    const int n_lambda = static_cast<int>(lambda_count);

    if (y.size() != n) {
        Rcpp::stop("length(y) must equal nrow(x).");
    }
    if (n_lambda < 1) {
        Rcpp::stop("lambda must contain at least one value.");
    }
    if (!std::isfinite(tol) || tol <= 0.0 || max_iter < 1 ||
        inner_loops < 1 || !std::isfinite(beta_ls) || beta_ls <= 0.0 ||
        beta_ls >= 1.0 || !std::isfinite(sigma_ls) || sigma_ls <= 0.0 ||
        sigma_ls >= 1.0) {
        Rcpp::stop("invalid group-lasso control parameters.");
    }
    std::size_t interrupt_counter = 0;
    for (R_xlen_t i = 0; i < y.size(); ++i) {
        check_user_interrupt_periodically(interrupt_counter);
        if (!std::isfinite(y[i])) {
            Rcpp::stop("y must contain only finite values.");
        }
    }
    for (R_xlen_t i = 0; i < lambda.size(); ++i) {
        check_user_interrupt_periodically(interrupt_counter);
        if (!std::isfinite(lambda[i]) || lambda[i] < 0.0) {
            Rcpp::stop("lambda must contain finite non-negative values.");
        }
    }

    std::vector<GroupSpec> groups = build_groups(x, group);
    const int n_groups = static_cast<int>(groups.size());
    const double sqrt_tol = std::sqrt(tol);
    const double min_scale = 1e-30;
    const double* x_ptr = x.begin();
    int max_group_size = 0;
    for (const GroupSpec& spec : groups) {
        check_user_interrupt_periodically(interrupt_counter);
        max_group_size = std::max(max_group_size, static_cast<int>(spec.cols.size()));
    }

    if (p > 0 && static_cast<std::size_t>(p) >
            static_cast<std::size_t>(R_XLEN_T_MAX) /
                static_cast<std::size_t>(n_lambda)) {
        Rcpp::stop("requested group-lasso solution path is too large.");
    }
    Rcpp::NumericMatrix coefficients(p, n_lambda);
    Rcpp::LogicalVector converged(n_lambda, true);

    std::vector<double> coef(static_cast<std::size_t>(p), 0.0);
    std::vector<double> norms(static_cast<std::size_t>(n_groups), 0.0);
    std::vector<double> residual(y.begin(), y.end());
    std::vector<double> xjd(static_cast<std::size_t>(n), 0.0);
    std::vector<double> coef_group(static_cast<std::size_t>(max_group_size), 0.0);
    std::vector<double> d(static_cast<std::size_t>(max_group_size), 0.0);
    std::vector<double> cond(static_cast<std::size_t>(max_group_size), 0.0);
    std::vector<double> gradient(static_cast<std::size_t>(max_group_size), 0.0);
    std::vector<double> coef_test(static_cast<std::size_t>(max_group_size), 0.0);
    std::vector<double> full_step(static_cast<std::size_t>(max_group_size), 0.0);
    std::vector<int> active_groups;
    active_groups.reserve(static_cast<std::size_t>(n_groups));

    double rss = squared_norm(residual, interrupt_counter);

    for (int pos = 0; pos < n_lambda; ++pos) {
        check_user_interrupt_periodically(interrupt_counter);
        const double lambda_value = lambda[pos];
        double fn_value = checked_double(
            static_cast<long double>(rss) +
                static_cast<long double>(
                    penalty_value(
                        groups, norms, lambda_value, interrupt_counter
                    )
                ),
            "objective"
        );
        double d_fn = 1.0;
        double d_par = 1.0;
        bool do_all = false;
        bool stalled = false;
        int counter = 1;
        int iter_count = 0;

        while (d_fn > tol || d_par > sqrt_tol || !do_all) {
            check_user_interrupt_periodically(interrupt_counter);
            if (iter_count >= max_iter) {
                converged[pos] = false;
                break;
            }

            const double fn_old = fn_value;
            d_par = 0.0;
            active_groups.clear();
            if (counter == 0 || counter > inner_loops) {
                do_all = true;
                active_groups.resize(static_cast<std::size_t>(n_groups));
                std::iota(active_groups.begin(), active_groups.end(), 0);
                counter = 1;
            } else {
                do_all = false;
                for (int g = 0; g < n_groups; ++g) {
                    check_user_interrupt_periodically(interrupt_counter);
                    if (norms[static_cast<std::size_t>(g)] != 0.0) {
                        active_groups.push_back(g);
                    }
                }
                if (active_groups.empty()) {
                    do_all = true;
                    active_groups.resize(static_cast<std::size_t>(n_groups));
                    std::iota(active_groups.begin(), active_groups.end(), 0);
                } else {
                    counter += 1;
                }
            }

            if (do_all) {
                iter_count += 1;
            }

            for (int g_index : active_groups) {
                check_user_interrupt_periodically(interrupt_counter);
                const GroupSpec& spec =
                    groups[static_cast<std::size_t>(g_index)];
                const int group_size = static_cast<int>(spec.cols.size());
                const double border = checked_double(
                    static_cast<long double>(spec.sqrt_size) *
                        static_cast<long double>(lambda_value),
                    "penalty threshold"
                );

                for (int k = 0; k < group_size; ++k) {
                    check_user_interrupt_periodically(interrupt_counter);
                    const int col =
                        spec.cols[static_cast<std::size_t>(k)];
                    const double* x_col = x_ptr + static_cast<std::size_t>(n) * col;
                    const double grad = checked_double(
                        -2.0L * static_cast<long double>(
                            stable_dot_ptr(
                                x_col, residual, n, interrupt_counter
                            )
                        ),
                        "gradient"
                    );
                    coef_group[static_cast<std::size_t>(k)] =
                        coef[static_cast<std::size_t>(col)];
                    gradient[static_cast<std::size_t>(k)] = grad;
                    cond[static_cast<std::size_t>(k)] = checked_double(
                        -static_cast<long double>(grad) +
                            static_cast<long double>(spec.hessian) *
                                static_cast<long double>(
                                    coef[static_cast<std::size_t>(col)]
                                ),
                        "proximal condition"
                    );
                }

                const double cond_norm =
                    l2_norm_n(cond, group_size, interrupt_counter);
                if (cond_norm > border) {
                    const double scale = checked_double(
                        (1.0L -
                            static_cast<long double>(border) /
                                static_cast<long double>(cond_norm)) /
                            static_cast<long double>(spec.hessian),
                        "proximal scale"
                    );
                    for (int k = 0; k < group_size; ++k) {
                        check_user_interrupt_periodically(interrupt_counter);
                        d[static_cast<std::size_t>(k)] =
                            checked_double(
                                static_cast<long double>(scale) *
                                    static_cast<long double>(
                                        cond[static_cast<std::size_t>(k)]
                                    ) -
                                    static_cast<long double>(
                                        coef_group[
                                            static_cast<std::size_t>(k)
                                        ]
                                    ),
                                "coefficient step"
                            );
                    }
                } else {
                    for (int k = 0; k < group_size; ++k) {
                        check_user_interrupt_periodically(interrupt_counter);
                        d[static_cast<std::size_t>(k)] = -coef_group[static_cast<std::size_t>(k)];
                    }
                }

                const bool has_step = cond_norm > border ||
                    norms[static_cast<std::size_t>(g_index)] != 0.0;
                if (!has_step) {
                    norms[static_cast<std::size_t>(g_index)] =
                        group_norm(coef, spec, interrupt_counter);
                    continue;
                }

                bool numerically_zero_step = true;
                const double step_roundoff =
                    16.0 * std::numeric_limits<double>::epsilon();
                for (int k = 0; k < group_size; ++k) {
                    check_user_interrupt_periodically(interrupt_counter);
                    if (std::fabs(d[static_cast<std::size_t>(k)]) >
                            step_roundoff *
                                (1.0 +
                                    std::fabs(
                                        coef_group[
                                            static_cast<std::size_t>(k)
                                        ]
                                    ))) {
                        numerically_zero_step = false;
                        break;
                    }
                }
                if (numerically_zero_step) {
                    norms[static_cast<std::size_t>(g_index)] =
                        group_norm(coef, spec, interrupt_counter);
                    continue;
                }

                double step_scale = 1.0;
                std::fill(xjd.begin(), xjd.end(), 0.0);

                for (int k = 0; k < group_size; ++k) {
                    check_user_interrupt_periodically(interrupt_counter);
                    if (d[static_cast<std::size_t>(k)] == 0.0) {
                        continue;
                    }
                    const int col =
                        spec.cols[static_cast<std::size_t>(k)];
                    const double* x_col = x_ptr + static_cast<std::size_t>(n) * col;
                    for (int i = 0; i < n; ++i) {
                        check_user_interrupt_periodically(interrupt_counter);
                        const std::size_t index = static_cast<std::size_t>(i);
                        xjd[index] = checked_multiply_add(
                            xjd[index], x_col[i],
                            d[static_cast<std::size_t>(k)],
                            "linear predictor update"
                        );
                    }
                }

                const double dot_res_xjd = stable_vector_dot(
                    residual, xjd, n, interrupt_counter
                );
                ScaledSumSquares xjd_norm;
                for (int i = 0; i < n; ++i) {
                    check_user_interrupt_periodically(interrupt_counter);
                    xjd_norm.add(xjd[static_cast<std::size_t>(i)]);
                }
                const double dot_xjd_xjd =
                    xjd_norm.squared_norm("step sum of squares");

                const double coef_norm_old =
                    norms[static_cast<std::size_t>(g_index)];
                for (int k = 0; k < group_size; ++k) {
                    check_user_interrupt_periodically(interrupt_counter);
                    full_step[static_cast<std::size_t>(k)] =
                        checked_double(
                            static_cast<long double>(
                                coef_group[static_cast<std::size_t>(k)]
                            ) +
                                static_cast<long double>(
                                    d[static_cast<std::size_t>(k)]
                                ),
                            "full coefficient step"
                        );
                }
                const double gradient_step = stable_vector_dot(
                    gradient, d, group_size, interrupt_counter
                );
                const double full_step_norm =
                    l2_norm_n(full_step, group_size, interrupt_counter);
                const double qh = checked_double(
                    static_cast<long double>(gradient_step) +
                        static_cast<long double>(border) *
                            (static_cast<long double>(full_step_norm) -
                                static_cast<long double>(coef_norm_old)),
                    "line-search directional derivative"
                );

                double rss_test = candidate_rss(
                    rss, step_scale, dot_res_xjd, dot_xjd_xjd
                );
                for (int k = 0; k < group_size; ++k) {
                    check_user_interrupt_periodically(interrupt_counter);
                    coef_test[static_cast<std::size_t>(k)] =
                        checked_double(
                            static_cast<long double>(
                                coef_group[static_cast<std::size_t>(k)]
                            ) +
                                static_cast<long double>(step_scale) *
                                    static_cast<long double>(
                                        d[static_cast<std::size_t>(k)]
                                    ),
                            "candidate coefficient"
                        );
                }
                double left = checked_double(
                    static_cast<long double>(rss_test) +
                        static_cast<long double>(border) *
                            static_cast<long double>(
                                l2_norm_n(
                                    coef_test, group_size,
                                    interrupt_counter
                                )
                            ),
                    "line-search objective"
                );
                double right = checked_double(
                    static_cast<long double>(rss) +
                        static_cast<long double>(border) *
                            static_cast<long double>(coef_norm_old) +
                        static_cast<long double>(sigma_ls) *
                            static_cast<long double>(step_scale) *
                            static_cast<long double>(qh),
                    "line-search bound"
                );

                if (line_search) {
                    while (left > right && step_scale > min_scale) {
                        check_user_interrupt_periodically(interrupt_counter);
                        step_scale = checked_double(
                            static_cast<long double>(step_scale) *
                                static_cast<long double>(beta_ls),
                            "line-search scale"
                        );
                        rss_test = candidate_rss(
                            rss, step_scale, dot_res_xjd, dot_xjd_xjd
                        );
                        for (int k = 0; k < group_size; ++k) {
                            check_user_interrupt_periodically(
                                interrupt_counter
                            );
                            coef_test[static_cast<std::size_t>(k)] =
                                checked_double(
                                    static_cast<long double>(
                                        coef_group[
                                            static_cast<std::size_t>(k)
                                        ]
                                    ) +
                                        static_cast<long double>(
                                            step_scale
                                        ) *
                                            static_cast<long double>(
                                                d[
                                                    static_cast<std::size_t>(
                                                        k
                                                    )
                                                ]
                                            ),
                                    "candidate coefficient"
                                );
                        }
                        left = checked_double(
                            static_cast<long double>(rss_test) +
                                static_cast<long double>(border) *
                                    static_cast<long double>(
                                        l2_norm_n(
                                            coef_test, group_size,
                                            interrupt_counter
                                        )
                                    ),
                            "line-search objective"
                        );
                        right = checked_double(
                            static_cast<long double>(rss) +
                                static_cast<long double>(border) *
                                    static_cast<long double>(coef_norm_old) +
                                static_cast<long double>(sigma_ls) *
                                    static_cast<long double>(step_scale) *
                                    static_cast<long double>(qh),
                            "line-search bound"
                        );
                    }
                }

                if (step_scale <= min_scale) {
                    converged[pos] = false;
                    stalled = true;
                    break;
                }

                for (int k = 0; k < group_size; ++k) {
                    check_user_interrupt_periodically(interrupt_counter);
                    const int col =
                        spec.cols[static_cast<std::size_t>(k)];
                    const double new_coef = coef_test[static_cast<std::size_t>(k)];
                    const double scaled =
                        std::fabs(
                            new_coef -
                            coef[static_cast<std::size_t>(col)]
                        ) /
                        (1.0 + std::fabs(new_coef));
                    if (!std::isfinite(scaled)) {
                        Rcpp::stop(
                            "non-finite parameter change in group-lasso "
                            "optimization"
                        );
                    }
                    d_par = std::max(d_par, scaled);
                    coef[static_cast<std::size_t>(col)] = new_coef;
                }
                ScaledSumSquares updated_residual_norm;
                for (int i = 0; i < n; ++i) {
                    check_user_interrupt_periodically(interrupt_counter);
                    const std::size_t index = static_cast<std::size_t>(i);
                    residual[index] = checked_double(
                        static_cast<long double>(residual[index]) -
                            static_cast<long double>(step_scale) *
                                static_cast<long double>(xjd[index]),
                        "residual update"
                    );
                    updated_residual_norm.add(residual[index]);
                }
                rss = updated_residual_norm.squared_norm(
                    "residual sum of squares"
                );
                norms[static_cast<std::size_t>(g_index)] =
                    l2_norm_n(
                        coef_test, group_size, interrupt_counter
                    );
            }

            if (stalled) {
                break;
            }
            fn_value = checked_double(
                static_cast<long double>(rss) +
                    static_cast<long double>(
                        penalty_value(
                            groups, norms, lambda_value,
                            interrupt_counter
                        )
                    ),
                "objective"
            );
            d_fn = checked_double(
                std::fabs(
                    static_cast<long double>(fn_old) -
                    static_cast<long double>(fn_value)
                ) /
                    (1.0L +
                        std::fabs(static_cast<long double>(fn_value))),
                "objective change"
            );

            if (d_fn <= tol && d_par <= sqrt_tol) {
                counter = 0;
            }
        }

        for (int col = 0; col < p; ++col) {
            check_user_interrupt_periodically(interrupt_counter);
            coefficients(col, pos) =
                coef[static_cast<std::size_t>(col)];
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("coefficients") = coefficients,
        Rcpp::Named("converged") = converged
    );
}
