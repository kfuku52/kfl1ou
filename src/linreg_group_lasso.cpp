#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

namespace {

struct GroupSpec {
    std::vector<int> cols;
    double hessian;
    double sqrt_size;
};

inline double dot_ptr(const double* x, const std::vector<double>& y, int n) {
    double out = 0.0;
    for (int i = 0; i < n; ++i) {
        out += x[i] * y[i];
    }
    return out;
}

inline double l2_norm(const std::vector<double>& x) {
    double out = 0.0;
    for (double value : x) {
        out += value * value;
    }
    return std::sqrt(out);
}

inline double l2_norm_n(const std::vector<double>& x, int size) {
    double out = 0.0;
    for (int i = 0; i < size; ++i) {
        out += x[static_cast<std::size_t>(i)] * x[static_cast<std::size_t>(i)];
    }
    return std::sqrt(out);
}

inline double group_norm(const std::vector<double>& coef, const GroupSpec& group) {
    double out = 0.0;
    for (int col : group.cols) {
        out += coef[col] * coef[col];
    }
    return std::sqrt(out);
}

inline double penalty_value(const std::vector<GroupSpec>& groups,
                            const std::vector<double>& norms,
                            double lambda) {
    double out = 0.0;
    for (std::size_t g = 0; g < groups.size(); ++g) {
        out += lambda * groups[g].sqrt_size * norms[g];
    }
    return out;
}

std::vector<GroupSpec> build_groups(const Rcpp::NumericMatrix& x,
                                    const Rcpp::IntegerVector& group) {
    const int n = x.nrow();
    const int p = x.ncol();

    if (group.size() != p) {
        Rcpp::stop("length(group) must equal ncol(x).");
    }

    std::unordered_map<int, std::size_t> positions;
    positions.reserve(static_cast<std::size_t>(p));

    std::vector<GroupSpec> groups;
    groups.reserve(static_cast<std::size_t>(p));
    const double* x_ptr = x.begin();

    for (int col = 0; col < p; ++col) {
        if (Rcpp::IntegerVector::is_na(group[col])) {
            Rcpp::stop("group indices must not contain NA values.");
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
            const double* x_col = x_ptr + static_cast<std::size_t>(n) * col;
            double sumsq = 0.0;
            for (int i = 0; i < n; ++i) {
                sumsq += x_col[i] * x_col[i];
            }
            max_diag = std::max(max_diag, 2.0 * sumsq);
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

    std::vector<GroupSpec> groups = build_groups(x, group);
    const double* x_ptr = x.begin();
    std::vector<double> residual(y.begin(), y.end());
    double max_value = 0.0;

    for (const GroupSpec& spec : groups) {
        double grad_sq = 0.0;
        for (int col : spec.cols) {
            const double* x_col = x_ptr + static_cast<std::size_t>(n) * col;
            const double grad = -2.0 * dot_ptr(x_col, residual, n);
            grad_sq += grad * grad;
        }
        max_value = std::max(max_value, std::sqrt(grad_sq) / spec.sqrt_size);
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
    const int n_lambda = lambda.size();

    if (y.size() != n) {
        Rcpp::stop("length(y) must equal nrow(x).");
    }

    std::vector<GroupSpec> groups = build_groups(x, group);
    const int n_groups = static_cast<int>(groups.size());
    const double sqrt_tol = std::sqrt(tol);
    const double min_scale = 1e-30;
    const double* x_ptr = x.begin();
    int max_group_size = 0;
    for (const GroupSpec& spec : groups) {
        max_group_size = std::max(max_group_size, static_cast<int>(spec.cols.size()));
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

    double rss = std::inner_product(residual.begin(), residual.end(), residual.begin(), 0.0);

    for (int pos = 0; pos < n_lambda; ++pos) {
        const double lambda_value = lambda[pos];
        double fn_value = rss + penalty_value(groups, norms, lambda_value);
        double d_fn = 1.0;
        double d_par = 1.0;
        bool do_all = false;
        int counter = 1;
        int iter_count = 0;

        while (d_fn > tol || d_par > sqrt_tol || !do_all) {
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
                    if (norms[g] != 0.0) {
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
                const GroupSpec& spec = groups[g_index];
                const int group_size = static_cast<int>(spec.cols.size());
                const double border = spec.sqrt_size * lambda_value;

                for (int k = 0; k < group_size; ++k) {
                    const int col = spec.cols[k];
                    const double* x_col = x_ptr + static_cast<std::size_t>(n) * col;
                    const double grad = -2.0 * dot_ptr(x_col, residual, n);
                    coef_group[static_cast<std::size_t>(k)] = coef[col];
                    gradient[static_cast<std::size_t>(k)] = grad;
                    cond[static_cast<std::size_t>(k)] = -grad + spec.hessian * coef[col];
                }

                const double cond_norm = l2_norm_n(cond, group_size);
                if (cond_norm > border) {
                    const double scale = (1.0 - border / cond_norm) / spec.hessian;
                    for (int k = 0; k < group_size; ++k) {
                        d[static_cast<std::size_t>(k)] =
                            scale * cond[static_cast<std::size_t>(k)] -
                            coef_group[static_cast<std::size_t>(k)];
                    }
                } else {
                    for (int k = 0; k < group_size; ++k) {
                        d[static_cast<std::size_t>(k)] = -coef_group[static_cast<std::size_t>(k)];
                    }
                }

                const bool has_step = cond_norm > border || norms[g_index] != 0.0;
                if (!has_step) {
                    norms[g_index] = group_norm(coef, spec);
                    continue;
                }

                double step_scale = 1.0;
                std::fill(xjd.begin(), xjd.end(), 0.0);

                for (int k = 0; k < group_size; ++k) {
                    if (d[static_cast<std::size_t>(k)] == 0.0) {
                        continue;
                    }
                    const int col = spec.cols[k];
                    const double* x_col = x_ptr + static_cast<std::size_t>(n) * col;
                    for (int i = 0; i < n; ++i) {
                        xjd[static_cast<std::size_t>(i)] +=
                            x_col[i] * d[static_cast<std::size_t>(k)];
                    }
                }

                double dot_res_xjd = 0.0;
                double dot_xjd_xjd = 0.0;
                for (int i = 0; i < n; ++i) {
                    dot_res_xjd += residual[i] * xjd[i];
                    dot_xjd_xjd += xjd[i] * xjd[i];
                }

                const double coef_norm_old = norms[g_index];
                double qh = 0.0;
                for (int k = 0; k < group_size; ++k) {
                    qh += gradient[static_cast<std::size_t>(k)] * d[static_cast<std::size_t>(k)];
                    full_step[static_cast<std::size_t>(k)] =
                        coef_group[static_cast<std::size_t>(k)] + d[static_cast<std::size_t>(k)];
                }
                qh += border * (l2_norm_n(full_step, group_size) - coef_norm_old);

                double rss_test = rss - 2.0 * step_scale * dot_res_xjd +
                    step_scale * step_scale * dot_xjd_xjd;
                for (int k = 0; k < group_size; ++k) {
                    coef_test[static_cast<std::size_t>(k)] =
                        coef_group[static_cast<std::size_t>(k)] +
                        step_scale * d[static_cast<std::size_t>(k)];
                }
                double left = rss_test + border * l2_norm_n(coef_test, group_size);
                double right = rss + border * coef_norm_old + sigma_ls * step_scale * qh;

                if (line_search) {
                    while (left > right && step_scale > min_scale) {
                        step_scale *= beta_ls;
                        rss_test = rss - 2.0 * step_scale * dot_res_xjd +
                            step_scale * step_scale * dot_xjd_xjd;
                        for (int k = 0; k < group_size; ++k) {
                            coef_test[static_cast<std::size_t>(k)] =
                                coef_group[static_cast<std::size_t>(k)] +
                                step_scale * d[static_cast<std::size_t>(k)];
                        }
                        left = rss_test + border * l2_norm_n(coef_test, group_size);
                        right = rss + border * coef_norm_old + sigma_ls * step_scale * qh;
                    }
                }

                if (step_scale <= min_scale) {
                    continue;
                }

                for (int k = 0; k < group_size; ++k) {
                    const int col = spec.cols[k];
                    const double new_coef = coef_test[static_cast<std::size_t>(k)];
                    const double scaled = std::fabs(new_coef - coef[col]) /
                        (1.0 + std::fabs(new_coef));
                    d_par = std::max(d_par, scaled);
                    coef[col] = new_coef;
                }
                for (int i = 0; i < n; ++i) {
                    residual[static_cast<std::size_t>(i)] -=
                        step_scale * xjd[static_cast<std::size_t>(i)];
                }
                rss = rss_test;
                norms[g_index] = l2_norm_n(coef_test, group_size);
            }

            fn_value = rss + penalty_value(groups, norms, lambda_value);
            d_fn = std::fabs(fn_old - fn_value) / (1.0 + std::fabs(fn_value));

            if (d_fn <= tol && d_par <= sqrt_tol) {
                counter = 0;
            }
        }

        for (int col = 0; col < p; ++col) {
            coefficients(col, pos) = coef[col];
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("coefficients") = coefficients,
        Rcpp::Named("converged") = converged
    );
}
