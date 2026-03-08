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

    Rcpp::NumericMatrix coefficients(p, n_lambda);
    Rcpp::LogicalVector converged(n_lambda, true);

    std::vector<double> coef(static_cast<std::size_t>(p), 0.0);
    std::vector<double> coef_old(static_cast<std::size_t>(p), 0.0);
    std::vector<double> norms(static_cast<std::size_t>(n_groups), 0.0);
    std::vector<double> residual(y.begin(), y.end());
    std::vector<double> xjd(static_cast<std::size_t>(n), 0.0);
    std::vector<double> coef_group;
    std::vector<double> d;
    std::vector<double> cond;
    std::vector<double> gradient;

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
            coef_old = coef;

            std::vector<int> active_groups;
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

            std::vector<double> start_pen(static_cast<std::size_t>(n_groups), 1.0);

            for (int g_index : active_groups) {
                const GroupSpec& spec = groups[g_index];
                const int group_size = static_cast<int>(spec.cols.size());
                const double border = spec.sqrt_size * lambda_value;

                coef_group.assign(static_cast<std::size_t>(group_size), 0.0);
                gradient.assign(static_cast<std::size_t>(group_size), 0.0);
                cond.assign(static_cast<std::size_t>(group_size), 0.0);
                d.assign(static_cast<std::size_t>(group_size), 0.0);

                for (int k = 0; k < group_size; ++k) {
                    const int col = spec.cols[k];
                    const double* x_col = x_ptr + static_cast<std::size_t>(n) * col;
                    const double grad = -2.0 * dot_ptr(x_col, residual, n);
                    coef_group[k] = coef[col];
                    gradient[k] = grad;
                    cond[k] = -grad + spec.hessian * coef[col];
                }

                const double cond_norm = l2_norm(cond);
                if (cond_norm > border) {
                    const double scale = (1.0 - border / cond_norm) / spec.hessian;
                    for (int k = 0; k < group_size; ++k) {
                        d[k] = scale * cond[k] - coef_group[k];
                    }
                } else {
                    for (int k = 0; k < group_size; ++k) {
                        d[k] = -coef_group[k];
                    }
                }

                bool has_step = false;
                for (double value : d) {
                    if (value != 0.0) {
                        has_step = true;
                        break;
                    }
                }
                if (!has_step) {
                    norms[g_index] = group_norm(coef, spec);
                    continue;
                }

                double step_scale = std::min(start_pen[g_index] / beta_ls, 1.0);
                std::fill(xjd.begin(), xjd.end(), 0.0);

                for (int k = 0; k < group_size; ++k) {
                    if (d[k] == 0.0) {
                        continue;
                    }
                    const int col = spec.cols[k];
                    const double* x_col = x_ptr + static_cast<std::size_t>(n) * col;
                    for (int i = 0; i < n; ++i) {
                        xjd[i] += x_col[i] * d[k];
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
                    qh += gradient[k] * d[k];
                }

                std::vector<double> full_step(static_cast<std::size_t>(group_size), 0.0);
                for (int k = 0; k < group_size; ++k) {
                    full_step[k] = coef_group[k] + d[k];
                }
                qh += border * (l2_norm(full_step) - coef_norm_old);

                double rss_test = rss - 2.0 * step_scale * dot_res_xjd +
                    step_scale * step_scale * dot_xjd_xjd;
                std::vector<double> coef_test(static_cast<std::size_t>(group_size), 0.0);
                for (int k = 0; k < group_size; ++k) {
                    coef_test[k] = coef_group[k] + step_scale * d[k];
                }
                double left = rss_test + border * l2_norm(coef_test);
                double right = rss + border * coef_norm_old + sigma_ls * step_scale * qh;

                if (line_search) {
                    while (left > right && step_scale > min_scale) {
                        step_scale *= beta_ls;
                        rss_test = rss - 2.0 * step_scale * dot_res_xjd +
                            step_scale * step_scale * dot_xjd_xjd;
                        for (int k = 0; k < group_size; ++k) {
                            coef_test[k] = coef_group[k] + step_scale * d[k];
                        }
                        left = rss_test + border * l2_norm(coef_test);
                        right = rss + border * coef_norm_old + sigma_ls * step_scale * qh;
                    }
                }

                if (step_scale <= min_scale) {
                    start_pen[g_index] = 1.0;
                    continue;
                }

                for (int k = 0; k < group_size; ++k) {
                    coef[spec.cols[k]] = coef_test[k];
                }
                for (int i = 0; i < n; ++i) {
                    residual[i] -= step_scale * xjd[i];
                }
                rss = rss_test;
                start_pen[g_index] = step_scale;
                norms[g_index] = group_norm(coef, spec);
            }

            fn_value = rss + penalty_value(groups, norms, lambda_value);
            d_par = 0.0;
            for (int col = 0; col < p; ++col) {
                const double scaled = std::fabs(coef[col] - coef_old[col]) / (1.0 + std::fabs(coef[col]));
                d_par = std::max(d_par, scaled);
            }
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
