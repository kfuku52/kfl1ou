#include <Rcpp.h>

#include <cmath>
#include <limits>
#include <vector>

extern "C" {
void effectiveSampleSize(int *Npo, int *npo, int *pNpo, int *rootpo,
                         double *transa, double *transb, int *des, int *anc,
                         int *edge, double *output);
void threepoint_l1ou(int *Npo, int *npo, int *pNpo, int *dYpo, int *dXpo,
                     int *rootpo, double *transa, double *transb, int *des,
                     int *anc, double *y, double *X, double *output);
}

// [[Rcpp::export]]
Rcpp::NumericVector effective_sample_size_c(int N, int n, int pN, int root,
                                            double transa,
                                            Rcpp::NumericVector transb,
                                            Rcpp::IntegerVector des,
                                            Rcpp::IntegerVector anc,
                                            Rcpp::IntegerVector edge) {
    if (N < 1 || n < 2 || pN < 1) {
        Rcpp::stop("invalid tree dimensions supplied to effective_sample_size_c");
    }
    const long long node_count = static_cast<long long>(n) + pN;
    if (node_count > std::numeric_limits<int>::max() ||
        N != node_count - 1 || root <= n || root > node_count) {
        Rcpp::stop(
            "tree dimensions or root are inconsistent in effective_sample_size_c"
        );
    }
    if (transb.size() != N || des.size() != N || anc.size() != N) {
        Rcpp::stop("tree edge arrays have inconsistent lengths");
    }
    if (!std::isfinite(transa) || transa < 0.0) {
        Rcpp::stop("tree edge lengths are invalid in effective_sample_size_c");
    }
    if (edge.size() < 1 || edge[edge.size() - 1] != N + 1) {
        Rcpp::stop("cut edges must end with the root-edge sentinel");
    }
    for (R_xlen_t i = 0; i < edge.size(); ++i) {
        if (edge[i] < 1 || edge[i] > N + 1 ||
            (i > 0 && edge[i] <= edge[i - 1])) {
            Rcpp::stop("cut edges must be strictly increasing valid indices");
        }
    }

    std::vector<int> parent(static_cast<std::size_t>(node_count + 1), 0);
    std::vector<int> incoming_edge(
        static_cast<std::size_t>(node_count + 1), -1
    );
    std::vector<int> last_child_edge(
        static_cast<std::size_t>(node_count + 1), -1
    );
    for (int i = 0; i < N; ++i) {
        if (!std::isfinite(transb[i]) || transb[i] < 0.0 ||
            anc[i] <= n || anc[i] > node_count || des[i] < 1 ||
            des[i] > node_count || des[i] == root || des[i] == anc[i] ||
            incoming_edge[des[i]] >= 0) {
            Rcpp::stop("tree topology or edge lengths are invalid in effective_sample_size_c");
        }
        parent[des[i]] = anc[i];
        incoming_edge[des[i]] = i;
        last_child_edge[anc[i]] = i;
    }
    for (int node = 1; node <= node_count; ++node) {
        if (node != root && incoming_edge[node] < 0) {
            Rcpp::stop("tree topology is disconnected in effective_sample_size_c");
        }
        if (node > n && last_child_edge[node] < 0) {
            Rcpp::stop(
                "internal tree node has no children in effective_sample_size_c"
            );
        }
        if (node > n && node != root &&
            last_child_edge[node] >= incoming_edge[node]) {
            Rcpp::stop(
                "tree edges must be supplied in postorder to effective_sample_size_c"
            );
        }
    }
    std::vector<unsigned char> state(
        static_cast<std::size_t>(node_count + 1), 0
    );
    state[root] = 2;
    for (int node = 1; node <= node_count; ++node) {
        int current = node;
        while (state[current] == 0) {
            state[current] = 1;
            current = parent[current];
        }
        if (state[current] == 1) {
            Rcpp::stop("tree topology contains a cycle in effective_sample_size_c");
        }
        current = node;
        while (state[current] == 1) {
            state[current] = 2;
            current = parent[current];
        }
    }

    Rcpp::NumericVector output(edge.size());
    effectiveSampleSize(&N, &n, &pN, &root, &transa, transb.begin(), des.begin(),
                        anc.begin(), edge.begin(), output.begin());
    return output;
}

// [[Rcpp::export]]
Rcpp::NumericVector threepoint_l1ou_c(int N, int n, int pN, int dY, int dX,
                                      int root, double transa,
                                      Rcpp::NumericVector transb,
                                      Rcpp::IntegerVector des,
                                      Rcpp::IntegerVector anc,
                                      Rcpp::NumericVector y,
                                      Rcpp::NumericVector X) {
    if (N < 1 || n < 2 || pN < 1 || dY < 1 || dX < 1) {
        Rcpp::stop("invalid dimensions supplied to threepoint_l1ou_c");
    }
    const long long node_count = static_cast<long long>(n) + pN;
    if (node_count > std::numeric_limits<int>::max() ||
        N != node_count - 1 || root <= n || root > node_count) {
        Rcpp::stop("tree dimensions or root are inconsistent in threepoint_l1ou_c");
    }
    const std::size_t nodes = static_cast<std::size_t>(node_count);
    const std::size_t max_size = std::numeric_limits<std::size_t>::max();
    const std::size_t max_double_elements = max_size / sizeof(double);
    auto checked_product = [](std::size_t left, std::size_t right) {
        if (right != 0 && left > max_size / right) {
            Rcpp::stop("requested working memory is too large in threepoint_l1ou_c");
        }
        return left * right;
    };
    const std::size_t dy = static_cast<std::size_t>(dY);
    const std::size_t dx = static_cast<std::size_t>(dX);
    const std::size_t y_workspace = checked_product(nodes, dy);
    const std::size_t x_workspace = checked_product(nodes, dx);
    const std::size_t yy_workspace = checked_product(y_workspace, dy);
    const std::size_t xx_workspace = checked_product(x_workspace, dx);
    const std::size_t xy_workspace = checked_product(x_workspace, dy);
    if (yy_workspace > max_double_elements ||
        xx_workspace > max_double_elements ||
        xy_workspace > max_double_elements) {
        Rcpp::stop("requested working memory is too large in threepoint_l1ou_c");
    }
    if (!std::isfinite(transa) || transa < 0.0 || transb.size() != N ||
        des.size() != N || anc.size() != N) {
        Rcpp::stop("tree edge arrays are invalid in threepoint_l1ou_c");
    }
    const long long expected_y = static_cast<long long>(n) * dY;
    const long long expected_x = static_cast<long long>(n) * dX;
    if (y.size() != expected_y || X.size() != expected_x) {
        Rcpp::stop("response or design dimensions are inconsistent in threepoint_l1ou_c");
    }

    std::vector<int> parent(static_cast<std::size_t>(node_count + 1), 0);
    std::vector<int> incoming_edge(static_cast<std::size_t>(node_count + 1), -1);
    std::vector<int> last_child_edge(static_cast<std::size_t>(node_count + 1), -1);
    for (int i = 0; i < N; ++i) {
        if (!std::isfinite(transb[i]) || transb[i] < 0.0 ||
            anc[i] <= n || anc[i] > node_count || des[i] < 1 ||
            des[i] > node_count || des[i] == root || des[i] == anc[i] ||
            incoming_edge[des[i]] >= 0) {
            Rcpp::stop("tree topology is invalid in threepoint_l1ou_c");
        }
        parent[des[i]] = anc[i];
        incoming_edge[des[i]] = i;
        last_child_edge[anc[i]] = i;
    }
    for (int node = 1; node <= node_count; ++node) {
        if (node != root && incoming_edge[node] < 0) {
            Rcpp::stop("tree topology is disconnected in threepoint_l1ou_c");
        }
        if (node > n && last_child_edge[node] < 0) {
            Rcpp::stop("internal tree node has no children in threepoint_l1ou_c");
        }
        if (node > n && node != root &&
            last_child_edge[node] >= incoming_edge[node]) {
            Rcpp::stop("tree edges must be supplied in postorder to threepoint_l1ou_c");
        }
    }
    std::vector<unsigned char> state(
        static_cast<std::size_t>(node_count + 1), 0
    );
    state[root] = 2;
    for (int node = 1; node <= node_count; ++node) {
        int current = node;
        while (state[current] == 0) {
            state[current] = 1;
            current = parent[current];
        }
        if (state[current] == 1) {
            Rcpp::stop("tree topology contains a cycle in threepoint_l1ou_c");
        }
        current = node;
        while (state[current] == 1) {
            state[current] = 2;
            current = parent[current];
        }
    }
    for (R_xlen_t i = 0; i < y.size(); ++i) {
        if (!std::isfinite(y[i])) {
            Rcpp::stop("response contains non-finite values in threepoint_l1ou_c");
        }
    }
    for (R_xlen_t i = 0; i < X.size(); ++i) {
        if (!std::isfinite(X[i])) {
            Rcpp::stop("design contains non-finite values in threepoint_l1ou_c");
        }
    }

    const long long max_output = std::numeric_limits<int>::max();
    const long long dy_square = static_cast<long long>(dY) * dY;
    const long long dx_square = static_cast<long long>(dX) * dX;
    if (dy_square > max_output || dx_square > max_output) {
        Rcpp::stop("requested output is too large in threepoint_l1ou_c");
    }
    const long long output_size = 2LL + dY + dy_square + dX + dx_square +
        static_cast<long long>(dX) * dY;
    if (output_size > max_output) {
        Rcpp::stop("requested output is too large in threepoint_l1ou_c");
    }
    Rcpp::NumericVector output(static_cast<R_xlen_t>(output_size));
    threepoint_l1ou(&N, &n, &pN, &dY, &dX, &root, &transa, transb.begin(),
                    des.begin(), anc.begin(), y.begin(), X.begin(),
                    output.begin());
    return output;
}
