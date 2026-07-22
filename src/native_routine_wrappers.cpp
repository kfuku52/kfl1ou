#include <Rcpp.h>

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
    if (N < 1 || n < 2 || pN < 1 || root < 1 || root > n + pN) {
        Rcpp::stop("invalid tree dimensions supplied to effective_sample_size_c");
    }
    if (transb.size() != N || des.size() != N || anc.size() != N) {
        Rcpp::stop("tree edge arrays have inconsistent lengths");
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
    Rcpp::NumericVector output(2 + dY + dY * dY + dX + dX * dX + dX * dY);
    threepoint_l1ou(&N, &n, &pN, &dY, &dX, &root, &transa, transb.begin(),
                    des.begin(), anc.begin(), y.begin(), X.begin(),
                    output.begin());
    return output;
}
