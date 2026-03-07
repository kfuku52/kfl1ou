#include <Rcpp.h>
#include <dlfcn.h>

namespace {

using get_threads_fn = int (*)();
using set_threads_fn = void (*)(int);

template <typename Fn>
Fn lookup_symbol(const char* name) {
    return reinterpret_cast<Fn>(dlsym(RTLD_DEFAULT, name));
}

int get_threads(const char* name) {
    get_threads_fn fn = lookup_symbol<get_threads_fn>(name);
    if (fn == nullptr) {
        return NA_INTEGER;
    }
    return fn();
}

bool set_threads(const char* name, int n_threads) {
    set_threads_fn fn = lookup_symbol<set_threads_fn>(name);
    if (fn == nullptr) {
        return false;
    }
    fn(n_threads);
    return true;
}

}  // namespace

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List get_runtime_thread_settings() {
    return Rcpp::List::create(
        Rcpp::Named("openblas_threads") = get_threads("openblas_get_num_threads"),
        Rcpp::Named("openmp_threads") = get_threads("omp_get_max_threads")
    );
}

// [[Rcpp::export]]
Rcpp::List set_runtime_thread_settings(int blas_threads = -1, int openmp_threads = -1) {
    bool openblas_set = false;
    bool openmp_set = false;

    if (blas_threads > 0) {
        openblas_set = set_threads("openblas_set_num_threads", blas_threads) ||
            set_threads("goto_set_num_threads", blas_threads);
    }
    if (openmp_threads > 0) {
        openmp_set = set_threads("omp_set_num_threads", openmp_threads);
    }

    return Rcpp::List::create(
        Rcpp::Named("openblas_set") = openblas_set,
        Rcpp::Named("openmp_set") = openmp_set
    );
}
