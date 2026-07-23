#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <Rcpp.h>

using namespace Rcpp;
typedef Matrix<REALSXP> NumericMatrix;

#define RASSERT(condition){if(!(condition)){throw std::range_error(std::string("internal error!@")+__FILE__);}}

// [[Rcpp::plugins(cpp11)]]
void one_step(const int i1, const int i2, const int e1, const int e2,
        const int counter,
        Rcpp::NumericMatrix &edgeList, //the third column contains lengths
        Rcpp::NumericVector &tips, 
        Rcpp::NumericMatrix &F, Rcpp::NumericMatrix &G, 
        Rcpp::NumericMatrix &D, Rcpp::NumericMatrix &B,
        double &rootEdge){

    const int nEdges = edgeList.nrow();

    //int e1=-1, e2=-1;
    //for (int i=0; i<nEdges; ++i)
    //    if (edgeList(i,1) == i1)
    //        e1 = i;
    //    else if (edgeList(i,1) == i2)
    //        e2 = i;

    RASSERT( e1!=-1 && e2!=-1 );

    const int i3  = edgeList(e1,0);
    RASSERT( edgeList(e2,0) == i3 );
    
    const double t1 = edgeList(e1,2);
    const double t2 = edgeList(e2,2);

    int e3 = -1; // -1 means root index
    for(int i=0; i<nEdges; ++i) 
        if( edgeList(i,1) == i3){
            e3 = i;
            break;
        }

    double t3;
    if (e3 == -1)
        t3 = rootEdge>0?rootEdge:0;
    else
        t3 = edgeList(e3,2);
    
    const double  u = t1+t2;
    if (!std::isfinite(u) || u <= 0.0) {
        Rcpp::stop("encountered sibling branches with non-positive total length");
    }
    const double us = std::sqrt(u);

    D(_,counter) = (F(_,i1) - F(_,i2))/us;
    B(_,counter) = (G(_,i1)*t1 - G(_,i2)*t2)/us;

    F(_,i3) = (F(_,i1)*t2 + F(_,i2)*t1)/u;
    G(_,i3) = G(_,i1) + G(_,i2);

    if( e3>=0 )
        edgeList(e3,2) = t3 + 1/(1/t1+1/t2);
    else
        rootEdge = t3 + 1/(1/t1+1/t2);

    // erase-remove idiom
    tips.erase( std::remove( tips.begin(), tips.end(), i1 ), tips.end() );
    tips.erase( std::remove( tips.begin(), tips.end(), i2 ), tips.end() );
    tips.push_back(i3);
}



// [[Rcpp::export]]
Rcpp::List cmp_sqrt_OU_covariance(Rcpp::NumericMatrix edgeList, int nTips, double rootEdge){
    if (edgeList.ncol() != 3 || nTips < 2) {
        Rcpp::stop("edgeList must describe a rooted, strictly bifurcating tree");
    }
    const long long node_count_ll = 2LL * nTips - 1LL;
    const long long edge_count_ll = 2LL * (nTips - 1LL);
    if (node_count_ll > std::numeric_limits<int>::max() ||
        edge_count_ll > std::numeric_limits<int>::max() ||
        edgeList.nrow() != edge_count_ll) {
        Rcpp::stop("tree dimensions are too large or inconsistent");
    }
    const int node_count = static_cast<int>(node_count_ll);
    const int edge_count = static_cast<int>(edge_count_ll);
    const long long workspace_elements =
        static_cast<long long>(nTips) * node_count;
    const long long square_elements =
        static_cast<long long>(nTips) * nTips;
    const long long max_r_length =
        static_cast<long long>(std::numeric_limits<R_xlen_t>::max());
    if (workspace_elements > max_r_length || square_elements > max_r_length) {
        Rcpp::stop("requested covariance workspace is too large");
    }
    if (!std::isfinite(rootEdge) || rootEdge < 0.0) {
        Rcpp::stop("root edge must be finite and non-negative");
    }

    std::vector<int> parent(static_cast<std::size_t>(node_count), -1);
    std::vector<int> incoming_edge(static_cast<std::size_t>(node_count), -1);
    std::vector<int> child_count(static_cast<std::size_t>(node_count), 0);
    std::vector<int> last_child_edge(static_cast<std::size_t>(node_count), -1);
    for (int edge = 0; edge < edge_count; ++edge) {
        const double parent_value = edgeList(edge, 0);
        const double child_value = edgeList(edge, 1);
        const double length = edgeList(edge, 2);
        if (!std::isfinite(parent_value) || !std::isfinite(child_value) ||
            !std::isfinite(length) || length < 0.0 ||
            parent_value != std::floor(parent_value) ||
            child_value != std::floor(child_value) ||
            parent_value < nTips || parent_value >= node_count ||
            child_value < 0.0 || child_value >= node_count) {
            Rcpp::stop("edgeList contains invalid node identifiers or lengths");
        }
        const int ancestor = static_cast<int>(parent_value);
        const int descendant = static_cast<int>(child_value);
        if (descendant == ancestor || incoming_edge[descendant] >= 0) {
            Rcpp::stop("edgeList contains invalid tree topology");
        }
        parent[descendant] = ancestor;
        incoming_edge[descendant] = edge;
        child_count[ancestor] += 1;
        last_child_edge[ancestor] = edge;
    }

    int root = -1;
    for (int node = 0; node < node_count; ++node) {
        if (incoming_edge[node] < 0) {
            if (node < nTips || root >= 0) {
                Rcpp::stop("edgeList must contain exactly one internal root");
            }
            root = node;
        }
        if (node >= nTips && child_count[node] != 2) {
            Rcpp::stop("every internal node must have exactly two children");
        }
        if (node >= nTips && incoming_edge[node] >= 0 &&
            last_child_edge[node] >= incoming_edge[node]) {
            Rcpp::stop("edgeList must be supplied in postorder");
        }
    }
    if (root < 0) {
        Rcpp::stop("edgeList must contain exactly one internal root");
    }
    for (int edge = 0; edge < edge_count; edge += 2) {
        if (edgeList(edge, 0) != edgeList(edge + 1, 0)) {
            Rcpp::stop("sibling edges must be adjacent in postorder");
        }
    }

    std::vector<unsigned char> state(
        static_cast<std::size_t>(node_count), 0
    );
    state[root] = 2;
    for (int node = 0; node < node_count; ++node) {
        int current = node;
        while (state[current] == 0) {
            state[current] = 1;
            current = parent[current];
        }
        if (state[current] == 1) {
            Rcpp::stop("edgeList contains a cycle");
        }
        current = node;
        while (state[current] == 1) {
            state[current] = 2;
            current = parent[current];
        }
    }

    Rcpp::NumericMatrix F(nTips, node_count);
    Rcpp::NumericMatrix G(nTips, node_count);
    for(int i=0; i<nTips; ++i)
        F(i,i) = G(i,i) = 1;

    Rcpp::NumericMatrix D(nTips,nTips);
    Rcpp::NumericMatrix B(nTips,nTips);

    Rcpp::NumericVector tips(nTips);
    for(int i=0; i<tips.size(); ++i)
        tips[i] = i;

    int counter = 0;
    for(int i=0; i + 1 < edgeList.nrow() && tips.size() > 1; i+=2)
        one_step( edgeList(i,1), edgeList(i+1,1), i, i+1, counter++, edgeList, tips, F, G, D, B, rootEdge);

    if (tips.size() != 1) {
        Rcpp::stop("failed to reduce the tree to a single root state");
    }
    if (counter != nTips - 1) {
        Rcpp::stop("unexpected number of independent tree contrasts");
    }
    if (!std::isfinite(rootEdge) || rootEdge <= 0.0) {
        Rcpp::stop("the reduced root contrast has non-positive length");
    }
    
    for(int i=0; i<F.nrow(); ++i){
        D(i,counter) = F(i,tips[0])/std::sqrt(rootEdge);
        B(i,counter) = G(i,tips[0])*std::sqrt(rootEdge);
    }

    //sqrtInvSigma # instead of D
    //sqrtSigma    # instead of B
    
    return( Rcpp::List::create( Rcpp::Named("sqrtInvSigma") = D, Rcpp::Named("sqrtSigma") = B) );
}
