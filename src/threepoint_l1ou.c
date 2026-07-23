#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>

static size_t checked_size_product(size_t left, size_t right,
                                   const char *label) {
  if (right != 0 && left > ((size_t)-1) / right) {
    error("%s allocation size overflows size_t", label);
  }
  return left * right;
}

static void *zero_alloc(size_t count, size_t size, const char *label) {
  size_t bytes = checked_size_product(count, size, label);
  void *memory = (void *)R_alloc(count, (int)size);
  memset(memory, 0, bytes);
  return memory;
}

void threepoint_l1ou(int *Npo, int *npo, int *pNpo, int *dYpo, int *dXpo,
                     int *rootpo, double *transa, double *transb, int *des,
                     int *anc, double *y, double *X, double *output) {

  int N = *Npo;
  int n = *npo;
  int pN = *pNpo;
  int dY = *dYpo;
  int dX = *dXpo;
  if (N < 1 || n < 2 || pN < 1 || dY < 1 || dX < 1 ||
      n > INT_MAX - pN) {
    error("invalid dimensions supplied to threepoint_l1ou");
  }
  int npN = n + pN;
  if (N != npN - 1 || dY > INT_MAX / dY || dX > INT_MAX / dX ||
      dX > INT_MAX / dY) {
    error("dimension products are too large in threepoint_l1ou");
  }
  int dY2 = dY * dY;
  int dX2 = dX * dX;
  int dXY = dX * dY;
  long long output_count =
      2LL + dY + dY2 + dX + dX2 + dXY;
  if (output_count > INT_MAX) {
    error("requested output is too large in threepoint_l1ou");
  }
  int r = *rootpo;
  r--;
  double rootEdge = *transa;

  size_t node_count = (size_t)npN;
  size_t y1_count =
      checked_size_product(node_count, (size_t)dY, "response workspace");
  size_t yy_count =
      checked_size_product(y1_count, (size_t)dY, "response cross-product");
  size_t X1_count =
      checked_size_product(node_count, (size_t)dX, "design workspace");
  size_t XX_count =
      checked_size_product(X1_count, (size_t)dX, "design cross-product");
  size_t Xy_count =
      checked_size_product(X1_count, (size_t)dY, "design-response cross-product");
  double *logd =
      (double *)zero_alloc(node_count, sizeof(double), "node workspace");
  double *vec11 =
      (double *)zero_alloc(node_count, sizeof(double), "node workspace");
  double *yy =
      (double *)zero_alloc(yy_count, sizeof(double), "response cross-product");
  double *y1 =
      (double *)zero_alloc(y1_count, sizeof(double), "response workspace");
  double *Xy = (double *)zero_alloc(
      Xy_count, sizeof(double), "design-response cross-product");
  double *X1 =
      (double *)zero_alloc(X1_count, sizeof(double), "design workspace");
  double *XX =
      (double *)zero_alloc(XX_count, sizeof(double), "design cross-product");
  int *zero =
      (int *)zero_alloc(node_count, sizeof(int), "zero-edge workspace");

  for (int iedge = 0; iedge < N + 1; iedge++) {
    zero[iedge] = -1;
  }

  for (int iedge = 0; iedge < N + 1; iedge++) {
    double el;
    double invel;
    int di;
    int anci = 0;
    if (iedge < N) {
      el = transb[iedge];
      di = des[iedge] - 1;
      anci = anc[iedge] - 1;
    } else {
      el = rootEdge;
      di = r;
    }

    int iY1 = di;
    int iX1 = di;
    int iYY = di;
    int iXY = di;
    int iXX = di;

    if (di < n) {
      if (el > 0) {
        invel = 1 / el;
      } else {
        invel = 1.0;
      }
      if (el > 0) {
        logd[di] = log(el);
        vec11[di] = invel;
      } else {
        if (zero[anci] >= 0) {
          error("two or more sister external edges have length 0, V is singular\n (node %d in pruning-wise-ordered tree)", anci + 1);
        } else {
          zero[anci] = di;
        }
      }
      int jY = di;
      for (int j = 0; j < dY; j++) {
        y1[iY1] = y[jY] * invel;
        int kY = di;
        for (int k = 0; k < dY; k++) {
          yy[iYY] = y1[iY1] * y[kY];
          iYY += npN;
          kY += n;
        }
        int kX = di;
        for (int k = 0; k < dX; k++) {
          Xy[iXY] = y1[iY1] * X[kX];
          iXY += npN;
          kX += n;
        }
        iY1 += npN;
        jY += n;
      }
      int jX = di;
      for (int j = 0; j < dX; j++) {
        X1[iX1] = X[jX] * invel;
        int kX = di;
        for (int k = 0; k < dX; k++) {
          XX[iXX] = X1[iX1] * X[kX];
          iXX += npN;
          kX += n;
        }
        iX1 += npN;
        jX += n;
      }
    } else {
      int goodchildren = 1;
      double ev;
      double ev2;
      if (zero[di] >= 0) {
        if (el <= 0) {
          error("One external edge and its parent both have length 0\n (node %d in pruning-wise-ordered tree). To avoid this situation,\n please make a polytomy or resolve it differently ", di + 1);
        }
        goodchildren = 0;
      }
      if (goodchildren) {
        logd[di] += log(1 + el * vec11[di]);
        ev = el / (1 + el * vec11[di]);
        ev2 = 1 / (1 + el * vec11[di]);
        for (int j = 0; j < dY; j++) {
          double tmp1 = ev * y1[iY1];
          int kY1 = di;
          for (int k = 0; k < dY; k++) {
            yy[iYY] -= tmp1 * y1[kY1];
            iYY += npN;
            kY1 += npN;
          }
          int kX1 = di;
          for (int k = 0; k < dX; k++) {
            Xy[iXY] -= tmp1 * X1[kX1];
            iXY += npN;
            kX1 += npN;
          }
          iY1 += npN;
        }
        for (int j = 0; j < dX; j++) {
          double tmp1 = ev * X1[iX1];
          int kX1 = di;
          for (int k = 0; k < dX; k++) {
            XX[iXX] -= tmp1 * X1[kX1];
            iXX += npN;
            kX1 += npN;
          }
          iX1 += npN;
        }
      } else {
        logd[di] += log(el);
        int d0 = zero[di];
        double fac = 1 / el + vec11[di];
        iY1 = 0;
        iX1 = 0;
        for (int j = 0; j < dY; j++) {
          double tmp1 = fac * y1[d0 + iY1];
          int kY1 = 0;
          for (int k = 0; k < dY; k++) {
            yy[iYY] += tmp1 * y1[d0 + kY1] - y1[d0 + iY1] * y1[di + kY1] - y1[di + iY1] * y1[d0 + kY1];
            iYY += npN;
            kY1 += npN;
          }
          int kX1 = 0;
          for (int k = 0; k < dX; k++) {
            Xy[iXY] += tmp1 * X1[d0 + kX1] - y1[d0 + iY1] * X1[di + kX1] - y1[di + iY1] * X1[d0 + kX1];
            iXY += npN;
            kX1 += npN;
          }
          iY1 += npN;
        }
        for (int j = 0; j < dX; j++) {
          double tmp1 = fac * X1[d0 + iX1];
          int kX1 = 0;
          for (int k = 0; k < dX; k++) {
            XX[iXX] += tmp1 * X1[d0 + kX1] - X1[d0 + iX1] * X1[di + kX1] - X1[di + iX1] * X1[d0 + kX1];
            iXX += npN;
            kX1 += npN;
          }
          iX1 += npN;
        }
      }
      if (goodchildren) {
        iY1 = di;
        for (int j = 0; j < dY; j++) {
          y1[iY1] *= ev2;
          iY1 += npN;
        }
        iX1 = di;
        for (int j = 0; j < dX; j++) {
          X1[iX1] *= ev2;
          iX1 += npN;
        }
        vec11[di] *= ev2;
      } else {
        invel = 1 / el;
        iY1 = 0;
        for (int j = 0; j < dY; j++) {
          y1[di + iY1] = y1[zero[di] + iY1] * invel;
          iY1 += npN;
        }
        iX1 = 0;
        for (int j = 0; j < dX; j++) {
          X1[di + iX1] = X1[zero[di] + iX1] * invel;
          iX1 += npN;
        }
        vec11[di] = invel;
      }
    }

    if ((iedge < N) && ((di >= n) || (el > 0))) {
      logd[anci] += logd[di];
      iY1 = 0;
      iX1 = 0;
      iYY = 0;
      iXX = 0;
      iXY = 0;
      for (int j = 0; j < dY; j++) {
        y1[anci + iY1] += y1[di + iY1];
        for (int k = 0; k < dY; k++) {
          yy[anci + iYY] += yy[di + iYY];
          iYY += npN;
        }
        for (int k = 0; k < dX; k++) {
          Xy[anci + iXY] += Xy[di + iXY];
          iXY += npN;
        }
        iY1 += npN;
      }
      for (int j = 0; j < dX; j++) {
        X1[anci + iX1] += X1[di + iX1];
        for (int k = 0; k < dX; k++) {
          XX[anci + iXX] += XX[di + iXX];
          iXX += npN;
        }
        iX1 += npN;
      }
      vec11[anci] += vec11[di];
    }
  }

  output[0] = logd[r];
  output[1] = vec11[r];
  int p = 2;
  int ikXY = r;
  for (int j = 0; j < dY; j++) {
    output[p + j] = y1[ikXY];
    ikXY += npN;
  }
  p += dY;
  ikXY = r;
  for (int j = 0; j < dY2; j++) {
    output[p + j] = yy[ikXY];
    ikXY += npN;
  }
  p += dY2;
  ikXY = r;
  for (int j = 0; j < dX; j++) {
    output[p + j] = X1[ikXY];
    ikXY += npN;
  }
  p += dX;
  ikXY = r;
  for (int j = 0; j < dX2; j++) {
    output[p + j] = XX[ikXY];
    ikXY += npN;
  }
  p += dX2;
  ikXY = r;
  for (int j = 0; j < dXY; j++) {
    output[p + j] = Xy[ikXY];
    ikXY += npN;
  }

}
