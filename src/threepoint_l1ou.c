#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <R_ext/Utils.h>

static size_t checked_size_product(size_t left, size_t right,
                                   const char *label) {
  if (right != 0 && left > ((size_t)-1) / right) {
    error("%s allocation size overflows size_t", label);
  }
  return left * right;
}

static size_t checked_size_sum(size_t left, size_t right,
                               const char *label) {
  if (left > ((size_t)-1) - right) {
    error("%s size overflows size_t", label);
  }
  return left + right;
}

static void *zero_alloc(size_t count, size_t size, const char *label) {
  size_t bytes = checked_size_product(count, size, label);
  void *memory = (void *)R_alloc(count, (int)size);
  memset(memory, 0, bytes);
  return memory;
}

static void check_interrupt_periodically(size_t *counter) {
  ++(*counter);
  if (((*counter) & 0x3fffU) == 0U) {
    R_CheckUserInterrupt();
  }
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
  size_t node_count = checked_size_sum((size_t)n, (size_t)pN, "node");
  if (node_count > (size_t)INT_MAX || (size_t)N != node_count - 1U) {
    error("tree dimensions are inconsistent in threepoint_l1ou");
  }
  size_t dy = (size_t)dY;
  size_t dx = (size_t)dX;
  size_t dY2 = checked_size_product(dy, dy, "response output");
  size_t dX2 = checked_size_product(dx, dx, "design output");
  size_t dXY = checked_size_product(dx, dy, "cross-product output");
  if (dY2 > (size_t)INT_MAX || dX2 > (size_t)INT_MAX ||
      dXY > (size_t)INT_MAX) {
    error("dimension products are too large in threepoint_l1ou");
  }
  size_t output_count = 2U;
  output_count = checked_size_sum(output_count, dy, "output");
  output_count = checked_size_sum(output_count, dY2, "output");
  output_count = checked_size_sum(output_count, dx, "output");
  output_count = checked_size_sum(output_count, dX2, "output");
  output_count = checked_size_sum(output_count, dXY, "output");
  if (output_count > (size_t)INT_MAX) {
    error("requested output is too large in threepoint_l1ou");
  }
  int r = *rootpo;
  r--;
  double rootEdge = *transa;

  size_t y1_count =
      checked_size_product(node_count, dy, "response workspace");
  size_t yy_count =
      checked_size_product(y1_count, dy, "response cross-product");
  size_t X1_count =
      checked_size_product(node_count, dx, "design workspace");
  size_t XX_count =
      checked_size_product(X1_count, dx, "design cross-product");
  size_t Xy_count =
      checked_size_product(X1_count, dy, "design-response cross-product");
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

  for (size_t inode = 0; inode < node_count; inode++) {
    zero[inode] = -1;
  }

  size_t interrupt_counter = 0;
  for (int iedge = 0; iedge < N + 1; iedge++) {
    check_interrupt_periodically(&interrupt_counter);
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

    size_t iY1 = (size_t)di;
    size_t iX1 = (size_t)di;
    size_t iYY = (size_t)di;
    size_t iXY = (size_t)di;
    size_t iXX = (size_t)di;

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
      size_t jY = (size_t)di;
      for (int j = 0; j < dY; j++) {
        y1[iY1] = y[jY] * invel;
        size_t kY = (size_t)di;
        for (int k = 0; k < dY; k++) {
          check_interrupt_periodically(&interrupt_counter);
          yy[iYY] = y1[iY1] * y[kY];
          iYY += node_count;
          kY += (size_t)n;
        }
        size_t kX = (size_t)di;
        for (int k = 0; k < dX; k++) {
          check_interrupt_periodically(&interrupt_counter);
          Xy[iXY] = y1[iY1] * X[kX];
          iXY += node_count;
          kX += (size_t)n;
        }
        iY1 += node_count;
        jY += (size_t)n;
      }
      size_t jX = (size_t)di;
      for (int j = 0; j < dX; j++) {
        X1[iX1] = X[jX] * invel;
        size_t kX = (size_t)di;
        for (int k = 0; k < dX; k++) {
          check_interrupt_periodically(&interrupt_counter);
          XX[iXX] = X1[iX1] * X[kX];
          iXX += node_count;
          kX += (size_t)n;
        }
        iX1 += node_count;
        jX += (size_t)n;
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
        logd[di] += log1p(el * vec11[di]);
        ev = el / (1 + el * vec11[di]);
        ev2 = 1 / (1 + el * vec11[di]);
        for (int j = 0; j < dY; j++) {
          double tmp1 = ev * y1[iY1];
          size_t kY1 = (size_t)di;
          for (int k = 0; k < dY; k++) {
            check_interrupt_periodically(&interrupt_counter);
            yy[iYY] -= tmp1 * y1[kY1];
            iYY += node_count;
            kY1 += node_count;
          }
          size_t kX1 = (size_t)di;
          for (int k = 0; k < dX; k++) {
            check_interrupt_periodically(&interrupt_counter);
            Xy[iXY] -= tmp1 * X1[kX1];
            iXY += node_count;
            kX1 += node_count;
          }
          iY1 += node_count;
        }
        for (int j = 0; j < dX; j++) {
          double tmp1 = ev * X1[iX1];
          size_t kX1 = (size_t)di;
          for (int k = 0; k < dX; k++) {
            check_interrupt_periodically(&interrupt_counter);
            XX[iXX] -= tmp1 * X1[kX1];
            iXX += node_count;
            kX1 += node_count;
          }
          iX1 += node_count;
        }
      } else {
        logd[di] += log(el);
        int d0 = zero[di];
        double fac = 1 / el + vec11[di];
        iY1 = 0;
        iX1 = 0;
        for (int j = 0; j < dY; j++) {
          double tmp1 = fac * y1[(size_t)d0 + iY1];
          size_t kY1 = 0;
          for (int k = 0; k < dY; k++) {
            check_interrupt_periodically(&interrupt_counter);
            yy[iYY] += tmp1 * y1[(size_t)d0 + kY1] -
                y1[(size_t)d0 + iY1] * y1[(size_t)di + kY1] -
                y1[(size_t)di + iY1] * y1[(size_t)d0 + kY1];
            iYY += node_count;
            kY1 += node_count;
          }
          size_t kX1 = 0;
          for (int k = 0; k < dX; k++) {
            check_interrupt_periodically(&interrupt_counter);
            Xy[iXY] += tmp1 * X1[(size_t)d0 + kX1] -
                y1[(size_t)d0 + iY1] * X1[(size_t)di + kX1] -
                y1[(size_t)di + iY1] * X1[(size_t)d0 + kX1];
            iXY += node_count;
            kX1 += node_count;
          }
          iY1 += node_count;
        }
        for (int j = 0; j < dX; j++) {
          double tmp1 = fac * X1[(size_t)d0 + iX1];
          size_t kX1 = 0;
          for (int k = 0; k < dX; k++) {
            check_interrupt_periodically(&interrupt_counter);
            XX[iXX] += tmp1 * X1[(size_t)d0 + kX1] -
                X1[(size_t)d0 + iX1] * X1[(size_t)di + kX1] -
                X1[(size_t)di + iX1] * X1[(size_t)d0 + kX1];
            iXX += node_count;
            kX1 += node_count;
          }
          iX1 += node_count;
        }
      }
      if (goodchildren) {
        iY1 = (size_t)di;
        for (int j = 0; j < dY; j++) {
          y1[iY1] *= ev2;
          iY1 += node_count;
        }
        iX1 = (size_t)di;
        for (int j = 0; j < dX; j++) {
          X1[iX1] *= ev2;
          iX1 += node_count;
        }
        vec11[di] *= ev2;
      } else {
        invel = 1 / el;
        iY1 = 0;
        for (int j = 0; j < dY; j++) {
          y1[(size_t)di + iY1] =
              y1[(size_t)zero[di] + iY1] * invel;
          iY1 += node_count;
        }
        iX1 = 0;
        for (int j = 0; j < dX; j++) {
          X1[(size_t)di + iX1] =
              X1[(size_t)zero[di] + iX1] * invel;
          iX1 += node_count;
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
        y1[(size_t)anci + iY1] += y1[(size_t)di + iY1];
        for (int k = 0; k < dY; k++) {
          check_interrupt_periodically(&interrupt_counter);
          yy[(size_t)anci + iYY] += yy[(size_t)di + iYY];
          iYY += node_count;
        }
        for (int k = 0; k < dX; k++) {
          check_interrupt_periodically(&interrupt_counter);
          Xy[(size_t)anci + iXY] += Xy[(size_t)di + iXY];
          iXY += node_count;
        }
        iY1 += node_count;
      }
      for (int j = 0; j < dX; j++) {
        X1[(size_t)anci + iX1] += X1[(size_t)di + iX1];
        for (int k = 0; k < dX; k++) {
          check_interrupt_periodically(&interrupt_counter);
          XX[(size_t)anci + iXX] += XX[(size_t)di + iXX];
          iXX += node_count;
        }
        iX1 += node_count;
      }
      vec11[anci] += vec11[di];
    }
  }

  output[0] = logd[r];
  output[1] = vec11[r];
  size_t p = 2U;
  size_t ikXY = (size_t)r;
  for (int j = 0; j < dY; j++) {
    check_interrupt_periodically(&interrupt_counter);
    output[p + (size_t)j] = y1[ikXY];
    ikXY += node_count;
  }
  p += dy;
  ikXY = (size_t)r;
  for (size_t j = 0; j < dY2; j++) {
    check_interrupt_periodically(&interrupt_counter);
    output[p + j] = yy[ikXY];
    ikXY += node_count;
  }
  p += dY2;
  ikXY = (size_t)r;
  for (int j = 0; j < dX; j++) {
    check_interrupt_periodically(&interrupt_counter);
    output[p + (size_t)j] = X1[ikXY];
    ikXY += node_count;
  }
  p += dx;
  ikXY = (size_t)r;
  for (size_t j = 0; j < dX2; j++) {
    check_interrupt_periodically(&interrupt_counter);
    output[p + j] = XX[ikXY];
    ikXY += node_count;
  }
  p += dX2;
  ikXY = (size_t)r;
  for (size_t j = 0; j < dXY; j++) {
    check_interrupt_periodically(&interrupt_counter);
    output[p + j] = Xy[ikXY];
    ikXY += node_count;
  }

}
