#include <math.h>
#include <stdlib.h>

#include <R.h>

void threepoint_l1ou(int *Npo, int *npo, int *pNpo, int *dYpo, int *dXpo,
                     int *rootpo, double *transa, double *transb, int *des,
                     int *anc, double *y, double *X, double *output) {

  int N = *Npo;
  int n = *npo;
  int pN = *pNpo;
  int npN = n + pN;
  int dY = *dYpo;
  int dX = *dXpo;
  int r = *rootpo;
  r--;
  double rootEdge = *transa;

  double *logd = (double *)calloc(npN, sizeof(double));
  double *vec11 = (double *)calloc(npN, sizeof(double));
  double *yy = (double *)calloc(npN * dY * dY, sizeof(double));
  double *y1 = (double *)calloc(npN * dY, sizeof(double));
  double *Xy = (double *)calloc(npN * dX * dY, sizeof(double));
  double *X1 = (double *)calloc(npN * dX, sizeof(double));
  double *XX = (double *)calloc(npN * dX * dX, sizeof(double));
  int *zero = (int *)calloc(npN, sizeof(int));

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
  for (int j = 0; j < dY * dY; j++) {
    output[p + j] = yy[ikXY];
    ikXY += npN;
  }
  p += dY * dY;
  ikXY = r;
  for (int j = 0; j < dX; j++) {
    output[p + j] = X1[ikXY];
    ikXY += npN;
  }
  p += dX;
  ikXY = r;
  for (int j = 0; j < dX * dX; j++) {
    output[p + j] = XX[ikXY];
    ikXY += npN;
  }
  p += dX * dX;
  ikXY = r;
  for (int j = 0; j < dX * dY; j++) {
    output[p + j] = Xy[ikXY];
    ikXY += npN;
  }

  free(logd);
  free(vec11);
  free(y1);
  free(yy);
  free(X1);
  free(XX);
  free(Xy);
  free(zero);
}
