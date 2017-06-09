//
// Created by Moskovskaya Elizaveta on 19/05/2017.
//
#include <NTL/ZZXFactoring.h>
#include <NTL/mat_GF2.h>
#include <vector>

#include "gauss.h"
#include "norm.h"

using namespace NTL;

void fillExpMatrix(mat_GF2& matrix, const Vec<Pair<long, long>>& answer,
                   const Vec<Pair<ZZ, ZZ>>& RFB, const Vec<Pair<ZZ, ZZ>>& AFB, const Vec<Pair<ZZ, ZZ>>& QCFB,
                   const ZZX& poly, const ZZ& m, long d) {
    for (long i = 0L; i < answer.length(); ++i) {
        ZZ s = answer[i].a + m * answer[i].b;
        ZZ n = norm(answer[i].a, answer[i].b, m, poly, d);
        if (s < 0L) {
            matrix[i][0] = GF2(1L);
        }
        for (long j = 0L; j < RFB.length(); ++j) {
            long k = 0L;
            while (s != 0 && s % RFB[j].b == 0) {
                ++k;
                s /= RFB[j].b;
            }
            matrix[i][j + 1] = GF2(k);
        }
        for (long j = 0L; j < AFB.length(); ++j) {
            long k = 0L;
            if ((answer[i].a + AFB[j].a * answer[i].b) % AFB[j].b == 0) {
                while (n != 0 && n % AFB[j].b == 0) {
                    ++k;
                    n /= AFB[j].b;
                }
            }
            matrix[i][j + RFB.length() + 1] = GF2(k);
        }
        for (int j = 0L; j < QCFB.length(); ++j) {
            ZZ_p::init(QCFB[j].b);
            if (power(conv<ZZ_p>(answer[i].a + QCFB[j].a * answer[i].b), (QCFB[j].b - 1) / 2) != 1) {
                matrix[i][j + AFB.length() + RFB.length() + 1] = GF2(1L);
            }
        }
    }
}

bool gauss (mat_GF2& matrix, vec_GF2 & ans) {
    long m = matrix.NumCols();
    long n = matrix.NumRows();
    for (long col=0, row=0; col<m && row<n; ++col) {
        for (long i=row; i<n; ++i) {
            if (!IsZero(matrix[i][col])) {
                ans[i] = 1;
                swap (matrix[i], matrix[row]);
                break;
            }
        }

        if (IsZero(matrix[row][col])) {
            continue;
        }
        for (int i=0; i<n; ++i) {
            if (i != row && !IsZero(matrix[i][col])) {
                matrix[i] += matrix[row];
            }

        }
        ++row;
    }
}