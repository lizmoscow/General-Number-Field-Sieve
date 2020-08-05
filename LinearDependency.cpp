//
// Created by Moskovskaya Elizaveta on 19/05/2017.
//
#include <NTL/ZZXFactoring.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_GF2.h>
#include <vector>

#include "gauss.h"
#include "norm.h"

using namespace NTL;

typedef mat_GF2 mat_;
typedef vec_GF2 vec_;
typedef GF2 type_;
const ZZ F(8L);

void fillExpMatrix(mat_& matrix, const Vec<Pair<long, long>>& answer,
                   const Vec<Pair<ZZ, ZZ>>& RFB, const Vec<Pair<ZZ, ZZ>>& AFB, const Vec<Pair<ZZ, ZZ>>& QCFB,
                   const ZZX& poly, const ZZ& m, long d) {
    for (long i = 0L; i < answer.length(); ++i) {
        ZZ_pPush push();
        ZZ_p::init(F);
        ZZ s = answer[i].a + m * answer[i].b;
        ZZ n = norm(answer[i].a, answer[i].b, m, poly, d);
        if (s < 0L) {
            matrix[0][i] = type_(1L);
        }
        for (long j = 0L; j < RFB.length(); ++j) {
            long k = 0L;
            while (s != 0 && s % RFB[j].b == 0) {
                ++k;
                s /= RFB[j].b;
            }
            matrix[j + 1][i] = type_(k % 2);
        }
        for (long j = 0L; j < AFB.length(); ++j) {
            long k = 0L;
            if ((answer[i].a + AFB[j].a * answer[i].b) % AFB[j].b == 0) {
                while (n != 0 && n % AFB[j].b == 0) {
                    ++k;
                    n /= AFB[j].b;
                }
            }
            matrix[j + RFB.length() + 1][i] = type_(k % 2);
        }
        for (int j = 0L; j < QCFB.length(); ++j) {
            ZZ_p::init(QCFB[j].b);
            if (power(conv<ZZ_p>(answer[i].a + QCFB[j].a * answer[i].b), (QCFB[j].b - 1) / 2) != 1) {
                matrix[j + AFB.length() + RFB.length() + 1][i] = type_(1L);
            }
        }
    }
}

void gauss (mat_GF2& matrix, Vec<vec_GF2> &vecs) {
    long m = matrix.NumCols();
    long n = matrix.NumRows();
    std::vector<long> where(m, -1);
    for (long col = 0, row = 0; col < m && row < n; ++col) {
        for (long i = row; i < n; ++i) {
            if (!IsZero(matrix[i][col])) {
                swap (matrix[i], matrix[row]);
                break;
            }
        }

        if (IsZero(matrix[row][col])) {
            continue;
        }
        where[col] = row;

        for (long i = 0; i < n; ++i) {
            if (i != row && !IsZero(matrix[i][col])) {
                matrix[i] += matrix[row];
            }

        }
        ++row;
    }

    vec_GF2 vec;
    vec.SetLength(m);
    for (long i = 0; i < m; ++i) {
        if (where[i] == -1) {
            clear(vec);
            for (long j = 0; j < m; ++j) {
                if (where[j] != -1) {
                    vec[j] = matrix[where[j]][i];
                }
                else {
                    if (j == i) {
                        vec[j] = GF2(1);
                    }
                }
            }
            vecs.append(vec);
        }
    }

}


void randComb(const Vec<vec_GF2> &vecs, vec_GF2 &ans) {
    vec_GF2 randomCombination;
    random(randomCombination, vecs.length());

    for (long i = 0; i < vecs[0].length(); ++i) {
        ans[i] = GF2(0);
        for (long j = 0; j < vecs.length(); ++j) {
            if (!IsZero(randomCombination[j])) {
                ans[i] += randomCombination[j] * vecs[j][i];
            }
        }
    }
}


void mulDiag(mat_ &res, const mat_ &matrix, const vec_ & diag) {
    if (matrix.NumCols() != diag.length()) {
        throw "Given vector does not have an appropriate length!";
    }
    for (long i = 0; i < matrix.NumRows(); ++i) {
        for (long j = 0; j < matrix.NumCols(); ++j) {
            res[i][j] = matrix[i][j] * diag[j];
        }
    }
}

bool helpLanczos(vec_ &x, const mat_ &a, const vec_ &b) {
    Vec<vec_> w;
    w.append(b);
    //std::cerr << w[0] << "\n";
    long count = 0;


    vec_ wTemp;
    wTemp.SetLength(a.NumRows());
    vec_ alpha1, alpha2;
    type_ numerator, denominator;
    wTemp = a * w[count];
    alpha1 =  a * w[count];
    alpha2 = a * w[count];
    numerator = alpha1 * alpha2;
    alpha1 = w[count];
    alpha2 = a * w[count];
    denominator = alpha1 * alpha2;
    if (IsZero(denominator)) {
        return false;
    }
    wTemp -= numerator / denominator * w[count];
    w.append(wTemp);
    //std::cerr << w[1] << "\n";
    ++count;

    type_ res = w[count] * a * w[count];
    while (!IsZero(res)) {
        mul(wTemp, a, w[count]);
        for (long i = 0; i <= 1; ++i) {
            alpha1 = a * w[count - i];
            alpha2 = a * w[count - i];
            numerator = alpha1 * alpha2;
            alpha1 = w[count - i];
            alpha2 = a * w[count - i];
            denominator = alpha1 * alpha2;
            if (IsZero(denominator)) {
                return false;
            }
            wTemp -= numerator / denominator * w[count - i];
        }
        w.append(wTemp);
        ++count;
        res = w[count] * a * w[count];
    }
    if (!IsZero(w[count])) {
        return false;
    }
    std::cerr << x.length() << " " << w[0].length() << " " << w[1].length() << '\n';
    for (long i = 0; i < count - 1; ++i) {
        x += (w[i] * b) / (w[i] * w[i]) * w[i];
    }
    return true;
}

void lanczos(const mat_ &matrix, vec_ &x) {
    ZZ_pPush push();
    ZZ_p::init(F);
    long n = matrix.NumRows();
    long m = matrix.NumCols();
    vec_ b0, b;
    b0.SetLength(n);
    for (long i = 0; i < n; ++i) {
        b0[i] =  - matrix[i][m - 1];
    }
    vec_ xx;
    xx.SetLength(m);

    while (true) {
        mat_ a = transpose(matrix);
        b = b0;
        vec_ d;
        d.SetLength(n);
        /*for (long i = 0; i < d.length(); ++i) {
            d[i] = ZZ_p(rand() % 2);
        }*/
        for(long i = 0; i < n; ++i)
        {
            do
            {
                d[i] = random_GF2();
            }
            while(IsZero(d[i]));
            //sqr(D[i], D[i]);
        }
        mulDiag(a, a, d);
        //std::cerr << a.NumRows() << " x " << a.NumCols() << " " << d.length() << " x " << d.length() << " " << b.length();
        b = a * b;
        mul(a, a, matrix);
        if (helpLanczos(xx, a, b)) {
            std::cout << matrix << "\n" << xx << "\n" << b0 << "\n" << a << "\n" << b << "\n";
            vec_ cmp;
            cmp = matrix * xx;
            if (cmp == b0) {
                x = xx;
                return;
            }
        }
    }
}
