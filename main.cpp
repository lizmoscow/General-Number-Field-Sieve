//
// Created by Moskovskaya Elizaveta on 13/05/2017.
//

#include <iostream>
#include <fstream>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZX.h>
#include <NTL/ZZXFactoring.h>
#include <vector>
#include <NTL/mat_GF2.h>
//#include <linbox/matrix/sparse-matrix.h>
//#include <linbox/solutions/solve.h>
//#include <linbox/solutions/methods.h>

#include "delta.h"
#include "degreeInit.h"
#include "polyInit.h"
#include "RFBGen.h"
#include "AFBGen.h"
#include "QCFBGen.h"
#include "sieving.h"
#include "norm.h"
#include "mulDiag.h"
#include "fillExpMatrix.h"
#include "lanczos.h"
#include "gauss.h"
#include "fieldsFinder.h"
#include "ShanksTonelli.h"
#include "ChineeseRemainder.h"
#include "computeY.h"

using namespace NTL;

/*typedef Givaro::Modular<long> Field;
void fillExpSpMatrix(LinBox::SparseMatrix<Field>& matrix, const Vec<Pair<long, long>> &answer,
                     const Vec<Pair<ZZ, ZZ>>& RFB, const Vec<Pair<ZZ, ZZ>>& AFB,
                     const Vec<Pair<ZZ, ZZ>>& QCFB, const ZZX& poly, const ZZ& m, long d);*/

int main( ) {
    try {
        ZZ n;
        std::cin >> n;
        long s = 0;
        //std::cin >> s;

        auto start = std::chrono::system_clock::now();

        long degree = degreeInit(n);
        std::cout << degree << '\n';
        ZZ m = TruncToZZ(pow(conv<RR>(n), RR(1) / conv<RR>(degree)));
        std::cout << "m = " << m << '\n';
        ZZX poly = polyInit(n, m, degree);
        std::cout << "polynomial: " << poly << '\n';
        ZZ z = NextPrime(TruncToZZ(exp(pow(0.5 * log(n) * log(conv<RR>(log(n))),RR(1) / 2))));
        std::cout << "z = " << z << '\n';
        Vec<Pair<ZZ, ZZ>> RFB = RFBGen(m, z);
        std::cout << "rational factor base: " << RFB << '\n';
        Vec<Pair<ZZ, ZZ>> AFB = AFBGen(m, z, poly);
        std::cout << "algebraic factor base: " << AFB << '\n';
        Vec<Pair<ZZ, ZZ>> QCFB = QCFBGen(m, z, poly);
        std::cout << "quadratic character factor base: " << QCFB << '\n';
        long minPairAmount = RFB.length() + AFB.length() + QCFB.length() + 2L;
        std::cout << "minimal amount of pairs required = " << minPairAmount << '\n';
        long B = conv<long>(TruncToZZ(exp(RR(delta()) * pow(conv<RR>(log(n))
                                                            * pow(log(conv<RR>(log(n))), RR(2.0)), RR(1.0) / RR(3.0)))));
        std::cout << "sieving border = " << B << '\n';
        Vec<Pair<long, long>> answer = sieving(RFB, AFB, B, m, poly, degree, minPairAmount);
        std::cout << "smooth pairs: " << answer << '\n';


//-------------- Old linear algebra part -------------------------
        mat_GF2 matrix, mCopy;
        matrix.SetDims(minPairAmount - 1, answer.length());
        fillExpMatrix(matrix, answer, RFB, AFB, QCFB, poly, m, degree);
        mCopy = matrix;
        vec_GF2 res;
        res.SetLength(answer.length());
        gauss(matrix, res, s);
        std::cout << matrix << "\n";
        std::cout << res << "\n";
        vec_GF2 vec;
        vec.SetLength(matrix.NumRows());
        vec_GF2 colCopy;
        colCopy.SetLength(matrix.NumRows());
        for (long i = 0; i < matrix.NumCols(); ++i) {
            for (long j = 0; j < matrix.NumRows(); ++j) {
                colCopy[j] = mCopy[j][i];
            }
            vec += (res[i]) * colCopy;
        }
        //std::cout << "vec = " << vec << "\n";
        Vec<Pair<long, long>> result;
        for (long i = 0; i < res.length(); ++i) {
            if (!IsZero(res[i])) {
                result.append(answer[i]);
            }
        }
        std::cout << "pairs, multiplication of which gives a square: " << result << '\n';
//-----------------------------------------------------------------


//-------------------- New linear algebra part --------------------
        /*Field f(2);
        LinBox::SparseMatrix<Field> matrix(f, minPairAmount - 1, answer.length());
        fillExpSpMatrix(matrix, answer, RFB, AFB, QCFB, poly, m, degree);*/

        /*std::fstream out;
        out.open("/Users/liza_moskovskaya/General-Number-Field-Sieve/matrix.txt");
        matrix.write(out);
        out.close();
        std::cout << matrix.rowdim() << " * " << matrix.coldim() << '\n';*/

/*    LinBox::DenseMatrix<Field> NullSpace(f,(size_t)(matrix.rowdim()),(size_t)(matrix.coldim()));
    LinBox::GaussDomain<Field> GD(f);
    GD.nullspacebasisin(NullSpace, matrix);
    std::cerr << "NullsSpace dimensions:" << NullSpace.rowdim() << 'x' << NullSpace.coldim() << std::endl;*/

        /*int64_t p = (int64_t)(matrix.coldim());
        int64_t t = (int64_t)(matrix.rowdim());
        LinBox::DenseVector<Field> v(f, p), b (f, t);
        LinBox::solve(v, matrix, b, LinBox::Method::Lanczos());
        std::cout << "(Lanczos) Solution is [";
        for(auto it=v.begin();it != v.end(); ++it)
            f.write(std::cout, *it) << " ";
        std::cout << "]" << std::endl;*/
//-----------------------------------------------------------------


        Vec<ZZ> fields;
        Vec<ZZ_pX> roots;
        fieldsFinder(fields, roots, poly, result, n, degree);
        std::cout << "fields: " << fields << '\n';
        std::cout << "roots: " << roots << '\n';

/*    fields.append(conv<ZZ>(9851L));
    fields.append(conv<ZZ>(9907L));
    fields.append(conv<ZZ>(9929L));

    ZZ_pPush push();
    ZZ_p::init(fields[0]);
    ZZ_pX polynimial1;
    SetCoeff(polynimial1, 0, conv<ZZ_p>(4716L));
    SetCoeff(polynimial1, 1, conv<ZZ_p>(6752L));
    SetCoeff(polynimial1, 2, conv<ZZ_p>(4458L));
    polynimial1.normalize();
    roots.append(polynimial1);
    ZZ_pPush push();
    ZZ_p::init(fields[1]);
    ZZ_pX polynimial2;
    SetCoeff(polynimial2, 0, conv<ZZ_p>(6847L));
    SetCoeff(polynimial2, 1, conv<ZZ_p>(3088L));
    SetCoeff(polynimial2, 2, conv<ZZ_p>(1218L));
    polynimial2.normalize();
    roots.append(polynimial2);
    ZZ_pPush push();
    ZZ_p::init(fields[2]);
    ZZ_pX polynimial3;
    SetCoeff(polynimial3, 0, conv<ZZ_p>(1257L));
    SetCoeff(polynimial3, 1, conv<ZZ_p>(1168L));
    SetCoeff(polynimial3, 2, conv<ZZ_p>(2141L));
    polynimial3.normalize();
    roots.append(polynimial3);

    m = 31L;
    poly[0] = 8L;
    poly[1] = 29L;
    poly[2] = 15L;
    poly[3] = 1L;
    result.append(Pair<long, long>(-8, 3));
    result.append(Pair<long, long>(-8, 7));
    result.append(Pair<long, long>(-5, 2));
    result.append(Pair<long, long>(-2, 1));
    result.append(Pair<long, long>(2, 1));
    result.append(Pair<long, long>(2, 3));
    result.append(Pair<long, long>(4, 1));
    result.append(Pair<long, long>(4, 11));
    result.append(Pair<long, long>(11, 7));
    result.append(Pair<long, long>(13, 1));
    result.append(Pair<long, long>(19, 4));
    result.append(Pair<long, long>(24, 55));
    result.append(Pair<long, long>(25, 2));
    result.append(Pair<long, long>(44, 9));
    result.append(Pair<long, long>(104, 1));
    result.append(Pair<long, long>(119, 11));
    std::cout << "pairs, multiplication of which gives a square: " << result << '\n';
    std::cout << "fields: " << fields << '\n';
    roots.append(ShanksTonelli(fields[0], result, poly, degree));
    roots.append(ShanksTonelli(fields[1], result, poly, degree));
    roots.append(ShanksTonelli(fields[2], result, poly, degree));
    std::cout << "roots: " << roots << '\n';*/

        ZZ x = ChineeseRemainder(fields, roots, m, n);
        ZZ y = computeY(result, poly, m, n);
        std::cout << "x = " << x << " y = " << y << '\n';
        ZZ x_ = GCD(x - y, n);
        ZZ y_ = GCD(x + y, n);
        std:: cout << "x_ = " << x_ << " y_ = " << y_ << "\n";

        auto diff = std::chrono::system_clock::now() - start;
        auto sec = std::chrono::duration_cast<std::chrono::seconds>(diff);
        std::cout << "\nduration: " << sec.count();
    }
    catch(const char *e) {
        std::cerr << e;
    }

    return 0;
}


/*void fillExpSpMatrix(LinBox::SparseMatrix<Field>& matrix, const Vec<Pair<long, long>>& answer,
                     const Vec<Pair<ZZ, ZZ>>& RFB, const Vec<Pair<ZZ, ZZ>>& AFB, const Vec<Pair<ZZ, ZZ>>& QCFB,
                     const ZZX& poly, const ZZ& m, long degree) {

    Givaro::Modular<short> f(2);
    for (long i = 0L; i < answer.length(); ++i) {
        ZZ s = answer[i].a + m * answer[i].b;
        ZZ n = norm(answer[i].a, answer[i].b, m, poly, degree);
        if (s < 0L) {
            matrix.setEntry(0, i, 1);
        }
        for (long j = 0L; j < RFB.length(); ++j) {
            long k = 0L;
            while (s != 0 && s % RFB[j].b == 0) {
                ++k;
                s /= RFB[j].b;
            }
            matrix.setEntry(j + 1, i, k % 2);
        }
        for (long j = 0L; j < AFB.length(); ++j) {
            long k = 0L;
            if ((answer[i].a + AFB[j].a * answer[i].b) % AFB[j].b == 0) {
                while (n != 0 && n % AFB[j].b == 0) {
                    ++k;
                    n /= AFB[j].b;
                }
            }
            matrix.setEntry(j + RFB.length() + 1, i, k % 2);
        }
        for (int j = 0L; j < QCFB.length(); ++j) {
            ZZ_p::init(QCFB[j].b);
            if (power(conv<ZZ_p>(answer[i].a + QCFB[j].a * answer[i].b), (QCFB[j].b - 1) / 2) != 1) {
                matrix.setEntry(j + AFB.length() + RFB.length() + 1, i, 1);
            }
        }
    }
}*/
