//
// Created by Moskovskaya Elizaveta on 13/05/2017.
//

#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZX.h>
#include <NTL/ZZXFactoring.h>
#include <vector>
#include <NTL/mat_GF2.h>

#include "delta.h"
#include "degreeInit.h"
#include "polyInit.h"
#include "RFBGen.h"
#include "AFBGen.h"
#include "QCFBGen.h"
#include "sieving.h"
#include "fillExpMatrix.h"
#include "gauss.h"
#include "fieldsFinder.h"
#include "ShanksTonelli.h"
#include "ChineeseRemainder.h"
#include "computeY.h"

using namespace NTL;

int main( ) {
    ZZ n;
    std::cin >> n;

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
    mat_GF2 matrix;
    matrix.SetDims(answer.length(), minPairAmount - 1);
    fillExpMatrix(matrix, answer, RFB, AFB, QCFB, poly, m, degree);
    vec_GF2 res;
    res.SetLength(answer.length());
    gauss(matrix, res);
    std::cout << "answer: " << res << '\n';
    Vec<Pair<long, long>> result;
    for (long i = 0; i < res.length(); ++i) {
        if (!IsZero(res[i])) {
            result.append(answer[i]);
        }
    }
    std::cout << "pairs, multiplication of which gives a square: " << result << '\n';

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

    auto diff = std::chrono::system_clock::now() - start;
    auto sec = std::chrono::duration_cast<std::chrono::seconds>(diff);
    std::cout << "\nduration: " << sec.count();

    return 0;
}
