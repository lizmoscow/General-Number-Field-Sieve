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
    std::cout << res << '\n';
    Vec<Pair<long, long>> resanswer;
    for (long i = 0; i < res.length(); ++i) {
        if (!IsZero(res[i])) {
            resanswer.append(answer[i]);
        }
    }
    std::cout << resanswer << '\n';


    auto diff = std::chrono::system_clock::now() - start;
    auto sec = std::chrono::duration_cast<std::chrono::seconds>(diff);
    std::cout << "\nduration: " << sec.count();

    return 0;
}