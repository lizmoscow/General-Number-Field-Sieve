//
// Created by Moskovskaya Elizaveta on 19/06/2017.
//

#include <iostream>
#include <fstream>
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
#include "randComb.h"
#include "fieldsFinder.h"
#include "ShanksTonelli.h"
#include "ChineeseRemainder.h"
#include "computeY.h"

using namespace NTL;

void initialization(ZZ n, long &degree, ZZ &m, ZZX &poly, Vec<Pair<long, long>> &answer, Vec<vec_GF2> &vecs,
                    Vec<ZZ> &fields) {
    degree = degreeInit(n);
    //std::cout << degree << '\n';
    m = TruncToZZ(pow(conv<RR>(n), RR(1) / conv<RR>(degree)));
    //std::cout << "m = " << m << '\n';
    poly = polyInit(n, m, degree);
    //std::cout << "polynomial: " << poly << '\n';
    ZZ z = NextPrime(TruncToZZ(exp(pow(0.5 * log(n) * log(conv<RR>(log(n))),RR(1) / 2))));
    //std::cout << "z = " << z << '\n';
    Vec<Pair<ZZ, ZZ>> RFB = RFBGen(m, z);
    //std::cout << "rational factor base: " << RFB << '\n';
    Vec<Pair<ZZ, ZZ>> AFB = AFBGen(m, z, poly);
    //std::cout << "algebraic factor base: " << AFB << '\n';
    Vec<Pair<ZZ, ZZ>> QCFB = QCFBGen(m, z, poly);
    //std::cout << "quadratic character factor base: " << QCFB << '\n';
    long minPairAmount = RFB.length() + AFB.length() + QCFB.length() + 2L;
    //std::cout << "minimal amount of pairs required = " << minPairAmount << '\n';
    long B = conv<long>(TruncToZZ(exp(RR(delta()) * pow(conv<RR>(log(n))
                                                        * pow(log(conv<RR>(log(n))), RR(2.0)), RR(1.0) / RR(3.0)))));
    //std::cout << "sieving border = " << B << '\n';
    sieving(answer, RFB, AFB, B, m, poly, degree, minPairAmount);
    //std::cout << "smooth pairs: " << answer << '\n';

    mat_GF2 matrix;
    matrix.SetDims(minPairAmount - 1, answer.length());
    fillExpMatrix(matrix, answer, RFB, AFB, QCFB, poly, m, degree);
    gauss(matrix, vecs);
    fieldsFinder(fields, poly, n);
    //std::cout << "fields: " << fields << '\n';
}


Pair<ZZ, ZZ> solving(const ZZ &n, const long &degree, const ZZ &m, const ZZX &poly,
                     const Vec<Pair<long, long>> &answer, Vec<vec_GF2> &vecs, const Vec<ZZ> &fields) {
    vec_GF2 res;
    res.SetLength(vecs[0].length());
    randComb(vecs, res);
    Vec<Pair<long, long>> result;
    for (long i = 0; i < res.length(); ++i) {
        if (!IsZero(res[i])) {
            result.append(answer[i]);
        }
    }
    //std::cout << "pairs, multiplication of which gives a square: " << result << '\n';

    Vec<ZZ_pX> roots;
    for (long i = 0; i < fields.length(); ++i) {
        roots.append(ShanksTonelli(fields[i], result, poly, degree));
    }
    //std::cout << "roots: " << roots << '\n';

    ZZ x = ChineeseRemainder(fields, roots, m, n);
    ZZ y = computeY(result, poly, m, n);
    //std::cout << "x = " << x << " y = " << y << '\n';
    ZZ x_ = GCD(x - y, n);
    ZZ y_ = GCD(x + y, n);
    //std:: cout << "x_ = " << x_ << " y_ = " << y_ << "\n";
    return Pair<ZZ, ZZ>(x_, y_);
}


ZZ GNFS(const ZZ &n) {
    long degree;
    ZZ m;
    ZZX poly;
    Vec<Pair<long, long>> answer;
    Vec<vec_GF2> vecs;
    Vec<ZZ> fields;

    initialization(n, degree, m, poly, answer, vecs, fields);

    Pair<ZZ, ZZ> res(ZZ(1), ZZ(1));
    while ((IsOne(res.a) || res.a == n) && (IsOne(res.b)) || res.b == n) {
        res = solving(n, degree, m, poly, answer, vecs, fields);
    }
    return (IsOne(res.a)) ? res.b : res.a;
}