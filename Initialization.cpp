//
// Created by Moskovskaya Elizaveta on 14/05/2017.
//

#include <NTL/RR.h>
#include <NTL/ZZX.h>
#include <NTL/ZZXFactoring.h>

using namespace NTL;

inline double epsilon( ) {
    static double epsilon = 0.000001;
    return epsilon;
}

inline double delta( ) {
    static double delta = pow(8.0 / 9, 1.0 / 3) + epsilon();
    return delta;
}

inline ZZ eval(const ZZX& poly, const ZZ& point) {
    ZZ value = ZZ(0);
    for (long i = 0L; i <= deg(poly); ++i) {
        value += coeff(poly, i) * power(point, i);
    }
    return value;
}

long degreeInit(const ZZ& n) {
    long degree = conv<long>(TruncToZZ(pow(2 / RR(delta()), RR(1) / 2) *
                                       pow(conv<RR>(log(n)) / log(conv<RR>(log(n))), RR(1) / 3)));
    return degree % 2 == 1 ? degree : degree + 1;
}

ZZX polyInit(const ZZ& n, const ZZ& m, long d) {
    ZZX poly;
    ZZ temp = ZZ(n);
    ZZ coef;
    for (long i = 0L; i < d; ++i) {
        coef = TruncToZZ(conv<RR>(temp) / conv<RR>(power(m, d - i)));
        temp -= power(m, d - i) * coef;
        SetCoeff(poly, d - i, coef);
    }
    SetCoeff(poly, 0, temp);
    return poly;
}

Vec<Pair<ZZ, ZZ>> RFBGen(const ZZ &m, const ZZ &z) {
    Vec<Pair<ZZ, ZZ>> RFB;
    ZZ prime = ZZ(2);
    while (prime <= z) {
        RFB.append(Pair<ZZ,ZZ>(m % prime, prime));
        prime = NextPrime(prime + 1);
    }
    return RFB;
}

Vec<Pair<ZZ, ZZ>> AFBGen(const ZZ &m, const ZZ &z, const ZZX &poly) {
    Vec<Pair<ZZ, ZZ>> AFB;
    ZZ border = NextPrime(3 * z);
    ZZ prime = ZZ(2);
    while (prime < border) {
        for (ZZ q = ZZ(0); q < prime; ++q) {
            if (IsZero(eval(poly, q) % prime)) {
                AFB.append(Pair<ZZ, ZZ>(q % prime, prime));
            }
        }
        prime = NextPrime(prime + 1);
    }
    return AFB;
}

Vec<Pair<ZZ, ZZ>> QCFBGen(const ZZ& m, const ZZ& z, const ZZX& poly) {
    Vec<Pair<ZZ, ZZ>> QCFB;
    ZZ prime = NextPrime(3.4 * z);
    ZZ border = prime + z / 2;
    while (prime < border) {
        for (ZZ q = ZZ(0); q < prime; ++q) {
            if (IsZero(eval(poly, q) % prime)) {
                QCFB.append(Pair<ZZ, ZZ>(q % prime, prime));
            }
        }
        prime = NextPrime(prime + 1);
    }
    return QCFB;
}