//
// Created by Moskovskaya Elizaveta on 14/05/2017.
//

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZXFactoring.h>

using namespace NTL;

inline bool compare(const Pair<long, long>& first, const Pair<long, long>& second) {
    if (first.b > second.b || (first.b == second.b && first.a > second.a)) {
        return true;
    }
    return false;
}

void intersection(Vec<Pair<long,long>>::iterator begin1, Vec<Pair<long,long>>::iterator end1,
                  Vec<Pair<long,long>>::iterator begin2, Vec<Pair<long,long>>::iterator end2,
                  Vec<Pair<long,long>>& answer, bool (*cmp)(const Pair<long, long>&, const Pair<long, long>&)) {
    while (begin1 != end1 && begin2 != end2) {
        if (cmp(*begin2, *begin1)) {
            ++begin1;
        }
        else if (cmp(*begin1, *begin2)) {
            ++begin2;
        }
        else {
            answer.append(*begin1);
            ++begin1;
            ++begin2;
        }
    }
}

inline ZZ rSieveArrayInit(long a, long b, const ZZ& m, const ZZX& poly, long d) {
    return abs(a + b * m);
}

inline ZZ aSieveArrayInit(long a, long b, const ZZ& m, const ZZX& poly, long d) {
    ZZ res;
    for (long i = 0L; i <= d; ++i) {
        res += coeff(poly, i) * pow(a, i) * pow(-b, d - i);
    }
    return abs(res);
}

void arrayInit(Vec<ZZ>& sievingArray, long B, const ZZ& m, const ZZX& poly, long d, long b,
               ZZ (*pSieveArrayInit) (long, long, const ZZ&, const ZZX&, long)) {
    for (long a = -B; a <= B; ++a) {
        sievingArray.at(B + a) = pSieveArrayInit(a, b, m, poly, d);
    }
}

inline long stepFunction(const Pair<ZZ, ZZ>& ideal, long b, long k, long B) {
    static ZZ border = power(conv<ZZ>(2L), (long)sizeof(long) * 8 - 1) - 1;
    ZZ step = -b * ideal.a + k * ideal.b;
    return step <= border ? conv<long>(step) + B : B + 1;
}

void sieve(Vec<Pair<long, long>>& pairs, const Vec<Pair<ZZ, ZZ>>& FB, Vec<ZZ>& sievingArray, long B, long b) {
    long k;
    long pos;
    for (long i = 0L; i < FB.length(); ++i) {
        k = conv<long>(CeilToZZ((-B + b * conv<RR>(FB.at(i).a)) / conv<RR>(FB.at(i).b)));
        while (stepFunction(FB.at(i), b, k, B) < 0L) {
            ++k;
        }
        pos = stepFunction(FB.at(i), b, k, B);
        while (pos <= 2 * B) {
            while (sievingArray[pos] % FB.at(i).b == 0L and sievingArray[pos] != 0L) {
                sievingArray[pos] /= FB.at(i).b;
            }
            pos = stepFunction(FB.at(i), b, ++k, B);
        }
    }
    for (long i = 0L; i <= 2 * B; ++i) {
        if (sievingArray[i] == 1L) {
            pairs.append(Pair<long, long>(-B + i, b));
        }
    }
}

Vec<Pair<long, long>> sieving(const Vec<Pair<ZZ, ZZ>>& RFB, const Vec<Pair<ZZ, ZZ>>& AFB,
                          long B, const ZZ& m, const ZZX& poly, long degree, long minPairAmount) {
    ZZ (*prSieveArrayInit) (long, long, const ZZ&, const ZZX&, long) = &rSieveArrayInit;
    ZZ (*paSieveArrayInit) (long, long, const ZZ&, const ZZX&, long) = &aSieveArrayInit;
    bool (*cmp)(const Pair<long, long>&, const Pair<long, long>&) = &compare;

    Vec<Pair<long, long>> answer;

    Vec<ZZ> sievingArray;
    sievingArray.SetLength(2 * B + 1);

    long b = 1L;

    while (answer.length() < minPairAmount) {
        Vec<Pair<long, long>> rationalRes;
        arrayInit(sievingArray, B, m, poly, degree, b, prSieveArrayInit);
        sieve(rationalRes, RFB, sievingArray, B, b);

        Vec<Pair<long, long>> algebraicRes;
        arrayInit(sievingArray, B, m, poly, degree, b, paSieveArrayInit);
        sieve(algebraicRes, AFB, sievingArray, B, b);

        intersection(rationalRes.begin(), rationalRes.end(), algebraicRes.begin(), algebraicRes.end(), answer, cmp);

        ++b;
    }

    return answer;
}

