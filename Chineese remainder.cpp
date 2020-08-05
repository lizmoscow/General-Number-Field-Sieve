//
// Created by Moskovskaya Elizaveta on 14/06/2017.
//
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/RR.h>

using namespace NTL;

ZZ fieldF(const ZZ &prev, const ZZX &poly) {
    ZZ myPrime(NextPrime(prev + 1));

    ZZ_pX one;
    SetCoeff(one, 0);
    one.normalize();

    ZZX a;
    SetCoeff(a, conv<long>(myPrime));
    SetCoeff(a, 1, -1);
    a.normalize();
    ZZX remainder;
    rem(remainder, a, poly);

    ZZ_pPush push;
    ZZ_p::init(myPrime);
    ZZ_pX r = GCD(conv<ZZ_pX>(remainder), conv<ZZ_pX>(poly));

    while (r != one) {
        SetCoeff(a, conv<long>(myPrime), 0);
        myPrime = NextPrime(myPrime + 1);
        SetCoeff(a, conv<long>(myPrime));
        a.normalize();
        rem(remainder, a, poly);

        ZZ_pPush push;
        ZZ_p::init(myPrime);
        r = GCD(conv<ZZ_pX>(remainder), conv<ZZ_pX>(poly));
    }
    return myPrime;
}


long logOrder(const ZZ& field, const ZZX& poly, const ZZ_pE& lambda) {
    ZZ_pX uno;
    SetCoeff(uno, 0);
    uno.normalize();
    ZZ_pE one = conv<ZZ_pE>(uno);
    uno.kill();

    if (lambda == conv<ZZ_pE>(one)) {
        return 0;
    }
    ZZ pow = power(field, deg(poly)) - 1;
    long count = 1;
    for (long i = 2; i < TruncToZZ(sqrt(conv<RR>(pow))); i *= 2, ++count) {
        if (pow % i == 0 && power(lambda, i) == one) {
            return count;
        }
    }
    return 0;
}

long order(const ZZ& field, const ZZX& poly, const ZZ_pE& lambda) {
    ZZ_pX uno;
    SetCoeff(uno, 0);
    uno.normalize();
    ZZ_pE one = conv<ZZ_pE>(uno);
    uno.kill();

    if (lambda == conv<ZZ_pE>(one)) {
        return 0;
    }
    ZZ pow = power(field, deg(poly)) - 1;
    for (long i = 2; i < TruncToZZ(sqrt(conv<RR>(pow))); i++) {
        if (pow % i == 0 && power(lambda, i) == one) {
            return i;
        }
    }
    return 0;
}



ZZ_pX ShanksTonelli(const ZZ& field, const Vec<Pair<long, long>>& pairs, const ZZX& poly, long d) {
    ZZ_pPush push;
    ZZ_p::init(field);
    ZZX dpoly2 = conv<ZZX>(sqr(diff(poly)));
    for (long i = 0L; i < pairs.length(); ++i) {
        ZZX a;
        SetCoeff(a, 0L, pairs[i].a);
        SetCoeff(a, 1L, pairs[i].b);
        a.normalize();
        rem(dpoly2, dpoly2 * a, poly);
        //dpoly2 *= a;
    }
    ZZ_pX delta0(conv<ZZ_pX>(dpoly2));
    ZZ q(power(field, d));
    ZZ q1(q - 1L);
    long r = 0L;
    while (q1 % 2L == 0L) {
        q1 /= 2L;
        ++r;
    }
    long s = conv<long>(q1);
    q1.kill();

    ZZ_pEPush pushE;
    ZZ_pE::init(conv<ZZ_pX>(poly));
    ZZ_pE delta(conv<ZZ_pE>(delta0));
    ZZ_pE deltaS = power(delta, s);
    ZZ_pX theta;
    SetCoeff(theta, 0);
    SetCoeff(theta, 1);
    theta.normalize();
    ZZ_pE zeta(power(conv<ZZ_pE>(theta), s));
    theta.kill();
    ZZ_pE lambda(deltaS);
    ZZ_pE omega(power(delta, (s + 1L) / 2L));
    ZZ_pX one;
    SetCoeff(one, 0);
    one.normalize();
    long count = 0;
    while (lambda != conv<ZZ_pE>(one) && count < r) {
        long m = logOrder(field, poly, lambda);
        //long m = (long)(log(conv<ZZ>(order(field, poly, lambda))) / log(conv<ZZ>(2L)));
        if (m == 0) {
            throw "Problems with ordering!";
        }
        if (m > r - 1) {
            return ZZ_pX();
        }
        lambda *= power(zeta, conv<long>(power(conv<ZZ>(2L), r - m)));
        omega *= power(zeta, conv<long>(power(conv<ZZ>(2L), r - m - 1)));
        ++count;
    }
    if (count >= r) {
        return ZZ_pX();
    }
    ZZ_pX nu = conv<ZZ_pX>(omega);

    ZZ n0 = conv<ZZ>(rep(power(omega, (conv<long>(power(field, d)) - 1L) / (conv<long>(field) - 1L)))[0]);
    ZZ n1 = conv<ZZ>(rep(power(delta, (conv<long>(power(field, d)) - 1L) / (conv<long>(field) - 1L)))[0]);

    if (n0 > (field - 1) / 2) {
        n0 -= field;
    }
    if (n1 > (field - 1) / 2) {
        n1 -= field;
    }
    if (n0 * n1 < 0) {
        nu = conv<ZZ_pX>(-omega);
    }
    return nu;
}


void fieldsFinder(Vec<ZZ> &fields, const ZZX &poly, const ZZ &n) {
    ZZ field(9800L);
    ZZ mul(1L);
    while (mul < n) {
        field = fieldF(field, poly);
        fields.append(field);
        mul *= field;
    }
}

ZZ eval(const ZZX &poly, const ZZ &x) {
    ZZ res(0L);
    for (long i = 0; i <= deg(poly); ++i) {
        res += poly[i] * power(x, i);
    }
    return res;
}

ZZ ChineeseRemainder(const Vec<ZZ> &fields, const Vec<ZZ_pX> &roots, const ZZ &m, const ZZ&n) {
    ZZ p(1L);
    for (long i = 0; i < fields.length(); ++i) {
        p *= fields[i];
    }
    Vec<ZZ> ps;
    Vec<ZZ> as;
    Vec<ZZ> xs;
    ZZ s, t, d;
    RR rr;
    for (long i = 0; i < fields.length(); ++i) {
        ps.append(p / fields[i]);
        XGCD(d, s, t, p / fields[i], fields[i]);
        as.append(s % fields[i]);
        ZZ_pPush push;
        ZZ_p::init(fields[i]);
        ZZX root = conv<ZZX>(roots[i]);
        xs.append((eval(root, m)) % fields[i]);
        rr += conv<RR>(as[i] * xs[i]) / conv<RR>(fields[i]);
    }
    ZZ r = TruncToZZ(rr);
    ZZ x = (-r * p) % n;
    for (long i = 0; i < fields.length(); ++i) {
        x = (x + as[i] * xs[i] * ps[i]) % n;
    }
    return x;
}


ZZ computeY(const Vec<Pair<long, long>> &pairs, const ZZX &poly, const ZZ &m, const ZZ &n) {
    ZZ y1(1L);
    for (long i = 0; i < pairs.length(); ++i) {
        y1 *= pairs[i].a + pairs[i].b * m;
    }
    ZZX dpoly = diff(poly);
    ZZ y2 = eval(dpoly, m);
    return (TruncToZZ(sqrt(conv<RR>(y1))) * y2) % n;
}