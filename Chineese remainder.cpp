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

ZZ fieldF(ZZ &prev, ZZX &poly) {
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


long logOrder(ZZ& field, ZZX& poly, ZZ_pE& lambda) {
    ZZ_pX one;
    SetCoeff(one, 0);
    one.normalize();
    if (lambda == conv<ZZ_pE>(one)) {
        return 0;
    }
    ZZ pow = power(field, deg(poly)) - 1;
    std::cerr << " pow = " << pow;
    long count = 1;
    std::cerr << " count = " << count;
    for (long i = 2; i < TruncToZZ(sqrt(conv<RR>(pow))); i *= 2) {
        if (pow % i == 0 && power(lambda, i) == conv<ZZ_pE>(one)) {
            return count;
        }
        ++count;
        std::cerr << " lambda = " << power(lambda, i) << " count = " << count;
    }
    return 0;
}

ZZ_pX ShanksTonelli(ZZ& field, Vec<Pair<long, long>>& pairs, ZZX& poly, long d) {
    ZZ_pPush push;
    ZZ_p::init(field);
    ZZ_pX dpoly2 = conv<ZZ_pX>(sqr(diff(poly)));
    for (int i = 0; i < pairs.length(); ++i) {
        ZZ_pX a;
        SetCoeff(a, 0, pairs[i].a);
        SetCoeff(a, 1, pairs[1].b);
        a.normalize();
        dpoly2 *= a;
    }
    ZZ q(power(field, d));
    std::cerr << " q = " << q;
    ZZ q1(q - 1L);
    long r = 0L;
    while (q1 % 2L == 0L) {
        q1 /= 2L;
        ++r;
    }
    long s = conv<long>(q1);
    q1.kill();
    std::cerr << " s = " << s << " r = " << r;


    ZZ_pEPush pushE;
    ZZ_pE::init(conv<ZZ_pX>(poly));
    ZZ_pE delta(conv<ZZ_pE>(dpoly2));
    /*for (int i = 0; i < deg(dpoly2); ++i) {
        SetCoeff(delta, i, conv<ZZ_p>(dpoly2[i]));
    }*/
    std::cerr << " delta = " << delta;
    ZZ_pE deltaS = power(delta, s);
    std::cerr << " deltaS = " << deltaS;
    ZZ_pX theta;
    SetCoeff(theta, 0);
    SetCoeff(theta, 1);
    theta.normalize();
    ZZ_pE zeta(power(conv<ZZ_pE>(theta), s));
    theta.kill();
    std::cerr << " zeta = " << zeta;
    ZZ_pE lambda(deltaS);
    std::cerr << " lambda = " << lambda;
    ZZ_pE omega(power(delta, (s + 1L) / 2L));
    std::cerr << " omega = " << omega << "\n";
    ZZ_pX one;
    SetCoeff(one, 0);
    one.normalize();
    long count = 0;
    while (lambda != conv<ZZ_pE>(one) && count < r) {
        long m = logOrder(field, poly, lambda);
        std::cerr << " m = " << m;
        if (m == 0) {
            throw "Problems with ordering!";
        }
        if (m > r - 1) {
            std::cout << '\n';
            return ZZ_pX();
        }
        lambda *= power(zeta, conv<long>(power(conv<ZZ>(2L), r - m)));
        std::cerr << " lambda = " << lambda;
        omega *= power(zeta, conv<long>(power(conv<ZZ>(2L), r - m - 1)));
        std::cerr << " omega = " << omega << '\n';
        ++count;
    }
    if (count >= r) {
        std::cout << '\n';
        return ZZ_pX();
    }
    ZZ_pX nu = conv<ZZ_pX>(omega);
    std::cerr << " nu = " << nu;

    ZZ n0 = conv<ZZ>(rep(power(omega, (conv<long>(power(field, d)) - 1L) / (conv<long>(field) - 1L)))[0]);
    std::cerr << " n0 = " << n0;
    ZZ n1 = conv<ZZ>(rep(power(delta, (conv<long>(power(field, d)) - 1L) / (conv<long>(field) - 1L)))[0]);
    std::cerr << " n1 = " << n1 << '\n';
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


void fieldsFinder(Vec<ZZ> &fields, Vec<ZZ_pX> &roots, ZZX &poly, Vec<Pair<long, long>> &pairs, ZZ &n, long d) {
    ZZ field(9800L);
    ZZ mul(1L);
    while (mul < n && fields.length() < 3) {
        field = fieldF(field, poly);
        ZZ_pX nu = ShanksTonelli(field, pairs, poly, d);
        if (!IsZero(nu)) {
            fields.append(field);
            roots.append(nu);
            mul *= field;
        }
    }
}




