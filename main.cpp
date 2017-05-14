//
// Created by Moskovskaya Elizaveta on 13/05/2017.
//

#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZX.h>
#include <NTL/ZZXFactoring.h>
#include <vector>

#include "degreeInit.h"
#include "polyInit.h"
#include "RFBGen.h"
#include "AFBGen.h"
#include "QCFBGen.h"

using namespace NTL;

int main( ) {
    ZZ n;
    std::cin >> n;
    long degree = degreeInit(n);
    std::cout << degree << '\n';
    ZZ m = TruncToZZ(pow(conv<RR>(n), 1 / conv<RR>(degree)));
    std::cout << m << '\n';
    ZZX poly = polyInit(n, m, degree);
    std::cout << poly << '\n';
    ZZ z = NextPrime(TruncToZZ(exp(pow(0.5 * log(n) * log(conv<RR>(log(n))),RR(1) / 2))));
    std::cout << z << '\n';
    Vec<Pair<ZZ, ZZ>> RFB = RFBGen(m, z);
    std::cout << RFB.length() << '\n';
    Vec<Pair<ZZ, ZZ>> AFB = AFBGen(m, z, poly);
    std::cout << AFB.length() << '\n';
    Vec<Pair<ZZ, ZZ>> QCFB = QCFBGen(m, z, poly);
    std::cout << QCFB.length() << '\n';
    long minPairAmount = RFB.length() + AFB.length() + QCFB.length() + 2L;
    std::cout << minPairAmount;
    return 0;
}