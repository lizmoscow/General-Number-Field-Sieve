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
#include <fstream>

#include "GNFS.h"

#define TEST

using namespace NTL;

int main(int argc, char* argv[]) {
    try {

#ifdef TEST

        std::string filename = "/Users/liza_moskovskaya/General-Number-Field-Sieve/test.csv";
        std::ofstream myfile;
        myfile.open(filename);
        myfile << "exp,time,\n";
        ZZ x, y, factor;
        for (long i = 1L; i <= 10; ++i) {
            if (i == 1L) {
                x = 50;
            }
            else {
                x = NextPrime(power(ZZ(10L), i));
            }
            y = NextPrime(10 * power(ZZ(10L), i));
            auto start = std::chrono::system_clock::now();
            for (long j = 0; j < 5; ++j) {
                x = NextPrime(x + 1);
                y = NextPrime(y + 1);
                factor = GNFS(x * y);
                std::cout << "**********************************************************************************\n" <<
                          factor << " * " << x * y / factor << " = " << x * y << "\n" <<
                          "**********************************************************************************\n";
            }
            auto diff = std::chrono::system_clock::now() - start;
            auto sec = std::chrono::duration_cast<std::chrono::seconds>(diff);
            myfile << i << "," << sec.count() << ",\n";
            std::cout << i << " : " << sec.count() << "\n\n";
        }
        myfile.close();

#endif

#ifndef TEST

        if (argc > 1) {
            ZZ n(atoi(argv[1])), factor;
            auto start = std::chrono::system_clock::now();

            factor = GNFS(n);

            std::cout << "**********************************************************************************\n" <<
                      factor << " * " << n / factor << " = " << n << "\n" <<
                      "**********************************************************************************\n";

            auto diff = std::chrono::system_clock::now() - start;
            auto sec = std::chrono::duration_cast<std::chrono::seconds>(diff);
            std::cout << "\nduration\n: " << sec.count();
        }

#endif
    }
    catch(const char *e) {
        std::cerr << e;
    }

    return 0;
}