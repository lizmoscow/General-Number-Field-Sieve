//
// Created by Moskovskaya Elizaveta on 15/05/2017.
//

#ifndef GENERAL_NUMBER_FIELD_SIEVE_SIEVING_H
#define GENERAL_NUMBER_FIELD_SIEVE_SIEVING_H

void sieving(NTL::Vec<NTL::Pair<long, long>> &answer, const NTL::Vec<NTL::Pair<NTL::ZZ, NTL::ZZ>>& RFB,
                                        const NTL::Vec<NTL::Pair<NTL::ZZ, NTL::ZZ>>& AFB,
                                        long B, const NTL::ZZ& m, const NTL::ZZX& poly,
                                        long degree, long minPairAmount);

#endif //GENERAL_NUMBER_FIELD_SIEVE_SIEVING_H
