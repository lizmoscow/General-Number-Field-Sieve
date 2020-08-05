//
// Created by Moskovskaya Elizaveta on 19/06/2017.
//

#ifndef GENERAL_NUMBER_FIELD_SIEVE_INITIALIZATION_H
#define GENERAL_NUMBER_FIELD_SIEVE_INITIALIZATION_H

void initialization(NTL::ZZ n, long &degree, NTL::ZZ &m, NTL::ZZX &poly,
                    NTL::Vec<NTL::Pair<long, long>> &answer, NTL::mat_GF2 & matrix);

#endif //GENERAL_NUMBER_FIELD_SIEVE_INITIALIZATION_H
