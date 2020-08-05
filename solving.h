//
// Created by Moskovskaya Elizaveta on 19/06/2017.
//

#ifndef GENERAL_NUMBER_FIELD_SIEVE_SOLVING_H
#define GENERAL_NUMBER_FIELD_SIEVE_SOLVING_H

NTL::Pair<NTL::ZZ, NTL::ZZ> solving(const NTL::ZZ &n, const long &degree, const NTL::ZZ &m, const NTL::ZZX &poly,
                                    const NTL::Vec<NTL::Pair<long, long>> &answer, NTL::mat_GF2 &matrix);

#endif //GENERAL_NUMBER_FIELD_SIEVE_SOLVING_H
