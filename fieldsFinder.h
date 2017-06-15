//
// Created by Moskovskaya Elizaveta on 14/06/2017.
//

#ifndef GENERAL_NUMBER_FIELD_SIEVE_FIELDSFINDER_H
#define GENERAL_NUMBER_FIELD_SIEVE_FIELDSFINDER_H

void fieldsFinder(NTL::Vec<NTL::ZZ> &fields, NTL::Vec<NTL::ZZ_pX> &roots,
                  NTL::ZZX &poly, NTL::Vec<NTL::Pair<long, long>> &pairs, NTL::ZZ &n, long d);

#endif //GENERAL_NUMBER_FIELD_SIEVE_FIELDSFINDER_H
