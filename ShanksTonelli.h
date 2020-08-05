//
// Created by Moskovskaya Elizaveta on 17/06/2017.
//

#ifndef GENERAL_NUMBER_FIELD_SIEVE_SHANKSTONELLI_H
#define GENERAL_NUMBER_FIELD_SIEVE_SHANKSTONELLI_H

NTL::ZZ_pX ShanksTonelli(const NTL::ZZ& field, const NTL::Vec<NTL::Pair<long, long>>& pairs, const NTL::ZZX& poly, long d);

#endif //GENERAL_NUMBER_FIELD_SIEVE_SHANKSTONELLI_H
