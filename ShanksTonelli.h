//
// Created by Moskovskaya Elizaveta on 15/06/2017.
//

#ifndef GENERAL_NUMBER_FIELD_SIEVE_SHANKSTONELLI_H
#define GENERAL_NUMBER_FIELD_SIEVE_SHANKSTONELLI_H

NTL::ZZ_pX ShanksTonelli(NTL::ZZ& field, NTL::Vec<NTL::Pair<long, long>>& pairs, NTL::ZZX& polynomial, long degree);

#endif //GENERAL_NUMBER_FIELD_SIEVE_SHANKSTONELLI_H
