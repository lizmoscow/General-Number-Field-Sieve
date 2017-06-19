//
// Created by Moskovskaya Elizaveta on 18/06/2017.
//

#ifndef GENERAL_NUMBER_FIELD_SIEVE_LANCZOS_H
#define GENERAL_NUMBER_FIELD_SIEVE_LANCZOS_H

void lanczos(const NTL::mat_GF2 &matrix, NTL::vec_GF2 &x);
bool block_lanczos(const NTL::mat_GF2 & A, NTL::vec_GF2 & x, const NTL::vec_GF2 & y);

#endif //GENERAL_NUMBER_FIELD_SIEVE_LANCZOS_H
