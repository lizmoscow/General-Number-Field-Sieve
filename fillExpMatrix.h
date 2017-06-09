//
// Created by Moskovskaya Elizaveta on 20/05/2017.
//

#ifndef GENERAL_NUMBER_FIELD_SIEVE_FILLEXPMATRIX_H
#define GENERAL_NUMBER_FIELD_SIEVE_FILLEXPMATRIX_H

void fillExpMatrix(NTL::mat_GF2&, const NTL::Vec<NTL::Pair<long, long>>&,
                   const NTL::Vec<NTL::Pair<NTL::ZZ, NTL::ZZ>>&, const NTL::Vec<NTL::Pair<NTL::ZZ, NTL::ZZ>>&,
                   const NTL::Vec<NTL::Pair<NTL::ZZ, NTL::ZZ>>&, const NTL::ZZX&, const NTL::ZZ&, long);

#endif //GENERAL_NUMBER_FIELD_SIEVE_FILLEXPMATRIX_H
