//
// Created by Moskovskaya Elizaveta on 17/06/2017.
//

#ifndef GENERAL_NUMBER_FIELD_SIEVE_FILLEXPSPMATRIX_H
#define GENERAL_NUMBER_FIELD_SIEVE_FILLEXPSPMATRIX_H
#include <linbox/matrix/sparse-matrix.h>

void fillExpSpMatrix(LinBox::SparseMatrix<Givaro::Modular<short>>& matrix, const NTL::Vec<NTL::Pair<long, long>> &answer,
                     const NTL::Vec<NTL::Pair<NTL::ZZ, NTL::ZZ>>& RFB, const NTL::Vec<NTL::Pair<NTL::ZZ, NTL::ZZ>>& AFB,
                     const NTL::Vec<NTL::Pair<NTL::ZZ, NTL::ZZ>>& QCFB, const NTL::ZZX& poly, const NTL::ZZ& m, long d);

#endif //GENERAL_NUMBER_FIELD_SIEVE_FILLEXPSPMATRIX_H
