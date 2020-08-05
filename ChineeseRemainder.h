//
// Created by Moskovskaya Elizaveta on 16/06/2017.
//

#ifndef GENERAL_NUMBER_FIELD_SIEVE_CHINEESEREMAINDER_H
#define GENERAL_NUMBER_FIELD_SIEVE_CHINEESEREMAINDER_H

NTL::ZZ ChineeseRemainder(const NTL::Vec<NTL::ZZ> &fields, const NTL::Vec<NTL::ZZ_pX> &roots,
                          const NTL::ZZ &m, const NTL::ZZ&n);

#endif //GENERAL_NUMBER_FIELD_SIEVE_CHINEESEREMAINDER_H
