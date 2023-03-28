/***********************************************************
 *                                                         *
 *     Basic include for all GMP polynomial functions.     *
 *                                                         *
 **********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void mod_add(mpz_t a, mpz_t b, mpz_t c, mpz_t n);
void mod_sub(mpz_t a, mpz_t b, mpz_t c, mpz_t n);
void mod_mul(mpz_t a, mpz_t b, mpz_t c, mpz_t n);
void mod_div(mpz_t a, mpz_t b, mpz_t c, mpz_t n);
void mod_neg(mpz_t a, mpz_t b, mpz_t n);
void minit(mpz_t m);
void mget(mpz_t mod);
void mset(mpz_t prm);
void madd(mpz_t a, mpz_t b, mpz_t c);
void msub(mpz_t a, mpz_t b, mpz_t c);
void mmul(mpz_t a, mpz_t b, mpz_t c);
void mdiv(mpz_t a, mpz_t b, mpz_t c);
void minv(mpz_t a, mpz_t b);
void mneg(mpz_t a, mpz_t b);
int msqrt(mpz_t x, mpz_t a);
void mrand(mpz_t rand);
int msqr(mpz_t x);
void mpowi(mpz_t a, mpz_t b, long i);
