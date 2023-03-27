/**************************************************************
 *                                                            *
 *     polynomials in GMP for elliptic curve extensions.      *
 *                                                            *
 *************************************************************/

#include "modulo.h"

#ifndef POLY_H
#define POLY_H

#define MAXDEGREE 32

typedef struct
{
  unsigned long deg;
  mpz_t coef[MAXDEGREE];
}POLY;

void poly_init(POLY *p);
void poly_clear(POLY *p);
void poly_print(POLY a);
void poly_printf(char *string, POLY a);
void poly_copy(POLY *a, POLY b);
void poly_irrd_set(POLY i);
void poly_irrd_get(POLY *i);
void poly_q_get(mpz_t pk);
void poly_add(POLY *c, POLY a, POLY b);
void poly_sub(POLY *c, POLY a, POLY b);
void poly_normal(POLY *a);
void poly_mulprep(POLY f);
void poly_mul(POLY *rslt, POLY a, POLY b);
void poly_euclid(POLY *q, POLY *r, POLY a, POLY b);
void poly_gcd(POLY *d, POLY a, POLY b);
void poly_sqm(POLY *x2, POLY x, POLY a, int flag);
void poly_xp(POLY *xp, POLY x);
int poly_cmp(POLY a, POLY b);
int poly_irreducible(POLY *f, long n);
void poly_pow(POLY *h, POLY g, mpz_t k);
void gpow_p2(POLY *h, POLY g);
void poly_pseudo_div(POLY *Q, POLY *R, POLY A, POLY B);
void poly_cont(mpz_t cont, POLY A);
void poly_resltnt(mpz_t rsltnt, POLY A, POLY B);
int poly_sqr(POLY x);
void poly_sqrt(POLY *sqt, POLY a);
void poly_div(POLY *a, POLY b, POLY c);
void poly_rand(POLY *rnd);

#endif
