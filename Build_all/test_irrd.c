/***********************************************************
 *                                                         *
 *   Test polynomial routines to find irreducible          *
 *   polynomials. Expand to test poly-eliptic routines.    *
 *                                                         *
 **********************************************************/

#include "poly_eliptic.h"
#include "eliptic.h"

//#define PRIME 1223
#define PRIME 43

int main(int argc, char *argv[])
{
  POLY r, tst, pow, tst2;
  mpz_t n, pk2;
  POLY_CURVE Ex;
  CURVE E;
  int ck;
  POLY_POINT P1, P2;
  
  mpz_init_set_str(n, "394574950729", 10);
//  mpz_init_set_ui(n, PRIME);
  minit(n);
  poly_init(&r);
  if(poly_irreducible(&r, 2))
    poly_printf("Found!\n", r);
  else
    printf("no irreducible polynomial found...\n");
  poly_irrd_set(r);
  poly_mulprep(r);
  poly_init(&tst);
  poly_init(&pow);
  poly_init(&tst2);
  tst.deg = 1;
  mpz_set_ui(tst.coef[1], 2);
  mpz_set_ui(tst.coef[0], 0);
  poly_printf("starting x:\n", tst);
  curve_init(&E);
  mpz_set_ui(E.a4, 23);
  mpz_set_ui(E.a6, 42);
  poly_curve_init(&Ex);
  Ex.a4.deg = 0;
  mpz_set_ui(Ex.a4.coef[0], 23);
  mpz_set_ui(Ex.a6.coef[0], 42);
  poly_curve_printf("Extended curve coefficients:\n", Ex);
  poly_point_init(&P1);
  poly_point_init(&P2);
  poly_elptic_embed(&P1, &P2, tst, Ex);

  poly_point_printf("P1:\n", P1);
  poly_point_printf("P2:\n", P2);
  poly_fofx(&tst2, tst, Ex);
  poly_printf("fofx:\n", tst2);
  ck = poly_sqr(tst2);
  printf("residue: %d\n", ck);

  if(ck)
    poly_sqrt(&pow, tst2);
  poly_printf("sqrt = \n", pow);

/* test p^k-1 / 2  */

  mpz_init(pk2);
  mpz_pow_ui(pk2, n, 2);
  mpz_sub_ui(pk2, pk2, 1);
  while(!mpz_tstbit(pk2, 0))
    mpz_divexact_ui(pk2, pk2, 2);
  gmp_printf("q = %Zd\n", pk2);
  poly_pow(&pow, tst2, pk2);
  poly_printf("fofx^(q)\n", pow);
}
