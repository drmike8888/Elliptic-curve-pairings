/***********************************************************
 *                                                         *
 *    For very simple polynomial show powers modulo a      *
 *    fixed prime polynomial.                              *
 *                                                         *
 **********************************************************/

#include "poly.h"

//#define PRIME 43
#define PRIME 1847

int main(int argc, char *argv[])
{
  POLY r, tst, pow;
  mpz_t n, pk2;
  int ck;

  if(argc < 2)
  {
    printf("Use: ./poly_exp_test <exponent>\n");
    exit(-1);
  }
  mpz_init_set_ui(n, PRIME);
  minit(n);
  mpz_init_set_str(pk2, argv[1], 10);
  gmp_printf("pk2= %Zd\n", pk2);
  poly_init(&r);
  r.deg = 2;
//  mpz_set_ui(r.coef[0], 3);
    mpz_set_ui(r.coef[0], 2);
  mpz_set_ui(r.coef[1], 1);
  mpz_set_ui(r.coef[2], 1);
  poly_mulprep(r);
  poly_printf("r = ", r);
  poly_init(&tst);
  poly_init(&pow);
  tst.deg = 1;
  mrand(tst.coef[1]);
  mrand(tst.coef[0]);
  poly_pow(&pow, tst, pk2);
  poly_printf("taking ", tst);
  printf("to power %s\n", argv[1]);
  poly_printf("gives ", pow);
}
