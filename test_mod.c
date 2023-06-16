/*************************************************
 *                                               *
 *   Make sure modular base functions work.      *
 *                                               *
 ************************************************/

#include "modulo.h"
#include "poly.h"

int main(int argc, char *argv[])
{
  mpz_t a, b, c, n;
  POLY A, B, C, D;
  int i, j;
  POLY q, r, x, xp, xp2, xp3;
  
  mpz_inits(a, b, c, n, NULL);
  mpz_set_str(n, "54881133298817149", 10);
  mpz_set_str(b, "54881122831868268", 10);
  mpz_set_str(c, "10022120677842", 10);

  mod_add(a, b, c, n);
  gmp_printf("b+c = %Zd\n", a);
  mod_sub(a, b, c, n);
  gmp_printf("b-c = %Zd\n", a);
  mod_mul(a, b, c, n);
  gmp_printf("b*c = %Zd\n", a);
  mod_div(a, b, c, n);
  gmp_printf("b/c = %Zd\n", a);

  minit(n);
  poly_init(&A);
  poly_init(&B);
  poly_init(&C);
  A.deg = 1;
  mpz_set_ui(A.coef[0], 1);
  mpz_set_ui(A.coef[1], 1);
  poly_print(A);
  B.deg = 2;
  mpz_set(B.coef[0], a);
  mpz_set(B.coef[1], b);
  mpz_set(B.coef[2], c);
  poly_print(B);
/*
  printf("\n");
  poly_add(&C, A, B);
  poly_print(C);
  printf("\n");
  poly_normal(C.coef, C.deg+1);
  poly_print(C);
  printf("\n");
*/
  mpz_set_ui(n, 7);
  minit(n);
  C.deg = 4;
  mpz_set_ui(C.coef[4], 1);
  mpz_set_ui(C.coef[3], 2);
  mpz_set_ui(C.coef[2], 1);
  mpz_set_ui(C.coef[1], 3);
  mpz_set_ui(C.coef[0], 5);
  poly_mulprep(C);
  printf("\n");
  A.deg = 3;
  mpz_set_ui(A.coef[1], 3);
  mpz_set_ui(A.coef[2], 2);
  mpz_set_ui(A.coef[3], 1);
  mpz_set_ui(B.coef[0], 3);
  mpz_set_ui(B.coef[1], 5);
  mpz_set_ui(B.coef[2], 1);
  poly_print(A);
  printf("\n");
  poly_print(B);
  printf("\n");
  poly_init(&D);
  poly_mul(&D, A, B);
  poly_print(D);

  poly_init(&q);
  poly_init(&r);
  printf("D:\n");
  poly_print(D);
  printf("A:\n");
  poly_print(A);
  poly_euclid(&q, &r, A, D);
  printf("q:\n");
  poly_print(q);
  printf("r:\n");
  poly_print(r);

  poly_gcd(&C, A, D);
  printf("gcd(A, D):\n");
  poly_print(C);

  B.deg = 4;
  mpz_set_ui(B.coef[4], 1);
  mpz_set_ui(B.coef[3], 0);
  mpz_set_ui(B.coef[2], 0);
  mpz_set_ui(B.coef[1], 0);
  mpz_set_ui(B.coef[0], 6);
  A.deg = 2;
  mpz_set_ui(A.coef[2], 1);
  mpz_set_ui(A.coef[1], 0);
  mpz_set_ui(A.coef[0], 1);
  printf("A:\n");
  poly_print(A);
  printf("B:\n");
  poly_print(B);
  poly_gcd(&C, A, B);
  printf("gcd(A, B):\n");
  poly_print(C);
  printf("\n");
  
  poly_init(&x);
  x.deg = 1;
  mpz_set_ui(x.coef[1], 1);
  poly_init(&xp);
  poly_xp(&xp, x);
  printf("x^p:\n");
  poly_print(xp);
  printf("\n");
  poly_init(&xp2);
  poly_xp(&xp2, xp);
  poly_printf("x^p^2:\n", xp2);

  poly_init(&xp3);
  poly_xp(&xp3, xp2);
  poly_printf("x^p^3:\n", xp3);
}

