/***********************************************************
 *                                                         *
 *    Verify elliptic routines work as required.           *
 *                                                         *
 **********************************************************/

#include "modulo.h"
#include "eliptic.h"

int main(int argc, char *argv[])
{
  int i, j;
  POINT P, Q, R;
  CURVE E;
  mpz_t r, p;

  if(argc<4)
  {
    printf("Use: ./eliptic_test <prime> <a4> <a6>\n");
    exit(-1);
  }
  if(mpz_init_set_str(p, argv[1], 10) < 0)
  {
    printf("invalid prime input\n");
    exit(-2);
  }
  minit(p);

  curve_init(&E);
  if(mpz_init_set_str(E.a4, argv[2], 10) < 0)
  {
    printf("invalid a4 input\n");
    exit(-3);
  }
  if(mpz_init_set_str(E.a6, argv[3], 10) < 0)
  {
    printf("invalid a6 input\n");
    exit(-4);
  }
  gmp_printf("curve: y^2 = x^3 + %Zd*x + %Zd\n", E.a4, E.a6);

  mpz_init(r);
  mrand(r);
  point_init(&P);
  point_init(&Q);
  point_init(&R);
  elptic_embed(&P, &Q, r, E);
  point_printf("P1: ", P);
  point_printf("P2: ", Q);

  elptic_sum(&R, Q, Q, E);
  point_printf("2 P2: ", R);
  mpz_set_ui(r, 13);
  elptic_mul(&R, P, r, E);
  point_printf("13*P1: ", R);
}
  
  
