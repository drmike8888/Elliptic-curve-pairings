/*******************************************************************
 *                                                                 *
 *   Understand how quotient group has same number of points as    *
 *   torsion group.                                                *
 *   Especially in cyclic only curve.                              *
 *                                                                 *
 ******************************************************************/

#include "poly_eliptic.h"
#include "eliptic.h"
#include "pairing.h"

#define M 43
#define BASE 55
#define XTEND 1815

int main(int argc, char *argv[])
{
  FILE *dat;
  mpz_t prm, x, *bseordr, *xtndordr, factors[8], tor;
  POLY irrd, xtnd, t1, t2, t3, t4;
  POINT *basepnt;
  POLY_POINT *xtndpnt;
  int i, j, k, quotnt_base, quotnt_xtnd;
  POLY_CURVE Ex;
  CURVE E;
  int numbsepnts[3], numxtdpnts[7];
  
/*  initialize base prime and polynomial  */
  
  mpz_init_set_ui(prm, M);
  minit(prm);
  poly_init(&irrd);
  if(poly_irreducible(&irrd, 2))
    poly_printf("Found irreducible polynomial:\n", irrd);
  else
    printf("no irreducible polynomial found...\n");
  poly_irrd_set(irrd);
  poly_mulprep(irrd);

/*  find every point on base curve and order of that point */

  basepnt = (POINT*)malloc(sizeof(POINT)*BASE);
  for(i=0; i<BASE; i++)
    point_init(&basepnt[i]);
  bseordr = (mpz_t*)malloc(sizeof(mpz_t)*BASE);
  for(i=0; i<BASE; i++)
    mpz_init(bseordr[i]);
  curve_init(&E);
  mpz_set_ui(E.a4, 23);
  mpz_set_ui(E.a6, 42);
  for(i=0; i<8; i++)
    mpz_init(factors[i]);
  mpz_set_ui(factors[0], 5);
  mpz_set_ui(factors[1], 11);
  mpz_set_ui(factors[2], 55);
  mpz_set_ui(factors[3], 3);
  mpz_set_ui(factors[4], 33);
  mpz_set_ui(factors[5], 15);
  mpz_set_ui(factors[6], 165);

  mpz_init_set(tor, factors[1]);  // torsion order
  for(i=0; i<3; i++)
    numbsepnts[i] = 0;
  for(i=0; i<7; i++)
    numxtdpnts[i] = 0;
  i = 0;
  j = 0;
  quotnt_base = 0;
  mpz_init(x);
  
  while(i<M-1)
  {
    printf("j: %d ", j);
    mpz_set_ui(x, i);
    elptic_embed(&basepnt[j], &basepnt[j+1], x, E);
    point_printf("base point ", basepnt[j]);
    k = get_order(bseordr[j], basepnt[j], E,  factors, 3);
    numbsepnts[k] += 2;
    mpz_set(bseordr[j+1], bseordr[j]);
    if(mpz_cmp(tor, bseordr[j]))
       quotnt_base += 2;            // if not order 11, member of E/nE
    while(mpz_cmp_ui(basepnt[j].x, i) >= 0) i++;
    j += 2;
  }

  printf("number of points in base E/nE group: %d\n", quotnt_base);
  for(i=0; i<3; i++)
    gmp_printf("order %Zd has %d points\n", factors[i], numbsepnts[i]);
  printf("\n");
  
/*  find every point on extension curve and order of that point  */

  poly_curve_init(&Ex);
  Ex.a4.deg = 0;
  mpz_set_ui(Ex.a4.coef[0], 23);
  mpz_set_ui(Ex.a6.coef[0], 42);
  xtndpnt = (POLY_POINT*)malloc(sizeof(POLY_POINT)*XTEND);
  for(i=0; i<XTEND; i++)
    poly_point_init(&xtndpnt[i]);
  xtndordr = (mpz_t*)malloc(sizeof(mpz_t)*XTEND);
  for(i=0; i<XTEND; i++)
    mpz_init(xtndordr[i]);
  j = 0;
  poly_init(&xtnd);
  quotnt_xtnd = 0;
  while(j < XTEND-1)
  {
    poly_elptic_embed(&xtndpnt[j], &xtndpnt[j+1], xtnd, Ex);
    k = poly_get_order(xtndordr[j], xtndpnt[j], Ex, factors, 7);
    numxtdpnts[k] += 2;
    mpz_set(xtndordr[j + 1], xtndordr[j]);
    if(mpz_cmp(tor, xtndordr[j]))
      quotnt_xtnd += 2;
    poly_copy(&xtnd, xtndpnt[j].x);
    FF_bump(&xtnd);
    j += 2;
  }
  printf("number of points on extended E/nE group: %d\n", quotnt_xtnd);  
  for(i=0; i<7; i++)
    gmp_printf("order %Zd has %d points\n", factors[i], numxtdpnts[i]);
}
  
