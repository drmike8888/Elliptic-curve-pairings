/*****************************************************************
 *                                                               *
 *   Compute Tate pairing using tiny field so all points can be  *
 *   examined.                                                   *
 *                                                               *
 ****************************************************************/

#include "poly_eliptic.h"
#include "eliptic.h"
#include "pairing.h"
#include <time.h>

#define M 43
#define XTEND 1815

/*  choose random point with specific order.
    Input grp array, point_list array, length of point_list for this order, 
    and 0 for G_1 or 1 for G_2 point.
    Returns index of point from point_list.
    NOTE: point_list is one dimensional for returned index.
*/

int rndselect(long *grp, int *point_list, int nmpnt, int type)
{
  int j, k, r;

  r = -1;
  while(r < 0)
  {
    k = rand() % nmpnt;
    j = point_list[k];
    if(!type && (grp[2*j] == 1))
      r = j;
    else if(type && (grp[2*j] > 1))
      r = j;
  }
  return j;
}

int main(int argc, char *argv[])
{
  FILE *pnts;
  POLY_CURVE Ex;
  CURVE E;
  POLY irrd, xtnd, t1, t2, t3, t4;
  POLY_POINT *xtndpnt;
  mpz_t prm, x, factors[7], tor;
  mpz_t *xtndordr;
  int i, j, k, m, which, skip;
  POLY_POINT T, P, Q, S, TpQ;
  long *grp;
  int numxtdpnts[7];
  int *point_list, pdex[7];
  struct timespec ts;
  
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

/*  set up factors for tiny curve */

  for(i=0; i<7; i++)
    mpz_init(factors[i]);
  mpz_set_ui(factors[0], 3);
  mpz_set_ui(factors[1], 5);
  mpz_set_ui(factors[2], 11);
  mpz_set_ui(factors[3], 15);
  mpz_set_ui(factors[4], 33);
  mpz_set_ui(factors[5], 55);
  mpz_set_ui(factors[6], 165);
  for(i=0; i<7; i++)
    numxtdpnts[i] = 0;

  mpz_init_set(tor, factors[2]);  // torsion order
  i = 0;
  j = 0;
  mpz_init(x);
  grp = (long*)malloc(sizeof(long)*2*1815);
  
/*  find every point on extension curve and order of that point  */

  poly_curve_init(&Ex);
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
  while(j < XTEND-1)
  {
    poly_elptic_embed(&xtndpnt[j], &xtndpnt[j+1], xtnd, Ex);
    k = poly_get_order(xtndordr[j], xtndpnt[j], Ex, factors, 7);
    grp[2*j] = g1g2(xtndpnt[j]);
    grp[2*j + 1] = mpz_get_ui(xtndordr[j]);
    numxtdpnts[k] += 2;
    mpz_set(xtndordr[j + 1], xtndordr[j]);
    poly_copy(&xtnd, xtndpnt[j].x);
    j++;
    grp[2*j] = grp[2*j - 2];
    grp[2*j + 1] = grp[2*j - 1];
    FF_bump(&xtnd);
    j++;
  }
  for(i=0; i<7; i++)
  {
    gmp_printf("order %Zd has %d points\n", factors[i], numxtdpnts[i]);
    pdex[i] = 0;
  }

/* use nanosecond clock to initialize random generator */

  clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
  srand(ts.tv_nsec);
  
/* create 2D array with index for each order.
   use point list to find random points of given order */

  point_list = (int*)malloc(sizeof(int)*7*1000);
  for(j=0; j<XTEND; j++)
  {
    for(i=0; i<7; i++)
    {
      if(!mpz_cmp_ui(factors[i], grp[2*j + 1]))
      {
	k = i*1000 + pdex[i];
	point_list[k] = j;
	pdex[i]++;
      }
    }
  }

/* point S is reference. Choose point of order 3  */
  
  poly_point_init(&S);
  poly_point_init(&P);
  poly_point_init(&Q);
  poly_point_init(&T);
  poly_init(&t1);
  poly_init(&t2);
  poly_init(&t3);
  poly_init(&t4);

  k = rndselect(grp, &point_list[000], numxtdpnts[0], 1);
  poly_point_copy(&S, xtndpnt[k]);
  poly_point_printf("order 3 S:\n", S);

/*  Do groups of 4 tests with 
    G1 x G1
    G1 x G2
    G2 x G1
    G2 x G2
    using same and different orders */
  
/* Do 4 tests of order 11 for all points  */
  
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);

  printf("===========================================\n");  
  printf("Tate G1 x G1 (order 11)\n\n");
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);

  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_point_init(&TpQ);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  printf("===========================================\n");  
  printf("Tate G1 x G2* (order 11)\n\n");
  
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, Q, P, S, tor, Ex);
  poly_printf("tate(Q, P): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("t(P, Q)* t(Q, P): ", t3);

  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  printf("===========================================\n");  
  printf("Tate G2 x G1 (order 11)\n\n");
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);
  
  printf("===========================================\n");  
  printf("Tate G2 x G2 (order 11)\n\n");
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  printf("===========================================\n");  
  printf("Tate G1 x G1* (order 55)\n\n");

  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  printf("===========================================\n");  
  printf("Tate G1 x G2* (order 55)\n\n");
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  printf("===========================================\n");  
  printf("Tate G2 x G1* (order 55)\n\n");
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);
  
  printf("===========================================\n");  
  printf("Tate G2 x G2* (order 55)\n\n");
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 2);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  k = rndselect(grp, &point_list[1000], numxtdpnts[1], 0);
  poly_point_copy(&S, xtndpnt[k]);
  poly_point_printf("order 5 S:\n", S);
  printf("===========================================\n");  
  printf("Tate G1 x G2* (order 11x33)\n\n");
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[4000], numxtdpnts[4], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[4000], numxtdpnts[4], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);
  
  printf("===========================================\n");  
  printf("Tate G2 x G2* (order 11x33)\n\n");
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[4000], numxtdpnts[4], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[4000], numxtdpnts[4], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  printf("===========================================\n");  
  printf("Tate G2 x G1* (order 33x11)\n\n");
  k = rndselect(grp, &point_list[4000], numxtdpnts[4], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);
  
  printf("===========================================\n");  
  printf("Tate G2 x G2* (order 33x11)\n\n");
  k = rndselect(grp, &point_list[4000], numxtdpnts[4], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  /*  fails catastrophically!*/
  k = rndselect(grp, &point_list[000], numxtdpnts[0], 1);
  poly_point_copy(&S, xtndpnt[k]);
  poly_point_printf("order 3 S:\n", S);
  printf("===========================================\n");  
  printf("Tate G1 x G1* (order 11x55)\n\n");
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);
  /**/  
  printf("===========================================\n");  
  printf("Tate G1 x G2* (order 11x55)\n\n");
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  printf("===========================================\n");  
  printf("Tate G2 x G1* (order 11x55)\n\n");
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);
  
  printf("===========================================\n");  
  printf("Tate G2 x G2* (order 11x55)\n\n");
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  /*  fails catastrophically! */
  printf("===========================================\n");  
  printf("Tate G1 x G1* (order 55x11)\n\n");
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  printf("===========================================\n");  
  printf("Tate G1 x G2* (order 55x11)\n\n");
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 0);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);
  
  printf("===========================================\n");  
  printf("Tate G2 x G1* (order 55x11)\n\n");
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);
  
  printf("===========================================\n");  
  printf("Tate G2 x G2* (order 55x11)\n\n");
  k = rndselect(grp, &point_list[5000], numxtdpnts[5], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);
  /**/
  k = rndselect(grp, &point_list[1000], numxtdpnts[1], 0);
  poly_point_copy(&S, xtndpnt[k]);
  poly_point_printf("order 5 S:\n", S);
  printf("===========================================\n");  
  printf("Tate G1 x G2* (order 11x165)\n\n");
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[6000], numxtdpnts[6], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[6000], numxtdpnts[6], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  printf("===========================================\n");  
  printf("Tate G2 x G2* (order 11x165)\n\n");
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[6000], numxtdpnts[6], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[6000], numxtdpnts[6], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  printf("===========================================\n");  
  printf("Tate G2 x G1* (order 165x11)\n\n");
  k = rndselect(grp, &point_list[6000], numxtdpnts[6], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 0);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);

  printf("===========================================\n");  
  printf("Tate G2 x G2* (order 165x11)\n\n");
  k = rndselect(grp, &point_list[6000], numxtdpnts[6], 1);
  poly_point_copy(&P, xtndpnt[k]);
  poly_point_printf("P:\n", P);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&Q, xtndpnt[k]);
  poly_point_printf("Q:\n", Q);
  k = rndselect(grp, &point_list[2000], numxtdpnts[2], 1);
  poly_point_copy(&T, xtndpnt[k]);
  poly_point_printf("T:\n", T);
  tate(&t1, P, Q, S, tor, Ex);
  poly_printf("tate(P, Q): ", t1);
  tate(&t2, P, T, S, tor, Ex);
  poly_printf("tate(P, T): ", t2);
  poly_mul(&t3, t1, t2);
  poly_printf("(P,Q)*(P,T): ", t3);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  tate(&t4, P, TpQ, S, tor, Ex);
  poly_printf("(P, T+Q): ", t4);
}
