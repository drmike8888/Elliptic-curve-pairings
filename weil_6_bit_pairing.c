/*********************************************************************
 *                                                                   *
 *    Tiny example to show how field extension works.  Curve was     *
 *    found with brute force search and luck.                        *
 *    base prime is 43, number of points on base curve is 11 and     *
 *    Weil pairing works with points of order 11 not on G1.          *
 *                                                                   *
 ********************************************************************/

#include "poly_eliptic.h"
#include "eliptic.h"
#include "pairing.h"

#define M 43

int main(int argc, char *argv[])
{
  FILE *pnts;
  POLY_CURVE Ex;
  CURVE E;
  POLY_POINT Px1, Px2, Qx, Tx;
  POINT P1, P2;
  mpz_t prm, x, ordr, factors[8], tor;
  POLY irrd, xtnd, t1, t2, t3, t4;
  int i, j, k, m, which, skip;
  POLY_POINT table[1815], T, P, Q, S, TpQ;
  long *grp;
  
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

/*  find every point on base curve  */

  pnts = fopen("all_points_prm_43.dat", "w");
  curve_init(&E);
  mpz_set_ui(E.a4, 23);
  mpz_set_ui(E.a6, 42);
  poly_curve_init(&Ex);
  point_init(&P1);
  point_init(&P2);
  mpz_inits(x, ordr, NULL);
  for(i=0; i<8; i++)
    mpz_init(factors[i]);
  mpz_set_ui(factors[0], 5);
  mpz_set_ui(factors[1], 11);
  mpz_set_ui(factors[2], 55);
  mpz_set_ui(factors[3], 3);
  mpz_set_ui(factors[4], 33);
  mpz_set_ui(factors[5], 15);
  mpz_set_ui(factors[6], 165);
  mpz_set_ui(factors[7], 1815);
  
  i = 1;
  j = 0;
  while(i<M-1)
  {
    mpz_set_ui(x, i);
    elptic_embed(&P1, &P2, x, E);
    get_order(ordr, P1, E, factors, 3);
    gmp_fprintf(pnts, "%2d: (%Zd, %Zd) order: %Zd\n", j, P1.x, P1.y, ordr);
    j++;
    gmp_fprintf(pnts, "%2d: (%Zd, %Zd)\n", j, P2.x, P2.y);
    j++;
    while(mpz_cmp_ui(P1.x, i) >= 0) i++;
  }
  fprintf(pnts, "\n");

/*  find every point on extension curve  */

  Ex.a4.deg = 0;
  mpz_set_ui(Ex.a4.coef[0], 23);
  mpz_set_ui(Ex.a6.coef[0], 42);
  poly_point_init(&Px1);
  poly_point_init(&Px2);
  poly_init(&xtnd);
  j = 0;
  grp = (long*)malloc(sizeof(long)*2*1815);
  while(j < 1814)
  {
    poly_elptic_embed(&Px1, &Px2, xtnd, Ex);
    poly_get_order(ordr, Px1, Ex, factors, 8);
    gmp_fprintf(pnts, "%2d: x = ", j);
    if(Px1.x.deg)
      gmp_fprintf(pnts, "%Zd*x + ", Px1.x.coef[1]);
    gmp_fprintf(pnts, "%Zd  y = ", Px1.x.coef[0]);
    if(Px1.y.deg)
      gmp_fprintf(pnts, "%Zd*x + ", Px1.y.coef[1]);
    gmp_fprintf(pnts, "%Zd  order: %Zd | %ld\n", Px1.y.coef[0], ordr, g1g2(Px1));
    poly_point_init(&table[j]);
    poly_point_copy(&table[j], Px1);
    grp[2*j] = g1g2(Px1);
    grp[2*j + 1] = mpz_get_ui(ordr);
    j++;
    gmp_fprintf(pnts, "%2d: x = ", j);
    if(Px2.x.deg)
      gmp_fprintf(pnts, "%Zd*x + ", Px2.x.coef[1]);
    gmp_fprintf(pnts, "%Zd  y = ", Px2.x.coef[0]);
    if(Px2.y.deg)
      gmp_fprintf(pnts, "%Zd*x + ", Px2.y.coef[1]);
    gmp_fprintf(pnts, "%Zd | %ld\n", Px2.y.coef[0], g1g2(Px2));
    poly_point_init(&table[j]);
    poly_point_copy(&table[j], Px2);
    grp[2*j + 1] = g1g2(Px2);
    grp[2*j] = mpz_get_ui(ordr);
    j++;
    poly_copy(&xtnd, Px1.x);
    FF_bump(&xtnd);
  }
  fclose(pnts);

/*  check cardinality routine  */

  mpz_set_si(x, -11);
  cardinality(ordr, x, 2);
  gmp_printf("cardinality of extension: %Zd\n", ordr);
  k = mpz_get_ui(ordr);
  printf("k=%d\n", k);
  
/*  pick 3 independant points and 1 not in same order
    not general - fixed for this problem.

    This is G1 x G1 test  */
  printf("G1 x G1 test--------------------\n");
  which = 0;
  for(i=1; i<k; i++)
  {
    if(!which)
    {
      if((grp[2*i] == 11) && (grp[2*i + 1] == 1))
      {
	printf("%d %d\n", which, i);
	poly_point_copy(&P, table[i]);
	which++;
	i++;
      }
      continue;
    }
    if(which == 1)
    {
      if((grp[2*i] == 11) && (grp[2*i + 1] == 1))
      {
	printf("%d %d\n", which, i);
	poly_point_copy(&Q, table[i - 1]);
	which++;
	i++;
      }
      continue;
    }
    if(which == 2)
    {
      if((grp[2*i] == 11) && (grp[2*i + 1] == 1))
      {
	printf("%d %d\n", which, i);
	poly_point_copy(&T, table[i]);
	which++;
	i++;
      }
      continue;
    }
    if(which == 3)
    {
      i=1;
      while(i<k)
      {
	if((grp[2*i] == 55) && (grp[2*i + 1] == 4))
        {
	  printf("%d %d\n", which, i);
	  poly_point_copy(&S, table[i]);
	  which++;
	  break;
	}
	i++;
      }
      continue;
    }
    if(which == 4)
      break;
  }
  poly_point_printf("P:\n", P);
  poly_point_printf("Q:\n", Q);
  poly_point_printf("T:\n", T);
  poly_point_printf("S:\n", S);

/* compute (P,Q)  */

  mpz_init_set_ui(tor, 11);
  poly_init(&t1);
  weil(&t1, P, Q, S, tor, Ex);
  poly_printf("G1 x G1 weil (P, Q): ", t1);

/* compute (P,T)  */

  poly_init(&t2);
  weil(&t2, P, T, S, tor, Ex);
  poly_printf("G1 x G1 weil (P, T): ", t2);

/* compute (P,T+Q)  */

  poly_point_init(&TpQ);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  poly_point_printf("T + Q: \n", TpQ);
  weil(&t3, P, TpQ, S, tor, Ex);
  poly_printf("G1 x G1 weil (P, T+Q): ", t3);

  poly_init(&t4);
  poly_mul(&t4, t1, t2);
  poly_printf("G1 x G1 (P, T)*(P, Q): ", t4);
  printf("--------------------------------------------\n\n");
  
/*    This is G1 x G2 test  */

  printf(" G1 x G2 test-----------------------------\n");
  which = 0;
  for(i=1; i<k; i++)
  {
    if(!which)
    {
      if((grp[2*i] == 11) && (grp[2*i + 1] == 1))
      {
	printf("%d %d\n", which, i);
	poly_point_copy(&P, table[i]);
	which++;
	i++;
      }
      continue;
    }
    if(which == 1)
    {
      if((grp[2*i] == 11) && (grp[2*i + 1] == 4))
      {
	printf("%d %d\n", which, i);
	poly_point_copy(&Q, table[i - 1]);
	which++;
	i+=40;
      }
      continue;
    }
    if(which == 2)
    {
      if((grp[2*i] == 11) && (grp[2*i + 1] == 4))
      {
	printf("%d %d\n", which, i);
	poly_point_copy(&T, table[i]);
	which++;
	i++;
      }
      continue;
    }
    if(which == 3)
    {
      i=1;
      while(i<k)
      {
	if((grp[2*i] == 55) && (grp[2*i + 1] == 4))
        {
	  printf("%d %d\n", which, i);
	  poly_point_copy(&S, table[i]);
	  which++;
	  break;
	}
	i++;
      }
      continue;
    }
    if(which == 4)
      break;
  }
  if(which != 4) exit(0);
  poly_point_printf("P:\n", P);
  poly_point_printf("Q:\n", Q);
  poly_point_printf("T:\n", T);
  poly_point_printf("S:\n", S);

/* compute (P,Q)  */

  mpz_init_set_ui(tor, 11);
  poly_init(&t1);
  weil(&t1, P, Q, S, tor, Ex);
  poly_printf("G1 x G2 weil (P, Q): ", t1);

/* compute (P,T)  */

  poly_init(&t2);
  weil(&t2, P, T, S, tor, Ex);
  poly_printf("G1 x G2 weil (P, T): ", t2);

/* compute (P,T+Q)  */

  poly_point_init(&TpQ);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  poly_point_printf("T + Q: \n", TpQ);
  weil(&t3, P, TpQ, S, tor, Ex);
  poly_printf("G1 x G2 weil (P, T+Q): ", t3);

  poly_init(&t4);
  poly_mul(&t4, t1, t2);
  poly_printf("G1 x G2 weil (P, T)*(P, Q): ", t4);

/* take pairing to torsion power */

  poly_pow(&t2, t4, factors[1]);
  poly_printf("Result to torsion power: ", t2);
  printf("-----------------------------------------------\n\n");

/*    This is G2 x G2 test  */

  printf("G2 x G2 test ------------------------\n");
  which = 0;
  skip = 0;
  for(i=0; i<k; i++)
  {
    if(!which)
    {
      if((grp[2*i] == 11) && (grp[2*i + 1] == 4))
      {
	printf("%d %d\n", which, i);
	poly_point_copy(&P, table[i]);
	if(!skip)
	  skip++;
	else
	{
	  which++;
	  i++;
	}
      }
      continue;
    }
    if(which == 1)
    {
      if((grp[2*i] == 11) && (grp[2*i + 1] == 4))
      {
	printf("%d %d\n", which, i);
	poly_point_copy(&Q, table[i - 1]);
	which++;
	i++;
      }
      continue;
    }
    if(which == 2)
    {
      if((grp[2*i] == 11) && (grp[2*i + 1] == 4))
      {
	printf("%d %d\n", which, i);
	poly_point_copy(&T, table[i]);
	which++;
	i++;
      }
      continue;
    }
    if(which == 3)
    {
      i=1;
      while(i<k)
      {
	if((grp[2*i] == 55) && (grp[2*i + 1] == 4))
        {
	  printf("%d %d\n", which, i);
	  poly_point_copy(&S, table[i]);
	  which++;
	  break;
	}
	i++;
      }
      continue;
    }
    if(which == 4)
      break;
  }
  if(which != 4) exit(0);
  poly_point_printf("P:\n", P);
  poly_point_printf("Q:\n", Q);
  poly_point_printf("T:\n", T);
  poly_point_printf("S:\n", S);

/* compute (P,Q)  */

  mpz_init_set_ui(tor, 11);
  poly_init(&t1);
  weil(&t1, P, Q, S, tor, Ex);
  poly_printf("G2 x G2 weil (P, Q): ", t1);

/* compute (P,T)  */

  poly_init(&t2);
  weil(&t2, P, T, S, tor, Ex);
  poly_printf("G2 x G2 weil (P, T): ", t2);

/* compute (P,T+Q)  */

  poly_point_init(&TpQ);
  poly_elptic_sum(&TpQ, T, Q, Ex);
  poly_point_printf("T + Q: \n", TpQ);
  weil(&t3, P, TpQ, S, tor, Ex);
  poly_printf("G2 x G2 weil (P, T+Q): ", t3);

  poly_init(&t4);
  poly_mul(&t4, t1, t2);
  poly_printf("G2 x G2 weil (P, T)*(P, Q): ", t4);

/* take pairing to torsion power */

  poly_pow(&t2, t4, factors[1]);
  poly_printf("Result to torsion power: ", t2);
  printf("\n\n");
/*
  pnts = fopen("double.check", "w");
  for(i=0; i<1814; i++)
  {
    poly_get_order(ordr, table[i], Ex, factors, 8);
    gmp_fprintf(pnts, "%2d: x = ", i);
    if(table[i].x.deg)
      gmp_fprintf(pnts, "%Zd*x + ", table[i].x.coef[1]);
    gmp_fprintf(pnts, "%Zd  y = ", table[i].x.coef[0]);
    if(table[i].y.deg)
      gmp_fprintf(pnts, "%Zd*x + ", table[i].y.coef[1]);
    gmp_fprintf(pnts, "%Zd  order: %Zd | %ld\n", table[i].y.coef[0], ordr, g1g2(table[i]));
  }
  fclose(pnts);
*/
}
