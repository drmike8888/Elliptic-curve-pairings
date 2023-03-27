/******************************************************************
 *                                                                *
 *   Subroutines to compute pairings and associated utilities.    *
 *                                                                *
 *****************************************************************/

#include "modulo.h"
#include "eliptic.h"
#include "poly_eliptic.h"
#include "pairing.h"

/* determine if point is in G1 or G2 group  
   returns 1 for x in group 1, y in group 1
           2 for x in group 2, y in group 1
           3 for x in group 1, y in group 2
           4 for x in group 2, y in group 2  */

long g1g2(POLY_POINT P)
{
  if((P.x.deg > 0) && (P.y.deg > 0))
    return 4;
  if(P.x.deg > 0)
    return 2;
  if(P.y.deg > 0)
    return 3;
  return 1;
}

/*  compute cardinality of extension field.
    Enter with trace of Frobenius and
    extension degree.
*/

void cardinality(mpz_t crd, mpz_t t, long k)
{
  mpz_t *v, t1, t2, p, pk;
  int i;
  
  v = (mpz_t *)malloc(sizeof(mpz_t)*(k+1));
  for(i=0; i<=k; i++)
    mpz_init(v[i]);
  mpz_set_ui(v[0], 2);
  mpz_set(v[1], t);
  mpz_inits(t1, t2, NULL);
  mget(p);
  for(i=2; i<=k; i++)
  {
    mpz_mul(t1, t, v[i - 1]);
    mpz_mul(t2, p, v[i - 2]);
    mpz_sub(v[i], t1, t2);
  }
  poly_q_get(pk);
  mpz_set(crd, pk);
  mpz_add_ui(crd, crd, 1);
  mpz_sub(crd, crd, v[k]);
  for(i=0; i<=k; i++)
    mpz_clear(v[i]);
  mpz_clears(t1, t2, p, pk, NULL);
}

/*  determine slope of line between two points
    and return miller's h value for 3rd point.
    inputs are curve E, points P, Q and R. Also
    need field prime z. AEC pg 394
*/

void hpq(POLY *h, POLY_POINT P, POLY_POINT Q, POLY_POINT R, POLY_CURVE E)
{
  POLY t, lmbda, b, tx, t1;

  if((poly_test_point(P)) || poly_test_point(Q))
  {
    h->deg = 0;    // if P or Q is infinity, return 1
    mpz_set_ui(h->coef[0], 1);
    return;
  }
  poly_init(&t);
  poly_init(&b);
  poly_add(&b, P.y, Q.y);  // divide by zero?
  if(!b.deg && !mpz_cmp_ui(b.coef[0], 0))
  {
    poly_sub(&b, P.x, Q.x);      // x_p - x_q
    if(!b.deg && !mpz_cmp_ui(b.coef[0], 0)) // Really P == -Q?
    {
      poly_sub(h, R.x, P.x);   // slope=infinity, return x_r - x_p
      poly_clear(&t);
      poly_clear(&b);
      return;
    }
    poly_sub(&t, P.y, Q.y);
  }
  else
  {
/*  compute lambda (slope between P and Q) using secure form  */

    poly_init(&t1);
    poly_init(&tx);
    poly_mul(&t, P.x, P.x);
    poly_mul(&t1, P.x, Q.x);
    poly_mul(&tx, Q.x, Q.x);
    poly_add(&t, t, t1);
    poly_add(&t, t, tx);
    poly_add(&t, t, E.a4);
    poly_add(&b, Q.y, P.y);
  }
  poly_init(&lmbda);
  poly_div(&lmbda, t, b);

/* finally compute h  */
  
  poly_sub(&t, R.y, P.y);
  poly_sub(&tx, R.x, P.x);
  poly_mul(&tx, tx, lmbda);
  poly_sub(&t, t, tx);
  poly_mul(&tx, lmbda, lmbda);
  poly_sub(&b, R.x, tx);
  poly_add(&b, b, P.x);
  poly_add(&b, b, Q.x);
  poly_div(h, t, b);
  poly_clear(&t);
  poly_clear(&b);
  poly_clear(&tx);
  poly_clear(&t1);
  poly_clear(&lmbda);
}

/* compute Miller's algorithm for Weil Pairing 
   input curve E, points of same order P & Q
   and point of different order R. Also need
   torsion order m.  
*/

void miller(POLY *f, POLY_POINT P, POLY_POINT R, mpz_t m, POLY_CURVE E)
{
  POLY_POINT T;
  POLY h;
  long mask;

  mask = mpz_sizeinbase(m, 2) - 2;
  poly_point_init(&T);
  poly_point_copy(&T, P);
  f->deg = 0;
  mpz_set_ui(f->coef[0], 1);
  poly_init(&h);
  while(mask>=0)
  {
    hpq(&h, T, T, R, E);
    poly_mul(f, *f, *f);
    poly_mul(f, *f, h);
    poly_elptic_sum(&T, T, T, E);
    if(mpz_tstbit(m, mask))
    {
      hpq(&h, T, P, R, E);
      poly_mul(f, *f, h);
      poly_elptic_sum(&T, T, P, E);
    }
    mask--;
  }
  poly_point_clear(&T);
  poly_clear(&h);
}

/* compute Weil pairing  based on Silverman's book page 396.
   Need curve E, points P, Q of same order and S of different
   order, m is order of P, Q.  
   Returns w which is mth root of 1.
*/

void weil(POLY *w, POLY_POINT P, POLY_POINT Q, POLY_POINT S, mpz_t m, POLY_CURVE E)
{
  POLY_POINT QpS, mS, PmS;
  POLY t1, t2, t3, t4, w1, w2;

  poly_point_init(&QpS);
  poly_elptic_sum(&QpS, Q, S, E);
  poly_point_init(&mS);
  poly_copy(&mS.x, S.x);
  poly_point_init(&PmS);
  poly_sub(&mS.y, PmS.y, S.y);   // -y = 0 - y
  poly_elptic_sum(&PmS, P, mS, E);
  poly_init(&t1);
  miller(&t1, P, QpS, m, E);
  poly_init(&t2);
  miller(&t2, P, S, m, E);
  poly_init(&t3);
  miller(&t3, Q, PmS, m, E);
  poly_init(&t4);
  miller(&t4, Q, mS, m, E);
  poly_init(&w1);
  poly_div(&w1, t1, t2);
  poly_init(&w2);
  poly_div(&w2, t3, t4);
  poly_div(w, w1, w2);

  poly_clear(&w1);
  poly_clear(&w2);
  poly_clear(&t1);
  poly_clear(&t2);
  poly_clear(&t3);
  poly_clear(&t4);
  poly_point_clear(&QpS);
  poly_point_clear(&PmS);
  poly_point_clear(&mS);
}

/* compute Tate pairing  based on Silverman's book page 396.
   Need curve E, points P, Q with Q in torsion group and S of different
   order, m is order of Q.
*/

void tate(POLY *t, POLY_POINT P, POLY_POINT Q, POLY_POINT S, mpz_t m, POLY_CURVE E)
{
  POLY_POINT QpS;
  POLY t1, t2;
  mpz_t pw;
  
  poly_point_init(&QpS);
  poly_init(&t1);
  poly_init(&t2);
  poly_elptic_sum(&QpS, Q, S, E);
  miller(&t1, P, QpS, m, E);
  miller(&t2, P, S, m, E);
  poly_div(t, t1, t2);
  poly_q_get(pw);
  mpz_sub_ui(pw, pw, 1);
  mpz_divexact(pw, pw, m);
  poly_pow(t, *t, pw);
  poly_point_clear(&QpS);
  poly_clear(&t1);
  poly_clear(&t2);
  mpz_clear(pw);
}

/*  find order of a point.  Requires list of factors,
    curve parameters and point to test.  
    Return index found in factors list.   */

int get_order(mpz_t order, POINT P, CURVE E, mpz_t *factors, int n)
{
  int i;
  POINT R;

  point_init(&R);
  for(i=0; i<n; i++)
  {
    elptic_mul(&R, P, factors[i], E);
    if(test_point(R))
      break;
  }
  if(i<n)
    mpz_set(order, factors[i]);
  else
  {
    printf("missing order in base!!\n");
    exit(-3);
  }
  point_clear(&R);
  return i;
}

/* same thing but for polynomials  */

int poly_get_order(mpz_t order, POLY_POINT P, POLY_CURVE E, mpz_t *factors, int n)
{
  int i;
  POLY_POINT R;

  poly_point_init(&R);
  for(i=0; i<n; i++)
  {
    poly_elptic_mul(&R, P, factors[i], E);
    if(poly_test_point(R))
      break;
  }
  if(i<n)
    mpz_set(order, factors[i]);
  else
  {
    printf("missing order in xtended!!\n");
    exit(-4);
  }
  poly_point_clear(&R);
  return i;
}

