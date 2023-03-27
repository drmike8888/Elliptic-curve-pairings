/****************************************************************
 *                                                              *
 *    Routines to compute point addition and multiplication     *
 *    over an elliptic curve on finite field.  Embedding points *
 *    on curve and other utilities included.  Routines are      *
 *    security focused rather than efficiency focused.          *
 *                                                              *
 ***************************************************************/

#include "modulo.h"
#include "eliptic.h"

/*  structure initialization routines  */

void point_init(POINT *P)
{
  mpz_inits(P->x, P->y, NULL);
}

void point_clear(POINT *P)
{
  mpz_clears(P->x, P->y, NULL);
}

void curve_init(CURVE *E)
{
  mpz_inits(E->a4, E->a6, NULL);
}

void curve_clear(CURVE *E)
{
  mpz_clears(E->a4, E->a6, NULL);
}

/*  duplicate a point  */

void point_copy(POINT *R, POINT P)
{
  mpz_set(R->x, P.x);
  mpz_set(R->y, P.y);
}

/*  test if point is at infinity.  Since all our curves
    have a6 != 0, the point (0, 0) never happens - use
    that as point at infinity. Return 1 if at infinity,
    0 otherwise.
*/

int test_point(POINT P)
{
  if(!mpz_cmp_ui(P.x, 0) && !mpz_cmp_ui(P.y, 0))
    return 1;
  return 0;
}

/*  output x, y values of a point.  Enter with string for label.  */

void point_printf(char *str, POINT P)
{
  printf("%s", str);
  gmp_printf("(%Zd,  %Zd)\n", P.x, P.y);
}

/*  compute RHS of elliptic curve equation
      f(x) = x^3 + a4*x + a6
    Enter with x and curve parameters
    Return f(x) value (which should already be initialzied)
*/

void fofx(mpz_t f, mpz_t x, CURVE E)
{
  mpz_t t1, t2;

  mpz_inits(t1, t2, NULL);
  mmul(t1, x, x);
  mmul(t1, t1, x);
  mmul(t2, E.a4, x);
  madd(f, t1, t2);
  madd(f, f, E.a6);
  mpz_clears(t1, t2, NULL);
}

/*  embed random value onto a curve.
    Enter with value x, and curve E.  Increments x (mod p)
    by 1 until point found on curve.
    Returns 2 points, first one with smaller y value.
*/

void elptic_embed(POINT *P1, POINT *P2, mpz_t x, CURVE E)
{
  mpz_t f, one;
  int done;

  mpz_init(f);
  mpz_init_set_ui(one, 1);
  mpz_set(P1->x, x);
  done = 0;
  while(!done)
  {
    fofx(f, P1->x, E);
    if(msqr(f) > 0)
      done = 1;
    else
      madd(P1->x, P1->x, one);
  }
  mpz_set(P2->x, P1->x);
  msqrt(P1->y, f);
  mneg(P2->y, P1->y);
  done = mpz_cmp(P2->y, P1->y);
  if(done < 0)
    mpz_swap(P2->y, P1->y);
  mpz_clears(f, one, NULL);
}

/*  perform R = P + Q using "side channel proof" form.
    This routine is used for both addition and doubling.
    Enter with points P and Q.  
    Returns point R (space should already be initialized)
*/

void elptic_sum(POINT *R, POINT P, POINT Q, CURVE E)
{
  mpz_t lmbda, ltp, lbt, t1, t2, t3;
  POINT rslt;
  
/*  see if either input is point at infinity  */

  if(test_point(P))
  {
    point_copy(R, Q);
    return;
  }
  if(test_point(Q))
  {
    point_copy(R, P);
    return;
  }

/*  compute lambda using general form  */

  mpz_inits(t1, t2, t3, ltp, lbt, lmbda, NULL);
  mmul(t1, P.x, P.x);      // x1*x1
  mmul(t2, P.x, Q.x);      // x1*x2
  mmul(t3, Q.x, Q.x);      // x2*x2
  madd(ltp, t1, t2);
  madd(ltp, ltp, t3);
  madd(ltp, ltp, E.a4);    // summed
  madd(lbt, P.y, Q.y);     // y1 + y2
  if(!mpz_cmp_ui(lbt, 0))  // if(P = - Q)
  {
    msub(lbt, Q.x, P.x);   // x2 - x1
    if(!mpz_cmp_ui(lbt, 0))// Really P == -Q?
    {
      mpz_set_ui(R->x, 0);   // return point at infinity
      mpz_set_ui(R->y, 0);
      mpz_clears(t1, t2, t3, ltp, lbt, lmbda, NULL);
      return;
    }
    msub(ltp, Q.y, P.y);   // y2 - y1
  }
  mdiv(lmbda, ltp, lbt);

/*  finally compute resulting point  */

  point_init(&rslt);
  mmul(t1, lmbda, lmbda);   // lambda^2
  madd(t2, P.x, Q.x);       // x1 + x2
  msub(rslt.x, t1, t2);       // x3 = difference
  msub(t1, P.x, rslt.x);      // x1 - x3
  mmul(t2, t1, lmbda);      // *lambda
  msub(rslt.y, t2, P.y);      // - y1
  point_copy(R, rslt);
  mpz_clears(t1, t2, t3, ltp, lbt, lmbda, NULL);
  point_clear(&rslt);
}

/*  Elliptic curve point multiplication Q = k*P  
    Enter with pre-initialzed storage Q
    point P, constant k and curve E.
*/

void elptic_mul(POINT *Q, POINT P, mpz_t k, CURVE E)
{
  int bit, j;
  POINT R;

  point_init(&R);
  point_copy(&R, P);
  j = mpz_sizeinbase(k, 2) - 2;
  while(j >= 0)
  {
    elptic_sum(&R, R, R, E);
    bit = mpz_tstbit(k, j);
    if(bit)
      elptic_sum(&R, R, P, E);
    j--;
  }
  point_copy(Q, R);
  point_clear(&R);
}

/*  create a random point on a curve.
    Enter with place for point to go 
    and curve to put it on.
*/

void point_rand(POINT *P, CURVE E)
{
  mpz_t r;
  POINT mP;
  
  mpz_init(r);
  mrand(r);
  point_init(&mP);
  if(mpz_tstbit(r, 0))
    elptic_embed(P, &mP, r, E);
  else
    elptic_embed(&mP, P, r, E);
  mpz_clear(r);
  point_clear(&mP);
}
