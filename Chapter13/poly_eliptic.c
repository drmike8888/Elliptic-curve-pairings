/***************************************************************
 *                                                             *
 *   Utility routines for eliptic curves on extension fields.  *
 *   Same as eliptic.c, but for polynomial arguments.          *
 *                                                             *
 **************************************************************/

#include "poly_eliptic.h"

/*  structure initialization routines  */

void poly_point_init(POLY_POINT *P)
{
  int i;

  P->x.deg = 0;
  P->y.deg = 0;
  for(i=0; i<MAXDEGREE; i++)
    mpz_inits(P->x.coef[i], P->y.coef[i], NULL);
}

void poly_point_clear(POLY_POINT *P)
{
  int i;

  for(i=0; i<MAXDEGREE; i++)
    mpz_clears(P->x.coef[i], P->y.coef[i], NULL);
}

void poly_curve_init(POLY_CURVE *E)
{
  int i;

  E->a4.deg = 0;
  E->a6.deg = 0;
  for(i=0; i<MAXDEGREE; i++)
    mpz_inits(E->a4.coef[i], E->a6.coef[i], NULL);
}

void poly_curve_clear(POLY_CURVE *E)
{
  int i;

  for(i=0; i<MAXDEGREE; i++)
    mpz_clears(E->a4.coef[i], E->a6.coef[i], NULL);
}

/*  duplicate a point  */

void poly_point_copy(POLY_POINT *R, POLY_POINT P)
{
  int i;

  poly_copy(&R->x, P.x);
  poly_copy(&R->y, P.y);
}

/*  test if point is at infinity.  Since all our curves
    have a6 != 0, the point (0, 0) never happens - use
    that as point at infinity. Return 1 if at inifinity,
    0 otherwise.
*/

int poly_test_point(POLY_POINT P)
{
  int i;

  if(P.x.deg || P.y.deg)
    return 0;
  if(!mpz_cmp_ui(P.x.coef[0], 0) && !mpz_cmp_ui(P.y.coef[0], 0))
    return 1;
  return 0;
}

/*  output x, y values of a poly point.  Enter with string for label.  */

void poly_point_printf(char *str, POLY_POINT P)
{
  printf("%s", str);
  printf("x: ");
  poly_print(P.x);
  poly_printf("y: ", P.y);
}

/*  output a4, a6 values of a poly curve.  Enter with string for label.  */

void poly_curve_printf(char *str, POLY_CURVE E)
{
  printf("%s", str);
  poly_printf("a4: ", E.a4);
  poly_printf("a6: ", E.a6);
}

/*  compute RHS of elliptic curve equation
      f(x) = x^3 + a4*x + a6
    Enter with x and curve parameters
    Return f(x) value (which should already be initialzied)
    NOTE: mulprep routine should be done with prime polynomial!
*/

void poly_fofx(POLY *f, POLY x, POLY_CURVE E)
{
  POLY t1, t2;

  poly_init(&t1);
  poly_init(&t2);
  poly_mul(&t1, x, x);
  poly_mul(&t1, t1, x);
  poly_mul(&t2, E.a4, x);
  poly_add(f, t1, t2);
  poly_add(f, *f, E.a6);
  poly_clear(&t1);
  poly_clear(&t2);
}

/*  bump polynomial finite field.  
    increments lowest value until overflow.  Works
    up chain.  Operates in place.
*/

void FF_bump(POLY *x)
{
  int i;
  mpz_t one;
  POLY ird;
  
  mpz_init_set_ui(one, 1);
  poly_init(&ird);
  poly_irrd_get(&ird);
  i = 0;
  while(i < ird.deg)
  {
    madd(x->coef[i], x->coef[i], one);
    if(mpz_cmp_ui(x->coef[i], 0))
      return;
    i++;
    if((i > x->deg) && (x->deg < ird.deg))
      x->deg++;
  }
  mpz_clear(one);
  poly_clear(&ird);
}

/*  embed random value onto a curve.
    Enter with value x, and curve E.  
    Increments x (mod irreducible polynomial) by 1 until point found on curve.
    Returns 2 points, first one with smaller y value (from highest degree down).
*/

void poly_elptic_embed(POLY_POINT *P1, POLY_POINT *P2, POLY x, POLY_CURVE E)
{
  POLY f;
  int done, i;
  mpz_t tmp;
  
  poly_init(&f);
  poly_copy(&P1->x, x);
  done = 0;
  while(!done)
  {
    poly_fofx(&f, P1->x, E);
    if(poly_sqr(f) > 0)
      done = 1;
    else
      FF_bump(&(P1->x));
  }
  poly_copy(&(P2->x), P1->x);
  poly_sqrt(&(P1->y), f);
  for(i=0; i<=P1->y.deg; i++)
    mneg(P2->y.coef[i], P1->y.coef[i]);
  P2->y.deg = P1->y.deg;
  done = mpz_cmp(P2->y.coef[P2->y.deg], P1->y.coef[P1->y.deg]);
  if(done < 0)
  {
    mpz_init(tmp);
    for(i=0; i<=P1->y.deg; i++)
    {
      mpz_set(tmp, P1->y.coef[i]);
      mpz_set(P1->y.coef[i], P2->y.coef[i]);
      mpz_set(P2->y.coef[i], tmp);
    }
    mpz_clear(tmp);
  }
  poly_clear(&f);
}

/*  perform R = P + Q using "side channel proof" form.
    This routine is used for both addition and doubling.
    Enter with points P and Q.  
    Returns point R (space should already be initialized)
*/

void poly_elptic_sum(POLY_POINT *R, POLY_POINT P, POLY_POINT Q, POLY_CURVE E)
{
  POLY lmbda, ltp, lbt, t1, t2, t3;
  POLY_POINT rslt;
  
/*  see if either input is point at infinity  */

  if(poly_test_point(P))
  {
    poly_point_copy(R, Q);
    return;
  }
  if(poly_test_point(Q))
  {
    poly_point_copy(R, P);
    return;
  }

/*  compute lambda using general form  */

  poly_init(&t1);
  poly_init(&t2);
  poly_init(&t3);
  poly_init(&ltp);
  poly_init(&lbt);
  poly_init(&lmbda);
  poly_mul(&t1, P.x, P.x);
  poly_mul(&t2, P.x, Q.x);
  poly_mul(&t3, Q.x, Q.x);
  poly_add(&ltp, t1, t2);
  poly_add(&ltp, ltp, t3);
  poly_add(&ltp, ltp, E.a4);
  poly_add(&lbt, P.y, Q.y);

  if(!lbt.deg && !mpz_cmp_ui(lbt.coef[0], 0))  // if(P = - Q)
  {
    poly_sub(&lbt, Q.x, P.x);      // x2 - x1
    if(!lbt.deg && !mpz_cmp_ui(lbt.coef[0], 0)) // Really P == -Q?
    {
      R->x.deg = 0;
      mpz_set_ui(R->x.coef[0], 0);   // return point at infinity
      R->y.deg = 0;
      mpz_set_ui(R->y.coef[0], 0);
      poly_clear(&t1);
      poly_clear(&t2);
      poly_clear(&t3);
      poly_clear(&ltp);
      poly_clear(&lbt);
      poly_clear(&lmbda);
      return;
    }
    poly_sub(&ltp, Q.y, P.y);
  }

/*  finally compute resulting point  */

  poly_div(&lmbda, ltp, lbt);
  poly_point_init(&rslt);
  poly_mul(&t1, lmbda, lmbda);
  poly_add(&t2, P.x, Q.x);
  poly_sub(&rslt.x, t1, t2);
  poly_sub(&t1, P.x, rslt.x);
  poly_mul(&t2, t1, lmbda);
  poly_sub(&rslt.y, t2, P.y);
  poly_point_copy(R, rslt);

  poly_clear(&t1);
  poly_clear(&t2);
  poly_clear(&t3);
  poly_clear(&ltp);
  poly_clear(&lbt);
  poly_clear(&lmbda);
  poly_point_clear(&rslt);
}

/*  Elliptic curve point multiplication Q = k*P  
    Enter with pre-initialzed storage Q
    point P, constant k and curve E.
*/

void poly_elptic_mul(POLY_POINT *Q, POLY_POINT P, mpz_t k, POLY_CURVE E)
{
  int bit, j;
  POLY_POINT R;

  poly_point_init(&R);
  poly_point_copy(&R, P);
  j = mpz_sizeinbase(k, 2) - 2;
  while(j >= 0)
  {
    poly_elptic_sum(&R, R, R, E);
    bit = mpz_tstbit(k, j);
    if(bit)
      poly_elptic_sum(&R, R, P, E);
    j--;
  }
  poly_point_copy(Q, R);
  poly_point_clear(&R);
}

/*  create a random point on a curve.
    Enter with place for point to go 
    and curve to put it on.
*/

void poly_point_rand(POLY_POINT *P, POLY_CURVE E)
{
  POLY r;
  POLY_POINT mP;
  
  poly_init(&r);
  poly_rand(&r);
  poly_point_init(&mP);
  if(mpz_tstbit(r.coef[0], 0))
    poly_elptic_embed(P, &mP, r, E);
  else
    poly_elptic_embed(&mP, P, r, E);
  poly_clear(&r);
  poly_point_clear(&mP);
}



