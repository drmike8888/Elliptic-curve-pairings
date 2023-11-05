/****************************************************
 *                                                  *
 *    Simple subroutines for polynomial arithmetic. *
 *                                                  *
 ***************************************************/

#include "poly.h"

//#define DEBUG

static mpz_t *table = NULL;          // used for multiplication
static long poly_degree;
static POLY irrd;                    // irreducible basis polynomial (degree k)
static mpz_t ptok;                   // p^k

/*  initialize polynomial
    Enter with pointer to polynomial storage place and modulo prime n.
    Uses GMP init to initialize coefficents and sets degree to zero
    as well as all coefficients to zero.
*/

void poly_init(POLY *p)
{
  int i;

  for(i=0; i<MAXDEGREE; i++)
    mpz_init(p->coef[i]);
  p->deg = 0;
}

/* clear out polynomial coefficient storage  */

void poly_clear(POLY *p)
{
  int i;

  for(i=0; i<MAXDEGREE; i++)
    mpz_clear(p->coef[i]);
}

/*  print a polynomial to stdout.
    format same as pari/gp or sagemath input
*/

void poly_print(POLY a)
{
  int i;
  mpz_t prm;

  mget(prm);
  for(i=a.deg; i>0; i--)
  {
    if(mpz_cmp_ui(a.coef[i], 0))
      gmp_printf("Mod(%Zd, %Zd)*x^%d + ", a.coef[i], prm, i);
  }
  gmp_printf("Mod(%Zd, %Zd) ", a.coef[0], prm);
  printf("\n");
  mpz_clear(prm);
}

/* with string and new line after  */

void poly_printf(char *string, POLY a)
{
  printf("%s", string);
  poly_print(a);
  printf("\n");
}

/*  copy polynomial data from one structure to another.  */

void poly_copy(POLY *a, POLY b)
{
  int i;

  a->deg = b.deg;
  for(i=0; i<=b.deg; i++)
    mpz_set(a->coef[i], b.coef[i]);
}

/*  setup irreducible polynomial used for extension field.  */

void poly_irrd_set(POLY i)
{
  poly_init(&irrd);
  poly_copy(&irrd, i);
  mget(ptok);
  mpz_pow_ui(ptok, ptok, i.deg);
}

/*  return irreducible polynomial used for extension field.  */

void poly_irrd_get(POLY *i)
{
  poly_copy(i, irrd);
}

/*  return cardinality of extension field  */

void poly_q_get(mpz_t pk)
{
  mpz_init(pk);
  mpz_set(pk, ptok);
}

/*  Add two polynomials.  c = a + b  
    Over write operations allowed.   
*/
void poly_add(POLY *c, POLY a, POLY b)
{
  int i, dc;
  POLY rslt;

  poly_init(&rslt);
  if(a.deg > b.deg)
  {
    rslt.deg = a.deg;
    for(dc=a.deg; dc>b.deg; dc--)
      mpz_set(rslt.coef[dc], a.coef[dc]);
  }
  else if(b.deg > a.deg)
  {
    rslt.deg = b.deg;
    for(dc=b.deg; dc>a.deg; dc--)
      mpz_set(rslt.coef[dc], b.coef[dc]);
  }
  else
  {
    dc = a.deg;
    rslt.deg = a.deg;
  }
  while(dc >= 0)
  {
    madd(rslt.coef[dc], a.coef[dc], b.coef[dc]);
    dc--;
  }
  i = rslt.deg;
  while((i > 0) && (!mpz_cmp_ui(rslt.coef[i], 0)))
  {
    rslt.deg--;
    i--;
  }
  for(i=0; i<=rslt.deg; i++)
    mpz_set(c->coef[i], rslt.coef[i]);
  c->deg = rslt.deg;
  poly_clear(&rslt);
}

/*  Subtract two polynomials.  c = a - b  
    Over write operations allowed.   
*/
void poly_sub(POLY *c, POLY a, POLY b)
{
  int i;
  POLY bneg;

  poly_init(&bneg);
  bneg.deg = b.deg;
  for(i=0; i<=b.deg; i++)
  {
    mpz_init_set(bneg.coef[i], b.coef[i]);
    mneg(bneg.coef[i], bneg.coef[i]);
  }
  poly_add(c, a, bneg);
  poly_clear(&bneg);
}

/*  normalize a list of coefficients modulo modn.
    Enter with pointer to coefficients and degree n, 
    order 0 .. n. 
    operates in place.
*/

void poly_normal(POLY *a)
{
  int i, tst;
  mpz_t c;

  mpz_init(c);
  for(i=a->deg; i>=0; i--)
  {
    tst = mpz_cmp_ui(a->coef[i], 1);
    if(tst < 0)
      continue;
    if(!tst)
      return;                     // already normal
    minv(c, a->coef[i]);
    break;
  }
  if(i < 0) return;  // all zeros!
  while(i >= 0)
  {
    mmul(a->coef[i], a->coef[i], c);
    i--;
  }
}

#ifdef DEBUG
void poly_debug(int n)
{
  int i, j;

  for(i=0; i<2*n; i++)
  {
    for(j=n-1; j>=0; j--)
      gmp_printf("%Zd  ", table[i*n + j]);
    printf("\n");
  }
}
#endif

/*  Create lookup table for polynomial modular multiplication.
    The table will contain all powers of x from 0 to 2*n modulo
    the input f(x) polynomial (assumed monic).  The first 0 to 
    n-1 are trivial but it makes coding easier.  The rest of the 
    table will be polynomials of degree n-1 or less.  
*/

void poly_mulprep(POLY f)
{
  int i, j, tst;
  mpz_t tmp;
  POLY fnrml;

  poly_init(&fnrml);
  poly_copy(&fnrml, f);
  poly_degree = f.deg;
  if(table)
    free(table);
  table = (mpz_t*)malloc(sizeof(mpz_t)*poly_degree*poly_degree*2);
  for(i=0; i<2*poly_degree; i++)
    for(j=0; j<poly_degree; j++)
      mpz_init(table[poly_degree*i + j]);
  mpz_init(tmp);
  
/*  set lowest degree terms to x^j  */
  
  for(i=0; i<poly_degree; i++)
    mpz_set_ui(table[poly_degree*i + i], 1);

/*  set x^n entry to -f(x)  */

  poly_normal(&fnrml);
  for(j=0; j<poly_degree; j++)
  {
    mpz_neg(tmp, fnrml.coef[j]);
    mpz_set(table[poly_degree*poly_degree + j], tmp);
  }
  for(i=poly_degree+1; i<2*poly_degree; i++)
  {
/*   add x^n entry to rotated coefficients.  */

    for(j=1; j<poly_degree; j++)
    {
      mmul(tmp, table[poly_degree*poly_degree + j], table[i*poly_degree - 1]);
      madd(table[i*poly_degree + j], table[(i-1)*poly_degree + j - 1], tmp);
    }
    mmul(tmp, table[poly_degree*poly_degree], table[i*poly_degree - 1]);
    mpz_set(table[i*poly_degree], tmp);
  }
  mpz_clear(tmp);
  poly_clear(&fnrml);
#ifdef DEBUG
  poly_debug(poly_degree);
#endif
}

/* Multiply two polynomials modulo prime polynomial. 
   Assumes prime polynomial table already created with 
   mulprep routine.  Also assumes result space has 
   already been initialized.
*/

void poly_mul(POLY *rslt, POLY a, POLY b)
{
  int i, j, m, n;
  mpz_t coef[2*MAXDEGREE], acf[2*MAXDEGREE], bcf[2*MAXDEGREE];
  mpz_t tmp;

/*  initialize space */
  
  m = a.deg;
  n = b.deg;
  for(i=0; i<=n+m; i++)
    mpz_inits(coef[i], acf[i], bcf[i], NULL);
  for(i=0; i<=m; i++)
    mpz_set(acf[i], a.coef[i]);
  for(i=0; i<=n; i++)
    mpz_set(bcf[i], b.coef[i]);
  mpz_init(tmp);

/*  create double length list of coefficients  */
  
  for(i=0; i<=n+m; i++)
  {
    for(j=0; j<=i; j++)
    {
      mmul(tmp, acf[j], bcf[i - j]);
      madd(coef[i], coef[i], tmp);
    }
  }

/* combine upper powers with lower using table  */

  if(n+m < poly_degree)
  {
    rslt->deg = n+m;
    for(i=0; i<=n+m; i++)
      mpz_set(rslt->coef[i], coef[i]);
  }
  else
  {
    rslt->deg = poly_degree - 1;
    for(i=0; i<poly_degree; i++)
      mpz_set(rslt->coef[i], coef[i]);
    for(i=poly_degree; i<=n+m; i++)
    {
      for(j=0; j<poly_degree; j++)
      {
	mmul(tmp, coef[i], table[i*poly_degree + j]);
	madd(rslt->coef[j], rslt->coef[j], tmp);
      }
    }
  }

/*  check max degree has non zero coefficient  */

  while((mpz_cmp_ui(rslt->coef[rslt->deg], 0) <= 0) && (rslt->deg > 0))
    rslt->deg--;
  for(i=0; i<=n+m; i++)
    mpz_clears(coef[i], acf[i], bcf[i], NULL);
  mpz_clear(tmp);
}

/* polynomial version of Euclids division algorithm.
   Enter with polynomials a and b (for a/b).
   Returns q and r such that a = q*b + r.
   Assumes modulus set for coefficients, q and r initialized.
*/

void poly_euclid(POLY *q, POLY *r, POLY a, POLY b)
{
  int i, j;
  mpz_t s;
  POLY tmp;
  int k;
  
/*  set r = a, q = 0  */
  
  poly_copy(r, a);
  q->deg = 0;
  mpz_set_ui(q->coef[0], 0);
  if(b.deg > a.deg) return;
  
/*  perform division of each coefficient  */

  mpz_init(s);
  poly_init(&tmp);
  q->deg = a.deg - b.deg;
  while((r->deg >= b.deg) && r->deg)
  {
    j = r->deg - b.deg;
    mdiv(s, r->coef[r->deg], b.coef[b.deg]);
    mpz_set(q->coef[j], s);
    for(i=0; i<=b.deg; i++)
      mmul(tmp.coef[i + j], s, b.coef[i]);
    tmp.deg = r->deg;
    poly_sub(r, *r, tmp);
  }
  mpz_clear(s);
  poly_clear(&tmp);
}

/*  polynomial GCD routine.
    d = gcd(a, b) where a and b are polynomials
    and d has been initialzed.
*/

void poly_gcd(POLY *d, POLY a, POLY b)
{
  POLY aw, bw, q, r;

  poly_init(&aw);
  if(poly_cmp(a, aw))
  {
    poly_copy(d, b);
    poly_clear(&aw);
    return;
  }
  poly_init(&bw);
  if(poly_cmp(b, bw))
  {
    poly_copy(d, a);
    poly_clear(&aw);
    poly_clear(&bw);
    return;
  }
  if(a.deg >= b.deg)
  {
    poly_copy(&aw, a);
    poly_copy(&bw, b);
  }
  else
  {
    poly_copy(&aw, b);
    poly_copy(&bw, a);
  }
  poly_init(&q);
  poly_init(&r);
  while(bw.deg > 0)
  {
    poly_euclid(&q, &r, aw, bw);
    poly_copy(&aw, bw);
    poly_copy(&bw, r);
  }
  if(!mpz_cmp_ui(bw.coef[0], 0))
    poly_copy(d, aw);
  else
    poly_copy(d, bw);
  poly_clear(&r);
  poly_clear(&q);
  poly_clear(&bw);
  poly_clear(&aw);
}

/*  polynomial square and multiply with flag.
    input is x, output is x^2 flag = 0 and
    a*x^2 if flag is set.  
*/

void poly_sqm(POLY *x2, POLY x, POLY a, int flag)
{
  POLY tmp;

  poly_init(&tmp);
  poly_mul(&tmp, x, x);
  if(flag)
    poly_mul(x2, tmp, a);
  else
    poly_copy(x2, tmp);
  poly_clear(&tmp);
}

/*  compute x^p where p is the field prime.  
    both minit and poly_mulprep must have been run.
    input is previous run of poly_xp to get x^(p^i)
*/

void poly_xp(POLY *xp, POLY x)
{
  int i, bitcnt, bit;
  mpz_t prm;
  
  mget(prm);
  bitcnt = mpz_sizeinbase(prm, 2) - 2;
  poly_copy(xp, x);
  while(bitcnt >= 0)
  {
    bit = mpz_tstbit(prm, bitcnt);
    poly_sqm(xp, *xp, x, bit);
    bitcnt--;
  }
  mpz_clear(prm);
}

/*  compare 2 polynomials.  return 1 if equal, 0 if not  */

int poly_cmp(POLY a, POLY b)
{
  int i;

  if(a.deg != b.deg)
    return 0;
  for(i=a.deg; i>=0; i--)
    if(mpz_cmp(a.coef[i], b.coef[i]))
      return 0;
  return 1;
}

/*  Attempt to find irreducible trinomial of degree n.  
    Scans x^n + x + j for p-1 >= j >= 1.
    Clobbers mulprep table.  
    Returns 1 and sets f to first found irreducible trinomial, 
    returns 0 if no irreducible trinomial found, or n > MAXDEGREE. 
    Use Ben-Or from Modern Computer Algebra 14.40.
    Lidl & Niederreiter pg. 130: number of x^n + x + a being 
    irreducible is p/n, so for p >> n should be quick.
*/

int poly_irreducible(POLY *f, long n)
{
  POLY q, r, x, xp[MAXDEGREE/2], xpm1;
  mpz_t j, prime;
  long i, mlimt, done;
  
  if(n > MAXDEGREE)
    return 0;
  poly_init(&r);
  r.deg = n;
  mpz_set_ui(r.coef[n], 1);
  mpz_set_ui(r.coef[1], 1);
  poly_init(&x);
  x.deg = 1;
  mpz_set_ui(x.coef[1], 1);
  poly_init(&q);
  poly_init(&xpm1);
  mlimt = n/2;
  for(i=0; i<mlimt; i++)
    poly_init(&xp[i]);
  mget(prime);
  mpz_init_set_ui(j, 2);
  done = 0;
  while((mpz_cmp(j, prime) < 0) && !done)
  {
    mpz_set(r.coef[0], j);
    mpz_add_ui(j, j, 1);
    poly_mulprep(r);
    i = 0;
    q.deg = 0;
    while((i<mlimt) && !q.deg)
    {
      if(!i)
	poly_xp(&xp[0], x);
      else
	poly_xp(&xp[i], xp[i-1]);
      poly_sub(&xpm1, xp[i], x);
      poly_gcd(&q, xpm1, r);
      i++;
    }
    if(q.deg)
      continue;         // r divides x^p^i

/*  if we get this far, it is irreducible  */

    done = 1;
    poly_copy(f, r);
  }
  poly_clear(&r);
  poly_clear(&x);
  poly_clear(&q);
  for(i=0; i<mlimt; i++)
    poly_clear(&xp[i]);
  poly_clear(&xpm1);
  mpz_clears(j, prime, NULL);
  if(done)
    return 1;
  return 0;
}


/*  compute h(x) = g(x)^k mod f(x) 
    assumes mulprep was run, f(x) = irreducible polynomial
    and result space is initialized.
*/

void poly_pow(POLY *h, POLY g, mpz_t k)
{
  int bitcnt, bit;
  POLY a;

  poly_init(&a);
  poly_copy(&a, g);
  bitcnt = mpz_sizeinbase(k, 2) - 2;
  while(bitcnt >= 0)
  {
    bit = mpz_tstbit(k, bitcnt);
    poly_sqm(&a, a, g, bit);
    bitcnt--;
  }
  poly_copy(h, a);
  poly_clear(&a);
}

/*  compute h(x) = g(x)^{(p-1)/2} mod A(x)
    Ensure mulprep table setup.  Assumes result
    space already initialized.
*/

void gpow_p2(POLY *h, POLY g)
{
  mpz_t prm;

  mget(prm);
  mpz_sub_ui(prm, prm, 1);
  mpz_divexact_ui(prm, prm, 2);
  poly_pow(h, g, prm);
  mpz_clear(prm);
}

/*  Pseudo-division used in resultant algorithm.
    This version based on Cohen pg 112.
    Enter with two polynomials A(deg m), B(deg n)
    Returns Q, R with d^(m-n+1)*A = B*Q + R
*/

void poly_pseudo_div(POLY *Q, POLY *R, POLY A, POLY B)
{
  long e, i, k;
  mpz_t d;
  POLY S, T;
  
  poly_init(&S);
  poly_init(&T);
  poly_copy(R, A);
  Q->deg = 0;
  mpz_set_ui(Q->coef[0], 0);
  e = A.deg - B.deg + 1;
  mpz_init_set(d, B.coef[B.deg]);
  while(R->deg >= B.deg)
  {
    S.deg = R->deg - B.deg;
    mpz_set(S.coef[S.deg], R->coef[R->deg]);
    for(i=0; i<=Q->deg; i++)
      mmul(Q->coef[i], Q->coef[i], d);
    poly_add(Q, *Q, S);
    for(i=0; i<=R->deg; i++)
      mmul(R->coef[i], R->coef[i], d);
    k = S.deg;
    for(i=0; i<=B.deg; i++)
      mmul(T.coef[i + k], B.coef[i], S.coef[k]);
    T.deg = B.deg + k;
    poly_sub(R, *R, T);
    e--;
    mpz_set_ui(S.coef[k], 0);
  }
  if(e >= 1)
  {
    mpowi(d, d, e);
    for(i=0; i<=Q->deg; i++)
      mmul(Q->coef[i], Q->coef[i], d);
    for(i=0; i<=R->deg; i++)
      mmul(R->coef[i], R->coef[i], d);
  }
  mpz_clear(d);
  poly_clear(&S);
  poly_clear(&T);
}

/*  compute the content of a polynomial.
    this is just the gcd of all the coefficients.
*/

void poly_cont(mpz_t cont, POLY A)
{
  int i;
  mpz_t rslt;

  mpz_init_set(rslt, A.coef[A.deg]);
  for(i=A.deg-1; i>=0; i--)
  {
    mpz_gcd(rslt, rslt, A.coef[i]);
    if(!mpz_cmp_ui(rslt, 1))
      break;
  }
  mpz_set(cont, rslt);
  mpz_clear(rslt);
}

/*  compute the resultant of two polynomials A and B.
    Output is a mod p value.
*/
void poly_resltnt(mpz_t rsltnt, POLY A, POLY B)
{
  POLY Aa, Bb, Q, R;
  mpz_t g, h, ta, tb, a, b;
  long dlta, s, i;

/*  if either A or B is zero, resultant is zero  */
  
  if(!A.deg && !A.coef[0])
  {
    mpz_set_ui(rsltnt, 0);
    return;
  }
  if(!B.deg && !B.coef[0])
  {
    mpz_set_ui(rsltnt, 0);
    return;
  }

/* initialize local variables  */

  mpz_inits(g, h, ta, tb, a, b, NULL);
  poly_cont(a, A);
  poly_cont(b, B);
  poly_init(&Aa);
  if(mpz_cmp_ui(a, 1))  // content != 1
  {
    Aa.deg = A.deg;
    for(i=0; i<=A.deg; i++)
      mdiv(Aa.coef[i], A.coef[i], a);
    mpowi(ta, a, B.deg);
  }
  else
  {
    poly_copy(&Aa, A);
    mpz_set_ui(ta, 1);
  }
  poly_init(&Bb);
  if(mpz_cmp_ui(b, 1))  // content != 1
  {
    Bb.deg = B.deg;
    for(i=0; i<=B.deg; i++)
      mdiv(Bb.coef[i], B.coef[i], b);
    mpowi(tb, b, A.deg);
  }
  else
  {
    poly_copy(&Bb, B);
    mpz_set_ui(tb, 1);
  }
  mpz_set_ui(g, 1);
  mpz_set_ui(h, 1);
  poly_init(&Q);
  s = 1;
  if(A.deg < B.deg)  // make A larger than B
  {
    poly_copy(&Q, Aa);
    poly_copy(&Aa, Bb);
    poly_copy(&Bb, Q);
  }
  poly_init(&R);
  
/*  perform pseudo division  */
  
  while(Bb.deg > 0)
  {
    dlta = Aa.deg - Bb.deg;
    if((Aa.deg & 1) && (Bb.deg & 1))
      s = -s;
    poly_pseudo_div(&Q, &R, Aa, Bb);
    poly_copy(&Aa, Bb);
    mpowi(a, h, dlta);
    mmul(b, a, g);
    Bb.deg = R.deg;
    for(i=0; i<=R.deg; i++)
      mdiv(Bb.coef[i], R.coef[i], b);
    
/*  save new h, g values  */

    mpz_set(g, Aa.coef[Aa.deg]);
    i = 1 - dlta;
    mpowi(a, h, i);
    mpowi(b, g, dlta);
    mmul(h, a, b);
  }

/*  finished, compute final resultant  */
  
  i = 1 - Aa.deg;
  mpowi(a, h, i);
  mpowi(b, Bb.coef[0], Aa.deg);
  mmul(h, a, b);
  mmul(rsltnt, h, ta);
  mmul(rsltnt, rsltnt, tb);
  if(s < 0)
    mneg(rsltnt, rsltnt);
  poly_clear(&Aa);
  poly_clear(&Bb);
  poly_clear(&Q);
  poly_clear(&R);
  mpz_clears(g, h, ta, tb, a, b, NULL);
}

/*  test if polynomial is quadratic residue.
    Assumes irrd has already been set.  Will crash hard if not.
    return 1 if it is, 0 if non-residue.
*/

int poly_sqr(POLY x)
{
  mpz_t res, p, q;
  int k;
  
  mpz_inits(res, q, NULL);
  poly_resltnt(res, x, irrd);
  mget(p);
  mpz_sub_ui(q, p, 1);
  mpz_div_ui(q, q, 2);
  mpz_powm(res, res, q, p);
  k = mpz_cmp_ui(res, 1);
  mpz_clears(res, p, q, NULL);
  if(!k)
    return 1;
  return 0;
}

/*  compute square root of a polynomial modulo irreducible polynomial.
    Since this is over a field extension, (p^k + 1)/4 works if p^k
    is congruent to 3 mod 4.  Otherwise we use Tonelli-Shanks using
    polynomials instead of GF(p).
                sqt = sqrt(a)  if a is quadratic residue
		garbage otherwise
*/

void poly_sqrt(POLY *sqt, POLY a)
{
  long r, i, ck, m, m2;
  mpz_t pk, q;
  POLY x, b, y, one, bpw, t;

/*  see if p^k = 3 mod 4  */
  
  poly_q_get(pk);
  mpz_init(q);
  if(mpz_tstbit(pk, 0) && mpz_tstbit(pk, 1))
  {
    poly_init(&y);
    mpz_add_ui(q, pk, 1);
    mpz_divexact_ui(q, q, 4);
    poly_pow(&y, a, q);
    poly_copy(sqt, y);
    poly_clear(&y);
    mpz_clears(q, pk, NULL);
    return;
  }

/*  p^k = 1 mod 4, so do Tonelli-Shanks  */

  mpz_sub_ui(q, pk, 1);
  r = 0;
  while(!mpz_tstbit(q, 0))
  {
    mpz_divexact_ui(q, q, 2);
    r++;
  }
  poly_init(&y);
  ck = 1;
  while(ck)
  {
    poly_rand(&y);
    ck = poly_sqr(y);
  }
  poly_pow(&y, y, q);
  poly_init(&b);
  poly_pow(&b, a, q);
  poly_init(&x);
  mpz_add_ui(pk, q, 1);
  mpz_divexact_ui(pk, pk, 2);
  poly_pow(&x, a, pk);
  poly_init(&one);
  one.deg = 0;
  mpz_set_ui(one.coef[0], 1);
  poly_init(&bpw);
  poly_init(&t);
  while(!poly_cmp(b, one))
  {
    m = 0;
    while(!poly_cmp(bpw, one))
    {
      m++;
      m2 = 1 << m;
      mpz_set_ui(pk, m2);            // b^(2^m)
      poly_pow(&bpw, b, pk);
      if(m == r)
      {
	printf("square root failed\n");
	return;
      }
    }
    mpz_set_ui(bpw.coef[0], 0);
    i = r - m - 1;
    m2 = 1 << i;
    mpz_set_ui(pk, m2);            // y^(2^i)
    poly_pow(&t, y, pk);
    poly_mul(&y, t, t);           // t^2
    r = m;
    poly_mul(&x, x, t);
    poly_mul(&b, b, y);
    if(i < 0) exit(0);
  }
  poly_copy(sqt, x);
  poly_clear(&x);
  poly_clear(&b);
  poly_clear(&y);
  poly_clear(&one);
  poly_clear(&bpw);
  poly_clear(&t);
  mpz_clears(pk, q, NULL);
}

/*  compute polynomial inverse modulo irrd.
    a = 1/b mod irrd
*/

void poly_invert(POLY *a, POLY b)
{
  mpz_t rho;
  POLY r, u, q, one, t, w, tmp, y, v;
  int done, i;

  mpz_init(rho);
  poly_init(&v);
  poly_copy(&v, irrd);
  poly_init(&r);
  poly_init(&u);
  poly_copy(&u, b);
  poly_normal(&u);
  poly_init(&w);
  poly_init(&q);
  poly_init(&one);
  mpz_set_ui(one.coef[0], 1);
  poly_init(&t);
  poly_init(&y);
  minv(y.coef[0], b.coef[b.deg]);
  poly_init(&tmp);
  done = 0;
  while(!done)
  {
    poly_euclid(&q, &r, v, u);
    minv(rho, r.coef[r.deg]);
    poly_normal(&r);
    poly_mul(&tmp, q, y);
    poly_sub(&t, w, tmp);
    for(i=0; i<=t.deg; i++)
      mmul(t.coef[i], t.coef[i], rho);
    poly_copy(&w, y);
    poly_copy(&y, t);
    poly_copy(&v, u);
    poly_copy(&u, r);
    done = poly_cmp(r, one);
  }
  poly_copy(a, t);
  mpz_clear(rho);
  poly_clear(&r);
  poly_clear(&u);
  poly_clear(&q);
  poly_clear(&one);
  poly_clear(&t);
  poly_clear(&w);
  poly_clear(&tmp);
  poly_clear(&y);
  poly_clear(&v);
}

/*  divide two polynomials which are modulo irrd.
    Computes a = b / c mod irrd.
*/

void poly_div(POLY *a, POLY b, POLY c)
{
  POLY q;
  int i;
  
  poly_init(&q);
  poly_invert(&q, c);
  poly_mul(a, b, q);
  poly_clear(&q);
}

/*  create random polynomial with coefficients
    one degree less than irreducible polynomial.
*/

void poly_rand(POLY *rnd)
{
  int i;

  rnd->deg = irrd.deg - 1;
  for(i=0; i<irrd.deg; i++)
    mrand(rnd->coef[i]);
}
