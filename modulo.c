/**************************************************************
 *                                                            *
 *     Fundamental modular arithmetic fuctions using GMP.     *
 *                                                            *
 *************************************************************/

#include "modulo.h"

static mpz_t modulus;
static gmp_randstate_t state;

/* Generic modulo operations  */

/*  a = b + c mod n  */

void mod_add(mpz_t a, mpz_t b, mpz_t c, mpz_t n)
{
  mpz_t rslt;

  mpz_init(rslt);
  mpz_add(rslt, b, c);
  mpz_mod(a, rslt, n);
  mpz_clear(rslt);
}

/*  a = b - c mod n  */

void mod_sub(mpz_t a, mpz_t b, mpz_t c, mpz_t n)
{
  mpz_t rslt;

  mpz_init(rslt);
  mpz_sub(rslt, b, c);
  mpz_mod(a, rslt, n);
  mpz_clear(rslt);
}

/*  a = b * c mod n  */

void mod_mul(mpz_t a, mpz_t b, mpz_t c, mpz_t n)
{
  mpz_t rslt;

  mpz_init(rslt);
  mpz_mul(rslt, b, c);
  mpz_mod(a, rslt, n);
  mpz_clear(rslt);
}

/*  a = b / c mod n  */

void mod_div(mpz_t a, mpz_t b, mpz_t c, mpz_t n)
{
  mpz_t rslt;

  mpz_init(rslt);
  if(!mpz_invert(rslt, c, n))
  {
    printf("division by zero in div_mod!\n");
    mpz_clear(rslt);
    exit(-1);
  }
  mpz_mul(rslt, b, rslt);
  mpz_mod(a, rslt, n);
  mpz_clear(rslt);
}

/*  a = - b mod n  */

void mod_neg(mpz_t a, mpz_t b, mpz_t n)
{
  mpz_t rslt;

  mpz_init(rslt);
  mpz_neg(rslt, b);
  mpz_mod(a, rslt, n);
  mpz_clear(rslt);
}
  
/*  Fixed modulus operations  */

void minit(mpz_t m)
{
  mpz_init_set(modulus, m);
  gmp_randinit_mt(state);
}

/*  return local modulus  */

void mget(mpz_t mod)
{
  mpz_init_set(mod, modulus);
}

/* reset modulus (used in snark stuff) */

void mset(mpz_t prm)
{
  mpz_set(modulus, prm);
}

/*  a = b + c mod n  */

void madd(mpz_t a, mpz_t b, mpz_t c)
{
  mod_add(a, b, c, modulus);
}

/*  a = b - c mod n  */

void msub(mpz_t a, mpz_t b, mpz_t c)
{
  mod_sub(a, b, c, modulus);
}

/*  a = b * c mod n  */

void mmul(mpz_t a, mpz_t b, mpz_t c)
{
  mod_mul(a, b, c, modulus);
}

/*  a = b / c mod n  */

void mdiv(mpz_t a, mpz_t b, mpz_t c)
{
  mod_div(a, b, c, modulus);
}

/*  a = 1 / b mod n  */

void minv(mpz_t a, mpz_t b)
{
  mpz_invert(a, b, modulus);
}

/*  a = - b mod n  */

void mneg(mpz_t a, mpz_t b)
{
  mod_neg(a, b, modulus);
}

/* return a random value in range
   0 to modulus - 1.  minit() must
   have been run first.
*/

void mrand(mpz_t rand)
{
  mpz_urandomm(rand, state, modulus);
  while(!mpz_cmp_ui(rand, 0) || !mpz_cmp_ui(rand, 1))
	mpz_urandomm(rand, state, modulus);
}

/*  return legendre symbol of input.  */

int msqr(mpz_t x)
{
  return mpz_legendre(x, modulus);
}

/*  a = b^i mod n  */

void mpowi(mpz_t a, mpz_t b, long i)
{
  if(i < 0)
  {
    minv(a, b);
    mpz_powm_ui(a, a, -i, modulus);
  }
  else if(!i)
    mpz_set_ui(a, 1);
  else
    mpz_powm_ui(a, b, i, modulus);
}

/*  compute square root mod n
    x = sqrt(a) mod n  
    Return 1 on success, 0 if a is not
    a quadratic residue mod n.  */

int msqrt(mpz_t x, mpz_t a)
{
  mpz_t n, q, y, b, t, t1;
  long e, i, r, cmp, m;

  if(!msqr(a))  // no point if not a quadratic residue
    return 0;
  
/* is modulus ~ 3 mod 4?  */

  if(mpz_tstbit(modulus, 0) && mpz_tstbit(modulus, 1))
  {
    mpz_init_set(q, modulus);
    mpz_add_ui(q, q, 1);
    mpz_divexact_ui(q, q, 4);
    mpz_powm(x, a, q, modulus);
    mpz_clear(q);
    return 1;
  }
  mpz_inits(n, q, y, b, t, t1, NULL);

/*  p - 1 = 2^e * q  */
  
  mpz_sub_ui(q, modulus, 1);
  e = mpz_scan1(q, 0);
  i = e;
  while(i)
  {
    mpz_divexact_ui(q, q, 2);
    i--;
  }

/*  find a generator  */
  
  i = 1;
  while(i >= 0)
  {
    mrand(n);
    i = mpz_legendre(n, modulus);
  }
 

/*  initialize working components  */

  mpz_powm(y, n, q, modulus);
  r = e;
  mpz_sub_ui(q, q, 1);
  mpz_divexact_ui(q, q, 2);
  mpz_powm(x, a, q, modulus);
  mmul(b, x, x);
  mmul(b, b, a);
  mmul(x, x, a);
  
/*  loop on algorithm until finished or failure.  */
  
  cmp = mpz_cmp_ui(b, 1);
  while(cmp)
  {
    m = 1;
    mpz_set(t1, b);
    while(m < r)
    {
      mpowi(t1, t1, 2);
      if(!mpz_cmp_ui(t1, 1))
	break;
      m++;
    }
    if(r == m)
    {
       mpz_clears(n, q, y, b, t, t1, NULL);
       return 0;
    }
    i = r - m - 1;
    mpz_set(t, y);
    while(i)
    {
      mpowi(t, t, 2);
      i--;
    }
    mmul(y, t, t);
    r = m;
    mmul(x, x, t);
    mmul(b, b, y);
    cmp = mpz_cmp_ui(b, 1);
  }
  mpz_clears(n, q, y, b, t, t1, NULL);
  return 1;
}
