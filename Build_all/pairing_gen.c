/*****************************************************************
 *                                                               *
 *     Create pairing friendly curves using algorithms 6.19,     *
 *     6.20 in "A Taxonomy of Pairing-Friendly Elliptic Curves"  *
 *     D. Freeman, M. Scott, and E Teske (2010)                  *
 *                                                               *
 ****************************************************************/

#include "modulo.h"
#include "poly.h"
#include <math.h>

/*  Compute cyclotomic polynomial Phi_4k(z).
    Enter with k, alpha and x.
    Checks k is prime, and returns -1 if not.
    z = alpha*x^2  
    Assumes k <= 37.  returns -1 if over that.
*/

int phi4k(mpz_t r, long k, mpz_t alpha, mpz_t x)
{
  int i;
  mpz_t z[36], ck;

  if(k > 37)
    return -1;
  mpz_init_set_ui(ck, k);
  if(mpz_probab_prime_p(ck, 5) != 2)
    return -1;
  mpz_set_ui(ck, 1);
  mpz_init(z[0]);
  mpz_mul(z[0], x, x);
  mpz_mul(z[0], alpha, z[0]);
  mpz_sub(ck, ck, z[0]);
  for(i=1; i<k-1; i++)
  {
    mpz_init(z[i]);
    mpz_mul(z[i], z[i - 1], z[0]);
    if(i & 1)
      mpz_add(ck, ck, z[i]);
    else
      mpz_sub(ck, ck, z[i]);
  }
  mpz_set(r, ck);
  mpz_clear(ck);
  for(i=0; i<k-1; i++)
    mpz_clear(z[i]);
  return 1;
}

/*  Choose k and alpha based on log2(r).
    Values based on table 1.1 in Taxonomy.
    Enter with log2r, returns k, alpha and
    x starting point.
*/

void xstart(long lg2r, long *k, long *max, long *alphabase)
{
  double lgu, lgx, lga, ud, xd, w;

/*  choose k and alpha base from limits on table  */
  
  if(lg2r < 160)
  {
    *k = 7;
    *alphabase = 3;
  }
  else if(lg2r < 192)
  {
    *k = 7;
    *alphabase = 3;
  }
  else if(lg2r < 224)
  {
    *k = 11;
    *alphabase = 3;
  }
  else if(lg2r < 256)
  {
    *k = 19;
    *alphabase = 3;
  }
  else if(lg2r < 320)
  {
    *k = 19;
    *alphabase = 3;
  }
  else if(lg2r < 384)
  {
    *k = 19;
    *alphabase = 43;
  }
  else if(lg2r < 448)
  {
    *k = 23;
    *alphabase = 67;
  }
  else
  {
    *k = 31;
    *alphabase = 3;
  }

/*  alpha = u^2 * alphabase.  r is order (alpha*x^2)^(k-1)
    so compute u and x with assumption alpha and x are about
    the same bit size to start with.
*/
  lga = log2(*alphabase)/2.0;
  w = (double)lg2r/(*k - 1.0)/2.0 - lga;
  *max = exp2(w + .5);
}

/*  compute q(z) using equation 6.20 and 6.19 from
    Taxonomy paper.
    Enter with k, alpha and x. 
    If these parameters give 4q exactly, 
    Returns q (which should have already been initialized)
    and value 1, otherwise returns 0 with q untouched.
*/

int qofz(mpz_t q, long k, mpz_t alpha, mpz_t x)
{
  long k1, k2, prm;
  mpz_t  t1, t2, t3, t4, z;

  mpz_inits(z, t1, t2, t3, t4, NULL);
  k1 = k + 1;
  k2 = k1/2;
  mpz_mul(z, x, x);
  mpz_mul(z, z, alpha);
  mpz_pow_ui(t1, z, k1);
  mpz_pow_ui(t2, z, k);
  mpz_add_ui(t4, z, 1);
  mpz_pow_ui(t3, z, k2);
  mpz_mul_ui(t3, t3, 4);
  if(k2 & 1)
    mpz_sub(t1, t1, t3);
  else
    mpz_add(t1, t1, t3);
  mpz_add(t1, t1, t2);
  mpz_add(t1, t1, t4);
  prm = 0;
  if(mpz_divisible_2exp_p(t1, 2))
  {
    mpz_fdiv_q_ui(q, t1, 4);
    prm = 1;
  }
  mpz_clears(z, t1, t2, t3, t4, NULL);
  if(prm)
    return 1;
  return 0;
}

/*  create alpha from u and alphabase as mpz value. */

void mkalpha(mpz_t alpha, long u, long alphabase)
{
  long alph;

  alph = u*u*alphabase;
  mpz_set_ui(alpha, alph);
}

/*  compute t(z) using equation 6.20 and 6.19 from
    Taxonomy paper.  (trace of Frobenius)
    Enter with k, alpha and x. 
*/

void tofz(mpz_t t, long k, mpz_t alpha, mpz_t x)
{
  mpz_t z, t1, t2;
  long k1;
  
  mpz_inits(z, t1, t2, NULL);
  mpz_mul(z, x, x);
  mpz_mul(z, z, alpha);
  k1 = (k + 1)/2;
  mpz_pow_ui(t1, z, k1);
  mpz_set_ui(t2, 1);
  if(k1 & 1)
    mpz_sub(t, t2, t1);
  else
    mpz_add(t, t2, t1);
  mpz_clears(z, t1, t2, NULL);
}

int main(int argc, char *argv[])
{
  FILE *pair;
  mpz_t r, alpha, x, q, x0, t;
  long k, u, lg2r, alphabase, rpm, qpm, rsz, qsz;
  long max;
  double rho;
  int j, jlo, jhi, m;
  double xsz, asz;
  char filename[128];
  
  if(argc < 2)
  {
    printf("Use: ./pairing_gen <log2(r)>\n");
    printf("   where log2(r) is 2 to 512\n");
    exit(-1);
  }
  lg2r = atol(argv[1]);
  sprintf(filename, "pair.%03ld", lg2r);
  pair = fopen(filename, "w");
  if(lg2r < 2)
  {
    printf("Come on, that's too small!\n");
    exit(-1);
  }
  if(lg2r > 576)
    printf("OK, but security will be questionable.\n");

  mpz_inits(r, alpha, x, x0, q, NULL);
  xstart(lg2r, &k, &max, &alphabase);
  fprintf(pair, "k= %ld alphabase = %ld max = %ld\n", k, alphabase, max);
  for(m = 1; m<max; m++)
  {
    mpz_set_ui(x0, m);
    jlo = max/(m+1);
    jhi = max/m;
    if(jlo == jhi) jhi++;
    fprintf(pair, "%d %d %d\n", m, jlo, jhi);
    for(j=jlo; j<jhi; j++)
    {
      mkalpha(alpha, j, alphabase);
//      gmp_printf("k= %ld alpha = %Zd  x = %Zd\n", k, alpha, x0);
      phi4k(r, k, alpha, x0);
      rpm = mpz_probab_prime_p(r, 25);
      rsz = mpz_sizeinbase(r, 2);
//      if(rsz > lg2r + 5)
//	continue;
//      gmp_printf("r = %Zd  numbits: %ld\n", r, rsz);
      if(rpm)
      {
	if(qofz(q, k, alpha, x0))
        {
	  qsz = mpz_sizeinbase(q, 2);
//	  gmp_printf("q: %Zd  bits: %ld\n", q, qsz);
	  qpm = mpz_probab_prime_p(q, 25);
	  if(qpm)
	  {
	    tofz(t, k, alpha, x0);
	    gmp_fprintf(pair, "k= %ld alpha = %Zd  x = %Zd\n", k, alpha, x0);
	    gmp_fprintf(pair, "r = %Zd numbits: %ld\n", r, rsz);
	    gmp_fprintf(pair, "q = %Zd  numbits: %ld\n", q, qsz);
	    fprintf(pair, "rho = %lf\n", (double)qsz/(double)rsz);
	    gmp_fprintf(pair, "t = %Zd\n", t);
	  }
	}
      }
    }
  }
  fclose(pair);
  mpz_clears(r, alpha, x, x0, q, NULL);
}

  
