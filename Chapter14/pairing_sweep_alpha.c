/*****************************************************************
 *                                                               *
 *     Create pairing friendly curves using algorithms 6.19,     *
 *     6.20, and 6.24 in "A Taxonomy of Pairing-Friendly         *
 *     Elliptic Curves" D. Freeman, M. Scott, and E Teske (2010) *
 *                                                               *
 ****************************************************************/

#include "modulo.h"
#include "poly.h"
#include <math.h>

/*  Compute cyclotomic polynomial Phi_4k(z).
    Enter with k, alpha and x.
    z = alpha*x^2  
*/
void phi4k(mpz_t r, long k, mpz_t alpha, mpz_t x)
{
  int i;
  mpz_t z[36], ck;

  mpz_init_set_ui(ck, 1);
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
}

/*  compute q(z) using equation 6.20 and 6.19 from
    Taxonomy paper.
    Enter with k, alpha and x. 
    If these parameters give 4q exactly, 
    Returns q (which should have already been initialized)
*/

void qofz_20(mpz_t q, long k, mpz_t alpha, mpz_t x)
{
  long k1, k2;
  mpz_t  t1, t2, t3, t4, z;

  mpz_inits(z, t1, t2, t3, t4, NULL);
  k1 = k + 1;
  k2 = k1/2;
  mpz_mul(z, x, x);
  mpz_mul(z, z, alpha);
  mpz_pow_ui(t1, z, k1);
  mpz_pow_ui(t2, z, k);
  mpz_pow_ui(t3, z, k2);
  mpz_mul_ui(t3, t3, 4);
  mpz_add_ui(t4, z, 1);
  mpz_add(t1, t1, t2);
  mpz_add(t1, t1, t3);
  mpz_add(t1, t1, t4);
  mpz_fdiv_q_ui(q, t1, 4);
  mpz_clears(z, t1, t2, t3, t4, NULL);
}

/*  compute q(z) using equation 6.2 and 6.19 from
    Taxonomy paper.
    Enter with k, alpha and x. 
    If these parameters give 4q exactly, 
    Returns q (which should have already been initialized)
*/

void qofz_2(mpz_t q, long k, mpz_t alpha, mpz_t x)
{
  long k1, k2;
  mpz_t  t1, t2, t3, t4, z;

  mpz_inits(z, t1, t2, t3, t4, NULL);
  k1 = k + 1;
  k2 = k + 2;
  mpz_mul(z, x, x);
  mpz_mul(z, z, alpha);
  mpz_pow_ui(t1, z, k2);
  mpz_pow_ui(t2, z, k1);
  mpz_mul_ui(t2, t2, 2);
  mpz_pow_ui(t3, z, k);
  mpz_sub_ui(t4, z, 1);
  mpz_mul(t4, t4, t4);
  mpz_add(t1, t1, t2);
  mpz_add(t1, t1, t3);
  mpz_add(t1, t1, t4);
  mpz_fdiv_q_ui(q, t1, 4);
  mpz_clears(z, t1, t2, t3, t4, NULL);
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

void tofz_20(mpz_t t, long k, mpz_t alpha, mpz_t x)
{
  mpz_t z, t1;
  long k1;

  mpz_inits(z, t1, NULL);
  mpz_mul(z, x, x);
  mpz_mul(z, z, alpha);
  k1 = (k + 1)/2;
  mpz_pow_ui(t1, z, k1);
  mpz_add_ui(t, t1, 1);
  mpz_clears(z, t1, NULL);
}

/*  compute t(z) using equation 6.2 and 6.19 from
    Taxonomy paper.  (trace of Frobenius)
    Enter with k, alpha and x. 
*/

void tofz_2(mpz_t t, long k, mpz_t alpha, mpz_t x)
{
  mpz_t z, one;

  mpz_init(z);
  mpz_mul(z, x, x);
  mpz_mul(z, z, alpha);
  mpz_init_set_ui(one, 1);
  mpz_sub(t, one, z);
  mpz_clears(z, one, NULL);
}

int main(int argc, char *argv[])
{
  FILE *pair;
  mpz_t r, alpha, x, q, t;
  long k, u, lg2r, alphabase, rpm, qpm, rsz, qsz;
  long max, atab[5] = {7, 11, 15, 19, 23};
  double rho;
  int j, m, nmrpm, a, kdex;
  int ktab[9] = {5, 7, 11, 13, 17, 19, 23, 29, 31};
  int algt[9] = {1, 0,  0,  1,  1,  0,  0,  1,  0}; // 1 == 6.2  0 == 6.20
  double xsz, asz;
  char filename[128];

  if(argc < 2)
  {
    printf("Use: ./pairing_sweep_alpha  <log2(r) max range>\n");
    exit(-1);
  }
  lg2r = atol(argv[1]);
  if(lg2r < 3)
  {
    printf("need more bits to work with\n");
    exit(-2);
  }
  printf("choose embedding degree k from 5, 7, 11, 13, 17, 19, 23, 29, 31: ");
  fflush(stdout);
  scanf("%ld", &k);
  for(kdex=0; kdex<9; kdex++)
    if(k == ktab[kdex])
      break;
  if(kdex > 8)
  {
    printf("k must be from list.\n");
    exit(-1);
  }
  sprintf(filename, "pairings_alpha.%02ld", k);
  pair = fopen(filename, "w");
  nmrpm = 0;
  mpz_inits(r, alpha, x, q, t, NULL);
  max = exp2(lg2r/2/(k - 1));
  for(a=0; a<8; a++)
  {
    for(m=0; m<5; m++)
    {
      alphabase = atab[m] + a*20;
      printf("processing alpha = %ld\n", alphabase);
      mpz_set_ui(alpha, alphabase);
      for(j=1; j<max; j+=2)
      {
	mpz_set_ui(x, j);
//      gmp_printf("k= %ld alpha = %Zd  x = %Zd\n", k, alpha, x);
	phi4k(r, k, alpha, x);
	rsz = mpz_sizeinbase(r, 2);
	if(rsz > lg2r)
	  break;
	rpm = mpz_probab_prime_p(r, 25);
	if(rpm)
        {
	  nmrpm++;
	  if(algt[kdex])
	  {
	    qofz_2(q, k, alpha, x);
	    tofz_2(t, k, alpha, x);
	  }
	  else
	  {
	    qofz_20(q, k, alpha, x);
	    tofz_20(t, k, alpha, x);
	  }
	  qsz = mpz_sizeinbase(q, 2);
//	  gmp_printf("q: %Zd  bits: %ld\n", q, qsz);
	  qpm = mpz_probab_prime_p(q, 25);
	  if(qpm)
	  {
	    gmp_fprintf(pair, "k= %ld alpha = %Zd  x = %Zd\n", k, alpha, x);
	    gmp_fprintf(pair, "r = %Zd numbits: %ld\n", r, rsz);
	    gmp_fprintf(pair, "q = %Zd  numbits: %ld\n", q, qsz);
	    fprintf(pair, "rho = %lf\n", (double)qsz/(double)rsz);
	    gmp_fprintf(pair, "t = %Zd\n\n", t);
	  }
	}
      }
    }
  }
  fclose(pair);
  mpz_clears(r, alpha, x, q, t, NULL);
  printf("found %d r primes\n", nmrpm);
}

  
