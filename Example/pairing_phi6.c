/**************************************************************
 *                                                            *
 *    Look at Algorithm 6.6 in"A Taxonomy of Pairing-Friendly *
 *    Elliptic Curves" D. Freeman, M. Scott, and              *
 *    E Teske (2010)                                          *
 *                                                            *
 *************************************************************/

#include "modulo.h"

/*  table of coeficients for Cyclotomic polynomials of order
    6*prime.  Since only interested in a few fuctions, this is
    a simple way to create polynomials.
*/

static int p5cf[9] = {1, 1, 0, -1, -1, -1, 0, 1, 1};
static int p7cf[13] = {1, 1, 0, -1, -1, 0, 1, 0, -1, -1, 0, 1, 1};
static int p11cf[21] = {1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, -1,
  0, 1, 1, 0, -1, -1, 0, 1, 1};
static int p13cf[25] = {1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0,
  1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1, 1};
static int p17cf[33] = {1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0,
  1, 1, 0, -1, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1, 1};
static int p19cf[37] = {1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0,
  1, 1, 0, -1, -1, 0, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1, 1,
  0, -1, -1, 0, 1, 1};
static int p23cf[45] = {1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0,
  1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, -1, 0, 1, 1, 0, -1, -1,
  0, 1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1, 1};
static int p29cf[57] = {1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0,
  1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, -1,
  0, 1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1,
  0, 1, 1, 0, -1, -1, 0, 1, 1};
static int p31cf[61] = {1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0,
  1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0,
  1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1, 1,
  0, -1, -1, 0, 1, 1, 0, -1, -1, 0, 1, 1};

/*  create pointer and get number of coeficients from input value.
    returns -1 count and NULL pointer if k is a bad value.
*/

int phicoef(int k, int *coef[])
{
  int j;

  j = 2*k - 1;
  *coef = NULL;
  if((k < 5) || (k > 31))
    return -1;
  switch(k)
  {
    case 5:
      *coef = p5cf;
      break;
    case 7:
      *coef = p7cf;
      break;
    case 11:
      *coef = p11cf;
      break;
    case 13:
      *coef = p13cf;
      break;
    case 17:
      *coef = p17cf;
      break;
    case 19:
      *coef = p19cf;
      break;
    case 23:
      *coef = p23cf;
      break;
    case 29:
      *coef = p29cf;
      break;
    case 31:
      *coef = p31cf;
      break;
    default:
      j = -1;
      break;
  }
  return j;
}

/* compute Phi_6k(x) for input k and x.
   k must be prime in range 5 to 31.
   returns 1 on success, 0 on failure.
   r is computed response on success.
*/

int phi6k(mpz_t r, long k, mpz_t x)
{
  int i;
  int *cf, len;
  mpz_t z;

  len = phicoef(k, &cf);
  if(len < 0)
    return 0;
  mpz_init_set_ui(z, 1);
  mpz_set_ui(r, 0);
  for(i=0; i<len; i++)
  {
    if(cf[i] > 0)
      mpz_add(r, r, z);
    else if(cf[i] < 0)
      mpz_sub(r, r, z);
    mpz_mul(z, z, x);
  }
  mpz_clear(z);
  return 1;
}

/*  compute t(x) for k ~ 1 mod 6  */

void tofx1(mpz_t t, long k, mpz_t x)
{
  mpz_t xk;

  mpz_init(xk);
  mpz_pow_ui(xk, x, k+1);
  mpz_set_ui(t, 1);
  mpz_add(t, t, x);
  mpz_sub(t, t, xk);
  mpz_clear(xk);
}

/*  compute t(x) for k ~ 5 mod 6  */

void tofx5(mpz_t t, long k, mpz_t x)
{
  mpz_t x3;

  mpz_init(x3);
  mpz_pow_ui(x3, x, 3);
  mpz_add_ui(t, x3, 1);
  mpz_clear(x3);
}

/*  compute q(x) for k ~ 1 mod 6  */

void qofx1(mpz_t q, long k, mpz_t x)
{
  mpz_t xp, t1, t2, t3;

  mpz_inits(xp, t1, t2, t3, NULL);
  mpz_pow_ui(xp, x, k);
  mpz_mul(t1, xp, xp);         // x^2k
  mpz_sub(t1, t1, xp);         // + x^k 
  mpz_add_ui(t1, t1, 1);       // + 1
  mpz_add_ui(t2, x, 1);
  mpz_mul(t2, t2, t2);         // (x + 1)^2
  mpz_pow_ui(t3, x, 2*k+1);    // x^(2k+1)
  mpz_mul(t1, t1, t2);
  if(mpz_divisible_ui_p(t1, 3))
    mpz_divexact_ui(t1, t1, 3);
  else
  {
    mpz_set_ui(t2, 3);
    mpz_fdiv_q(t1, t1, t2);
  }
  mpz_sub(q, t1, t3);
  mpz_clears(xp, t1, t2, t3, NULL);
}

/*  compute q(x) for k ~ 5 mod 6  */

void qofx5(mpz_t q, long k, mpz_t x)
{
  mpz_t xp, t1, t2, t3;

  mpz_inits(xp, t1, t2, t3, NULL);
  mpz_pow_ui(xp, x, k);
  mpz_mul(t1, xp, xp);         // x^2k
  mpz_sub(t1, t1, xp);         // + x^k 
  mpz_add_ui(t1, t1, 1);       // + 1
  mpz_mul(t2, x, x);           // x^2
  mpz_sub(t2, t2, x);          // -x
  mpz_add_ui(t2, t2, 1);       // + 1
  mpz_mul(t3, xp, x);          // x^(k+1)
  mpz_mul(t1, t1, t2);
  if(mpz_divisible_ui_p(t1, 3))
    mpz_divexact_ui(t1, t1, 3);
  else
  {
    mpz_set_ui(t2, 3);
    mpz_fdiv_q(t1, t1, t2);
  }
  mpz_add(q, t1, t3);
  mpz_clears(xp, t1, t2, t3, NULL);
}

int main(int argc, char *argv[])
{
  FILE *pair;
  int i, max;
  long k, rpm, qpm, rsz, qsz;
  mpz_t r, x, t, q;
  char filename[128];

  if(argc < 3)
  {
    printf("use: ./pairing_phi6 <embedding degree> <search limit>\n");
    printf("primes in 5 to 31 allowed for embedding degree\n");
    exit(-1);
  }

  k = atol(argv[1]);
  if((k < 5) || (k > 31))
  {
    printf("bad embedding degree.\n");
    exit(-2);
  }
  
  mpz_inits(r, x, t, q, NULL);
  mpz_set_ui(x, 1);
  if(!phi6k(r, k, x))
  {
    printf("bad embedding degree!\n");
    exit(-3);
  }
  max = atoi(argv[2]);
  if(max < 2)
  {
    printf("bad limit!!\n");
    exit(-4);
  }
  sprintf(filename, "pairing_phi6.%02ld", k);
  pair = fopen(filename, "w");
  for(i=1; i<=max; i++)
  {
    mpz_set_ui(x, i);
    phi6k(r, k, x);
    rpm = mpz_probab_prime_p(r, 25);
    rsz = mpz_sizeinbase(r, 2);
    if(rpm)
    {
      if(k%6 == 1)
      {
	qofx1(q, k, x);
	tofx1(t, k, x);
      }
      else    // k%6 == 5
      {
	qofx5(q, k, x);
	tofx5(t, k, x);
      }
      qpm = mpz_probab_prime_p(q, 25);
      qsz = mpz_sizeinbase(q, 2);
      if(qpm)
      {
	gmp_fprintf(pair, "k= %ld  x = %Zd\n", k, x);
	gmp_fprintf(pair, "r = %Zd numbits: %ld\n", r, rsz);
	gmp_fprintf(pair, "q = %Zd  numbits: %ld\n", q, qsz);
	fprintf(pair, "rho = %lf\n", (double)qsz/(double)rsz);
	gmp_fprintf(pair, "t = %Zd\n\n", t);
      }
    }
  }
  fclose(pair);
  mpz_clears(r, x, t, q, NULL);
}
