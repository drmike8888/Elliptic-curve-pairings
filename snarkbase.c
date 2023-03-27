/************************************************************
 *                                                          *
 *   Create simple specific Quadratic Arithmetic Program    *
 *   to build for example SNARK.                            *
 *                                                          *
 ***********************************************************/

#include "signature.h"

/* compute p_i = (r_i - r_j)(r_i - r_k)...
   Enter with list of primes for each r_i, number of them
   and index for which one is positive.
*/

void p_i(mpz_t p, int i, mpz_t *list, int n)
{
  int j;
  mpz_t term, rslt;

  mpz_inits(term, rslt, NULL);
  mpz_set_ui(rslt, 1);
  for(j=0; j<n; j++)
  {
    if(j == i)
      continue;
    msub(term, list[i], list[j]);
    mmul(rslt, rslt, term);
  }
  mpz_set(p, rslt);
  mpz_clears(term, rslt, NULL);
}

/* create list of start and end indexes 
   for given number of rows
   array k should have room for up to 2*n ints.
*/

void startk(int *k, int rows, int limt)
{
  int m, index;

  index = 0;
  for(m=0; m<rows; m++)
  {
    k[2*m] = index;
    index++;
  }
  index = limt - 1;
  for(m=rows-1; m>=0; m--)
  {
    k[2*m + 1] = index;
    index--;
  }
}

/* Compute Lagrange interpolant coefficients. If one function,
   then input i = j. Routine returns n-1 coefficients for li(x)
   x^(n-1), ... ,x^0.

   if i != j, 
   Compute product of two Lagrange interpolants divided by t(x)
   Enter with list of primes and number of gates and
   the indexes i and j for the cross product.
   Returns array of coefficients from x^(n-2) ... x^0.
   Assumption is array is pre-initialized.
*/

void li_lj(mpz_t *coef, int i, int j, mpz_t *list, int n)
{
  int *k, m, cfdex, done, bmp, r, cflimt;
  mpz_t sum, term, *sublst;

  if(i == j)
    cflimt = n - 1;
  else
    cflimt = n - 2;
  sublst = (mpz_t*)malloc(sizeof(mpz_t)*cflimt);
  r = 0;
  for(m=0; m<n; m++)
  {
    if((m != i) && (m != j))
    {
      mpz_init(sublst[r]);
      mpz_neg(sublst[r], list[m]);
      r++;
    }
  }
  k = (int*)malloc(sizeof(int)*n*2);
  mpz_set_ui(coef[0], 1);
  mpz_inits(sum, term, NULL);
  for(cfdex = 1; cfdex < cflimt+1; cfdex++)
  {
    mpz_set_ui(sum, 0);
    startk(k, cfdex, cflimt);
    done = 0;
    while(!done)
    {
      mpz_set(term, sublst[k[0]]);
      for(m=1; m<cfdex; m++)
	mmul(term, term, sublst[k[2*m]]);
      madd(sum, sum, term);
      m = cfdex - 1;
      bmp = 0;
      while(!bmp && !done)
      {
	if(k[2*m] != k[2*m + 1])
	{
	  bmp = 1;
	  k[2*m]++;
	}
	else
	{
	  while(!bmp)
	  {
	    m--;
	    if(m < 0)
	    {
	      done = 1;
	      bmp = 2;
	    }
	    else
	    {
	      if(k[2*m] != k[2*m + 1])
	      {
		k[2*m]++;
		while(m < cfdex-1)
		{
		  m++;
		  k[2*m] = k[2*(m-1)] + 1;
		}
		bmp = 1;
	      }
	    }
	  }
	}
      }
    }
    mpz_set(coef[cfdex], sum);
  }
  mpz_clears(sum, term, NULL);
  for(m=0; m<cflimt; m++)
    mpz_clear(sublst[m]);
}

/* Compute coefficients for li(x).
   input list of primes for gates, number of gates
   and index i. Returns n coefficients for degree n-1
   polynomials in order x^(n-1), ..., x^0
*/

void liofx(mpz_t *coef, int i, mpz_t *list, int n)
{
  mpz_t pi;
  int j;
  
  li_lj(coef, i, i, list, n);
  mpz_init(pi);
  p_i(pi, i, list, n);
  minv(pi, pi);
  for(j=0; j<n; j++)
    mmul(coef[j], coef[j], pi);
  mpz_clear(pi);
}

/* Compute coefficients for li(x)lj(x)/t(x).
   input list of primes for gates, number of gates
   and indecies i and j. Returns n-1 coefficients for
   degree n-2 polynomial in order x^(n-2), ..., x^0
*/

void liljofx(mpz_t *coef, int i, int j, mpz_t *list, int n)
{
  mpz_t pi, pj;
  int k;

  li_lj(coef, i, j, list, n);
  mpz_inits(pi, pj, NULL);
  p_i(pi, i, list, n);
  p_i(pj, j, list, n);
  mmul(pi, pi, pj);
  minv(pi, pi);
  for(k=0; k<n-1; k++)
    mmul(coef[k], coef[k], pi);
  mpz_clears(pi, pj, NULL);
}

/*  Create table of li_lj coefficients for all combinations
    of cross terms. Enter with list of primes for gates and
    number of gates n. Creates space for all coefficients,
    then initializes them. Total number of elements to free
    will be (n-1)^2 * n / 2. 
*/

mpz_t* all_lilj(mpz_t *list, int n)
{
  mpz_t *table;
  int i, j, k;

  k = n*(n-1)*(n-1)/2;
  table = (mpz_t*)malloc(sizeof(mpz_t)*k);
  for(i=0; i<k; i++)
    mpz_init(table[i]);
  k = 0;
  for(i=0; i<n-1; i++)
  {
    for(j=i+1; j<n; j++)
    {
      liljofx(&table[k], i, j, list, n);
      k += n - 1;
    }
  }
  return table;
}

/* use coefficients of Lagrange polynomials to compute value.
   Enter with z value, pointer to coefficients l(z) and degree
   of coefficients. Assumes result spot is pre-initialized.
   coefficients are in decreasing power: 
       c_0*x^deg + c_1*x^(deg-1) +...+ c_0
*/

void lcalc(mpz_t rslt, mpz_t z, mpz_t *coef, int deg)
{
  int i;

  mpz_set(rslt, coef[0]);
  for(i=1; i<=deg; i++)
  {
    mmul(rslt, rslt, z);
    madd(rslt, rslt, coef[i]);
  }
}

/*  compute t(z) for given z and list of n primes for gates. */

void tofzgrth(mpz_t t, mpz_t z, mpz_t *list, int n)
{
  int i;
  mpz_t tmp;

  mpz_init(tmp);
  mpz_set_ui(t, 1);
  for(i=0; i<n; i++)
  {
    msub(tmp, z, list[i]);
    mmul(t, t, tmp);
  }
}

/*  multiply matrix by column vector and sum all columns.
    Not a standard matrix operation, but required for
    combining CRS values. Matrix is length rows and width
    columns.
*/

void matflat(mpz_t *vector, mpz_t *mat, int width, mpz_t *coef, int length)
{
  int i, j;
  mpz_t *cmt;

  cmt = (mpz_t*)malloc(sizeof(mpz_t)*width*length);
  for(i=0; i<length; i++)
  {
    for(j=0; j<width; j++)
    {
      mpz_init(cmt[i*width + j]);
      mmul(cmt[i*width + j], mat[i*width + j], coef[i]);
    }
  }
  for(j=0; j<width; j++)
  {
    mpz_set_ui(vector[j], 0);
    for(i=0; i<length; i++)
      madd(vector[j], vector[j], cmt[i*width + j]);
  }
  for(i=0; i<length*width; i++)
    mpz_clear(cmt[i]);
}

