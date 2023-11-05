/******************************************************************
 *                                                                *
 *   Parse *.output files from base_curve*.c (created using       *
 *   ./base_curve* > file.output). Sort all curves found in       *
 *   decreasing order.                                            *
 *                                                                *
 *****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <ctype.h>
#include <string.h>

/*  find next good line to work with  */

int nxtgd(char *bufr, int n, int end)
{
  int done = 0;
  int k = n;

  while((!done) && (k < end))
  {
    while(bufr[k] != '2') k++;
    k++;
    if(bufr[k] == '^') done = 1;
  }
  return k;
}

/*  convert ' ' to null so string converstion will work  */

void endstr(char *bufr, int n)
{
  int i = n;
  
  while(bufr[i] != '\n') i++;
  bufr[i] = 0;
}

/*  from end of a6, get a6 as string  */

void get_a6(char *bufr, int n, char *a6)
{
  int i;

  while(bufr[n] == ' ') n--;
  i = 4;
  a6[i] = 0;
  i--;
  while(bufr[n] != ' ')
  {
    a6[i] = bufr[n];
    i--;
    n--;
  }
  while(i >= 0)
  {
    a6[i] = ' ';
    i--;
  }
}

/*  find number of * on line  */

int numstar(char *bufr, int n)
{
  int i, j;

  i = n;
  j = 0;
  while(bufr[i] != '\n')
  {
    if(bufr[i] == '*') j++;
    i++;
  }
  return j;
}

int main(int argc, char *argv[])
{
  FILE *otpt, *sv;
  int i, j, k, xp, m, *cof;
  mpz_t *prime, test;
  char *bufr, *ptr, filename[1024];
  char *a6, *strt;
  int  *used, s;
  
/* read in whole file  */
  
  if(argc < 2)
  {
    printf("USE: ./pull_curves <filename>\n");
    exit(-1);
  }
  otpt = fopen(argv[1], "r");
  if(!otpt)
  {
    printf("can't find file %s\n", argv[1]);
    exit(-2);
  }
  k = 0;
  bufr = (char *)malloc(4*1024*1024);
  while(!feof(otpt))
  {
    bufr[k] = fgetc(otpt);
    k++;
  }
  k -= 2;
  fclose(otpt);

  sprintf(filename, "%s", argv[1]);
  ptr = strstr(filename, ".out");
  if(!ptr)
    strcat(filename, ".saved");
  else
    sprintf(ptr, ".saved");
  sv = fopen(filename, "w");

/*  all a4 = 1, so no need to save this.  
    create space based on file for 160  */

  a6 = (char*)malloc(32768*5);   // 5 bytes per entry
  prime = (mpz_t*)malloc(sizeof(mpz_t)*32768);
  cof = (int*)malloc(sizeof(int)*32768);
  used = (int*)malloc(sizeof(int)*32768);
  
/*  look for 2^, pull factors and save to array  */

  i = 0;
  j = 0;
  while(i < k)
  {
    i = nxtgd(bufr, i, k);  // find next curve
    if(i >= k) break;
    s = i - 2;               // end of a6
    get_a6(bufr, s, &a6[j*5]);
    cof[j] = 1;
    i++;
    m = bufr[i] & 0x0f;           // power of 2
    i++;
    while(bufr[i] != ' ')
    {
      m *= 10;
      m += bufr[i] & 0x0f;
      i++;
    }
    if(m)
    {
      while(m)
      {
	cof[j] *= 2;
	m--;
      }
    }
    m = numstar(bufr, i);             //  more cofactors?
    while(m > 1)
    {
      while(!isdigit(bufr[i])) i++;
      xp = bufr[i] & 0x0f;
      i++;
      while(bufr[i] != ' ')
      {
	xp *= 10;
	xp += bufr[i] & 0x0f;
	i++;
      }
      cof[j] *= xp;
      m--;
      while(bufr[i] != '*') i++;
    }
    while(!isdigit(bufr[i])) i++;
    endstr(bufr, i);
    mpz_init_set_str(prime[j], &bufr[i], 10);
    used[j] = 0;
    j++;
  }

/*  find largest primes and output them in decreasing order  */

  mpz_init(test);
  for(i=0; i<j; i++)
  {
    for(k=0; k<j; k++)
    {
      if(used[k]) continue;
      if(mpz_cmp(prime[k], test) > 0)
      {
	mpz_set(test, prime[k]);
	m = k;
      }
    }
    used[m] = 1;
    gmp_fprintf(sv, "1 %s %d %Zx\n", &a6[5*m], cof[m], prime[m]);
    mpz_set_ui(test, 0);
  }
  fclose(sv);
}
