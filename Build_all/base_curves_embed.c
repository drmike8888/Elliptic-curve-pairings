/**********************************************************************
 *                                                                    *
 *    Find out if curves with good primes have low embedding degree.  *
 *    Make "low" as large as possible (> 1000?)                       *
 *                                                                    *
 *********************************************************************/

#include "modulo.h"
#include <string.h>

int dgt(char *fld)
{
  int i, k;

  k = 0;
  for(i=0; i<3; i++)
  {
    k *= 10;
    k += fld[i] & 0x0f;
  }
  return k;
}
  
int main(int argc, char *argv[])
{
  FILE *csv, *dat;
  int i, j, k, fnd;
  long a4, r;
  char filename[1024], *fld, dummy[256], a6[8];
  mpz_t prm, p, t, tpw;

  if(argc < 2)
  {
    printf("Use: ./base_curves_embed <filename>\n");
    printf("     where <filename> is *_full_sort.csv\n");
    exit(-1);
  }

  csv = fopen(argv[1], "r");
  if(!csv)
  {
    printf("can't find file %s\n", argv[1]);
    exit(-2);
  }

  sprintf(filename, "%s", argv[1]);
  fld = strstr(filename, "_full_sort.csv");
  if(!fld)
  {
    printf("wrong kind of file - must be fully processed\n");
    exit(-3);
  }
  fld -= 3;
  k = dgt(fld);
  switch(k)
  {
    case 160:
      mpz_init_set_str(prm, "ac000000000000000000000000000000000000001", 16);
      break;
    case 256:
      mpz_init_set_str(prm, "2b000000000000000000000000000000000000000000000000000000000000001", 16);
      break;
    case 384:
      mpz_init_set_str(prm, "2e00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", 16);
      break;
    case 512:
      mpz_init_set_str(prm, "e20000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", 16);
      break;
    default:
      printf("don't have prime for field size %d\n", k);
      exit(-4);
  }
  fld += 3;
  sprintf(fld, "%s", "_embed.dat");
  dat = fopen(filename, "w");
  mpz_inits(p, t, tpw, NULL);
  
/*  for each line, compute t-1 and then see if (t-1)^k = 1 mod p  */

  while(!feof(csv))
  {
    gmp_fscanf(csv, "%ld %s %ld %Zd %s", &a4, a6, &r, p, dummy);
    minit(p);
    mpz_mul_si(t, p, r);
    mpz_sub(t, prm, t);         // now have t - 1
    mpz_set(tpw, t);
    fnd = 0;
    k = 2;
    while((!fnd) && (k < 257))
    {
      mmul(tpw, tpw, t);
      j = mpz_cmp_ui(tpw, 1);
      if(!j)
	fnd = 1;
      else
	k++;
    }
    gmp_fprintf(dat, "%ld %s %ld %Zd", a4, a6, r, p);
    if(fnd)
      fprintf(dat, "  %d\n", k);
    else
      fprintf(dat, " > 256\n");
  }
  fclose(dat);
  fclose(csv);
}

