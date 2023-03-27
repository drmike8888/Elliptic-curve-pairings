/**********************************************************************
 *                                                                    *
 *    Find out if curves with good primes have low embedding degree.  *
 *    Make "low" as large as possible (> 1000?)                       *
 *                                                                    *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>

int dgt(char *fld)
{
  int i, k;

  k = 0;
  for(i=0; i<3; i++)
  {
    k += fld[i] & 0x0f;
    k *= 10;
  }
  return k;
}
  
int main(int argc, char *argv[])
{
  FILE *csv, *dat;
  int i, j, k;
  long a4, r;
  char filename[1024], *fld, dummy[256], a6[8];
  mpz_t prm, p, k;

  if(argc < 2)
  {
    printf("Use: ./base_curves_embed <filename>\n");
    printf("     where <filename> is *_full_sort.csv\n");
    exit(-1);
  }

  csv = fopen(argv[1]);
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
      mpz_init_set_str(prm, "0xac000000000000000000000000000000000000001", 16);
      break;
    case 256:
      mpz_init_set_str(prm, "0x2b000000000000000000000000000000000000000000000000000000000000001", 16);
      break;
    case 384:
      mpz_init_set_str(prm, "0x2e00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", 16);
      break;
    case 512:
      mpz_init_set_str(prm, "0xe20000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", 16);
      break;
    defalult:
      printf("don't have prime for field size %d\n", k);
      exit(-4);
  }
  fld += 3;
  sprintf(fld, "%s", "_embed.dat");
  dat = fopen(filename, "w");
  mpz_init(p);
  while(!feof(csv))
  {
    gmp_fscanf(csv, "%ld %s %ld %Zd %s", &a4, a6, &r, p);
    
