/***********************************************************
 *                                                         *
 *    Create SNARK QAP parameters v, w, y, h and t.        *
 *                                                         *
 **********************************************************/

#include "signature.h"
#include "snarkbase.h"

int main(int argc, char *argv[])
{
  FILE *qap;
  mpz_t list[11];
  int i, j, sw[10], l;
  mpz_t *htable;
  mpz_t *vj, *wj, *yj;
  SIG_SYSTEM sig;

/*  read in system parameters from previously generated file  */

  get_system("curve_11_parameters.bin", &sig);

/* "working in exponent" so use torsion for prime modulus */
  
  minit(sig.tor);

  for(i=0; i<11; i++)
    mpz_init(list[i]);

  mpz_set_ui(list[0], 31);
  mpz_set_ui(list[1], 37);
  mpz_set_ui(list[2], 41);
  mpz_set_ui(list[3], 43);
  mpz_set_ui(list[4], 47);
  mpz_set_ui(list[5], 3);
  mpz_set_ui(list[6], 5);
  mpz_set_ui(list[7], 7);
  mpz_set_ui(list[8], 11);
  mpz_set_ui(list[9], 13);

/*  create left, right and output coefficient 
    matricies for each wire */

  vj = (mpz_t*)malloc(sizeof(mpz_t)*10*5);
  wj = (mpz_t*)malloc(sizeof(mpz_t)*10*5);
  yj = (mpz_t*)malloc(sizeof(mpz_t)*10*5);
  for(i=0; i<50; i++)
  {
    mpz_init(vj[i]);
    mpz_init(wj[i]);
    mpz_init(yj[i]);
  }
  liofx(vj, 0, list, 5);                // v0 = l0
  liofx(&wj[5], 2, list, 5);
  for(j=0; j<5; j++)
    madd(wj[5 + j], wj[5 + j], vj[j]);  // w1 = l0 + l2
  liofx(&vj[10], 1, list, 5);           // v2 = l1
  liofx(&wj[10], 2, list, 5);           // w2 = l2
  liofx(&wj[20], 3, list, 5);           // w4 = l3
  liofx(&wj[15], 1, list, 5);           // w3 = l1 + l3
  for(j=0; j<5; j++)
    madd(wj[15 + j], wj[15 + j], wj[20 + j]); 
  liofx(&vj[25], 2, list, 5);           // v5 = l2
  liofx(&yj[25], 0, list, 5);           // y5 = l0
  liofx(&vj[30], 3, list, 5);           // v6 = l3
  liofx(&yj[30], 1, list, 5);           // y6 = l1
  liofx(&vj[35], 4, list, 5);           // v7 = l4
  liofx(&yj[35], 2, list, 5);           // y7 = l2
  liofx(&wj[40], 4, list, 5);           // w8 = l4
  liofx(&yj[40], 3, list, 5);           // y8 = l3
  liofx(&yj[45], 4, list, 5);           // y9 = l4

/* build up h(x) table  */
  
  htable = all_lilj(list, 5);

/* permutation of indexing to allow statement indexes
   to preceed witness indexesl.  */

  sw[0] = 0;
  sw[1] = 1;
  sw[2] = 2;
  sw[3] = 5;
  sw[4] = 7;
  sw[5] = 3;
  sw[6] = 4;
  sw[7] = 6;
  sw[8] = 8;
  sw[9] = 9;
  l = 4;
  
  qap = fopen("snark.qap", "w");
  i = 5;
  j = 10;
  fwrite(&i, sizeof(int), 1, qap);  // n
  fwrite(&j, sizeof(int), 1, qap);  // m
  fwrite(sw, sizeof(int), 10, qap);
  fwrite(&l, sizeof(int), 1, qap);
  for(i=0; i<5; i++)
    mpz_out_raw(qap, list[i]);
  for(i=0; i<50; i++)
    mpz_out_raw(qap, vj[i]);
  for(i=0; i<50; i++)
    mpz_out_raw(qap, wj[i]);
  for(i=0; i<50; i++)
    mpz_out_raw(qap, yj[i]);
  for(i=0; i<40; i++)
    mpz_out_raw(qap, htable[i]);
  fclose(qap);
/*
  printf("n = %d\n", 5);
  printf("m = 10\n");
  for(i=0; i<10; i++)
    printf("sw[%d] = %d\n", i, sw[i]);
  printf("list:\n");
  for(i=0; i<5; i++)
    gmp_printf("%Zd ", list[5+i]);
  printf("\n");
  printf("vj: \n");
  for(i=0; i<10; i++)
  {
    for(j=0; j<5; j++)
      gmp_printf("%Zd ", vj[i*5 + j]);
    printf("\n");
  }
  printf("wj: \n");
  for(i=0; i<10; i++)
  {
    for(j=0; j<5; j++)
      gmp_printf("%Zd ", wj[i*5 + j]);
    printf("\n");
  }
  printf("yj: \n");
  for(i=0; i<10; i++)
  {
    for(j=0; j<5; j++)
      gmp_printf("%Zd ", yj[i*5 + j]);
    printf("\n");
  }
  printf("h table:\n");
  for(i=0; i<10; i++)
  {
    for(j=0; j<4; j++)
      gmp_printf("%Zd ", htable[i*4 + j]);
    printf("\n");
  }
*/  
  for(i=0; i<11; i++)
    mpz_clear(list[i]);
  for(i=0; i<40; i++)
    mpz_clear(htable[i]);
  free(htable);
  for(i=0; i<50; i++)
  {
    mpz_clear(vj[i]);
    mpz_clear(wj[i]);
    mpz_clear(yj[i]);
  }
  free(vj);
  free(wj);
  free(yj);
}
