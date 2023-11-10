/***********************************************************
 *                                                         *
 *    Verify one record of SNARK using CRS.                *
 *                                                         *
 **********************************************************/

#include "signature.h"
#include "snarkbase.h"

int main(int argc, char *argv[])
{
  FILE *crs, *qap;
  int i, j, k, m, n, l;
  mpz_t *aj;
  SIG_SYSTEM sig;
  POINT aG, bG, dG, *zG, *thtaG, *ztG;
  POLY_POINT bH, gH, dH, *zH;
  POINT A, C, Tmp, V;
  POLY_POINT B, P, R;
  POLY eab, eVH, eCH, eAB;

/*  read in system parameters from previously generated file  */

  get_system("curve_11_parameters.bin", &sig);
  minit(sig.prime);
  poly_irrd_set(sig.irrd);
  poly_mulprep(sig.irrd);

/* read in CRS values next  */
  
  crs = fopen("snark.crs", "r");
  fread(&n, sizeof(int), 1, crs);
  fread(&m, sizeof(int), 1, crs);
  fread(&l, sizeof(int), 1, crs);
  point_init(&aG);
  point_init(&bG);
  point_init(&dG);
  point_read(&aG, crs);
  point_read(&bG, crs);
  point_read(&dG, crs);
  zG = (POINT*)malloc(sizeof(POINT)*n);
  for(i=0; i<n; i++)
  {
    point_init(&zG[i]);
    point_read(&zG[i], crs);
  }
  thtaG = (POINT*)malloc(sizeof(POINT)*m);
  for(i=0; i<m; i++)
  {
    point_init(&thtaG[i]);
    point_read(&thtaG[i], crs);
  }
  j = n - 1;
  ztG = (POINT*)malloc(sizeof(POINT)*j);
  for(i=0; i<j; i++)
  {
    point_init(&ztG[i]);
    point_read(&ztG[i], crs);
  }
  poly_point_init(&bH);
  poly_point_read(&bH, crs);
  poly_point_init(&dH);
  poly_point_read(&dH, crs);
  poly_point_init(&gH);
  poly_point_read(&gH, crs);
  zH = (POLY_POINT*)malloc(sizeof(POLY_POINT)*n);
  for(i=0; i<n; i++)
  {
    poly_point_init(&zH[i]);
    poly_point_read(&zH[i], crs);
  }
  fclose(crs);
/*
  printf("\nreading crs file\n");
  printf("n = %d m = %d l = %d\n", n, m, l);
  point_printf("aG: ", aG);
  point_printf("bG: ", bG);
  point_printf("dG: ", dG);
  for(i=0; i<n; i++)
  {
    printf("zG[%d]: ", i);
    point_printf(" ", zG[i]);
  }
  for(i=0; i<m; i++)
  {
    printf("thtaG[%d]:", i);
    point_printf(" ", thtaG[i]);
  }
  poly_point_printf("bH: ", bH);
  poly_point_printf("dH: ", dH);
  poly_point_printf("gH: ", gH);
  for(i=0; i<n; i++)
  {
    printf("zH[%d]", i);
    poly_point_printf(" ", zH[i]);
  }
*/
/* read in record data with proof values */

  qap = fopen("snark_record.0", "r");
  fread(&n, sizeof(int), 1, qap);
  fread(&m, sizeof(int), 1, qap);
  fread(&l, sizeof(int), 1, qap);
  aj = (mpz_t*)malloc(sizeof(mpz_t)*(l+1));
  for(i=0; i<=l; i++)
  {
    mpz_init(aj[i]);
    mpz_inp_raw(aj[i], qap);
  }
  point_init(&A);
  point_read(&A, qap);
  point_init(&C);
  point_read(&C, qap);
  poly_point_init(&B);
  poly_point_read(&B, qap);
  fclose(qap);
  printf("reading in record 0\n");
  printf("n = %d m = %d l = %d\n", n, m, l);
  for(i=0; i<=l; i++)
    gmp_printf("aj[%d]: %Zd\n", i, aj[i]);
  point_printf("A: ", A);
  point_printf("C: ", C);
  poly_point_printf("B: ", B);
  
  /*  COMPUTE VERIFIER VALUES  */
/* statement values times groth points */

  point_init(&V);
  for(i=0; i<=l; i++)
  {
    elptic_mul(&Tmp, thtaG[i], aj[i], sig.E);
    printf("aj[%d]:", i);
    point_printf(" ", Tmp);
    elptic_sum(&V, V, Tmp, sig.E);
    point_printf("V: ", V);
  }
  
  poly_point_init(&R);
  poly_point_rand(&R, sig.Ex);
  poly_init(&eab);
  poly_point_init(&P);
  tog2(&P, aG);
  tate(&eab, P, bH, R, sig.tor, sig.Ex);
  poly_printf("e(a*G, b*H): ", eab);
  
  poly_init(&eVH);
  tog2(&P, V);
  tate(&eVH, P, gH, R, sig.tor, sig.Ex);
  poly_printf("e(V, gH): ", eVH);
  poly_init(&eCH);
  tog2(&P, C);
  tate(&eCH, P, dH, R, sig.tor, sig.Ex);
  poly_printf("e(C, dH): ", eCH);
  poly_init(&eAB);
  tog2(&P, A);
  tate(&eAB, P, B, R, sig.tor, sig.Ex);
  poly_printf("e(A, B): ", eAB);
  poly_mul(&eab, eab, eVH);
  poly_mul(&eab, eab, eCH);
  poly_printf("e()*e()*e(): ", eab);

  if(poly_cmp(eab, eAB))
    printf("Record verifies!\n");
  else
    printf("Record falsified!\n");
  
  for(i=0; i<=l; i++)
    mpz_clear(aj[i]);
  free(aj);
  point_clear(&aG);
  point_clear(&bG);
  point_clear(&dG);

  for(i=0; i<n; i++)
    point_clear(&zG[i]);
  for(i=0; i<m; i++)
    point_clear(&thtaG[i]);
  poly_point_clear(&bH);
  poly_point_clear(&gH);
  poly_point_clear(&dH);
  for(i=0; i<n; i++)
    poly_point_clear(&zH[i]);
  for(i=0; i<n-1; i++)
    point_clear(&ztG[i]);

  point_clear(&A);
  point_clear(&Tmp);
  poly_point_clear(&B);
  poly_point_clear(&P);
  point_clear(&C);
  point_clear(&V);
  poly_clear(&eab);
  poly_clear(&eVH);
  poly_clear(&eCH);
  poly_clear(&eAB);
  poly_point_clear(&R);
}
