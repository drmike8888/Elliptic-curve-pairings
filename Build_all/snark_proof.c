/***********************************************************
 *                                                         *
 *    Create proof for one record of data using CRS and    *
 *    QAP data.                                            *
 *                                                         *
 **********************************************************/

#include "signature.h"
#include "snarkbase.h"

/*  create combined coefficients for specific example problem
    with 5 gates and 4 inputs. Enter with values for
    a_0 through a_9 where a_0 = 1, a_1, a_2, a_5 and a_7 are 
    public (statement) and a_3, a_4, a_6, a_8 and a_9 are private (witness).
    The 10 l_i*l_j coefficient result space should already be initialized.
*/

void crossterms(mpz_t *crscoef, mpz_t *a)
{
  mpz_t t0, t1, t2, t3, t4, t5;

  mpz_inits(t0, t1, t2, t3, t4, t5, NULL);

  // l0-l1
  mmul(t1, a[1], a[2]);
  mmul(t2, a[0], a[3]);
  madd(crscoef[0], t1, t2);
  // l0-l2
  madd(t0, a[1], a[2]);
  mmul(t1, t0, a[0]);
  mmul(t2, a[1], a[5]);
  madd(crscoef[1], t1, t2);
  // l0-l3
  madd(t1, a[3], a[4]);
  mmul(t2, t1, a[0]);
  mmul(t3, a[1], a[6]);
  madd(crscoef[2], t2, t3);
  // l0-l4
  mmul(t4, a[0], a[8]);
  mmul(t5, a[1], a[7]);
  madd(crscoef[3], t4, t5);
  // l1-l2
  mmul(t2, t0, a[2]);
  mmul(t3, a[3], a[5]);
  madd(crscoef[4], t2, t3);
  // l1-l3
  mmul(t2, t1, a[2]);
  mmul(t3, a[3], a[6]);
  madd(crscoef[5], t2, t3);
  // l1-l4
  mmul(t4, a[3], a[7]);
  mmul(t5, a[2], a[8]);
  madd(crscoef[6], t4, t5);
  // l2-l3
  mmul(t3, t0, a[6]);
  mmul(t4, t1, a[5]);
  madd(crscoef[7], t3, t4);
  // l2-l4
  mmul(t2, t0, a[7]);
  mmul(t3, a[5], a[8]);
  madd(crscoef[8], t2, t3);
  // l3-l4
  mmul(t2, t1, a[7]);
  mmul(t3, a[6], a[8]);
  madd(crscoef[9], t2, t3);
  mpz_clears(t0, t1, t2, t3, t4, t5, NULL);
}

/*  create values for all a_j for wires in QAP. 
    Enter with a2 = dose, a3 = patient and a1 = medicine numbers.
    picks random a3 and computes all other wires by form of QAP.
*/

void wires(mpz_t *aj, mpz_t a1, mpz_t a2, mpz_t a3)
{
  mpz_t tmp;

  mpz_set(aj[1], a1);
  mpz_set(aj[2], a2);
  mpz_set(aj[3], a3);
  mpz_set_ui(aj[0], 1);
  mpz_set(aj[5], a1);
  mrand(aj[4]);
  mmul(aj[6], a2, a3);
  mpz_init_set(tmp, a1);
  madd(tmp, tmp, a2);
  mmul(aj[7], tmp, a1);
  mpz_set(tmp, a3);
  madd(tmp, tmp, aj[4]);
  mmul(aj[8], tmp, aj[6]);
  mmul(aj[9], aj[8], aj[7]);
  mpz_clear(tmp);
}

/* create random values in range of torsion order.
   Reads in file from random.org and calls hash
   function to return the two values. Extremely 
   specialized for this example.
*/

void rand_rs_get(mpz_t r, mpz_t s, mpz_t tor)
{
  FILE *rnd;
  char *bufr;
  int i;

  rnd = fopen("random.rs_data", "r");
  if(!rnd)
  {
    printf("can't find file random.rs_data\n");
    exit(-3);
  }
  bufr = (char*)malloc(256);
  for(i=0; i<128; i++)
    fscanf(rnd, "%hhd", &bufr[i]);
  fclose(rnd);
  hash1(r, bufr, 48, tor);
  hash1(s, &bufr[64], 48, tor);
}

int main(int argc, char *argv[])
{
  FILE *qap, *crs;
  mpz_t *list, *coef;
  int i, j, k, m, n, l, *sw;
  mpz_t *htable, *hxcoef;
  mpz_t *vj, *wj, *yj, *aj;
  SIG_SYSTEM sig;
  POINT aG, bG, dG, *zG, *thtaG, *ztG;
  POLY_POINT bH, gH, dH, *zH;
  mpz_t *av, *aw, r, s;
  POINT A, C, Tmp, V;
  POLY_POINT B, P, R;

/*  read in system parameters from previously generated file  */

  get_system("curve_11_parameters.bin", &sig);
  minit(sig.prime);
  poly_irrd_set(sig.irrd);
  poly_mulprep(sig.irrd);

/*  read in QAP parameters */

  printf("open snark.qap file\n");
  qap = fopen("snark.qap", "r");
  if(!qap)
  {
    printf("can't find file snark.qap\n");
    exit(-2);
  }
  fread(&n, sizeof(int), 1, qap);
  fread(&m, sizeof(int), 1, qap);
  sw = (int*)malloc(sizeof(int)*m);
  fread(sw, sizeof(int), m, qap);
  list = (mpz_t*)malloc(sizeof(mpz_t)*n);
  fread(&l, sizeof(int), 1, qap);
  for(i=0; i<n; i++)
  {
    mpz_init(list[i]);
    mpz_inp_raw(list[i], qap);
  }
  
/*  read in left, right and output coefficient 
    matricies for each wire */

  k = n*m;
  vj = (mpz_t*)malloc(sizeof(mpz_t)*k);
  wj = (mpz_t*)malloc(sizeof(mpz_t)*k);
  yj = (mpz_t*)malloc(sizeof(mpz_t)*k);
  for(i=0; i<k; i++)
  {
    mpz_init(vj[i]);
    mpz_inp_raw(vj[i], qap);
  }
  for(i=0; i<k; i++)
  {
    mpz_init(wj[i]);
    mpz_inp_raw(wj[i], qap);
  }
  for(i=0; i<k; i++)
  {
    mpz_init(yj[i]);
    mpz_inp_raw(yj[i], qap);
  }
  j = n - 1;
  k = n*j*j/2;
  htable = (mpz_t*)malloc(sizeof(mpz_t)*k);
  for(i=0; i<k; i++)
  {
    mpz_init(htable[i]);
    mpz_inp_raw(htable[i], qap);
  }
  fclose(qap);
/*
  printf("n = %d\n", n);
  printf("m = %d\n", m);
  for(i=0; i<m; i++)
    printf("sw[%d] = %d\n", i, sw[i]);
  printf("list:\n");
  for(i=0; i<n; i++)
    gmp_printf("%Zd ", list[i]);
  printf("\n");
  printf("vj: \n");
  for(i=0; i<m; i++)
  {
    for(j=0; j<n; j++)
      gmp_printf("%Zd ", vj[i*n + j]);
    printf("\n");
  }
  printf("wj: \n");
  for(i=0; i<m; i++)
  {
    for(j=0; j<n; j++)
      gmp_printf("%Zd ", wj[i*n + j]);
    printf("\n");
  }
  printf("yj: \n");
  for(i=0; i<m; i++)
  {
    for(j=0; j<n; j++)
      gmp_printf("%Zd ", yj[i*n + j]);
    printf("\n");
  }
  printf("h table:\n");
  for(i=0; i<m; i++)
  {
    for(j=0; j<n-1; j++)
      gmp_printf("%Zd ", htable[i*(n-1) + j]);
    printf("\n");
  }
*/
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
    point_printf(" ", thtaG[sw[i]]);
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
  /* COMPUTE PROVER VALUES  */
/* create statement and witness values  */

  mset(sig.tor);                              // back to torsion value
  aj = (mpz_t*)malloc(sizeof(mpz_t)*m);
  for(i=0; i<m; i++)
    mpz_init(aj[i]);
  mpz_set_ui(aj[1], 2036);                    // statement
  mpz_set_ui(aj[2], 1700000);                 // statement
  mpz_set_ui(aj[3], 49);                      // witness
/*
  mpz_set_ui(aj[1], 2);
  mpz_set_ui(aj[2], 3);
  mpz_set_ui(aj[3], 5);
*/
  wires(aj, aj[1], aj[2], aj[3]);
/*  for(i=0; i<m; i++)
    gmp_printf("aj[%d]: %Zd\n", i, aj[i]);
*/  
/* build up h(x) values  */
  
  coef = (mpz_t*)malloc(sizeof(mpz_t)*m);
  for(i=0; i<m; i++)
    mpz_init(coef[i]);
  crossterms(coef, aj);
  hxcoef = (mpz_t*)malloc(sizeof(mpz_t)*j);
  for(i=0; i<j; i++)
    mpz_init(hxcoef[i]);
  matflat(hxcoef, htable, j, coef, m);

/* choose random values r and s */

  mpz_inits(r, s, NULL);
  rand_rs_get(r, s, sig.tor);

  av = (mpz_t*)malloc(sizeof(mpz_t)*n);
  aw = (mpz_t*)malloc(sizeof(mpz_t)*n); 
  for(j=0; j<n; j++)
  {
    mpz_init(av[j]);
    mpz_init(aw[j]);
  }
  matflat(av, vj, n, aj, m);
  matflat(aw, wj, n, aj, m);

/* Compute point A  */

  mset(sig.prime);                              // assume irrd stays ok!
  point_init(&A);
  point_init(&Tmp);
  point_copy(&A, aG);
  j = n - 1;
  for(i=0; i<n; i++)
  {
    elptic_mul(&Tmp, zG[j - i], av[i], sig.E);
    elptic_sum(&A, A, Tmp, sig.E);
  }
  elptic_mul(&Tmp, dG, r, sig.E);
  elptic_sum(&A, A, Tmp, sig.E);
  point_printf("A: ", A);
  printf("\n");
  
/* compute point B */

  poly_point_init(&B);
  poly_point_init(&P);
  poly_point_copy(&B, bH);
  for(i=0; i<n; i++)
  {
    poly_elptic_mul(&P, zH[j - i], aw[i], sig.Ex);
    poly_elptic_sum(&B, B, P, sig.Ex);
  }
  poly_elptic_mul(&P, dH, s, sig.Ex);
  poly_elptic_sum(&B, B, P, sig.Ex);
  poly_point_printf("B: ", B);
  
/* compute point C */

  point_init(&C);
  point_copy(&C, bG);            //  r(betaG + awG)
  for(i=0; i<n; i++)
  {
    elptic_mul(&Tmp, zG[j - i], aw[i], sig.E);
    elptic_sum(&C, C, Tmp, sig.E);
  }
  elptic_mul(&C, C, r, sig.E);
  
  point_copy(&Tmp, A);
  elptic_mul(&Tmp, Tmp, s, sig.E);   // + s A
  elptic_sum(&C, C, Tmp, sig.E);

  for(i=l+1; i<m; i++)
  {
    elptic_mul(&Tmp, thtaG[i], aj[sw[i]], sig.E);
    elptic_sum(&C, C, Tmp, sig.E);
  }
  j = n - 1;
  for(i=0; i<j; i++)                          // + h(z)*t(z)/delta G
  {
    elptic_mul(&Tmp, ztG[j - 1 - i], hxcoef[i], sig.E);
    elptic_sum(&C, C, Tmp, sig.E);
  }
  point_printf("C: ", C);

/*  write out proof for this record */

  printf("write out record 0\n");
  qap = fopen("snark_record.0", "w");
  fwrite(&n, sizeof(int), 1, qap);
  fwrite(&m, sizeof(int), 1, qap);
  fwrite(&l, sizeof(int), 1, qap);
  for(i=0; i<=l; i++)
    mpz_out_raw(qap, aj[sw[i]]);
  point_write(&A, qap);
  point_write(&C, qap);
  poly_point_write(&B, qap);
  fclose(qap);

  for(i=0; i<n; i++)
    mpz_clear(list[i]);
  for(i=0; i<m; i++)
    mpz_clear(coef[i]);
  j = n - 1;
  k = m*j;
  for(i=0; i<k; i++)
    mpz_clear(htable[i]);
  free(htable);
  for(i=0; i<j; i++)
    mpz_clear(hxcoef[i]);
  k = n * m;
  for(i=0; i<k; i++)
  {
    mpz_clear(vj[i]);
    mpz_clear(wj[i]);
    mpz_clear(yj[i]);
  }
  free(vj);
  free(wj);
  free(yj);
  for(i=0; i<m; i++)
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
  for(i=0; i<j; i++)
    point_clear(&ztG[i]);
  for(j=0; j<n; j++)
  {
    mpz_clear(av[j]);
    mpz_clear(aw[j]);
  }
  point_clear(&A);
  point_clear(&Tmp);
  mpz_clears(r, s, NULL);
  poly_point_clear(&B);
  poly_point_clear(&P);
  point_clear(&C);
}
