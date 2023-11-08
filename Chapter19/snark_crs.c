/***********************************************************
 *                                                         *
 *    Create CRS for QAP data                              *
 *                                                         *
 **********************************************************/

#include "signature.h"
#include "snarkbase.h"

int main(int argc, char *argv[])
{
  FILE *qap, *crs;
  mpz_t *list;
  int i, j, k, m, n, *sw, l;
  mpz_t *htable;
  mpz_t *vj, *wj, *yj;
  mpz_t alpha, beta, gamma, delta, z;
  mpz_t *zpow, *theta, tmp, tdlta;
  SIG_SYSTEM sig;
  POINT aG, bG, dG, zG[5], thtaG[10], ztG[4];
  POLY_POINT bH, gH, dH, zH[5];
  mpz_t r, s;
  POINT A, C, Tmp, V;
  POLY_POINT B, P, R;
  POLY eab, eVH, eCH, eAB;

/*  read in system parameters from previously generated file  */

  get_system("curve_11_parameters.bin", &sig);

/* "working in exponent" so use torsion for prime modulus */
  
  minit(sig.tor);

/*  read in QAP parameters */

  printf("read in QAP parameters\n");
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
/* choose random values for common reference string.
   These should be erased so no one can re-use them
   to create invalid proofs.
*/
  mpz_inits(alpha, beta, gamma, delta, z, tmp, NULL);
  mrand(alpha);
  mrand(beta);
  mrand(gamma);
  mrand(delta);
  mrand(z);

  /*  create all powers of z for each aj coefficient */

  printf("compute powers of z\n");
  zpow = (mpz_t*)malloc(sizeof(mpz_t)*n);
  for(i=0; i<n; i++)
    mpz_init(zpow[i]);
  mpz_set_ui(zpow[0], 1);
  for(i=1; i<n; i++)
    mmul(zpow[i], zpow[i - 1], z);

/* create beta*vj + alpha*wj + yj for each aj coefficient */

  printf("compute theta terms\n");
  theta = (mpz_t*)malloc(sizeof(mpz_t)*m);
  for(i=0; i<m; i++)
    mpz_init(theta[i]);

  for(i=0; i<m; i++)
  {
    lcalc(tmp, z, &vj[i*n], j);     // j = n-1
    mmul(theta[i], tmp, beta);      // beta * vj
    lcalc(tmp, z, &wj[i*n], j);
    mmul(tmp, tmp, alpha);
    madd(theta[i], theta[i], tmp);  // + alpha * wj
    lcalc(tmp, z, &yj[i*n], j);
    mmul(tmp, tmp, tmp);     
    madd(theta[i], theta[i], tmp);  // + yj^2
  }

/* statement indexes are 0 to l
   witness indexes are l+1 to m-1  */

  for(i=0; i<=l; i++)
    mdiv(theta[sw[i]], theta[sw[i]], gamma);
  for(i=l+1; i<m; i++)
    mdiv(theta[sw[i]], theta[sw[i]], delta);

  mpz_init(tdlta);
  tofzgrth(tmp, z, list, n);
  mdiv(tdlta, tmp, delta);

/* create CRS points on G1 
   start by initializing field prime and extension field.  */

  printf("create CRS points on G1\n");
  mset(sig.prime);
  poly_irrd_set(sig.irrd);
  poly_mulprep(sig.irrd);
  point_init(&aG);
  elptic_mul(&aG, sig.G1, alpha, sig.E);
  point_init(&bG);
  elptic_mul(&bG, sig.G1, beta, sig.E);
  point_init(&dG);
  elptic_mul(&dG, sig.G1, delta, sig.E);
  for(i=0; i<n; i++)
  {
    point_init(&zG[i]);
    elptic_mul(&zG[i], sig.G1, zpow[i], sig.E);
  }

  for(i=0; i<m; i++)
  {
    point_init(&thtaG[i]);
    elptic_mul(&thtaG[i], sig.G1, theta[i], sig.E);
  }

/*  terms for z^i*t(z)/delta * G */

  point_init(&ztG[0]);
  elptic_mul(&ztG[0], sig.G1, tdlta, sig.E);
  j = n - 1;
  for(i=1; i<j; i++)
  {
    point_init(&ztG[i]);
    elptic_mul(&ztG[i], ztG[i-1], z, sig.E);
  }

/* create CRS points on G2 */

  printf("create CRS points on G2\n");
  poly_point_init(&bH);
  poly_elptic_mul(&bH, sig.G2, beta, sig.Ex);
  poly_point_init(&gH);
  poly_elptic_mul(&gH, sig.G2, gamma, sig.Ex);
  poly_point_init(&dH);
  poly_elptic_mul(&dH, sig.G2, delta, sig.Ex);
  for(i=0; i<n; i++)
  {
    poly_point_init(&zH[i]);
    poly_elptic_mul(&zH[i], sig.G2, zpow[i], sig.Ex);
  }

/* now that CRS computed, write all data to disk */

  printf("\nwriting crs file\n");
  crs = fopen("snark.crs", "w");
  fwrite(&n, sizeof(int), 1, crs);
  fwrite(&m, sizeof(int), 1, crs);
  fwrite(&l, sizeof(int), 1, crs);
  point_write(&aG, crs);
  point_write(&bG, crs);
  point_write(&dG, crs);
  for(i=0; i<n; i++)
    point_write(&zG[i], crs);
  for(i=0; i<m; i++)
    point_write(&thtaG[sw[i]], crs);
  j = n - 1;
  for(i=0; i<j; i++)
    point_write(&ztG[i], crs);

  poly_point_write(&bH, crs);
  poly_point_write(&dH, crs);
  poly_point_write(&gH, crs);
  for(i=0; i<n; i++)
    poly_point_write(&zH[i], crs);
  fclose(crs);
/*
  printf("\nwriting crs file\n");
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
  mpz_clears(alpha, beta, gamma, delta, z, tmp, NULL);
  for(i=0; i<n; i++)
    mpz_clear(zpow[i]);
  for(i=0; i<m; i++)
    mpz_clear(theta[i]);
  for(i=0; i<n; i++)
    mpz_clear(list[i]);
  free(sw);
  j = n - 1;
  for(i=0; i<m*j; i++)
    mpz_clear(htable[i]);
  free(htable);
  for(i=0; i<m*n; i++)
  {
    mpz_clear(vj[i]);
    mpz_clear(wj[i]);
    mpz_clear(yj[i]);
  }
  free(vj);
  free(wj);
  free(yj);
  mpz_clear(tdlta);
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
}
