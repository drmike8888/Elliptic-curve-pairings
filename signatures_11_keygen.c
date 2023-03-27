/********************************************************************
 *                                                                  *
 *    Create a file of public and private keys for use with key     *
 *    aggregation schemes.                                          *
 *                                                                  *
 *******************************************************************/

#include "signature.h"
#include <string.h>

#define EMBED_DEGREE 11

int main(int argc, char *argv[])
{
  FILE *rnd, *key;
  POLY_CURVE Ex;
  CURVE E;
  POLY_POINT G2, Tst2, PK;
  POINT G1, Tst1;
  POLY irrd;
  mpz_t prm, t, cardE, cardEx, tor;
  mpz_t cobse, coxtd, sk;
  int i, j, k;
  unsigned char *rndat;

/*  read in random data  */

  rnd = fopen("random.data", "r");
  if(!rnd)
  {
    printf("can't find file random.data\n");
    exit(-1);
  }
  rndat = (unsigned char*)malloc(2048);
  i = 0;
  while(!feof(rnd))
  {
    fscanf(rnd, "%hhd", &rndat[i]);
    i++;
  }
  fclose(rnd);

/*  parameter generation:
    prime, base curve, torsion value,
    irreducible polynomial, extension curve  */
  
  mpz_init_set_str(prm,
      "3252011917820513804209601668228184687614933135935522716335534348924920411135269", 10);
  minit(prm);

  curve_init(&E);
  mpz_set_str(E.a4,
      "2176336228457045175187176003692720598604044414672361175116793345654679844960925", 10);
  mpz_set_str(E.a6,
      "1870980220760335983335580434512960352657871567376133141577832331267457576039896", 10);
  
  poly_init(&irrd);
  if(poly_irreducible(&irrd, EMBED_DEGREE))
    poly_printf("Found irreducible polynomial:\n", irrd);
  else
    printf("no irreducible polynomial found...\n");
  poly_irrd_set(irrd);
  poly_mulprep(irrd);
  poly_curve_init(&Ex);    // this sets .deg to 0
  mpz_set(Ex.a4.coef[0], E.a4);
  mpz_set(Ex.a6.coef[0], E.a6);
  mpz_init_set_str(tor,
       "848222711223348251273816343240321778327803280959919975351835864551", 10);

/* compute #E base curve = p + 1 - t  */
  
  mpz_init_set_str(t, "3606666653515050472962077101037353515626", 10);
  mpz_inits(cardE, cardEx, cobse, coxtd, NULL);
  mpz_add_ui(cardE, prm, 1);
  mpz_sub(cardE, cardE, t);
  gmp_printf("base has %Zd points\n", cardE);
  mpz_div(cobse, cardE, tor);
  
/*  compute #E extension curves  */

  cardinality(cardEx, t, EMBED_DEGREE);
  gmp_printf("extension has %Zd points\n", cardEx);
  mpz_div(coxtd, cardEx, tor);
  mpz_div(coxtd, coxtd, tor);   // remove tor^2 for cofactor

/*  create generator points from random points and co-factors  */

  point_init(&G1);
  point_rand(&G1, E);
  elptic_mul(&G1, G1, cobse, E);
  point_printf("G1 generator: ", G1);

  poly_point_init(&G2);
  poly_point_rand(&G2, Ex);
//  gmp_printf("coxtd: %Zd\n", coxtd);
  poly_elptic_mul(&G2, G2, coxtd, Ex);
  poly_point_printf("G2 generator:\n", G2);

/*  save system parameters to disk  */

  key = fopen("curve_11_parameters.bin", "w");
  mpz_out_raw(key, prm);
  mpz_out_raw(key, E.a4);
  mpz_out_raw(key, E.a6);
  mpz_out_raw(key, cardE);
  mpz_out_raw(key, tor);
  mpz_out_raw(key, cobse);
  mpz_out_raw(key, G1.x);
  mpz_out_raw(key, G1.y);
  poly_write(&irrd, key);
  mpz_out_raw(key, cardEx);
  mpz_out_raw(key, coxtd);
  poly_point_write(&G2, key);
  fclose(key);
  
/*  create sets of keys using random data  */
  
  key = fopen("key_data_11.skpk", "w");
  poly_point_init(&PK);
  mpz_init(sk);
  k = 20;
  fwrite(&k, sizeof(int), 1, key);
  for(i=0; i<k; i++)
  {
    keygen(sk, &PK, G2, Ex, &rndat[i*32], 32);
    mpz_out_raw(key, sk);
    poly_point_write(&PK, key);
  }
  fclose(key);

/* debug output and input back is correct  */

  gmp_printf("prime: %Zd\n", prm);
  gmp_printf("E.a4: %Zd\n", E.a4);
  gmp_printf("E.a6: %Zd\n", E.a6);
  gmp_printf("cardE: %Zd\n", cardE);
  gmp_printf("tor: %Zd\n", tor);
  gmp_printf("cobse: %Zd\n", cobse);
  gmp_printf("G1.x: %Zd\n", G1.x);
  gmp_printf("G1.y: %Zd\n", G1.y);
  poly_printf("irrd:\n", irrd);
  gmp_printf("cardEx: %Zd\n", cardEx);
  gmp_printf("coxtd: %Zd\n", coxtd);
  poly_point_printf("G2:\n", G2);
}
