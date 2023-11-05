/**********************************************************************
 *                                                                    *
 *   Create signatures using pairings.  Use paper "Compact            *
 *   Multi-Signatures for Smaller Blockchains" for algorithms.        *
 *   Example is for GF(p^11) found from pairing_sweep_alpha and       *
 *   get_curve.                                                       *
 *                                                                    *
 *********************************************************************/

#include "signature.h"
#include <string.h>

#define EMBED_DEGREE 11

int main(int argc, char *argv[])
{
  FILE *key;
  int i, j, k, m;
  SIG_SYSTEM sig;
  POLY_POINT PK[32], APK;
  mpz_t  sk[32], aj[32];
  unsigned short *msg;
  POINT S[32], Sigma;
  long list[32], nmlst;
  POINT *Sgrp, Ssum;
  POINT *Muvec, *Mumatrix, *Memkey;
  POLY_POINT Pgrpsum;

/*  read in system parameters from previously generated file  */

  get_system("curve_11_parameters.bin", &sig);
  minit(sig.prime);
  poly_irrd_set(sig.irrd);
  poly_mulprep(sig.irrd);
  
/*  read in public and private keys for testing signature aggregation  */
  
  key = fopen("key_data_11.skpk", "r");
  if(!key)
  {
    printf("can't find file key_data_11.skpk\n");
    exit(-2);
  }
  fread(&k, sizeof(int), 1, key);
  if(k > 32)
  {
    printf("need more space for testing: %d\n", k);
    exit(-3);
  }
  for(i=0; i<k; i++)
  {
    poly_point_init(&PK[i]);
    mpz_init(sk[i]);
    mpz_inp_raw(sk[i], key);
    poly_point_read(&PK[i], key);
  }
  fclose(key);
  printf("computing aj_hash:\n");
  aj_hash(aj, sig, PK, k);
  
/*  combine the a_j*PK_j into one aggregate point of all public key * hashes  */

  poly_point_init(&APK);  // automatically point at infinity
  aj_sum(&APK, sig, PK, aj, k);
  poly_point_printf("APK:\n", APK);

/* read in random data (from random.org) for message to sign  */

  key = fopen("message.dat", "r");
  if(!key)
  {
    printf("can't find message.dat\n");
    exit(-4);
  }
  msg = (unsigned short*)malloc(sizeof(unsigned short)*102);
  for(i=0; i<100; i++)
    fscanf(key, "%hd", &msg[i]);
  fclose(key);

/*  compute signature for each key  */

  for(i=0; i<k; i++)
  {
    point_init(&S[i]);
    sign(&S[i], sig, sk[i], aj[i], (unsigned char*)msg, 200);
  }

/*  combiner takes all signatures and creates aggregate  */

  point_init(&Sigma);
  agregat_sig(&Sigma, sig, S, k);
  point_printf("Sigma:\n", Sigma);

/*  compute Weil pairings of aggregate signatures and keys with
    base G2 and message: e(Sigma, G2) = e(H0, apk)??            */

  if(multisig_verify(sig, Sigma, APK, (unsigned char*)msg, 200))
    printf("e(sigma, g2) matches e(H0, apk)\n");
  else
    printf("Signature FAILS verification\n");

/*  Begin subgroup testing section */

/* allocate space for column vector and matrix for membership key */

  Muvec = (POINT*)malloc(k*sizeof(POINT));
  Mumatrix = (POINT*)malloc(k*k*sizeof(POINT));
  for(i=0; i<k; i++)
  {
    point_init(&Muvec[i]);
    for(j=0; j<k; j++)
      point_init(&Mumatrix[i*k + j]);
  }

/* Fill mu matrix one column at a time. This simulates each node
   computing their column data and transmitting all values except
   diagonal to other nodes.
*/

  for(j=0; j<k; j++)
  {
    mu_column(Muvec, APK, sig, aj[j], sk[j], k);
    for(i=0; i<k; i++)
      point_copy(&Mumatrix[i*k + j], Muvec[i]);
  }
  
/*  with all components of a row collected including their own
    diagonal terms each node computes its membership key.
*/

  Memkey = (POINT*)malloc(k*sizeof(POINT));
  for(i=0; i<k; i++)
  {
    point_init(&Memkey[i]);
    membership_key(&Memkey[i], &Mumatrix[i*k], k, sig.E);
  }

/*  from random.org:
22	46	5	28	42
2	53	62	59	58
21	36	32	17	33
36	32	62	63	53
61	48	55	10	24
54	41	27	50	36
19	53	16	7	27
42	37	8	59	34
28	8	54	5	43
14	55	57	16	8
4	26	56	51	27
39	26	41	58	21
40	33	50	6	26
57	51	23	3	19
52	28	47	54	54
31	33	19	17	58
36	38	56	7	40
28	25	30	5	1
16	60	60	25	48
56	31	57	31	53
*/
  
// for this test choose random keys from list

  list[0] = 5;
  list[1] = 2;
  list[2] = 17;
  list[3] = 10;
  list[4] = 16;
  list[5] = 19;
  list[6] = 7;
  list[7] = 8;
  list[8] = 14;
  list[9] = 4;
  nmlst = 10;

/*  use list to create subgroup public key 
    and signature of message */

  Sgrp = (POINT*)malloc(nmlst*sizeof(POINT));
  for(i=0; i<nmlst; i++)
  {
    point_init(&Sgrp[i]);
    j = list[i];
    subgrp_sign(&Sgrp[i], APK, sig, Memkey[j], sk[j], (unsigned char*)msg, 200);
  }

/*  send all signatures to combiner to create final signature data block */

  poly_point_init(&Pgrpsum);
  point_init(&Ssum);
  subgrp_combine(&Pgrpsum, &Ssum, sig, Sgrp, PK, list, nmlst);

  printf("subgroup index list:\n");
  for(i=0; i<nmlst; i++)
    printf("%ld ", list[i]);
  printf("\n");
  printf("subgroup signature:\n");
  poly_point_printf("Public Key Aggregation:\n", Pgrpsum);
  point_printf("Aggregate Signature: ", Ssum);
  printf("\n");
  
/*  compute three pairings and check that signatures match  */

  j = subgrp_verify(sig, APK, (unsigned char*)msg, 200,
		    list, nmlst, Pgrpsum, Ssum);
  if(j)
    printf("Subgroup Aggregate Signature Verifies!\n");
  else
    printf("Subgroup Aggregate Signature FAILS!!\n");
  free(msg);
  free(Muvec);
  free(Mumatrix);
  free(Memkey);
  free(Sgrp);
  for(i=0; i<k; i++)
  {
    poly_point_clear(&PK[i]);
    point_clear(&S[i]);
    mpz_clear(aj[i]);
  }
  poly_point_clear(&APK);
  point_clear(&Sigma);
}
