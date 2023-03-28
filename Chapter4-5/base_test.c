/****************************************************************************
 *                                                                          *
 *    Test basic protocols.  Start with key exchange.  Assume transmiison   *
 *    of public key points and setup of curves is uniform and accomplished. *
 *                                                                          *
 ***************************************************************************/

#include "base_protocols.h"

/*  skip n lines of text.  returns start of next line.  */

int skipln(char *bfr, int strt, int skp)
{
  int i;

  i = strt;
  while(skp)
  {
    while(bfr[i] != '\n') i++;
    i++;
    skp--;
  }
  return i;
}

int main(int argc, char *argv[])
{
  FILE *parm;
  mpz_t prm, sk, sok;
  long lvl;
  BASE_SYSTEM base;
  POINT Pk, Pok, Rk, Rok;
  char filename[256], *bufr, *ptr;
  int i, k;
  mpz_t myshare, theirshare, myrand, theirrand;
  ECDSA sig;
  SCHNORR snr;
  
/*  read in parameters based on level from command line  */
  
  if(argc < 2)
  {
    printf("Use: ./base_test <level>\n");
    exit(-1);
  }
  lvl = atol(argv[1]);
  sprintf(filename, "Curve_%ld_params.dat", lvl);
  mpz_inits(prm, base.order, NULL);
  parm = fopen(filename, "r");
  if(!parm)
  {
    printf("can't find file %s\n", filename);
    exit(-2);
  }
  bufr = (char*)malloc(1024*4);
  i = 0;
  while((!feof(parm)) && (i < 1024))
  {
    bufr[i] = fgetc(parm);
    i++;
  }
  fclose(parm);

/*  convert text to big numbers  */
  
  i = skipln(bufr, 0, 1);
  gmp_sscanf(&bufr[i], "%Zx", prm);
  gmp_printf("%Zd\n", prm);
  minit(prm);
  i = skipln(bufr, i, 2);
  gmp_sscanf(&bufr[i], "%Zx", base.order);
  gmp_printf("%Zd\n", base.order);
  i = skipln(bufr, i, 2);
  sscanf(&bufr[i], "%ld", &base.cofactor);
  i = skipln(bufr, i, 2);
  curve_init(&base.E);
  gmp_sscanf(&bufr[i], "%Zx %Zx", base.E.a4, base.E.a6);
  gmp_printf("E: %Zx  %Zx\n", base.E.a4, base.E.a6);
  i = skipln(bufr, i, 3);
  point_init(&base.Base);
  gmp_sscanf(&bufr[i], "%Zx %Zx", base.Base.x, base.Base.y);
  point_printf("Base point: ", base.Base);
/*
  point_init(&Rk);
  elptic_mul(&Rk, Base, ordr, E);
  point_printf("check: ", Rk);
  exit(0);
*/
/*  Generate secret key and public key from input 
    phrase.
*/

  printf("Input pass phrase to generate secret key: ");
  fflush(stdout);
  i = 1024;
  getline(&bufr, (size_t*)&i, stdin);
  mpz_init(sk);
  point_init(&Pk);
  gen_key(sk, &Pk, bufr, base);
  gmp_printf("secret key: %Zx\n", sk);
  point_printf("Public key: ", Pk);

/*  create another secret and public key, this time
    using a fixed phrase to represent "the other side".
*/

  mpz_init(sok);
  point_init(&Pok);
  sprintf(bufr, "Secret Key Test For Other Side 157 164 218 149 124 108 253 26 40 ");
  gen_key(sok, &Pok, bufr, base);
  point_printf("Other side Public key: ", Pok);

/*  Test Diffie-Hellman routine for both sides  */

  mpz_inits(myshare, theirshare, NULL);
  diffie_hellman(myshare, sk, Pok, base.E);
  diffie_hellman(theirshare, sok, Pk, base.E);
  if(!mpz_cmp(myshare, theirshare))
    printf("Keys match.\n");
  else
    printf("Keys DON'T match!\n");
  gmp_printf("my keyshare:    %Zx\n", myshare);
  gmp_printf("their keyshare: %Zx\n", theirshare);

/*  Next generate random secret and public keys for MQV test  */

  mpz_inits(myrand, theirrand, NULL);
  point_init(&Rk);
  point_init(&Rok);
  mqv_ephem(myrand, &Rk, base);
  mqv_ephem(theirrand, &Rok, base);

/*  Each side sends the other the public key, and then computes the shared secret  */

  mqv_share(myshare, sk, Pk, myrand, Rk, Pok, Rok, base);

  mqv_share(theirshare, sok, Pok, theirrand, Rok, Pk, Rk, base);
  if(!mpz_cmp(myshare, theirshare))
    printf("MQV Keys match.\n");
  else
    printf("MQV Keys DON'T match!\n");
  gmp_printf("my keyshare:    %Zx\n", myshare);
  gmp_printf("their keyshare: %Zx\n", theirshare);

/*  now create a test for a digital signature  */

  parm = fopen("sign_test.txt", "r");
  if(!parm)
  {
    printf("sign_test.txt not found??\n");
    exit(-7);
  }
  k=0;
  while((!feof(parm)) && (k < 4*1024))
  {
    bufr[k] = fgetc(parm);
    k++;
  }
  k -= 2;
  sig_init(&sig);
  ecdsa_sign(&sig, sk, Pk, bufr, k, base);

/*  and check that the verify routine works  */

  if(ecdsa_verify(sig, Pk, bufr, k, base))
    printf("message verifies!\n");
  else
    printf("message does not match original signed.\n");

/*  next perform Schnorr signature  */

  snr_init(&snr);
  schnorr_sign(&snr, sk, Pk, bufr, k, base);
  
/*  and check that the verify routine works  */

  if(schnorr_verify(snr, Pk, bufr, k, base))
    printf("Schnorr message verifies!\n");
  else
    printf("Schnorr message does not match original signed.\n");

  free(bufr);
  mpz_clears(sk, sok, theirshare, myshare, prm, base.order, NULL);
  mpz_clears(myrand, theirrand, NULL);
  point_clear(&base.Base);
  point_clear(&Pk);
  point_clear(&Pok);
  curve_clear(&base.E);
  snr_clear(&snr);
  sig_clear(&sig);
}
