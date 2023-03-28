/**********************************************************************
 *                                                                    *
 *    Subroutines to compute Aggregated Multiple Signature algorithm. *
 *                                                                    *
 *********************************************************************/

#include "signature.h"

/*  Low level key Generation routine.
    Input variable ptr is a string of len bytes,
    Extension Curve Ex and G2 base point.  
    Output is private key sk and public key PK.
*/

void keygen(mpz_t sk, POLY_POINT *PK, POLY_POINT G2, POLY_CURVE Ex,
	    unsigned char *ptr, int len)
{

  mpz_import(sk, len, -1, 1, 0, 0, ptr);
  poly_elptic_mul(PK, G2, sk, Ex);
}

/*  add number of bits of security to bits in prime to find number of
    bytes required from extensible hash. Used in hash0 and hash1.
*/

long secbyte(mpz_t p)
{
  int m, k, b;
  
/*  find security level and add to bit size of prime  */
  
  m = mpz_sizeinbase(p, 2);
  if(m < 208)
    k = 80;
  else if(m < 320)
    k = 128;
  else if(m < 448)
    k = 192;
  else
    k = 256;
  b = m + k + 7;
  b >>= 3;
  return b;
}

/*  Hash to GF(r) using KangarooTwelve.
    Enter with pointer to data as bytes, number of bytes and torsion
    value r.
    Returns value mod r.
    NOTE: we assume p is twice the protection level.  Protection
    levels are 80, 128, 192 and 256 bits, so we add 10, 16, 24 or
    32 bytes to ensure result mod p is actually random.  See
    IETF irtf-cfrg-hash-to-curve-11 for details section 5.
*/

void hash1(mpz_t hsh, unsigned char *dat, long len, mpz_t r)
{
  unsigned char *outp, *dst;
  mpz_t zero;
  long b;
  
  dst = (char*)malloc(24);
  sprintf(dst, "Hash_1 pring&sig");

/*  find security level and add to bit size of prime  */

  b = secbyte(r);
  outp = (unsigned char*)malloc(b + 2);
  KangarooTwelve(dat, len, outp, b, dst, 16);

/*  convert hash bytes into an integer  */

  mpz_import(hsh, b, -1, 1, 0, 0, outp);
  mpz_init(zero);
  mod_add(hsh, hsh, zero, r);      // force to be mod r
  mpz_clear(zero);
  free(outp);
  free(dst);
}

/*  Hash to G1 using elptic_embed. Same issue as hash1
    with expansion of output and conversion to integer.
*/

void hash0(POINT *H, SIG_SYSTEM sig, unsigned char *dat, long len)
{
  unsigned char *outp, *dst;
  mpz_t hsh, zero, prm;
  POINT mP;
  long b;

  dst = (char*)malloc(24);
  sprintf(dst, "signat H_0 parng");
  mget(prm);
  b = secbyte(prm);
  outp = (unsigned char*)malloc(b+2);
  KangarooTwelve(dat, len, outp, b, dst, 16);

/*  convert hash bytes into an integer  */

  mpz_inits(hsh, zero, NULL);
  mpz_import(hsh, b , -1, 1, 0, 0, outp);
  madd(hsh, hsh, zero);  
  point_init(&mP);
  elptic_embed(H, &mP, hsh, sig.E);
  elptic_mul(H, *H, sig.cobse, sig.E);
  mpz_clears(hsh, zero, prm, NULL);
  point_clear(&mP);
  free(outp);
  free(dst);
}
  
/*  Hash to G1 using elptic_embed. Only difference with
    hash0 is domain separation tag.
*/

void hash2(POINT *H, SIG_SYSTEM sig, unsigned char *dat, long len)
{
  unsigned char *outp, *dst;
  mpz_t hsh, zero, prm;
  POINT mP;
  long b;

  dst = (char*)malloc(24);
  sprintf(dst, "signat H_2 parng");
  mget(prm);
  b = secbyte(prm);
  outp = (unsigned char*)malloc(b+2);
  KangarooTwelve(dat, len, outp, b, dst, 16);

/*  convert hash bytes into an integer  */

  mpz_inits(hsh, zero, NULL);
  mpz_import(hsh, b , -1, 1, 0, 0, outp);
  madd(hsh, hsh, zero);  
  point_init(&mP);
  elptic_embed(H, &mP, hsh, sig.E);
  elptic_mul(H, *H, sig.cobse, sig.E);
  mpz_clears(hsh, zero, prm, NULL);
  point_clear(&mP);
  free(outp);
  free(dst);
}
  
/* read in signature curve parameters.
   Enter with filename and pointer to system structure.
   This routine will intialize all spaces in the
   struct.
*/

void get_system(char *filename, SIG_SYSTEM *sig)
{
  FILE *sys;
  int i;

  sys = fopen(filename, "r");
  if(!sys)
  {
    printf("can't find file %s\n", filename);
    exit(-1);
  }
  // initialize struct space
  mpz_inits(sig->prime, sig->cardE, sig->tor, sig->cobse, sig->cardEx, sig->coxtd, NULL);
  poly_init(&sig->irrd);
  curve_init(&sig->E);
  point_init(&sig->G1);
  poly_point_init(&sig->G2);
  poly_curve_init(&sig->Ex);

  // read in raw data
  mpz_inp_raw(sig->prime, sys);
  mpz_inp_raw(sig->E.a4, sys);
  mpz_inp_raw(sig->E.a6, sys);
  mpz_inp_raw(sig->cardE, sys);
  mpz_inp_raw(sig->tor, sys);
  mpz_inp_raw(sig->cobse, sys);
  mpz_inp_raw(sig->G1.x, sys);
  mpz_inp_raw(sig->G1.y, sys);

  poly_read(&sig->irrd, sys);
  mpz_inp_raw(sig->cardEx, sys);
  mpz_inp_raw(sig->coxtd, sys);
  poly_point_read(&sig->G2, sys);
  mpz_set(sig->Ex.a4.coef[0], sig->E.a4);
  mpz_set(sig->Ex.a6.coef[0], sig->E.a6);
}

/*  convert G1 point to G2 point.  
    Enter with space for G2 initalized.
*/

void tog2(POLY_POINT *G2, POINT G1)
{
  G2->x.deg = 0;
  mpz_set(G2->x.coef[0], G1.x);
  G2->y.deg = 0;
  mpz_set(G2->y.coef[0], G1.y);
}

/*  copy a poly_point to character buffer.
    Enter with pointer to buffer location, poly_point and
    size of each coefficient in bytes and number of coefficients.
*/

void point2text(unsigned char *ptr, POLY_POINT P, long msz, long deg)
{
  long j;

  for(j=0; j<2*deg*msz; j++)
    ptr[j] = 0;
  for(j=0; j<deg; j++)
    mpz_export(&ptr[msz*j], NULL, -1, 1, 0, 0, P.x.coef[j]);
  for(j=0; j<deg; j++)
    mpz_export(&ptr[msz*(deg + j)], NULL, -1, 1, 0, 0, P.y.coef[j]);
}

/*  compute self aggregate hash for every key. 
    Enter with number of public keys, array of public keys and
    system paramters structure.
    Each ajhsh[j] is initialized here.  Should be numkey long.  */

void aj_hash(mpz_t *ajhsh, SIG_SYSTEM sig, POLY_POINT *PK,
	     long numkey)
{
  long i, j, k, m, bfsz;
  unsigned char *srcbfr;

  m = (mpz_sizeinbase(sig.prime, 16) + 1)/2;
  bfsz = m*sig.irrd.deg*2*(numkey + 1);
  srcbfr = (unsigned char*)malloc(bfsz);
  for(i=0; i<numkey; i++)
    point2text(&srcbfr[m*(sig.irrd.deg*2*(i + 1))], PK[i], m, sig.irrd.deg);

  for(i=0; i<numkey; i++)
  {
    mpz_init(ajhsh[i]);
    point2text(srcbfr, PK[i], m, sig.irrd.deg);
    hash1(ajhsh[i], srcbfr, bfsz, sig.tor);
  }
  free(srcbfr);
}

/*  compute individual signature of a message.
    Enter with pointer to message, it's length, self aggregate public key,
    secret key and sig struct.  Returns S_j which must be pre-initialized.
*/

void sign(POINT *S, SIG_SYSTEM sig, mpz_t sk, mpz_t aj, 
	  unsigned char *msg, long msgln)
{
  mpz_t xpnt;
  POINT H0;

  point_init(&H0);
  hash0(&H0, sig, msg, msgln);
  mpz_init(xpnt);
  mpz_mul(xpnt, aj, sk);
  elptic_mul(S, H0, xpnt, sig.E);
  mpz_clear(xpnt);
  point_clear(&H0);
}

/*  combine all individual signatures into a final aggregate.
    Enter with sig struct, array of signatures and number of them.
    Return Sigma must be pre-initialized.
*/

void agregat_sig(POINT *Sigma, SIG_SYSTEM sig, POINT *S, long nmsig)
{
  int i;
  
  point_copy(Sigma, S[0]);
  if(nmsig == 1) return;
  for(i=1; i<nmsig; i++)     
    elptic_sum(Sigma, *Sigma, S[i], sig.E);
}

/*  compute apk = sum(a_j * PK_j).
    Enter with sig struct, public key list and 
    self aggregate hash list.
*/

void aj_sum(POLY_POINT *APK, SIG_SYSTEM sig, POLY_POINT *PK,
	    mpz_t *ajhsh, long nmkey)
{
  int i;
  POLY_POINT Tmp;

  APK->x.deg = 0;
  mpz_set_ui(APK->x.coef[0], 0);
  APK->y.deg = 0;
  mpz_set_ui(APK->y.coef[0], 0);
  poly_point_init(&Tmp);
  for(i=0; i<nmkey; i++)
  {
    poly_elptic_mul(&Tmp, PK[i], ajhsh[i], sig.Ex);
    poly_elptic_sum(APK, *APK, Tmp, sig.Ex);
  }
  poly_point_clear(&Tmp);
}
  
/*  Verify aggregate signature.
    Enter with G1 points Sigma (aggregate signature), message hash H0,
    sig struct and G2 point APK.  Returns 1 if good, 0 if bad.  */

int multisig_verify(SIG_SYSTEM sig, POINT Sigma, POLY_POINT APK,
		    unsigned char *msg, long msgln)
{
  POLY_POINT S2, H20, R;
  POINT H0;
  POLY w1, w2;
  int eql;
  
  poly_point_init(&S2);
  tog2(&S2, Sigma);

/*  compute Weil pairings of aggregate signatures and keys with
    base G2 and message: e(Sigma, G2) = e(H0, apk)??            */

  point_init(&H0);
  hash0(&H0, sig, msg, msgln);
  poly_point_init(&H20);
  poly_point_init(&R);
  tog2(&H20, H0);
  poly_init(&w1);
  poly_init(&w2);
  poly_point_rand(&R, sig.Ex);
  weil(&w1, S2, sig.G2, R, sig.tor, sig.Ex);
  poly_printf("e(sigma, g2): ", w1);
  weil(&w2, H20, APK, R, sig.tor, sig.Ex);
  poly_printf("e(H0, apk): ", w2);
  eql = poly_cmp(w1, w2);
  poly_clear(&w1);
  poly_clear(&w2);
  poly_point_clear(&H20);
  poly_point_clear(&S2);
  poly_point_clear(&R);
  point_clear(&H0);
  return eql;
}

/*  compute column of membership key matrix.
    mu_(i,j) = a_j*sk*H_2(APK, i)
    Enter with private key, this nodes a_j, 
    the full signature parameter set, 
    APK = sum of all a_j*P_j and maximum i. 
    Returns vector of G_1 points. 
    Assumes muvec intialized to nmkey depth.
*/

void mu_column(POINT *muvec, POLY_POINT APK, SIG_SYSTEM sig,
	       mpz_t aj, mpz_t sk, long nmkey)
{
  long i, j, bfsz, m;
  mpz_t apsi;
  unsigned char *hshbfr;
  
  mpz_init(apsi);
  mod_mul(apsi, sk, aj, sig.tor);
  m = (mpz_sizeinbase(sig.prime, 16) + 1)/2;
  bfsz = m*sig.irrd.deg*2 + sizeof(long);
  hshbfr = (unsigned char*)malloc(bfsz + 8);
  point2text(hshbfr, APK, m, sig.irrd.deg);
  for(i=0; i<nmkey; i++)
  {
    *(long*)(&hshbfr[m*sig.irrd.deg*2]) = i;
    hash2(&muvec[i], sig, hshbfr, bfsz);
    elptic_mul(&muvec[i], muvec[i], apsi, sig.E);
  }
  free(hshbfr);
  mpz_clear(apsi);
}

/* Every node computes their membership key as the sum
   of a row from the mu matrix. Each node will call this
   routine once.
*/

void membership_key(POINT *Memkey, POINT *Murow, long numkey, CURVE E)
{
  long j;
  POINT Tmp;

  point_init(&Tmp);
  for(j=0; j<numkey; j++)
      elptic_sum(&Tmp, Tmp, Murow[j], E);
  point_copy(Memkey, Tmp);
  point_clear(&Tmp);
}

/*  subgroup signature subroutine
    Enter with signature system parameters, APK,
    membership key, private key, message to sign and
    length of message.
    Returns signature as point on G_1.
 */
void subgrp_sign(POINT *S, POLY_POINT APK, SIG_SYSTEM sig,
	   POINT Mk, mpz_t sk, unsigned char *msg, long msgln)
{
  long j, m, bfsz;
  unsigned char *apkmsg;
  POINT H0;

  m = (mpz_sizeinbase(sig.prime, 16) + 1)/2;
  bfsz = m*sig.irrd.deg*2 + msgln;
  apkmsg = (unsigned char*)malloc(bfsz + 8);
  point2text(apkmsg, APK, m, sig.irrd.deg);
  for(j=0; j<msgln; j++)
    apkmsg[m*2*sig.irrd.deg + j] = msg[j];
  point_init(&H0);
  hash0(&H0, sig, apkmsg, bfsz);
  elptic_mul(S, H0, sk, sig.E);
  elptic_sum(S, *S, Mk, sig.E);
  point_clear(&H0);
  free(apkmsg);
}

/*  subgroup combine subroutine
    Enter with index list of signers, number of them,
    vector of signatures, vector of all public keys and
    signature system parameters.
    Returns sum of public keys and signatures in list order.
    Assumes result spaces already intialized.
*/

void subgrp_combine(POLY_POINT *PK, POINT *Ssum, SIG_SYSTEM sig,
		    POINT *Svec, POLY_POINT *Pubky, long *list, long nmlst)
{
  long i, j;

  mpz_set_ui(Ssum->x, 0);
  mpz_set_ui(Ssum->y, 0);
  PK->x.deg = 0;
  mpz_set_ui(PK->x.coef[0], 0);
  PK->y.deg = 0;
  mpz_set_ui(PK->y.coef[0], 0);
  for(i=0; i<nmlst; i++)
  {
    j = list[i];
    elptic_sum(Ssum, *Ssum, Svec[i], sig.E);
    poly_elptic_sum(PK, *PK, Pubky[j], sig.Ex);
  }
}

/*  Verify subgroup signature.
    Enter with signature system parameters,
    APK = sum of a_j*P_j,
    message being signed and its length,
    ordered list of keys used to create signature
      and public key sum and its length,
    aggregate sum of signatures and public keys.
    Returns 1 on pass, 0 on fail.
*/

int subgrp_verify(SIG_SYSTEM sig, POLY_POINT APK, unsigned char *msg, long msgln,
		  long *list, long nmlst, POLY_POINT PK, POINT Ssum)
{
  POLY w1, w2, w3, ck;
  POINT *H2, Hsum, H0;
  POLY_POINT H02, R;
  long i, j, bfsz, m;
  unsigned char *bufr;
  int eql;

/* compute hash of APK and message  */

  m = (mpz_sizeinbase(sig.prime, 16) + 1)/2;
  bfsz = m*sig.irrd.deg*2 + msgln;
  bufr = (unsigned char*)malloc(bfsz + sizeof(long));
  point2text(bufr, APK, m, sig.irrd.deg);
  for(j=0; j<msgln; j++)
    bufr[m*2*sig.irrd.deg + j] = msg[j];

/* compute pairing of hash point with aggregate sum of keys */
  
  point_init(&H0);
  hash0(&H0, sig, bufr, bfsz);
  poly_point_init(&H02);
  tog2(&H02, H0);
  poly_point_init(&R);
  poly_point_rand(&R, sig.Ex);
  poly_init(&w1);
  weil(&w1, H02, PK, R,  sig.tor, sig.Ex);
  poly_printf("e(H0(A, m), PK)= ", w1);
  
/* compute array of points H2(APK, j) in list order.
   note that APK is already in buffer. */

  H2 = (POINT*)malloc(sizeof(POINT)*nmlst);
  bfsz = m*sig.irrd.deg*2 + sizeof(long);
  for(i=0; i<nmlst; i++)
  {
    j = list[i];
    point_init(&H2[i]);
    *(long*)(&bufr[m*2*sig.irrd.deg]) = j;
    hash2(&H2[i], sig, bufr, bfsz);
  }
  
/* sum all hashes to a single point and compute pairing with APK */

  point_init(&Hsum);
  for(i=0; i<nmlst; i++)
    elptic_sum(&Hsum, Hsum, H2[i], sig.E);

  tog2(&H02, Hsum);
  poly_init(&w2);
  weil(&w2, H02, APK, R, sig.tor, sig.Ex);
  poly_printf("e(H2(APK,j), APK) = ", w2);
  
/* compute pairing of subgroup aggregate signature with base point  */

  tog2(&H02, Ssum);
  poly_init(&w3);
  weil(&w3, H02, sig.G2, R, sig.tor, sig.Ex);

  poly_printf("e(s, g_2): ", w3);
  poly_init(&ck);
  poly_mul(&ck, w1, w2);
  eql = poly_cmp(ck, w3);
  poly_clear(&w1);
  poly_clear(&w2);
  poly_clear(&w3);
  poly_clear(&ck);
  for(i=0; i<nmlst; i++)
    point_clear(&H2[i]);
  point_clear(&Hsum);
  point_clear(&H0);
  free(bufr);
  poly_point_clear(&H02);
  poly_point_clear(&R);
  return eql;
}

/* Write and read points, polynomials and poly_points to
   a disk file. Enter with pointer to object, and FILE pointer.
   On reads, object must be pre-initialized.
*/

void point_write(POINT *P, FILE *f)
{
  mpz_out_raw(f, P->x);
  mpz_out_raw(f, P->y);
}

void point_read(POINT *P, FILE *f)
{
  mpz_inp_raw(P->x, f);
  mpz_inp_raw(P->y, f);
}

void poly_write(POLY *p, FILE *f)
{
  int i, k;

  k = p->deg;
  fwrite(&k, sizeof(long), 1, f);
  for(i=0; i<=k; i++)
    mpz_out_raw(f, p->coef[i]);
}

void poly_read(POLY *p, FILE *f)
{
  int i, k;

  fread(&k, sizeof(long), 1, f);
  p->deg = k;
  for(i=0; i<=k; i++)
    mpz_inp_raw(p->coef[i], f);
}

void poly_point_write(POLY_POINT *P, FILE *f)
{
  poly_write(&P->x, f);
  poly_write(&P->y, f);
}

void poly_point_read(POLY_POINT *P, FILE *f)
{
  poly_read(&P->x, f);
  poly_read(&P->y, f);
}
