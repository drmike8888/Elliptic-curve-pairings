/********************************************************************
 *                                                                  *
 *   Subroutines to operate basic elliptic curve protocols.         *
 *   Includes public/private key generation, secret key sharing     *
 *   using public keys and digital signatures.                      *
 *                                                                  *
 *******************************************************************/

#include "base_protocols.h"

/*  Hash to specified order using KangarooTwelve.
    Enter with pointer to data as bytes, and number of bytes.
    Returns value mod n.
    NOTE: we assume n is twice the protection level.  Protection
    levels are 80, 128, 192 and 256 bits, so we add 10, 16, 24 or
    32 bytes to ensure result mod p is actually random.  See
    IETF irtf-cfrg-hash-to-curve-11 section 5 for details.
*/

void hash(mpz_t hsh, unsigned char *dat, long len, mpz_t prm)
{
  unsigned char *outp, *dst;
  int m, k, b;
  
  dst = (char*)malloc(24);
  sprintf(dst, "Hash_b pring&sig");

/*  find security level and add to bit size of prime  */
  
  m = mpz_sizeinbase(prm, 2);
  if(m < 208)
    k = 80;
  else if(m < 320)
    k = 128;
  else if(m < 448)
    k = 192;
  else
    k = 256;
  b = m + k;
  if(b & 7)
    b = (b >> 3) + 1;
  else
    b >>= 3;
  outp = (unsigned char*)malloc(b + 2);
  KangarooTwelve(dat, len, outp, b, dst, 16);

/*  convert hash bytes into an integer  */

  mpz_import(hsh, b, -1, 1, 0, 0, outp);
  mpz_mod(hsh, hsh, prm);      // force to be mod n
  free(outp);
  free(dst);
}

/*  convert input phrase into private and public key on given curve.
    pre-initialze sk and PK before calling.  phrase is assumed to
    be null terminated C string.  Input curve and base point should
    be same for all users.  */

void gen_key(mpz_t sk, POINT *PK, unsigned char *phrase, BASE_SYSTEM bse)
{
  long np;
  
  np = 0;
  while(phrase[np]) np++;
  hash(sk, phrase, np, bse.order);
  elptic_mul(PK, bse.Base, sk, bse.E);
}

/*  Compute Diffie-Hellman key share using my private key and other sides
    public key.  Input curve to operate on.  Returns x value of point as
    shared key.
*/

void diffie_hellman(mpz_t keyshare, mpz_t my_key, POINT Their_key, CURVE E)
{
  POINT Tmp;

  point_init(&Tmp);
  elptic_mul(&Tmp, Their_key, my_key, E);
  mpz_set(keyshare, Tmp.x);
  point_clear(&Tmp);
}

/*  create random value and use that as ephemeral key to
    generate random point for transmission to other side.
    Enter with base system.
    Returns key and point (which were pre-initialized).
*/

void mqv_ephem(mpz_t ephm, POINT *Eph, BASE_SYSTEM bse)
{
  mrand(ephm);
  elptic_mul(Eph, bse.Base, ephm, bse.E);
}

/*  NIST avf(x) function. Take bottom half of x value ored with
    2^(f/2) where f is number of bits in point order.
*/

void avf(mpz_t z, mpz_t x, BASE_SYSTEM bse)
{
  long f, i;
  mpz_t mask;

  mpz_init(mask);
  f = (mpz_sizeinbase(bse.order, 2) >> 1) + 1;
  for(i=0; i<f; i++)
    mpz_setbit(mask, i);
  mpz_and(z, x, mask);
  mpz_setbit(z, f);
  mpz_clear(mask);
}

/*  Compute MQV key share using my private and public key, 
    my ephemeral private and public key, their public key
    and ephemeral public key along with base system.
*/

void mqv_share(mpz_t keyshare, mpz_t my_key, POINT MY_KEY,
	       mpz_t my_ephm, POINT My_Ephm,
	       POINT Their_key, POINT Their_Ephm,
	       BASE_SYSTEM bse)
{
  POINT U;
  mpz_t s, z;

/*  compute s = k + R.x * sk  */

  mpz_inits(s, z, NULL);
  avf(z, My_Ephm.x, bse);
  mod_mul(s, z, my_key, bse.order);
  mod_add(s, s, my_ephm, bse.order);

/*  compute U = R' + R.x' * P'  */

  point_init(&U);
  avf(z, Their_Ephm.x, bse);
  elptic_mul(&U, Their_key, z, bse.E);
  elptic_sum(&U, U, Their_Ephm, bse.E);

/*  compute key share value = r*s*U (x component)  */

  elptic_mul(&U, U, s, bse.E);
  if(bse.cofactor > 1)
  {
    mpz_set_ui(z, bse.cofactor);
    elptic_mul(&U, U, z, bse.E);
  }
  mpz_set(keyshare, U.x);
  point_clear(&U);
  mpz_clears(s, z, NULL);
}

/*  initialize and clear signature struct.  */

void sig_init(ECDSA *sig)
{
  mpz_inits(sig->c, sig->d, NULL);
}

void sig_clear(ECDSA *sig)
{
  mpz_clears(sig->c, sig->d, NULL);
}

/* Elliptic Curve Digital Signature Algorithm (ECDSA)  
   Enter with public key, secret key, message, and
   base system.  Returns signature struct (which must be
   pre-initialized).
*/

void ecdsa_sign(ECDSA *sig, mpz_t sk, POINT Pk,
	       unsigned char *msg, long len, BASE_SYSTEM base)
{
  mpz_t k, e, tmp;
  POINT R;
  
  mpz_inits(e, k, tmp, NULL);
  hash(e, msg, len, base.order);        // hash message to prime
  mrand(k);
  point_init(&R);
  elptic_mul(&R, base.Base, k, base.E);
  mpz_mod(sig->c, R.x, base.order);  // translate R.x to base point order
  mod_mul(tmp, sk, sig->c, base.order);
  mod_add(tmp, tmp, e, base.order);
  mod_div(sig->d, tmp, k, base.order);   // d = (e + sk*c)/k
  point_clear(&R);
  mpz_clears(e, k, tmp, NULL);
}

/* ECDSA verify routine.
   Enter with signature struct, base system and signers public key
   along with the message that was signed.
   returns 1 on success, 0 on failure.
*/

int ecdsa_verify(ECDSA sig, POINT Pk, unsigned char *msg, long len, 
		BASE_SYSTEM base)
{
  mpz_t h, h1, h2, e;
  POINT R, S, T;
  int rtn;
  
  mpz_inits(h, h1, h2, e, NULL);
  hash(e, msg, len, base.order);        // hash message to prime
  mpz_invert(h, sig.d, base.order);
  mod_mul(h1, e, h, base.order);
  mod_mul(h2, sig.c, h, base.order);
  point_init(&T);
  elptic_mul(&T, Pk, h2, base.E);
  point_init(&S);
  elptic_mul(&S, base.Base, h1, base.E);
  point_init(&R);
  elptic_sum(&R, T, S, base.E);
  mpz_mod(h, R.x, base.order);  // convert R.x mod n
  if(!mpz_cmp(h, sig.c))      // check it's good
    rtn = 1;
  else
    rtn = 0;
  point_clear(&R);
  point_clear(&S);
  point_clear(&T);
  mpz_clears(h, h1, h2, e, NULL);
  return rtn;
}

/*  initialize and clear signature struct.  */

void snr_init(SCHNORR *sig)
{
  point_init(&sig->Q);
  mpz_init(sig->s);
}

void snr_clear(SCHNORR *sig)
{
  point_clear(&sig->Q);
  mpz_clear(sig->s);
}

/* Schnorr signature, one of the versions.
   Similar to ECDSA, but output is a point + value.
   Enter with public key, secret key, message, and
   base system.  Returns signature struct (which must be
   pre-initialized).
*/

void schnorr_sign(SCHNORR *sig, mpz_t sk, POINT Pk,
	       unsigned char *msg, long len, BASE_SYSTEM base)
{
  mpz_t k, e, tmp;
  unsigned char *cat;
  int xsz, i;
  
  mpz_inits(e, k, tmp, NULL);
  mrand(k);
  elptic_mul(&sig->Q, base.Base, k, base.E);
  xsz = (mpz_sizeinbase(sig->Q.x, 16) + 1)/2;
  cat = (unsigned char*)malloc(len + 2*xsz + 2);
  mpz_export(cat, NULL, -1, 1, 0, 0, sig->Q.x);
  mpz_export(&cat[xsz], NULL, -1, 1, 0, 0, sig->Q.y);
  for(i=0; i<len; i++)
    cat[2*xsz + i] = msg[i];
  hash(e, cat, len+2*xsz, base.order);        // hash message to prime
  mod_mul(tmp, sk, e, base.order);
  mod_sub(sig->s, k, tmp, base.order);
  mpz_clears(e, k, tmp, NULL);
}

/*  Schnorr verify.
    Enter with signature, public key, message, length and
    system.  Returns 1 if pass, 0 if fail.
*/

int schnorr_verify(SCHNORR sig, POINT Pk, unsigned char *msg,
		   long len, BASE_SYSTEM base)
{
  mpz_t e;
  unsigned char *cat;
  int xsz, i;
  POINT U, V, Qck;
  
  mpz_init(e);
  xsz = (mpz_sizeinbase(sig.Q.x, 16) + 1)/2;
  cat = (unsigned char*)malloc(len + 2*xsz + 2);
  mpz_export(cat, NULL, -1, 1, 0, 0, sig.Q.x);
  mpz_export(&cat[xsz], NULL, -1, 1, 0, 0, sig.Q.y);
  for(i=0; i<len; i++)
    cat[2*xsz + i] = msg[i];
  hash(e, cat, len+2*xsz, base.order);        // hash message to prime
  point_init(&U);
  point_init(&V);
  point_init(&Qck);
  elptic_mul(&U, base.Base, sig.s, base.E);
  elptic_mul(&V, Pk, e, base.E);
  elptic_sum(&Qck, U, V, base.E);  // compute Q from signature

/*  verify computed Q matches signature Q  */
  
  if((!mpz_cmp(sig.Q.x, Qck.x)) && (!mpz_cmp(sig.Q.y, Qck.y)))
    i = 1;
  else
    i = 0;
  mpz_clear(e);
  point_clear(&U);
  point_clear(&V);
  point_clear(&Qck);
  return i;
}
