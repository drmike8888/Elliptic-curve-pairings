#include "modulo.h"
#include "eliptic.h"
#include "k12.header/KangarooTwelve.h"

typedef struct
{
  mpz_t c;
  mpz_t d;
}ECDSA;

typedef struct
{
  POINT Q;
  mpz_t s;
}SCHNORR;
  
typedef struct
{
  long  cofactor;
  mpz_t order;
  POINT Base;
  CURVE E;
}BASE_SYSTEM;

void hash(mpz_t hsh, unsigned char *dat, long len, mpz_t prm);
void gen_key(mpz_t sk, POINT *PK, unsigned char *phrase, BASE_SYSTEM bse);
void diffie_hellman(mpz_t keyshare, mpz_t my_key, POINT Their_key, CURVE E);
void mqv_ephem(mpz_t ephm, POINT *Eph, BASE_SYSTEM bse);
void mqv_share(mpz_t keyshare, mpz_t my_key, POINT MY_KEY,
	       mpz_t my_ephm, POINT My_Ephm,
	       POINT Their_key, POINT Their_Ephm,
	       BASE_SYSTEM bse);
void sig_init(ECDSA *sig);
void sig_clear(ECDSA *sig);
void ecdsa_sign(ECDSA *sig, mpz_t sk, POINT Pk,
	       unsigned char *msg, long len, BASE_SYSTEM base);
int ecdsa_verify(ECDSA sig, POINT Pk, unsigned char *msg, long len, 
		BASE_SYSTEM base);
void snr_init(SCHNORR *sig);
void snr_clear(SCHNORR *sig);
void schnorr_sign(SCHNORR *sig, mpz_t sk, POINT Pk,
		  unsigned char *msg, long len, BASE_SYSTEM base);
int schnorr_verify(SCHNORR sig, POINT Pk, unsigned char *msg,
		   long len, BASE_SYSTEM base);
