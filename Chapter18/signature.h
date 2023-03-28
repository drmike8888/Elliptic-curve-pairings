#include "poly_eliptic.h"
#include "eliptic.h"
#include "pairing.h"
#include "k12.header/KangarooTwelve.h"

void keygen(mpz_t sk, POLY_POINT *PK, POLY_POINT G2, POLY_CURVE Ex,
	    unsigned char *ptr, int len);
void hash1(mpz_t hsh, unsigned char *dat, long len, mpz_t r);
void hash0(POINT *H, SIG_SYSTEM sig, unsigned char *dat, long len);
void hash2(POINT *H, SIG_SYSTEM sig, unsigned char *dat, long len);
void point2text(unsigned char *ptr, POLY_POINT P, long msz, long deg);
void aj_hash(mpz_t *ajhsh, SIG_SYSTEM sig, POLY_POINT *PK,
	     long numkey);
void get_system(char *filename, SIG_SYSTEM *sig);
void tog2(POLY_POINT *G2, POINT G1);
void sign(POINT *S, SIG_SYSTEM sig, mpz_t sk, mpz_t aj,
	   unsigned char *msg, long msgln);
void agregat_sig(POINT *Sigma, SIG_SYSTEM sig, POINT *S, long nmsig);
void aj_sum(POLY_POINT *APK, SIG_SYSTEM sig, POLY_POINT *PK,
	    mpz_t *ajhsh, long nmkey);
int multisig_verify(SIG_SYSTEM sig, POINT Sigma, POLY_POINT APK,
		    unsigned char *msg, long msgln);
void mu_column(POINT *muvec, POLY_POINT APK, SIG_SYSTEM sig,
	       mpz_t aj, mpz_t sk, long nmkey);
void subgrp_sign(POINT *S, POLY_POINT APK, SIG_SYSTEM sig,
		 POINT Mk, mpz_t sk, unsigned char *msg, long msgln);
void subgrp_combine(POLY_POINT *PK, POINT *Ssum, SIG_SYSTEM sig,
		    POINT *Svec, POLY_POINT *Pubky, long *list, long nmlst);
int subgrp_verify(SIG_SYSTEM sig, POLY_POINT APK, unsigned char *msg, long msgln,
		  long *list, long nmlst, POLY_POINT PK, POINT Ssum);
void membership_key(POINT *Memkey, POINT *Murow, long numkey, CURVE E);
void point_write(POINT *P, FILE *f);
void point_read(POINT *P, FILE *f);
void poly_write(POLY *p, FILE *f);
void poly_read(POLY *p, FILE *f);
void poly_point_write(POLY_POINT *P, FILE *f);
void poly_point_read(POLY_POINT *P, FILE *f);
