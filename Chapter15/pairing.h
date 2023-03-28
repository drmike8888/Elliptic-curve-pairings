/*  pairing related function and structures  */

typedef struct
{
  mpz_t prime;
  CURVE E;
  mpz_t cardE;
  mpz_t tor;
  mpz_t cobse;
  POINT G1;
  POLY irrd;
  mpz_t cardEx;
  mpz_t coxtd;
  POLY_POINT G2;
  POLY_CURVE Ex;
} SIG_SYSTEM;

long g1g2(POLY_POINT P);
void cardinality(mpz_t crd, mpz_t t, long k);
void weil(POLY *w, POLY_POINT P, POLY_POINT Q, POLY_POINT S, mpz_t m, POLY_CURVE E);
int get_order(mpz_t order, POINT P, CURVE E, mpz_t *factors, int n);
int poly_get_order(mpz_t order, POLY_POINT P, POLY_CURVE E, mpz_t *factors, int n);
void tate(POLY *t, POLY_POINT P, POLY_POINT Q, POLY_POINT S, mpz_t m, POLY_CURVE E);

