/********************************************************
 *                                                      *
 *    Use GMP to compute elliptic curve addition and    *
 *    multiplication.  Set up basic structures to make  *
 *    this easy to follow.                              *
 *                                                      *
 *******************************************************/

typedef struct
{
  mpz_t   x;
  mpz_t   y;
}POINT;

typedef struct
{
  mpz_t  a4;
  mpz_t  a6;
}CURVE;

void point_init(POINT *P);
void point_clear(POINT *P);
void curve_init(CURVE *E);
void curve_clear(CURVE *E);
void point_copy(POINT *R, POINT P);
int test_point(POINT P);
void fofx(mpz_t f, mpz_t x, CURVE E);
void elptic_sum(POINT *R, POINT P, POINT Q, CURVE E);
void elptic_embed(POINT *P1, POINT *P2, mpz_t x, CURVE E);
void point_printf(char *str, POINT P);
void elptic_mul(POINT *Q, POINT P, mpz_t k, CURVE E);
void point_rand(POINT *P, CURVE E);
