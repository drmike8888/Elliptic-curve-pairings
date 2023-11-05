/************************************************************************
 *                                                                      *
 *     Eliptic curve utilities header for field extension curves.       *
 *                                                                      *
 ***********************************************************************/

#include "poly.h"

typedef struct
{
  POLY   x;
  POLY   y;
}POLY_POINT;

typedef struct
{
  POLY  a4;
  POLY  a6;
}POLY_CURVE;

void poly_point_init(POLY_POINT *P);
void poly_point_clear(POLY_POINT *P);
void poly_curve_init(POLY_CURVE *E);
void poly_curve_clear(POLY_CURVE *E);
void poly_point_copy(POLY_POINT *R, POLY_POINT P);
int poly_test_point(POLY_POINT P);
void poly_point_printf(char *str, POLY_POINT P);
void poly_curve_printf(char *str, POLY_CURVE E);
void poly_fofx(POLY *f, POLY x, POLY_CURVE E);
void FF_bump(POLY *x);
void poly_elptic_embed(POLY_POINT *P1, POLY_POINT *P2, POLY x, POLY_CURVE E);
void poly_elptic_sum(POLY_POINT *R, POLY_POINT P, POLY_POINT Q, POLY_CURVE E);
void poly_elptic_mul(POLY_POINT *Q, POLY_POINT P, mpz_t k, POLY_CURVE E);
void poly_point_rand(POLY_POINT *P, POLY_CURVE E);
