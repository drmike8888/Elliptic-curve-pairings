/**************************************************************
 *                                                            *
 *   Using data from curve search program, create a file for  *
 *   each field size for use as system level parameters.      *
 *                                                            *
 *************************************************************/

#include "modulo.h"
#include "eliptic.h"

/*  look up prime value and return it.  return 0 on error  
      */

void getprm(mpz_t prm, int size)
{
  switch(size)
  {
    case 160:
      mpz_set_str(prm, "ac000000000000000000000000000000000000001", 16);
      break;
    case 256:
      mpz_set_str(prm, "2b000000000000000000000000000000000000000000000000000000000000001", 16);
      break;
    case 384:
      mpz_set_str(prm, "2e00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", 16);
      break;
    case 512:
      mpz_set_str(prm, "e20000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", 16);
      break;
    default:
      printf("don't have prime for field size %d\n", size);
      mpz_set_ui(prm, 0);
      break;
  }
}

/*  look up curve parameters and associated cardinality.  Return
    all zeros if not valid size.  This initializes E and carde.    */

void getcurve(CURVE *E, mpz_t order, mpz_t cofac, int size)
{
  curve_init(E);
  switch(size)
  {
    case 160:
//1 782e 1 ac0000000000000000006543ba11adf8eb6345c77
      mpz_set_str(order, "ac0000000000000000006543ba11adf8eb6345c77", 16);
      mpz_set_ui(cofac, 1);
      mpz_set_ui(E->a4, 1);
      mpz_set_ui(E->a6, 0x782e);
      break;
    case 256:
//1  a87 1 2b0000000000000000000000000000002e7f521c85bba055a6e2161b956a47f69
      mpz_set_str(order, "2b0000000000000000000000000000002e7f521c85bba055a6e2161b956a47f69", 16);
      mpz_set_ui(cofac, 1);
      mpz_set_ui(E->a4, 1);
      mpz_set_ui(E->a6, 0xa87);
      break;
    case 384:
//1  310 1 2e00000000000000000000000000000000000000000000002275cc5f2f7fcc15352a2c993900a851b3a75365a9ac54733
      mpz_set_str(order, "2e00000000000000000000000000000000000000000000002275cc5f2f7fcc15352a2c993900a851b3a75365a9ac54733", 16);
      mpz_set_ui(cofac, 1);
      mpz_set_ui(E->a4, 1);
      mpz_set_ui(E->a6, 0x310);
      break;
    case 512:
//1   41 1 e2000000000000000000000000000000000000000000000000000000000000007788830d091dc57e3af7d7bbd15386ee9414602d88d1e6489cd056336922bbf4d

      mpz_set_str(order, "e2000000000000000000000000000000000000000000000000000000000000007788830d091dc57e3af7d7bbd15386ee9414602d88d1e6489cd056336922bbf4d", 16);
      mpz_set_ui(cofac, 1);
      mpz_set_ui(E->a4, 1);
      mpz_set_ui(E->a6, 0x41);
      break;
    default:
      printf("don't have prime for field size %d\n", size);
      mpz_set_ui(order, 0);
      mpz_set_ui(E->a4, 0);
      mpz_set_ui(E->a6, 0);
      break;
  }
}

int main(int argc, char *argv[])
{
  FILE *param;
  mpz_t prm, order, cofac;
  CURVE E;
  POINT Base;
  int i, list[4] = {160, 256, 384, 512};
  char filename[256];

  mpz_inits(prm, order, cofac, NULL);
  curve_init(&E);
  point_init(&Base);
  for(i=0; i<4; i++)
  {
    getprm(prm, list[i]);
    minit(prm);
    getcurve(&E, order, cofac, list[i]);
    point_rand(&Base, E);
    elptic_mul(&Base, Base, cofac, E);
    sprintf(filename, "Curve_%d_params.dat", list[i]);
    param = fopen(filename, "w");
    fprintf(param, "prime\n");
    gmp_fprintf(param, "%Zx\n", prm);
    fprintf(param, "order\n");
    gmp_fprintf(param, "%Zx\n", order);
    fprintf(param, "cofactor\n");
    gmp_fprintf(param, "%Zd\n", cofac);
    fprintf(param, "curve(a4   a6)\n");
    gmp_fprintf(param, "%Zx\n", E.a4);
    gmp_fprintf(param, "%Zx\n", E.a6);
    fprintf(param, "basepoint(x   y)\n");
    gmp_fprintf(param, "%Zx\n", Base.x);
    gmp_fprintf(param, "%Zx\n", Base.y);
    fclose(param);
  }
}  
  
