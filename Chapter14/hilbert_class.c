/************************************************************
 *                                                          *
 *     Print out Hilbert Class Polynomials of j(tau) using  *
 *   pari/gp function polclass().  The values of D come     *
 *   from the values used as alpha in the GMP routines      *
 *   which look for pairing friendly curves.                *
 *                                                          *
 ***********************************************************/

#include <pari/pari.h>

int main(int argc, char *argv[])
{
  FILE *hlbrt;
  GEN D;
  long atab[5] = {7, 11, 15, 19, 23};
  int i, j;
  pari_sp av;
   
  pari_init(2*1024*1024*512, 5*1024*512);

  hlbrt = fopen("Hilbert_Polynomials.list", "w");
  for(j=0; j<8; j++)
  {
    for(i=0; i<5; i++)
    {
      D = stoi(-(atab[i] + j*20));
      pari_fprintf(hlbrt, "%Ps : %Ps\n", D, polclass(D, 0 , -1));
    }
  }
  fclose(hlbrt);
}
