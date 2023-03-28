/**********************************************************************
 *                                                                    *
 *    For primes out of Reisel book, find curves with large primes    *
 *    and high embedding degree.                                      *
 *                                                                    *
 *********************************************************************/

#include <pari/pari.h>

/*
  primes from Reisel for fields of 160, 256, 384 and 512:
  43*2^158+1   0xac000000000000000000000000000000000000001
  43*2^252+1   0x2b000000000000000000000000000000000000000000000000000000000000001
  23*2^381+1   0x2e00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
  113*2^509+1  0xe20000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
*/

int main()
{
  GEN y, E, f, z, a1, a3, a2, ell5, num, tmp;
  int k, pw2, pw[4];
  unsigned long a4coef, a6coef;
  pari_sp av;

/* initialize pari and finite field  */

  pari_init(1024*1024*1024, 5*1024*512);
  z = gp_read_str("0xac000000000000000000000000000000000000001");
  y = ffgen(z, -1);

  a1 = gen_0;
  a3 = gen_0;
  a2 = gen_0;
  ell5 = zerovec( 5);
  gel(ell5, 1) = a1;
  gel(ell5, 2) = a2;
  gel(ell5, 3) = a3;
  a4coef = 1;
  a6coef = 0x40;

  num = muluu(3, 5);
  num = mului(7, num);
  num = mului(11, num);
  while(a6coef < 0xfffff)
  {
    printf("%0lx  %0lx", a4coef, a6coef);
    fflush(stdout);
    av = avma;
    gel(ell5, 4) = stoi(a4coef);
    gel(ell5, 5) = stoi(a6coef);
    E = ellinit(ell5, y, 0);

    f = ellcard(E, NULL);
    k = bittest(f, 0);
    pw2 = 0;
    while(!k)
    {
      pw2++;
      f = gdivexact(f, gen_2);
      k = bittest(f, 0);
    }
    tmp = gcdii(f, num);
    if(isint1(tmp))
    {
      if(isprime(f))
	pari_printf(" 2^%d * %Ps\n", pw2, f);
      else
	printf("\n");
    }
    else
    {
      pw[0] = itos(tmp);
      k = 1;
      f = gdivexact(f, tmp);
      while((k < 4) && !isint1(tmp))
      {
	tmp = gcdii(f, num);
	if(isint1(tmp))
	{
	  if(isprime(f))
	  {
	    printf(" 2^%d * %d ", pw2, pw[0]);
	    k--;
	    while(k)
	    {
	      printf("* %d  ", pw[k]);
	      k--;
	    }
	    pari_printf("* %Ps\n", f);
	  }
	  else
	    printf("\n");
	}
	else
	{
	  pw[k] = itos(tmp);
	  f=gdivexact(f, tmp);
	  k++;
	}
      }
      if(k >= 4)
	printf("\n");
    }
    gerepileall( av, 1, &ell5);
    a6coef++;
  }
  printf("all done\n");
  return 0;
}
