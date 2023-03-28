/*******************************************************************
 *                                                                 *
 *      Use values from pairings_sweep_alpha to find Hilbert Class *
 *    Polynomial and j value which gives curve of correct order.   *
 *    Input is alpha, q and t.  Outputs j roots of HCP and curve   *
 *    parameters a, b.                                             *
 *                                                                 *
 ******************************************************************/

#include "poly.h"
#include "eliptic.h"
#include <string.h>

/*  compute roots of second order polynomial.
    Enter with 2nd order polynomial hc
    Returns 2 roots j1 and j2. (which should be
    pre-initialized)  */

void tworoots(mpz_t j1, mpz_t j2, POLY hc)
{
  mpz_t d, e, f;

  if(hc.deg != 2)
    return;
  if(mpz_cmp_ui(hc.coef[2], 1))
  {
    mdiv(hc.coef[1], hc.coef[1], hc.coef[2]);
    mdiv(hc.coef[0], hc.coef[0], hc.coef[2]);
  }
  mpz_inits(d, e, NULL);
  mmul(d, hc.coef[1], hc.coef[1]);
  mpz_init_set_ui(f, 4);
  mmul(e, hc.coef[0], f);
  msub(d, d, e);
  msqrt(e, d);
  mpz_neg(j1, hc.coef[1]);
  mpz_neg(j2, hc.coef[1]);
  madd(j1, j1, e);
  msub(j2, j2, e);
  mpz_set_ui(f, 2);
  mdiv(j1, j1, f);
  mdiv(j2, j2, f);
  mpz_clears(d, e, f, NULL);
}

int main(int argc, char *argv[])
{
  FILE *hilb;
  POLY hc, stack[16], x1, hofx, Aofx, Bofx, rem;
  int i, j, k, alpha, xs, xe, sign, stckp, done;
  mpz_t root[16], p, c, a4[16], a6[16], tp, b, j0;
  mpz_t two, tri, t, crde;
  char *hcdat;
  POINT R, P0;
  CURVE E;
  
  if(argc < 4)
  {
    printf("Use: ./get_curve <discriminant> <prime> <t>\n");
    printf("   values from output of parings_alpha.k\n");
    exit(-1);
  }

/*  read in Hilbert Class Polynomials  */
  
  hilb = fopen("Hilbert_Polynomials.list", "r");
//    hilb = fopen("Hilbert_Polynomials.test", "r");
  if(!hilb)
  {
    printf("can't find file Hilbert_Polynomials.list\n");
    exit(-2);
  }
  hcdat = (char*)malloc(6*1024);
  k = 0;
  while(!feof(hilb))
  {
    hcdat[k] = fgetc(hilb);
    k++;
  }
  fclose(hilb);
  k -= 2;
  
/*  see if discriminat in list  */

  alpha = atol(argv[1]);
  if(alpha < 0)
    alpha *= -1;
  i = 0;
  j = 0;
  while((j < alpha) && (i < k))
  {
    sscanf(&hcdat[i], "%d", &j);
    j *= -1;
    if(j == alpha)
      break;
    while(hcdat[i] != '\n')
      i++;
    i++;
  }
  if(j != alpha)
  {
    printf("invalid discriminant %d\n", alpha);
    exit(-3);
  }

/*  Set up mod q operations  */

  if(mpz_init_set_str(p, argv[2], 10) < 0)
  {
    printf("invalid prime string\n");
    exit(-4);
  }
  minit(p);

/*  input factor t = p + 1 - #E  */

  if(mpz_init_set_str(t, argv[3], 10) < 0)
  {
    printf("invalid t string\n");
    exit(-5);
  }

/*  read in hilbert polynomial  */

  poly_init(&hc);
  xs = i;
  while(hcdat[xs] != 'x') xs++;
  if(hcdat[xs + 1] == '^')
  {
    xs += 2;
    sscanf(&hcdat[xs], "%ld", &hc.deg);
  }
  else
    hc.deg = 1;
  mpz_set_ui(hc.coef[hc.deg], 1);
  j = hc.deg - 1;
  while(j >= 0)
  {
    while((hcdat[xs] != '+') && (hcdat[xs] != '-'))
      xs++;
    if(hcdat[xs] == '-')
      sign = -1;
    else
      sign = 1;
    xs++;
    xe = xs+1;
    while((hcdat[xe] != '*') && (hcdat[xe] != '\n'))
      xe++;
    hcdat[xe] = 0;
    mpz_set_str(hc.coef[j], &hcdat[xs], 10);
    if(sign < 0)
      mpz_neg(hc.coef[j], hc.coef[j]);
    mpz_mod(hc.coef[j], hc.coef[j], p);   // ensure size mod p
    j--;
  }

/*  for order 1 and 2, output result directly  */

  stckp = 0;
  if(hc.deg < 3)
  {
    if(hc.deg == 1)
    {
      mpz_init_set(root[0], hc.coef[0]);
      mpz_neg(root[0], root[0]);
    }
    else
    {
      mpz_inits(root[0], root[1], NULL);
      tworoots(root[0], root[1], hc);
    }
  }
  else
  {

/* find h(x)  = x^p mod hc(x)  */

    poly_mulprep(hc);
    poly_init(&x1);
    x1.deg = 1;
    mpz_set_ui(x1.coef[1], 1);
    poly_init(&hofx);
    poly_xp(&hofx, x1);
    poly_sub(&hofx, hofx, x1);
    
/*  A(x) = gcd(x^p - x, hc(x))  */
    
    poly_init(&Aofx);
    poly_gcd(&Aofx, hc, hofx);
    if(!Aofx.deg)
    {
      printf("no roots found for this combination:\n");
      poly_print(hc);
      gmp_printf("prime: %Zd\n", p);
      exit(-5);
    }
    poly_init(&Bofx);
    poly_init(&rem);
    
/*  push first A(x) on stack  */

    for(i=0; i<16; i++)
      poly_init(&stack[i]);
    poly_copy(&stack[stckp], Aofx);

/* initialize all roots to zero  */
  
    for(i=0; i<hc.deg; i++)
      mpz_init(root[i]);
    j = 0;                 // root index counter
    stckp++;
  }
  while(stckp)
  {
/*  pop next A(x) off stack */

    stckp--;
    poly_copy(&Aofx, stack[stckp]);
    done = 0;
    while(!done)
    {
      if(!Aofx.coef[0])         // is A(0) = 0?
      {
	j++;                    // this root is zero
	for(i=1; i<=Aofx.deg; i++)
	  mpz_set(Aofx.coef[i - 1], Aofx.coef[i]);  // A(x)/x
	Aofx.deg--;
	if(!Aofx.deg)
	  done = 1;
      }
      if(Aofx.deg == 1)
      {
	if(!mpz_cmp_ui(Aofx.coef[1], 1))
	  mpz_set(root[j], Aofx.coef[0]);
	else
	  mdiv(root[j], Aofx.coef[0], Aofx.coef[1]);
	mneg(root[j], root[j]);
	j++;
	done = 1;
      }
      else if(Aofx.deg == 2)
      {
	tworoots(root[j+1], root[j], Aofx);
	j += 2;
	done = 1;
      }
      if(!done)
      {
	Bofx.deg = 0;
	poly_mulprep(Aofx);
	while(!Bofx.deg || (Bofx.deg == Aofx.deg))
	{
	  mrand(x1.coef[0]);
	  gpow_p2(&hofx, x1);
	  mpz_sub_ui(hofx.coef[0], hofx.coef[0], 1);
	  poly_gcd(&Bofx, hofx, Aofx);
	}
	poly_copy(&stack[stckp], Bofx);
	stckp++;
	poly_euclid(&stack[stckp], &rem, Aofx, Bofx);
	stckp++;
	done = 1;
      }
    }
    if(stckp > 15) exit(0);
  }

/*  output all found roots  */

  for(i=0; i<hc.deg; i++)
    gmp_printf("%d:  %Zd\n", i, root[i]);

/*  Compute a4 and a6 for each root  */

  mpz_inits(c, tp, b, j0, two, tri, NULL);
  mpz_set_ui(j0, 1728);
  mpz_set_ui(two, 2);
  mpz_set_ui(tri, 3);
  for(i=0; i<hc.deg; i++)
  {
    mpz_inits(a4[i], a6[i], NULL);
    mpz_set(tp, root[i]);
    msub(b, j0, tp);
    mdiv(c, tp, b);
    mmul(a4[i], c, tri);
    mmul(a6[i], c, two);
    gmp_printf("%d: a4= %Zd  a6= %Zd\n", i, a4[i], a6[i]);
  }

/*  compute #E * random point on each curve
    if we get point at infinity, this the curve we want.  */

  mpz_init_set(crde, p);
  mpz_add_ui(crde, crde, 1);
  mpz_sub(crde, crde, t);
  gmp_printf("#E = %Zd\n", crde);
  point_init(&R);
  curve_init(&E);
  point_init(&P0);
  done = 0;
  for(i=0; i<hc.deg; i++)
  {
    mpz_set(E.a4, a4[i]);
    mpz_set(E.a6, a6[i]);
    point_rand(&R, E);
    elptic_mul(&P0, R, crde, E);
    if(test_point(P0))
    {
      printf("curve %d is right curve!\n", i);
      done = 1;
    }
    else
      printf("curve %d is not right curve.\n", i);
  }

/*  if no curve is right, compute twist of first one and try again.  */
  
  if(!done)
  {
    k = 1;
    while(k >= 0)
    {
      mrand(c);
      k = msqr(c);
    }
    mmul(b, c, c);
    mmul(E.a4, a4[0], b);  // c^2 * a4
    mmul(b, b, c);
    mmul(E.a6, a6[0], b);  // c^3 * a6
    point_rand(&R, E);
    elptic_mul(&P0, R, crde, E);
    if(test_point(P0))
    {
      gmp_printf("a4 = %Zd  a6 = %Zd\n", E.a4, E.a6);
      printf("twist is right curve!\n");
    }
    else
      printf("twist is not right curve.\n");
  }
  mpz_clears(c, tp, b, j0, two, tri, NULL);
  mpz_clear(p);
  poly_clear(&x1);
  poly_clear(&Aofx);
  poly_clear(&Bofx);
  poly_clear(&hofx);
  for(i=0; i<=hc.deg; i++)
    mpz_clears(root[i], a4[i], a6[i], NULL);
  point_clear(&R);
  curve_clear(&E);
}

