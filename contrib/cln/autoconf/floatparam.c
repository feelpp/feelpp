/* Determine some float parameters, much like gcc's "enquire.c". */
/* Bruno Haible 24.8.1996 */

/* This program expects to be compiled by an ANSI C or C++ compiler. */

#include <stdio.h>

typedef int boolean;
#define TRUE  1
#define FALSE 0

#ifdef HAVE_LONGDOUBLE
typedef long double ldouble;
#endif

static void header (void)
{
  printf("/* Rounding modes, for use below */\n");
  printf("#define rounds_to_nearest        0  /* 0.5 ulp */\n");
  printf("#define rounds_to_zero           1  /* 1 ulp */\n");
  printf("#define rounds_to_infinity       2  /* 1 ulp */\n");
  printf("#define rounds_to_minus_infinity 3  /* 1 ulp */\n");
  printf("\n");
}

#define check(type,typeprefix,typestr,equalfn,mainfn)  \
static boolean equalfn (volatile type* x, volatile type* y);		\
static void mainfn (void)						\
{									\
  int mant_bits;							\
  int epsilon_bits = -1;						\
  int negepsilon_bits = -1;						\
  { type x = 1.0; type y; type z;					\
    for (y = 1.0; ; y = 0.5*y)						\
      { z = x + y; if (equalfn(&x,&z)) break;				\
        z = z - x; if (!equalfn(&y,&z)) break;				\
        epsilon_bits++;							\
  }   }									\
  { type x = 1.0; type y; type z;					\
    for (y = -1.0; ; y = 0.5*y)						\
      { z = x + y; if (equalfn(&x,&z)) break;				\
        z = z - x; if (!equalfn(&y,&z)) break;				\
        negepsilon_bits++;						\
  }   }									\
  printf("/* Properties of type `%s': */\n",typestr);			\
  printf("/* Largest n for which 1+2^(-n) is exactly represented is %d. */\n",epsilon_bits); \
  printf("/* Largest n for which 1-2^(-n) is exactly represented is %d. */\n",negepsilon_bits); \
  if (negepsilon_bits <= epsilon_bits)					\
    { printf("#error \"No exponent jump at 1.0 for type %s!\"\n",typestr); \
      mant_bits = -1;							\
    }									\
  else									\
    { if (negepsilon_bits > epsilon_bits+1)				\
        printf("/* Base for type `%s' is 2^%d\n",typestr,negepsilon_bits-epsilon_bits); \
      mant_bits = epsilon_bits+1;					\
      printf("#define %s_mant_bits %d\n",typeprefix,mant_bits);		\
    }									\
  { int i; type x, y1, y2, ys1, ys2, z1, z2, zs1, zs2;			\
    x = 1.0; for (i = 0; i < epsilon_bits; i++) { x = 0.5*x; }		\
    y1 = 1.0 + 5.0*x; y2 = 1.0 + 6.0*x;					\
    ys1 = 1.0 + 5.4*x; ys2 = 1.0 + 5.6*x;				\
    z1 = -1.0 + (-5.0)*x; z2 = -1.0 + (-6.0)*x;				\
    zs1 = -1.0 + (-5.4)*x; zs2 = -1.0 + (-5.6)*x;			\
    if (equalfn(&ys1,&y1) && equalfn(&ys2,&y2) && equalfn(&zs1,&z1) && equalfn(&zs2,&z2)) \
      printf("#define %s_rounds rounds_to_nearest\n",typeprefix);	\
    else if (equalfn(&ys1,&y1) && equalfn(&ys2,&y1) && equalfn(&zs1,&z1) && equalfn(&zs2,&z1)) \
      printf("#define %s_rounds rounds_to_zero\n",typeprefix);		\
    else if (equalfn(&ys1,&y2) && equalfn(&ys2,&y2) && equalfn(&zs1,&z1) && equalfn(&zs2,&z1)) \
      printf("#define %s_rounds rounds_to_infinity\n",typeprefix);	\
    else if (equalfn(&ys1,&y1) && equalfn(&ys2,&y1) && equalfn(&zs1,&z2) && equalfn(&zs2,&z2)) \
      printf("#define %s_rounds rounds_to_minus_infinity\n",typeprefix); \
    else								\
      printf("#error \"Unknown rounding mode for type %s!\"\n",typestr); \
  }									\
  printf("\n");								\
}									\
static boolean equalfn (volatile type* x, volatile type* y)		\
{									\
  return *x == *y;							\
}									\

check(float,"float","float",equal_float,main_float)
check(double,"double","double",equal_double,main_double)
#ifdef HAVE_LONGDOUBLE
check(ldouble,"long_double","long double",equal_ldouble,main_ldouble)
#endif

/* Some systems (arm/linux) store doubles as little endian but with higher
 * and lower word reversed. */
static void flipped_double (void)
{
  typedef struct { unsigned lo, hi; } dfloat;
  union { dfloat eksplicit; double machine_double; } x;
  x.machine_double = 2;
  dfloat test = x.eksplicit;
  if (test.lo==0 && test.hi!=0) {
    printf("#define double_wordorder_bigendian_p 0\n");
  } else if (test.lo!=0 && test.hi==0) {
    printf("#define double_wordorder_bigendian_p 1\n");
  } else {
    /* Dazed and confused!  Better not define anything.
	 * Code should rely on CL_CPU_BIG_ENDIAN_P instead. */
  }
  printf("\n");
}
	 
int main()
{
  header();
  main_float();
  main_double();
#ifdef HAVE_LONGDOUBLE
  main_ldouble();
#endif
  flipped_double();

  if (ferror(stdout) || fclose(stdout)) return 1;
  return 0;
}
