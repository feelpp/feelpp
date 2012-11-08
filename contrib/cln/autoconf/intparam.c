/* Bestimmung einiger Maschinen-Parameter und -Abhängigkeiten */
/* und Ausgabe in ein Include-File */
/* Bruno Haible 10.9.1991, 12.10.1992, 6.12.1992, 24.10.1993 */

/* Auf einigen Systemen werden in <sys/types.h> die Typen uchar, ushort, uint, */
/* ulong definiert. Normalerweise reicht _POSIX_SOURCE aus, dies zu verhindern, */
/* bei AIX 3.2.5 (rs6000-ibm-aix3.2.5) jedoch nicht. Wir müssen Gewalt anwenden. */
#define _POSIX_SOURCE
#define uchar  os_uchar
#define ushort os_ushort
#define uint   os_uint
#define ulong  os_ulong
#include <stdio.h>
#undef ulong
#undef uint
#undef ushort
#undef uchar

#define loop  while(1)

typedef int boolean;
#define TRUE  1
#define FALSE 0

#ifdef __CHAR_UNSIGNED__
typedef signed char  schar;
#else
typedef char  schar;
#endif
typedef unsigned char  uchar;
typedef unsigned short  ushort;
typedef unsigned /* int */  uint;
typedef unsigned long  ulong;
#ifdef HAVE_LONGLONG
typedef long long  longlong;
typedef unsigned long long  ulonglong;
#endif
typedef int (function)();

static int random_table[256] = /* 2048 zufällige Bits, hier von pi */
  { 0xC9,0x0F,0xDA,0xA2,0x21,0x68,0xC2,0x34,0xC4,0xC6,0x62,0x8B,
    0x80,0xDC,0x1C,0xD1,0x29,0x02,0x4E,0x08,0x8A,0x67,0xCC,0x74,
    0x02,0x0B,0xBE,0xA6,0x3B,0x13,0x9B,0x22,0x51,0x4A,0x08,0x79,
    0x8E,0x34,0x04,0xDD,0xEF,0x95,0x19,0xB3,0xCD,0x3A,0x43,0x1B,
    0x30,0x2B,0x0A,0x6D,0xF2,0x5F,0x14,0x37,0x4F,0xE1,0x35,0x6D,
    0x6D,0x51,0xC2,0x45,0xE4,0x85,0xB5,0x76,0x62,0x5E,0x7E,0xC6,
    0xF4,0x4C,0x42,0xE9,0xA6,0x37,0xED,0x6B,0x0B,0xFF,0x5C,0xB6,
    0xF4,0x06,0xB7,0xED,0xEE,0x38,0x6B,0xFB,0x5A,0x89,0x9F,0xA5,
    0xAE,0x9F,0x24,0x11,0x7C,0x4B,0x1F,0xE6,0x49,0x28,0x66,0x51,
    0xEC,0xE4,0x5B,0x3D,0xC2,0x00,0x7C,0xB8,0xA1,0x63,0xBF,0x05,
    0x98,0xDA,0x48,0x36,0x1C,0x55,0xD3,0x9A,0x69,0x16,0x3F,0xA8,
    0xFD,0x24,0xCF,0x5F,0x83,0x65,0x5D,0x23,0xDC,0xA3,0xAD,0x96,
    0x1C,0x62,0xF3,0x56,0x20,0x85,0x52,0xBB,0x9E,0xD5,0x29,0x07,
    0x70,0x96,0x96,0x6D,0x67,0x0C,0x35,0x4E,0x4A,0xBC,0x98,0x04,
    0xF1,0x74,0x6C,0x08,0xCA,0x18,0x21,0x7C,0x32,0x90,0x5E,0x46,
    0x2E,0x36,0xCE,0x3B,0xE3,0x9E,0x77,0x2C,0x18,0x0E,0x86,0x03,
    0x9B,0x27,0x83,0xA2,0xEC,0x07,0xA2,0x8F,0xB5,0xC5,0x5D,0xF0,
    0x6F,0x4C,0x52,0xC9,0xDE,0x2B,0xCB,0xF6,0x95,0x58,0x17,0x18,
    0x39,0x95,0x49,0x7C,0xEA,0x95,0x6A,0xE5,0x15,0xD2,0x26,0x18,
    0x98,0xFA,0x05,0x10,0x15,0x72,0x8E,0x5A,0x8A,0xAA,0xC4,0x2D,
    0xAD,0x33,0x17,0x0D,0x04,0x50,0x7A,0x33,0xA8,0x55,0x21,0xAB,
    0xDF,0x1C,0xBA,0x65,
  };
#define random_table_length  (8*256)
static int random_position = -1;
int next_random_bit(void)
  { random_position++;
    if (random_position==random_table_length) random_position = 0;
    return (random_table[random_position/8] >> (random_position % 8)) & 1;
  }

void printf_underscored (const char* string)
  { char c;
    while (!((c = *string++) == '\0')) { printf("%c",(c==' ' ? '_' : c)); }
  }

/* string_length(string) is the same as strlen(string). */
/* Better avoid to depend on <string.h>. */
int string_length (char* string)
  { int count = 0;
    while (!(*string++ == '\0')) { count++; }
    return count;
  }

static int char_bitsize, short_bitsize, int_bitsize, long_bitsize;
static int uchar_bitsize, ushort_bitsize, uint_bitsize, ulong_bitsize;
static boolean char_uchar_same, short_ushort_same, int_uint_same, long_ulong_same;
static int pointer_bitsize;
#ifdef HAVE_LONGLONG
static int longlong_bitsize, ulonglong_bitsize;
static boolean longlong_ulonglong_same;
#endif

void main1(void) {
#define get_unsigned_integer_bitsize(type,where)  \
  { type x = 1;                               \
    int bits = 0;                             \
    loop {                                    \
      if (x==0) break;                        \
      x = x+x;                                \
      bits++;                                 \
      if (bits==1000) { bits = -1; break; }   \
    }                                         \
    where = bits;                             \
  }
#define get_signed_integer_bitsize(type,unsigned_type,where)  \
  { /* Signed integer overflow is "undefined behaviour" in C99, and gcc-4.3    \
       (without -fwrapv option) actually does weird things when signed integer \
       overflow occurs. Therefore perform the addition on the unsigned type.   \
       Drawback: This will not detect cases where the signed type has more bits\
       than the unsigned type but the same size according to sizeof. Blech. */ \
    type x = 1;                               \
    int bits = 0;                             \
    loop {                                    \
      if (x==0) break;                        \
      x = (unsigned_type)x + (unsigned_type)x;\
      bits++;                                 \
      if (bits==1000) { bits = -1; break; }   \
    }                                         \
    where = bits;                             \
  }
#define print_integer_bitsize(type,typestr,where)  \
  { if (where >= 0)                                                  \
      { printf("/* Integers of type %s have %ld bits. */\n",typestr,(long)where); \
        if (!(typestr[0] == 'u'))                                    \
          { printf("#define "); printf_underscored(typestr); printf("_bitsize %ld\n",(long)where); } \
        printf("\n");                                                \
      }                                                              \
      else                                                           \
      { printf("#error \"Integers of type %s have no binary representation!!\"\n",typestr); } \
    if (!(where == char_bitsize * sizeof(type)))                     \
      { printf("#error \"Formula BITSIZE(T) = SIZEOF(T) * BITSPERBYTE does not hold for type %s!!\"\n",typestr); } \
  }
  get_signed_integer_bitsize(schar,uchar,char_bitsize);
  get_signed_integer_bitsize(short,ushort,short_bitsize);
  get_signed_integer_bitsize(int,uint,int_bitsize);
  get_signed_integer_bitsize(long,ulong,long_bitsize);
  print_integer_bitsize(schar,"char",char_bitsize);
  print_integer_bitsize(short,"short",short_bitsize);
  print_integer_bitsize(int,"int",int_bitsize);
  print_integer_bitsize(long,"long",long_bitsize);
#ifdef HAVE_LONGLONG
  get_signed_integer_bitsize(longlong,ulonglong,longlong_bitsize);
  print_integer_bitsize(longlong,"long long",longlong_bitsize);
#endif
  get_unsigned_integer_bitsize(uchar,uchar_bitsize);
  get_unsigned_integer_bitsize(ushort,ushort_bitsize);
  get_unsigned_integer_bitsize(uint,uint_bitsize);
  get_unsigned_integer_bitsize(ulong,ulong_bitsize);
  print_integer_bitsize(uchar,"unsigned char",uchar_bitsize);
  print_integer_bitsize(ushort,"unsigned short",ushort_bitsize);
  print_integer_bitsize(uint,"unsigned int",uint_bitsize);
  print_integer_bitsize(ulong,"unsigned long",ulong_bitsize);
#ifdef HAVE_LONGLONG
  get_unsigned_integer_bitsize(ulonglong,ulonglong_bitsize);
  print_integer_bitsize(ulonglong,"unsigned long long",ulonglong_bitsize);
#endif
}

void main2(void) {
#define compare_integer_bitsizes(typestr1,typestr2,type1_bitsize,type2_bitsize)  \
  { if (!(type1_bitsize==type2_bitsize))                                                       \
      printf("#error \"Integer types %s and %s have different sizes!!\"\n",typestr1,typestr2); \
  }
  compare_integer_bitsizes("char","unsigned char",char_bitsize,uchar_bitsize);
  compare_integer_bitsizes("short","unsigned short",short_bitsize,ushort_bitsize);
  compare_integer_bitsizes("int","unsigned int",int_bitsize,uint_bitsize);
  compare_integer_bitsizes("long","unsigned long",long_bitsize,ulong_bitsize);
#ifdef HAVE_LONGLONG
  compare_integer_bitsizes("long long","unsigned long long",longlong_bitsize,ulonglong_bitsize);
#endif
}

#define get_a_random(type,bitsize,where)  \
  { type x = 0;                                          \
    int i = bitsize;                                     \
    while (i>0) { x = (x<<1) + next_random_bit(); i--; } \
    where = x;                                           \
  }
#define get_a_random_twice(type1,type2,bitsize,where1,where2)  \
  { type1 x1 = 0; type2 x2 = 0;                 \
    int i = bitsize;                            \
    while (i>0)                                 \
      { type1 b = next_random_bit();            \
        x1 = ((x1<<1) + b); x2 = ((x2<<1) + b); \
        i--;                                    \
      }                                         \
    where1 = x1; where2 = x2;                   \
  }

void main3(void) {
#define compare_integer_representation(type1,type2,typestr1,typestr2,type1_bitsize,type2_bitsize,where)  \
  { if ((type1_bitsize>=0) && (type2_bitsize>=0) && (type1_bitsize==type2_bitsize))           \
      { int i,j;                                                                              \
        type1 sample1; type2 sample2;                                                         \
        where = TRUE;                                                                         \
        for (i = 0; i<100; i++)                                                               \
          { get_a_random_twice(type1,type2,type1_bitsize,sample1,sample2);                    \
            if (!(sample1 == (type1)(sample2))) { where = FALSE; }                            \
            if (!(sample2 == (type2)(sample1))) { where = FALSE; }                            \
          }                                                                                   \
        for (i = 0; i<100; i++)                                                               \
          { get_a_random(type1,type1_bitsize,sample1);                                        \
            sample2 = (type2)(sample1);                                                       \
            for (j = 0; j < type1_bitsize; j++)                                               \
              if (!( ((sample1 & ((type1)1<<j)) == 0)                                         \
                     == ((sample2 & ((type2)1<<j)) == 0)                                      \
                 ) )                                                                          \
                { where = FALSE; }                                                            \
          }                                                                                   \
        if (where)                                                                            \
          { printf("/* Integer types %s and %s have the same binary representation. */\n",typestr1,typestr2); } \
          else                                                                                \
          { printf("#error \"Integer types %s and %s have different binary representations!!\"\n",typestr1,typestr2); } \
      }                                                                                       \
      else                                                                                    \
      { where = FALSE; }                                                                      \
  }
  compare_integer_representation(schar,uchar,"char","unsigned char",char_bitsize,uchar_bitsize,char_uchar_same);
  compare_integer_representation(short,ushort,"short","unsigned short",short_bitsize,ushort_bitsize,short_ushort_same);
  compare_integer_representation(int,uint,"int","unsigned int",int_bitsize,uint_bitsize,int_uint_same);
  compare_integer_representation(long,ulong,"long","unsigned long",long_bitsize,ulong_bitsize,long_ulong_same);
#ifdef HAVE_LONGLONG
  compare_integer_representation(longlong,ulonglong,"long long","unsigned long long",longlong_bitsize,ulonglong_bitsize,longlong_ulonglong_same);
#endif
  printf("\n");
}

void main4(void) {
#define test_integer_ushift(type,typestr,type_bitsize)  \
  if (type_bitsize >= 0)                                                                        \
    { int i,j,shc;                                                                              \
      type sample1,sample2;                                                                     \
      boolean left_works = TRUE, right_works = TRUE;                                            \
      for (i = 0; i<100; i++)                                                                   \
        { get_a_random(type,type_bitsize,sample1);                                              \
          for (shc = 0; shc < type_bitsize; shc++)                                              \
            { sample2 = sample1 << shc;                                                         \
              for (j=0; j < type_bitsize; j++)                                                  \
                { if (!( ((sample2 & ((type)1<<j)) == 0)                                        \
                         ==                                                                     \
                         (j < shc ? TRUE : ((sample1 & ((type)1<<(j-shc))) == 0))               \
                     ) )                                                                        \
                    { left_works = FALSE; }                                                     \
        }   }   }                                                                               \
      for (i = 0; i<100; i++)                                                                   \
        { get_a_random(type,type_bitsize,sample1);                                              \
          for (shc = 0; shc < type_bitsize; shc++)                                              \
            { sample2 = sample1 >> shc;                                                         \
              for (j=0; j < type_bitsize; j++)                                                  \
                { if (!( ((sample2 & ((type)1<<j)) == 0)                                        \
                         ==                                                                     \
                         (j >= type_bitsize-shc ? TRUE : ((sample1 & ((type)1<<(j+shc))) == 0)) \
                     ) )                                                                        \
                    { right_works = FALSE; }                                                    \
        }   }   }                                                                               \
      if (!left_works)                                                                          \
        { printf("#error \"Left shift of integers of type %s does not work!!\"\n",typestr); }   \
      if (!right_works)                                                                         \
        { printf("#error \"Right shift of integers of type %s does not work!!\"\n",typestr); }  \
    }
#define test_integer_sshift(type,typestr,type_bitsize)  \
  if (type_bitsize >= 0)                                                                       \
    { int i,j,shc;                                                                             \
      type sample1,sample2;                                                                    \
      boolean left_works = TRUE, right_works = TRUE;                                           \
      for (i = 0; i<100; i++)                                                                  \
        { get_a_random(type,type_bitsize,sample1);                                             \
          for (shc = 0; shc < type_bitsize; shc++)                                             \
            { sample2 = sample1 << shc;                                                        \
              for (j=0; j < type_bitsize; j++)                                                 \
                { if (!( ((sample2 & ((type)1<<j)) == 0)                                       \
                         ==                                                                    \
                         (j < shc ? TRUE : ((sample1 & ((type)1<<(j-shc))) == 0))              \
                     ) )                                                                       \
                    { left_works = FALSE; }                                                    \
        }   }   }                                                                              \
      for (i = 0; i<100; i++)                                                                  \
        { get_a_random(type,type_bitsize,sample1);                                             \
          for (shc = 0; shc < type_bitsize; shc++)                                             \
            { sample2 = sample1 >> shc;                                                        \
              for (j=0; j < type_bitsize; j++)                                                 \
                { if (!( ((sample2 & ((type)1<<j)) == 0)                                       \
                         ==                                                                    \
                         ((sample1 & ((type)1<< (j+shc>=type_bitsize ? type_bitsize-1 : j+shc))) == 0) \
                     ) )                                                                       \
                    { right_works = FALSE; }                                                   \
        }   }   }                                                                              \
      if (!left_works)                                                                         \
        { printf("#error \"Left shift of integers of type %s does not work!!\"\n",typestr); }  \
      if (!right_works)                                                                        \
        { printf("#error \"Right shift of integers of type %s does not work!!\"\n",typestr); } \
    }
  test_integer_ushift(uchar,"unsigned char",uchar_bitsize);
  test_integer_ushift(ushort,"unsigned short",ushort_bitsize);
  test_integer_ushift(uint,"unsigned int",uint_bitsize);
  test_integer_ushift(ulong,"unsigned long",ulong_bitsize);
#ifdef HAVE_LONGLONG
  test_integer_ushift(ulonglong,"unsigned long long",ulonglong_bitsize);
#endif
  test_integer_sshift(schar,"char",char_bitsize);
  test_integer_sshift(short,"short",short_bitsize);
  test_integer_sshift(int,"int",int_bitsize);
  test_integer_sshift(long,"long",long_bitsize);
#ifdef HAVE_LONGLONG
  test_integer_sshift(longlong,"long long",longlong_bitsize);
#endif
}

void main5(void) {
#define test_integer_casts(type1,type2,typestr1,typestr2,type1_bitsize,type2_bitsize,want)  \
  if (type1_bitsize <= type2_bitsize)                                                                      \
    { int i,j;                                                                                             \
      boolean modifies = FALSE;                                                                            \
      boolean zero_extends = TRUE;                                                                         \
      boolean sign_extends = TRUE;                                                                         \
      for (i = 0; i<100; i++)                                                                              \
        { type1 sample1;                                                                                   \
          type2 sample2;                                                                                   \
          get_a_random(type1,type1_bitsize,sample1);                                                       \
          sample2 = (type2)sample1;                                                                        \
          if (!(sample1 == (type1)sample2)) { modifies = TRUE; }                                           \
          for (j = 0; j<type1_bitsize; j++)                                                                \
            if (!( ((sample1 & ((type1)1<<j)) == 0) == ((sample2 & ((type2)1<<j)) == 0) ))                 \
              { zero_extends = FALSE; sign_extends = FALSE; }                                              \
          for (j = type1_bitsize; j<type2_bitsize; j++)                                                    \
            if (!((sample2 & ((type2)1<<j)) == 0))                                                         \
              { zero_extends = FALSE; }                                                                    \
          for (j = type1_bitsize; j<type2_bitsize; j++)                                                    \
            if (!( ((sample1 & ((type1)1<<(type1_bitsize-1))) == 0) == ((sample2 & ((type2)1<<j)) == 0) )) \
              { sign_extends = FALSE; }                                                                    \
        }                                                                                                  \
      if (modifies)                                                                                        \
        printf("#error \"Casts: (%s)(%s)(x) == x does not hold for every %s x !!\"\n",typestr1,typestr2,typestr1); \
      if (zero_extends && sign_extends)                                                                    \
        { if (!(type1_bitsize == type2_bitsize))                                                           \
            printf("#error \"Casts from %s to %s works by identity!!\"\n",typestr1,typestr2);              \
        }                                                                                                  \
      if (zero_extends && !sign_extends)                                                                   \
        { if ((type1_bitsize == type2_bitsize) || !(typestr1[0] == 'u') || !(want==1))                     \
            printf("#error \"Casts from %s to %s works by zero-extend!!\"\n",typestr1,typestr2);           \
        }                                                                                                  \
      if (sign_extends && !zero_extends)                                                                   \
        { if ((type1_bitsize == type2_bitsize) || (typestr1[0] == 'u') || !(want==2))                      \
            printf("#error \"Casts from %s to %s works by sign-extend!!\"\n",typestr1,typestr2);           \
        }                                                                                                  \
      if (!sign_extends && !zero_extends)                                                                  \
        printf("#error \"Casts from %s to %s works in an unknown manner!!\"\n",typestr1,typestr2);         \
    }
  /* erst Casts zwischen Integers vermutlich gleicher Größe: */
  test_integer_casts(schar,uchar,"char","unsigned char",char_bitsize,uchar_bitsize,0);
  test_integer_casts(short,ushort,"short","unsigned short",short_bitsize,ushort_bitsize,0);
  test_integer_casts(int,uint,"int","unsigned int",int_bitsize,uint_bitsize,0);
  test_integer_casts(long,ulong,"long","unsigned long",long_bitsize,ulong_bitsize,0);
  test_integer_casts(uchar,schar,"unsigned char","char",uchar_bitsize,char_bitsize,0);
  test_integer_casts(ushort,short,"unsigned short","short",ushort_bitsize,short_bitsize,0);
  test_integer_casts(uint,int,"unsigned int","int",uint_bitsize,int_bitsize,0);
  test_integer_casts(ulong,long,"unsigned long","long",ulong_bitsize,long_bitsize,0);
#ifdef HAVE_LONGLONG
  test_integer_casts(longlong,ulonglong,"long long","unsigned long long",longlong_bitsize,ulonglong_bitsize,0);
  test_integer_casts(ulonglong,longlong,"unsigned long long","long long",ulonglong_bitsize,longlong_bitsize,0);
#endif
  /* dann Casts zwischen Integers unterschiedlicher Größe, aber gleichen Vorzeichens: */
  test_integer_casts(uchar,ushort,"unsigned char","unsigned short",uchar_bitsize,ushort_bitsize,1);
  test_integer_casts(uchar,uint,"unsigned char","unsigned int",uchar_bitsize,uint_bitsize,1);
  test_integer_casts(uchar,ulong,"unsigned char","unsigned long",uchar_bitsize,ulong_bitsize,1);
  test_integer_casts(ushort,uint,"unsigned short","unsigned int",ushort_bitsize,uint_bitsize,1);
  test_integer_casts(ushort,ulong,"unsigned short","unsigned long",ushort_bitsize,ulong_bitsize,1);
  test_integer_casts(uint,ulong,"unsigned int","unsigned long",uint_bitsize,ulong_bitsize,1);
#ifdef HAVE_LONGLONG
  test_integer_casts(uchar,ulonglong,"unsigned char","unsigned long long",uchar_bitsize,ulonglong_bitsize,1);
  test_integer_casts(ushort,ulonglong,"unsigned short","unsigned long long",ushort_bitsize,ulonglong_bitsize,1);
  test_integer_casts(uint,ulonglong,"unsigned int","unsigned long long",uint_bitsize,ulonglong_bitsize,1);
  test_integer_casts(ulong,ulonglong,"unsigned long","unsigned long long",ulong_bitsize,ulonglong_bitsize,1);
#endif
  test_integer_casts(schar,short,"char","short",char_bitsize,short_bitsize,2);
  test_integer_casts(schar,int,"char","int",char_bitsize,int_bitsize,2);
  test_integer_casts(schar,long,"char","long",char_bitsize,long_bitsize,2);
  test_integer_casts(short,int,"short","int",short_bitsize,int_bitsize,2);
  test_integer_casts(short,long,"short","long",short_bitsize,long_bitsize,2);
  test_integer_casts(int,long,"int","long",int_bitsize,long_bitsize,2);
#ifdef HAVE_LONGLONG
  test_integer_casts(schar,longlong,"char","long long",char_bitsize,longlong_bitsize,2);
  test_integer_casts(short,longlong,"short","long long",short_bitsize,longlong_bitsize,2);
  test_integer_casts(int,longlong,"int","long long",int_bitsize,longlong_bitsize,2);
  test_integer_casts(long,longlong,"long","long long",long_bitsize,longlong_bitsize,2);
#endif
  /* dann Casts zwischen Integers unterschiedlicher Größe und unterschiedlichen Vorzeichens: */
  test_integer_casts(uchar,short,"unsigned char","short",uchar_bitsize,short_bitsize,1);
  test_integer_casts(uchar,int,"unsigned char","int",uchar_bitsize,int_bitsize,1);
  test_integer_casts(uchar,long,"unsigned char","long",uchar_bitsize,long_bitsize,1);
  test_integer_casts(ushort,int,"unsigned short","int",ushort_bitsize,int_bitsize,1);
  test_integer_casts(ushort,long,"unsigned short","long",ushort_bitsize,long_bitsize,1);
  test_integer_casts(uint,long,"unsigned int","long",uint_bitsize,long_bitsize,1);
#ifdef HAVE_LONGLONG
  test_integer_casts(uchar,longlong,"unsigned char","long long",uchar_bitsize,longlong_bitsize,1);
  test_integer_casts(ushort,longlong,"unsigned short","long long",ushort_bitsize,longlong_bitsize,1);
  test_integer_casts(uint,longlong,"unsigned int","long long",uint_bitsize,longlong_bitsize,1);
  test_integer_casts(ulong,longlong,"unsigned long","long long",ulong_bitsize,longlong_bitsize,1);
#endif
  test_integer_casts(schar,ushort,"char","unsigned short",char_bitsize,ushort_bitsize,2);
  test_integer_casts(schar,uint,"char","unsigned int",char_bitsize,uint_bitsize,2);
  test_integer_casts(schar,ulong,"char","unsigned long",char_bitsize,ulong_bitsize,2);
  test_integer_casts(short,uint,"short","unsigned int",short_bitsize,uint_bitsize,2);
  test_integer_casts(short,ulong,"short","unsigned long",short_bitsize,ulong_bitsize,2);
  test_integer_casts(int,ulong,"int","unsigned long",int_bitsize,ulong_bitsize,2);
#ifdef HAVE_LONGLONG
  test_integer_casts(schar,ulonglong,"char","unsigned long long",char_bitsize,ulonglong_bitsize,2);
  test_integer_casts(short,ulonglong,"short","unsigned long long",short_bitsize,ulonglong_bitsize,2);
  test_integer_casts(int,ulonglong,"int","unsigned long long",int_bitsize,ulonglong_bitsize,2);
  test_integer_casts(long,ulonglong,"long","unsigned long long",long_bitsize,ulonglong_bitsize,2);
#endif
}

void main6(void) {
#define check_sizeof_pointer(type,typestr)  \
  { if (!(sizeof(type) <= sizeof(long)))                                 \
      printf("#error \"Type %s does not fit into a long!!\"\n",typestr); \
  }
  check_sizeof_pointer(char*,"char *");
  check_sizeof_pointer(long*,"long *");
  check_sizeof_pointer(function*,"function *");
  pointer_bitsize = char_bitsize * sizeof(char*);
  printf("/* Pointers of type %s have %ld bits. */\n","char *",(long)pointer_bitsize);
  printf("#define pointer_bitsize %ld\n",(long)pointer_bitsize);
  printf("\n");
}

void main7(void) {
#define test_pointer_casts(type1,type2,typestr1,typestr2)  \
  if (!(sizeof(type1) == sizeof(type2)))                                                               \
    { printf("#error \"Pointer types %s and %s have different sizes!!\"\n",typestr1,typestr2); }       \
    else                                                                                               \
    { int i;                                                                                           \
      ulong differences1 = 0, differences2 = 0;                                                        \
      for (i = 0; i<100; i++)                                                                          \
        { ulong sample;                                                                                \
          type1 sample1;                                                                               \
          type2 sample2;                                                                               \
          get_a_random(ulong,ulong_bitsize,sample);                                                    \
          sample1 = (type1)sample;                                                                     \
          sample2 = (type2)sample;                                                                     \
          differences1 |= ((ulong)sample1 ^ (ulong)(type1)(sample2));                                  \
          differences2 |= ((ulong)sample2 ^ (ulong)(type2)(sample1));                                  \
        }                                                                                              \
      if (differences1==0)                                                                             \
        printf("/* Casts from %s to %s is OK (does nothing). */\n",typestr2,typestr1);                 \
      else                                                                                             \
      if (differences1 == ~(ulong)0)                                                                   \
        printf("#error \"Casts from %s to %s work in an unknown way!!\"\n",typestr2,typestr1);         \
      else                                                                                             \
        printf("#error \"Casts from %s to %s modify part 0x%8lX of pointer!!\"\n",typestr2,typestr1,differences1); \
      if (differences2==0)                                                                             \
        printf("/* Casts from %s to %s is OK (does nothing). */\n",typestr1,typestr2);                 \
      else                                                                                             \
      if (differences2 == ~(ulong)0)                                                                   \
        printf("#error \"Casts from %s to %s work in an unknown way!!\"\n",typestr1,typestr2);         \
      else                                                                                             \
        printf("#error \"Casts from %s to %s modify part 0x%8lX of pointer!!\"\n",typestr1,typestr2,differences2); \
    }
  test_pointer_casts(char*,long*,"char *","long *");
  test_pointer_casts(char*,function*,"char *","function *");
  printf("\n");
}

void main8(void) {
/* The following macro works only in C, not in C++, because C++ restricts the
   use of NULL pointers and also because C++ forbids defining types within a
   cast. */
#define alignmentof(type)  \
  (int)(&((struct { char dummy1; type dummy2; } *)0)->dummy2)
#define get_alignment(type,typestr)  \
  { struct { char dummy1; type dummy2; } dummy;                                                           \
    long alignment = (char*)&dummy.dummy2 - (char*)&dummy;                                                \
    printf("/* Type %s has sizeof = %ld and alignment = %ld. */\n",typestr,(long)sizeof(type),alignment); \
    if (!(typestr[0] == 'u') && !(typestr[string_length(typestr)-1] == '*'))                              \
      { printf("#define sizeof_"); printf_underscored(typestr); printf(" %ld\n",(long)sizeof(type));      \
        printf("#define alignment_"); printf_underscored(typestr); printf(" %ld\n",alignment);            \
      }                                                                                                   \
    if (!((alignment & (alignment-1)) == 0))                                                              \
      printf("#error \"The alignment %ld of type %s is not a power of two!!\"\n",alignment,typestr);      \
    printf("\n");                                                                                         \
  }
  get_alignment(char,"char"); get_alignment(uchar,"unsigned char");
  get_alignment(short,"short"); get_alignment(ushort,"unsigned short");
  get_alignment(int,"int"); get_alignment(uint,"unsigned int");
  get_alignment(long,"long"); get_alignment(ulong,"unsigned long");
#ifdef HAVE_LONGLONG
  get_alignment(longlong,"long long"); get_alignment(ulonglong,"unsigned long long");
#endif
  get_alignment(float,"float");
  get_alignment(double,"double");
  get_alignment(char*,"char *");
  get_alignment(long*,"long *");
  get_alignment(function*,"function *");
}

void main9(void) {
#define get_endian(type,typestr,type_bitsize)  \
  { if (type_bitsize == uchar_bitsize * sizeof(type))                                            \
      { union { uchar einzeln[sizeof(type)]; type gesamt; } x;                              \
        int i,j;                                                                                 \
        boolean big_endian = TRUE;                                                               \
        boolean little_endian = TRUE;                                                            \
        for (i = 0; i<100; i++)                                                                  \
          { type sample;                                                                         \
            get_a_random(type,type_bitsize,sample);                                              \
            x.gesamt = sample;                                                                   \
            for (j = 0; j<sizeof(type); j++, sample >>= uchar_bitsize)                           \
              { if (!( (sample & (((type)1<<uchar_bitsize)-1)) == x.einzeln[j] ))                \
                  { little_endian = FALSE; }                                                     \
                if (!( (sample & (((type)1<<uchar_bitsize)-1)) == x.einzeln[sizeof(type)-1-j] )) \
                  { big_endian = FALSE; }                                                        \
          }   }                                                                                  \
        if (big_endian && little_endian)                                                         \
          { if (!(sizeof(type) == 1))                                                            \
              printf("#error \"Endianness of type %s in memory doesn't matter.\"\n",typestr); }  \
        if (big_endian && !little_endian)                                                        \
          { printf("/* Type %s is stored BIG-ENDIAN in memory (i.e. like mc68000 or sparc). */\n",typestr); \
            printf("#define "); printf_underscored(&typestr[9]); printf("_big_endian\n");        \
          }                                                                                      \
        if (little_endian && !big_endian)                                                        \
          { printf("/* Type %s is stored LITTLE-ENDIAN in memory (i.e. like Z80 or VAX). */\n",typestr); \
            printf("#define "); printf_underscored(&typestr[9]); printf("_little_endian\n");     \
          }                                                                                      \
        if (!big_endian && !little_endian)                                                       \
          { printf("#error \"Type %s is stored in memory in an obscure manner!!\"\n",typestr); } \
      }                                                                                          \
      else                                                                                       \
      { printf("#error \"Endianness makes no sense for type %s !!\"\n",typestr); }               \
  }
  get_endian(uchar,"unsigned char",uchar_bitsize);
  get_endian(ushort,"unsigned short",ushort_bitsize);
  get_endian(uint,"unsigned int",uint_bitsize);
  get_endian(ulong,"unsigned long",ulong_bitsize);
#ifdef HAVE_LONGLONG
  get_endian(ulonglong,"unsigned long long",ulonglong_bitsize);
#endif
  printf("\n");
}

long get_stack_direction(void)
  { char dummy;
    static char* dummyaddr = (char*)0;
    if (!(dummyaddr == (char*)0))
      { return (&dummy) - dummyaddr; }
    else
      { dummyaddr = &dummy;
        { long result = get_stack_direction();
          /* The next assignment avoids tail recursion elimination (IRIX 6.4 CC). */
          dummyaddr = (char*)0;
          return result;
      } }
  }

void main10(void)
  { long stack_direction = get_stack_direction();
    if (stack_direction > 0)
      { printf("/* Stack grows up, ca. %ld bytes per function call. */\n",(long)stack_direction);
        printf("#define stack_grows_up\n");
      }
    else if (stack_direction < 0)
      { printf("/* Stack grows down, ca. %ld bytes per function call. */\n",-(long)stack_direction);
        printf("#define stack_grows_down\n");
      }
    else
      printf("#error \"Unknown stack model -- incorrect C semantics!!\"\n");
  }

int main()
{ main1();
  main2();
  main3();
  main4();
  main5();
  main6();
  main7();
  main8();
  main9();
  main10();
  if (ferror(stdout) || fclose(stdout)) return 1;
  return 0;
}
