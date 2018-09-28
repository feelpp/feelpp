/*
 * The following program was used to compute the first 100,000,000 decimal
 * digits of e = exp(1), on December 18-20, 1998.
 * Timings on a Sun UltraSparc-II (296 MHz), running Solaris 2.6, equipped
 * with 512 MB RAM and 2 GB swap:
 *
 * 100 digits:
 *   computation of e:      real time:     0.002 s, run time:     0.000 s
 *   conversion to decimal: real time:     0.003 s, run time:     0.000 s
 * 1000 digits:
 *   computation of e:      real time:     0.018 s, run time:     0.020 s
 *   conversion to decimal: real time:     0.028 s, run time:     0.020 s
 * 10000 digits:
 *   computation of e:      real time:     0.488 s, run time:     0.480 s
 *   conversion to decimal: real time:     1.059 s, run time:     1.060 s
 * 100000 digits:
 *   computation of e:      real time:     8.139 s, run time:     8.010 s
 *   conversion to decimal: real time:    16.593 s, run time:    16.540 s
 * 1000000 digits:
 *   computation of e:      real time:   122.383 s, run time:   121.020 s
 *   conversion to decimal: real time:   252.524 s, run time:   250.760 s
 * 10000000 digits:
 *   computation of e:      real time:  2152.061 s, run time:  2056.430 s
 *   conversion to decimal: real time:  3579.670 s, run time:  3388.990 s
 * 100000000 digits:
 *   computation of e:      real time: 40061.367 s, run time: 30449.630 s
 *   conversion to decimal: real time: 54507.003 s, run time: 40063.510 s
 */

#include <cln/number.h>
#include <cln/io.h>
#include <cln/integer.h>
#include <cln/integer_io.h>
#include <cln/float.h>
#include <cln/float_io.h>
#include <cln/real.h>
#include <cln/complex.h>
#include <cstring>
#include <cln/timing.h>
#include <cmath>

using namespace std;
using namespace cln;

void
sum_exp1 (uintC a, uintC b, cl_I & first, cl_I & second)
{
  switch (b - a)
    {
    case 1:
      first = second = b;
      break;
    case 2:
        {
          cl_I s = (a + b) >> 1;
          second = s * b;
          first = second + b;
        }
      break;
    default:
      {
        cl_I lp, lq, rp, rq, tmp;
        uintC mid = (a + b) >> 1;
        sum_exp1 (a, mid, lp, lq);
        sum_exp1 (mid, b, rp, rq);
        tmp = lp * rq;
        first = tmp + rp; 
        second = lq * rq; 
      }
      break;
    }
}

namespace cln {
  extern const cl_LF cl_I_to_LF(const cl_I&, uintC);
}

void
const_exp1 (cl_LF & result, uintC dec)
{
  uintC c = (uintC) (dec * ::log (10.0));
  uintC n = dec;
  uintC actuallen = (uintC)(3.321928094 * dec / intDsize);
  n = (uintC) ((n + c) / ::log ((double)n));
  n = (uintC) ((n + c) / ::log ((double)n));
  n = (uintC) ((n + c) / ::log ((double)n));

  n += 2;
  actuallen += 2;

  cout << "n = " << n << endl;
  cout << "actuallen = " << actuallen << endl;
  cl_I p, q;
  sum_exp1 (0, n, p, q);
  cout << "sum_exp1 ends ok" << endl;
  result = The(cl_LF)(cl_I_to_LF (p, actuallen) / cl_I_to_LF (q, actuallen));
  cout << "const_exp1 returns ok" << endl;
}

int
main (int argc, char *argv[])
{
  long digits = 100;
  while (argc >= 3) {
        if (!strcmp(argv[1],"-n")) {
            digits = atol(argv[2]);
            argc -= 2; argv += 2;
            continue;
        }
        break;
    }
    if (argc < 1)
        return(1);

  cl_LF c1;
  long l = digits;
  cout << "\nCalculating exp1 to " << l << " decimals" << endl;
  { CL_TIMING;
    const_exp1 (c1, l);
  }
  { CL_TIMING;
    cout << "@" << endl;
    cout << c1 << endl;
    cout << "@" << endl;
  }
  return(0);
}
