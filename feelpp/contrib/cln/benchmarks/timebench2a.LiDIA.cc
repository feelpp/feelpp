#include <LiDIA/bigint.h>
#include <LiDIA/bigfloat.h>
#include <LiDIA/timer.h>
#include <cstdlib>
#include <cstring>

int main (int argc, char * argv[])
{
	int digits = 100;
	int repetitions = 1;
	while (argc >= 3) {
		if (!strcmp(argv[1],"-r")) {
			repetitions = atoi(argv[2]);
			argc -= 2; argv += 2;
			continue;
		}
		if (!strcmp(argv[1],"-n")) {
			digits = atoi(argv[2]);
			argc -= 2; argv += 2;
			continue;
		}
		break;
	}
	if (argc < 1)
		exit(1);

	cerr << "Number of digits: " << digits << "\n";
	cerr << "Number of repetitions: " << repetitions << "\n";

	bigint pow; power(pow, (bigint)10,digits);
	bigfloat::precision(digits*2);
	bigint x1; truncate(x1, ((sqrt((bigfloat)5)+1)/2) * (pow*pow));
	bigfloat::precision(digits);
	bigint x2; truncate(x2, sqrt((bigfloat)3) * pow);
	bigint x3 = pow+1;

	cerr << "multiplication\n";
	{ bigint r = x1*x2;
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigint r = x1*x2; }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "division\n";
	{ bigint q; bigint r; div_rem(q,r, x1,x2);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigint q; bigint r; div_rem(q,r, x1,x2); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << q << endl << r << endl << endl;
	}

	cerr << "isqrt\n";
	{ bigint r; sqrt(r, x3);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigint r; sqrt(r, x3); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "gcd\n";
	{ bigint r = gcd(x1,x2);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigint r = gcd(x1,x2); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

}
