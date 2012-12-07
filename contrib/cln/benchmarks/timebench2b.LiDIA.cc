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
	cerr << "Number of repetitions (except for pi,euler,e): " << repetitions << "\n";

	bigfloat::precision(digits);
	bigfloat x1 = sqrt((bigfloat)2);
	bigfloat x2 = sqrt((bigfloat)3);
	bigfloat x3 = log((bigfloat)2);

	cerr << "multiplication\n";
	{ bigfloat r = x1*x2;
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = x1*x2; }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "sqrt\n";
	{ bigfloat r = sqrt(x3);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = sqrt(x3); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "pi\n";
	{ bigfloat r;
	  { timer t; t.set_print_mode(0); t.start_timer();
	    r = Pi();
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "eulerconst\n";
	{ bigfloat r;
	  { timer t; t.set_print_mode(0); t.start_timer();
	    r = Euler();
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "e\n";
	{ bigfloat r;
	  { timer t; t.set_print_mode(0); t.start_timer();
	    r = E();
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "exp\n";
	{ bigfloat r = exp(-x1);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = exp(-x1); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "log\n";
	{ bigfloat r = log(x2);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = log(x2); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "sin\n";
	{ bigfloat r = sin(5*x1);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = sin(5*x1); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "cos\n";
	{ bigfloat r = cos(5*x1);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = cos(5*x1); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "asin\n";
	{ bigfloat r = asin(x3);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = asin(x3); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "acos\n";
	{ bigfloat r = acos(x3);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = acos(x3); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "atan\n";
	{ bigfloat r = atan(x3);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = atan(x3); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "sinh\n";
	{ bigfloat r = sinh(x2);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = sinh(x2); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "cosh\n";
	{ bigfloat r = cosh(x2);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = cosh(x2); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "asinh\n";
	{ bigfloat r = asinh(x3);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = asinh(x3); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "acosh\n";
	{ bigfloat r = acosh(1+x3);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = acosh(1+x3); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

	cerr << "atanh\n";
	{ bigfloat r = atanh(x3);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { r = atanh(x3); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << r << endl << endl;
	}

}
