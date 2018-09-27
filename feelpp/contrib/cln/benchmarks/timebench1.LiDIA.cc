#include <LiDIA/bigfloat.h>
#include <LiDIA/timer.h>

int main (int argc, char * argv[])
{
	int repetitions = 1;
	if ((argc >= 3) && !strcmp(argv[1],"-r")) {
		repetitions = atoi(argv[2]);
		argc -= 2; argv += 2;
	}
	if (argc < 1)
		exit(1);

	bigfloat::precision(1000);

	cerr << "pi\n";
	{ bigfloat p;
	  { timer t; t.set_print_mode(0); t.start_timer();
	    p = Pi();
	    t.stop_timer(); cerr << t << endl;
	  }
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = Pi(); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "gamma\n";
	{ bigfloat p;
	  { timer t; t.set_print_mode(0); t.start_timer();
	    p = Euler();
	    t.stop_timer(); cerr << t << endl;
	  }
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = Euler(); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "e\n";
	{ bigfloat p = E();
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = E(); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "sqrt(3)\n";
	{ bigfloat p = sqrt((bigfloat)3);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = sqrt((bigfloat)3); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "exp(log(2))\n";
	{ bigfloat p = exp(log((bigfloat)2));
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = exp(log((bigfloat)2)); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "log(exp(2))\n";
	{ bigfloat p = log(exp((bigfloat)2));
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = log(exp((bigfloat)2)); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "sin(pi/3)\n";
	{ bigfloat p = sin(Pi()/3);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = sin(Pi()/3); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "cos(pi/3)\n";
	{ bigfloat p = cos(Pi()/3);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = cos(Pi()/3); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "arcsin(sqrt(3)/2)\n";
	{ bigfloat p = asin(sqrt((bigfloat)3)/2);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = asin(sqrt((bigfloat)3)/2); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "arccos(sqrt(3)/2)\n";
	{ bigfloat p = acos(sqrt((bigfloat)3)/2);
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = acos(sqrt((bigfloat)3)/2); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "sinh(log(2))\n";
	{ bigfloat p = sinh(log((bigfloat)2));
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = sinh(log((bigfloat)2)); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "cosh(log(2))\n";
	{ bigfloat p = cosh(log((bigfloat)2));
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = cosh(log((bigfloat)2)); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "arsinh(pi)\n";
	{ bigfloat p = asinh(Pi());
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = asinh(Pi()); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

	cerr << "arcosh(pi)\n";
	{ bigfloat p = acosh(Pi());
	  { timer t; t.set_print_mode(0); t.start_timer();
	    for (int rep = repetitions; rep > 0; rep--)
	      { bigfloat p = acosh(Pi()); }
	    t.stop_timer(); cerr << t << endl;
	  }
	  cout << p << endl << endl;
	}

}



