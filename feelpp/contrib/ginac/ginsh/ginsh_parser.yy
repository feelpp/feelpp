/** @file ginsh_parser.yy
 *
 *  Input grammar definition for ginsh.
 *  This file must be processed with yacc/bison. */

/*
 *  GiNaC Copyright (C) 1999-2011 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


/*
 *  Definitions
 */

%{
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_RUSAGE
#include <sys/resource.h>
#else
#include <ctime>
#endif

#ifdef HAVE_UNISTD_H
#include <sys/types.h>
#include <unistd.h>
#endif

#include <stdexcept>

#include "ginsh.h"

#define YYERROR_VERBOSE 1

#ifdef HAVE_LIBREADLINE
// Original readline settings
static int orig_completion_append_character;
static const char *orig_basic_word_break_characters;

#if (RL_VERSION_MAJOR >= 5)
#define GINAC_RL_COMPLETER_CAST(a) const_cast<char *>((a))
#else
#define GINAC_RL_COMPLETER_CAST(a) (a)
#endif
#endif // HAVE_LIBREADLINE

// Expression stack for %, %% and %%%
static void push(const ex &e);
static ex exstack[3];
// Assigned symbols
static exmap assigned_symbol_table;

// Start and end time for the time() function
#ifdef HAVE_RUSAGE
static struct rusage start_time, end_time;
#define START_TIMER getrusage(RUSAGE_SELF, &start_time);
#define STOP_TIMER getrusage(RUSAGE_SELF, &end_time);
#define PRINT_TIME_USED cout << \
   (end_time.ru_utime.tv_sec - start_time.ru_utime.tv_sec) + \
       (end_time.ru_stime.tv_sec - start_time.ru_stime.tv_sec) + \
       double(end_time.ru_utime.tv_usec - start_time.ru_utime.tv_usec) / 1e6 + \
       double(end_time.ru_stime.tv_usec - start_time.ru_stime.tv_usec) / 1e6 \
                       << 's' << endl;
#else
static std::clock_t start_time, end_time;
#define START_TIMER start_time = std::clock();
#define STOP_TIMER end_time = std::clock();
#define PRINT_TIME_USED \
  cout << double(end_time - start_time)/CLOCKS_PER_SEC << 's' << endl;
#endif

// Table of functions (a multimap, because one function may appear with different
// numbers of parameters)
typedef ex (*fcnp)(const exprseq &e);
typedef ex (*fcnp2)(const exprseq &e, int serial);

struct fcn_desc {
	fcn_desc() : p(NULL), num_params(0), is_ginac(false), serial(0) {}
	fcn_desc(fcnp func, int num) : p(func), num_params(num), is_ginac(false), serial(0) {}
	fcn_desc(fcnp2 func, int num, int ser) : p((fcnp)func), num_params(num), is_ginac(true), serial(ser) {}

	fcnp p;		// Pointer to function
	int num_params;	// Number of parameters (0 = arbitrary)
	bool is_ginac;	// Flag: function is GiNaC function
	int serial;	// GiNaC function serial number (if is_ginac == true)
};

typedef multimap<string, fcn_desc> fcn_tab;
static fcn_tab fcns;

static fcn_tab::const_iterator find_function(const ex &sym, int req_params);

// Table to map help topics to help strings
typedef multimap<string, string> help_tab;
static help_tab help;

static void insert_fcn_help(const char *name, const char *str);
static void print_help(const string &topic);
static void print_help_topics(void);
%}

/* Tokens (T_LITERAL means a literal value returned by the parser, but not
   of class numeric or symbol (e.g. a constant or the FAIL object)) */
%token T_NUMBER T_SYMBOL T_LITERAL T_DIGITS T_QUOTE T_QUOTE2 T_QUOTE3
%token T_EQUAL T_NOTEQ T_LESSEQ T_GREATEREQ

%token T_QUIT T_WARRANTY T_PRINT T_IPRINT T_PRINTLATEX T_PRINTCSRC T_TIME
%token T_XYZZY T_INVENTORY T_LOOK T_SCORE T_COMPLEX_SYMBOLS T_REAL_SYMBOLS

/* Operator precedence and associativity */
%right '='
%left T_EQUAL T_NOTEQ
%left '<' '>' T_LESSEQ T_GREATEREQ
%left '+' '-'
%left '*' '/'
%nonassoc NEG
%right '^'
%nonassoc '!'

%start input


/*
 *  Grammar rules
 */

%%
input	: /* empty */
	| input line
	;

line	: ';'
	| exp ';' {
		try {
			cout << $1 << endl;
			push($1);
		} catch (exception &e) {
			cerr << e.what() << endl;
			YYERROR;
		}
	}
	| exp ':' {
		try {
			push($1);
		} catch (exception &e) {
			std::cerr << e.what() << endl;
			YYERROR;
		}
	}
	| T_PRINT '(' exp ')' ';' {
		try {
			$3.print(print_tree(std::cout));
		} catch (exception &e) {
			std::cerr << e.what() << endl;
			YYERROR;
		}
	}
	| T_IPRINT '(' exp ')' ';' {
		try {
			ex e = $3;
			if (!e.info(info_flags::integer))
				throw (std::invalid_argument("argument to iprint() must be an integer"));
			long i = ex_to<numeric>(e).to_long();
			cout << i << endl;
			cout << "#o" << oct << i << endl;
			cout << "#x" << hex << i << dec << endl;
		} catch (exception &e) {
			cerr << e.what() << endl;
			YYERROR;
		}
	}
	| T_PRINTLATEX '(' exp ')' ';' {
		try {
			$3.print(print_latex(std::cout)); cout << endl;
		} catch (exception &e) {
			std::cerr << e.what() << endl;
			YYERROR;
		}
	}
	| T_PRINTCSRC '(' exp ')' ';' {
		try {
			$3.print(print_csrc_double(std::cout)); cout << endl;
		} catch (exception &e) {
			std::cerr << e.what() << endl;
			YYERROR;
		}
	}
	| '?' T_SYMBOL 		{print_help(ex_to<symbol>($2).get_name());}
	| '?' T_TIME		{print_help("time");}
	| '?' T_PRINT		{print_help("print");}
	| '?' T_IPRINT		{print_help("iprint");}
	| '?' T_PRINTLATEX	{print_help("print_latex");}
	| '?' T_PRINTCSRC	{print_help("print_csrc");}
	| '?' '?'		{print_help_topics();}
	| T_QUIT		{YYACCEPT;}
	| T_WARRANTY {
		cout << "This program is free software; you can redistribute it and/or modify it under\n";
		cout << "the terms of the GNU General Public License as published by the Free Software\n";
		cout << "Foundation; either version 2 of the License, or (at your option) any later\n";
		cout << "version.\n";
		cout << "This program is distributed in the hope that it will be useful, but WITHOUT\n";
		cout << "ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS\n";
		cout << "FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more\n";
		cout << "details.\n";
		cout << "You should have received a copy of the GNU General Public License along with\n";
		cout << "this program. If not, write to the Free Software Foundation, 675 Mass Ave,\n";
		cout << "Cambridge, MA 02139, USA.\n";
	}
	| T_XYZZY		{cout << "Nothing happens.\n";}
	| T_INVENTORY		{cout << "You're not carrying anything.\n";}
	| T_LOOK		{cout << "You're in a twisty little maze of passages, all alike.\n";}
	| T_SCORE {
		cout << "If you were to quit now, you would score ";
		cout << (syms.size() > 350 ? 350 : syms.size());
		cout << " out of a possible 350.\n";
	}
	| T_REAL_SYMBOLS { symboltype = domain::real; }
	| T_COMPLEX_SYMBOLS { symboltype = domain::complex; }
	| T_TIME { START_TIMER } '(' exp ')' { STOP_TIMER PRINT_TIME_USED }
	| error ';'		{yyclearin; yyerrok;}
	| error ':'		{yyclearin; yyerrok;}
	;

exp	: T_NUMBER		{$$ = $1;}
	| T_SYMBOL		{
		exmap::const_iterator i = assigned_symbol_table.find($1);
		if (i == assigned_symbol_table.end())
			$$ = $1;
		else
			$$ = i->second.eval();
	}
	| '\'' T_SYMBOL '\''	{$$ = $2;}
	| T_LITERAL		{$$ = $1;}
	| T_DIGITS		{$$ = $1;}
	| T_QUOTE		{$$ = exstack[0];}
	| T_QUOTE2		{$$ = exstack[1];}
	| T_QUOTE3		{$$ = exstack[2];}
	| T_SYMBOL '(' exprseq ')' {
		fcn_tab::const_iterator i = find_function($1, $3.nops());
		if (i->second.is_ginac) {
			$$ = ((fcnp2)(i->second.p))(ex_to<exprseq>($3), i->second.serial);
		} else {
			$$ = (i->second.p)(ex_to<exprseq>($3));
		}
	}
	| T_DIGITS '=' T_NUMBER	{$$ = $3; Digits = ex_to<numeric>($3).to_int();}
	| T_SYMBOL '=' exp	{$$ = $3; assigned_symbol_table[$1] = $3; }
	| exp T_EQUAL exp	{$$ = $1 == $3;}
	| exp T_NOTEQ exp	{$$ = $1 != $3;}
	| exp '<' exp		{$$ = $1 < $3;}
	| exp T_LESSEQ exp	{$$ = $1 <= $3;}
	| exp '>' exp		{$$ = $1 > $3;}
	| exp T_GREATEREQ exp	{$$ = $1 >= $3;}
	| exp '+' exp		{$$ = $1 + $3;}
	| exp '-' exp		{$$ = $1 - $3;}
	| exp '*' exp		{$$ = $1 * $3;}
	| exp '/' exp		{$$ = $1 / $3;}
	| '-' exp %prec NEG	{$$ = -$2;}
	| '+' exp %prec NEG	{$$ = $2;}
	| exp '^' exp		{$$ = power($1, $3);}
	| exp '!'		{$$ = factorial($1);}
	| '(' exp ')'		{$$ = $2;}
	| '{' list_or_empty '}'	{$$ = $2;}
	| '[' matrix ']'	{$$ = lst_to_matrix(ex_to<lst>($2));}
	;

exprseq	: exp			{$$ = exprseq($1);}
	| exprseq ',' exp	{exprseq es(ex_to<exprseq>($1)); $$ = es.append($3);}
	;

list_or_empty: /* empty */	{$$ = *new lst;}
	| list			{$$ = $1;}
	;

list	: exp			{$$ = lst($1);}
	| list ',' exp		{lst l(ex_to<lst>($1)); $$ = l.append($3);}
	;

matrix	: '[' row ']'		{$$ = lst($2);}
	| matrix ',' '[' row ']' {lst l(ex_to<lst>($1)); $$ = l.append($4);}
	;

row	: exp			{$$ = lst($1);}
	| row ',' exp		{lst l(ex_to<lst>($1)); $$ = l.append($3);}
	;


/*
 *  Routines
 */

%%
// Error print routine
int yyerror(const char *s)
{
	cerr << s << " at " << yytext << endl;
	return 0;
}

// Push expression "e" onto the expression stack (for ", "" and """)
static void push(const ex &e)
{
	exstack[2] = exstack[1];
	exstack[1] = exstack[0];
	exstack[0] = e;
}


/*
 *  Built-in functions
 */

static ex f_collect(const exprseq &e) {return e[0].collect(e[1]);}
static ex f_collect_distributed(const exprseq &e) {return e[0].collect(e[1], true);}
static ex f_collect_common_factors(const exprseq &e) {return collect_common_factors(e[0]);}
static ex f_convert_H_to_Li(const exprseq &e) {return convert_H_to_Li(e[0], e[1]);}
static ex f_degree(const exprseq &e) {return e[0].degree(e[1]);}
static ex f_denom(const exprseq &e) {return e[0].denom();}
static ex f_eval1(const exprseq &e) {return e[0].eval();}
static ex f_evalf1(const exprseq &e) {return e[0].evalf();}
static ex f_evalm(const exprseq &e) {return e[0].evalm();}
static ex f_eval_integ(const exprseq &e) {return e[0].eval_integ();}
static ex f_expand(const exprseq &e) {return e[0].expand();}
static ex f_factor(const exprseq &e) {return factor(e[0]);}
static ex f_gcd(const exprseq &e) {return gcd(e[0], e[1]);}
static ex f_has(const exprseq &e) {return e[0].has(e[1]) ? ex(1) : ex(0);}
static ex f_lcm(const exprseq &e) {return lcm(e[0], e[1]);}
static ex f_lcoeff(const exprseq &e) {return e[0].lcoeff(e[1]);}
static ex f_ldegree(const exprseq &e) {return e[0].ldegree(e[1]);}
static ex f_lsolve(const exprseq &e) {return lsolve(e[0], e[1]);}
static ex f_nops(const exprseq &e) {return e[0].nops();}
static ex f_normal1(const exprseq &e) {return e[0].normal();}
static ex f_numer(const exprseq &e) {return e[0].numer();}
static ex f_numer_denom(const exprseq &e) {return e[0].numer_denom();}
static ex f_pow(const exprseq &e) {return pow(e[0], e[1]);}
static ex f_sqrt(const exprseq &e) {return sqrt(e[0]);}
static ex f_sqrfree1(const exprseq &e) {return sqrfree(e[0]);}
static ex f_subs2(const exprseq &e) {return e[0].subs(e[1]);}
static ex f_tcoeff(const exprseq &e) {return e[0].tcoeff(e[1]);}

#define CHECK_ARG(num, type, fcn) if (!is_a<type>(e[num])) throw(std::invalid_argument("argument " #num " to " #fcn "() must be a " #type))

static ex f_charpoly(const exprseq &e)
{
	CHECK_ARG(0, matrix, charpoly);
	return ex_to<matrix>(e[0]).charpoly(e[1]);
}

static ex f_coeff(const exprseq &e)
{
	CHECK_ARG(2, numeric, coeff);
	return e[0].coeff(e[1], ex_to<numeric>(e[2]).to_int());
}

static ex f_content(const exprseq &e)
{
	return e[0].content(e[1]);
}

static ex f_decomp_rational(const exprseq &e)
{
	return decomp_rational(e[0], e[1]);
}

static ex f_determinant(const exprseq &e)
{
	CHECK_ARG(0, matrix, determinant);
	return ex_to<matrix>(e[0]).determinant();
}

static ex f_diag(const exprseq &e)
{
	size_t dim = e.nops();
	matrix &m = *new matrix(dim, dim);
	for (size_t i=0; i<dim; i++)
		m.set(i, i, e.op(i));
	return m;
}

static ex f_diff2(const exprseq &e)
{
	CHECK_ARG(1, symbol, diff);
	return e[0].diff(ex_to<symbol>(e[1]));
}

static ex f_diff3(const exprseq &e)
{
	CHECK_ARG(1, symbol, diff);
	CHECK_ARG(2, numeric, diff);
	return e[0].diff(ex_to<symbol>(e[1]), ex_to<numeric>(e[2]).to_int());
}

static ex f_divide(const exprseq &e)
{
	ex q;
	if (divide(e[0], e[1], q))
		return q;
	else
		return fail();
}

static ex f_eval2(const exprseq &e)
{
	CHECK_ARG(1, numeric, eval);
	return e[0].eval(ex_to<numeric>(e[1]).to_int());
}

static ex f_evalf2(const exprseq &e)
{
	CHECK_ARG(1, numeric, evalf);
	return e[0].evalf(ex_to<numeric>(e[1]).to_int());
}

static ex f_find(const exprseq &e)
{
	exset found;
	e[0].find(e[1], found);
	lst l;
	for (exset::const_iterator i = found.begin(); i != found.end(); ++i)
		l.append(*i);
	return l;
}

static ex f_fsolve(const exprseq &e)
{
	CHECK_ARG(1, symbol, fsolve);
	CHECK_ARG(2, numeric, fsolve);
	CHECK_ARG(3, numeric, fsolve);
	return fsolve(e[0], ex_to<symbol>(e[1]), ex_to<numeric>(e[2]), ex_to<numeric>(e[3]));
}

static ex f_integer_content(const exprseq &e)
{
	return e[0].expand().integer_content();
}

static ex f_integral(const exprseq &e)
{
	CHECK_ARG(0, symbol, integral);
	return integral(e[0], e[1], e[2], e[3]);
}

static ex f_inverse(const exprseq &e)
{
	CHECK_ARG(0, matrix, inverse);
	return ex_to<matrix>(e[0]).inverse();
}

static ex f_is(const exprseq &e)
{
	CHECK_ARG(0, relational, is);
	return (bool)ex_to<relational>(e[0]) ? ex(1) : ex(0);
}

class apply_map_function : public map_function {
	ex apply;
public:
	apply_map_function(const ex & a) : apply(a) {}
	virtual ~apply_map_function() {}
	ex operator()(const ex & e) { return apply.subs(wild() == e, true); }
};

static ex f_map(const exprseq &e)
{
	apply_map_function fcn(e[1]);
	return e[0].map(fcn);
}

static ex f_match(const exprseq &e)
{
	exmap repls;
	if (e[0].match(e[1], repls)) {
		lst repl_lst;
		for (exmap::const_iterator i = repls.begin(); i != repls.end(); ++i)
			repl_lst.append(relational(i->first, i->second, relational::equal));
		return repl_lst;
	}
	throw std::runtime_error("FAIL");
}

static ex f_normal2(const exprseq &e)
{
	CHECK_ARG(1, numeric, normal);
	return e[0].normal(ex_to<numeric>(e[1]).to_int());
}

static ex f_op(const exprseq &e)
{
	CHECK_ARG(1, numeric, op);
	int n = ex_to<numeric>(e[1]).to_int();
	if (n < 0 || n >= (int)e[0].nops())
		throw(std::out_of_range("second argument to op() is out of range"));
	return e[0].op(n);
}

static ex f_prem(const exprseq &e)
{
	return prem(e[0], e[1], e[2]);
}

static ex f_primpart(const exprseq &e)
{
	return e[0].primpart(e[1]);
}

static ex f_quo(const exprseq &e)
{
	return quo(e[0], e[1], e[2]);
}

static ex f_rank(const exprseq &e)
{
	CHECK_ARG(0, matrix, rank);
	return ex_to<matrix>(e[0]).rank();
}

static ex f_rem(const exprseq &e)
{
	return rem(e[0], e[1], e[2]);
}

static ex f_resultant(const exprseq &e)
{
	CHECK_ARG(2, symbol, resultant);
	return resultant(e[0], e[1], ex_to<symbol>(e[2]));
}

static ex f_series(const exprseq &e)
{
	CHECK_ARG(2, numeric, series);
	return e[0].series(e[1], ex_to<numeric>(e[2]).to_int());
}

static ex f_sprem(const exprseq &e)
{
	return sprem(e[0], e[1], e[2]);
}

static ex f_sqrfree2(const exprseq &e)
{
	CHECK_ARG(1, lst, sqrfree);
	return sqrfree(e[0], ex_to<lst>(e[1]));
}

static ex f_subs3(const exprseq &e)
{
	CHECK_ARG(1, lst, subs);
	CHECK_ARG(2, lst, subs);
	return e[0].subs(ex_to<lst>(e[1]), ex_to<lst>(e[2]));
}

static ex f_trace(const exprseq &e)
{
	CHECK_ARG(0, matrix, trace);
	return ex_to<matrix>(e[0]).trace();
}

static ex f_transpose(const exprseq &e)
{
	CHECK_ARG(0, matrix, transpose);
	return ex_to<matrix>(e[0]).transpose();
}

static ex f_unassign(const exprseq &e)
{
	CHECK_ARG(0, symbol, unassign);
	exmap::iterator i = assigned_symbol_table.find(e[0]);
	if (i != assigned_symbol_table.end())
		assigned_symbol_table.erase(i);
	return e[0];
}

static ex f_unit(const exprseq &e)
{
	return e[0].unit(e[1]);
}

static ex f_dummy(const exprseq &e)
{
	throw(std::logic_error("dummy function called (shouldn't happen)"));
}

// Tables for initializing the "fcns" map and the function help topics
struct fcn_init {
	const char *name;
	fcnp p;
	int num_params;
};

static const fcn_init builtin_fcns[] = {
	{"charpoly", f_charpoly, 2},
	{"coeff", f_coeff, 3},
	{"collect", f_collect, 2},
	{"collect_common_factors", f_collect_common_factors, 1},
	{"collect_distributed", f_collect_distributed, 2},
	{"content", f_content, 2},
	{"convert_H_to_Li", f_convert_H_to_Li, 2},
	{"decomp_rational", f_decomp_rational, 2},
	{"degree", f_degree, 2},
	{"denom", f_denom, 1},
	{"determinant", f_determinant, 1},
	{"diag", f_diag, 0},
	{"diff", f_diff2, 2},
	{"diff", f_diff3, 3},
	{"divide", f_divide, 2},
	{"eval", f_eval1, 1},
	{"eval", f_eval2, 2},
	{"evalf", f_evalf1, 1},
	{"evalf", f_evalf2, 2},
	{"evalm", f_evalm, 1},
	{"eval_integ", f_eval_integ, 1},
	{"expand", f_expand, 1},
	{"factor", f_factor, 1},
	{"find", f_find, 2},
	{"fsolve", f_fsolve, 4},
	{"gcd", f_gcd, 2},
	{"has", f_has, 2},
	{"integer_content", f_integer_content, 1},
	{"integral", f_integral, 4},
	{"inverse", f_inverse, 1},
	{"iprint", f_dummy, 0},      // for Tab-completion
	{"is", f_is, 1},
	{"lcm", f_lcm, 2},
	{"lcoeff", f_lcoeff, 2},
	{"ldegree", f_ldegree, 2},
	{"lsolve", f_lsolve, 2},
	{"map", f_map, 2},
	{"match", f_match, 2},
	{"nops", f_nops, 1},
	{"normal", f_normal1, 1},
	{"normal", f_normal2, 2},
	{"numer", f_numer, 1},
	{"numer_denom", f_numer_denom, 1},
	{"op", f_op, 2},
	{"pow", f_pow, 2},
	{"prem", f_prem, 3},
	{"primpart", f_primpart, 2},
	{"print", f_dummy, 0},       // for Tab-completion
	{"print_csrc", f_dummy, 0},  // for Tab-completion
	{"print_latex", f_dummy, 0}, // for Tab-completion
	{"quo", f_quo, 3},
	{"rank", f_rank, 1},
	{"rem", f_rem, 3},
	{"resultant", f_resultant, 3},
	{"series", f_series, 3},
	{"sprem", f_sprem, 3},
	{"sqrfree", f_sqrfree1, 1},
	{"sqrfree", f_sqrfree2, 2},
	{"sqrt", f_sqrt, 1},
	{"subs", f_subs2, 2},
	{"subs", f_subs3, 3},
	{"tcoeff", f_tcoeff, 2},
	{"time", f_dummy, 0},        // for Tab-completion
	{"trace", f_trace, 1},
	{"transpose", f_transpose, 1},
	{"unassign", f_unassign, 1},
	{"unit", f_unit, 2},
	{NULL, f_dummy, 0}           // End marker
};

struct fcn_help_init {
	const char *name;
	const char *help;
};

static const fcn_help_init builtin_help[] = {
	{"acos", "inverse cosine function"},
	{"acosh", "inverse hyperbolic cosine function"},
	{"asin", "inverse sine function"},
	{"asinh", "inverse hyperbolic sine function"},
	{"atan", "inverse tangent function"},
	{"atan2", "inverse tangent function with two arguments"},
	{"atanh", "inverse hyperbolic tangent function"},
	{"beta", "Beta function"},
	{"binomial", "binomial function"},
	{"cos", "cosine function"},
	{"cosh", "hyperbolic cosine function"},
	{"exp", "exponential function"},
	{"factorial", "factorial function"},
	{"lgamma", "natural logarithm of Gamma function"},
	{"tgamma", "Gamma function"},
	{"log", "natural logarithm"},
	{"psi", "psi function\npsi(x) is the digamma function, psi(n,x) the nth polygamma function"},
	{"sin", "sine function"},
	{"sinh", "hyperbolic sine function"},
	{"tan", "tangent function"},
	{"tanh", "hyperbolic tangent function"},
	{"zeta", "zeta function\nzeta(x) is Riemann's zeta function, zetaderiv(n,x) its nth derivative.\nIf x is a GiNaC::lst, it is a multiple zeta value\nzeta(x,s) is an alternating Euler sum"},
	{"Li2", "dilogarithm"},
	{"Li3", "trilogarithm"},
	{"Li", "(multiple) polylogarithm"},
	{"S", "Nielsen's generalized polylogarithm"},
	{"H", "harmonic polylogarithm"},
	{"Order", "order term function (for truncated power series)"},
	{"Derivative", "inert differential operator"},
	{NULL, NULL}	// End marker
};

#include "ginsh_extensions.h"


/*
 *  Add functions to ginsh
 */

// Functions from fcn_init array
static void insert_fcns(const fcn_init *p)
{
	while (p->name) {
		fcns.insert(make_pair(string(p->name), fcn_desc(p->p, p->num_params)));
		p++;
	}
}

static ex f_ginac_function(const exprseq &es, int serial)
{
	return GiNaC::function(serial, es).eval(1);
}

// All registered GiNaC functions
namespace GiNaC {
static void ginsh_get_ginac_functions(void)
{
	vector<function_options> gfv = function::get_registered_functions();
	vector<function_options>::const_iterator i = gfv.begin(), end = gfv.end();
	unsigned serial = 0;
	while (i != end) {
		fcns.insert(make_pair(i->get_name(), fcn_desc(f_ginac_function, i->get_nparams(), serial)));
		++i;
		serial++;
	}
}
}


/*
 *  Find a function given a name and number of parameters. Throw exceptions on error.
 */

static fcn_tab::const_iterator find_function(const ex &sym, int req_params)
{
	const string &name = ex_to<symbol>(sym).get_name();
	typedef fcn_tab::const_iterator I;
	pair<I, I> b = fcns.equal_range(name);
	if (b.first == b.second)
		throw(std::logic_error("unknown function '" + name + "'"));
	else {
		for (I i=b.first; i!=b.second; i++)
			if ((i->second.num_params == 0) || (i->second.num_params == req_params))
				return i;
	}
	throw(std::logic_error("invalid number of arguments to " + name + "()"));
}


/*
 *  Insert help strings
 */

// Normal help string
static void insert_help(const char *topic, const char *str)
{
	help.insert(make_pair(string(topic), string(str)));
}

// Help string for functions, automatically generates synopsis
static void insert_fcn_help(const char *name, const char *str)
{
	typedef fcn_tab::const_iterator I;
	pair<I, I> b = fcns.equal_range(name);
	if (b.first != b.second) {
		string help_str = string(name) + "(";
		for (int i=0; i<b.first->second.num_params; i++) {
			if (i)
				help_str += ", ";
			help_str += "expression";
		}
		help_str += ") - ";
		help_str += str;
		help.insert(make_pair(string(name), help_str));
	}
}

// Help strings for functions from fcn_help_init array
static void insert_help(const fcn_help_init *p)
{
	while (p->name) {
		insert_fcn_help(p->name, p->help);
		p++;
	}
}


/*
 *  Print help to cout
 */

// Help for a given topic
static void print_help(const string &topic)
{
	typedef help_tab::const_iterator I;
	pair<I, I> b = help.equal_range(topic);
	if (b.first == b.second)
		cout << "no help for '" << topic << "'\n";
	else {
		for (I i=b.first; i!=b.second; i++)
			cout << i->second << endl;
	}
}

// List of help topics
static void print_help_topics(void)
{
	cout << "Available help topics:\n";
	help_tab::const_iterator i;
	string last_name = string("*");
	int num = 0;
	for (i=help.begin(); i!=help.end(); i++) {
		// Don't print duplicates
		if (i->first != last_name) {
			if (num)
				cout << ", ";
			num++;
			cout << i->first;
			last_name = i->first;
		}
	}
	cout << "\nTo get help for a certain topic, type ?topic\n";
}


/*
 *  Function name completion functions for readline
 */

static char *fcn_generator(const char *text, int state)
{
	static int len;				// Length of word to complete
	static fcn_tab::const_iterator index;	// Iterator to function being currently considered

	// If this is a new word to complete, initialize now
	if (state == 0) {
		index = fcns.begin();
		len = strlen(text);
	}

	// Return the next function which partially matches
	while (index != fcns.end()) {
		const char *fcn_name = index->first.c_str();
		++index;
		if (strncmp(fcn_name, text, len) == 0)
			return strdup(fcn_name);
	}
	return NULL;
}

#ifdef HAVE_LIBREADLINE
static char **fcn_completion(const char *text, int start, int end)
{
	if (rl_line_buffer[0] == '!') {
		// For shell commands, revert back to filename completion
		rl_completion_append_character = orig_completion_append_character;
		rl_basic_word_break_characters = const_cast<char*>(orig_basic_word_break_characters);
		rl_completer_word_break_characters = GINAC_RL_COMPLETER_CAST(rl_basic_word_break_characters);
		return rl_completion_matches(text, rl_filename_completion_function);
	} else {
		// Otherwise, complete function names
		rl_completion_append_character = '(';
		rl_basic_word_break_characters = " \t\n\"#$%&'()*+,-./:;<=>?@[\\]^`{|}~";
		rl_completer_word_break_characters = GINAC_RL_COMPLETER_CAST(rl_basic_word_break_characters);
		return rl_completion_matches(text, fcn_generator);
	}
}
#endif // HAVE_LIBREADLINE

static void ginsh_readline_init(char* name)
{
#ifdef HAVE_LIBREADLINE
	// Init readline completer
	rl_readline_name = name;
	rl_attempted_completion_function = fcn_completion;
	orig_completion_append_character = rl_completion_append_character;
	orig_basic_word_break_characters = rl_basic_word_break_characters;
#endif // HAVE_LIBREADLINE
}

void greeting(void)
{
    cout << "ginsh - GiNaC Interactive Shell (GiNaC V" << GINACLIB_VERSION << ")" << endl;
    cout << "  __,  _______  Copyright (C) 1999-2011 Johannes Gutenberg University Mainz,\n"
         << " (__) *       | Germany.  This is free software with ABSOLUTELY NO WARRANTY.\n"
         << "  ._) i N a C | You are welcome to redistribute it under certain conditions.\n"
         << "<-------------' For details type `warranty;'.\n" << endl;
    cout << "Type ?? for a list of help topics." << endl;
}

/*
 *  Main program
 */

int main(int argc, char **argv)
{
	// Print banner in interactive mode
	if (isatty(0)) 
		greeting();
	assigned_symbol_table = exmap();

	// Init function table
	insert_fcns(builtin_fcns);
	insert_fcns(extended_fcns);
	ginsh_get_ginac_functions();

	// Init help for operators (automatically generated from man page)
	insert_help("operators", "Operators in falling order of precedence:");
#include "ginsh_op_help.h"

	// Init help for built-in functions (automatically generated from man page)
#include "ginsh_fcn_help.h"

	// Help for GiNaC functions is added manually
	insert_help(builtin_help);
	insert_help(extended_help);

	// Help for other keywords
	insert_help("print", "print(expression) - dumps the internal structure of the given expression (for debugging)");
	insert_help("iprint", "iprint(expression) - prints the given integer expression in decimal, octal, and hexadecimal bases");
	insert_help("print_latex", "print_latex(expression) - prints a LaTeX representation of the given expression");
	insert_help("print_csrc", "print_csrc(expression) - prints a C source code representation of the given expression");

	ginsh_readline_init(argv[0]);

	// Init input file list, open first file
	num_files = argc - 1;
	file_list = argv + 1;
	if (num_files) {
		yyin = fopen(*file_list, "r");
		if (yyin == NULL) {
			cerr << "Can't open " << *file_list << endl;
			exit(1);
		}
		num_files--;
		file_list++;
	}

	// Parse input, catch all remaining exceptions
	int result;
again:	try {
		result = yyparse();
	} catch (exception &e) {
		cerr << e.what() << endl;
		goto again;
	}
	return result;
}
