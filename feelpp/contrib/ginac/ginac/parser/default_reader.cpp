/** @file default_reader.cpp
 *
 *  Implementation of the default and builtin readers (part of GiNaC's parser).
 **/

/*
 *  GiNaC Copyright (C) 1999-2016 Johannes Gutenberg University Mainz, Germany
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
#include <cmath>
#include <random>
#include <algorithm>
#include "parse_context.h"
#include "power.h"
#include "lst.h"
#include "operators.h"
#include "inifcns.h"
#include "function.h"
#include "constant.h"
#include <cstdint> // for uintptr_t

namespace GiNaC
{
DECLARE_FUNCTION_2P(rand);
static ex rand_eval(const ex & a,const ex & b) 
{
	return rand(a,b).hold();
}

static void rand_print_latex( const ex& arg1, const ex& arg2, const print_context& c )
{
    c.s << "{"; arg1.print(c); c.s << "|}";
}

static void rand_print_csrc_float( const ex& arg1, const ex& arg2, const print_context& c )
{
    c.s << "[](){ auto [dis, gen] = uniformDistribution("; 
	arg1.print(c); 
	c.s << ",";
	arg2.print( c );
	c.s << "); return dis(gen); }()";
}
REGISTER_FUNCTION(rand, eval_func( rand_eval).
                       evalf_func( rand_eval).
                       print_func<print_latex>( rand_print_latex).
                       print_func<print_csrc_float>( rand_print_csrc_float).
                       print_func<print_csrc_double>( rand_print_csrc_float));

DECLARE_FUNCTION_2P( uniform );
static ex uniform_eval( const ex& a, const ex& b )
{
    return uniform( a,b ).hold();
}

static void uniform_print_latex( const ex& arg1, const ex& arg2, const print_context& c )
{
    c.s << "{";
    arg1.print( c );
    c.s << "|}";
}

static void uniform_print_csrc_float( const ex& arg1, const ex& arg2, const print_context& c )
{
    c.s << "[](){ auto [dis, gen] = uniformDistribution(";
    arg1.print( c );
    c.s << ",";
    arg2.print( c );
    c.s << "); return dis(gen); }()";
}
REGISTER_FUNCTION( uniform, eval_func( uniform_eval ).
                         evalf_func( uniform_eval ).
						 print_func<print_latex>( uniform_print_latex ).
						 print_func<print_csrc_float>( uniform_print_csrc_float ).
						 print_func<print_csrc_double>( uniform_print_csrc_float ) );

DECLARE_FUNCTION_2P( loguniform );
static ex loguniform_eval( const ex& a, const ex& b )
{
    return loguniform( a,b ).hold();
}

static void loguniform_print_latex( const ex& arg1, const ex& arg2, const print_context& c )
{
    c.s << "{";
    arg1.print( c );
    c.s << "|}";
}

static void loguniform_print_csrc_float( const ex& arg1, const ex& arg2, const print_context& c )
{
    c.s << "[](){ auto [dis, gen] = uniformDistribution( std::log(";
    arg1.print( c );
    c.s << "),std::log(";
    arg2.print( c );
    c.s << ")); return exp(dis(gen)); }()";
}
REGISTER_FUNCTION( loguniform, eval_func( loguniform_eval ).
                         evalf_func( loguniform_eval ).
						 print_func<print_latex>( loguniform_print_latex ).
						 print_func<print_csrc_float>( loguniform_print_csrc_float ).
						 print_func<print_csrc_double>( loguniform_print_csrc_float ) );

DECLARE_FUNCTION_2P( normal );
static ex normal_eval( const ex& a, const ex& b )
{
    return normal( a,b ).hold();
}

static void normal_print_latex( const ex& arg1, const ex& arg2, const print_context& c )
{
    c.s << "{";
    arg1.print( c );
    c.s << "|}";
}

static void normal_print_csrc_float( const ex& arg1, const ex& arg2, const print_context& c )
{
    c.s << "[](){ auto [dis, gen] = normalDistribution(";
    arg1.print( c );
    c.s << ",";
    arg2.print( c );
    c.s << "); return dis(gen); }()";
}
REGISTER_FUNCTION( normal, eval_func( normal_eval ).
                         evalf_func( normal_eval ).
						 print_func<print_latex>( normal_print_latex ).
						 print_func<print_csrc_float>( normal_print_csrc_float ).
						 print_func<print_csrc_double>( normal_print_csrc_float ) );

DECLARE_FUNCTION_2P( lognormal );
static ex lognormal_eval( const ex& a, const ex& b )
{
    return lognormal( a,b ).hold();
}

static void lognormal_print_latex( const ex& arg1, const ex& arg2, const print_context& c )
{
    c.s << "{";
    arg1.print( c );
    c.s << "|}";
}

static void lognormal_print_csrc_float( const ex& arg1, const ex& arg2, const print_context& c )
{
    c.s << "[](){ auto [dis, gen] = lognormalDistribution(";
    arg1.print( c );
    c.s << ",";
    arg2.print( c );
    c.s << "); return dis(gen); }()";
}
REGISTER_FUNCTION( lognormal, eval_func( lognormal_eval ).
                         evalf_func( lognormal_eval ).
						 print_func<print_latex>( lognormal_print_latex ).
						 print_func<print_csrc_float>( lognormal_print_csrc_float ).
						 print_func<print_csrc_double>( lognormal_print_csrc_float ) );

DECLARE_FUNCTION_2P(mod);
static ex mod_eval(const ex & x, const ex& y ) 
{ 
	if (is_exactly_a<numeric>(y) && is_exactly_a<numeric>(x))
		return mod(ex_to<numeric>(x), ex_to<numeric>(y));
	
	return mod(x, y).hold();
}
static ex mod_evalf(const ex & x, const ex& y ) 
{ 
	if (is_exactly_a<numeric>(y) && is_exactly_a<numeric>(x))
		return mod(ex_to<numeric>(x), ex_to<numeric>(y));
	
	return mod(x, y).hold();	
}
                                                                                
static void mod_print_latex(const ex & arg1, const ex&  arg2,const print_context & c)
{
    c.s << "{"; arg1.print(c); c.s << "\%"; arg2.print(c); c.s << "|}";
}
                                                                                
static void mod_print_csrc_float(const ex & arg1, const ex & arg2, const print_context & c)
{
    c.s << "fmod("; arg1.print(c); c.s << ","; arg2.print(c); c.s << ")";
}

                                                                                
REGISTER_FUNCTION(mod, eval_func( mod_eval).
                       evalf_func( mod_evalf).
                       print_func<print_latex>( mod_print_latex).
                       print_func<print_csrc_float>( mod_print_csrc_float).
                       print_func<print_csrc_double>( mod_print_csrc_float));

DECLARE_FUNCTION_1P(floor);
static ex floor_eval(const ex & x) 
{ 
	if (is_exactly_a<numeric>(x))
		return floor(ex_to<numeric>(x));
	
	return floor(x).hold();
}
static ex floor_evalf(const ex & x) 
{ 
	if (is_exactly_a<numeric>(x))
		return floor(ex_to<numeric>(x));
	
	return floor(x).hold();	
}
                                                                                
static void floor_print_latex(const ex & arg1, const print_context & c)
{
    c.s << "{"; arg1.print(c); c.s << "|}";
}
                                                                                
static void floor_print_csrc_float(const ex & arg1, const print_context & c)
{
    c.s << "std::floor("; arg1.print(c); c.s << ")";
}
REGISTER_FUNCTION(floor, eval_func( floor_eval).
                       evalf_func( floor_evalf).
                       print_func<print_latex>( floor_print_latex).
                       print_func<print_csrc_float>( floor_print_csrc_float).
                       print_func<print_csrc_double>( floor_print_csrc_float));
DECLARE_FUNCTION_1P(fract);
static ex fract_eval(const ex & x) 
{ 
	if (is_exactly_a<numeric>(x))
		return ex_to<numeric>(x)-floor(ex_to<numeric>(x));
	
	return fract(x).hold();
}

                                                                                
static void fract_print_latex(const ex & arg1, const print_context & c)
{
    c.s << "{"; arg1.print(c); c.s << "|}";
}
                                                                                
static void fract_print_csrc_float(const ex & arg1, const print_context & c)
{
    arg1.print(c); c.s << " - std::floor(";arg1.print(c); c.s << ")";
}
REGISTER_FUNCTION(fract, eval_func( fract_eval).
                       evalf_func( fract_eval).
                       print_func<print_latex>( fract_print_latex).
                       print_func<print_csrc_float>( fract_print_csrc_float).
                       print_func<print_csrc_double>( fract_print_csrc_float));
DECLARE_FUNCTION_1P(ceil);
static ex ceil_eval(const ex & x) 
{ 
	if (is_exactly_a<numeric>(x))
		return ceil(ex_to<numeric>(x));
	
	return ceil(x).hold();
}
static ex ceil_evalf(const ex & x) 
{ 
	if (is_exactly_a<numeric>(x))
		return ceil(ex_to<numeric>(x));
	
	return ceil(x).hold();	
}
                                                                                
static void ceil_print_latex(const ex & arg1, const print_context & c)
{
    c.s << "{"; arg1.print(c); c.s << "|}";
}
                                                                                
static void ceil_print_csrc_float(const ex & arg1, const print_context & c)
{
    c.s << "std::ceil("; arg1.print(c); c.s << ")";
}
                                                                                
REGISTER_FUNCTION(ceil, eval_func( ceil_eval).
                       evalf_func( ceil_evalf).
                       print_func<print_latex>( ceil_print_latex).
                       print_func<print_csrc_float>( ceil_print_csrc_float).
                       print_func<print_csrc_double>( ceil_print_csrc_float));

DECLARE_FUNCTION_1P(sign);
static ex sign_eval(const ex & x) 
{ 
	if (is_exactly_a<numeric>(x))
		return (numeric(0.) < ex_to<numeric>(x)) - (ex_to<numeric>(x) < numeric(0.));
	
	return sign(x).hold();
}
static ex sign_evalf(const ex & x) 
{ 
	if (is_exactly_a<numeric>(x))
		return (numeric(0.) < ex_to<numeric>(x)) - (ex_to<numeric>(x) < numeric(0.));
	
	return sign(x).hold();	
}
                                                                                
static void sign_print_latex(const ex & arg1,const print_context & c)
{
    c.s << "{"; arg1.print(c); c.s << "|}";
}
                                                                                
static void sign_print_csrc_float(const ex & arg1, const print_context & c)
{
    c.s << "(double(0) < "; arg1.print(c); c.s << ")-(";arg1.print(c); c.s << " < double(0))"; 
}
                                                                                
REGISTER_FUNCTION(sign, eval_func( sign_eval).
                       evalf_func( sign_evalf).
                       print_func<print_latex>( sign_print_latex).
                       print_func<print_csrc_float>( sign_print_csrc_float).
                       print_func<print_csrc_double>( sign_print_csrc_float));

DECLARE_FUNCTION_3P(clamp)
static ex clamp_eval( const ex & x, const ex & lo, const ex& hi ) 
{ 
	if (is_exactly_a<numeric>(x) && is_exactly_a<numeric>(lo) && is_exactly_a<numeric>(hi))
	{
		return std::clamp((ex_to<numeric>(x) - ex_to<numeric>(lo)) / (ex_to<numeric>(hi) - ex_to<numeric>(lo)), numeric(0.0), numeric(1.0));
	}
	
	return clamp(x, lo, hi).hold();
}
                                                                                
static void clamp_print_latex(const ex & x, const ex&  lo, const ex& hi,const print_context & c)
{
    c.s << "{"; lo.print(c); c.s << "\%"; hi.print(c); c.s << "|}";
}

static void clamp_print_csrc_float(const ex & x, const ex&  lo, const ex& hi,const print_context & c)												                                
{
 	c.s << "std::clamp( "; x.print(c); c.s << ", "; lo.print(c); c.s << ","; hi.print(c); c.s << ")";
}

                                                                                
REGISTER_FUNCTION(clamp, eval_func(clamp_eval).
                       evalf_func(clamp_eval).
                       print_func<print_latex>(clamp_print_latex).
                       print_func<print_csrc_float>(clamp_print_csrc_float).
                       print_func<print_csrc_double>(clamp_print_csrc_float));

DECLARE_FUNCTION_2P(step1);
static ex step1_eval( const ex & x, const ex & edge ) 
{ 
	if (is_exactly_a<numeric>(x) && is_exactly_a<numeric>(edge) )
	{
		if ( ex_to<numeric>(x) < ex_to<numeric>(edge) )
			return 0.;
		return 1;
	}
	
	return step1(x, edge).hold();
}
                                                                                
static void step1_print_latex(const ex & arg1, const ex&  arg2,const print_context & c)
{
    c.s << "{ step1("; arg1.print(c); c.s << ","; arg2.print(c); c.s << ")}";
}
                                                                                
static void step1_print_csrc_float(const ex & arg1, const ex & arg2, const print_context & c)
{
    c.s << "("; arg1.print(c); c.s << " < "; arg2.print(c); c.s << ") ? 0. : 1.";
}

                                                                                
REGISTER_FUNCTION(step1, eval_func(step1_eval).
                       evalf_func(step1_eval).
                       print_func<print_latex>(step1_print_latex).
                       print_func<print_csrc_float>(step1_print_csrc_float).
                       print_func<print_csrc_double>(step1_print_csrc_float));

DECLARE_FUNCTION_3P(smoothstep)
static ex smoothstep_eval( const ex & x, const ex & lo, const ex& hi ) 
{ 
	if (is_exactly_a<numeric>(x) && is_exactly_a<numeric>(lo) && is_exactly_a<numeric>(hi))
	{
		auto t = std::clamp((ex_to<numeric>(x) - ex_to<numeric>(lo)) / (ex_to<numeric>(hi) - ex_to<numeric>(lo)), numeric(0.0), numeric(1.0));
		return t * t * (3.0 - 2.0 * t);;
	}
	
	return smoothstep(x, lo, hi).hold();
}
                                                                                
static void smoothstep_print_latex(const ex & x, const ex&  lo, const ex& hi,const print_context & c)
{
    c.s << "{"; lo.print(c); c.s << "\%"; hi.print(c); c.s << "|}";
}

static void smoothstep_print_csrc_float(const ex & x, const ex&  lo, const ex& hi,const print_context & c)												                                
{
 	c.s << "[]( const double& t ){ return t * t * (3.0 - 2.0 * t); }( std::clamp( "; 
	c.s << "("; x.print(c); c.s << "-"; lo.print(c); c.s << ")/"; 
	c.s << "("; hi.print(c); c.s << "-"; lo.print(c); c.s << ")"; 
	c.s << ", 0.0, 1.0) )";
}

                                                                                
REGISTER_FUNCTION(smoothstep, eval_func(smoothstep_eval).
                       evalf_func(smoothstep_eval).
                       print_func<print_latex>(smoothstep_print_latex).
                       print_func<print_csrc_float>(smoothstep_print_csrc_float).
                       print_func<print_csrc_double>(smoothstep_print_csrc_float));


DECLARE_FUNCTION_3P(rectangle);
static ex rectangle_eval( const ex & x, const ex & lo, const ex& hi ) 
{ 
	if (is_exactly_a<numeric>(x) && is_exactly_a<numeric>(lo) && is_exactly_a<numeric>(hi) )
	{
		if ( ex_to<numeric>(x) < ex_to<numeric>(lo) )
			return 0.;
		else if ( ex_to<numeric>(x) > ex_to<numeric>(hi) )
			return 0.;
		return 1;
	}
	
	return rectangle(x, lo, hi).hold();
}
                                                                                
static void rectangle_print_latex(const ex & arg1, const ex&  arg2,const ex&  arg3,const print_context & c)
{
    c.s << "{ rectangle("; arg1.print(c); c.s << ","; arg2.print(c); c.s << ")}";
}
                                                                                
static void rectangle_print_csrc_float(const ex & x, const ex & lo, const ex & hi, const print_context & c)
{
    c.s << "(("; x.print(c); c.s << " < "; lo.print(c); c.s << ") ? 0. :";
	c.s << "("; hi.print(c); c.s << " < "; x.print(c); c.s << ") ? 0. : 1.)";
}

                                                                                
REGISTER_FUNCTION(rectangle, eval_func(rectangle_eval).
                       evalf_func(rectangle_eval).
                       print_func<print_latex>(rectangle_print_latex).
                       print_func<print_csrc_float>(rectangle_print_csrc_float).
                       print_func<print_csrc_double>(rectangle_print_csrc_float));

DECLARE_FUNCTION_3P(triangle);
static ex triangle_eval( const ex & x, const ex & lo, const ex& hi ) 
{ 
	if (is_exactly_a<numeric>(x) && is_exactly_a<numeric>(lo) && is_exactly_a<numeric>(hi) )
	{
		auto t = std::clamp( ( ex_to<numeric>(x)-ex_to<numeric>(lo))/(ex_to<numeric>(hi)-ex_to<numeric>(lo)), numeric(0), numeric(0) );
		return 1-abs(ex_to<numeric>(t));
	}
	return triangle(x, lo, hi).hold();
}
                                                                                
static void triangle_print_latex(const ex & arg1, const ex&  arg2, const ex&  arg3,const print_context & c)
{
    c.s << "{ triangle("; arg1.print(c); c.s << ","; arg2.print(c); c.s << ")}";
}
                                                                                
static void triangle_print_csrc_float(const ex & x, const ex & lo, const ex & hi, const print_context & c)
{
	c.s << "1-std::abs(";
	c.s << "std::clamp( ";
	c.s << " -1 + 2*( "; x.print(c); c.s << "-"; lo.print(c); c.s << ")/"; 
	c.s << "("; hi.print(c); c.s << "-"; lo.print(c); c.s << ")"; 
	c.s << ", -1.0, 1.0)";
	c.s << ")";
	
}

                                                                                
REGISTER_FUNCTION(triangle, eval_func(triangle_eval).
                       evalf_func(triangle_eval).
                       print_func<print_latex>(triangle_print_latex).
                       print_func<print_csrc_float>(triangle_print_csrc_float).
                       print_func<print_csrc_double>(triangle_print_csrc_float));


DECLARE_FUNCTION_5P(mapabcd);
static ex mapabcd_eval( const ex & x, const ex & a, const ex& b, const ex & c, const ex& d ) 
{ 
	if (is_exactly_a<numeric>(x) && is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b)&& is_exactly_a<numeric>(c) && is_exactly_a<numeric>(d) )
	{
		return std::clamp( ex_to<numeric>(c) + ( ex_to<numeric>(d)-ex_to<numeric>(c))*( ex_to<numeric>(x)-ex_to<numeric>(a))/(ex_to<numeric>(b)-ex_to<numeric>(a)), 
				     	   ex_to<numeric>(a), ex_to<numeric>(b) );
	}
	return mapabcd( x, a, b, c, d ).hold();
}
                                                                                
static void mapabcd_print_latex(const ex & x, const ex & a, const ex& b, const ex & c, const ex& d, const print_context & co)
{
    co.s << "{ mapabcd("; x.print(co); co.s << ","; a.print(co); co.s << ","; b.print(co);co.s << ","; c.print(co);co.s << ","; d.print(co);co.s << ")}";
}
                                                                                
static void mapabcd_print_csrc_float(const ex & x, const ex & a, const ex & b, const ex & c, const ex & d, const print_context & co)
{
	co.s << "std::clamp( ";
	c.print(co); co.s << "+ "; 
	co.s << "("; d.print(co); co.s << "-("; c.print(co); co.s << "))"; 
	co.s << "*( "; x.print(co); co.s << "-("; a.print(co); co.s << "))/"; 
	co.s << "("; b.print(co); co.s << "-("; a.print(co); co.s << "))"; 
	co.s << ","; c.print(co); co.s << ","; d.print(co); 
	co.s << ")";
	
}

                                                                                
REGISTER_FUNCTION(mapabcd, eval_func(mapabcd_eval).
                       evalf_func(mapabcd_eval).
                       print_func<print_latex>(mapabcd_print_latex).
                       print_func<print_csrc_float>(mapabcd_print_csrc_float).
                       print_func<print_csrc_double>(mapabcd_print_csrc_float));


DECLARE_FUNCTION_4P(pulse);
static ex pulse_eval( const ex & x, const ex & a, const ex& b, const ex & p ) 
{ 
	if (is_exactly_a<numeric>(x) && is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b)&& is_exactly_a<numeric>(p) )
	{
		if ( mod( ex_to<numeric>(x), ex_to<numeric>(p) ) < ex_to<numeric>(a) )
			return 0.;
		else if ( mod( ex_to<numeric>(x), ex_to<numeric>(p)  ) > ex_to<numeric>(b)   )
			return 0.;
		return mod( ex_to<numeric>(x), ex_to<numeric>(p) );
	}
	return pulse( x, a, b, p ).hold();
}
                                                                                
static void pulse_print_latex(const ex & x, const ex & a, const ex& b, const ex & p, const print_context & co)
{
    co.s << "{ pulse("; x.print(co); co.s << ","; a.print(co); co.s << ","; b.print(co);co.s << ","; p.print(co); co.s << ")}";
}
                                                                                
static void pulse_print_csrc_float(const ex & x, const ex & a, const ex & b, const ex & p, const print_context & co)
{
	co.s << "(( std::fmod("; x.print(co); co.s << ","; p.print(co); co.s << ") < "; a.print(co); co.s << ") ? 0. :";
	co.s << "("; b.print(co); co.s << " < ";  co.s << "( std::fmod("; x.print(co); co.s << ","; p.print(co); co.s << ")"; co.s << ") ? 0. : 1.))";	
}

/**
 * @brief rectangle pulse of period @p
 * @return the ginac expression of the pulse
 */                                                                                
REGISTER_FUNCTION(pulse, eval_func(pulse_eval).
                       evalf_func(pulse_eval).
                       print_func<print_latex>(pulse_print_latex).
                       print_func<print_csrc_float>(pulse_print_csrc_float).
                       print_func<print_csrc_double>(pulse_print_csrc_float));

DECLARE_FUNCTION_3P(sinewave);
static ex sinewave_eval( const ex & x, const ex & f, const ex & phi ) 
{ 
	if (is_exactly_a<numeric>(x) && is_exactly_a<numeric>(f)&& is_exactly_a<numeric>(phi) )
	{
		return sin( 2*Pi*ex_to<numeric>(f)*ex_to<numeric>(x) + ex_to<numeric>(phi) );
	}
	return sinewave( x, f, phi ).hold();
}
                                                                                
static void sinewave_print_latex(const ex & x, const ex & f, const ex & phi,  const print_context & co)
{
    co.s << "{sin(2*pi*"; f.print(co); co.s << "*"; x.print(co); co.s << "+"; phi.print(co); co.s << ")}";
}
                                                                                
static void sinewave_print_csrc_float(const ex & x, const ex & f, const ex & phi, const print_context & co)
{
	co.s << "("; co.s << "sin(2*pi*"; f.print(co); co.s << "*"; x.print(co); co.s << "+"; phi.print(co); co.s << "))";
}

/**
 * @brief rectangle sinewave of period @p
 * @return the ginac expression of the sinewave
 */                                                                                
REGISTER_FUNCTION(sinewave, eval_func(sinewave_eval).
                       evalf_func(sinewave_eval).
                       print_func<print_latex>(sinewave_print_latex).
                       print_func<print_csrc_float>(sinewave_print_csrc_float).
                       print_func<print_csrc_double>(sinewave_print_csrc_float));
					   
static ex mod_reader(const exvector& ev)
{
	return GiNaC::mod(ev[0],ev[1]);
}

static ex sqrt_reader(const exvector& ev)
{
	return GiNaC::sqrt(ev[0]);
}

static ex pow_reader(const exvector& ev)
{
	return GiNaC::pow(ev[0], ev[1]);
}
static ex min_reader(const exvector& ev)
{
	return (ev[0]+ev[1])/2-GiNaC::abs(ev[0]-ev[1])/2;//GiNaC::min(ev[0], ev[1]);
}
static ex max_reader(const exvector& ev)
{
	return (ev[0]+ev[1])/2+GiNaC::abs(ev[0]-ev[1])/2; //GiNaC::max(ev[0], ev[1]);
}

static ex power_reader(const exvector& ev)
{
	return GiNaC::power(ev[0], ev[1]);
}

static ex lst_reader(const exvector& ev)
{
	return GiNaC::lst(ev.begin(), ev.end());
}


// function::registered_functions() is protected, but we need to access it
// TODO: add a proper const method to the `function' class, so we don't
// need this silly hack any more.
class registered_functions_hack : public function
{
public:
	static const std::vector<function_options>& get_registered_functions()
	{
		return function::registered_functions();
	}
private:
	registered_functions_hack();
	registered_functions_hack(const registered_functions_hack&);
	registered_functions_hack& operator=(const registered_functions_hack&);
};

// Encode an integer into a pointer to a function. Since functions
// are aligned (the minimal alignment depends on CPU architecture)
// we can distinguish between pointers and integers.
static reader_func encode_serial_as_reader_func(unsigned serial)
{
	uintptr_t u = (uintptr_t)serial;
	u = (u << 1) | (uintptr_t)1;
	reader_func ptr = (reader_func)((void *)u);
	return ptr;
}

const prototype_table& get_default_reader()
{
	using std::make_pair;
	static bool initialized = false;
	static prototype_table reader;
	if (!initialized) {
		
		reader[make_pair("sqrt", 1)] = sqrt_reader;
		reader[make_pair("pow", 2)] = pow_reader;
		reader[make_pair("mod", 2)] = mod_reader;
		reader[make_pair("power", 2)] = power_reader;
        reader[make_pair("min", 2)] = min_reader;
        reader[make_pair("max", 2)] = max_reader;
		reader[make_pair("lst", 0)] = lst_reader;
		unsigned serial = 0;
		for (auto & it : registered_functions_hack::get_registered_functions()) {
			prototype proto = make_pair(it.get_name(), it.get_nparams());
			reader[proto] = encode_serial_as_reader_func(serial);
			++serial;
		}
		initialized = true;
	}
	return reader;
}

const prototype_table& get_builtin_reader()
{
	using std::make_pair;
	static bool initialized = false;
	static prototype_table reader;
	if (!initialized) {
		
		reader[make_pair("sqrt", 1)] = sqrt_reader;
		reader[make_pair("pow", 2)] = pow_reader;
		reader[make_pair("mod", 2)] = mod_reader;
		reader[make_pair("power", 2)] = power_reader;
        reader[make_pair("min", 2)] = min_reader;
        reader[make_pair("max", 2)] = max_reader;
		reader[make_pair("lst", 0)] = lst_reader;
		enum {
			log,
			exp,
			sin,
			cos,
			tan,
			asin,
			acos,
			atan,
			sinh,
			cosh,
			tanh,
			asinh,
			acosh,
			atanh,
			atan2,
			Li2,
			Li3,
			zetaderiv,
			Li,
			S,
			H,
			lgamma,
			tgamma,
			beta,
			factorial,
			binomial,
			Order,
			NFUNCTIONS
		};
		auto it = registered_functions_hack::get_registered_functions().begin();
		unsigned serial = 0;
		for ( ; serial<NFUNCTIONS; ++it, ++serial ) {
			prototype proto = make_pair(it->get_name(), it->get_nparams());
			reader[proto] = encode_serial_as_reader_func(serial);
		}
		initialized = true;
	}
	return reader;
}

} // namespace GiNaC
extern template GiNaC::registered_class_info GiNaC::container<std::vector>::reg_info;
extern template GiNaC::registered_class_info GiNaC::container<std::list>::reg_info;
