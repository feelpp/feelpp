/** @file function.cpp
 *
 *  Implementation of class of symbolic functions. */

/*
 *  This file was generated automatically by function.py.
 *  Please do not modify it directly, edit function.cppy instead!
 *  function.py options: maxargs=14
 *
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

#include "function.h"
#include "operators.h"
#include "fderivative.h"
#include "ex.h"
#include "lst.h"
#include "symmetry.h"
#include "print.h"
#include "power.h"
#include "archive.h"
#include "inifcns.h"
#include "utils.h"
#include "hash_seed.h"
#include "remember.h"

#include <iostream>
#include <limits>
#include <list>
#include <stdexcept>
#include <string>

namespace GiNaC {

//////////
// helper class function_options
//////////

function_options::function_options()
{
	initialize();
}

function_options::function_options(std::string const & n, std::string const & tn)
{
	initialize();
	set_name(n, tn);
}

function_options::function_options(std::string const & n, unsigned np)
{
	initialize();
	set_name(n, std::string());
	nparams = np;
}

function_options::~function_options()
{
	// nothing to clean up at the moment
}

void function_options::initialize()
{
	set_name("unnamed_function", "\\\\mbox{unnamed}");
	nparams = 0;
	eval_f = evalf_f = real_part_f = imag_part_f = conjugate_f = expand_f
		= derivative_f = expl_derivative_f = power_f = series_f = nullptr;
	info_f = nullptr;
	evalf_params_first = true;
	use_return_type = false;
	eval_use_exvector_args = false;
	evalf_use_exvector_args = false;
	conjugate_use_exvector_args = false;
	real_part_use_exvector_args = false;
	imag_part_use_exvector_args = false;
	expand_use_exvector_args = false;
	derivative_use_exvector_args = false;
	expl_derivative_use_exvector_args = false;
	power_use_exvector_args = false;
	series_use_exvector_args = false;
	print_use_exvector_args = false;
	info_use_exvector_args = false;
	use_remember = false;
	functions_with_same_name = 1;
	symtree = 0;
}

function_options & function_options::set_name(std::string const & n,
                                              std::string const & tn)
{
	name = n;
	if (tn==std::string())
		TeX_name = "\\\\mbox{"+name+"}";
	else
		TeX_name = tn;
	return *this;
}

function_options & function_options::latex_name(std::string const & tn)
{
	TeX_name = tn;
	return *this;
}

// the following lines have been generated for max. 14 parameters
function_options & function_options::eval_func(eval_funcp_1 e) 
{
	test_and_set_nparams(1);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_2 e) 
{
	test_and_set_nparams(2);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_3 e) 
{
	test_and_set_nparams(3);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_4 e) 
{
	test_and_set_nparams(4);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_5 e) 
{
	test_and_set_nparams(5);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_6 e) 
{
	test_and_set_nparams(6);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_7 e) 
{
	test_and_set_nparams(7);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_8 e) 
{
	test_and_set_nparams(8);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_9 e) 
{
	test_and_set_nparams(9);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_10 e) 
{
	test_and_set_nparams(10);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_11 e) 
{
	test_and_set_nparams(11);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_12 e) 
{
	test_and_set_nparams(12);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_13 e) 
{
	test_and_set_nparams(13);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::eval_func(eval_funcp_14 e) 
{
	test_and_set_nparams(14);
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_1 e) 
{
	test_and_set_nparams(1);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_2 e) 
{
	test_and_set_nparams(2);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_3 e) 
{
	test_and_set_nparams(3);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_4 e) 
{
	test_and_set_nparams(4);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_5 e) 
{
	test_and_set_nparams(5);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_6 e) 
{
	test_and_set_nparams(6);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_7 e) 
{
	test_and_set_nparams(7);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_8 e) 
{
	test_and_set_nparams(8);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_9 e) 
{
	test_and_set_nparams(9);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_10 e) 
{
	test_and_set_nparams(10);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_11 e) 
{
	test_and_set_nparams(11);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_12 e) 
{
	test_and_set_nparams(12);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_13 e) 
{
	test_and_set_nparams(13);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_14 e) 
{
	test_and_set_nparams(14);
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_1 e) 
{
	test_and_set_nparams(1);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_2 e) 
{
	test_and_set_nparams(2);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_3 e) 
{
	test_and_set_nparams(3);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_4 e) 
{
	test_and_set_nparams(4);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_5 e) 
{
	test_and_set_nparams(5);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_6 e) 
{
	test_and_set_nparams(6);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_7 e) 
{
	test_and_set_nparams(7);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_8 e) 
{
	test_and_set_nparams(8);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_9 e) 
{
	test_and_set_nparams(9);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_10 e) 
{
	test_and_set_nparams(10);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_11 e) 
{
	test_and_set_nparams(11);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_12 e) 
{
	test_and_set_nparams(12);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_13 e) 
{
	test_and_set_nparams(13);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_14 e) 
{
	test_and_set_nparams(14);
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_1 e) 
{
	test_and_set_nparams(1);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_2 e) 
{
	test_and_set_nparams(2);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_3 e) 
{
	test_and_set_nparams(3);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_4 e) 
{
	test_and_set_nparams(4);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_5 e) 
{
	test_and_set_nparams(5);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_6 e) 
{
	test_and_set_nparams(6);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_7 e) 
{
	test_and_set_nparams(7);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_8 e) 
{
	test_and_set_nparams(8);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_9 e) 
{
	test_and_set_nparams(9);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_10 e) 
{
	test_and_set_nparams(10);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_11 e) 
{
	test_and_set_nparams(11);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_12 e) 
{
	test_and_set_nparams(12);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_13 e) 
{
	test_and_set_nparams(13);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_14 e) 
{
	test_and_set_nparams(14);
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_1 e) 
{
	test_and_set_nparams(1);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_2 e) 
{
	test_and_set_nparams(2);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_3 e) 
{
	test_and_set_nparams(3);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_4 e) 
{
	test_and_set_nparams(4);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_5 e) 
{
	test_and_set_nparams(5);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_6 e) 
{
	test_and_set_nparams(6);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_7 e) 
{
	test_and_set_nparams(7);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_8 e) 
{
	test_and_set_nparams(8);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_9 e) 
{
	test_and_set_nparams(9);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_10 e) 
{
	test_and_set_nparams(10);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_11 e) 
{
	test_and_set_nparams(11);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_12 e) 
{
	test_and_set_nparams(12);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_13 e) 
{
	test_and_set_nparams(13);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_14 e) 
{
	test_and_set_nparams(14);
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_1 e) 
{
	test_and_set_nparams(1);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_2 e) 
{
	test_and_set_nparams(2);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_3 e) 
{
	test_and_set_nparams(3);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_4 e) 
{
	test_and_set_nparams(4);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_5 e) 
{
	test_and_set_nparams(5);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_6 e) 
{
	test_and_set_nparams(6);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_7 e) 
{
	test_and_set_nparams(7);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_8 e) 
{
	test_and_set_nparams(8);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_9 e) 
{
	test_and_set_nparams(9);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_10 e) 
{
	test_and_set_nparams(10);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_11 e) 
{
	test_and_set_nparams(11);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_12 e) 
{
	test_and_set_nparams(12);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_13 e) 
{
	test_and_set_nparams(13);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_14 e) 
{
	test_and_set_nparams(14);
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_1 e) 
{
	test_and_set_nparams(1);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_2 e) 
{
	test_and_set_nparams(2);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_3 e) 
{
	test_and_set_nparams(3);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_4 e) 
{
	test_and_set_nparams(4);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_5 e) 
{
	test_and_set_nparams(5);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_6 e) 
{
	test_and_set_nparams(6);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_7 e) 
{
	test_and_set_nparams(7);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_8 e) 
{
	test_and_set_nparams(8);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_9 e) 
{
	test_and_set_nparams(9);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_10 e) 
{
	test_and_set_nparams(10);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_11 e) 
{
	test_and_set_nparams(11);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_12 e) 
{
	test_and_set_nparams(12);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_13 e) 
{
	test_and_set_nparams(13);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_14 e) 
{
	test_and_set_nparams(14);
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_1 e) 
{
	test_and_set_nparams(1);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_2 e) 
{
	test_and_set_nparams(2);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_3 e) 
{
	test_and_set_nparams(3);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_4 e) 
{
	test_and_set_nparams(4);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_5 e) 
{
	test_and_set_nparams(5);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_6 e) 
{
	test_and_set_nparams(6);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_7 e) 
{
	test_and_set_nparams(7);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_8 e) 
{
	test_and_set_nparams(8);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_9 e) 
{
	test_and_set_nparams(9);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_10 e) 
{
	test_and_set_nparams(10);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_11 e) 
{
	test_and_set_nparams(11);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_12 e) 
{
	test_and_set_nparams(12);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_13 e) 
{
	test_and_set_nparams(13);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_14 e) 
{
	test_and_set_nparams(14);
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_1 e) 
{
	test_and_set_nparams(1);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_2 e) 
{
	test_and_set_nparams(2);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_3 e) 
{
	test_and_set_nparams(3);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_4 e) 
{
	test_and_set_nparams(4);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_5 e) 
{
	test_and_set_nparams(5);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_6 e) 
{
	test_and_set_nparams(6);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_7 e) 
{
	test_and_set_nparams(7);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_8 e) 
{
	test_and_set_nparams(8);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_9 e) 
{
	test_and_set_nparams(9);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_10 e) 
{
	test_and_set_nparams(10);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_11 e) 
{
	test_and_set_nparams(11);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_12 e) 
{
	test_and_set_nparams(12);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_13 e) 
{
	test_and_set_nparams(13);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_14 e) 
{
	test_and_set_nparams(14);
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_1 e) 
{
	test_and_set_nparams(1);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_2 e) 
{
	test_and_set_nparams(2);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_3 e) 
{
	test_and_set_nparams(3);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_4 e) 
{
	test_and_set_nparams(4);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_5 e) 
{
	test_and_set_nparams(5);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_6 e) 
{
	test_and_set_nparams(6);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_7 e) 
{
	test_and_set_nparams(7);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_8 e) 
{
	test_and_set_nparams(8);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_9 e) 
{
	test_and_set_nparams(9);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_10 e) 
{
	test_and_set_nparams(10);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_11 e) 
{
	test_and_set_nparams(11);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_12 e) 
{
	test_and_set_nparams(12);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_13 e) 
{
	test_and_set_nparams(13);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_14 e) 
{
	test_and_set_nparams(14);
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_1 e) 
{
	test_and_set_nparams(1);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_2 e) 
{
	test_and_set_nparams(2);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_3 e) 
{
	test_and_set_nparams(3);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_4 e) 
{
	test_and_set_nparams(4);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_5 e) 
{
	test_and_set_nparams(5);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_6 e) 
{
	test_and_set_nparams(6);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_7 e) 
{
	test_and_set_nparams(7);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_8 e) 
{
	test_and_set_nparams(8);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_9 e) 
{
	test_and_set_nparams(9);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_10 e) 
{
	test_and_set_nparams(10);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_11 e) 
{
	test_and_set_nparams(11);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_12 e) 
{
	test_and_set_nparams(12);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_13 e) 
{
	test_and_set_nparams(13);
	info_f = info_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_14 e) 
{
	test_and_set_nparams(14);
	info_f = info_funcp(e);
	return *this;
}
// end of generated lines

function_options & function_options::eval_func(eval_funcp_exvector e)
{
	eval_use_exvector_args = true;
	eval_f = eval_funcp(e);
	return *this;
}
function_options & function_options::evalf_func(evalf_funcp_exvector e)
{
	evalf_use_exvector_args = true;
	evalf_f = evalf_funcp(e);
	return *this;
}
function_options & function_options::conjugate_func(conjugate_funcp_exvector e)
{
	conjugate_use_exvector_args = true;
	conjugate_f = conjugate_funcp(e);
	return *this;
}
function_options & function_options::real_part_func(real_part_funcp_exvector e)
{
	real_part_use_exvector_args = true;
	real_part_f = real_part_funcp(e);
	return *this;
}
function_options & function_options::imag_part_func(imag_part_funcp_exvector e)
{
	imag_part_use_exvector_args = true;
	imag_part_f = imag_part_funcp(e);
	return *this;
}
function_options & function_options::expand_func(expand_funcp_exvector e)
{
	expand_use_exvector_args = true;
	expand_f = expand_funcp(e);
	return *this;
}
function_options & function_options::derivative_func(derivative_funcp_exvector e)
{
	derivative_use_exvector_args = true;
	derivative_f = derivative_funcp(e);
	return *this;
}
function_options & function_options::expl_derivative_func(expl_derivative_funcp_exvector e)
{
	expl_derivative_use_exvector_args = true;
	expl_derivative_f = expl_derivative_funcp(e);
	return *this;
}
function_options & function_options::power_func(power_funcp_exvector e)
{
	power_use_exvector_args = true;
	power_f = power_funcp(e);
	return *this;
}
function_options & function_options::series_func(series_funcp_exvector e)
{
	series_use_exvector_args = true;
	series_f = series_funcp(e);
	return *this;
}
function_options & function_options::info_func(info_funcp_exvector e)
{
	info_use_exvector_args = true;
	info_f = info_funcp(e);
	return *this;
}

// end of generated lines

function_options & function_options::set_return_type(unsigned rt, const return_type_t* rtt)
{
	use_return_type = true;
	return_type = rt;
	if (rtt != nullptr)
		return_type_tinfo = *rtt;
	else
		return_type_tinfo = make_return_type_t<function>();
	return *this;
}

function_options & function_options::do_not_evalf_params()
{
	evalf_params_first = false;
	return *this;
}

function_options & function_options::remember(unsigned size,
                                              unsigned assoc_size,
                                              unsigned strategy)
{
	use_remember = true;
	remember_size = size;
	remember_assoc_size = assoc_size;
	remember_strategy = strategy;
	return *this;
}

function_options & function_options::overloaded(unsigned o)
{
	functions_with_same_name = o;
	return *this;
}

function_options & function_options::set_symmetry(const symmetry & s)
{
	symtree = s;
	return *this;
}
	
void function_options::test_and_set_nparams(unsigned n)
{
	if (nparams==0) {
		nparams = n;
	} else if (nparams!=n) {
		// we do not throw an exception here because this code is
		// usually executed before main(), so the exception could not
		// be caught anyhow
		std::cerr << "WARNING: " << name << "(): number of parameters ("
		          << n << ") differs from number set before (" 
		          << nparams << ")" << std::endl;
	}
}

void function_options::set_print_func(unsigned id, print_funcp f)
{
	if (id >= print_dispatch_table.size())
		print_dispatch_table.resize(id + 1);
	print_dispatch_table[id] = f;
}

/** This can be used as a hook for external applications. */
unsigned function::current_serial = 0;


GINAC_IMPLEMENT_REGISTERED_CLASS(function, exprseq)

//////////
// default constructor
//////////

// public

function::function() : serial(0)
{
}

//////////
// other constructors
//////////

// public

function::function(unsigned ser) : serial(ser)
{
}

// the following lines have been generated for max. 14 parameters
function::function(unsigned ser, const ex & param1)
	: exprseq{param1}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2)
	: exprseq{param1, param2}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3)
	: exprseq{param1, param2, param3}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4)
	: exprseq{param1, param2, param3, param4}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5)
	: exprseq{param1, param2, param3, param4, param5}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6)
	: exprseq{param1, param2, param3, param4, param5, param6}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7)
	: exprseq{param1, param2, param3, param4, param5, param6, param7}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8)
	: exprseq{param1, param2, param3, param4, param5, param6, param7, param8}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9)
	: exprseq{param1, param2, param3, param4, param5, param6, param7, param8, param9}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10)
	: exprseq{param1, param2, param3, param4, param5, param6, param7, param8, param9, param10}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11)
	: exprseq{param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, param11}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11, const ex & param12)
	: exprseq{param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, param11, param12}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11, const ex & param12, const ex & param13)
	: exprseq{param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, param11, param12, param13}, serial(ser)
{
}
function::function(unsigned ser, const ex & param1, const ex & param2, const ex & param3, const ex & param4, const ex & param5, const ex & param6, const ex & param7, const ex & param8, const ex & param9, const ex & param10, const ex & param11, const ex & param12, const ex & param13, const ex & param14)
	: exprseq{param1, param2, param3, param4, param5, param6, param7, param8, param9, param10, param11, param12, param13, param14}, serial(ser)
{
}

function::function(unsigned ser, const exprseq & es) : exprseq(es), serial(ser)
{

	// Force re-evaluation even if the exprseq was already evaluated
	// (the exprseq copy constructor copies the flags)
	clearflag(status_flags::evaluated);
}

function::function(unsigned ser, const exvector & v)
  : exprseq(v), serial(ser)
{
}

function::function(unsigned ser, exvector && v)
  : exprseq(std::move(v)), serial(ser)
{
}

//////////
// archiving
//////////

/** Construct object from archive_node. */
void function::read_archive(const archive_node& n, lst& sym_lst)
{
	inherited::read_archive(n, sym_lst);
	// Find serial number by function name
	std::string s;
	if (n.find_string("name", s)) {
		unsigned int ser = 0;
		for (auto & it : registered_functions()) {
			if (s == it.name) {
				serial = ser;
				return;
			}
			++ser;
		}
		throw (std::runtime_error("unknown function '" + s + "' in archive"));
	} else
		throw (std::runtime_error("unnamed function in archive"));
}

/** Archive the object. */
void function::archive(archive_node &n) const
{
	inherited::archive(n);
	GINAC_ASSERT(serial < registered_functions().size());
	n.add_string("name", registered_functions()[serial].name);
}

GINAC_BIND_UNARCHIVER(function);

//////////
// functions overriding virtual functions from base classes
//////////

// public

void function::print(const print_context & c, unsigned level) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];
	const std::vector<print_funcp> &pdt = opt.print_dispatch_table;

	// Dynamically dispatch on print_context type
	const print_context_class_info *pc_info = &c.get_class_info();

next_context:
	unsigned id = pc_info->options.get_id();
	if (id >= pdt.size() || pdt[id] == nullptr) {

		// Method not found, try parent print_context class
		const print_context_class_info *parent_pc_info = pc_info->get_parent();
		if (parent_pc_info) {
			pc_info = parent_pc_info;
			goto next_context;
		}

		// Method still not found, use default output
		if (is_a<print_tree>(c)) {

			c.s << std::string(level, ' ') << class_name() << " "
			    << opt.name << " @" << this
			    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
			    << ", nops=" << nops()
			    << std::endl;
			unsigned delta_indent = static_cast<const print_tree &>(c).delta_indent;
			for (size_t i=0; i<seq.size(); ++i)
				seq[i].print(c, level + delta_indent);
			c.s << std::string(level + delta_indent, ' ') << "=====" << std::endl;

		} else if (is_a<print_csrc>(c)) {

			// Print function name in lowercase
			std::string lname = opt.name;
			size_t num = lname.size();
			for (size_t i=0; i<num; i++)
				lname[i] = tolower(lname[i]);
			c.s << lname;
			printseq(c, '(', ',', ')', exprseq::precedence(), function::precedence());

		} else if (is_a<print_latex>(c)) {
			c.s << opt.TeX_name;
			printseq(c, '(', ',', ')', exprseq::precedence(), function::precedence());
		} else {
			c.s << opt.name;
			printseq(c, '(', ',', ')', exprseq::precedence(), function::precedence());
		}

	} else {

		// Method found, call it
		current_serial = serial;
		if (opt.print_use_exvector_args)
			((print_funcp_exvector)pdt[id])(seq, c);
		else switch (opt.nparams) {
			// the following lines have been generated for max. 14 parameters
			case 1:
				((print_funcp_1)(pdt[id]))(seq[0], c);
				break;
			case 2:
				((print_funcp_2)(pdt[id]))(seq[0], seq[1], c);
				break;
			case 3:
				((print_funcp_3)(pdt[id]))(seq[0], seq[1], seq[2], c);
				break;
			case 4:
				((print_funcp_4)(pdt[id]))(seq[0], seq[1], seq[2], seq[3], c);
				break;
			case 5:
				((print_funcp_5)(pdt[id]))(seq[0], seq[1], seq[2], seq[3], seq[4], c);
				break;
			case 6:
				((print_funcp_6)(pdt[id]))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], c);
				break;
			case 7:
				((print_funcp_7)(pdt[id]))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], c);
				break;
			case 8:
				((print_funcp_8)(pdt[id]))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], c);
				break;
			case 9:
				((print_funcp_9)(pdt[id]))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], c);
				break;
			case 10:
				((print_funcp_10)(pdt[id]))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], c);
				break;
			case 11:
				((print_funcp_11)(pdt[id]))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], c);
				break;
			case 12:
				((print_funcp_12)(pdt[id]))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], c);
				break;
			case 13:
				((print_funcp_13)(pdt[id]))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], c);
				break;
			case 14:
				((print_funcp_14)(pdt[id]))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], seq[13], c);
				break;
			// end of generated lines
		default:
			throw(std::logic_error("function::print(): invalid nparams"));
		}
	}
}

ex function::eval() const
{
	if (flags & status_flags::evaluated) {
		return *this;
	}

	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	// Canonicalize argument order according to the symmetry properties
	if (seq.size() > 1 && !(opt.symtree.is_zero())) {
		exvector v = seq;
		GINAC_ASSERT(is_a<symmetry>(opt.symtree));
		int sig = canonicalize(v.begin(), ex_to<symmetry>(opt.symtree));
		if (sig != std::numeric_limits<int>::max()) {
			// Something has changed while sorting arguments, more evaluations later
			if (sig == 0)
				return _ex0;
			return ex(sig) * thiscontainer(std::move(v));
		}
	}

	if (opt.eval_f==nullptr) {
		return this->hold();
	}

	bool use_remember = opt.use_remember;
	ex eval_result;
	if (use_remember && lookup_remember_table(eval_result)) {
		return eval_result;
	}
	current_serial = serial;
	if (opt.eval_use_exvector_args)
		eval_result = ((eval_funcp_exvector)(opt.eval_f))(seq);
	else
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
		case 1:
			eval_result = ((eval_funcp_1)(opt.eval_f))(seq[0]);
			break;
		case 2:
			eval_result = ((eval_funcp_2)(opt.eval_f))(seq[0], seq[1]);
			break;
		case 3:
			eval_result = ((eval_funcp_3)(opt.eval_f))(seq[0], seq[1], seq[2]);
			break;
		case 4:
			eval_result = ((eval_funcp_4)(opt.eval_f))(seq[0], seq[1], seq[2], seq[3]);
			break;
		case 5:
			eval_result = ((eval_funcp_5)(opt.eval_f))(seq[0], seq[1], seq[2], seq[3], seq[4]);
			break;
		case 6:
			eval_result = ((eval_funcp_6)(opt.eval_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5]);
			break;
		case 7:
			eval_result = ((eval_funcp_7)(opt.eval_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6]);
			break;
		case 8:
			eval_result = ((eval_funcp_8)(opt.eval_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7]);
			break;
		case 9:
			eval_result = ((eval_funcp_9)(opt.eval_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8]);
			break;
		case 10:
			eval_result = ((eval_funcp_10)(opt.eval_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9]);
			break;
		case 11:
			eval_result = ((eval_funcp_11)(opt.eval_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10]);
			break;
		case 12:
			eval_result = ((eval_funcp_12)(opt.eval_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11]);
			break;
		case 13:
			eval_result = ((eval_funcp_13)(opt.eval_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12]);
			break;
		case 14:
			eval_result = ((eval_funcp_14)(opt.eval_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], seq[13]);
			break;
		// end of generated lines
	default:
		throw(std::logic_error("function::eval(): invalid nparams"));
	}
	if (use_remember) {
		store_remember_table(eval_result);
	}
	return eval_result;
}

ex function::evalf() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	// Evaluate children first
	exvector eseq;
	if (!opt.evalf_params_first)
		eseq = seq;
	else {
		eseq.reserve(seq.size());
		for (auto & it : seq) {
			eseq.push_back(it.evalf());
		}
	}

	if (opt.evalf_f==nullptr) {
		return function(serial,eseq).hold();
	}
	current_serial = serial;
	if (opt.evalf_use_exvector_args)
		return ((evalf_funcp_exvector)(opt.evalf_f))(seq);
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
		case 1:
			return ((evalf_funcp_1)(opt.evalf_f))(eseq[0]);
		case 2:
			return ((evalf_funcp_2)(opt.evalf_f))(eseq[0], eseq[1]);
		case 3:
			return ((evalf_funcp_3)(opt.evalf_f))(eseq[0], eseq[1], eseq[2]);
		case 4:
			return ((evalf_funcp_4)(opt.evalf_f))(eseq[0], eseq[1], eseq[2], eseq[3]);
		case 5:
			return ((evalf_funcp_5)(opt.evalf_f))(eseq[0], eseq[1], eseq[2], eseq[3], eseq[4]);
		case 6:
			return ((evalf_funcp_6)(opt.evalf_f))(eseq[0], eseq[1], eseq[2], eseq[3], eseq[4], eseq[5]);
		case 7:
			return ((evalf_funcp_7)(opt.evalf_f))(eseq[0], eseq[1], eseq[2], eseq[3], eseq[4], eseq[5], eseq[6]);
		case 8:
			return ((evalf_funcp_8)(opt.evalf_f))(eseq[0], eseq[1], eseq[2], eseq[3], eseq[4], eseq[5], eseq[6], eseq[7]);
		case 9:
			return ((evalf_funcp_9)(opt.evalf_f))(eseq[0], eseq[1], eseq[2], eseq[3], eseq[4], eseq[5], eseq[6], eseq[7], eseq[8]);
		case 10:
			return ((evalf_funcp_10)(opt.evalf_f))(eseq[0], eseq[1], eseq[2], eseq[3], eseq[4], eseq[5], eseq[6], eseq[7], eseq[8], eseq[9]);
		case 11:
			return ((evalf_funcp_11)(opt.evalf_f))(eseq[0], eseq[1], eseq[2], eseq[3], eseq[4], eseq[5], eseq[6], eseq[7], eseq[8], eseq[9], eseq[10]);
		case 12:
			return ((evalf_funcp_12)(opt.evalf_f))(eseq[0], eseq[1], eseq[2], eseq[3], eseq[4], eseq[5], eseq[6], eseq[7], eseq[8], eseq[9], eseq[10], eseq[11]);
		case 13:
			return ((evalf_funcp_13)(opt.evalf_f))(eseq[0], eseq[1], eseq[2], eseq[3], eseq[4], eseq[5], eseq[6], eseq[7], eseq[8], eseq[9], eseq[10], eseq[11], eseq[12]);
		case 14:
			return ((evalf_funcp_14)(opt.evalf_f))(eseq[0], eseq[1], eseq[2], eseq[3], eseq[4], eseq[5], eseq[6], eseq[7], eseq[8], eseq[9], eseq[10], eseq[11], eseq[12], eseq[13]);
		// end of generated lines
	}
	throw(std::logic_error("function::evalf(): invalid nparams"));
}

/**
 *  This method is defined to be in line with behavior of function::return_type()
 */
ex function::eval_ncmul(const exvector & v) const
{
	// If this function is called then the list of arguments is non-empty
	// and the first argument is non-commutative, see  function::return_type()
	return seq.begin()->eval_ncmul(v);
}

unsigned function::calchash() const
{
	unsigned v = golden_ratio_hash(make_hash_seed(typeid(*this)) ^ serial);
	for (size_t i=0; i<nops(); i++) {
		v = rotate_left(v);
		v ^= this->op(i).gethash();
	}

	if (flags & status_flags::evaluated) {
		setflag(status_flags::hash_calculated);
		hashvalue = v;
	}
	return v;
}

ex function::thiscontainer(const exvector & v) const
{
	return function(serial, v);
}

ex function::thiscontainer(exvector && v) const
{
	return function(serial, std::move(v));
}

/** Implementation of ex::series for functions.
 *  \@see ex::series */
ex function::series(const relational & r, int order, unsigned options) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.series_f==nullptr) {
		return basic::series(r, order);
	}
	ex res;
	current_serial = serial;
	if (opt.series_use_exvector_args) {
		try {
			res = ((series_funcp_exvector)(opt.series_f))(seq, r, order, options);
		} catch (do_taylor) {
			res = basic::series(r, order, options);
		}
		return res;
	}
	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
		case 1:
			try {
				res = ((series_funcp_1)(opt.series_f))(seq[0], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 2:
			try {
				res = ((series_funcp_2)(opt.series_f))(seq[0], seq[1], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 3:
			try {
				res = ((series_funcp_3)(opt.series_f))(seq[0], seq[1], seq[2], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 4:
			try {
				res = ((series_funcp_4)(opt.series_f))(seq[0], seq[1], seq[2], seq[3], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 5:
			try {
				res = ((series_funcp_5)(opt.series_f))(seq[0], seq[1], seq[2], seq[3], seq[4], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 6:
			try {
				res = ((series_funcp_6)(opt.series_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 7:
			try {
				res = ((series_funcp_7)(opt.series_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 8:
			try {
				res = ((series_funcp_8)(opt.series_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 9:
			try {
				res = ((series_funcp_9)(opt.series_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 10:
			try {
				res = ((series_funcp_10)(opt.series_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 11:
			try {
				res = ((series_funcp_11)(opt.series_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 12:
			try {
				res = ((series_funcp_12)(opt.series_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 13:
			try {
				res = ((series_funcp_13)(opt.series_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		case 14:
			try {
				res = ((series_funcp_14)(opt.series_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], seq[13], r, order, options);
			} catch (do_taylor) {
				res = basic::series(r, order, options);
			}
			return res;
		// end of generated lines
	}
	throw(std::logic_error("function::series(): invalid nparams"));
}

/** Implementation of ex::conjugate for functions. */
ex function::conjugate() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (opt.conjugate_f==nullptr) {
		return conjugate_function(*this).hold();
	}

	if (opt.conjugate_use_exvector_args) {
		return ((conjugate_funcp_exvector)(opt.conjugate_f))(seq);
	}

	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
		case 1:
			return ((conjugate_funcp_1)(opt.conjugate_f))(seq[0]);
		case 2:
			return ((conjugate_funcp_2)(opt.conjugate_f))(seq[0], seq[1]);
		case 3:
			return ((conjugate_funcp_3)(opt.conjugate_f))(seq[0], seq[1], seq[2]);
		case 4:
			return ((conjugate_funcp_4)(opt.conjugate_f))(seq[0], seq[1], seq[2], seq[3]);
		case 5:
			return ((conjugate_funcp_5)(opt.conjugate_f))(seq[0], seq[1], seq[2], seq[3], seq[4]);
		case 6:
			return ((conjugate_funcp_6)(opt.conjugate_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5]);
		case 7:
			return ((conjugate_funcp_7)(opt.conjugate_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6]);
		case 8:
			return ((conjugate_funcp_8)(opt.conjugate_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7]);
		case 9:
			return ((conjugate_funcp_9)(opt.conjugate_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8]);
		case 10:
			return ((conjugate_funcp_10)(opt.conjugate_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9]);
		case 11:
			return ((conjugate_funcp_11)(opt.conjugate_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10]);
		case 12:
			return ((conjugate_funcp_12)(opt.conjugate_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11]);
		case 13:
			return ((conjugate_funcp_13)(opt.conjugate_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12]);
		case 14:
			return ((conjugate_funcp_14)(opt.conjugate_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], seq[13]);
		// end of generated lines
	}
	throw(std::logic_error("function::conjugate(): invalid nparams"));
}

/** Implementation of ex::real_part for functions. */
ex function::real_part() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (opt.real_part_f==nullptr)
		return basic::real_part();

	if (opt.real_part_use_exvector_args)
		return ((real_part_funcp_exvector)(opt.real_part_f))(seq);

	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
		case 1:
			return ((real_part_funcp_1)(opt.real_part_f))(seq[0]);
		case 2:
			return ((real_part_funcp_2)(opt.real_part_f))(seq[0], seq[1]);
		case 3:
			return ((real_part_funcp_3)(opt.real_part_f))(seq[0], seq[1], seq[2]);
		case 4:
			return ((real_part_funcp_4)(opt.real_part_f))(seq[0], seq[1], seq[2], seq[3]);
		case 5:
			return ((real_part_funcp_5)(opt.real_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4]);
		case 6:
			return ((real_part_funcp_6)(opt.real_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5]);
		case 7:
			return ((real_part_funcp_7)(opt.real_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6]);
		case 8:
			return ((real_part_funcp_8)(opt.real_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7]);
		case 9:
			return ((real_part_funcp_9)(opt.real_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8]);
		case 10:
			return ((real_part_funcp_10)(opt.real_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9]);
		case 11:
			return ((real_part_funcp_11)(opt.real_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10]);
		case 12:
			return ((real_part_funcp_12)(opt.real_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11]);
		case 13:
			return ((real_part_funcp_13)(opt.real_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12]);
		case 14:
			return ((real_part_funcp_14)(opt.real_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], seq[13]);
		// end of generated lines
	}
	throw(std::logic_error("function::real_part(): invalid nparams"));
}

/** Implementation of ex::imag_part for functions. */
ex function::imag_part() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (opt.imag_part_f==nullptr)
		return basic::imag_part();

	if (opt.imag_part_use_exvector_args)
		return ((imag_part_funcp_exvector)(opt.imag_part_f))(seq);

	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
		case 1:
			return ((imag_part_funcp_1)(opt.imag_part_f))(seq[0]);
		case 2:
			return ((imag_part_funcp_2)(opt.imag_part_f))(seq[0], seq[1]);
		case 3:
			return ((imag_part_funcp_3)(opt.imag_part_f))(seq[0], seq[1], seq[2]);
		case 4:
			return ((imag_part_funcp_4)(opt.imag_part_f))(seq[0], seq[1], seq[2], seq[3]);
		case 5:
			return ((imag_part_funcp_5)(opt.imag_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4]);
		case 6:
			return ((imag_part_funcp_6)(opt.imag_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5]);
		case 7:
			return ((imag_part_funcp_7)(opt.imag_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6]);
		case 8:
			return ((imag_part_funcp_8)(opt.imag_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7]);
		case 9:
			return ((imag_part_funcp_9)(opt.imag_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8]);
		case 10:
			return ((imag_part_funcp_10)(opt.imag_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9]);
		case 11:
			return ((imag_part_funcp_11)(opt.imag_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10]);
		case 12:
			return ((imag_part_funcp_12)(opt.imag_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11]);
		case 13:
			return ((imag_part_funcp_13)(opt.imag_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12]);
		case 14:
			return ((imag_part_funcp_14)(opt.imag_part_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], seq[13]);
		// end of generated lines
	}
	throw(std::logic_error("function::imag_part(): invalid nparams"));
}

/** Implementation of ex::info for functions. */
bool function::info(unsigned inf) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (opt.info_f==nullptr) {
		return basic::info(inf);
	}

	if (opt.info_use_exvector_args) {
		return ((info_funcp_exvector)(opt.info_f))(seq, inf);
	}

	switch (opt.nparams) {
		// the following lines have been generated for max. 14 parameters
		case 1:
			return ((info_funcp_1)(opt.info_f))(seq[0], inf);
		case 2:
			return ((info_funcp_2)(opt.info_f))(seq[0], seq[1], inf);
		case 3:
			return ((info_funcp_3)(opt.info_f))(seq[0], seq[1], seq[2], inf);
		case 4:
			return ((info_funcp_4)(opt.info_f))(seq[0], seq[1], seq[2], seq[3], inf);
		case 5:
			return ((info_funcp_5)(opt.info_f))(seq[0], seq[1], seq[2], seq[3], seq[4], inf);
		case 6:
			return ((info_funcp_6)(opt.info_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], inf);
		case 7:
			return ((info_funcp_7)(opt.info_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], inf);
		case 8:
			return ((info_funcp_8)(opt.info_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], inf);
		case 9:
			return ((info_funcp_9)(opt.info_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], inf);
		case 10:
			return ((info_funcp_10)(opt.info_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], inf);
		case 11:
			return ((info_funcp_11)(opt.info_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], inf);
		case 12:
			return ((info_funcp_12)(opt.info_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], inf);
		case 13:
			return ((info_funcp_13)(opt.info_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], inf);
		case 14:
			return ((info_funcp_14)(opt.info_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], seq[13], inf);
		// end of generated lines
	}
	throw(std::logic_error("function::info(): invalid nparams"));
}

// protected

/** Implementation of ex::diff() for functions. It applies the chain rule,
 *  except for the Order term function.
 *  \@see ex::diff */
ex function::derivative(const symbol & s) const
{
	ex result;

	try {
		// Explicit derivation
		result = expl_derivative(s);
	} catch (...) {
		// Chain rule
		ex arg_diff;
		size_t num = seq.size();
		for (size_t i=0; i<num; i++) {
			arg_diff = seq[i].diff(s);
			// We apply the chain rule only when it makes sense.  This is not
			// just for performance reasons but also to allow functions to
			// throw when differentiated with respect to one of its arguments
			// without running into trouble with our automatic full
			// differentiation:
			if (!arg_diff.is_zero())
				result += pderivative(i)*arg_diff;
		}
	}
	return result;
}

int function::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	if (serial != o.serial)
		return serial < o.serial ? -1 : 1;
	else
		return exprseq::compare_same_type(o);
}

bool function::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	if (serial != o.serial)
		return false;
	else
		return exprseq::is_equal_same_type(o);
}

bool function::match_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	return serial == o.serial;
}

unsigned function::return_type() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.use_return_type) {
		// Return type was explicitly specified
		return opt.return_type;
	} else {
		// Default behavior is to use the return type of the first
		// argument. Thus, exp() of a matrix behaves like a matrix, etc.
		if (seq.empty())
			return return_types::commutative;
		else
			return seq.begin()->return_type();
	}
}

return_type_t function::return_type_tinfo() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.use_return_type) {
		// Return type was explicitly specified
		return opt.return_type_tinfo;
	} else {
		// Default behavior is to use the return type of the first
		// argument. Thus, exp() of a matrix behaves like a matrix, etc.
		if (seq.empty())
			return make_return_type_t<function>();
		else
			return seq.begin()->return_type_tinfo();
	}
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

// protected

ex function::pderivative(unsigned diff_param) const // partial differentiation
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];
	
	if (opt.derivative_f) {
		// Invoke the defined derivative function.
		current_serial = serial;
		if (opt.derivative_use_exvector_args)
			return ((derivative_funcp_exvector)(opt.derivative_f))(seq, diff_param);
		switch (opt.nparams) {
			// the following lines have been generated for max. 14 parameters
			case 1:
				return ((derivative_funcp_1)(opt.derivative_f))(seq[0], diff_param);
			case 2:
				return ((derivative_funcp_2)(opt.derivative_f))(seq[0], seq[1], diff_param);
			case 3:
				return ((derivative_funcp_3)(opt.derivative_f))(seq[0], seq[1], seq[2], diff_param);
			case 4:
				return ((derivative_funcp_4)(opt.derivative_f))(seq[0], seq[1], seq[2], seq[3], diff_param);
			case 5:
				return ((derivative_funcp_5)(opt.derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], diff_param);
			case 6:
				return ((derivative_funcp_6)(opt.derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], diff_param);
			case 7:
				return ((derivative_funcp_7)(opt.derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], diff_param);
			case 8:
				return ((derivative_funcp_8)(opt.derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], diff_param);
			case 9:
				return ((derivative_funcp_9)(opt.derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], diff_param);
			case 10:
				return ((derivative_funcp_10)(opt.derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], diff_param);
			case 11:
				return ((derivative_funcp_11)(opt.derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], diff_param);
			case 12:
				return ((derivative_funcp_12)(opt.derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], diff_param);
			case 13:
				return ((derivative_funcp_13)(opt.derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], diff_param);
			case 14:
				return ((derivative_funcp_14)(opt.derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], seq[13], diff_param);
			// end of generated lines
		}
	}
	// No derivative defined? Fall back to abstract derivative object.
	return fderivative(serial, diff_param, seq);
}

ex function::expl_derivative(const symbol & s) const // explicit differentiation
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.expl_derivative_f) {
		// Invoke the defined explicit derivative function.
		current_serial = serial;
		if (opt.expl_derivative_use_exvector_args)
			return ((expl_derivative_funcp_exvector)(opt.expl_derivative_f))(seq, s);
		switch (opt.nparams) {
			// the following lines have been generated for max. 14 parameters
			case 1:
				return ((expl_derivative_funcp_1)(opt.expl_derivative_f))(seq[0], s);
			case 2:
				return ((expl_derivative_funcp_2)(opt.expl_derivative_f))(seq[0], seq[1], s);
			case 3:
				return ((expl_derivative_funcp_3)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], s);
			case 4:
				return ((expl_derivative_funcp_4)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], seq[3], s);
			case 5:
				return ((expl_derivative_funcp_5)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], s);
			case 6:
				return ((expl_derivative_funcp_6)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], s);
			case 7:
				return ((expl_derivative_funcp_7)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], s);
			case 8:
				return ((expl_derivative_funcp_8)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], s);
			case 9:
				return ((expl_derivative_funcp_9)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], s);
			case 10:
				return ((expl_derivative_funcp_10)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], s);
			case 11:
				return ((expl_derivative_funcp_11)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], s);
			case 12:
				return ((expl_derivative_funcp_12)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], s);
			case 13:
				return ((expl_derivative_funcp_13)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], s);
			case 14:
				return ((expl_derivative_funcp_14)(opt.expl_derivative_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], seq[13], s);
			// end of generated lines
		}
	}
	// There is no fallback for explicit derivative.
	throw(std::logic_error("function::expl_derivative(): explicit derivation is called, but no such function defined"));
}

ex function::power(const ex & power_param) const // power of function
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];
	
	if (opt.power_f) {
		// Invoke the defined power function.
		current_serial = serial;
		if (opt.power_use_exvector_args)
			return ((power_funcp_exvector)(opt.power_f))(seq,  power_param);
		switch (opt.nparams) {
			// the following lines have been generated for max. 14 parameters
			case 1:
				return ((power_funcp_1)(opt.power_f))(seq[0], power_param);
			case 2:
				return ((power_funcp_2)(opt.power_f))(seq[0], seq[1], power_param);
			case 3:
				return ((power_funcp_3)(opt.power_f))(seq[0], seq[1], seq[2], power_param);
			case 4:
				return ((power_funcp_4)(opt.power_f))(seq[0], seq[1], seq[2], seq[3], power_param);
			case 5:
				return ((power_funcp_5)(opt.power_f))(seq[0], seq[1], seq[2], seq[3], seq[4], power_param);
			case 6:
				return ((power_funcp_6)(opt.power_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], power_param);
			case 7:
				return ((power_funcp_7)(opt.power_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], power_param);
			case 8:
				return ((power_funcp_8)(opt.power_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], power_param);
			case 9:
				return ((power_funcp_9)(opt.power_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], power_param);
			case 10:
				return ((power_funcp_10)(opt.power_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], power_param);
			case 11:
				return ((power_funcp_11)(opt.power_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], power_param);
			case 12:
				return ((power_funcp_12)(opt.power_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], power_param);
			case 13:
				return ((power_funcp_13)(opt.power_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], power_param);
			case 14:
				return ((power_funcp_14)(opt.power_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], seq[13], power_param);
			// end of generated lines
		}
	}
	// No power function defined? Fall back to returning a power object.
	return dynallocate<GiNaC::power>(*this, power_param).setflag(status_flags::evaluated);
}

ex function::expand(unsigned options) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.expand_f) {
		// Invoke the defined expand function.
		current_serial = serial;
		if (opt.expand_use_exvector_args)
			return ((expand_funcp_exvector)(opt.expand_f))(seq,  options);
		switch (opt.nparams) {
			// the following lines have been generated for max. 14 parameters
			case 1:
				return ((expand_funcp_1)(opt.expand_f))(seq[0], options);
			case 2:
				return ((expand_funcp_2)(opt.expand_f))(seq[0], seq[1], options);
			case 3:
				return ((expand_funcp_3)(opt.expand_f))(seq[0], seq[1], seq[2], options);
			case 4:
				return ((expand_funcp_4)(opt.expand_f))(seq[0], seq[1], seq[2], seq[3], options);
			case 5:
				return ((expand_funcp_5)(opt.expand_f))(seq[0], seq[1], seq[2], seq[3], seq[4], options);
			case 6:
				return ((expand_funcp_6)(opt.expand_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], options);
			case 7:
				return ((expand_funcp_7)(opt.expand_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], options);
			case 8:
				return ((expand_funcp_8)(opt.expand_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], options);
			case 9:
				return ((expand_funcp_9)(opt.expand_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], options);
			case 10:
				return ((expand_funcp_10)(opt.expand_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], options);
			case 11:
				return ((expand_funcp_11)(opt.expand_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], options);
			case 12:
				return ((expand_funcp_12)(opt.expand_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], options);
			case 13:
				return ((expand_funcp_13)(opt.expand_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], options);
			case 14:
				return ((expand_funcp_14)(opt.expand_f))(seq[0], seq[1], seq[2], seq[3], seq[4], seq[5], seq[6], seq[7], seq[8], seq[9], seq[10], seq[11], seq[12], seq[13], options);
			// end of generated lines
		}
	}
	// No expand function defined? Return the same function with expanded arguments (if required)
	if (options & expand_options::expand_function_args)
		return inherited::expand(options);
	else
		return (options == 0) ? setflag(status_flags::expanded) : *this;
}

std::vector<function_options> & function::registered_functions()
{
	static std::vector<function_options> rf = std::vector<function_options>();
	return rf;
}

bool function::lookup_remember_table(ex & result) const
{
	return remember_table::remember_tables()[this->serial].lookup_entry(*this,result);
}

void function::store_remember_table(ex const & result) const
{
	remember_table::remember_tables()[this->serial].add_entry(*this,result);
}

// public

unsigned function::register_new(function_options const & opt)
{
	size_t same_name = 0;
	for (auto & i : registered_functions()) {
		if (i.name==opt.name) {
			++same_name;
		}
	}
	if (same_name>=opt.functions_with_same_name) {
		// we do not throw an exception here because this code is
		// usually executed before main(), so the exception could not
		// caught anyhow
		std::cerr << "WARNING: function name " << opt.name
		          << " already in use!" << std::endl;
	}
	registered_functions().push_back(opt);
	if (opt.use_remember) {
		remember_table::remember_tables().
			push_back(remember_table(opt.remember_size,
			                         opt.remember_assoc_size,
			                         opt.remember_strategy));
	} else {
		remember_table::remember_tables().push_back(remember_table());
	}
	return registered_functions().size()-1;
}

/** Find serial number of function by name and number of parameters.
 *  Throws exception if function was not found. */
unsigned function::find_function(const std::string &name, unsigned nparams)
{
	unsigned serial = 0;
	for (auto & it : function::registered_functions()) {
		if (it.get_name() == name && it.get_nparams() == nparams)
			return serial;
		++serial;
	}
	throw (std::runtime_error("no function '" + name + "' with " + std::to_string(nparams) + " parameters defined"));
}

/** Return the print name of the function. */
std::string function::get_name() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	return registered_functions()[serial].name;
}

} // namespace GiNaC

extern template GiNaC::registered_class_info GiNaC::container<std::vector>::reg_info;
