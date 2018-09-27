//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 25 Jul 2018
//! @copyright 2018 Feel++ Consortium
//!
#include <iostream>


#include <pybind11/pybind11.h>

namespace py = pybind11;

class C : public std::enable_shared_from_this<C>
{
public:
    C() = default;
    C( C const& c ) = default;
    C( C && c ) = default;
    C& operator=( C && c ) = default;
    virtual ~C() = default;
    std::vector<double> a;
};
class A {

public:
    A() = default;
    A( A const& a ) = default;
    void hello() const { std::cout << "Hello\n"; }
    class B : public C
    
    {
    public:
        B()=default;
        B( B const& c )  = default;
        B( B && b ) = default;
        B& operator=( B&& b ) = default;
        void hi() const { std::cout << "Hi\n"; }
    };
};    
PYBIND11_MODULE(discr, m )
{
    py::class_<A,std::shared_ptr<A>>(m,"A")
        .def( py::init<>() )
        .def( "hello", &A::hello, "hello" );
    py::class_<A::B,std::shared_ptr<A::B>>(m,"B_A")
        .def( py::init<>() )
        .def( "Hi", &A::B::hi, "hi" );
}

