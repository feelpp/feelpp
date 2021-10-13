/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-10-04

  Copyright (C) 2004 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_factory.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-10-04
 */
#include <iostream>

#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

class A
{
public:
    virtual ~A() {}
    virtual const char* hello() const = 0;
};



class B
    : public A
{
public:
    ~B() override {}
    const char* hello() const override
    {
        return "hello";
    }
};

class C
    : public A
{
public:
    ~C() override {}
    const char* hello() const override
    {
        return "hie";
    }
};

class D
    : public B
{
public:
    const char* hello() const override
    {
        return "Yo";
    }
};

class E
    : public C
{
public:
    const char* hello() const override
    {
        return "Ciao";
    }
};
class F
    : public D
{
public:
    F()
        : D(), str( "salut" )
    {}
    F( const F& f )
        : D(), str( f.str )
    {
        VLOG(1) << "calling F::copy constructor\n";
    }
    const char* hello() const override
    {
        return str.c_str();
    }
    std::string str;
};

typedef Feel::Singleton< Feel::Factory< A, std::string > > AFactory;
typedef Feel::Singleton< Feel::FactoryClone< A > > AFactoryClone;

namespace
{
std::unique_ptr<A> createB()
{
    return std::make_unique< B>();
}
std::unique_ptr<A> createC()
{
    return std::make_unique< C>();
}
std::unique_ptr<A> createD()
{
    return std::make_unique< D>();
}
std::unique_ptr<A> createE()
{
    return std::make_unique< E>();
}
std::unique_ptr<A> createF()
{
    return std::make_unique< F>();
}
const bool regB = AFactory::instance().registerProduct( "B", &createB );
const bool regC = AFactory::instance().registerProduct( "C", &createC );
const bool regD = AFactory::instance().registerProduct( "D", &createD );
const bool regE = AFactory::instance().registerProduct( "E", &createE );
const bool regF = AFactory::instance().registerProduct( "F", &createF );

// cloning: dolly is not far away ;)
std::unique_ptr<A> createFc( A const* f )
{
    return std::make_unique< F>( ( F const& )*f );
}
const bool regFc = AFactoryClone::instance().registerProduct( typeid( F ), &createFc );
}

int
main( int /*argc*/, char** /*argv*/ )
{
    std::cerr << "B hello must be hello : " << AFactory::instance().createObject( "B" )->hello() << "\n";
    std::cerr << "C hello must be hie   : " << AFactory::instance().createObject( "C" )->hello() << "\n";
    std::cerr << "D hello must be Yo    : " << AFactory::instance().createObject( "D" )->hello() << "\n";
    std::cerr << "E hello must be Ciao  : " << AFactory::instance().createObject( "E" )->hello() << "\n";
    std::cerr << "F hello must be salut : " << AFactory::instance().createObject( "F" )->hello() << "\n";

    F f;
    std::cerr << "Clone F hello must be salut : " << AFactoryClone::instance().createObject( &f )->hello() << "\n";
}
