//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! Date: 03 Jan 2017
//!
//! Copyright (C) 2017 Feel++ Consortium
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
#include <feel/feel.hpp>

struct f_t2
{
    int i;
    f_t2( int _i ) : i (_i) {}
    f_t2( f_t2 const& f ) : i( f.i ) {}
    f_t2( f_t2 && f ) : i( std::move( f.i) ) {}
    f_t2& operator=( f_t2&& f ) { i = std::move( f.i ); return *this; }
};
int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description laplacianoptions( "Laplacian options" );
	laplacianoptions.add_options()
        ( "n", po::value<int>()->default_value( 1.0 ), "coeff" )
        ( "wr", po::value<bool>()->default_value( true ), "coeff" )
        ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
		;

	Environment env( _argc=argc, _argv=argv,
                     _desc=laplacianoptions,
                     _about=about(_name="qs_laplacian",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));


    using f_t = Mesh<Simplex<FEELPP_DIM,1>>::face_type;
    f_t f,f1;
    f = std::move(f1);

    tic();
    std::vector<f_t2> vref;
    if ( boption( "wr" ) )
        vref.reserve( (int)doption( "mu" ) );
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t2 f(i);
        vref.push_back( std::move(f) );
    }
    using namespace std::string_literals;
    
    toc( "inserting "s + std::to_string( (int)doption( "mu" ) ) + " vref" );

    tic();
    std::vector<f_t2> vref2;
    if ( boption( "wr" ) )
        vref2.reserve( (int)doption( "mu" ) );
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        vref2.emplace_back( i );
    }
    using namespace std::string_literals;
    
    toc( "inserting "s + std::to_string( (int)doption( "mu" ) ) + " vref2" );
    
    
    tic();
    std::unordered_map<int, f_t2> fsetref;
    if ( boption( "wr" ) )
        fsetref.reserve( (int)doption( "mu" ) );
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t2 f(i);
        fsetref.insert( std::make_pair( i, f ) );
    }
    using namespace std::string_literals;
    
    toc( "inserting "s + std::to_string( (int)doption( "mu" ) ) + " f_t2" );

    tic();
    std::unordered_map<int, f_t2> fsetref2;
    if ( boption( "wr" ) )
        fsetref2.reserve( (int)doption( "mu" ) );
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t2 f(i);
        fsetref.emplace( i, std::move(f) );
    }
    using namespace std::string_literals;
    
    toc( "emplacing 1: "s + std::to_string( (int)doption( "mu" ) ) + " f_t2" );

    tic();
    std::unordered_map<int, f_t2> fsetref3;
    if ( boption( "wr" ) )
        fsetref3.reserve( (int)doption( "mu" ) );
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        fsetref.emplace( i, i );
    }
    using namespace std::string_literals;
    
    toc( "emplacing 2: "s + std::to_string( (int)doption( "mu" ) ) + " f_t2" );

    tic();
    std::vector<f_t> vfset;
    if ( boption( "wr" ) )
        vfset.reserve( (int)doption( "mu" ) );
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t f;
        f.setId( i );
        vfset.push_back( std::move(f) );
    }
    using namespace std::string_literals;
    
    toc( "move inserting "s + std::to_string( (int)doption( "mu" ) ) + " faces (vfset)" );

    tic();
    std::vector<f_t> vfset2;
    if ( boption( "wr" ) )
        vfset2.reserve( (int)doption( "mu" ) );
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t f;
        f.setId( i );
        vfset2.emplace_back( std::move(f) );
    }
    using namespace std::string_literals;
    
    toc( "move emplaecing "s + std::to_string( (int)doption( "mu" ) ) + " faces (vfset2)" );
    
    tic();
    std::unordered_map<int, f_t> fset;
    if ( boption( "wr" ) )
        fset.reserve( (int)doption( "mu" ) );
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t f;
        f.setId( i );
        fset.insert( std::make_pair( i, f ) );
    }
    using namespace std::string_literals;
    
    toc( "copy inserting "s + std::to_string( (int)doption( "mu" ) ) + " faces" );

    {
    tic();
    std::unordered_map<int, f_t> fset;
    if ( boption( "wr" ) )
        fset.reserve( (int)doption( "mu" ) );
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t f;
        f.setId( i );
        fset.insert( std::piecewise_construct, std::forward_as_tuple( i ),
                     std::forward_as_tuple( std::move(f) ) );
    }
    using namespace std::string_literals;
    
    toc( "move inserting "s + std::to_string( (int)doption( "mu" ) ) + " faces" );
    }

    {
    tic();
    std::unordered_map<int, f_t> fset;
    if ( boption( "wr" ) )
        fset.reserve( (int)doption( "mu" ) );
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t f;
        f.setId( i );
        fset.emplace( std::piecewise_construct, std::forward_as_tuple( i ),
                      std::forward_as_tuple( i ) );
    }
    using namespace std::string_literals;
    
    toc( "emplace "s + std::to_string( (int)doption( "mu" ) ) + " faces" );
    }
    tic();
    std::unordered_map<int, f_t> fset1;
    if ( boption( "wr" ) )
        fset1.reserve( (int)doption( "mu" ) );
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t f;
        f.setId( i );
        fset1.emplace( i, std::move(f) );
    }
    using namespace std::string_literals;
    
    toc( "emplaceing "s + std::to_string( (int)doption( "mu" ) ) + " faces" );

    tic();
    auto m  = new Mesh<Simplex<FEELPP_DIM,1>>;
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t f;
        f.setId( i );
        m->addFace( m->endFace(), std::move(f) );
    }
    toc( "move insert "s + std::to_string( (int)doption( "mu" ) ) + " faces in mesh" );
    tic();
    auto m1  = new Mesh<Simplex<FEELPP_DIM,1>>;
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t f;
        f.setId( i );
        m1->addFace( m1->endFace(), f );
    }
    toc( "copy insert "s + std::to_string( (int)doption( "mu" ) ) + " faces in mesh" );

    tic();
    auto m2  = new Mesh<Simplex<FEELPP_DIM,1>>;
    for( int i = 0; i < (int)doption( "mu" ); ++i )
    {
        f_t f;
        f.setId( i );
        m2->emplaceFace( m->endFace(), std::move(f) );
    }
    toc( "move emplace 1: "s + std::to_string( (int)doption( "mu" ) ) + " faces in mesh" );
    tic();
    {
        auto m  = new Mesh<Simplex<FEELPP_DIM,1>>;
        for( int i = 0; i < (int)doption( "mu" ); ++i )
        {
            m->emplaceFace( m->endFace(), i );
        }
    }
    toc( "move emplace 2: "s + std::to_string( (int)doption( "mu" ) ) + " faces in mesh" );

}
