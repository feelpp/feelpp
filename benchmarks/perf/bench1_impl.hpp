/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-06-14

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file bench1_impl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-06-14
 */
#include <bench1.hpp>
#include <boost/algorithm/string/replace.hpp>

#define DO_BENCH( TheExpr, str1, str  )                             \
    {                                                               \
        const int Qorder = ExpressionOrder<decltype(elements(Xh->mesh())),decltype(TheExpr)>::value; \
    typename _Q<Qorder>::template apply<FSType::fe_type::nDim,double,Simplex>::type quad; \
    timer.restart();                                                    \
    form2(Xh,Xh,_matrix=M) += integrate( _range=elements(Xh->mesh()),  _expr=TheExpr );      \
    ofs << std::setw(20) << std::string(str1)+"_"+::boost::algorithm::replace_all_copy( std::string(str), " ", "_" ) << " " \
        << std::setw(4) << FSType::fe_type::nOrder << " "               \
        << std::setw(8) << std::setprecision( 5 ) << timer.elapsed() << " " \
        << std::setw(5) << Xh->mesh()->numElements() << " "             \
        << std::setw(5) << Xh->nLocalDof() << " "                       \
        << std::setw(5) << Qorder << " "                                \
        << std::setw(5) << quad.nPoints() << std::endl;                 \
    }
/**/

#define DO_BENCH2( TheExpr, TheExprQ, str1, str  )                  \
    {                                                               \
        const int Qorder = ExpressionOrder<decltype(elements(Xh->mesh())),decltype(TheExprQ)>::value; \
    typename _Q<Qorder>::template apply<FSType::fe_type::nDim,double,Simplex>::type quad; \
    timer.restart();                                                \
    form2(Xh,Xh,_matrix=M) += integrate( _range=elements(Xh->mesh()), _expr=TheExpr, _quad=_Q<Qorder>() );      \
    ofs << std::setw(20) << std::string(str1)+"_"+::boost::algorithm::replace_all_copy( std::string(str), " ", "_" ) << " " \
        << std::setw(4) << FSType::fe_type::nOrder << " "               \
        << std::setw(8) << std::setprecision( 5 ) << timer.elapsed() << " " \
        << std::setw(5) << Xh->mesh()->numElements() << " "             \
        << std::setw(5) << Xh->nLocalDof() << " "                       \
        << std::setw(5) << Qorder << " "                                \
        << std::setw(5) << quad.nPoints() << std::endl;                 \
    }
/**/

namespace Feel
{

template<typename FSType>
void
Bench1::R( boost::shared_ptr<FSType> const& Xh  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );


    typedef fusion::vector<Lagrange<0, Scalar> > dp0_basis_type;
    typedef FunctionSpace<typename FSType::mesh_type, dp0_basis_type, Discontinuous> dp0_space_type;
    typename dp0_space_type::pointer_type P0h = dp0_space_type::New( Xh->mesh() );
    typename dp0_space_type::element_type w( P0h );
    w = vf::project( P0h, elements( Xh->mesh() ), Px() );
    form2( Xh,Xh,_matrix=M,_init=true );
    boost::timer timer;

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    std::ofstream ofs( ( boost::format( "perf-R-%1%-%2%.dat" ) % FSType::fe_type::nDim % FSType::fe_type::nOrder ).str().c_str(), std::ios::app );
    ofs.precision( 4 );
    ofs.width( 6 );
    DO_BENCH( idt( u )*id( v ), "R", "const" );
    DO_BENCH( idv( w )*idt( u )*id( v ), "R", "p0" );
    DO_BENCH( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*idt( u )*id( v ), "R", "xyz" );
    DO_BENCH2( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*idt( u )*id( v ),
               idt( u )*id( v ),
               "R", "xyz(const)" );
    DO_BENCH( idv( w )*( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*idt( u )*id( v ), "R", "p0 xyz" );
    DO_BENCH2( idv( w )*( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*idt( u )*id( v ),
               idt( u )*id( v ),
               "R", "p0 xyz(const)" );
    DO_BENCH( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*idt( u )*id( v ), "R", "val(xyz)" );
    DO_BENCH2( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*idt( u )*id( v ),
               idt( u )*id( v ),
               "R", "val(xyz)(const)" );
    DO_BENCH( idv( w )*val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*idt( u )*id( v ),
              "R", "p0 val(xyz)" );
    DO_BENCH2( idv( w )*val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*idt( u )*id( v ),
               idt( u )*id( v ),
               "R", "p0 val(xyz)(const)" );

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}

template<typename FSType>
void
Bench1::D( boost::shared_ptr<FSType> const& Xh  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    form2( Xh,Xh,_matrix=M,_init=true );
    boost::timer timer;

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    std::ofstream ofs( ( boost::format( "perf-D-%1%-%2%.dat" ) % FSType::fe_type::nDim % FSType::fe_type::nOrder ).str().c_str(), std::ios::app );
    ofs.precision( 4 );
    ofs.width( 6 );
    DO_BENCH( gradt( u )*trans( grad( v ) ), "D", "const" );
    DO_BENCH2( gradt( u )*trans( grad( v ) ), idt( u )*id( v ), "D", "const (R)" );
    DO_BENCH( dxt( u )*dx( v )+dyt( u )*dy( v ), "D", "const2" );
    DO_BENCH( dxt( u )*dx( v )+dyt( u )*dy( v )+dzt( u )*dz( v ), "D", "const3" );
    DO_BENCH( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*gradt( u )*trans( grad( v ) ), "D", "xyz" );
    DO_BENCH2( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*gradt( u )*trans( grad( v ) ),
               gradt( u )*trans( grad( v ) ),"D", "xyz(const)" );
    DO_BENCH( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*gradt( u )*trans( grad( v ) ), "D", "val(xyz)" );
    DO_BENCH2( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*gradt( u )*trans( grad( v ) ),
               gradt( u )*trans( grad( v ) ),
               "D", "val(xyz)(const)" );

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}
template<typename FSType>
void
Bench1::A( boost::shared_ptr<FSType> const& Xh, mpl::int_<1>  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    typedef FunctionSpace<typename FSType::mesh_type, bases<Lagrange<0,Vectorial,Discontinuous> > > P0h_type;
    boost::shared_ptr<P0h_type> P0h = P0h_type::New( Xh->mesh() );
    auto beta0 = P0h->element( "beta" );

    typedef FunctionSpace<typename FSType::mesh_type, bases<Lagrange<FSType::fe_type::nOrder,Vectorial> > > Pkh_type;
    boost::shared_ptr<Pkh_type> Pkh = Pkh_type::New( Xh->mesh() );
    auto betak = Pkh->element( "betak" );

    form2( Xh,Xh,_matrix=M,_init=true );
    boost::timer timer;

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    std::ofstream ofs( ( boost::format( "perf-A-%1%-%2%.dat" ) % FSType::fe_type::nDim % FSType::fe_type::nOrder ).str().c_str(), std::ios::app );
    ofs.precision( 4 );
    ofs.width( 6 );

    DO_BENCH( ( gradt( u )*vec( constant( 1.0 ) ) )*id( v ), "A", "const" );
    DO_BENCH2( ( gradt( u )*vec( constant( 1.0 ) ) )*id( v ),idt( u )*id( v ),"A", "const (R)" );
    DO_BENCH( ( gradt( u )*vec( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ) ) )*id( v ), "A", "xyz" );
    DO_BENCH( ( gradt( u )*vec( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ) ) )*id( v ), "A", "val(xyz)" );
    DO_BENCH2( ( gradt( u )*vec( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ) ) )*id( v ),
               ( gradt( u )*vec( constant( 1.0 ) ) )*id( v ),
               "A", "xyz(const)" );
    DO_BENCH2( ( gradt( u )*vec( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ) ) )*id( v ),
               ( gradt( u )*vec( constant( 1.0 ) ) )*id( v ),
               "A", "val(xyz)(const)" );


#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}

template<typename FSType>
void
Bench1::A( boost::shared_ptr<FSType> const& Xh, mpl::int_<2>  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    typedef FunctionSpace<typename FSType::mesh_type, bases<Lagrange<0,Vectorial,Discontinuous> > > P0h_type;
    boost::shared_ptr<P0h_type> P0h = P0h_type::New( Xh->mesh() );
    auto beta0 = P0h->element( "beta" );

    typedef FunctionSpace<typename FSType::mesh_type, bases<Lagrange<FSType::fe_type::nOrder,Vectorial> > > Pkh_type;
    boost::shared_ptr<Pkh_type> Pkh = Pkh_type::New( Xh->mesh() );
    auto betak = Pkh->element( "betak" );

    form2( Xh,Xh,_matrix=M,_init=true );
    boost::timer timer;

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    std::ofstream ofs( ( boost::format( "perf-A-%1%-%2%.dat" ) % FSType::fe_type::nDim % FSType::fe_type::nOrder ).str().c_str(), std::ios::app );
    ofs.precision( 4 );
    ofs.width( 6 );

    DO_BENCH( ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ) ) )*id( v ), "A", "const" );
    DO_BENCH2( ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ) ) )*id( v ),idt( u )*id( v ),"A", "const (R)" );
    DO_BENCH( ( gradt( u )*vec( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ),( ( Px()*Px()*Px() )+( Py()*Py() ) ) ) )*id( v ), "A", "xyz" );
    DO_BENCH( ( gradt( u )*vec( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ),val( ( Px()*Px()*Px() )+( Py()*Py() ) ) ) )*id( v ), "A", "val(xyz)" );
    DO_BENCH2( ( gradt( u )*vec( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ),( ( Px()*Px()*Px() )+( Py()*Py() ) ) ) )*id( v ),
               ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ) ) )*id( v ),
               "A", "xyz(const)" );
    DO_BENCH2( ( gradt( u )*vec( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ),val( ( Px()*Px()*Px() )+( Py()*Py() ) ) ) )*id( v ),
               ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ) ) )*id( v ),
               "A", "val(xyz)(const)" );

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}

template<typename FSType>
void
Bench1::A( boost::shared_ptr<FSType> const& Xh, mpl::int_<3>  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    typedef FunctionSpace<typename FSType::mesh_type, bases<Lagrange<0,Vectorial,Discontinuous> > > P0h_type;
    boost::shared_ptr<P0h_type> P0h = P0h_type::New( Xh->mesh() );
    auto beta0 = P0h->element( "beta" );

    typedef FunctionSpace<typename FSType::mesh_type, bases<Lagrange<FSType::fe_type::nOrder,Vectorial> > > Pkh_type;
    boost::shared_ptr<Pkh_type> Pkh = Pkh_type::New( Xh->mesh() );
    auto betak = Pkh->element( "betak" );

    form2( Xh,Xh,_matrix=M,_init=true );
    boost::timer timer;

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    std::ofstream ofs( ( boost::format( "perf-A-%1%-%2%.dat" ) % FSType::fe_type::nDim % FSType::fe_type::nOrder ).str().c_str(), std::ios::app );
    ofs.precision( 4 );
    ofs.width( 6 );

    DO_BENCH( ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ),constant( 1.0 ) ) )*id( v ), "A", "const" );
    DO_BENCH2( ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ),constant( 1.0 ) ) )*id( v ),idt( u )*id( v ),"A", "const (R)" );
    DO_BENCH( ( gradt( u )*vec( ( ( ( Px()*Px()*Px() )+( Py()*Py() ) )*Pz() ),( ( Px()*Px()*Px() )+( Py()*Py() ) ),( Px()*Px()*Px() ) ) )*id( v ), "A", "xyz" );
    DO_BENCH( ( gradt( u )*vec( val( ( ( Px()*Px()*Px() )+( Py()*Py() ) )*Pz() ),val( ( Px()*Px()*Px() )+( Py()*Py() ) ),val( Px()*Px()*Px() ) ) )*id( v ), "A", "val(xyz)" );
    DO_BENCH2( ( gradt( u )*vec( ( ( ( Px()*Px()*Px() )+( Py()*Py() ) )*Pz() ),( ( Px()*Px()*Px() )+( Py()*Py() ) ),( Px()*Px()*Px() ) ) )*id( v ),
               ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ),constant( 1.0 ) ) )*id( v ),
               "A", "xyz(const)" );
    DO_BENCH2( ( gradt( u )*vec( val( ( ( Px()*Px()*Px() )+( Py()*Py() ) )*Pz() ),val( ( Px()*Px()*Px() )+( Py()*Py() ) ),val( Px()*Px()*Px() ) ) )*id( v ),
               ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ),constant( 1.0 ) ) )*id( v ),
               "A", "val(xyz)(const)" );

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}

template<typename FSType>
void
Bench1::DR( boost::shared_ptr<FSType> const& Xh )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    form2( Xh,Xh,_matrix=M,_init=true );
    boost::timer timer;

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    std::ofstream ofs( ( boost::format( "perf-DR-%1%-%2%.dat" ) % FSType::fe_type::nDim % FSType::fe_type::nOrder ).str().c_str(), std::ios::app );
    ofs.precision( 4 );
    ofs.width( 6 );
    DO_BENCH( gradt( u )*trans( grad( v ) )+idt( u )*id( v ), "DR", "const" );
    DO_BENCH( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ), "DR", "xyz" );
    DO_BENCH( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ), "DR", "val(xyz)" );
    DO_BENCH2( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ),
               gradt( u )*trans( grad( v ) )+idt( u )*id( v ),
               "DR", "xyz(const)" );
    DO_BENCH2( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ),
               gradt( u )*trans( grad( v ) )+idt( u )*id( v ),
               "DR", "val(xyz)(const)" );

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}
template<typename FSType>
void
Bench1::ADR( boost::shared_ptr<FSType> const& Xh, mpl::int_<1>  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    form2( Xh,Xh,_matrix=M,_init=true );
    boost::timer timer;

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    std::ofstream ofs( ( boost::format( "perf-ADR-%1%-%2%.dat" ) % FSType::fe_type::nDim % FSType::fe_type::nOrder ).str().c_str(), std::ios::app );
    ofs.precision( 4 );
    ofs.width( 6 );
    DO_BENCH( gradt( u )*trans( grad( v ) )+idt( u )*id( v )+( gradt( u )*vec( constant( 1.0 ) ) )*id( v ), "ADR", "const" );
    DO_BENCH( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
              ( gradt( u )*vec( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ) ) )*id( v ), "ADR", "xyz" );
    DO_BENCH( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
              ( gradt( u )*vec( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ) ) )*id( v ), "ADR", "val(xyz)" );
    DO_BENCH2( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
               ( gradt( u )*vec( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ) ) )*id( v ),
               gradt( u )*trans( grad( v ) )+idt( u )*id( v )+( gradt( u )*vec( constant( 1.0 ) ) )*id( v ),
               "ADR", "xyz(const)" );
    DO_BENCH2( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
               ( gradt( u )*vec( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ) ) )*id( v ),
               gradt( u )*trans( grad( v ) )+idt( u )*id( v )+( gradt( u )*vec( constant( 1.0 ) ) )*id( v ),
               "ADR", "val(xyz)(const)" );
#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}

template<typename FSType>
void
Bench1::ADR( boost::shared_ptr<FSType> const& Xh, mpl::int_<2>  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    form2( Xh,Xh,_matrix=M,_init=true );
    boost::timer timer;

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    std::ofstream ofs( ( boost::format( "perf-ADR-%1%-%2%.dat" ) % FSType::fe_type::nDim % FSType::fe_type::nOrder ).str().c_str(), std::ios::app );
    ofs.precision( 4 );
    ofs.width( 6 );
    DO_BENCH( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) +
              ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ) ) )*id( v ), "ADR", "const" );
    DO_BENCH( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
              ( gradt( u )*vec( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ),( ( Px()*Px()*Px() )+( Py()*Py() ) ) ) )*id( v ), "ADR", "xyz" );
    DO_BENCH( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
              ( gradt( u )*vec( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ),val( ( Px()*Px()*Px() )+( Py()*Py() ) ) ) )*id( v ), "ADR", "val(xyz)" );
    DO_BENCH2( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
               ( gradt( u )*vec( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ),( ( Px()*Px()*Px() )+( Py()*Py() ) ) ) )*id( v ),
               gradt( u )*trans( grad( v ) )+idt( u )*id( v ) +
               ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ) ) )*id( v ),
               "ADR", "xyz(const)" );
    DO_BENCH2( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
               ( gradt( u )*vec( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() ),val( ( Px()*Px()*Px() )+( Py()*Py() ) ) ) )*id( v ),
               gradt( u )*trans( grad( v ) )+idt( u )*id( v ) +
               ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ) ) )*id( v ),
               "ADR", "val(xyz)(const)" );

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}
template<typename FSType>
void
Bench1::ADR( boost::shared_ptr<FSType> const& Xh, mpl::int_<3>  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype  M( M_backend->newMatrix( Xh, Xh ) );

    form2( Xh,Xh,_matrix=M,_init=true );
    boost::timer timer;

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStart( "perf" );
#endif
    std::ofstream ofs( ( boost::format( "perf-ADR-%1%-%2%.dat" ) % FSType::fe_type::nDim % FSType::fe_type::nOrder ).str().c_str(), std::ios::app );
    ofs.precision( 4 );
    ofs.width( 6 );
    DO_BENCH( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) +
              ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ),constant( 1.0 ) ) )*id( v ), "ADR", "const" );
    DO_BENCH( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
              ( gradt( u )*vec( ( ( ( Px()*Px()*Px() )+( Py()*Py() ) )*Pz() ),( ( Px()*Px()*Px() )+( Py()*Py() ) ),( Px()*Px()*Px() ) ) )*id( v ), "ADR", "xyz" );
    DO_BENCH( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
              ( gradt( u )*vec( val( ( ( Px()*Px()*Px() )+( Py()*Py() ) )*Pz() ),val( ( Px()*Px()*Px() )+( Py()*Py() ) ),val( Px()*Px()*Px() ) ) )*id( v ), "ADR", "val(xyz)" );
    DO_BENCH2( ( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
               ( gradt( u )*vec( ( ( ( Px()*Px()*Px() )+( Py()*Py() ) )*Pz() ),( ( Px()*Px()*Px() )+( Py()*Py() ) ),( Px()*Px()*Px() ) ) )*id( v ),
               gradt( u )*trans( grad( v ) )+idt( u )*id( v ) +
               ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ),constant( 1.0 ) ) )*id( v ),
               "ADR", "xyz(const)" );
    DO_BENCH2( val( ( Px()*Px()*Px() )+( Py()*Py() )*Pz() )*( gradt( u )*trans( grad( v ) )+idt( u )*id( v ) ) +
               ( gradt( u )*vec( val( ( ( Px()*Px()*Px() )+( Py()*Py() ) )*Pz() ),val( ( Px()*Px()*Px() )+( Py()*Py() ) ),val( Px()*Px()*Px() ) ) )*id( v ),
               gradt( u )*trans( grad( v ) )+idt( u )*id( v ) +
               ( gradt( u )*vec( constant( 1.0 ),constant( 1.0 ),constant( 1.0 ) ) )*id( v ),
               "ADR", "val(xyz)(const)" );

#if defined(FEELPP_HAS_GOOGLE_PROFILER_H)
    ProfilerStop();
#endif
}
template<typename MeshType, int Order>
void
Bench1::bench1( boost::shared_ptr<MeshType> & mesh )
{


    const int nDim = MeshType::nDim;
    std::ostringstream femstr;
    femstr << "FEM_PK(" << nDim << "," << Order << ")";
    std::ostringstream geostr;
    geostr << "GT_PK(" << nDim << "," << 1 << ")";

    LOG(INFO) << "------------------------------------------------------------" << "\n";
    LOG(INFO) << "dimension : " << nDim << "\n";
    LOG(INFO) << "      fem : " << femstr.str() << "\n";
    LOG(INFO) << "      geo : " << geostr.str() << "\n";
    LOG(INFO) << "++++++++++++++++++++++++++++++" << "\n";

    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;
    typedef FunctionSpace<MeshType, basis_type> space_type;
    typename space_type::pointer_type Xh = space_type::New( mesh );



    //IM<nDim, 2*Order, double, Simplex> im;
    //IM<nDim, 2*Order,double> im;
    //IM_PK<nDim, 4> im;
    //IMSimplex<2,2*Order> im;

    //v.space()->gm()->setCacheInformation( QDR, mesh->numElements(), QDR_NPTS);
    //v.space()->fe()->setCacheInformation( QDR, mesh->numElements(), QDR_NPTS);

    LOG(INFO) << "dof : " << Xh->nDof() << "\n"
          << "elt : " << Xh->mesh()->numElements() << "\n";


    boost::timer timer;

    R( Xh );
    D( Xh );
    A( Xh, mpl::int_<nDim>() );
    //DR( Xh );
    //ADR( Xh, mpl::int_<nDim>() );
    LOG(INFO) << "------------------------------------------------------------" << "\n";

}



}
