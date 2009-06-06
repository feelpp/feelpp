/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-04-05

  Copyright (C) 2005,2006 EPFL

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
   \file test_mixed.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-04-05
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/profiler.hpp>

#include <lifeconfig.h>

#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/debug.hpp>
#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/importer.hpp>
#include <life/lifefilters/exporterEnsight.hpp>
#include <life/lifefilters/exporterGMSH.hpp>
#include <life/lifediscr/geoMap.hpp>
#include <life/lifediscr/quadRule.hpp>
#include <life/lifediscr/refFE.hpp>
#include <life/lifediscr/FEFactory.hpp>
#include <life/lifemesh/RegionTree.hpp>
#include <life/lifediscr/FESpace.hpp>
#include <life/lifediscr/mixedFESpace.hpp>

#include <life/lifealg/SolverPETSC.hpp>
#include <life/lifealg/gmres.hpp>
#include <life/lifealg/preconditioner.hpp>
#include <life/lifealg/solverpardiso.hpp>

#include <life/lifevf/vf.hpp>

using namespace Life;
using namespace Life::vf;

typedef ublas::matrix_range<csr_matrix_type> matrix_range_type;

std::string exporter_str;

struct Simulation
{
    double h;
    double T0;
    double T;
    double dt;
    double nu;
    bool check;

    double gammabc;
    double gammau;
    double gammap;
};

node_type zero( node_type const& __n )
{
    return ublas::scalar_vector<double>( __n.size(), 0 );
}
node_type one( node_type const& __n )
{
    return ublas::unit_vector<double>( __n.size(), 0 );
}
node_type x2y2( node_type const& __n )
{
    node_type __x(__n.size() );
    for( id_type __i = 0; __i < __x.size(); ++__i )
    {
        __x[__i] = __n[__i]*__n[__i];
    }
    return __x;
}
node_type xm2( node_type const& __n )
{
    node_type __x(__n.size() );
    __x[0] = -2+__n[1];
    __x[1] = -2+__n[0];
    return __x;
}
double xy( node_type const& __n )
{
    return __n[0]*__n[1];
}
double xyp( node_type const& __n )
{
    return 2*(__n[0]*__n[1])+__n[0]*__n[1];
}

template<typename Elem>
void test( csr_matrix_type const& m, Elem& __U, Elem& __F,
           Elem& __V, Elem& __S )
{

    __U.element1() = __U.feSpace1()->project( x2y2 );
    __U.element2() = __U.feSpace2()->project( xy );
    __S.element1() = __U.feSpace1()->project( xm2 );
    __S.element2() = __U.feSpace2()->project( xyp );

    ublas::axpy_prod( m, __U, __V );

    std::cerr << "v1  = " << __V.element1() << "\n"
            << "v2 = " << __V.element2() << "\n"
            << "s1  = " << __S.element1() << "\n"
            << "s2 = " << __S.element2() << "\n";

    __V-=__S;
    std::cerr << "error = " << ublas::norm_2( __V ) << "\n";
    std::cerr << "error1 = " << ublas::norm_2( __V.element1() ) << "\n";
    std::cerr << "error2 = " << ublas::norm_2( __V.element2() ) << "\n";

    std::cerr << "error 1 = " << __V.element1() << "\n"
            << "error 2 = " << __V.element2() << "\n";
}

template<typename FEType>
class StokesBlock
{
public:

    typedef FEType fespace_type;
    typedef double value_type;
    typedef typename fespace_type::fespace_1_type velocity_space_type;
    typedef typename fespace_type::fespace_2_type pressure_space_type;
    typedef typename velocity_space_type::element_type velocity_type;
    typedef typename pressure_space_type::element_type pressure_type;

    StokesBlock( fespace_type const& __u )
        :
        _M_A( __u.feSpace1()->nDof(), __u.feSpace1()->nDof() ),
        _M_D( __u.feSpace2()->nDof(), __u.feSpace1()->nDof() ),
        _M_G( __u.feSpace1()->nDof(), __u.feSpace2()->nDof() ),
        _M_P( __u.feSpace2()->nDof(), __u.feSpace2()->nDof() ),
        _M_u( __u.feSpace1()->newElement( "velocity" ) ),
        _M_v( __u.feSpace1()->newElement( "velocity_trial" ) ),
        _M_p( __u.feSpace2()->newElement( "pressure" ) ),
        _M_q( __u.feSpace2()->newElement( "pressure_trial" ) )
        {
            double nu = 1.0;
            std::cerr << "Stokes block A\n";
            //
            // stiffness block
            //
            BilinearForm<velocity_space_type> __vf11( _M_u, _M_v, _M_A );

            LIFE_ASSERT( _M_A.size2() == _M_u.size() )
                        ( _M_A.size2() )( _M_u.size() ).error( "invalid dimensions" );
            LIFE_ASSERT( _M_A.size1() == _M_v.size() )
                        ( _M_A.size1() )( _M_v.size() ).error( "invalid dimensions" );

            __vf11 = ( int2d( "IM_TRIANGLE(4)", nu*dot(grad( _M_u ), grad( _M_v ) ) ) +
                       int1d( "IM_TRIANGLE(4)", -nu*dot( grad( _M_u ), N() )*id( _M_v )
                              + -nu*dot( grad( _M_v ), N() )*id( _M_u ) +
                              1000*nu*id( _M_u )*id( _M_v ) +
                              1000*dot( id( _M_u ), N() )*dot( id( _M_v ), N() ) ) );

            std::cerr << "Stokes block A done\n";
            //std::cout << "A = " << _M_A <<  "\n";
            //
            // gradient block
            //
            std::cerr << "Stokes block G\n";
            BilinearForm<pressure_space_type, velocity_space_type> __vf12( _M_p, _M_v, _M_G );

            LIFE_ASSERT( _M_G.size2() == _M_p.size() )
                        ( _M_G.size2() )( _M_p.size() ).error( "invalid dimensions" );
            LIFE_ASSERT( _M_G.size1() == _M_v.size() )
                        ( _M_G.size1() )( _M_v.size() ).error( "invalid dimensions" );


            __vf12 = ( int2d( "IM_TRIANGLE(4)", -1*id( _M_p )*div( _M_v ) )+
                       int1d( "IM_TRIANGLE(4)", -id( _M_p )*N()*id( _M_v ) ) );

            std::cerr << "Stokes block G done\n";
            //std::cout << "G = " << _M_G <<  "\n";
            //
            // divergence block
            //
            BilinearForm<velocity_space_type, pressure_space_type> __vf21( _M_u, _M_q, _M_D );

            LIFE_ASSERT( _M_G.size1() == _M_q.size() )
                        ( _M_G.size1() )( _M_q.size() ).error( "invalid dimensions" );
            LIFE_ASSERT( _M_G.size2() == _M_u.size() )
                        ( _M_G.size2() )( _M_u.size() ).error( "invalid dimensions" );

            __vf21 = ( int2d( "IM_TRIANGLE(4)", id( _M_q )*div( _M_u ) )+
                       int1d( "IM_TRIANGLE(4)", -id( _M_q )*N()*id( _M_u ) ) );

            //std::cout << "D = " << _M_D <<  "\n";

            std::cout << "D + G^T = " << ublas::norm_frobenius( _M_D + ublas::trans( _M_G ) ) << "\n";
            //
            // pressure/mass block
            //
            BilinearForm<pressure_space_type> __vf22( _M_p, _M_q, _M_P );
            __vf22 = int2d( "IM_TRIANGLE(3)", 1e-6*id( _M_p )*id( _M_q ) );

        }
    template<typename Elem>
    void mult( Elem const& __x, Elem& __y ) const
        {
            __y.element1() = ublas::prod( _M_A, __x.element1() )+ublas::prod( _M_G, __x.element2() );
            __y.element2() = ublas::prod( _M_D, __x.element1() )+ublas::prod( _M_P, __x.element2() );
        }
private:

    csr_matrix_type _M_A;
    csr_matrix_type _M_D;
    csr_matrix_type _M_G;
    csr_matrix_type _M_P;

    velocity_type _M_u;
    velocity_type _M_v;

    pressure_type _M_p;
    pressure_type _M_q;
};
template<typename FEType, typename Elem>
void axpy_prod( StokesBlock<FEType> const& __stokes, Elem const& __x, Elem& __y, bool init = true )
{
    __stokes.mult( __x, __y );
}
template<typename MeshType, typename MatrixType, typename Elem>
void
        solve( std::string const& name, GetPot /*gpdata*/,
       boost::shared_ptr<MeshType> const& __m,
       MatrixType const& __M, Elem & __u, Elem const& __F,
       Elem const& V, Elem const& S )
{
#if 0
    iteration_ptrtype __it = iteration_ptrtype( iteration_type::New() );
    __it->setMaximumNumberOfIterations( 100 );
    __it->setRelativePrecision( 1e-6 );
    __it->setInitialResidual( ublas::norm_2( __F ) );
    PreconditionerIdentity<ublas::vector<double> > prec_id;

    gmres( __M, __u, __F, prec_id, 20, *__it );
#else
    spy(__M,"stokes");
    SolverPardiso sp;
    sp.setMatrix( __M );
    sp.solve( __u.data().begin(), __F.data().begin() );
//     std::cout << "u = " << __u << "\n";
#endif
    //SolverPETSC __solver;
    //__solver.setOptionsFromGetPot( gpdata );
    //__solver.setMaxNonZeroEntriesInRow( P1->maxNonZeroEntriesPerRow() );

    //__solver.setMatribx( __M );
    //__solver.solve( __u, __F );

    //
    // save the results
    //
    typedef exporter<MeshType> export_type;
    typedef typename exporter<MeshType>::timeset_type timeset_type;

    boost::shared_ptr<export_type> __ensight;
    if ( exporter_str == "ensight" )
        __ensight = boost::shared_ptr<export_type>(new exporterEnsight<MeshType>( name ) );
    else if ( exporter_str == "gmsh" )
        __ensight = boost::shared_ptr<export_type>(new exporterGMSH<MeshType>( name ) );
    else
        return;

    typename export_type::timeset_ptrtype __ts( new timeset_type( name ) );
    __ts->setTimeIncrement( 0 );
    __ensight->addTimeSet( __ts );

    typename timeset_type::step_ptrtype __step = __ts->step( 1.0 );
    __step->setMesh( __m  );
    __step->addNodalVector( "velocity", __u.element1().size(), __u.element1().begin(), __u.element1().end()  );
    __step->addNodalVector( "rhs_velocity", __F.element1().size(), __F.element1().begin(), __F.element1().end()  );
    __step->addNodalScalar( "pressure", __u.element2().size(), __u.element2().begin(), __u.element2().end()  );
    __step->addNodalScalar( "rhs_pressure", __F.element2().size(), __F.element2().begin(), __F.element2().end()  );
    __step->addNodalScalar( "error_velocity", V.element1().size(), V.element1().begin(), V.element1().end()  );
    __step->addNodalScalar( "error_pressure", V.element2().size(), V.element2().begin(), V.element2().end()  );
    __step->addNodalScalar( "exact_velocity", S.element1().size(), S.element1().begin(), S.element1().end()  );
    __step->addNodalScalar( "exact_pressure", S.element2().size(), S.element2().begin(), S.element2().end()  );
    __ensight->save();
}
template<typename Elem>
void
getA( csr_matrix_type& __A, Elem& __u, Elem& __v )
{
    BilinearForm<typename Elem::fespace_type> __vf_A( __u, __v, __A );

    __vf_A = int2d( "IM_TRIANGLE(4)", dot(grad( __u ), grad( __v ) ) );
}
template<typename Elem1,typename Elem2>
void
getDt( csr_matrix_type& __A, Elem1& __u, Elem2& __v )
{
    BilinearForm<typename Elem1::fespace_type,typename Elem2::fespace_type> __vf_A( __u, __v, __A );

    __vf_A = int2d( "IM_TRIANGLE(4)", -id(__u)*div(__v));
}
template<typename Elem1,typename Elem2>
void
getD( csr_matrix_type& __A, Elem1& __u, Elem2& __v )
{
    BilinearForm<typename Elem1::fespace_type,typename Elem2::fespace_type> __vf_A( __u, __v, __A );

    __vf_A = int2d( "IM_TRIANGLE(4)", id(__v)*div(__u));
}
template<typename Elem>
void
getF( Elem& __A, Elem& __v )
{
    LinearForm<typename Elem::fespace_type> __vf_A( __v, __A );

    __vf_A = int2d( "IM_TRIANGLE(4)", id(__v) );
}
template<typename Elem>
void
getP( csr_matrix_type& __A, Elem& __u, Elem& __v )
{
    BilinearForm<typename Elem::fespace_type> __vf_A( __u, __v, __A );

    __vf_A = int2d( "IM_TRIANGLE(4)", 1e-6*id ( __u )*id( __v ) );
}

void test_mixed_fv2d( Simulation const& /*simdata*/, GetPot /*gpdata*/ )
{
#if 0
    typedef RegionMesh2D<LinearTriangle> mesh_type;
    //typedef RegionMesh2D<QuadraticTriangle> mesh_type;
    boost::shared_ptr<mesh_type> aMesh( new mesh_type );

    Gmsh __gmsh;
    __gmsh.setOrder( GMSH_ORDER_ONE );
    std::string fname = __gmsh.generateSquare( "cavity", simdata.h );

    importer import( fname, GMSH );
    import.import( *aMesh, 1 );

    typedef MixedFESpace<mesh_type, 2, 1> fespace_type;
    fespace_type::pointer_type P2P1 = fespace_type::New( aMesh, "FEM_PK(2,2)", "FEM_PK(2,1)", "GT_PK(2,1)" );
    fespace_type::element_type U = ( P2P1->newElement( "u_square" ) );
    fespace_type::element_type V = ( P2P1->newElement( "v_square" ) );
    fespace_type::element_type S = ( P2P1->newElement( "s_square" ) );
    fespace_type::element_type __F = ( P2P1->newElement( "f_square" ) );

    fespace_type::element_1_type u = U.element1();
    fespace_type::element_1_type v = U.element1();
    fespace_type::element_2_type p = V.element2();
    fespace_type::element_2_type q = V.element2();

    double nu = 1.0;


    MixedLinearForm<mesh_type, 2, 1> __vf_F( U, __F );
//       __vf_F( v ) = int2d( "IM_TRIANGLE(4)", id( v ) );
//      __vf_F( q ) = int2d( "IM_TRIANGLE(4)", id( q ) );
    //__vf_F( v ) = int1d( "IM_TRIANGLE(4)", 20, ( - ( id( q )*N() - nu*dot( grad( v ), N() ) )*unit( 1 )
    //+ 1000*nu*unit( 1 )*id( v ) ) );
    std::cerr << "F = " << __F  << "\n";
    std::cerr << "1^T F = " << ublas::inner_prod( ScalarVector( v.size(), 1 ), __F.element1() ) << "\n";
    std::cerr << "1^T F = " << ublas::inner_prod( ScalarVector( __F.size(), 1 ), __F ) << "\n";
#if 1
    typedef FESpace<mesh_type,2> P1_type;
    P1_type::pointer_type P1 = P1_type::New(aMesh,"FEM_PK(2,2)", "GT_PK(2,1)" );
    P1_type::element_type u_p1 = ( P1->newElement( "u_p1_square" ) );
    P1_type::element_type v_p1 = ( P1->newElement( "v_p1_square" ) );
    typedef FESpace<mesh_type,1> P1p_type;
    P1p_type::pointer_type P1p = P1p_type::New(aMesh,"FEM_PK(2,1)", "GT_PK(2,1)" );
    P1p_type::element_type p_p1 = ( P1p->newElement( "p_p1_square" ) );
    P1p_type::element_type q_p1 = ( P1p->newElement( "q_p1_square" ) );
    csr_matrix_type __A,__P, __Dt, __D;
    getA( __A, u_p1, v_p1 );
    getDt( __Dt, p_p1, v_p1 );
    getD( __D, u_p1, q_p1 );
    getP( __P, q_p1, q_p1 );

    getF( u_p1, v_p1);
    typedef ublas::vector_range<ublas::vector<double> > vector_range_type;
    vector_range_type vr( __F, ublas::range(0,u.feSpace()->nDof() ));
    std::cerr << "error F = " << ublas::norm_2( vr - u_p1 ) << "\n";

    csr_matrix_type __M1;
    MixedBilinearForm<fespace_type> __vf_M1( U, V, __M1 );


#define TEST_P 1
#define TEST_A 1
#define TEST_Dt 1
#define TEST_D 1

    V = __F;
#if defined(TEST_A)
    __vf_M1( u, v ) = int2d( "IM_TRIANGLE(4)", nu*dot(grad( u ), grad( v ) ) );

    spy(__M1,"after_A");
    std::cerr << "******************************************************8\n";

    std::cerr << "error A : " << ublas::norm_frobenius( matrix_range_type( __M1,
                                  ublas::range( 0, u.feSpace()->nDof() ),
                                  ublas::range( 0, v.feSpace()->nDof() ) ) - __A ) << "\n";
    std::cerr << "******************************************************8\n";
#endif
            //+ on( 20, u, __F, ::one )+on( 10, u, __F, ::zero );
#if defined(TEST_Dt)
     __vf_M1( p, v ) = int2d( "IM_TRIANGLE(4)", -id( p )*div( v ) );

     std::cerr << "******************************************************8\n";
     std::cerr << "error Dt : " << ublas::norm_frobenius( matrix_range_type( __M1,
             ublas::range( 0, u.feSpace()->nDof()),
             ublas::range( v.feSpace()->nDof(), u.feSpace()->nDof()+q.feSpace()->nDof() ) ) - __Dt ) << "\n";
     std::cerr << "******************************************************8\n";
#endif

#if defined(TEST_D)
     __vf_M1( u, q ) = int2d( "IM_TRIANGLE(4)", id( q )*div( u ) );//+on( 20, u, __F, ::one )+on( 10, u, __F, ::zero );

     std::cerr << "******************************************************8\n";
     std::cerr << "error D : " << ublas::norm_frobenius( matrix_range_type( __M1,
             ublas::range( u.feSpace()->nDof(), u.feSpace()->nDof()+q.feSpace()->nDof() ),
             ublas::range( 0, u.feSpace()->nDof() ) ) - __D ) << "\n";

//      std::cerr << "D1 = " << matrix_range_type( __M1,
//              ublas::range( u.feSpace()->nDof(), u.feSpace()->nDof()+q.feSpace()->nDof() ),
//              ublas::range( 0, u.feSpace()->nDof() ) ) << "\n";
//      std::cerr << "D2 = " << __D << "\n";
     std::cerr << "******************************************************8\n";
#endif
    spy(__M1,"before_P");
#if defined(TEST_P)
    __vf_M1( p, q ) = int2d( "IM_TRIANGLE(4)", 1e-6*id( p )*id( q ) )+on( 20, u, __F, ::one )+on( 10, u, __F, ::zero);
    spy(__M1,"after_P");
    spy(__P,"P");
    std::cerr << "******************************************************8\n";
    matrix_range_type M1_P( __M1,
                       ublas::range( u.feSpace()->nDof(), u.feSpace()->nDof()+p.feSpace()->nDof() ),
                       ublas::range( v.feSpace()->nDof(), u.feSpace()->nDof()+q.feSpace()->nDof() ) );
    std::cerr << "error P : " << ublas::norm_frobenius( M1_P - __P ) << "\n";
    ublas::vector<double> __vp( __P.size1() );
    __vp = ublas::scalar_vector<double>(__vp.size(),1);
    std::cerr << "1 P 1 = " << ublas::inner_prod( __vp,ublas::prod( __P,__vp ) ) << "\n";
    std::cerr << "1 MP 1 = " << ublas::inner_prod( __vp,ublas::prod( M1_P,__vp ) ) << "\n";
//     std::cerr << "P1 = " << M1_P << "\n";
//     std::cerr << "P2 = " << __P << "\n";

    std::cerr << "******************************************************8\n";
#endif

    std::cerr << "F = " << __F  << "\n";
    std::cerr << "F1 = " << __F.element1()  << "\n";
    std::cerr << "F2 = " << __F.element2()  << "\n";
    if ( simdata.h >= 0.9 )
        std::cerr << "M1 = " << __M1 << "\n";

    U = ScalarVector( U.size(), 1 );
    std::cerr << "1^T M1 1 = " <<  ublas::inner_prod( U,ublas::prod( __M1, U ) ) << "\n";
    solve( "cavity", gpdata, aMesh, __M1, U, __F, V, S );
#else
#if 0
    StokesBlock<fespace_type> Stokes( *P2P1 );
    solve( "cavity", gpdata, aMesh, Stokes, U, __F );
#else
    csr_matrix_type __M1;
    MixedBilinearForm<fespace_type> __vf_M1( U, V, __M1 );
      __vf_M1( u, v ) = int2d( "IM_TRIANGLE(4)", nu*dot(grad( u ), grad( v ) ) );
      __vf_M1( p, v ) = int2d( "IM_TRIANGLE(4)", -id( p )*div( v ) );
      __vf_M1( u, q ) = int2d( "IM_TRIANGLE(4)", id( q )*div( u ) );
      __vf_M1( p, q ) = int2d( "IM_TRIANGLE(3)", id( p )*id( q ) );

    test( __M1, U, __F, V, S );
    if ( __h >= 0.9 )
        std::cerr << "M1 = " << __M1 << "\n";

    U = ScalarVector( U.size(), 1 );
    std::cerr << "1^T M1 1 = " << ublas::inner_prod( U,ublas::prod( __M1, U ) ) << "\n";
    solve( "cavity", gpdata, aMesh, __M1, U, __F, V, S );
#endif
#endif
#endif
}

#if 0
void test_mixed_fv3d( double __h, GetPot gpdata )
{
    typedef RegionMesh3D<LinearTetra> mesh_type;
    //typedef RegionMesh3D<QuadraticTetra> mesh_type;
    boost::shared_ptr<mesh_type> aMesh( new mesh_type );

    Gmsh __gmsh;
    __gmsh.setOrder( GMSH_ORDER_ONE );
    std::string fname = __gmsh.generateCube( "cavity3d", __h );

    importer import( fname, GMSH );
    import.import( *aMesh, 1 );

    typedef MixedFESpace<mesh_type, 3, 1> fespace_type;
    fespace_type::pointer_type P2P1 = fespace_type::New( aMesh, "FEM_PK(3,2)", "FEM_PK(3,1)", "GT_PK(3,1)" );
    fespace_type::element_type U = ( P2P1->newElement( "u_square" ) );
    fespace_type::element_type V = ( P2P1->newElement( "v_square" ) );
    fespace_type::element_type __F = ( P2P1->newElement( "f_square" ) );

    fespace_type::element_1_type u = U.element1();
    fespace_type::element_1_type v = U.element1();
    fespace_type::element_2_type p = V.element2();
    fespace_type::element_2_type q = V.element2();

    double nu = 1.0;


    MixedLinearForm<mesh_type, 3, 1> __vf_F( U, __F );
//     __vf_F( v ) = int2d( "IM_TRIANGLE(4)", id( v ) );
    //__vf_F( v ) = int1d( "IM_TRIANGLE(4)", 20, ( - ( id( q )*N() - nu*dot( grad( v ), N() ) )*unit( 1 )
    //+ 1000*nu*unit( 1 )*id( v ) ) );
//     std::cerr << "F = " << __F  << "\n";
    std::cerr << "1^T F = " << ublas::inner_prod( ScalarVector( v.size(), 1 ), __F.element1() ) << "\n";

    boost::timer __timer;
    csr_matrix_type __M1;
    MixedBilinearForm<fespace_type> __vf_M1( U, V, __M1 );
    __vf_M1( u, v ) = int3d( "IM_TETRAHEDRON(15)", nu*dot(grad( u ), grad( v ) ) );
    std::cerr << "C has been built in " << __timer.elapsed() << "\n";
    __timer.restart();

    __vf_M1( p, v ) = int3d( "IM_TETRAHEDRON(15)", -id( p )*div( v ) );
    std::cerr << "Dt has been built in " << __timer.elapsed() << "\n";
    __timer.restart();

    __vf_M1( u, q ) = int3d( "IM_TETRAHEDRON(15)", id( q )*div( u ) );
    std::cerr << "D has been built in " << __timer.elapsed() << "\n";
    __timer.restart();

    __vf_M1( p, q ) = int3d( "IM_TETRAHEDRON(4)", 1e-6*id( p )*id( q ) )+
            on( 20, u, __F, ::one )+on( 10, u, __F, ::zero );
    std::cerr << "P has been built in " << __timer.elapsed() << "\n";

//     std::cerr << "F = " << __F  << "\n";
//     std::cerr << "F1 = " << __F.element1()  << "\n";
//     std::cerr << "F2 = " << __F.element2()  << "\n";
    if ( __h >= 0.9 )
        std::cerr << "M1 = " << __M1 << "\n";

    __timer.restart();
    U = ScalarVector( U.size(), 1 );
    std::cerr << "1^T M1 1 = " << ublas::inner_prod( U,ublas::prod( __M1, U ) ) << "\n";
    solve( "cavity3d", gpdata, aMesh, __M1, U, __F, V, V );
    std::cerr << "system solved in " << __timer.elapsed() << "\n";
}
#endif
void test_mixed_ns2d( Simulation const& simdata, GetPot gpdata )
{
    using namespace Life::vf;

    typedef RegionMesh2D<LinearTriangle> mesh_type;
    //typedef RegionMesh2D<QuadraticTriangle> mesh_type;
    boost::shared_ptr<mesh_type> aMesh( new mesh_type );

    Gmsh __gmsh;
    __gmsh.setOrder( GMSH_ORDER_ONE );
    std::string fname = __gmsh.generateSquare( "cavity2d", simdata.h );

    importer import( fname, GMSH );
    import.import( *aMesh, 1 );

//     const std::string QDR=QDR;
//     const int QD_nbp = 4;
    const std::string QDR="IM_TRIANGLE(3)";
    const int QD_nbp = 3;
    typedef MixedFESpace<mesh_type, VECTORIAL, SCALAR> fespace_type;
    fespace_type::pointer_type P2P1 = fespace_type::New( aMesh, "FEM_PK(2,1)", "FEM_PK(2,1)", "GT_PK(2,1)" );
    fespace_type::element_type U = ( P2P1->newElement( "u_cavity2d" ) );
    fespace_type::element_type V = ( P2P1->newElement( "v_cavity2d" ) );
    ublas::vector<double> __F( U.size() );

    fespace_type::element_1_type u = U.element1();
    fespace_type::element_1_type v = U.element1();
    fespace_type::element_2_type p = V.element2();
    fespace_type::element_2_type q = V.element2();



    MixedLinearForm<fespace_type> __vf_F( U, __F );
//     __vf_F( v ) = int2d( QDR, (1./simdata.dt)*id( v ) );
    //__vf_F( v ) = int1d( QDR, 20, ( - ( id( q )*N() - nu*dot( grad( v ), N() ) )*unit( 1 )
    //+ 1000*nu*unit( 1 )*id( v ) ) );
//     std::cerr << "F = " << __F  << "\n";
    std::cerr << "1^T F = " << ublas::inner_prod( ScalarVector( v.size(), 1 ),
                                                  ublas::project(__F,ublas::range(0,U.element1().size())) ) << "\n";


    csr_matrix_type __NS;
    MixedBilinearForm<fespace_type> __vf_NS( U, V, __NS );
    //
    // save the results
    //
    typedef exporter<mesh_type> export_type;
    typedef exporter<mesh_type>::timeset_type timeset_type;

    boost::shared_ptr<export_type> __ensight;
    if ( exporter_str == "ensight" )
        __ensight = boost::shared_ptr<export_type>(new exporterEnsight<mesh_type>( "cavity2d" ) );
    else if ( exporter_str == "gmsh" )
        __ensight = boost::shared_ptr<export_type>(new exporterGMSH<mesh_type>( "cavity2d") );
    else
        return;

    export_type::timeset_ptrtype __ts( new timeset_type( "cavity2d" ) );
    __ts->setTimeIncrement( simdata.dt );
    __ensight->addTimeSet( __ts );

    v.feSpace()->gm()->setCacheInformation( QDR, aMesh->numElements(), QD_nbp);
    p.feSpace()->gm()->setCacheInformation( QDR, aMesh->numElements(), QD_nbp);
    v.feSpace()->fe()->setCacheInformation( QDR, aMesh->numElements(), QD_nbp);
    p.feSpace()->fe()->setCacheInformation( QDR, aMesh->numElements(), QD_nbp);

    fespace_type::element_1_type::component_type ux = u.comp(X);
    fespace_type::element_1_type::component_type uy = u.comp(Y);
    fespace_type::element_1_type::component_type vx = v.comp(X);
    fespace_type::element_1_type::component_type vy = v.comp(Y);

    U.element1() = ublas::scalar_vector<double>(U.element1().size(),1);

    __vf_F( v ) = int2d( QDR, vfprod(vec(U.element1()),id( v ) ) );

    std::cout << "1^T F = " << ublas::inner_prod( ScalarVector( v.size(), 1 ),
                                                  ublas::project(__F,ublas::range(0,U.element1().size())));

    if ( simdata.T0 == 0 )
        U.element1() = ublas::scalar_vector<double>(U.element1().size(),0);
    else
    {
        // reload previous time steps
        // it should depend on the time step
        timeset_type::step_ptrtype __step = __ts->step( simdata.T0 );
        __step->load();
        __step->setMesh( aMesh );
        U.element1() = __step->nodalVector( "velocity" );
        U.element2() = __step->nodalScalar( "pressure" );
    }
/*
    +  int1d(Th)((p*N.x  - dn(u1))*v1+(p*N.y - dn(u2))*v2)
  +  int1d(Th)((-q*N.x  - dn(v1))*u1+(-q*N.y - dn(v2))*u2)
  +  int1d(Th)(gammabc/lenEdge*(u1*v1+u2*v2))
  +  int1d(Th)(gammabc/lenEdge*(u1*N.x+u2*N.y)*(v1*N.x+v2*N.y))
//  Right hand side imposed boundary conditions
  -  int1d(Th,3)((-q*N.x  - dn(v1))*1.e0)
  -  int1d(Th,3)(gammabc/lenEdge*1.0*v1)
  -  int1d(Th,3)(gammabc/lenEdge*(1.0*N.x)*v1*N.x)
 */
    bool build_pattern = true;
    for( double t = simdata.T0+simdata.dt; t <= simdata.T ; t += simdata.dt )
    {
        std::cout << "************************************************************\n"
                  << "Time = " << t << " || time step = " << simdata.dt << " || Final Time = " << simdata.T << "\n";


        boost::timer __timer;
        size_type pf = PATTERN_FULL;
        size_type pd = PATTERN_DIAGONAL;

         __vf_F( v ) = int2d( QDR, vfprod(vec(U.element1()),id( v ) )/simdata.dt );
#if 0
        __vf_F( vx ) = int1d( QDR, 20, simdata.gammabc*(simdata.nu/hFace())*id(vx)*1.0-
                                       simdata.nu*dn(vx)*1.0+
                                       simdata.gammabc*id(vx)*Nx()*1.0*Nx() );
         __vf_F( vy ) = int1d( QDR, 20, simdata.gammabc*id(vy)*Ny()*1.0*Nx());
         __vf_F( q ) = int1d( QDR, 20, -1.0*id( q )*Nx() );
#else
         __vf_F( q ) = int1d( QDR, 20, 0*id(q) );
#endif



        std::cout << "rhs has been built in " << __timer.elapsed() << "\n";
        __timer.restart();


        fespace_type::element_1_type::component_type Ux = U.element1().comp(X);
        fespace_type::element_1_type::component_type Uy = U.element1().comp(Y);

        __vf_NS( u, v, build_pattern, pd  ) = int2d( QDR,
                                                      id(u)*id(v)/simdata.dt+
                                                      dot(vec(U.element1()),grad(u))*id(v)+
                                                              simdata.nu*dot(grad( u ), grad( v ) ) );
//                 intif( QDR, simdata.gammau*vf::abs(dot(N(),vec(U.element1()))/(hFace()^2.0))*jump(dn(u))*jump(dn(v)));
#if 0
        __vf_NS( ux, vx, false, pd  ) += intif( QDR, simdata.gammau*((((id(U.feSpace1(),Ux)^2.0)+(id(U.feSpace1(),Uy)^2.0))^0.5)*(hFace()^2.0))*jump(dx(ux))*jump(dx(vx)));
        __vf_NS( uy, vy, false, pd  ) += intif( QDR, simdata.gammau*((((id(U.feSpace1(),Ux)^2.0)+(id(U.feSpace1(),Uy)^2.0))^0.5)*(hFace()^2.0))*jump(dy(uy))*jump(dy(vy)));
        __vf_NS( ux, vy, build_pattern, pd  ) += intif( QDR, simdata.gammau*((((id(U.feSpace1(),Ux)^2.0)+(id(U.feSpace1(),Uy)^2.0))^0.5)*(hFace()^2.0))*jump(dx(ux))*jump(dy(vy)));
        __vf_NS( uy, vx, build_pattern, pd  ) += intif( QDR,
 simdata.gammau*((((id(U.feSpace1(),Ux)^2.0)+(id(U.feSpace1(),Uy)^2.0))^0.5)*(hFace()^2.0))*jump(dy(uy))*jump(dx(vx)));
#endif
#if 0
                int1d( QDR, 10, -simdata.nu*(dn(u)*id(v)+dn(v)*id(u))+
                            simdata.gammabc*(simdata.nu/hFace())*id(u)*id(v)+
                                    simdata.gammabc*vfprod(N(),id(u))*vfprod(N(),id(v)))+
                int1d( QDR, 20, -simdata.nu*(dn(u)*id(v)+dn(v)*id(u))+
                                    simdata.gammabc*(simdata.nu/hFace())*id(u)*id(v)+
                                            simdata.gammabc*vfprod(N(),id(u))*vfprod(N(),id(v)));
#endif
        std::cout << "C has been built in " << __timer.elapsed() << "\n";
        __timer.restart();

        __vf_NS( p, v, build_pattern, pf  ) = int2d( QDR, -1*id( p )*div( v ) );
//                 int1d(QDR, 10, vfprod(N(),id(p))*id(v) )+
//                 int1d(QDR, 20, vfprod(N(),id(p))*id(v) );
        // weak dirichlet condition: \int p.N v
//         __vf_NS( p, vx, false  ) += int1d(QDR, id(p)*Nx()*id(vx) );
//         __vf_NS( p, vy, false  ) += int1d(QDR, id(p)*Ny()*id(vy) );

        std::cout << "Dt has been built in " << __timer.elapsed() << "\n";
        __timer.restart();

        __vf_NS( u, q, build_pattern, pf  ) = int2d( QDR, id( q )*div( u ) );
//                 int1d(QDR, 10, -1.0*vfprod(N(),id(q))*id(u) )+
//                 int1d(QDR, 20, -1.0*vfprod(N(),id(q))*id(u) );
        // weak dirichlet condition: \int - q.N u
//         __vf_NS( ux, q, false ) += int1d(QDR, -id(q)*Nx()*id(ux) );
//         __vf_NS( uy, q, false ) += int1d(QDR, -id(q)*Ny()*id(uy) );

        std::cout << "D has been built in " << __timer.elapsed() << "\n";
        __timer.restart();

        __vf_NS( p, q, build_pattern, pd  ) = int2d( QDR, 1e-6*id( p )*id( q ) )+
                                              intif(  QDR,
                                                      simdata.gammap*(hFace()^3)*jump(dn(p))*jump(dn(q))/simdata.nu )+
                on( 20, u, __F, oneX() )+on( 10, u, __F, 0 );
        std::cout << "P has been built in " << __timer.elapsed() << "\n";

        spy(__NS,"pattern_cavity2d");
        // from now on don't build pattern anymore
        build_pattern = false;
        __timer.restart();

        SolverPardiso sp;
        sp.setMatrix( __NS );
        sp.solve( U.data().begin(), __F.data().begin() );
        std::cout << "Problem  has been solved in " << __timer.elapsed() << "\n";

        timeset_type::step_ptrtype __step = __ts->step( t );
        __step->setMesh( aMesh );
        __step->addNodalVector( "velocity", U.element1().size(), U.element1().begin(), U.element1().end()  );
//         __step->addNodalVector( "rhs_velocity", __F.element1().size(), __F.element1().begin(), __F.element1().end()  );
        __step->addNodalScalar( "pressure", U.element2().size(), U.element2().begin(), U.element2().end()  );
//         __step->addNodalScalar( "rhs_pressure", __F.element2().size(), __F.element2().begin(), __F.element2().end()  );
        __ensight->save();

        std::cout << "||velocity||_inf = " << ublas::norm_inf( U.element1() ) << "\n";
        if ( simdata.check && ublas::norm_inf( U.element1() ) > 1 )
        {
            break;
        }
    }
}
int
main( int argc,  char** argv )
{
    Life::Assert::setLog( "assertions.log");

    std::string data_file_name;

    Simulation S;
    bool __stokes = false;
    bool __navier_stokes = false;

    po::options_description desc("Specific options");
    desc.add_options()
            ("file,f", Life::po::value<std::string>(&data_file_name)->default_value( "data" ), "data file name")
            ("export,x", Life::po::value<std::string>(&exporter_str)->default_value( "ensight" ), "export type")
            ("stokes,s", Life::po::value<bool>(&__stokes)->default_value( false ), "stokes solve")
            ("ns,n", Life::po::value<bool>(&__navier_stokes)->default_value( true ), "navier stokes solve")
            ("h", Life::po::value<double>(&S.h)->default_value( 0.1 ), "h value")
            ("dt", Life::po::value<double>(&S.dt)->default_value( 0.1 ), "time step value")
            ("T0", Life::po::value<double>(&S.T0)->default_value( 0 ), "start time value")
            ("Tf", Life::po::value<double>(&S.T)->default_value( 1 ), "final time value")
            ("nu", Life::po::value<double>(&S.nu)->default_value( 1 ), "viscosity value")
            ("gammabc", Life::po::value<double>(&S.gammabc)->default_value( 10 ), "weak dirichlet coef")
            // stabilization
            ("gammau", Life::po::value<double>(&S.gammau)->default_value( 2.5e-2 ), "u stabilization coef")
            ("gammap", Life::po::value<double>(&S.gammap)->default_value( 2.5e-2 ), "p stabilization coef")
            ("check,c", Life::po::value<bool>(&S.check)->default_value( false ), "check that ||u||_inf <= 1")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    GetPot data_file(data_file_name.c_str());

    if ( __stokes )
        test_mixed_fv2d( S, data_file );
    if ( __navier_stokes )
        test_mixed_ns2d( S, data_file );
}
