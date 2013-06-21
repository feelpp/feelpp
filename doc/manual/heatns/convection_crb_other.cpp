/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-03-04

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)

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
   \file convection_other.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Elisa Schenone
   \date 2012-08-13
 */
#include <boost/lexical_cast.hpp>

#include "convection_crb.hpp"

// Gmsh geometry/mesh generator
#include <feel/feelfilters/gmsh.hpp>

// gmsh importer
#include <feel/feelfilters/gmsh.hpp>


// ****** CONSTRUCTEURS ****** //

ConvectionCrb::ConvectionCrb( )
:
M_backend( backend_type::build( BACKEND_PETSC ) ),
exporter( Exporter<mesh_type>::New( "ensight" ) ),
M_Dmu( new parameterspace_type )
{
}

ConvectionCrb::ConvectionCrb( po::variables_map const& vm )
:
M_vm( vm ),
M_backend( backend_type::build( vm ) ),
exporter( Exporter<mesh_type>::New( vm, "convection" ) ),
M_Dmu( new parameterspace_type )
{
}

// <int Order_s, int Order_p, int Order_t>
Feel::gmsh_ptrtype
ConvectionCrb::createMesh()
{

    timers["mesh"].first.restart();
    gmsh_ptrtype gmshp( new gmsh_type );
    gmshp->setWorldComm( Environment::worldComm() );

    double h = this->vm()["hsize"]. as<double>();
    double l = this->vm()["length"]. as<double>();

    std::ostringstream ostr;

    ostr << gmshp->preamble()
         << "a=" << 0 << ";\n"
         << "b=" << l << ";\n"
         << "c=" << 0 << ";\n"
         << "d=" << 1 << ";\n"
         << "hBis=" << h << ";\n"
         << "Point(1) = {a,c,0.0,hBis};\n"
         << "Point(2) = {b,c,0.0,hBis};\n"
         << "Point(3) = {b,d,0.0,hBis};\n"
         << "Point(4) = {a,d,0.0,hBis};\n"
         << "Point(5) = {b/2,d,0.0,hBis};\n"
         << "Point(6) = {b/2,d/2,0.0,hBis};\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Line(3) = {3,5};\n"
         << "Line(4) = {5,4};\n"
         << "Line(5) = {4,1};\n"
         << "Line(6) = {5,6};\n"
         << "Line Loop(7) = {1,2,3,4,5,6};\n"
         //<< "Line Loop(7) = {1,2,3,4,5};\n"
         << "Plane Surface(8) = {7};\n";
#if CONVECTION_DIM == 2
    ostr << "Physical Line(\"Tinsulated\") = {1,3,4};\n"
         << "Physical Line(\"Tfixed\") = {5};\n"
         << "Physical Line(\"Tflux\") = {2};\n"
        // << "Physical Line(\"Fflux\") = {6};\n"
         << "Physical Line(\"F.wall\") = {3, 4, 5, 1, 2};\n"
         << "Physical Surface(\"domain\") = {8};\n";
#else
    ostr << "Extrude {0, 0, 1} {\n"
         << "   Surface{8};\n"
         << "}\n"
         << "Physical Surface(\"Tfixed\") = {35};\n"
         << "Physical Surface(\"Tflux\") = {23};\n"
         << "Physical Surface(\"Tinsulated\") = {19, 40, 8, 31, 27};\n"
        //<< "Physical Surface(\"Fflux\") = {39};\n"
         << "Physical Surface(\"F.wall\") = {31, 27, 23, 19, 35, 40, 8};\n"
         << "Physical Volume(\"domain\") = {1};\n";
#endif
    std::ostringstream fname;
    fname << "domain";

    gmshp->setPrefix( fname.str() );
    gmshp->setDescription( ostr.str() );
    LOG(INFO) << "[timer] createMesh(): " << timers["mesh"].second << "\n";

    return gmshp;

}



void
ConvectionCrb::exportResults( element_type& U )
{
    exporter->step( 0 )->setMesh( P1h->mesh() );
    exporter->step( 0 )->add( "u", U. element<0>() );
    exporter->step( 0 )->add( "p", U. element<1>() );
    exporter->step( 0 )->add( "T", U. element<2>() );
    exporter->save();

}

// <int Order_s, int Order_p, int Order_t>
void ConvectionCrb ::exportResults( element_ptrtype& U, int i )
{
    exporter->step( i )->setMesh( P1h->mesh() );
    exporter->step( i )->add( "u", U-> element<0>() );
    exporter->step( i )->add( "p", U-> element<1>() );
    exporter->step( i )->add( "T", U-> element<2>() );
    exporter->save();
}

// <int Order_s, int Order_p, int Order_t>
void ConvectionCrb ::exportResults( element_type& U, double t )
{
    exporter->step( t )->setMesh( P1h->mesh() );
    exporter->step( t )->add( "u", U. element<0>() );
    exporter->step( t )->add( "p", U. element<1>() );
    exporter->step( t )->add( "T", U. element<2>() );
    exporter->save();
}

void
ConvectionCrb::solve( sparse_matrix_ptrtype& D,
              element_type& u,
              vector_ptrtype& F )
{

    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    M_backend->solve( D, D, U, F );
    u = *U;
}

typename ConvectionCrb::element_type
ConvectionCrb::solve( parameter_type const& mu )
{
    this->solve( mu, pT );
    return *pT;
}

void
ConvectionCrb::solve( parameter_type const& mu, element_ptrtype& T )
{
    using namespace vf;
    Feel::ParameterSpace<2>::Element M_current_mu( mu );
    static int i = 1;

    backend(_rebuild=true)->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );
    //M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobianWithoutAffineDecomposition, boost::ref( *this ), _1, _2 );
    backend()->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    vector_ptrtype R( M_backend->newVector( Xh ) );
    sparse_matrix_ptrtype J( M_backend->newMatrix( Xh,Xh ) );

    double gr = mu( 0 );
    double pr = mu( 1 );
    T->zero();
    bool use_continuity = this->vm()["use_continuity"].as<bool>();
    int N=1;


    LOG(INFO) << "nonlinear solve " << i << "  for gr=" << gr << " pr=" << pr << "\n";

    if( use_continuity )
        N=std::max( 1.0,std::max( std::ceil( std::log( gr ) ),std::ceil( std::log( pr )-std::log( 1.e-2 ) ) ) );

    for ( int i = 0; i < N; ++i )
    {

        if( use_continuity )
        {
            int denom = ( N==1 )?1:N-1;
            M_current_Grashofs = math::exp( math::log( 1. )+i*( math::log( gr )-math::log( 1. ) )/denom );
            M_current_Prandtl = math::exp( math::log( 1.e-2 )+i*( math::log( pr )-math::log( 1.e-2 ) )/denom );
        }
        else
        {
            M_current_Grashofs = gr;
            M_current_Prandtl = pr;
        }


        //std::cout << "i/N = " << i+1 << "/" << N <<std::endl;
        //std::cout << " intermediary Grashof = " << M_current_Grashofs<<std::endl;
        //std::cout<< " and Prandtl = " << M_current_Prandtl << "\n"<<std::endl;

        M_current_mu << M_current_Grashofs, M_current_Prandtl;
        this->computeBetaQm( M_current_mu );
        this->update( M_current_mu );
//        T->print(std::cout);
        //M_backend->nlSolve(_solution=T);
        //M_backend->nlSolver()->setReuse( 1, 1 );
        //std::ostringstream ostr;
        //ostr << M_backend->prefix() << "-" << i++;
        i++;
        auto p = preconditioner( _prefix=M_backend->prefix(),_pc=M_backend->pcEnumType()/*LU_PRECOND*/,
                                 _backend=M_backend,
                                 _pcfactormatsolverpackage=M_backend->matSolverPackageEnumType(),
                                 _rebuild=true );
        backend()->nlSolve(_jacobian=J , _solution=T , _residual=R, _prec=p, _reuse_jac=false, _reuse_prec=false );
#if 0
        if ( exporter->doExport() )
        {
//            T->print(std::cout);
            LOG(INFO) << "exportResults starts\n";
            this->exportResults( T, i );
            LOG(INFO) << "exportResults done\n";

        }
#endif
    }

}

void
ConvectionCrb::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    M_backend->solve( _matrix=M,  _solution=u, _rhs=f );
}

double
ConvectionCrb::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}
double
ConvectionCrb::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}

void
ConvectionCrb::run( const double * X, unsigned long N, double * Y, unsigned long P )
{

/*    using namespace vf;
    Feel::ParameterSpace<2>::Element mu( M_Dmu );
    mu << X[0], X[1];
    static int do_init = true;

    if ( do_init )
    {
//        double meshSize = X[2];
        this->init();
        do_init = false;
    }

    this->solve( mu, pT );
 */
//    double mean = integrate( elements( mesh ), chi( ( Px() >= -0.1 ) && ( Px() <= 0.1 ) )*idv( *pT ) ).evaluate()( 0,0 )/0.2;
//    Y[0]=mean;

}



double
ConvectionCrb::output( int output_index, parameter_type const& mu , element_type& unknown, bool need_to_solve)
{
    using namespace vf;
    //this->solve( mu, pT );

    auto mesh = Xh->mesh();
    auto U = Xh->element( "u" );
    U = *pT;

    element_0_type u = U. element<0>(); // fonction vitesse
    element_1_type p = U. element<1>(); // fonction pression
    element_2_type t = U. element<2>(); // fonction temperature
#if defined( FEELPP_USE_LM )
    element_3_type xi = U. element<3>(); // fonction multipliers
#endif

    // value mean-pressure
    double meas = integrate( elements( mesh ),constant( 1.0 )  ).evaluate()( 0, 0 );
    //std::cout << "measure(Omega)=" << meas << " (should be equal to 1)\n";
    //std::cout << "mean pressure = "
    //<< integrate( elements( mesh ) ,idv( p ) ).evaluate()( 0,0 )/meas << "\n";

#if defined( FEELPP_USE_LM )
    LOG(INFO) << "value of the Lagrange multiplier xi= " << xi( 0 ) << "\n";
    // std::cout << "value of the Lagrange multiplier xi= " << xi( 0 ) << "\n";
#endif

    double mean_div_u = integrate( elements( mesh ),
                                  divv( u ) ).evaluate()( 0, 0 );
    //std::cout << "mean_div(u)=" << mean_div_u << "\n";

    double div_u_error_L2 = integrate( elements( mesh ),
                                      divv( u )*divv( u ) ).evaluate()( 0, 0 );
    //std::cout << "||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";

    double AverageTdomain = integrate( elements( mesh ) , idv( t ) ).evaluate()( 0,0 ) ;
    //std::cout << "AverageTdomain = " << AverageTdomain << std::endl;


    double output = 0.0;

    // right hand side (compliant)
    if ( output_index == 0 )
    {
        output = 0.0 ;

    }

    // output
    if ( output_index == 1 )
    {
        // calcul le nombre de Nusselt
        double AverageT = integrate( markedfaces( mesh,"Tflux" ) ,
                                    idv( t ) ).evaluate()( 0,0 ) ;
        //std::cout << "AverageT = " << AverageT << std::endl;

        output = AverageT;

    }

    if( output_index == 2 )
    {
        //auto ux = u.comp(X);
        auto ux = u.comp<X>();
        //output = integrate( elements(mesh) ,  trans( idv( u ) )*idv( u )  ).evaluate()( 0,0 ) ;
        //double output2 = integrate( elements(mesh) ,  idv( u.comp(X) )  ).evaluate()( 0,0 ) ;
        //double output4 = integrate( elements(mesh) ,  idv( ux )  ).evaluate()( 0,0 ) ;
        output = integrate( elements(mesh) ,  idv( u )(0)  ).evaluate()( 0,0 ) ;
        //double output3 = integrate( elements(mesh) ,  trans( idv( u ) ) * vec( cst(1) , cst(0) ) ).evaluate()( 0,0 ) ;
        //std::cout<<"output = "<<output<<" et output2 = "<<output2<<" et output3 = "<<output3<<" et output4 = "<<output4<<std::endl;
    }

    return output;

}

// instantiation
// class ConvectionCrb<2,1,2>;
