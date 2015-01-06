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
#include <feel/feelcrb/eim.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

typedef Eigen::MatrixXd matrixN_type;
typedef Eigen::VectorXd vectorN_type;

typedef std::vector< std::vector< double > > beta_vector_type;

void ConvectionCrb::initModel()
{
    mesh_ptrtype mesh;

    if (this->vm()["readMesh"]. as<int>()){
        std::string repository = this->vm()["input_dir"]. as<std::string>() ;
        std::string file_mesh = this->vm()["mesh_name"]. as<std::string>() ;;
        std::string complete_name = repository + file_mesh;
        LOG(INFO) << "Meshes read in file : " << complete_name <<std::endl;

        mesh  =  loadGMSHMesh( _mesh=new mesh_type,
                              _filename=complete_name,
                              _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
    }
    else{
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=createMesh(),
                               _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
    }

    Xh = space_type::New( mesh );
    RbXh = rbfunctionspace_type::New( _model=this->shared_from_this() , _mesh=mesh );

    LOG(INFO)<<"number of dofs : "<<Xh->nLocalDof()<<"\n";
    pT = element_ptrtype( new element_type( Xh ) );

    element_type U( Xh, "u" );
    element_type Un( Xh, "un" );
    element_type V( Xh, "v" );
    element_0_type u = U. element<0>(); // fonction vitesse
    element_0_type un = Un. element<0>(); // fonction vitesse
    element_0_type v = V. element<0>(); // fonction test vitesse
    element_1_type p = U. element<1>(); // fonction pression
    element_1_type pn = Un. element<1>(); // fonction pression
    element_1_type q = V. element<1>(); // fonction test pression
    element_2_type t = U. element<2>(); // fonction temperature
    element_2_type tn = Un. element<2>(); // fonction temperature
    element_2_type s = V. element<2>(); // fonction test temperature
#if defined( FEELPP_USE_LM )
    element_3_type xi = U. element<3>(); // fonction multipliers
    element_3_type eta = V. element<3>(); // fonction test multipliers
#endif

    using namespace Feel::vf;
    Feel::ParameterSpace<2>::Element mu_min( M_Dmu );
    double Grmin( this->vm()["Grmin"]. as<double>() );
    double Prmin( this->vm()["Prmin"]. as<double>() );
    mu_min << Grmin, Prmin;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<2>::Element mu_max( M_Dmu );
    double Grmax( this->vm()["Grmax"]. as<double>() );
    double Prmax( this->vm()["Prmax"]. as<double>() );
    mu_max << Grmax, Prmax;
    M_Dmu->setMax( mu_max );


#if 0
    //==== EIM
    auto Pset = M_Dmu->sampling();
    int sampling_size = M_vm["eim.sampling-size"].template as<int>();
    std::string file_name = ( boost::format("eim_Pset_%1%") % sampling_size ).str();
    std::ifstream file ( file_name );
    if( ! file )
    {
        Pset->randomize( sampling_size );
        Pset->writeOnFile(file_name);
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(file_name);
    }

    M_mu = M_Dmu->element();
    M_unknown = Xh->element();
    auto unknown_u = M_unknown.element<0>();

    auto Uspace = U_space_type::New( mesh );
    auto Tspace = T_space_type::New( mesh );
    auto uu = trans( idv( unknown_u ) )*idv( unknown_u );
    auto eim_uu = eim( _model=this,
                       _element=M_unknown,
                       _parameter=M_mu,
                       _expr=uu,
                       //_expr=sin(cst_ref(M_mu(0))),
                       //_expr=sin(cst_ref(M_mu(0))*trans( idv(u) )*idv( u )),
                       //_expr=Feel::vf::cst_ref(1),
                       _space=Uspace,
                       _name="eim_uu",
                       _options=this->vm(),
                       _sampling=Pset );
    M_funs.push_back( eim_uu );
    //=====
#endif


    //  initialisation de A0, A1, A2
    M_Aqm.resize( Qa() );
    for (int ii=0; ii<Qa(); ii++)
    {
        M_Aqm[ii].resize(1);
        M_Aqm[ii][0] = M_backend->newMatrix( Xh, Xh );
    }

    M_Fqm.resize( Nl() );
    for (int ii=0; ii<Nl(); ii++)
    {
        M_Fqm[ii].resize( Ql(ii) );
        //if( ii < 2 )
        {
            for (int jj=0; jj<Ql(ii); jj++)
            {
                M_Fqm[ii][jj].resize(1);
                M_Fqm[ii][jj][0] = M_backend->newVector( Xh );
            }
        }
        //for the 2nd output, we use EIM
#if 0
        if( ii == 2 )
        {
            M_Fqm[2][0].resize( eim_uu->mMax() );
            for(int m=0; m<eim_uu->mMax(); m++)
                M_Fqm[2][0][m] = M_backend->newVector( Xh );
        }
#endif
    }

    D = M_backend->newMatrix( Xh, Xh );
    F = M_backend->newVector( Xh );

    M_A_tril = M_backend->newMatrix( Xh , Xh );


    double gamma( this->vm()["penalbc"]. as<double>() );
    double k=this->vm()["k"]. as<double>();
    double nu=this->vm()["nu"]. as<double>();
    double rho=this->vm()["rho"]. as<double>();
    int weakdir=this->vm()["weakdir"]. as<int>();


    int adim=this->vm()["adim"]. as<int>();

    // Fluid
    // diffusion

    // M_Aqm[0] = C = grad(u)*grav(v)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0], _init=true ) = integrate( _range=elements( mesh ),
                                                                          _expr=trace( gradt( u )*trans( grad( v ) ) )  );

    // heat diffusion: M_Aqm[1] = G = grad(t)*grad(s)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[1][0], _init=true ) = integrate( _range=elements( mesh ),
                                                                          _expr=gradt( t )*trans( grad( s ) ) );

    // pressure-velocity terms
    // M_Aqm[2] = -B = -p*div(v)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0], _init=true ) = integrate ( _range=elements( mesh ), _expr=- idt( p ) * div( v ) );

    // M_Aqm[2] += B^t = q*div(u)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] ) += integrate ( _range=elements( mesh ), _expr=divt( u ) * id( q ) );

    if ( weakdir == 1 )
    {
        // weak Dirichlet condition at the walls (u=0)
        form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] ) += integrate ( marked2faces( mesh, "F.wall" ), -trans( gradt( u )*N() )*id( v ) );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] ) += integrate ( marked2faces( mesh, "F.wall" ), -trans( -idt( p )*N() )*id( v ) );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )  += integrate ( marked2faces( mesh, "F.wall" ), -trans( grad( v )*N() )*idt( u ) );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )  += integrate ( marked2faces( mesh, "F.wall" ), -trans( -id( q )*N() )*idt( u ) );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] )  += integrate ( marked2faces( mesh, "F.wall" ), +gamma*trans( idt( u ) )*id( v )/hFace() );
    }

#if defined( FEELPP_USE_LM )
    // multipliers for zero-mean pressure
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] ) += integrate ( _range=elements( mesh ), _expr=id( q )*idt( xi ) );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] ) += integrate ( _range=elements( mesh ), _expr=idt( p )*id( eta ) );
#endif


    double expansion = 1;

    if ( adim == 0 ) expansion=3.7e-3;

    // Temperature


#if CONVECTION_DIM==2
    // buyoancy forces c(theta,v)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] ) +=integrate( _range=elements( mesh ),
                                                              _expr=-expansion*idt( t )*( trans( vec( constant( 0. ),constant( 1.0 ) ) )*id( v ) ) );
#else
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] ) +=integrate( _range=elements( mesh ),
                                                              _expr=-expansion*idt( t )*( trans( vec( cst(0.), constant( 0. ),constant( 1.0 ) ) )*id( v ) ) );
#endif

    if ( weakdir == 1 )
    {
        // weak Dirichlet on temperature (T=0|left wall)
        form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[1][0] ) += integrate ( markedfaces( mesh, "Tfixed" ),
                                                                         - gradt( t )*N()*id( s ) );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[1][0] ) += integrate ( markedfaces( mesh,"Tfixed" ),
                                                                    - grad( s )*N()*idt( t ) );
        form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[2][0] ) += integrate ( markedfaces( mesh, "Tfixed"  ),
                                                                    gamma*idt( t )*id( s )/hFace() );
    }

    M = M_backend->newMatrix( Xh, Xh );

    form2( Xh, Xh, _matrix=M, _init=true ) = integrate( elements( mesh ), trans( id( v ) )*idt( u ) + trace( grad( v )*trans( gradt( u ) ) ) );

    form2( Xh, Xh, _matrix=M ) += integrate( _range=elements( mesh ), _expr=id( q )*idt( p ) + grad( q )*trans( gradt( p ) ) );

    form2( Xh, Xh, _matrix=M ) += integrate( _range=elements( mesh ), _expr=id( s )*idt( t ) + grad( s )*trans( gradt( t ) ) );

    M->close();

    form1( Xh, _vector=M_Fqm[0][0][0] ) =
        integrate ( markedfaces( mesh, "Tflux"),
                    // heat flux on the right side
                    - id( s )  );
    M_Fqm[0][0][0]->close();
    //output : \int_{\Omega} T
    form1( _test=Xh, _vector=M_Fqm[1][0][0] ) = integrate( _range=markedfaces( mesh,"Tflux" ), _expr= id(s) );//*(1.0/area) ) ;
    M_Fqm[1][0][0]->close();


    //int mMax = eim_uu->mMax();
    //for(int m=0; m<mMax; m++)
    {
        form1( _test=Xh, _vector=M_Fqm[2][0][0] ) = integrate( _range=elements(mesh), _expr=trans( id( v ) )*id( v) );//random ( false )
        //form1( _test=Xh, _vector=M_Fqm[2][0][0] ) = integrate( _range=elements(mesh), _expr=id( v )( 0 ) );//compilation failed
        //form1( _test=Xh, _vector=M_Fqm[2][0][m] ) = integrate( _range=elements(mesh), _expr=idv(eim_uu->q(m)) );
        //M_Fqm[2][0][m]->close();
    }
    //form1( _test=Xh, _vector=M_Fqm[2][0][0] ) = integrate( _range=elements(mesh), _expr=id( v )( 0 ) );
    //    form1( _test=Xh, _vector=M_Fqm[2][0][0] ) = integrate( _range=elements(mesh), _expr= trans( id(v.comp(X)) )  ) ;
    M_Fqm[2][0][0]->close();

}

// \return the number of terms in affine decomposition of left hand
// side bilinear form
int ConvectionCrb::Qa() const
{
    return 4;
}
int ConvectionCrb::QaTri() const
{
    return 1;
}
int ConvectionCrb::mMaxA( int q )
{
    if ( q < 4 )
        return 1;
    else
        throw std::logic_error( "[Model] ERROR : try to acces to mMaxA(q) with a bad value of q");
}

/**
 * there is at least one output which is the right hand side of the
 * primal problem
 *
 * \return number of outputs associated to the model
 */
int ConvectionCrb::Nl() const
{
    return 3;
}

int ConvectionCrb::mMaxF( int output_index, int q)
{
    int dumy=0;
    if( output_index < 2 )
    {
        if ( q < Ql(output_index) )
            return 1;
        else
            throw std::logic_error( "[Model] ERROR : try to acces to mMaxF(output_index,q) with a bad value of q");
    }
    if( output_index == 2 )
    {
        auto eim_uu = M_funs[0];
        return eim_uu->mMax();
    }
    if( output_index > 2 )
            throw std::logic_error( "[Model] ERROR : try to acces to mMaxF(output_index,q) with a bad value of output_index");
    return dumy;
}

/**
 * \param l the index of output
 * \return number of terms  in affine decomposition of the \p q th output term
 */
int ConvectionCrb::Ql( int l ) const
{
    if ( l == 0 ) return 1;
    return 1;
}


/**
 * \brief compute the beta coefficient for both bilinear and linear form
 * \param mu parameter to evaluate the coefficients
 */
boost::tuple<beta_vector_type, std::vector<beta_vector_type> >
ConvectionCrb::computeBetaQm( parameter_type const& mu, double time)
{
#if 0
    auto eim_uu = M_funs[0];
    int mMax = eim_uu->mMax();
    vectorN_type beta_uu = eim_uu->beta( mu );
#endif

    M_betaAqm.resize( Qa() );
    M_betaAqm[0].resize(1);
    M_betaAqm[1].resize(1);
    M_betaAqm[2].resize(1);
    M_betaAqm[3].resize(1);
    M_betaAqm[0][0] = 1/math::sqrt( mu( 0 ) );
    M_betaAqm[1][0] = 1/( mu( 1 )*math::sqrt( mu ( 0 ) ) );
    M_betaAqm[2][0] = 1;
    M_betaAqm[3][0] = 1;

    M_betaFqm.resize( Nl() );
    M_betaFqm[0].resize( Ql( 0 ) );
    M_betaFqm[0][0].resize(1);

    //M_betaFqm[0][1].resize(1);
    M_betaFqm[0][0][0] = 1/( mu( 1 )*math::sqrt( mu ( 0 ) ) );
    //M_betaFqm[0][1][0] = 1/( mu( 1 )*math::sqrt( mu ( 0 ) ) );

    M_betaFqm[1].resize( Ql( 1 ) );
    M_betaFqm[1][0].resize(1);
    M_betaFqm[1][0][0] = 1;

    double area = integrate( elements(Xh->mesh()), vf::cst(1.) ).evaluate()(0,0);
    M_betaFqm[2].resize( Ql( 2 ) );
    M_betaFqm[2][0].resize(1);
    //M_betaFqm[2][0].resize(mMax);
    //for(int m=0; m<mMax; m++)
    {
        //M_betaFqm[2][0][m] = beta_uu(m);
        M_betaFqm[2][0][0] = 1;
    }
        //M_betaFqm[2][0].resize(1);
        //M_betaFqm[2][0][0] = 1;

    return boost::make_tuple( M_betaAqm, M_betaFqm );
}

void
ConvectionCrb::update( parameter_type const& mu )
{
    D->zero();

    for ( size_type q = 0; q < (M_betaAqm.size()-1); ++q )
    {
        for ( size_type m = 0; m < mMaxA(q); ++m )
        {
            D->addMatrix( M_betaAqm[q][m] , M_Aqm[q][m] );
        }
        //D->addMatrix( M_betaAqm[q], M_Aqm[q] );
    }

    D->close();
    F->zero();

    for ( size_type q = 0; q < M_Fqm[0].size(); ++q )
    {
        for ( size_type m = 0; m < mMaxF(0,q); ++m )
        {
            F->add( M_betaFqm[0][q][m], M_Fqm[0][q][m] );
        }
        //std::cout << "[affine decomp] scale q=" << q << " with " << M_betaFqm[0][q] << "\n";
        //F->add( M_betaFqm[0][q], M_Fqm[0][q] );
    }
    F->close();

}



void ConvectionCrb ::updateJacobianWithoutAffineDecomposition( const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    boost::timer ti;
    mesh_ptrtype mesh = Xh->mesh();
    element_type U( Xh, "u" );
    element_type Un( Xh, "un" );
    element_type V( Xh, "v" );
    element_0_type u = U. element<0>(); // fonction vitesse
    element_0_type un = Un. element<0>(); // fonction vitesse
    element_0_type v = V. element<0>(); // fonction test vitesse
    element_1_type p = U. element<1>(); // fonction pression
    element_1_type pn = Un. element<1>(); // fonction pression
    element_1_type q = V. element<1>(); // fonction test pression
    element_2_type t = U. element<2>(); // fonction temperature
    element_2_type tn = Un. element<2>(); // fonction temperature
    element_2_type s = V. element<2>(); // fonction test temperature
#if defined( FEELPP_USE_LM )
    element_3_type xi = U. element<3>(); // fonction multipliers
    element_3_type eta = V. element<3>(); // fonction test multipliers
#endif

    //M_oplin = oplin_ptrtype( new oplin_type( Xh, Xh, M_backend ) );
    //M_oplin->mat().zero();

    double gr= M_current_Grashofs;
    double sqgr( 1/math::sqrt( gr ) );
    double pr = M_current_Prandtl;
    double sqgrpr( 1/( pr*math::sqrt( gr ) ) );
    double gamma( this->vm()["penalbc"]. as<double>() );
    double k=this->vm()["k"]. as<double>();
    double nu=this->vm()["nu"]. as<double>();
    double rho=this->vm()["rho"]. as<double>();
    int weakdir=this->vm()["weakdir"]. as<int>();
    int adim=this->vm()["adim"]. as<int>();
    //choix de la valeur des paramètres dimensionnés ou adimensionnés
    double a=0.0,b=0.0,c=0.0;
    double pC=1;

    if ( adim == 0 ) pC = this->vm()["pC"]. as<double>();

    if ( adim==1 )
    {
        a=1;
        b=sqgr;
        c=sqgrpr;
    }

    else
    {
        a=rho;
        b=nu;
        c=k;
    }


    //if first time
    if( ! J )
    {
        J =  M_backend->newMatrix( _test=Xh, _trial=Xh );
        M_D =  M_backend->newMatrix( _test=Xh, _trial=Xh );
        M_L =  M_backend->newMatrix( _test=Xh, _trial=Xh );

        // Fluid
        // diffusion
        auto bf = form2( _test=Xh, _trial=Xh, _matrix=M_L );
        bf =integrate( _range=elements( mesh ),
                       _expr=cst( b )*trace( gradt( u )*trans( grad( v ) ) )  );
        LOG(INFO) << "[initLinearOperator] Fluid Diffusion terms done\n";

        // pressure-velocity terms
        bf += integrate ( _range=elements( mesh ), _expr=- idt( p ) * div( v ) );
        bf += integrate ( _range=elements( mesh ), _expr=divt( u ) * id( q ) );

        if ( weakdir == 1 )
        {
            // weak Dirichlet condition at the walls (u=0)
            auto SigmaNt = ( -idt( p )*N()+cst_ref( sqgr )*gradt( u )*N() );
            auto SigmaN = ( -id( q )*N()+cst_ref( sqgr )*grad( v )*N() );
            bf  += integrate ( marked2faces( mesh, "F.wall" ), -trans( SigmaNt )*id( v ) );
            bf  += integrate ( marked2faces( mesh, "F.wall" ), -trans( SigmaN )*idt( u ) );
            bf  += integrate ( marked2faces( mesh, "F.wall" ), +gamma*trans( idt( u ) )*id( v )/hFace() );
        }

        LOG(INFO) << "[initLinearOperator] Fluid Pressure-Velocity terms done\n";

#if defined( FEELPP_USE_LM )
        // multipliers for zero-mean pressure
        bf+= integrate ( _range=elements( mesh ), _expr=id( q )*idt( xi ) );
        bf+= integrate ( _range=elements( mesh ), _expr=idt( p )*id( eta ) );
        LOG(INFO) << "[initLinearOperator] Fluid Pressure-Multipliers terms done\n";
#endif

        LOG(INFO) << "[initLinearOperator] done in " << ti.elapsed() << "s\n";
        ti.restart();

        double expansion = 1;
        if ( adim == 0 ) expansion=3.7e-3;

        // Temperature
#if CONVECTION_DIM==2
        // buyoancy forces c(theta,v)
        bf +=integrate( _range=elements( mesh ),
                        _expr=-expansion*idt( t )*( trans( vec( constant( 0. ),constant( 1.0 ) ) )*id( v ) ) );
#else
        bf +=integrate( _range=elements( mesh ),
                        _expr=-expansion*idt( t )*( trans( vec( cst(0.), constant( 0. ),constant( 1.0 ) ) )*id( v ) ) );
#endif

        LOG(INFO) << "[initLinearOperator] temperature Force terms done\n";
        // heat conduction/diffusion: e(beta1,theta,chi)+f(theta,chi)
        bf  += integrate( _range=elements( mesh ),
                          _expr=cst( c )*gradt( t )*trans( grad( s ) ) );
        LOG(INFO) << "[initLinearOperator] Temperature Diffusion terms done\n";

        if ( weakdir == 1 )
        {
            // weak Dirichlet on temperature (T=0|left wall)
            bf  += integrate ( markedfaces( mesh, "Tfixed" ),
                               - gradt( t )*N()*id( s )*cst_ref( sqgrpr ) );
            bf  += integrate ( markedfaces( mesh, "Tfixed" ),
                               - grad( s )*N()*idt( t )*cst_ref( sqgrpr ) );
            bf  += integrate ( markedfaces( mesh,"Tfixed" ),
                               gamma*idt( t )*id( s )/hFace() );
        }

        LOG(INFO) << "[initLinearOperator2] done in " << ti.elapsed() << "s\n";
        ti.restart();
    }// end of if ( ! J )
    else
    {
        M_D->zero();
        J->zero();
    }


    // Fluid-NS
    // fluid convection derivatives: attention 2 terms
    form2( _test=Xh,_trial=Xh, _matrix=M_D )  += integrate ( _range=elements( mesh ),_expr=cst( a )*trans( id( v ) )*( gradv( u ) )*idt( u ) );
    form2( _test=Xh,_trial=Xh, _matrix=M_D )  +=integrate ( _range=elements( mesh ), _expr=cst( a )*trans( id( v ) )*( gradt( u )*idv( u ) ) );

    LOG(INFO) << "[updateJacobian1] Convection terms done\n";
    std::cout << "[updateJacobian1] Convection terms done"<<std::endl;

    //
    // temperature derivatives
    //
    // heat convection by the fluid: attention 2 terms
    form2( _test=Xh, _trial=Xh, _matrix=M_D ) +=
        integrate ( elements( mesh ),
                    pC*grad( s )*( idv( t )*idt( u ) ) );

    form2( _test=Xh, _trial=Xh, _matrix=M_D ) +=
        integrate ( elements( mesh ),
                    pC*grad( s )*( idt( t )*idv( u ) ) );

    form2( _test=Xh, _trial=Xh, _matrix=M_D ) +=
        integrate ( boundaryfaces( mesh ),
                    pC*( trans( idv( u ) )*N() )*id( s )*idt( t ) );

    form2( _test=Xh, _trial=Xh, _matrix=M_D ) +=
        integrate ( boundaryfaces( mesh ),
                    pC*( trans( idt( u ) )*N() )*id( s )*idv( t ) );

    LOG(INFO) << "[updateJacobian2] Temperature convection terms done\n";


    LOG(INFO) << "[updateJacobian2] Temperature weak Dirichlet BC terms done\n";

    M_D->close();
    M_L->close();
    J->close();

    //conditions fortes de dir
    auto Rtemp =  M_backend->newVector( Xh );
    double T0 = this->vm()["T0"]. as<double>();

    if ( weakdir == 0 )
    {
        //vitesse
        form2( Xh, Xh, _matrix=M_D )  += on( boundaryfaces( mesh ),u, Rtemp,one()*0. );

        if ( adim==1 )
            //temperature
            form2( Xh, Xh, _matrix=M_D )  += on ( markedfaces( mesh, "Tfixed" ),t,Rtemp,cst( 0.0 ) );

        else
            form2( Xh, Xh, _matrix=M_D )  += on ( markedfaces( mesh, "Tfixed" ),t,Rtemp,cst( T0 ) );
    }
    LOG(INFO) << "[updateJacobian] done in " << ti.elapsed() << "s\n";

    J->addMatrix( 1.0, M_L);
    J->addMatrix( 1.0, M_D);


}

void ConvectionCrb ::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{

    LOG(INFO) << "[updateJacobian] start\n";

    auto mesh = Xh->mesh();

    auto U = Xh->element( "u" );
    U = *X;
    auto V = Xh->element( "v" );
    auto W = Xh->element( "v" );
    auto u = U. element<0>(); // fonction vitesse
    auto v = V. element<0>(); // fonction test vitesse
    auto p = U. element<1>(); // fonction pression
    auto q = V. element<1>(); // fonction test pression
    auto t = U. element<2>(); // fonction temperature
    auto s = V. element<2>(); // fonction test temperature
#if defined( FEELPP_USE_LM )
    auto xi = U. element<3>(); // fonction multipliers
    auto eta = V. element<3>(); // fonction test multipliers
#endif


    double gamma( this->vm()["penalbc"]. as<double>() );
    int adim=this->vm()["adim"]. as<int>();

    //conditions fortes de dir
    auto Rtemp =  M_backend->newVector( Xh );
    int weakdir( this->vm()["weakdir"]. as<int>() );
    double T0 = this->vm()["T0"]. as<double>();

    double pC = 1.0;
    M_Aqm[3][0]->zero();

    // Fluid-NS
    // fluid convection derivatives: attention 2 terms

    form2( _test=Xh,_trial=Xh, _matrix=M_Aqm[3][0] )  += integrate ( _range=elements( mesh ), _expr=cst( 1. )*trans( id( v ) )*( gradv( u ) )*idt( u ) );
    form2( _test=Xh,_trial=Xh, _matrix=M_Aqm[3][0] )  += integrate ( _range=elements( mesh ), _expr=cst( 1. )*trans( id( v ) )*( gradt( u ) )*idv( u ) );

    //    // temperature derivatives
    //
    // heat convection by the fluid: attention 2 terms
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[3][0] ) +=
    integrate ( elements( mesh ),
               pC*grad( s )*( idv( t )*idt( u ) ) );

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[3][0] ) +=
    integrate ( elements( mesh ),
               pC*grad( s )*( idt( t )*idv( u ) ) );


    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[3][0] ) +=
    integrate ( boundaryfaces( mesh ),
               pC*( trans( idv( u ) )*N() )*id( s )*idt( t ) );

    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[3][0] ) +=
    integrate ( boundaryfaces( mesh ),
               pC*( trans( idt( u ) )*N() )*id( s )*idv( t ) );

    if ( weakdir == 0 )
    {

        //vitesse
        form2( Xh, Xh, _matrix=M_Aqm[3][0] )  += on( boundaryfaces( mesh ),u, Rtemp,one()*0. );

        if ( adim==1 )
            //temperature
            form2( Xh, Xh, _matrix=M_Aqm[3][0] )  += on ( markedfaces( mesh, "Tfixed" ),t,Rtemp,cst( 0.0 ) );

        else
            form2( Xh, Xh, _matrix=M_Aqm[3][0] )  += on ( markedfaces( mesh, "Tfixed" ),t,Rtemp,cst( T0 ) );
    }

    M_Aqm[3][0]->close();

    J->zero();
    J->addMatrix(1.,D);
    if ( this->vm()["enable-convection-terms"]. as<bool>() )
        J->addMatrix(1.,M_Aqm[3][0]);

}


//return the jacobian matrix evaluated at X
typename ConvectionCrb::sparse_matrix_ptrtype
ConvectionCrb::jacobian( const element_type& X )
{
    sparse_matrix_ptrtype J;
    J = M_backend->newMatrix( _test=Xh, _trial=Xh );

    vector_ptrtype XX( M_backend->newVector( Xh ) );
    *XX = X;

    updateJacobian( XX,  J);

    return J;
}
//return the residual vector evaluated at X
typename ConvectionCrb::vector_ptrtype
ConvectionCrb::residual( const element_type& X )
{
    vector_ptrtype R ( M_backend->newVector( Xh ) );

    vector_ptrtype XX( M_backend->newVector( Xh ) );
    *XX = X;

    updateResidual( XX, R );

    return R;
}

typename ConvectionCrb ::sparse_matrix_ptrtype
ConvectionCrb::computeTrilinearForm( const element_type& X )
{
    auto mesh = Xh->mesh();
    auto U = Xh->element( "U" );
    U = X;
    auto V = Xh->element( "v" );
    auto W = Xh->element( "w" );
    auto u = U. element<0>(); // fonction vitesse
    auto v = V. element<0>(); // fonction test vitesse
    auto p = U. element<1>(); // fonction pression
    auto q = V. element<1>(); // fonction test pression
    auto t = U. element<2>(); // fonction temperature
    auto s = V. element<2>(); // fonction test temperature
#if defined( FEELPP_USE_LM )
    auto xi = U. element<3>(); // fonction multipliers
    auto eta = V. element<3>(); // fonction test multipliers
#endif

    //sparse_matrix_ptrtype A_tril;
    //A_tril = M_backend->newMatrix( _test=Xh, _trial=Xh );

    double gamma( this->vm()["penalbc"]. as<double>() );

    int adim=this->vm()["adim"]. as<int>();

    //conditions fortes de dir
    auto Rtemp =  M_backend->newVector( Xh );
    int weakdir( this->vm()["weakdir"]. as<int>() );
    double T0 = this->vm()["T0"]. as<double>();

    // Fluid-NS
    form2( _test=Xh,_trial=Xh, _matrix=M_A_tril ) = integrate ( _range=elements( mesh ), _expr=trans( id( v ) )*( gradt( v ) )*idv( u )  );

    //    // temperature derivatives
    //
    // heat convection by the fluid: attention 2 terms
    form2( _test=Xh, _trial=Xh, _matrix=M_A_tril ) += integrate ( elements( mesh ), gradv( t )*( idt( s )*id( v ) ) );

    form2( _test=Xh, _trial=Xh, _matrix=M_A_tril ) += integrate ( boundaryfaces( mesh ), trans( id( v ) )*N() *idt( s )*idv( t ) );

    if ( weakdir == 0 )
    {
        //vitesse
        //form2( Xh, Xh, _matrix=M_A_tril )  += on( boundaryfaces( mesh ),u, Rtemp,one()*0. );
        form2( Xh, Xh, _matrix=M_A_tril )  += on( boundaryfaces( mesh ),v, Rtemp, one()*0 );
        if ( adim==1 )
            //temperature
            //form2( Xh, Xh, _matrix=M_A_tril )  += on ( markedfaces( mesh, "Tfixed" ),t,Rtemp,cst( 0.0 ) );
            form2( Xh, Xh, _matrix=M_A_tril )  += on ( markedfaces( mesh, "Tfixed" ),s,Rtemp,cst( 0.0 ) );
        else
            //form2( Xh, Xh, _matrix=M_A_tril )  += on ( markedfaces( mesh, "Tfixed" ),t,Rtemp,cst( T0 ) );
            form2( Xh, Xh, _matrix=M_A_tril )  += on ( markedfaces( mesh, "Tfixed" ),s,Rtemp,cst( T0 ) );
    }
    //std::cout<<"proc "<< Environment::worldComm().globalRank() << "total time for Model computeTrilinearForm : "<<ti.elapsed()<<"s end"<<std::endl;
    return M_A_tril;

}

// instantiation
// class ConvectionCrb<2,1,2>;
