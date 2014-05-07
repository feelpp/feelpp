/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-12-10

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file opusmodel.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-12-10
 */
#if !defined(OPUSMODELRB_IMPL_HPP_)
#define OPUSMODELRB_IMPL_HPP_ 1

#include <feel/feel.hpp>

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/LU>


#include<feel/feelcore/debugeigen.hpp>
#include<opusmodelrb.hpp>
//#include <opusoffline.hpp>
//#include <opusonlines1.hpp>

namespace Feel
{
using namespace Eigen;
using namespace vf;

template<int OrderU, int OrderP, int OrderT>
OpusModelRB<OrderU,OrderP,OrderT>::OpusModelRB( OpusModelRB const& om )
    :
    super( om ),
    M_is_initialized( om.M_is_initialized ),
    M_mesh( om.M_mesh ),
    M_Dmu( new parameterspace_type )
{
    initParametrization();
}
template<int OrderU, int OrderP, int OrderT>
OpusModelRB<OrderU,OrderP,OrderT>::OpusModelRB( po::variables_map const& vm )
    :
    //super( 2, vm ),
    super( vm ),
    backend( backend_type::build( vm, "backend.crb.fem" ) ),
    backendM( backend_type::build( vm, "backend.crb.norm" ) ),
    M_meshSize( vm["hsize"].template as<double>() ),
    M_is_steady( vm["steady"].template as<bool>() ),
    M_is_initialized( false ),
    M_mesh( new mesh_type ),
    M_mesh_air( new mesh_type ),
    M_mesh_line( new mesh12_type ),
    M_mesh_cross_section_2( new mesh12_type ),
    //M_exporter( Exporter<mesh_type>::New( vm, "opus" ) ),
    M_Dmu( new parameterspace_type )

{

    LOG(INFO) << "[constructor::vm] constructor, build backend done" << "\n";
    initParametrization();
}
template<int OrderU, int OrderP, int OrderT>
OpusModelRB<OrderU,OrderP,OrderT>::OpusModelRB(  )
    :
    super(),
    backend(),
    backendM(),
    M_meshSize( 1 ),
    M_is_steady( true ),
    M_is_initialized( false ),
    M_mesh( new mesh_type ),
    M_mesh_air( new mesh_type ),
    M_mesh_line( new mesh12_type ),
    M_mesh_cross_section_2( new mesh12_type ),
    //M_exporter( Exporter<mesh_type>::New( "ensight", "opus" ) ),
    M_Dmu( new parameterspace_type )

{
    LOG(INFO) << "[default] constructor, build backend" << "\n";
    backend = backend_type::build( BACKEND_PETSC );
    backendM = backend_type::build( BACKEND_PETSC );
    LOG(INFO) << "[default] constructor, build backend done" << "\n";
    initParametrization();
    LOG(INFO) << "[default] init done" << "\n";
}
template<int OrderU, int OrderP, int OrderT>
void
OpusModelRB<OrderU,OrderP,OrderT>::initParametrization()
{
    parameter_type mu_min( M_Dmu );
    //mu_min << 0.2, 1e-5, 1e6, 0.1, 5e-2;
    mu_min << 0.2, 1e-5, 1e6, 0.1, 4e-3;
    //mu_min << 0.2, 1e-5, 0, 0.1, 4e-3;
    M_Dmu->setMin( mu_min );
    parameter_type mu_max( M_Dmu );
    mu_max << 150, 1e-2, 1e6, 1e2, 5e-2;
    M_Dmu->setMax( mu_max );

    std::cout << "  -- Dmu min : "  << M_Dmu->min() << "\n";
    std::cout << "  -- Dmu max : "  << M_Dmu->max() << "\n";

}

template<int OrderU, int OrderP, int OrderT>
void
OpusModelRB<OrderU, OrderP, OrderT>::initializationField( element_ptrtype& initial_field,parameter_type const& mu )
{
    initial_field->setOnes();
    initial_field->scale( M_T0 );
}

template<int OrderU, int OrderP, int OrderT>
void
OpusModelRB<OrderU,OrderP,OrderT>::initModel()
{
    LOG(INFO) << " -- OpusModelRB::init\n";
    LOG(INFO) << "   - initialized: " << M_is_initialized << "\n";

    if ( M_is_initialized ) return;

    M_is_initialized = true;

    double e_AIR_ref = 5e-2; // m
#if 0
    typedef Gmsh gmsh_type;
    typedef boost::shared_ptr<gmsh_type> gmsh_ptrtype;

    std::string mesh_name, mesh_desc;
    gmsh_type gmsh;
    gmsh_ptrtype Gmsh_ptrtype;

    //boost::tie( mesh_name, mesh_desc, Gmsh_ptrtype ) = this->data()->createMesh( M_meshSize );
    Gmsh_ptrtype  = this->data()->createMesh( M_meshSize );
    //warning : now mesh_name is empty and fname = .msh
    std::string fname = gmsh.generate( mesh_name, mesh_desc );

    LOG(INFO) << "Generated mesh thermal\n";
    ImporterGmsh<mesh_type> import( fname );
    M_mesh->accept( import );
    LOG(INFO) << "Imported mesh thermal\n";

    Gmsh_ptrtype  = this->data()->createMeshLine( 1 );
    fname = gmsh.generate( mesh_name, mesh_desc );
    ImporterGmsh<mesh12_type> import12( fname );
    M_mesh_line->accept( import12 );
    LOG(INFO) << "Imported mesh line\n";

    Gmsh_ptrtype  = this->data()->createMeshCrossSection2( 0.2 );
    fname = gmsh.generate( mesh_name, mesh_desc );
    ImporterGmsh<mesh12_type> import_cross_section2( fname );
    M_mesh_cross_section_2->accept( import_cross_section2 );
    LOG(INFO) << "Imported mesh cross section 2\n";
#else
    LOG(INFO) << "   - Loading mesh thermal h=" << M_meshSize << "\n";
    M_mesh = createGMSHMesh( _mesh=new mesh_type,
                             _desc =  this->data()->createMesh( M_meshSize, true ),
                             _update = MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
    LOG(INFO) << "   - Imported mesh thermal\n";
    M_mesh_line = createGMSHMesh( _mesh=new mesh12_type,
                                  _desc =  this->data()->createMeshLine( 1 ),
                                  _update = MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
    LOG(INFO) << "   - Imported mesh line\n";
    M_mesh_cross_section_2 = createGMSHMesh( _mesh=new mesh12_type,
                             _desc =  this->data()->createMeshCrossSection2( 0.2 ),
                             _update = MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
    LOG(INFO) << "   - [init] Imported mesh cross section 2\n";
#endif // 0

    M_P1h = p1_functionspace_type::New( M_mesh_line );
    LOG(INFO) << "   - P1h built\n";
    M_P0h = p0_space_type::New( M_mesh );
    LOG(INFO) << "   - P0h built\n";
    typedef typename node<double>::type node_type;
    node_type period( 2 );
    //period[0]=this->data()->component("PCB").e()+this->data()->component("AIR").e();
    period[0]=this->data()->component( "PCB" ).e()+e_AIR_ref;
    period[1]=0;
    LOG(INFO) << "   - period built\n";

    //M_Th = temp_functionspace_type::New( _mesh=M_mesh,
    //                                     _periodicity=Periodic<1,2,value_type>( period ) );
    M_Th = temp_functionspace_type::New( _mesh=M_mesh,
                                         _periodicity=periodicity(Periodic<>( 1, 2 , period ) ) );
    M_RbTh = temp_rbfunctionspace_type::New(_model=this->shared_from_this(),
                                            _mesh=M_mesh,
                                            _periodicity=periodicity(Periodic<>( 1, 2 , period ) ) );
    LOG(INFO) << "   - Th built\n";
    //M_grad_Th = grad_temp_functionspace_type::New( _mesh=M_mesh );
    //LOG(INFO) << "grad Th built\n";

    pT = element_ptrtype( new element_type( M_Th ) );
    pV = element_ptrtype( new element_type( M_Th ) );

    M_bdf_poly = element_ptrtype(  new element_type( M_Th ) ) ;

    LOG(INFO) << "   - pT  built\n";

    LOG(INFO) << "   - Generated function space\n";
    LOG(INFO) << "   -  o        number of elements :  " << M_mesh->numElements() << "\n";
    LOG(INFO) << "   -  o          number of points :  " << M_mesh->numPoints() << "\n";
    LOG(INFO) << "   -  o number of local dof in Th :  " << M_Th->nLocalDof() << "\n";
    LOG(INFO) << "   -  o       number of dof in Th :  " << M_Th->nDof() << "\n";
    LOG(INFO) << "   -  o       number of dof in Th :  " << M_Th->dof()->nDof() << "\n";

    domains = p0_element_ptrtype( new p0_element_type( M_P0h, "domains" ) );
    *domains = vf::project( M_P0h, elements( M_P0h->mesh() ),
                            chi( emarker() == M_Th->mesh()->markerName( "PCB" ) )*M_Th->mesh()->markerName( "PCB" )+
                            chi( emarker() == M_Th->mesh()->markerName( "AIR123" ) )*M_Th->mesh()->markerName( "AIR123" )+
                            chi( emarker() == M_Th->mesh()->markerName( "AIR4" ) )*M_Th->mesh()->markerName( "AIR4" )+
                            chi( emarker() == M_Th->mesh()->markerName( "IC1" ) )*M_Th->mesh()->markerName( "IC1" ) +
                            chi( emarker() == M_Th->mesh()->markerName( "IC2" ) )*M_Th->mesh()->markerName( "IC2" ) );

    k = p0_element_ptrtype( new p0_element_type( M_P0h, "k" ) );
    *k = vf::project( M_P0h, elements( M_P0h->mesh() ),
                      chi( emarker() == M_Th->mesh()->markerName( "PCB" ) )*this->data()->component( "PCB" ).k()+
                      chi( emarker() == M_Th->mesh()->markerName( "AIR123" ) )*this->data()->component( "AIR" ).k()+
                      chi( emarker() == M_Th->mesh()->markerName( "AIR4" ) )*this->data()->component( "AIR" ).k()+
                      chi( emarker() == M_Th->mesh()->markerName( "IC1" ) )*this->data()->component( "IC1" ).k()+
                      chi( emarker() == M_Th->mesh()->markerName( "IC2" ) )*this->data()->component( "IC2" ).k() );
    rhoC = p0_element_ptrtype( new p0_element_type( M_P0h, "rhoC" ) );
    *rhoC = vf::project( M_P0h, elements( M_P0h->mesh() ),
                         chi( emarker() == M_Th->mesh()->markerName( "PCB" ) )*this->data()->component( "PCB" ).rhoC()+
                         chi( emarker() == M_Th->mesh()->markerName( "AIR123" ) )*this->data()->component( "AIR" ).rhoC()+
                         chi( emarker() == M_Th->mesh()->markerName( "AIR4" ) )*this->data()->component( "AIR" ).rhoC()+
                         chi( emarker() == M_Th->mesh()->markerName( "IC1" ) )*this->data()->component( "IC1" ).rhoC() +
                         chi( emarker() == M_Th->mesh()->markerName( "IC2" ) )*this->data()->component( "IC2" ).rhoC() );

    Q = p0_element_ptrtype( new p0_element_type( M_P0h, "Q" ) );
    *Q = vf::project( M_P0h, elements( M_P0h->mesh() ),
                      chi( emarker() == M_Th->mesh()->markerName( "IC1" ) )*this->data()->component( "IC1" ).Q()
                      +chi( emarker() == M_Th->mesh()->markerName( "IC2" ) )*this->data()->component( "IC2" ).Q() );
    LOG(INFO) << "   - [OpusModel::OpusModel] P0 functions allocated\n";



    double e_AIR = this->data()->component( "AIR" ).e();

    double e_PCB = this->data()->component( "PCB" ).e();
    double e_IC = this->data()->component( "IC1" ).e();
    //double L_IC = this->data()->component("IC1").h();



    auto chi_AIR = chi( Px() >= e_PCB+e_IC );
    //AUTO( ft, (constant(1.0-math::exp(-M_time/3.0 ) ) ) );
    auto ft = ( constant( 1.0 ) );
    //AUTO( vy, (constant(3.)/(2.*(e_AIR-e_IC)))*M_flow_rate*(1.-vf::pow((Px()-((e_AIR+e_IC)/2+e_PCB))/((e_AIR-e_IC)/2),2)) );
    auto vy = ( constant( 3. )/( 2.*( e_AIR-e_IC ) ) )*( 1.-vf::pow( ( Px()-( ( e_AIR+e_IC )/2+e_PCB ) )/( ( e_AIR-e_IC )/2 ),2 ) );
    //double x_mid = e_PCB+(e_IC+e_AIR)/2;
    //AUTO( vy, (constant(3)/(2*e_AIR))*M_flow_rate*(1-vf::pow((Px()-(x_mid))/(e_AIR/2),2))*ft*chi_AIR );
    //auto conv_coeff = vec( constant(0.), vy*ft*chi_AIR );
    auto conv_coeff = vec( constant( 0. ), vy );


    auto k_AIR = this->data()->component( "AIR" ).k();
    auto detJ44 = ( e_AIR - e_IC )/( e_AIR_ref - e_IC );
    auto detJinv44 = ( e_AIR_ref - e_IC )/( e_AIR - e_IC );
    auto J44 = mat<2,2>( cst( detJ44 ), cst( 0. ), cst( 0. ), cst( 1. ) );
    auto Jinv44 = mat<2,2>( cst( detJinv44 ), cst( 0. ), cst( 0. ), cst( 1. ) );
    auto K44 = k_AIR * Jinv44;

    size_type pattern = Pattern::COUPLED | Pattern::EXTENDED;

    //  initialisation de A1 et A2
    M_Aqm.resize( Qa() );
    for(int i=0; i< Qa(); i++)
        M_Aqm[i].resize(1);
    M_Aqm[0][0] = backend->newMatrix( _test=M_Th, _trial=M_Th , _pattern=pattern );
    //form2( _test=M_Th, _trial=M_Th, _matrix=M_Aqm[0] );

    for ( int q = 1; q < Qa(); ++q )
    {
        M_Aqm[q][0] = backend->newMatrix( _test=M_Th, _trial=M_Th , _pattern=pattern );
    }

    // mass matrix
    M_Mqm.resize( Qm() );

    for ( int q = 0; q < Qm(); ++q )
    {
        M_Mqm[q].resize( 1 );
        M_Mqm[q][0] = backend->newMatrix( _test=M_Th, _trial=M_Th , _pattern=pattern );
    }

    // outputs
    M_L.resize( Nl() );

    for ( int l = 0; l < Nl(); ++l )
    {
        M_L[l].resize( Ql( l ) );

        for ( int q = 0; q < Ql( l ); ++q )
        {
            M_L[l][q].resize( 1 );
            M_L[l][q][0] = backend->newVector( M_Th );
        }
    }


    D = backend->newMatrix( _test=M_Th, _trial=M_Th , _pattern=pattern );
    L.resize( Nl() );

    for ( int l = 0; l < Nl(); ++l )
    {
        L[l] = backend->newVector( M_Th );
    }

    //Mass = backend->newMatrix( M_Th, M_Th );

    M_temp_bdf = bdf( _space=M_Th, _vm=this->vm(), _name="temperature" , _prefix="temperature" );

    using namespace Feel::vf;

    element_ptrtype bdf_poly ( new element_type ( M_Th ) );

    element_type u( M_Th, "u" );
    element_type v( M_Th, "v" );
    element_type w( M_Th, "w" );

    LOG(INFO) << "   - Number of dof " << M_Th->nLocalDof() << "\n";

    M_T0 = 300;

    std::vector<std::string> markers;
    markers.push_back( "Gamma_4_AIR1" );
    markers.push_back( "Gamma_4_AIR4" );
    markers.push_back( "Gamma_4_PCB" );

    LOG(INFO) << "   - Dirichlet T0=" << M_T0 << "\n";
    double surf = integrate( markedelements( M_mesh,"IC2" ), constant( 1. ) ).evaluate()( 0, 0 );
    //
    // output 0

    form1( M_Th, M_L[0][0][0], _init=true ) =
        integrate( markedelements( M_mesh,"IC1" ),
                   id( v ) );
    form1( M_Th, M_L[0][0][0] ) +=
        integrate( markedelements( M_mesh,"IC2" ),
                   id( v ) );
    M_L[0][0][0]->close();
    form1( M_Th, M_L[0][1][0], _init=true ) =
        integrate( markedfaces( M_mesh,"Gamma_4_AIR1" ),
                   constant( M_T0 )*idv( *k )*( -grad( w )*N()+
                           this->data()->gammaBc()*id( w )/hFace() )
                 );
    form1( M_Th, M_L[0][1][0] ) +=
        integrate( markedfaces( M_mesh,"Gamma_4_PCB" ),
                   constant( M_T0 )*idv( *k )*( -grad( w )*N()+this->data()->gammaBc()*id( w )/hFace() )
                 );
    M_L[0][1][0]->close();

    // grad terms in dirichlet condition on AIR4
    // x normal derivative term
    form1( M_Th, M_L[0][2][0], _init=true ) =
        integrate( markedfaces( M_mesh,"Gamma_4_AIR4" ),
                   -constant( M_T0 )* k_AIR*dx( w )*Nx() );
    M_L[0][2][0]->close();
    // y normal derivative term
    form1( M_Th, M_L[0][3][0], _init=true ) =
        integrate( markedfaces( M_mesh,"Gamma_4_AIR4" ),
                   -constant( M_T0 )* k_AIR*dy( w )*Ny() );
    M_L[0][3][0]->close();
    // penalisation term in dirichlet constant on AIR4
    form1( M_Th, M_L[0][4][0], _init=true ) =
        integrate( markedfaces( M_mesh,"Gamma_4_AIR4" ),
                   M_T0*k_AIR*this->data()->gammaBc()*id( w )/hFace() );
    M_L[0][4][0]->close();
    LOG(INFO) << "   - rhs 0 done\n";

    // output 1
    form1( M_Th, M_L[1][0][0], _init=true ) =
        integrate( markedelements( M_mesh,"IC2" ),
                   id( v )/surf
                 );
    M_L[1][0][0]->close();
    LOG(INFO) << "   - rhs 1 done\n";
    // output 2
    // term associated with AIR3 : mult by 1/ea
    form1( M_Th, M_L[2][0][0], _init=true ) =
        integrate( markedfaces( M_mesh,"Gamma_3_AIR3" ),
                   id( v ) );
    M_L[2][0][0]->close();
    // term associated with AIR4 : mult J44/ea
    form1( M_Th, M_L[2][1][0], _init=true ) =
        integrate( markedfaces( M_mesh,"Gamma_3_AIR4" ),
                   id( v )
                 );
    M_L[2][1][0]->close();
    LOG(INFO) << "   - rhs 2 done\n";

    form1( M_Th, M_L[3][0][0], _init=true ) =
        integrate( markedfaces( M_mesh,"Gamma_3_AIR3" ),
                   id( v ) );
    M_L[3][0][0]->close();
    // term associated with AIR4 : mult J44/ea
    form1( M_Th, M_L[3][1][0], _init=true ) =
        integrate( markedfaces( M_mesh,"Gamma_3_AIR4" ),
                   id( v )
                 );
    M_L[3][1][0]->close();
    LOG(INFO) << "   - rhs 3 done\n";

    //
    // left hand side terms
    //
    //size_type pattern = Pattern::COUPLED | Pattern::EXTENDED;
    // matrix to merge all Aq
    form2( M_Th, M_Th, D, _init=true, _pattern=pattern ) =
        integrate( elements( M_mesh ), 0*idt( u )*id( v ) )+
        integrate( markedfaces( M_Th->mesh(), M_Th->mesh()->markerName( "AIR4" ) ),
                   this->data()->gammaTemp()*( vf::pow( hFace(),2.0 )/constant( std::pow( OrderT,3.5 ) ) )*
                   0.*( leftfacet( dyt( u )*Ny() ) * leftface( dy( w )*Ny() )+
                        rightfacet( dyt( u )*Ny() ) * rightface( dy( w )*Ny() )+
                        leftfacet( dyt( u )*Ny() ) * rightface( dy( w )*Ny() )+
                        rightfacet( dyt( u )*Ny() ) * leftface( dy( w )*Ny() ) ) );
    D->close();
    LOG(INFO) << "   - D  done\n";

    int AqIndex = 0;
    //test diffusion coeff
    double surfpcb = integrate( markedelements( M_mesh,"PCB" ),constant( 1.0 ) ).evaluate()( 0, 0 );
    LOG(INFO) << "   - k_PCB " << this->data()->component( "PCB" ).k() << " " <<
          integrate( markedelements( M_mesh,"PCB" ),idv( *k ) ).evaluate()( 0, 0 )/surfpcb << "\n";

    //
    // Conduction terms
    //
    // PCB + AIR123
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedelements( M_mesh,"PCB" ),
                   idv( *k )*( gradt( u )*trans( grad( v ) ) ) );
    form2( M_Th, M_Th, M_Aqm[AqIndex][0] ) +=
        integrate( markedelements( M_mesh,"AIR123" ),
                   idv( *k )*( gradt( u )*trans( grad( v ) ) ) );
    // boundary conditions (diffusion terms)
    form2( M_Th, M_Th, M_Aqm[AqIndex][0] ) +=
        integrate( markedfaces( M_mesh,"Gamma_4_AIR1" ),
                   idv( *k )*( -gradt( u )*N()*id( w )
                               -grad( w )*N()*idt( u )
                               +this->data()->gammaBc()*idt( u )*id( w )/hFace() ) );

    form2( M_Th, M_Th, M_Aqm[AqIndex][0] ) +=
        integrate( markedfaces( M_mesh,"Gamma_4_PCB" ),
                   idv( *k )*( -gradt( u )*N()*id( w )
                               -grad( w )*N()*idt( u )
                               +this->data()->gammaBc()*idt( u )*id( w )/hFace() ) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();

    // boundary condition for AIR4 (depends on ea and D)
    // x normal derivative term
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_mesh,"Gamma_4_AIR4" ),
                   k_AIR*( -dxt( u )*Nx()*id( w ) -dx( w )*Nx()*idt( u ) ) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();
    // y normal derivative term
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_mesh,"Gamma_4_AIR4" ),
                   k_AIR*( -dyt( u )*Ny()*id( w ) -dy( w )*Ny()*idt( u ) ) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();
    // penalisation term
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_mesh,"Gamma_4_AIR4" ),
                   k_AIR*this->data()->gammaBc()*idt( u )*id( w )/hFace() );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();

    //
    // IC{1,2} terms
    //
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedelements( M_mesh,"IC1" ),
                   ( gradt( u )*trans( grad( v ) ) ) );
    form2( M_Th, M_Th, M_Aqm[AqIndex][0] ) +=
        integrate( markedelements( M_mesh,"IC2" ),
                   ( gradt( u )*trans( grad( v ) ) ) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();


    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedelements( M_mesh,"AIR4" ),
                   k_AIR*dxt( u )*trans( dx( w ) ) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedelements( M_mesh,"AIR4" ),
                   k_AIR*dyt( u )*trans( dy( w ) ) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();

    //
    // Convection terms : only y derivative since v=(0,vy) and take vy = 1
    //
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedelements( M_mesh,"AIR4" ),
                   idv( *rhoC )*dyt( u )*id( w ) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedelements( M_mesh,"AIR4" ),
                   Px()*idv( *rhoC )*dyt( u )*id( w ) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedelements( M_mesh,"AIR4" ),
                   Px()*Px()*idv( *rhoC )*dyt( u )*id( w ) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();

    LOG(INFO) << "   - Aq[5]  done\n";
    // no convection AIR123
    //form2( M_Th, M_Th, M_Aqm[2] ) +=
    //integrate( markedelements(M_mesh,"AIR123"),
    //idv(*rhoC)*(gradt(u)*(conv_coeff))*id(w) );

    //form2( M_Th, M_Th, M_Aqm[2] ) +=
    //integrate( markedfaces(M_mesh,"Gamma_4_AIR1"),
    //idv(*rhoC)*(trans(N())*(conv_coeff))*id(w)*idt(u) );
#if 0
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_mesh,"Gamma_4_AIR4" ),
                   idv( *rhoC )*( Ny()*id( w )*idt( u ) ) );
    //idv(*rhoC)*(trans(N())*(conv_coeff))*id(w)*idt(u) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();

    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_mesh,"Gamma_4_AIR4" ),
                   Px()*idv( *rhoC )*( Ny()*id( w )*idt( u ) ) );
    //idv(*rhoC)*(trans(N())*(conv_coeff))*id(w)*idt(u) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();

    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_mesh,"Gamma_4_AIR4" ),
                   Px()*Px()*idv( *rhoC )*( Ny()*id( w )*idt( u ) ) );
    //idv(*rhoC)*(trans(N())*(conv_coeff))*id(w)*idt(u) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();
#endif
    //
    // Conductance terms
    //
    AUTO( N_IC_PCB,vec( constant( -1. ),constant( 0. ) ) );
    LOG(INFO) << "[add discontinuous interface at boundary " << M_Th->mesh()->markerName( "Gamma_IC1_PCB" ) << "\n";

    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_mesh, "Gamma_IC1_PCB" ),
                   ( trans( jump( id( w ) ) )*N_IC_PCB )*( trans( jumpt( idt( u ) ) )*N_IC_PCB ) );
    LOG(INFO) << "[add discontinuous interface at boundary " << M_Th->mesh()->markerName( "Gamma_IC2_PCB" ) << "\n";
    form2( M_Th, M_Th, M_Aqm[AqIndex][0] ) +=
        integrate( markedfaces( M_mesh, "Gamma_IC2_PCB" ),
                   ( trans( jump( id( w ) ) )*N_IC_PCB )*( trans( jumpt( idt( u ) ) )*N_IC_PCB ) );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();


    // stabilisation terms : AIR4 (in AIR123 velocity is zero)
    // coefficient is Jinv44(1,1) for x terms
    // coefficient is Jinv44(2,2)=1 for y terms

    // x terms
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_Th->mesh(), M_Th->mesh()->markerName( "AIR4" ) ),
                   this->data()->gammaTemp()*( vf::pow( hFace(),2.0 )/constant( std::pow( OrderT,3.5 ) ) )*leftfacev( vf::abs( Ny() ) )*
                   ( leftfacet( dxt( u )*Nx() ) * leftface( dx( w )*Nx() ) +
                     rightfacet( dxt( u )*Nx() ) * rightface( dx( w )*Nx() )+
                     leftfacet( dxt( u )*Nx() ) * rightface( dx( w )*Nx() )+
                     rightfacet( dxt( u )*Nx() ) * leftface( dx( w )*Nx() ) )
                 );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();

    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_Th->mesh(), M_Th->mesh()->markerName( "AIR4" ) ),
                   this->data()->gammaTemp()*( vf::pow( hFace(),2.0 )/constant( std::pow( OrderT,3.5 ) ) )*leftfacev( Px()*vf::abs( Ny() ) )*
                   ( leftfacet( dxt( u )*Nx() ) * leftface( dx( w )*Nx() ) +
                     rightfacet( dxt( u )*Nx() ) * rightface( dx( w )*Nx() )+
                     leftfacet( dxt( u )*Nx() ) * rightface( dx( w )*Nx() )+
                     rightfacet( dxt( u )*Nx() ) * leftface( dx( w )*Nx() ) )
                 );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_Th->mesh(), M_Th->mesh()->markerName( "AIR4" ) ),
                   this->data()->gammaTemp()*( vf::pow( hFace(),2.0 )/constant( std::pow( OrderT,3.5 ) ) )*leftfacev( Px()*Px()*vf::abs( Ny() ) )*
                   ( leftfacet( dxt( u )*Nx() ) * leftface( dx( w )*Nx() ) +
                     rightfacet( dxt( u )*Nx() ) * rightface( dx( w )*Nx() )+
                     leftfacet( dxt( u )*Nx() ) * rightface( dx( w )*Nx() )+
                     rightfacet( dxt( u )*Nx() ) * leftface( dx( w )*Nx() ) )
                 );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();


    // y terms
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_Th->mesh(), M_Th->mesh()->markerName( "AIR4" ) ),
                   this->data()->gammaTemp()*( vf::pow( hFace(),2.0 )/constant( std::pow( OrderT,3.5 ) ) )*leftfacev( vf::abs( Ny() ) )*
                   ( leftfacet( dyt( u )*Ny() ) * leftface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * rightface( dy( w )*Ny() )+
                     leftfacet( dyt( u )*Ny() ) * rightface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * leftface( dy( w )*Ny() ) )

                 );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();

    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_Th->mesh(), M_Th->mesh()->markerName( "AIR4" ) ),
                   this->data()->gammaTemp()*( vf::pow( hFace(),2.0 )/constant( std::pow( OrderT,3.5 ) ) )*leftfacev( Px()*vf::abs( Ny() ) )*
                   ( leftfacet( dyt( u )*Ny() ) * leftface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * rightface( dy( w )*Ny() )+
                     leftfacet( dyt( u )*Ny() ) * rightface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * leftface( dy( w )*Ny() ) )

                 );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();

    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_Th->mesh(), M_Th->mesh()->markerName( "AIR4" ) ),
                   this->data()->gammaTemp()*( vf::pow( hFace(),2.0 )/constant( std::pow( OrderT,3.5 ) ) )*leftfacev( Px()*Px()*vf::abs( Ny() ) )*
                   ( leftfacet( dyt( u )*Ny() ) * leftface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * rightface( dy( w )*Ny() )+
                     leftfacet( dyt( u )*Ny() ) * rightface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * leftface( dy( w )*Ny() ) )

                 );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();

    // xy terms
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_Th->mesh(), M_Th->mesh()->markerName( "AIR4" ) ),
                   this->data()->gammaTemp()*( vf::pow( hFace(),2.0 )/constant( std::pow( OrderT,3.5 ) ) )*leftfacev( vf::abs( Ny() ) )*
                   ( leftfacet( dxt( u )*Nx() ) * leftface( dy( w )*Ny() )+
                     leftfacet( dyt( u )*Ny() ) * leftface( dx( w )*Nx() )+
                     leftfacet( dxt( u )*Nx() ) * rightface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * rightface( dx( w )*Nx() )+

                     leftfacet( dyt( u )*Ny() ) * rightface( dx( w )*Nx() )+
                     rightfacet( dxt( u )*Nx() ) * leftface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * leftface( dx( w )*Nx() )+
                     rightfacet( dxt( u )*Nx() ) * rightface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * rightface( dx( w )*Nx() ) )
                 );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();

    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_Th->mesh(), M_Th->mesh()->markerName( "AIR4" ) ),
                   this->data()->gammaTemp()*( vf::pow( hFace(),2.0 )/constant( std::pow( OrderT,3.5 ) ) )*leftfacev( Px()*vf::abs( Ny() ) )*
                   ( leftfacet( dxt( u )*Nx() ) * leftface( dy( w )*Ny() )+
                     leftfacet( dyt( u )*Ny() ) * leftface( dx( w )*Nx() )+
                     leftfacet( dxt( u )*Nx() ) * rightface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * rightface( dx( w )*Nx() ) +

                     leftfacet( dyt( u )*Ny() ) * rightface( dx( w )*Nx() )+
                     rightfacet( dxt( u )*Nx() ) * leftface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * leftface( dx( w )*Nx() )+
                     rightfacet( dxt( u )*Nx() ) * rightface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * rightface( dx( w )*Nx() ) )
                 );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();
    form2( M_Th, M_Th, M_Aqm[AqIndex][0], _init=true, _pattern=pattern ) =
        integrate( markedfaces( M_Th->mesh(), M_Th->mesh()->markerName( "AIR4" ) ),
                   this->data()->gammaTemp()*( vf::pow( hFace(),2.0 )/constant( std::pow( OrderT,3.5 ) ) )*leftfacev( Px()*Px()*vf::abs( Ny() ) )*
                   ( leftfacet( dxt( u )*Nx() ) * leftface( dy( w )*Ny() )+
                     leftfacet( dyt( u )*Ny() ) * leftface( dx( w )*Nx() )+
                     leftfacet( dxt( u )*Nx() ) * rightface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * rightface( dx( w )*Nx() ) +

                     leftfacet( dyt( u )*Ny() ) * rightface( dx( w )*Nx() )+
                     rightfacet( dxt( u )*Nx() ) * leftface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * leftface( dx( w )*Nx() )+
                     rightfacet( dxt( u )*Nx() ) * rightface( dy( w )*Ny() )+
                     rightfacet( dyt( u )*Ny() ) * rightface( dx( w )*Nx() ) )
                 );
    LOG(INFO) << "   - Aq[" << AqIndex << "]  done\n";
    M_Aqm[AqIndex++][0]->close();


    //mas matrix
    //form2( M_Th, M_Th, Mass, _init=true, _pattern=pattern );


    form2( M_Th, M_Th, M_Mqm[0][0], _init=true, _pattern=pattern ) =
        integrate ( markedelements( M_mesh,"PCB" ) ,    idv( rhoC )*idt( u )*id( w ) ) +
        integrate ( markedelements( M_mesh,"IC1" ) ,    idv( rhoC )*idt( u )*id( w ) ) +
        integrate ( markedelements( M_mesh,"IC2" ) ,    idv( rhoC )*idt( u )*id( w ) ) +
        integrate ( markedelements( M_mesh,"AIR123" ) , idv( rhoC )*idt( u )*id( w ) ) ;

    form2( M_Th, M_Th, M_Mqm[1][0], _init=true, _pattern=pattern ) =
        integrate ( markedelements( M_mesh,"AIR4" ) , idv( rhoC )*idt( u )*id( w ) ) ;

    M_Mqm[0][0]->close();
    M_Mqm[1][0]->close();

    //
    // H_1 scalar product
    //
    M = backendM->newMatrix( M_Th, M_Th );

    form2( M_Th, M_Th, M, _init=true ) =
        integrate( elements( M_mesh ),
                   id( u )*idt( v )
                   +grad( u )*trans( gradt( u ) )
                 );
    M->close();




    //
    // L_2 scalar product
    //
    Mpod = backendM->newMatrix( M_Th, M_Th );

    form2( M_Th, M_Th, Mpod, _init=true ) =
        integrate( elements( M_mesh ),
                   id( u )*idt( v )
                 );
    Mpod->close();


    LOG(INFO) << "   - M and Mpod  done\n";
    LOG(INFO) << "OpusModelRB::init done\n";
}

template<int OrderU, int OrderP, int OrderT>
int
OpusModelRB<OrderU,OrderP,OrderT>::Qa() const
{
    //return 17;
    return 20;
}
template<int OrderU, int OrderP, int OrderT>
int
OpusModelRB<OrderU,OrderP,OrderT>::Qm() const
{
    return 2;
}

template<int OrderU, int OrderP, int OrderT>
int
OpusModelRB<OrderU,OrderP,OrderT>::Nl() const
{
    return 4;
}

template<int OrderU, int OrderP, int OrderT>
int
OpusModelRB<OrderU,OrderP,OrderT>::Ql( int l ) const
{
    switch ( l )
    {
    case 1 :
        return 1;

    case 2 :
        return 2;

    case 3 :
        return 2;

    default:
    case 0 :
        return 5;
    }
}

template<int OrderU, int OrderP, int OrderT>
int
OpusModelRB<OrderU,OrderP,OrderT>::mMaxA( int q )
{
    if ( q < Qa() )
        return 1;
    else
        throw std::logic_error( "[Model OpusModelRb] ERROR : try to acces to mMaxA(q) with a bad value of q");

}

template<int OrderU, int OrderP, int OrderT>
int
OpusModelRB<OrderU,OrderP,OrderT>::mMaxM( int q )
{
    if ( q < Qm() )
        return 1;
    else
        throw std::logic_error( "[Model OpusModelRb] ERROR : try to acces to mMaxM(q) with a bad value of q");
}

template<int OrderU, int OrderP, int OrderT>
int
OpusModelRB<OrderU,OrderP,OrderT>::mMaxF( int output_index, int q)
{
    int max = 0;
    if( output_index < Nl() )
    {
        if ( q < Ql( output_index ) )
            max = 1;
    }
    else
        throw std::logic_error( "[Model OpusModelRb] ERROR : try to acces to mMaxF(output_index,q) with a bad value of q");

    return max;
}


template<int OrderU, int OrderP, int OrderT>
typename OpusModelRB<OrderU,OrderP,OrderT>::beta_vectors_type
OpusModelRB<OrderU,OrderP,OrderT>::computeBetaQm( element_type const&  T, parameter_type const& mu, double time , bool only_terms_time_dependent=false )
{
    return computeBetaQm( mu , time , only_terms_time_dependent );
}

template<int OrderU, int OrderP, int OrderT>
typename OpusModelRB<OrderU,OrderP,OrderT>::beta_vectors_type
OpusModelRB<OrderU,OrderP,OrderT>::computeBetaQm( parameter_type const& mu, double time , bool only_terms_time_dependent=false )
{
    //LOG(INFO) << "[OpusModelRB::computeBetaQm] mu = " << mu << "\n";
    double kIC = mu( 0 );
    double D = mu( 1 );
    double Q = mu( 2 );
    double r = mu( 3 );
    double e_AIR = mu( 4 );

    double e_AIR_ref = 5e-2; // m
    double e_PCB = this->data()->component( "PCB" ).e();
    double e_IC = this->data()->component( "IC1" ).e();


    double c_1 = 3./( 2.*( e_AIR-e_IC ) );
    double c_2 = ( e_AIR-e_IC )/2;
    double c_3 = ( ( e_AIR+e_IC )/2+e_PCB );
    double TJ = ( e_AIR-e_IC )/( e_AIR_ref-e_IC );
    double Tb = ( e_AIR_ref-e_AIR )*( e_PCB+e_IC )/( e_AIR_ref-e_IC );

    //auto ft = (constant(1.0));
    double ft = 1.0-math::exp( -time/3 );

    //auto vy = (1.-vf::pow((Px()-((e_AIR+e_IC)/2+e_PCB))/((e_AIR-e_IC)/2),2));
    //auto conv_coeff = D*vy;
    double denom = ( ( e_IC - e_AIR )*( e_IC - e_AIR_ref )*( e_IC - e_AIR_ref ) );
    double conv1 = 6*( e_PCB + e_AIR_ref )*( e_PCB + e_IC )/denom;
    double conv2 = -  6*( 2*e_PCB + e_IC + e_AIR_ref )/denom;
    double conv3 = 6/denom;

    double k_AIR = this->data()->component( "AIR" ).k();
    double detJ44 = ( e_AIR - ( e_IC ) )/( e_AIR_ref - ( e_IC ) );
    double detJinv44 = ( e_AIR_ref - ( e_IC ) )/( e_AIR - ( e_IC ) );
    double J44xx = detJ44;
    double J44yy = 1.;
    double Jinv44xx = detJinv44;
    double Jinv44yy = 1.;

    //LOG(INFO) << "detJ44 = " << detJ44 << "\n";
    //LOG(INFO) << "D= " << D << "\n";
#if 0
    double Dnum = integrate( markedfaces( M_Th->mesh(),M_Th->mesh()->markerName( "Gamma_4_AIR4" ) ),
                             -D*( conv1+conv2*Px()+conv3*Px()*Px() )*Ny()*detJ44 ).evaluate()( 0, 0 );
#endif
    //LOG(INFO) << "Dnum= " << Dnum << "\n";

    int AqIndex = 0;
    M_betaAqm.resize( Qa() );
    for(int i=0; i<this->Qa(); i++)
        M_betaAqm[i].resize(1);

    M_betaAqm[AqIndex++][0] = 1;
    M_betaAqm[AqIndex++][0] = Jinv44xx*detJ44; //
    M_betaAqm[AqIndex++][0] = Jinv44yy*detJ44; //
    M_betaAqm[AqIndex++][0] = detJ44; //
    M_betaAqm[AqIndex++][0] = kIC; //
    M_betaAqm[AqIndex++][0] = Jinv44xx*Jinv44xx*detJ44; // AIR4  diffusion
    M_betaAqm[AqIndex++][0] = Jinv44yy*Jinv44yy*detJ44; // AIR4  diffusion
    M_betaAqm[AqIndex++][0] = ft*D*conv1*Jinv44yy*detJ44; //
    M_betaAqm[AqIndex++][0] = ft*D*conv2*Jinv44yy*detJ44; //
    M_betaAqm[AqIndex++][0] = ft*D*conv3*Jinv44yy*detJ44; //
    //M_betaAqm[AqIndex++][0] = 0*conv1*detJ44; //
    //M_betaAqm[AqIndex++][0] = 0*conv2*detJ44; //
    //M_betaAqm[AqIndex++][0] = 0*conv3*detJ44; //
    M_betaAqm[AqIndex++][0] = r; //
    M_betaAqm[AqIndex++][0] = ft*D*conv1*Jinv44xx*Jinv44xx*detJ44; // x
    M_betaAqm[AqIndex++][0] = ft*D*conv2*Jinv44xx*Jinv44xx*detJ44; // x
    M_betaAqm[AqIndex++][0] = ft*D*conv3*Jinv44xx*Jinv44xx*detJ44; // x
    M_betaAqm[AqIndex++][0] = ft*D*conv1*Jinv44yy*Jinv44yy*detJ44; // y
    M_betaAqm[AqIndex++][0] = ft*D*conv2*Jinv44yy*Jinv44yy*detJ44; // y
    M_betaAqm[AqIndex++][0] = ft*D*conv3*Jinv44yy*Jinv44yy*detJ44; // y
    M_betaAqm[AqIndex++][0] = ft*D*conv1*Jinv44xx*Jinv44yy*detJ44; // x y
    M_betaAqm[AqIndex++][0] = ft*D*conv2*Jinv44xx*Jinv44yy*detJ44; // x y
    M_betaAqm[AqIndex++][0] = ft*D*conv3*Jinv44xx*Jinv44yy*detJ44; // x y

    //LOG(INFO) << "BetaQ = " << M_betaAqm << "\n";
    M_betaL.resize( Nl() );
    // l = 0
    M_betaL[ 0 ].resize( Ql( 0 ) );
    for(int i=0; i < Ql( 0 ); i++)
        M_betaL[ 0 ][ i ].resize( 1 );
    M_betaL[ 0 ][ 0 ][ 0 ] = Q * ( 1.0-math::exp( -time ) ); //
    M_betaL[ 0 ][ 1 ][ 0 ] = 1; // start Dirichlet terms
    M_betaL[ 0 ][ 2 ][ 0 ] = Jinv44xx*detJ44; // ea : dx Nx term dirichlet
    M_betaL[ 0 ][ 3 ][ 0 ] = Jinv44yy*detJ44; // ea : dy Ny term dirichlet
    M_betaL[ 0 ][ 4 ][ 0 ] = detJ44; // ea : penalisation term dirichlet

    //LOG(INFO) << "BetaL[0] = " << M_betaL[0] << "\n";

    // l =1
    M_betaL[ 1 ].resize( Ql( 1 ) );
    for(int i=0; i < Ql( 1 ); i++)
        M_betaL[ 1 ][ i ].resize( 1 );
    M_betaL[ 1 ][ 0 ][ 0 ] = 1; //
    //LOG(INFO) << "BetaL[1] = " << M_betaL[1] << "\n";

    // l = 2
    M_betaL[ 2 ].resize( Ql( 2 ) );
    for(int i=0; i < Ql( 2 ); i++)
        M_betaL[ 2 ][ i ].resize( 1 );
    M_betaL[ 2 ][ 0 ][ 0 ]= 1./e_AIR;//1/ea; // AIR3
    M_betaL[ 2 ][ 1 ][ 0 ] = detJ44/e_AIR;//J44/ea; // AIR4
    //LOG(INFO) << "betaL[2] = " << M_betaL[2] << "\n";

    M_betaL[ 3 ].resize( Ql( 3 ) );
    for(int i=0; i < Ql( 3 ); i++)
        M_betaL[ 3 ][ i ].resize( 1 );
    M_betaL[ 3 ][ 0 ][ 0 ]= 1.;
    M_betaL[ 3 ][ 1 ][ 0 ]= detJ44;
    //LOG(INFO) << "BetaL[3] = " << M_betaL[3] << "\n";


    M_betaMqm.resize( Qm() );
    for(int i=0; i<Qm(); i++)
        M_betaMqm[i].resize( 1 );
    M_betaMqm[ 0 ][ 0 ] = 1 ;
    M_betaMqm[ 1 ][ 0 ] = detJ44;

    return boost::make_tuple( M_betaMqm, M_betaAqm, M_betaL );
}

template<int OrderU, int OrderP, int OrderT>
OpusModelRB<OrderU,OrderP,OrderT>::~OpusModelRB()
{}
template<int OrderU, int OrderP, int OrderT>
typename OpusModelRB<OrderU,OrderP,OrderT>::sparse_matrix_ptrtype
OpusModelRB<OrderU,OrderP,OrderT>::newMatrix() const
{
    auto Dnew = backend->newMatrix( M_Th, M_Th );
    *Dnew  = *D;
    Dnew->zero();
    return Dnew;
}
template<int OrderU, int OrderP, int OrderT>
typename OpusModelRB<OrderU,OrderP,OrderT>::vector_ptrtype
OpusModelRB<OrderU,OrderP,OrderT>::newVector() const
{
    return backend->newVector( M_Th );
}

template<int OrderU, int OrderP, int OrderT>
void
OpusModelRB<OrderU,OrderP,OrderT>::update( parameter_type const& mu , double time )
{
    this->computeBetaQm( mu , time );

    double Fr = mu( 1 );
    double e_AIR = mu( 4 );

    double e_AIR_ref = 5e-2; // m
    double e_PCB = this->data()->component( "PCB" ).e();
    double e_IC = this->data()->component( "IC1" ).e();

    double c_1 = 3./( 2.*( e_AIR-e_IC ) );
    double c_2 = ( e_AIR-e_IC )/2;
    double c_3 = ( ( e_AIR+e_IC )/2+e_PCB );
    double TJ = ( e_AIR-e_IC )/( e_AIR_ref-e_IC );
    double Tb = ( e_AIR_ref-e_AIR )*( e_PCB+e_IC )/( e_AIR_ref-e_IC );

    //auto ft = (constant(1.0));
    double ft = 1.0-math::exp( -time/3 );
    //auto vy = (1.-vf::pow((Px()-((e_AIR+e_IC)/2+e_PCB))/((e_AIR-e_IC)/2),2));
    //auto conv_coeff = D*vy;
    double denom = ( ( e_IC - e_AIR )*( e_IC - e_AIR_ref )*( e_IC - e_AIR_ref ) );
    double conv1 = 6*( e_PCB + e_AIR_ref )*( e_PCB + e_IC )/denom;
    double conv2 = -  6*( 2*e_PCB + e_IC + e_AIR_ref )/denom;
    double conv3 = 6/denom;

    *pV = vf::project( M_Th, markedelements( M_Th->mesh(), M_Th->mesh()->markerName( "AIR4" ) ),
                       ft*Fr*( conv1+conv2*Px()+conv3*Px()*Px() ) );
    LOG(INFO) << "[update(mu)] pV done\n";
    boost::timer ti;
    D->zero();
    //*D = M_Aqm[0][0];
    //D->scale( M_betaAqm[0][0] ) );

    for ( size_type q = 0; q < Qa(); ++q )
    {
        for ( size_type m = 0; m < mMaxA(q); ++m )
        {
            //LOG(INFO) << "[affine decomp] scale q=" << q << " with " << M_betaAqm[q][m] << "\n";
            D->addMatrix( M_betaAqm[q][m] , M_Aqm[q][m] );
        }
    }

    LOG(INFO) << "[update(mu,"<<time<<")] D assembled in " << ti.elapsed() << "s\n";
    ti.restart();

    for ( int l = 0; l < Nl(); ++l )
    {

        L[l]->zero();
        for ( size_type q = 0; q < Ql( l ); ++q )
        {
            for ( size_type m = 0; m < mMaxF( l , q ); ++m )
            {
                L[l]->add( M_betaL[ l ][ q ][ m ], M_L[ l ][ q ][ m ] );
            }
        }

        LOG(INFO) << "[update(mu,"<<time<<")] L[" << l << "] assembled in " << ti.elapsed() << "s\n";
        ti.restart();
    }


    //mass matrix contribution
    auto vec_bdf_poly = backend->newVector( M_Th );

    for ( size_type q = 0; q < Qm(); ++q )
    {
        for ( size_type m = 0; m < mMaxM(q); ++m )
        {
            //left hand side
            D->addMatrix( M_betaMqm[q][m]*M_bdf_coeff, M_Mqm[q][m] );
            //right hand side
            *vec_bdf_poly = *M_bdf_poly;
            vec_bdf_poly->scale( M_betaMqm[q][m] );
            L[0]->addVector( *vec_bdf_poly, *M_Mqm[q][m] );
        }
    }

    LOG(INFO) << "[update(mu,"<<time<<")] add mass matrix contributions in " << ti.elapsed() << "s\n";
    ti.restart();

}

template<int OrderU, int OrderP, int OrderT>
typename OpusModelRB<OrderU,OrderP,OrderT>::element_type
OpusModelRB<OrderU,OrderP,OrderT>::solve( parameter_type const& mu )
{
    //element_ptrtype T( new element_type( M_Th ) );

    this->solve( mu, pT );
    //this->exportResults( *pT );

    std::vector<double> LT( this->Nl() );

    for ( int l = 0; l < this->Nl()-1; ++l )
    {
        LT[l] = inner_product( *L[l], *pT );
        LOG(INFO) << "LT(" << l << ")=" << LT[l] << "\n";
    }

    LT[3] = inner_product( *L[3], *pV );
    LOG(INFO) << "LT(" << 3 << ")=" << LT[3] << "\n";

    return *pT;
}

template<int OrderU, int OrderP, int OrderT>
void
OpusModelRB<OrderU,OrderP,OrderT>::solve( parameter_type const& mu, element_ptrtype& T )
{
    boost::timer ti;
    //LOG(INFO) << "solve(mu,T) for parameter " << mu << "\n";
    using namespace Feel::vf;

    //initialization of temperature
    *T = vf::project( M_Th, elements( M_Th->mesh() ), constant( M_T0 ) );
    M_temp_bdf->initialize( *T );

    if ( M_is_steady )
    {
        M_temp_bdf->setSteady();
    }

    M_temp_bdf->start();
    M_bdf_coeff = M_temp_bdf->polyDerivCoefficient( 0 );

    for ( M_temp_bdf->start(); !M_temp_bdf->isFinished(); M_temp_bdf->next() )
    {
        *M_bdf_poly = M_temp_bdf->polyDeriv();
        this->update( mu , M_temp_bdf->time() );
        LOG(INFO) << "[solve(mu)] : time = "<<M_temp_bdf->time()<<"\n";
        LOG(INFO) << "[solve(mu)] update(mu) done in " << ti.elapsed() << "s\n";
        ti.restart();
        LOG(INFO) << "[solve(mu)] start solve\n";
        backend->solve( _matrix=D,  _solution=*T, _rhs=L[0] );

#if(0)
        auto ret = backend->solve( _matrix=D,  _solution=*T, _rhs=L[0], _reuse_prec=( M_temp_bdf->iteration() >=2 ) );

        if ( !ret.template get<0>() )
        {
            LOG(INFO)<<"WARNING : we have not converged ( nb_it : "<<ret.template get<1>()<<" and residual : "<<ret.template get<2>() <<" ) \n";
        }

#endif

        LOG(INFO) << "[solve(mu)] solve done in " << ti.elapsed() << "s\n";
        ti.restart();
        this->exportResults( M_temp_bdf->time(), *T , mu );
        LOG(INFO) << "[solve(mu)] export done in " << ti.elapsed() << "s\n";
        ti.restart();

        M_temp_bdf->shiftRight( *T );

    }


}

template<int OrderU, int OrderP, int OrderT>
void
OpusModelRB<OrderU,OrderP,OrderT>::solve( parameter_type const& mu, element_ptrtype& T, vector_ptrtype const& rhs, bool transpose )
{
    //LOG(INFO) << "solve(mu,T) for parameter " << mu << "\n";
    using namespace Feel::vf;


    *T = vf::project( M_Th, elements( M_Th->mesh() ), constant( M_T0 ) );
    M_temp_bdf->initialize( *T );


    if ( M_is_steady )
    {
        M_temp_bdf->setSteady();
    }

    M_temp_bdf->start();
    M_bdf_coeff = M_temp_bdf->polyDerivCoefficient( 0 );

    for ( M_temp_bdf->start(); !M_temp_bdf->isFinished(); M_temp_bdf->next() )
    {
        *M_bdf_poly = M_temp_bdf->polyDeriv();
        this->update( mu , M_temp_bdf->time() );

        if ( transpose )
        {
            auto ret = backend->solve( _matrix=D->transpose(),  _solution=*T, _rhs=rhs , _reuse_prec=( M_temp_bdf->iteration() >=2 ) );

            if ( !ret.template get<0>() )
            {
                LOG(INFO)<<"WARNING : we have not converged ( nb_it : "<<ret.template get<1>()<<" and residual : "<<ret.template get<2>() <<" ) \n";
            }

        }

        else
        {
            auto ret = backend->solve( _matrix=D,  _solution=*T, _rhs=rhs , _reuse_prec=( M_temp_bdf->iteration() >=2 ) );

            if ( !ret.template get<0>() )
            {
                LOG(INFO)<<"WARNING : we have not converged ( nb_it : "<<ret.template get<1>()<<" and residual : "<<ret.template get<2>() <<" ) \n";
            }
        }

        this->exportResults( M_temp_bdf->time(), *T ,mu );


        M_temp_bdf->shiftRight( *T );
    }




}

template<int OrderU, int OrderP, int OrderT>
void
OpusModelRB<OrderU,OrderP,OrderT>::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //LOG(INFO) << "l2solve(u,f)\n";
    //backendM->solve( _matrix=M,  _solution=u, _rhs=f, _prec=M );
    //backendM = backend_type::build( BACKEND_PETSC );
    backendM->solve( _matrix=M, _solution=u, _rhs=f );
    //LOG(INFO) << "l2solve(u,f) done\n";
}


/**
 * H1 scalar product
 */
template<int OrderU, int OrderP, int OrderT>
typename OpusModelRB<OrderU,OrderP,OrderT>::sparse_matrix_ptrtype
OpusModelRB<OrderU,OrderP,OrderT>::energyMatrix ( void )
{
    return M;
}

template<int OrderU, int OrderP, int OrderT>
double
OpusModelRB<OrderU,OrderP,OrderT>::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}
template<int OrderU, int OrderP, int OrderT>
double
OpusModelRB<OrderU,OrderP,OrderT>::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );

}

template<int OrderU, int OrderP, int OrderT>
double
OpusModelRB<OrderU,OrderP,OrderT>::scalarProductForPod( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return Mpod->energy( x, y );
}
template<int OrderU, int OrderP, int OrderT>
double
OpusModelRB<OrderU,OrderP,OrderT>::scalarProductForPod( vector_type const& x, vector_type const& y )
{
    return Mpod->energy( x, y );

}

template<int OrderU, int OrderP, int OrderT>
double
OpusModelRB<OrderU,OrderP,OrderT>::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{
    need_to_solve=true;
    if( need_to_solve )
        this->solve( mu, pT );
    else
    {
        if ( M_is_steady )
            M_temp_bdf->setSteady();
        this->computeBetaQm( mu , M_temp_bdf->timeFinal() );
        this->update( mu , M_temp_bdf->timeFinal() );
        *pT = u;
    }

    vector_ptrtype U( backend->newVector( M_Th ) );
    *U = *pT;
    LOG(INFO) << "S1 = " << inner_product( *L[1], *U ) << "\n S2 = " << inner_product( *L[2], *U ) << "\n";
    return inner_product( L[output_index], U );
}

template<int OrderU, int OrderP, int OrderT>
void
OpusModelRB<OrderU,OrderP,OrderT>::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    LOG(INFO) << "[OpusModel::run] input/output relationship\n";

    parameter_type mu( M_Dmu );

    mu << /*kIC*/X[0],/*D*/X[1], /*Q*/X[2], /*r*/X[3], /*ea*/X[4];

#if 0
    M_is_steady = X[0];

    if ( M_is_steady )
    {
        mu << /*kIC*/X[1],/*D*/X[2], /*Q*/X[3], /*r*/X[4], /*ea*/X[5];
    }

    else
    {
        mu << /*kIC*/X[1],/*D*/X[2], /*Q*/X[3], /*r*/X[4], /*ea*/X[5], /*dt*/X[6], /*Tf*/X[7];
    }

#endif

    for ( unsigned long i = 0; i < N; ++i )
        LOG(INFO) << "[OpusModelRB::run] X[" << i << "]=" << X[i] << "\n";

    this->data()->component( "IC1" ).setK( X[0] );
    this->data()->component( "IC2" ).setK( X[0] );
    this->data()->component( "AIR" ).setFlowRate( X[1] );
    this->data()->component( "IC1" ).setQ( X[2] );
    this->data()->component( "IC2" ).setQ( X[2] );

    for ( unsigned long i = 0; i < N; ++i )
        LOG(INFO) << "[OpusModel::run] X[" << i << "]=" << X[i] << "\n";

    this->data()->component( "AIR" ).setE( X[4] );
    M_meshSize = X[5];
    LOG(INFO) << "[OpusModelRB::run] parameters set\n";

    this->data()->print();

    LOG(INFO) << "[OpusModelRB::run] parameters print done\n";

    LOG(INFO) << "[OpusModelRB::run] init\n";
    this->initModel();
    LOG(INFO) << "[OpusModelRB::run] init done\n";

    *pT = vf::project( M_Th, elements( M_Th->mesh() ), constant( M_T0 ) );
    M_temp_bdf->initialize( *pT );


    this->solve( mu, pT );
    LOG(INFO) << "[OpusModelRB::run] solve done\n";

    vector_ptrtype U( backend->newVector( M_Th ) );
    *U = *pT;
    Y[0] = inner_product( *L[1], *U );
    Y[1] = inner_product( *L[2], *U );
    LOG(INFO) << "[OpusModel::run] run done, set outputs\n";

    for ( unsigned long i = 0; i < P; ++i )
        LOG(INFO) << "[OpusModel::run] Y[" << i << "]=" << Y[i] << "\n";

}
template<int OrderU, int OrderP, int OrderT>
void
OpusModelRB<OrderU,OrderP,OrderT>::run()
{}

template<int OrderU, int OrderP, int OrderT>
void
OpusModelRB<OrderU,OrderP,OrderT>::exportResults( double time, temp_element_type& T , parameter_type const& mu )
{
    std::ostringstream osstr ;

    std::string exp_name;
    export_ptrtype exporter;
    std::string mu_str;

    for ( int i=0; i<mu.size(); i++ )
    {
        mu_str= mu_str + ( boost::format( "_%1%" ) %mu[i] ).str() ;
    }

    exp_name = "T" + ( boost::format( "%1%" ) %time ).str()+"_mu_"+mu_str;

    exporter = export_ptrtype( Exporter<mesh_type>::New( "ensight", exp_name  ) );
    exporter->step( time )->setMesh( T.functionSpace()->mesh() );



    int j = time;
    osstr<<j;
    //LOG(INFO) << "exportresults : " << this->data()->doExport() << "\n";
    //if ( this->data()->doExport() )
    {
        LOG(INFO) << "exporting...\n";
        exporter->step( time )->setMesh( T.functionSpace()->mesh() );
        exporter->step( time )->add( "Domains", *domains );
        exporter->step( time )->add( "k", *k );
        exporter->step( time )->add( "rhoC", *rhoC );
        exporter->step( time )->add( "Q", *Q );
        exporter->step( time )->add( "T", T );
        exporter->step( time )->add( "V", *pV );
        //M_exporter->step(0)->add( "Velocity",  U.template element<0>() );
        //M_exporter->step(0)->add( "Pressure",  U.template element<1>() );

        using namespace  vf;
        //typename grad_temp_functionspace_type::element_type g( M_grad_Th, "k*grad(T)" );
        //g = vf::project( M_grad_Th, elements( M_grad_Th->mesh() ), trans(idv(*k)*gradv(T)) );
        //M_exporter->step(0)->add( "k_grad(T)", g );
        exporter->save();
        LOG(INFO) << "exporting done.\n";
    }
}

} // Feel

#endif

