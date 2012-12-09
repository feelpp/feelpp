/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-11-15

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
   \file opusmodelthermal.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-11-15
 */
#ifndef _OPUSMODELTHERMAL_HPP_
#define _OPUSMODELTHERMAL_HPP_ 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/foreach.hpp>

#include <feel/feel.hpp>
//#include <feel/options.hpp>
//#include <feel/feelalg/backend.hpp>
//#include <feel/feelalg/solvereigen.hpp>
//#include <feel/feeldiscr/functionspace.hpp>
//#include <feel/feeldiscr/region.hpp>
//#include <feel/feeldiscr/operatorlagrangep1.hpp>
//#include <feel/feelpoly/im.hpp>
//#include <feel/feeldiscr/bdf2.hpp>
//#include <feel/feelvf/vf.hpp>


#include <opusdefs.hpp>

namespace Feel
{
Feel::po::options_description opusModelThermalOptions();
/**
 * \addtogroup Models
 * \\@{
 */
/**
 * \class OpusModelThermal
 * \brief Opus heat transfer Model
 *
 * This class implements the heat transfert model. It is a template
 * class parametrized by the function space type associated with the
 * temperature field.
 *
 */
template<typename SpaceType>
class OpusModelThermal : public OpusModelBase
{
    typedef OpusModelBase super;
public:
#define Entity Simplex
    /**
     * Typedefs  and Constants
     */
    static const uint16_type Dim = SpaceType::nDim;
    static const uint16_type Order = SpaceType::basis_type::nOrder;

    static const uint16_type imOrder = 2*Order;
    static const uint16_type GeoOrder = 1;
    typedef OpusModelThermal<SpaceType> self_type;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef SpaceType functionspace_type;
    typedef boost::shared_ptr<SpaceType> functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;

    /**
     * constructor: Xh space and some space functions are initialized
     */
    OpusModelThermal( po::variables_map const& _vm, functionspace_ptrtype const& Xh );

    ~OpusModelThermal() {}

    // \return initial temperature and on Gamma_4
    double T0() const
    {
        return M_T0;
    }

    // \return stabilization parameter
    double gammaTemp() const
    {
        return M_gamma_temp;
    }

    template< typename MassExpr,
              typename DiffExpr,
              typename ConvExpr,
              typename RhsExpr >
    void update( double time, MassExpr const& mass, DiffExpr const& diff, ConvExpr const& conv, RhsExpr const& rhsExpr );




    void solve( element_type& T );
    void solve( boost::shared_ptr<element_type>& T )
    {
        solve( *T );
    }

    //! load from the property tree
    //void load( const boost::property_tree::ptree& pt );

    //! save the property tree
    //void save( boost::property_tree::ptree& pt );

    //! get the thermal conductance between PCB and IC{1,2}
    void thermalConductance() const
    {
        return M_c;
    }

    //! set the thermal conductance between PCB and IC{1,2}
    void setThermalConductance( double c )
    {
        M_c = c;
    }

private:

    po::variables_map vm;

    //! initial and Gamma_4 boundary temperature
    double M_T0;

    //! penalisation coefficient temperature problem stabilisation
    double M_gamma_temp;
    double M_c;
    bool M_stab;

    backend_ptrtype M_backend;

    sparse_matrix_ptrtype M_D;
    sparse_matrix_ptrtype M_Dt;
    sparse_matrix_ptrtype M_M;
    vector_ptrtype M_F;

    functionspace_ptrtype M_Xh;
}; // OpusModelThermal

template<typename SpaceType>
OpusModelThermal<SpaceType>::OpusModelThermal( po::variables_map const& _vm,
        functionspace_ptrtype const& Xh )
    :
    super(),
    vm( _vm ),
    M_T0( vm["thermal.T0"].template as<double>() ),
    M_gamma_temp( vm["thermal.gamma-temp"].template as<double>() ),
    M_c( vm["thermal.c"].template as<double>() ),
    M_stab( vm["thermal.stab"].template as<bool>() ),

    M_backend( backend_type::build( vm, "thermal" ) ),
    M_D(),
    M_Dt(),
    M_M(),
    M_F(),

    M_Xh( Xh )
{
    LOG(INFO) << "[OpusModelThermal] constructor start\n";
    FEELPP_ASSERT( M_backend != 0 ).error( "[OpusModelThermal] invalid backend" );
    FEELPP_ASSERT( M_Xh != 0 ).error( "[OpusModelThermal] invalid functionspace_ptrtype" );


    M_D = M_backend->newMatrix( M_Xh, M_Xh );
    FEELPP_ASSERT( M_D != 0 ).error( "invalid CDR matrix" );
    LOG(INFO) << "[OpusModelThermal] D allocated\n";

    M_Dt = M_backend->newMatrix( M_Xh, M_Xh );
    FEELPP_ASSERT( M_D != 0 ).error( "invalid CDR matrix" );
    LOG(INFO) << "[OpusModelThermal] Dt allocated\n";

    M_M = M_backend->newMatrix( M_Xh, M_Xh );
    FEELPP_ASSERT( M_M != 0 ).error( "invalid mass matrix" );
    LOG(INFO) << "[OpusModelThermal] M allocated\n";

    M_F = M_backend->newVector( M_Xh );
    FEELPP_ASSERT( M_F != 0 ).error( "invalid vector" );
    LOG(INFO) << "[OpusModelThermal] F allocated\n";
    sparse_matrix_ptrtype M( M_backend->newMatrix( M_Xh, M_Xh ) );

}


template<typename SpaceType>
template< typename MassExpr,
          typename DiffExpr,
          typename ConvExpr,
          typename RhsExpr >
void
OpusModelThermal<SpaceType>::update( double time,
                                     MassExpr const& mass_coeff, // rank-0
                                     DiffExpr const& diff_coeff, // rank-0
                                     ConvExpr const& conv_coeff, // rank-1 (column vector)
                                     RhsExpr const& rhs_coeff    // rank-0
                                   )
{
    using namespace vf;
    element_type u( M_Xh, "trial" );
    element_type w( M_Xh, "test" );
    static bool do_init = true;
    // reinitialize completely if
    //  - steady case
    //  -
    do_init = ( vm["steady"].template as<bool>() ||
                ( math::abs( time - ( vm["bdf.time-initial"].template as<double>() + vm["bdf.time-step"].template as<double>() ) ) < 1e-10 ) );

    if ( do_init )
        std::cout << "  -- initialize thermal model\n";

    std::vector<std::string> markers;
    markers.push_back( "Gamma_4_AIR1" );
    markers.push_back( "Gamma_4_AIR4" );
    markers.push_back( "Gamma_4_PCB" );

    size_type pattern = Pattern::COUPLED | Pattern::EXTENDED;

    if ( do_init )
        //if ( 1 )
    {
        //M_D->zero();
        form2( _test=M_Xh, _trial=M_Xh, _matrix=M_D, _init=do_init, _pattern=pattern )
            = integrate( _range=elements( M_Xh->mesh() ),
                         _expr=( mass_coeff )*idt( u )*id( w )+( diff_coeff )*gradt( u )*trans( grad( w ) ) );

#if defined( OPUS_WITH_THERMAL_DISCONTINUITY )
        auto N_IC_PCB = vec( constant( -1. ),constant( 0. ) );
        LOG(INFO) << "[add discontinuous interface at boundary " << M_Xh->mesh()->markerName( "Gamma_IC1_PCB" ) << "\n";

        form2( M_Xh, M_Xh, M_D ) += integrate( _range=markedfaces( M_Xh->mesh(), "Gamma_IC1_PCB" ),
                                               _expr=M_c*( trans( jump( id( w ) ) )*N_IC_PCB )*( trans( jumpt( idt( u ) ) )*N_IC_PCB ) );
        LOG(INFO) << "[add discontinuous interface at boundary " << M_Xh->mesh()->markerName( "Gamma_IC2_PCB" ) << "\n";
        form2( M_Xh, M_Xh, M_D ) += integrate( _range=markedfaces( M_Xh->mesh(), "Gamma_IC2_PCB" ),
                                               _expr=M_c*( trans( jump( id( w ) ) )*N_IC_PCB )*( trans( jumpt( idt( u ) ) )*N_IC_PCB ) );
#endif

#if 1

        if ( M_stab )
            form2( M_Xh, M_Xh, M_D ) +=
                integrate( internalfaces( M_Xh->mesh() ),
                           //markedfaces(M_Xh->mesh(), M_Xh->mesh()->markerName( "AIR4" )),
                           //val(this->gammaTemp()*(vf::pow(hFace(),2.0)/constant(std::pow(Order,3.5)))*
                           //     //sqrt(trans(conv_coeff)*conv_coeff))*
                           //     (abs(trans(N())*conv_coeff)))*
                           0.*
                           ( jumpt( gradt( u ) ) * jump( grad( w ) ) )
                           ,_Q<2*Order>() );

#endif
        BOOST_FOREACH( std::string marker, markers )
        {
            LOG(INFO) << "[add weakbc boundary terms velocity] boundary " << marker << " id : "
                  << M_Xh->mesh()->markerName( marker ) << "\n";
            form2( _test=M_Xh, _trial=M_Xh, _matrix=M_D ) +=
                integrate( markedfaces( M_Xh->mesh(),marker ),
                           _expr=( diff_coeff )*( -gradt( u )*N()*id( w )
                                                  -grad( w )*N()*idt( u )
                                                  +this->data()->gammaBc()*idt( u )*id( w )/hFace() ) );

        }
        M_D->close();

        form2( M_Xh, M_Xh, M_M, _init=do_init ) =
            integrate( elements( M_Xh->mesh() ),
                       idt( u )*id( w ) + gradt( u )*trans( grad( w ) ),_Q<2*Order>() );
        M_M->close();
    }

#if 1

    if ( ( math::exp( -time/3 ) > 1e-14 ) || this->data()->isSteady() )
    {
        form2( M_Xh, M_Xh, M_Dt, _init=do_init, _pattern=pattern );
        {
            boost::timer ti;
            // add convection
            std::vector<std::string> air_regions = boost::assign::list_of( "AIR4" )( "AIR123" );
            BOOST_FOREACH( auto region, air_regions )
            {
                form2( M_Xh, M_Xh, M_Dt, _pattern=pattern )+=
                    integrate( markedelements( M_Xh->mesh(),region ),
                               ( gradt( u )*( conv_coeff ) )*id( w ) );
                LOG(INFO) << "[conv terms] region:" << region << " assembly time= " << ti.elapsed() << "\n";
                ti.restart();

            }
        }

#if defined( OPUS_WITH_THERMAL_DISCONTINUITY )
        auto N_IC_PCB = vec( constant( -1. ),constant( 0. ) );
        LOG(INFO) << "[add discontinuous interface at boundary "
              << M_Xh->mesh()->markerName( "Gamma_IC1_PCB" ) << "\n";

        form2( M_Xh, M_Xh, M_Dt ) +=
            integrate( markedfaces( M_Xh->mesh(), M_Xh->mesh()->markerName( "Gamma_IC1_PCB" ) ),
                       0.*( trans( jump( id( w ) ) )*N_IC_PCB )*( trans( jumpt( idt( u ) ) )*N_IC_PCB ),_Q<2*( Order )>() );
        LOG(INFO) << "[add discontinuous interface at boundary "
              << M_Xh->mesh()->markerName( "Gamma_IC2_PCB" ) << "\n";
        form2( M_Xh, M_Xh, M_Dt ) +=
            integrate( markedfaces( M_Xh->mesh(), M_Xh->mesh()->markerName( "Gamma_IC2_PCB" ) ),
                       0.*( trans( jump( id( w ) ) )*N_IC_PCB )*( trans( jumpt( idt( u ) ) )*N_IC_PCB ), _Q<2*( Order )>() );
#endif

        // add stabilisation
        if ( M_stab )
        {
            size_type n =
                std::distance( boost::get<1>( markedfaces( M_Xh->mesh(),
                                              M_Xh->mesh()->markerName( "AIR4" ) ) ),
                               boost::get<2>( markedfaces( M_Xh->mesh(),
                                              M_Xh->mesh()->markerName( "AIR4" ) ) ) );
            size_type nt = std::distance( boost::get<1>( internalfaces( M_Xh->mesh() ) ),
                                          boost::get<2>( internalfaces( M_Xh->mesh() ) ) );

            LOG(INFO) << "n AIR faces = " << n << "\n";
            LOG(INFO) << "n total faces = " << nt << "\n";
            boost::timer ti;
#if defined(FEELPP_HAS_GOOGLE_PROFILER_H )
            //ProfilerStart( "faces_integration_profile" );
#endif

            form2( M_Xh, M_Xh, M_Dt ) +=
                integrate( // markedfaces(M_Xh->mesh(), M_Xh->mesh()->markerName( "AIR4" )),
                    internalfaces( M_Xh->mesh() ),
                    val( this->gammaTemp()*( vf::pow( hFace(),2.0 )/constant( std::pow( Order,3.5 ) ) )*
                         ( abs( trans( N() )*conv_coeff ) ) )*
                    ( jumpt( gradt( u ) ) * jump( grad( w ) ) )
                    ,_Q<2*( Order )>() );
#if defined(FEELPP_HAS_GOOGLE_PROFILER_H )
            //ProfilerStop();
#endif

            LOG(INFO) << "[stab terms] time= " << ti.elapsed() << "\n";
        }

        M_Dt->close();
        // add constant wrt time contribution
        M_Dt->addMatrix( 1.0, M_D );
    }

#endif
    // right hand side
    form1( M_Xh, M_F, _init=do_init ) = integrate( _range=elements( M_Xh->mesh() ), _expr=( rhs_coeff )*id( w ) );
    BOOST_FOREACH( std::string marker, markers )
    {
        LOG(INFO) << "[add weakbc boundary terms velocity] boundary "
              << marker << " id : " << M_Xh->mesh()->markerName( marker ) << "\n";
        form1( _test=M_Xh, _vector=M_F ) +=
            integrate( _range=markedfaces( M_Xh->mesh(), marker ),
                       _expr=constant( M_T0 )*( diff_coeff )*( -grad( w )*N()+
                               this->data()->gammaBc()*id( w )/hFace() ) );
        LOG(INFO) << "[add weakbc boundary terms velocity] boundary "
              << marker << " id : " << M_Xh->mesh()->markerName( marker ) << "done \n";
    }
    M_F->close();


}


template<typename SpaceType>
void
OpusModelThermal<SpaceType>::solve( element_type& T )
{
    M_backend->solve( _matrix=M_Dt,  _solution=T, _rhs=M_F );
}

#if 0
template<typename SpaceType>
void
OpusModelThermal<SpaceType>::load( const boost::property_tree::ptree& pt )
{
    // load c and stab
}

template<typename SpaceType>
void
OpusModelThermal<SpaceType>::save( boost::property_tree::ptree& pt )
{
    // save c and stab
}
#endif // 0

} // Feel
/** \\@} */
#endif



