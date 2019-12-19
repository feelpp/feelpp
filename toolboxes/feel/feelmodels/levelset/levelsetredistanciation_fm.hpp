//! -*- mode: c++; coding: utf-9; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
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
//! @file levelsetredistanciation_fm.hpp
//! @author Thibaut Metivet <thibaut.metivet@gmail.com>
//! @date 25 Jul 2019
//! @copyright 2019 Feel++ Consortium
//!

#ifndef _LEVELSET_REDISTANCIATION_FM_HPP
#define _LEVELSET_REDISTANCIATION_FM_HPP 1

#include <feel/feelmodels/levelset/levelsetredistanciation.hpp>
#include <feel/feelmodels/levelset/levelsetfilters.hpp>
#include <feel/feelmodels/modelcore/utils.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feells/fastmarching.hpp>
#include <feel/feells/distancepointtoface.hpp>

namespace Feel
{

template<typename FunctionSpaceType>
class LevelSetRedistanciationFM : 
    public LevelSetRedistanciation<FunctionSpaceType>
{
    public:
        typedef LevelSetRedistanciation<FunctionSpaceType> super_type;
        typedef LevelSetRedistanciationFM<FunctionSpaceType> self_type;
        typedef std::shared_ptr<self_type> self_ptrtype;
        //--------------------------------------------------------------------//
        // Space
        typedef FunctionSpaceType functionspace_type;
        typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;
        typedef typename functionspace_type::template Basis<0>::type basis_type;
        static const uint16_type functionSpaceOrder = functionspace_type::fe_type::nOrder;
        static const uint16_type nOrder = functionspace_type::fe_type::nOrder;
        static constexpr uint16_type nDim = functionspace_type::nDim;
        typedef typename functionspace_type::value_type value_type;
        // Element
        typedef typename functionspace_type::element_type element_type;
        typedef std::shared_ptr<element_type> element_ptrtype;
        // Mesh
        typedef typename functionspace_type::convex_type convex_type;
        typedef typename functionspace_type::mesh_type mesh_type;
        // Periodicity
        typedef typename functionspace_type::periodicity_0_type periodicity_type;
        static const bool functionSpaceIsPeriodic = functionspace_type::is_periodic;

        // Scalar projectors
        typedef Projector<functionspace_type, functionspace_type> projector_type;
        typedef std::shared_ptr<projector_type> projector_ptrtype;

        // Redist P1 space
        typedef Lagrange<1, Scalar> basis_P1_type;
        typedef FunctionSpace<
            mesh_type, 
            bases<basis_P1_type>, 
            value_type,
            Periodicity<NoPeriodicity>,
            mortars<NoMortar>
                > functionspace_P1_type;
        typedef std::shared_ptr<functionspace_P1_type> functionspace_P1_ptrtype;

        // Gradient discontinuous space
        static const uint16_type functionSpaceDiscontinuousOrder = ( functionSpaceOrder == 1 ) ? 0 : functionSpaceOrder;
        typedef typename FeelModels::detail::ChangeBasisContinuity<Discontinuous, typename FeelModels::detail::ChangeBasisOrder<functionSpaceDiscontinuousOrder, basis_type>::type>::type basis_discontinuous_type;
        //typedef Lagrange<nOrder, Scalar, Discontinuous> basis_discontinuous_type;
        typedef FunctionSpace<
            mesh_type, 
            bases<basis_discontinuous_type>, 
            value_type,
            Periodicity<NoPeriodicity>,
            mortars<NoMortar>
                > functionspace_discontinuous_type;
        typedef std::shared_ptr<functionspace_discontinuous_type> functionspace_discontinuous_ptrtype;

        // Fast-marching space
        static constexpr bool UseRedistP1Space = !( functionSpaceOrder == 1 && !functionSpaceIsPeriodic );
        typedef typename mpl::if_< 
            mpl::bool_<UseRedistP1Space>,
            functionspace_P1_type,
            functionspace_type
                >::type functionspace_FM_type;
        typedef std::shared_ptr<functionspace_FM_type> functionspace_FM_ptrtype;

        // Fast-marching
        typedef FastMarching<functionspace_FM_type> fastmarching_type;
        typedef std::shared_ptr< fastmarching_type > fastmarching_ptrtype;

        // Range
        typedef elements_reference_wrapper_t<mesh_type> range_elements_type;

        // Initialisation method
        enum class FastMarchingInitialisationMethod { 
            NONE=0, 
            ILP_NODAL, ILP_L2, ILP_SMOOTH, 
            ILDIST,
            HJ_EQ, IL_HJ_EQ
        };

        typedef boost::bimap<std::string, FastMarchingInitialisationMethod> fastmarchinginitialisationmethodidmap_type;
        static const fastmarchinginitialisationmethodidmap_type FastMarchingInitialisationMethodMap;

        //--------------------------------------------------------------------//
        //--------------------------------------------------------------------//
        //--------------------------------------------------------------------//
        // Constructor
        LevelSetRedistanciationFM( 
                functionspace_ptrtype const& space,
                std::string const& prefix = "" );

        // Projectors
        projector_ptrtype const& projectorL2( bool buildOnTheFly = true ) const;
        void setProjectorL2( projector_ptrtype const& p ) { M_projectorL2 = p; }
        projector_ptrtype const& projectorSM( bool buildOnTheFly = true ) const;
        void setProjectorSM( projector_ptrtype const& p ) { M_projectorSM = p; }

        // discontinuous space
        functionspace_discontinuous_ptrtype const& functionSpaceDiscontinuous() const { return M_spaceDiscontinuous; }
        // Fast-marching space
        functionspace_FM_ptrtype const& functionSpaceFM() const { return M_spaceFM; }

        // Fast-marching
        fastmarching_ptrtype const& fastMarching() const { return M_fastMarching; }

        // Fast-marching initialisation
        FastMarchingInitialisationMethod fastMarchingInitialisationMethod() const { return M_fastMarchingInitialisationMethod; }
        void setFastMarchingInitialisationMethod( FastMarchingInitialisationMethod m ) { M_fastMarchingInitialisationMethod = m; }

        // Redistanciation
        element_type initFastMarching( element_type const& phi, range_elements_type const& rangeInitialElts ) const;
        element_type runFastMarching( element_type const& phi, range_elements_type const& rangeInitialElts ) const;
        element_type run( element_type const& phi, range_elements_type const& rangeInitialElts ) const;
        element_type run( element_type const& phi ) const;

        // Parameters
        void loadParametersFromOptionsVm();

    private:
        //--------------------------------------------------------------------//
        // Projectors
        projector_ptrtype M_projectorL2;
        projector_ptrtype M_projectorSM;
        //--------------------------------------------------------------------//
        //--------------------------------------------------------------------//
        // Reinit P1 operators
        typedef OperatorInterpolation<
            functionspace_type, // from space
            functionspace_P1_type // to space
                > op_interpolation_to_P1_type;

        typedef OperatorInterpolation<
            functionspace_P1_type, // from space
            functionspace_type // to space
                > op_interpolation_from_P1_type;

        typedef std::shared_ptr<op_interpolation_to_P1_type> op_interpolation_to_P1_ptrtype;
        typedef std::shared_ptr<op_interpolation_from_P1_type> op_interpolation_from_P1_ptrtype;

        typedef OperatorLagrangeP1<functionspace_type> op_lagrangeP1_type;
        typedef std::shared_ptr<op_lagrangeP1_type> op_lagrangeP1_ptrtype;

        op_interpolation_to_P1_ptrtype M_opInterpolationToP1;
        op_interpolation_from_P1_ptrtype M_opInterpolationFromP1;
        op_lagrangeP1_ptrtype M_opLagrangeP1;
        //--------------------------------------------------------------------//
        // PN-1 discontinuous space
        functionspace_discontinuous_ptrtype M_spaceDiscontinuous;
        // Fast-marching space
        functionspace_FM_ptrtype M_spaceFM;

        //--------------------------------------------------------------------//
        // Fast-marching
        fastmarching_ptrtype M_fastMarching;
        FastMarchingInitialisationMethod M_fastMarchingInitialisationMethod;
        //--------------------------------------------------------------------//
        //--------------------------------------------------------------------//
        // Options
};

template<typename FunctionSpaceType>
const typename LevelSetRedistanciationFM<FunctionSpaceType>::fastmarchinginitialisationmethodidmap_type
LevelSetRedistanciationFM<FunctionSpaceType>::FastMarchingInitialisationMethodMap = 
boost::assign::list_of< typename LevelSetRedistanciationFM<FunctionSpaceType>::fastmarchinginitialisationmethodidmap_type::relation >
    ( "none", FastMarchingInitialisationMethod::NONE )
    ( "ilp-l2", FastMarchingInitialisationMethod::ILP_L2 )
    ( "ilp-smooth", FastMarchingInitialisationMethod::ILP_SMOOTH )
    ( "ilp-nodal", FastMarchingInitialisationMethod::ILP_NODAL )
    ( "ildist", FastMarchingInitialisationMethod::ILDIST )
    ( "hj", FastMarchingInitialisationMethod::HJ_EQ )
    ( "il-hj", FastMarchingInitialisationMethod::IL_HJ_EQ )
;

template<typename FunctionSpaceType>
LevelSetRedistanciationFM<FunctionSpaceType>::LevelSetRedistanciationFM( 
        functionspace_ptrtype const& space,
        std::string const& prefix )
    : super_type( space, prefix )
{
    // Load parameters
    this->loadParametersFromOptionsVm();
    // Init
    M_spaceDiscontinuous = functionspace_discontinuous_type::New(
            _mesh=space->mesh()
            );

    if constexpr( UseRedistP1Space )
    {
        M_opLagrangeP1 = lagrangeP1( 
                _space=this->functionSpace(), 
                _update=MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES 
                );
        M_spaceFM = functionspace_P1_type::New(
                _mesh=M_opLagrangeP1->mesh(),
                _periodicity=periodicity( NoPeriodicity() )
                );

        if( functionSpaceOrder > 1 )
        {
            M_opInterpolationToP1 = opInterpolation(
                    _domainSpace = this->functionSpace(),
                    _imageSpace = this->functionSpaceFM(),
                    _type = InterpolationNonConforme(false)
                    );
            M_opInterpolationFromP1 = opInterpolation(
                    _domainSpace = this->functionSpaceFM(),
                    _imageSpace = this->functionSpace(),
                    _type = InterpolationNonConforme(false)
                    );
        }
    }
    else
    {
        M_spaceFM = this->functionSpace();
    }

    M_fastMarching.reset( new fastmarching_type( M_spaceFM ) );
}

template<typename FunctionSpaceType>
typename LevelSetRedistanciationFM<FunctionSpaceType>::projector_ptrtype const& 
LevelSetRedistanciationFM<FunctionSpaceType>::projectorL2( bool buildOnTheFly ) const
{
    if( !M_projectorL2 && buildOnTheFly )
    {
        auto backendName = prefixvm( this->prefix(), "projector-l2" );
        auto backendProjectorL2 = Backend<value_type>::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpace()->worldCommPtr()
                );
        const_cast<self_type*>(this)->M_projectorL2 = projector(
                this->functionSpace(),
                this->functionSpace(),
                backendProjectorL2
                );
    }

    return M_projectorL2;
}

template<typename FunctionSpaceType>
typename LevelSetRedistanciationFM<FunctionSpaceType>::projector_ptrtype const& 
LevelSetRedistanciationFM<FunctionSpaceType>::projectorSM( bool buildOnTheFly ) const
{
    if( !M_projectorSM && buildOnTheFly )
    {
        double projectorSMCoeff = this->functionSpace()->mesh()->hAverage() / functionSpaceOrder * doption( _name="smooth-coeff", _prefix=prefixvm(this->prefix(),"projector-sm") );
        auto backendName = prefixvm( this->prefix(), "projector-sm" );
        auto backendProjectorSM = Backend<value_type>::build(
                soption( _prefix=backendName, _name="backend" ),
                backendName,
                this->functionSpace()->worldCommPtr()
                );
        const_cast<self_type*>(this)->M_projectorSM = projector(
                this->functionSpace(),
                this->functionSpace(),
                backendProjectorSM,
                DIFF, projectorSMCoeff, 30
                );
    }

    return M_projectorSM;
}

template<typename FunctionSpaceType>
void
LevelSetRedistanciationFM<FunctionSpaceType>::loadParametersFromOptionsVm()
{
    const std::string fm_init_method = soption( _name="fm-init-method", _prefix=this->prefix() );
    CHECK(FastMarchingInitialisationMethodMap.left.count(fm_init_method)) << fm_init_method <<" is not in the list of possible fast-marching initialisation methods\n";
    M_fastMarchingInitialisationMethod = FastMarchingInitialisationMethodMap.left.at(fm_init_method);
}

template<typename FunctionSpaceType>
typename LevelSetRedistanciationFM<FunctionSpaceType>::element_type 
LevelSetRedistanciationFM<FunctionSpaceType>::initFastMarching( element_type const& phi, range_elements_type const& rangeInitialElts ) const
{
    auto phiRedist = this->functionSpace()->elementPtr();
    switch( M_fastMarchingInitialisationMethod )
    {
        case FastMarchingInitialisationMethod::ILP_NODAL :
        {
            phiRedist->setConstant( 1e8 );

            auto const modGradPhi = this->functionSpaceDiscontinuous()->elementPtr();
            modGradPhi->on( _range=rangeInitialElts, _expr=norm2(gradv(phi)) );

            static const uint16_type nDofPerElt = functionspace_type::fe_type::nDof;
            auto itElt = boost::get<1>( rangeInitialElts );
            auto enElt = boost::get<2>( rangeInitialElts );
            for( ; itElt != enElt; ++itElt )
            {
                auto const elt = boost::unwrap_ref( *itElt );
                size_type const eltId = elt.id();

                for( uint16_type j = 0; j < nDofPerElt; ++j )
                {
                    size_type dofId = phiRedist->functionSpace()->dof()->localToGlobalId( eltId, j );
                    value_type modgradphi;
                    if constexpr( functionSpaceDiscontinuousOrder == 0 )
                        modgradphi = modGradPhi->localToGlobal( eltId, 0, 0 );
                    else
                        modgradphi = modGradPhi->localToGlobal( eltId, j, 0 );

                    value_type sdist = phi.localToGlobal( eltId, j, 0 ) / modgradphi;
                    if( std::abs( sdist ) < std::abs( (*phiRedist)(dofId) ) )
                        (*phiRedist)(dofId) = sdist;
                }
            }
        }
        break;

        case FastMarchingInitialisationMethod::ILP_L2 :
        {
            auto const modGradPhi = this->projectorL2()->project( norm2( gradv(phi) ) );
            phiRedist->on( 
                    _range=rangeInitialElts, 
                    _expr=idv(phi)/idv(modGradPhi)
                    );
        }
        break;

        case FastMarchingInitialisationMethod::ILP_SMOOTH :
        {
            auto const modGradPhi = this->projectorSM()->project( norm2( gradv(phi) ) );
            phiRedist->on( 
                    _range=rangeInitialElts, 
                    _expr=idv(phi)/idv(modGradPhi) 
                    );
        }
        break;

        case FastMarchingInitialisationMethod::ILDIST :
        {
            CHECK( convex_type::is_simplex ) << "ILDIST is only implemented for simplex" << std::endl;
            CHECK( functionSpaceOrder == 1 ) << "ILDIST is only implemented at order 1" << std::endl;

            phiRedist->setConstant( 1e8 );

            static const uint16_type nDofPerElt = functionspace_type::fe_type::nDof;
            auto itElt = boost::get<1>( rangeInitialElts );
            auto enElt = boost::get<2>( rangeInitialElts );
            for( ; itElt != enElt; ++itElt )
            {
                auto const elt = boost::unwrap_ref( *itElt );
                size_type const eltId = elt.id();

                std::vector< size_type > positiveDofIds, negativeDofIds;
                for( uint16_type j = 0; j < nDofPerElt; ++j )
                {
                    size_type dofId = phi.functionSpace()->dof()->localToGlobalId( eltId, j );
                    if( phi.localToGlobal( eltId, j, 0 ) < 0. )
                        negativeDofIds.push_back( dofId );
                    else
                        positiveDofIds.push_back( dofId );
                }
                CHECK( negativeDofIds.size() != 0 && positiveDofIds.size() != 0 ) << "ILPDIST only works with initial elements lying on the interface\n";

                if constexpr( nDim == 2)
                {
                    typedef Eigen::Matrix<value_type, nDim, 1, Eigen::ColMajor> pt_type;
                    std::vector< pt_type > segmentPts;
                    for( size_type const dofId1: negativeDofIds )
                    {
                        for( size_type const dofId2: positiveDofIds )
                        {
                            auto const& pt1 = boost::get<0>( phi.functionSpace()->dof()->dofPoint( dofId1 ) );
                            auto P1 = eigenMap<nDim>( pt1 );
                            value_type phi1 = phi(dofId1);

                            auto const& pt2 = boost::get<0>( phi.functionSpace()->dof()->dofPoint( dofId2 ) );
                            auto P2 = eigenMap<nDim>( pt2 );
                            value_type phi2 = phi(dofId2);

                            segmentPts.emplace_back( P1 + ( phi1 / (phi1-phi2) ) * (P2-P1) );
                        }
                    }
                    CHECK( segmentPts.size() == 2 ) << "should have 2 intersection points (have " << segmentPts.size() << ")" << std::endl;

                    for( uint16_type j = 0; j < nDofPerElt; ++j )
                    {
                        size_type dofId = phi.functionSpace()->dof()->localToGlobalId( eltId, j );
                        auto const& pt = boost::get<0>( phi.functionSpace()->dof()->dofPoint( dofId ) );
                        auto P = eigenMap<nDim>( pt );

                        value_type dist = Feel::detail::geometry::distancePointToSegment( P, segmentPts[0], segmentPts[1] );
                        if( dist < std::abs( (*phiRedist)(dofId) ) )
                            (*phiRedist)(dofId) = ( phi(dofId) > 0. ) ? dist : -dist;
                    }
                }
                else if constexpr( nDim == 3)
                {
                    //TODO
                }
            }
        }
        break;

        case FastMarchingInitialisationMethod::HJ_EQ :
        {
            CHECK(false) << "TODO\n";
            //*phi = *explicitHJ(max_iter, dtau, tol);
        }
        break;
        case FastMarchingInitialisationMethod::IL_HJ_EQ :
        {
            CHECK(false) << "TODO\n";
            //*phi = *explicitHJ(max_iter, dtau, tol);
        }
        break;
        case FastMarchingInitialisationMethod::NONE :
        {
            // Nothing to do.
        }
        break;
    } 

    return *phiRedist;
}

template<typename FunctionSpaceType>
typename LevelSetRedistanciationFM<FunctionSpaceType>::element_type 
LevelSetRedistanciationFM<FunctionSpaceType>::runFastMarching( element_type const& phi, range_elements_type const& rangeInitialElts ) const
{
    if constexpr( UseRedistP1Space )
    {
        auto phi_redist = this->functionSpace()->element();
        auto phi_FM = this->functionSpaceFM()->element();

        // Project onto P1 space
        if( functionSpaceOrder > 1 )
        {
            M_opInterpolationToP1->apply( phi, phi_FM );
        }
        else
        {
            phi_FM = vf::project( this->functionSpaceFM(), elements(this->mesh()), idv(phi) );
        }

        // Run fast-marching
        phi_FM = this->fastMarching()->run( phi_FM, rangeInitialElts );

        // Project back onto function-space
        if( functionSpaceOrder > 1 )
        {
            M_opInterpolationFromP1->apply( phi_FM, phi_redist );
        }
        else
        {
            phi_redist = vf::project( this->functionSpace(), elements(this->mesh()), idv(phi_FM) );
        }

        return phi_redist;
    }
    else
    {
        // Directly run fast-marching
        return this->fastMarching()->run( phi, rangeInitialElts );
    }
}

template<typename FunctionSpaceType>
typename LevelSetRedistanciationFM<FunctionSpaceType>::element_type 
LevelSetRedistanciationFM<FunctionSpaceType>::run( element_type const& phi, range_elements_type const& rangeInitialElts ) const
{
    return this->runFastMarching( 
            this->initFastMarching( phi, rangeInitialElts ), 
            rangeInitialElts 
            );
}

template<typename FunctionSpaceType>
typename LevelSetRedistanciationFM<FunctionSpaceType>::element_type 
LevelSetRedistanciationFM<FunctionSpaceType>::run( element_type const& phi ) const
{
    auto rangeInitialElements = levelsetInterfaceElements( phi );
    return this->run( phi, rangeInitialElements );
}

} // namespace Feel

#endif //_LEVELSET_REDISTANCIATION_FM_HPP
