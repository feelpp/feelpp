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
        static inline const uint16_type functionSpaceOrder = functionspace_type::fe_type::nOrder;
        static inline const uint16_type nOrder = functionspace_type::fe_type::nOrder;
        static constexpr uint16_type nDim = functionspace_type::nDim;
        using size_type = typename functionspace_type::size_type;
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

        // Marker P0d space
        typedef Lagrange<0, Scalar, Discontinuous> basis_P0d_type;
        typedef FunctionSpace<
            mesh_type, 
            bases<basis_P0d_type>, 
            value_type,
            Periodicity<NoPeriodicity>,
            mortars<NoMortar>
                > functionspace_P0d_type;
        typedef std::shared_ptr<functionspace_P0d_type> functionspace_P0d_ptrtype;

        // Gradient discontinuous space
        static inline const uint16_type functionSpaceDiscontinuousOrder = ( functionSpaceOrder == 1 ) ? 0 : functionSpaceOrder;
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
        // Fast-marching space
        functionspace_FM_ptrtype M_spaceFM;
        // Marker P0d spaces
        functionspace_P0d_ptrtype M_spaceP0d;
        functionspace_P0d_ptrtype M_spaceP0dIsoPN;

        //--------------------------------------------------------------------//
        // Fast-marching
        fastmarching_ptrtype M_fastMarching;
        FastMarchingInitialisationMethod M_fastMarchingInitialisationMethod;
        //--------------------------------------------------------------------//
        //--------------------------------------------------------------------//
        // Options
};

} // namespace Feel

#endif //_LEVELSET_REDISTANCIATION_FM_HPP
