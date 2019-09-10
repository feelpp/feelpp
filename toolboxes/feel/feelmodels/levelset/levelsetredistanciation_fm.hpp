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
        static const uint16_type functionSpaceOrder = functionspace_type::fe_type::nOrder;
        typedef typename functionspace_type::value_type value_type;
        // Element
        typedef typename functionspace_type::element_type element_type;
        typedef std::shared_ptr<element_type> element_ptrtype;
        // Mesh
        typedef typename functionspace_type::mesh_type mesh_type;
        // Periodicity
        typedef typename functionspace_type::periodicity_0_type periodicity_type;
        static const bool functionSpaceIsPeriodic = functionspace_type::is_periodic;

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

        //--------------------------------------------------------------------//
        //--------------------------------------------------------------------//
        //--------------------------------------------------------------------//
        // Constructor
        LevelSetRedistanciationFM( 
                functionspace_ptrtype const& space,
                std::string const& prefix = "" );

        // Fast-marching space
        functionspace_FM_ptrtype const& functionSpaceFM() const { return M_spaceFM; }

        // Fast-marching
        fastmarching_ptrtype const& fastMarching() const { return M_fastMarching; }

        // Redistanciation
        element_type run( element_type const& phi, range_elements_type const& rangeInitialElts ) const;
        element_type run( element_type const& phi ) const;

        // Parameters
        void loadParametersFromOptionsVm();

    private:
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

        //--------------------------------------------------------------------//
        // Fast-marching
        fastmarching_ptrtype M_fastMarching;
        //--------------------------------------------------------------------//
        //--------------------------------------------------------------------//
        // Options
};

template<typename FunctionSpaceType>
LevelSetRedistanciationFM<FunctionSpaceType>::LevelSetRedistanciationFM( 
        functionspace_ptrtype const& space,
        std::string const& prefix )
    : super_type( space, prefix )
{
    // Load parameters
    this->loadParametersFromOptionsVm();
    // Init
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
void
LevelSetRedistanciationFM<FunctionSpaceType>::loadParametersFromOptionsVm()
{
}

template<typename FunctionSpaceType>
typename LevelSetRedistanciationFM<FunctionSpaceType>::element_type 
LevelSetRedistanciationFM<FunctionSpaceType>::run( element_type const& phi, range_elements_type const& rangeInitialElts ) const
{
    if constexpr( UseRedistP1Space )
    {
        auto phi_reinit = this->functionSpace()->element();
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
            M_opInterpolationFromP1->apply( phi_FM, phi_reinit );
        }
        else
        {
            phi_reinit = vf::project( this->functionSpace(), elements(this->mesh()), idv(phi_FM) );
        }

        return phi_reinit;
    }
    else
    {
        // Directly run fast-marching
        return this->fastMarching()->run( phi, rangeInitialElts );
    }
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
