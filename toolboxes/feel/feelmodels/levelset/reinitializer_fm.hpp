/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 Date: 2016-05-20

 Copyright (C) 2016 Université Joseph Fourier (Grenoble I)

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
 \file reinitializer_fm.hpp
 \author Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 \date 2016-05-20
 */
#ifndef _REINITIALIZER_FM_HPP
#define _REINITIALIZER_FM_HPP 1

#include <feel/feelmodels/levelset/reinitializer.hpp>
#include <feel/feells/reinit_fms_impl.hpp>

namespace Feel
{

template<typename FunctionSpaceType>
class ReinitializerFM
: public Reinitializer<FunctionSpaceType>
{
public:
    typedef Reinitializer<FunctionSpaceType> super_type;
    typedef ReinitializerFM<FunctionSpaceType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;
    //--------------------------------------------------------------------//
    // Function spaces
    typedef FunctionSpaceType functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;

    static inline const uint16_type nOrder = functionspace_type::fe_type::nOrder;
    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename functionspace_type::value_type value_type;

    typedef typename functionspace_type::periodicity_0_type periodicity_type;
    static const bool is_periodic = functionspace_type::is_periodic;

    template< template<uint16_type Dim> class PointSetT >
    struct BasisReinitP1
    {
        typedef Lagrange<1, PointSetT> type;
    };
    typedef typename 
        mpl::if_<mpl::bool_<functionspace_type::is_scalar>,
            mpl::identity<typename BasisReinitP1<Scalar>::type>,
            typename mpl::if_<mpl::bool_<functionspace_type::is_vectorial>,
                mpl::identity<typename BasisReinitP1<Vectorial>::type>,
                mpl::identity<typename BasisReinitP1<Tensor2>::type> 
            >::type
        >::type::type basis_reinitP1_type;
    typedef FunctionSpace<
        mesh_type, 
        bases<basis_reinitP1_type>, 
        value_type,
        Periodicity<NoPeriodicity>,
        mortars<NoMortar>
            > functionspace_reinitP1_type;
    typedef std::shared_ptr<functionspace_reinitP1_type> functionspace_reinitP1_ptrtype;

    //--------------------------------------------------------------------//
    // Range
    typedef Range<mesh_type,MESH_ELEMENTS> range_elements_type;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Constructor
    ReinitializerFM( 
            functionspace_ptrtype const& space,
            std::string const& prefix = "" );
    //--------------------------------------------------------------------//
    // Run reinitialization
    element_type run( element_type const& phi, range_elements_type const& rangeInitialElts );
    element_type run( element_type const& phi );
    //--------------------------------------------------------------------//
    // Parameters
    void loadParametersFromOptionsVm();

    void setUseMarker2AsMarkerDone( bool val = true ) { M_useMarker2AsMarkerDone = val; }
    bool useMarker2AsMarkerDone() const { return M_useMarker2AsMarkerDone; }

private:
    template<typename FST>
    struct UseReinitP1Space
    {
        static constexpr bool value = !( FST::fe_type::nOrder == 1 && !FST::is_periodic );
    };
    static constexpr bool use_reinitP1_space = UseReinitP1Space<functionspace_type>::value;

    //template<typename FST>
    //struct ReinitFunctionSpace
    //{
        //typedef typename mpl::if_< 
            //use_reinitP1_space,
            //functionspace_reinitP1_type,
            //functionspace_type
                //>::type type;
    //};

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    template<typename FST>
    void init( functionspace_ptrtype const& space, 
            typename std::enable_if< UseReinitP1Space<FST>::value >::type* =0);
    template<typename FST>
    void init( functionspace_ptrtype const& space, 
            typename std::enable_if< !UseReinitP1Space<FST>::value >::type* =0);
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    template<typename FST>
    element_type run_impl( element_type const& phi, range_elements_type const& rangeInitialElts,
            typename std::enable_if< UseReinitP1Space<FST>::value >::type* =0);
    template<typename FST>
    element_type run_impl( element_type const& phi, 
            typename std::enable_if< UseReinitP1Space<FST>::value >::type* =0);
    template<typename FST>
    element_type run_impl( element_type const& phi, range_elements_type const& rangeInitialElts,
            typename std::enable_if< !UseReinitP1Space<FST>::value >::type* =0);
    template<typename FST>
    element_type run_impl( element_type const& phi, 
            typename std::enable_if< !UseReinitP1Space<FST>::value >::type* =0);

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Reinit P1 operators
    typedef OperatorInterpolation<
        functionspace_type, // from space
        functionspace_reinitP1_type // to space
        > op_interpolation_to_P1_type;

    typedef OperatorInterpolation<
        functionspace_reinitP1_type, // from space
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
    // Reinit P1 space
    functionspace_reinitP1_ptrtype M_spaceReinitP1;

    //--------------------------------------------------------------------//
    // ReinitializerFMS (takes only P1 non-periodic space)
    typedef typename mpl::if_< 
        mpl::bool_<use_reinitP1_space>,
        functionspace_reinitP1_type,
        functionspace_type
            >::type reinitializerFMS_functionspace_type;
    //typedef typename ReinitFunctionSpace<functionspace_type>::type reinitializerFMS_functionspace_type;
    
    typedef ReinitializerFMS< reinitializerFMS_functionspace_type, periodicity_type > reinitializerFMS_type;
    typedef std::shared_ptr<reinitializerFMS_type> reinitializerFMS_ptrtype;

    reinitializerFMS_ptrtype M_reinitializerFMS;
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Options
    bool M_useMarker2AsMarkerDone;
};

template<typename FunctionSpaceType>
ReinitializerFM<FunctionSpaceType>::ReinitializerFM( 
        functionspace_ptrtype const& space,
        std::string const& prefix )
    : super_type( space, prefix )
{
    this->loadParametersFromOptionsVm();
    this->init<FunctionSpaceType>( space );
}

template<typename FunctionSpaceType>
typename ReinitializerFM<FunctionSpaceType>::element_type 
ReinitializerFM<FunctionSpaceType>::run( element_type const& phi, range_elements_type const& rangeInitialElts )
{
    return run_impl<FunctionSpaceType>( phi, rangeInitialElts );
}

template<typename FunctionSpaceType>
typename ReinitializerFM<FunctionSpaceType>::element_type 
ReinitializerFM<FunctionSpaceType>::run( element_type const& phi )
{
    return run_impl<FunctionSpaceType>( phi );
}

template<typename FunctionSpaceType>
void
ReinitializerFM<FunctionSpaceType>::loadParametersFromOptionsVm()
{
    M_useMarker2AsMarkerDone = boption( _name="use-marker2-as-done", _prefix=this->prefix() );
}

// Init if using own reinit P1 space
template<typename FunctionSpaceType>
template<typename FST>
void 
ReinitializerFM<FunctionSpaceType>::init(functionspace_ptrtype const& space, 
        typename std::enable_if< UseReinitP1Space<FST>::value >::type*)
{
    M_opLagrangeP1 = lagrangeP1( _space=space, _update=MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES );
    M_spaceReinitP1 = functionspace_reinitP1_type::New(
            _mesh=M_opLagrangeP1->mesh(),
            _periodicity=periodicity(NoPeriodicity())
            );
    M_reinitializerFMS.reset(
            new reinitializerFMS_type( M_spaceReinitP1, this->M_periodicity )
            );

    if( nOrder > 1 )
    {
        M_opInterpolationToP1 = opInterpolation(
                _domainSpace = space,
                _imageSpace = M_spaceReinitP1,
                _type = InterpolationNonConforme(false)
                );
        M_opInterpolationFromP1 = opInterpolation(
                _domainSpace = M_spaceReinitP1,
                _imageSpace = space,
                _type = InterpolationNonConforme(false)
                );
    }
}

// Init if using provided space
template<typename FunctionSpaceType>
template<typename FST>
void 
ReinitializerFM<FunctionSpaceType>::init(functionspace_ptrtype const& space, 
        typename std::enable_if< !UseReinitP1Space<FST>::value >::type*)
{
    M_reinitializerFMS.reset(
            new reinitializerFMS_type( space, this->M_periodicity )
            );
}

// Run if using own reinit P1 space
template<typename FunctionSpaceType>
template<typename FST>
typename ReinitializerFM<FunctionSpaceType>::element_type 
ReinitializerFM<FunctionSpaceType>::run_impl( element_type const& phi, range_elements_type const& rangeInitialElts,
        typename std::enable_if< UseReinitP1Space<FST>::value >::type*)
{
    auto phi_reinit = this->functionSpace()->element();
    auto phi_reinitP1 = M_spaceReinitP1->element();

    if( nOrder > 1 )
    {
        M_opInterpolationToP1->apply( phi, phi_reinitP1 );
    }
    else
    {
        phi_reinitP1 = vf::project( M_spaceReinitP1, elements(this->mesh()), idv(phi) );
    }

    //phi_reinitP1 = M_reinitializerFMS->march( phi_reinitP1, this->useMarker2AsMarkerDone() );
    phi_reinitP1 = M_reinitializerFMS->march( phi_reinitP1 );
    // TODO: provide possible rangeInitialElts

    if( nOrder > 1 )
    {
        M_opInterpolationFromP1->apply( phi_reinitP1, phi_reinit );
    }
    else
    {
        phi_reinit = vf::project( this->functionSpace(), elements(this->mesh()), idv(phi_reinitP1) );
    }
    
    return phi_reinit;
}
template<typename FunctionSpaceType>
template<typename FST>
typename ReinitializerFM<FunctionSpaceType>::element_type 
ReinitializerFM<FunctionSpaceType>::run_impl( element_type const& phi, 
        typename std::enable_if< UseReinitP1Space<FST>::value >::type*)
{
    auto phi_reinit = this->functionSpace()->element();
    auto phi_reinitP1 = M_spaceReinitP1->element();

    if( nOrder > 1 )
    {
        M_opInterpolationToP1->apply( phi, phi_reinitP1 );
    }
    else
    {
        phi_reinitP1 = vf::project( _space=M_spaceReinitP1, _range=elements(this->mesh()), _expr=idv(phi) );
    }

    phi_reinitP1 = M_reinitializerFMS->march( phi_reinitP1 );

    if( nOrder > 1 )
    {
        M_opInterpolationFromP1->apply( phi_reinitP1, phi_reinit );
    }
    else
    {
        phi_reinit = vf::project( _space=this->functionSpace(), _range=elements(this->mesh()), _expr=idv(phi_reinitP1) );
    }
    
    return phi_reinit;
}

// Run if using provided space
template<typename FunctionSpaceType>
template<typename FST>
typename ReinitializerFM<FunctionSpaceType>::element_type 
ReinitializerFM<FunctionSpaceType>::run_impl( element_type const& phi, range_elements_type const& rangeInitialElts, 
        typename std::enable_if< !UseReinitP1Space<FST>::value >::type*)
{
    return M_reinitializerFMS->march( phi, rangeInitialElts );
}
template<typename FunctionSpaceType>
template<typename FST>
typename ReinitializerFM<FunctionSpaceType>::element_type 
ReinitializerFM<FunctionSpaceType>::run_impl( element_type const& phi, 
        typename std::enable_if< !UseReinitP1Space<FST>::value >::type*)
{
    return M_reinitializerFMS->march( phi );
}

} // namespace Feel

#endif
