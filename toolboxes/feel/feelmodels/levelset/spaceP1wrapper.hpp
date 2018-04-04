/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 Date: 2016-09-08

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
 \file spaceP1wrapper.hpp
 \author Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 \date 2016-09-08
 */
#ifndef _SPACEP1WRAPPER_HPP
#define _SPACEP1WRAPPER_HPP 1

#include <feel/feeldiscr/operatorlagrangep1.hpp>

namespace Feel
{

template<typename FunctionSpaceType>
class SpaceP1Wrapper
{
public:
    typedef SpaceP1Wrapper<FunctionSpaceType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
    //--------------------------------------------------------------------//
    // Function spaces
    typedef FunctionSpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    static const uint16_type nOrder = functionspace_type::fe_type::nOrder;
    typedef typename functionspace_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

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
            > functionspace_P1_type;
    typedef boost::shared_ptr<functionspace_P1_type> functionspace_P1_ptrtype;

    typedef typename functionspace_P1_type::element_type element_P1_type;
    typedef boost::shared_ptr<element_P1_type> element_P1_ptrtype;

    template<typename FST>
    struct UseP1Space
    {
        static constexpr bool value = !( FST::fe_type::nOrder == 1 && !FST::is_periodic );
    };
    static constexpr bool use_P1_space = UseP1Space<functionspace_type>::value;

    typedef typename mpl::if_< 
        mpl::bool_<use_P1_space>,
        functionspace_P1_type,
        functionspace_type
            >::type functionspace_P1wrapped_type;
    typedef typename mpl::if_< 
        mpl::bool_<use_P1_space>,
        functionspace_P1_ptrtype,
        functionspace_ptrtype
            >::type functionspace_P1wrapped_ptrtype;

    typedef typename mpl::if_< 
        mpl::bool_<use_P1_space>,
        element_P1_type,
        element_type
            >::type element_P1wrapped_type;
    typedef typename mpl::if_< 
        mpl::bool_<use_P1_space>,
        element_P1_ptrtype,
        element_ptrtype
            >::type element_P1wrapped_ptrtype;

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

    typedef boost::shared_ptr<op_interpolation_to_P1_type> op_interpolation_to_P1_ptrtype;
    typedef boost::shared_ptr<op_interpolation_from_P1_type> op_interpolation_from_P1_ptrtype;

    typedef OperatorLagrangeP1<functionspace_type> op_lagrangeP1_type;
    typedef boost::shared_ptr<op_lagrangeP1_type> op_lagrangeP1_ptrtype;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Constructor
    SpaceP1Wrapper( functionspace_ptrtype const& space );

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Function spaces
    functionspace_ptrtype const& functionSpace() const { return M_space; }
    mesh_ptrtype const& mesh() const { return M_space->mesh(); }

    functionspace_P1wrapped_ptrtype const& functionSpaceP1Wrapped() { return functionSpaceP1Wrapped_impl<functionspace_type>(); }

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // toP1
    element_P1wrapped_type toP1( element_type const& phi ) { return toP1_impl<functionspace_type>( phi ); }
    element_P1wrapped_ptrtype toP1( element_ptrtype const& phi ) 
    { 
        auto phiP1 = this->functionSpaceP1Wrapped()->elementPtr();
        *phiP1 = toP1_impl<functionspace_type>( *phi );
        return phiP1;
    }
    // fromP1
    element_type fromP1( element_P1wrapped_type const& phi ) { return fromP1_impl<functionspace_type>( phi ); }
    element_ptrtype fromP1( element_P1wrapped_ptrtype const& phi ) 
    {
        auto phiP1 = this->functionSpace()->elementPtr();
        *phiP1 = fromP1_impl<functionspace_type>( *phi );
        return phiP1;
    }

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
private:
    template<typename FST>
    void init( functionspace_ptrtype const& space, 
            typename std::enable_if< UseP1Space<FST>::value >::type* =0);
    template<typename FST>
    void init( functionspace_ptrtype const& space, 
            typename std::enable_if< !UseP1Space<FST>::value >::type* =0);

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Function space P1 wrapped
    template<typename FST>
    functionspace_P1_ptrtype const& functionSpaceP1Wrapped_impl( 
            typename std::enable_if< UseP1Space<FST>::value >::type* =0)
    {
        return M_spaceP1;
    }
    template<typename FST>
    functionspace_ptrtype const& functionSpaceP1Wrapped_impl( 
            typename std::enable_if< !UseP1Space<FST>::value >::type* =0)
    {
        return M_space;
    }
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // toP1
    template<typename FST>
    element_P1_type toP1_impl( element_type const& phi, 
            typename std::enable_if< UseP1Space<FST>::value >::type* =0);
    template<typename FST>
    element_type toP1_impl( element_type const& phi, 
            typename std::enable_if< !UseP1Space<FST>::value >::type* =0);
    // fromP1
    template<typename FST>
    element_type fromP1_impl( element_P1_type const& phi, 
            typename std::enable_if< UseP1Space<FST>::value >::type* =0);
    template<typename FST>
    element_type fromP1_impl( element_type const & phi, 
            typename std::enable_if< !UseP1Space<FST>::value >::type* =0);

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    op_interpolation_to_P1_ptrtype M_opInterpolationToP1;
    op_interpolation_from_P1_ptrtype M_opInterpolationFromP1;
    op_lagrangeP1_ptrtype M_opLagrangeP1;
    //--------------------------------------------------------------------//
    // Function spaces
    functionspace_ptrtype M_space;
    periodicity_type M_periodicity;
    functionspace_P1_ptrtype M_spaceP1;
};

template<typename FunctionSpaceType>
SpaceP1Wrapper<FunctionSpaceType>::SpaceP1Wrapper( 
        functionspace_ptrtype const& space ) : 
    M_space(space),
    M_periodicity(boost::fusion::at_c<0>(space->periodicity()))
{
    this->init<FunctionSpaceType>( space );
}

// Init if using own reinit P1 space
template<typename FunctionSpaceType>
template<typename FST>
void 
SpaceP1Wrapper<FunctionSpaceType>::init(functionspace_ptrtype const& space, 
        typename std::enable_if< UseP1Space<FST>::value >::type*)
{
    M_opLagrangeP1 = lagrangeP1( space );
    M_spaceP1 = functionspace_P1_type::New(
            _mesh=M_opLagrangeP1->mesh(),
            _periodicity=periodicity(NoPeriodicity())
            );

    if( nOrder > 1 )
    {
        M_opInterpolationToP1 = opInterpolation(
                _domainSpace = space,
                _imageSpace = M_spaceP1,
                _type = InterpolationNonConforme(false)
                );
        M_opInterpolationFromP1 = opInterpolation(
                _domainSpace = M_spaceP1,
                _imageSpace = space,
                _type = InterpolationNonConforme(false)
                );
    }
}

// Init if using provided space
template<typename FunctionSpaceType>
template<typename FST>
void 
SpaceP1Wrapper<FunctionSpaceType>::init(functionspace_ptrtype const& space, 
        typename std::enable_if< !UseP1Space<FST>::value >::type*)
{
}

// toP1 if using own reinit P1 space
template<typename FunctionSpaceType>
template<typename FST>
typename SpaceP1Wrapper<FunctionSpaceType>::element_P1_type 
SpaceP1Wrapper<FunctionSpaceType>::toP1_impl( element_type const& phi, 
        typename std::enable_if< UseP1Space<FST>::value >::type*)
{

    if( nOrder > 1 )
    {
        auto phi_P1 = M_spaceP1->element();
        M_opInterpolationToP1->apply( phi, phi_P1 );
        return phi_P1;
    }
    else
    {
        return vf::project( M_spaceP1, elements(this->mesh()), idv(phi) );
    }

}

// toP1 if using provided space
template<typename FunctionSpaceType>
template<typename FST>
typename SpaceP1Wrapper<FunctionSpaceType>::element_type 
SpaceP1Wrapper<FunctionSpaceType>::toP1_impl( element_type const& phi, 
        typename std::enable_if< !UseP1Space<FST>::value >::type*)
{
    return phi;
}

// fromP1 if using own reinit P1 space
template<typename FunctionSpaceType>
template<typename FST>
typename SpaceP1Wrapper<FunctionSpaceType>::element_type 
SpaceP1Wrapper<FunctionSpaceType>::fromP1_impl( element_P1_type const& phi_P1, 
        typename std::enable_if< UseP1Space<FST>::value >::type*)
{
    if( nOrder > 1 )
    {
        auto phi = this->functionSpace()->element();
        M_opInterpolationFromP1->apply( phi_P1, phi );
        return phi;
    }
    else
    {
        return vf::project( this->functionSpace(), elements(this->mesh()), idv(phi_P1) );
    }
}

// fromP1 if using provided space
template<typename FunctionSpaceType>
template<typename FST>
typename SpaceP1Wrapper<FunctionSpaceType>::element_type 
SpaceP1Wrapper<FunctionSpaceType>::fromP1_impl( element_type const& phi, 
        typename std::enable_if< !UseP1Space<FST>::value >::type*)
{
    return phi;
}

} // namespace Feel

#endif
