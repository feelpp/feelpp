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
 \file distlabels_fm.hpp
 \author Thibaut Metivet <thibaut.metivet@univ-grenoble-alpes.fr>
 \date 2016-09-08
 */
#ifndef _DISTLABELS_FM_HPP
#define _DISTLABELS_FM_HPP 1

#include <feel/feelmodels/levelset/reinitializer.hpp>
#include <feel/feells/distlabels_fms.hpp>

namespace Feel
{

template<typename FunctionSpaceType>
class DistLabelsFM
: public Reinitializer<FunctionSpaceType>
{
public:
    typedef Reinitializer<FunctionSpaceType> super_type;
    typedef DistLabelsFM<FunctionSpaceType> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;
    //--------------------------------------------------------------------//
    // Function spaces
    typedef FunctionSpaceType functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    static const uint16_type nOrder = functionspace_type::fe_type::nOrder;
    typedef typename functionspace_type::mesh_type mesh_type;
    typedef typename functionspace_type::value_type value_type;

    typedef typename functionspace_type::periodicity_0_type periodicity_type;
    static const bool is_periodic = functionspace_type::is_periodic;

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Constructor
    DistLabelsFM( 
            functionspace_ptrtype const& space,
            std::string const& prefix = "" );
    //--------------------------------------------------------------------//
    // Run reinitialization
    self_type& run( element_type const& phi );
    //--------------------------------------------------------------------//
    // Accessors
    element_type distance() const;
    element_type nearestNeighbourLabel() const;
    element_type nextNearestNeighbourDistance() const;
    element_type nextNearestNeighbourLabel() const;

    void setSelfLabel( element_ptrtype const& label ); 
    //--------------------------------------------------------------------//
    // Parameters
    void setUseMarker2AsMarkerDone( bool val = true ) { M_useMarker2AsMarkerDone = val; }
    bool useMarker2AsMarkerDone() const { return M_useMarker2AsMarkerDone; }

private:
    //--------------------------------------------------------------------//
    // Space P1 wrapper
    typedef SpaceP1Wrapper<functionspace_type> spaceP1wrapper_type;
    typedef boost::shared_ptr<spaceP1wrapper_type> spaceP1wrapper_ptrtype;

    spaceP1wrapper_ptrtype M_spaceP1Wrapper;

    //--------------------------------------------------------------------//
    // LabelDistanceFMS (takes only P1 non-periodic space)
    typedef typename mpl::if_< 
        mpl::bool_<spaceP1wrapper_type::use_P1_space>,
        typename spaceP1wrapper_type::functionspace_P1_type,
        functionspace_type
            >::type labeldistanceFMS_functionspace_type;
    //typedef typename ReinitFunctionSpace<functionspace_type>::type labeldistanceFMS_functionspace_type;
    
    typedef LabelDistanceFMS< labeldistanceFMS_functionspace_type, periodicity_type > labeldistanceFMS_type;
    typedef boost::shared_ptr<labeldistanceFMS_type> labeldistanceFMS_ptrtype;

    labeldistanceFMS_ptrtype M_labeldistanceFMS;
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    // Options
    bool M_useMarker2AsMarkerDone;
};

template<typename FunctionSpaceType>
DistLabelsFM<FunctionSpaceType>::DistLabelsFM( 
        functionspace_ptrtype const& space,
        std::string const& prefix ) :
    super_type( space, prefix ),
    M_useMarker2AsMarkerDone( false )
{
    M_spaceP1Wrapper.reset(
            new spaceP1wrapper_type( space )
            );
    M_labeldistanceFMS.reset(
            new labeldistanceFMS_type( M_spaceP1Wrapper->functionSpaceP1Wrapped(), this->M_periodicity )
            );
}

template<typename FunctionSpaceType>
typename DistLabelsFM<FunctionSpaceType>::self_type&
DistLabelsFM<FunctionSpaceType>::run( element_type const& phi )
{
    M_labeldistanceFMS->run( M_spaceP1Wrapper->toP1(phi), this->useMarker2AsMarkerDone() );
    return *this;
}

template<typename FunctionSpaceType>
typename DistLabelsFM<FunctionSpaceType>::element_type
DistLabelsFM<FunctionSpaceType>::distance() const
{
    return M_spaceP1Wrapper->fromP1( *(M_labeldistanceFMS->getNearestNeighbourDistance()) );
}

template<typename FunctionSpaceType>
typename DistLabelsFM<FunctionSpaceType>::element_type
DistLabelsFM<FunctionSpaceType>::nearestNeighbourLabel() const
{
    return M_spaceP1Wrapper->fromP1( *(M_labeldistanceFMS->getNearestNeighbourLabel()) );
}

template<typename FunctionSpaceType>
typename DistLabelsFM<FunctionSpaceType>::element_type
DistLabelsFM<FunctionSpaceType>::nextNearestNeighbourDistance() const
{
    return M_spaceP1Wrapper->fromP1( *(M_labeldistanceFMS->getNextNearestNeighbourDistance()) );
}

template<typename FunctionSpaceType>
typename DistLabelsFM<FunctionSpaceType>::element_type
DistLabelsFM<FunctionSpaceType>::nextNearestNeighbourLabel() const
{
    return M_spaceP1Wrapper->fromP1( *(M_labeldistanceFMS->getNextNearestNeighbourLabel()) );
}

template<typename FunctionSpaceType>
void
DistLabelsFM<FunctionSpaceType>::setSelfLabel( element_ptrtype const& label )
{ 
    M_labeldistanceFMS->setSelfLabel( M_spaceP1Wrapper->toP1(label) );
}

} // namespace Feel

#endif
