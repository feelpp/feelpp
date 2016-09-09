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
#include <feel/feelmodels/levelset/spaceP1wrapper.hpp>
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
    ReinitializerFM( 
            functionspace_ptrtype const& space,
            std::string const& prefix = "" );
    //--------------------------------------------------------------------//
    // Run reinitialization
    void run( element_type const& phi );
    //--------------------------------------------------------------------//
    // Accessors
    element_type distance() const;
    //--------------------------------------------------------------------//
    // Parameters
    void loadParametersFromOptionsVm();

    void setUseMarker2AsMarkerDone( bool val = true ) { M_useMarker2AsMarkerDone = val; }
    bool useMarker2AsMarkerDone() const { return M_useMarker2AsMarkerDone; }

private:
    //--------------------------------------------------------------------//
    // Space P1 wrapper
    typedef SpaceP1Wrapper<functionspace_type> spaceP1wrapper_type;
    typedef boost::shared_ptr<spaceP1wrapper_type> spaceP1wrapper_ptrtype;

    spaceP1wrapper_ptrtype M_spaceP1Wrapper;
    //--------------------------------------------------------------------//
    // ReinitializerFMS (takes only P1 non-periodic space)
    typedef typename mpl::if_< 
        mpl::bool_<spaceP1wrapper_type::use_P1_space>,
        typename spaceP1wrapper_type::functionspace_P1_type,
        functionspace_type
            >::type reinitializerFMS_functionspace_type;
    
    typedef ReinitializerFMS< reinitializerFMS_functionspace_type, periodicity_type > reinitializerFMS_type;
    typedef boost::shared_ptr<reinitializerFMS_type> reinitializerFMS_ptrtype;

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
    M_spaceP1Wrapper.reset(
            new spaceP1wrapper_type( space )
            );
    M_reinitializerFMS.reset(
            new reinitializerFMS_type( M_spaceP1Wrapper->functionSpaceP1Wrapped(), this->M_periodicity )
            );
}

template<typename FunctionSpaceType>
void
ReinitializerFM<FunctionSpaceType>::run( element_type const& phi )
{
    M_reinitializerFMS->run( M_spaceP1Wrapper->toP1(phi), this->useMarker2AsMarkerDone() ); 
}

template<typename FunctionSpaceType>
typename ReinitializerFM<FunctionSpaceType>::element_type
ReinitializerFM<FunctionSpaceType>::distance() const
{
    return M_spaceP1Wrapper->fromP1( *(M_reinitializerFMS->getDistance()) );
}

template<typename FunctionSpaceType>
void
ReinitializerFM<FunctionSpaceType>::loadParametersFromOptionsVm()
{
    M_useMarker2AsMarkerDone = boption( _name="use-marker2-as-done", _prefix=this->prefix() );
}

} // namespace Feel

#endif
