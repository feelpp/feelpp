/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Copyright (C) 2010 University of Coimbra

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
 \file ale.cpp
 \author Goncalo Pena <gpena@mat.uc.pt>
 \date 2010-10-12
 */

#include <feel/feelmodels/modelmesh/ale_impl.hpp>

#include <boost/preprocessor/comparison/greater_equal.hpp>

#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/unary.hpp>
#include <feel/feelvf/one.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/trans.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/trace.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/val.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/on.hpp>

#include <feel/feelmodels/modelmesh/metricmeshadaptation.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>

namespace Feel
{
namespace FeelModels
{
namespace ALE_IMPL
{
#if 0
namespace detailALE
{
template < typename SpaceLowType,typename SpaceHighType >
std::shared_ptr<SpaceHighType>
buildSpaceHigh(std::shared_ptr<SpaceLowType> spaceLow, mpl::bool_<false> /**/ )
{
    return SpaceHighType::New( _mesh=spaceLow->mesh() );
}

template < typename SpaceLowType,typename SpaceHighType >
std::shared_ptr<SpaceHighType>
buildSpaceHigh(std::shared_ptr<SpaceLowType> spaceLow, mpl::bool_<true> /**/ )
{
    return spaceLow;
}
}
#endif

template < class Convex, int Order >
ALE<Convex,Order>::ALE( std::string const& prefix, worldcomm_ptr_t const& worldcomm,
                        ModelBaseRepository const& modelRep )
    :
    super_type( prefix,worldcomm,modelRep ),
    M_verboseSolverTimer(boption(_prefix=this->prefix(),_name="verbose_solvertimer")),
    M_verboseSolverTimerAllProc(boption(_prefix=this->prefix(),_name="verbose_solvertimer_allproc")),
    M_alemeshTypeName( soption( _name="type",_prefix=this->prefix() ) ),
    M_doHoCorrection( boption(_prefix=this->prefix(),_name="apply-ho-correction") ),
    M_isInitHarmonicExtension( false ),
    M_isInitWinslow( false )
{}

template < class Convex, int Order >
ALE<Convex,Order>::ALE( mesh_ptrtype mesh, std::string const& prefix, worldcomm_ptr_t const& worldcomm,
                        ModelBaseRepository const& modelRep )
    :
    ALE( prefix,worldcomm,modelRep )
{
    this->log( "ALE", "constructor 1", "start" );

    M_reference_mesh = mesh;

    this->createALE();

    this->preCompute();

    this->log( "ALE", "constructor 1", "finish" );
}
template < class Convex, int Order >
ALE<Convex,Order>::ALE( mesh_ptrtype mesh, range_elements_type const& rangeElt,
                        std::string const& prefix, worldcomm_ptr_t const& worldcomm,
                        ModelBaseRepository const& modelRep )
    :
    ALE( prefix,worldcomm,modelRep )
{
    this->log( "ALE", "constructor 2", "start" );

    M_reference_mesh = mesh;

    this->createALE( rangeElt );

    this->preCompute();

    this->log( "ALE", "constructor 2", "finish" );
}

//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
std::shared_ptr<std::ostringstream>
ALE<Convex,Order>::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n   Physical Markers";
    for ( std::string bctype : { "moving","fixed","free" } )
    {
        if ( this->bcToMarkers().find(bctype) != this->bcToMarkers().end() )
        {
            *_ostr << "\n     -- " << bctype << " : ";
            auto it = this->bcToMarkers().find(bctype)->second.begin();
            auto en = this->bcToMarkers().find(bctype)->second.end();
            for ( int cptMark = 0 ; it!=en ; ++it,++cptMark )
            {
                if ( cptMark > 0 ) *_ostr << " , ";
                *_ostr << *it;
            }
        }
    }
        // if ( this->flagSet().find("moving") != this->flagSet().end() )
    // {
    //     *_ostr << "\n     -- moving : ";
    //     auto it = this->flagSet().find("moving")->second.begin();
    //     auto en = this->flagSet().find("moving")->second.end();
    //     for ( int cptMark = 0 ; it!=en ; ++it,++cptMark )
    //     {
    //         if ( cptMark > 0 ) *_ostr << " , ";
    //         *_ostr << *it;
    //     }
    // }
    // if ( this->flagSet().find("free") != this->flagSet().end() )
    // {
    //     *_ostr << "\n     -- free : ";
    //     auto it = this->flagSet().find("free")->second.begin();
    //     auto en = this->flagSet().find("free")->second.end();
    //     for ( int cptMark = 0 ; it!=en ; ++it,++cptMark )
    //     {
    //         if ( cptMark > 0 ) *_ostr << " , ";
    //         *_ostr << *it;
    //     }
    // }


#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
    if ( M_alemeshTypeName == "harmonic" && M_isInitHarmonicExtension )
        *_ostr << M_harmonicextensionFactory->getInfo()->str();
#endif
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
    if ( M_alemeshTypeName == "winslow" && M_isInitWinslow )
        *_ostr << M_winslowFactory->getInfo()->str();
#endif
    return _ostr;
}

//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
typename ALE<Convex,Order>::mesh_ptrtype
ALE<Convex,Order>::referenceMesh()
{
    return M_reference_mesh;
}

//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
void
ALE<Convex,Order>::init()
{
    if ( M_alemeshTypeName == "harmonic" )
    {
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
        this->createHarmonicExtension();
#else
        CHECK( false ) << " FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION is turned to OFF";
#endif
    }
    else if ( M_alemeshTypeName == "winslow" )
    {
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
        this->createWinslow();
#else
        CHECK( false ) << " FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW is turned to OFF";
#endif
    }
    else
        CHECK( false ) << "invalid alemesh type : " << M_alemeshTypeName;

    M_dofsHighOnBoundary.clear();
    if ( M_fspaceHigh )
    {
        for ( std::string const& bctype : std::vector<std::string>( { "moving","fixed" } ) )
        {
            auto & dofsHighOnBoundary = M_dofsHighOnBoundary[bctype];
            for ( auto const& faceWrap : markedfaces(M_fspaceHigh->mesh(),this->markers(bctype) ) )
            {
                auto const& face = unwrap_ref( faceWrap );
                auto facedof = M_fspaceHigh->dof()->faceLocalDof( face.id() );
                for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
                    dofsHighOnBoundary.insert( it->index() );
            }
        }
    }
}

//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
void
ALE<Convex,Order>::createALE( std::optional<range_elements_type> const& rangeElt )
{
    if ( rangeElt )
        M_fspaceLow = space_low_type::New( _mesh=M_reference_mesh,_range=*rangeElt );
    else
        M_fspaceLow = space_low_type::New( _mesh=M_reference_mesh );
    M_aleLow.reset( new element_low_type( M_fspaceLow, "low_order_ALE_map" ) );
    M_displacementLow.reset( new element_low_type( M_fspaceLow, "low_order_displacement" ) );
    M_identityLow.reset( new element_low_type( M_fspaceLow, "low_order_identity_map" ) );

    //this->createALEHO( mpl::bool_< ( Order > Order_low ) >() );
    if constexpr (  Order > Order_low )
    {
        M_bHigh = backend_type::build( soption( _name="backend" ), prefixvm(this->prefix(),"ho"), this->worldCommPtr() );
        //M_fspaceHigh = detailALE::buildSpaceHigh<space_low_type,space_high_type>(M_fspaceLow,mpl::bool_<isEqualOrderAndOrderLow>() );
        if constexpr ( std::is_same_v<space_low_type,space_high_type> )
            M_fspaceHigh = M_fspaceLow;
        else
            M_fspaceHigh = space_high_type::New( _mesh=M_fspaceLow->mesh(),_range=M_fspaceLow->template meshSupport<0>() );

#if ALE_WITH_BOUNDARYELEMENT
        M_fspaceHighLocal = space_high_type::New( createSubmesh( M_reference_mesh, boundaryelements(M_reference_mesh) ) );
#endif
        M_aleHigh.reset( new element_high_type( M_fspaceHigh, "high_order_ALE map" ) );
        M_displacementHigh.reset( new element_high_type( M_fspaceHigh, "high_order_displacement" ) );
        M_identityHigh.reset( new element_high_type( M_fspaceHigh, "high_order_identity map" ) );
#if ALE_WITH_BOUNDARYELEMENT
        M_harmonicHigh =  M_bHigh->newMatrix( M_fspaceHighLocal, M_fspaceHighLocal );
        M_rhsHigh = M_bHigh->newVector( M_fspaceHighLocal );
#else
        M_harmonicHigh = M_bHigh->newMatrix( M_fspaceHigh, M_fspaceHigh );
        M_rhsHigh = M_bHigh->newVector( M_fspaceHigh );
#endif

    }
}

//-------------------------------------------------------------------------------------------//

#if 0
template < class Convex, int Order >
void
ALE<Convex,Order>::restart( mesh_ptrtype mesh )
{
    /* Update reference mesh */
    M_reference_mesh = mesh;

    this->createALE();

    this->preCompute();
}
#endif

//-------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//

#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
template < class Convex, int Order >
void
ALE<Convex,Order>::createHarmonicExtension()
{
    this->log( "ALE", "createHarmonicExtension", "start" );

    M_harmonicextensionFactory.reset( new harmonicextension_type( M_fspaceLow,
                                                                  Feel::backend(_rebuild=true,_name=this->prefix() ),
                                                                  prefixvm(this->prefix(),"harmonic"),
                                                                  this->repository() ) );
    M_harmonicextensionFactory->setMarkersInBoundaryCondition(this->bcToMarkers());
    M_harmonicextensionFactory->init();
    M_isInitHarmonicExtension = true;

    this->log( "ALE", "createHarmonicExtension", "finish" );
}
#endif

//-------------------------------------------------------------------------------------------//

#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
template < class Convex, int Order >
void
ALE<Convex,Order>::createWinslow()
{
    this->log( "ALE", "createWinslow", "start" );

    M_winslowFactory.reset(new winslow_type( M_fspaceLow,//M_reference_mesh,
                                             prefixvm(this->prefix(),"winslow") ) );
    M_winslowFactory->setMarkersInBoundaryCondition(this->bcToMarkers());
    M_winslowFactory->init();
    M_isInitWinslow=true;

    this->log( "ALE", "createWinslow", "finish" );
}
#endif


//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
void
ALE<Convex,Order>::generateMap( ale_map_element_type const & dispOnBoundary,
                                ale_map_element_type const & oldDisp )
{
    this->log( "ALE", "generateMap", "start" );

    if ( M_alemeshTypeName == "harmonic" )
        generateLowOrderMap_HARMONIC( dispOnBoundary );
    else if ( M_alemeshTypeName == "winslow" )
        generateLowOrderMap_WINSLOW( dispOnBoundary, oldDisp );
    else
        CHECK( false ) << "wrong arg alemesh.type";

    interpolateLow2High( mpl::bool_< ( Order > Order_low ) >() );

    if (M_doHoCorrection)
    {
        // Run if Order_low is different of Order
        updateBoundaryElements( dispOnBoundary, mpl::bool_< ( Order > Order_low ) >() );
    }

    this->log( "ALE", "generateMap", "finish" );
}

//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
void
ALE<Convex,Order>::generateLowOrderMap_HARMONIC( ale_map_element_type const & dispOnBoundary )
{
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
    M_harmonicextensionFactory->setMarkersInBoundaryCondition(this->bcToMarkers());
    M_harmonicextensionFactory->generateALEMap(dispOnBoundary);

    *M_displacementLow = *M_harmonicextensionFactory->displacement();

    // update ALE map from displacement
    *M_aleLow = *M_identityLow;
    *M_aleLow += *M_displacementLow;
#endif
}

//-------------------------------------------------------------------------------------------//

/**
 * Creates the low order ALE map, given the boundary's displacement
 */
template < class Convex, int Order >
void
ALE<Convex,Order>::generateLowOrderMap_WINSLOW( ale_map_element_type const & dispOnBoundary,
                                                ale_map_element_type const & oldDisp )
{
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
    M_winslowFactory->setMarkersInBoundaryCondition(this->bcToMarkers());
    M_winslowFactory->generateALEMap(dispOnBoundary,oldDisp);

    // interpolate disp
    M_displacementLow->on( _range=elements(M_displacementLow->mesh()),
                           _expr=idv(M_winslowFactory->displacement()) );
    if ( M_displacementLow->functionSpace()->dof()->buildDofTableMPIExtended() )
        sync( *M_displacementLow, "=" );

    *M_aleLow = *M_identityLow;
    *M_aleLow += *M_displacementLow;
#endif
}

//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
void
ALE<Convex,Order>::preCompute()
{
    using namespace Feel::vf;

    M_aleLow->zero();
    M_displacementLow->zero();
    M_identityLow->on( _range=elements( /*M_reference_mesh*/support(M_fspaceLow) ), _expr=P() );

    this->preComputeHO( mpl::bool_< ( Order > Order_low ) >() );
}
//-------------------------------------------------------------------------------------------//
template < class Convex, int Order >
void
ALE<Convex,Order>::preComputeHO( mpl::false_ )
{}
template < class Convex, int Order >
void
ALE<Convex,Order>::preComputeHO( mpl::true_ )
{
    M_aleHigh->zero();
    M_displacementHigh->zero();
    *M_identityHigh = vf::project( _space=M_fspaceHigh, _range=elements( /*M_reference_mesh*/support(M_fspaceHigh) ), _expr=P() );

    using namespace Feel::vf;
#if ALE_WITH_BOUNDARYELEMENT
    auto z = fspaceHighLocal->element( "z" );
    auto w = fspaceHighLocal->element( "w" );
    //M_timer.restart();
    form2( _test=M_fspaceHighLocal, _trial=M_fspaceHighLocal, _matrix=M_harmonicHigh ) =
        integrate( elements(M_fspaceHighLocal->mesh()), trace( trans(gradt(z))*grad(w) ) );
    //form2( fspaceHighLocal, fspaceHighLocal, harmonicHigh ) +=
    //    integrate( boundaryfaces(fspaceHighLocal->mesh()), - trans((gradt(z)*N()))*id(w) );

#else
    auto z = M_fspaceHigh->element( "z" );
    auto w = M_fspaceHigh->element( "w" );
    //M_timer.restart();
    form2( _test=M_fspaceHigh, _trial=M_fspaceHigh, _matrix=M_harmonicHigh ) =
        integrate( elements(M_fspaceHigh->mesh()), trace( trans(gradt(z))*grad(w) ) );
    //form2( fspaceHigh, fspaceHigh, harmonicHigh ) +=
    //    integrate( boundaryfaces(fspaceHigh->mesh()), - trans((gradt(z)*N()))*id(w) );
#endif
    //LOG(INFO) << "[ALE] Time to generate high order harmonic extension matrix: " << M_timer.elapsed() << "\n";
    M_harmonicHigh->close();
}
//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
void
ALE<Convex,Order>::interpolateLow2High( mpl::bool_<true> )
{
    interpolate( M_fspaceHigh, *M_displacementLow, *M_displacementHigh );
    interpolate( M_fspaceHigh, *M_aleLow, *M_aleHigh );
}
template < class Convex, int Order >
void
ALE<Convex,Order>::interpolateLow2High( mpl::bool_<false> )
{}

//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
void
ALE<Convex,Order>::updateBoundaryElements( ale_map_element_type const & dispOnBoundary, mpl::bool_<true> )
{

    updateBoundaryElements(dispOnBoundary,mpl::bool_<true>(), mpl::int_<0>() );
}
//-------------------------------------------------------------------------------------------//
template < class Convex, int Order >
void
ALE<Convex,Order>::updateBoundaryElements( ale_map_element_type const & dispOnBoundary, mpl::bool_<true>, mpl::int_<0>  )
{

    static const bool hasDofInElt = mpl::if_< mpl::equal_to<mpl::int_<Dim>,mpl::int_<3> >,
                                              typename mpl::if_< mpl::greater< mpl::int_<space_high_type::fe_type::nDofPerVolume>,mpl::int_<0> >,
                                                                 mpl::bool_<true>,
                                                                 mpl::bool_<false> >::type,
                                              typename mpl::if_< mpl::greater< mpl::int_<space_high_type::fe_type::nDofPerFace>,mpl::int_<0> >,
                                                                 mpl::bool_<true>,
                                                                 mpl::bool_<false> >::type >::type::value;

    updateBoundaryElements(dispOnBoundary,mpl::bool_<true>(), mpl::int_<0>(), mpl::bool_<hasDofInElt>() );

}
//-------------------------------------------------------------------------------------------//
template < class Convex, int Order >
void
ALE<Convex,Order>::updateBoundaryElements( ale_map_element_type const & dispOnBoundary, mpl::bool_<true>, mpl::int_<0>, mpl::bool_<false> )
{
    using namespace Feel::vf;
    boost::mpi::timer timerHighOrderBoundary;

    //auto correction = idv(dispOnBoundary) - idv(M_displacementLow);
    auto correction = idv(dispOnBoundary) - idv(M_displacementHigh);
    auto zero = 0*one();

    auto dispHighAux = M_fspaceHigh->element();
    dispHighAux.on(_range=markedfaces(M_fspaceHigh->mesh(), this->markers("moving")),
                   _expr=correction);
    // also at fixed boundary (if boundary is ho, else is zero)
    dispHighAux.on(_range=markedfaces(M_fspaceHigh->mesh(), this->markers("fixed")),
                   _expr=correction);
    // synch values
    auto itFindDofsMoving = M_dofsHighOnBoundary.find( "moving" );
    CHECK( itFindDofsMoving != M_dofsHighOnBoundary.end() ) << "dofsHighOnBoundary moving not initialized";
    auto itFindDofsFixed = M_dofsHighOnBoundary.find( "fixed" );
    CHECK( itFindDofsFixed != M_dofsHighOnBoundary.end() ) << "dofsHighOnBoundary fixed not initialized";
    std::set<size_type> dofsHighOnBoundaryMovingAndFixed;
    dofsHighOnBoundaryMovingAndFixed.insert( itFindDofsMoving->second.begin(),itFindDofsMoving->second.end() );
    dofsHighOnBoundaryMovingAndFixed.insert( itFindDofsFixed->second.begin(),itFindDofsFixed->second.end() );
    sync( dispHighAux, "=", dofsHighOnBoundaryMovingAndFixed );
    *M_displacementHigh += dispHighAux;

    *M_aleHigh = *M_identityHigh;
    *M_aleHigh += *M_displacementHigh;

    double timerHighOrderBoundaryElapsed = timerHighOrderBoundary.elapsed();
    if (this->verboseSolverTimer()) Feel::FeelModels::Log(this->prefix()+".ALE","updateBoundaryElements<true,0,false>",
                                                   (boost::format("finish in %1%")%timerHighOrderBoundaryElapsed).str(),
                                                   this->worldComm(),this->verboseSolverTimerAllProc());
}

//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
void
ALE<Convex,Order>::updateBoundaryElements( ale_map_element_type const & dispOnBoundary, mpl::bool_<true>, mpl::int_<0>, mpl::bool_<true> )
{
    boost::mpi::timer M_timer;
#if !ALE_WITH_BOUNDARYELEMENT
    using namespace Feel::vf;

    boost::mpi::timer timerHighOrderBoundary;

    //auto correction = val(idv(dispOnBoundary) + idv(M_identityHigh) - idv(M_aleLow));
    auto correction = idv(dispOnBoundary) - idv(M_displacementLow);
    auto zero = 0*one();

    M_rhsHigh->zero();
    M_rhsHigh->close();

    auto v = M_fspaceHigh->element( "v" );
    // Impose homogeneous Dirichlet boundary condition on the points of the mesh
    form2( _test=M_fspaceHigh, _trial=M_fspaceHigh, _matrix=M_harmonicHigh ) +=
        on( _range=internalfaces(M_reference_mesh),
            _element=v,
            _rhs=M_rhsHigh,
            _expr=zero );
#if 0
    std::list<std::string> listMarkerNoDispImposed;
    for ( uint16_type i=0; i < this->flagSet("fixed").size(); ++i )
        listMarkerNoDispImposed.push_back( this->flagSet("fixed",i) );
    for ( uint16_type i=0; i < this->flagSet("free").size(); ++i )
        listMarkerNoDispImposed.push_back( this->flagSet("free",i) );
    // TODO integrator on with list of markedfaces
#endif

    /* Impose boundary conditions for the moving boundary = position - aleLow */
    if ( this->bcToMarkers().find("moving") != this->bcToMarkers().end() )
    {
        auto dmlose = Feel::FeelModels::detail::distributeMarkerListOnSubEntity(M_reference_mesh, this->markers("moving") );
        auto const& faceMarkers = std::get<0>( dmlose );
        if ( !faceMarkers.empty() )
        {
            form2( _test=M_fspaceHigh, _trial=M_fspaceHigh, _matrix=M_harmonicHigh ) +=
                on( _range=markedfaces( M_reference_mesh, faceMarkers ),
                    _element=v,
                    _rhs=M_rhsHigh,
                    _expr=correction );
        }
    }
    for ( std::string bctype : { "fixed","free" } )
    {
        if ( this->bcToMarkers().find(bctype) != this->bcToMarkers().end() )
        {
            auto dmlose = Feel::FeelModels::detail::distributeMarkerListOnSubEntity(M_reference_mesh, this->markers(bctype) );
            auto const& faceMarkers = std::get<0>( dmlose );
            if ( !faceMarkers.empty() )
            {
                form2( _test=M_fspaceHigh, _trial=M_fspaceHigh, _matrix=M_harmonicHigh ) +=
                    on( _range=markedfaces( M_reference_mesh, faceMarkers ),
                        _element=v,
                        _rhs=M_rhsHigh,
                        _expr=idv(dispOnBoundary)/*zero*/ );
            }
        }
    }

    //std::cout << "  -- ale : solving for high order correction" << std::endl;
    M_timer.restart();

    M_bHigh->solve(_matrix=M_harmonicHigh, _solution=v, _rhs=M_rhsHigh);

    LOG(INFO) << "[ALE] Time to solve high order harmonic extension operator: " << M_timer.elapsed() << "\n";

    //v.updateGlobalValues();

    // be careful accumulate=true works only because v is 0 on all dof in the intersection of element including on moving boundary
#if 0 //!!!!!!
    *M_displacementHigh = vf::project( _space=M_fspaceHigh, _range=boundaryelements(M_fspaceHigh->mesh()), _expr=idv(v), _accumulate=true);
#else

    //auto dispHighAux = vf::project( _space=M_fspaceHigh, _range=boundaryelements(M_fspaceHigh->mesh()), _expr=idv(v));

    auto dispHighAux = M_fspaceHigh->element();
    dispHighAux = vf::project( _space=M_fspaceHigh, _range=boundaryelements(M_fspaceHigh->mesh()), _expr=idv(v));

    *M_displacementHigh += dispHighAux;
    //displacementHigh += v;
#endif

    *M_aleHigh = *M_identityHigh;
    *M_aleHigh += *M_displacementHigh;

    if ( this->verbose() )
    {
        double normDisp = M_displacementHigh->l2Norm();
        double normDispImposed = dispOnBoundary.l2Norm();
        double normDispAux = dispHighAux.l2Norm();
        double normSol = v.l2Norm();
        Feel::FeelModels::Log(this->prefix()+".ALE","updateBoundaryElements<true>",
                       (boost::format("normDisp %1% normDispImposed %2% normDispAux %3% normSol%4%") %normDisp %normDispImposed %normDispAux %normSol ).str(),
                       this->worldComm(),this->verboseAllProc());
    }

#if 0
    std::vector<geo_element_type> it_elt( 1 );
    /* For in all the elements in contact with the moving boundary */
    for( auto iter = elementsOnMovingBoundary.begin(); iter != elementsOnMovingBoundary.end(); ++iter )
    {
        it_elt.clear();
        std::cout << "Treat element with id " << *iter << "\n";

        mesh_ptrtype oneElementMesh ( new mesh_type );
        it_elt.push_back( reference_mesh->element(*iter) );

        std::cout << "Create submesh from this element...\n";
        reference_mesh->createSubmesh( *oneElementMesh, it_elt.begin(), it_elt.end() );
        std::cout << "Create submesh from this element: done!\n";

        std::cout << "One Element mesh has " << oneElementMesh->numElements() << " elements\n";
        std::cout << "with id : " << oneElementMesh->element(0).id() << "\n";

        std::cout << "Face markers: \n";
        std::cout << "Face 0: " << oneElementMesh->element(0).edge(0).marker() << "\n";
        std::cout << "Face 1: " << oneElementMesh->element(0).edge(1).marker() << "\n";
        std::cout << "Face 2: " << oneElementMesh->element(0).edge(2).marker() << "\n";
    }
#endif

    double timerHighOrderBoundaryElapsed = timerHighOrderBoundary.elapsed();
    std::ostringstream ostr;ostr<< timerHighOrderBoundaryElapsed <<"s";
    if (this->verboseSolverTimer()) Feel::FeelModels::Log(this->prefix()+".ALE","updateBoundaryElements<true,0,true>", "finish in "+ostr.str(),
                                                   this->worldComm(),this->verboseSolverTimerAllProc());

#endif

}


template < class Convex, int Order >
void
ALE<Convex,Order>::updateBoundaryElements( ale_map_element_type const & dispOnBoundary, mpl::bool_<true>, mpl::int_<1> )
{
#if ALE_WITH_BOUNDARYELEMENT
    using namespace Feel::vf;

    boost::mpi::timer timerHighOrderBoundary;

    //auto correction = val(idv(dispOnBoundary) + idv(M_identityHigh) - idv(M_aleLow));
    auto correction = idv(dispOnBoundary) - idv(M_displacementLow);
    auto zero = 0*one();

    //auto rhs = b->newVector( fspaceHighLocal );
    M_rhsHigh->zero();
    M_rhsHigh->close();

    std::cout << "  -- ale : building high order correction\n";
    auto v = fspaceHighLocal->element( "v" );
    std::cout << "     + setting 0 displacement on all internal faces" << std::endl;
    /* Impose homogeneous Dirichlet boundary condition on the points of the mesh */
    form2( _test=M_fspaceHighLocal, _trial=M_fspaceHighLocal, _matrix=M_harmonicHigh ) +=
        on( _range=internalfaces(M_fspaceHighLocal->mesh()), _element=v, _rhs=M_rhsHigh, _expr=zero );

    /* Impose boundary conditions for the moving boundary = position - aleLow */
    std::cout << "     + projecting correction" << std::endl;
    auto corrLocal = vf::project( _space=M_fspaceHighLocal, _range=elements(M_fspaceHighLocal->mesh()), _expr=correction );
    std::cout << "     + setting correction displacement on all boundary faces" << std::endl;


#if 1
    form2( _test=M_fspaceHighLocal, _trial=M_fspaceHighLocal, _matrix=M_harmonicHigh ) +=
        on( _range=boundaryfaces( M_fspaceHighLocal->mesh()), _element=v, _rhs=M_rhsHigh, _expr=print(idv(corrLocal),"corrlocal=") );
#else

    std::cout << "\n ----------1--------\n";

    auto hola = boost::make_tuple( mpl::size_t<MESH_FACES>(),
                                   M_fspaceHighLocal->mesh()->beginFace(),
                                   M_fspaceHighLocal->mesh()->endFace() );


    form2( _test=M_fspaceHighLocal, _trial=M_fspaceHighLocal, _matrix=M_harmonicHigh ) +=
        on( _range=boundaryfaces( M_fspaceHighLocal->mesh()), _element=v, _rhs=M_rhsHigh, _expr=zero );

    std::cout << "\n ----------2--------\n";

    for ( uint16_type i=0; i < this->flagSet("moving").size(); ++i )
    {
        form2( _test=M_fspaceHighLocal, _trial=M_fspaceHighLocal, _matrix=M_harmonicHigh ) +=
            on( _range=markedfaces( M_fspaceHighLocal->mesh(), this->flagSet("moving",i) ), _element=v,
                _rhs=M_rhsHigh, _expr=correction );
        std::cout << "\n ----------3--------\n";
    }
#endif

    //std::cout << "  -- ale : solving for high order correction" << std::endl;
    M_timer.restart();

    M_bHigh->solve(_matrix=M_harmonicHigh, _solution=v, _rhs=M_rhsHigh);

    LOG(INFO) << "[ALE] Time to solve high order harmonic extension operator: " << M_timer.elapsed() << "\n";

    //v.updateGlobalValues();

    // be careful accumulate=true works only because v is 0 on all dof in the intersection of element including on moving boundary

    //auto dispHighAux = vf::project( _space=M_fspaceHigh, _range=boundaryelements(M_fspaceHigh->mesh()), _expr=idv(v));

    auto dispHighAux = M_fspaceHigh->element();
    dispHighAux = vf::project( _space=M_fspaceHigh, _range=boundaryelements(M_fspaceHigh->mesh()), _expr=idv(v));

#endif

}


template < class Convex, int Order >
void
ALE<Convex,Order>::updateBoundaryElements( ale_map_element_type const & dispOnBoundary,//std::vector<elem_type> const& polyDisplacementSet,
                                           mpl::bool_<false> )
{
    //std::cout << "\nWARNING!!!!!!updateBoundaryElements\n";
}

template < class Convex, int Order >
void
ALE<Convex,Order>::initMetricMeshAdaptation()
{
    this->M_metricMeshAdaptation = std::make_shared<typename super_type::metricmeshadaptation_type>( M_fspaceLow, "");
    this->M_metricMeshAdaptation->init();
}

template < class Convex, int Order >
void
ALE<Convex,Order>::updateMetricMeshAdaptation( Expr<GinacExVF<2>> const& e )
{
    if ( !this->M_metricMeshAdaptation )
        this->initMetricMeshAdaptation();
    this->M_metricMeshAdaptation->update( e );
    this->updateMetricMeshAdaptationForUse();
}

template < class Convex, int Order >
void
ALE<Convex,Order>::updateMetricMeshAdaptation( typename super_type::metricmeshadaptation_type::element_scalar_type const& u )
{
    if ( !this->M_metricMeshAdaptation )
        this->initMetricMeshAdaptation();
    this->M_metricMeshAdaptation->update( u );
    this->updateMetricMeshAdaptationForUse();
}

template < class Convex, int Order >
void
ALE<Convex,Order>::updateMetricMeshAdaptationForUse()
{
    if ( M_alemeshTypeName == "winslow" )
    {
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
        M_winslowFactory->setMetricMeshAdaptation( this->M_metricMeshAdaptation->scalarMetric() );
#else
        CHECK( false ) << "WINSLOW is disable";
#endif
    }
    else
        CHECK( false ) << "mesh adaptation only supported by Winslow model";
}

} // namespace ALE_IMPL
} // namespace FeelModels
} // namespace Feel

