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


namespace Feel
{
namespace FeelModels
{
namespace ALE_IMPL
{

namespace detailALE
{
template < typename SpaceLowType,typename SpaceHighType >
boost::shared_ptr<SpaceHighType>
buildSpaceHigh(boost::shared_ptr<SpaceLowType> spaceLow, bool moveGhostEltFromExtendedStencil, mpl::bool_<false> /**/ )
{
    return SpaceHighType::New( _mesh=spaceLow->mesh(),_worldscomm=spaceLow->worldsComm(),
                               _extended_doftable=std::vector<bool>(1,moveGhostEltFromExtendedStencil) );
}

template < typename SpaceLowType,typename SpaceHighType >
boost::shared_ptr<SpaceHighType>
buildSpaceHigh(boost::shared_ptr<SpaceLowType> spaceLow, bool moveGhostEltFromExtendedStencil, mpl::bool_<true> /**/ )
{
    return spaceLow;
}
}



template < class Convex, int Order >
ALE<Convex,Order>::ALE( mesh_ptrtype mesh, std::string prefix, WorldComm const& worldcomm, bool moveGhostEltFromExtendedStencil,
                        ModelBaseRepository const& modelRep )
    :
    super_type( mesh,prefix,worldcomm,moveGhostEltFromExtendedStencil,modelRep ),
    M_verboseSolverTimer(boption(_prefix=this->prefix(),_name="verbose_solvertimer")),
    M_verboseSolverTimerAllProc(boption(_prefix=this->prefix(),_name="verbose_solvertimer_allproc")),
    M_reference_mesh( mesh ),
    M_alemeshTypeName( soption( _name="type",_prefix=this->prefix() ) ),
    M_doHoCorrection( boption(_prefix=this->prefix(),_name="apply-ho-correction") ),
    M_isInitHarmonicExtension( false ),
    M_isInitWinslow( false ),
    M_moveGhostEltFromExtendedStencil( moveGhostEltFromExtendedStencil )
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix(),"constructor", "start",
                                        this->worldComm(),this->verboseAllProc());
    this->createALE();

    this->preCompute();
#if 0
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
#endif
    if (this->verbose()) Feel::FeelModels::Log(this->prefix(),"constructor", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//-------------------------------------------------------------------------------------------//
/**
 * copy constructor
 */
template < class Convex, int Order >
ALE<Convex,Order>::ALE( ALE const& tc )
    :
    super_type( tc ),
    M_verboseSolverTimer( tc.M_verboseSolverTimer ),
    M_verboseSolverTimerAllProc( tc.M_verboseSolverTimerAllProc ),
    M_reference_mesh( tc.M_reference_mesh ),
    M_fspaceLow( tc.M_fspaceLow ),
    M_fspaceHigh( tc.M_fspaceHigh ),
    M_fspaceHighLocal( tc.M_fspaceHighLocal ),
    M_aleLow( tc.M_aleLow ),
    M_displacementLow( tc.M_displacementLow ),
    M_identityLow( tc.M_identityLow ),
    M_aleHigh( tc.M_aleHigh ),
    M_displacementHigh( tc.M_displacementHigh ),
    M_identityHigh( tc.M_identityHigh ),
    M_bHigh( tc.M_bHigh),
    M_harmonicHigh( tc.M_harmonicHigh ),
    M_rhsHigh( tc.M_rhsHigh ),
    M_alemeshTypeName( tc.M_alemeshTypeName ),
    M_doHoCorrection( tc.M_doHoCorrection ),
#if defined( FSI_ENABLE_HARMONICEXTENSION )
    M_harmonicextensionFactory( tc.M_harmonicextensionFactory ),
#endif
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
    M_winslowFactory( tc.M_winslowFactory ),
#endif
    M_isInitHarmonicExtension( tc.M_isInitHarmonicExtension ),
    M_isInitWinslow( tc.M_isInitWinslow )
{
}

//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
ALE<Convex,Order>::~ALE() {}

//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
boost::shared_ptr<std::ostringstream>
ALE<Convex,Order>::getInfo() const
{
    boost::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );
    *_ostr << "\n   Physical Markers";
    if ( this->flagSet().find("fixed") != this->flagSet().end() )
    {
        *_ostr << "\n     -- fixed : ";
        auto it = this->flagSet().find("fixed")->second.begin();
        auto en = this->flagSet().find("fixed")->second.end();
        for ( int cptMark = 0 ; it!=en ; ++it,++cptMark )
        {
            if ( cptMark > 0 ) *_ostr << " , ";
            *_ostr << *it;
        }
    }
    if ( this->flagSet().find("moving") != this->flagSet().end() )
    {
        *_ostr << "\n     -- moving : ";
        auto it = this->flagSet().find("moving")->second.begin();
        auto en = this->flagSet().find("moving")->second.end();
        for ( int cptMark = 0 ; it!=en ; ++it,++cptMark )
        {
            if ( cptMark > 0 ) *_ostr << " , ";
            *_ostr << *it;
        }
    }
    if ( this->flagSet().find("free") != this->flagSet().end() )
    {
        *_ostr << "\n     -- free : ";
        auto it = this->flagSet().find("free")->second.begin();
        auto en = this->flagSet().find("free")->second.end();
        for ( int cptMark = 0 ; it!=en ; ++it,++cptMark )
        {
            if ( cptMark > 0 ) *_ostr << " , ";
            *_ostr << *it;
        }
    }


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
            for ( auto const& faceWrap : markedfaces(M_fspaceHigh->mesh(),this->flagSet(bctype) ) )
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
ALE<Convex,Order>::createALE()
{
    M_fspaceLow = space_low_type::New( _mesh=M_reference_mesh,_worldscomm=std::vector<WorldComm>(1,this->worldComm()),
                                       _extended_doftable=std::vector<bool>(1,M_moveGhostEltFromExtendedStencil) );
    M_aleLow.reset( new element_low_type( M_fspaceLow, "low_order_ALE_map" ) );
    M_displacementLow.reset( new element_low_type( M_fspaceLow, "low_order_displacement" ) );
    M_identityLow.reset( new element_low_type( M_fspaceLow, "low_order_identity_map" ) );

    //M_fspaceHigh = detailALE::buildSpaceHigh<space_low_type,space_high_type>(M_fspaceLow, M_moveGhostEltFromExtendedStencil,
    //                                                                         mpl::bool_<isEqualOrderAndOrderLow>() );
    this->createALEHO( mpl::bool_< ( Order > Order_low ) >() );
}
template < class Convex, int Order >
void
ALE<Convex,Order>::createALEHO( mpl::true_ )
{
    M_bHigh = backend_type::build( soption( _name="backend" ), prefixvm(this->prefix(),"ho"), this->worldComm() );
    M_fspaceHigh = detailALE::buildSpaceHigh<space_low_type,space_high_type>(M_fspaceLow, M_moveGhostEltFromExtendedStencil,
                                                                             mpl::bool_<isEqualOrderAndOrderLow>() );
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
template < class Convex, int Order >
void
ALE<Convex,Order>::createALEHO( mpl::false_ )
{

}

//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
void
ALE<Convex,Order>::restart( mesh_ptrtype mesh )
{
    /* Update reference mesh */
    M_reference_mesh = mesh;

    this->createALE();

    this->preCompute();
}

//-------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------//

#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
template < class Convex, int Order >
void
ALE<Convex,Order>::createHarmonicExtension()
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix(),"createHarmonicExtension", "start",
                                        this->worldComm(),this->verboseAllProc());
    /*bool useGhostEltFromExtendedStencil = M_fspaceLow->dof()->buildDofTableMPIExtended() && M_reference_mesh->worldComm().localSize()>1;
     M_harmonicextensionFactory.reset( new harmonicextension_type( M_reference_mesh, M_bLow, prefixvm(M_prefix,"alemesh.harmonic"),
     M_worldComm,useGhostEltFromExtendedStencil) );*/
    M_harmonicextensionFactory.reset( new harmonicextension_type( M_fspaceLow,
                                                                  Feel::backend(_rebuild=true,_name=this->prefix() ),
                                                                  prefixvm(this->prefix(),"harmonic"),
                                                                  this->repository() ) );
    M_harmonicextensionFactory->setflagSet(this->flagSet());
    M_harmonicextensionFactory->init();
    M_isInitHarmonicExtension = true;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix(),"createHarmonicExtension", "finish",
                                        this->worldComm(),this->verboseAllProc());
}
#endif

//-------------------------------------------------------------------------------------------//

#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_WINSLOW )
template < class Convex, int Order >
void
ALE<Convex,Order>::createWinslow()
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix(),"createWinslow", "start",
                                        this->worldComm(),this->verboseAllProc());

    bool useGhostEltFromExtendedStencil = M_fspaceLow->dof()->buildDofTableMPIExtended() && M_reference_mesh->worldComm().localSize()>1;
    M_winslowFactory.reset(new winslow_type( M_reference_mesh,
                                             prefixvm(this->prefix(),"winslow"),
                                             this->worldComm(),useGhostEltFromExtendedStencil) );
    //M_winslowFactory.reset( new winslow_type( M_fspaceLow,M_bLow,prefixvm(M_prefix,"alemesh.winslow")/*M_prefix*/) );
    M_winslowFactory->setflagSet(this->flagSet());
    M_winslowFactory->init();
    M_isInitWinslow=true;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix(),"createWinslow", "finish",
                                        this->worldComm(),this->verboseAllProc());
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
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MethodsNum","generateMap", "start",
                                        this->worldComm(),this->verboseAllProc());

    if ( M_alemeshTypeName == "harmonic" )
        generateLowOrderMap_HARMONIC( dispOnBoundary );
    else if ( M_alemeshTypeName == "winslow" )
        generateLowOrderMap_WINSLOW( dispOnBoundary, oldDisp );
    else FEELPP_ASSERT (false).error( "wrong arg alemesh.type" );

    interpolateLow2High( mpl::bool_< ( Order > Order_low ) >() );

    if (M_doHoCorrection)
    {
        // Run if Order_low is different of Order
        updateBoundaryElements( dispOnBoundary, mpl::bool_< ( Order > Order_low ) >() );
    }

#if 0
    CHECK( M_displacementLow->functionSpace()->nLocalDofWithGhost() == dispOnBoundary.functionSpace()->nLocalDofWithGhost() ) << "BLEM ALE 1\n";
    CHECK( M_displacementLow->functionSpace()->nLocalDofWithoutGhost() == dispOnBoundary.functionSpace()->nLocalDofWithoutGhost() ) << "BLEM ALE 1\n";
    for ( uint16_type i=0; i < this->flagSet("moving").size(); ++i )
    {
        for ( auto const& faces : markedfaces(dispOnBoundary.mesh(), this->flagSet("moving",i) ) )
        {
            //for ( uint16_type l =0; l < ale_map_functionspace_type::dof_type::fe_type::nLocalDof; ++l )
            for ( uint16_type l =0; l < M_displacementLow->functionSpace()->dof()->nLocalDofOnFace(true); ++l )
            {
                int ncdof  = ale_map_functionspace_type::dof_type::is_product?ale_map_functionspace_type::dof_type::nComponents:1;
                for ( uint16_type c1 = 0; c1 < ncdof; ++c1 )
                {
                    //const size_type thedof =  boost::get<0>(dispOnBoundary.functionSpace()->dof()->localToGlobal(faces,l,c1) );
                    const size_type thedof =  boost::get<0>(dispOnBoundary.functionSpace()->dof()->faceLocalToGlobal(faces.id(),l,c1) );

                    CHECK( thedof < M_displacementLow->functionSpace()->nLocalDofWithGhost() ) << "invalid dof " << thedof
                                                                                               << "faces.id() " <<faces.id()
                                                                                               <<" l "<<l<<" c1 " << c1 << "\n";
                    CHECK( std::abs( (*M_displacementLow)(thedof) - dispOnBoundary(thedof) ) < 1e-12 ) << "error ALE at dof "<< thedof << " :  "
                                                                                                       << (*M_displacementLow)(thedof) << " vs " << dispOnBoundary(thedof) << "\n";
                }
            }
        }
    }
#endif


    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".MethodsNum","generateMap", "finish",
                                        this->worldComm(),this->verboseAllProc());
}

//-------------------------------------------------------------------------------------------//

template < class Convex, int Order >
void
ALE<Convex,Order>::generateLowOrderMap_HARMONIC( ale_map_element_type const & dispOnBoundary )
{
#if defined( FEELPP_TOOLBOXES_HAS_MESHALE_HARMONICEXTENSION )
    using namespace Feel::vf;

    M_harmonicextensionFactory->setflagSet(this->flagSet());
    M_harmonicextensionFactory->generateALEMap(dispOnBoundary);

    // interpolate disp
#if 1
    //interpolate( M_displacementLow->functionSpace(), M_harmonicextensionFactory->displacement(), M_displacementLow, INTERPOLATE_SAME_MESH );
    *M_displacementLow = *M_harmonicextensionFactory->displacement();
#else
    bool useGhostEltFromExtendedStencil = M_displacementLow->functionSpace()->dof()->buildDofTableMPIExtended() &&
        M_harmonicextensionFactory->functionSpace()->dof()->buildDofTableMPIExtended() &&
        M_displacementLow->mesh()->worldComm().localSize()>1;
    *M_displacementLow = vf::project(_space=M_displacementLow->functionSpace(),
                                     _range=elements(M_displacementLow->mesh(),useGhostEltFromExtendedStencil),
                                     _expr=idv(M_harmonicextensionFactory->displacement()) );
#endif
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
    M_winslowFactory->setflagSet(this->flagSet());
    M_winslowFactory->generateALEMap(dispOnBoundary,oldDisp);

    // interpolate disp
    //interpolate( M_fspaceHigh, M_displacementLow, M_displacementHigh, INTERPOLATE_SAME_MESH );
    bool useGhostEltFromExtendedStencil = M_displacementLow->functionSpace()->dof()->buildDofTableMPIExtended() &&
        M_winslowFactory->functionSpace()->dof()->buildDofTableMPIExtended() &&
        M_displacementLow->mesh()->worldComm().localSize()>1;
    EntityProcessType entityProcess = (useGhostEltFromExtendedStencil)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    *M_displacementLow = vf::project(_space=M_displacementLow->functionSpace(),
                                     _range=elements(M_displacementLow->mesh(),entityProcess),
                                     _expr=idv(M_winslowFactory->displacement()) );
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
    //bool useGhostEltFromExtendedStencil = M_fspaceLow->dof()->buildDofTableMPIExtended() && M_reference_mesh->worldComm().localSize()>1;
    EntityProcessType entityProcess = (M_moveGhostEltFromExtendedStencil)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    *M_identityLow = vf::project( _space=M_fspaceLow, _range=elements( M_reference_mesh,entityProcess ), _expr=P() );

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
    EntityProcessType entityProcess = (M_moveGhostEltFromExtendedStencil)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    *M_identityHigh = vf::project( _space=M_fspaceHigh, _range=elements( M_reference_mesh,entityProcess ), _expr=P() );

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
    interpolate( M_fspaceHigh, *M_displacementLow, *M_displacementHigh, INTERPOLATE_SAME_MESH );
    interpolate( M_fspaceHigh, *M_aleLow, *M_aleHigh, INTERPOLATE_SAME_MESH );
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
    dispHighAux.on(_range=markedfaces(M_fspaceHigh->mesh(), this->flagSet("moving")),
                   _expr=correction);
    // also at fixed boundary (if boundary is ho, else is zero)
    dispHighAux.on(_range=markedfaces(M_fspaceHigh->mesh(), this->flagSet("fixed")),
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
    for ( uint16_type i=0; i < this->flagSet("moving").size(); ++i )
    {
        form2( _test=M_fspaceHigh, _trial=M_fspaceHigh, _matrix=M_harmonicHigh ) +=
            on( _range=markedfaces( M_reference_mesh, this->flagSet("moving",i) ),
                _element=v,
                _rhs=M_rhsHigh,
                _expr=correction );
    }

    for ( uint16_type i=0; i < this->flagSet("fixed").size(); ++i )
    {
        form2( _test=M_fspaceHigh, _trial=M_fspaceHigh, _matrix=M_harmonicHigh ) +=
            on( _range=markedfaces( M_reference_mesh, this->flagSet("fixed",i) ),
                _element=v,
                _rhs=M_rhsHigh,
                _expr=idv(dispOnBoundary)/*zero*/ );
    }
    for ( uint16_type i=0; i < this->flagSet("free").size(); ++i )
    {
        form2( _test=M_fspaceHigh, _trial=M_fspaceHigh, _matrix=M_harmonicHigh ) +=
            on( _range=markedfaces( M_reference_mesh, this->flagSet("free",i) ),
                _element=v,
                _rhs=M_rhsHigh,
                _expr=idv(dispOnBoundary)/*zero*/ );
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


} // namespace ALE_IMPL
} // namespace FeelModels
} // namespace Feel

