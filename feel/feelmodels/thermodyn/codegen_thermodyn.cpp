/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include "codegen_thermodyn.hpp"

#include <feel/feelfilters/geotool.hpp>

#undef THERMODYNAMICS_MESH
#undef THERMODYNAMICS_MESH1
#undef THERMODYNAMICS_MESH2
#undef THERMODYNAMICS_MESH3
#include "thermodyn.mesh"


namespace Feel {

namespace FeelModels {

THERMODYNAMICS_CLASS_NAME::THERMODYNAMICS_CLASS_NAME( std::string __prefix,
                                                      bool __buildMesh,
                                                      WorldComm const& __worldComm,
                                                      std::string __subPrefix,
                                                      std::string __appliShortRepository )
    :
    super_type( __prefix,__worldComm,__buildMesh,__subPrefix,__appliShortRepository)
{
    this->log("ThermoDynamics","constructor", "start");

    this->setFilenameSaveInfo( prefixvm(this->prefix(),"ThermoDynamics.info") );
    //-----------------------------------------------------------------------------//
    // load info from .bc file
    this->loadConfigBCFile();
    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    // build mesh, space, exporter,...
    if (__buildMesh) this->build();
    //-----------------------------------------------------------------------------//
    this->log("ThermoDynamics","constructor", "finish");
}

void
THERMODYNAMICS_CLASS_NAME::loadConfigBCFile()
{
    auto const bcDef = THERMODYNAMICS_BC(this/*->shared_from_this()*/);
    this->clearMarkerDirichletBC();
    this->clearMarkerNeumannBC();
    ForEachBC( bcDef, cl::dirichlet,
               this->addMarkerDirichletBC("elimination",PhysicalName) );
    ForEachBC( bcDef, cl::neumann_scal,
               this->addMarkerNeumannBC(NeumannBCShape::SCALAR,PhysicalName) );
}

void
THERMODYNAMICS_CLASS_NAME::loadConfigMeshFile(std::string const& geofilename)
{
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","loadConfigMeshFile", "start",
                                               this->worldComm(),this->verboseAllProc());

#if !defined(THERMODYNAMICS_MESH) && !defined(THERMODYNAMICS_MESH1) && !defined(THERMODYNAMICS_MESH2) && !defined(THERMODYNAMICS_MESH3)
    CHECK( false ) << " THERMODYNAMICS_MESH is not define : probably .mesh is wrong\n";
#endif

    switch ( this->geotoolMeshIndex() )
    {
    default :
    case 0 :
    {
#if defined(THERMODYNAMICS_MESH)
        THERMODYNAMICS_MESH(this->meshSize(), geofilename );
        M_mesh=mesh;
#endif
    }
    break;
    case 1 :
    {
#if defined(THERMODYNAMICS_MESH1)
        THERMODYNAMICS_MESH1(this->meshSize(), geofilename );
        M_mesh=mesh;
#endif
    }
    break;
    case 2 :
    {
#if defined(THERMODYNAMICS_MESH2)
        THERMODYNAMICS_MESH2(this->meshSize(), geofilename );
        M_mesh=mesh;
#endif
    }
    break;
    {
#if defined(THERMODYNAMICS_MESH3)
        THERMODYNAMICS_MESH3(this->meshSize(), geofilename );
        M_mesh=mesh;
#endif
    }
    break;

    } // switch

    this->log("ThermoDynamics","loadConfigMeshFile", "finish");
}


void
THERMODYNAMICS_CLASS_NAME::init(bool buildMethodNum)
{
    super_type::init( buildMethodNum, this->shared_from_this() );
}


void
THERMODYNAMICS_CLASS_NAME::updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    this->log("ThermoDynamics","updateBCStrongDirichletLinearPDE","start");

    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto const& u = *this->fieldTemperature();

    auto rowStartInMatrix = this->rowStartInMatrix();
    auto colStartInMatrix = this->colStartInMatrix();
    auto rowStartInVector = this->rowStartInVector();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );

    auto const bcDef = THERMODYNAMICS_BC(this->shared_from_this());

    ForEachBC( bcDef, cl::dirichlet,
               bilinearForm_PatternCoupled +=
               /**/ on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",PhysicalName ) /*PhysicalName*/),
                        _element=u,
                        _rhs=F,
                        _expr=Expression ) );

    this->log("ThermoDynamics","updateBCStrongDirichletLinearPDE","finish");
}

void
THERMODYNAMICS_CLASS_NAME::updateSourceTermLinearPDE( vector_ptrtype& F, bool buildCstPart ) const
{
    if ( M_overwritemethod_updateSourceTermLinearPDE != NULL )
    {
        M_overwritemethod_updateSourceTermLinearPDE(F,buildCstPart);
        return;
    }

#if defined(THERMODYNAMICS_VOLUME_FORCE)
    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto rowStartInVector = this->rowStartInVector();
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=rowStartInVector );
    auto const& v = *this->fieldTemperature();
    auto f = THERMODYNAMICS_VOLUME_FORCE(this->shared_from_this()) ;

    if ( !buildCstPart )
    {
        myLinearForm +=
            integrate( _range=elements(mesh),
                       _expr= f*id(v),
#if defined(THERMODYNAMICS_VOLUME_FORCE_QUADORDER)
                       _quad=_Q<THERMODYNAMICS_VOLUME_FORCE_QUADORDER>(),
                       _quad1=_Q<THERMODYNAMICS_VOLUME_FORCE_QUADORDER>(),
#endif
                       _geomap=this->geomap() );
    }
#endif
}

void
THERMODYNAMICS_CLASS_NAME::updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const
{
    auto mesh = this->mesh();
    auto Xh = this->spaceTemperature();
    auto rowStartInMatrix = this->rowStartInMatrix();
    auto colStartInMatrix = this->colStartInMatrix();
    auto rowStartInVector = this->rowStartInVector();
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );
    auto const& v = *this->fieldTemperature();
    //boundaries conditions
    auto const bcDef = THERMODYNAMICS_BC(this->shared_from_this());

    if ( !buildCstPart )
    {
        ForEachBC( bcDef, cl::neumann_scal,
                   bilinearForm_PatternCoupled +=
                   /**/ integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,PhysicalName) ),
                                   _expr= Expression*id(v),
                                   _geomap=this->geomap() ) );
    }

}


} // end namespace FeelModels
} // end namespace Feel
