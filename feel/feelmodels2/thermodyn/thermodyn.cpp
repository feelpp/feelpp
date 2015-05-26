/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

//#include "thermodyn.hpp"
#include <feel/feelmodels2/thermodyn/thermodyn.hpp>

#include <feel/feelvf/vf.hpp>
/*#include <feel/feelvf/form.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/operators.hpp>
 #include <feel/feelvf/operations.hpp>*/

namespace Feel {

namespace FeelModels {

    template< typename ConvexType, int OrderTemp>
    ThermoDynamics<ConvexType,OrderTemp>::ThermoDynamics( bool __isStationary,
                                                          std::string __prefix,
                                                          WorldComm const& __worldComm,
                                                          bool __buildMesh,
                                                          std::string __subPrefix,
                                                          std::string __appliShortRepository )
        :
        super_type( __isStationary,__prefix,__worldComm,__buildMesh,__subPrefix,__appliShortRepository)
    {
        if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","constructor", "start",
                                            this->worldComm(),this->verboseAllProc());

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
        if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","constructor", "finish",
                                            this->worldComm(),this->verboseAllProc());
    }

    template< typename ConvexType, int OrderTemp>
    void
    ThermoDynamics<ConvexType,OrderTemp>::loadConfigBCFile()
    {
        this->clearMarkerDirichletBC();
        this->clearMarkerNeumannBC();

        M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( this->prefix()/*"thermo"*/, "Dirichlet" );
        for( auto const& d : M_bcDirichlet )
            this->addMarkerDirichletBC("elimination", marker(d) );
        M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( this->prefix()/*"thermo"*/, "Neumann" );
        for( auto const& d : M_bcNeumann )
            this->addMarkerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d));
    }

    template< typename ConvexType, int OrderTemp>
    void
    ThermoDynamics<ConvexType,OrderTemp>::loadConfigMeshFile(std::string const& geofilename)
    {
        CHECK( false ) << "not allow";
    }


    template< typename ConvexType, int OrderTemp>
    void
    ThermoDynamics<ConvexType,OrderTemp>::init(bool buildMethodNum)
    {
        super_type::init( buildMethodNum, this->shared_from_this() );
    }


    template< typename ConvexType, int OrderTemp>
    void
    ThermoDynamics<ConvexType,OrderTemp>::updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
    {
        if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","updateBCStrongDirichletLinearPDE","start",
                                            this->worldComm(),this->verboseAllProc());

        auto mesh = this->mesh();
        auto Xh = this->spaceTemperature();
        auto const& u = *this->fieldTemperature();
        auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                                  _pattern=size_type(Pattern::COUPLED),
                                                  _rowstart=this->rowStartInMatrix(),
                                                  _colstart=this->colStartInMatrix() );

        for( auto const& d : M_bcDirichlet )
        {
            bilinearForm_PatternCoupled +=
                on( _range=markedfaces(mesh, this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                    _element=u,_rhs=F,_expr=expression(d) );
        }

        if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".ThermoDynamics","updateBCStrongDirichletLinearPDE","finish",
                                                   this->worldComm(),this->verboseAllProc());
    }

    template< typename ConvexType, int OrderTemp>
    void
    ThermoDynamics<ConvexType,OrderTemp>::updateSourceTermLinearPDE( vector_ptrtype& F, bool buildCstPart ) const
    {
#if 0 // TODO
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
#endif
    }

    template< typename ConvexType, int OrderTemp>
    void
    ThermoDynamics<ConvexType,OrderTemp>::updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const
    {
        auto mesh = this->mesh();
        auto Xh = this->spaceTemperature();
        auto const& v = *this->fieldTemperature();
        auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                                  _pattern=size_type(Pattern::COUPLED),
                                                  _rowstart=this->rowStartInMatrix(),
                                                  _colstart=this->colStartInMatrix() );

        if ( !buildCstPart )
        {
            for( auto const& d : M_bcNeumann )
            {
                bilinearForm_PatternCoupled +=
                    integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d)) ),
                               _expr= expression(d)*id(v),
                               _geomap=this->geomap() );
            }
        }
    }


} // end namespace FeelModels
} // end namespace Feel
