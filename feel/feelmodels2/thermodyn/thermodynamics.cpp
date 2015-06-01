/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelmodels2/thermodyn/thermodynamics.hpp>

#include <feel/feelvf/vf.hpp>
/*#include <feel/feelvf/form.hpp>
#include <feel/feelvf/on.hpp>
#include <feel/feelvf/operators.hpp>
 #include <feel/feelvf/operations.hpp>*/

namespace Feel {

namespace FeelModels {

THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
THERMODYNAMICS_CLASS_TEMPLATE_TYPE::ThermoDynamics( bool __isStationary,
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

    THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
    void
    THERMODYNAMICS_CLASS_TEMPLATE_TYPE::loadConfigBCFile()
    {
        this->clearMarkerDirichletBC();
        this->clearMarkerNeumannBC();

        M_bcDirichlet = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "Dirichlet" );
        for( auto const& d : M_bcDirichlet )
            this->addMarkerDirichletBC("elimination", marker(d) );
        M_bcNeumann = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "Neumann" );
        for( auto const& d : M_bcNeumann )
            this->addMarkerNeumannBC(super_type::NeumannBCShape::SCALAR,marker(d));

        M_volumicForcesProperties = this->modelProperties().boundaryConditions().getScalarFields( "temperature", "VolumicForces" );
    }

    THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
    void
    THERMODYNAMICS_CLASS_TEMPLATE_TYPE::loadConfigMeshFile(std::string const& geofilename)
    {
        CHECK( false ) << "not allow";
    }


    THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
    void
    THERMODYNAMICS_CLASS_TEMPLATE_TYPE::init(bool buildMethodNum)
    {
        super_type::init( buildMethodNum, this->shared_from_this() );
    }

    THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
    void
    THERMODYNAMICS_CLASS_TEMPLATE_TYPE::solve()
    {
        M_bcDirichlet.setParameterValues( this->modelProperties().parameters().toParameterValues() );
        M_bcNeumann.setParameterValues( this->modelProperties().parameters().toParameterValues() );
        M_volumicForcesProperties.setParameterValues( this->modelProperties().parameters().toParameterValues() );
        super_type::solve();
    }

    THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
    void
    THERMODYNAMICS_CLASS_TEMPLATE_TYPE::updateBCStrongDirichletLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
    {
        if ( M_bcDirichlet.empty() ) return;

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

    THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
    void
    THERMODYNAMICS_CLASS_TEMPLATE_TYPE::updateSourceTermLinearPDE( vector_ptrtype& F, bool buildCstPart ) const
    {
        if ( this->M_overwritemethod_updateSourceTermLinearPDE != NULL )
        {
            this->M_overwritemethod_updateSourceTermLinearPDE(F,buildCstPart);
            return;
        }

        if ( M_volumicForcesProperties.empty() ) return;

        if ( !buildCstPart )
        {
            auto myLinearForm = form1( _test=this->spaceTemperature(), _vector=F,
                                       _rowstart=this->rowStartInVector() );
            auto const& v = *this->fieldTemperature();

            for( auto const& d : M_volumicForcesProperties )
            {
                if ( marker(d).empty() )
                    myLinearForm +=
                        integrate( _range=elements(this->mesh()),
                                   _expr= expression(d)*id(v),
                                   _geomap=this->geomap() );
                else
                    myLinearForm +=
                        integrate( _range=markedelements(this->mesh(),marker(d)),
                                   _expr= expression(d)*id(v),
                                   _geomap=this->geomap() );
            }
        }
    }

    THERMODYNAMICS_CLASS_TEMPLATE_DECLARATIONS
    void
    THERMODYNAMICS_CLASS_TEMPLATE_TYPE::updateWeakBCLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F,bool buildCstPart) const
    {
        if ( M_bcNeumann.empty() ) return;

        if ( !buildCstPart )
        {
            auto mesh = this->mesh();
            auto Xh = this->spaceTemperature();
            auto const& v = *this->fieldTemperature();
            auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                                      _pattern=size_type(Pattern::COUPLED),
                                                      _rowstart=this->rowStartInMatrix(),
                                                      _colstart=this->colStartInMatrix() );
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
