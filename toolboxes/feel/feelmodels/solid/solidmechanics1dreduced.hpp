/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef FEELPP_TOOLBOXES_SOLIDMECHANICS_1DREDUCED_HPP
#define FEELPP_TOOLBOXES_SOLIDMECHANICS_1DREDUCED_HPP 1

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelts/newmark.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>
#include <feel/feelmodels/solid/solidmechanics1dreducedboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisDisplacementType>
class SolidMechanics1dReduced : public ModelNumerical,
                                public ModelPhysics<ConvexType::nRealDim>
{
    using super_physics_type = ModelPhysics<ConvexType::nRealDim>;
public :
    using super_type = ModelNumerical;
    using size_type = typename super_type::size_type;
    typedef SolidMechanics1dReduced<ConvexType,BasisDisplacementType> self_type;
    // TODO static assert only 1d
    typedef std::shared_ptr<self_type> self_ptrtype;
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static inline const uint16_type nDim = convex_type::nDim;
    static inline const uint16_type nOrderGeo = convex_type::nOrder;
    static inline const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    // materials properties
    typedef MaterialsProperties<nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;

    // basis
    using basis_displacement_type = BasisDisplacementType;
    static inline const uint16_type nOrderDisplacement = basis_displacement_type::nOrder;
    // function space displacement
    typedef FunctionSpace<mesh_type, bases<basis_displacement_type> > space_displacement_type;
    typedef std::shared_ptr<space_displacement_type> space_displacement_ptrtype;
    typedef typename space_displacement_type::element_type element_displacement_type;
    typedef std::shared_ptr<element_displacement_type> element_displacement_ptrtype;

    typedef typename space_displacement_type::component_functionspace_type space_displacement_component_type;
    typedef typename space_displacement_type::component_functionspace_ptrtype space_displacement_component_ptrtype;
    typedef typename space_displacement_component_type::element_type element_displacement_component_type;
    typedef std::shared_ptr<element_displacement_component_type> element_displacement_component_ptrtype;
    //___________________________________________________________________________________//
    // time step newmark
    typedef Newmark<space_displacement_component_type> newmark_type;
    typedef std::shared_ptr<newmark_type> newmark_ptrtype;
    //___________________________________________________________________________________//
    // exporter
    typedef Exporter<mesh_type,nOrderGeo> exporter_type;
    typedef std::shared_ptr<exporter_type> exporter_ptrtype;


#if 0
    typedef Range<mesh_type,MESH_ELEMENTS> range_elt_type;
    typedef OperatorInterpolation<space_stress_scal_type, space_1dreduced_type ,range_elt_type> op_interpolation2dTo1d_normalstress_type;
    typedef std::shared_ptr<op_interpolation2dTo1d_normalstress_type> op_interpolation2dTo1d_normalstress_ptrtype;
#endif

    struct FieldTag
    {
        static auto displacement( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
    };

    SolidMechanics1dReduced( std::string const& prefix,
                             std::string const& keyword,// = "cfpde",
                             worldcomm_ptr_t const& worldComm,// = Environment::worldCommPtr(),
                             std::string const& subPrefix,//  = "",
                             ModelBaseRepository const& modelRep = ModelBaseRepository() );

    std::shared_ptr<self_type> shared_from_this() { return std::dynamic_pointer_cast<self_type>( super_type::shared_from_this() ); }

    void init();

    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

    // physical parameters
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { CHECK( !this->isUpdatedForUse() ) << "setMaterialsProperties can be called only before called isUpdatedForUse";  M_materialsProperties = mp; }

    // mesh
    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }

    BlocksBaseGraphCSR buildBlockMatrixGraph() const override;

    space_displacement_component_ptrtype const& functionSpace1dReduced() const { return M_spaceDisp; }

    element_displacement_component_type & fieldDisplacementScal1dReduced() { return *M_fieldDisp; }
    element_displacement_component_type const & fieldDisplacementScal1dReduced() const { return *M_fieldDisp; }
    element_displacement_component_ptrtype fieldDisplacementScal1dReducedPtr() const { return M_fieldDisp; }

    element_displacement_type & fieldDisplacementVect1dReduced() { return *M_fieldDisp_vect; }
    element_displacement_type const & fieldDisplacementVect1dReduced() const { return *M_fieldDisp_vect; }
    element_displacement_ptrtype const & fieldDisplacementVect1dReducedPtr() const { return M_fieldDisp_vect; }

    element_displacement_component_type & fieldVelocityScal1dReduced() { return *M_fieldVelocity; }
    element_displacement_component_type const& fieldVelocityScal1dReduced() const { return *M_fieldVelocity; }
    element_displacement_type & fieldVelocityVect1dReduced() { return *M_fieldVelocity_vect; }
    element_displacement_type const& fieldVelocityVect1dReduced() const { return *M_fieldVelocity_vect; }
    element_displacement_ptrtype const& fieldVelocityVect1dReducedPtr() const { return M_fieldVelocity_vect; }

    newmark_ptrtype & timeStepNewmark() { return M_timeStepNewmark; }
    newmark_ptrtype const& timeStepNewmark() const { return M_timeStepNewmark; }

    //backend_ptrtype const& backend1dReduced() const { return M_backend_1dReduced; }
    //model_algebraic_factory_ptrtype algebraicFactory1dReduced() const { return M_algebraicFactory_1dReduced; }
    //BlocksBaseGraphCSR buildBlockMatrixGraph1dReduced() const;

    double thickness1dReduced() const { return M_thickness_1dReduced; }
    double radius1dReduced() const { return M_radius_1dReduced; }

    void startTimeStep();
    void updateTimeStep();
    std::shared_ptr<TSBase> timeStepBase() const { return M_timeStepNewmark; }

    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& paramValues );


    void predictorDispl();
    void updateVelocity();

    void updateInterfaceDispFrom1dDisp();
    void updateInterfaceVelocityFrom1dVelocity();
    element_displacement_ptrtype extendVelocity1dReducedVectorial( element_displacement_component_type const& vel1d ) const;


    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    auto modelFields( std::string const& prefix = "" ) const
        {
            return this->modelFields( this->fieldDisplacementScal1dReducedPtr(), prefix );
        }

    template <typename DisplacementFieldType>
    auto modelFields( DisplacementFieldType const& field_s, std::string const& prefix = "" ) const
        {
            auto mfield_disp = modelField<FieldCtx::ID,DisplacementFieldType>( FieldTag::displacement(this) );
            mfield_disp.add( FieldTag::displacement(this), prefix,"displacement", field_s, "s", this->keyword() );
            return Feel::FeelModels::modelFields( mfield_disp );
        }


    template <typename SymbExprType>
    auto exprPostProcessExports( SymbExprType const& se, std::string const& prefix = "" ) const
        {
            //return hana::concat( this->materialsProperties()->exprPostProcessExports( this->mesh(), this->physicsAvailable(),se ),
            //                     this->exprPostProcessExportsToolbox( se, prefix ) );
            return this->materialsProperties()->exprPostProcessExports( this->mesh(), this->physicsAvailable(),se );
        }

    template <typename ModelFieldsType,typename SymbolsExpr,typename ExportsExprType>
    void exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr );

    template <typename SymbolsExpr>
    void exportResults( double time, SymbolsExpr const& se )
        {
            this->exportResults( time, this->modelFields(), se, this->exprPostProcessExports( se ) );
        }


    // assembly methods for linear system
    //void updateLinearPDE( DataUpdateLinear & data ) const override;
    template <typename ModelContextType>
    void updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const;
    //void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;
    template <typename ModelContextType>
    void updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mctx ) const;


private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initMaterialProperties();
    void initFunctionSpaces();
    void initBoundaryConditions();
    void initTimeStep();
    void initInitialConditions();
    void initPostProcess() override;

private :

    // physical parameter
    materialsproperties_ptrtype M_materialsProperties;

    // mesh
    Range<mesh_type,MESH_ELEMENTS> M_rangeMeshElements;

    // function space
    space_displacement_component_ptrtype M_spaceDisp;
    //element disp,vel,acc
    element_displacement_component_ptrtype M_fieldDisp;
    element_displacement_component_ptrtype M_fieldVelocity;
    element_displacement_component_ptrtype M_fieldAcceleration;
    // vectorial 1d_reduced space
    space_displacement_ptrtype M_spaceDispVect;
    element_displacement_ptrtype M_fieldDisp_vect;
    element_displacement_ptrtype M_fieldVelocity_vect;
    // time discretisation
    std::string M_timeStepping;
    newmark_ptrtype M_timeStepNewmark;

    // boundary conditions
    using boundary_conditions_type = SolidMechanics1dReducedBoundaryConditions<nRealDim>;
    std::shared_ptr<boundary_conditions_type> M_boundaryConditions;
#if 0
    map_scalar_field<2> M_bcDirichlet;
    map_scalar_field<2> M_volumicForcesProperties;
#endif

    // exporter
    exporter_ptrtype M_exporter;

    // axi-sym properties
    double M_thickness_1dReduced, M_radius_1dReduced;

};


template< typename ConvexType, typename BasisDisplacementType>
template <typename ModelFieldsType, typename SymbolsExpr, typename ExportsExprType>
void
SolidMechanics1dReduced<ConvexType,BasisDisplacementType>::exportResults( double time, ModelFieldsType const& mfields, SymbolsExpr const& symbolsExpr, ExportsExprType const& exportsExpr )
{
    this->log("SolidMechanics1dReduced","exportResults", "start");
    this->timerTool("PostProcessing").start();

    this->executePostProcessExports( M_exporter, time, mfields, symbolsExpr, exportsExpr );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("SolidMechanics1dReduced","exportResults", "finish");
}



template< typename ConvexType, typename BasisDisplacementType>
template <typename ModelContextType>
void
SolidMechanics1dReduced<ConvexType,BasisDisplacementType>::updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    this->log( "SolidMechanics1dReduced","updateLinearPDE","start");

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _buildCstPart = data.buildCstPart();
    bool BuildNonCstPart = !_buildCstPart;
    bool BuildCstPart = _buildCstPart;

    //auto mesh = M_mesh_1dReduced;
    auto Xh1 = M_spaceDisp;
    //auto u = Xh1->element(), v = Xh1->element();
    auto const& u = *M_fieldDisp;
    auto const& v = u;

    auto const& se = mctx.symbolsExpr();

    double epp = this->thickness1dReduced();//0.1;
    double k=2.5;
    double G=1e5;
    double gammav=0.01;
    double R0 = this->radius1dReduced();//0.5;
    //---------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------//

    auto rowStartInMatrix = this->rowStartInMatrix();
    auto colStartInMatrix = this->colStartInMatrix();
    auto rowStartInVector = this->rowStartInVector();
    auto bilinearForm1dreduced = form2( _test=Xh1,_trial=Xh1,_matrix=A,
                                        _rowstart=rowStartInMatrix,
                                        _colstart=colStartInMatrix );
    auto linearFormDisplacement = form1( _test=Xh1, _vector=F,
                                      _rowstart=rowStartInVector );

    //---------------------------------------------------------------------------------------//
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nRealDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = this->materialsProperties()->materialProperties( matName );
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );

            auto const& densityProp = this->materialsProperties()->density( matName );
            auto densityExpr = expr( densityProp.expr(), se );

            auto E = expr( matProperties.property( "Young-modulus" ).exprScalar(), se );
            auto mu = expr( matProperties.property( "Poisson-ratio" ).exprScalar(), se );

            // acceleration term
            if (!this->isStationary())
            {
                if (BuildCstPart)
                {
                    bilinearForm1dreduced +=
                        integrate( _range=range,
                                   _expr=this->timeStepNewmark()->polySecondDerivCoefficient()*densityExpr*epp*idt(u)*id(v) );
                }

                if (BuildNonCstPart)
                {
                    linearFormDisplacement +=
                        integrate( _range=range,
                                   _expr=densityExpr*epp*idv(this->timeStepNewmark()->polyDeriv() )*id(v) );
                }

                // viscoealstic term
                if (BuildCstPart)
                {
                    bilinearForm1dreduced +=
                        integrate( _range=range,
                                   _expr=this->timeStepNewmark()->polyFirstDerivCoefficient()*gammav*dxt(u)*dx(v) );
                }
                if (BuildNonCstPart)
                {
                    linearFormDisplacement +=
                        integrate( _range=range,
                                   _expr=gammav*dxv(this->timeStepNewmark()->polyFirstDeriv())*dx(v) );
                }
            }
            //---------------------------------------------------------------------------------------//
            // reaction term
            if (BuildCstPart)
            {
                bilinearForm1dreduced +=
                    integrate( _range=range,
                               _expr=(E*epp/((1-mu*mu)*R0*R0))*idt(u)*id(v) );
            }
            //---------------------------------------------------------------------------------------//
            // diffusion term
            if (BuildCstPart)
            {
                bilinearForm1dreduced +=
                    integrate( _range=range,
                               _expr=k*G*epp*dxt(u)*dx(v) );
            }

            //---------------------------------------------------------------------------------------//
            // source term
            if (BuildNonCstPart)
            {
                // linearFormDisplacement +=
                //     integrate( _range=M_rangeMeshElements1dReduced,
                //                _expr=idv(*M_stress_1dReduced)*id(v) );
#if 0 // VINCENT
                for( auto const& d : M_volumicForcesProperties )
                {
                    //std::cout << "apply
                    auto rangeBodyForceUsed = markers(d).empty()? /*M_rangeMeshElements*/elements(this->mesh()) : markedelements(this->mesh(),markers(d));
                    linearFormDisplacement +=
                        integrate( _range=rangeBodyForceUsed,
                                   _expr= expression(d,se)*id(v),
                                   _geomap=this->geomap() );
                }
#endif
            }
        }
    }
    //---------------------------------------------------------------------------------------//

    this->log( "SolidMechanics1dReduced","updateLinearPDE","finish");

}

template< typename ConvexType, typename BasisDisplacementType>
template <typename ModelContextType>
void
SolidMechanics1dReduced<ConvexType,BasisDisplacementType>::updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    if ( !M_boundaryConditions->hasTypeDofElimination() )
        return;

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto const& se = mctx.symbolsExpr();

    auto Xh = M_spaceDisp;
    auto const& u = *M_fieldDisp;
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );

    M_boundaryConditions->applyDofEliminationLinear( bilinearForm, F, this->mesh(), u, se );
}


} // namespace FeelModels
} // namespace Feel

#endif
