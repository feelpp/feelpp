/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2011-07-17

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
 \file fluidmechanics.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-17
 */

#ifndef FEELPP_TOOLBOXES_FLUIDMECHANICS_HPP
#define FEELPP_TOOLBOXES_FLUIDMECHANICS_HPP 1


#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/projectors.hpp>

#include <feel/feelmodels/modelcore/traits.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>

#include <feel/feelmodels/fluid/fluidmechanicsmaterialproperties.hpp>
#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>

#if defined( FEELPP_MODELS_HAS_MESHALE )
#include <feel/feelmodels/modelmesh/meshale.hpp>
#endif

#include <feel/feelmodels/modelcore/stabilizationglsparameterbase.hpp>
#include <feel/feelmodels/modelcore/rangedistributionbymaterialname.hpp>


namespace Feel
{
namespace FeelModels
{

/** 
 * Fluid Mechanics Toolbox
 * \ingroup Toolboxes
 */
template< typename ConvexType, typename BasisVelocityType,
          typename BasisPressureType = Lagrange< (BasisVelocityType::nOrder>1)? (BasisVelocityType::nOrder-1):BasisVelocityType::nOrder, Scalar,Continuous,PointSetFekete>,
          typename BasisDVType=Lagrange<0, Scalar,Discontinuous/*,PointSetFekete*/> >
class FluidMechanics : public ModelNumerical,
                       public ModelPhysics<ConvexType::nDim>,
                       public std::enable_shared_from_this< FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType,BasisDVType> >,
                       public MarkerManagementDirichletBC,
                       public MarkerManagementNeumannBC,
                       public MarkerManagementALEMeshBC,
                       public MarkerManagementSlipBC,
                       public MarkerManagementPressureBC
{
public:
    using super_type = ModelNumerical;
    using super2_type = ModelPhysics<ConvexType::nDim>;
    using size_type = typename super_type::size_type;
    typedef FluidMechanics< ConvexType,BasisVelocityType,BasisPressureType,BasisDVType > self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // trace mesh
    typedef typename mesh_type::trace_mesh_type trace_mesh_type;
    typedef typename mesh_type::trace_mesh_ptrtype trace_mesh_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // basis fluid
    static const bool useMixedBasis = true;
    static const uint16_type nOrderVelocity = BasisVelocityType::nOrder;
    static const uint16_type nOrderPressure = BasisPressureType::nOrder;
    typedef BasisVelocityType basis_fluid_u_type;
    typedef BasisPressureType basis_fluid_p_type;
    typedef Lagrange<0, Scalar,Continuous> basis_l_type;
    //___________________________________________________________________________________//
    // mixed basis
    //typedef bases<basis_fluid_u_type,basis_fluid_p_type> basis_fluid_type;
    //___________________________________________________________________________________//
    // function space velocity
    typedef FunctionSpace<mesh_type, bases<basis_fluid_u_type> > space_velocity_type;
    typedef std::shared_ptr<space_velocity_type> space_velocity_ptrtype;
    typedef typename space_velocity_type::element_type element_velocity_type;
    typedef std::shared_ptr<element_velocity_type> element_velocity_ptrtype;
    typedef typename space_velocity_type::element_external_storage_type element_velocity_external_storage_type;
    // function space component of velocity
    typedef typename space_velocity_type::component_functionspace_type component_space_velocity_type;
    typedef std::shared_ptr<component_space_velocity_type> component_space_velocity_ptrtype;
    typedef typename component_space_velocity_type::element_type component_element_velocity_type;
    typedef std::shared_ptr<component_element_velocity_type> component_element_velocity_ptrtype;
    // function space pressure
    typedef FunctionSpace<mesh_type, bases<basis_fluid_p_type> > space_pressure_type;
    typedef std::shared_ptr<space_pressure_type> space_pressure_ptrtype;
    typedef typename space_pressure_type::element_type element_pressure_type;
    typedef std::shared_ptr<element_pressure_type> element_pressure_ptrtype;
    typedef typename space_pressure_type::element_external_storage_type element_pressure_external_storage_type;
    // function space for lagrange multiplier which impose the mean pressure
    typedef FunctionSpace<mesh_type, bases<basis_l_type> > space_meanpressurelm_type;
    typedef std::shared_ptr<space_meanpressurelm_type> space_meanpressurelm_ptrtype;
    // function space velocity on trace
    typedef FunctionSpace<trace_mesh_type, bases<basis_fluid_u_type> > space_trace_velocity_type;
    typedef std::shared_ptr<space_trace_velocity_type> space_trace_velocity_ptrtype;
    typedef typename space_trace_velocity_type::element_type element_trace_velocity_type;
    typedef std::shared_ptr<element_trace_velocity_type> element_trace_velocity_ptrtype;
    // function space component of velocity on trace
    typedef typename space_trace_velocity_type::component_functionspace_type space_trace_velocity_component_type;
    typedef std::shared_ptr<space_trace_velocity_component_type> space_trace_velocity_component_ptrtype;
    typedef typename space_trace_velocity_component_type::element_type element_trace_velocity_component_type;
    typedef std::shared_ptr<element_trace_velocity_component_type> element_trace_velocity_component_ptrtype;
    // function space P0 continuous vectorial on trace
    typedef FunctionSpace<trace_mesh_type, bases<Lagrange<0, Vectorial,Continuous>>> space_trace_p0c_vectorial_type;
    typedef std::shared_ptr<space_trace_p0c_vectorial_type> space_trace_p0c_vectorial_ptrtype;
    typedef typename space_trace_p0c_vectorial_type::element_type element_trace_p0c_vectorial_type;
    typedef std::shared_ptr<element_trace_p0c_vectorial_type> element_trace_p0c_vectorial_ptrtype;
    // function space P0 continuous scalar on trace
    typedef FunctionSpace<trace_mesh_type, bases<Lagrange<0, Scalar,Continuous>>> space_trace_p0c_scalar_type;
    typedef std::shared_ptr<space_trace_p0c_scalar_type> space_trace_p0c_scalar_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // ALE
#if defined( FEELPP_MODELS_HAS_MESHALE )
    typedef MeshALE<convex_type> mesh_ale_type;
    typedef std::shared_ptr<mesh_ale_type> mesh_ale_ptrtype;
    // ref ALE mesh
    typedef typename mesh_ale_type::mesh_ref_type mesh_ref_type;
    typedef typename mesh_ale_type::mesh_ref_ptrtype mesh_ref_ptrtype;
    // mesh disp
    typedef typename mesh_ale_type::ale_map_functionspace_type space_mesh_disp_type;
    typedef typename mesh_ale_type::ale_map_element_type element_mesh_disp_type;
    typedef std::shared_ptr<element_mesh_disp_type> element_mesh_disp_ptrtype;
    // mesh velocity (whole domain)
    typedef typename mesh_ale_type::ale_map_element_type element_meshvelocity_type;
    typedef std::shared_ptr<element_meshvelocity_type> element_meshvelocity_ptrtype;
    // case where structure displacement is scalar!
    typedef typename space_mesh_disp_type::component_functionspace_type space_mesh_disp_scalar_type;
    typedef std::shared_ptr<space_mesh_disp_scalar_type> space_mesh_disp_scalar_ptrtype;
    typedef typename space_mesh_disp_scalar_type::element_type element_mesh_disp_scalar_type;
    typedef std::shared_ptr<element_mesh_disp_scalar_type> element_mesh_disp_scalar_ptrtype;
#endif
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // function space stress
    //typedef bases<Lagrange<nOrderVelocity-1+space_alemapdisc_type::basis_type::nOrder, Vectorial,Discontinuous,PointSetFekete> > basis_stress_type;
    typedef Lagrange<nOrderVelocity-1+mesh_type::nOrder, Vectorial,Discontinuous,PointSetFekete> basis_normalstress_type;
    typedef FunctionSpace<trace_mesh_type, bases<basis_normalstress_type> > space_normalstress_type;
    typedef std::shared_ptr<space_normalstress_type> space_normalstress_ptrtype;
    typedef typename space_normalstress_type::element_type element_normalstress_type;
    typedef std::shared_ptr<element_normalstress_type> element_normalstress_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // function space vorticity
    typedef typename mpl::if_< mpl::equal_to<mpl::int_<nDim>,mpl::int_<2> >,
                               Lagrange<nOrderVelocity-1, Scalar,Discontinuous,PointSetFekete>,
                               Lagrange<nOrderVelocity-1, Vectorial,Discontinuous,PointSetFekete> >::type basis_vorticity_type;
    typedef FunctionSpace<mesh_type, bases<basis_vorticity_type> > space_vorticity_type;
    typedef std::shared_ptr<space_vorticity_type> space_vorticity_ptrtype;
    typedef typename space_vorticity_type::element_type element_vorticity_type;
    typedef std::shared_ptr<element_vorticity_type> element_vorticity_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // functionspace for rho, mu, nu
    typedef BasisDVType basis_densityviscosity_type;
    static const uint16_type nOrderDensityViscosity = BasisDVType::nOrder;
    typedef FunctionSpace<mesh_type, bases<basis_densityviscosity_type> > space_densityviscosity_type;
    // viscosity model desc
    typedef FluidMechanicsMaterialProperties<space_densityviscosity_type> material_properties_type; // TO REMOVE
    typedef std::shared_ptr<material_properties_type> material_properties_ptrtype; // TO REMOVE
    typedef MaterialsProperties<nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;


    typedef bases<Lagrange<nOrderVelocity, Vectorial,Continuous,PointSetFekete> > basis_vectorial_PN_type;
    typedef FunctionSpace<mesh_type, basis_vectorial_PN_type> space_vectorial_PN_type;
    typedef std::shared_ptr<space_vectorial_PN_type> space_vectorial_PN_ptrtype;
    typedef typename space_vectorial_PN_type::element_type element_vectorial_PN_type;
    typedef std::shared_ptr<element_vectorial_PN_type> element_vectorial_PN_ptrtype;
    //___________________________________________________________________________________//
    // stabilization
    typedef StabilizationGLSParameterBase<mesh_type> stab_gls_parameter_type;
    typedef std::shared_ptr<stab_gls_parameter_type> stab_gls_parameter_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // algebraic tools
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;
    typedef typename model_algebraic_factory_type::graph_type graph_type;
    typedef typename model_algebraic_factory_type::graph_ptrtype graph_ptrtype;
    typedef typename model_algebraic_factory_type::indexsplit_type indexsplit_type;
    typedef typename model_algebraic_factory_type::indexsplit_ptrtype indexsplit_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // time
    typedef Bdf<space_velocity_type> bdf_velocity_type;
    typedef std::shared_ptr<bdf_velocity_type> bdf_velocity_ptrtype;
    typedef Bdf<space_pressure_type> savets_pressure_type;
    typedef std::shared_ptr<savets_pressure_type> savets_pressure_ptrtype;
    typedef Bdf<space_trace_p0c_vectorial_type> bdf_trace_p0c_vectorial_type;
    typedef std::shared_ptr<bdf_trace_p0c_vectorial_type> bdf_trace_p0c_vectorial_ptrtype;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    typedef elements_reference_wrapper_t<mesh_type> range_elements_type;
    typedef faces_reference_wrapper_t<mesh_type> range_faces_type;
    //___________________________________________________________________________________//
    // fluid inlet
    typedef typename basis_fluid_u_type::component_basis_type basis_fluidinlet_type;
    typedef FunctionSpace<trace_mesh_type, bases<basis_fluidinlet_type> > space_fluidinlet_type;
    typedef std::shared_ptr<space_fluidinlet_type> space_fluidinlet_ptrtype;
    typedef typename space_fluidinlet_type::element_type element_fluidinlet_type;
    typedef std::shared_ptr<element_fluidinlet_type> element_fluidinlet_ptrtype;
    typedef OperatorInterpolation<space_fluidinlet_type, component_space_velocity_type,
                                  range_faces_type> op_interpolation_fluidinlet_type;
    typedef std::shared_ptr<op_interpolation_fluidinlet_type> op_interpolation_fluidinlet_ptrtype;
    //___________________________________________________________________________________//
    // windkessel model
    typedef bases<Lagrange<0, Scalar,Continuous>,Lagrange<0, Scalar,Continuous> > basis_fluidoutlet_windkessel_type;
    typedef FunctionSpace<trace_mesh_type, basis_fluidoutlet_windkessel_type > space_fluidoutlet_windkessel_type;
    typedef std::shared_ptr<space_fluidoutlet_windkessel_type> space_fluidoutlet_windkessel_ptrtype;
    typedef typename space_fluidoutlet_windkessel_type::element_type element_fluidoutlet_windkessel_type;
    typedef std::shared_ptr<element_fluidoutlet_windkessel_type> element_fluidoutlet_windkessel_ptrtype;
#if defined( FEELPP_MODELS_HAS_MESHALE )
    typedef typename MeshALE<typename trace_mesh_type::shape_type>::ale_map_functionspace_type space_fluidoutlet_windkessel_mesh_disp_type;
    typedef std::shared_ptr<space_fluidoutlet_windkessel_mesh_disp_type> space_fluidoutlet_windkessel_mesh_disp_ptrtype;
    typedef typename space_fluidoutlet_windkessel_mesh_disp_type::element_type element_fluidoutlet_windkessel_mesh_disp_type;
    typedef std::shared_ptr<element_fluidoutlet_windkessel_mesh_disp_type> element_fluidoutlet_windkessel_mesh_disp_ptrtype;
    // typedef boost::tuple<boost::mpl::size_t<MESH_ELEMENTS>,
    //                      typename MeshTraits<trace_mesh_type>::element_const_iterator,
    //                      typename MeshTraits<trace_mesh_type>::element_const_iterator> range_fluidoutlet_windkessel_type;
    typedef OperatorInterpolation<space_mesh_disp_type,
                                  space_fluidoutlet_windkessel_mesh_disp_type/*,
                                                                              range_fluidoutlet_windkessel_type*/> op_interpolation_fluidoutlet_windkessel_meshdisp_type;
    typedef std::shared_ptr<op_interpolation_fluidoutlet_windkessel_meshdisp_type> op_interpolation_fluidoutlet_windkessel_meshdisp_ptrtype;
#endif

    struct FieldTag
    {
        static auto velocity( self_type const* t ) { return ModelFieldTag<self_type,0>( t ); }
        static auto pressure( self_type const* t ) { return ModelFieldTag<self_type,1>( t ); }
        static auto mesh_displacement( self_type const* t ) { return ModelFieldTag<self_type,2>( t ); }
        //static auto body_translational_velocity( BodyBoundaryCondition const* t ) { return BodyBoundaryCondition::FieldTag::translational_velocity( t ); }
        //static auto body_angular_velocity( BodyBoundaryCondition const* t ) { return BodyBoundaryCondition::FieldTag::angular_velocity( t ); }
    };


    //___________________________________________________________________________________//

    class Body //: public ModelPhysics<nDim>,
    //  public std::enable_shared_from_this<Body>
    {
    public :
        using moment_of_inertia_type = typename mpl::if_< mpl::equal_to<mpl::int_<nDim>,mpl::int_<3> >,
                                                       eigen_matrix_type<nDim, nDim>,
                                                       eigen_matrix_type<1, 1> >::type;

        Body() = default;
            // :
            // ModelPhysics<nDim>( "body" )
            // {}
        explicit Body( std::shared_ptr<ModelPhysics<nRealDim>> const& mphysics )
            :
            M_modelPhysics( mphysics ),
            M_mass( 0 )
            {}
        Body( Body const& ) = default;
        Body( Body && ) = default;

        void setup( pt::ptree const& p, ModelMaterials const& mats, mesh_ptrtype mesh );

        void updateForUse();

        bool hasMaterialsProperties() const { return (M_materialsProperties? true : false); }

        void setMass( double m ) { M_mass = m; }
        void setMomentOfInertia( moment_of_inertia_type const& m ) { M_momentOfInertia = m; }
        void setMomentOfInertia( double val ) { M_momentOfInertia = val*moment_of_inertia_type::Identity(); }
        void setMassCenter( eigen_vector_type<nRealDim> const& massCenter ) { M_massCenter = massCenter; }
        double mass() const { return M_mass; }
        moment_of_inertia_type const& momentOfInertia() const { return M_momentOfInertia; }
        eigen_vector_type<nRealDim> const& massCenter() const { return M_massCenter; }

        auto massExpr() const { return cst( M_mass ); }
        auto momentOfInertiaExpr() const
            {
                if constexpr ( nDim == 2 )
                    return cst(M_momentOfInertia(0,0));
                else
                    return mat<3,3>( cst(M_momentOfInertia(0,0)),cst(M_momentOfInertia(0,1)),cst(M_momentOfInertia(0,2)),
                                     cst(M_momentOfInertia(1,0)),cst(M_momentOfInertia(1,1)),cst(M_momentOfInertia(1,2)),
                                     cst(M_momentOfInertia(2,0)),cst(M_momentOfInertia(2,1)),cst(M_momentOfInertia(2,2)) );
            }
        auto massCenterExpr() const
            {
                if constexpr ( nDim == 2 )
                    return vec( cst(M_massCenter(0)), cst(M_massCenter(1)) );
                else
                    return vec( cst(M_massCenter(0)), cst(M_massCenter(1)), cst(M_massCenter(2)) );
            }


        template <typename ExprType>
        double evaluateMassFromDensity( Expr<ExprType> const& densityExpr ) const
            {
                CHECK( M_materialsProperties ) << "no materialsProperties defined";
                auto mom = M_materialsProperties->materialsOnMesh(M_mesh);
                double mass = 0;
                for ( auto const& rangeData : mom->rangeMeshElementsByMaterial() )
                {
                    auto const& range = rangeData.second;
                    mass += integrate(_range=range,_expr=densityExpr).evaluate()(0,0);
                }
                return mass;
            }

        void setParameterValues( std::map<std::string,double> const& mp )
            {
                if ( M_materialsProperties )
                    M_materialsProperties->setParameterValues( mp );
            }

    private :
        std::shared_ptr<ModelPhysics<nRealDim>> M_modelPhysics;
        mesh_ptrtype M_mesh;
        materialsproperties_ptrtype M_materialsProperties;
        eigen_vector_type<nRealDim> M_massCenter;//, M_massCenterRef;
        double M_mass;
        moment_of_inertia_type M_momentOfInertia;
    };

    // bc body
    class BodyBoundaryCondition
    {
    public :
        typedef typename mpl::if_< mpl::equal_to<mpl::int_<nDim>,mpl::int_<2> >,
                                   space_trace_p0c_scalar_type,
                                   space_trace_p0c_vectorial_type >::type space_trace_angular_velocity_type;
        typedef std::shared_ptr<space_trace_angular_velocity_type> space_trace_angular_velocity_ptrtype;
        typedef typename space_trace_angular_velocity_type::element_type element_trace_angular_velocity_type;
        typedef std::shared_ptr<element_trace_angular_velocity_type> element_trace_angular_velocity_ptrtype;
        typedef Bdf<space_trace_angular_velocity_type> bdf_trace_angular_velocity_type;
        typedef std::shared_ptr<bdf_trace_angular_velocity_type> bdf_trace_angular_velocity_ptrtype;

        struct FieldTag
        {
            static auto translational_velocity( BodyBoundaryCondition const* t ) { return ModelFieldTag<BodyBoundaryCondition,0>( t ); }
            static auto angular_velocity( BodyBoundaryCondition const* t ) { return ModelFieldTag<BodyBoundaryCondition,1>( t ); }
        };

        BodyBoundaryCondition();
        BodyBoundaryCondition( BodyBoundaryCondition const& ) = default;
        BodyBoundaryCondition( BodyBoundaryCondition && ) = default;

        void setup( std::string const& bodyName, pt::ptree const& p, self_type const& fluidToolbox );
        void init( self_type const& fluidToolbox );
        void updateForUse( self_type const& fluidToolbox );


        void initTimeStep( self_type const& fluidToolbox, int bdfOrder, int nConsecutiveSave, std::string const& myFileFormat )
            {
                M_bdfTranslationalVelocity = fluidToolbox.createBdf( M_XhTranslationalVelocity, "body."+M_name+".translational-velocity", bdfOrder, nConsecutiveSave, myFileFormat );
                M_bdfAngularVelocity = fluidToolbox.createBdf( M_XhAngularVelocity, "body."+M_name+".angular-velocity", bdfOrder, nConsecutiveSave, myFileFormat );

                if ( fluidToolbox.doRestart() )
                {
                    M_bdfTranslationalVelocity->restart();
                    M_bdfAngularVelocity->restart();
                    *M_fieldTranslationalVelocity = M_bdfTranslationalVelocity->unknown(0);
                    *M_fieldAngularVelocity = M_bdfAngularVelocity->unknown(0);
                }
            }

        void startTimeStep()
            {
                M_bdfTranslationalVelocity->start( *M_fieldTranslationalVelocity );
                M_bdfAngularVelocity->start( *M_fieldAngularVelocity );
            }
        void updateTimeStep()
            {
                M_bdfTranslationalVelocity->next( *M_fieldTranslationalVelocity );
                M_bdfAngularVelocity->next( *M_fieldAngularVelocity );
            }

        std::string const& name() const { return M_name; }

        range_faces_type const& rangeMarkedFacesOnFluid() const { return M_rangeMarkedFacesOnFluid; }
        trace_mesh_ptrtype mesh() const { return M_mesh; }

        std::set<std::string>/*ModelMarkers*/ const& markers() const { return M_markers; }

        space_trace_p0c_vectorial_ptrtype spaceTranslationalVelocity() const { return M_XhTranslationalVelocity; }
        space_trace_angular_velocity_ptrtype spaceAngularVelocity() const { return M_XhAngularVelocity; }
        element_trace_p0c_vectorial_ptrtype fieldTranslationalVelocityPtr() const { return M_fieldTranslationalVelocity; }
        element_trace_angular_velocity_ptrtype fieldAngularVelocityPtr() const { return M_fieldAngularVelocity; }

        bdf_trace_p0c_vectorial_ptrtype bdfTranslationalVelocity() const { return M_bdfTranslationalVelocity; }
        bdf_trace_angular_velocity_ptrtype bdfAngularVelocity() const { return M_bdfAngularVelocity; }

        Body const& body() const { return *M_body; }
        auto massExpr() const { return M_body->massExpr(); }
        auto momentOfInertiaExpr() const { return M_body->momentOfInertiaExpr(); }
        auto massCenterExpr() const
            {
                return M_body->massCenterExpr();
            }

        bool hasTranslationalVelocityExpr() const { return M_translationalVelocityExpr.template hasExpr<nDim,1>(); }
        auto const& translationalVelocityExpr() const { return M_translationalVelocityExpr.template expr<nDim,1>(); }
        bool hasAngularVelocityExpr() const {
            if constexpr ( nDim == 2 )
                return M_angularVelocityExpr.template hasExpr<1,1>();
            else
                return M_angularVelocityExpr.template hasExpr<nDim,1>();
        }
        auto const& angularVelocityExpr() const
            {
                if constexpr ( nDim == 2 )
                    return M_angularVelocityExpr.expr<1,1>();
                else
                    return M_angularVelocityExpr.expr<nDim,1>();
            }

        auto rigidVelocityExpr() const
            {
                if constexpr ( nDim == 2 )
                    return this->translationalVelocityExpr() + this->angularVelocityExpr()*vec(-Py()+this->massCenterExpr()(1,0),Px()-this->massCenterExpr()(0,0) );
                else
                    return this->translationalVelocityExpr() + cross( this->angularVelocityExpr(), P()-this->massCenterExpr() );
            }
        auto rigidVelocityExprFromFields() const
            {
                if constexpr ( nDim == 2 )
                    return idv(M_fieldTranslationalVelocity) + idv(M_fieldAngularVelocity)*vec(-Py()+this->massCenterExpr()(1,0),Px()-this->massCenterExpr()(0,0) );
                else
                    return idv(M_fieldTranslationalVelocity) + cross( idv(M_fieldAngularVelocity), P()-this->massCenterExpr() );
            }

        sparse_matrix_ptrtype matrixPTilde_translational() const { return M_matrixPTilde_translational; }
        sparse_matrix_ptrtype matrixPTilde_angular() const { return M_matrixPTilde_angular; }


        //---------------------------------------------------------------------------//
        // elastic velocity
        //---------------------------------------------------------------------------//
        bool hasElasticVelocity() const { return ( M_fieldElasticVelocity? true : false ); }

        bool hasElasticVelocityFromExpr() const { return !M_elasticVelocityExprBC.empty(); }

        element_trace_velocity_ptrtype fieldElasticVelocityPtr() const { return M_fieldElasticVelocity; }

        auto elasticVelocityExpr() const { CHECK( this->hasElasticVelocity() ) << "no elastic velocity"; return idv(M_fieldElasticVelocity); }

        void updateElasticVelocityFromExpr( self_type const& fluidToolbox );

        //---------------------------------------------------------------------------//
        // gravity
        bool gravityForceEnabled() const { return M_gravityForceEnabled; }
        //double massOfFluid() const { return M_massOfFluid; }
        eigen_vector_type<nRealDim> const& gravityForceWithMass() const { return M_gravityForceWithMass; }

        //---------------------------------------------------------------------------//

        void setParameterValues( std::map<std::string,double> const& mp )
            {
                M_translationalVelocityExpr.setParameterValues( mp );
                M_angularVelocityExpr.setParameterValues( mp );
                for ( auto & [bcName,eve] : M_elasticVelocityExprBC )
                    std::get<0>( eve ).setParameterValues( mp );
                if ( M_body )
                    M_body->setParameterValues( mp );
            }

    private :
        std::string M_name;
        ModelMarkers M_markers;
        range_faces_type M_rangeMarkedFacesOnFluid;
        trace_mesh_ptrtype M_mesh;
        space_trace_p0c_vectorial_ptrtype M_XhTranslationalVelocity;
        space_trace_angular_velocity_ptrtype M_XhAngularVelocity;
        element_trace_p0c_vectorial_ptrtype M_fieldTranslationalVelocity;
        element_trace_angular_velocity_ptrtype M_fieldAngularVelocity;
        bdf_trace_p0c_vectorial_ptrtype M_bdfTranslationalVelocity;
        bdf_trace_angular_velocity_ptrtype M_bdfAngularVelocity;
        sparse_matrix_ptrtype M_matrixPTilde_translational, M_matrixPTilde_angular;
        ModelExpression M_translationalVelocityExpr, M_angularVelocityExpr;

        std::shared_ptr<Body> M_body;
        eigen_vector_type<nRealDim> M_massCenterRef;

        space_trace_velocity_ptrtype M_XhElasticVelocity;
        element_trace_velocity_ptrtype M_fieldElasticVelocity;
        std::map<std::string, std::tuple< ModelExpression, std::set<std::string>>> M_elasticVelocityExprBC;

        bool M_gravityForceEnabled;
        //double M_massOfFluid;
        eigen_vector_type<nRealDim> M_gravityForceWithMass;
    };

    class BodySetBoundaryCondition : public std::map<std::string,BodyBoundaryCondition>
    {
    public:
        void initTimeStep( self_type const& fluidToolbox, int bdfOrder, int nConsecutiveSave, std::string const& myFileFormat )
            {
                for ( auto & [name,bpbc] : *this )
                    bpbc.initTimeStep( fluidToolbox, bdfOrder, nConsecutiveSave, myFileFormat );
            }
        void startTimeStep()
            {
                for ( auto & [name,bpbc] : *this )
                    bpbc.startTimeStep();
            }
        void updateTimeStep()
            {
                for ( auto & [name,bpbc] : *this )
                    bpbc.updateTimeStep();
            }
        void updateForUse( self_type const& fluidToolbox );
        void updateAlgebraicFactoryForUse( self_type const& fluidToolbox, model_algebraic_factory_ptrtype algebraicFactory );

        void init( self_type const& fluidToolbox )
            {
                for ( auto & [name,bpbc] : *this )
                    bpbc.init( fluidToolbox );
            }
        void setParameterValues( std::map<std::string,double> const& mp )
            {
                for ( auto & [name,bpbc] : *this )
                    bpbc.setParameterValues( mp );
            }
        bool hasTranslationalVelocityExpr() const
            {
                for ( auto const& [name,bpbc] : *this )
                    if ( bpbc.hasTranslationalVelocityExpr() )
                        return true;
                return false;
            }
        bool hasAngularVelocityExpr() const
            {
                for ( auto const& [name,bpbc] : *this )
                    if ( bpbc.hasAngularVelocityExpr() )
                        return true;
                return false;
            }
        bool hasElasticVelocity() const
            {
                for ( auto const& [name,bpbc] : *this )
                    if ( bpbc.hasElasticVelocity() )
                        return true;
                return false;
            }
        bool hasElasticVelocityFromExpr() const
            {
                for ( auto const& [name,bpbc] : *this )
                    if ( bpbc.hasElasticVelocityFromExpr() )
                        return true;
                return false;
            }

        auto modelFields( self_type const& fluidToolbox, std::string const& prefix = "" ) const
            {
                using _field_translational_ptrtype = std::decay_t<decltype(this->begin()->second.fieldTranslationalVelocityPtr())>;
                using _field_angular_ptrtype = std::decay_t<decltype(this->begin()->second.fieldAngularVelocityPtr())>;

                std::map<std::string,std::tuple<_field_translational_ptrtype,_field_angular_ptrtype>> registerFields;
                for ( auto const& [name,bpbc] : *this )
                {
                    registerFields[name] = std::make_tuple( bpbc.fieldTranslationalVelocityPtr(), bpbc.fieldAngularVelocityPtr() );
                }
                return this->modelFieldsImpl( fluidToolbox,registerFields,prefix );
            }
        auto modelFields( self_type const& fluidToolbox, vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
            {
                using _field_translational_ptrtype = std::decay_t<decltype( this->begin()->second.spaceTranslationalVelocity()->elementPtr( *sol,rowStartInVector ) )>;
                using _field_angular_ptrtype = std::decay_t<decltype(this->begin()->second.spaceAngularVelocity()->elementPtr( *sol, rowStartInVector ) )>;

                std::map<std::string,std::tuple<_field_translational_ptrtype,_field_angular_ptrtype>> registerFields;
                for ( auto const& [name,bpbc] : *this )
                {
                    size_type startBlockIndexTranslationalVelocity = fluidToolbox.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".translational-velocity");
                    size_type startBlockIndexAngularVelocity = fluidToolbox.startSubBlockSpaceIndex("body-bc."+bpbc.name()+".angular-velocity");
                    registerFields[name] = std::make_tuple( bpbc.spaceTranslationalVelocity()->elementPtr( *sol, rowStartInVector+startBlockIndexTranslationalVelocity ),
                                                            bpbc.spaceAngularVelocity()->elementPtr( *sol, rowStartInVector+startBlockIndexAngularVelocity ) );
                }
                return this->modelFieldsImpl( fluidToolbox,registerFields,prefix );
            }
    private:

        template <typename _field_translational_ptrtype, typename _field_angular_ptrtype>
        auto modelFieldsImpl( self_type const& fluidToolbox, std::map<std::string,std::tuple<_field_translational_ptrtype,_field_angular_ptrtype>> const& registerFields, std::string const& prefix ) const
            {
                auto mfieldTranslational = modelField<FieldCtx::ID,_field_translational_ptrtype>( BodyBoundaryCondition::FieldTag::translational_velocity(nullptr) );
                auto mfieldAngular = modelField<FieldCtx::ID,_field_angular_ptrtype>( BodyBoundaryCondition::FieldTag::angular_velocity(nullptr) );
                for ( auto const& [name,bpbc] : *this )
                {
                    auto const& field_translational = std::get<0>( registerFields.find( name )->second );
                    auto const& field_angular = std::get<1>( registerFields.find( name )->second );
                    std::string prefixBase = prefixvm( prefix, (boost::format("body.%1%")%name).str() );
                    std::string prefix_symbol = prefixvm( fluidToolbox.keyword(), (boost::format("body_%1%")%name).str(), "_");
                    mfieldTranslational.add( BodyBoundaryCondition::FieldTag::translational_velocity(&bpbc), prefixBase, "translational-velocity", field_translational, "V",  prefix_symbol );
                    mfieldAngular.add( BodyBoundaryCondition::FieldTag::angular_velocity(&bpbc), prefixBase, "angular-velocity", field_angular, "W",  prefix_symbol );
                }
                return Feel::FeelModels::modelFields( mfieldTranslational, mfieldAngular );
            }
    };


    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    // export
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    typedef Exporter<trace_mesh_type,nOrderGeo> export_trace_type;
    typedef std::shared_ptr<export_trace_type> export_trace_ptrtype;
    //typedef Exporter<mesh_type,nOrderGeo> gmsh_export_type;
    //typedef std::shared_ptr<gmsh_export_type> gmsh_export_ptrtype;
    //___________________________________________________________________________________//
    // export ho
#if 1 //defined(FEELPP_HAS_VTK)
    //fais comme ca car bug dans opeartorlagrangeP1 pour les champs vectorielles
    typedef FunctionSpace<mesh_type,bases<Lagrange<nOrderVelocity,Scalar,Continuous,PointSetFekete> > > space_create_ho_type;
    // mesh
    typedef Mesh<Simplex<nDim,1,nDim> > mesh_visu_ho_type;
    //function space vectorial
    typedef FunctionSpace<mesh_visu_ho_type,bases<Lagrange<1,Vectorial,Continuous,PointSetFekete> > > space_vectorial_visu_ho_type;
    typedef std::shared_ptr<space_vectorial_visu_ho_type> space_vectorial_visu_ho_ptrtype;
    typedef typename space_vectorial_visu_ho_type::element_type element_vectorial_visu_ho_type;
    typedef std::shared_ptr<element_vectorial_visu_ho_type> element_vectorial_visu_ho_ptrtype;
    // function space scalar
    //typedef FunctionSpace<mesh_visu_ho_type,bases<Lagrange<1,Scalar,Continuous,PointSetFekete> > > space_scalar_visu_ho_type;
    typedef typename space_vectorial_visu_ho_type::component_functionspace_type space_scalar_visu_ho_type;
    typedef std::shared_ptr<space_scalar_visu_ho_type> space_scalar_visu_ho_ptrtype;
    typedef typename space_scalar_visu_ho_type::element_type element_scalar_visu_ho_type;
    typedef std::shared_ptr<element_scalar_visu_ho_type> element_scalar_visu_ho_ptrtype;
    // function space vectorial discontinuos
    typedef FunctionSpace<mesh_visu_ho_type,bases<Lagrange<1, Vectorial,Discontinuous,PointSetFekete> > > space_vectorialdisc_visu_ho_type;
    typedef std::shared_ptr<space_vectorialdisc_visu_ho_type> space_vectorialdisc_visu_ho_ptrtype;
    typedef typename space_vectorialdisc_visu_ho_type::element_type element_vectorialdisc_visu_ho_type;
    typedef std::shared_ptr<element_vectorialdisc_visu_ho_type> element_vectorialdisc_visu_ho_ptrtype;
    //___________________________________________________________________________________//
    //
    // typedef boost::tuple<boost::mpl::size_t<MESH_ELEMENTS>,
    //                      typename MeshTraits<mesh_visu_ho_type>::element_const_iterator,
    //                      typename MeshTraits<mesh_visu_ho_type>::element_const_iterator> range_visu_ho_type;
    //___________________________________________________________________________________//

    typedef OperatorInterpolation<space_velocity_type,
                                  space_vectorial_visu_ho_type > op_interpolation_visu_ho_vectorial_type;
    typedef std::shared_ptr<op_interpolation_visu_ho_vectorial_type> op_interpolation_visu_ho_vectorial_ptrtype;

    typedef OperatorInterpolation<space_pressure_type,
                                  space_scalar_visu_ho_type> op_interpolation_visu_ho_scalar_type;
    typedef std::shared_ptr<op_interpolation_visu_ho_scalar_type> op_interpolation_visu_ho_scalar_ptrtype;

#if defined( FEELPP_MODELS_HAS_MESHALE )
    typedef OperatorInterpolation<space_mesh_disp_type,
                                  space_vectorial_visu_ho_type> op_interpolation_visu_ho_meshdisp_type;
    typedef std::shared_ptr<op_interpolation_visu_ho_meshdisp_type> op_interpolation_visu_ho_meshdisp_ptrtype;
#endif

#if 0
    typedef OperatorInterpolation<space_normalstress_type,
                                  space_vectorialdisc_visu_ho_type/*,
                                                                   range_visu_ho_type*/> op_interpolation_visu_ho_vectorialdisc_type;
    typedef std::shared_ptr<op_interpolation_visu_ho_vectorialdisc_type> op_interpolation_visu_ho_vectorialdisc_ptrtype;
#endif
    //___________________________________________________________________________________//

    typedef Exporter<mesh_visu_ho_type> export_ho_type;
    typedef std::shared_ptr<export_ho_type> export_ho_ptrtype;
#endif

    // measure tools for points evaluation
    typedef MeasurePointsEvaluation<space_velocity_type,space_pressure_type> measure_points_evaluation_type;
    typedef std::shared_ptr<measure_points_evaluation_type> measure_points_evaluation_ptrtype;

    using force_type = Eigen::Matrix<typename super_type::value_type, nDim, 1, Eigen::ColMajor>;
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//
    //___________________________________________________________________________________//

    //___________________________________________________________________________________//
    // constructor
    explicit FluidMechanics( std::string const& prefix,
                             std::string const& keyword = "fluid",
                             worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr(),
                             std::string const& subPrefix = "",
                             ModelBaseRepository const& modelRep = ModelBaseRepository() );
    FluidMechanics( self_type const & M ) = default;

    static self_ptrtype New( std::string const& prefix,
                             std::string const& keyword = "fluid",
                             worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                             std::string const& subPrefix = "",
                             ModelBaseRepository const& modelRep = ModelBaseRepository() );
    //___________________________________________________________________________________//

    static std::string expandStringFromSpec( std::string const& expr );

private :
    void loadParameterFromOptionsVm();
    void initMesh();
    void createFunctionSpaces();
    void createALE();
    void initBoundaryConditions();
    void initFluidInlet();
    void initFluidOutlet();
    void initUserFunctions();
    void initPostProcess() override;
    void createPostProcessExporters();
public :
    void init( bool buildModelAlgebraicFactory=true );
    void initAlgebraicFactory();

    void createFunctionSpacesNormalStress();
    void createFunctionSpacesVorticity();
    void createFunctionSpacesSourceAdded();

    void loadMesh(mesh_ptrtype __mesh );
    void setMesh( mesh_ptrtype const& mesh ) { M_mesh = mesh; }

    void updateMarkedZonesInMesh();

    std::shared_ptr<std::ostringstream> getInfo() const override;
    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"FluidMechanicsMesh.path"); }
    //___________________________________________________________________________________//

    mesh_ptrtype const& mesh() const { return M_mesh; }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }
    std::shared_ptr<RangeDistributionByMaterialName<mesh_type> > rangeDistributionByMaterialName() const { return M_rangeDistributionByMaterialName; }

    space_velocity_ptrtype const& functionSpaceVelocity() const { return M_XhVelocity; }
    space_pressure_ptrtype const& functionSpacePressure() const { return M_XhPressure; }

    element_velocity_type & fieldVelocity() { return *M_fieldVelocity; }
    element_velocity_type const& fieldVelocity() const { return *M_fieldVelocity; }
    element_velocity_ptrtype & fieldVelocityPtr() { return M_fieldVelocity; }
    element_velocity_ptrtype const& fieldVelocityPtr() const { return M_fieldVelocity; }
    element_pressure_type & fieldPressure() { return *M_fieldPressure; }
    element_pressure_type const& fieldPressure() const { return *M_fieldPressure; }
    element_pressure_ptrtype const& fieldPressurePtr() const { return M_fieldPressure; }

    element_velocity_ptrtype const& fieldConvectionVelocityExtrapolatedPtr() const { return M_fieldConvectionVelocityExtrapolated; }

    element_normalstress_ptrtype & fieldNormalStressPtr() { return M_fieldNormalStress; }
    element_normalstress_ptrtype const& fieldNormalStressPtr() const { return M_fieldNormalStress; }
    element_normalstress_type const& fieldNormalStress() const { return *M_fieldNormalStress; }
    element_normalstress_ptrtype & fieldWallShearStressPtr() { return M_fieldWallShearStress; }
    element_normalstress_ptrtype const& fieldWallShearStressPtr() const { return M_fieldWallShearStress; }
    element_normalstress_type const& fieldWallShearStress() const { return *M_fieldWallShearStress; }

    element_vorticity_ptrtype const& fieldVorticityPtr() const { return M_fieldVorticity; }
    element_vorticity_ptrtype & fieldVorticityPtr() { return M_fieldVorticity; }
    element_vorticity_type const& fieldVorticity() const {  CHECK( M_fieldVorticity ) << "fieldVorticity not init"; return *M_fieldVorticity; }
    element_vorticity_type & fieldVorticity() {  CHECK( M_fieldVorticity ) << "fieldVorticity not init";return *M_fieldVorticity; }

    bool useExtendedDofTable() const;

    // fields defined by user (in json or external to this class)
    std::map<std::string,component_element_velocity_ptrtype> const& fieldsUserScalar() const { return M_fieldsUserScalar; }
    std::map<std::string,element_velocity_ptrtype> const& fieldsUserVectorial() const { return M_fieldsUserVectorial; }
    bool hasFieldUserScalar( std::string const& key ) const { return M_fieldsUserScalar.find( key ) != M_fieldsUserScalar.end(); }
    bool hasFieldUserVectorial( std::string const& key ) const { return M_fieldsUserVectorial.find( key ) != M_fieldsUserVectorial.end(); }
    component_element_velocity_ptrtype const& fieldUserScalarPtr( std::string const& key ) const {
        CHECK( this->hasFieldUserScalar( key ) ) << "field name " << key << " not registered"; return M_fieldsUserScalar.find( key )->second; }
    element_velocity_ptrtype const& fieldUserVectorialPtr( std::string const& key ) const {
        CHECK( this->hasFieldUserVectorial( key ) ) << "field name " << key << " not registered"; return M_fieldsUserVectorial.find( key )->second; }
    component_element_velocity_type const& fieldUserScalar( std::string const& key ) const { return *this->fieldUserScalarPtr( key ); }
    element_velocity_type const& fieldUserVectorial( std::string const& key ) const { return *this->fieldUserVectorialPtr( key ); }

    void registerCustomFieldScalar( std::string const& name )
        {
            if ( M_fieldsUserScalar.find( name ) == M_fieldsUserScalar.end() )
                M_fieldsUserScalar[name];
        }
    void registerCustomFieldVectorial( std::string const& name )
        {
            if ( M_fieldsUserVectorial.find( name ) == M_fieldsUserVectorial.end() )
                M_fieldsUserVectorial[name];
        }
    template <typename ExprT>
    void updateCustomField( std::string const& name, vf::Expr<ExprT> const& e )
        {
            this->updateCustomField( name, e, M_rangeMeshElements );
        }
    template <typename ExprT, typename OnRangeType>
    void updateCustomField( std::string const& name, vf::Expr<ExprT> const& e, OnRangeType const& range, std::enable_if_t< ExprTraits<OnRangeType, vf::Expr<ExprT>>::shape::is_scalar>* = nullptr )
        {
            if ( M_fieldsUserScalar.find( name ) == M_fieldsUserScalar.end() || !M_fieldsUserScalar[name] )
                M_fieldsUserScalar[name] = this->functionSpaceVelocity()->compSpace()->elementPtr();
             M_fieldsUserScalar[name]->on(_range=range,_expr=e );
        }
    template <typename ExprT, typename OnRangeType>
    void updateCustomField( std::string const& name, vf::Expr<ExprT> const& e, OnRangeType const& range, std::enable_if_t< ExprTraits<OnRangeType, vf::Expr<ExprT>>::shape::is_vectorial>* = nullptr )
        {
            if ( M_fieldsUserVectorial.find( name ) == M_fieldsUserVectorial.end() || !M_fieldsUserVectorial[name] )
                M_fieldsUserVectorial[name] = this->functionSpaceVelocity()->elementPtr();
            M_fieldsUserVectorial[name]->on(_range=range,_expr=e );
        }

    //___________________________________________________________________________________//
    // algebraic data
    backend_ptrtype backend() { return M_backend; }
    backend_ptrtype const& backend() const { return  M_backend; }
    typename super_type::block_pattern_type blockPattern() const override;
    virtual BlocksBaseGraphCSR buildBlockMatrixGraph() const override;
    graph_ptrtype buildMatrixGraph() const override;
    virtual int nBlockMatrixGraph() const;
    model_algebraic_factory_ptrtype algebraicFactory() { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    virtual size_type nLocalDof() const;
    void buildBlockVector();
    BlocksBaseVector<double> blockVectorSolution() { return M_blockVectorSolution; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    void updateBlockVectorSolution();

    //___________________________________________________________________________________//
    // time step scheme
    std::string const& timeStepping() const { return M_timeStepping; }
    bdf_velocity_ptrtype timeStepBDF() { return M_bdfVelocity; }
    bdf_velocity_ptrtype const& timeStepBDF() const { return M_bdfVelocity; }
    std::shared_ptr<TSBase> timeStepBase() { return this->timeStepBDF(); }
    std::shared_ptr<TSBase> timeStepBase() const { return this->timeStepBDF(); }
    void initTimeStep();
    void startTimeStep();
    void updateTimeStep();

    // init/update user functions defined in json
    void updateUserFunctions( bool onlyExprWithTimeSymbol = false );

    // post process
    void exportResults() { this->exportResults( this->currentTime() ); }
    void exportResults( double time );
    template <typename SymbolsExpr>
    void exportResults( double time, SymbolsExpr const& symbolsExpr );

    // void exportFields( double time );
    // bool updateExportedFields( export_ptrtype exporter, std::set<std::string> const& fields, double time );
    // bool updateExportedFieldsOnTrace( export_trace_ptrtype exporter, std::set<std::string> const& fields, double time );
    void setDoExport(bool b);
private :
    void executePostProcessMeasures( double time );
    template <typename TupleFieldsType,typename SymbolsExpr>
    void executePostProcessMeasures( double time, TupleFieldsType const& tupleFields, SymbolsExpr const& symbolsExpr );
    void updateConvectionVelocityExtrapolated();
    void updateTimeStepCurrentResidual();
    //void exportResultsImpl( double time );
    void exportResultsImplHO( double time );
public :
    //___________________________________________________________________________________//
    // ale mesh
#if defined( FEELPP_MODELS_HAS_MESHALE )
    mesh_ale_ptrtype meshALE() { return M_meshALE; }
    mesh_ale_ptrtype const& meshALE() const { return M_meshALE; }
    element_mesh_disp_ptrtype meshDisplacementOnInterface() { return M_meshDisplacementOnInterface; }
    element_meshvelocity_type & meshVelocity() { return *M_meshALE->velocity(); }
    element_meshvelocity_type const & meshVelocity() const { return *M_meshALE->velocity(); }
#endif
    //___________________________________________________________________________________//

    bool applyMovingMeshBeforeSolve() const { return M_applyMovingMeshBeforeSolve; }
    void setApplyMovingMeshBeforeSolve( bool b ) { M_applyMovingMeshBeforeSolve = b; }
    bool isMoveDomain() const { return M_isMoveDomain; }
#if 0
    std::string const& modelName() const;
    void setModelName( std::string const& type );
#endif
    std::string const& solverName() const;
    void setSolverName( std::string const& type );

    bool isStationaryModel() const;

    bool startBySolveNewtonian() const { return M_startBySolveNewtonian; }
    void startBySolveNewtonian( bool b ) { M_startBySolveNewtonian=b; }
    bool hasSolveNewtonianAtKickOff() const { return M_hasSolveNewtonianAtKickOff; }
    void hasSolveNewtonianAtKickOff( bool b ) { M_hasSolveNewtonianAtKickOff=b; }

    bool startBySolveStokesStationary() const { return M_startBySolveStokesStationary; }
    void startBySolveStokesStationary( bool b ) { M_startBySolveStokesStationary=b; }
    bool hasSolveStokesStationaryAtKickOff() const { return M_hasSolveStokesStationaryAtKickOff; }
    void hasSolveStokesStationaryAtKickOff( bool b ) { M_hasSolveStokesStationaryAtKickOff=b; }

    //___________________________________________________________________________________//
    // fsi parameters

    bool useFSISemiImplicitScheme() const { return M_useFSISemiImplicitScheme; }
    void useFSISemiImplicitScheme(bool b) { M_useFSISemiImplicitScheme=b; }
    /*FEELPP_DEPRECATED*/ std::string couplingFSIcondition() const { return M_couplingFSIcondition; }
    /*FEELPP_DEPRECATED*/ void couplingFSIcondition(std::string s) { M_couplingFSIcondition=s; }

    std::set<std::string> const& markersFSI() const { return M_markersFSI; }
    //___________________________________________________________________________________//
    // stabilization
    bool stabilizationGLS() const { return M_stabilizationGLS; }
    std::string const& stabilizationGLSType() const { return M_stabilizationGLSType; }
    stab_gls_parameter_ptrtype const& stabilizationGLSParameterConvectionDiffusion() const { return M_stabilizationGLSParameterConvectionDiffusion; }
    stab_gls_parameter_ptrtype const& stabilizationGLSParameterPressure() const { return M_stabilizationGLSParameterPressure; }
    range_elements_type const& stabilizationGLSEltRangeConvectionDiffusion( std::string const& matName ) const
        {
            auto itFind = M_stabilizationGLSEltRangeConvectionDiffusion.find( matName );
            CHECK( itFind != M_stabilizationGLSEltRangeConvectionDiffusion.end() ) << "not found with material matName";
            return itFind->second;
        }
    range_elements_type const& stabilizationGLSEltRangePressure( std::string const& matName ) const
        {
            auto itFind = M_stabilizationGLSEltRangePressure.find( matName );
            CHECK( itFind != M_stabilizationGLSEltRangePressure.end() ) << "not found with material matName";
            return itFind->second;
        }
    void setStabilizationGLSDoAssembly( bool b) { M_stabilizationGLSDoAssembly = b; }

    bool applyCIPStabOnlyOnBoundaryFaces() const { return M_applyCIPStabOnlyOnBoundaryFaces; }
    void applyCIPStabOnlyOnBoundaryFaces(bool b) { M_applyCIPStabOnlyOnBoundaryFaces=b; }
    bool doCIPStabConvection() const { return M_doCIPStabConvection; }
    void doCIPStabConvection(bool b) { M_doCIPStabConvection=b; }
    bool doCIPStabDivergence() const { return M_doCIPStabDivergence; }
    void doCIPStabDivergence(bool b) { M_doCIPStabDivergence=b; }
    bool doCIPStabPressure() const { return M_doCIPStabPressure; }
    void doCIPStabPressure(bool b) { M_doCIPStabPressure=b; }
    double stabCIPConvectionGamma() const { return M_stabCIPConvectionGamma; }
    double stabCIPDivergenceGamma() const { return M_stabCIPDivergenceGamma; }
    double stabCIPPressureGamma() const { return M_stabCIPPressureGamma; }

    bool doStabDivDiv() const { return M_doStabDivDiv; }
    void doStabDivDiv(bool b) { M_doStabDivDiv=b; }

    bool doStabConvectionEnergy() const { return M_doStabConvectionEnergy; }
    void doStabConvectionEnergy(bool b) { M_doStabConvectionEnergy=b; }

    bool definePressureCst() const { return M_definePressureCst; }
    void setDefinePressureCst(bool b) { M_definePressureCst = b; }
    std::string definePressureCstMethod() const { return M_definePressureCstMethod; }
    void setDefinePressureCstMethod(std::string s) { M_definePressureCstMethod = s; }
    double definePressureCstPenalisationBeta() const { return M_definePressureCstPenalisationBeta; }

    void updateDefinePressureCst();

    //___________________________________________________________________________________//
    // physical parameters rho,mu,nu,...
    material_properties_ptrtype & materialProperties() { return M_materialProperties; }
    material_properties_ptrtype const& materialProperties() const { return M_materialProperties; }

    void updateRho(double rho)
    {
        this->materialProperties()->setCstDensity(rho);
    }
    void updateMu(double mu)
    {
        this->materialProperties()->setCstDynamicViscosity(mu);
        M_pmmNeedUpdate = true;
    }
    template < typename ExprT >
    void updateRho(vf::Expr<ExprT> const& __expr)
    {
        this->materialProperties()->updateDensityField( __expr );
    }
    template < typename ExprT >
    void updateMu(vf::Expr<ExprT> const& __expr)
    {
        this->materialProperties()->updateDynamicViscosityField( __expr );
        M_pmmNeedUpdate = true;
    }

    //___________________________________________________________________________________//
    // toolbox fields
    //___________________________________________________________________________________//

    auto modelFields( std::string const& prefix = "" ) const
        {
            return this->modelFields( this->fieldVelocityPtr(), this->fieldPressurePtr(), M_bodySetBC.modelFields( *this, prefix ), prefix );
        }
    auto modelFields( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto field_u = this->fieldVelocity().functionSpace()->elementPtr( *sol, rowStartInVector+this->startSubBlockSpaceIndex("velocity") );
            auto field_p = this->fieldPressure().functionSpace()->elementPtr( *sol, rowStartInVector+this->startSubBlockSpaceIndex("pressure") );
            auto mfields_body = M_bodySetBC.modelFields( *this, sol, rowStartInVector, prefix );
            return this->modelFields( field_u, field_p, mfields_body, prefix );
        }
    template <typename VelocityFieldType,typename PressureFieldType,typename ModelFieldsBodyType>
    auto modelFields( VelocityFieldType const& field_u, PressureFieldType const& field_p, ModelFieldsBodyType const& mfields_body, std::string const& prefix = "" ) const
        {
            auto mfields_ale = this->modelFieldsMeshALE( prefix );
            return Feel::FeelModels::modelFields( modelField<FieldCtx::ID|FieldCtx::MAGNITUDE/*|FieldCtx::GRAD|FieldCtx::GRAD_NORMAL*/>( FieldTag::velocity(this), prefix, "velocity", field_u, "U", this->keyword() ),
                                                  modelField<FieldCtx::ID>( FieldTag::pressure(this), prefix, "pressure", field_p, "P", this->keyword() ),
                                                  mfields_body, mfields_ale
                                                  );
        }

    auto trialSelectorModelFields( size_type startBlockSpaceIndex = 0 ) const
        {
            return Feel::FeelModels::selectorModelFields( selectorModelField( FieldTag::velocity(this), "velocity", startBlockSpaceIndex + this->startSubBlockSpaceIndex("velocity") ),
                                                          selectorModelField( FieldTag::pressure(this), "pressure", startBlockSpaceIndex + this->startSubBlockSpaceIndex("pressure") )
                                                          );
        }

    //___________________________________________________________________________________//
    // symbols expression
    //___________________________________________________________________________________//

    template <typename ModelFieldsType>
    auto symbolsExpr( ModelFieldsType const& mfields ) const
        {
            auto seToolbox = this->symbolsExprToolbox( mfields );
            auto seParam = this->symbolsExprParameter();
            //auto seMat = this->materialsProperties()->symbolsExpr();
            auto seFields = mfields.symbolsExpr();
            return Feel::vf::symbolsExpr( seToolbox, seParam/*, seMat*/, seFields );
        }
    auto symbolsExpr( std::string const& prefix = "" ) const { return this->symbolsExpr( this->modelFields( prefix ) ); }

    template <typename ModelFieldsType>
    auto symbolsExprToolbox( ModelFieldsType const& mfields ) const
        {
            return symbols_expression_empty_t{};
        }

    //___________________________________________________________________________________//
    // model context helper
    //___________________________________________________________________________________//

    // template <typename ModelFieldsType>
    // auto modelContext( ModelFieldsType const& mfields, std::string const& prefix = "" ) const
    //     {
    //         return Feel::FeelModels::modelContext( mfields, this->symbolsExpr( mfields ) );
    //     }
    auto modelContext( std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( prefix );
            return Feel::FeelModels::modelContext( std::move( mfields ), this->symbolsExpr( mfields ) );
        }
    auto modelContext( vector_ptrtype sol, size_type rowStartInVector = 0, std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields( sol, rowStartInVector, prefix );
            return Feel::FeelModels::modelContext( std::move( mfields ), this->symbolsExpr( mfields ) );
        }


    //___________________________________________________________________________________//
    // fields
    template <typename SymbolsExpr>
    void updateFields( SymbolsExpr const& symbolsExpr )
        {
            //this->materialProperties()->updateFields( symbolsExpr );
            this->updateVorticity();
        }
#if 0
    auto allFields( std::string const& prefix = "" ) const
        {
            std::map<std::string,element_normalstress_ptrtype> fields_normalstress;
            fields_normalstress[prefixvm(prefix,"trace.normal-stress")] = this->fieldNormalStressPtr();
            fields_normalstress[prefixvm(prefix,"trace.wall-shear-stress")] = this->fieldWallShearStressPtr();
#if 0
            std::map<std::string,element_trace_p0c_vectorial_ptrtype> fields_NoSlipRigidTranslationalVelocity;
            std::map<std::string,typename BodyBoundaryCondition::element_trace_angular_velocity_ptrtype> fields_NoSlipRigidAngularVelocity;
            if ( !M_bodySetBC.empty() )
            {
                fields_NoSlipRigidTranslationalVelocity[prefixvm(prefix,"trace.body.translational-velocity")]= M_bodySetBC.begin()->second.fieldTranslationalVelocityPtr();
                fields_NoSlipRigidAngularVelocity[prefixvm(prefix,"trace.body.angular-velocity")]=M_bodySetBC.begin()->second.fieldAngularVelocityPtr();
            }
#endif
            std::map<std::string, typename mesh_ale_type::ale_map_element_ptrtype> fields_disp;
            if ( this->isMoveDomain() )
                fields_disp[prefixvm(prefix,"displacement")] = this->meshALE()->displacement();
            return hana::make_tuple( std::make_pair( prefixvm( prefix,"velocity"),this->fieldVelocityPtr() ),
                                     std::make_pair( prefixvm( prefix,"pressure"),this->fieldPressurePtr() ),
                                     std::make_pair( prefixvm( prefix,"vorticity"),this->fieldVorticityPtr() ),
                                     fields_disp,
                                     fields_normalstress
                                     //,fields_NoSlipRigidTranslationalVelocity,fields_NoSlipRigidAngularVelocity
                                     );
        }
#endif
    //___________________________________________________________________________________//
    template <typename SymbExprType>
    auto exprPostProcessExports( SymbExprType const& se, std::string const& prefix = "" ) const
        {
            return hana::make_tuple();
        }

    //___________________________________________________________________________________//
    // boundary conditions + body forces
    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& paramValues );

    map_vector_field<nDim,1,2> const& bcDirichlet() const { return M_bcDirichlet; }
    map_vector_field<nDim,1,2>& bcDirichlet() { return M_bcDirichlet; }
    std::map<ComponentType,map_scalar_field<2> > const& bcDirichletComponents() const { return M_bcDirichletComponents; }
    std::map<ComponentType,map_scalar_field<2> > & bcDirichletComponents() { return M_bcDirichletComponents; }
    map_scalar_field<2> const& bcNeumannScalar() const { return M_bcNeumannScalar; }
    map_scalar_field<2> const& bcPressure() const { return M_bcPressure; }
    map_vector_field<nDim,1,2> const& bcNeumannVectorial() const { return M_bcNeumannVectorial; }
    map_matrix_field<nDim,nDim,2> const& bcNeumannTensor2() const { return M_bcNeumannTensor2; }
    map_vector_field<nDim,1,2> const& bodyForces() const { return M_volumicForcesProperties; }

    bool hasDirichletBC() const
        {
            return ( !M_bcDirichlet.empty() ||
                     !M_bcDirichletComponents.find(Component::X)->second.empty() ||
                     !M_bcDirichletComponents.find(Component::Y)->second.empty() ||
                     !M_bcDirichletComponents.find(Component::Z)->second.empty() );
        }

    
    // boundary conditions
    double dirichletBCnitscheGamma() const { return M_dirichletBCnitscheGamma; }
    void setDirichletBCnitscheGamma( double val) { M_dirichletBCnitscheGamma=val; }

    std::set<std::string> const& markersNameMovingBoundary() const { return this->markerALEMeshBC("moving"); }
    //___________________________________________________________________________________//
    // dirichlet with Lagrange multiplier
    trace_mesh_ptrtype const& meshDirichletLM() const { return M_meshDirichletLM; }
    space_trace_velocity_ptrtype const& XhDirichletLM() const { return M_XhDirichletLM; }
    //___________________________________________________________________________________//
    // impose mean pressure with P0 Lagrange multiplier
    space_meanpressurelm_ptrtype const& XhMeanPressureLM( int k ) const { return M_XhMeanPressureLM[k]; }
    //___________________________________________________________________________________//
    // fluid inlet bc
    bool hasFluidInlet() const { return !M_fluidInletDesc.empty(); }
    bool hasFluidInlet( std::string const& type ) const
    {
        for (auto const& inletbc : M_fluidInletDesc )
            if ( std::get<1>( inletbc ) == type )
                return true;
        return false;
    }
    void updateFluidInletVelocity();
    //___________________________________________________________________________________//
    // fluid outlets bc
    bool hasFluidOutlet() const { return !M_fluidOutletsBCType.empty(); }
    bool hasFluidOutletFree() const { return this->hasFluidOutlet("free"); }
    bool hasFluidOutletWindkessel() const { return this->hasFluidOutlet("windkessel"); }
    bool hasFluidOutlet(std::string const& type) const
    {
        for (auto const& outletbc : M_fluidOutletsBCType )
            if ( std::get<1>( outletbc ) == type )
                return true;
        return false;
    }
    bool hasFluidOutletWindkesselImplicit() const
    {
        for (auto const& outletbc : M_fluidOutletsBCType )
            if ( std::get<1>( outletbc ) == "windkessel" && std::get<0>( std::get<2>( outletbc ) ) == "implicit" )
                return true;
        return false;
    }
    bool hasFluidOutletWindkesselExplicit() const
    {
        for (auto const& outletbc : M_fluidOutletsBCType )
            if ( std::get<1>( outletbc ) == "windkessel" && std::get<0>( std::get<2>( outletbc ) ) == "explicit" )
                return true;
        return false;
    }
    int nFluidOutlet() const { return M_fluidOutletsBCType.size(); }
    int nFluidOutletWindkesselImplicit() const
    {
        int res=0;
        for (auto const& outletbc : M_fluidOutletsBCType )
            if ( std::get<1>( outletbc ) == "windkessel" && std::get<0>( std::get<2>( outletbc ) ) == "implicit" )
                ++res;
        return res;
    }
    std::map<int,std::vector<double> > const& fluidOutletWindkesselPressureDistalOld() const { return M_fluidOutletWindkesselPressureDistal_old; }
    trace_mesh_ptrtype const& fluidOutletWindkesselMesh() const { return M_fluidOutletWindkesselMesh; }
    space_fluidoutlet_windkessel_ptrtype const& fluidOutletWindkesselSpace() { return M_fluidOutletWindkesselSpace; }


    bool hasStrongDirichletBC() const
        {
            bool hasStrongDirichletBC = this->hasMarkerDirichletBCelimination() || this->hasFluidInlet() || M_bcMarkersMovingBoundaryImposed.hasMarkerDirichletBCelimination() || this->hasMarkerPressureBC();
            return hasStrongDirichletBC;
        }

    //___________________________________________________________________________________//

    void updateRangeDistributionByMaterialName( std::string const& key, range_faces_type const& rangeFaces );
    //___________________________________________________________________________________//

    std::shared_ptr<typename space_pressure_type::element_type>/*element_fluid_pressure_ptrtype*/ const& velocityDiv() const { return M_velocityDiv; }
    std::shared_ptr<typename space_pressure_type::element_type>/*element_fluid_pressure_ptrtype*/ velocityDiv() { return M_velocityDiv; }
    bool velocityDivIsEqualToZero() const { return M_velocityDivIsEqualToZero; }

    //___________________________________________________________________________________//

    // update normal stress
    //void updateNormalStressOnCurrentMesh();
    void updateNormalStressOnCurrentMesh( std::string const& nameOfRange, element_normalstress_ptrtype & fieldToUpdate );
    // update normal stress in reference ALE mesh
    void updateNormalStressOnReferenceMesh( std::string const& nameOfRange, element_normalstress_ptrtype & fieldToUpdate );

    void updateWallShearStress( std::string const& nameOfRange, element_normalstress_ptrtype & fieldToUpdate );
    void updateVorticity();

    template < typename ExprT >
    void updateVelocity(vf::Expr<ExprT> const& __expr)
    {
        M_fieldVelocity->on(_range=M_rangeMeshElements,_expr=__expr );
    }
    template < typename ExprT >
    void updatePressure(vf::Expr<ExprT> const& __expr)
    {
        M_fieldPressure->on(_range=M_rangeMeshElements,_expr=__expr );
    }

    template < typename ExprT >
    void updateSourceAdded(vf::Expr<ExprT> const& __expr)
    {
        if (!M_XhSourceAdded) this->createFunctionSpacesSourceAdded();
        M_SourceAdded->on(_range=elements( this->mesh()),_expr=__expr );
        M_haveSourceAdded=true;
    }
    template < typename ExprT >
    void updateVelocityDiv(vf::Expr<ExprT> const& __expr)
    {
        //if (!M_velocityDiv) M_velocityDiv=M_Xh->template functionSpace<1>()->elementPtr();
        if (!M_velocityDiv)
            M_velocityDiv.reset(new typename space_pressure_type::element_type/*element_fluid_pressure_type*/(M_XhPressure,"velocityDiv") );
        //*M_velocityDiv = vf::project(_space=M_Xh->template functionSpace<1>(),_range=elements( this->mesh()),_expr=__expr);
        M_velocityDiv->on(_range=elements(this->mesh()),_expr=__expr);
        M_velocityDivIsEqualToZero=false;
    }

#if defined( FEELPP_MODELS_HAS_MESHALE )
    template <typename element_mecasol_ptrtype>
    void updateStructureDisplacement(element_mecasol_ptrtype const & structSol);

    //void updateStructureDisplacementBis( typename mesh_ale_type::ale_map_element_ptrtype disp );
    void updateALEmesh();

    template <typename element_vel_mecasol_ptrtype>
    void updateStructureVelocity(element_vel_mecasol_ptrtype velstruct);
#endif
    //___________________________________________________________________________________//

    double computeMeshArea( std::string const& marker = "" ) const;
    double computeMeshArea( std::set<std::string> const& markers ) const;

    // compute measures : drag,lift,flow rate, mean pressure, mean div, norm div
    force_type computeForce( std::string const& markerName ) const;
    double computeFlowRate( std::string const& marker, bool useExteriorNormal=true ) const;
    double computeFlowRate( std::list<std::string> const& markers, bool useExteriorNormal=true ) const;
    double computePressureSum() const;
    double computePressureMean() const;
    double computeVelocityDivergenceSum() const;
    double computeVelocityDivergenceMean() const;
    double computeVelocityDivergenceNormL2() const;

#if 0
#if 0
    // Averaged Preassure computed on a set of slice (false for compute on actual mesh)
    template <typename SetMeshSlicesType>
    std::vector<double> computeAveragedPreassure( SetMeshSlicesType const & setMeshSlices,mpl::bool_<false> /**/);
    // Averaged Preassure computed on a set of slice (true for compute on ref mesh)
    template <typename SetMeshSlicesType>
    std::vector<double> computeAveragedPreassure( SetMeshSlicesType const & setMeshSlices,mpl::bool_<true> /**/);
    // Flow rate computed on a set of slice (false for compute on actual mesh)
    template <typename SetMeshSlicesType>
    std::vector<double> computeFlowRate(SetMeshSlicesType const & setMeshSlices,mpl::bool_<false> /**/);
    // Flow rate computed on a set of slice (true for compute on ref mesh)
    template <typename SetMeshSlicesType>
    std::vector<double> computeFlowRate(SetMeshSlicesType const & setMeshSlices,mpl::bool_<true> /**/);
#else
    typedef Mesh<Simplex<1,1,nRealDim> > mesh_slice1d_type;
    typedef std::shared_ptr<mesh_slice1d_type> mesh_slice1d_ptrtype;
    typedef typename mpl::at_c<typename space_fluid_velocity_type::bases_list,0>::type basis_slice_velocity_type;
    typedef FunctionSpace<mesh_slice1d_type, bases<basis_slice_velocity_type> > space_slice_velocity_type;
    typedef OperatorInterpolation<space_fluid_velocity_type,space_slice_velocity_type> op_interp_velocity_type;
    typedef std::shared_ptr<op_interp_velocity_type> op_interp_velocity_ptrtype;

    typedef typename mpl::at_c<typename space_fluid_pressure_type::bases_list,0>::type basis_slice_pressure_type;
    typedef FunctionSpace<mesh_slice1d_type, bases<basis_slice_pressure_type> > space_slice_pressure_type;
    typedef OperatorInterpolation<space_fluid_pressure_type,space_slice_pressure_type> op_interp_pressure_type;
    typedef std::shared_ptr<op_interp_pressure_type> op_interp_pressure_ptrtype;

#if defined( FEELPP_MODELS_HAS_MESHALE )
    typedef typename mpl::at_c<typename space_mesh_disp_type::bases_list,0>::type basis_slice_meshdisp_type;
    typedef FunctionSpace<mesh_slice1d_type, bases<basis_slice_meshdisp_type> > space_slice_meshdisp_type;
    typedef OperatorInterpolation<space_mesh_disp_type,space_slice_meshdisp_type> op_interp_meshdisp_type;
    typedef std::shared_ptr<op_interp_meshdisp_type> op_interp_meshdisp_ptrtype;

    std::vector<double> computeAveragedPreassure( std::vector<mesh_slice1d_ptrtype> const& setMeshSlices,
                                                  std::vector<op_interp_pressure_ptrtype> const& opInterp,
                                                  bool computeOnRefMesh=false,
                                                  std::vector<op_interp_meshdisp_ptrtype> const& opInterpMeshDisp = std::vector<op_interp_meshdisp_ptrtype>() );
    std::vector<double> computeFlowRate(std::vector<mesh_slice1d_ptrtype> const& setMeshSlices,
                                        std::vector<op_interp_velocity_ptrtype> const& opInterp,
                                        bool computeOnRefMesh=false,
                                        std::vector<op_interp_meshdisp_ptrtype> const& opInterpMeshDisp = std::vector<op_interp_meshdisp_ptrtype>() );
#else
    //TODO?
#endif

#endif
#endif
    //___________________________________________________________________________________//

    void solve();
    //___________________________________________________________________________________//
    void preSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const override;
    void postSolveNewton( vector_ptrtype rhs, vector_ptrtype sol ) const override;
    void preSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const override;
    void postSolvePicard( vector_ptrtype rhs, vector_ptrtype sol ) const override;
    void preSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const override;
    void postSolveLinear( vector_ptrtype rhs, vector_ptrtype sol ) const override;
    //___________________________________________________________________________________//

    void initInHousePreconditioner();
    void updateInHousePreconditioner( DataUpdateLinear & data ) const override;
    void updateInHousePreconditioner( DataUpdateJacobian & data ) const override;
    typedef OperatorPCDBase<typename space_velocity_type::value_type> operatorpcdbase_type;
    //typedef std::shared_ptr<operatorpcdbase_type> operatorpcdbase_ptrtype;
    void addUpdateInHousePreconditionerPCD( std::string const& name, std::function<void(operatorpcdbase_type &)> const& init,
                                            std::function<void(operatorpcdbase_type &,DataUpdateBase &)> const& up = std::function<void(operatorpcdbase_type &,DataUpdateBase &)>() )
        {
            M_addUpdateInHousePreconditionerPCD[name] = std::make_pair(init,up);
        }

    std::shared_ptr<operatorpcdbase_type> operatorPCD() const { return M_operatorPCD; }
    bool hasOperatorPCD() const { return ( M_operatorPCD.use_count() > 0 ); }
private :
    void updateInHousePreconditionerPMM( sparse_matrix_ptrtype const& mat, vector_ptrtype const& vecSol ) const;
    void updateInHousePreconditionerPCD( sparse_matrix_ptrtype const& mat, vector_ptrtype const& vecSol, DataUpdateBase & data ) const;
public :

    //___________________________________________________________________________________//

    // non linear (newton)
    void updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const override;
    void updateJacobian( DataUpdateJacobian & data ) const override;
    void updateResidual( DataUpdateResidual & data ) const override;

    void updateResidualStabilisation( DataUpdateResidual & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p ) const;
    void updateJacobianStabilisation( DataUpdateJacobian & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p ) const;
    template<typename DensityExprType, typename ViscosityExprType, typename... ExprT>
    void updateResidualStabilisationGLS( DataUpdateResidual & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p,
                                         Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                         std::string const& matName, const ExprT&... exprs ) const;
    template<typename DensityExprType, typename ViscosityExprType, typename... ExprT>
    void updateJacobianStabilisationGLS( DataUpdateJacobian & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p,
                                         Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu,
                                         std::string const& matName, const ExprT&... exprs ) const;
    void updateJacobianWeakBC( DataUpdateJacobian & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p ) const;
    void updateResidualWeakBC( DataUpdateResidual & data, element_velocity_external_storage_type const& u, element_pressure_external_storage_type const& p ) const;
    void updateJacobianDofElimination( DataUpdateJacobian & data ) const override;
    void updateResidualDofElimination( DataUpdateResidual & data ) const override;

    // linear
    void updateLinearPDE( DataUpdateLinear & data ) const override;
    void updateLinearPDEWeakBC( DataUpdateLinear & data ) const;
    void updateLinearPDEStabilisation( DataUpdateLinear & data ) const;
    template<typename DensityExprType, typename ViscosityExprType, typename AdditionalRhsType = hana::tuple<>, typename AdditionalMatType = hana::tuple<> >
    void updateLinearPDEStabilisationGLS( DataUpdateLinear & data, Expr<DensityExprType> const& rho, Expr<ViscosityExprType> const& mu, std::string const& matName,
                                          AdditionalRhsType const& addRhsTuple = hana::make_tuple(), AdditionalMatType const& addMatTuple = hana::make_tuple() ) const;
    void updateLinearPDEDofElimination( DataUpdateLinear & data ) const override;

    void updatePicard( DataUpdateLinear & data ) const override;
    double updatePicardConvergence( vector_ptrtype const& Unew, vector_ptrtype const& Uold ) const override;

    //___________________________________________________________________________________//

private :
    void updateBoundaryConditionsForUse();

    //protected:
    virtual size_type initStartBlockIndexFieldsInMatrix();
    virtual int initBlockVector();


    auto modelFieldsMeshALE( std::string const& prefix = "" ) const
        {
            using _field_disp_ptrtype = typename mesh_ale_type::ale_map_element_ptrtype;
            auto mfieldDisp = modelField<FieldCtx::ID,_field_disp_ptrtype>( FieldTag::mesh_displacement(this) );
            if ( this->isMoveDomain() )
                mfieldDisp.add( FieldTag::mesh_displacement(this), prefix, "displacement", this->meshALE()->displacement(), "disp", this->keyword() );
            return Feel::FeelModels::modelFields( mfieldDisp );
        }

    //----------------------------------------------------
    // mesh
    mesh_ptrtype M_mesh;
    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;
    MeshMover<mesh_type> M_mesh_mover;
    trace_mesh_ptrtype M_meshTrace;
    // fluid space and solution
    space_velocity_ptrtype M_XhVelocity;
    space_pressure_ptrtype M_XhPressure;
    element_velocity_ptrtype M_fieldVelocity;
    element_pressure_ptrtype M_fieldPressure;
    element_velocity_ptrtype M_fieldConvectionVelocityExtrapolated; // with Oseen solver
    // lagrange multiplier space for mean pressure
    std::vector<space_meanpressurelm_ptrtype> M_XhMeanPressureLM;
    // trace mesh and space
    trace_mesh_ptrtype M_meshDirichletLM;
    space_trace_velocity_ptrtype M_XhDirichletLM;
    // lagrange multiplier for impose pressure bc
    trace_mesh_ptrtype M_meshLagrangeMultiplierPressureBC;
    space_trace_velocity_component_ptrtype M_spaceLagrangeMultiplierPressureBC;
    element_trace_velocity_component_ptrtype M_fieldLagrangeMultiplierPressureBC1, M_fieldLagrangeMultiplierPressureBC2;
    // body bc
    BodySetBoundaryCondition M_bodySetBC;
    // time discrtisation fluid
    std::string M_timeStepping;
    bdf_velocity_ptrtype M_bdfVelocity;
    savets_pressure_ptrtype M_savetsPressure;
    double M_timeStepThetaValue;
    vector_ptrtype M_timeStepThetaSchemePreviousContrib;
    //----------------------------------------------------
    // normak boundary stress ans WSS
    space_normalstress_ptrtype M_XhNormalBoundaryStress;
    element_normalstress_ptrtype M_fieldNormalStress;
    element_normalstress_ptrtype M_fieldWallShearStress;
    // vorticity space
    space_vorticity_ptrtype M_XhVorticity;
    element_vorticity_ptrtype M_fieldVorticity;
    // fields defined in json
    std::map<std::string,component_element_velocity_ptrtype> M_fieldsUserScalar;
    std::map<std::string,element_velocity_ptrtype> M_fieldsUserVectorial;
    //----------------------------------------------------
    // mesh ale tool and space
    bool M_isMoveDomain;
#if defined( FEELPP_MODELS_HAS_MESHALE )
    mesh_ale_ptrtype M_meshALE;
    element_mesh_disp_ptrtype M_meshDisplacementOnInterface;
#endif
    //----------------------------------------------------
    // physical properties/parameters and space
    material_properties_ptrtype M_materialProperties; // TO REMOVE
    materialsproperties_ptrtype M_materialsProperties;
    // boundary conditions + body forces
    map_vector_field<nDim,1,2> M_bcDirichlet;
    std::map<ComponentType,map_scalar_field<2> > M_bcDirichletComponents;
    map_scalar_field<2> M_bcNeumannScalar, M_bcPressure;
    map_vector_field<nDim,1,2> M_bcNeumannVectorial;
    map_matrix_field<nDim,nDim,2> M_bcNeumannTensor2;
    map_vector_field<nDim,1,2> M_bcMovingBoundaryImposed;
    MarkerManagementDirichletBC M_bcMarkersMovingBoundaryImposed;
    map_vector_field<nDim,1,2> M_volumicForcesProperties;
    //---------------------------------------------------
    std::shared_ptr<RangeDistributionByMaterialName<mesh_type> > M_rangeDistributionByMaterialName;
    // range of mesh faces by material : (type -> ( matName -> ( faces range ) )
    //std::map<std::string,std::map<std::string,faces_reference_wrapper_t<mesh_type>>> M_rangeMeshFacesByMaterial;
    //---------------------------------------------------
    space_vectorial_PN_ptrtype M_XhSourceAdded;
    element_vectorial_PN_ptrtype M_SourceAdded;
    bool M_haveSourceAdded;
    //----------------------------------------------------
    std::shared_ptr<typename space_pressure_type::element_type>/*element_fluid_pressure_ptrtype*/ M_velocityDiv;
    bool M_velocityDivIsEqualToZero;
    //----------------------------------------------------
    //std::string M_modelName;
    std::string M_solverName;

    double M_dirichletBCnitscheGamma;

    bool M_useFSISemiImplicitScheme;
    std::string M_couplingFSIcondition;
    std::set<std::string> M_markersFSI;

    bool M_startBySolveNewtonian, M_hasSolveNewtonianAtKickOff;
    bool M_startBySolveStokesStationary, M_hasSolveStokesStationaryAtKickOff;
    bool M_applyMovingMeshBeforeSolve;
    //----------------------------------------------------
    // stabilization
    bool M_stabilizationGLS, M_stabilizationGLSDoAssembly;
    std::string M_stabilizationGLSType;
    stab_gls_parameter_ptrtype M_stabilizationGLSParameterConvectionDiffusion;
    stab_gls_parameter_ptrtype M_stabilizationGLSParameterPressure;
    std::map<std::string,range_elements_type> M_stabilizationGLSEltRangeConvectionDiffusion;
    std::map<std::string,range_elements_type> M_stabilizationGLSEltRangePressure;

    bool M_applyCIPStabOnlyOnBoundaryFaces;
    // stabilisation available
    bool M_doCIPStabConvection,M_doCIPStabDivergence,M_doCIPStabPressure;
    double M_stabCIPConvectionGamma,M_stabCIPDivergenceGamma,M_stabCIPPressureGamma;
    element_velocity_ptrtype M_fieldMeshVelocityUsedWithStabCIP;
    bool M_doStabDivDiv;
    bool M_doStabConvectionEnergy; // see Nobile thesis
    //----------------------------------------------------
    bool M_definePressureCst;
    bool M_definePressureCstOnlyOneZoneAppliedOnWholeMesh;
    std::vector<std::set<std::string> > M_definePressureCstMarkers;
    std::vector<range_elements_type> M_definePressureCstMeshRanges;
    std::string M_definePressureCstMethod;
    double M_definePressureCstPenalisationBeta;
    std::vector<std::pair<vector_ptrtype,std::set<size_type> > > M_definePressureCstAlgebraicOperatorMeanPressure;
    //----------------------------------------------------
    // fluid inlet bc
    std::vector< std::tuple<std::string,std::string, scalar_field_expression<2> > > M_fluidInletDesc; // (marker,type,vmax expr)
    std::map<std::string,trace_mesh_ptrtype> M_fluidInletMesh;
    std::map<std::string,space_fluidinlet_ptrtype> M_fluidInletSpace;
    std::map<std::string,element_fluidinlet_ptrtype > M_fluidInletVelocity;
    std::map<std::string,std::tuple<component_element_velocity_ptrtype,
                                    op_interpolation_fluidinlet_ptrtype > > M_fluidInletVelocityInterpolated;
    std::map<std::string,std::tuple<element_fluidinlet_ptrtype,double,double> > M_fluidInletVelocityRef;//marker->(uRef,maxURef,flowRateRef)
    //----------------------------------------------------
    // fluid outlet 0d (free, windkessel)
    std::vector< std::tuple<std::string,std::string, std::tuple<std::string,double,double,double> > > M_fluidOutletsBCType;
    mutable std::map<int,double> M_fluidOutletWindkesselPressureDistal,M_fluidOutletWindkesselPressureProximal;
    std::map<int,std::vector<double> > M_fluidOutletWindkesselPressureDistal_old;
    trace_mesh_ptrtype M_fluidOutletWindkesselMesh;
    space_fluidoutlet_windkessel_ptrtype M_fluidOutletWindkesselSpace;
    //----------------------------------------------------
    // post-process field exported
    // std::set<std::string> M_postProcessFieldExported;
    // std::set<std::string> M_postProcessFieldOnTraceExported;
    // std::set<std::string> M_postProcessUserFieldExported;

    // exporter option
    bool M_isHOVisu;
    // exporter fluid
    export_ptrtype M_exporter;
    export_trace_ptrtype M_exporterTrace;
    export_trace_ptrtype M_exporterFluidOutlet;
    export_trace_ptrtype M_exporterLagrangeMultiplierPressureBC;
    // exporter fluid ho
#if 1 //defined(FEELPP_HAS_VTK)
    export_ho_ptrtype M_exporter_ho;
    space_vectorial_visu_ho_ptrtype M_XhVectorialVisuHO;
    space_scalar_visu_ho_ptrtype M_XhScalarVisuHO;
    space_vectorialdisc_visu_ho_ptrtype M_XhVectorialDiscVisuHO;

    element_vectorial_visu_ho_ptrtype M_velocityVisuHO;
    element_scalar_visu_ho_ptrtype M_pressureVisuHO;
    element_vectorial_visu_ho_ptrtype M_meshdispVisuHO;
    //element_vectorialdisc_visu_ho_ptrtype M_normalStressVisuHO;
    //element_vectorialdisc_visu_ho_ptrtype M_fieldWallShearStressVisuHO;

    op_interpolation_visu_ho_vectorial_ptrtype M_opIvelocity;
    op_interpolation_visu_ho_scalar_ptrtype M_opIpressure;
#if defined( FEELPP_MODELS_HAS_MESHALE )
    op_interpolation_visu_ho_meshdisp_ptrtype M_opImeshdisp;
    MeshMover<mesh_visu_ho_type> M_meshmover_visu_ho;
#endif
    //op_interpolation_visu_ho_vectorialdisc_ptrtype M_opIstress;
#endif
    // post-process measure at point
    measure_points_evaluation_ptrtype M_measurePointsEvaluation;
    // post-process measure forces (lift,drag) and flow rate
    std::vector< ModelMeasuresForces > M_postProcessMeasuresForces;
    std::vector< ModelMeasuresFlowRate > M_postProcessMeasuresFlowRate;
    // post-process measure fields
    std::map<std::string,std::string> M_postProcessMeasuresFields;
    //----------------------------------------------------
    //----------------------------------------------------
    // algebraic data/tools
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;
    //----------------------------------------------------
    // overwrite assembly process : source terms
    typedef boost::function<void ( vector_ptrtype& F, bool buildCstPart )> updateSourceTermLinearPDE_function_type;
    updateSourceTermLinearPDE_function_type M_overwritemethod_updateSourceTermLinearPDE;
    typedef boost::function<void ( vector_ptrtype& R )> updateSourceTermResidual_function_type;
    updateSourceTermResidual_function_type M_overwritemethod_updateSourceTermResidual;
    //----------------------------------------------------
    bool M_preconditionerAttachPMM, M_preconditionerAttachPCD;
    mutable bool M_pmmNeedUpdate;
    std::shared_ptr<operatorpcdbase_type> M_operatorPCD;
    std::map<std::string,std::pair<std::function<void(operatorpcdbase_type &)>,std::function<void(operatorpcdbase_type &, DataUpdateBase &)> > > M_addUpdateInHousePreconditionerPCD;

}; // FluidMechanics


template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType>
template <typename SymbolsExpr>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType,BasisDVType>::exportResults( double time, SymbolsExpr const& symbolsExpr )
{
    this->log("FluidMechanics","exportResults", (boost::format("start at time %1%")%time).str() );
    this->timerTool("PostProcessing").start();

    // this->modelProperties().parameters().updateParameterValues();
    // auto paramValues = this->modelProperties().parameters().toParameterValues();
    // this->modelProperties().postProcess().setParameterValues( paramValues );

    this->updateFields( symbolsExpr );

    auto fields = this->modelFields();
    if constexpr ( nOrderGeo <= 2 )
    {
        this->executePostProcessExports( M_exporter, time, fields, symbolsExpr );
        this->executePostProcessExports( M_exporterTrace, "trace_mesh", time, fields, symbolsExpr );
    }
    if ( M_isHOVisu )
        this->exportResultsImplHO( time );

    this->executePostProcessMeasures( time, fields, symbolsExpr );
    this->executePostProcessSave( (this->isStationary())? invalid_uint32_type_value : M_bdfVelocity->iteration(), fields );

    if ( this->isMoveDomain() && this->hasPostProcessExportsField( "alemesh" ) )
        this->meshALE()->exportResults( time );

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("FluidMechanics","exportResults", "finish" );
}

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType, typename BasisDVType>
template <typename TupleFieldsType, typename SymbolsExpr>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType,BasisDVType>::executePostProcessMeasures( double time, TupleFieldsType const& tupleFields, SymbolsExpr const& symbolsExpr )
{
    bool hasMeasure = false;

    // forces (lift,drag) measures
    for ( auto const& ppForces : M_postProcessMeasuresForces )
    {
        CHECK( ppForces.meshMarkers().size() == 1 ) << "TODO";
        auto measuredForce = this->computeForce( ppForces.meshMarkers().front() );
        std::string name = ppForces.name();
        this->postProcessMeasuresIO().setMeasure( "drag_"+name, measuredForce(0,0) );
        this->postProcessMeasuresIO().setMeasure( "lift_"+name, measuredForce(1,0) );
        hasMeasure = true;
    }
    // flow rate measures
    for ( auto const& ppFlowRate : M_postProcessMeasuresFlowRate )
    {
        double valFlowRate = this->computeFlowRate( ppFlowRate.meshMarkers(), ppFlowRate.useExteriorNormal() );
        this->postProcessMeasuresIO().setMeasure("flowrate_"+ppFlowRate.name(),valFlowRate);
        hasMeasure = true;
    }

    if ( true )
    {
        bool hasMeasuresPressure = M_postProcessMeasuresFields.find("pressure") != M_postProcessMeasuresFields.end();
        bool hasMeasuresVelocityDivergence = M_postProcessMeasuresFields.find( "velocity-divergence" ) != M_postProcessMeasuresFields.end();
        double area = 0;
        if ( hasMeasuresPressure || hasMeasuresVelocityDivergence )
            area = this->computeMeshArea();
        if ( hasMeasuresPressure )
        {
            double pressureSum = this->computePressureSum();
            double pressureMean = pressureSum/area;
            this->postProcessMeasuresIO().setMeasure("pressure_sum",pressureSum);
            this->postProcessMeasuresIO().setMeasure("pressure_mean",pressureMean);
            hasMeasure = true;
        }
        if ( hasMeasuresVelocityDivergence )
        {
            double velocityDivergenceSum = this->computeVelocityDivergenceSum();
            double velocityDivergenceMean = velocityDivergenceSum/area;
            double velocityDivergenceNormL2 = this->computeVelocityDivergenceNormL2();
            this->postProcessMeasuresIO().setMeasure("velocity_divergence_sum",velocityDivergenceNormL2);
            this->postProcessMeasuresIO().setMeasure("velocity_divergence_mean",velocityDivergenceMean);
            this->postProcessMeasuresIO().setMeasure("velocity_divergence_normL2",velocityDivergenceNormL2);
            hasMeasure = true;
        }
    }


    bool hasMeasureNorm = this->updatePostProcessMeasuresNorm( this->mesh(), M_rangeMeshElements, symbolsExpr, tupleFields );
    bool hasMeasureStatistics = this->updatePostProcessMeasuresStatistics( this->mesh(), M_rangeMeshElements, symbolsExpr, tupleFields );
    bool hasMeasurePoint = this->updatePostProcessMeasuresPoint( M_measurePointsEvaluation, tupleFields );
    if ( hasMeasureNorm || hasMeasureStatistics || hasMeasurePoint )
        hasMeasure = true;

    if ( hasMeasure )
    {
        if ( !this->isStationary() )
            this->postProcessMeasuresIO().setMeasure( "time", time );
        this->postProcessMeasuresIO().exportMeasures();
        this->upload( this->postProcessMeasuresIO().pathFile() );
    }
}

//---------------------------------------------------------------------------------------------------------//
#if 0
template <typename SetMeshSlicesType>
std::vector<double>
FLUIDMECHANICS_CLASS_NAME::computeAveragedPreassure( SetMeshSlicesType const & setMeshSlices,mpl::bool_<false> /**/)
{
    using namespace Feel::vf;

    auto solFluid = this->getSolution();

    auto p = solFluid->element<1>();

    auto setMS = setMeshSlices.getContainerOfMeshSlices();

    auto nbSlice = setMS.size();
    std::vector<double> res(nbSlice);

    for (uint16_type i = 0 ; i<nbSlice ; ++i)
    {
        auto meshSlice = setMS[i];

        double area = integrate(_range=elements(meshSlice),
                                _expr=cst(1.) ).evaluate()(0,0);

        res[i] = (1./area)*(integrate(_range=elements(meshSlice),
                                      _expr=idv(p) ).evaluate()(0,0));
    }

    return res;
}

//---------------------------------------------------------------------------------------------------------//

template <typename SetMeshSlicesType>
std::vector<double>
FLUIDMECHANICS_CLASS_NAME::computeAveragedPreassure( SetMeshSlicesType const & setMeshSlices,mpl::bool_<true> /**/)
{
    using namespace Feel::vf;

    this->meshALE()->revertReferenceMesh();

    auto solFluid = this->getSolution();

    auto p = solFluid->element<1>();

    // Identity matrix
    auto Id = BOOST_PP_IF(BOOST_PP_EQUAL(FLUIDMECHANICS_DIM,2),
                          oneX()*trans(oneX()) + oneY()*trans(oneY()),
                          oneX()*trans(oneX()) + oneY()*trans(oneY())+oneZ()*trans(oneZ()) );
    // Deformation tensor
    auto Fa = Id+gradv(*M_meshALE->displacement());
#if (FLUIDMECHANICS_DIM==2)
    auto Fa11 = trans(Fa*oneX())*oneX();
    auto Fa12 = trans(Fa*oneY())*oneX();
    auto Fa21 = trans(Fa*oneX())*oneY();
    auto Fa22 = trans(Fa*oneY())*oneY();
    auto detFa = Fa11*Fa22-Fa21*Fa12;
    //sans le determinant devant car il s annule avec un terme apres
    //auto InvFa = mat<2,2>( Fa22,-Fa12,-Fa21,Fa11);
#endif
#if (FLUIDMECHANICS_DIM==3)
    auto Fa11 = trans(Fa*oneX())*oneX();
    auto Fa12 = trans(Fa*oneY())*oneX();
    auto Fa13 = trans(Fa*oneZ())*oneX();
    auto Fa21 = trans(Fa*oneX())*oneY();
    auto Fa22 = trans(Fa*oneY())*oneY();
    auto Fa23 = trans(Fa*oneZ())*oneY();
    auto Fa31 = trans(Fa*oneX())*oneZ();
    auto Fa32 = trans(Fa*oneY())*oneZ();
    auto Fa33 = trans(Fa*oneZ())*oneZ();
    auto detFa = Fa11*(Fa22*Fa33-Fa23*Fa32) - Fa21*(Fa12*Fa33-Fa13*Fa32) + Fa31*(Fa12*Fa23 - Fa13*Fa22);
    //sans le determinant devant car il s annule avec un terme apres
    //auto InvFa = mat<3,3>( Fa22*Fa33-Fa23*Fa32 , Fa13*Fa32-Fa12*Fa33 , Fa12*Fa23-Fa13*Fa22,
    //                       Fa23*Fa31-Fa21*Fa33 , Fa11*Fa33-Fa13*Fa31 , Fa13*Fa21-Fa11*Fa23,
    //                       Fa21*Fa32-Fa22*Fa31 , Fa12*Fa31-Fa11*Fa32 , Fa11*Fa22-Fa12*Fa21
    //                        );
#endif


    auto setMS = setMeshSlices.getContainerOfMeshSlices();

    auto nbSlice = setMS.size();
    std::vector<double> res(nbSlice);

    for (uint16_type i = 0 ; i<nbSlice ; ++i)
    {
        auto meshSlice = setMS[i];

        double area = integrate(_range=elements(meshSlice),
                                _expr=cst(1.) ).evaluate()(0,0);

        res[i] = (1./area)*(integrate(_range=elements(meshSlice),
                                      _expr=idv(p)*detFa ).evaluate()(0,0));
    }


    this->meshALE()->revertMovingMesh();

    return res;
}

//---------------------------------------------------------------------------------------------------------//

// Flow rate computed on a set of slice
template <typename SetMeshSlicesType>
std::vector<double>
FLUIDMECHANICS_CLASS_NAME::computeFlowRate(SetMeshSlicesType const & setMeshSlices,mpl::bool_<false> /**/)
{
    using namespace Feel::vf;

    auto solFluid = this->getSolution();
    auto u = solFluid->element<0>();

    auto setMS = setMeshSlices.getContainerOfMeshSlices();
    auto nbSlice = setMS.size();
    std::vector<double> res(nbSlice);

    auto dirVelocity = vec(cst(1.0),cst(0.));
    for (uint16_type i = 0 ; i<nbSlice ; ++i)
    {
        auto meshSlice = setMS[i];
        res[i] = integrate(_range=elements(meshSlice),
                           _expr=trans(idv(u))*dirVelocity ).evaluate()(0,0);
    }

    return res;

}

//---------------------------------------------------------------------------------------------------------//

// Flow rate computed on a set of slice
template <typename SetMeshSlicesType>
std::vector<double>
FLUIDMECHANICS_CLASS_NAME::computeFlowRate(SetMeshSlicesType const & setMeshSlices,mpl::bool_<true> /**/)
{
    using namespace Feel::vf;

    auto solFluid = this->getSolution();
    auto u = solFluid->element<0>();

    // Identity matrix
    auto Id = BOOST_PP_IF(BOOST_PP_EQUAL(FLUIDMECHANICS_DIM,2),
                          oneX()*trans(oneX()) + oneY()*trans(oneY()),
                          oneX()*trans(oneX()) + oneY()*trans(oneY())+oneZ()*trans(oneZ()) );
    // Deformation tensor
    auto Fa = Id+gradv(*M_meshALE->displacement());
#if (FLUIDMECHANICS_DIM==2)
    auto Fa11 = trans(Fa*oneX())*oneX();
    auto Fa12 = trans(Fa*oneY())*oneX();
    auto Fa21 = trans(Fa*oneX())*oneY();
    auto Fa22 = trans(Fa*oneY())*oneY();
    auto detFa = Fa11*Fa22-Fa21*Fa12;
    //sans le determinant devant car il s annule avec un terme apres
    //auto InvFa = mat<2,2>( Fa22,-Fa12,-Fa21,Fa11);
#endif
#if (FLUIDMECHANICS_DIM==3)
    auto Fa11 = trans(Fa*oneX())*oneX();
    auto Fa12 = trans(Fa*oneY())*oneX();
    auto Fa13 = trans(Fa*oneZ())*oneX();
    auto Fa21 = trans(Fa*oneX())*oneY();
    auto Fa22 = trans(Fa*oneY())*oneY();
    auto Fa23 = trans(Fa*oneZ())*oneY();
    auto Fa31 = trans(Fa*oneX())*oneZ();
    auto Fa32 = trans(Fa*oneY())*oneZ();
    auto Fa33 = trans(Fa*oneZ())*oneZ();
    auto detFa = Fa11*(Fa22*Fa33-Fa23*Fa32) - Fa21*(Fa12*Fa33-Fa13*Fa32) + Fa31*(Fa12*Fa23 - Fa13*Fa22);
    //sans le determinant devant car il s annule avec un terme apres
    //auto InvFa = mat<3,3>( Fa22*Fa33-Fa23*Fa32 , Fa13*Fa32-Fa12*Fa33 , Fa12*Fa23-Fa13*Fa22,
    //                       Fa23*Fa31-Fa21*Fa33 , Fa11*Fa33-Fa13*Fa31 , Fa13*Fa21-Fa11*Fa23,
    //                       Fa21*Fa32-Fa22*Fa31 , Fa12*Fa31-Fa11*Fa32 , Fa11*Fa22-Fa12*Fa21
    //                        );
#endif



    auto setMS = setMeshSlices.getContainerOfMeshSlices();
    auto nbSlice = setMS.size();
    std::vector<double> res(nbSlice);

    auto dirVelocity = vec(cst(1.0),cst(0.));
    for (uint16_type i = 0 ; i<nbSlice ; ++i)
    {
        auto meshSlice = setMS[i];
        res[i] = integrate(_range=elements(meshSlice),
                           _expr=trans(idv(u))*dirVelocity ).evaluate()(0,0);
    }

    return res;

}
#endif
//---------------------------------------------------------------------------------------------------------//

} // namespace FeelModels
} // namespace Feel

#include <feel/feelmodels/fluid/fluidmechanicsupdatestabilisationgls.hpp>

#endif /* FEELPP_TOOLBOXES_FLUIDMECHANICS_HPP */


