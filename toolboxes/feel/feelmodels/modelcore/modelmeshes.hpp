/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */
#pragma once

#include <feel/feelmesh/meshbase.hpp>
#include <feel/feelmesh/enums.hpp>
#include <feel/feeldiscr/geometricspace.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/nedelec.hpp>
#include <feel/feelfilters/databymeshentity.hpp>
#include <feel/feelmodels/modelcore/modelbase.hpp>
#include <feel/feelmodels/modelpostprocess.hpp>
#include <feel/feelmodels/modelmarkers.hpp>
#include <feel/feelmodels/modelcore/modelfields.hpp>

#include <feel/feelmodels/modelcore/modelmeasurespointsevaluation.hpp>
#include <feel/feelmodels/modelmesh/meshale.hpp>


namespace Feel
{
namespace FeelModels
{

template <typename IndexType>
class ModelMeshes;

/**
 * @brief Common class for Mesh Models
 * @ingroup ModelCore
 * 
 * @tparam IndexType 
 */
template <typename IndexType>
class ModelMeshCommon
{
public:
    using index_type = IndexType;
    using mesh_base_type = MeshBase<index_type>;
    using mesh_base_ptrtype = std::shared_ptr<mesh_base_type>;

    class ImportConfig
    {
    public :
        ImportConfig()
            :
            M_generatePartitioning( false ),
            M_numberOfPartition( 1 ),
            M_straightenMesh( true ),
            M_meshComponents( MESH_UPDATE_FACES|MESH_UPDATE_EDGES ),
            M_loadByMasterRankOnly( false )
            {}
        ImportConfig( ImportConfig const& ) = default;
        ImportConfig( ImportConfig && ) = default;
        ImportConfig( ModelMeshes<IndexType> const& mMeshes );

        std::string const& inputFilename() const { return M_inputFilename; }
        std::string const& meshFilename() const { return M_meshFilename; }
        std::string const& geoFilename() const { return M_geoFilename; }
        bool generatePartitioning() const { return M_generatePartitioning; }
        int numberOfPartition() const { return M_numberOfPartition; }
        double meshSize() const { return M_meshSize; }
        bool straightenMesh() const { return M_straightenMesh; }
        size_type meshComponents() const { return M_meshComponents; }
        bool loadByMasterRankOnly() const { return M_loadByMasterRankOnly; }

        void setStraightenMesh( bool b ) { M_straightenMesh = b; }
        void setMeshComponents( size_type c ) { M_meshComponents = c; }

        bool hasMeshFilename() const { return !M_meshFilename.empty(); }
        bool hasGeoFilename() const { return !M_geoFilename.empty(); }

        void setup( nl::json const& jarg, ModelMeshes<IndexType> const& mMeshes );
        void updateForUse( ModelMeshes<IndexType> const& mMeshes );

        void setupInputMeshFilenameWithoutApplyPartitioning( std::string const& filename );
        void setupSequentialAndLoadByMasterRankOnly();

        void updateInformationObject( nl::json & p ) const;
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    private :
        std::string M_inputFilename, M_meshFilename, M_geoFilename;
        bool M_generatePartitioning;
        int M_numberOfPartition;
        double M_meshSize;
        bool M_straightenMesh;
        size_type M_meshComponents;
        bool M_loadByMasterRankOnly;
    };


    ModelMeshCommon() = default;
    ModelMeshCommon( ModelMeshes<IndexType> const& mMeshes ) : M_importConfig( mMeshes ) {}
    ModelMeshCommon( ModelMeshCommon const& ) = default;
    ModelMeshCommon( ModelMeshCommon && ) = default;

    ImportConfig & importConfig() { return M_importConfig; }
    ImportConfig const& importConfig() const { return M_importConfig; }

    bool hasMesh() const { return M_mesh? true : false; }

    template <typename MeshType = mesh_base_type>
    auto mesh() const
        {
            if constexpr( std::is_same_v<MeshType,mesh_base_type> )
                return M_mesh;
            else
                return std::dynamic_pointer_cast<MeshType>( M_mesh );
        }

    std::string const& meshFilename() const { return M_meshFilename; }

    void setMesh( mesh_base_ptrtype m, std::string const& filename = "" ) { M_mesh = m; M_meshFilename = filename;  }

    template <typename MeshType>
    void initMeasurePointsEvaluationTool()
        {
            if ( !M_measurePointsEvaluation )
            {
                auto m = this->mesh<MeshType>();
                auto geospace = std::make_shared<GeometricSpace<MeshType>>( m );
                M_measurePointsEvaluation = std::make_shared<MeasurePointsEvaluation<MeshType>>( geospace );
            }
        }

    template <typename MeshType>
    auto measurePointsEvaluationTool()
        {
            return std::dynamic_pointer_cast<MeasurePointsEvaluation<MeshType>>( M_measurePointsEvaluation );
        }

    bool hasMeshMotion() const { return M_meshMotionTool? true : false; }

    template <typename MeshType>
    auto meshMotionTool() const
        {
            return std::dynamic_pointer_cast<MeshALE<typename MeshType::shape_type>>( M_meshMotionTool );
        }

    void setMeshMotionTool( std::shared_ptr<ModelBase> meshMotionTool ) { M_meshMotionTool = meshMotionTool; }

    template <typename SpaceType>
    auto createFunctionSpace( std::string const& basis )
        {
            auto itFind = M_functionSpaces.find( basis );
            if ( itFind == M_functionSpaces.end() )
            {
                auto Vh = SpaceType::New(_mesh=this->mesh<typename SpaceType::mesh_type>());
                M_functionSpaces[basis] = Vh;
                return Vh;
            }
            else
                return std::dynamic_pointer_cast<SpaceType>( itFind->second );
        }

    template <typename SpaceType = FunctionSpaceBase>
    auto functionSpace( std::string const& basis ) const
        {
            auto itFind = M_functionSpaces.find( basis );
            CHECK( itFind!=  M_functionSpaces.end() ) << "no space with basis " << basis;
            if constexpr( std::is_same_v<SpaceType,FunctionSpaceBase> )
                return itFind->second;
            else
                return std::dynamic_pointer_cast<SpaceType>( itFind->second );
        }

    void clearFunctionSpaces() { M_functionSpaces.clear(); }

private:
    ImportConfig M_importConfig;
    mesh_base_ptrtype M_mesh;
    std::string M_meshFilename;
    std::map<std::string, std::shared_ptr<FunctionSpaceBase> > M_functionSpaces;
    std::shared_ptr<MeasurePointsEvaluationBase> M_measurePointsEvaluation;
    std::shared_ptr<ModelBase> M_meshMotionTool;
};


/**
 * @brief Mesh Model class
 * @ingroup ModelCore
 * 
 * @tparam IndexType
 */
template <typename IndexType>
class ModelMesh
{
public :
    using index_type = IndexType;
    using mesh_base_type = MeshBase<index_type>;
    using mesh_base_ptrtype = std::shared_ptr<mesh_base_type>;
    using collection_data_by_mesh_entity_type = CollectionOfDataByMeshEntity<index_type>;
    using import_config_type = typename ModelMeshCommon<IndexType>::ImportConfig;
private :
    /**
     * @brief A struct that represents the setup of fields.
     *
     * This struct contains information about the name, basis, and parts of a fields.
     */
    struct FieldsSetup
    {
          /**
           * @brief A struct that represents the setup of a part.
           *
           * This struct contains information about the filename, model expression, and markers of a part.
           */
          struct PartSetup
          {
              /**
               * @brief Default constructor.
               */
              PartSetup() = default;
  
              /**
               * @brief Returns the filename of the part.
               *
               * @return A reference to the filename string.
               */
              std::string const& filename() const { return M_filename; }
  
              /**
               * @brief Returns whether the part has a model expression.
               *
               * @return True if the part has a model expression, false otherwise.
               */
              bool hasExpr() const { return M_mexpr ? true : false; }
  
              /**
               * @brief Returns the model expression of the part.
               *
               * @return A reference to the model expression object.
               */
              ModelExpression const& mExpr() const { return *M_mexpr; }
  
              /**
               * @brief Returns the markers of the part.
               *
               * @return A reference to the model markers object.
               */
              ModelMarkers const& markers() const { return M_markers; }
  
              /**
               * @brief Creates a new PartSetup object from the given model meshes and JSON data.
               *
               * @param mMeshes The model meshes object.
               * @param jarg The JSON data.
               * @return An optional PartSetup object if successful, std::nullopt otherwise.
               */
              static std::optional<PartSetup> create( ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg );
  
              private:
              std::string M_filename;                 ///< The filename of the part.
              std::optional<ModelExpression> M_mexpr; ///< The model expression of the part.
              ModelMarkers M_markers;                 ///< The markers of the part.
          };
  
          /**
           * @brief Constructs a new FieldsSetup object from the given name, model meshes, and JSON data.
           *
           * @param name The name of the fields.
           * @param mMeshes The model meshes object.
           * @param jarg The JSON data.
           */
          FieldsSetup( std::string const& name, ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg );
  
          /**
           * @brief Returns the name of the fields.
           *
           * @return A reference to the name string.
           */
          std::string const& name() const { return M_name; }
  
          /**
           * @brief Returns the basis of the fields.
           *
           * @return A reference to the basis string.
           */
          std::string const& basis() const { return M_basis; }
  
          /**
           * @brief Returns the parts of the fields.
           *
           * @return A reference to the vector of part setup objects.
           */
          std::vector<PartSetup> const& parts() const { return M_parts; }
  
          private:
          std::string M_name;             ///< The name of the fields.
          std::string M_basis;            ///< The basis of the fields.
          std::vector<PartSetup> M_parts; ///< The parts of the fields.
    };

    struct DistanceToRangeSetup
    {
        struct Normalization {
            Normalization( std::string const& name, std::string const& type ) : M_name( name ), M_type( type ) {}
            Normalization( std::string const& name, std::string const& type, std::pair<ModelExpression,ModelExpression> const& range ) : M_name( name ), M_type( type ), M_range( range ) {}
            std::string const& name() const { return M_name; }
            std::string const& type() const { return M_type; }
            std::pair<ModelExpression,ModelExpression> const& range() const { return M_range.value(); }
            static std::optional<Normalization> create( ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg );
        private :
            std::string M_name;
            std::string M_type;
            std::optional<std::pair<ModelExpression,ModelExpression>> M_range;
        };

        DistanceToRangeSetup( std::string const& name, ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg );
        std::string const& name() const { return M_name; }
        std::set<std::string> const& markers() const { return M_markers; }
        std::vector<Normalization> const& normalizations() const { return M_normalizations; }
        double maxDistance() const { return M_maxDistance.template hasExpr<1,1>()? M_maxDistance.template expr<1,1>().evaluate()(0,0) : -1.; }

    private :
        std::string M_name;
        std::set<std::string> M_markers;
        std::vector<Normalization> M_normalizations;
        ModelExpression M_maxDistance;
    };
    struct MeshMotionSetup
    {
        MeshMotionSetup( ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg );

        std::set<std::string> const& computationalDomainMarkers() { return M_computationalDomainMarkers; }
        std::map<std::string,std::tuple<ModelExpression,std::set<std::string>>> const& displacementImposed() const { return M_displacementImposed; }
        std::set<std::string> const& displacementZeroMarkers() const { return M_displacementZeroMarkers; }
        std::set<std::string> const& displacementFreeMarkers() const { return M_displacementFreeMarkers; }

        void setParameterValues( std::map<std::string,double> const& paramValues )
            {
                for ( auto & [name,dispData] : M_displacementImposed )
                    std::get<0>( dispData ).setParameterValues( paramValues );
            }

    private :
        std::set<std::string> M_computationalDomainMarkers;
        std::map<std::string,std::tuple<ModelExpression,std::set<std::string>>> M_displacementImposed;
        std::set<std::string> M_displacementZeroMarkers, M_displacementFreeMarkers;
    };

public :

    struct MeshAdaptation
    {

        struct Event {
            enum class Type { never=0, after_import, after_init, each_time_step/*, mesh_quality*/ };
            virtual ~Event() {};
            virtual Type type() const { return Type::never; }
        };
        struct EventAfterImport : public Event
        {
            typename Event::Type type() const override { return Event::Type::after_import; }
        };
        struct EventAfterInit : public Event
        {
            typename Event::Type type() const override { return Event::Type::after_init; }
        };
        struct EventEachTimeStep : public Event
        {
            EventEachTimeStep( double t, size_type index ) : M_currentTime( t ), M_currentIndex( index ) {}
            typename Event::Type type() const override { return Event::Type::each_time_step; }
            double currentTime() const { return M_currentTime; }
            size_type currentIndex() const { return M_currentIndex; }
        private :
            double M_currentTime;
            size_type M_currentIndex;
        };

        template <typename Event::Type TT>
        static auto createEvent() { return std::make_shared<Event>(); }
        template <>
        static auto createEvent<Event::Type::after_import>() { return std::make_shared<EventAfterImport>(); }
        template <>
        static auto createEvent<Event::Type::after_init>() { return std::make_shared<EventAfterInit>(); }
        template <typename Event::Type TT>
        static auto createEvent( double t, size_type index ) { return std::make_shared<Event>(); }
        template <>
        static auto createEvent<Event::Type::each_time_step>( double t, size_type index ) { return std::make_shared<EventEachTimeStep>( t, index ); }

    private :

        template <typename TT>
        friend class ModelMesh;

        struct Execute;

        struct Setup
        {
            Setup( ModelMeshes<IndexType> const& mMeshes, nl::json const& jarg );

            ModelExpression const& metricExpr() const { return M_metric; }

            std::set<std::string> const& requiredMarkers() const { return M_requiredMarkers; }

            template <typename SymbolsExprType = symbols_expression_empty_t>
            bool isExecutedWhen( std::shared_ptr<Event> e, SymbolsExprType const& se = symbols_expression_empty_t{} ) const
                {
                    if ( M_executionEvents.find( e->type() ) == M_executionEvents.end() )
                        return false;

                    switch ( e->type() )
                    {
                    case Event::Type::each_time_step:
                    {
                        auto e_ets = std::dynamic_pointer_cast<EventEachTimeStep>( e );
                        if ( e_ets->currentIndex() == M_eventEachTimeStep_lastExecutionIndex )
                            return false;
                        if ( ( e_ets->currentIndex() % M_eventEachTimeStep_frequency ) != 0 )
                            return false;
                        if ( M_eventEachTimeStep_condition.template hasExpr<1,1>() )
                        {
                            auto conditionExpr = expr( M_eventEachTimeStep_condition.template expr<1,1>(), se );
                            if ( std::abs(conditionExpr.evaluate()(0,0)) < 1e-12 )
                                return false;
                        }
                        break;
                    }
                    default:
                        break;
                    }
                    return true;
                }

            void setParameterValues( std::map<std::string,double> const& paramValues )
                {
                    M_metric.setParameterValues( paramValues );
                    M_eventEachTimeStep_condition.setParameterValues( paramValues );
                }

            void updateInformationObject( nl::json & jsonInfo ) const;

            static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

            friend struct Execute;
        private:

            nl::json M_remesherSetup;
            std::set<std::string> M_requiredMarkers;
            std::set<typename Event::Type> M_executionEvents;
            size_type M_eventEachTimeStep_frequency;
            size_type M_eventEachTimeStep_lastExecutionIndex;
            ModelExpression M_eventEachTimeStep_condition;
            ModelExpression M_metric;
            fs::path M_tmpDir; // store tmp mesh file in parallel, should be removed when ParMmg work!
        };

        struct Execute
        {
            template <typename MeshType,typename SymbolsExprType>
            static
            std::shared_ptr<MeshType> execute( typename MeshAdaptation::Setup const& mas, std::shared_ptr<MeshType> inputMesh, std::shared_ptr<Event> event, SymbolsExprType const& se )
                {
                    if ( !mas.isExecutedWhen( event,se ) )
                        return {};
                    Execute mae{mas};
                    return mae.execute(inputMesh,se);
                }

            Execute( typename MeshAdaptation::Setup const& mas ) : M_mas( mas ) {}

            template <typename MeshType,typename SymbolsExprType>
            std::shared_ptr<MeshType> execute( std::shared_ptr<MeshType> inputMesh, SymbolsExprType const& se )
                {
                    auto Vh = Pch<1>( inputMesh );
                    //auto met = Vh->element();
                    //M_scalarMetricField = Vh->elementPtr();
                    auto metField = Vh->elementPtr();
                    M_scalarMetricField = metField;
                    auto metExpr = expr( M_mas.metricExpr().template expr<1,1>(), se );
                    metField->on( _range=elements(inputMesh),_expr=metExpr );

                    if  constexpr ( (MeshType::nDim == MeshType::nRealDim) && MeshType::nOrder == 1 )
                                      return this->executeImpl( inputMesh );
                    else
                        return std::shared_ptr<MeshType>{};
                }

        private:
            template <typename MeshType>
            std::shared_ptr<MeshType> executeImpl( std::shared_ptr<MeshType> inputMesh );

        private:
            typename MeshAdaptation::Setup const& M_mas;
            std::shared_ptr<Vector<double>> M_scalarMetricField;
        };
    };

    //public :

    ModelMesh( std::string const& name )
        :
        M_name( name ),
        M_mmeshCommon( std::make_shared<ModelMeshCommon<IndexType>>() )
        {}
    ModelMesh( std::string const& name, ModelMeshes<IndexType> const& mMeshes );
    ModelMesh( ModelMesh const& ) = default;
    ModelMesh( ModelMesh && ) = default;

    //! return json metadata
    nl::json const& metadata() const { return M_metadata; }

    void setup( nl::json const& jarg, ModelMeshes<IndexType> const& mMeshes );

    void setupRestart( ModelMeshes<IndexType> const& mMeshes );

    void setMesh( mesh_base_ptrtype m, std::string const& meshFilename = "" ) { M_mmeshCommon->setMesh( m, meshFilename ); }

    void setAsShared( ModelMesh<IndexType> const& m )
        {
            M_mmeshCommon = m.M_mmeshCommon;
        }

    import_config_type & importConfig() { return M_mmeshCommon->importConfig(); }
    import_config_type const& importConfig() const { return M_mmeshCommon->importConfig(); }

    template <typename MeshType>
    void updateForUse( ModelMeshes<IndexType> const& mMeshes );

    template <typename MeshType>
    void initMeasurePointsEvaluationTool()
        {
            M_mmeshCommon->template initMeasurePointsEvaluationTool<MeshType>();
        }
    template <typename MeshType>
    auto measurePointsEvaluationTool()
        {
            return M_mmeshCommon->template measurePointsEvaluationTool<MeshType>();
        }

    template <typename MeshType = mesh_base_type>
    auto mesh() const { return M_mmeshCommon->template mesh<MeshType>(); }

    bool hasMeshMotion() const { return M_mmeshCommon->hasMeshMotion(); }

    template <typename MeshType>
    auto meshMotionTool() const
        {
            return M_mmeshCommon->template meshMotionTool<MeshType>();
        }


    std::map<std::string,collection_data_by_mesh_entity_type> const& collectionOfDataByMeshEntity() const { return M_codbme; }

    /*!
     * \brief Returns a tuple of tuples containing the names and types of the basis fields supported by a mesh.
     *
     * The function template takes a template parameter `MeshType` which is a type that represents a mesh. 
     * The function returns a tuple of tuples, where each inner tuple contains a string and a type. 
     * The string represents the name of a basis field type supported by the mesh, and the type represents the type of the basis field.
     *
     * The function template is a static member function of the class that it is defined in. 
     * It is a constexpr function, which means that it can be evaluated at compile-time if the template parameter `MeshType` is a constant expression.
     *
     * The function template uses a C++17 feature called structured bindings to create the tuple of tuples. 
     * It uses the `hana` library to create the tuples and the `type_c` function to create the type objects.
     *
     * The function template is specialized for different types of meshes, depending on their dimension and the number of real dimensions. 
     * For example, for a 1D mesh with 1 real dimension, the function returns a tuple of tuples containing the names and types of the Lagrange basis fields with 1, 2, 0, and 1 DoFs. 
     * For a 2D mesh with 2 real dimensions, the function returns a tuple of tuples containing the names and types of the Lagrange basis fields (continuous order 1,2 and discontinuous 0,1 and vectorial continuous order 1,2 as well as Nedelec basis field 1st family of order 1)
     *
     * \tparam MeshType The type of the mesh.
     * \return A tuple of tuples containing the names and types of the basis fields supported by the mesh.
     */
    template <typename MeshType>
    static constexpr auto basisFieldTypeSupported()
    {
        if constexpr ( MeshType::nDim == 1 )
        {
            if constexpr( MeshType::nRealDim == 1 )
                return hana::make_tuple( hana::make_tuple( "Pch1", hana::type_c<Lagrange<1,Scalar,Continuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pch2", hana::type_c<Lagrange<2,Scalar,Continuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pdh0", hana::type_c<Lagrange<0,Scalar,Discontinuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pdh1", hana::type_c<Lagrange<1,Scalar,Discontinuous,PointSetFekete>> )
                                        );
            else
                return hana::make_tuple( hana::make_tuple( "Pch1", hana::type_c<Lagrange<1,Scalar,Continuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pch2", hana::type_c<Lagrange<2,Scalar,Continuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pdh0", hana::type_c<Lagrange<0,Scalar,Discontinuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pdh1", hana::type_c<Lagrange<1,Scalar,Discontinuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pchv1", hana::type_c<Lagrange<1,Vectorial,Continuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pchv2", hana::type_c<Lagrange<2,Vectorial,Continuous,PointSetFekete>> )
                                        );
        }
        else
            return hana::make_tuple( hana::make_tuple( "Pch1", hana::type_c<Lagrange<1,Scalar,Continuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pch2", hana::type_c<Lagrange<2,Scalar,Continuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pdh0", hana::type_c<Lagrange<0,Scalar,Discontinuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pdh1", hana::type_c<Lagrange<1,Scalar,Discontinuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pchv1", hana::type_c<Lagrange<1,Vectorial,Continuous,PointSetFekete>> ),
                                        hana::make_tuple( "Pchv2", hana::type_c<Lagrange<2,Vectorial,Continuous,PointSetFekete>> ),
                                        hana::make_tuple( "Ned1h0", hana::type_c<Nedelec<0,NedelecKind::NED1>> )
                                        );
    }
    /**
     * @brief Returns a tuple of model fields for a given mesh.
     *
     * @tparam MeshType The type of mesh for which to return the model fields.
     *
     * @param prefix_field A string prefix to use for the field names.
     * @param prefix_symbol A string prefix to use for the symbol names.
     *
     * @return A tuple of model fields for the given mesh.
     */
    template <typename MeshType>
    auto modelFields( std::string const& prefix_field = "", std::string const& prefix_symbol = "" ) const
    {
        static constexpr auto tuple_t_basis = ModelMesh<IndexType>::template basisFieldTypeSupported<MeshType>();

        auto mf_fields = hana::fold( tuple_t_basis, model_fields_empty_t{}, [this,prefix_field,prefix_symbol]( auto const& r, auto const& cur )
                                    {
                                        using basis_type = typename std::decay_t<decltype(hana::at_c<1>( cur ) )>::type;
                                        //using basis_type = typename std::decay_t<decltype(cur)>::type;
                                        using space_type = FunctionSpace<MeshType, bases<basis_type> >;
                                        using field_type = typename space_type::element_type;
                                        using field_ptrtype = std::shared_ptr<field_type>;

                                        auto mftag = ModelFieldTag< ModelMesh<index_type>,0>( this );
                                        auto mf = modelField<FieldCtx::FULL,field_ptrtype>( mftag );

                                        for ( auto const& [name,fieldBase] : M_fields )
                                        {
                                            if ( auto u = std::dynamic_pointer_cast<field_type>( fieldBase ) )
                                                mf.add( mftag, prefixvm( prefix_field,"fields" ), name, u, name, prefixvm( prefix_symbol, "fields", "_" ) );
                                        }
                                        return Feel::FeelModels::modelFields( r,mf );
                                    });

        // distance to range
        using distange_to_range_space_type = FunctionSpace<MeshType, bases<Lagrange<MeshType::nOrder,Scalar,Continuous,PointSetFekete>>>;
        using distange_to_range_field_type = typename distange_to_range_space_type::element_type;
        using distange_to_range_field_ptrtype = std::shared_ptr<distange_to_range_field_type>;
        auto mftag_distToRange = ModelFieldTag< ModelMesh<index_type>,1>( this );
        auto mf_distToRange = modelField<FieldCtx::ID,distange_to_range_field_ptrtype>( mftag_distToRange );
        for ( auto const& [name,fieldBase] : M_distanceToRanges )
        {
            if ( auto u = std::dynamic_pointer_cast<distange_to_range_field_type>( fieldBase ) )
                mf_distToRange.add( mftag_distToRange, prefixvm( prefix_field, "distanceToRange" ), name, u, name, prefixvm( prefix_symbol, "distanceToRange", "_" ) );
        }

        using mf_meshmotion_type = std::decay_t<decltype(this->meshMotionTool<MeshType>()->modelFields(""))>;
        mf_meshmotion_type mf_meshmotion;
        if ( auto mmt = this->meshMotionTool<MeshType>() )
            mf_meshmotion = mmt->modelFields( prefixvm( prefix_field,"meshmotion") );

        return Feel::FeelModels::modelFields( std::move(mf_fields), std::move(mf_distToRange), std::move( mf_meshmotion ) );
    }

    template <typename MeshType = mesh_base_type, bool AddFields = true>
    auto symbolsExpr( std::string const& prefix_symbol = "" ) const
    {
        using _se_type = std::decay_t< decltype( M_codbme.begin()->second.symbolsExpr() )>;
        _se_type res;
        for ( auto const& [name,data] : M_codbme )
        {
            res =  Feel::vf::symbolsExpr( res, data.symbolsExpr( prefixvm( prefix_symbol, "data_"+name, "_" ) ) );
        }

        if constexpr( AddFields && !std::is_same_v<MeshType,mesh_base_type> )
        {
            return Feel::vf::symbolsExpr( res, this->modelFields<MeshType>( "",prefix_symbol ).symbolsExpr() );
        }
        else
        {
            return res;
        }
    }

    void updateTime( double time )
    {
        for ( auto & [name,data] : M_codbme )
            data.updateTime( time );
    }


    void updateInformationObject( nl::json & p, std::string const& prefix_symbol ) const;

    static tabulate_informations_ptr_t tabulateInformations( nl::json const& p, TabulateInformationProperties const& tabInfoProp );


    void setParameterValues( std::map<std::string,double> const& paramValues )
    {
        if ( M_meshMotionSetup )
            M_meshMotionSetup->setParameterValues( paramValues );
        for ( typename MeshAdaptation::Setup & mas : M_meshAdaptationSetup )
                mas.setParameterValues( paramValues );
    }


    template <typename MeshType>
    void
    updateDistanceToRange();


    template <typename MeshType,typename SymbolsExprType>
    void
    updateMeshMotion( SymbolsExprType const& se )
        {
            auto meshALE = this->meshMotionTool<MeshType>();
            if ( !meshALE )
                return;
            if ( !M_meshMotionSetup )
                return;

            if ( !M_meshMotionSetup->displacementImposed().empty() )
            {
                bool meshIsOnRefAtBegin = meshALE->isOnReferenceMesh();
                if ( !meshIsOnRefAtBegin )
                    meshALE->revertReferenceMesh( false );
                meshALE->revertInitialDomain( false );

                auto mesh = this->mesh<MeshType>();
                for ( auto const& [name,dispData] : M_meshMotionSetup->displacementImposed() )
                {
                    auto const& [mexpr,markers] = dispData;
                    meshALE->updateDisplacementImposedOnInitialDomain( M_name/*this->keyword()*/,
                                                                       Feel::vf::expr( mexpr.template expr<MeshType::nRealDim,1>(), se ),
                                                                       markedfaces(mesh,markers) );
                }

                meshALE->revertReferenceMesh( false );
                if ( !meshIsOnRefAtBegin )
                    meshALE->revertMovingMesh( false );
            }

            meshALE->updateMovingMesh();
        }


    template <typename MeshType,typename SymbolsExprType>
    void
    updateMeshAdaptation( std::shared_ptr<typename MeshAdaptation::Event> event, SymbolsExprType const& se )
        {
             auto oldmesh = this->mesh<MeshType>();

             for ( typename MeshAdaptation::Setup & mas : M_meshAdaptationSetup )
             {
                 auto out = MeshAdaptation::Execute::execute( mas, oldmesh, event, se );
                 if ( !out )
                     continue;

                 if ( M_functionApplyRemesh )
                     std::invoke( M_functionApplyRemesh, oldmesh, out );
                 else
                     this->applyRemesh( out );

                 break;
             }
        }


    using function_apply_remesh = std::function<void ( mesh_base_ptrtype,mesh_base_ptrtype )>;
    void setFunctionApplyRemesh( function_apply_remesh f ) { M_functionApplyRemesh = f; }

    template <typename MeshType>
    void applyRemesh( std::shared_ptr<MeshType> const& newMesh );

    /**
     * @brief Updates a field in a model mesh based on an expression and a range of elements.
     * This function template updates a field in a model mesh based on an expression and a range of elements. The function takes several parameters:
     * @param name A string representing the name of the field to be updated.
     * @param expr An expression that represents the value to be assigned to the field.
     * @param range A range of elements over which the field should be updated.
     * @param basis A string representing the basis function space to be used for the field. If no basis is provided, the function will default to a scalar or vectorial function space based on the expression type.
     * @param initZero A boolean flag indicating whether the field should be initialized to zero before being updated.
     * 
     * The function template uses Hana library to iterate over a tuple of supported basis function spaces for the given mesh type. 
     * For each basis function space, it checks if it matches the provided basis parameter. 
     *  - If a match is found, it creates a function space of the appropriate type using the provided basis and checks if the field already exists in the mesh. 
     *  - If the field does not exist, it is created using the function space. 
     *  - If the field already exists but has a different basis, the field is rebuilt. 
     *  - Finally, the function applies the expression to the field over the specified range of elements.
     */
    template <typename MeshType, typename ExprTypr, typename RangeType, std::enable_if_t< is_filter_v<std::decay_t<RangeType>>, bool> = true>
    void updateField( std::string const& name, Expr<ExprTypr> const& expr, RangeType const& range, std::string const& basis = "", bool initZero = false )
    {
        using expr_shape_type = typename ExprTraits<RangeType,Expr<ExprTypr>>::shape;
        std::string basisUsed = basis;
        if ( basis.empty() )
        {
            if ( expr_shape_type::is_scalar )
                basisUsed = "Pch1";
            else if ( expr_shape_type::is_vectorial )
                    basisUsed = "Pchv1";
            // TODO tensor case
        }

        auto themesh = this->mesh<MeshType>();
        if ( !themesh )
            return;
        static constexpr auto tuple_t_basis = ModelMesh<IndexType>::template basisFieldTypeSupported<MeshType>();
        hana::for_each( tuple_t_basis, [this,&name,&expr,&basisUsed,&range,&initZero,&themesh]( auto const& b ) {
                                            if ( basisUsed != hana::at_c<0>( b ) )
                                                return;
                                            using basis_type = typename std::decay_t<decltype(hana::at_c<1>( b ) )>::type;
                                            using space_type = FunctionSpace<MeshType, bases<basis_type> >;
                                            if constexpr ( expr_shape_type::M == space_type::nComponents1 && expr_shape_type::N == space_type::nComponents2 )
                                            {
                                                bool doInit = false;
                                                if ( M_fields.find( name ) == M_fields.end() )
                                                    doInit = true;
                                                else if ( !dynamic_cast<typename space_type::element_type*>( M_fields.at( name ).get() ) )
                                                    doInit = true;
                                                // init field if necessary
                                                if ( doInit )
                                                {
                                                    auto Vh = M_mmeshCommon->template createFunctionSpace<space_type>( basisUsed );
                                                    auto u = Vh->elementPtr();
                                                    M_fields[name] = u;
                                                }
                                                // apply nodal interpolation
                                                auto u = dynamic_cast<typename space_type::element_type*>( M_fields.at( name ).get() );
                                                if ( initZero )
                                                    u->zero();
                                                u->on(_range=range,_expr=expr);
                                            }
                                            else
                                                throw std::runtime_error( fmt::format("ModelMesh<IndexType>::updateField expression and space incompatible : expr [{},{}] vs space [{},{}]",
                                                                                        expr_shape_type::M, expr_shape_type::N, space_type::nComponents1, space_type::nComponents2 ) );
                                        } );
    }

    /**
     * @brief Updates a field with a given name and expression on all mesh elements.
     *
     * This function template is specialized for a given MeshType and ExprType and is used to update a field with a given name on all mesh elements.
     * @param name The name of the field to update.
     * @param expr The expression used to update the field.
     * @param basis The basis used to update the field.
     * @tparam MeshType The type of mesh to update the field on.
     * @tparam ExprType The type of expression used to update the field.
     */
    template <typename MeshType, typename ExprTypr>
    void updateField( std::string const& name, Expr<ExprTypr> const& expr, std::string const& basis = "" )
    {
        auto themesh = this->mesh<MeshType>();
        if ( !themesh )
            return;
        this->updateField<MeshType>( name, expr, elements(themesh), basis );
    }

    /**
     * @brief Updates a field with data from another field.
     *
     * This function template is specialized for a function space element FieldType and a mesh type MeshType
     * The purpose of this function is to update a field with data from another field. The function takes in the name of the field to be updated, the field to be used for the update, and an optional basis string.
     *
     * @param name The name of the field to be updated.
     * @param u The field to be used for the update.
     * @param basis An optional basis string.
     *
     * \todo get a default basis if empty (iterate over all space, maybe one will match or use also continuity property)
     */
    template <typename MeshType, typename FieldType, std::enable_if_t<is_functionspace_element_v<std::decay_t<FieldType>>, bool> = true>
    void updateField( std::string const& name, FieldType const& u, std::string const& basis = "" )
    {
        this->updateField<MeshType>( name, idv(u), basis );
    }

    /**
     * @brief Removes a field given its name from the map.
     *
     * This function removes the field with the given name from the M_fields map and returns true if the field was successfully removed, or false if the field was not found in the map.
     *
     * @param name The name of the field to remove.
     * @return true if the field was successfully removed, or false if the field was not found in the map.
     */
    bool removeField( std::string const& name ) { return M_fields.erase( name ) > 0; }

private:
    std::string M_name;
    nl::json M_metadata;

    std::shared_ptr<ModelMeshCommon<IndexType>> M_mmeshCommon;

    std::map<std::string,collection_data_by_mesh_entity_type> M_codbme;

    std::vector<FieldsSetup> M_fieldsSetup;
    std::vector<DistanceToRangeSetup> M_distanceToRangeSetup;
    std::map<std::string, std::shared_ptr<Vector<double>> > M_fields;
    std::map<std::string, std::shared_ptr<Vector<double>> > M_distanceToRanges;
    std::optional<MeshMotionSetup> M_meshMotionSetup;
    std::vector<typename MeshAdaptation::Setup> M_meshAdaptationSetup;
    function_apply_remesh M_functionApplyRemesh;

};

/**
 * @brief Mesh Model Collection
 * @ingroup ModelCore
 * 
 * @tparam IndexType 
 */
template <typename IndexType>
class ModelMeshes : protected std::map<std::string,std::shared_ptr<ModelMesh<IndexType>>>,
                    virtual public ModelBase
{
public:
    using index_type = IndexType;
    using modelmesh_type = ModelMesh<IndexType>;
    using mesh_base_type = typename modelmesh_type::mesh_base_type;
    using mesh_base_ptrtype = typename modelmesh_type::mesh_base_ptrtype;
    using mesh_adaptation_type = typename modelmesh_type::MeshAdaptation;

    ModelMeshes() : ModelBase("")
    {
        auto me = std::make_shared<ModelMesh<IndexType>>( this->keyword(), *this );
        this->emplace( std::make_pair( this->keyword(), std::move( me ) ) );
    }
    ModelMeshes( ModelMeshes const& ) = default;
    ModelMeshes( ModelMeshes && ) = default;

    bool hasModelMesh( std::string const& meshName ) const
    {
        return this->find( meshName ) != this->end();
    }
    ModelMesh<IndexType> & modelMesh( std::string const& meshName )
    {
        auto itFindMesh = this->find( meshName );
        CHECK( itFindMesh != this->end() ) << "mesh not found : " << meshName;
        return *itFindMesh->second;
    }
    ModelMesh<IndexType> & modelMesh() { return this->modelMesh( this->keyword() ); }

    ModelMesh<IndexType> const& modelMesh( std::string const& meshName ) const
    {
        auto itFindMesh = this->find( meshName );
        CHECK( itFindMesh != this->end() ) << "mesh not found : " << meshName;
        return *itFindMesh->second;
    }
    ModelMesh<IndexType> const& modelMesh() const { return this->modelMesh( this->keyword() ); }

    void setup( nl::json const& jarg, std::set<std::string> const& keywordsToSetup );
    void setup( nl::json const& jarg ) { this->setup( jarg, {this->keyword()} ); }

    void setupRestart( std::string const& meshName )
    {
        if ( !this->hasModelMesh( meshName ) )
            this->emplace( std::make_pair( meshName, std::make_shared<ModelMesh<IndexType>>( meshName ) ) );
        this->modelMesh( meshName ).setupRestart( *this );
    }

    void setMesh( std::string const& meshName, mesh_base_ptrtype m )
    {
        if ( !this->hasModelMesh( meshName ) )
            this->emplace( std::make_pair( meshName, std::make_shared<ModelMesh<IndexType>>( meshName ) ) );
        this->modelMesh( meshName ).setMesh( m );
    }

    //! shared data (mesh,space,...) between several ModelMesh
    void setModelMeshAsShared( std::string const& meshName, ModelMesh<IndexType> const& m )
    {
        if ( !this->hasModelMesh( meshName ) )
            this->emplace( std::make_pair( meshName, std::make_shared<ModelMesh<IndexType>>( meshName ) ) );
        this->modelMesh( meshName ).setAsShared( m );
    }
    void setModelMeshAsShared( ModelMesh<IndexType> const& m ) { this->setModelMeshAsShared( this->keyword(), m ); }

    template <typename MeshType>
    void updateForUse( std::string const& meshName )
    {
        if ( this->hasModelMesh( meshName ) )
            this->modelMesh( meshName ).template updateForUse<MeshType>( *this );
    }

    template <typename MeshType>
    void initMeasurePointsEvaluationTool( std::string const& meshName )
        {
            if ( this->hasModelMesh( meshName ) )
                this->modelMesh( meshName ).template initMeasurePointsEvaluationTool<MeshType>();
        }
    template <typename MeshType>
    auto measurePointsEvaluationTool( std::string const& meshName )
        {
            return this->modelMesh( meshName ).template measurePointsEvaluationTool<MeshType>();
        }


    template <typename MeshType = mesh_base_type>
        auto mesh( std::string const& meshName ) const
    {
        return this->modelMesh( meshName ).template mesh<MeshType>();
    }

    template <typename MeshType>
    auto modelFields( std::string const& prefix_field, std::string const& prefix_symbol = "" ) const
        {
            using _mf_type = std::decay_t< decltype( this->begin()->second->template modelFields<MeshType>() )>;
            _mf_type res;
            for ( auto const& [meshName,mMesh] : *this )
                res = Feel::FeelModels::modelFields( res, mMesh->template modelFields<MeshType>( prefixvm( prefix_field, meshName ), prefixvm( prefix_symbol, meshName, "_" )  ) );
            return res;
        }

    template <typename MeshType = mesh_base_type, bool AddFields = true>
    auto symbolsExpr( std::string const& prefix_symbol = "meshes" ) const
    {
        using _se_type = std::decay_t< decltype( this->begin()->second->template symbolsExpr<MeshType,AddFields>() )>;
        _se_type res;
        for ( auto const& [meshName,mMesh] : *this )
            res =  Feel::vf::symbolsExpr( res, mMesh->template symbolsExpr<MeshType,AddFields>( prefixvm( prefix_symbol, meshName, "_" ) ) );
        return res;
    }

    void updateTime( double time )
    {
        for ( auto & [meshName,mMesh] : *this )
            mMesh->updateTime( time );
    }

    void updateInformationObject( nl::json & p ) const override { this->updateInformationObject( p, "meshes" ); }
    void updateInformationObject( nl::json & p, std::string const& prefix_symbol ) const;

    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

    void setParameterValues( std::map<std::string,double> const& paramValues )
        {
            for ( auto & [name,mmesh] : *this )
                mmesh->setParameterValues( paramValues );
        }

    bool hasMeshMotion( std::string const& meshName ) const
        {
            if ( !hasModelMesh( meshName ) )
                return false;
            return this->modelMesh( meshName ).hasMeshMotion();
        }

    template <typename MeshType>
    auto meshMotionTool( std::string const& meshName ) const
        {
            return this->modelMesh( meshName ).template meshMotionTool<MeshType>();
        }

    template <typename MeshType,typename SymbolsExprType>
    void
    updateMeshMotion( std::string const& meshName, SymbolsExprType const& se )
    {
        if ( this->hasModelMesh( meshName ) )
            this->modelMesh( meshName ).template updateMeshMotion<MeshType>( se );
    }

    template <typename MeshType,typename SymbolsExprType>
    void
    updateMeshAdaptation( std::string const& meshName, std::shared_ptr<typename mesh_adaptation_type::Event> event, SymbolsExprType const& se )
        {
            if ( this->hasModelMesh( meshName ) )
                this->modelMesh( meshName ).template updateMeshAdaptation<MeshType>( event, se );
        }


    template <typename MeshType>
    void applyRemesh( std::string const& meshName, std::shared_ptr<MeshType> const& newMesh )
        {
            this->modelMesh( meshName ).applyRemesh( newMesh );
        }

    std::string repository_meshes() const { return (fs::path(this->repository().root())/fmt::format("{}.meshes",this->keyword())).string(); }

    void saveMetadata() const
        {
            nl::json metadata;
            for ( auto const& [meshName,mMesh] : *this )
            {
                metadata[meshName] = mMesh->metadata();
            }


            if ( this->worldComm().isMasterRank() )
            {
                if ( !fs::exists( this->repository_meshes()) )
                    fs::create_directories( this->repository_meshes() );

                //M_tmpDir = (fs::path(mMeshes.repository().root())/"meshes")/"tmp";
                std::ofstream ofile( (fs::path(this->repository_meshes())/"metadata.json").string(), std::ios::trunc);
                if ( ofile )
                {
                    ofile << metadata.dump(1);
                    ofile.close();
                }
            }

        }
};


} // namespace FeelModels
} // namespace Feel
