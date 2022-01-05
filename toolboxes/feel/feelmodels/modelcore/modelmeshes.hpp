/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */
#pragma once

#include <feel/feelmesh/meshbase.hpp>
#include <feel/feelmesh/enums.hpp>
#include <feel/feeldiscr/geometricspace.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/databymeshentity.hpp>
#include <feel/feelmodels/modelcore/modelbase.hpp>
#include <feel/feelmodels/modelpostprocess.hpp>
#include <feel/feelmodels/modelmarkers.hpp>
#include <feel/feelmodels/modelcore/modelfields.hpp>

#include <feel/feelmodels/modelcore/modelmeasurespointsevaluation.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename IndexType>
class ModelMeshes;

template <typename IndexType>
class ModelMeshCommon
{
public:
    using index_type = IndexType;
    using mesh_base_type = MeshBase<index_type>;
    using mesh_base_ptrtype = std::shared_ptr<mesh_base_type>;

    ModelMeshCommon() = default;
    ModelMeshCommon( ModelMeshCommon const& ) = default;
    ModelMeshCommon( ModelMeshCommon && ) = default;

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

private:
    mesh_base_ptrtype M_mesh;
    std::string M_meshFilename;
    std::map<std::string, std::shared_ptr<FunctionSpaceBase> > M_functionSpaces;
    std::shared_ptr<MeasurePointsEvaluationBase> M_measurePointsEvaluation;
};

template <typename IndexType>
class ModelMesh
{
public :
    using index_type = IndexType;
    using mesh_base_type = MeshBase<index_type>;
    using mesh_base_ptrtype = std::shared_ptr<mesh_base_type>;
    using collection_data_by_mesh_entity_type = CollectionOfDataByMeshEntity<index_type>;

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

    ModelMesh( std::string const& name )
        :
        M_name( name ),
        M_mmeshCommon( std::make_shared<ModelMeshCommon<IndexType>>() )
        {}
    ModelMesh( std::string const& name, ModelMeshes<IndexType> const& mMeshes );
    ModelMesh( ModelMesh const& ) = default;
    ModelMesh( ModelMesh && ) = default;

    void setup( nl::json const& jarg, ModelMeshes<IndexType> const& mMeshes );

    void setupRestart( ModelMeshes<IndexType> const& mMeshes );

    void setMesh( mesh_base_ptrtype m, std::string const& meshFilename = "" ) { M_mmeshCommon->setMesh( m, meshFilename ); }

    void setAsShared( ModelMesh<IndexType> const& m )
        {
            M_mmeshCommon = m.M_mmeshCommon;
        }

    ImportConfig & importConfig() { return M_importConfig; }
    ImportConfig const& importConfig() const { return M_importConfig; }

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

    std::map<std::string,collection_data_by_mesh_entity_type> const& collectionOfDataByMeshEntity() const { return M_codbme; }

    template <typename MeshType>
    auto modelFields( std::string const& prefix_symbol = "" ) const
        {
            static constexpr auto tuple_t_basis = hana::to_tuple(hana::tuple_t<
                                                                 Lagrange<1,Scalar,Continuous,PointSetFekete>,
                                                                 Lagrange<2,Scalar,Continuous,PointSetFekete>,
                                                                 Lagrange<1,Vectorial,Continuous,PointSetFekete>,
                                                                 Lagrange<2,Vectorial,Continuous,PointSetFekete>
                                                                 >);

            std::string prefix = prefix_symbol;

            auto mf_fields = hana::fold( tuple_t_basis, model_fields_empty_t{}, [this,prefix,prefix_symbol]( auto const& r, auto const& cur )
                                       {
                                           using basis_type = typename std::decay_t<decltype(cur)>::type;
                                           using space_type = FunctionSpace<MeshType, bases<basis_type> >;
                                           using field_type = typename space_type::element_type;
                                           using field_ptrtype = std::shared_ptr<field_type>;

                                           auto mftag = ModelFieldTag< ModelMesh<index_type>,0>( this );
                                           auto mf = modelField<FieldCtx::ID,field_ptrtype>( mftag );

                                           for ( auto const& [name,fieldBase] : M_fields )
                                           {
                                               if ( auto u = std::dynamic_pointer_cast<field_type>( fieldBase ) )
                                                   mf.add( mftag, prefix, name, u, name, prefixvm( prefix_symbol, "fields", "_" ) );
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
                    mf_distToRange.add( mftag_distToRange, prefix, name, u, name, prefixvm( prefix_symbol, "distanceToRange", "_" ) );
            }

            return Feel::FeelModels::modelFields( std::move(mf_fields), std::move(mf_distToRange) );
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
                return Feel::vf::symbolsExpr( res, this->modelFields<MeshType>( prefix_symbol ).symbolsExpr() );
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


    void updateInformationObject( nl::json & p ) const;

    static tabulate_informations_ptr_t tabulateInformations( nl::json const& p, TabulateInformationProperties const& tabInfoProp );



private:
    std::string M_name;
    std::shared_ptr<ModelMeshCommon<IndexType>> M_mmeshCommon;
    ImportConfig M_importConfig;
    std::map<std::string,collection_data_by_mesh_entity_type> M_codbme;


    struct FieldsSetup
    {
        FieldsSetup( std::string const& name, nl::json const& jarg )
            :
            M_name( name )
            {
                if ( jarg.contains("filename") )
                    M_filename = Environment::expand( jarg.at("filename") );
                else
                    CHECK( false ) << "filename required";

                if ( jarg.contains("basis") )
                    M_basis = Environment::expand( jarg.at("basis") );
                else
                    CHECK( false ) << "basis required";
            }
        std::string const& name() const { return M_name; }
        std::string const& filename() const { return M_filename; }
        std::string const& basis() const { return M_basis; }
    private :
        std::string M_name;
        std::string M_filename;
        std::string M_basis;
    };

    struct DistanceToRangeSetup
    {
        DistanceToRangeSetup( std::string const& name, nl::json const& jarg )
            :
            M_name( name )
            {
                auto itFind = jarg.find("markers");
                if ( itFind == jarg.end() )
                    itFind = jarg.find("marker");
                if ( itFind != jarg.end() )
                {
                    ModelMarkers mm;
                    mm.setup( *itFind );
                    M_markers = mm;
                }
            }
        std::string const& name() const { return M_name; }
        std::set<std::string> const& markers() const { return M_markers; }
    private :
        std::string M_name;
        std::set<std::string> M_markers;
    };

    std::vector<FieldsSetup> M_fieldsSetup;
    std::vector<DistanceToRangeSetup> M_distanceToRangeSetup;
    std::map<std::string, std::shared_ptr<Vector<double>> > M_fields;
    std::map<std::string, std::shared_ptr<Vector<double>> > M_distanceToRanges;
};

template <typename IndexType>
class ModelMeshes : protected std::map<std::string,std::shared_ptr<ModelMesh<IndexType>>>,
                    virtual public ModelBase
{
    using index_type = IndexType;
    using mesh_base_type = MeshBase<IndexType>;
    using mesh_base_ptrtype = std::shared_ptr<mesh_base_type>;
public:
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

    void setup( pt::ptree const& pt );

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

    void updateInformationObject( nl::json & p ) const override;

    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

};


} // namespace FeelModels
} // namespace Feel
