/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */
#pragma once

#include <feel/feelmesh/meshbase.hpp>
#include <feel/feelfilters/databymeshentity.hpp>

namespace Feel
{
namespace FeelModels
{

template <typename IndexType>
class ModelMesh
{
public :
    using index_type = IndexType;
    using mesh_base_type = MeshBase<index_type>;
    using mesh_base_ptrtype = std::shared_ptr<mesh_base_type>;
    using collection_data_by_mesh_entity_type = CollectionOfDataByMeshEntity<index_type>;

    ModelMesh() = default;
    ModelMesh( ModelMesh const& ) = default;
    ModelMesh( ModelMesh && ) = default;

    void setup( pt::ptree const& pt );

    void setMesh( mesh_base_ptrtype m ) { M_mesh = m; }

    template <typename MeshType>
    void updateForUse()
        {
            for ( auto & [name,data] : M_codbme )
            {
                data.setMesh( M_mesh );
                data.template updateForUse<MeshType>();
            }
        }

    template <typename MeshType = mesh_base_type>
    auto mesh() const
        {
            if constexpr( std::is_same_v<MeshType,mesh_base_type> )
                            return M_mesh;
            else
            {
                auto m = std::dynamic_pointer_cast<MeshType>( M_mesh );
                CHECK( m ) << "dynamic_pointer_cast fail";
                return m;
            }
        }

    std::map<std::string,collection_data_by_mesh_entity_type> const& collectionOfDataByMeshEntity() const { return M_codbme; }

    auto symbolsExpr( std::string const& prefix_symbol = "" ) const
        {
            using _se_type = std::decay_t< decltype( M_codbme.begin()->second.symbolsExpr() )>;
            _se_type res;
            for ( auto const& [name,data] : M_codbme )
            {
                res =  Feel::vf::symbolsExpr( res, data.symbolsExpr( prefixvm( prefix_symbol, "data_"+name, "_" ) ) );
            }
            return res;
        }

    void updateTime( double time )
        {
            for ( auto & [name,data] : M_codbme )
                data.updateTime( time );
        }


private:
    mesh_base_ptrtype M_mesh;
    std::map<std::string,collection_data_by_mesh_entity_type> M_codbme;
};

template <typename IndexType>
class ModelMeshes : protected std::map<std::string,std::shared_ptr<ModelMesh<IndexType>>>,
                    virtual public ModelBase
{
    using index_type = IndexType;
    using mesh_base_type = MeshBase<IndexType>;
    using mesh_base_ptrtype = std::shared_ptr<mesh_base_type>;
public:
    ModelMeshes() : ModelBase("") {}
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
    ModelMesh<IndexType> const& modelMesh( std::string const& meshName ) const
    {
        auto itFindMesh = this->find( meshName );
        CHECK( itFindMesh != this->end() ) << "mesh not found : " << meshName;
        return *itFindMesh->second;
    }

    void setup( pt::ptree const& pt );

    void setMesh( std::string const& meshName, mesh_base_ptrtype m )
    {
        if ( !this->hasModelMesh( meshName ) )
            this->emplace( std::make_pair( meshName, std::make_shared<ModelMesh<IndexType>>() ) );
        this->modelMesh( meshName ).setMesh( m );
    }

    template <typename MeshType>
    void updateForUse( std::string const& meshName )
    {
        this->modelMesh( meshName ).template updateForUse<MeshType>();
    }

    template <typename MeshType = mesh_base_type>
        auto mesh( std::string const& meshName ) const
    {
        return this->modelMesh( meshName ).template mesh<MeshType>();
    }

    auto symbolsExpr( std::string const& prefix_symbol = "meshes" ) const
    {
        using _se_type = std::decay_t< decltype( this->begin()->second->symbolsExpr() )>;
        _se_type res;
        for ( auto const& [meshName,mMesh] : *this )
            res =  Feel::vf::symbolsExpr( res, mMesh->symbolsExpr( prefixvm( prefix_symbol, meshName, "_" ) ) );
        return res;
    }

    void updateTime( double time )
    {
        for ( auto & [meshName,mMesh] : *this )
            mMesh->updateTime( time );
    }
};

template <typename IndexType>
void
ModelMeshes<IndexType>::setup( pt::ptree const& pt )
{
    for ( auto const& item : pt )
    {
        std::string meshName = item.first;
        auto me = std::make_shared<ModelMesh<IndexType>>();
        me->setup( item.second );
        this->emplace( std::make_pair( meshName, std::move( me ) ) );
    }
}

template <typename IndexType>
void
ModelMesh<IndexType>::setup( pt::ptree const& pt )
{
    if ( auto importPtree = pt.get_child_optional("Import") )
    {
        // TODO
    }

    if ( auto dataPtree = pt.get_child_optional("Data") )
    {
        for ( auto const& item : *dataPtree )
        {
            std::string dataName = item.first;
            collection_data_by_mesh_entity_type c;
            c.setup( item.second );
            M_codbme.emplace( std::make_pair( dataName, std::move( c ) ) );
        }
    }

}


} // namespace FeelModels
} // namespace Feel
