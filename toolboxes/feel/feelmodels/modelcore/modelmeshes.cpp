/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/modelmeshes.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmodels/modelmesh/createmesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>

namespace Feel
{
namespace FeelModels
{


template <typename IndexType>
ModelMesh<IndexType>::ImportConfig::ImportConfig( ModelMeshes<IndexType> const& mMeshes )
    :
    M_generatePartitioning( boption(_prefix=mMeshes.prefix(),_name="gmsh.partition",_vm=mMeshes.clovm()) ),
    M_numberOfPartition( mMeshes.worldComm().localSize() ),
    M_meshSize( doption(_prefix=mMeshes.prefix(),_name="gmsh.hsize",_vm=mMeshes.clovm()) ),
    M_meshComponents( MESH_UPDATE_FACES|MESH_UPDATE_EDGES )
{
    if ( mMeshes.clovm().count( prefixvm(mMeshes.prefix(),"mesh.filename").c_str() ) )
        M_inputFilename = Environment::expand( soption(_prefix=mMeshes.prefix(),_name="mesh.filename",_vm=mMeshes.clovm()) );
}

template <typename IndexType>
void
ModelMesh<IndexType>::ImportConfig::setup( pt::ptree const& pt, ModelMeshes<IndexType> const& mMeshes )
{
    if ( auto filenameOpt = pt.template get_optional<std::string>( "filename" ) )
        M_inputFilename = Environment::expand( *filenameOpt );

    if ( auto generatePartitioningOpt = pt.template get_optional<bool>( "partition" ) )
        M_generatePartitioning = *generatePartitioningOpt;

    if ( auto nPartitionOpt = pt.template get_optional<int>( "number-of-partition" ) )
    {
        CHECK( *nPartitionOpt > 0 && *nPartitionOpt <= mMeshes.worldComm().localSize() ) << "invalid number of partition : " << *nPartitionOpt;
        M_numberOfPartition = *nPartitionOpt;
    }

    if ( auto hsizeOpt = pt.template get_optional<double>( "hsize" ) )
        M_meshSize = *hsizeOpt;
}

template <typename IndexType>
void
ModelMesh<IndexType>::ImportConfig::setupInputMeshFilenameWithoutApplyPartitioning( std::string const& filename )
{
    // TODO
}

template <typename IndexType>
void
ModelMesh<IndexType>::ImportConfig::updateForUse( ModelMeshes<IndexType> const& mMeshes )
{
    if ( M_inputFilename.empty() )
        return;
    std::string meshfile = M_inputFilename;
    RemoteData rdTool( M_inputFilename, mMeshes.worldCommPtr() );
    if ( rdTool.canDownload() )
    {
        auto dowloadedData = rdTool.download( (fs::path(Environment::downloadsRepository())/fs::path(mMeshes.prefix())/fs::path("meshes")).string() );
        CHECK( dowloadedData.size() > 0 ) << "no data download";
        meshfile = dowloadedData[0];
        if ( dowloadedData.size() == 2 )
        {
            if ( fs::path( dowloadedData[0] ).extension() == ".h5" && fs::path( dowloadedData[1] ).extension() == ".json" )
                meshfile = dowloadedData[1];
        }
    }

    if ( fs::path( meshfile ).extension() == ".geo" )
        M_geoFilename = meshfile;
    else
        M_meshFilename = meshfile;
}

template <typename IndexType>
ModelMesh<IndexType>::ModelMesh( std::string const& name, ModelMeshes<IndexType> const& mMeshes )
    :
    M_name( name ),
    M_importConfig( mMeshes )
{}


template <typename IndexType>
void
ModelMesh<IndexType>::setupRestart( ModelMeshes<IndexType> const& mMeshes )
{
    std::string meshFilename;
    if ( mMeshes.worldComm().isMasterRank() )
    {
        fs::path thedir = mMeshes.rootRepository();
        std::string fileNameMeshPath = (thedir/prefixvm(M_name,"mesh.path")).string();
        std::ifstream file( fileNameMeshPath.c_str() );
        if ( !file )
            CHECK( false ) << "Fail to open the txt file containing path of mesh file : " << fileNameMeshPath << "\n";

        std::string meshFilename;
        if ( ! ( file >> meshFilename ) )
            CHECK( false ) << "Fail to read the msh path in file : " << fileNameMeshPath << "\n";
        file.close();
    }
    mpi::broadcast( mMeshes.worldComm().localComm(), meshFilename, mMeshes.worldComm().masterRank() );

    M_importConfig.setupInputMeshFilenameWithoutApplyPartitioning( meshFilename );
}


template <typename IndexType>
template <typename MeshType>
void
ModelMesh<IndexType>::updateForUse( ModelMeshes<IndexType> const& mMeshes )
{
    using mesh_type = MeshType;

    if ( !M_mesh )
    {
        M_importConfig.updateForUse( mMeshes );

        if ( M_importConfig.hasMeshFilename() )
        {
            std::string const& inputMeshFilename = M_importConfig.meshFilename();
            mMeshes.log("createMeshModel","", "load mesh file : " + inputMeshFilename);

            std::string rootpath = mMeshes.rootRepository();
            std::string meshPartitionedFilename = (fs::path( rootpath ) / (mMeshes.prefix() + ".json")).string();
            std::string meshFileExt = fs::path( inputMeshFilename ).extension().string();
            bool generatePartitioning = M_importConfig.generatePartitioning();
            if ( generatePartitioning && meshFileExt != ".msh" )
                CHECK( false ) << "Can not rebuild at this time the mesh partitionining with other format than .msh : TODO";

            M_mesh = loadMesh(_mesh=new mesh_type( M_name, mMeshes.worldCommPtr() ),
                              _filename=inputMeshFilename,
                              _prefix=mMeshes.prefix(),
                              _vm=mMeshes.clovm(),
                              _worldcomm=mMeshes.worldCommPtr(),
                              _rebuild_partitions=generatePartitioning,
                              _rebuild_partitions_filename=meshPartitionedFilename,
                              _partitions=M_importConfig.numberOfPartition(),
                              _savehdf5=0,
                              _update= M_importConfig.meshComponents()/*MESH_UPDATE_EDGES|MESH_UPDATE_FACES*/);

            M_meshFilename = (generatePartitioning)? meshPartitionedFilename : M_importConfig.meshFilename();
        }
        else if ( M_importConfig.hasGeoFilename() )
        {
            std::string const& inputGeoFilename = M_importConfig.geoFilename();
            std::string path = mMeshes.rootRepository();

            std::string mshfile = (fs::path( path ) / mMeshes.prefix()).string();
            if ( M_importConfig.numberOfPartition() > 1 )
                mshfile += ".json";
            else
                mshfile += ".msh";

            gmsh_ptrtype geodesc = geo( _filename=inputGeoFilename,
                                        _prefix=mMeshes.prefix(),
                                        _vm=mMeshes.clovm(),
                                        _worldcomm=mMeshes.worldCommPtr() );
            // allow to have a geo and msh file with a filename equal to prefix
            geodesc->setPrefix(mMeshes.prefix());
            M_mesh = createGMSHMesh(_mesh=new mesh_type( M_name, mMeshes.worldCommPtr() ),
                                    _desc=geodesc,
                                    _prefix=mMeshes.prefix(),
                                    _vm=mMeshes.clovm(),
                                    _worldcomm=mMeshes.worldCommPtr(),
                                    _partitions=M_importConfig.numberOfPartition(),
                                    _update=M_importConfig.meshComponents(),
                                    _directory=mMeshes.rootRepository() );
            M_meshFilename = mshfile;
        }

        if ( M_mesh && !M_meshFilename.empty() && mMeshes.worldComm().isMasterRank() )
        {
            fs::path thedir = mMeshes.rootRepository();//fs::path( fileSavePath ).parent_path();
            std::string fileNameMeshPath = (thedir/prefixvm(M_name,"mesh.path")).string();
            if ( !fs::exists(thedir))
                fs::create_directories(thedir);
            std::ofstream file( fileNameMeshPath.c_str(), std::ios::out|std::ios::trunc);
            file << M_meshFilename;
            file.close();
        }
    } // if ( !M_mesh )

    // update data by mesh entity
    for ( auto & [name,data] : M_codbme )
    {
        data.setMesh( M_mesh );
        data.template updateForUse<MeshType>();
    }
}


template class ModelMesh<uint32_type>;

template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<2,1>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<3,1>>>( ModelMeshes<uint32_type> const& );
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 2 )
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<2,2>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<3,2>>>( ModelMeshes<uint32_type> const& );
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 3 )
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<2,3>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<3,3>>>( ModelMeshes<uint32_type> const& );
#endif
#if BOOST_PP_GREATER_EQUAL( FEELPP_MESH_MAX_ORDER, 4 )
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<2,4>>>( ModelMeshes<uint32_type> const& );
template void ModelMesh<uint32_type>::updateForUse<Mesh<Simplex<3,4>>>( ModelMeshes<uint32_type> const& );
#endif
} // namespace FeelModels
} // namespace Feel