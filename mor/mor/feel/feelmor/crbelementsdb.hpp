/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-06-16

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file crbelementsdb.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-06-16
 */
#ifndef __CRBElementsDB_H
#define __CRBElementsDB_H 1

#include <string>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/fusion/algorithm.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelmor/crbdb.hpp>
#include <feel/feelcore/hdf5.hpp>

namespace Feel
{

template<typename ModelType>
class CRBElementsDB : public CRBDB
{

    typedef  CRBDB super;

public :

    typedef ModelType model_type;
    typedef std::shared_ptr<model_type> model_ptrtype;

    //! element of the functionspace type
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;

    //! mesh type
    typedef typename model_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    //! function space type
    typedef typename model_type::space_type space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;

    typedef typename model_type::rbfunctionspace_type rbfunctionspace_type;
    typedef std::shared_ptr<rbfunctionspace_type> rbfunctionspace_ptrtype;

    typedef typename rbfunctionspace_type::rb_basis_type wn_type;

    //! constructors
    CRBElementsDB( std::string const& name = "defaultname_crbelementdb",
                   std::string const& ext = "elements",
                   worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
    :
        super( name, ext, worldComm ),
        M_fileFormat( soption(_name="crb.db.format") ),
        M_useMonolithicRbSpace( true ),
        M_N( 0 )
    {
#ifndef FEELPP_HAS_HDF5
        if ( M_fileFormat == "hdf5" )
        {
            LOG(INFO) << "CRB db format hdf5 unsupported. Switching to boost.";
            M_fileFormat = "boost";
        }
#endif
        if ( M_fileFormat == "hdf5" )
            this->setDBFilename( fmt::format( "{}.{}.h5",
                                   this->name(), ext ) );
        else
            this->setDBFilename( fmt::format( "{}.{}_p{}.crbdb",
                                   this->name(),
                                   ext,
                                   this->worldComm().globalRank()));
    }

    CRBElementsDB( std::string const& name,
                   std::string const& ext,
                   model_ptrtype const & model )
        :
        CRBElementsDB( name, ext )
        {
            M_model = model;
        }

    //! destructor
    ~CRBElementsDB() override
    {}

    void setup( boost::property_tree::ptree const& ptree, std::string const& dbDir )
        {
            CHECK( M_rbSpace ) << "no rbspace";
            if ( !M_rbSpace->mesh/*functionSpace*/() )
                M_rbSpace->setup( ptree, dbDir );
            size_type rbdim = ptree.template get<int>( "dimension" );
            this->setMN( rbdim );
            std::string dbname = ptree.template get<std::string>( "database-filename" );
            if ( this->worldComm().globalSize() > 1 && M_fileFormat != "hdf5" )
                dbname.replace(dbname.end()-std::string("_p0.crbdb").size(), dbname.end(),
                               (boost::format("_p%1%.crbdb")%this->worldComm().globalRank()).str() );
            fs::path dbnamePath = fs::path( dbname );
            this->setDBFilename( dbnamePath.filename().string() );
            
            // get function space name (needed if nproc offline != nproc online)
            std::string functionspaceName = ptree.template get<std::string>( "functionspace-filename" );
            fs::path fsnamePath = fs::path( functionspaceName );
            this->setFSFilename( fsnamePath.filename().string() );

            if ( dbnamePath.is_absolute() )
                this->setDBDirectory( dbnamePath.parent_path().string() );
            else if ( !dbDir.empty() )
                this->setDBDirectory( dbDir );
            this->setIsLoaded( false );
            CHECK( this->loadDB() ) << "loading of crb basis function fails";
        }


    /**
     * save the database
     */
    void saveDB() override;

    /**
     * save a new element to the database
     */
    void saveNewElementToDB() override;

    //!
    //! 
    //!
    bool loadDB() override;

    //!
    //! 
    //!
    void loadDB( std::string const& filename, crb::load l ) override {}
    
#ifdef FEELPP_HAS_HDF5
    /**
     * save the database in hdf5 format
     */
    void saveHDF5DB();

    /**
     * save a new element to the database in hdf5 format
     */
    void saveHDF5newElementToDB();

    /**
     * load the database in hdf5 format
     */
    void loadHDF5DB();
#endif

    boost::tuple<wn_type, wn_type> wn() const
    {
        //return boost::make_tuple( M_WN , M_WNdu );
        CHECK( M_rbSpace ) << "no reduced basis space defined";
        return boost::make_tuple( M_rbSpace->primalRB(), M_rbSpace->dualRB() );
    }

    bool useMonolithicRbSpace() const
    {
        return M_useMonolithicRbSpace;
    }

    void setUseMonolithicRbSpace( bool b )
    {
        M_useMonolithicRbSpace = (space_type::nSpaces>1)? b : true;
    }

    void setMN( size_type MN )
    {
        M_N = MN;
    }

    void setWn( boost::tuple< wn_type, wn_type > const& WN )
    {
        M_rbSpace->setPrimalBasis( WN.template get<0>() );
        M_rbSpace->setDualBasis( WN.template get<1>() );
    }

    void setModel( model_ptrtype const& model )
    {
        M_model = model;
        this->setDBDirectory( M_model->uuid() );
        M_rbSpace = model->rBFunctionSpace();
        M_useMonolithicRbSpace = model->useMonolithicRbSpace();
    }

    std::string dbFileFormat()
    {
        return M_fileFormat;
    }

private :

    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const;

    template<class Archive>
    void load( Archive & ar, const unsigned int version ) ;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    std::string M_fileFormat;

    bool M_useMonolithicRbSpace;

    size_type M_N;

    rbfunctionspace_ptrtype M_rbSpace;

    model_ptrtype M_model;


};//class CRBElementsDB

template<typename ModelType>
void
CRBElementsDB<ModelType>::saveNewElementToDB()
{    
    if ( M_fileFormat == "boost" )
    {        
        this->saveDB();     
    }
    else
    {        
        this->saveHDF5newElementToDB();     
    }
}



template<typename ModelType>
void
CRBElementsDB<ModelType>::saveDB()
{
    if ( M_fileFormat == "hdf5" )
    {
#ifdef FEELPP_HAS_HDF5
        this->saveHDF5DB();
#else
        CHECK(false) << "Feel++ not compiled with hdf5";
#endif
    }
    else
    /* save in boost format by default */
    {
        auto p = this->dbLocalPath() / this->dbFilename();
        if( this->worldComm().isMasterRank() )
            std::cout << "CRBElementsDB::saveDB : " << p << std::endl;
        fs::ofstream ofs( p );

        if ( ofs )
        {
            boost::archive::binary_oarchive oa( ofs );
            // write class instance to archive
            oa << *this;
            // archive and stream closed when destructors are called
        }
    }
}

template<typename ModelType>
bool
CRBElementsDB<ModelType>::loadDB()
{
#if 0
    bool rebuild_db = boption(_name="crb.rebuild-database");
    int Nrestart = ioption(_name="crb.restart-from-N");
    if ( rebuild_db && Nrestart < 1 )
        return false;
#endif

    if( this->isDBLoaded() )
        return true;

    fs::path db = this->lookForDB();

    // if ( db.empty() )
    //     return false;
        
    // Check if the first primal basis element exists
    std::ostringstream hdf5File_primal;
    fs::path p = this->dbLocalPath() / fs::path(this->dbFilename());
    p.replace_extension("");
    std::string additional_part_primal = std::string("-primal-") +  std::to_string(1);
    hdf5File_primal << p.string() + additional_part_primal << ".h5";
    LOG(INFO) << fmt::format("hdf5File_primal {}",hdf5File_primal.str());
    
    if( (!fs::exists( db ) || db.empty() ) && !fs::exists(hdf5File_primal.str()))   
        return false;    

    if ( M_fileFormat == "hdf5" )
    {
#ifdef FEELPP_HAS_HDF5
        this->loadHDF5DB();
#else
        CHECK(false) << "Feel++ not compiled with hdf5";
#endif
        std::cout << "Loading " << db << " done...\n";
        this->setIsLoaded( true );
        return true;
    }
    else
    {
        std::cout << "Loading " << db << "...\n";
        fs::ifstream ifs( db );

        if ( ifs )
        {
            boost::archive::binary_iarchive ia( ifs );
            // write class instance to archive
            ia >> *this;
            //std::cout << "Loading " << db << " done...\n";
            this->setIsLoaded( true );
            // archive and stream closed when destructors are called
            return true;
        }
    }

    return false;
}



namespace detail
{
template<class Archive,typename RbSpaceType>
struct SaveDatabaseCompositeByBlock
{
    SaveDatabaseCompositeByBlock( Archive & ar, RbSpaceType const& rbSpace )
        :
        M_ar( ar ),
        M_rbSpace( rbSpace )
    {}
    template <typename T>
    void operator()( T & x ) const
    {
        typedef typename T::first_type key_type;

        auto subRbSpace = M_rbSpace.template rbFunctionSpace<key_type::value>();
        auto const& WN = subRbSpace->primalRB();
        auto const& WNdu = subRbSpace->dualRB();
        size_type numberOfPrimalBasis = WN.size();
        size_type numberOfDualBasis = WNdu.size();
        M_ar & BOOST_SERIALIZATION_NVP( numberOfPrimalBasis );
        M_ar & BOOST_SERIALIZATION_NVP( numberOfDualBasis );
        for( size_type i=0; i<numberOfPrimalBasis; i++ )
            M_ar & BOOST_SERIALIZATION_NVP( unwrap_ptr( WN[i] ) );
        for( size_type i=0; i<numberOfDualBasis; i++ )
            M_ar & BOOST_SERIALIZATION_NVP( unwrap_ptr( WNdu[i] ) );
    }
    Archive & M_ar;
    RbSpaceType const& M_rbSpace;
};

template<class Archive,typename RbSpaceType>
void
saveDatabaseCompositeByBlock( Archive & ar, RbSpaceType const& rbSpace, typename std::enable_if<RbSpaceType::element_type::is_composite>::type* = nullptr )
{
    boost::fusion::for_each( rbSpace.rbfunctionspaces(), SaveDatabaseCompositeByBlock<Archive,RbSpaceType>( ar,rbSpace ) );
}
template<class Archive,typename RbSpaceType>
void
saveDatabaseCompositeByBlock( Archive & ar, RbSpaceType const& rbSpace, typename std::enable_if<!RbSpaceType::element_type::is_composite>::type* = nullptr )
{}

template<class Archive,typename RbSpaceType>
struct LoadDatabaseCompositeByBlock
{
    LoadDatabaseCompositeByBlock( Archive & ar, RbSpaceType const& rbSpace )
        :
        M_ar( ar ),
        M_rbSpace( rbSpace )
    {}
    template <typename T>
    void operator()( T & x ) const
    {
        typedef typename T::first_type key_type;

        auto subRbSpace = M_rbSpace.template rbFunctionSpace<key_type::value>();

        size_type numberOfPrimalBasis = 0, numberOfDualBasis = 0;
        M_ar & BOOST_SERIALIZATION_NVP( numberOfPrimalBasis );
        M_ar & BOOST_SERIALIZATION_NVP( numberOfDualBasis );

        size_type numberOfPrimalBasisLoaded = numberOfPrimalBasis;//std::min( numberOfPrimalBasis,M_N );
        size_type numberOfDualBasisLoaded = numberOfDualBasis;//std::min( numberOfDualBasis,M_N );

        if ( subRbSpace->dimension() < numberOfPrimalBasisLoaded )
            subRbSpace->setDimension( numberOfPrimalBasisLoaded );
        auto & WN = subRbSpace->primalRB();
        auto & WNdu = subRbSpace->dualRB();
        if ( WN.size() < numberOfPrimalBasisLoaded )
            WN.resize( numberOfPrimalBasisLoaded );
        if ( WNdu.size() < numberOfDualBasisLoaded )
            WNdu.resize( numberOfDualBasisLoaded );

        CHECK( subRbSpace->functionSpace() ) << "rbspace does not defined a fespace";
        auto Xh = subRbSpace->functionSpace();

        for( size_type i = 0 ; i < numberOfPrimalBasisLoaded ; i++ )
        {
            auto & wni = WN[i];
            if ( !wni )
                wni = Xh->elementPtr( (boost::format( "fem-primal-%1%" ) % ( i ) ).str() );
            M_ar & BOOST_SERIALIZATION_NVP( unwrap_ptr( wni ) );
        }
        for( size_type i = 0 ; i < numberOfDualBasisLoaded ; i++ )
        {
            auto & wndui = WNdu[i];
            if ( !wndui )
                wndui = Xh->elementPtr( (boost::format( "fem-dual-%1%" ) % ( i ) ).str() );
            M_ar & BOOST_SERIALIZATION_NVP( unwrap_ptr( wndui ) );
        }
    }
    Archive & M_ar;
    RbSpaceType const& M_rbSpace;
};

template<class Archive,typename RbSpaceType>
void
loadDatabaseCompositeByBlock( Archive & ar, RbSpaceType const& rbSpace, typename std::enable_if<RbSpaceType::element_type::is_composite>::type* = nullptr )
{
    boost::fusion::for_each( rbSpace.rbfunctionspaces(), LoadDatabaseCompositeByBlock<Archive,RbSpaceType>( ar,rbSpace ) );
}
template<class Archive,typename RbSpaceType>
void
loadDatabaseCompositeByBlock( Archive & ar, RbSpaceType const& rbSpace, typename std::enable_if<!RbSpaceType::element_type::is_composite>::type* = nullptr )
{}

} // namespace detail

template<typename ModelType>
template<class Archive>
void
CRBElementsDB<ModelType>::save( Archive & ar, const unsigned int version ) const
{
#if 0
    auto mesh = mesh_type::New();
    auto is_mesh_loaded = mesh->load( _name="mymesh",_path=this->dbLocalPath(),_type="binary" );

    if ( ! is_mesh_loaded )
    {
        auto first_element = M_WN[0];
        mesh = first_element.functionSpace()->mesh() ;
        mesh->save( _name="mymesh",_path=this->dbLocalPath(),_type="binary" );
    }
#endif
    ar & BOOST_SERIALIZATION_NVP( M_useMonolithicRbSpace );
    if ( M_useMonolithicRbSpace )
    {
        auto const& WN = M_rbSpace->primalRB();
        auto const& WNdu = M_rbSpace->dualRB();
        size_type numberOfPrimalBasis = WN.size();
        size_type numberOfDualBasis = WNdu.size();
        LOG( INFO ) << "saving Elements DB";
        ar & BOOST_SERIALIZATION_NVP( numberOfPrimalBasis );
        ar & BOOST_SERIALIZATION_NVP( numberOfDualBasis );
        for( size_type i=0; i<numberOfPrimalBasis; i++ )
            ar & BOOST_SERIALIZATION_NVP( unwrap_ptr( WN[i] ) );
        for( size_type i=0; i<numberOfDualBasis; i++ )
            ar & BOOST_SERIALIZATION_NVP( unwrap_ptr( WNdu[i] ) );
    }
    else
    {
        Feel::detail::saveDatabaseCompositeByBlock( ar,*M_rbSpace );
    }
    LOG( INFO ) << "Elements DB saved";
}

#ifdef FEELPP_HAS_HDF5

template<typename ModelType>
void
CRBElementsDB<ModelType>::saveHDF5newElementToDB()
{
    auto & WN = M_rbSpace->primalRB();
    auto & WNdu = M_rbSpace->dualRB();

    int size = WN.size();

    std::ostringstream hdf5File_primal,hdf5File_dual; 
    std::string additional_part_primal = std::string("-primal-") +  std::to_string(size);
    std::string additional_part_dual = std::string("-dual-") +  std::to_string(size);
    fs::path p = this->dbLocalPath() / fs::path(this->dbFilename());
    p.replace_extension("");
    hdf5File_primal << p.string() + additional_part_primal << ".h5";
    hdf5File_dual << p.string() + additional_part_dual << ".h5";

    LOG(INFO) << fmt::format("String hdf 5 file primal {}",p.string() + additional_part_primal);
    LOG(INFO) << fmt::format("String hdf 5 file dual {}",p.string() + additional_part_dual);

    HDF5 hdf5_primal, hdf5_dual;

    hsize_t dims[1];
    hsize_t offset[1];

    dims[0] = 1;
    offset[0] = 0;

    /* only do this on proc 0 */
    if(this->worldComm().isMasterRank())
    {
        /* If a previous db already exists, we remove it */
        if( boost::filesystem::exists( hdf5File_primal.str() ) )
        {
            boost::filesystem::remove( hdf5File_primal.str() );
        }
        if( boost::filesystem::exists( hdf5File_dual.str() ) )
        {
            boost::filesystem::remove( hdf5File_dual.str() );
        }

        /* Select the correct size for data */
        hid_t memDataType;
        if(sizeof(int) == 4)
        { memDataType = H5T_NATIVE_INT; }
        else
        { memDataType = H5T_NATIVE_LLONG; }

        hdf5_primal.openFile( hdf5File_primal.str(), this->worldComm().subWorldCommSeq(), true, true );
        hdf5_primal.createTable( "dbSize", memDataType, dims, 1 );
        hdf5_primal.write( "dbSize", memDataType, dims, offset, &size, 1 );
        hdf5_primal.closeTable( "dbSize" );     
        hdf5_primal.closeFile();

        hdf5_dual.openFile( hdf5File_dual.str(), this->worldComm().subWorldCommSeq(), true, true );
        hdf5_dual.createTable( "dbSize", memDataType, dims, 1 );
        hdf5_dual.write( "dbSize", memDataType, dims, offset, &size, 1 );
        hdf5_dual.closeTable( "dbSize" );        
        hdf5_dual.closeFile();
    }
    this->worldComm().barrier();

    LOG( INFO ) << "saving HDF5 Elements DB";

    std::ostringstream tableName;
    LOG( INFO ) << hdf5File_primal.str();
    tableName << "WN[" << size-1 << "]";
    WN[size-1]->saveHDF5(hdf5File_primal.str(), tableName.str(), true);

    tableName.str("");
    tableName << "WNdu[" << size-1 << "]";
    WNdu[size-1]->saveHDF5(hdf5File_dual.str(), tableName.str(), true);

    LOG( INFO ) << "Elements DB saved in hdf5";



}

template<typename ModelType>
void
CRBElementsDB<ModelType>::saveHDF5DB()
{
    auto & WN = M_rbSpace->primalRB();
    auto & WNdu = M_rbSpace->dualRB();

    int size = WN.size();

    std::ostringstream hdf5File;
    fs::path p = this->dbLocalPath() / fs::path(this->dbFilename());
    p.replace_extension("");
    hdf5File << p.string() << ".h5";

    HDF5 hdf5;
    hsize_t dims[1];
    hsize_t offset[1];

    dims[0] = 1;
    offset[0] = 0;

    /* only do this on proc 0 */
    if(this->worldComm().isMasterRank())
    {
        /* If a previous db already exists, we remove it */
        if( boost::filesystem::exists( hdf5File.str() ) )
        {
            boost::filesystem::remove( hdf5File.str() );
        }

        /* Select the correct size for data */
        hid_t memDataType;
        if(sizeof(int) == 4)
        { memDataType = H5T_NATIVE_INT; }
        else
        { memDataType = H5T_NATIVE_LLONG; }

        hdf5.openFile( hdf5File.str(), this->worldComm().subWorldCommSeq(), true, true );
        hdf5.createTable( "dbSize", memDataType, dims, 1 );
        hdf5.write( "dbSize", memDataType, dims, offset, &size, 1 );
        hdf5.closeTable( "dbSize" );
        hdf5.closeFile();
    }
    this->worldComm().barrier();

    LOG( INFO ) << "saving HDF5 Elements DB";
    for(int i=0; i < size; i++)
    {
        std::ostringstream tableName;
        LOG( INFO ) << hdf5File.str();
        tableName << "WN[" << i << "]";
        WN[i]->saveHDF5(hdf5File.str(), tableName.str(), true);

        tableName.str("");
        tableName << "WNdu[" << i << "]";
        WNdu[i]->saveHDF5(hdf5File.str(), tableName.str(), true);
    }
    LOG( INFO ) << "Elements DB saved in hdf5";
}
#endif

template<typename ModelType>
template<class Archive>
void
CRBElementsDB<ModelType>::load( Archive & ar, const unsigned int version )
{
    LOG( INFO ) << " loading Elements DB ... ";

    //mesh_ptrtype mesh;
    space_ptrtype Xh;

    if ( M_rbSpace && M_rbSpace->functionSpace() )
        Xh = M_rbSpace->functionSpace();
    else
    {
        if ( !M_model )
        {
            LOG(INFO) << "[load] model not initialized, loading fdb files...\n";
            auto mesh = mesh_type::New();
            bool is_mesh_loaded = mesh->load( _name="mymesh",_path=this->dbLocalPath(),_type="binary" );
            Xh = space_type::New( mesh );
            LOG(INFO) << "[load] loading fdb files done.\n";
        }
        else
        {
            LOG(INFO) << "[load] get mesh/Xh from model...\n";
            //mesh = M_model->functionSpace()->mesh();
            Xh = M_model->functionSpace();
            LOG(INFO) << "[load] get mesh/Xh from model done.\n";
        }
    }

    ar & BOOST_SERIALIZATION_NVP( M_useMonolithicRbSpace );

    if ( M_useMonolithicRbSpace )
    {
        size_type numberOfPrimalBasis = 0, numberOfDualBasis = 0;
        ar & BOOST_SERIALIZATION_NVP( numberOfPrimalBasis );
        ar & BOOST_SERIALIZATION_NVP( numberOfDualBasis );

        size_type numberOfPrimalBasisLoaded = std::min( numberOfPrimalBasis,M_N );
        size_type numberOfDualBasisLoaded = std::min( numberOfDualBasis,M_N );

        if ( M_rbSpace->dimension() < numberOfPrimalBasisLoaded )
            M_rbSpace->setDimension( numberOfPrimalBasisLoaded );
        auto & WN = M_rbSpace->primalRB();
        auto & WNdu = M_rbSpace->dualRB();

        LOG( INFO ) << "loading Elements DB (boost)";
        for( size_type i = 0 ; i < numberOfPrimalBasisLoaded ; i++ )
        {
            auto & wni = WN[i];
            if ( !wni )
                wni = Xh->elementPtr( (boost::format( "fem-primal-%1%" ) % ( i ) ).str() );
            ar & BOOST_SERIALIZATION_NVP( unwrap_ptr( wni ) );
        }

        for( size_type i = 0 ; i < numberOfDualBasisLoaded ; i++ )
        {
            auto & wndui = WNdu[i];
            if ( !wndui )
                wndui = Xh->elementPtr( (boost::format( "fem-dual-%1%" ) % ( i ) ).str() );
            ar & BOOST_SERIALIZATION_NVP( unwrap_ptr( wndui ) );
        }
        M_rbSpace->updatePrimalBasisForUse();
        M_rbSpace->updateDualBasisForUse();
    }
    else
    {
        Feel::detail::loadDatabaseCompositeByBlock( ar,*M_rbSpace );
    }
    LOG( INFO ) << "Elements DB loaded";
}

#ifdef FEELPP_HAS_HDF5
template<typename ModelType>
void
CRBElementsDB<ModelType>::loadHDF5DB()
{
    LOG( INFO ) << " loading HDF5 Elements DB ... ";

    /* If we are loading a hdf5 database */
    /* We passed a dummy databased in the load archive */
    /* so we have to find the right path where they are located */

    /* build the filename of db 0 */
    std::ostringstream hdf5File;
    fs::path dbpath = this->dbLocalPath();
    fs::path p = dbpath / fs::path(this->dbFilename());
    p.replace_extension("");
    hdf5File << p.string() << ".h5";

    if (false)//boost::filesystem::exists( hdf5File.str() ) )
    {
        HDF5 hdf5;
        hsize_t dims[1];
        hsize_t offset[1];

        dims[0] = 1;
        offset[0] = 0;

        int size;
        /* Select the correct size for data */
        hid_t memDataType;
        if(sizeof(int) == 4)
        { memDataType = H5T_NATIVE_INT; }
        else
        { memDataType = H5T_NATIVE_LLONG; }

        hdf5.openFile( hdf5File.str(), this->worldComm(), true, false );
        hdf5.openTable( "dbSize", dims);
        hdf5.read( "dbSize", memDataType, dims, offset, &size, 1 );
        hdf5.closeTable( "dbSize");
        hdf5.closeFile();

        this->setMN(size);

        auto & WN = M_rbSpace->primalRB();
        auto & WNdu = M_rbSpace->dualRB();

        WN.resize( M_N );
        WNdu.resize( M_N );

        space_ptrtype Xh;
        if ( M_rbSpace && M_rbSpace->functionSpace() )
            Xh = M_rbSpace->functionSpace();
        else
        {
            if ( !M_model )
            {
                LOG(INFO) << "[load] model not initialized, loading fdb files...\n";
                auto mesh = mesh_type::New();
                bool is_mesh_loaded = mesh->load( _name="mymesh",_path=this->dbLocalPath(),_type="binary" );
                Xh = space_type::New( mesh );
                LOG(INFO) << "[load] loading fdb files done.\n";
            }
            else
            {
                LOG(INFO) << "[load] get mesh/Xh from model...\n";
                auto mesh = M_model->functionSpace()->mesh();
                Xh = M_model->functionSpace();
                LOG(INFO) << "[load] get mesh/Xh from model done.\n";
            }
        }
        // read the functionspace file 
        fs::path functionspacePath = dbpath / fs::path(this->fsFilename());
        std::string spacefileName = functionspacePath.string();
        LOG(INFO) << fmt::format("Reading functionspace file {}",spacefileName);    
        std::optional<std::vector<index_type>> spacesRelation;    
        auto ch = (spacefileName.substr(0, spacefileName.find_last_of("."))).back();
        std::string str_nbproc(1,ch);
        int nproc_db = std::stoi(str_nbproc);    
        LOG(INFO) << fmt::format("Number of processors for the saved database {}",nproc_db);
        // load the relation between current dof table and on database file
        if (  nproc_db !=this->worldComm().globalSize() )
        {
            if(auto *file = fopen(spacefileName.c_str(), "r"))
            {
                fclose(file);
                spacesRelation = M_rbSpace->functionSpace()->relationFromFile( spacefileName );                
            }
        }

        LOG( INFO ) << "loading Elements DB (hdf5)";
        for(int i=0; i<M_N; i++)
        {
            std::ostringstream tableName;
            tableName << "WN[" << i << "]";
            auto & wni = WN[i];
            if ( !wni )
                wni = Xh->elementPtr( (boost::format( "fem-primal-%1%" ) % ( i ) ).str() );
            
            // if nproc offline != nproc online
            if(this->worldComm().globalSize() != nproc_db)
                wni->loadHDF5(hdf5File.str(),spacesRelation, tableName.str());
            else
                wni->loadHDF5(hdf5File.str(), tableName.str());

            tableName.str("");
            tableName << "WNdu[" << i << "]";
            auto & wndui = WNdu[i];
            if ( !wndui )
                wndui = Xh->elementPtr( (boost::format( "fem-dual-%1%" ) % ( i ) ).str() );
            
            // if nproc offline != nproc online
            if(this->worldComm().globalSize() != nproc_db)
                wndui->loadHDF5(hdf5File.str(),spacesRelation, tableName.str());
            else
                wndui->loadHDF5(hdf5File.str(), tableName.str());
        }
        M_rbSpace->updatePrimalBasisForUse();
        M_rbSpace->updateDualBasisForUse();
        LOG( INFO ) << "Elements DB loaded";
    }
    else
    {
        std::ostringstream hdf5File;
        fs::path dbpath = this->dbLocalPath();
        fs::path p = dbpath / fs::path(this->dbFilename());
        p.replace_extension("");
        hdf5File << p.string() << ".h5";

        
        HDF5 hdf5_primal, hdf5_dual;
        hsize_t dims[1];
        hsize_t offset[1];

        dims[0] = 1;
        offset[0] = 0;

        LOG(INFO) << fmt::format("Size of the reduced basis {}",M_N);

        auto & WN = M_rbSpace->primalRB();
        auto & WNdu = M_rbSpace->dualRB();

        WN.resize( M_N );
        WNdu.resize( M_N );

        space_ptrtype Xh;
        if ( M_rbSpace && M_rbSpace->functionSpace() )
            Xh = M_rbSpace->functionSpace();
        else
        {
            if ( !M_model )
            {
                LOG(INFO) << "[load] model not initialized, loading fdb files...\n";
                auto mesh = mesh_type::New();
                bool is_mesh_loaded = mesh->load( _name="mymesh",_path=this->dbLocalPath(),_type="binary" );
                Xh = space_type::New( mesh );
                LOG(INFO) << "[load] loading fdb files done.\n";
            }
            else
            {
                LOG(INFO) << "[load] get mesh/Xh from model...\n";
                auto mesh = M_model->functionSpace()->mesh();
                Xh = M_model->functionSpace();
                LOG(INFO) << "[load] get mesh/Xh from model done.\n";
            }
        }

        // read the functionspace file 
        fs::path functionspacePath = dbpath / fs::path(this->fsFilename());
        std::string spacefileName = functionspacePath.string();
        LOG(INFO) << fmt::format("Reading functionspace file {}",spacefileName);    
        std::optional<std::vector<index_type>> spacesRelation;    
        auto ch = (spacefileName.substr(0, spacefileName.find_last_of("."))).back();
        std::string str_nbproc(1,ch);
        int nproc_db = std::stoi(str_nbproc);    
        LOG(INFO) << fmt::format("Number of processors for the saved database {}",nproc_db);
        // load the relation between current dof table and on database file
        if (  nproc_db !=this->worldComm().globalSize() )
        {
            if(auto *file = fopen(spacefileName.c_str(), "r"))
            {
                fclose(file);
                spacesRelation = M_rbSpace->functionSpace()->relationFromFile( spacefileName );                
            }
        }

        for(int i=0; i<M_N; i++)
        {         
            std::ostringstream hdf5File_primal,hdf5File_dual;    
            std::string additional_part_primal = std::string("-primal-") +  std::to_string(i+1);
            std::string additional_part_dual = std::string("-dual-") +  std::to_string(i+1);
            fs::path p = this->dbLocalPath() / fs::path(this->dbFilename());
            p.replace_extension("");
            hdf5File_primal << p.string() + additional_part_primal << ".h5";
            hdf5File_dual << p.string() + additional_part_dual << ".h5";

            LOG(INFO) << fmt::format("String hdf 5 reading file primal {}",p.string() + additional_part_primal);
            LOG(INFO) << fmt::format("String hdf 5 reading file dual {}",p.string() + additional_part_dual);

            std::ostringstream tableName;
            tableName << "WN[" << i << "]";
            auto & wni = WN[i];
            if ( !wni )
                wni = Xh->elementPtr( (boost::format( "fem-primal-%1%" ) % ( i ) ).str() );
            
            // if nproc offline != nproc online
            if(this->worldComm().globalSize() != nproc_db)
                wni->loadHDF5(hdf5File_primal.str(),spacesRelation, tableName.str());
            else
                wni->loadHDF5(hdf5File_primal.str(), tableName.str());

            tableName.str("");
            tableName << "WNdu[" << i << "]";
            auto & wndui = WNdu[i];
            if ( !wndui )
                wndui = Xh->elementPtr( (boost::format( "fem-dual-%1%" ) % ( i ) ).str() );
            
            // if nproc offline != nproc online
            if(this->worldComm().globalSize() != nproc_db)
                wndui->loadHDF5(hdf5File_dual.str(),spacesRelation, tableName.str());
            else
                wndui->loadHDF5(hdf5File_dual.str(), tableName.str());

        }
        M_rbSpace->updatePrimalBasisForUse();
        M_rbSpace->updateDualBasisForUse();
        LOG( INFO ) << "Elements DB split loaded";
    }

}
#endif

}//Feel



namespace boost
{
namespace serialization
{
template< typename T>
struct version< Feel::CRBElementsDB<T> >
{
    // at the moment the version of the CRB DB is 0. if any changes is done
    // to the format it is mandatory to increase the version number below
    // and use the new version number of identify the new entries in the DB
    typedef mpl::int_<0> type;
    typedef mpl::integral_c_tag tag;
    static const unsigned int value = version::type::value;
};
template<typename T> const unsigned int version<Feel::CRBElementsDB<T> >::value;
}
}

#endif /* __CRBElementsDB_H */
