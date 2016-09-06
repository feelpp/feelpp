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

#include <feel/feelcore/feel.hpp>

namespace Feel
{

template<typename ModelType>
class CRBElementsDB : public CRBDB
{

    typedef  CRBDB super;

public :

    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;

    //! element of the functionspace type
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;

    //! mesh type
    typedef typename model_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! function space type
    typedef typename model_type::space_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef typename model_type::rbfunctionspace_type rbfunctionspace_type;
    typedef typename model_type::rbfunctionspace_ptrtype rbfunctionspace_ptrtype;

    typedef std::vector<element_type> wn_type;

    //! constructors
    CRBElementsDB( std::string const& name = "defaultname_crbelementdb",
                   WorldComm const& worldComm = Environment::worldComm() )
    :
        super( name, worldComm ),
        M_fileFormat( soption(_name="crb.db.format") ),
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
            this->setDBFilename( ( boost::format( "%1%.h5" )
                                   %this->name() ).str() );
        else
            this->setDBFilename( ( boost::format( "%1%_p%2%.crbdb" )
                                   %this->name()
                                   %this->worldComm().globalRank()
                                   ).str() );
    }

    CRBElementsDB( std::string const& name,
                   model_ptrtype const & model )
        :
        CRBElementsDB( name )
        {
            M_model = model;
        }

    //! destructor
    ~CRBElementsDB()
    {}

    void setup( boost::property_tree::ptree const& ptree, std::string const& dbDir )
        {
            CHECK( M_rbSpace ) << "no rbspace";
            if ( !M_rbSpace->mesh/*functionSpace*/() )
                M_rbSpace->setup( ptree, dbDir );
            size_type rbdim = ptree.template get<int>( "dimension" );
            this->setMN( rbdim );
            std::string dbname = ptree.template get<std::string>( "database-filename" );
            this->setDBFilename( fs::path( dbname ).filename().string() );
            if ( dbDir.empty() )
                this->setDBDirectory( fs::path( dbname ).parent_path().string() );
            else
                this->setDBDirectory( dbDir );
            this->setIsLoaded( false );
            CHECK( this->loadDB() ) << "loading of crb basis function fails";
        }


    /**
     * save the database
     */
    void saveDB();

    /**
     * load the database
     */
    bool loadDB();

#ifdef FEELPP_HAS_HDF5
    /**
     * save the database in hdf5 format
     */
    void saveHDF5DB();

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
        M_rbSpace = model->rBFunctionSpace();
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

    size_type M_N;

    rbfunctionspace_ptrtype M_rbSpace;

    model_ptrtype M_model;


};//class CRBElementsDB



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
        fs::ofstream ofs( this->dbLocalPath() / this->dbFilename() );

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

    if ( db.empty() )
        return false;

    if ( !fs::exists( db ) )
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
    auto const& M_WN = M_rbSpace->primalRB();
    auto const& M_WNdu = M_rbSpace->dualRB();
    int sizepr = M_WN.size();
    int sizedu = M_WNdu.size();

    LOG( INFO ) << "saving Elements DB";
    for(int i=0; i<sizepr; i++)
        ar & BOOST_SERIALIZATION_NVP( M_WN[i] );
    for(int i=0; i<sizedu; i++)
        ar & BOOST_SERIALIZATION_NVP( M_WNdu[i] );
    LOG( INFO ) << "Elements DB saved";
}

#ifdef FEELPP_HAS_HDF5
template<typename ModelType>
void
CRBElementsDB<ModelType>::saveHDF5DB()
{
    auto & M_WN = M_rbSpace->primalRB();
    auto & M_WNdu = M_rbSpace->dualRB();

    int size = M_WN.size();

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
        tableName << "M_WN[" << i << "]";
        M_WN[i].saveHDF5(hdf5File.str(), tableName.str(), true);

        tableName.str("");
        tableName << "M_WNdu[" << i << "]";
        M_WNdu[i].saveHDF5(hdf5File.str(), tableName.str(), true);
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

    auto & M_WN = M_rbSpace->primalRB();
    auto & M_WNdu = M_rbSpace->dualRB();

    M_WN.resize( M_N );
    M_WNdu.resize( M_N );

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
    element_type temp = Xh->element();

    LOG( INFO ) << "loading Elements DB (boost)";
    for( int i = 0 ; i < M_N ; i++ )
    {
        temp.setName( (boost::format( "fem-primal-%1%" ) % ( i ) ).str() );
        ar & BOOST_SERIALIZATION_NVP( temp );
        M_WN[i] = temp;
    }

    for( int i = 0 ; i < M_N ; i++ )
    {
        temp.setName( (boost::format( "fem-dual-%1%" ) % ( i ) ).str() );
        ar & BOOST_SERIALIZATION_NVP( temp );
        M_WNdu[i] = temp;
    }
    M_rbSpace->updatePrimalBasisForUse();
    M_rbSpace->updateDualBasisForUse();
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
    /* If the path does not exist then the db are in the system path */
    if ( ! fs::exists( hdf5File.str() ) )
    {
        dbpath = this->dbSystemPath();
        p = dbpath / fs::path(this->dbFilename());
        p.replace_extension("");
        hdf5File.str("");
        hdf5File << p.string() << ".h5";
    }

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

    auto & M_WN = M_rbSpace->primalRB();
    auto & M_WNdu = M_rbSpace->dualRB();

    M_WN.resize( M_N );
    M_WNdu.resize( M_N );

    mesh_ptrtype mesh;
    space_ptrtype Xh;

    if ( !M_model )
    {
        LOG(INFO) << "[load] model not initialized, loading fdb files...\n";
        mesh = mesh_type::New();
        bool is_mesh_loaded = mesh->load( _name="mymesh",_path=this->dbLocalPath(),_type="binary" );
        Xh = space_type::New( mesh );
        LOG(INFO) << "[load] loading fdb files done.\n";
    }
    else
    {
        LOG(INFO) << "[load] get mesh/Xh from model...\n";
        mesh = M_model->functionSpace()->mesh();
        Xh = M_model->functionSpace();
        LOG(INFO) << "[load] get mesh/Xh from model done.\n";
    }

    element_type temp = Xh->element();

    LOG( INFO ) << "loading Elements DB (hdf5)";
    for(int i=0; i<M_N; i++)
    {
        std::ostringstream tableName;
        tableName << "M_WN[" << i << "]";
        temp.setName( (boost::format( "fem-primal-%1%" ) % ( i ) ).str() );
        temp.loadHDF5(hdf5File.str(), tableName.str());
        M_WN[i] = temp;

        tableName.str("");
        tableName << "M_WNdu[" << i << "]";
        temp.setName( (boost::format( "fem-dual-%1%" ) % ( i ) ).str() );
        temp.loadHDF5(hdf5File.str(), tableName.str());
        M_WNdu[i] = temp;
    }
    M_rbSpace->updatePrimalBasisForUse();
    M_rbSpace->updateDualBasisForUse();
    LOG( INFO ) << "Elements DB loaded";

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
