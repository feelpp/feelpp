/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-06-15

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
   \file crbdb.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-06-15
 */
#ifndef __CRBDB_H
#define __CRBDB_H 1

#include <string>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/serialization/version.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/commobject.hpp>
#include <feel/feelcrb/crbenums.hpp>

namespace Feel
{
/**
 * \class CRBDB
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class CRBDB : public CommObject
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{
    using super = CommObject;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    CRBDB() = delete;
    
    //! default constructor
    CRBDB( std::string const& name,
           std::string const& ext, 
           worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

    CRBDB( std::string const& name,
           std::string const& ext,
           uuids::uuid const& i,
           worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    
    //! copy constructor
    CRBDB( CRBDB const & ) = default;
    //! destructor
    ~CRBDB() override;

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    CRBDB& operator=( CRBDB const & o ) = default;

    //@}

    /** @name Accessors
     */
    //@{

    //! \return name
    std::string const& name() const
    {
        return M_name;
    }

    //!
    //! return extension
    //!
    std::string const& extension() const { return M_ext; }
    
    //!
    //! @return the uuid of the CRBDB
    //!
    uuids::uuid const& id() const
    {
        return M_uuid;
    }

    //!
    //! @return the uuid from a path to the json file
    //!
    uuids::uuid id( fs::path p ) const;

    //!
    //! @return the string representation of the uuid
    //!
    std::string idStr() const { return boost::lexical_cast<std::string>( M_uuid ); }
    
    //!
    //! set UUID 
    //!
    void setId( uuids::uuid const& i )
        {
            M_uuid = i;
            this->setDBDirectory( M_uuid );
        } 
        
    //! \return the DB json filename
    std::string jsonFilename() const
        {
            return fs::path( M_dbfilename ).stem().string()+".json";
        }
    //! \return the abslute DB json filename
    std::string absoluteJsonFilename() const
    {
        return (this->dbLocalPath() / fs::path(this->jsonFilename()) ).string();
    }
    //! \return the DB filename
    std::string const& dbFilename() const
    {
        return M_dbfilename;
    }
    //! \return the abslute DB filename
    std::string absoluteDbFilename() const
    {
        return (this->dbLocalPath() / fs::path(this->dbFilename()) ).string();
    }

    //! \return the DB filename with the proc number
    std::string dbFilenameProc() const
    {
        auto dbf = fs::path(M_dbfilename);
        return (boost::format("%1%_p%2%%3%")
                % dbf.stem().string()
                % this->worldCommPtr()->globalRank()
                % dbf.extension().string() ).str();
    }

    //! \return the abslute DB filename with the proc number
    std::string absoluteDbFilenameProc() const
    {
        return (this->dbLocalPath() / fs::path(this->dbFilenameProc()) ).string();
    }
    //! \return the mesh filename
    std::string meshFilename() const
    {
        return (boost::format("%1%_mesh_p%2%.json")
                % fs::path( M_dbfilename ).stem().string()
                % this->worldCommPtr()->size()).str();
    }
    //! \return the abslute mesh filename
    std::string absoluteMeshFilename() const
    {
        return (this->dbLocalPath() / fs::path(this->meshFilename()) ).string();
    }

    //! \return prefix directory
    std::string dbDirectory( bool withSubDir = true ) const
    {
        if ( withSubDir && !M_dbSubDirectory.empty() )
            return (fs::path( M_dbDirectory )/ M_dbSubDirectory ).string();
        return M_dbDirectory;
    }

    //!
    //! set the DB directory according to a name and a UUID
    //!
    void setDBDirectory( std::string const& name, uuids::uuid const& i );

    //!
    //! set the DB directory according to a UUID
    //!
    void setDBDirectory( uuids::uuid const& i );

    //! \return sub directory
    std::string const& dbSubDirectory() const
    {
        return M_dbSubDirectory;
    }

    //! \return relative path
    std::string dbRelativePath() const
    {
        return ( M_dbSubDirectory.empty() )? this->dbFilename() :
            (fs::path( M_dbSubDirectory ) / fs::path(this->dbFilename())).string();
    }

    //! \return the db local path
    virtual fs::path dbLocalPath() const;

    //! \return path to database, empty path if not found
    virtual fs::path lookForDB() const;

    /**
     * find Uuid corresponding to load mode l and set the DB directory
     * @param[in] l load mode to find the DB
     * @param[in] uid_file uuid or filename in case of load mode needing them
     * @return true if DB found, false otherwise
     */
    bool findDBUuid(int l, std::string const& uid_file = "");

    //! \return \c variables_map
    po::variables_map const& vm() const
    {
        return Environment::vm();
    }

    //! \return true if the DB has been loaded, false otherwise
    bool isDBLoaded() const
    {
        return M_isloaded;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //!
    //! set the name of the model for the DB
    //! @code
    //! this->setName( "heat" );
    //! @endcode
    //!
    void setName( std::string const& name )
    {
        M_name = name;
    }

    //!
    //! set the name of the model for the DB
    //! @code
    //! this->setDBFilename( "heat.crbdb" );
    //! @endcode
    //!
    void setDBFilename( std::string const& filename )
    {
        M_dbfilename = filename;
    }

    //! set DB directory
    void setDBDirectory( std::string const& directory )
    {
        M_dbDirectory = directory;
    }

    //! add a subdirectory to the database directory
    void addDBSubDirectory( std::string const& subdirectory )
    {
        if ( M_dbSubDirectory.empty() )
            M_dbSubDirectory = subdirectory;
        else
            M_dbSubDirectory = (fs::path( M_dbSubDirectory ) / fs::path(subdirectory )).string();
        //M_dbDirectory = ( fs::path( M_dbDirectory )/fs::path( subdirectory ) ).string();
    }

    //! set DB sub directory
    void setDBSubDirectory( std::string const& directory )
        {
            M_dbSubDirectory = directory;
        }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * save the CRB database
     */
    virtual void saveDB();

    /**
     * load the CRB database
     */
    virtual bool loadDB();

    //!
    //! 
    //!
    virtual void loadDB( std::string const& filename, crb::load l ) = 0;

    //!
    //! @return path to file \p f and check json extension
    //!
    fs::path db( std::string const& ) const;

    /**
     * check existence of db file
     * @param[in] filename The path of the db file
     * @return The uuid of the corresponding db
     */
    uuids::uuid idFromDBFilename( std::string const& filename ) const;

    //!
    //!
    //!
    virtual void loadDBFromId( std::string const& id, crb::load l = crb::load::rb, std::string const& root = Environment::rootRepository() ) ;

    /**
     * check existence of db file
     * @param[in] id The uuid of the db
     * @return The uuid of the corresponding db
     */
    uuids::uuid idFromId( std::string const& id ) const;

    //!
    //! @return fs::path from DB \p id
    //! check existence of json metadata file
    //! @param id db id
    //! @param root root repository for CRB DB
    //!
    fs::path dbFromId( std::string const& id, std::string const& root = Environment::rootRepository() ) const;
    
    /**
     * check existence of db file
     * @param[in] last What last file type
     * @return The uuid of the corresponding db or uuid_nil if nothing found
     */
    uuids::uuid idFromDBLast( crb::last last = crb::last::modified ) const;

    //!
    //! 
    //!
    virtual void loadDBLast( crb::last last = crb::last::modified, crb::load l = crb::load::rb, std::string const& root = Environment::rootRepository() );

    //!
    //! @return fs::path from DB \p last
    //! check existence of json metadata file
    //! @param last what last file type 
    //! @param root root repository for CRB DB
    //!
    fs::path dbLast( crb::last last = crb::last::modified, std::string const& root = Environment::rootRepository() ) const;
    
    //@}

protected:
    void setIsLoaded( bool isloaded )
    {
        M_isloaded = isloaded;
    }

    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const
    {}

    template<class Archive>
    void load( Archive & ar, const unsigned int version )
    {}

    BOOST_SERIALIZATION_SPLIT_MEMBER()

private:

    std::string M_name;
    std::string M_ext;
    uuids::uuid M_uuid;
    std::string M_dbfilename;
    std::string M_dbDirectory;
    std::string M_dbSubDirectory;
    bool M_isloaded;


};
}
#endif /* __CRBDB_H */
