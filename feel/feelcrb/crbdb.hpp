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

namespace Feel
{
/**
 * \class CRBDB
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class CRBDB
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    CRBDB( std::string const& name = "defaultname_crbdb", WorldComm const& worldComm = Environment::worldComm() );

    //! copy constructor
    CRBDB( CRBDB const & ) = default;
    //! destructor
    virtual ~CRBDB();

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

    //! \return the mpi communicators
    WorldComm const& worldComm() const { return M_worldComm; }

    //! \return name
    std::string const& name() const
    {
        return M_name;
    }

    //! \return the DB filename
    std::string const& dbFilename() const
    {
        return M_dbfilename;
    }

    //! \return prefix directory
    std::string const& dbDirectory() const
    {
        return M_dbDirectory;
    }

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

    //! \return the db system path
    fs::path dbSystemPath() const;

    //! \return path to database, empty path if not found
    virtual fs::path lookForDB() const;

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

    //! set the DB filename
    void setName( std::string const& name )
    {
        M_name = name;
    }

    //! set the DB filename
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
        M_dbDirectory = ( fs::path( M_dbDirectory )/fs::path( subdirectory ) ).string();
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
    //! mpi communicators
    WorldComm const& M_worldComm;

    std::string M_name;
    std::string M_dbfilename;
    std::string M_dbDirectory;
    std::string M_dbSubDirectory;
    bool M_isloaded;


};
}
#endif /* __CRBDB_H */
