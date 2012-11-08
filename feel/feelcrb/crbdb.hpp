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
    CRBDB();
    CRBDB( std::string prefixdir,
           std::string name,
           std::string dbprefix,
           po::variables_map const& vm );
    //! copy constructor
    CRBDB( CRBDB const & );
    //! destructor
    virtual ~CRBDB();

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    CRBDB& operator=( CRBDB const & o )
    {
        if ( this != &o )
        {
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    //! \return prefix directory
    std::string const& prefixDirectory() const
    {
        return M_prefixdir;
    }

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

    //! \return the db local path
    fs::path dbLocalPath() const;

    //! \return the db system path
    fs::path dbSystemPath() const;

    //! \return path to database, empty path if not found
    fs::path lookForDB() const;

    //! \return \c variables_map
    po::variables_map vm()
    {
        return M_vm;
    }

    //! \return \c variables_map
    po::variables_map vm() const
    {
        return M_vm;
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
    void setDBFilename( std::string const& filename )
    {
        M_dbfilename = filename;
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

    std::string M_prefixdir;
    std::string M_name;
    std::string M_dbfilename;
    po::variables_map M_vm;
    bool M_isloaded;


};
}
#endif /* __CRBDB_H */
