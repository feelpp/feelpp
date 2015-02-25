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

template<typename SpaceType>
class CRBElementsDB : public CRBDB
{

    typedef  CRBDB super;

public :
    typedef SpaceType space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    //! element of the functionspace type
    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;

    //! mesh type
    typedef typename space_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef std::vector<element_type> wn_type;

    //! constructors
    CRBElementsDB()
    :
        super()
    {
    }

    CRBElementsDB( std::string prefixdir,
                std::string name,
                std::string dbprefix )
    :
        super( prefixdir,
               name,
               dbprefix),
        M_N( 0 )
    {}


    //! destructor
    ~CRBElementsDB()
    {}


    /**
     * save the database
     */
    void saveDB();

    /**
     * load the database
     */
    bool loadDB();


    boost::tuple<wn_type, wn_type> wn()
    {
        return boost::make_tuple( M_WN , M_WNdu );
    }

    void setMN( size_type MN )
    {
        M_N = MN;
    }

    void setWn( boost::tuple< wn_type, wn_type > WN )
    {
        auto primal = WN.template get<0>();
        auto dual = WN.template get<1>();
        M_WN = primal;
        M_WNdu = dual;
    }

private :

    friend class boost::serialization::access;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    size_type M_N;

    wn_type M_WN;
    wn_type M_WNdu;
};//class CRBElementsDB

template<typename SpaceType>
void
CRBElementsDB<SpaceType>::saveDB()
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

template<typename SpaceType>
bool
CRBElementsDB<SpaceType>::loadDB()
{

    bool rebuild_db = boption(_name="crb.rebuild-database");
    int Nrestart = ioption(_name="crb.restart-from-N");
    if ( rebuild_db && Nrestart < 1 )
        return false;

    if( this->isDBLoaded() )
        return true;

    fs::path db = this->lookForDB();

    if ( db.empty() )
        return false;

    if ( !fs::exists( db ) )
        return false;

    //std::cout << "Loading " << db << "...\n";
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

    return false;
}


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
