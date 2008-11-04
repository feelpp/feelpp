/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2004-11-09

  Copyright (C) 2004 EPFL
  Copyright (C) 2007-2008 Université Joseph Fourier (Grenoble I)

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
   \file Exporter.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2004-11-09
 */
#ifndef __Exporter_H
#define __Exporter_H 1

#include <life/lifecore/life.hpp>
#include <life/lifecore/visitor.hpp>
#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>
#include <life/lifediscr/timeset.hpp>

namespace Life
{
/**
 * \enum
 */
enum file_type
    {
        ASCII  = 0,
        BINARY = 1
    };

/**
 * \class Exporter
 * \brief export Life generated data to some file formats
 *
 * use the visitor pattern
 *
 * \ingroup Exporter
 * @author Christophe Prud'homme
 */
template<typename MeshType>
class Exporter
    :
        public VisitorBase,
        public Visitor<MeshType>
{
public:


    /** @name Typedefs
     */
    //@{
    typedef VisitorBase super1;
    typedef Visitor<MeshType> super2;

    typedef TimeSet<MeshType> timeset_type;
    typedef boost::shared_ptr<timeset_type> timeset_ptrtype;
    typedef std::set<timeset_ptrtype> timeset_set_type;
    typedef typename timeset_set_type::iterator timeset_iterator;
    typedef typename timeset_set_type::const_iterator timeset_const_iterator;

    struct Factory
    {
        typedef Life::Singleton< Life::Factory< Exporter<MeshType>, std::string > > type;
    };
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Exporter( std::string const& type, std::string const& prefix = "", int freq = 1 );

    Exporter( po::variables_map const& vm, std::string const& exporter_prefix = "" );

    Exporter( Exporter const & exporter );

    virtual ~Exporter();

    static Exporter<MeshType>* New( std::string const& exporter )
    {
        return Factory::type::instance().createObject( exporter );
    }

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the type of exporter
     */
    std::string const& type() const { return M_type; }

    /**
     * \return the prefix of the exported file names
     */
    std::string const& prefix() const { return M_prefix; }

    /**
     * \return the frequency at which the results are saved
     */
    int freq() const { return M_freq; }

    /**
     * \return the file type format (ASCII or BINARY)
     */
    file_type fileType() const { return M_ft; }

    //@}

    /** @name  Mutators
     */
    //@{

    virtual Exporter<MeshType>* setOptions( po::variables_map const& vm, std::string const& exp_prefix = "" );

    Exporter<MeshType>* setType( std::string const& __type )
    {
        M_type = __type;
        return this;
    }
    Exporter<MeshType>* setPrefix( std::string const& __prefix )
    {
        M_prefix = __prefix;
        return this;
    }

    Exporter<MeshType>* setFreq( int __freq )
    {
        M_freq = __freq;
        return this;
    }

    Exporter<MeshType>* setFileType( file_type __ft )
    {
        M_ft = __ft;
        return this;
    }

    timeset_iterator beginTimeSet() { return M_ts_set.begin(); }
    timeset_iterator endTimeSet() { return M_ts_set.end(); }

    timeset_const_iterator beginTimeSet() const { return M_ts_set.begin(); }
    timeset_const_iterator endTimeSet() const { return M_ts_set.end(); }


    //@}

    /** @name  Methods
     */
    //@{

    void addTimeSet( timeset_ptrtype const& __ts )
    {
        if ( __ts )
            M_ts_set.insert( __ts );
    }


    virtual void save() const = 0;



    //@}
protected:

    std::string M_type;
    std::string M_prefix;
    int M_freq;
    file_type M_ft;

    timeset_set_type M_ts_set;
};



po::options_description exporter_options( std::string const& prefix = "" );

}
#endif /* __Exporter_H */
