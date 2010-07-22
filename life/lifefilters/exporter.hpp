/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2004-11-09

  Copyright (C) 2004 EPFL
  Copyright (C) 2007-2008 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

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
 * \ingroup Exporter
 *
 * Use the visitor and factory pattern.
 *
 * Here is a snippet on how to use the Exporter class
 * \code
 * #include <life/lifefilters/exporter.hpp>
 * typedef Exporter<mesh_type> export_type;
 * typedef boost::shared_ptr<export_type> export_ptrtype;
 * // vm is a po::variables_map to get the command lines options
 * export_ptrtype exporter( export_type::New( vm );
 * // U is an element of a function space which we want to visualise
 * exporter->step(0)->setMesh( U.functionSpace()->mesh() );
 * exporter->step(0)->add( "u", U );
 * \endcode
 *
 * \sa Laplacian
 *
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
    typedef std::vector<timeset_ptrtype> timeset_set_type;
    typedef typename timeset_set_type::iterator timeset_iterator;
    typedef typename timeset_set_type::const_iterator timeset_const_iterator;
    typedef typename timeset_type::step_type step_type;
    typedef typename timeset_type::step_ptrtype step_ptrtype;
    struct Factory
    {
        typedef Life::Singleton< Life::Factory< Exporter<MeshType>, std::string > > type;
    };
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * Constructor
     * \param type string containing the type of exporter (gmsh, ensight,...)
     * \param prefix the prefix for the file names of the exported data
     * \param freq an integer giving the frequency at which the data should be saved
     */
    Exporter( std::string const& type, std::string const& prefix = "", int freq = 1 );

    /**
     * Constructor
     * \param vm \p variables_map containing the type of exporter and other exporter options
     * \param prefix the prefix for the file names of the exported data
     * \param freq an integer giving the frequency at which the data should be saved
     */
    Exporter( po::variables_map const& vm, std::string const& exporter_prefix = "" );

    /**
     * copy constructor
     */
    Exporter( Exporter const & exporter );

    /**
     * destructor
     */
    virtual ~Exporter();

    /**
     * Static function instantiating from the Exporter Factory an exporter out
     * of the \p exportername and using \p prefix for the prefix of the data
     * files.
     */
    static Exporter<MeshType>* New( std::string const& exportername, std::string prefix = "export" );

    /**
     * Static function instantiating from the Exporter Factory an exporter out
     * of the variables_map \p vm and using \p prefix for the prefix of the data
     * files.
     */
    static Exporter<MeshType>* New( po::variables_map const& vm, std::string prefix = "export" );

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return true if doing the export, false otherwise
     */
    bool doExport() const { return M_do_export; }

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

    /**
     * \return the path to the saved files
     */
    std::string path() const { return M_path; }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the doExport to \p do_export
     */
    void setDoExport( bool do_export ) {  M_do_export = do_export; }

    /**
     * set the options from the \p variables_map \p vm as well as the prefix \p
     * exp_prefix
     */
    virtual Exporter<MeshType>* setOptions( po::variables_map const& vm, std::string const& exp_prefix = "" );

    /**
     * set to \p __type the type of exporter (gmsh, ensight...)
     */
    Exporter<MeshType>* setType( std::string const& __type )
    {
        M_type = __type;
        return this;
    }

    /**
     * add an extra path to the current directory to save the data using the \p
     * boost::format object \p fmt
     */
    Exporter<MeshType>* addPath( boost::format fmt );

    /**
     * set the prefix to \p __prefix
     */
    Exporter<MeshType>* setPrefix( std::string const& __prefix )
    {
        M_prefix = __prefix;
        return this;
    }

    /**
     * set the save frequency to \p __freq
     */
    Exporter<MeshType>* setFreq( int __freq )
    {
        M_freq = __freq;
        return this;
    }

    /**
     * set the \p file type to \p __ft (binary or ascii)
     */
    Exporter<MeshType>* setFileType( file_type __ft )
    {
        M_ft = __ft;
        return this;
    }

    timeset_iterator beginTimeSet() { return M_ts_set.begin(); }
    timeset_iterator endTimeSet() { return M_ts_set.end(); }

    timeset_const_iterator beginTimeSet() const { return M_ts_set.begin(); }
    timeset_const_iterator endTimeSet() const { return M_ts_set.end(); }


    timeset_ptrtype defaultTimeSet() { return M_ts_set.front(); }
    timeset_ptrtype timeSet( int ts ) { return M_ts_set[ts]; }

    step_ptrtype step( double time ) { return M_ts_set.back()->step( time ); }
    step_ptrtype step( double time, int s ) { return M_ts_set[s]->step( time ); }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * add the timeset \p __ts to the Exporter
     */
    void addTimeSet( timeset_ptrtype const& __ts )
    {
        if ( __ts )
            M_ts_set.push_back( __ts );
    }


    /**
     * this p save function is defined by the Exporter subclasses and implement
     * saving the data to files
     */
    virtual void save() const = 0;



    //@}
protected:

    bool M_do_export;
    std::string M_type;
    std::string M_prefix;
    int M_freq;
    file_type M_ft;
    std::string M_path;

    mutable timeset_set_type M_ts_set;
};



po::options_description exporter_options( std::string const& prefix = "" );

} // Life

//#if !defined( LIFE_INSTANTIATION_MODE )
//# include <life/lifefilters/exporter.cpp>
//#endif // LIFE_INSTANTIATION_MODE

#endif /* __Exporter_H */
