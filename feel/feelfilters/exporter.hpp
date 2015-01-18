/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-11-09

  Copyright (C) 2004 EPFL
  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2015 Feel++ Consortium

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-11-09
 */
#ifndef __Exporter_H
#define __Exporter_H 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/visitor.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

#include <feel/feeldiscr/timeset.hpp>

namespace Feel
{
enum ExporterGeometry
{
    EXPORTER_GEOMETRY_STATIC = 0,
    EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY = 1,
    EXPORTER_GEOMETRY_CHANGE = 2

};
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
 * \brief export Feel generated data to some file formats
 * \ingroup Exporter
 *
 * Use the visitor and factory pattern.
 *
 * Here is a snippet on how to use the Exporter class
 * \code
 * #include <feel/feelfilters/exporter.hpp>
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
template<typename MeshType, int N = 1>
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

    typedef Exporter<MeshType,N> etype;
    typedef boost::shared_ptr<etype> ptrtype;
    typedef TimeSet<MeshType,N> timeset_type;
    typedef typename timeset_type::mesh_type mesh_type;
    typedef typename timeset_type::mesh_ptrtype mesh_ptrtype;
    typedef boost::shared_ptr<timeset_type> timeset_ptrtype;
    typedef std::vector<timeset_ptrtype> timeset_set_type;
    typedef typename timeset_set_type::iterator timeset_iterator;
    typedef typename timeset_set_type::const_iterator timeset_const_iterator;
    typedef typename timeset_type::step_type step_type;
    typedef typename timeset_type::step_ptrtype step_ptrtype;
    
    struct Factory
    {
        typedef Feel::Singleton< Feel::Factory< Exporter<MeshType,N>, std::string > > type;
    };
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     */
    Exporter( WorldComm const& worldComm = Environment::worldComm() );

    /**
     * Constructor
     * \param type string containing the type of exporter (gmsh, ensight,...)
     * \param prefix the prefix for the file names of the exported data
     * \param freq an integer giving the frequency at which the data should be saved
     */
    Exporter( std::string const& type,
              std::string const& prefix = "",
              int freq = 1,
              WorldComm const& worldComm = Environment::worldComm() );

    /**
     * Constructor
     * \param vm \p variables_map containing the type of exporter and other exporter options
     * \param prefix the prefix for the file names of the exported data
     * \param freq an integer giving the frequency at which the data should be saved
     */
    Exporter( po::variables_map const& vm,
              std::string const& exporter_prefix = "",
              WorldComm const& worldComm = Environment::worldComm() ) FEELPP_DEPRECATED;

    /**
     * Constructor
     * \param prefix the prefix for the file names of the exported data
     * \param freq an integer giving the frequency at which the data should be saved
     */
    Exporter( std::string const& exporter_prefix,
              WorldComm const& worldComm = Environment::worldComm() );

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
    static boost::shared_ptr<Exporter<MeshType,N> > New( std::string const& exportername,
                                                         std::string prefix,
                                                         WorldComm const& worldComm = Environment::worldComm() );

    /**
     * Static function instantiating from the Exporter Factory an exporter out
     * of the variables_map \p vm and using \p prefix for the prefix of the data
     * files.
     */
    static boost::shared_ptr<Exporter<MeshType,N> > New( po::variables_map const& vm = Environment::vm(),
                                                         std::string prefix = Environment::about().appName(),
                                                         WorldComm const& worldComm = Environment::worldComm() ) FEELPP_DEPRECATED;
    /**
     * Static function instantiating from the Exporter Factory an exporter using
     * \p prefix for the prefix of the data files.
     */    
    static boost::shared_ptr<Exporter<MeshType,N> > New( std::string prefix,
                                                         WorldComm const& worldComm = Environment::worldComm() );

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
    bool doExport() const
    {
        return M_do_export;
    }

    /**
     * \return the type of exporter
     */
    std::string const& type() const
    {
        return M_type;
    }

    /**
     * \return the prefix of the exported file names
     */
    std::string const& prefix() const
    {
        return M_prefix;
    }

    /**
     * \return the frequency at which the results are saved
     */
    int freq() const
    {
        return M_freq;
    }


    /**
     * \return the frequency at which the results are saved
     */
    int cptOfSave() const
    {
        return M_cptOfSave;
    }

    /**
     * \return the file type format (ASCII or BINARY)
     */
    file_type fileType() const
    {
        return M_ft;
    }

    /**
     * \return the path to the saved files
     */
    std::string path() const
    {
        return M_path;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the doExport to \p do_export
     */
    void setDoExport( bool do_export )
    {
        M_do_export = do_export;
    }

    /**
     * set the options from the \p variables_map \p vm as well as the prefix \p
     * exp_prefix
     */
    virtual Exporter<MeshType,N>* setOptions( std::string const& exp_prefix = "" );

    /**
     * set the options from the \p variables_map \p vm as well as the prefix \p
     * exp_prefix
     */
    virtual Exporter<MeshType,N>* setOptions( po::variables_map const& vm, std::string const& exp_prefix = "" ) FEELPP_DEPRECATED { return setOptions( exp_prefix ); }


    /**
     * set to \p __type the type of exporter (gmsh, ensight...)
     */
    Exporter<MeshType,N>* setType( std::string const& __type )
    {
        M_type = __type;
        return this;
    }

    /**
     * add an extra path to the current directory to save the data using the \p
     * boost::format object \p fmt
     */
    Exporter<MeshType,N>* addPath( boost::format fmt );

    /**
     * set the path
     */
    void setPath( std::string path )
    {
        if ( !fs::exists(path) && this->worldComm().isMasterRank() )
            fs::create_directories( path );
        // be sure that all process can find the path after
        this->worldComm().barrier();

        M_path = path;
    }

    /**
     * set the prefix to \p __prefix
     */
    Exporter<MeshType,N>* setPrefix( std::string const& __prefix )
    {
        M_prefix = __prefix;
        return this;
    }

    /**
     * set the save frequency to \p __freq
     */
    Exporter<MeshType,N>* setFreq( int __freq )
    {
        M_freq = __freq;
        return this;
    }

    /**
     * set the \p file type to \p __ft (binary or ascii)
     */
    Exporter<MeshType,N>* setFileType( file_type __ft )
    {
        M_ft = __ft;
        return this;
    }

    timeset_iterator beginTimeSet()
    {
        return M_ts_set.begin();
    }
    timeset_iterator endTimeSet()
    {
        return M_ts_set.end();
    }

    timeset_const_iterator beginTimeSet() const
    {
        return M_ts_set.begin();
    }
    timeset_const_iterator endTimeSet() const
    {
        return M_ts_set.end();
    }


    timeset_ptrtype defaultTimeSet()
    {
        return M_ts_set.front();
    }
    timeset_ptrtype timeSet( int ts )
        {
            return M_ts_set[ts];
        }

    bool useSingleTransientFile() const
        {
            return M_use_single_transient_file;
        }
    void setUseSingleTransientFile( bool s )
        {
            M_use_single_transient_file = s;
        }
    void
    setMesh( mesh_ptrtype mesh, ExporterGeometry exgeo = EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY )
        {
            M_ex_geometry = exgeo;
            M_ts_set.back()->setMesh( mesh );
            //this->step( 0 )->setMesh( mesh );
        }
    template<typename F>
    void
    add( std::string const& name, F const& u )
        {
            this->step( 0 )->add( this->prefix()+"."+name, u );
        }

    void
    addRegions()
        {
            this->step( 0 )->addRegions( this->prefix() );
        }
    step_ptrtype step( double time )
    {
        if ( this->cptOfSave() % this->freq()  )
            return M_ts_set.back()->step( time, true );

        else
            return M_ts_set.back()->step( time, false );
    }
    step_ptrtype step( double time, int s )
    {
        if ( M_cptOfSave % this->freq()  )
            return M_ts_set[s]->step( time, true );

        else
            return M_ts_set[s]->step( time, false );
    }

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

    /**
     * save in a file set of time which are been exported
     */
    void saveTimeSet() const
    {
        if ( ! this->worldComm().isMasterRank() )
            {
                ++M_cptOfSave;
                return;
            }

        auto __ts_it = this->beginTimeSet();
        auto __ts_en = this->endTimeSet();

        for ( ; __ts_it != __ts_en ; ++__ts_it )
        {
            auto filename = this->path()+"/"+prefix()+".timeset";
            ( *__ts_it )->save( filename );
        }

        ++M_cptOfSave;
    }

    /**
     * reload from file set of time which are been exported
     */
    void restart( double __time )
    {
        auto __ts_it = this->beginTimeSet();
        auto __ts_en = this->endTimeSet();

        for ( ; __ts_it != __ts_en ; ++__ts_it )
        {
            auto filename = this->path()+"/"+prefix()+".timeset";
            ( *__ts_it )->load( filename,__time );

            M_cptOfSave = ( *__ts_it )->numberOfTotalSteps()+1;
        }
    }

    void setWorldComm( WorldComm const& wc )
        {
            M_worldComm = wc;
        }
    WorldComm const& worldComm() const
        {
            return M_worldComm;
        }

    ExporterGeometry exporterGeometry() const { return M_ex_geometry; }
    //@}
protected:

    WorldComm M_worldComm;

    bool M_do_export;
    bool M_use_single_transient_file;
    std::string M_type;
    std::string M_prefix;
    int M_freq;
    mutable int M_cptOfSave;
    file_type M_ft;
    std::string M_path;
    ExporterGeometry M_ex_geometry;

    mutable timeset_set_type M_ts_set;
};




namespace detail
{
template<typename Args>
struct compute_exporter_return
{
    typedef typename boost::remove_pointer<
        typename boost::remove_const<
            typename boost::remove_reference<
                typename parameter::binding<Args, tag::mesh>::type
                >::type
            >::type
        >::type::element_type mesh_type;
    //typename Feel::vf::detail::clean_type<Args, tag::mesh>::type::element_type mesh_type;
    //typedef typename parameter::value_type<Args, tag::order>::type order_type;
    //typedef boost::shared_ptr<Exporter<mesh_type,order_type::value> > type;
    typedef Exporter<mesh_type,mesh_type::nOrder> type;
    typedef boost::shared_ptr<type> ptrtype;
    //typedef boost::shared_ptr<Exporter<Mesh<Simplex<2> >,1> > type;

};
}
BOOST_PARAMETER_FUNCTION( ( typename Feel::detail::compute_exporter_return<Args>::ptrtype ),
                          exporter,                                       // 2. name of the function template
                          tag,                                        // 3. namespace of tag types
                          ( required                                  // 4. one required parameter, and
                            ( mesh, * )
                          ) // required
                          ( optional                                  // 4. one required parameter, and
                            ( fileset, *, boption(_name="exporter.fileset") )
                            ( order, *, mpl::int_<1>() )
                            ( name,  *, Environment::about().appName() )
                            ( geo,   *, soption(_name="exporter.geometry") )
                            ( path, *( boost::is_convertible<mpl::_,std::string> ), std::string(".") )
                          ) )
{
    typedef typename Feel::detail::compute_exporter_return<Args>::type exporter_type;
    auto e =  exporter_type::New( name,mesh->worldComm() );
    e->setPrefix( name );
    e->setUseSingleTransientFile( fileset );
    if ( geo == "change_coords_only" )
        e->setMesh( mesh, EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY );
    else if ( geo == "change" )
        e->setMesh( mesh, EXPORTER_GEOMETRY_CHANGE );
    else // default
        e->setMesh( mesh, EXPORTER_GEOMETRY_STATIC );
    e->setPath( path );
    // addRegions not work with transient simulation!
    //e->addRegions();
    return e;
    //return Exporter<Mesh<Simplex<2> >,1>::New();
}

namespace meta
{
template<typename MeshType, int N = 1>
struct Exporter
{
    typedef Feel::Exporter<MeshType,N> type;
    typedef boost::shared_ptr<type> ptrtype;
};

}
} // Feel

//#if !defined( FEELPP_INSTANTIATION_MODE )
# include <feel/feelfilters/exporterimpl.hpp>
//#endif // FEELPP_INSTANTIATION_MODE


#endif /* __Exporter_H */
