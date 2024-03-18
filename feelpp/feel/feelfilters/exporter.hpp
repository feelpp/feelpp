/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-11-09

  Copyright (C) 2004 EPFL
  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2016 Feel++ Consortium

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
#ifndef FEELPP_FILTERS_EXPORTER_H
#define FEELPP_FILTERS_EXPORTER_H

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/visitor.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

#include <feel/feeldiscr/timeset.hpp>
#include <feel/feelfilters/enums.hpp>

namespace Feel
{

/**
 * \class Exporter
 * \brief export Feel generated data to some file formats
 * \ingroup Exporter
 *
 * \tparam MeshType     mesh type
 * \tparam N            mesh geometrical order
 *
 * Use the visitor and factory pattern.
 *
 * Here is a snippet on how to use the Exporter class
 * \code
 * #include <feel/feelfilters/exporter.hpp>
 * typedef Exporter<mesh_type> export_type;
 * typedef std::shared_ptr<export_type> export_ptrtype;
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
class FEELPP_EXPORT Exporter
    :
        public CommObject,
        public VisitorBase,
        public Visitor<MeshType>
{
public:


    /** @name Typedefs
     */
    //@{
    using super = CommObject;
    typedef VisitorBase super1;
    typedef Visitor<MeshType> super2;

    typedef Exporter<MeshType,N> etype;
    typedef std::shared_ptr<etype> ptrtype;
    typedef Feel::detail::TimeSet<MeshType,N> timeset_type;
    typedef typename timeset_type::mesh_type mesh_type;
    typedef typename timeset_type::mesh_ptrtype mesh_ptrtype;
    typedef std::shared_ptr<timeset_type> timeset_ptrtype;
    typedef std::vector<timeset_ptrtype> timeset_set_type;
    typedef typename timeset_set_type::iterator timeset_iterator;
    typedef typename timeset_set_type::const_iterator timeset_const_iterator;
    typedef typename timeset_type::step_type step_type;
    typedef typename timeset_type::step_ptrtype step_ptrtype;
    typedef typename timeset_type::step_set_type step_set_type;

    typedef typename mesh_type::index_type index_type;
    struct Factory
    {
        typedef Feel::Singleton< Feel::Factory< Exporter<MeshType,N>, std::string > > type;
    };

  protected :
    typedef std::vector<std::pair<timeset_ptrtype,step_set_type>> steps_write_on_disk_type;
  public :

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     */
    Exporter( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

    /**
     * Constructor
     * \param type string containing the type of exporter (gmsh, ensight,...)
     * \param prefix the prefix for the file names of the exported data
     * \param freq an integer giving the frequency at which the data should be saved
     */
    Exporter( std::string const& type,
              std::string const& prefix = "",
              int freq = 1,
              worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

    /**
     * Constructor
     * \param vm \p variables_map containing the type of exporter and other exporter options
     * \param prefix the prefix for the file names of the exported data
     * \param freq an integer giving the frequency at which the data should be saved
     */
    Exporter( po::variables_map const& vm,
              std::string const& exporter_prefix = "",
              worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() ) FEELPP_DEPRECATED;

    /**
     * Constructor
     * \param prefix the prefix for the file names of the exported data
     * \param freq an integer giving the frequency at which the data should be saved
     */
    Exporter( std::string const& exporter_prefix,
              worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

    /**
     * copy constructor
     */
    Exporter( Exporter const & exporter );

    /**
     * destructor
     */
    ~Exporter() override;

    /**
     * Static function instantiating from the Exporter Factory an exporter out
     * of the \p exportername and using \p prefix for the prefix of the data
     * files.
     */
    static std::shared_ptr<Exporter<MeshType,N> > New( std::string const& exportername,
                                                         std::string prefix,
                                                         worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

    /**
     * Static function instantiating from the Exporter Factory an exporter out
     * of the variables_map \p vm and using \p prefix for the prefix of the data
     * files.
     */
    static std::shared_ptr<Exporter<MeshType,N> > New( po::variables_map const& vm = Environment::vm(),
                                                         std::string prefix = Environment::about().appName(),
                                                         worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() ) FEELPP_DEPRECATED;
    /**
     * Static function instantiating from the Exporter Factory an exporter using
     * \p prefix for the prefix of the data files.
     */
    static std::shared_ptr<Exporter<MeshType,N> > New( std::string prefix,
                                                         worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );

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
     * \return the file type format (ASCII or BINARY)
     */
    FileType fileType() const
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

        if(M_ts_set.size() > 0)
        {
            M_ts_set.back()->setName( M_prefix );
        }

        return this;
    }

    /**
     * get the prefix to \p __prefix
     */
    std::string getPrefix()
    {
        return M_prefix;
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
    Exporter<MeshType,N>* setFileType( FileType __ft )
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
        CHECK( !M_ts_set.empty() ) << "time set is empty";
        return M_ts_set.front();
    }
    timeset_ptrtype timeSet( int ts )
        {
            CHECK( 0 <= ts && ts < M_ts_set.size() ) << "invalid time set index " << ts;
            return M_ts_set[ts];
        }
    /**
     * @return true if export by parts if true, false otherwise
     */
    bool exportByParts() const
    {
        return M_export_by_parts;
    }
    /**
     * set the export by parts to \p b
     */
    void setExportByParts( bool b )
    {
        M_export_by_parts = b;
    }
    /**
     * @return true if use single transient file if true, false otherwise
     */
    bool useSingleTransientFile() const
        {
            return M_use_single_transient_file;
        }
    /**
     * set the use single transient file to \p s
     */
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


    //! export a scalar quantity \p u with name \p name
    //! \param name name of the scalar quantity
    //! \param u scalar quantity to be exported
    //! \param cst true if the scalar is constant over time, false otherwise
    template<typename T>
    void
    add( std::string const& name, T const& u, bool cst = false,
         typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr )
        {
            this->step( 0 )->add( name, u, cst );
        }

    //! export field \p u with name \p name. If both nodal and element representations are present
    //! in \p reps, the suffixes _n and _e are added to the name.
    //! \param name name of field exported
    //! \param u the field (element of function space)
    //! \param reps representation of the field exported. It can be a string (nodal or element) or set of string
    template<typename F>
    void
    add( std::string const& name, F const& u, typename step_type::variant_representation_arg_type reps = "",
         typename std::enable_if<is_functionspace_element_v<F>>::type* = nullptr )
        {
            this->step( 0 )->add( name, u, reps );
        }

    //! export expression \p expr with name \p name. If both nodal and element representations are present
    //! in \p reps, the suffixes _n and _e are added to the name.
    //! \param name name of field exported
    //! \param expr a Feel++ expression (scalar, vectorial of dim 2 or 3,  square matrix 2*2 or 3*3 )
    //! \param reps representation of the field exported. It can be a string (nodal or element) or set of string
    template<typename ExprT>
    void
    add( std::string const& name, ExprT const& expr, typename step_type::variant_representation_arg_type reps = "",
         typename std::enable_if_t<std::is_base_of_v<ExprBase,ExprT> >* = nullptr )
    {
        this->step( 0 )->add( name, expr, reps );
    }

    //! export expression \p expr with name \p name which is defined on \p rangeElt. If both nodal and element
    //! representations are present in \p reps, the suffixes _n and _e are added to the name.
    //! \param name name of field exported
    //! \param expr a Feel++ expression (scalar, vectorial of dim 2 or 3,  square matrix 2*2 or 3*3 )
    //! \param rangeElt collection of mesh element
    //! \param reps representation of the field exported. It can be a string (nodal or element) or set of string
    template<typename ExprT, typename EltWrapperT = elements_reference_wrapper_t<mesh_type>>
    void
    add( std::string const& name, ExprT const& expr, EltWrapperT const& rangeElt, typename step_type::variant_representation_arg_type reps = "",
         typename std::enable_if_t<std::is_base_of_v<ExprBase,ExprT> && is_filter_v<EltWrapperT> >* = nullptr )
    {
        this->step( 0 )->add( name, expr, rangeElt, reps );
    }

    //! export the mesh partitioning
    void
    addRegions()
        {
            this->step( 0 )->addRegions( "" );
        }

    /**
     * @return the step shared_ptr at time \p time
     */
    step_ptrtype step( double time )
    {
        CHECK( !M_ts_set.empty() ) << "timeset is empty";
        return this->step( time, M_ts_set.size() -1 );
    }
    step_ptrtype step( double time, int s )
    {
        CHECK( s >= 0 && s < M_ts_set.size() ) << "invalid timeset index " << s;
        timeset_ptrtype __ts = M_ts_set[s];
        return __ts->step( time, this->freq() );
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * add the timeset \p __ts to the Exporter
     */
    uint16_type addTimeSet( timeset_ptrtype const& __ts )
    {
        if ( __ts )
        {
            M_ts_set.push_back( __ts );
            return M_ts_set.size()-1;
        }
        return invalid_v<uint16_type>;
    }
    //! add the timeset with name \p __tsname to the Exporter
    //! if the name is empty, use the prefix exporter
    uint16_type addTimeSet( std::string const& tsname = "" )
    {
        std::string tsnameUsed = tsname.empty()? this->prefix() : tsname;
        return this->addTimeSet( timeset_ptrtype( new timeset_type( tsnameUsed ) ) );
    }

    //! save timeset in memory on disk
    void save() const
    {
        if ( !this->worldComm().isActive() )
            return;

        bool hasStepToWrite = false;
        steps_write_on_disk_type stepsToWriteOnDisk;
        timeset_const_iterator __ts_it = this->beginTimeSet();
        timeset_const_iterator __ts_en = this->endTimeSet();
        for ( ; __ts_it != __ts_en ; ++__ts_it )
        {
            timeset_ptrtype __ts = *__ts_it;
            auto steps = __ts->stepsToWriteOnDisk();
            stepsToWriteOnDisk.push_back( std::make_pair( __ts, steps ) );
            if ( !steps.empty() || (__ts->numberOfSteps() == 0 && __ts->hasMesh() ) )
                hasStepToWrite = true;
        }

        if ( hasStepToWrite )
            this->save( stepsToWriteOnDisk );

        for ( auto & [__ts, steps ] : stepsToWriteOnDisk )
        {
            for ( auto & step : steps )
            {
                step->setState( STEP_ON_DISK );
                step->cleanup();
            }
            // save metadata file (sould be done even the if step is not write on the disk
            std::string filename = (fs::path(this->path()) / (prefix()+".timeset")).string();
            __ts->save( filename, this->worldComm() );
        }
    }

    //!
    //! serve results though a webserver
    //!
    virtual void serve() const;

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
            if ( !fs::exists( filename ) )
                return;
            ( *__ts_it )->load( filename,__time );
        }
    }

    ExporterGeometry exporterGeometry() const { return M_ex_geometry; }
    //@}
protected:

    /**
     * this p save function is defined by the Exporter subclasses and implement
     * saving the data to files
     */
    virtual void save( steps_write_on_disk_type const& stepsToWriteOnDisk ) const = 0;


    bool M_do_export;
    bool M_export_by_parts;
    bool M_use_single_transient_file;
    std::string M_type;
    std::string M_prefix;
    int M_freq;
    FileType M_ft;
    std::string M_path;
    ExporterGeometry M_ex_geometry;

    mutable timeset_set_type M_ts_set;
};



template <typename ... Ts>
auto exporter( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && mesh = args.get(_mesh);
    bool fileset = args.get_else_invocable(_fileset,[](){ return boption(_name="exporter.fileset"); } );
    bool byparts = args.get_else_invocable(_byparts,[](){ return boption(_name="exporter.byparts"); } );
    std::string const& name = args.get_else_invocable(_name,[](){ return Environment::about().appName(); } );
    std::string const& geo = args.get_else_invocable(_geo, [](){ return soption(_name="exporter.geometry"); } );
    auto && path = args.get_else_invocable(_path, [&name](){ return std::string((fs::path(Environment::exportsRepository())/fs::path(soption("exporter.format"))/name).string()); } );

    using mesh_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(mesh)>>>;
    using exporter_type = Exporter<mesh_type,mesh_type::nOrder>;

    auto e =  exporter_type::New( name,mesh->worldCommPtr() );
    e->setPrefix( name );
    e->setUseSingleTransientFile( fileset );
    e->setExportByParts( byparts );
    if ( std::string(geo).compare("change_coords_only") == 0 )
        e->setMesh( mesh, EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY );
    else if ( std::string(geo).compare("change") == 0 )
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
    typedef std::shared_ptr<type> ptrtype;
};

}
} // Feel

//#if !defined( FEELPP_INSTANTIATION_MODE )
# include <feel/feelfilters/exporterimpl.hpp>
//#endif // FEELPP_INSTANTIATION_MODE


#endif /* __Exporter_H */
