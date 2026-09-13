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

#include <typeinfo>
#include <cmath>
#include <iomanip>
#include <limits>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/visitor.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

#include <feel/feeldiscr/timeset.hpp>
#include <feel/feelfilters/enums.hpp>
#include <feel/feelfilters/exporterio.hpp>
#include <feel/feelmesh/meshfragmentation.hpp>

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
    using field_set_type = Feel::detail::ExportFieldSet<MeshType,N>;
    using field_set_ptrtype = std::shared_ptr<field_set_type>;
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

    /** @return Whether the storage strategy writes dataset fields once.
     * Other backends materialize the same immutable snapshot in temporal records.
     * This is a storage capability, not a restriction on exporter-level add().
     */
    virtual bool supportsNativeStaticFields() const { return false; }

    /** @brief Select payload aggregation on backends supporting this policy.
     * All ranks must select the same policy. Automatic preserves the backend
     * default; unsupported explicit policies throw rather than being ignored.
     */
    virtual void setIOPolicy( ExporterIOPolicy policy )
    {
        if ( policy != ExporterIOPolicy::Automatic )
            throw std::invalid_argument( "this exporter does not support an explicit I/O policy" );
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
        if (hasStaticFields() && __type!=M_type) throw std::logic_error("set backend before dataset fields");
        M_type = __type;
        return this;
    }

    /**
     * add an extra path to the current directory to save the data using the \p
     * boost::format object \p fmt
     */
    Exporter<MeshType,N>* addPath( boost::format fmt );

    /**
     * @brief Collectively create and validate the shared output directory.
     * Every rank performs the directory operation: an MPI barrier alone does
     * not invalidate another node's negative filesystem lookup cache.
     * @throws std::runtime_error if any rank cannot access the directory.
     */
    void setPath( std::string path )
    {
        validateStaticOutputPath(path);
        std::string error;
        try
        {
            fs::create_directories( path );
            if ( !fs::is_directory( path ) )
                error = "not a directory: " + path;
        }
        catch ( std::exception const& e ) { error = e.what(); }
        auto message=collectiveError(error);
        if (!message.empty()) throw std::runtime_error("Exporter::setPath: "+message);
        M_path = path;
    }

    /**
     * set the prefix to \p __prefix
     */
    Exporter<MeshType,N>* setPrefix( std::string const& __prefix )
    {
        if (hasStaticFields() && __prefix!=M_prefix) throw std::logic_error("cannot rename an exporter containing dataset fields");
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
     * @return mesh fragmentation used as mesh parts of exporter
     */
    MeshFragmentation<mesh_type> const& meshFragmentation() const
    {
        return M_meshFragmentation;
    }

    /**
     * set the mesh fragmentation used as mesh parts of exporter
     */
    template <typename MFT>
    void setMeshFragmentation( MFT && mf )
    {
        if (hasStaticFields()) throw std::logic_error("set mesh fragmentation before dataset fields");
        M_meshFragmentation = std::forward<MFT>( mf );
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
            if (hasStaticFields() && (mesh!=M_staticMesh || exgeo!=M_ex_geometry))
                throw std::logic_error("cannot replace dataset mesh or geometry policy");
            M_mesh = mesh;
            M_ex_geometry = exgeo;
            M_ts_set.back()->setMesh( mesh );
            //this->step( 0 )->setMesh( mesh );
        }


    /** @brief Snapshot a time-independent per-case scalar on the exporter.
     * @param name Dataset-level name, reserved against temporal fields.
     * @param u Value copied now.
     * @param cst Compatibility flag: must be true; use step(t)->add for temporal scalars.
     */
    template<typename T>
    void add(std::string const& name,T const& u,bool cst=true,
             typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr)
    {
        std::ostringstream schema;
        schema << "scalar:" << std::setprecision(std::numeric_limits<T>::max_digits10) << u;
        std::string error=cst ? "" : "use step(t)->add for time-dependent scalars";
        if (!std::isfinite(static_cast<scalar_type>(u))) error="dataset constants must be finite";
        registerDatasetField(name,"",schema.str(),error,
                             [&](auto& fields) { fields.add(sanitize(name),u,true); });
    }

    /** @brief Snapshot an FE field independently of all temporal time sets.
     * @param name Dataset-level name, sanitized as for Step::add.
     * @param u Field copied now; later source mutation cannot change this snapshot.
     * @param reps Nodal/element representation(s), with the same defaults as Step::add.
     *
     * Collective registration precedes all steps and saves. No temporal step is
     * manufactured by this call. Native Gold writes once; packed Gold and other
     * backends materialize this same snapshot in each output record. Mesh/layout
     * changes require a new exporter. See step(t)->add for transient fields.
     */
    template<typename F>
    void add(std::string const& name,F const& u,typename step_type::variant_representation_arg_type reps="",
             typename std::enable_if<is_functionspace_element_v<decay_type<Feel::remove_shared_ptr_type<std::remove_pointer_t<F>>>>>::type* = nullptr)
    {
        std::string error;
        if constexpr (is_ptr_or_shared_ptr<F>::value)
            if (!u) error="null dataset field";
        validateDatasetError(error);
        auto const& field=unwrap_ptr(u);
        registerDatasetField(name,reps,typeid(decay_type<decltype(field)>).name(),validateDatasetFunction(field),
                             [&](auto& fields) { fields.add(sanitize(name),field,reps); });
    }

    /** @brief Snapshot an expression evaluated once on the exporter mesh.
     * @param name Dataset field name.
     * @param expr Expression evaluated at registration, not at later save calls.
     * @param reps Nodal/element representation(s), as for Step::add.
     */
    template<typename ExprT>
    void add(std::string const& name,ExprT const& expr,typename step_type::variant_representation_arg_type reps="",
             typename std::enable_if_t<std::is_base_of_v<ExprBase,ExprT>>* = nullptr)
    {
        registerDatasetField(name,reps,typeid(ExprT).name(),"",
                             [&](auto& fields) { fields.add(sanitize(name),expr,reps); });
    }

    /** @brief Snapshot an expression over an exporter-mesh range.
     * @param name Dataset field name.
     * @param expr Expression evaluated once; values outside the range remain zero.
     * @param rangeElt Element/boundary range belonging to the exporter mesh.
     * @param reps Nodal/element representation(s), as for Step::add.
     */
    template<typename ExprT,typename EltWrapperT=Range<mesh_type,MESH_ELEMENTS>>
    void add(std::string const& name,ExprT const& expr,EltWrapperT const& rangeElt,
             typename step_type::variant_representation_arg_type reps="",
             typename std::enable_if_t<std::is_base_of_v<ExprBase,ExprT> && is_filter_v<EltWrapperT>>* = nullptr)
    {
        std::string error=rangeElt.mesh()==M_mesh.get() ? "" : "expression range must use the exporter mesh";
        registerDatasetField(name,reps,typeid(ExprT).name(),error,
                             [&](auto& fields) { fields.add(sanitize(name),expr,rangeElt,reps); });
    }

    /** @brief Snapshot mesh partition IDs as the dataset field "pid". */
    void addRegions()
    {
        registerDatasetField("pid","element","partition-id","",[](auto& fields) { fields.addRegions(""); });
    }

    /** @return Whether this exporter owns time-independent field snapshots. */
    bool hasStaticFields() const { return bool(M_staticFields); }

    /** @return Time-neutral dataset storage; it has no Step or TimeSet reference. */
    field_set_ptrtype const& staticFields() const { return M_staticFields; }

    /** @return Whether the dataset snapshot has been successfully materialized. */
    bool staticFieldsWritten() const { return M_staticWritten; }


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
        M_staticOnlyAdapter=false;
        CHECK( s >= 0 && s < M_ts_set.size() ) << "invalid timeset index " << s;
        timeset_ptrtype __ts = M_ts_set[s];
        auto result = __ts->step(time,this->freq());
        if (hasStaticFields()) result->fieldSet()->reserveNames(M_staticNames);
        return result;
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

        if (hasStaticFields())
        {
            validateStaticFields();
            validateStaticOutputPath(this->path());
            // Temporal-only formats need a compatibility output record, but
            // registration itself never constructs a Step or references a TimeSet.
            if (!supportsNativeStaticFields() && M_ts_set.front()->numberOfSteps()==0)
            {
                M_ts_set.front()->step(0,this->freq());
                M_staticOnlyAdapter=true;
            }
        }
        bool hasStepToWrite=false;
        steps_write_on_disk_type stepsToWriteOnDisk;
        for (auto const& ts : M_ts_set)
        {
            auto steps=ts->stepsToWriteOnDisk();
            stepsToWriteOnDisk.emplace_back(ts,steps);
            if (!steps.empty() || (ts->numberOfSteps()==0 && ts->hasMesh())) hasStepToWrite=true;
        }
        if (hasStepToWrite)
        {
            bool const materialize=hasStaticFields() && !supportsNativeStaticFields();
            auto exposeSnapshot=[&](bool attach) {
                if (materialize)
                    for (auto const& [ts,steps] : stepsToWriteOnDisk)
                        for (auto const& step : steps) step->fieldSet()->materialize(*M_staticFields,attach);
            };
            exposeSnapshot(true);
            try { this->save(stepsToWriteOnDisk); }
            catch (...) { exposeSnapshot(false); throw; }
            exposeSnapshot(false);
            if (materialize) markStaticFieldsWritten(false);
            M_hasSaved=true;
        }
        for (auto& [ts,steps] : stepsToWriteOnDisk)
        {
            for (auto const& step : steps)
            {
                step->setState(STEP_ON_DISK);
                step->cleanup();
            }
            ts->save((fs::path(this->path())/(ts->name()+".timeset")).string(),this->worldComm());
        }
        if (hasStaticFields() && M_staticWritten) writeStaticFieldGuard();
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
        if (hasStaticFields() || fs::exists(fs::path(this->path())/(this->prefix()+".static-fields")))
            throw std::logic_error("restart with dataset fields requires a new exporter/output directory");
        auto __ts_it = this->beginTimeSet();
        auto __ts_en = this->endTimeSet();

        for ( ; __ts_it != __ts_en ; ++__ts_it )
        {
            auto filename = this->path()+"/"+(*__ts_it)->name()+".timeset";
            if ( !fs::exists( filename ) )
                return;
            ( *__ts_it )->load( filename,__time );
            M_hasSaved=true;
        }
    }

    ExporterGeometry exporterGeometry() const { return M_ex_geometry; }
    //@}
protected:

    /** @return Whether the only temporal record is a stationary-format adapter. */
    bool isStationaryDataset() const
    { return M_staticOnlyAdapter && M_ts_set.size()==1 && M_ts_set.front()->numberOfSteps()==1; }

    /**
     * this p save function is defined by the Exporter subclasses and implement
     * saving the data to files
     */
    virtual void save( steps_write_on_disk_type const& stepsToWriteOnDisk ) const = 0;


    /** @brief Mark successful output; only native storage may release copied values. */
    void markStaticFieldsWritten(bool releasePayload=true) const
    {
        M_staticOutputPath=fs::absolute(this->path()).string();
        M_staticWritten=true;
        if (releasePayload) M_staticFields->cleanup();
    }

private:
    /** @brief Return one bounded error on every rank, without all-gathering P strings.
     * Success reduces one rank identifier. Failure broadcasts at most 4 KiB
     * from the first failing rank, preserving collective exception behavior.
     */
    std::string collectiveError(std::string const& error) const
    {
        auto const& comm=this->worldComm().comm();
        int source=error.empty()?comm.size():comm.rank();
        int first=mpi::all_reduce(comm,source,mpi::minimum<int>());
        if (first==comm.size()) return {};
        std::string message=comm.rank()==first?error.substr(0,4096):std::string{};
        mpi::broadcast(comm,message,first);
        return message;
    }

    /** @brief Propagate validation errors before collective conversion or I/O. */
    void validateDatasetError(std::string const& error) const
    {
        auto message=collectiveError(error);
        if (!message.empty()) throw std::invalid_argument("Exporter dataset fields: "+message);
    }

    /** @brief Validate mesh, communicator and support, including mixed-space subfields. */
    template<typename F>
    std::string validateDatasetFunction(F const& field) const
    {
        if constexpr (F::functionspace_type::nSpaces > 1)
        {
            std::string error;
            hana::for_each(hana::make_range(hana::int_c<0>,hana::int_c<F::functionspace_type::nSpaces>),
                          [&](auto i) {
                              auto e=validateDatasetFunction(field.template element<decltype(i)::value>());
                              if (!e.empty()) error=e;
                          });
            return error;
        }
        else
        {
            if (!field.worldComm().isActive()) return "field communicator is inactive";
            int relation=MPI_UNEQUAL;
            MPI_Comm_compare(this->worldComm().comm(),field.worldComm().comm(),&relation);
            if (!field.worldComm().isActive() || (relation!=MPI_IDENT && relation!=MPI_CONGRUENT))
                return "field communicator differs from exporter communicator";
            if (!M_mesh || !M_mesh->isSameMesh(field.functionSpace()->mesh()))
                return "field must use the exporter mesh";
            if (field.functionSpace()->dof()->hasMeshSupport() && field.functionSpace()->dof()->meshSupport()->isPartialSupport())
                return "partial-support FE snapshots are unsupported; use an explicit expression range";
            return "";
        }
    }

    /** @brief Validate a collective schema and populate an exporter-owned FieldSet. */
    template<typename Populate>
    void registerDatasetField(std::string const& name,typename field_set_type::variant_representation_arg_type const& reps,
                              std::string const& shape,std::string error,Populate&& populate)
    {
        std::set<std::string> representations;
        try { representations=field_set_type::representationType(reps); }
        catch (std::invalid_argument const& e) { error=e.what(); }
        std::string key=sanitize(name);
        std::ostringstream schema;
        schema << key << ':' << shape;
        for (auto const& rep : representations) schema << ':' << rep;
        auto expectedSchema=schema.str();
        mpi::broadcast(this->worldComm().comm(),expectedSchema,this->worldComm().masterRank());
        if (expectedSchema!=schema.str()) error="dataset field schema differs between MPI ranks";
        if (!M_mesh) error="set the exporter mesh before registering dataset fields";
        else if (M_ts_set.size()!=1) error="dataset fields currently require one output mesh/time sequence";
        else if (M_hasSaved || M_ts_set.front()->numberOfSteps()) error="register dataset fields before creating steps or saving";
        else if (key.empty() || key.find_first_of("/\\")!=std::string::npos) error="invalid dataset field name";
        else if (M_staticNames->count(key) || M_staticNames->count(key+"_n") || M_staticNames->count(key+"_e"))
            error="duplicate dataset field name: "+key;
        validateDatasetError(error);
        // Build transactionally: a rejected mixed-field subname must not leave
        // earlier subfields in the dataset or overwrite an existing snapshot.
        if (M_staticFields) validateStaticFields();
        auto candidate=std::make_shared<field_set_type>(M_staticFields ? M_staticFields->spaceCache() : nullptr);
        candidate->setMesh(M_mesh);
        candidate->reserveNames(M_staticNames);
        populate(*candidate);
        if (!M_staticFields)
        {
            M_staticFields=candidate;
            M_staticMesh=M_mesh;
            M_staticStructureRevision=M_mesh->functionSpaceStructuralRevision();
            M_staticGeometryRevision=M_mesh->geometryRevision();
            M_staticSequenceName=M_ts_set.front()->name();
        }
        else M_staticFields->materialize(*candidate,true);
        M_staticNames->insert(key);
        auto names=M_staticFields->names();
        M_staticNames->insert(names.begin(),names.end());
    }

    /** @brief Reject stale geometry/layout and name collisions before opening files. */
    void validateStaticFields() const
    {
        std::string error;
        if (M_ts_set.size()!=1) error="dataset fields currently require one output mesh/time sequence";
        if (M_mesh!=M_staticMesh ||
            M_staticMesh->functionSpaceStructuralRevision()!=M_staticStructureRevision ||
            M_staticMesh->geometryRevision()!=M_staticGeometryRevision)
            error="mesh changed after dataset registration; start a new exporter";
        for (auto const& ts : M_ts_set)
        {
            if (!ts->hasMesh() || ts->mesh()!=M_staticMesh) error="time sequence uses a different dataset mesh";
            if (ts->name()!=M_staticSequenceName) error="cannot rename the output sequence after dataset registration";
            for (auto it=ts->beginStep();it!=ts->endStep();++it)
            {
                auto const& step=*it;
                if (step->isIgnored()) continue;
                if (!step->hasMesh() || step->mesh()!=M_staticMesh) error="step uses a different dataset mesh";
                for (auto const& name : step->fieldSet()->names())
                    if (M_staticNames->count(name)) error="temporal field collides with dataset field: "+name;
            }
        }
        validateDatasetError(error);
    }

    /** @brief Reject moving already-materialized snapshots to another destination. */
    void validateStaticOutputPath(std::string const& path) const
    {
        if (M_staticWritten && M_staticOutputPath!=fs::absolute(path).string())
            throw std::logic_error("cannot change the output path after writing dataset fields");
    }

    /** @brief Publish an exporter-level restart guard, not a temporal field manifest. */
    void writeStaticFieldGuard() const
    {
        std::string error;
        if (this->worldComm().isMasterRank())
        {
            std::ofstream out((fs::path(this->path())/(this->prefix()+".static-fields")).string());
            out << "Feel++ dataset field snapshots; restart unsupported\n";
            for (auto const& name : *M_staticNames) out << name << '\n';
            out.close();
            if (!out) error="cannot write dataset restart guard";
        }
        validateDatasetError(error);
    }

    //! Exporter-owned dataset geometry and time-neutral field storage.
    mesh_ptrtype M_mesh, M_staticMesh;
    field_set_ptrtype M_staticFields;
    //! Dataset names; temporal fields only borrow weak reservations.
    std::shared_ptr<std::set<std::string>> M_staticNames=std::make_shared<std::set<std::string>>();

    //! Mesh-layout revision binding the snapshot to exported ordering.
    uint64_type M_staticStructureRevision=0, M_staticGeometryRevision=0;
    //! Materialization/registration lifecycle, independent of temporal steps.
    mutable bool M_staticWritten=false, M_hasSaved=false, M_staticOnlyAdapter=false;
    //! Fixed destination after the first successful snapshot write.
    mutable std::string M_staticOutputPath;
    //! Filename-layout guard only; static field ownership remains time-neutral.
    std::string M_staticSequenceName;

protected:

    bool M_do_export;
    MeshFragmentation<mesh_type> M_meshFragmentation;
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
    using mesh_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(mesh)>>>;
    using mesh_fragmentation_type = MeshFragmentation<mesh_type>;
    auto && meshFragmentation = args.get_else_invocable(_byparts,[](){ return mesh_fragmentation_type(boption(_name="exporter.byparts")? mesh_fragmentation_type::Strategy::AllMarkedElements : mesh_fragmentation_type::Strategy::None ); } );

    std::string const& name = args.get_else_invocable(_name,[](){ return Environment::about().appName(); } );
    std::string const& geo = args.get_else_invocable(_geo, [](){ return soption(_name="exporter.geometry"); } );
    auto && path = args.get_else_invocable(_path, [&name](){ return std::string((fs::path(Environment::exportsRepository())/fs::path(soption("exporter.format"))/name).string()); } );

    using exporter_type = Exporter<mesh_type,mesh_type::nOrder>;

    auto e =  exporter_type::New( name,mesh->worldCommPtr() );
    e->setPrefix( name );
    e->setUseSingleTransientFile( fileset );
    e->setMeshFragmentation( meshFragmentation );
    if ( std::string(geo).compare("change_coords_only") == 0 )
        e->setMesh( mesh, EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY );
    else if ( std::string(geo).compare("change") == 0 )
        e->setMesh( mesh, EXPORTER_GEOMETRY_CHANGE );
    else if ( geo == "static" )
        e->setMesh( mesh, EXPORTER_GEOMETRY_STATIC );
    else
        throw std::invalid_argument( "unknown exporter geometry: " + geo );
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
    typedef Feel::Exporter<decay_type<MeshType>,N> type;
    typedef std::shared_ptr<type> ptrtype;
};

}

template<typename MeshType, int N = 1>
using exporter_t = typename meta::Exporter<MeshType,N>::type;

template<typename MeshType, int N = 1>
using exporter_ptr_t = typename meta::Exporter<MeshType,N>::ptrtype;

} // Feel

//#if !defined( FEELPP_INSTANTIATION_MODE )
# include <feel/feelfilters/exporterimpl.hpp>
//#endif // FEELPP_INSTANTIATION_MODE


#endif /* __Exporter_H */
