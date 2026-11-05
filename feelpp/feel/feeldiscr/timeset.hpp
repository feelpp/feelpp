/*

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-12

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2012 Universite Joseph Fourier

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
//! \file timeSet.hpp
//! \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! \date 2005-01-12
#ifndef FEELPP_DISCR_TIMESET_HPP
#define FEELPP_DISCR_TIMESET_HPP 1

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <variant>

#include <boost/operators.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/serialization/map.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>

// #include <boost/numeric/ublas/vector_serialize.hpp>

#include <feel/feelcore/context.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelcomplex.hpp>

#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/lagrange.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/interpolate.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pdh.hpp>

#include <feel/feeldiscr/elementdiv.hpp>
#include <feel/feeldiscr/exportfieldset.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelvf/vf.hpp>

#define TS_INITIAL_INDEX 1

namespace Feel
{

//! \brief Legacy serialized step-state flags, retained for API/archive compatibility.
//! Scheduling and publication state belong to TimeSet::Step, not ExportFieldSet.
enum
{
    STEP_NEW = ( 1 << 0 ),       //!< Sample has not been explicitly marked non-new.
    STEP_HAS_DATA = ( 1 << 1 ),  //!< Geometry or field metadata is present.
    STEP_ON_DISK = ( 1 << 2 ),   //!< Current field revision has been published.
    STEP_IN_MEMORY = ( 1 << 3 ), //!< Field buffers are marked resident.
    STEP_IGNORED = ( 1 << 4 ),   //!< Output frequency suppresses this sample.
    STEP_OVERWRITE = ( 1 << 10 ) //!< Legacy overwrite policy.
};

template <typename A0, typename A1, typename A2, typename A3, typename A4> class FunctionSpace;
template <typename MeshType, int N> class Exporter;
namespace detail
{
//! \class TimeSet
//! \ingroup SpaceTime
//! \brief data TimeSet
//!
//! \tparam MeshType     Mesh type
//! \tparam N            Mesh geometrical order
//!
//! \author Christophe Prud'homme
template <typename MeshType, int N = 1> class TimeSet
{
  public:
    //! \name Subclasses
    //@{
    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    using size_type = typename mesh_type::size_type;

    typedef Pdh_type<MeshType, 0> scalar_p0_space_type;
    typedef Pch_type<MeshType, N /*1*/> scalar_p1_space_type;
    typedef std::shared_ptr<scalar_p0_space_type> scalar_p0_space_ptrtype;
    typedef std::shared_ptr<scalar_p1_space_type> scalar_p1_space_ptrtype;

    typedef typename scalar_p0_space_type::element_type element_scalar_type;
    typedef typename scalar_p1_space_type::element_type nodal_scalar_type;

    //! \brief A temporal sample composed of a time-neutral field set.
    //! Time/index metadata belongs here; field conversion and payload storage
    //! belong to ExportFieldSet. The latter can also be owned directly by an
    //! exporter without constructing a Step or referencing a TimeSet.
    //! Step owns scheduling and publication policy. Field revisions invalidate the
    //! published state even when registration uses fieldSet() directly.
    //! Mutation is not thread-safe; owners must honor collective registration.
    class Step : public boost::equality_comparable<Step>, public boost::less_than_comparable<Step>
    {
      public:
        using field_set_type = ExportFieldSet<MeshType, N>;
        using field_set_ptrtype = std::shared_ptr<field_set_type>;
        using step_type = Step;
        using mesh_type = MeshType;
        using mesh_ptrtype = std::shared_ptr<mesh_type>;
        using size_type = typename mesh_type::size_type;
        using scalar_p0_space_type = typename field_set_type::scalar_p0_space_type;
        using scalar_p1_space_type = typename field_set_type::scalar_p1_space_type;
        using scalar_p0_space_ptrtype = typename field_set_type::scalar_p0_space_ptrtype;
        using scalar_p1_space_ptrtype = typename field_set_type::scalar_p1_space_ptrtype;
        using element_scalar_type = typename field_set_type::element_scalar_type;
        using nodal_scalar_type = typename field_set_type::nodal_scalar_type;
        using nodal_scalar_ptrtype = typename field_set_type::nodal_scalar_ptrtype;
        using element_scalar_ptrtype = typename field_set_type::element_scalar_ptrtype;
        using nodal_field_type = typename field_set_type::nodal_field_type;
        using element_field_type = typename field_set_type::element_field_type;
        using map_scalar_type = typename field_set_type::map_scalar_type;
        using map_complex_type = typename field_set_type::map_complex_type;
        using map_nodal_type = typename field_set_type::map_nodal_type;
        using map_element_type = typename field_set_type::map_element_type;
        using scalar_iterator = typename field_set_type::scalar_iterator;
        using scalar_const_iterator = typename field_set_type::scalar_const_iterator;
        using complex_iterator = typename field_set_type::complex_iterator;
        using complex_const_iterator = typename field_set_type::complex_const_iterator;
        using nodal_iterator = typename field_set_type::nodal_iterator;
        using nodal_const_iterator = typename field_set_type::nodal_const_iterator;
        using element_iterator = typename field_set_type::element_iterator;
        using element_const_iterator = typename field_set_type::element_const_iterator;
        using variant_representation_arg_type =
            typename field_set_type::variant_representation_arg_type;

        //! \brief Release this sample's field-set ownership.
        ~Step() = default;

        //! \return Field data composed into this temporal sample.
        field_set_ptrtype const &fieldSet() const { return M_fields; }

        //! \return Physical time of this sample.
        Real time() const { return M_time; }

        //! \return Temporal sequence index.
        size_type index() const { return M_index; }

        //! \return Index among non-ignored output samples.
        size_type activeIndex() const { return M_activeIndex; }

        //! \return Whether the legacy new-sample flag is set.
        bool isNew() const { return M_state.test( STEP_NEW ); }

        //! \return Whether this sample contains geometry or values.
        bool hasData() const { return M_fields->hasData(); }

        //! \return Whether this sample has an on-disk representation.
        bool isOnDisk() const
        {
            return M_state.test( STEP_ON_DISK ) && M_writtenRevision == M_fields->revision();
        }

        //! \return Whether component buffers are resident.
        bool isInMemory() const { return M_fields->isInMemory(); }

        //! \return Whether output frequency suppresses this sample.
        bool isIgnored() const { return M_state.test( STEP_IGNORED ); }

        //! \return Legacy state bits assembled from temporal policy and field residency.
        size_type state() const
        {
            return ( M_state.context() & ~( STEP_HAS_DATA | STEP_IN_MEMORY | STEP_ON_DISK ) ) |
                   ( hasData() ? STEP_HAS_DATA : 0 ) | ( isInMemory() ? STEP_IN_MEMORY : 0 ) |
                   ( isOnDisk() ? STEP_ON_DISK : 0 );
        }

        //! \brief Add legacy state flags, recording the published revision when requested.
        void setState( size_type state )
        {
            M_state.set( state & ~( STEP_HAS_DATA | STEP_IN_MEMORY ) );
            M_fields->restoreDataState( hasData() || ( state & STEP_HAS_DATA ),
                                        isInMemory() || ( state & STEP_IN_MEMORY ) );
            if ( state & STEP_ON_DISK )
            {
                M_writtenRevision = M_fields->revision();
            }
        }

        //! \return Whether geometry has been associated.
        bool hasMesh() const { return M_fields->hasMesh(); }

        //! \return Mesh associated with this sample.
        mesh_ptrtype mesh() const { return M_fields->mesh(); }

        //! \brief Associate geometry with this sample's field storage.
        void setMesh( mesh_ptrtype const &mesh ) { M_fields->setMesh( mesh ); }

        //! \brief Add a scalar, FE field or expression to this temporal sample.
        //! All existing representation/range overloads are forwarded unchanged.
        template <typename... Args> void add( Args &&...args )
        {
            if ( isIgnored() )
            {
                return;
            }
            M_fields->add( std::forward<Args>( args )... );
        }

        //! \brief Preserve initializer-list names for mixed FE fields.
        template <typename F> void add( std::initializer_list<std::string> names, F const &field )
        {
            if ( isIgnored() )
            {
                return;
            }
            M_fields->add( names, field );
        }

        //! \brief Compatibility wrapper for per-case scalar registration.
        FEELPP_DEPRECATED void addScalar( std::string const &name, scalar_type value,
                                          bool cst = false )
        {
            this->add( name, value, cst );
        }

        //! \brief Register a complex quantity using the existing backend convention.
        void addComplex( std::string const &name, complex_type const &value, bool cst = false )
        {
            if ( isIgnored() )
            {
                return;
            }
            M_fields->addComplex( name, value, cst );
        }

        //! \brief Add partition IDs to this temporal sample.
        void addRegions( std::string const &prefix = "" ) { this->addRegions( prefix, prefix ); }

        //! \brief Add partition IDs with a distinct output filename prefix.
        void addRegions( std::string const &prefix, std::string const &filename )
        {
            if ( isIgnored() )
            {
                return;
            }
            M_fields->addRegions( prefix, filename );
        }

        //! \return First per-case scalar.
        scalar_const_iterator beginScalar() const { return M_fields->beginScalar(); }

        //! \return End of per-case scalars.
        scalar_const_iterator endScalar() const { return M_fields->endScalar(); }

        //! \return Named per-case scalar.
        scalar_type scalar( std::string const &name ) const { return M_fields->scalar( name ); }

        //! \return First nodal field.
        nodal_const_iterator beginNodal() const { return M_fields->beginNodal(); }

        //! \return End of nodal fields.
        nodal_const_iterator endNodal() const { return M_fields->endNodal(); }

        //! \return Named nodal field.
        nodal_field_type const &nodal( std::string const &name ) const
        {
            return M_fields->nodal( name );
        }

        //! \return Range of nodal fields.
        auto nodal() const { return M_fields->nodal(); }

        //! \return First element field.
        element_const_iterator beginElement() const { return M_fields->beginElement(); }

        //! \return End of element fields.
        element_const_iterator endElement() const { return M_fields->endElement(); }

        //! \return Named element field.
        element_field_type const &element( std::string const &name ) const
        {
            return M_fields->element( name );
        }

        //! \return Persistent display name after buffer cleanup.
        std::string const &fieldName( std::string const &key, bool nodal ) const
        {
            return M_fields->fieldName( key, nodal );
        }

        //! \brief Release values while retaining the schema.
        void cleanup() { M_fields->cleanup(); }

        //! \brief Restore the legacy in-memory flag; this does not read field values.
        void load() { M_fields->restoreDataState( hasData(), true ); }

        //! \brief Emit temporal and field-state diagnostics at verbose level.
        void showMe( std::string const &text ) const
        {
            DVLOG( 2 ) << text << " index: " << M_index << " time: " << M_time
                       << " isNew: " << isNew() << " isIgnored: " << isIgnored()
                       << " isOnDisk: " << isOnDisk();
            M_fields->showMe( text );
        }

        //! \return Equality by sample index.
        bool operator==( Step const &other ) const { return index() == other.index(); }

        //! \return Ordering by sample index.
        bool operator<( Step const &other ) const { return index() < other.index(); }

      private:
        friend class TimeSet;
        //! \brief Construct an empty sample for temporal bookkeeping.
        Step() = default;

        //! \brief Construct a temporal sample sharing its sequence's conversion workspace.
        Step( TimeSet *ts, Real time, size_type index, size_type activeIndex,
              size_type state = STEP_NEW | STEP_OVERWRITE )
            : M_time( time ), M_index( index ), M_activeIndex( activeIndex ), M_state( 0 ),
              M_fields( std::make_shared<field_set_type>( ts->M_fieldSpaceCache ) )
        {
            setState( state );
        }

        //! Temporal metadata exists only in Step, never in ExportFieldSet.
        Real M_time = 0;
        size_type M_index = 0;                        //!< Sequence index.
        size_type M_activeIndex = 0;                  //!< Index excluding ignored samples.
        Context M_state{ STEP_NEW | STEP_OVERWRITE }; //!< Temporal/output policy flags.
        std::uint64_t M_writtenRevision = 0;          //!< Field revision recorded at publication.
        //! Independently owned payload/schema with a shared conversion workspace.
        field_set_ptrtype M_fields = std::make_shared<field_set_type>();
    };

    //! \brief Order sample pointers by sequence index.
    struct ltstep
    {
        //! \return Whether the first sample precedes the second.
        bool operator()( std::shared_ptr<Step> const &s1, std::shared_ptr<Step> const &s2 ) const
        {
            return *s1 < *s2;
        }
    };

    //! \name Typedefs
    //@{

    typedef Step step_type;
    typedef std::shared_ptr<Step> step_ptrtype;
    typedef std::set<step_ptrtype, ltstep> step_set_type;
    typedef typename step_set_type::iterator step_iterator;
    typedef typename step_set_type::const_iterator step_const_iterator;
    typedef typename step_set_type::reverse_iterator step_reverse_iterator;
    typedef typename step_set_type::const_reverse_iterator step_const_reverse_iterator;
    //@}

    //@}

    //! \name Constructors, destructor
    //@{

    //! constructor for a time set
    //!
    //! \param filename name of the file that stores the timeset information
    //! \param init if true, remove the file that stores the timset info if it exists
    TimeSet( std::string filename = "undefined", bool init = false );

    //! \brief Copy time-set metadata without copying samples.
    TimeSet( TimeSet const & );

    //! \brief Release the time-set index.
    ~TimeSet();

    //@}

    //! \name Operator overloads
    //@{

    //! \brief Assign time-set metadata.
    TimeSet &operator=( TimeSet const & );

    //@}

    //! \name Accessors
    //@{

    //! \return the name of the time set
    std::string const &name() const { return M_name; }

    //! \brief Update the display name of the time set.
    void setName( std::string const &name ) const { M_name = name; }

    //! \return the index of the time set
    uint32_type index() const { return M_index; }

    //! \return the number of steps already stored
    size_type numberOfSteps() const { return M_step_set.size(); }

    //! \return the number of active steps
    size_type numberOfActiveSteps() const
    {
        size_type res = 0;
        auto it = this->beginStep();
        auto end = this->endStep();
        for ( ; it != end; ++it )
        {
            step_ptrtype sample = *it;
            if ( !sample->isIgnored() )
                ++res;
        }
        return res;
    }

    //! return the first active step
    step_ptrtype firstActiveStep() const
    {
        auto it = this->beginStep();
        auto end = this->endStep();
        for ( ; it != end; ++it )
        {
            step_ptrtype sample = *it;
            if ( !sample->isIgnored() )
                return sample;
        }
        step_ptrtype sample;
        return sample;
    }

    //! \return the time increment between two steps
    Real timeIncrement() const { return M_time_increment; }

    //! \return Resident, non-ignored samples whose current data is not published.
    step_set_type stepsToWriteOnDisk() const
    {
        auto it = this->beginStep();
        auto end = this->endStep();
        step_set_type stepToWriteOnDisk;
        for ( ; it != end; ++it )
        {
            step_ptrtype sample = *it;
            if ( sample->isOnDisk() || sample->isIgnored() )
                continue;
            if ( sample->hasData() && sample->isInMemory() )
                stepToWriteOnDisk.insert( sample );
        }
        return stepToWriteOnDisk;
    }

    //! \name  Mutators
    //@{

    //! set the time increment
    void setTimeIncrement( Real increment ) { M_time_increment = increment; }

    //! \brief Set the legacy in-memory step-retention count.
    void setNumberOfStepsInMemory( uint16_type i ) { M_keep_steps = i; }

    //@}

    //! \name  Methods
    //@{

    //! \param time time at which we want to get the step
    //! \param freq active step by frequence
    //! \return a step defined at time \c time if not found then generate a new one
    step_ptrtype step( Real time, int freq = 1 );

    //! \return Iterator to the first sample.
    step_iterator beginStep() { return M_step_set.begin(); }

    //! \return Const iterator to the first sample.
    step_const_iterator beginStep() const { return M_step_set.begin(); }

    //! \return Reverse iterator to the last sample.
    step_reverse_iterator rbeginStep() { return M_step_set.rbegin(); }

    //! \return Const reverse iterator to the last sample.
    step_const_reverse_iterator rbeginStep() const { return M_step_set.rbegin(); }

    //! \return End iterator for samples.
    step_iterator endStep() { return M_step_set.end(); }

    //! \return Const end iterator for samples.
    step_const_iterator endStep() const { return M_step_set.end(); }

    //! \return Reverse end iterator for samples.
    step_reverse_iterator rendStep() { return M_step_set.rend(); }

    //! \return Const reverse end iterator for samples.
    step_const_reverse_iterator rendStep() const { return M_step_set.rend(); }

    //! \return Insertion position and whether the sample index was new.
    std::pair<step_iterator, bool> insertStep( step_ptrtype sample )
    {
        return M_step_set.insert( sample );
    }

    //! \brief Release all samples, retaining sequence metadata and conversion spaces.
    void clear() { M_step_set.clear(); }

    //! \brief Restore temporal metadata and discard samples later than time.
    //! Field payloads are not deserialized; static-field restart remains unsupported.
    void load( std::string const &filename, Real time )
    {
        if ( fs::exists( filename + ".static-fields" ) )
            throw std::logic_error( "static field restart is not supported" );
        std::ifstream ifs( filename );
        // load data from archive
        boost::archive::text_iarchive ia( ifs );
        ia >> *this;

        resetPreviousTime( time );
    }

    //! \brief Save temporal metadata on the communicator's master rank.
    void save( std::string const &filename, WorldComm const &worldComm )
    {
        if ( worldComm.isMasterRank() )
        {
            std::ofstream ofs( filename );
            // save data from archive
            boost::archive::text_oarchive oa( ofs );
            oa << *this;
        }
    }

    //@}

  protected:
    //! name of the time set
    mutable std::string M_name;

    //! index of the time set
    uint32_type M_index;

    //! steps
    step_set_type M_step_set;

    //! time increment
    Real M_time_increment;

    uint16_type M_keep_steps;

  private:
    friend class boost::serialization::access;

    //! \brief Preserve the legacy metadata archive layout and step-state bitset.
    template <class Archive>
    FEELPP_NO_EXPORT void serialize( Archive &ar, const unsigned int /*version*/ )
    {
        ar &boost::serialization::make_nvp( "name", M_name );
        ar &boost::serialization::make_nvp( "index", M_index );
        ar &boost::serialization::make_nvp( "time_increment", M_time_increment );
        ar &boost::serialization::make_nvp( "keep_steps", M_keep_steps );

        if ( Archive::is_saving::value )
        {
            size_type s = M_step_set.size();
            ar &boost::serialization::make_nvp( "number_of_steps", s );

            step_iterator it = M_step_set.begin();
            step_iterator end = M_step_set.end();

            for ( ; it != end; ++it )
            {
                double t = ( *it )->time();
                ar &boost::serialization::make_nvp( "time", t );
                size_type ind = ( *it )->index();
                ar &boost::serialization::make_nvp( "index", ind );
                size_type state = ( *it )->state();
                ar &boost::serialization::make_nvp( "state", state );
            }
        }

        if ( Archive::is_loading::value )
        {
            size_type s( 0 );
            size_type nActiveIndex = 0;
            size_type activeIndex = invalid_v<size_type>;
            ar &boost::serialization::make_nvp( "number_of_steps", s );

            for ( size_type i = 0; i < s; ++i )
            {
                double t = 0;
                ar &boost::serialization::make_nvp( "time", t );

                size_type ind = 0;
                ar &boost::serialization::make_nvp( "index", ind );

                size_type state = 0;
                ar &boost::serialization::make_nvp( "state", state );
                Context ctxState( state );
                if ( !ctxState.test( STEP_IGNORED ) )
                {
                    ++nActiveIndex;
                    activeIndex = nActiveIndex;
                }
                else
                    activeIndex = invalid_v<size_type>;

                step_iterator stepIt;
                bool inserted;
                boost::tie( stepIt, inserted ) = M_step_set.insert(
                    step_ptrtype( new Step( this, t, ind, activeIndex, state ) ) );

                CHECK( inserted ) << "insertion failed at t=" << t << " and ind=" << ind;
            }
        }
    }

    //! delete all time next this time (after a restart by example)
    FEELPP_NO_EXPORT void resetPreviousTime( Real time );

  public:
    //! \brief Associate geometry with this temporal sequence.
    void setMesh( mesh_ptrtype m ) { M_mesh = m; }

    //! \return Whether geometry has been associated with this sequence.
    bool hasMesh() const { return M_mesh != std::nullopt; }

    //! \return Geometry associated with this sequence.
    mesh_ptrtype mesh() const
    {
        DLOG_IF( WARNING, hasMesh() ) << "Time Set has no mesh data structure associated\n";
        return M_mesh.value();
    }

  public:
    std::optional<mesh_ptrtype> M_mesh;

    //! Conversion workspace shared by this sequence's independently owned field sets.
    std::shared_ptr<typename Step::field_set_type::SpaceCache> M_fieldSpaceCache =
        std::make_shared<typename Step::field_set_type::SpaceCache>();

  private:
    //! keep track of the current index of the TimeSets to ensure they get
    static uint32_type S_currentIndex;
};

template <typename MeshType, int N>
inline FEELPP_EXPORT bool operator<( TimeSet<MeshType, N> const &ts1,
                                     TimeSet<MeshType, N> const &ts2 )
{
    return ts1.index() < ts2.index();
}

template <typename MeshType, int N>
inline FEELPP_EXPORT bool operator<( std::shared_ptr<TimeSet<MeshType, N>> const &ts1,
                                     std::shared_ptr<TimeSet<MeshType, N>> const &ts2 )
{
    return ts1->index() < ts2->index();
}

template <typename MeshType, int N>
FEELPP_NO_EXPORT bool operator<( typename TimeSet<MeshType, N>::step_type const &s1,
                                 typename TimeSet<MeshType, N>::step_type const &s2 )
{
    // return s1->index() < s2->index();
    return s1->time() < s2->time();
}

template <typename MeshType, int N>
FEELPP_NO_EXPORT bool operator<( typename TimeSet<MeshType, N>::step_ptrtype const &s1,
                                 typename TimeSet<MeshType, N>::step_ptrtype const &s2 )
{
    return s1->index() < s2->index();
}

template <typename MeshType, int N>
TimeSet<MeshType, N>::TimeSet( std::string name, bool init )
    : M_name( name ), M_index( S_currentIndex++ ), M_time_increment( 0.1 ), M_keep_steps( 1 )
{
}

template <typename MeshType, int N>
TimeSet<MeshType, N>::TimeSet( TimeSet const &ts )
    : M_name( ts.M_name ), M_index( ts.M_index ), M_time_increment( ts.M_time_increment ),
      M_keep_steps( ts.M_keep_steps )
{
}

template <typename MeshType, int N> TimeSet<MeshType, N>::~TimeSet()
{
    --S_currentIndex;
}

template <typename MeshType, int N> void TimeSet<MeshType, N>::resetPreviousTime( Real time )
{
    step_iterator it = M_step_set.begin();
    step_iterator end = M_step_set.end();
    bool find = false;

    while ( !find && it != end )
    {
        double t = ( *it )->time();
        double eps = 1e-10;
        if ( ( t - eps ) <= time )
            ++it;
        else
            find = true;
    }

    if ( find )
        M_step_set.erase( it, end );
}

template <typename MeshType, int N>
typename TimeSet<MeshType, N>::step_ptrtype TimeSet<MeshType, N>::step( Real time, int freq )
{
    size_type nStepBefore = 0, nActiveStepBefore = 0;
    step_iterator stepIt = beginStep();
    for ( ; stepIt != endStep(); ++stepIt )
    {
        if ( math::abs( ( *stepIt )->time() - time ) < 1e-10 )
            break;
        ++nStepBefore;
        if ( !( *stepIt )->isIgnored() )
            ++nActiveStepBefore;
    }

    bool ignoreStep = ( nStepBefore % freq ) > 0;

    if ( stepIt == endStep() )
    {
        bool inserted;
        DVLOG( 2 ) << "[TimeSet<MeshType, N>::step] Inserting new step at time " << time
                   << " with index " << numberOfSteps();

        size_type theState =
            ignoreStep ? STEP_NEW | STEP_OVERWRITE | STEP_IGNORED : STEP_NEW | STEP_OVERWRITE;
        size_type activeIndex = ignoreStep ? invalid_v<size_type> : nActiveStepBefore + 1;
        step_ptrtype thestep(
            new Step( this, time, nStepBefore /*numberOfSteps()*/ + 1, activeIndex, theState ) );
        if ( this->hasMesh() && !ignoreStep )
            thestep->setMesh( this->mesh() );
        boost::tie( stepIt, inserted ) = insertStep( thestep );

        DVLOG( 2 ) << "[TimeSet<MeshType, N>::step] step was inserted properly? "
                   << ( inserted ? "yes" : "no" ) << "\n";
        DVLOG( 2 ) << "[TimeSet<MeshType, N>::step] step index : " << ( *stepIt )->index()
                   << " time : " << ( *stepIt )->time() << "\n";
        // cleanup();
    }

    else
    {
        DVLOG( 2 ) << "[TimeSet<MeshType, N>::step] found step at time " << time << " with index "
                   << ( *stepIt )->index() << "\n";
        ( *stepIt )->showMe( "TimesSet::step(t)" );
    }

    return *stepIt;
}

template <typename MeshType, int N>
uint32_type TimeSet<MeshType, N>::S_currentIndex = TS_INITIAL_INDEX;

} // namespace detail
} // namespace Feel
#endif // FEELPP_DISCR_TIMESET_HPP
