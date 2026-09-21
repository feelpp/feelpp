/*

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-11-09

  Copyright (C) 2004,2005 EPFL
  Copyright (C) 2007-2012 Universite Joseph Fourier (Grenoble I)

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
//! \file ExporterEnsightGold.hpp
//! \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
//! \date 2006-11-26
#ifndef FEELPP_FILTERS_EXPORTERENSIGHTGOLD_HPP
#define FEELPP_FILTERS_EXPORTERENSIGHTGOLD_HPP 1

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#ifdef __linux__
#include <linux/magic.h>
#include <sys/vfs.h>
#endif

#include <feel/feelfilters/detail/fileindex.hpp>
#include <feel/feelfilters/detail/meshcontiguousnumberingmapping.hpp>
#include <feel/feelmesh/filters.hpp>

namespace Feel
{

//! \class ExporterEnsightGold
//! \brief exporter to EnsightGold format
//!
//! \ingroup Exporter
//! @author Christophe Prud'homme
//! @author Alexandre Ancel
template <typename MeshType, int N> class ExporterEnsightGold : public Exporter<MeshType, N>
{
    typedef Exporter<MeshType, N> super;

  public:
    //! @name Typedefs
    //@{

    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    using index_type = typename mesh_type::index_type;
    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_iterator timeset_iterator;
    typedef typename super::timeset_const_iterator timeset_const_iterator;
    //! \return False for packed/merged files, whose static-variable combination
    //! is not supported correctly by the verified ParaView legacy reader.
    bool supportsNativeStaticFields() const override { return !M_mergeTimeSteps; }

    //! \brief Select automatic, independent or bounded-memory root writes.
    //! All ranks must select the same policy before the first save. Automatic
    //! uses Root on Linux NFS; Independent requires a validated filesystem.
    void setIOPolicy( ExporterIOPolicy policy ) override
    {
        if ( !M_resolvedIOPath.empty() )
            throw std::logic_error( "EnSight I/O policy is frozen after the first save" );
        M_ioPolicy = policy;
    }

  protected:
    using mesh_contiguous_numbering_mapping_type =
        Feel::detail::MeshContiguousNumberingMapping<mesh_type, float>;
    using mesh_contiguous_numbering_mapping_ptrtype =
        std::shared_ptr<mesh_contiguous_numbering_mapping_type>;
    using step_ptrtype = typename super::step_ptrtype;
    using steps_write_on_disk_type = typename super::steps_write_on_disk_type;

  public:
    //@}

    //! @name Constructors, destructor
    //@{
    //! The elements that are supported by the EnSight6 format are:
    //!
    //! \htmlonly
    //! <pre>
    //! 1                 1------------------2        1----------2--------3
    //! point                   two node bar                three node bar
    //!
    //!
    //! 7
    //! 4-------------3          4-------------3
    //! 3                 |             |          |             |
    //! 3                        /\                |             |          |             |
    //! /\                      /  \               |             |        8 |             | 6
    //! /  \               6    /    \  5           |             |          |             |
    //! /    \                  /      \             |             |          |             |
    //! /      \                /        \            |             |          |             |
    //! /        \              /          \           |             |          |      5      |
    //! /          \            /    4       \          1-------------2          1-------------2
    //! 1------------2           1------------2
    //! three node triangle       six node triangle       four node quadrangle     eight node
    //! quadrangle
    //!
    //!
    //! /\
//! / |\
//! /  |4\
//! /   |  \
//! /    |   \
//! /     |    \
//! 1------|-----\
//! \     |    3/
    //! \    |    /
    //! \  2|   /
    //! \  |  /
    //! \ | /
    //! \\2/
    //!
    //! four node tetrahedron
    //! </pre>
    //! \endhtmlonly
    explicit ExporterEnsightGold( worldcomm_ptr_t const &worldComm = Environment::worldCommPtr() );

    ExporterEnsightGold( std::string const &__p = "default", int freq = 1,
                         worldcomm_ptr_t const &worldComm = Environment::worldCommPtr() );

    explicit ExporterEnsightGold( std::string const &exp_prefix,
                                  worldcomm_ptr_t const &worldComm = Environment::worldCommPtr() );

    //! \brief Copy configuration but start with fresh I/O and numbering caches.
    ExporterEnsightGold( ExporterEnsightGold const &ex ) : super( ex )
    {
        init();
        M_ioPolicy = ex.M_ioPolicy;
    }

    ~ExporterEnsightGold() override;

    //@}

    //! @name Operator overloads
    //@{

    //@}

    //! @name Accessors
    //@{

    //! \return the ensight element type
    std::string const &elementType() const { return M_element_type; }

    //@}

    //! @name  Mutators
    //@{

    //@}

    //! @name  Methods
    //@{

    void visit( mesh_type *mesh ) override;

    //@}

  protected:
    //! save the timeset
    void save( steps_write_on_disk_type const &stepsToWriteOnDisk ) const override;

  private:
    //! \brief Collectively resolve and cache the I/O policy for this destination.
    void resolveIOPolicy() const
    {
        if ( M_resolvedIOPath == this->path() )
            return;
        int requested = static_cast<int>( M_ioPolicy );
        int minimum = 0;
        int maximum = 0;
        checkMpiIo(
            MPI_Allreduce( &requested, &minimum, 1, MPI_INT, MPI_MIN, this->worldComm().comm() ),
            "minimum I/O policy" );
        checkMpiIo(
            MPI_Allreduce( &requested, &maximum, 1, MPI_INT, MPI_MAX, this->worldComm().comm() ),
            "maximum I/O policy" );
        if ( minimum != maximum || minimum < static_cast<int>( ExporterIOPolicy::Automatic ) ||
             maximum > static_cast<int>( ExporterIOPolicy::Root ) )
            throw std::invalid_argument(
                "EnSight I/O policy must be valid and identical on all ranks" );
        int useRoot = M_ioPolicy == ExporterIOPolicy::Root;
#ifdef __linux__
        if ( M_ioPolicy == ExporterIOPolicy::Automatic )
        {
            struct statfs info;
            if ( ::statfs( this->path().c_str(), &info ) != 0 )
                checkMpiIo( MPI_ERR_IO, "cannot identify output filesystem" );
            useRoot = info.f_type == NFS_SUPER_MAGIC;
        }
#endif
        int anyRoot = 0;
        checkMpiIo(
            MPI_Allreduce( &useRoot, &anyRoot, 1, MPI_INT, MPI_MAX, this->worldComm().comm() ),
            "resolve I/O policy" );
        M_useRootIO = anyRoot;
        if ( M_useRootIO && !M_payloadComm )
            M_payloadComm.emplace( this->worldComm().comm(), mpi::comm_duplicate );
        if ( M_useRootIO && this->worldComm().isMasterRank() )
            M_payloadBuffer.resize( 1024 * 1024 );
        M_resolvedIOPath = this->path();
    }

    //! \brief Collectively write a distributed payload, including empty ranks.
    //! Independent uses the original disjoint MPI writes. Root exchanges byte
    //! ranges on a private communicator and streams at most 1 MiB at a time,
    //! avoiding both cross-client NFS writes and a full-field gather allocation.
    //! \param offset Absolute byte offset of this rank's range.
    //! \param count Number of datatype entries in this rank's input buffer.
    void writePayload( MPI_File file, MPI_Offset offset, void const *data, std::size_t count,
                       MPI_Datatype datatype, MPI_Status *status ) const
    {
        int const width = validateWriteRange( offset, data, count, datatype );
        if ( !M_useRootIO )
        {
            writeAt( file, offset, data, count, datatype, status );
            return;
        }
        int rank = M_payloadComm->rank();
        int size = M_payloadComm->size();
        if ( count > std::size_t( std::numeric_limits<MPI_Offset>::max() / width ) )
            checkMpiIo( MPI_ERR_COUNT, "payload byte count overflow" );
        MPI_Offset descriptor[2] = { offset, MPI_Offset( count ) * width };
        int const root = this->worldComm().masterRank();

        std::vector<MPI_Offset> ranges( rank == root ? 2 * size : 0 );
        checkMpiIo( MPI_Gather( descriptor, 2, MPI_OFFSET, ranges.data(), 2, MPI_OFFSET, root,
                                *M_payloadComm ),
                    "gather payload ranges" );
        constexpr MPI_Offset chunkSize = 1024 * 1024;
        auto const *bytes = static_cast<char const *>( data );
        if ( rank != root )
        {
            for ( MPI_Offset sent = 0; sent < descriptor[1]; sent += chunkSize )
            {
                int n = int( std::min( chunkSize, descriptor[1] - sent ) );
                checkMpiIo( MPI_Send( bytes + sent, n, MPI_BYTE, root, 0, *M_payloadComm ),
                            "send payload" );
            }
        }
        else
        {
            for ( int peer = 0; peer < size; ++peer )
                for ( MPI_Offset received = 0; received < ranges[2 * peer + 1];
                      received += chunkSize )
                {
                    int n = int( std::min( chunkSize, ranges[2 * peer + 1] - received ) );
                    char const *buffer = nullptr;
                    if ( peer == root )
                        buffer = bytes + received;
                    else
                    {
                        checkMpiIo( MPI_Recv( M_payloadBuffer.data(), n, MPI_BYTE, peer, 0,
                                              *M_payloadComm, status ),
                                    "receive payload" );
                        int receivedCount = 0;
                        checkMpiIo( MPI_Get_count( status, MPI_BYTE, &receivedCount ),
                                    "payload receive count" );
                        if ( receivedCount != n )
                            checkMpiIo( MPI_ERR_IO, "short payload receive" );
                        buffer = M_payloadBuffer.data();
                    }
                    writeAt( file, ranges[2 * peer] + received, buffer, n, MPI_BYTE, status );
                }
        }
    }

    //! \brief Abort the exporter communicator on an MPI-I/O failure.
    //! Independent writes must not throw on one rank while peers enter a
    //! collective close. MPI_Abort is the fail-fast policy for unrecoverable I/O.
    //! \param error MPI return code.
    //! \param operation Operation/path included in the diagnostic.
    void checkMpiIo( int error, std::string const &operation ) const
    {
        if ( error == MPI_SUCCESS )
            return;
        char message[MPI_MAX_ERROR_STRING];
        int length = 0;
        MPI_Error_string( error, message, &length );
        std::cerr << "EnSight Gold " << operation << ": " << std::string( message, length )
                  << std::endl;
        MPI_Abort( this->worldComm().comm(), error );
        std::terminate();
    }

    //! \brief Validate contiguous basic-datatype byte ranges without overflowing offsets.
    //! \return Datatype width in bytes; invalid input aborts the exporter communicator.
    int validateWriteRange( MPI_Offset offset, void const *data, std::size_t count,
                            MPI_Datatype datatype ) const
    {
        int width = 0;
        MPI_Aint lower = 0;
        MPI_Aint extent = 0;
        checkMpiIo( MPI_Type_size( datatype, &width ), "datatype size" );
        checkMpiIo( MPI_Type_get_extent( datatype, &lower, &extent ), "datatype extent" );
        if ( width <= 0 || lower != 0 || extent != width || offset < 0 ||
             count > std::size_t( ( std::numeric_limits<MPI_Offset>::max() - offset ) / width ) ||
             ( count && !data ) )
            checkMpiIo( MPI_ERR_COUNT, "invalid or overflowing write range" );
        return width;
    }

    //! \brief Checked independent write, including MPI's 32-bit count limit.
    void writeAt( MPI_File file, MPI_Offset offset, void const *data, std::size_t count,
                  MPI_Datatype datatype, MPI_Status *status ) const
    {
        validateWriteRange( offset, data, count, datatype );
        if ( count > std::size_t( std::numeric_limits<int>::max() ) || offset < 0 )
            checkMpiIo( MPI_ERR_COUNT, "write count/offset exceeds supported range" );
        checkMpiIo( MPI_File_write_at( file, offset, data, int( count ), datatype, status ),
                    "MPI_File_write_at" );
        int written = 0;
        checkMpiIo( MPI_Get_count( status, datatype, &written ), "MPI_Get_count" );
        if ( written != int( count ) )
            checkMpiIo( MPI_ERR_IO, "short MPI_File_write_at" );
    }

    //! \brief Checked collective open, reporting the failing filename.
    void openFile( MPI_Comm comm, char const *filename, int mode, MPI_Info info,
                   MPI_File *file ) const
    {
        checkMpiIo( MPI_File_open( comm, filename, mode, info, file ),
                    "MPI_File_open " + std::string( filename ) );
    }

    //! \brief Checked collective close.
    void closeFile( MPI_File *file ) const
    {
        checkMpiIo( MPI_File_close( file ), "MPI_File_close" );
    }

    //! init the ensight exporter
    FEELPP_NO_EXPORT void init();

    //! write the '' file for ensight
    FEELPP_NO_EXPORT void writeSoSFile() const;

    //! write case file variables
    template <typename Iterator, typename TSt>
    void writeCaseFileVariables( Iterator it, Iterator end, std::string const &loc, TSt const &__ts,
                                 std::ostream &__out, typename super::field_set_ptrtype const &step,
                                 bool isStatic = false ) const;

    //! write the 'case' file for ensight
    FEELPP_NO_EXPORT void writeCaseFile() const;

    //! write the 'geo' file for ensight
    FEELPP_NO_EXPORT void writeGeoFiles( timeset_ptrtype __ts, mesh_ptrtype mesh, int timeIndex,
                                         bool isFirstStep ) const;

    FEELPP_NO_EXPORT void writeGeoMarkers( MPI_File fh,
                                           mesh_contiguous_numbering_mapping_type const &mp,
                                           bool writeHeaderBeginFile, bool writeBeginEndTimeSet,
                                           Feel::detail::FileIndex &index ) const;

    FEELPP_NO_EXPORT void
    writeGeoMarkedFaces( MPI_File fh, mesh_ptrtype mesh,
                         std::pair<const std::string, std::vector<index_type>> &m ) const;

    FEELPP_NO_EXPORT void writeGeoMarkedElements( MPI_File fh,
                                                  mesh_contiguous_numbering_mapping_type const &mp,
                                                  int part ) const;

    //! write the variables file for ensight
    FEELPP_NO_EXPORT void writeVariableFiles( timeset_ptrtype __ts, step_ptrtype step ) const;

    template <bool IsNodal, typename Iterator>
    FEELPP_NO_EXPORT void saveFields( timeset_ptrtype __ts,
                                      typename super::field_set_ptrtype __step, bool writeNewFile,
                                      std::string const &filenameStepIndex, bool isFirstStep,
                                      Iterator __var, Iterator en, bool isStatic = false ) const;

  private:
    //! Requested policy and destination-specific resolution; metadata always uses the master.
    ExporterIOPolicy M_ioPolicy = ExporterIOPolicy::Automatic;
    mutable std::string M_resolvedIOPath;
    mutable bool M_useRootIO = false;
    //! Private message context prevents collisions with application communication.
    mutable std::optional<mpi::communicator> M_payloadComm;
    //! Bounded master-only staging storage, independent of global mesh size/rank count.
    mutable std::vector<char> M_payloadBuffer;
    mutable std::string M_filename;
    std::string M_element_type;
    std::string M_face_type;
    bool M_mergeTimeSteps;
    int M_packTimeSteps;

    // mapping allow to get ordering between Feel++ and Ensight format with curve element
    std::map<std::string, std::vector<uint16_type>> M_nodesOrderingInElementToEnsight;

    /* Number of digits used in timesteps */
    /* Set to 4 by default: range [0000; 9999] for timesteps */
    mutable int M_timeExponent;
    // file position for explicit pointers
    mutable MPI_Offset posInFile;
    mutable std::map<std::string, mesh_contiguous_numbering_mapping_ptrtype> M_cache_mp;
    //! Array-to-DOF maps belong to a time set and a part, not to a part alone.
    mutable std::map<std::string, std::map<int, std::vector<size_type>>> M_mapNodalArrayToDofId;
    mutable std::map<std::string, std::map<int, std::vector<size_type>>> M_mapElementArrayToDofId;
};

} // namespace Feel

// #if !defined( FEELPP_INSTANTIATION_MODE )
#include <feel/feelfilters/exporterensightgold_impl.hpp>
// #endif // FEELPP_INSTANTIATION_MODE

#endif /* __ExporterEnsightGold_H */
