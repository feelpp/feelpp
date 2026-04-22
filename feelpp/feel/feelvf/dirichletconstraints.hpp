#ifndef FEELPP_VF_DIRICHLETCONSTRAINTS_HPP
#define FEELPP_VF_DIRICHLETCONSTRAINTS_HPP 1

#include <algorithm>
#include <compare>
#include <cstdint>
#include <functional>
#include <iterator>
#include <map>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <set>
#include <utility>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include <feel/feelalg/vector.hpp>
#include <feel/feelcore/context.hpp>
#include <feel/feelcore/typetraits.hpp>
#include <feel/feelmesh/enums.hpp>

namespace Feel::vf
{

inline bool isEliminationDirichlet( Feel::Context const& onContext ) noexcept
{
    return onContext.test( ContextOn::ELIMINATION );
}

inline bool isPenalisationDirichlet( Feel::Context const& onContext ) noexcept
{
    return !onContext.test( ContextOn::ELIMINATION | ContextOn::SYMMETRIC );
}

inline bool supportsDeferredDirichlet( Feel::Context const& onContext ) noexcept
{
    return isEliminationDirichlet( onContext ) || isPenalisationDirichlet( onContext );
}

enum class DeferredDirichletPolicy : std::uint8_t
{
    automatic = 0,
    deferred = 1,
    immediate = 2
};

constexpr bool usesDeferredDirichlet( DeferredDirichletPolicy policy ) noexcept
{
    return policy != DeferredDirichletPolicy::immediate;
}

inline bool shouldDeferDirichlet( DeferredDirichletPolicy policy,
                                  Feel::Context const& onContext,
                                  bool allowDeferred = true ) noexcept
{
    return allowDeferred &&
           supportsDeferredDirichlet( onContext ) &&
           usesDeferredDirichlet( policy );
}

template<typename T>
constexpr T deferredDirichletPenalty()
{
    return static_cast<T>( 1e30 );
}

enum class DeferredDirichletEntity : std::uint8_t
{
    unspecified = 0,
    element = 1,
    face = 2,
    edge = 3,
    point = 4
};

constexpr std::uint8_t deferredDirichletEntityPriority( DeferredDirichletEntity entity ) noexcept
{
    return static_cast<std::uint8_t>( entity );
}

template<ElementsType ET>
constexpr std::uint8_t deferredDirichletEntityPriority() noexcept
{
    if constexpr ( ET == MESH_ELEMENTS )
        return deferredDirichletEntityPriority( DeferredDirichletEntity::element );
    else if constexpr ( ET == MESH_FACES )
        return deferredDirichletEntityPriority( DeferredDirichletEntity::face );
    else if constexpr ( ET == MESH_EDGES )
        return deferredDirichletEntityPriority( DeferredDirichletEntity::edge );
    else if constexpr ( ET == MESH_POINTS )
        return deferredDirichletEntityPriority( DeferredDirichletEntity::point );
    else
        return deferredDirichletEntityPriority( DeferredDirichletEntity::unspecified );
}

template<typename T>
class DeferredDirichletSet
{
public:
    using value_type = T;
    using real_type = typename type_traits<value_type>::real_type;

    struct Entry
    {
        std::vector<int> dofs;
        std::vector<value_type> values;
        Feel::Context onContext;
        double valueOnDiagonal = 1.0;
        std::uint8_t entityPriority = deferredDirichletEntityPriority( DeferredDirichletEntity::unspecified );
    };

    struct CompatibilityKey
    {
        Feel::size_type onContext = 0;
        double valueOnDiagonal = 1.0;

        auto operator<=>( CompatibilityKey const& ) const = default;
    };

    using entry_type = Entry;
    using entry_list_type = std::vector<entry_type>;

    void append( entry_type const& entry )
    {
        this->append( entry.dofs, entry.values, entry.onContext, entry.valueOnDiagonal, entry.entityPriority );
    }

    void append( std::vector<int> dofs,
                 std::vector<value_type> values,
                 Feel::Context const& onContext,
                 double valueOnDiagonal,
                 std::uint8_t entityPriority = deferredDirichletEntityPriority( DeferredDirichletEntity::unspecified ) )
    {
        M_entries.push_back( normalizeEntry( std::move( dofs ), std::move( values ), onContext, valueOnDiagonal, entityPriority ) );
    }

    void append( DeferredDirichletSet const& other )
    {
        for ( auto const& entry : other.entries() )
            this->append( entry );
    }

    [[nodiscard]] bool empty() const noexcept { return M_entries.empty(); }
    [[nodiscard]] std::size_t size() const noexcept { return M_entries.size(); }

    void clear() noexcept { M_entries.clear(); }

    [[nodiscard]] entry_list_type const& entries() const noexcept { return M_entries; }

    [[nodiscard]] entry_list_type mergedEntries() const
    {
        struct Candidate
        {
            value_type value = value_type( 0 );
            std::uint8_t entityPriority = deferredDirichletEntityPriority( DeferredDirichletEntity::unspecified );
            std::size_t order = 0;
            bool initialized = false;
        };

        struct GroupState
        {
            CompatibilityKey key;
            std::map<int, Candidate> candidates;
        };

        std::vector<GroupState> groupedEntries;
        groupedEntries.reserve( M_entries.size() );

        auto const isBetterCandidate = []( Candidate const& current,
                                           std::uint8_t entityPriority,
                                           std::size_t order,
                                           value_type const& value )
        {
            if ( !current.initialized )
                return true;
            if ( entityPriority != current.entityPriority )
                return entityPriority > current.entityPriority;
            if ( order != current.order )
                return order > current.order;
            return !equivalentValues( current.value, value );
        };

        for ( std::size_t entryIndex = 0; entryIndex < M_entries.size(); ++entryIndex )
        {
            auto const& entry = M_entries[entryIndex];
            CompatibilityKey const key = compatibilityKey( entry );
            auto it = std::ranges::find( groupedEntries, key, []( auto const& groupedEntry ) { return groupedEntry.key; } );
            if ( it == groupedEntries.end() )
            {
                groupedEntries.push_back( GroupState{ .key = key } );
                it = std::prev( groupedEntries.end() );
            }

            for ( std::size_t k = 0; k < entry.dofs.size(); ++k )
            {
                auto& candidate = it->candidates[entry.dofs[k]];
                if ( isBetterCandidate( candidate, entry.entityPriority, entryIndex, entry.values[k] ) )
                {
                    candidate.value = entry.values[k];
                    candidate.entityPriority = entry.entityPriority;
                    candidate.order = entryIndex;
                    candidate.initialized = true;
                }
            }
        }

        entry_list_type merged;
        merged.reserve( groupedEntries.size() );
        for ( auto& groupedEntry : groupedEntries )
        {
            entry_type mergedEntry = makeEmptyEntry( groupedEntry.key );
            mergedEntry.dofs.reserve( groupedEntry.candidates.size() );
            mergedEntry.values.reserve( groupedEntry.candidates.size() );
            for ( auto const& [dof, candidate] : groupedEntry.candidates )
            {
                mergedEntry.dofs.push_back( dof );
                mergedEntry.values.push_back( candidate.value );
            }
            merged.push_back( std::move( mergedEntry ) );
        }
        return merged;
    }

private:
    using dof_value_type = std::pair<int, value_type>;

    static CompatibilityKey compatibilityKey( entry_type const& entry )
    {
        return CompatibilityKey{
            .onContext = entry.onContext.context(),
            .valueOnDiagonal = entry.valueOnDiagonal
        };
    }

    static entry_type makeEmptyEntry( CompatibilityKey const& key )
    {
        return entry_type{
            .dofs = {},
            .values = {},
            .onContext = Feel::Context( key.onContext ),
            .valueOnDiagonal = key.valueOnDiagonal,
            .entityPriority = deferredDirichletEntityPriority( DeferredDirichletEntity::unspecified )
        };
    }

    static entry_type normalizeEntry( std::vector<int> dofs,
                                      std::vector<value_type> values,
                                      Feel::Context const& onContext,
                                      double valueOnDiagonal,
                                      std::uint8_t entityPriority )
    {
        if ( dofs.size() != values.size() )
        {
            std::ostringstream msg;
            msg << "invalid deferred Dirichlet data: dofs=" << dofs.size()
                << " values=" << values.size();
            throw std::invalid_argument( msg.str() );
        }

        std::vector<dof_value_type> dofValues;
        dofValues.reserve( dofs.size() );
        for ( std::size_t k = 0; k < dofs.size(); ++k )
            dofValues.emplace_back( dofs[k], std::move( values[k] ) );

        std::ranges::sort( dofValues, std::less<>{}, &dof_value_type::first );

        entry_type normalized = makeEmptyEntry( CompatibilityKey{
            .onContext = onContext.context(),
            .valueOnDiagonal = valueOnDiagonal
        } );
        normalized.entityPriority = entityPriority;
        normalized.dofs.reserve( dofValues.size() );
        normalized.values.reserve( dofValues.size() );

        for ( auto const& [dof, value] : dofValues )
        {
            if ( !normalized.dofs.empty() && normalized.dofs.back() == dof )
            {
                if ( !equivalentValues( normalized.values.back(), value ) )
                    throw std::invalid_argument( incompatibleValueMessage( dof, normalized.values.back(), value ) );
                continue;
            }

            normalized.dofs.push_back( dof );
            normalized.values.push_back( value );
        }

        return normalized;
    }

    static bool equivalentValues( value_type const& lhs, value_type const& rhs )
    {
        real_type const scale = std::max( { real_type( 1 ), math::abs( lhs ), math::abs( rhs ) } );
        real_type const tolerance = real_type( 100 ) * type_traits<value_type>::epsilon() * scale;
        return math::abs( lhs - rhs ) <= tolerance;
    }

    static std::string incompatibleValueMessage( int dof, value_type const& lhs, value_type const& rhs )
    {
        std::ostringstream msg;
        msg << "conflicting deferred Dirichlet values for dof " << dof
            << ": " << lhs << " vs " << rhs;
        return msg.str();
    }

private:
    entry_list_type M_entries;
};

template<typename VectorPtrType>
auto cloneVectorWithValues( VectorPtrType const& source )
{
    auto copy = source->clone();
    copy->zero();
    copy->add( typename std::decay_t<decltype( *copy )>::value_type( 1 ), *source );
    if ( !copy->closed() )
        copy->close();
    return copy;
}

template<typename VectorPtrType, typename SourceVectorPtrType>
void copyVectorValues( VectorPtrType const& destination,
                      SourceVectorPtrType const& source )
{
    destination->zero();
    destination->add( typename std::decay_t<decltype( *destination )>::value_type( 1 ), *source );
    if ( !destination->closed() )
        destination->close();
}

template<typename VectorPtrType, typename CandidateMapType>
void synchronizeDeferredDirichletCandidates( VectorPtrType const& values,
                                             CandidateMapType const& localCandidates )
{
    using vector_type = std::decay_t<decltype( *values )>;
    using value_type = typename vector_type::value_type;
    using size_type = typename vector_type::size_type;

    struct Candidate
    {
        std::uint8_t entityPriority = deferredDirichletEntityPriority( DeferredDirichletEntity::unspecified );
        std::size_t order = 0;
        value_type value = value_type( 0 );
        rank_type rank = invalid_rank_type_value;
        bool initialized = false;
    };

    auto const dataMap = values->mapPtr();
    if ( !dataMap )
    {
        for ( auto const& [gpdof, candidate] : localCandidates )
            values->set( gpdof, candidate.value );
        return;
    }

    if ( dataMap->worldComm().localSize() == 1 )
    {
        for ( auto const& [gpdof, candidate] : localCandidates )
            values->set( gpdof, candidate.value );
        return;
    }

    auto const isBetterCandidate = []( Candidate const& current,
                                       std::uint8_t entityPriority,
                                       std::size_t order,
                                       value_type const& value,
                                       rank_type rank )
    {
        if ( !current.initialized )
            return true;
        if ( entityPriority != current.entityPriority )
            return entityPriority > current.entityPriority;
        if ( order != current.order )
            return order > current.order;
        if ( rank != current.rank )
            return rank < current.rank;
        return math::abs( current.value - value ) > type_traits<value_type>::epsilon();
    };

    using owner_payload_type = boost::tuple<size_type, std::uint8_t, std::uint64_t, value_type>;
    using value_payload_type = boost::tuple<size_type, value_type>;

    rank_type const localRank = dataMap->worldComm().localRank();
    std::map<size_type, Candidate> ownerCandidates;
    std::map<rank_type, std::vector<owner_payload_type>> dataToOwner, dataFromGhost;

    for ( auto const& [gpdof, candidate] : localCandidates )
    {
        size_type const gcdof = dataMap->mapGlobalProcessToGlobalCluster( gpdof );
        if ( dataMap->dofGlobalProcessIsGhost( gpdof ) )
        {
            rank_type const ownerRank = dataMap->procOnGlobalCluster( gcdof );
            CHECK( ownerRank != invalid_rank_type_value ) << "proc not found for gcdof: " << gcdof;
            dataToOwner[ownerRank].push_back( boost::make_tuple( gcdof, candidate.entityPriority, static_cast<std::uint64_t>( candidate.order ), candidate.value ) );
            continue;
        }

        auto& ownerCandidate = ownerCandidates[gcdof];
        if ( isBetterCandidate( ownerCandidate, candidate.entityPriority, candidate.order, candidate.value, localRank ) )
        {
            ownerCandidate.entityPriority = candidate.entityPriority;
            ownerCandidate.order = candidate.order;
            ownerCandidate.value = candidate.value;
            ownerCandidate.rank = localRank;
            ownerCandidate.initialized = true;
        }
    }

    int const nbRequest = 2 * dataMap->neighborSubdomains().size();
    auto* reqs = new mpi::request[nbRequest];
    std::map<rank_type, std::size_t> sizeRecv, sizeSend;
    int cptRequest = 0;

    for ( rank_type neighborRank : dataMap->neighborSubdomains() )
    {
        sizeSend[neighborRank] = dataToOwner[neighborRank].size();
        reqs[cptRequest++] = dataMap->worldComm().localComm().isend( neighborRank, 0, sizeSend[neighborRank] );
        reqs[cptRequest++] = dataMap->worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
    }
    mpi::wait_all( reqs, reqs + cptRequest );

    cptRequest = 0;
    for ( rank_type neighborRank : dataMap->neighborSubdomains() )
    {
        std::size_t const nSendData = dataToOwner[neighborRank].size();
        if ( nSendData > 0 )
            reqs[cptRequest++] = dataMap->worldComm().localComm().isend( neighborRank, 0, dataToOwner[neighborRank].data(), nSendData );

        std::size_t const nRecvData = sizeRecv[neighborRank];
        dataFromGhost[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[cptRequest++] = dataMap->worldComm().localComm().irecv( neighborRank, 0, dataFromGhost[neighborRank].data(), nRecvData );
    }
    mpi::wait_all( reqs, reqs + cptRequest );

    for ( auto const& [remoteRank, payloads] : dataFromGhost )
    {
        for ( auto const& payload : payloads )
        {
            size_type const gcdof = boost::get<0>( payload );
            auto const entityPriority = boost::get<1>( payload );
            auto const order = static_cast<std::size_t>( boost::get<2>( payload ) );
            auto const& value = boost::get<3>( payload );
            auto& ownerCandidate = ownerCandidates[gcdof];
            if ( isBetterCandidate( ownerCandidate, entityPriority, order, value, remoteRank ) )
            {
                ownerCandidate.entityPriority = entityPriority;
                ownerCandidate.order = order;
                ownerCandidate.value = value;
                ownerCandidate.rank = remoteRank;
                ownerCandidate.initialized = true;
            }
        }
    }

    std::map<rank_type, std::vector<value_payload_type>> dataToGhost, dataFromOwner;
    for ( auto const& [gcdof, ownerCandidate] : ownerCandidates )
    {
        if ( !ownerCandidate.initialized )
            continue;

        size_type const gpdof = dataMap->worldIndexToProcessIndex( gcdof );
        DCHECK( !dataMap->dofGlobalProcessIsGhost( gpdof ) ) << "gpdof " << gpdof << " must be active";
        values->set( gpdof, ownerCandidate.value );

        auto itShared = dataMap->activeDofSharedOnCluster().find( gpdof );
        if ( itShared == dataMap->activeDofSharedOnCluster().end() )
            continue;

        for ( rank_type neighborRank : itShared->second )
        {
            if ( neighborRank == localRank )
                continue;
            dataToGhost[neighborRank].push_back( boost::make_tuple( gcdof, ownerCandidate.value ) );
        }
    }

    cptRequest = 0;
    sizeRecv.clear();
    sizeSend.clear();
    for ( rank_type neighborRank : dataMap->neighborSubdomains() )
    {
        sizeSend[neighborRank] = dataToGhost[neighborRank].size();
        reqs[cptRequest++] = dataMap->worldComm().localComm().isend( neighborRank, 0, sizeSend[neighborRank] );
        reqs[cptRequest++] = dataMap->worldComm().localComm().irecv( neighborRank, 0, sizeRecv[neighborRank] );
    }
    mpi::wait_all( reqs, reqs + cptRequest );

    cptRequest = 0;
    for ( rank_type neighborRank : dataMap->neighborSubdomains() )
    {
        std::size_t const nSendData = dataToGhost[neighborRank].size();
        if ( nSendData > 0 )
            reqs[cptRequest++] = dataMap->worldComm().localComm().isend( neighborRank, 0, dataToGhost[neighborRank].data(), nSendData );

        std::size_t const nRecvData = sizeRecv[neighborRank];
        dataFromOwner[neighborRank].resize( nRecvData );
        if ( nRecvData > 0 )
            reqs[cptRequest++] = dataMap->worldComm().localComm().irecv( neighborRank, 0, dataFromOwner[neighborRank].data(), nRecvData );
    }
    mpi::wait_all( reqs, reqs + cptRequest );
    delete [] reqs;

    for ( auto const& [remoteRank, payloads] : dataFromOwner )
    {
        (void)remoteRank;
        for ( auto const& payload : payloads )
        {
            size_type const gcdof = boost::get<0>( payload );
            size_type const gpdof = dataMap->worldIndexToProcessIndex( gcdof );
            DCHECK( dataMap->dofGlobalProcessIsGhost( gpdof ) ) << "gpdof " << gpdof << " must be ghost";
            values->set( gpdof, boost::get<1>( payload ) );
        }
    }
}

template<typename EntryRange, typename MatrixPtrType, typename VectorPtrType>
void applyDeferredDirichletEntries( EntryRange const& entries,
                                    MatrixPtrType const& matrix,
                                    VectorPtrType const& rhsVector )
{
    using entry_type = std::ranges::range_value_t<EntryRange>;
    using key_type = std::pair<Feel::size_type, double>;
    using vector_type = std::decay_t<decltype( *rhsVector )>;
    using value_type = typename vector_type::value_type;
    using size_type = typename vector_type::size_type;

    struct GroupState
    {
        key_type key;
        VectorPtrType values;
        std::set<int> dofSet;
        struct LocalCandidate
        {
            std::uint8_t entityPriority = deferredDirichletEntityPriority( DeferredDirichletEntity::unspecified );
            std::size_t order = 0;
            value_type value = value_type( 0 );
        };
        std::map<size_type, LocalCandidate> localCandidates;
    };

    std::vector<GroupState> groupedStates;
    std::size_t entryIndex = 0;

    auto findOrCreateGroupState = [&]( key_type const& key ) -> GroupState&
    {
        auto it = std::ranges::find( groupedStates, key, []( auto const& groupedState ) { return groupedState.key; } );
        if ( it == groupedStates.end() )
        {
            groupedStates.push_back( GroupState{
                .key = key,
                .values = rhsVector->clone(),
                .dofSet = {}
            } );
            groupedStates.back().values->zero();
            it = std::prev( groupedStates.end() );
        }
        return *it;
    };

    for ( auto const& bc : entries )
    {
        if ( bc.dofs.size() != bc.values.size() )
        {
            std::ostringstream msg;
            msg << "invalid deferred Dirichlet data: dofs=" << bc.dofs.size()
                << " values=" << bc.values.size();
            throw std::invalid_argument( msg.str() );
        }

        key_type const key{ bc.onContext.context(), bc.valueOnDiagonal };
        auto& groupState = findOrCreateGroupState( key );

        if ( !bc.dofs.empty() )
        {
            groupState.dofSet.insert( bc.dofs.begin(), bc.dofs.end() );
            for ( std::size_t k = 0; k < bc.dofs.size(); ++k )
            {
                auto& candidate = groupState.localCandidates[static_cast<size_type>( bc.dofs[k] )];
                if ( bc.entityPriority > candidate.entityPriority ||
                     ( bc.entityPriority == candidate.entityPriority && entryIndex >= candidate.order ) )
                {
                    candidate.entityPriority = bc.entityPriority;
                    candidate.order = entryIndex;
                    candidate.value = bc.values[k];
                }
            }
        }
        ++entryIndex;
    }

    for ( auto& groupState : groupedStates )
    {
        if ( groupState.dofSet.empty() )
            continue;

        synchronizeDeferredDirichletCandidates( groupState.values, groupState.localCandidates );

        std::vector<int> dofs( groupState.dofSet.begin(), groupState.dofSet.end() );
        Feel::Context const onContext( groupState.key.first );

        if ( isPenalisationDirichlet( onContext ) )
        {
            auto const penalty = deferredDirichletPenalty<value_type>();
            for ( int dof : dofs )
            {
                matrix->set( dof, dof, penalty );
                rhsVector->set( dof, groupState.values->operator()( dof ) * penalty );
            }
            continue;
        }

        if ( !groupState.values->closed() )
            groupState.values->close();
        matrix->zeroRows( dofs, *groupState.values, *rhsVector, onContext, groupState.key.second );
    }
}

} // namespace Feel::vf

#endif /* FEELPP_VF_DIRICHLETCONSTRAINTS_HPP */
