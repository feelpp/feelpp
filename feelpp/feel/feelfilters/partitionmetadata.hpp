/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
#pragma once

#include <feel/feelcore/json.hpp>

#include <algorithm>
#include <functional>
#include <istream>
#include <map>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Feel::detail
{
/**
 * @brief Stores marker values for sparse mesh-fragment identifiers.
 *
 * The table is immutable after finalize(). It stores all marker values in a
 * contiguous buffer and exposes non-owning spans, avoiding an allocation for
 * each fragment lookup. Sparse and legacy physical identifiers are supported.
 */
template <typename Value>
class PartitionMarkerTable
{
    /** @brief Describes one fragment range in the contiguous value buffer. */
    struct Entry
    {
        int id;
        std::size_t offset;
        std::size_t count;
    };

public:
    /**
     * @brief Adds marker values associated with a fragment identifier.
     * @param id Fragment identifier.
     * @param values Marker values to copy into the table.
     */
    void append( int id, std::vector<Value> const& values )
    {
        M_entries.push_back( { id, M_values.size(), values.size() } );
        M_values.insert( M_values.end(), values.begin(), values.end() );
    }

    /**
     * @brief Sorts entries and retains the first entry for every identifier.
     *
     * This must be called before at(). Duplicate identifiers preserve the
     * first occurrence, matching the legacy map::emplace behavior.
     */
    void finalize()
    {
        // Stable sorting preserves the first entry, like the previous map::emplace.
        // Legacy files can have several physical names for the same marker ID.
        std::stable_sort( M_entries.begin(), M_entries.end(), []( auto const& a, auto const& b ) { return a.id < b.id; } );
        M_entries.erase( std::unique( M_entries.begin(), M_entries.end(), []( auto const& a, auto const& b ) { return a.id == b.id; } ), M_entries.end() );
    }

    /**
     * @brief Returns the marker values for a fragment.
     * @param id Fragment identifier.
     * @return A non-owning view into the table storage.
     * @throw std::out_of_range if @p id is absent.
     */
    std::span<Value const> at( int id ) const
    {
        auto it = std::lower_bound( M_entries.begin(), M_entries.end(), id, []( auto const& entry, int value ) { return entry.id < value; } );
        if ( it == M_entries.end() || it->id != id )
            throw std::out_of_range( "unknown mesh fragment ID " + std::to_string( id ) );
        return std::span<Value const>( M_values ).subspan( it->offset, it->count );
    }

    /** @brief Returns the number of distinct fragment identifiers. */
    std::size_t size() const { return M_entries.size(); }

    /** @brief Returns the capacity reserved by the entry and value buffers. */
    std::size_t storageBytes() const { return M_entries.capacity() * sizeof( Entry ) + M_values.capacity() * sizeof( Value ); }

private:
    std::vector<Entry> M_entries;
    std::vector<Value> M_values;
};

/**
 * @brief Streams partition metadata without materializing a JSON document.
 *
 * The reader retains one physical or fragment record at a time. Unknown
 * extensions are ignored, while supported records are type-checked and invalid
 * input raises std::runtime_error.
 *
 * @tparam Value Marker value type.
 */
template <typename Value>
class PartitionMetadataReader : public nl::json_sax<nl::json>
{
    /** @brief Identifies the metadata object currently being parsed. */
    enum class State { Root, Mesh, Partition, Physicals, Fragmentation, Codimension, Physical, Fragment, Ignore };

    /** @brief Tracks the state and key of one nested JSON object or array. */
    struct Frame
    {
        State state;
        std::string key;
        bool array;
        int id = 0;
        int codimension = 0;
    };

public:
    /** @brief Callback invoked for every physical name, identifier and dimension. */
    using PhysicalCallback = std::function<void( std::string const&, int, int )>;

    /**
     * @brief Creates a reader that reports physical metadata through @p physical.
     * @param physical Callback receiving physical name, identifier and dimension.
     */
    explicit PartitionMetadataReader( PhysicalCallback physical ) : M_physical( std::move( physical ) ) {}

    /**
     * @brief Parses mesh partition metadata from a JSON stream.
     * @param input Metadata stream containing a top-level mesh object.
     * @throw std::runtime_error if the stream is malformed or does not describe a mesh.
     */
    void read( std::istream& input )
    {
        if ( !nl::json::sax_parse( input, this ) || !M_hasMesh )
            throw std::runtime_error( "invalid mesh partition metadata: expected mesh object" );
        for ( auto& [codimension, table] : M_fragments )
            table.finalize();
    }

    /** @brief HDF5 filename declared by the metadata, when present. */
    std::optional<std::string> M_h5Filename;
    /** @brief Declared number of partitions, when present. */
    std::optional<int> M_partitions;
    /** @brief True when the metadata explicitly contains fragment mappings. */
    bool M_hasFragmentation = false;
    /** @brief Marker tables indexed by mesh codimension. */
    std::map<int, PartitionMarkerTable<Value>> M_fragments;

    bool null() override { return scalar( nullptr ); }
    bool boolean( bool value ) override { return scalar( value ); }
    bool number_integer( number_integer_t value ) override { return scalar( value ); }
    bool number_unsigned( number_unsigned_t value ) override { return scalar( value ); }
    bool number_float( number_float_t value, string_t const& ) override { return scalar( value ); }
    bool string( string_t& value ) override { return scalar( value ); }
    bool binary( binary_t& ) override { throw std::runtime_error( "binary mesh metadata is unsupported" ); }
    bool start_object( std::size_t ) override { return start( false ); }
    bool start_array( std::size_t ) override { return start( true ); }
    bool key( string_t& value ) override
    {
        M_stack.back().key = value;
        return true;
    }
    bool end_object() override { return end(); }
    bool end_array() override { return end(); }
    bool parse_error( std::size_t position, std::string const&, nl::detail::exception const& error ) override
    {
        throw std::runtime_error( "invalid mesh partition metadata at byte " + std::to_string( position ) + ": " + error.what() );
    }

private:
    static void require( bool valid )
    {
        if ( !valid )
            throw std::runtime_error( "invalid mesh partition metadata record" );
    }

    static int integer( nl::json const& value )
    {
        require( value.is_number_integer() || value.is_string() );
        return value.is_string() ? std::stoi( value.template get<std::string>() ) : value.template get<int>();
    }

    bool scalar( nl::json const& value )
    {
        require( !M_stack.empty() );
        auto const& frame = M_stack.back();
        switch ( frame.state )
        {
        case State::Root:
            require( frame.key != "mesh" );
            break;
        case State::Mesh:
            if ( frame.key == "h5" )
                M_h5Filename = value.template get<std::string>();
            else if ( frame.key == "partition" || frame.key == "physicals" || frame.key == "fragmentation" )
                require( false );
            break;
        case State::Partition:
            if ( frame.key == "n" )
                M_partitions = integer( value );
            break;
        case State::Physicals:
        {
            require( value.is_string() );
            std::istringstream physical( value.template get<std::string>() );
            int id, dimension;
            std::string extra;
            require( bool( physical >> id >> dimension ) && !( physical >> extra ) );
            M_physical( frame.key, id, dimension );
            break;
        }
        case State::Physical:
            M_values.push_back( integer( value ) );
            require( M_values.size() <= 2 );
            break;
        case State::Fragment:
            require( value.is_number_integer() );
            M_values.push_back( value.template get<Value>() );
            break;
        case State::Fragmentation:
        case State::Codimension:
            require( false );
            break;
        case State::Ignore:
            break;
        }
        return true;
    }

    bool start( bool array )
    {
        Frame next{ State::Ignore, {}, array };
        if ( M_stack.empty() )
        {
            require( !array );
            next.state = State::Root;
        }
        else
        {
            auto const& parent = M_stack.back();
            if ( parent.state == State::Root && parent.key == "mesh" )
            {
                require( !array && !M_hasMesh );
                M_hasMesh = true;
                next.state = State::Mesh;
            }
            else if ( parent.state == State::Mesh )
            {
                if ( parent.key == "h5" )
                    require( false );
                if ( parent.key == "partition" )
                    next.state = State::Partition;
                else if ( parent.key == "physicals" )
                    next.state = State::Physicals;
                else if ( parent.key == "fragmentation" )
                {
                    next.state = State::Fragmentation;
                    M_hasFragmentation = true;
                }
                require( next.state == State::Ignore || !array );
            }
            else if ( parent.state == State::Physicals )
            {
                require( array );
                next.state = State::Physical;
                M_values.clear();
            }
            else if ( parent.state == State::Fragmentation )
            {
                require( !array );
                next.state = State::Codimension;
                next.codimension = std::stoi( parent.key );
                M_fragments[next.codimension];
            }
            else if ( parent.state == State::Codimension )
            {
                require( array );
                next.state = State::Fragment;
                next.codimension = parent.codimension;
                next.id = std::stoi( parent.key );
                M_values.clear();
            }
            else if ( parent.state == State::Physical || parent.state == State::Fragment ||
                      ( parent.state == State::Partition && parent.key == "n" ) )
                require( false );
        }
        M_stack.push_back( std::move( next ) );
        return true;
    }

    bool end()
    {
        auto const& frame = M_stack.back();
        if ( frame.state == State::Physical )
        {
            require( M_values.size() == 2 );
            M_physical( M_stack[M_stack.size() - 2].key, M_values[0], M_values[1] );
        }
        else if ( frame.state == State::Fragment )
            M_fragments[frame.codimension].append( frame.id, M_values );
        M_stack.pop_back();
        return true;
    }

    PhysicalCallback M_physical;
    bool M_hasMesh = false;
    std::vector<Frame> M_stack;
    std::vector<Value> M_values;
};
}
