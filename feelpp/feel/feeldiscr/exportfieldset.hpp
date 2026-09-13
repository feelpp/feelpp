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
//! \file exportfieldset.hpp
//! \brief Time-neutral export fields shared by exporters and temporal steps.
//!
//! Composition separates field lifetime from file encoding. This container owns
//! conversion buffers and persistent schemas, but has no time, step index,
//! TimeSet pointer or exporter pointer. A writer may store a dataset snapshot
//! once or materialize the same values in temporal records.
//! Owners have distinct field sets but may share compatible conversion spaces.
//! Registration and mutation are collective where required and are not thread-safe.
#ifndef FEELPP_DISCR_EXPORTFIELDSET_HPP
#define FEELPP_DISCR_EXPORTFIELDSET_HPP 1
#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <fmt/format.h>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelcomplex.hpp>

#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/lagrange.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/interpolate.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>

#include <feel/feeldiscr/elementdiv.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

namespace detail
{
//! \brief Time-neutral scalar/vector/tensor export storage and conversion.
//! \tparam MeshType Mesh defining the field layout.
//! \tparam N Geometric/export nodal order.
template <typename MeshType, int N = 1> class ExportFieldSet
{
  public:
    //! \name Mesh, conversion-space, field-schema and iterator types
    //@{
    using mesh_type = MeshType;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;
    using size_type = typename mesh_type::size_type;

    using scalar_p0_space_type = Pdh_type<MeshType, 0>;
    using scalar_p1_space_type = Pch_type<MeshType, N /*1*/>;
    using scalar_p0_space_ptrtype = std::shared_ptr<scalar_p0_space_type>;
    using scalar_p1_space_ptrtype = std::shared_ptr<scalar_p1_space_type>;

    using element_scalar_type = typename scalar_p0_space_type::element_type;
    using nodal_scalar_type = typename scalar_p1_space_type::element_type;

    using nodal_scalar_ptrtype = std::shared_ptr<nodal_scalar_type>;
    using element_scalar_ptrtype = std::shared_ptr<element_scalar_type>;
    using nodal_field_type =
        std::pair<FunctionSpaceType, std::vector<std::vector<nodal_scalar_ptrtype>>>;
    using element_field_type =
        std::pair<FunctionSpaceType, std::vector<std::vector<element_scalar_ptrtype>>>;

    using map_scalar_type = std::map<std::string, std::pair<scalar_type, bool>>;
    using map_complex_type = std::map<std::string, std::pair<complex_type, bool>>;
    using map_nodal_type = std::map<std::string, nodal_field_type>;
    using map_element_type = std::map<std::string, element_field_type>;

    using scalar_iterator = typename map_scalar_type::iterator;
    using scalar_const_iterator = typename map_scalar_type::const_iterator;
    using complex_iterator = typename map_complex_type::iterator;
    using complex_const_iterator = typename map_complex_type::const_iterator;
    using nodal_iterator = typename map_nodal_type::iterator;
    using nodal_const_iterator = typename map_nodal_type::const_iterator;
    using element_iterator = typename map_element_type::iterator;
    using element_const_iterator = typename map_element_type::const_iterator;

    using variant_representation_arg_type = std::variant<std::string, std::set<std::string>>;

    //@}

    //! \brief Reusable scalar conversion spaces, independent of field lifetime.
    //! Local references avoid repeated manager lookups within a sequence or snapshot.
    struct SpaceCache
    {
        scalar_p0_space_ptrtype M_p0; //!< Element conversion space.
        scalar_p1_space_ptrtype M_p1; //!< Nodal conversion space.
    };

    //! \brief Construct empty field storage with an optional shared conversion cache.
    //! \param cache Compatible conversion spaces, shared without temporal ownership.
    explicit ExportFieldSet( std::shared_ptr<SpaceCache> cache = {} )
        : M_spaceCache( cache ? std::move( cache ) : std::make_shared<SpaceCache>() )
    {
    }

    //! \name Accessors
    //@{

    //! \return Whether geometry or field data has been associated.
    [[nodiscard]] bool hasData() const { return M_hasData; }

    //! \return Whether component buffers are marked resident.
    [[nodiscard]] bool isInMemory() const { return M_isInMemory; }

    //! \return Content revision, independent of output time and publication state.
    //! Registration and mesh changes advance this value; cleanup does not.
    [[nodiscard]] std::uint64_t revision() const { return M_revision; }

    //! \brief Restore residency metadata without reading or allocating field values.
    //! \param hasData Whether geometry or field data was recorded.
    //! \param isInMemory Whether resident component buffers were recorded.
    //! Used for legacy metadata serialization; this is not payload loading.
    void restoreDataState( bool hasData, bool isInMemory )
    {
        M_hasData = hasData;
        M_isInMemory = isInMemory;
    }

    //! \return true if the mesh is available
    bool hasMesh() const { return M_mesh != std::nullopt; }

    //! \return a mesh
    mesh_ptrtype mesh() const { return M_mesh.value(); }

    //! \return the begin iterator for scalars
    scalar_const_iterator beginScalar() const { return M_scalar.begin(); }

    //! \return the end iterator for scalars
    scalar_const_iterator endScalar() const { return M_scalar.end(); }

    //! get the scalar with name n
    //! \param name name of the nodal scalar field
    //! \return the scalar value
    scalar_type scalar( std::string const &name ) const
    {
        if ( M_scalar.find( sanitize( name ) ) == M_scalar.end() )
        {
            std::ostringstream error;
            error << "invalid scalar value name " << sanitize( name );
            throw std::logic_error( error.str() );
        }

        return M_scalar.find( sanitize( name ) )->second.first;
    }

    //! get the nodal field with name n
    //! \param name name of the nodal field
    //! \return the nodal field
    nodal_field_type const &nodal( std::string const &name ) const
    {
        auto itFindNodalField = M_nodal.find( sanitize( name ) );
        CHECK( itFindNodalField != M_nodal.end() )
            << "invalid nodal field name " << sanitize( name );
        return itFindNodalField->second;
    }

    //! get the element field with name n
    //! \param name name of the element field
    //! \return the element field
    element_field_type const &element( std::string const &name ) const
    {
        auto itFindElementField = M_element.find( sanitize( name ) );
        CHECK( itFindElementField != M_element.end() )
            << "invalid nodal field name " << sanitize( name );
        return itFindElementField->second;
    }

    //@}

    //! \name  Mutators
    //@{

    //! \brief Associate the mesh used for conversion and output ordering.
    void setMesh( mesh_ptrtype const &mesh )
    {
        DVLOG( 2 ) << "[ExportFieldSet::setMesh] setMesh start\n";
        M_mesh = mesh;

        markModified();
        DVLOG( 2 ) << "[ExportFieldSet::setMesh] setMesh done\n";
    }

    //! \brief Legacy scalar registration; prefer the floating-point add overload.
    FEELPP_DEPRECATED void addScalar( std::string const &name, scalar_type const &value,
                                      bool cst = false )
    {
        validateFieldName( sanitize( name ), sanitize( name ) );
        M_scalar[sanitize( name )] = std::make_pair( value, cst );
        markModified();
    }

    //! \brief Store a scalar; cst selects the backend's per-case constant encoding.
    template <typename T>
    void add( std::string const &name, T const &value, bool cst = false,
              typename std::enable_if<std::is_floating_point<T>::value>::type * = nullptr )
    {
        validateFieldName( sanitize( name ), sanitize( name ) );
        M_scalar[sanitize( name )] = std::make_pair( value, cst );
        markModified();
    }

    //! \brief Store a complex quantity using the existing backend convention.
    void addComplex( std::string const &name, complex_type const &value, bool cst = false )
    {
        M_complex[sanitize( name )] = std::pair{ value, cst };
        markModified();
    }

    //! \brief Add generated mesh-region fields to this container.
    //! \details some regions can be automatically generated such as the
    //! process id map
    //!
    //! \param  prefix prefix string for the region
    void addRegions( std::string const &prefix = "" ) { this->addRegions( prefix, prefix ); }

    //! \brief Add partition IDs with separate display-name and filename prefixes.
    void addRegions( std::string const &prefix, std::string const &prefixfname )
    {

        VLOG( 1 ) << "[ExportFieldSet] Adding regions...\n";

        auto scalarElementSpace = this->scalarFunctionSpace<false>();

        VLOG( 1 ) << "[ExportFieldSet] adding pid...\n";
        this->add( prefixvm( prefix, "pid" ), prefixvm( prefixfname, "pid" ),
                   regionProcess( scalarElementSpace ) );
    }

    //! \brief Convert an FE field using explicit names for mixed-space components.
    template <typename FunctionType>
    void
    add( std::initializer_list<std::string> name, FunctionType const &func,
         typename std::enable_if<is_functionspace_element_v<decay_type<FunctionType>>>::type * =
             nullptr )
    {

        std::vector<std::string> str( name );
        add( str, func );
    }

    //! \brief Convert single or mixed FE fields; one mixed name expands to name_0, name_1, etc.
    template <typename FunctionType>
    void
    add( std::variant<std::string, std::vector<std::string>> const &name, FunctionType const &func,
         variant_representation_arg_type const &reps = "",
         typename std::enable_if<is_functionspace_element_v<decay_type<FunctionType>>>::type * =
             nullptr )
    {

        std::vector<std::string> names;
        if ( auto nameStringPtr = std::get_if<std::string>( &name ) )
        {
            if ( !nameStringPtr->empty() )
                names.push_back( *nameStringPtr );
        }

        else if ( auto nameVectorPtr = std::get_if<std::vector<std::string>>( &name ) )
        {
            names = *nameVectorPtr;
        }
        CHECK( !names.empty() ) << "no field name given";

        constexpr int nSpaces = decay_type<FunctionType>::functionspace_type::nSpaces;
        if constexpr ( nSpaces > 1 )
        {
            if ( names.size() == 1 )
            {
                std::string nameGiven = names[0];
                names.resize( nSpaces );
                for ( int i = 0; i < nSpaces; i++ )
                    names[i] = fmt::format( "{}_{}", nameGiven, i );
            }

            hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nSpaces> ),
                            [this, &names, &func, &reps]( auto const &e )
                            {
                                constexpr int subId = std::decay_t<decltype( e )>::value;
                                auto const &subFunc = func.template element<subId>();
                                this->add( names[subId], names[subId], subFunc, reps );
                            } );
        }
        else
        {
            add( names[0], names[0], func, reps );
        }
    }

    //! \brief Convert an FE field into the requested nodal/element representations.
    template <typename FunctionType>
    void
    add( std::string const &name, std::string const &filename, FunctionType const &func,
         variant_representation_arg_type const &requestedRepresentations = "",
         typename std::enable_if<is_functionspace_element_v<decay_type<FunctionType>>>::type * =
             nullptr )
    {

        using FunctionDecayType = decay_type<FunctionType>;
        constexpr bool funcIsNodal =
            ( FunctionDecayType::is_continuous ||
              FunctionDecayType::functionspace_type::continuity_type::is_discontinuous_locally ) &&
            ( !FunctionDecayType::is_hcurl_conforming ) &&
            ( !FunctionDecayType::is_hdiv_conforming );

        std::string defaultRepr = funcIsNodal ? "nodal" : "element";
        std::set<std::string> reps = representationType( requestedRepresentations, defaultRepr );

        std::map<std::string, std::string> repToSuffix = { { "nodal", "_n" }, { "element", "_e" } };

        for ( std::string const &rep : reps )
        {
            std::string nameUsed =
                ( reps.size() > 1 ) ? sanitize( name + repToSuffix[rep] ) : sanitize( name );
            std::string fnameUsed = ( reps.size() > 1 ) ? filename + repToSuffix[rep] : filename;
            if ( rep == "element" )
                add<false, false>( nameUsed, fnameUsed, unwrap_ptr( func ) );
            else if ( rep == "nodal" )
            {
                if constexpr ( funcIsNodal )
                    add<true, false>( nameUsed, fnameUsed, unwrap_ptr( func ) );
                else
                    add<true, true>( nameUsed, fnameUsed, unwrap_ptr( func ) );
            }
        }
    }

    //! \brief Copy/interpolate FE components into scalar export buffers.
    template <bool isNodal, bool isElementToNodal, typename FunctionType>
    FEELPP_NO_EXPORT void add( std::string const &name, std::string const &filename,
                               FunctionType const &func )
    {
        if ( !func.worldComm().isActive() )
            return;

        validateFieldName( name, filename );
        ( isNodal ? M_nodalNames : M_elementNames )[filename] = name;

        std::string reprType = isNodal ? "nodal" : "element";
        tic();
        auto scalarSpace = this->scalarFunctionSpace<isNodal>( func );
        toc( fmt::format( "ExportFieldSet::add get scalar space {}", reprType ),
             Environment::logVerbosityLevel() > 0 );

        tic();
        auto &fieldsMap = this->fields<isNodal>();
        std::vector<ComponentType> mapIndicesToComponent = { ComponentType::X, ComponentType::Y,
                                                             ComponentType::Z };

        // TODO: move in migrate filter
        auto rangeInterpWithRelatedMesh = [&scalarSpace, &func]()
        {
            auto meshDomain = func.functionSpace()->mesh();
            auto meshImage = scalarSpace->mesh();
            auto spaceDomain = func.functionSpace();
            auto spaceImage = scalarSpace;
            bool hasMeshSupportPartialDomain =
                spaceDomain->dof()->hasMeshSupport() &&
                spaceDomain->dof()->meshSupport()->isPartialSupport();
            bool hasMeshSupportPartialImage = spaceImage->dof()->hasMeshSupport() &&
                                              spaceImage->dof()->meshSupport()->isPartialSupport();
            if ( meshImage->isSameMesh( meshDomain ) )
            {
                if ( hasMeshSupportPartialDomain && hasMeshSupportPartialImage )
                    return intersect( elements( support( spaceDomain ) ),
                                      elements( support( spaceImage ) ) );
                else if ( hasMeshSupportPartialDomain )
                    return elements( support( spaceDomain ) );
                else if ( hasMeshSupportPartialImage )
                    return elements( support( spaceImage ) );
                else
                    return elements( meshImage );
            }

            else if ( meshImage->isSubMeshFrom( meshDomain ) )
            {
                if ( hasMeshSupportPartialDomain )
                {
                    // WARNING just for test, need to use migrate() filter
                    return elements( support( spaceImage ) );
                }

                else if ( hasMeshSupportPartialImage )
                    return elements( support( spaceImage ) );
                else
                    return elements( meshImage );
            }

            else if ( meshImage->isParentMeshOf( meshDomain ) )
            {
                if ( hasMeshSupportPartialImage )
                {
                    CHECK( false ) << "TODO : range intersection with related meshes";
                    return intersect( elements( support( spaceDomain ) ),
                                      elements( support( spaceImage ) ) );
                }

                else if ( hasMeshSupportPartialDomain )
                    return elements( support( spaceDomain ) );
                else
                    return elements( meshDomain );
            }

            return elements( meshImage );
        };

        if constexpr ( FunctionType::is_scalar )
        {
            fieldsMap[filename].first = FunctionSpaceType::SCALAR;
            fieldsMap[filename].second.resize(
                1, { scalarSpace->elementPtr( name, func.description() ) } );
            if constexpr ( !isElementToNodal )
            {
                if ( func.functionSpace()->mesh()->isRelatedTo( scalarSpace->mesh() ) )
                {
                    auto rangeElt = rangeInterpWithRelatedMesh();
                    fieldsMap[filename].second[0][0]->on( _range = rangeElt, _expr = idv( func ),
                                                          _close = true );
                }
                else
                    interpolate( scalarSpace, func, *fieldsMap[filename].second[0][0] );
            }
            else
            {
                *fieldsMap[filename].second[0][0] =
                    div( sum( scalarSpace, idv( func ) * meas() ),
                         sum( scalarSpace, meas() ) ); // TODO optimize!!!
                fieldsMap[filename].second[0][0]->setName( name );
            }
        }

        else if constexpr ( FunctionType::is_vectorial )
        {
            fieldsMap[filename].first = FunctionSpaceType::VECTORIAL;
            fieldsMap[filename].second.resize( FunctionType::nComponents );
            for ( int c1 = 0; c1 < FunctionType::nComponents; ++c1 )
                fieldsMap[filename].second[c1] = {
                    scalarSpace->elementPtr( name, func.description() ) };

            if constexpr ( !isElementToNodal )
            {
                if ( func.functionSpace()->mesh()->isRelatedTo( scalarSpace->mesh() ) )
                {
                    auto rangeElt = rangeInterpWithRelatedMesh();
                    for ( int c1 = 0; c1 < FunctionType::nComponents; ++c1 )
                        fieldsMap[filename].second[c1][0]->on(
                            _range = rangeElt,
                            _expr = idv( func.comp( mapIndicesToComponent.at( c1 ) ) ),
                            _close = true );
                }
                else
                    interpolate( scalarSpace, func, fieldsMap[filename].second );
            }
            else
            {
                for ( int c1 = 0; c1 < FunctionType::nComponents; ++c1 )
                {
                    *fieldsMap[filename].second[c1][0] =
                        div( sum( scalarSpace, idv( func )( c1, 0 ) * meas() ),
                             sum( scalarSpace, meas() ) ); // TODO optimize!!!
                    fieldsMap[filename].second[c1][0]->setName( name );
                }
            }
        }

        else if constexpr ( FunctionType::is_tensor2 )
        {
            fieldsMap[filename].first = FunctionSpaceType::TENSOR2;
            fieldsMap[filename].second.resize( FunctionType::nComponents1 );
            for ( int c1 = 0; c1 < FunctionType::nComponents1; ++c1 )
            {
                fieldsMap[filename].second[c1].resize( FunctionType::nComponents2 );
                for ( int c2 = 0; c2 < FunctionType::nComponents2; ++c2 )
                    fieldsMap[filename].second[c1][c2] =
                        scalarSpace->elementPtr( name, func.description() );
            }

            if constexpr ( !isElementToNodal )
            {
                if ( func.functionSpace()->mesh()->isRelatedTo( scalarSpace->mesh() ) )
                {
                    auto rangeElt = rangeInterpWithRelatedMesh();
                    for ( int c1 = 0; c1 < FunctionType::nComponents1; ++c1 )
                        for ( int c2 = 0; c2 < FunctionType::nComponents2; ++c2 )
                            fieldsMap[filename].second[c1][c2]->on(
                                _range = rangeElt,
                                _expr = idv( func.comp( mapIndicesToComponent.at( c1 ),
                                                        mapIndicesToComponent.at( c2 ) ) ),
                                _close = true );
                }
                else
                    interpolate( scalarSpace, func, fieldsMap[filename].second );
            }
            else
            {
                for ( int c1 = 0; c1 < FunctionType::nComponents1; ++c1 )
                    for ( int c2 = 0; c2 < FunctionType::nComponents2; ++c2 )
                    {
                        *fieldsMap[filename].second[c1][c2] =
                            div( sum( scalarSpace, idv( func )( c1, c2 ) * meas() ),
                                 sum( scalarSpace, meas() ) ); // TODO optimize!!!
                        fieldsMap[filename].second[c1][c2]->setName( name );
                    }
            }
        }

        else if constexpr ( FunctionType::is_tensor2symm )
        {
            fieldsMap[filename].first = FunctionSpaceType::TENSOR2_SYMM;
            fieldsMap[filename].second.resize( FunctionType::nComponents1 );
            for ( int c1 = 0; c1 < FunctionType::nComponents1; ++c1 )
            {
                fieldsMap[filename].second[c1].resize( c1 + 1 );
                for ( int c2 = 0; c2 <= c1; ++c2 )
                    fieldsMap[filename].second[c1][c2] =
                        scalarSpace->elementPtr( name, func.description() );
            }

            if constexpr ( !isElementToNodal )
            {
                if ( func.functionSpace()->mesh()->isRelatedTo( scalarSpace->mesh() ) )
                {
                    auto rangeElt = rangeInterpWithRelatedMesh();
                    for ( int c1 = 0; c1 < FunctionType::nComponents1; ++c1 )
                        for ( int c2 = 0; c2 <= c1; ++c2 )
                            fieldsMap[filename].second[c1][c2]->on(
                                _range = rangeElt,
                                _expr = idv( func.comp( mapIndicesToComponent.at( c1 ),
                                                        mapIndicesToComponent.at( c2 ) ) ),
                                _close = true );
                }
                else
                    interpolate( scalarSpace, func, fieldsMap[filename].second );
            }
            else
            {
                for ( int c1 = 0; c1 < FunctionType::nComponents1; ++c1 )
                    for ( int c2 = 0; c2 <= c1; ++c2 )
                    {
                        *fieldsMap[filename].second[c1][c2] =
                            div( sum( scalarSpace, idv( func )( c1, c2 ) * meas() ),
                                 sum( scalarSpace, meas() ) ); // TODO optimize!!!
                        fieldsMap[filename].second[c1][c2]->setName( name );
                    }
            }
        }
        else
            CHECK( false ) << "invalid FunctionType";

        markModified();

        showMe( "ExportFieldSet::add" );
        toc( fmt::format( "ExportFieldSet::add functionspace element {}", name ),
             Environment::logVerbosityLevel() > 0 );
    }

    //! \brief Evaluate an expression on the whole associated mesh.
    template <typename ExprT>
        requires std::is_base_of_v<ExprBase, ExprT>
    void add( std::string const &name, ExprT const &expr,
              variant_representation_arg_type const &rep = "" )
    {
        this->add( name, name, expr, rep );
    }

    //! \brief Evaluate an expression on an explicit mesh range.
    template <typename ExprT, typename EltWrapperT = Range<mesh_type, MESH_ELEMENTS>>
        requires( std::is_base_of_v<ExprBase, ExprT> && is_filter_v<EltWrapperT> )
    void add( std::string const &name, ExprT const &expr, EltWrapperT const &rangElt,
              variant_representation_arg_type const &rep = "" )
    {
        this->add( name, name, expr, rangElt, rep );
    }

    //! \brief Evaluate an expression using distinct display and filename keys.
    template <typename ExprT>
        requires std::is_base_of_v<ExprBase, ExprT>
    void add( std::string const &name, std::string const &filename, ExprT const &expr,
              variant_representation_arg_type const &rep = "" )
    {

        CHECK( this->hasMesh() ) << "no mesh provided";
        this->add( name, filename, expr, elements( this->mesh() ), rep );
    }

    //! \brief Dispatch ranged expression conversion by representation.
    template <typename ExprT, typename EltWrapperT = Range<mesh_type, MESH_ELEMENTS>>
        requires( std::is_base_of_v<ExprBase, ExprT> && is_filter_v<EltWrapperT> )
    void add( std::string const &name, std::string const &filename, ExprT const &expr,
              EltWrapperT const &rangElt,
              variant_representation_arg_type const &requestedRepresentation = "" )
    {

        std::set<std::string> reps = representationType( requestedRepresentation, "nodal" );

        std::map<std::string, std::string> repToSuffix = { { "nodal", "_n" }, { "element", "_e" } };

        for ( std::string const &rep : reps )
        {
            std::string nameUsed =
                ( reps.size() > 1 ) ? sanitize( name + repToSuffix[rep] ) : sanitize( name );
            std::string fnameUsed = ( reps.size() > 1 ) ? filename + repToSuffix[rep] : filename;
            if ( rep == "element" )
                this->addExpr<false>( nameUsed, fnameUsed, expr, rangElt );
            else if ( rep == "nodal" )
                this->addExpr<true>( nameUsed, fnameUsed, expr, rangElt );
        }
    }

    //! \brief Evaluate each scalar/vector/tensor expression component into its own buffer.
    template <bool isNodal, typename ExprT, typename EltWrapperT = Range<mesh_type, MESH_ELEMENTS>>
    FEELPP_NO_EXPORT void addExpr( std::string const &name, std::string const &filename,
                                   ExprT const &expr, EltWrapperT const &rangeElt )
    {
        validateFieldName( name, filename );
        ( isNodal ? M_nodalNames : M_elementNames )[filename] = name;
        tic();
        auto scalarSpace = this->scalarFunctionSpace<isNodal>();
        auto &fieldsMap = this->fields<isNodal>();

        using ExprShapeType =
            typename ExprT::template evaluator_t<typename mesh_type::element_type>::shape;

        if constexpr ( ExprShapeType::is_scalar )
        {
            fieldsMap[filename].first = FunctionSpaceType::SCALAR;
            fieldsMap[filename].second.resize( 1 );
            fieldsMap[filename].second[0].resize( 1 );
            if ( !fieldsMap[filename].second[0][0] )
                fieldsMap[filename].second[0][0] = scalarSpace->elementPtr( name );
            fieldsMap[filename].second[0][0]->on( _range = rangeElt, _expr = expr, _close = true );
        }

        else if constexpr ( ExprShapeType::is_vectorial )
        {
            fieldsMap[filename].first = FunctionSpaceType::VECTORIAL;
            constexpr uint16_type nCompExpr =
                ( ExprShapeType::is_transposed ) ? ExprShapeType::N : ExprShapeType::M;
            fieldsMap[filename].second.resize( nCompExpr );
            for ( int c1 = 0; c1 < nCompExpr; ++c1 )
            {
                fieldsMap[filename].second[c1].resize( 1 );
                if ( !fieldsMap[filename].second[c1][0] )
                    fieldsMap[filename].second[c1][0] = scalarSpace->elementPtr( name );
                if constexpr ( ExprShapeType::is_transposed )
                    fieldsMap[filename].second[c1][0]->on( _range = rangeElt, _expr = expr( 0, c1 ),
                                                           _close = true ); // TODO : optimize
                else
                    fieldsMap[filename].second[c1][0]->on( _range = rangeElt, _expr = expr( c1, 0 ),
                                                           _close = true ); // TODO : optimize
            }
        }

        else if constexpr ( ExprShapeType::is_tensor2 )
        {
            fieldsMap[filename].first = FunctionSpaceType::TENSOR2;
            fieldsMap[filename].second.resize( ExprShapeType::M );
            for ( int c1 = 0; c1 < ExprShapeType::M; ++c1 )
            {
                fieldsMap[filename].second[c1].resize( ExprShapeType::N );
                for ( int c2 = 0; c2 < ExprShapeType::N; ++c2 )
                {
                    if ( !fieldsMap[filename].second[c1][c2] )
                        fieldsMap[filename].second[c1][c2] = scalarSpace->elementPtr( name );
                    fieldsMap[filename].second[c1][c2]->on( _range = rangeElt,
                                                            _expr = expr( c1, c2 ),
                                                            _close = true ); // TODO : optimize
                }
            }
        }
        else
        {
            CHECK( false ) << "expression shape not supported";
        }

        markModified();
        toc( fmt::format( "ExportFieldSet::add expression {}", name ),
             Environment::logVerbosityLevel() > 0 );
    }

    //@}

    //! \name  Methods
    //@{

    //! \return First nodal field schema and payload.
    nodal_const_iterator beginNodal() const { return M_nodal.begin(); }

    //! \return End of nodal fields.
    nodal_const_iterator endNodal() const { return M_nodal.end(); }

    //! \return Range of nodal fields.
    std::pair<nodal_const_iterator, nodal_const_iterator> nodal() const
    {
        return std::pair{ M_nodal.begin(), M_nodal.end() };
    }

    //! \return First element field schema and payload.
    element_const_iterator beginElement() const { return M_element.begin(); }

    //! \return End of element fields.
    element_const_iterator endElement() const { return M_element.end(); }

    //! \brief Emit residency-state diagnostics at verbose level.
    void showMe( std::string const &str ) const
    {

        DVLOG( 2 ) << str << " hasData() " << hasData() << "\n";

        DVLOG( 2 ) << str << " isInMemory() " << isInMemory() << "\n";
    }

    //! \brief Display name retained after component buffers are released.
    //! \param key Output filename key.
    //! \param nodal True for a nodal field, false for an element field.
    std::string const &fieldName( std::string const &key, bool nodal ) const
    {
        return ( nodal ? M_nodalNames : M_elementNames ).at( key );
    }

    //! \brief Release field values, retaining the lightweight field schema.
    void cleanup()
    {
        auto clear = []( auto &c ) { c.second.second.clear(); };
        std::for_each( M_nodal.begin(), M_nodal.end(), clear );
        std::for_each( M_element.begin(), M_element.end(), clear );
        M_isInMemory = false;
    }
    //@}

    //! \brief Borrow reserved dataset names for collision checks on temporal fields.
    void reserveNames( std::shared_ptr<std::set<std::string> const> const &names )
    {
        M_reservedFieldNames = names;
    }

    //! \return Conversion workspace shared by compatible field sets.
    std::shared_ptr<SpaceCache> const &spaceCache() const { return M_spaceCache; }

    //! \return All output keys currently present, including cleaned schemas.
    std::set<std::string> names() const
    {
        std::set<std::string> result;
        for ( auto const &[name, value] : M_scalar )
            result.insert( name );
        for ( auto const &[name, value] : M_nodal )
            result.insert( name );
        for ( auto const &[name, value] : M_element )
            result.insert( name );
        return result;
    }

    //! \brief Materialize or detach an immutable snapshot for a temporal-only writer.
    //! Component pointers are shared temporarily; field values are not recopied.
    //! Call only after checking that names do not collide with temporal fields.
    void materialize( ExportFieldSet const &snapshot, bool attach )
    {
        auto merge = [attach]( auto &destination, auto const &fields )
        {
            for ( auto const &[key, value] : fields )
                if ( attach )
                    destination.emplace( key, value );
                else
                    destination.erase( key );
        };
        merge( M_scalar, snapshot.M_scalar );
        merge( M_nodal, snapshot.M_nodal );
        merge( M_element, snapshot.M_element );
        merge( M_nodalNames, snapshot.M_nodalNames );
        merge( M_elementNames, snapshot.M_elementNames );
        if ( attach )
            markModified();
        else
            ++M_revision;
    }

  private:
    //! \brief Record a successful geometry/field registration.
    //! Owners compare revisions to detect changes after publishing their data.
    void markModified()
    {
        M_hasData = true;
        M_isInMemory = true;
        ++M_revision;
    }

    //! \return Mutable nodal payload map selected at compile time.
    template <bool isNodal>
    map_nodal_type &fields( typename std::enable_if<isNodal>::type * = nullptr )
    {
        return M_nodal;
    }

    //! \return Mutable element payload map selected at compile time.
    template <bool isNodal>
    map_element_type &fields( typename std::enable_if<!isNodal>::type * = nullptr )
    {
        return M_element;
    }

    //! \brief Reuse scalar export spaces within a time set or a detached dataset snapshot.
    //! Cache misses use Pch<N>'s automatic manager policy: compatible whole-mesh
    //! spaces can be shared across exporters when functionspace.manager.enable is true.
    //! Only the space is shared; each exported component owns its field values.
    template <bool isNodal, typename FunctionType = std::nullopt_t>
    FEELPP_NO_EXPORT scalar_p1_space_ptrtype
    scalarFunctionSpace( FunctionType const &func = std::nullopt,
                         typename std::enable_if<isNodal>::type * = nullptr )
    {
        auto &sharedSpace = M_spaceCache->M_p1;
        if ( !sharedSpace )
        {
            if constexpr ( is_functionspace_element_v<FunctionType> )
                if constexpr ( std::is_same_v<scalar_p1_space_type,
                                              typename FunctionType::functionspace_type> )
                    if ( func.mesh() == M_mesh &&
                         !func.functionSpace()->dof()->hasDofTableExtended() &&
                         support( func.functionSpace() )->isFullSupport() )
                        sharedSpace = func.functionSpace();
            if ( !sharedSpace )
                sharedSpace = Pch<N>( M_mesh.value() );
        }
        if ( M_mesh.value() == sharedSpace->mesh() )
            M_scalarP1 = sharedSpace;
        else if ( !M_scalarP1 || M_mesh.value() != M_scalarP1->mesh() )
            M_scalarP1 = Pch<N>( M_mesh.value() );
        return M_scalarP1;
    }

    //! \brief Element-space counterpart using Pdh<0>'s automatic manager policy.
    //! Detached snapshots retain their spaces without sharing mutable field values.
    template <bool isNodal, typename FunctionType = std::nullopt_t>
    FEELPP_NO_EXPORT scalar_p0_space_ptrtype
    scalarFunctionSpace( FunctionType const &func = std::nullopt,
                         typename std::enable_if<!isNodal>::type * = nullptr )
    {
        auto &sharedSpace = M_spaceCache->M_p0;
        if ( !sharedSpace )
        {
            if constexpr ( is_functionspace_element_v<FunctionType> )
                if constexpr ( std::is_same_v<scalar_p0_space_type,
                                              typename FunctionType::functionspace_type> )
                    if ( func.mesh() == M_mesh &&
                         !func.functionSpace()->dof()->hasDofTableExtended() &&
                         support( func.functionSpace() )->isFullSupport() )
                        sharedSpace = func.functionSpace();
            if ( !sharedSpace )
                sharedSpace = Pdh<0>( M_mesh.value() );
        }
        if ( M_mesh.value()->isSameMesh( sharedSpace->mesh() ) )
            M_scalarP0 = sharedSpace;
        else if ( !M_scalarP0 || !M_mesh.value()->isSameMesh( M_scalarP0->mesh() ) )
            M_scalarP0 = Pdh<0>( M_mesh.value() );
        return M_scalarP0;
    }

    //! \brief Reject names reserved by an exporter-owned dataset snapshot.
    //! The weak reference neither owns static values nor ties them to a time set.
    //! Exporter::save also validates direct time-set access before any output.
    void validateFieldName( std::string const &name, std::string const &key ) const
    {
        auto reserved = M_reservedFieldNames.lock();
        if ( !reserved )
            return;
        int invalid = reserved->count( name ) || reserved->count( key );
        int anyInvalid = 0;
        MPI_Allreduce( &invalid, &anyInvalid, 1, MPI_INT, MPI_MAX, mesh()->worldComm().comm() );
        if ( anyInvalid )
            throw std::invalid_argument( "a field is already registered on the exporter" );
    }

  public:
    //! \brief Normalize and validate nodal/element representation arguments.
    static std::set<std::string>
    representationType( variant_representation_arg_type const &requestedRepresentation,
                        std::string const &valueIfEmpty = "" )
    {
        std::set<std::string> reps;
        if ( auto repStringPtr = std::get_if<std::string>( &requestedRepresentation ) )
        {
            if ( !repStringPtr->empty() )
                reps.insert( *repStringPtr );
        }

        else if ( auto repSetPtr = std::get_if<std::set<std::string>>( &requestedRepresentation ) )
        {
            reps = *repSetPtr;
        }
        if ( reps.empty() && !valueIfEmpty.empty() )
            reps.insert( valueIfEmpty );

        for ( std::string const &rep : reps )
            if ( rep != "nodal" && rep != "element" )
                throw std::invalid_argument( "invalid representation: " + rep +
                                             "; use nodal or element" );

        return reps;
    }

  private:
    //! Weak name reservations; never owns another field set or a time set.
    std::weak_ptr<std::set<std::string> const> M_reservedFieldNames;
    //! Mesh defining the exported field ordering.
    std::optional<mesh_ptrtype> M_mesh;

    map_scalar_type M_scalar;   //!< Named per-case scalar values and constant flags.
    map_complex_type M_complex; //!< Legacy complex values and constant flags.
    map_nodal_type M_nodal;     //!< Nodal component buffers and shape descriptors.
    map_element_type M_element; //!< Element component buffers and shape descriptors.
    //! Persistent display names, independent of the component-vector lifetime.
    std::map<std::string, std::string> M_nodalNames;
    std::map<std::string, std::string> M_elementNames;

    bool M_hasData = false;       //!< Geometry/field metadata is present.
    bool M_isInMemory = false;    //!< Component buffers are marked resident.
    std::uint64_t M_revision = 0; //!< Content revision, not a temporal index.

    scalar_p0_space_ptrtype M_scalarP0; //!< Element space for this mesh.
    scalar_p1_space_ptrtype M_scalarP1; //!< Nodal space for this mesh.
    //! Shared conversion workspace; it has no temporal or exporter ownership.
    std::shared_ptr<SpaceCache> M_spaceCache;
};

} // namespace detail
} // namespace Feel
#endif // FEELPP_DISCR_EXPORTFIELDSET_HPP
