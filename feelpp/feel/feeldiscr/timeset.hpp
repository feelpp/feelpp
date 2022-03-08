/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
/**
   \file timeSet.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-12
 */
#ifndef FEELPP_DISCR_TIMESET_HPP
#define FEELPP_DISCR_TIMESET_HPP 1

#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <variant>
#include <optional>

#include <boost/operators.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/split_member.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>

// #include <boost/numeric/ublas/vector_serialize.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/feelcomplex.hpp>
#include <feel/feelcore/context.hpp>

#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/lagrange.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/interpolate.hpp>

#include <feel/feelmesh/filters.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/elementdiv.hpp>

#define TS_INITIAL_INDEX 1

namespace Feel
{

enum
{
    STEP_NEW       = ( 1<<0 ),
    STEP_HAS_DATA  = ( 1<<1 ),
    STEP_ON_DISK   = ( 1<<2 ),
    STEP_IN_MEMORY = ( 1<<3 ),
    STEP_IGNORED   = ( 1<<4 ),
    STEP_OVERWRITE = ( 1<<10 )
};
template<typename A0,typename A1,typename A2,typename A3,typename A4> class FunctionSpace;
namespace detail
{
/**
 * \class TimeSet
 * \ingroup SpaceTime
 * \brief data TimeSet
 *
 * \tparam MeshType     Mesh type
 * \tparam N            Mesh geometrical order
 *
 * \author Christophe Prud'homme
 */
template<typename MeshType, int N = 1>
class TimeSet
{
public:

    /** @name Subclasses
     */
    //@{
    typedef MeshType mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    using size_type = typename mesh_type::size_type;

    typedef Pdh_type<MeshType,0> scalar_p0_space_type;
    typedef Pch_type<MeshType,N/*1*/> scalar_p1_space_type;
    typedef std::shared_ptr<scalar_p0_space_type> scalar_p0_space_ptrtype;
    typedef std::shared_ptr<scalar_p1_space_type> scalar_p1_space_ptrtype;


    typedef typename scalar_p0_space_type::element_type element_scalar_type;
    typedef typename scalar_p1_space_type::element_type nodal_scalar_type;


    /**
       \brief a step in a time set
    */
    class Step
        :
        boost::equality_comparable<Step>,
        boost::less_than_comparable<Step>
    {
    public:

        /**
         */
        //@{
        typedef Step step_type;
        typedef MeshType mesh_type;
        typedef std::shared_ptr<mesh_type> mesh_ptrtype;
        using size_type = typename mesh_type::size_type;

        typedef Pdh_type<MeshType,0> scalar_p0_space_type;
        typedef Pch_type<MeshType,N/*1*/> scalar_p1_space_type;
        typedef std::shared_ptr<scalar_p0_space_type> scalar_p0_space_ptrtype;
        typedef std::shared_ptr<scalar_p1_space_type> scalar_p1_space_ptrtype;

        typedef typename scalar_p0_space_type::element_type element_scalar_type;
        typedef typename scalar_p1_space_type::element_type nodal_scalar_type;

        typedef std::shared_ptr<nodal_scalar_type> nodal_scalar_ptrtype;
        typedef std::shared_ptr<element_scalar_type> element_scalar_ptrtype;
        typedef std::pair<FunctionSpaceType,std::vector<std::vector<nodal_scalar_ptrtype>>> nodal_field_type;
        typedef std::pair<FunctionSpaceType,std::vector<std::vector<element_scalar_ptrtype>>> element_field_type;

        typedef std::map<std::string, std::pair<scalar_type,bool> > map_scalar_type;
        typedef std::map<std::string, std::pair<complex_type,bool> > map_complex_type;
        typedef std::map<std::string, nodal_field_type> map_nodal_type;
        typedef std::map<std::string, element_field_type> map_element_type;

        typedef typename map_scalar_type::iterator scalar_iterator;
        typedef typename map_scalar_type::const_iterator scalar_const_iterator;
        typedef typename map_complex_type::iterator complex_iterator;
        typedef typename map_complex_type::const_iterator complex_const_iterator;
        typedef typename map_nodal_type::iterator nodal_iterator;
        typedef typename map_nodal_type::const_iterator nodal_const_iterator;
        typedef typename map_element_type::iterator element_iterator;
        typedef typename map_element_type::const_iterator element_const_iterator;

        typedef std::variant<std::string, std::set<std::string> > variant_representation_arg_type;

        //@}

        /** @name Constructors, destructor
         */
        //@{

        /**
           destructor
        */
        FEELPP_NO_EXPORT ~Step();

        //@}

        /** @name Accessors
         */
        //@{

        /**
         *
         * @return
         */
        bool isNew() const
        {
            return M_state.test( STEP_NEW );
        }

        /**
         *
         * @return
         */
        bool hasData() const
        {
            return M_state.test( STEP_HAS_DATA );
        }

        /**
         *
         * @return
         */
        bool isOnDisk() const
        {
            return M_state.test( STEP_ON_DISK );
        }

        /**
         *
         * @return
         */
        bool isInMemory() const
        {
            return M_state.test( STEP_IN_MEMORY );
        }

        bool isIgnored() const
        {
            return M_state.test( STEP_IGNORED );
        }

        /**
         * get the current state of the step
         * @return the size encoded into a size_type
         */
        size_type state() const
        {
            return M_state.context();
        }

        /**
           \return the time associated with the Step
        */
        Real time() const
        {
            return M_time;
        }

        /**
           \return the index of the time set
        */
        size_type index() const
        {
            return M_index;
        }

        //! \return the active index of the time set
        size_type activeIndex() const { return M_activeIndex; }

        /**
           \return true if the mesh is available
        */
        bool hasMesh() const
        {
            return M_mesh != std::nullopt;
        }

        /**
           \return a mesh
        */
        mesh_ptrtype mesh()
        {
            return M_mesh.value();
        }

        /**
         * @return the begin iterator for scalars
         */
        scalar_const_iterator beginScalar() const { return M_scalar.begin(); }

        /**
         * @return the end iterator for scalars
         */
        scalar_const_iterator endScalar() const { return M_scalar.end(); }

        /**
         * get the scalar with name n
         * @param __n name of the nodal scalar field
         * @return the scalar value
         */
        scalar_type scalar( std::string const& __n ) const
        {
            if ( M_scalar.find( sanitize(__n) ) == M_scalar.end() )
            {
                std::ostringstream __err;
                __err << "invalid scalar value name " << sanitize(__n);
                throw std::logic_error( __err.str() );
            }

            return M_scalar.find( sanitize(__n) )->second.first;
        }


        /**
         * get the nodal field with name n
         * @param __n name of the nodal field
         * @return the nodal field
         */
        nodal_field_type const& nodal( std::string const& __n ) const
        {
            auto itFindNodalField =  M_nodal.find( sanitize(__n) );
            CHECK( itFindNodalField != M_nodal.end() ) << "invalid nodal field name " << sanitize(__n);
            return itFindNodalField->second;
        }
        /**
         * get the element field with name n
         * @param __n name of the element field
         * @return the element field
         */
        element_field_type const& element( std::string const& __n ) const
        {
            auto itFindElementField =  M_element.find( sanitize(__n) );
            CHECK( itFindElementField != M_element.end() ) << "invalid nodal field name " << sanitize(__n);
            return itFindElementField->second;
        }

        //@}

        /** @name  Mutators
         */
        //@{

        void setState( size_type __st )
        {
            M_state.set( __st );
            showMe( "Step::setState" );
        }
        void setMesh( mesh_ptrtype const& __m )
        {
            DVLOG(2) << "[TimeSet::setMesh] setMesh start\n";
            M_mesh = __m;

            M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            M_state.clear( STEP_ON_DISK );
            DVLOG(2) << "[TimeSet::setMesh] setMesh start\n";
        }

        FEELPP_DEPRECATED void addScalar( std::string const& name, scalar_type const& __s, bool cst = false )
        {
            if ( this->isIgnored() )
                return;
            M_scalar[sanitize(name)] =  std::make_pair( __s, cst );
            M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            M_state.clear( STEP_ON_DISK );
        }
        template<typename T>
        void add( std::string const& name, T const& __s, bool cst = false,
                  typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr )
            {
                if ( this->isIgnored() )
                    return;
                M_scalar[sanitize(name)] =  std::make_pair( __s, cst );
                M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
                M_state.clear( STEP_ON_DISK );
            }

        void addComplex( std::string const& name, complex_type const& __s, bool cst = false )
        {
            if ( this->isIgnored() )
                return;
            M_complex[sanitize(name)] =  std::pair{__s,cst};
            M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            M_state.clear( STEP_ON_DISK );
        }

        /**
         * @brief add regions to timeset
         * @details some regions can be automatically generated such as the
         * process id map
         *
         * @param  prefix prefix string for the region
         */
        void
        addRegions( std::string const& prefix = "" )
        {
            if ( this->isIgnored() )
                return;
            this->addRegions( prefix,prefix );
        }
        void
        addRegions( std::string const& prefix, std::string const& prefixfname )
        {
            if ( this->isIgnored() )
                return;

            VLOG(1) << "[timeset] Adding regions...\n";
#if 0
            worldcomm_ptr_t meshComm = (M_mesh.value()->worldComm().numberOfSubWorlds() > 1)?
                M_mesh.value()->worldComm().subWorld(M_mesh.value()->worldComm().numberOfSubWorlds()) :
                M_mesh.value()->worldCommPtr();
            this->updateScalarP0( makeWorldsComm(1,meshComm) );
#endif

            auto scalarElementSpace = this->scalarFunctionSpace<false>();

            VLOG(1) << "[timeset] adding pid...\n";
            this->add( prefixvm(prefix,"pid"), prefixvm(prefixfname,"pid"),  regionProcess( scalarElementSpace ) );
        }

        template<typename FunctionType>
        void add( std::initializer_list<std::string>  __n, FunctionType const& func,
                  typename std::enable_if<is_functionspace_element_v<decay_type<FunctionType>>>::type* = nullptr )
        {
            if ( this->isIgnored() )
                return;

            std::vector<std::string> str( __n );
            add( str, func );
        }

        template<typename FunctionType>
        void add( std::variant<std::string,std::vector<std::string>> const& __n, FunctionType const& func, variant_representation_arg_type const& reps = "",
                  typename std::enable_if<is_functionspace_element_v<decay_type<FunctionType>>>::type* = nullptr )
        {
            if ( this->isIgnored() )
                return;

            std::vector<std::string> names;
            if( auto  nameStringPtr = std::get_if<std::string>(&__n))
            {
                if ( !nameStringPtr->empty() )
                    names.push_back( *nameStringPtr );
            }
            else if ( auto nameVectorPtr = std::get_if<std::vector<std::string>>(&__n))
            {
                names = *nameVectorPtr;
            }
            CHECK( !names.empty() ) << "no field name given";

            constexpr int nSpaces = decay_type<FunctionType>::functionspace_type::nSpaces;
            if constexpr( nSpaces > 1 )
            {
                if ( names.size() == 1 )
                {
                    std::string nameGiven = names[0];
                    names.resize( nSpaces );
                    for(int i=0; i<nSpaces; i++)
                        names[i] = (boost::format("%1%_%2%")%nameGiven %i).str();
                }

                hana::for_each( hana::make_range( hana::int_c<0>, hana::int_c<nSpaces> ), [this,&names,&func,&reps]( auto const& e )
                                 {
                                     constexpr int subId = std::decay_t<decltype(e)>::value;
                                     auto const& subFunc = func.template element<subId>();
                                     this->add( names[subId], names[subId], subFunc, reps );
                                 } );
            }
            else
            {
                add( names[0], names[0], func, reps );
            }
        }

        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, variant_representation_arg_type const& _reps = "",
                  typename std::enable_if<is_functionspace_element_v<decay_type<FunctionType>>>::type* = nullptr )
        {
            if ( this->isIgnored() )
                return;

            using FunctionDecayType = decay_type<FunctionType>;
            constexpr bool funcIsNodal = ( FunctionDecayType::is_continuous || FunctionDecayType::functionspace_type::continuity_type::is_discontinuous_locally ) && (!FunctionDecayType::is_hcurl_conforming) && (!FunctionDecayType::is_hdiv_conforming);

            std::string defaultRepr = funcIsNodal? "nodal" : "element";
#if 0
            // if Lagrange discontinuous with order > 0, put nodal representation by default;
            if constexpr ( is_lagrange_polynomialset_v<typename FunctionDecayType::functionspace_type::fe_type> && !funcIsNodal )
            {
                if (  FunctionDecayType::functionspace_type::fe_type::nOrder > 0 )
                    defaultRepr = "nodal";
            }
#endif
            std::set<std::string> reps = representationType( _reps , defaultRepr );

            std::map<std::string,std::string> repToSuffix = { { "nodal", "_n" }, { "element", "_e" } };

            for ( std::string const& rep : reps )
            {
                std::string nameUsed = (reps.size() > 1)? sanitize(__n+repToSuffix[rep]) : sanitize(__n);
                std::string fnameUsed =  (reps.size() > 1)? __fname+repToSuffix[rep] : __fname;
                if ( rep == "element" )
                    add<false,false>( nameUsed,fnameUsed,unwrap_ptr(func) );
                else if ( rep == "nodal" )
                {
                    if constexpr ( funcIsNodal )
                        add<true,false>( nameUsed,fnameUsed,unwrap_ptr(func) );
                    else
                        add<true,true>( nameUsed,fnameUsed,unwrap_ptr(func) );
                }
            }
        }

        template<bool IsNodal,bool IsElementToNodal,typename FunctionType>
        FEELPP_NO_EXPORT void add( std::string const& __n, std::string const& __fname, FunctionType const& func )
        {
            if ( !func.worldComm().isActive() ) return;

            std::string reprType = IsNodal? "nodal":"element";
            tic();
            auto scalarSpace = this->scalarFunctionSpace<IsNodal>( func );
            toc( (boost::format("Timeset::add get scalar space %1%")%reprType).str(),FLAGS_v>0);

            tic();
            auto & fieldsMap = this->fields<IsNodal>();
            //std::vector<ComponentType> mapIndicesToComponent = { ComponentType::X, ComponentType::Y, ComponentType::Z };

            if constexpr ( FunctionType::is_scalar )
                {
                    fieldsMap[ __fname].first = FunctionSpaceType::SCALAR;
                    fieldsMap[ __fname].second.resize( 1, { scalarSpace->elementPtr( __n, func.description() ) } );
                    if constexpr ( !IsElementToNodal )
                        {
                            interpolate( scalarSpace, func, *fieldsMap[__fname].second[0][0] );
                        }
                    else
                    {
                        *fieldsMap[__fname].second[0][0] =  div( sum( scalarSpace, idv(func)*meas() ), sum( scalarSpace, meas() ) ); // TODO optimize!!!
                        fieldsMap[__fname].second[0][0]->setName( __n );
                    }
                }
            else if constexpr ( FunctionType::is_vectorial )
                {
                    fieldsMap[ __fname].first = FunctionSpaceType::VECTORIAL;
                    fieldsMap[ __fname].second.resize(  FunctionType::nComponents );
                    for ( int c1 = 0; c1 < FunctionType::nComponents ;++c1 )
                        fieldsMap[ __fname].second[c1] = { scalarSpace->elementPtr( __n, func.description() ) };
                    if constexpr ( !IsElementToNodal )
                         {
                             interpolate( scalarSpace, func, fieldsMap[__fname].second );
                         }
                    else
                    {
                        for ( int c1 = 0; c1 < FunctionType::nComponents ;++c1 )
                        {
                            *fieldsMap[__fname].second[c1][0] =  div( sum( scalarSpace, idv(func)(c1,0)*meas() ), sum( scalarSpace, meas() ) ); // TODO optimize!!!
                            fieldsMap[__fname].second[c1][0]->setName( __n );
                        }
                    }
                }
            else if constexpr ( FunctionType::is_tensor2 )
                {
                    fieldsMap[ __fname].first = FunctionSpaceType::TENSOR2;
                    fieldsMap[ __fname].second.resize(  FunctionType::nComponents1 );
                    for ( int c1 = 0; c1 < FunctionType::nComponents1 ;++c1 )
                    {
                        fieldsMap[ __fname].second[c1].resize( FunctionType::nComponents2 );
                        for ( int c2 = 0; c2 < FunctionType::nComponents2 ;++c2 )
                            fieldsMap[ __fname].second[c1][c2] = scalarSpace->elementPtr( __n, func.description() );
                    }
                    if constexpr ( !IsElementToNodal )
                    {
                        interpolate( scalarSpace, func, fieldsMap[__fname].second );
                    }
                    else
                    {
                        for ( int c1 = 0; c1 < FunctionType::nComponents1 ;++c1 )
                            for ( int c2 = 0; c2 < FunctionType::nComponents2 ;++c2 )
                            {
                                *fieldsMap[__fname].second[c1][c2] =  div( sum( scalarSpace, idv(func)(c1,c2)*meas() ), sum( scalarSpace, meas() ) ); // TODO optimize!!!
                                fieldsMap[__fname].second[c1][c2]->setName( __n );
                            }
                    }

                }
            else if constexpr ( FunctionType::is_tensor2symm )
                {
                    fieldsMap[ __fname].first = FunctionSpaceType::TENSOR2_SYMM;
                    fieldsMap[ __fname].second.resize( FunctionType::nComponents1 );
                    for ( int c1 = 0; c1 < FunctionType::nComponents1 ;++c1 )
                    {
                        fieldsMap[ __fname].second[c1].resize( c1+1 );
                        for ( int c2 = 0; c2 <= c1 ;++c2 )
                            fieldsMap[ __fname].second[c1][c2] = scalarSpace->elementPtr( __n, func.description() );
                    }
                    if constexpr ( !IsElementToNodal )
                    {
                        interpolate( scalarSpace, func, fieldsMap[__fname].second );
                    }
                    else
                    {
                        for ( int c1 = 0; c1 < FunctionType::nComponents1 ;++c1 )
                            for ( int c2 = 0; c2 <= c1 ;++c2 )
                            {
                                *fieldsMap[__fname].second[c1][c2] =  div( sum( scalarSpace, idv(func)(c1,c2)*meas() ), sum( scalarSpace, meas() ) ); // TODO optimize!!!
                                fieldsMap[__fname].second[c1][c2]->setName( __n );
                            }
                    }

                }
            else
                CHECK( false ) << "invalid FunctionType";

            M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            M_state.clear( STEP_ON_DISK );

            showMe( "Step::add" );
            toc((boost::format("Timeset::add functionspace element %1%")%__n).str(),FLAGS_v>0);
        }



        template<typename ExprT>
        void add( std::string const& __n, ExprT const& expr, variant_representation_arg_type const& rep = "",
                  typename std::enable_if_t<std::is_base_of_v<ExprBase,ExprT> >* = nullptr )
            {
                this->add( __n, __n, expr, rep );
            }
        template<typename ExprT, typename EltWrapperT = elements_reference_wrapper_t<mesh_type>>
        void add( std::string const& __n, ExprT const& expr,  EltWrapperT const& rangElt, variant_representation_arg_type const& rep = "",
                  typename std::enable_if_t<std::is_base_of_v<ExprBase,ExprT> >* = nullptr )
            {
                this->add( __n, __n, expr, rangElt, rep );
            }
        template<typename ExprT>
        void add( std::string const& __n, std::string const& __fname, ExprT const& expr, variant_representation_arg_type const& rep = "",
                  typename std::enable_if_t<std::is_base_of_v<ExprBase,ExprT> >* = nullptr )
            {
                if ( this->isIgnored() )
                    return;

                CHECK( this->hasMesh() ) << "no mesh provided";
                this->add( __n, __fname, expr, elements(this->mesh()), rep );
            }

        template<typename ExprT, typename EltWrapperT = elements_reference_wrapper_t<mesh_type>>
        void add( std::string const& __n, std::string const& __fname, ExprT const& expr, EltWrapperT const& rangElt, variant_representation_arg_type const& _rep = "",
                  typename std::enable_if_t<std::is_base_of_v<ExprBase,ExprT> >* = nullptr )
            {
                if ( this->isIgnored() )
                    return;

                std::set<std::string> reps = representationType( _rep , "nodal" );

                std::map<std::string,std::string> repToSuffix = { { "nodal", "_n" }, { "element", "_e" } };

                for ( std::string const& rep : reps )
                {
                    std::string nameUsed = (reps.size() > 1)? sanitize(__n+repToSuffix[rep]) : sanitize(__n);
                    std::string fnameUsed =  (reps.size() > 1)? __fname+repToSuffix[rep] : __fname;
                    if ( rep == "element" )
                        this->addExpr<false>( nameUsed,fnameUsed,expr,rangElt );
                    else if ( rep == "nodal" )
                        this->addExpr<true>( nameUsed,fnameUsed,expr,rangElt );
                }
            }
        template<bool IsNodal,typename ExprT, typename EltWrapperT = elements_reference_wrapper_t<mesh_type>>
        FEELPP_NO_EXPORT void addExpr( std::string const& __n, std::string const& __fname, ExprT const& expr, EltWrapperT const& rangeElt )
            {
                tic();
                auto scalarSpace = this->scalarFunctionSpace<IsNodal>();
                auto & fieldsMap = this->fields<IsNodal>();

                using ExprShapeType = typename ExprT::template evaluator_t<typename mesh_type::element_type>::shape;

                if constexpr ( ExprShapeType::is_scalar )
                    {
                        fieldsMap[ __fname].first = FunctionSpaceType::SCALAR;
                        //fieldsMap[ __fname].second.resize( 1, { scalarSpace->elementPtr( __n ) } );
                        fieldsMap[ __fname].second.resize( 1 );
                        fieldsMap[ __fname].second[0].resize( 1 );
                        if ( !fieldsMap[__fname].second[0][0] )
                            fieldsMap[__fname].second[0][0] = scalarSpace->elementPtr( __n );
                        fieldsMap[__fname].second[0][0]->on(_range=rangeElt,_expr=expr);
                    }
                else if constexpr ( ExprShapeType::is_vectorial )
                {
                    fieldsMap[ __fname].first = FunctionSpaceType::VECTORIAL;
                    constexpr uint16_type nCompExpr = ( ExprShapeType::is_transposed )? ExprShapeType::N : ExprShapeType::M;
                    fieldsMap[ __fname].second.resize( nCompExpr );
                    for ( int c1 = 0; c1 < nCompExpr ;++c1 )
                    {
                        fieldsMap[ __fname].second[c1].resize( 1 );
                        if ( !fieldsMap[ __fname].second[c1][0] )
                            fieldsMap[ __fname].second[c1][0] = scalarSpace->elementPtr( __n );
                        if constexpr ( ExprShapeType::is_transposed )
                            fieldsMap[__fname].second[c1][0]->on(_range=rangeElt,_expr=expr(0,c1)); // TODO : optimize
                        else
                            fieldsMap[__fname].second[c1][0]->on(_range=rangeElt,_expr=expr(c1,0) ); // TODO : optimize
                    }

                }
                else if constexpr ( ExprShapeType::is_tensor2 )
                {
                    fieldsMap[ __fname].first = FunctionSpaceType::TENSOR2;
                    fieldsMap[ __fname].second.resize( ExprShapeType::M );
                    for ( int c1 = 0; c1 < ExprShapeType::M ;++c1 )
                    {
                        fieldsMap[ __fname].second[c1].resize( ExprShapeType::N );
                        for ( int c2 = 0; c2 < ExprShapeType::N ;++c2 )
                        {
                            if ( !fieldsMap[ __fname].second[c1][c2] )
                                fieldsMap[ __fname].second[c1][c2] = scalarSpace->elementPtr( __n );
                            fieldsMap[__fname].second[c1][0]->on(_range=rangeElt,_expr=expr(c1,c2));  // TODO : optimize
                        }
                    }
                }
                else
                {
                    CHECK( false ) << "expression shape not supported";
                }

                M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
                M_state.clear( STEP_ON_DISK );
                toc((boost::format("Timeset::add expression %1%")%__n).str(),FLAGS_v>0);
            }

        //@}

        /** @name  Methods
         */
        //@{

        nodal_const_iterator beginNodal() const
        {
            return M_nodal.begin();
        }
        nodal_const_iterator endNodal() const
        {
            return M_nodal.end();
        }
        std::pair<nodal_const_iterator,nodal_const_iterator> nodal() const
        {
            return std::pair{M_nodal.begin(),M_nodal.end()};
        }
        element_const_iterator beginElement() const
        {
            return M_element.begin();
        }
        element_const_iterator endElement() const
        {
            return M_element.end();
        }
        bool operator==( Step const& __s ) const
        {
            return this->index() == __s.index();
        }
        bool operator<( Step const& __s ) const
        {
            return this->index() < __s.index();
        }

        void load()
        {
            this->setState( STEP_IN_MEMORY );
        }
        void showMe( std::string const& str ) const
        {
            DVLOG(2) << str << " index :  " << M_index << "\n";
            DVLOG(2) << str << " time :  " << M_time << "\n";
            DVLOG(2) << str << " isNew() " << isNew() << "\n";
            DVLOG(2) << str << " hasData() " << hasData() << "\n";
            DVLOG(2) << str << " isOnDisk() " << isOnDisk() << "\n";
            DVLOG(2) << str << " isInMemory() " << isInMemory() << "\n";
        }

        void cleanup()
            {
                auto clear = []( auto& c ) {
                    c.second.second.clear();
                    //for ( auto & c1 : c.second.second )
                    //    for ( auto & c2 : c1 )
                    //       c2->clear(); //????
                };
                std::for_each( M_nodal.begin(),
                               M_nodal.end(),
                               clear );
                std::for_each( M_element.begin(),
                               M_element.end(),
                               clear );
                M_state.clear( STEP_IN_MEMORY );
            }
        //@}


    private:

        /** TimeSet is a good friend */
        friend class TimeSet;

        /** @name Constructors, destructor
         */
        //@{

        /**
           only TimeSet can use the default constructor
        */
        FEELPP_NO_EXPORT Step();

        /**
           construct a step with index \c index at a specified time \c time

           \param time time associated with the Step
           \param index index of the step
           * \param state the state of the step
           */
        FEELPP_NO_EXPORT Step( TimeSet* ts, Real time, size_type index, size_type activeIndex, size_type __state = STEP_NEW|STEP_OVERWRITE );


        /**
           only TimeSet can use the copy constructor
        */
        FEELPP_NO_EXPORT Step( Step const & ) = default;



        //@}

        /** @name Operator overloads
         */
        //@{

        FEELPP_NO_EXPORT Step& operator=( Step const& ) = default;

        //@}

        template <bool IsNodal>
        map_nodal_type & fields( typename std::enable_if<IsNodal>::type* = nullptr ) { return M_nodal; }
        template <bool IsNodal>
        map_element_type & fields( typename std::enable_if<!IsNodal>::type* = nullptr ) { return M_element; }

        template<bool IsNodal,typename FunctionType =  std::nullopt_t>
        FEELPP_NO_EXPORT scalar_p1_space_ptrtype scalarFunctionSpace( FunctionType const& func = std::nullopt, typename std::enable_if<IsNodal>::type* = nullptr )
            {
                if ( !M_ts->M_scalar_p1 )
                {
                    /*if ( M_ts->M_vector_p1 && M_ts->M_vector_p1->hasCompSpace() )
                     {
                     M_ts->M_scalar_p1 = std::make_shared<scalar_p1_space_type>();
                     M_ts->M_scalar_p1->shallowCopy( M_ts->M_vector_p1->compSpace() );
                     }*/
                    if constexpr ( is_functionspace_element_v<FunctionType> )
                        {
                            if constexpr (std::is_same_v<scalar_p1_space_type,typename FunctionType::functionspace_type> )
                                {
                                    if ( ( func.mesh() == M_mesh ) && !func.functionSpace()->extendedDofTable() && support( func.functionSpace() )->isFullSupport() )
                                        M_ts->M_scalar_p1 = func.functionSpace();
                                }
                        }
                    if ( !M_ts->M_scalar_p1 )
                        M_ts->M_scalar_p1 = scalar_p1_space_type::New(_mesh=M_mesh.value() );
                    M_scalar_p1 = M_ts->M_scalar_p1;
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p1 created\n";
                }
                else if ( M_mesh.value() == M_ts->M_scalar_p1->mesh() )
                {
                    M_scalar_p1 = M_ts->M_scalar_p1;
                }

                if ( M_mesh.value() != M_ts->M_scalar_p1->mesh() && !M_scalar_p1 )
                {
                    M_scalar_p1 = scalar_p1_space_type::New(_mesh=M_mesh.value() );
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p1 created\n";
                }
                return M_scalar_p1;
            }
        template<bool IsNodal,typename FunctionType = std::nullopt_t>
        FEELPP_NO_EXPORT scalar_p0_space_ptrtype scalarFunctionSpace( FunctionType const& func = std::nullopt, typename std::enable_if<!IsNodal>::type* = nullptr )
        {
            if ( !M_ts->M_scalar_p0 )
            {
                if constexpr ( is_functionspace_element_v<FunctionType> )
                    {
                        if constexpr ( std::is_same_v<scalar_p0_space_type,typename FunctionType::functionspace_type> )
                        {
                            if ( ( func.mesh() == M_mesh ) && !func.functionSpace()->extendedDofTable() && support( func.functionSpace() )->isFullSupport() )
                                M_ts->M_scalar_p0 = func.functionSpace();
                        }
                    }
                if ( !M_ts->M_scalar_p0 )
                    M_ts->M_scalar_p0 = scalar_p0_space_type::New(_mesh=M_mesh.value() );
                M_scalar_p0 = M_ts->M_scalar_p0;
                DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p0 created\n";
            }
            else if ( M_mesh.value()->isSameMesh( M_ts->M_scalar_p0->mesh() ) )
            {
                M_scalar_p0 = M_ts->M_scalar_p0;
            }

            if ( !M_mesh.value()->isSameMesh( M_ts->M_scalar_p0->mesh() ) && !M_scalar_p0 )
            {
                M_scalar_p0 = scalar_p0_space_type::New(_mesh=M_mesh.value() );
                DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p0 created\n";
            }
            return M_scalar_p0;
        }


        friend class boost::serialization::access;

        template<class Archive>
        FEELPP_NO_EXPORT void save( Archive & ar, const unsigned int /*version*/ ) const
        {
#if 0
            size_type s;

            // values
            s = M_scalar.size();
            ar & boost::serialization::make_nvp( "map_scalar_size", s );

            scalar_const_iterator __its= M_scalar.begin();
            scalar_const_iterator __ens= M_scalar.end();

            for ( ; __its!= __ens; ++__its )
            {
                ar & boost::serialization::make_nvp( "scalar", *__its );
            }

            s = M_complex.size();
            ar & boost::serialization::make_nvp( "map_complex_size", s );

            complex_const_iterator __itc = M_complex.begin();
            complex_const_iterator __enc = M_complex.end();

            for ( ; __itc!= __enc; ++__itc )
            {
                ar & boost::serialization::make_nvp( "complex", ( std::string const& )__itc->first );
                ar & boost::serialization::make_nvp( "real", ( scalar_type const& )__itc->second.first.real() );
                ar & boost::serialization::make_nvp( "imaginary", ( scalar_type const& )__itc->second.first.imag() );
            }

            s = M_nodal_scalar.size();
            ar & boost::serialization::make_nvp( "map_nodal_scalar_size", s );
            DVLOG(2) << "(saving) serialized size of nodal scalar (" << s << ")\n";
            nodal_scalar_const_iterator __itscalar = M_nodal_scalar.begin();
            nodal_scalar_const_iterator __enscalar = M_nodal_scalar.end();

            for ( ; __itscalar != __enscalar; ++__itscalar )
            {
                ar & ( nodal_scalar_type const& ) __itscalar->second;
            }

            s = M_nodal_vector.size();
            ar & boost::serialization::make_nvp( "map_nodal_vector_size", s );

            nodal_vector_const_iterator __itvector= M_nodal_vector.begin();
            nodal_vector_const_iterator __envector= M_nodal_vector.end();

            for ( ; __itvector!= __envector; ++__itvector )
            {
                ar & ( nodal_vector_type const& ) __itvector->second;
            }

            s = M_nodal_tensor2.size();
            ar & boost::serialization::make_nvp( "map_nodal_tensor2_size", s );

            nodal_tensor2_const_iterator __ittensor2= M_nodal_tensor2.begin();
            nodal_tensor2_const_iterator __entensor2= M_nodal_tensor2.end();

            for ( ; __ittensor2!= __entensor2; ++__ittensor2 )
            {
                ar & ( nodal_tensor2_type const& ) __ittensor2->second;
            }

            s = M_nodal_tensor2symm.size();
            ar & boost::serialization::make_nvp( "map_nodal_tensor2symm_size", s );

            nodal_tensor2symm_const_iterator __ittensor2symm= M_nodal_tensor2symm.begin();
            nodal_tensor2symm_const_iterator __entensor2symm= M_nodal_tensor2symm.end();

            for ( ; __ittensor2symm!= __entensor2symm; ++__ittensor2symm )
            {
                ar & ( nodal_tensor2symm_type const& ) __ittensor2symm->second;
            }

            s = M_element_scalar.size();
            ar & boost::serialization::make_nvp( "map_element_scalar_size", s );
            DVLOG(2) << "(saving) serialized size of element scalar (" << s << ")\n";
            element_scalar_const_iterator __eitscalar = M_element_scalar.begin();
            element_scalar_const_iterator __eenscalar = M_element_scalar.end();

            for ( ; __eitscalar != __eenscalar; ++__eitscalar )
            {
                ar & ( element_scalar_type const& ) __eitscalar->second;
            }

            s = M_element_vector.size();
            ar & boost::serialization::make_nvp( "map_element_vector_size", s );

            element_vector_const_iterator __eitvector= M_element_vector.begin();
            element_vector_const_iterator __eenvector= M_element_vector.end();

            for ( ; __eitvector!= __eenvector; ++__eitvector )
            {
                ar & ( element_vector_type const& ) __eitvector->second;
            }

            s = M_element_tensor2.size();
            ar & boost::serialization::make_nvp( "map_element_tensor2_size", s );

            element_tensor2_const_iterator __eittensor2= M_element_tensor2.begin();
            element_tensor2_const_iterator __eentensor2= M_element_tensor2.end();

            for ( ; __eittensor2!= __eentensor2; ++__eittensor2 )
            {
                ar & ( element_tensor2_type const& ) __eittensor2->second;
            }

            s = M_element_tensor2symm.size();
            ar & boost::serialization::make_nvp( "map_element_tensor2symm_size", s );

            element_tensor2symm_const_iterator __eittensor2symm= M_element_tensor2symm.begin();
            element_tensor2symm_const_iterator __eentensor2symm= M_element_tensor2symm.end();

            for ( ; __eittensor2symm!= __eentensor2symm; ++__eittensor2symm )
            {
                ar & ( element_tensor2symm_type const& ) __eittensor2symm->second;
            }
#endif
        }

        template<class Archive>
        FEELPP_NO_EXPORT void load( Archive & ar, const unsigned int /*version*/ )
        {
#if 0
            // loading
            size_type s( 0 );

            ar & boost::serialization::make_nvp( "map_scalar_size", s );
            DVLOG(2) << "(loading) serialized size of  scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                std::pair<std::string, std::pair<scalar_type,bool> > v;
                ar & v;

                DVLOG(2) << "(loading) dserialized scalar  " << v.first
                              << " with value " << v.second.first << "\n";

                std::pair<scalar_iterator,bool> __it = M_scalar.insert( v );

                if ( __it.second )
                    DVLOG(2) << v.first << " loaded and inserted (value: " << v.second.first << ")\n";

                else
                    DVLOG(2) << v.first << " was loaded but not inserted (value: " << v.second.first << ")\n";
            }

            ar & boost::serialization::make_nvp( "map_complex_size", s );
            DVLOG(2) << "(loading) serialized size of  complex (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                std::pair<std::string, std::pair<complex_type,bool> > v;
                scalar_type v_real;
                scalar_type v_imag;

                ar & boost::serialization::make_nvp( "complex", v.first );
                ar & boost::serialization::make_nvp( "real", v_real );
                ar & boost::serialization::make_nvp( "imaginary", v_imag );

                v.second.first = complex_type( v_real, v_imag );

                //                     DVLOG(2) << "(loading) dserialized complex  " << v.first
                //                                   << " with value " << v.second << "\n";

                std::pair<complex_iterator,bool> __it = M_complex.insert( v );

                if ( __it.second )
                    DVLOG(2) << v.first << " loaded and inserted\n";

                else
                    DVLOG(2) << v.first << " was loaded but not inserted\n";
            }

            // nodal scalar
            s = 0;
            ar & boost::serialization::make_nvp( "map_nodal_scalar_size", s );
            DVLOG(2) << "(loading) serialized size of nodal scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                nodal_scalar_type v;
                ar & v;

                DVLOG(2) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<nodal_scalar_iterator,bool> __it = M_nodal_scalar.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    DVLOG(2) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    DVLOG(2) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }

            s = 0;
            ar & boost::serialization::make_nvp( "map_nodal_vector_size", s );
            DVLOG(2) << "(loading) serialized size of nodal scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                nodal_vector_type v;
                ar & v;

                DVLOG(2) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<nodal_vector_iterator,bool> __it =  M_nodal_vector.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    DVLOG(2) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    DVLOG(2) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }

            s = 0;
            ar & boost::serialization::make_nvp( "map_nodal_tensor2_size", s );
            DVLOG(2) << "(loading) serialized size of nodal scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                nodal_tensor2_type v;
                ar & v;

                DVLOG(2) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<nodal_tensor2_iterator,bool> __it =  M_nodal_tensor2.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    DVLOG(2) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    DVLOG(2) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }

            s = 0;
            ar & boost::serialization::make_nvp( "map_nodal_tensor2symm_size", s );
            DVLOG(2) << "(loading) serialized size of nodal scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                nodal_tensor2symm_type v;
                ar & v;

                DVLOG(2) << "(loading) dserialized scalar field " << v.name()
                         << " of size " << v.size() << "\n";

                std::pair<nodal_tensor2symm_iterator,bool> __it =  M_nodal_tensor2symm.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    DVLOG(2) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    DVLOG(2) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }

            // element scalar
            s = 0;
            ar & boost::serialization::make_nvp( "map_element_scalar_size", s );
            DVLOG(2) << "(loading) serialized size of element scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                element_scalar_type v;
                ar & v;

                DVLOG(2) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<element_scalar_iterator,bool> __it = M_element_scalar.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    DVLOG(2) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    DVLOG(2) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }

            s = 0;
            ar & boost::serialization::make_nvp( "map_element_vector_size", s );
            DVLOG(2) << "(loading) serialized size of element scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                element_vector_type v;
                ar & v;

                DVLOG(2) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<element_vector_iterator,bool> __it =  M_element_vector.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    DVLOG(2) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    DVLOG(2) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }

            s = 0;
            ar & boost::serialization::make_nvp( "map_element_tensor2_size", s );
            DVLOG(2) << "(loading) serialized size of element scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                element_tensor2_type v;
                ar & v;

                DVLOG(2) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<element_tensor2_iterator,bool> __it =  M_element_tensor2.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    DVLOG(2) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    DVLOG(2) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }
            s = 0;
            ar & boost::serialization::make_nvp( "map_element_tensor2symm_size", s );
            DVLOG(2) << "(loading) serialized size of element scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                element_tensor2symm_type v;
                ar & v;

                DVLOG(2) << "(loading) dserialized scalar field " << v.name()
                         << " of size " << v.size() << "\n";

                std::pair<element_tensor2symm_iterator,bool> __it =  M_element_tensor2symm.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    DVLOG(2) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    DVLOG(2) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }
#endif
        }

        template<class Archive>
        FEELPP_NO_EXPORT void serialize( Archive & ar, const unsigned int file_version )
        {
            boost::serialization::split_member( ar, *this, file_version );
        }

        static
        std::set<std::string>
        representationType( variant_representation_arg_type const& _rep, std::string const& valueIfEmpty = "" )
            {
                std::set<std::string> reps;
                if( auto repStringPtr = std::get_if<std::string>(&_rep))
                {
                    if ( !repStringPtr->empty() )
                        reps.insert( *repStringPtr );
                }
                else if ( auto repSetPtr = std::get_if<std::set<std::string>>(&_rep))
                {
                    reps = *repSetPtr;
                }
                if ( reps.empty() && !valueIfEmpty.empty() )
                    reps.insert( valueIfEmpty );

                for ( std::string const& rep : reps )
                    CHECK( rep == "nodal" || rep == "element" ) << "invalid represation " << rep << ": should be nodal or element";

                return reps;
            }


    private:

        TimeSet* M_ts;

        /**
           Time associated with the step
        */
        Real M_time;

        /**
           keep track of the step index
        */
        size_type M_index;
        //! active index
        size_type M_activeIndex;

        std::optional<mesh_ptrtype> M_mesh;

        map_scalar_type M_scalar;
        map_complex_type M_complex;
        map_nodal_type M_nodal;
        map_element_type M_element;

        Context M_state;

        scalar_p0_space_ptrtype M_scalar_p0;
        scalar_p1_space_ptrtype M_scalar_p1;
    };

    struct ltstep
    {
        bool operator()( std::shared_ptr<Step> const& s1,
                         std::shared_ptr<Step> const& s2 ) const
        {
            return *s1 < *s2;
        }
    };

    /** @name Typedefs
     */
    //@{

    typedef Step step_type;
    typedef std::shared_ptr<Step> step_ptrtype;
    typedef std::set<step_ptrtype,ltstep> step_set_type;
    typedef typename step_set_type::iterator step_iterator;
    typedef typename step_set_type::const_iterator step_const_iterator;
    typedef typename step_set_type::reverse_iterator step_reverse_iterator;
    typedef typename step_set_type::const_reverse_iterator step_const_reverse_iterator;
    //@}




    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * constructor for a time set
     *
     * @param filename name of the file that stores the timeset information
     * @param init if true, remove the file that stores the timset info if it exists
     */
    TimeSet( std::string filename = "undefined", bool init = false );
    TimeSet( TimeSet const & );
    ~TimeSet();

    //@}

    /** @name Operator overloads
     */
    //@{

    TimeSet& operator=( TimeSet const& );

    //@}

    /** @name Accessors
     */
    //@{

    /**
       \return the name of the time set
    */
    std::string const&  name() const
    {
        return M_name;
    }

    /**
       \return the name of the time set
    */
    void setName( std::string const& name ) const
    {
        M_name = name;
    }

    /**
       \return the index of the time set
    */
    uint32_type index() const
    {
        return M_index;
    }

    /**
       \return the number of steps already stored
    */
    size_type numberOfSteps() const
    {
        return M_step_set.size();
    }

    //! \return the number of active steps
    size_type numberOfActiveSteps() const
    {
        size_type res = 0;
        auto __it = this->beginStep();
        auto __end = this->endStep();
        for ( ; __it != __end ; ++__it )
        {
            step_ptrtype __step = *__it;
            if ( !__step->isIgnored() )
                ++res;
        }
        return res;
    }

    //! return the first active step
    step_ptrtype firstActiveStep() const
        {
            auto __it = this->beginStep();
            auto __end = this->endStep();
            for ( ; __it != __end ; ++__it )
            {
                step_ptrtype __step = *__it;
                if ( !__step->isIgnored() )
                    return __step;
            }
            step_ptrtype __step;
            return __step;
        }

    /**
       \return the time increment between two steps
    */
    Real timeIncrement() const
    {
        return M_time_increment;
    }

    step_set_type stepsToWriteOnDisk() const
        {
            auto __it = this->beginStep();
            auto __end = this->endStep();
            step_set_type stepToWriteOnDisk;
            for ( ; __it != __end ; ++__it )
            {
                step_ptrtype __step = *__it;
                if ( __step->isOnDisk() || __step->isIgnored() )
                    continue;
                if (  __step->hasData() && __step->isInMemory() )
                    stepToWriteOnDisk.insert( __step );
            }
            return stepToWriteOnDisk;
        }

    /** @name  Mutators
     */
    //@{

    /**
       set the time increment
    */
    void setTimeIncrement( Real __incr )
    {
        M_time_increment = __incr;
    }


    void setNumberOfStepsInMemory( uint16_type __i )
    {
        M_keep_steps = __i;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * \param __time time at which we want to get the step
     * \param freq active step by frequence
     * \return a step defined at time \c __time if not found then generate a new one
     */
    step_ptrtype step( Real __time, int freq = 1 );


    step_iterator beginStep()
    {
        return M_step_set.begin();
    }
    step_const_iterator beginStep() const
    {
        return M_step_set.begin();
    }

    step_reverse_iterator rbeginStep()
    {
        return M_step_set.rbegin();
    }
    step_const_reverse_iterator rbeginStep() const
    {
        return M_step_set.rbegin();
    }

    step_iterator endStep()
    {
        return M_step_set.end();
    }
    step_const_iterator endStep() const
    {
        return M_step_set.end();
    }

    step_reverse_iterator rendStep()
    {
        return M_step_set.rend();
    }
    step_const_reverse_iterator rendStep() const
    {
        return M_step_set.rend();
    }


    std::pair<step_iterator,bool> insertStep( step_ptrtype __step )
    {
        return M_step_set.insert( __step );
    }

    void clear()
    {
        M_step_set.clear();
    }

    void load( std::string const& _nameFile, Real __time )
    {
        std::ifstream ifs( _nameFile );
        // load data from archive
        boost::archive::text_iarchive ia( ifs );
        ia >> *this;

        resetPreviousTime( __time );
    }

    void save( std::string const& _nameFile, WorldComm const& worldComm )
    {
        if ( worldComm.isMasterRank() )
        {
            std::ofstream ofs( _nameFile );
            // save data from archive
            boost::archive::text_oarchive oa( ofs );
            oa << *this;
        }
    }


    //@}

protected:

    /**
       name of the time set
    */
    mutable std::string M_name;

    /**
       index of the time set
    */
    uint32_type M_index;

    /**
       steps
    */
    step_set_type M_step_set;

    /**
       time increment
    */
    Real M_time_increment;

    uint16_type M_keep_steps;
private:

    friend class boost::serialization::access;

    template<class Archive>
    FEELPP_NO_EXPORT void serialize( Archive & ar, const unsigned int /*version*/ )
    {
        ar & boost::serialization::make_nvp( "name", M_name );
        ar & boost::serialization::make_nvp( "index", M_index );
        ar & boost::serialization::make_nvp( "time_increment", M_time_increment );
        ar & boost::serialization::make_nvp( "keep_steps", M_keep_steps );

        if ( Archive::is_saving::value )
        {
            size_type s = M_step_set.size();
            ar & boost::serialization::make_nvp( "number_of_steps", s );

            step_iterator __it = M_step_set.begin();
            step_iterator __en = M_step_set.end();

            for ( ; __it != __en; ++__it )
            {
                double t = ( *__it )->time();
                ar & boost::serialization::make_nvp( "time", t );
                size_type ind = ( *__it )->index();
                ar & boost::serialization::make_nvp( "index", ind );
                size_type state = ( *__it )->state();
                ar & boost::serialization::make_nvp( "state", state );
            }

        }

        if ( Archive::is_loading::value )
        {
            size_type s( 0 );
            size_type nActiveIndex = 0;
            size_type activeIndex = invalid_v<size_type>;
            ar & boost::serialization::make_nvp( "number_of_steps", s );

            for ( size_type __i = 0; __i < s; ++__i )
            {
                double t = 0;
                ar & boost::serialization::make_nvp( "time", t );

                size_type ind = 0;
                ar & boost::serialization::make_nvp( "index", ind );

                size_type __state = 0;
                ar & boost::serialization::make_nvp( "state", __state );
                Context ctxState(__state);
                if ( !ctxState.test( STEP_IGNORED ) )
                {
                    ++nActiveIndex;
                    activeIndex = nActiveIndex;
                }
                else
                    activeIndex = invalid_v<size_type>;

                step_iterator __sit;
                bool __inserted;
                boost::tie( __sit, __inserted ) = M_step_set.insert( step_ptrtype( new Step( this, t, ind, activeIndex, __state ) ) );

                CHECK( __inserted ) <<  "insertion failed at t="<< t << " and ind="<< ind;
            }

        }
    }

    //! delete all time next this __time (after a restart by exemple)
    FEELPP_NO_EXPORT void resetPreviousTime( Real __time );

public:
    void setMesh( mesh_ptrtype m ) { M_mesh = m; }
    bool hasMesh() const { return M_mesh != std::nullopt; }
    mesh_ptrtype mesh() const
        {
            DLOG_IF( WARNING, hasMesh() ) << "Time Set has no mesh data structure associated\n";
            return M_mesh.value();
        }


public:
    std::optional<mesh_ptrtype> M_mesh;


    scalar_p0_space_ptrtype M_scalar_p0;
    scalar_p1_space_ptrtype M_scalar_p1;
private:

    /**
       keep track of the current index of the TimeSets to ensure they get
    */
    static uint32_type _S_current_index;


};
template<typename MeshType, int N>
inline FEELPP_EXPORT  bool operator<( TimeSet<MeshType, N> const& __ts1, TimeSet<MeshType, N> const& __ts2 )
{
    return __ts1.index() < __ts2.index();
}

template<typename MeshType, int N>
inline FEELPP_EXPORT bool operator<( std::shared_ptr<TimeSet<MeshType, N> > const& __ts1,
                                     std::shared_ptr<TimeSet<MeshType, N> > const& __ts2 )
{
    return __ts1->index() < __ts2->index();
}

template<typename MeshType, int N>
FEELPP_NO_EXPORT bool operator<( typename TimeSet<MeshType, N>::step_type const& __s1,
                                 typename TimeSet<MeshType, N>::step_type const& __s2 )
{
    //return __s1->index() < __s2->index();
    return __s1->time() < __s2->time();
}

template<typename MeshType, int N>
FEELPP_NO_EXPORT bool operator<( typename TimeSet<MeshType, N>::step_ptrtype const& __s1,
                                 typename TimeSet<MeshType, N>::step_ptrtype const& __s2 )
{
    return __s1->index() < __s2->index();
}

template<typename MeshType, int N>
TimeSet<MeshType, N>::TimeSet( std::string __name, bool init )
    :
    M_name( __name ),
    M_index( _S_current_index++ ),
    M_time_increment( 0.1 ),
    M_keep_steps( 1 )
{}

template<typename MeshType, int N>
TimeSet<MeshType, N>::TimeSet( TimeSet const& __ts )
    :
    M_name( __ts.M_name ),
    M_index( __ts.M_index ),
    M_time_increment( __ts.M_time_increment ),
    M_keep_steps( __ts.M_keep_steps )
{
}

template<typename MeshType, int N>
TimeSet<MeshType, N>::~TimeSet()
{
    --_S_current_index;
}

template<typename MeshType, int N>
void
TimeSet<MeshType, N>::resetPreviousTime( Real __time )
{
    step_iterator __it = M_step_set.begin();
    step_iterator __en = M_step_set.end();
    bool find=false;

    while ( !find &&  __it != __en )
    {
        double t = ( *__it )->time();
        double eps = 1e-10;
        if ( ( t-eps ) <= __time )
            ++__it;
        else
            find=true;
    }

    if ( find ) M_step_set.erase( __it,__en );

}


template<typename MeshType, int N>
typename TimeSet<MeshType, N>::step_ptrtype
TimeSet<MeshType, N>::step( Real __time,  int freq )
{
    size_type nStepBefore = 0, nActiveStepBefore = 0;
    step_iterator __sit = beginStep();
    for ( ; __sit != endStep(); ++__sit )
    {
        if ( math::abs( ( *__sit )->time() - __time ) < 1e-10 )
            break;
        ++nStepBefore;
        if ( !( *__sit )->isIgnored() )
            ++nActiveStepBefore;
    }

    bool __ignoreStep = (nStepBefore % freq) > 0;

    if ( __sit == endStep() )
    {
        bool __inserted;
        DVLOG(2) << "[TimeSet<MeshType, N>::step] Inserting new step at time " << __time << " with index " << numberOfSteps();

        size_type theState = __ignoreStep? STEP_NEW|STEP_OVERWRITE|STEP_IGNORED : STEP_NEW|STEP_OVERWRITE;
        size_type activeIndex = __ignoreStep? invalid_v<size_type> : nActiveStepBefore+1;
        step_ptrtype thestep( new Step( this, __time, nStepBefore/*numberOfSteps()*/ + 1, activeIndex, theState ) );
        if ( this->hasMesh() && !__ignoreStep )
            thestep->setMesh( this->mesh() );
        boost::tie( __sit, __inserted ) = insertStep( thestep );

        DVLOG(2) << "[TimeSet<MeshType, N>::step] step was inserted properly? " << ( __inserted?"yes":"no" ) << "\n";
        DVLOG(2) << "[TimeSet<MeshType, N>::step] step index : " << ( *__sit )->index() << " time : " << ( *__sit )->time() << "\n";
        //cleanup();
    }

    else
    {
        DVLOG(2) << "[TimeSet<MeshType, N>::step] found step at time " << __time << " with index "
                      << ( *__sit )->index() << "\n";
        ( *__sit )->showMe( "TimesSet::step(t)" );
    }

    return *__sit;
}

template<typename MeshType, int N>
uint32_type TimeSet<MeshType, N>::_S_current_index = TS_INITIAL_INDEX;

//
// TimeSet<MeshType, N>::Step
//

template<typename MeshType, int N>
TimeSet<MeshType, N>::Step::Step()
    :
    M_time( 0 ),
    M_index( 0 ),
    M_activeIndex( 0 ),
    M_state( STEP_NEW|STEP_OVERWRITE )
{
}

template<typename MeshType, int N>
TimeSet<MeshType, N>::Step::Step( TimeSet* ts, Real __t, size_type __index, size_type activeIndex, size_type __state )
    :
    M_ts( ts ),
    M_time( __t ),
    M_index( __index ),
    M_activeIndex( activeIndex ),
    M_state( __state )
{
    showMe( "Step::Step()" );
}
template<typename MeshType, int N>
TimeSet<MeshType, N>::Step::~Step()
{
    this->setState( STEP_ON_DISK );
    showMe( "Step::~Step()" );
}

}
}
#endif /* __TimeSet_H */
