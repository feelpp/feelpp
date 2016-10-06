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

#ifndef __timeSet_H
#define __timeSet_H 1

#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>


#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
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
#include <feel/feeldiscr/interpolate.hpp>
#include <feel/feeldiscr/subelements.hpp>

#include <feel/feelmesh/filters.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/elementdiv.hpp>

#define TS_INITIAL_INDEX 1

namespace Feel
{
namespace fs = boost::filesystem;

enum
{
    STEP_NEW       = ( 1<<0 ),
    STEP_HAS_DATA  = ( 1<<1 ),
    STEP_ON_DISK   = ( 1<<2 ),
    STEP_IN_MEMORY = ( 1<<3 ),
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
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    typedef FunctionSpace<MeshType, Feel::bases<Lagrange<0,Scalar,Discontinuous> >, Discontinuous > scalar_p0_space_type;
    typedef FunctionSpace<MeshType, Feel::bases<Lagrange<0,Vectorial,Discontinuous> >, Discontinuous > vector_p0_space_type;
    typedef FunctionSpace<MeshType, Feel::bases<Lagrange<0,Tensor2,Discontinuous> >, Discontinuous > tensor2_p0_space_type;
    typedef FunctionSpace<MeshType, Feel::bases<Lagrange<0,Tensor2Symm,Discontinuous> >, Discontinuous > tensor2symm_p0_space_type;
    typedef FunctionSpace<MeshType, Feel::bases<Lagrange<N,Scalar> > > scalar_p1_space_type;
    typedef FunctionSpace<MeshType, Feel::bases<Lagrange<N,Vectorial> > > vector_p1_space_type;
    typedef FunctionSpace<MeshType, Feel::bases<Lagrange<N,Tensor2> > > tensor2_p1_space_type;
    typedef FunctionSpace<MeshType, Feel::bases<Lagrange<N,Tensor2Symm> > > tensor2symm_p1_space_type;
    typedef boost::shared_ptr<scalar_p0_space_type> scalar_p0_space_ptrtype;
    typedef boost::shared_ptr<vector_p0_space_type> vector_p0_space_ptrtype;
    typedef boost::shared_ptr<tensor2_p0_space_type> tensor2_p0_space_ptrtype;
    typedef boost::shared_ptr<tensor2symm_p0_space_type> tensor2symm_p0_space_ptrtype;
    typedef boost::shared_ptr<scalar_p1_space_type> scalar_p1_space_ptrtype;
    typedef boost::shared_ptr<vector_p1_space_type> vector_p1_space_ptrtype;
    typedef boost::shared_ptr<tensor2_p1_space_type> tensor2_p1_space_ptrtype;
    typedef boost::shared_ptr<tensor2symm_p1_space_type> tensor2symm_p1_space_ptrtype;


    typedef typename scalar_p0_space_type::element_type element_scalar_type;
    typedef typename vector_p0_space_type::element_type element_vector_type;
    typedef typename tensor2_p0_space_type::element_type element_tensor2_type;
    typedef typename tensor2symm_p0_space_type::element_type element_tensor2symm_type;
    typedef typename scalar_p1_space_type::element_type nodal_scalar_type;
    typedef typename vector_p1_space_type::element_type nodal_vector_type;
    typedef typename tensor2_p1_space_type::element_type nodal_tensor2_type;
    typedef typename tensor2symm_p1_space_type::element_type nodal_tensor2symm_type;


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
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


        typedef FunctionSpace<MeshType, Feel::bases<Lagrange<0,Scalar,Discontinuous> >, Discontinuous > scalar_p0_space_type;
        typedef FunctionSpace<MeshType, Feel::bases<Lagrange<0,Vectorial,Discontinuous> >, Discontinuous > vector_p0_space_type;
        typedef FunctionSpace<MeshType, Feel::bases<Lagrange<0,Tensor2,Discontinuous> >, Discontinuous > tensor2_p0_space_type;
        typedef FunctionSpace<MeshType, Feel::bases<Lagrange<0,Tensor2Symm,Discontinuous> >, Discontinuous > tensor2symm_p0_space_type;
        typedef FunctionSpace<MeshType, Feel::bases<Lagrange<N,Scalar> > > scalar_p1_space_type;
        typedef FunctionSpace<MeshType, Feel::bases<Lagrange<N,Vectorial> > > vector_p1_space_type;
        typedef FunctionSpace<MeshType, Feel::bases<Lagrange<N,Tensor2> > > tensor2_p1_space_type;
        typedef FunctionSpace<MeshType, Feel::bases<Lagrange<N,Tensor2Symm> > > tensor2symm_p1_space_type;
        typedef boost::shared_ptr<scalar_p0_space_type> scalar_p0_space_ptrtype;
        typedef boost::shared_ptr<vector_p0_space_type> vector_p0_space_ptrtype;
        typedef boost::shared_ptr<tensor2_p0_space_type> tensor2_p0_space_ptrtype;
        typedef boost::shared_ptr<tensor2symm_p0_space_type> tensor2symm_p0_space_ptrtype;
        typedef boost::shared_ptr<scalar_p1_space_type> scalar_p1_space_ptrtype;
        typedef boost::shared_ptr<vector_p1_space_type> vector_p1_space_ptrtype;
        typedef boost::shared_ptr<tensor2_p1_space_type> tensor2_p1_space_ptrtype;
        typedef boost::shared_ptr<tensor2symm_p1_space_type> tensor2symm_p1_space_ptrtype;


        typedef typename scalar_p0_space_type::element_type element_scalar_type;
        typedef typename vector_p0_space_type::element_type element_vector_type;
        typedef typename tensor2_p0_space_type::element_type element_tensor2_type;
        typedef typename tensor2symm_p0_space_type::element_type element_tensor2symm_type;
        typedef typename scalar_p1_space_type::element_type nodal_scalar_type;
        typedef typename vector_p1_space_type::element_type nodal_vector_type;
        typedef typename tensor2_p1_space_type::element_type nodal_tensor2_type;
        typedef typename tensor2symm_p1_space_type::element_type nodal_tensor2symm_type;

        typedef std::map<std::string, std::pair<scalar_type,bool> > map_scalar_type;
        typedef std::map<std::string, std::pair<complex_type,bool> > map_complex_type;
        typedef std::map<std::string, nodal_scalar_type> map_nodal_scalar_type;
        typedef std::map<std::string, nodal_vector_type> map_nodal_vector_type;
        typedef std::map<std::string, nodal_tensor2_type> map_nodal_tensor2_type;
        typedef std::map<std::string, nodal_tensor2symm_type> map_nodal_tensor2symm_type;
        typedef std::map<std::string, element_scalar_type> map_element_scalar_type;
        typedef std::map<std::string, element_vector_type> map_element_vector_type;
        typedef std::map<std::string, element_tensor2_type> map_element_tensor2_type;
        typedef std::map<std::string, element_tensor2symm_type> map_element_tensor2symm_type;

        typedef typename map_scalar_type::iterator scalar_iterator;
        typedef typename map_scalar_type::const_iterator scalar_const_iterator;

        typedef typename map_complex_type::iterator complex_iterator;
        typedef typename map_complex_type::const_iterator complex_const_iterator;

        typedef typename map_nodal_scalar_type::iterator nodal_scalar_iterator;
        typedef typename map_nodal_scalar_type::const_iterator nodal_scalar_const_iterator;

        typedef typename map_nodal_vector_type::iterator nodal_vector_iterator;
        typedef typename map_nodal_vector_type::const_iterator nodal_vector_const_iterator;

        typedef typename map_nodal_tensor2_type::iterator nodal_tensor2_iterator;
        typedef typename map_nodal_tensor2_type::const_iterator nodal_tensor2_const_iterator;
        typedef typename map_nodal_tensor2symm_type::iterator nodal_tensor2symm_iterator;
        typedef typename map_nodal_tensor2symm_type::const_iterator nodal_tensor2symm_const_iterator;

        typedef typename map_element_scalar_type::iterator element_scalar_iterator;
        typedef typename map_element_scalar_type::const_iterator element_scalar_const_iterator;

        typedef typename map_element_vector_type::iterator element_vector_iterator;
        typedef typename map_element_vector_type::const_iterator element_vector_const_iterator;

        typedef typename map_element_tensor2_type::iterator element_tensor2_iterator;
        typedef typename map_element_tensor2_type::const_iterator element_tensor2_const_iterator;
        typedef typename map_element_tensor2symm_type::iterator element_tensor2symm_iterator;
        typedef typename map_element_tensor2symm_type::const_iterator element_tensor2symm_const_iterator;

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

        /**
           \return true if the mesh is available
        */
        bool hasMesh() const
        {
            return M_mesh != boost::none;
        }

        /**
           \return a mesh
        */
        mesh_ptrtype mesh()
        {
            return M_mesh.get();
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
         * get the nodal scalar field with name n
         * @param __n name of the nodal scalar field
         * @return the nodal scalar field
         */
        nodal_scalar_type nodalScalar( std::string const& __n ) const
        {
            if ( M_nodal_scalar.find( sanitize(__n) ) == M_nodal_scalar.end() )
            {
                std::ostringstream __err;
                __err << "invalid nodal scalar field name " << sanitize(__n);
                throw std::logic_error( __err.str() );
            }

            return M_nodal_scalar.find( sanitize(__n) )->second.first;
        }

        /**
         * get the nodal vector field with name n
         * @param __n name of the field
         * @return the nodal vector field
         */
        nodal_vector_type nodalVector( std::string const& __n ) const
        {
            if ( M_nodal_vector.find( sanitize(__n) ) == M_nodal_vector.end() )
            {
                std::ostringstream __err;
                __err << "invalid nodal vector field name " << sanitize(__n);
                throw std::logic_error( __err.str() );
            }

            return M_nodal_vector.find( sanitize(__n) )->second;
        }

        /**
         * get the nodal tensor2 field with name n
         * @param __n name of the field
         * @return the nodal tensor2 field
         */
        nodal_tensor2_type nodalTensor2( std::string const& __n ) const
        {
            if ( M_nodal_tensor2.find( sanitize(__n) ) == M_nodal_tensor2.end() )
            {
                std::ostringstream __err;
                __err << "invalid nodal tensor2 field name " << sanitize(__n);
                throw std::logic_error( __err.str() );
            }

            return M_nodal_tensor2.find( sanitize(__n) )->second;
        }

        /**
         * get the nodal tensor2 field with name n
         * @param __n name of the field
         * @return the nodal tensor2 field
         */
        nodal_tensor2symm_type nodalTensor2Symm( std::string const& __n ) const
            {
                if ( M_nodal_tensor2symm.find( sanitize(__n) ) == M_nodal_tensor2symm.end() )
                {
                    std::ostringstream __err;
                    __err << "invalid nodal tensor2 field name " << sanitize(__n);
                    throw std::logic_error( __err.str() );
                }

                return M_nodal_tensor2symm.find( sanitize(__n) )->second;
            }

        /**
         * get the element scalar field with name n
         * @param __n name of the element scalar field
         * @return the element scalar field
         */
        element_scalar_type elementScalar( std::string const& __n ) const
        {
            if ( M_element_scalar.find( sanitize(__n) ) == M_element_scalar.end() )
            {
                std::ostringstream __err;
                __err << "invalid element scalar field name " << sanitize(__n);
                throw std::logic_error( __err.str() );
            }

            return M_element_scalar.find( sanitize(__n) )->second;
        }

        /**
         * get the element tensor2 field with name n
         * @param __n name of the field
         * @return the element tensor2 field
         */
        element_tensor2_type elementTensor2( std::string const& __n ) const
            {
                if ( M_element_tensor2.find( sanitize(__n) ) == M_element_tensor2.end() )
                {
                    std::ostringstream __err;
                    __err << "invalid element tensor2 field name " << sanitize(__n);
                    throw std::logic_error( __err.str() );
                }

                return M_element_tensor2.find( sanitize(__n) )->second;
            }

        /**
         * get the element vector field with name n
         * @param __n name of the field
         * @return the element vector field
         */
        element_vector_type elementVector( std::string const& __n ) const
        {
            if ( M_element_vector.find( sanitize(__n) ) == M_element_vector.end() )
            {
                std::ostringstream __err;
                __err << "invalid element vector field name " << sanitize(__n);
                throw std::logic_error( __err.str() );
            }

            return M_element_vector.find( sanitize(__n) )->second;
        }

        /**
         * get the element tensor2symm field with name n
         * @param __n name of the field
         * @return the element tensor2symm field
         */
        element_tensor2symm_type elementTensor2Symm( std::string const& __n ) const
            {
                if ( M_element_tensor2symm.find( sanitize(__n) ) == M_element_tensor2symm.end() )
                {
                    std::ostringstream __err;
                    __err << "invalid element tensor2symm field name " << sanitize(__n);
                    throw std::logic_error( __err.str() );
                }

                return M_element_tensor2symm.find( sanitize(__n) )->second;
            }


        //@}

        /** @name  Mutators
         */
        //@{

        void setState( size_type __st )
        {
            executeState( __st );
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
            M_scalar[sanitize(name)] =  std::make_pair( __s, cst );
            M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            M_state.clear( STEP_ON_DISK );
        }
        template<typename T>
        void add( std::string const& name, T const& __s, bool cst = false,
                  typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr )
            {
                M_scalar[sanitize(name)] =  std::make_pair( __s, cst );
                M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
                M_state.clear( STEP_ON_DISK );
            }

        void addComplex( std::string const& name, complex_type const& __s )
        {
            M_complex[sanitize(name)] =  __s;
            M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            M_state.clear( STEP_ON_DISK );
        }

        void addNodal( nodal_scalar_type const& __s )
        {
            M_nodal_scalar[sanitize(__s.name())] =  __s;
            M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            M_state.clear( STEP_ON_DISK );
        }

        void addNodal( nodal_vector_type const& __s )
        {
            M_nodal_vector[sanitize(__s.name())] =  __s;
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
            this->addRegions( prefix,prefix );
        }
        void
        addRegions( std::string const& prefix, std::string const& prefixfname )
        {
            VLOG(1) << "[timeset] Adding regions...\n";
            WorldComm const& meshComm = (M_mesh.get()->worldComm().numberOfSubWorlds() > 1)?
                M_mesh.get()->worldComm().subWorld(M_mesh.get()->worldComm().numberOfSubWorlds()) :
                M_mesh.get()->worldComm();
            this->updateScalarP0( std::vector<WorldComm>(1,meshComm) );

            VLOG(1) << "[timeset] adding pid...\n";
            this->add( prefixvm(prefix,"pid"), prefixvm(prefixfname,"pid"),  regionProcess( M_scalar_p0 ) );

            //add( "marker", regionMarker( M_scalar_p0 ) );
            //add( "marker2", regionMarker2( M_scalar_p0 ) );
            //add( "marker3", regionMarker3( M_scalar_p0 ) );
        }

        template<typename FunctionType>
        void add( std::initializer_list<std::string>  __n, FunctionType const& func,
                  typename std::enable_if<is_functionspace_element_v<FunctionType>>::type* = nullptr )
        {
            std::vector<std::string> str( sanitize( __n ) );
            add_( str, func, mpl::bool_<(FunctionType::functionspace_type::nSpaces>1)>() );
        }

        template<typename FunctionType>
        void add( std::vector<std::string> const& __n, FunctionType const& func,
                  typename std::enable_if<is_functionspace_element_v<FunctionType>>::type* = nullptr )
        {
            add_( __n, func, mpl::bool_<(FunctionType::functionspace_type::nSpaces>1)>() );
        }
        template<typename FunctionType>
        void add( std::string const& __n, FunctionType const& func,
                  typename std::enable_if<is_functionspace_element_v<FunctionType>>::type* = nullptr )
        {
            tic();
            add_( sanitize(__n), func, mpl::bool_<(FunctionType::functionspace_type::nSpaces>1)>() );
            toc((boost::format("Timeset::add %1%")%__n).str(),FLAGS_v>0);
        }

        template<typename TSet>
        struct AddFunctionProduct
        {
            AddFunctionProduct( TSet& tset ) : M_tset( tset ) {}
            TSet& M_tset;

            template<typename T>
            void
            operator()( T const& fun ) const
                {
                    LOG(INFO) << "export "  << fun.name() << " ...\n";
                    M_tset.add_( sanitize(fun.name()), fun, mpl::bool_<false>() );
                }
        };
        template<typename FunctionType>
        void add_( std::vector<std::string> const& __n, FunctionType const& func, mpl::bool_<true> )
        {
            std::vector<std::string> s = __n;
            FEELPP_ASSERT( s.size() == FunctionType::functionspace_type::nSpaces );
            // implement elements() which returns a fusion vector of the components
            fusion::for_each( subelements(func,s), AddFunctionProduct<step_type>( *this ) );
        }
        template<typename FunctionType>
        void add_( std::string const& __n, FunctionType const& func, mpl::bool_<true> )
        {
            std::vector<std::string> s;
            // s.push_back(__n);
            for(int i=0; i<FunctionType::functionspace_type::nSpaces; i++)
                {
                    std::ostringstream i_str;
                    i_str << "_" << i;
                    s.push_back(__n + i_str.str());
                }
            FEELPP_ASSERT( s.size() == FunctionType::functionspace_type::nSpaces );
            // implement elements() which returns a fusion vector of the components
            fusion::for_each( subelements(func,s), AddFunctionProduct<step_type>( *this ) );
        }
        template<typename FunctionType>
        void add_( std::string const& __n, FunctionType const& func, mpl::bool_<false> )
        {
            typedef typename mpl::or_<is_shared_ptr<FunctionType>, boost::is_pointer<FunctionType> >::type is_ptr_or_shared_ptr;
            add( __n,func,is_ptr_or_shared_ptr() );
        }
        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func )
        {
            typedef typename mpl::or_<is_shared_ptr<FunctionType>, boost::is_pointer<FunctionType> >::type is_ptr_or_shared_ptr;
            add( sanitize(__n),__fname,func,is_ptr_or_shared_ptr() );
        }

        template<typename FunctionType>
        void add( std::string const& __n, FunctionType const& func, mpl::bool_<true> )
        {
            add( __n,*func, mpl::bool_<false>() );
        }
        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<true> )
        {
            add( __n,__fname,*func, mpl::bool_<false>() );
        }

        template<typename FunctionType>
        void add( std::string const& __n, FunctionType const& func, mpl::bool_<false> )
        {
            constexpr bool nodal = (FunctionType::is_continuous || FunctionType::functionspace_type::continuity_type::is_discontinuous_locally)&& (!FunctionType::is_hcurl_conforming) && (!FunctionType::is_hdiv_conforming);
            add( __n, __n, func, mpl::bool_<false>(), mpl::bool_<nodal>() );
        }
        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false> )
        {
            constexpr bool nodal = (FunctionType::is_continuous || FunctionType::functionspace_type::continuity_type::is_discontinuous_locally)&& (!FunctionType::is_hcurl_conforming) && (!FunctionType::is_hdiv_conforming);
            add( __n, __fname, func, mpl::bool_<false>(), mpl::bool_<nodal>() );
        }
        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false>, mpl::bool_<true>,
                  typename std::enable_if<FunctionType::is_scalar>::type* = nullptr)
        {
            bool extendeddof = (soption(_name="exporter.format") == "ensightgold");
            boost::timer t;

                if ( !M_ts->M_scalar_p1 )
                {
                    if ( M_ts->M_vector_p1 && M_ts->M_vector_p1->hasCompSpace() )
                    {
                        M_ts->M_scalar_p1 = boost::make_shared<scalar_p1_space_type>();
                        M_ts->M_scalar_p1->shallowCopy( M_ts->M_vector_p1->compSpace() );
                    }
                    else
                    {
                        M_ts->M_scalar_p1 = boost::make_shared<scalar_p1_space_type> ( M_mesh.get(),
                                                                                       MESH_RENUMBER | MESH_CHECK,
                                                                                       typename scalar_p1_space_type::periodicity_type(),
                                                                                       func.worldsComm(),
                                                                                       std::vector<bool>(1,extendeddof) );
                    }
                    M_scalar_p1 = M_ts->M_scalar_p1;
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p1 created\n";
                }

                else if ( M_mesh.get() == M_ts->M_scalar_p1->mesh() )
                {
                    M_scalar_p1 = M_ts->M_scalar_p1;
                }

                if ( M_mesh.get() != M_ts->M_scalar_p1->mesh() && !M_scalar_p1 )
                {
                    M_scalar_p1 = boost::make_shared<scalar_p1_space_type> ( M_mesh.get(),
                                                                             MESH_RENUMBER | MESH_CHECK,
                                                                             typename scalar_p1_space_type::periodicity_type(),
                                                                             func.worldsComm(),
                                                                             std::vector<bool>(1,extendeddof));
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p1 created\n";
                }

                M_nodal_scalar.insert( std::make_pair( __fname, M_scalar_p1->element( __n, func.description() ) ) );

                //M_nodal_scalar[__fname].setName( __n );
                //M_nodal_scalar[__fname].setFunctionSpace( M_scalar_p1 );
                if ( func.worldComm().isActive() )
                {
                    // interpolate field on visualisation space
                    interpolate( M_scalar_p1, func, M_nodal_scalar[__fname] );
#if 0
                    // if exporter use extended dof table but the field to export is not define on this part
                    // we put to local min value for these ghosts dofs (else can disturbe visualisation range)
                    if ( M_scalar_p1->dof()->buildDofTableMPIExtended() && !func.functionSpace()->dof()->buildDofTableMPIExtended() )
                    {
                        double thelocalmin = func.min( false );
                        size_type nGhostDofAddedInExtendedDofTable = M_scalar_p1->dof()->nGhostDofAddedInExtendedDofTable();
                        size_type startDof = (M_scalar_p1->dof()->lastDof()-M_scalar_p1->dof()->firstDof()+1)-nGhostDofAddedInExtendedDofTable;
                        for ( int k=0;k<nGhostDofAddedInExtendedDofTable;++k )
                            M_nodal_scalar[__fname].set( startDof+k, thelocalmin );
                    }
#endif
                }
                M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
                M_state.clear( STEP_ON_DISK );

                showMe( "Step::add" );
                DVLOG(2) << "[timset::add] scalar time : " << t.elapsed() << "\n";
        }
        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false>, mpl::bool_<true>,
                  typename std::enable_if<FunctionType::is_vectorial>::type* = nullptr)
            {
                bool extendeddof = (soption(_name="exporter.format") == "ensightgold");
                boost::timer t;

                if ( !M_ts->M_vector_p1 )
                {
                    M_ts->M_vector_p1 = boost::make_shared<vector_p1_space_type> ( M_mesh.get(),
                                                                                    MESH_RENUMBER | MESH_CHECK,
                                                                                    typename vector_p1_space_type::periodicity_type(),
                                                                                    func.worldsComm(),
                                                                                    std::vector<bool>(1,extendeddof) );
                    M_vector_p1 = M_ts->M_vector_p1;
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p1 created\n";
                }

                else if ( M_mesh.get() == M_ts->M_vector_p1->mesh() )
                {
                    M_vector_p1 = M_ts->M_vector_p1;
                }

                if ( M_mesh.get() != M_ts->M_vector_p1->mesh() && !M_vector_p1 )
                {
                    M_vector_p1 = boost::make_shared<vector_p1_space_type> ( M_mesh.get(),
                                                                             MESH_RENUMBER | MESH_CHECK,
                                                                             typename vector_p1_space_type::periodicity_type(),
                                                                             func.worldsComm(),
                                                                             std::vector<bool>(1,extendeddof) );
                    DVLOG(2) << "[timeset::add] setmesh :  " << t.elapsed() << "\n";
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space vector p1 created\n";
                }

                t.restart();
                M_nodal_vector.insert( std::make_pair( __fname,M_vector_p1->element( __n, func.description() ) ) );
                //M_nodal_vector[__fname].setName( __n );
                //M_nodal_vector[__fname].setFunctionSpace( M_vector_p1 );
                DVLOG(2) << "[timeset::add] setmesh :  " << t.elapsed() << "\n";
                t.restart();

                if ( func.worldComm().isActive() )
                {
                    // interpolate field on visualisation space
                    interpolate( M_vector_p1, func, M_nodal_vector[__fname] );
#if 0
                    // if exporter use extended dof table but the field to export is not define on this part
                    // we put to a random vectorial value for these ghosts dofs (else can disturbe visualisation range)
                    size_type nGhostDofAddedInExtendedDofTable = M_vector_p1->dof()->nGhostDofAddedInExtendedDofTable();
                    if ( nGhostDofAddedInExtendedDofTable>0 && !func.functionSpace()->dof()->buildDofTableMPIExtended() )
                    {
                        // get a random value (in the first elt valid)
                        std::vector<typename vector_p1_space_type::value_type> randomEval( vector_p1_space_type::nComponents1 );
                        for ( auto const& elt : elements(M_vector_p1->mesh()) )
                        {
                            if ( M_vector_p1->dof()->isElementDone( elt.id() ) && !elt.isGhostCell() )
                            {
                                for ( uint16_type comp=0; comp<vector_p1_space_type::nComponents1 ; ++comp )
                                {
                                    const size_type index = M_vector_p1->dof()->localToGlobal( elt.id(), 0, comp ).index();
                                    randomEval[comp] = M_nodal_vector[__fname]( index );
                                }
                                break;
                            }
                        }
                        // update extended elt with this vectorial value
                        for ( auto const& eltWrap : elements(M_vector_p1->mesh(),EntityProcessType::GHOST_ONLY) )
                        {
                            auto const& elt = boost::unwrap_ref(eltWrap);
                            for ( uint16_type locDof=0; locDof < vector_p1_space_type::fe_type::nLocalDof; ++locDof )
                                for ( uint16_type comp=0; comp<vector_p1_space_type::nComponents1 ; ++comp )
                                {
                                    const size_type index = M_vector_p1->dof()->localToGlobal( elt.id(), locDof, comp ).index();
                                    if ( index <= ( M_vector_p1->dof()->lastDof()-nGhostDofAddedInExtendedDofTable) )
                                        continue;
                                    M_nodal_vector[__fname].set( index, randomEval[comp] );
                                }
                        }
                    } // if ( nGhostDofAddedInExtendedDofTable>0 && !func.functionSpace()->dof()->buildDofTableMPIExtended() )
#endif
                } // if ( func.worldComm().isActive() )

                M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
                M_state.clear( STEP_ON_DISK );

                showMe( "Step::add" );
                DVLOG(2) << "[timset::add] scalar time : " << t.elapsed() << "\n";
            }
        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false>, mpl::bool_<true>,
                  typename std::enable_if<FunctionType::is_tensor2>::type* = nullptr)
            {
                bool extendeddof = (soption(_name="exporter.format") == "ensightgold");
                boost::timer t;

                if ( !M_ts->M_tensor2_p1 )
                {
                    M_ts->M_tensor2_p1 = boost::make_shared<tensor2_p1_space_type> ( M_mesh.get(),
                                                                                     MESH_RENUMBER | MESH_CHECK,
                                                                                     typename tensor2_p1_space_type::periodicity_type(),
                                                                                     func.worldsComm(),
                                                                                     std::vector<bool>(1,extendeddof) );
                    M_tensor2_p1 = M_ts->M_tensor2_p1;
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p1 created\n";
                }

                else if ( M_mesh.get() == M_ts->M_tensor2_p1->mesh() )
                {
                    M_tensor2_p1 = M_ts->M_tensor2_p1;
                }

            if ( M_mesh.get() != M_ts->M_tensor2_p1->mesh() && !M_tensor2_p1 )
                {
                    M_tensor2_p1 = boost::make_shared<tensor2_p1_space_type> ( M_mesh.get(),
                                                                               MESH_RENUMBER | MESH_CHECK,
                                                                               typename tensor2_p1_space_type::periodicity_type(),
                                                                               func.worldsComm(),
                                                                               std::vector<bool>(1,extendeddof) );
                    DVLOG(2) << "[timeset::add] setmesh :  " << t.elapsed() << "\n";
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space tensor2 p1 created\n";
                }

                t.restart();
                M_nodal_tensor2.insert( std::make_pair( __fname,M_tensor2_p1->element( __n, func.description() ) ) );
                //M_nodal_tensor2[__fname].setName( __n );
                //M_nodal_tensor2[__fname].setFunctionSpace( M_tensor2_p1 );
                DVLOG(2) << "[timeset::add] setmesh :  " << t.elapsed() << "\n";
                t.restart();

                if ( func.worldComm().isActive() )
                {
                    // interpolate field on visualisation space
                    interpolate( M_tensor2_p1, func, M_nodal_tensor2[__fname] );
#if 0
                    // if exporter use extended dof table but the field to export is not define on this part
                    // we put to a random vectorial value for these ghosts dofs (else can disturbe visualisation range)
                    size_type nGhostDofAddedInExtendedDofTable = M_tensor2_p1->dof()->nGhostDofAddedInExtendedDofTable();
                    if ( nGhostDofAddedInExtendedDofTable>0 && !func.functionSpace()->dof()->buildDofTableMPIExtended() )
                    {
                        // get a random value (in the first elt valid)
                        std::vector<typename tensor2_p1_space_type::value_type> randomEval( tensor2_p1_space_type::nComponents );
                        for ( auto const& elt : elements(M_tensor2_p1->mesh()) )
                        {
                            if ( M_tensor2_p1->dof()->isElementDone( elt.id() ) && !elt.isGhostCell() )
                            {
                                for ( uint16_type comp=0; comp<tensor2_p1_space_type::nComponents ; ++comp )
                                {
                                    const size_type index = M_tensor2_p1->dof()->localToGlobal( elt.id(), 0, comp ).index();
                                    randomEval[comp] = M_nodal_tensor2[__fname]( index );
                                }
                                break;
                            }
                        }
                        // update extended elt with this tensor2ial value
                        for ( auto const& eltWrap : elements(M_tensor2_p1->mesh(),EntityProcessType::GHOST_ONLY) )
                        {
                            auto const& elt = boost::unwrap_ref(eltWrap);
                            for ( uint16_type locDof=0; locDof < tensor2_p1_space_type::fe_type::nLocalDof; ++locDof )
                                for ( uint16_type comp=0; comp<tensor2_p1_space_type::nComponents1 ; ++comp )
                                {
                                    const size_type index = M_tensor2_p1->dof()->localToGlobal( elt.id(), locDof, comp ).index();
                                    if ( index <= ( M_tensor2_p1->dof()->lastDof()-nGhostDofAddedInExtendedDofTable) )
                                        continue;
                                    M_nodal_tensor2[__fname].set( index, randomEval[comp] );
                                }
                        }
                    } // if ( nGhostDofAddedInExtendedDofTable>0 && !func.functionSpace()->dof()->buildDofTableMPIExtended() )
#endif
                } // if ( func.worldComm().isActive() )

                M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
                M_state.clear( STEP_ON_DISK );

                showMe( "Step::add" );
                DVLOG(2) << "[timset::add] scalar time : " << t.elapsed() << "\n";
            }

        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false>, mpl::bool_<true>,
                  typename std::enable_if<FunctionType::is_tensor2symm>::type* = nullptr)
            {
                bool extendeddof = (soption(_name="exporter.format") == "ensightgold");
                boost::timer t;

                if ( !M_ts->M_tensor2symm_p1 )
                {
                    M_ts->M_tensor2symm_p1 = boost::make_shared<tensor2symm_p1_space_type> ( M_mesh.get(),
                                                                                     MESH_RENUMBER | MESH_CHECK,
                                                                                     typename tensor2symm_p1_space_type::periodicity_type(),
                                                                                     func.worldsComm(),
                                                                                     std::vector<bool>(1,extendeddof) );
                    M_tensor2symm_p1 = M_ts->M_tensor2symm_p1;
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p1 created\n";
                }

                else if ( M_mesh.get() == M_ts->M_tensor2symm_p1->mesh() )
                {
                    M_tensor2symm_p1 = M_ts->M_tensor2symm_p1;
                }

            if ( M_mesh.get() != M_ts->M_tensor2symm_p1->mesh() && !M_tensor2symm_p1 )
                {
                    M_tensor2symm_p1 = boost::make_shared<tensor2symm_p1_space_type> ( M_mesh.get(),
                                                                               MESH_RENUMBER | MESH_CHECK,
                                                                               typename tensor2symm_p1_space_type::periodicity_type(),
                                                                               func.worldsComm(),
                                                                               std::vector<bool>(1,extendeddof) );
                    DVLOG(2) << "[timeset::add] setmesh :  " << t.elapsed() << "\n";
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space tensor2symm p1 created\n";
                }

                t.restart();
                M_nodal_tensor2symm.insert( std::make_pair( __fname,M_tensor2symm_p1->element( __n, func.description() ) ) );
                //M_nodal_tensor2symm[__fname].setName( __n );
                //M_nodal_tensor2symm[__fname].setFunctionSpace( M_tensor2symm_p1 );
                DVLOG(2) << "[timeset::add] setmesh :  " << t.elapsed() << "\n";
                t.restart();

                if ( func.worldComm().isActive() )
                {
                    // interpolate field on visualisation space
                    interpolate( M_tensor2symm_p1, func, M_nodal_tensor2symm[__fname] );
#if 0
                    // if exporter use extended dof table but the field to export is not define on this part
                    // we put to a random vectorial value for these ghosts dofs (else can disturbe visualisation range)
                    size_type nGhostDofAddedInExtendedDofTable = M_tensor2symm_p1->dof()->nGhostDofAddedInExtendedDofTable();
                    if ( nGhostDofAddedInExtendedDofTable>0 && !func.functionSpace()->dof()->buildDofTableMPIExtended() )
                    {
                        // get a random value (in the first elt valid)
                        std::vector<typename tensor2symm_p1_space_type::value_type> randomEval( tensor2symm_p1_space_type::nComponents );
                        for ( auto const& elt : elements(M_tensor2symm_p1->mesh()) )
                        {
                            if ( M_tensor2symm_p1->dof()->isElementDone( elt.id() ) && !elt.isGhostCell() )
                            {
                                for ( uint16_type comp=0; comp<tensor2symm_p1_space_type::nComponents ; ++comp )
                                {
                                    const size_type index = M_tensor2symm_p1->dof()->localToGlobal( elt.id(), 0, comp ).index();
                                    randomEval[comp] = M_nodal_tensor2symm[__fname]( index );
                                }
                                break;
                            }
                        }
                        // update extended elt with this tensor2symmial value
                        for ( auto const& eltWrap : elements(M_tensor2symm_p1->mesh(),EntityProcessType::GHOST_ONLY) )
                        {
                            auto const& elt = boost::unwrap_ref(eltWrap);
                            for ( uint16_type locDof=0; locDof < tensor2symm_p1_space_type::fe_type::nLocalDof; ++locDof )
                                for ( uint16_type comp=0; comp<tensor2symm_p1_space_type::nComponents1 ; ++comp )
                                {
                                    const size_type index = M_tensor2symm_p1->dof()->localToGlobal( elt.id(), locDof, comp ).index();
                                    if ( index <= ( M_tensor2symm_p1->dof()->lastDof()-nGhostDofAddedInExtendedDofTable) )
                                        continue;
                                    M_nodal_tensor2symm[__fname].set( index, randomEval[comp] );
                                }
                        }
                    } // if ( nGhostDofAddedInExtendedDofTable>0 && !func.functionSpace()->dof()->buildDofTableMPIExtended() )
#endif
                } // if ( func.worldComm().isActive() )

                M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
                M_state.clear( STEP_ON_DISK );

                showMe( "Step::add" );
                DVLOG(2) << "[timset::add] scalar time : " << t.elapsed() << "\n";
            }

        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false>, mpl::bool_<false>,
                  typename std::enable_if<FunctionType::is_scalar>::type* = nullptr )
            {
                if ( !func.worldComm().isActive() ) return;

                this->updateScalarP0( func.worldsComm() );

                // interpolate field on visualisation space
                // If user asks for P0 visu [exporter.element-spaces options contains P0]
                if( soption(_name="exporter.element-spaces").find( "P0" ) != std::string::npos )
                {
                    M_element_scalar.insert( std::make_pair( __fname,M_scalar_p0->element( __n, func.description() ) ) );
                    interpolate( M_scalar_p0, func, M_element_scalar[__fname] );
                    DVLOG(2) << "[TimeSet::add] scalar p0 function " << __n << " added to exporter with filename " << __fname <<  "\n";
                }

                // If user asks for smoothed P0 (P1) visu [exporter.element-spaces options contains P1]
                if ( soption(_name="exporter.element-spaces").find( "P1" ) != std::string::npos )
                {
                    if ( !M_ts->M_scalar_p1 )
                    {
                        bool extendeddof = (soption(_name="exporter.format") == "ensightgold");
                        M_ts->M_scalar_p1 = boost::make_shared<scalar_p1_space_type> ( M_mesh.get(),
                                                                                       MESH_RENUMBER | MESH_CHECK,
                                                                                       typename scalar_p1_space_type::periodicity_type(),
                                                                                       func.worldsComm(),
                                                                                       std::vector<bool>(1,extendeddof) );
                        M_scalar_p1 = M_ts->M_scalar_p1;
                        DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p1 created\n";
                    }

                    auto scalar_P0toP1 = M_scalar_p1->element();
                    scalar_P0toP1 = div( sum( M_scalar_p1, idv(func)*meas() ), sum( M_scalar_p1, meas() ) );
                    this->add( __n+"_P1", __fname+"_P1", scalar_P0toP1, mpl::bool_<false>(), mpl::bool_<true>());
                }

                M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
                M_state.clear( STEP_ON_DISK );
                showMe( "Step::addElement" );
                DVLOG(2) << "[TimeSet::add] p0 function done\n";
            }

        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false>, mpl::bool_<false>,
                  typename std::enable_if<FunctionType::is_vectorial>::type* = nullptr )
            {
                bool extendeddof = (soption(_name="exporter.format") == "ensightgold");
                if ( !func.worldComm().isActive() ) return;

                if ( !M_ts->M_vector_p0 )
                {
                    M_ts->M_vector_p0 = boost::make_shared<vector_p0_space_type> ( M_mesh.get(),
                                                                                   MESH_RENUMBER | MESH_CHECK,
                                                                                   typename vector_p0_space_type::periodicity_type(),
                                                                                   func.worldsComm(),
                                                                                   std::vector<bool>(1,extendeddof) );
                    M_vector_p0 = M_ts->M_vector_p0;
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space vector p0 created\n";
                }

                else if ( M_mesh.get() == M_ts->M_vector_p0->mesh() )
                {
                    M_vector_p0 = M_ts->M_vector_p0;
                }

                if (  M_mesh.get() != M_ts->M_vector_p0->mesh() && !M_vector_p0 )
                {
                    M_vector_p0 = boost::make_shared<vector_p0_space_type> ( M_mesh.get(),
                                                                             MESH_RENUMBER | MESH_CHECK,
                                                                             typename vector_p0_space_type::periodicity_type(),
                                                                             func.worldsComm(),
                                                                             std::vector<bool>(1,extendeddof) );
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space vector p0 created\n";
                }

                // interpolate field on visualisation space
                // If user asks for P0 visu [exporter.element-spaces options contains P0]
                if( soption(_name="exporter.element-spaces").find( "P0" ) != std::string::npos )
                {
                    M_element_vector.insert( std::make_pair( __fname,M_vector_p0->element( __n, func.description() ) ) );
                    interpolate( M_vector_p0, func, M_element_vector[__fname] );
                    DVLOG(2) << "[TimeSet::add] vector p0 function " << __n << " added to exporter with filename " << __fname <<  "\n";
                }

                // If user asks for smoothed P0 (P1) visu [exporter.element-spaces options contains P1]
                if ( soption(_name="exporter.element-spaces").find( "P1" ) != std::string::npos )
                {
                    if ( !M_ts->M_vector_p1 )
                    {
                        M_ts->M_vector_p1 = boost::make_shared<vector_p1_space_type> ( M_mesh.get(),
                                                                                       MESH_RENUMBER | MESH_CHECK,
                                                                                       typename vector_p1_space_type::periodicity_type(),
                                                                                       func.worldsComm(),
                                                                                       std::vector<bool>(1,extendeddof) );
                        M_vector_p1 = M_ts->M_vector_p1;
                        DVLOG(2) << "[TimeSet::setMesh] setMesh space vector p1 created\n";
                    }

                    auto vector_P0toP1 = M_vector_p1->element();
                    //vector_P0toP1 = div( sum( M_vector_p1, trans(idv(func))*meas() ), sum( M_vector_p1, meas()*one() ) );
                    vector_P0toP1 = div( sum( M_vector_p1, idv(func)*meas() ), sum( M_vector_p1, meas()*one() ) );
                    this->add( __n+"_P1", __fname+"_P1", vector_P0toP1, mpl::bool_<false>(), mpl::bool_<true>());
                }

                M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
                M_state.clear( STEP_ON_DISK );
                showMe( "Step::addElement" );
                DVLOG(2) << "[TimeSet::add] p0 function done\n";
            }

        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false>, mpl::bool_<false>,
                  typename std::enable_if<FunctionType::is_tensor2>::type* = nullptr )
            {
                bool extendeddof = (soption(_name="exporter.format") == "ensightgold");
                if ( !func.worldComm().isActive() ) return;
                if ( !M_ts->M_tensor2_p0 )
                {
                    M_ts->M_tensor2_p0 = boost::make_shared<tensor2_p0_space_type> ( M_mesh.get(),
                                                                                   MESH_RENUMBER | MESH_CHECK,
                                                                                   typename tensor2_p0_space_type::periodicity_type(),
                                                                                   func.worldsComm(),
                                                                                   std::vector<bool>(1,extendeddof) );
                    M_tensor2_p0 = M_ts->M_tensor2_p0;
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space vector p0 created\n";
                }

                else if ( M_mesh.get() == M_ts->M_tensor2_p0->mesh() )
                {
                    M_tensor2_p0 = M_ts->M_tensor2_p0;
                }

                if (  M_mesh.get() != M_ts->M_tensor2_p0->mesh() && !M_tensor2_p0 )
                {
                    M_tensor2_p0 = boost::make_shared<tensor2_p0_space_type> ( M_mesh.get(),
                                                                             MESH_RENUMBER | MESH_CHECK,
                                                                             typename tensor2_p0_space_type::periodicity_type(),
                                                                             func.worldsComm(),
                                                                             std::vector<bool>(1,extendeddof) );
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space vector p0 created\n";
                }


                M_element_tensor2[__fname].setName( __n );
                M_element_tensor2[__fname].setFunctionSpace( M_tensor2_p0 );

                // interpolate field on visualisation space
                interpolate( M_tensor2_p0, func, M_element_tensor2[__fname] );

                M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
                M_state.clear( STEP_ON_DISK );
                showMe( "Step::addElement" );
                DVLOG(2) << "[TimeSet::add] p0 function done\n";
            }

        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false>, mpl::bool_<false>,
                  typename std::enable_if<FunctionType::is_tensor2symm>::type* = nullptr )
            {
                bool extendeddof = (soption(_name="exporter.format") == "ensightgold");
                if ( !func.worldComm().isActive() ) return;

                if ( !M_ts->M_tensor2symm_p0 )
                {
                    M_ts->M_tensor2symm_p0 = boost::make_shared<tensor2symm_p0_space_type> ( M_mesh.get(),
                                                                                   MESH_RENUMBER | MESH_CHECK,
                                                                                   typename tensor2symm_p0_space_type::periodicity_type(),
                                                                                   func.worldsComm(),
                                                                                   std::vector<bool>(1,extendeddof) );
                    M_tensor2symm_p0 = M_ts->M_tensor2symm_p0;
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space vector p0 created\n";
                }

                else if ( M_mesh.get() == M_ts->M_tensor2symm_p0->mesh() )
                {
                    M_tensor2symm_p0 = M_ts->M_tensor2symm_p0;
                }

                if (  M_mesh.get() != M_ts->M_tensor2symm_p0->mesh() && !M_tensor2symm_p0 )
                {
                    M_tensor2symm_p0 = boost::make_shared<tensor2symm_p0_space_type> ( M_mesh.get(),
                                                                             MESH_RENUMBER | MESH_CHECK,
                                                                             typename tensor2symm_p0_space_type::periodicity_type(),
                                                                             func.worldsComm(),
                                                                             std::vector<bool>(1,extendeddof) );
                    DVLOG(2) << "[TimeSet::setMesh] setMesh space vector p0 created\n";
                }

                M_element_tensor2symm[__fname].setName( __n );
                M_element_tensor2symm[__fname].setFunctionSpace( M_tensor2symm_p0 );

                // interpolate field on visualisation space
                interpolate( M_tensor2symm_p0, func, M_element_tensor2symm[__fname] );

                M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
                M_state.clear( STEP_ON_DISK );
                showMe( "Step::addElement" );
                DVLOG(2) << "[TimeSet::add] p0 function done\n";
            }

        //@}

        /** @name  Methods
         */
        //@{

        nodal_scalar_const_iterator beginNodalScalar() const
        {
            return M_nodal_scalar.begin();
        }
        nodal_scalar_const_iterator endNodalScalar() const
        {
            return M_nodal_scalar.end();
        }
        nodal_vector_const_iterator beginNodalVector() const
        {
            return M_nodal_vector.begin();
        }
        nodal_vector_const_iterator endNodalVector() const
        {
            return M_nodal_vector.end();
        }
        nodal_tensor2_const_iterator beginNodalTensor2() const
        {
            return M_nodal_tensor2.begin();
        }
        nodal_tensor2_const_iterator endNodalTensor2() const
        {
            return M_nodal_tensor2.end();
        }
        nodal_tensor2symm_const_iterator beginNodalTensor2Symm() const
        {
            return M_nodal_tensor2symm.begin();
        }
        nodal_tensor2symm_const_iterator endNodalTensor2Symm() const
        {
            return M_nodal_tensor2symm.end();
        }

        element_scalar_const_iterator beginElementScalar() const
        {
            return M_element_scalar.begin();
        }
        element_scalar_const_iterator endElementScalar() const
        {
            return M_element_scalar.end();
        }
        element_vector_const_iterator beginElementVector() const
        {
            return M_element_vector.begin();
        }
        element_vector_const_iterator endElementVector() const
        {
            return M_element_vector.end();
        }
        element_tensor2_const_iterator beginElementTensor2() const
        {
            return M_element_tensor2.begin();
        }
        element_tensor2_const_iterator endElementTensor2() const
        {
            return M_element_tensor2.end();
        }
        element_tensor2symm_const_iterator beginElementTensor2Symm() const
        {
            return M_element_tensor2symm.begin();
        }
        element_tensor2symm_const_iterator endElementTensor2Symm() const
        {
            return M_element_tensor2symm.end();
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
        void showMe( std::string str ) const
        {
            DVLOG(2) << str << " index :  " << M_index << "\n";
            DVLOG(2) << str << " time :  " << M_time << "\n";
            DVLOG(2) << str << " isNew() " << isNew() << "\n";
            DVLOG(2) << str << " hasData() " << hasData() << "\n";
            DVLOG(2) << str << " isOnDisk() " << isOnDisk() << "\n";
            DVLOG(2) << str << " isInMemory() " << isInMemory() << "\n";
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
        Step();

        /**
           construct a step with index \c index at a specified time \c time

           \param time time associated with the Step
           \param index index of the step
           * \param state the state of the step
           */
        Step( TimeSet* ts, Real time, size_type index, size_type __state = STEP_NEW|STEP_OVERWRITE );


        /**
           only TimeSet can use the copy constructor
        */
        Step( Step const & );



        //@}

        /** @name Operator overloads
         */
        //@{

        Step& operator=( Step const& );

        //@}

        void updateScalarP0( std::vector<WorldComm> const& worldsComm )
        {
            if ( !M_ts->M_scalar_p0 )
            {
                bool extendeddof = (soption(_name="exporter.format") == "ensightgold");
                M_ts->M_scalar_p0 = boost::make_shared<scalar_p0_space_type> ( M_mesh.get(),
                                                                               MESH_RENUMBER | MESH_CHECK,
                                                                               typename scalar_p0_space_type::periodicity_type(),
                                                                               worldsComm,
                                                                               std::vector<bool>(1,extendeddof) );
                M_scalar_p0 = M_ts->M_scalar_p0;
                DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p0 created\n";
            }
            else if ( M_mesh.get() == M_ts->M_scalar_p0->mesh() )
            {
                M_scalar_p0 = M_ts->M_scalar_p0;
            }

            if ( M_mesh.get() != M_ts->M_scalar_p0->mesh() && !M_scalar_p0 )
            {
                bool extendeddof = (soption(_name="exporter.format") == "ensightgold");
                M_scalar_p0 = boost::make_shared<scalar_p0_space_type> ( M_mesh.get(),
                                                                         MESH_RENUMBER | MESH_CHECK,
                                                                         typename scalar_p0_space_type::periodicity_type(),
                                                                         worldsComm,
                                                                         std::vector<bool>(1,extendeddof) );
                DVLOG(2) << "[TimeSet::setMesh] setMesh space scalar p0 created\n";
            }

        }

        friend class boost::serialization::access;

        template<class Archive>
        void save( Archive & ar, const unsigned int /*version*/ ) const
        {
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
        }

        template<class Archive>
        void load( Archive & ar, const unsigned int /*version*/ )
        {

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
        }

        template<class Archive>
        void serialize( Archive & ar, const unsigned int file_version )
        {
            boost::serialization::split_member( ar, *this, file_version );
        }

        /**
         * execute the steps needed to have a coherent state
         * @param state new state to execute
         */
        void executeState( size_type state );

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

        boost::optional<mesh_ptrtype> M_mesh;

        map_scalar_type M_scalar;
        map_complex_type M_complex;
        map_nodal_scalar_type M_nodal_scalar;
        map_nodal_vector_type M_nodal_vector;
        map_nodal_tensor2_type M_nodal_tensor2;
        map_nodal_tensor2symm_type M_nodal_tensor2symm;
        map_element_scalar_type M_element_scalar;
        map_element_vector_type M_element_vector;
        map_element_tensor2_type M_element_tensor2;
        map_element_tensor2symm_type M_element_tensor2symm;

        Context M_state;


        scalar_p0_space_ptrtype M_scalar_p0;
        vector_p0_space_ptrtype M_vector_p0;
        tensor2_p0_space_ptrtype M_tensor2_p0;
        tensor2symm_p0_space_ptrtype M_tensor2symm_p0;
        scalar_p1_space_ptrtype M_scalar_p1;
        vector_p1_space_ptrtype M_vector_p1;
        tensor2_p1_space_ptrtype M_tensor2_p1;
        tensor2symm_p1_space_ptrtype M_tensor2symm_p1;


    };

    struct ltstep
    {
        bool operator()( boost::shared_ptr<Step> const& s1,
                         boost::shared_ptr<Step> const& s2 ) const
        {
            return *s1 < *s2;
        }
    };

    /** @name Typedefs
     */
    //@{

    typedef Step step_type;
    typedef boost::shared_ptr<Step> step_ptrtype;
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

    /**
       \return the number of steps ignored
    */
    size_type numberOfStepsIgnored() const
    {
        return M_stepIgnored_set.size();
    }

    /**
       \return the number of steps already stored
    */
    size_type numberOfSteps( mpl::bool_<false> /**/ ) const
    {
        return M_step_set.size();
    }

    /**
       \return the number of steps ignored
    */
    size_type numberOfSteps( mpl::bool_<true> /**/ ) const
    {
        return M_stepIgnored_set.size();
    }

    /**
       \return the number of steps already stored + steps ignored
    */
    size_type numberOfTotalSteps() const
    {
        return M_step_set.size()+M_stepIgnored_set.size();
    }

    /**
       \return the time increment between two steps
    */
    Real timeIncrement() const
    {
        return M_time_increment;
    }

    /**
       \return the step set container
    */
    //step_set_type & step_set( mpl::bool_<false> /**/ ) { return M_step_set; }

    /**
       \return the step ignored set container
    */
    //step_set_type & step_set( mpl::bool_<true> /**/ ) { return M_stepIgnored_set; }

    //@}

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
     *        __ignoreStep : exporter don't save
     * \return a step defined at time \c __time if not found then generate a new one
     */
    step_ptrtype step( Real __time, bool __ignoreStep=false )
    {
        if ( __ignoreStep )
            return step< mpl::bool_<true> >( __time );

        else
            return step< mpl::bool_<false> >( __time );
    }

    template <typename IgnoreStepType>
    step_ptrtype step( Real __time );


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


    step_iterator beginStep( mpl::bool_<false> /**/ )
    {
        return M_step_set.begin();
    }
    step_const_iterator beginStep(  mpl::bool_<false> /**/ ) const
    {
        return M_step_set.begin();
    }
    step_iterator endStep( mpl::bool_<false> /**/ )
    {
        return M_step_set.end();
    }
    step_const_iterator endStep( mpl::bool_<false> /**/ ) const
    {
        return M_step_set.end();
    }

    step_iterator beginStep( mpl::bool_<true> /**/ )
    {
        return M_stepIgnored_set.begin();
    }
    step_const_iterator beginStep(  mpl::bool_<true> /**/ ) const
    {
        return M_stepIgnored_set.begin();
    }
    step_iterator endStep( mpl::bool_<true> /**/ )
    {
        return M_stepIgnored_set.end();
    }
    step_const_iterator endStep( mpl::bool_<true> /**/ ) const
    {
        return M_stepIgnored_set.end();
    }

    std::pair<step_iterator,bool> insertStep( step_ptrtype __step, mpl::bool_<false> /**/ )
    {
        return M_step_set.insert( __step );
    }
    std::pair<step_iterator,bool> insertStep( step_ptrtype __step, mpl::bool_<true> /**/ )
    {
        return M_stepIgnored_set.insert( __step );
    }

    void clear()
    {
        M_step_set.clear();
        M_stepIgnored_set.clear();
    }

    void load( std::string _nameFile, Real __time )
    {
        fs::ifstream ifs( _nameFile );
        // load data from archive
        boost::archive::text_iarchive ia( ifs );
        ia >> *this;

        resetPreviousTime( __time );
    }

    void save( std::string _nameFile )
    {
        setNumberOfStepsInMemory( 0 );
        cleanup();
        setNumberOfStepsInMemory( 1 );

        fs::ofstream ofs( _nameFile );
        // save data from archive
        boost::archive::text_oarchive oa( ofs );
        oa << *this;
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
       steps ignored because of frequence > 1 in exporter
    */
    step_set_type M_stepIgnored_set;

    /**
       time increment
    */
    Real M_time_increment;

    uint16_type M_keep_steps;
private:

    friend class boost::serialization::access;

    template<class Archive>
    void serialize( Archive & ar, const unsigned int /*version*/ )
    {
        ar & boost::serialization::make_nvp( "name", M_name );
        ar & boost::serialization::make_nvp( "index", M_index );
        ar & boost::serialization::make_nvp( "time_increment", M_time_increment );
        ar & boost::serialization::make_nvp( "keep_steps", M_keep_steps );

        if ( Archive::is_saving::value )
        {
            size_type s = M_step_set.size();
            ar & boost::serialization::make_nvp( "number_of_steps", s );
            size_type s2 = M_stepIgnored_set.size();
            ar & boost::serialization::make_nvp( "number_of_steps_ignored", s2 );

            step_iterator __it = M_step_set.begin();
            step_iterator __en = M_step_set.end();

            for ( ; __it != __en; ++__it )
            {
                if ( !( *__it )->isOnDisk() )
                {
                    DVLOG(2) << "not including step " << ( *__it )->index()
                                  << " at time " << ( *__it )->time() << "\n";
                    ( *__it )->showMe( "TimeSet::serialize" );
                    continue;
                }

                double t = ( *__it )->time();
                ar & boost::serialization::make_nvp( "time", t );
                size_type ind = ( *__it )->index();
                ar & boost::serialization::make_nvp( "index", ind );
                size_type state = ( *__it )->state();
                ar & boost::serialization::make_nvp( "state", state );
            }

            __it = M_stepIgnored_set.begin();
            __en = M_stepIgnored_set.end();

            for ( ; __it != __en; ++__it )
            {
                if ( !( *__it )->isOnDisk() )
                {
                    DVLOG(2) << "not including step " << ( *__it )->index()
                                  << " at time " << ( *__it )->time() << "\n";
                    ( *__it )->showMe( "TimeSet::serialize" );
                    continue;
                }

                double t = ( *__it )->time();
                ar & boost::serialization::make_nvp( "time_ignored", t );
                size_type ind = ( *__it )->index();
                ar & boost::serialization::make_nvp( "index_ignored", ind );
                size_type state = ( *__it )->state();
                ar & boost::serialization::make_nvp( "state_ignored", state );
            }
        }

        if ( Archive::is_loading::value )
        {
            size_type s( 0 );
            ar & boost::serialization::make_nvp( "number_of_steps", s );
            size_type s2( 0 );
            ar & boost::serialization::make_nvp( "number_of_steps_ignored", s2 );

            for ( size_type __i = 0; __i < s; ++__i )
            {
                double t = 0;
                ar & boost::serialization::make_nvp( "time", t );

                size_type ind = 0;
                ar & boost::serialization::make_nvp( "index", ind );

                size_type __state = 0;
                ar & boost::serialization::make_nvp( "state", __state );

                step_iterator __sit;
                bool __inserted;
                boost::tie( __sit, __inserted ) = M_step_set.insert( step_ptrtype( new Step( this, t,ind, __state ) ) );

                FEELPP_ASSERT( __inserted )( t )( ind ).error ( "insertion failed" );
            }

            for ( size_type __i = 0; __i < s2; ++__i )
            {
                double t = 0;
                ar & boost::serialization::make_nvp( "time_ignored", t );

                size_type ind = 0;
                ar & boost::serialization::make_nvp( "index_ignored", ind );

                size_type __state = 0;
                ar & boost::serialization::make_nvp( "state_ignored", __state );

                step_iterator __sit;
                bool __inserted;
                boost::tie( __sit, __inserted ) = M_stepIgnored_set.insert( step_ptrtype( new Step( this, t,ind, __state ) ) );

                FEELPP_ASSERT( __inserted )( t )( ind ).error ( "insertion failed" );
            }

        }
    }

    //! cleanup steps states
    void cleanup();

    //! delete all time next this __time (after a restart by exemple)
    void resetPreviousTime( Real __time );

public:
    void setMesh( mesh_ptrtype m ) { M_mesh = m; }
    bool hasMesh() const { return M_mesh != boost::none; }
    mesh_ptrtype mesh() const
        {
            DLOG_IF( WARNING, hasMesh() ) << "Time Set has no mesh data structure associated\n";
            return M_mesh.get();
        }


public:
    boost::optional<mesh_ptrtype> M_mesh;


    scalar_p0_space_ptrtype M_scalar_p0;
    vector_p0_space_ptrtype M_vector_p0;
    tensor2_p0_space_ptrtype M_tensor2_p0;
    tensor2symm_p0_space_ptrtype M_tensor2symm_p0;
    scalar_p1_space_ptrtype M_scalar_p1;
    vector_p1_space_ptrtype M_vector_p1;
    tensor2_p1_space_ptrtype M_tensor2_p1;
    tensor2symm_p1_space_ptrtype M_tensor2symm_p1;
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
inline FEELPP_EXPORT bool operator<( boost::shared_ptr<TimeSet<MeshType, N> > const& __ts1,
                                     boost::shared_ptr<TimeSet<MeshType, N> > const& __ts2 )
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
{
    std::ostringstream __str;
    __str << M_name << ".ts";

    if ( fs::exists( __str.str() ) && !init )
    {
#if 0
        std::ifstream ifs( __str.str().c_str() );
        boost::archive::binary_iarchive ia( ifs );

        ia >> *this;

        ifs.close();
#endif
    }

    else if ( fs::exists( __str.str() ) )
        fs::remove( __str.str() );
}

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
#if 0
    step_iterator __it = M_step_set.begin();
    step_iterator __en = M_step_set.end();

    for ( ; __it != __en; ++__it )
    {
        ( *__it )->setState( STEP_ON_DISK );
    }


    // FIXME: saving the timeset is disabled for now

    std::ostringstream __str;
    __str << M_name << ".ts";

    std::ofstream ofs( __str.str().c_str() );

    if ( !ofs.fail() )
    {
        boost::archive::binary_oarchive oa( ofs );

        oa << const_cast<TimeSet<MeshType, N>const&>( *this );

        ofs.close();
    }

#endif

    --_S_current_index;


}

template<typename MeshType, int N>
TimeSet<MeshType, N>&
TimeSet<MeshType, N>::operator=( TimeSet const& __ts )
{
    if ( this != &__ts )
    {
        M_name = __ts.M_name;
        M_index = __ts.M_index;
        M_time_increment = __ts.M_time_increment;
    }

    return *this;
}

template<typename MeshType, int N>
void
TimeSet<MeshType, N>::cleanup()
{
    step_iterator __it = M_step_set.begin();
    step_iterator __en = M_step_set.end();

    for ( ; __it != __en; ++__it )
    {
        if ( ( *__it )->index() <= M_step_set.size() - M_keep_steps )
            ( *__it )->setState( STEP_ON_DISK );
    }

    __it = M_stepIgnored_set.begin();
    __en = M_stepIgnored_set.end();

    for ( ; __it != __en; ++__it )
    {
        if ( ( *__it )->index() <= M_stepIgnored_set.size() - M_keep_steps )
            ( *__it )->setState( STEP_ON_DISK );
    }


#if 0
    std::ostringstream __str;
    __str << M_name << ".ts";

    std::ofstream ofs( __str.str().c_str() );

    if ( !ofs.fail() )
    {
        boost::archive::binary_oarchive oa( ofs );

        oa << const_cast<TimeSet<MeshType, N>const&>( *this );

        ofs.close();
    }

#endif
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
        if ( !( *__it )->isOnDisk() )
        {
            DVLOG(2) << "not including step " << ( *__it )->index()
                          << " at time " << ( *__it )->time() << "\n";
            ( *__it )->showMe( "TimeSet::resetPreviousTime" );
            ++__it;
        }

        else
        {
            double t = ( *__it )->time();
            double eps = 1e-10;

            if ( ( t-eps ) <= __time ) ++__it;

            else find=true;
        }
    }

    if ( find ) M_step_set.erase( __it,__en );

    __it = M_stepIgnored_set.begin();
    __en = M_stepIgnored_set.end();
    find=false;

    while ( !find &&  __it != __en )
    {
        if ( !( *__it )->isOnDisk() )
        {
            DVLOG(2) << "not including step " << ( *__it )->index()
                          << " at time " << ( *__it )->time() << "\n";
            ( *__it )->showMe( "TimeSet::resetPreviousTime" );
            ++__it;
        }

        else
        {
            double t = ( *__it )->time();
            double eps = 1e-10;

            if ( ( t-eps ) <= __time ) ++__it;

            else find=true;
        }
    }

    if ( find ) M_stepIgnored_set.erase( __it,__en );

}


template<typename MeshType, int N>
template <typename IgnoreStepType>
typename TimeSet<MeshType, N>::step_ptrtype
TimeSet<MeshType, N>::step( Real __time )
{
    namespace lambda = boost::lambda;
    /*
    step_iterator __sit = std::find_if( M_step_set.begin(), M_step_set.end(),
                                        lambda::bind( &Step::time, lambda::_1 ) - __time < 1e-10 &&
                                        lambda::bind( &Step::time, lambda::_1 ) - __time > -1e-10 );
    */
    step_iterator __sit = beginStep( mpl::bool_<IgnoreStepType::value>() );

    for ( ; __sit != endStep( mpl::bool_<IgnoreStepType::value>() ); ++__sit )
    {
        if ( math::abs( ( *__sit )->time() - __time ) < 1e-10 )
            break;
    }

    if ( __sit == endStep( mpl::bool_<IgnoreStepType::value>() ) )
    {
        bool __inserted;

        DVLOG(2) << "[TimeSet<MeshType, N>::step] Inserting new step at time " << __time << " with index "
                      << numberOfSteps( mpl::bool_<IgnoreStepType::value>() ) << "\n";

        step_ptrtype thestep( new Step( this, __time, numberOfSteps( mpl::bool_<IgnoreStepType::value>() ) + 1 ) );
        if ( this->hasMesh() )
            thestep->setMesh( this->mesh() );
        boost::tie( __sit, __inserted ) = insertStep( thestep,mpl::bool_<IgnoreStepType::value>() );

        //boost::tie( __sit, __inserted ) = M_step_set.insert( step_ptrtype( new Step( this, __time,
        //                                                                              M_step_set.size() + 1 ) ) );

        DVLOG(2) << "[TimeSet<MeshType, N>::step] step was inserted properly? " << ( __inserted?"yes":"no" ) << "\n";
        DVLOG(2) << "[TimeSet<MeshType, N>::step] step index : " << ( *__sit )->index() << " time : " << ( *__sit )->time() << "\n";
        namespace lambda = boost::lambda;
        //std::cerr << "time values:";
        //std::for_each( M_step_set.begin(), M_step_set.end(),
        //               std::cerr << lambda::bind( &step_type::time, *lambda::_1 ) << boost::lambda::constant( ' ' ) );
        //std::cerr << "\n";

        cleanup();
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
    M_state( STEP_NEW|STEP_OVERWRITE )
{
}

template<typename MeshType, int N>
TimeSet<MeshType, N>::Step::Step( TimeSet* ts, Real __t, size_type __index, size_type __state )
    :
    M_ts( ts ),
    M_time( __t ),
    M_index( __index ),
    M_state( __state )
{
    showMe( "Step::Step()" );
}

template<typename MeshType, int N>
TimeSet<MeshType, N>::Step::Step( Step const& __step )
    :
    M_time( __step.M_time ),
    M_index( __step.M_index ),
    M_state( __step.M_state )
{
}

template<typename MeshType, int N>
TimeSet<MeshType, N>::Step::~Step()
{
    this->setState( STEP_ON_DISK );
    showMe( "Step::~Step()" );
}

template<typename MeshType, int N>
void
TimeSet<MeshType, N>::Step::executeState( size_type __st )
{
    Context __state( __st );
    DVLOG(2) << "executeState: isOnDisk() :  " << __state.test( STEP_ON_DISK ) << "\n";
    DVLOG(2) << "executeState: isInMemory() :  " << __state.test( STEP_IN_MEMORY ) << "\n";

    if ( hasData() && isInMemory()  && __state.test( STEP_ON_DISK ) )
    {
#if 0
        DVLOG(2) << "saving step " << this->index() << " at time " << this->time() << " on disk\n";
        std::ostringstream __str;
        __str << "step-" << this->index();

        if ( !fs::exists( __str.str() ) || M_state.test( STEP_OVERWRITE ) )
        {
            std::ofstream ofs( __str.str().c_str() );
            boost::archive::binary_oarchive oa( ofs );

            oa << const_cast<typename TimeSet<MeshType, N>::Step const&>( *this );

            ofs.close();
        }

#endif

        if ( hasData() )
        {
            DVLOG(2) << "releasing step " << M_index << " at time " << M_time << " allocated memory for this step\n";
            auto clear = []( auto& c ) { c.second.clear(); };
            {
                std::for_each( M_nodal_scalar.begin(),
                              M_nodal_scalar.end(),
                              clear );
                std::for_each( M_nodal_vector.begin(),
                               M_nodal_vector.end(),
                              clear );
                std::for_each( M_nodal_tensor2.begin(),
                              M_nodal_tensor2.end(),
                              clear );
                std::for_each( M_nodal_tensor2symm.begin(),
                              M_nodal_tensor2symm.end(),
                              clear );
            }
            {
                std::for_each( M_element_scalar.begin(),
                              M_element_scalar.end(),
                              clear );
                std::for_each( M_element_vector.begin(),
                              M_element_vector.end(),
                              clear );
                std::for_each( M_element_tensor2.begin(),
                              M_element_tensor2.end(),
                              clear );
                std::for_each( M_element_tensor2symm.begin(),
                              M_element_tensor2symm.end(),
                              clear );
            }
            M_state.clear( STEP_IN_MEMORY );
        }

        M_state.set( STEP_ON_DISK );
    }

    if ( hasData() && isOnDisk()  && __state.test( STEP_IN_MEMORY ) )
    {
        DVLOG(2) << "loading step " << this->index() << " at time " << this->time() << " in memory\n";
#if 0
        std::ostringstream __str;
        __str << "step-" << this->index();

        if ( fs::exists( __str.str() ) )
        {
            std::ifstream ifs( __str.str().c_str() );
            boost::archive::binary_iarchive ia( ifs );

            ia >> *this;

            ifs.close();
        }

#endif
        M_state.set( STEP_IN_MEMORY|STEP_HAS_DATA );
        M_state.clear( STEP_ON_DISK );
    }
}
}
}
#endif /* __TimeSet_H */
