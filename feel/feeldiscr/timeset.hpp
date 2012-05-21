/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-01-12

  Copyright (C) 2005,2006 EPFL, INRIA, Politecnico di Milano

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
#include <feel/feelcore/context.hpp>

#include <feel/feelalg/glas.hpp>
#include <feel/feelpoly/lagrange.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/interpolate.hpp>


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
/**
 * \class TimeSet
 * \ingroup SpaceTime
 * \brief data TimeSet
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


    typedef FunctionSpace<MeshType, bases<Lagrange<0,Scalar,Discontinuous> >, Discontinuous > scalar_p0_space_type;
    typedef FunctionSpace<MeshType, bases<Lagrange<0,Vectorial,Discontinuous> >, Discontinuous > vector_p0_space_type;
    typedef FunctionSpace<MeshType, bases<Lagrange<0,Tensor2,Discontinuous> >, Discontinuous > tensor2_p0_space_type;
    typedef FunctionSpace<MeshType, bases<Lagrange<N,Scalar> > > scalar_p1_space_type;
    typedef FunctionSpace<MeshType, bases<Lagrange<N,Vectorial> > > vector_p1_space_type;
    typedef FunctionSpace<MeshType, bases<Lagrange<N,Tensor2> > > tensor2_p1_space_type;
    typedef boost::shared_ptr<scalar_p0_space_type> scalar_p0_space_ptrtype;
    typedef boost::shared_ptr<vector_p0_space_type> vector_p0_space_ptrtype;
    typedef boost::shared_ptr<tensor2_p0_space_type> tensor2_p0_space_ptrtype;
    typedef boost::shared_ptr<scalar_p1_space_type> scalar_p1_space_ptrtype;
    typedef boost::shared_ptr<vector_p1_space_type> vector_p1_space_ptrtype;
    typedef boost::shared_ptr<tensor2_p1_space_type> tensor2_p1_space_ptrtype;


    typedef typename scalar_p0_space_type::element_type element_scalar_type;
    typedef typename vector_p0_space_type::element_type element_vector_type;
    typedef typename tensor2_p0_space_type::element_type element_tensor2_type;
    typedef typename scalar_p1_space_type::element_type nodal_scalar_type;
    typedef typename vector_p1_space_type::element_type nodal_vector_type;
    typedef typename tensor2_p1_space_type::element_type nodal_tensor2_type;


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

        typedef MeshType mesh_type;
        typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


        typedef FunctionSpace<MeshType, bases<Lagrange<0,Scalar,Discontinuous> >, Discontinuous > scalar_p0_space_type;
        typedef FunctionSpace<MeshType, bases<Lagrange<0,Vectorial,Discontinuous> >, Discontinuous > vector_p0_space_type;
        typedef FunctionSpace<MeshType, bases<Lagrange<0,Tensor2,Discontinuous> >, Discontinuous > tensor2_p0_space_type;
        typedef FunctionSpace<MeshType, bases<Lagrange<N,Scalar> > > scalar_p1_space_type;
        typedef FunctionSpace<MeshType, bases<Lagrange<N,Vectorial> > > vector_p1_space_type;
        typedef FunctionSpace<MeshType, bases<Lagrange<N,Tensor2> > > tensor2_p1_space_type;
        typedef boost::shared_ptr<scalar_p0_space_type> scalar_p0_space_ptrtype;
        typedef boost::shared_ptr<vector_p0_space_type> vector_p0_space_ptrtype;
        typedef boost::shared_ptr<tensor2_p0_space_type> tensor2_p0_space_ptrtype;
        typedef boost::shared_ptr<scalar_p1_space_type> scalar_p1_space_ptrtype;
        typedef boost::shared_ptr<vector_p1_space_type> vector_p1_space_ptrtype;
        typedef boost::shared_ptr<tensor2_p1_space_type> tensor2_p1_space_ptrtype;


        typedef typename scalar_p0_space_type::element_type element_scalar_type;
        typedef typename vector_p0_space_type::element_type element_vector_type;
        typedef typename tensor2_p0_space_type::element_type element_tensor2_type;
        typedef typename scalar_p1_space_type::element_type nodal_scalar_type;
        typedef typename vector_p1_space_type::element_type nodal_vector_type;
        typedef typename tensor2_p1_space_type::element_type nodal_tensor2_type;

        typedef std::map<std::string, scalar_type> map_scalar_type;
        typedef std::map<std::string, complex_type> map_complex_type;
        typedef std::map<std::string, nodal_scalar_type> map_nodal_scalar_type;
        typedef std::map<std::string, nodal_vector_type> map_nodal_vector_type;
        typedef std::map<std::string, nodal_tensor2_type> map_nodal_tensor2_type;
        typedef std::map<std::string, element_scalar_type> map_element_scalar_type;
        typedef std::map<std::string, element_vector_type> map_element_vector_type;
        typedef std::map<std::string, element_tensor2_type> map_element_tensor2_type;

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

        typedef typename map_element_scalar_type::iterator element_scalar_iterator;
        typedef typename map_element_scalar_type::const_iterator element_scalar_const_iterator;

        typedef typename map_element_vector_type::iterator element_vector_iterator;
        typedef typename map_element_vector_type::const_iterator element_vector_const_iterator;

        typedef typename map_element_tensor2_type::iterator element_tensor2_iterator;
        typedef typename map_element_tensor2_type::const_iterator element_tensor2_const_iterator;

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
            return _M_state.test( STEP_NEW );
        }

        /**
         *
         * @return
         */
        bool hasData() const
        {
            return _M_state.test( STEP_HAS_DATA );
        }

        /**
         *
         * @return
         */
        bool isOnDisk() const
        {
            return _M_state.test( STEP_ON_DISK );
        }

        /**
         *
         * @return
         */
        bool isInMemory() const
        {
            return _M_state.test( STEP_IN_MEMORY );
        }

        /**
         * get the current state of the step
         * @return the size encoded into a size_type
         */
        size_type state() const
        {
            return _M_state.context();
        }

        /**
           \return the time associated with the Step
        */
        Real time() const
        {
            return _M_time;
        }

        /**
           \return the index of the time set
        */
        size_type index() const
        {
            return _M_index;
        }

        /**
           \return true if the mesh is available
        */
        bool hasMesh() const
        {
            return _M_mesh;
        }

        /**
           \return a mesh
        */
        mesh_ptrtype mesh()
        {
            return _M_mesh.get();
        }

        /**
         * get the scalar with name n
         * @param __n name of the nodal scalar field
         * @return the scalar value
         */
        scalar_type scalar( std::string const& __n ) const
        {
            if ( _M_scalar.find( __n ) == _M_scalar.end() )
            {
                std::ostringstream __err;
                __err << "invalid scalar value name " << __n;
                throw std::logic_error( __err.str() );
            }

            return _M_scalar.find( __n )->second;
        }


        /**
         * get the nodal scalar field with name n
         * @param __n name of the nodal scalar field
         * @return the nodal scalar field
         */
        nodal_scalar_type nodalScalar( std::string const& __n ) const
        {
            if ( _M_nodal_scalar.find( __n ) == _M_nodal_scalar.end() )
            {
                std::ostringstream __err;
                __err << "invalid nodal scalar field name " << __n;
                throw std::logic_error( __err.str() );
            }

            return _M_nodal_scalar.find( __n )->second;
        }

        /**
         * get the nodal vector field with name n
         * @param __n name of the field
         * @return the nodal vector field
         */
        nodal_vector_type nodalVector( std::string const& __n ) const
        {
            if ( _M_nodal_vector.find( __n ) == _M_nodal_vector.end() )
            {
                std::ostringstream __err;
                __err << "invalid nodal vector field name " << __n;
                throw std::logic_error( __err.str() );
            }

            return _M_nodal_vector.find( __n )->second;
        }

        /**
         * get the nodal tensor2 field with name n
         * @param __n name of the field
         * @return the nodal tensor2 field
         */
        nodal_tensor2_type nodalTensor2( std::string const& __n ) const
        {
            if ( _M_nodal_tensor2.find( __n ) == _M_nodal_tensor2.end() )
            {
                std::ostringstream __err;
                __err << "invalid nodal tensor2 field name " << __n;
                throw std::logic_error( __err.str() );
            }

            return _M_nodal_tensor2.find( __n )->second;
        }

        /**
         * get the element scalar field with name n
         * @param __n name of the element scalar field
         * @return the element scalar field
         */
        element_scalar_type elementScalar( std::string const& __n ) const
        {
            if ( _M_element_scalar.find( __n ) == _M_element_scalar.end() )
            {
                std::ostringstream __err;
                __err << "invalid element scalar field name " << __n;
                throw std::logic_error( __err.str() );
            }

            return _M_element_scalar.find( __n )->second;
        }

        /**
         * get the element vector field with name n
         * @param __n name of the field
         * @return the element vector field
         */
        element_vector_type elementVector( std::string const& __n ) const
        {
            if ( _M_element_vector.find( __n ) == _M_element_vector.end() )
            {
                std::ostringstream __err;
                __err << "invalid element vector field name " << __n;
                throw std::logic_error( __err.str() );
            }

            return _M_element_vector.find( __n )->second;
        }

        /**
         * get the element tensor2 field with name n
         * @param __n name of the field
         * @return the element tensor2 field
         */
        element_tensor2_type elementTensor2( std::string const& __n ) const
        {
            if ( _M_element_tensor2.find( __n ) == _M_element_tensor2.end() )
            {
                std::ostringstream __err;
                __err << "invalid element tensor2 field name " << __n;
                throw std::logic_error( __err.str() );
            }

            return _M_element_tensor2.find( __n )->second;
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
            Debug( 8000 ) << "[TimeSet::setMesh] setMesh start\n";
            _M_mesh = __m;

            _M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            _M_state.clear( STEP_ON_DISK );
            Debug( 8000 ) << "[TimeSet::setMesh] setMesh start\n";
        }

        void addScalar( std::string const& name, scalar_type const& __s )
        {
            _M_scalar[name] =  __s;
            _M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            _M_state.clear( STEP_ON_DISK );
        }

        void addComplex( std::string const& name, complex_type const& __s )
        {
            _M_complex[name] =  __s;
            _M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            _M_state.clear( STEP_ON_DISK );
        }

        void addNodal( nodal_scalar_type const& __s )
        {
            _M_nodal_scalar[__s.name()] =  __s;
            _M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            _M_state.clear( STEP_ON_DISK );
        }

        void addNodal( nodal_vector_type const& __s )
        {
            _M_nodal_vector[__s.name()] =  __s;
            _M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            _M_state.clear( STEP_ON_DISK );
        }

        void
        addRegions()
        {
            Debug() << "[timeset] Adding regions...\n";
            if ( !M_ts->_M_scalar_p0 )
            {
                Debug() << "[timeset] creating space...\n";
                M_ts->_M_scalar_p0 = scalar_p0_space_ptrtype( new scalar_p0_space_type ( _M_mesh.get() ) );
                _M_scalar_p0 = M_ts->_M_scalar_p0;
            }
            Debug() << "[timeset] adding pid...\n";
            add( "pid", regionProcess( _M_scalar_p0 ) );
            //add( "marker", regionMarker( _M_scalar_p0 ) );
            //add( "marker2", regionMarker2( _M_scalar_p0 ) );
            //add( "marker3", regionMarker3( _M_scalar_p0 ) );
        }
        template<typename FunctionType>
        void add( std::string const& __n, FunctionType const& func )
        {
            typedef typename mpl::or_<is_shared_ptr<FunctionType>, boost::is_pointer<FunctionType> >::type is_ptr_or_shared_ptr;
            add( __n,func,is_ptr_or_shared_ptr() );
        }
        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func )
        {
            typedef typename mpl::or_<is_shared_ptr<FunctionType>, boost::is_pointer<FunctionType> >::type is_ptr_or_shared_ptr;
            add( __n,__fname,func,is_ptr_or_shared_ptr() );
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
            add( __n, __n, func, mpl::bool_<false>(), mpl::bool_<FunctionType::is_continuous || FunctionType::functionspace_type::continuity_type::is_discontinuous_locally>() );
        }
        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false> )
        {
            add( __n, __fname, func, mpl::bool_<false>(), mpl::bool_<FunctionType::is_continuous || FunctionType::functionspace_type::continuity_type::is_discontinuous_locally>() );
        }
        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false>, mpl::bool_<true> )
        {

            if ( FunctionType::is_scalar )
            {
                boost::timer t;

                if ( !M_ts->_M_scalar_p1 )
                {
                    M_ts->_M_scalar_p1 = scalar_p1_space_ptrtype( new scalar_p1_space_type ( _M_mesh.get(),
                                         MESH_RENUMBER | MESH_CHECK,
                                         typename scalar_p1_space_type::periodicity_type(),
                                         func.worldsComm() ) );
                    _M_scalar_p1 = M_ts->_M_scalar_p1;
                    Debug( 8000 ) << "[TimeSet::setMesh] setMesh space scalar p1 created\n";
                }

                else if ( _M_mesh.get() == M_ts->_M_scalar_p1->mesh() )
                {
                    _M_scalar_p1 = M_ts->_M_scalar_p1;
                }

                if ( _M_mesh.get() != M_ts->_M_scalar_p1->mesh() && !_M_scalar_p1 )
                {
                    _M_scalar_p1 = scalar_p1_space_ptrtype( new scalar_p1_space_type ( _M_mesh.get(),
                                                            MESH_RENUMBER | MESH_CHECK,
                                                            typename scalar_p1_space_type::periodicity_type(),
                                                            func.worldsComm() ) );
                    Debug( 8000 ) << "[TimeSet::setMesh] setMesh space scalar p1 created\n";
                }

                _M_nodal_scalar.insert( std::make_pair( __fname, _M_scalar_p1->element( __n ) ) );

                //_M_nodal_scalar[__fname].setName( __n );
                //_M_nodal_scalar[__fname].setFunctionSpace( _M_scalar_p1 );
                if ( func.worldComm().isActive() )
                    interpolate( _M_scalar_p1, func, _M_nodal_scalar[__fname] );

                Debug( 8000 ) << "[timset::add] scalar time : " << t.elapsed() << "\n";
            }

            else if ( FunctionType::is_vectorial )
            {
                boost::timer t;

                if ( !M_ts->_M_vector_p1 )
                {
                    M_ts->_M_vector_p1 = vector_p1_space_ptrtype( new vector_p1_space_type ( _M_mesh.get(),
                                         MESH_RENUMBER | MESH_CHECK,
                                         typename vector_p1_space_type::periodicity_type(),
                                         func.worldsComm() ) );
                    _M_vector_p1 = M_ts->_M_vector_p1;
                    Debug( 8000 ) << "[TimeSet::setMesh] setMesh space scalar p1 created\n";
                }

                else if ( _M_mesh.get() == M_ts->_M_vector_p1->mesh() )
                {
                    _M_vector_p1 = M_ts->_M_vector_p1;
                }

                if ( _M_mesh.get() != M_ts->_M_vector_p1->mesh() && !_M_vector_p1 )
                {
                    _M_vector_p1 = vector_p1_space_ptrtype( new vector_p1_space_type ( _M_mesh.get(),
                                                            MESH_RENUMBER | MESH_CHECK,
                                                            typename vector_p1_space_type::periodicity_type(),
                                                            func.worldsComm() ) );
                    Debug( 8000 ) << "[timeset::add] setmesh :  " << t.elapsed() << "\n";
                    Debug( 8000 ) << "[TimeSet::setMesh] setMesh space vector p1 created\n";
                }

                t.restart();
                _M_nodal_vector.insert( std::make_pair( __fname,_M_vector_p1->element( __n ) ) );
                //_M_nodal_vector[__fname].setName( __n );
                //_M_nodal_vector[__fname].setFunctionSpace( _M_vector_p1 );
                Debug( 8000 ) << "[timeset::add] setmesh :  " << t.elapsed() << "\n";
                t.restart();

                if ( func.worldComm().isActive() )
                    interpolate( _M_vector_p1, func, _M_nodal_vector[__fname] );

                Debug( 8000 ) << "[timset::add] scalar time : " << t.elapsed() << "\n";
            }

            _M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            _M_state.clear( STEP_ON_DISK );

            showMe( "Step::add" );
        }

        template<typename FunctionType>
        void add( std::string const& __n, std::string const& __fname, FunctionType const& func, mpl::bool_<false>, mpl::bool_<false> )
        {
            if ( !func.worldComm().isActive() ) return;

            if ( FunctionType::is_scalar )
            {
                if ( !M_ts->_M_scalar_p0 )
                {
                    M_ts->_M_scalar_p0 = scalar_p0_space_ptrtype( new scalar_p0_space_type ( _M_mesh.get(),
                                         MESH_RENUMBER | MESH_CHECK,
                                         typename scalar_p0_space_type::periodicity_type(),
                                         func.worldsComm() ) );
                    _M_scalar_p0 = M_ts->_M_scalar_p0;
                    Debug( 8000 ) << "[TimeSet::setMesh] setMesh space scalar p0 created\n";
                }

                else if ( _M_mesh.get() == M_ts->_M_scalar_p0->mesh() )
                {
                    _M_scalar_p0 = M_ts->_M_scalar_p0;
                }

                if ( _M_mesh.get() != M_ts->_M_scalar_p0->mesh() && !_M_scalar_p0 )
                {
                    _M_scalar_p0 = scalar_p0_space_ptrtype( new scalar_p0_space_type ( _M_mesh.get(),
                                                            MESH_RENUMBER | MESH_CHECK,
                                                            typename scalar_p0_space_type::periodicity_type(),
                                                            func.worldsComm() ) );
                    Debug( 8000 ) << "[TimeSet::setMesh] setMesh space scalar p0 created\n";
                }

                _M_element_scalar.insert( std::make_pair( __fname,_M_scalar_p0->element( __n ) ) );
                //_M_element_scalar[__fname].setName( __n );
                //_M_element_scalar[__fname].setFunctionSpace( _M_scalar_p0 );

                if ( func.worldComm().isActive() )
                    interpolate( _M_scalar_p0, func, _M_element_scalar[__fname] );

                //std::copy( func.begin(), func.end(), _M_element_scalar[__fname].begin() );
                Debug( 8000 ) << "[TimeSet::add] scalar p0 function " << __n << " added to exporter with filename " << __fname <<  "\n";
            }

            else if ( FunctionType::is_vectorial )
            {
                if ( !M_ts->_M_vector_p0 )
                {
                    M_ts->_M_vector_p0 = vector_p0_space_ptrtype( new vector_p0_space_type ( _M_mesh.get(),
                                         MESH_RENUMBER | MESH_CHECK,
                                         typename vector_p0_space_type::periodicity_type(),
                                         func.worldsComm() ) );
                    _M_vector_p0 = M_ts->_M_vector_p0;
                    Debug( 8000 ) << "[TimeSet::setMesh] setMesh space vector p0 created\n";
                }

                else if ( _M_mesh.get() == M_ts->_M_vector_p0->mesh() )
                {
                    _M_vector_p0 = M_ts->_M_vector_p0;
                }

                if (  _M_mesh.get() != M_ts->_M_vector_p0->mesh() && !_M_vector_p0 )
                {
                    _M_vector_p0 = vector_p0_space_ptrtype( new vector_p0_space_type ( _M_mesh.get(),
                                                            MESH_RENUMBER | MESH_CHECK,
                                                            typename vector_p0_space_type::periodicity_type(),
                                                            func.worldsComm() ) );
                    Debug( 8000 ) << "[TimeSet::setMesh] setMesh space vector p0 created\n";
                    //_M_tensor2_p0 = tensor2_p0_space_ptrtype( new tensor2_p0_space_type ( _M_mesh.get() ) );
                }

                _M_element_vector.insert( std::make_pair( __fname,_M_vector_p0->element( __n ) ) );
                //_M_element_vector[__fname].setName( __n );
                //_M_element_vector[__fname].setFunctionSpace( _M_vector_p0 );

                if ( func.worldComm().isActive() )
                    interpolate( _M_vector_p0, func, _M_element_vector[__fname] );

                //std::copy( func.begin(), func.end(), _M_element_vector[__fname].begin() );
            }

            else if ( FunctionType::is_tensor2 )
            {
                _M_element_tensor2[__fname].setName( __n );
                _M_element_tensor2[__fname].setFunctionSpace( _M_tensor2_p0 );

                if ( func.worldComm().isActive() )
                    interpolate( _M_tensor2_p0, func, _M_element_tensor2[__fname] );

                //std::copy( func.begin(), func.end(), _M_element_tensor2[__fname].begin() );
            }

            _M_state.set( STEP_HAS_DATA|STEP_IN_MEMORY );
            _M_state.clear( STEP_ON_DISK );
            showMe( "Step::addElement" );
            Debug( 8000 ) << "[TimeSet::add] p0 function done\n";
        }

        //@}

        /** @name  Methods
         */
        //@{

        nodal_scalar_const_iterator beginNodalScalar() const
        {
            return _M_nodal_scalar.begin();
        }
        nodal_scalar_const_iterator endNodalScalar() const
        {
            return _M_nodal_scalar.end();
        }
        nodal_vector_const_iterator beginNodalVector() const
        {
            return _M_nodal_vector.begin();
        }
        nodal_vector_const_iterator endNodalVector() const
        {
            return _M_nodal_vector.end();
        }
        nodal_tensor2_const_iterator beginNodalTensor2() const
        {
            return _M_nodal_tensor2.begin();
        }
        nodal_tensor2_const_iterator endNodalTensor2() const
        {
            return _M_nodal_tensor2.end();
        }

        element_scalar_const_iterator beginElementScalar() const
        {
            return _M_element_scalar.begin();
        }
        element_scalar_const_iterator endElementScalar() const
        {
            return _M_element_scalar.end();
        }
        element_vector_const_iterator beginElementVector() const
        {
            return _M_element_vector.begin();
        }
        element_vector_const_iterator endElementVector() const
        {
            return _M_element_vector.end();
        }
        element_tensor2_const_iterator beginElementTensor2() const
        {
            return _M_element_tensor2.begin();
        }
        element_tensor2_const_iterator endElementTensor2() const
        {
            return _M_element_tensor2.end();
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
            Debug( 8005 ) << str << " index :  " << _M_index << "\n";
            Debug( 8005 ) << str << " time :  " << _M_time << "\n";
            Debug( 8005 ) << str << " isNew() " << isNew() << "\n";
            Debug( 8005 ) << str << " hasData() " << hasData() << "\n";
            Debug( 8005 ) << str << " isOnDisk() " << isOnDisk() << "\n";
            Debug( 8005 ) << str << " isInMemory() " << isInMemory() << "\n";
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

        friend class boost::serialization::access;

        template<class Archive>
        void save( Archive & ar, const unsigned int /*version*/ ) const
        {
            size_type s;

            // values
            s = _M_scalar.size();
            ar & boost::serialization::make_nvp( "map_scalar_size", s );

            scalar_const_iterator __its= _M_scalar.begin();
            scalar_const_iterator __ens= _M_scalar.end();

            for ( ; __its!= __ens; ++__its )
            {
                ar & boost::serialization::make_nvp( "scalar", *__its );
            }

            s = _M_complex.size();
            ar & boost::serialization::make_nvp( "map_complex_size", s );

            complex_const_iterator __itc = _M_complex.begin();
            complex_const_iterator __enc = _M_complex.end();

            for ( ; __itc!= __enc; ++__itc )
            {
                ar & boost::serialization::make_nvp( "complex", ( std::string const& )__itc->first );
                ar & boost::serialization::make_nvp( "real", ( scalar_type const& )__itc->second.real() );
                ar & boost::serialization::make_nvp( "imaginary", ( scalar_type const& )__itc->second.imag() );
            }

            s = _M_nodal_scalar.size();
            ar & boost::serialization::make_nvp( "map_nodal_scalar_size", s );
            Debug( 8005 ) << "(saving) serialized size of nodal scalar (" << s << ")\n";
            nodal_scalar_const_iterator __itscalar = _M_nodal_scalar.begin();
            nodal_scalar_const_iterator __enscalar = _M_nodal_scalar.end();

            for ( ; __itscalar != __enscalar; ++__itscalar )
            {
                ar & ( nodal_scalar_type const& ) __itscalar->second;
            }

            s = _M_nodal_vector.size();
            ar & boost::serialization::make_nvp( "map_nodal_vector_size", s );

            nodal_vector_const_iterator __itvector= _M_nodal_vector.begin();
            nodal_vector_const_iterator __envector= _M_nodal_vector.end();

            for ( ; __itvector!= __envector; ++__itvector )
            {
                ar & ( nodal_vector_type const& ) __itvector->second;
            }

            s = _M_nodal_tensor2.size();
            ar & boost::serialization::make_nvp( "map_nodal_tensor2_size", s );

            nodal_tensor2_const_iterator __ittensor2= _M_nodal_tensor2.begin();
            nodal_tensor2_const_iterator __entensor2= _M_nodal_tensor2.end();

            for ( ; __ittensor2!= __entensor2; ++__ittensor2 )
            {
                ar & ( nodal_tensor2_type const& ) __ittensor2->second;
            }

            s = _M_element_scalar.size();
            ar & boost::serialization::make_nvp( "map_element_scalar_size", s );
            Debug( 8005 ) << "(saving) serialized size of element scalar (" << s << ")\n";
            element_scalar_const_iterator __eitscalar = _M_element_scalar.begin();
            element_scalar_const_iterator __eenscalar = _M_element_scalar.end();

            for ( ; __eitscalar != __eenscalar; ++__eitscalar )
            {
                ar & ( element_scalar_type const& ) __eitscalar->second;
            }

            s = _M_element_vector.size();
            ar & boost::serialization::make_nvp( "map_element_vector_size", s );

            element_vector_const_iterator __eitvector= _M_element_vector.begin();
            element_vector_const_iterator __eenvector= _M_element_vector.end();

            for ( ; __eitvector!= __eenvector; ++__eitvector )
            {
                ar & ( element_vector_type const& ) __eitvector->second;
            }

            s = _M_element_tensor2.size();
            ar & boost::serialization::make_nvp( "map_element_tensor2_size", s );

            element_tensor2_const_iterator __eittensor2= _M_element_tensor2.begin();
            element_tensor2_const_iterator __eentensor2= _M_element_tensor2.end();

            for ( ; __eittensor2!= __eentensor2; ++__eittensor2 )
            {
                ar & ( element_tensor2_type const& ) __eittensor2->second;
            }
        }

        template<class Archive>
        void load( Archive & ar, const unsigned int /*version*/ )
        {

            // loading
            size_type s( 0 );

            ar & boost::serialization::make_nvp( "map_scalar_size", s );
            Debug( 8005 ) << "(loading) serialized size of  scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                std::pair<std::string, scalar_type> v;
                ar & v;

                Debug( 8005 ) << "(loading) dserialized scalar  " << v.first
                              << " with value " << v.second << "\n";

                std::pair<scalar_iterator,bool> __it = _M_scalar.insert( v );

                if ( __it.second )
                    Debug( 8005 ) << v.first << " loaded and inserted (value: " << v.second << ")\n";

                else
                    Debug( 8005 ) << v.first << " was loaded but not inserted (value: " << v.second << ")\n";
            }

            ar & boost::serialization::make_nvp( "map_complex_size", s );
            Debug( 8005 ) << "(loading) serialized size of  complex (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                std::pair<std::string, complex_type> v;
                scalar_type v_real;
                scalar_type v_imag;

                ar & boost::serialization::make_nvp( "complex", v.first );
                ar & boost::serialization::make_nvp( "real", v_real );
                ar & boost::serialization::make_nvp( "imaginary", v_imag );

                v.second = complex_type( v_real, v_imag );

                //                     Debug( 8005 ) << "(loading) dserialized complex  " << v.first
                //                                   << " with value " << v.second << "\n";

                std::pair<complex_iterator,bool> __it = _M_complex.insert( v );

                if ( __it.second )
                    Debug( 8005 ) << v.first << " loaded and inserted\n";

                else
                    Debug( 8005 ) << v.first << " was loaded but not inserted\n";
            }

            // nodal scalar
            s = 0;
            ar & boost::serialization::make_nvp( "map_nodal_scalar_size", s );
            Debug( 8005 ) << "(loading) serialized size of nodal scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                nodal_scalar_type v;
                ar & v;

                Debug( 8005 ) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<nodal_scalar_iterator,bool> __it = _M_nodal_scalar.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    Debug( 8005 ) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    Debug( 8005 ) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }

            s = 0;
            ar & boost::serialization::make_nvp( "map_nodal_vector_size", s );
            Debug( 8005 ) << "(loading) serialized size of nodal scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                nodal_vector_type v;
                ar & v;

                Debug( 8005 ) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<nodal_vector_iterator,bool> __it =  _M_nodal_vector.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    Debug( 8005 ) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    Debug( 8005 ) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }

            s = 0;
            ar & boost::serialization::make_nvp( "map_nodal_tensor2_size", s );
            Debug( 8005 ) << "(loading) serialized size of nodal scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                nodal_tensor2_type v;
                ar & v;

                Debug( 8005 ) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<nodal_tensor2_iterator,bool> __it =  _M_nodal_tensor2.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    Debug( 8005 ) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    Debug( 8005 ) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }

            // element scalar
            s = 0;
            ar & boost::serialization::make_nvp( "map_element_scalar_size", s );
            Debug( 8005 ) << "(loading) serialized size of element scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                element_scalar_type v;
                ar & v;

                Debug( 8005 ) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<element_scalar_iterator,bool> __it = _M_element_scalar.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    Debug( 8005 ) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    Debug( 8005 ) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }

            s = 0;
            ar & boost::serialization::make_nvp( "map_element_vector_size", s );
            Debug( 8005 ) << "(loading) serialized size of element scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                element_vector_type v;
                ar & v;

                Debug( 8005 ) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<element_vector_iterator,bool> __it =  _M_element_vector.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    Debug( 8005 ) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    Debug( 8005 ) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
            }

            s = 0;
            ar & boost::serialization::make_nvp( "map_element_tensor2_size", s );
            Debug( 8005 ) << "(loading) serialized size of element scalar (" << s << ")\n";

            for ( size_type __i = 0; __i < s; ++__i )
            {
                element_tensor2_type v;
                ar & v;

                Debug( 8005 ) << "(loading) dserialized scalar field " << v.name()
                              << " of size " << v.size() << "\n";

                std::pair<element_tensor2_iterator,bool> __it =  _M_element_tensor2.insert( std::make_pair( v.name(), v ) );

                if ( __it.second )
                    Debug( 8005 ) << v.name() << " loaded and inserted (size: " << v.size() << ")\n";

                else
                    Debug( 8005 ) << v.name() << " was loaded but not inserted (size: " << v.size() << ")\n";
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
        Real _M_time;

        /**
           keep track of the step index
        */
        size_type _M_index;

        boost::optional<mesh_ptrtype> _M_mesh;

        map_scalar_type _M_scalar;
        map_complex_type _M_complex;
        map_nodal_scalar_type _M_nodal_scalar;
        map_nodal_vector_type _M_nodal_vector;
        map_nodal_tensor2_type _M_nodal_tensor2;
        map_element_scalar_type _M_element_scalar;
        map_element_vector_type _M_element_vector;
        map_element_tensor2_type _M_element_tensor2;

        Context _M_state;


        scalar_p0_space_ptrtype _M_scalar_p0;
        vector_p0_space_ptrtype _M_vector_p0;
        tensor2_p0_space_ptrtype _M_tensor2_p0;
        scalar_p1_space_ptrtype _M_scalar_p1;
        vector_p1_space_ptrtype _M_vector_p1;
        tensor2_p1_space_ptrtype _M_tensor2_p1;


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
        return _M_name;
    }

    /**
       \return the index of the time set
    */
    uint32_type index() const
    {
        return _M_index;
    }

    /**
       \return the number of steps already stored
    */
    size_type numberOfSteps() const
    {
        return _M_step_set.size();
    }

    /**
       \return the number of steps ignored
    */
    size_type numberOfStepsIgnored() const
    {
        return _M_stepIgnored_set.size();
    }

    /**
       \return the number of steps already stored
    */
    size_type numberOfSteps( mpl::bool_<false> /**/ ) const
    {
        return _M_step_set.size();
    }

    /**
       \return the number of steps ignored
    */
    size_type numberOfSteps( mpl::bool_<true> /**/ ) const
    {
        return _M_stepIgnored_set.size();
    }

    /**
       \return the number of steps already stored + steps ignored
    */
    size_type numberOfTotalSteps() const
    {
        return _M_step_set.size()+_M_stepIgnored_set.size();
    }

    /**
       \return the time increment between two steps
    */
    Real timeIncrement() const
    {
        return _M_time_increment;
    }

    /**
       \return the step set container
    */
    //step_set_type & step_set( mpl::bool_<false> /**/ ) { return _M_step_set; }

    /**
       \return the step ignored set container
    */
    //step_set_type & step_set( mpl::bool_<true> /**/ ) { return _M_stepIgnored_set; }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
       set the time increment
    */
    void setTimeIncrement( Real __incr )
    {
        _M_time_increment = __incr;
    }


    void setNumberOfStepsInMemory( uint16_type __i )
    {
        _M_keep_steps = __i;
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
        return _M_step_set.begin();
    }
    step_const_iterator beginStep() const
    {
        return _M_step_set.begin();
    }

    step_reverse_iterator rbeginStep()
    {
        return _M_step_set.rbegin();
    }
    step_const_reverse_iterator rbeginStep() const
    {
        return _M_step_set.rbegin();
    }

    step_iterator endStep()
    {
        return _M_step_set.end();
    }
    step_const_iterator endStep() const
    {
        return _M_step_set.end();
    }

    step_reverse_iterator rendStep()
    {
        return _M_step_set.rend();
    }
    step_const_reverse_iterator rendStep() const
    {
        return _M_step_set.rend();
    }


    step_iterator beginStep( mpl::bool_<false> /**/ )
    {
        return _M_step_set.begin();
    }
    step_const_iterator beginStep(  mpl::bool_<false> /**/ ) const
    {
        return _M_step_set.begin();
    }
    step_iterator endStep( mpl::bool_<false> /**/ )
    {
        return _M_step_set.end();
    }
    step_const_iterator endStep( mpl::bool_<false> /**/ ) const
    {
        return _M_step_set.end();
    }

    step_iterator beginStep( mpl::bool_<true> /**/ )
    {
        return _M_stepIgnored_set.begin();
    }
    step_const_iterator beginStep(  mpl::bool_<true> /**/ ) const
    {
        return _M_stepIgnored_set.begin();
    }
    step_iterator endStep( mpl::bool_<true> /**/ )
    {
        return _M_stepIgnored_set.end();
    }
    step_const_iterator endStep( mpl::bool_<true> /**/ ) const
    {
        return _M_stepIgnored_set.end();
    }

    std::pair<step_iterator,bool> insertStep( step_ptrtype __step, mpl::bool_<false> /**/ )
    {
        return _M_step_set.insert( __step );
    }
    std::pair<step_iterator,bool> insertStep( step_ptrtype __step, mpl::bool_<true> /**/ )
    {
        return _M_stepIgnored_set.insert( __step );
    }

    void clear()
    {
        _M_step_set.clear();
        _M_stepIgnored_set.clear();
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
    std::string _M_name;

    /**
       index of the time set
    */
    uint32_type _M_index;

    /**
       steps
    */
    step_set_type _M_step_set;

    /**
       steps ignored because of frequence > 1 in exporter
    */
    step_set_type _M_stepIgnored_set;

    /**
       time increment
    */
    Real _M_time_increment;

    uint16_type _M_keep_steps;
private:

    friend class boost::serialization::access;

    template<class Archive>
    void serialize( Archive & ar, const unsigned int /*version*/ )
    {
        ar & boost::serialization::make_nvp( "name", _M_name );
        ar & boost::serialization::make_nvp( "index", _M_index );
        ar & boost::serialization::make_nvp( "time_increment", _M_time_increment );
        ar & boost::serialization::make_nvp( "keep_steps", _M_keep_steps );

        if ( Archive::is_saving::value )
        {
            size_type s = _M_step_set.size();
            ar & boost::serialization::make_nvp( "number_of_steps", s );
            size_type s2 = _M_stepIgnored_set.size();
            ar & boost::serialization::make_nvp( "number_of_steps_ignored", s2 );

            step_iterator __it = _M_step_set.begin();
            step_iterator __en = _M_step_set.end();

            for ( ; __it != __en; ++__it )
            {
                if ( !( *__it )->isOnDisk() )
                {
                    Debug( 8000 ) << "not including step " << ( *__it )->index()
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

            __it = _M_stepIgnored_set.begin();
            __en = _M_stepIgnored_set.end();

            for ( ; __it != __en; ++__it )
            {
                if ( !( *__it )->isOnDisk() )
                {
                    Debug( 8000 ) << "not including step " << ( *__it )->index()
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
                boost::tie( __sit, __inserted ) = _M_step_set.insert( step_ptrtype( new Step( this, t,ind, __state ) ) );

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
                boost::tie( __sit, __inserted ) = _M_stepIgnored_set.insert( step_ptrtype( new Step( this, t,ind, __state ) ) );

                FEELPP_ASSERT( __inserted )( t )( ind ).error ( "insertion failed" );
            }

        }
    }

    //! cleanup steps states
    void cleanup();

    //! delete all time next this __time (after a restart by exemple)
    void resetPreviousTime( Real __time );

public:
    boost::optional<mesh_ptrtype> _M_mesh;


    scalar_p0_space_ptrtype _M_scalar_p0;
    vector_p0_space_ptrtype _M_vector_p0;
    tensor2_p0_space_ptrtype _M_tensor2_p0;
    scalar_p1_space_ptrtype _M_scalar_p1;
    vector_p1_space_ptrtype _M_vector_p1;
    tensor2_p1_space_ptrtype _M_tensor2_p1;
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
    _M_name( __name ),
    _M_index( _S_current_index++ ),
    _M_time_increment( 0.1 ),
    _M_keep_steps( 1 )
{
    std::ostringstream __str;
    __str << _M_name << ".ts";

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
    _M_name( __ts._M_name ),
    _M_index( __ts._M_index ),
    _M_time_increment( __ts._M_time_increment ),
    _M_keep_steps( __ts._M_keep_steps )
{
}

template<typename MeshType, int N>
TimeSet<MeshType, N>::~TimeSet()
{
#if 0
    step_iterator __it = _M_step_set.begin();
    step_iterator __en = _M_step_set.end();

    for ( ; __it != __en; ++__it )
    {
        ( *__it )->setState( STEP_ON_DISK );
    }


    // FIXME: saving the timeset is disabled for now

    std::ostringstream __str;
    __str << _M_name << ".ts";

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
        _M_name = __ts._M_name;
        _M_index = __ts._M_index;
        _M_time_increment = __ts._M_time_increment;
    }

    return *this;
}

template<typename MeshType, int N>
void
TimeSet<MeshType, N>::cleanup()
{
    step_iterator __it = _M_step_set.begin();
    step_iterator __en = _M_step_set.end();

    for ( ; __it != __en; ++__it )
    {
        if ( ( *__it )->index() <= _M_step_set.size() - _M_keep_steps )
            ( *__it )->setState( STEP_ON_DISK );
    }

    __it = _M_stepIgnored_set.begin();
    __en = _M_stepIgnored_set.end();

    for ( ; __it != __en; ++__it )
    {
        if ( ( *__it )->index() <= _M_stepIgnored_set.size() - _M_keep_steps )
            ( *__it )->setState( STEP_ON_DISK );
    }


#if 0
    std::ostringstream __str;
    __str << _M_name << ".ts";

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
    step_iterator __it = _M_step_set.begin();
    step_iterator __en = _M_step_set.end();
    bool find=false;

    while ( !find &&  __it != __en )
    {
        if ( !( *__it )->isOnDisk() )
        {
            Debug( 8000 ) << "not including step " << ( *__it )->index()
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

    if ( find ) _M_step_set.erase( __it,__en );

    __it = _M_stepIgnored_set.begin();
    __en = _M_stepIgnored_set.end();
    find=false;

    while ( !find &&  __it != __en )
    {
        if ( !( *__it )->isOnDisk() )
        {
            Debug( 8000 ) << "not including step " << ( *__it )->index()
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

    if ( find ) _M_stepIgnored_set.erase( __it,__en );

}


template<typename MeshType, int N>
template <typename IgnoreStepType>
typename TimeSet<MeshType, N>::step_ptrtype
TimeSet<MeshType, N>::step( Real __time )
{
    namespace lambda = boost::lambda;
    /*
    step_iterator __sit = std::find_if( _M_step_set.begin(), _M_step_set.end(),
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

        Debug( 8000 ) << "[TimeSet<MeshType, N>::step] Inserting new step at time " << __time << " with index "
                      << numberOfSteps( mpl::bool_<IgnoreStepType::value>() ) << "\n";

        step_ptrtype thestep( new Step( this, __time, numberOfSteps( mpl::bool_<IgnoreStepType::value>() ) + 1 ) );

        boost::tie( __sit, __inserted ) = insertStep( thestep,mpl::bool_<IgnoreStepType::value>() );

        //boost::tie( __sit, __inserted ) = _M_step_set.insert( step_ptrtype( new Step( this, __time,
        //                                                                              _M_step_set.size() + 1 ) ) );

        Debug( 8000 ) << "[TimeSet<MeshType, N>::step] step was inserted properly? " << ( __inserted?"yes":"no" ) << "\n";
        Debug( 8000 ) << "[TimeSet<MeshType, N>::step] step index : " << ( *__sit )->index() << " time : " << ( *__sit )->time() << "\n";
        namespace lambda = boost::lambda;
        //std::cerr << "time values:";
        //std::for_each( _M_step_set.begin(), _M_step_set.end(),
        //               std::cerr << lambda::bind( &step_type::time, *lambda::_1 ) << boost::lambda::constant( ' ' ) );
        //std::cerr << "\n";

        cleanup();
    }

    else
    {
        Debug( 8000 ) << "[TimeSet<MeshType, N>::step] found step at time " << __time << " with index "
                      << ( *__sit )->index() << "\n";
        ( *__sit )->showMe( "TimesSet::step(t)" );
    }

    return *__sit;
}

template<typename MeshType, int N>
uint32_type TimeSet<MeshType, N>::_S_current_index = 1;

//
// TimeSet<MeshType, N>::Step
//

template<typename MeshType, int N>
TimeSet<MeshType, N>::Step::Step()
    :
    _M_time( 0 ),
    _M_index( 0 ),
    _M_state( STEP_NEW|STEP_OVERWRITE )
{
}

template<typename MeshType, int N>
TimeSet<MeshType, N>::Step::Step( TimeSet* ts, Real __t, size_type __index, size_type __state )
    :
    M_ts( ts ),
    _M_time( __t ),
    _M_index( __index ),
    _M_state( __state )
{
    showMe( "Step::Step()" );
}

template<typename MeshType, int N>
TimeSet<MeshType, N>::Step::Step( Step const& __step )
    :
    _M_time( __step._M_time ),
    _M_index( __step._M_index ),
    _M_state( __step._M_state )
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
    Debug( 8005 ) << "executeState: isOnDisk() :  " << __state.test( STEP_ON_DISK ) << "\n";
    Debug( 8005 ) << "executeState: isInMemory() :  " << __state.test( STEP_IN_MEMORY ) << "\n";

    if ( hasData() && isInMemory()  && __state.test( STEP_ON_DISK ) )
    {
#if 0
        Debug( 8005 ) << "saving step " << this->index() << " at time " << this->time() << " on disk\n";
        std::ostringstream __str;
        __str << "step-" << this->index();

        if ( !fs::exists( __str.str() ) || _M_state.test( STEP_OVERWRITE ) )
        {
            std::ofstream ofs( __str.str().c_str() );
            boost::archive::binary_oarchive oa( ofs );

            oa << const_cast<typename TimeSet<MeshType, N>::Step const&>( *this );

            ofs.close();
        }

#endif

        if ( hasData() )
        {
            Debug( 8005 ) << "releasing step " << _M_index << " at time " << _M_time << " allocated memory for this step\n";

            {
                nodal_scalar_iterator __itscalar = _M_nodal_scalar.begin();
                nodal_scalar_iterator __enscalar = _M_nodal_scalar.end();

                for ( ; __itscalar != __enscalar; ++__itscalar )
                {
                    __itscalar->second.clear();
                }

                nodal_vector_iterator __itvector= _M_nodal_vector.begin();
                nodal_vector_iterator __envector= _M_nodal_vector.end();

                for ( ; __itvector!= __envector; ++__itvector )
                {
                    __itvector->second.clear();
                }

                nodal_tensor2_iterator __ittensor2= _M_nodal_tensor2.begin();
                nodal_tensor2_iterator __entensor2= _M_nodal_tensor2.end();

                for ( ; __ittensor2!= __entensor2; ++__ittensor2 )
                {
                    __ittensor2->second.clear();
                }
            }
            {
                element_scalar_iterator __itscalar = _M_element_scalar.begin();
                element_scalar_iterator __enscalar = _M_element_scalar.end();

                for ( ; __itscalar != __enscalar; ++__itscalar )
                {
                    __itscalar->second.clear();
                }

                element_vector_iterator __itvector= _M_element_vector.begin();
                element_vector_iterator __envector= _M_element_vector.end();

                for ( ; __itvector!= __envector; ++__itvector )
                {
                    __itvector->second.clear();
                }

                element_tensor2_iterator __ittensor2= _M_element_tensor2.begin();
                element_tensor2_iterator __entensor2= _M_element_tensor2.end();

                for ( ; __ittensor2!= __entensor2; ++__ittensor2 )
                {
                    __ittensor2->second.clear();
                }
            }
            _M_state.clear( STEP_IN_MEMORY );
        }

        _M_state.set( STEP_ON_DISK );
    }

    if ( hasData() && isOnDisk()  && __state.test( STEP_IN_MEMORY ) )
    {
        Debug( 8005 ) << "loading step " << this->index() << " at time " << this->time() << " in memory\n";
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
        _M_state.set( STEP_IN_MEMORY|STEP_HAS_DATA );
        _M_state.clear( STEP_ON_DISK );
    }
}

}
#endif /* __TimeSet_H */
