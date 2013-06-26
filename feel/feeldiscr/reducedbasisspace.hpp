/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2013-04-07

   Copyright (C) 2004 EPFL
   Copyright (C) 2006-2012 Université Joseph Fourier (Grenoble I)

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
   \file ReducedBasisSpace.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-07
*/
#ifndef __ReducedBasisSpace_H
#define __ReducedBasisSpace_H 1

#include <boost/shared_ptr.hpp>

#include <feel/feel.hpp>
#include <feel/feelcrb/crb.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>


namespace Feel
{

/**
 * \class ReducedBasisSpace
 * \brief Reduced Basis Space class
 *
 * @author Christophe Prud'homme
 * @see
 */

template<typename ModelType,
         typename MeshType,
         typename A1 = parameter::void_,
         typename A2 = parameter::void_,
         typename A3 = parameter::void_,
         typename A4 = parameter::void_>
class ReducedBasisSpace : public FunctionSpace<MeshType,A1,A2,A3,A4>
{

    typedef FunctionSpace<MeshType,A1,A2,A3,A4> super;
    typedef boost::shared_ptr<super> super_ptrtype;

public :

    typedef double value_type;

    typedef ModelType model_type;
    typedef boost::shared_ptr<model_type> model_ptrtype;

    typedef MeshType mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename super::element_type space_element_type;
    typedef boost::shared_ptr<space_element_type> space_element_ptrtype;

    typedef std::vector< space_element_type > basis_type;
    typedef boost::shared_ptr<basis_type> basis_ptrtype;

    typedef ReducedBasisSpace<ModelType, MeshType, A1, A2, A3, A4> this_type;
    typedef boost::shared_ptr<this_type> this_ptrtype;

    typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> eigen_vector_type;


    ReducedBasisSpace( boost::shared_ptr<ModelType> model , boost::shared_ptr<MeshType> mesh )
        :
        super ( mesh ),
        M_mesh( mesh ),
        M_model( model )
    {
        LOG( INFO ) <<" ReducedBasisSpace constructor " ;
        this->init();
    }

    //copy constructor
    ReducedBasisSpace( ReducedBasisSpace const& rb )
        :
        super( rb.M_mesh ),
        M_basis( rb.M_basis ),
        M_mesh( rb.M_mesh )
    {}

    ReducedBasisSpace( this_ptrtype const& rb )
        :
        super ( rb->M_mesh ),
        M_basis( rb->M_basis ),
        M_mesh( rb->M_mesh )
    {}

    void init()
    {
        //M_basis = basis_ptrtype( new basis_type() );
    }

    /*
     * Get the mesh
     */
    mesh_ptrtype mesh()
    {
        return M_mesh;
    }

    /*
     * add a new basis
     */
    void addBasisElement( space_element_type const & e )
    {
        M_basis.push_back( e );
    }

    void addBasisElement( space_element_ptrtype const & e )
    {
        M_basis.push_back( *e );
    }


    /*
     * Get basis of the reduced basis space
     */
    basis_type basis()
    {
        return M_basis;
    }


    /*
     * return true if the function space is composite
     */
    bool isComposite()
    {
        return super::is_composite;
    }
    /*
     * size of the reduced basis space : number of basis functions
     */
    int size()
    {
        return M_basis.size();
    }

    /*
     * visualize basis of rb space
     */
    void visualize ()
    {
    }

    static this_ptrtype New ( boost::shared_ptr<ModelType>  const& model, mesh_ptrtype const& mesh)
    {
        return this_ptrtype ( new this_type( model ,  mesh) );
    }

    model_ptrtype model()
    {
        return M_model;
    }

    /*
     * class Context
     * stock the evaluation of basis functions in given points
     */
    class ContextRb : public FunctionSpace<MeshType,A1,A2,A3,A4>::Context
    {

        typedef typename FunctionSpace<MeshType,A1,A2,A3,A4>::Context super;

    public :

        typedef ReducedBasisSpace<ModelType,MeshType,A1,A2,A3,A4> rbspace_type;
        typedef boost::shared_ptr<rbspace_type> rbspace_ptrtype;

        typedef rbspace_type::super_ptrtype functionspace_ptrtype;

        typedef rbspace_type::eigen_vector_type eigen_vector_type;

        typedef rbspace_type::basis_type basis_type;
        typedef boost::shared_ptr<basis_type> basis_ptrtype;

        typedef Eigen::MatrixXd eigen_matrix_type;


        ~ContextRb() {}

        ContextRb( rbspace_ptrtype rbspace )
            :
            super( rbspace->model()->functionSpace() ),
            M_rbspace( rbspace )
        {}

        //only because the function functionSpace()
        //returns an object : functionspace_ptrtype
        ContextRb( functionspace_ptrtype functionspace )
            :
            super( functionspace )
        {}

        void update ( )
        {
            int nb_pts = this->nPoints();
            int N = M_rbspace->size();
            M_phi.resize( N , nb_pts );
            for( int i=0; i<N; i++)
            {
                //evaluation of the i^th basis function
                //to all nodes in the FEM context ctx
                auto evaluation = M_rbspace->evaluateBasis( i , *this );
                for(int j=0; j<evaluation.size(); j++)
                    M_phi( i , j ) = evaluation( j );
            }
            DVLOG( 2 ) << "matrix M_phi : \n"<<M_phi;
        }

        /*
         * for given element coefficients, evaluate the element at node given in context_fem
         */
        eigen_vector_type id( eigen_vector_type coeffs , bool need_to_update=true)
        {
            if( need_to_update )
                this->update();
            int npts = super::nPoints();
            auto prod = coeffs.transpose()*M_phi;
            Eigen::Matrix<value_type, Eigen::Dynamic, 1> result( npts );
            result = prod.transpose();
            DVLOG( 2 ) << "coeffs : \n"<<coeffs;
            DVLOG( 2 ) << "prod : \n"<<prod;
            DVLOG( 2 ) << "result : \n"<<result;
            return result;
        }

        virtual functionspace_ptrtype functionSpace() const
        {
            return M_rbspace;
        }

    private :
        rbspace_ptrtype M_rbspace;
        eigen_matrix_type M_phi;
    };

    /**
     * \return function space context
     */
    ContextRb context() { return ContextRb( boost::dynamic_pointer_cast< ReducedBasisSpace<ModelType,MeshType,A1,A2,A3,A4> >( this->shared_from_this() ) ); }

    /*
    * evaluate the i^th basis function at nodes
    * given by the FEM context ctx
    */
    eigen_vector_type evaluateBasis( int i , ContextRb const& ctx )
    {
        return evaluateFromContext( _context=ctx , _expr=idv(M_basis[i]) );
    }


    /*
     *  contains coefficients in the RB expression
     *  i.e. contains coefficients u_i in
     *  u_N ( \mu ) = \sum_i^N u_i ( \mu ) \xi_i( x )
     *  where \xi_i( x ) are basis function
     */
    template<typename T = double,  typename Cont = Eigen::Matrix<T,Eigen::Dynamic,1> >
    class Element
        :
        public Cont,boost::addable<Element<T,Cont> >, boost::subtractable<Element<T,Cont> >
    {
    public:
        typedef ReducedBasisSpace<ModelType,MeshType,A1,A2,A3,A4> rbspace_type;
        typedef rbspace_type functionspace_type;
        typedef boost::shared_ptr<rbspace_type> rbspace_ptrtype;

        typedef Eigen::MatrixXd eigen_matrix_type;
        typedef rbspace_type::eigen_vector_type eigen_vector_type;

        //element of the FEM function space
        typedef rbspace_type::space_element_type space_element_type;

        //RB context
        typedef rbspace_type::ContextRb ctxrb_type;

        typedef T value_type;

        typedef Cont super;
        typedef Cont container_type;


        friend class ReducedBasisSpace<ModelType,MeshType,A1,A2,A3,A4>;

        Element()
        {}

        Element( rbspace_ptrtype const& rbspace , std::string const& name="u")
            :
            M_rbspace( rbspace ),
            M_name( name )
        {
            this->resize( M_rbspace.size() );
        }

        void setReducedBasisSpace( rbspace_ptrtype rbspace )
        {
            M_rbspace = rbspace;
        }

        int size() const
        {
            return super::size();
        }

        /**
         * \return the container read-only
         */
        super const& container() const
        {
            return *this;
        }

        /**
         * \return the container read-write
         */
        super & container()
        {
            return *this;
        }


        void setCoefficient( int index , double value )
        {
            int size=super::size();
            FEELPP_ASSERT( index < size )(index)(size).error("invalid index");
            this->operator()( index ) = value;
        }

        value_type  operator()( size_t i ) const
        {
            int size=super::size();
            FEELPP_ASSERT( i < size )(i)(size).error("invalid index");
            return super::operator()( i );
        }

        value_type& operator()( size_t i )
        {
            int size = super::size();
            FEELPP_ASSERT( i < size )(i)(size).error("invalid index");
            return super::operator()( i );
        }

        Element& operator+=( Element const& _e )
        {
            int size1=super::size();
            int size2=_e.size();
            FEELPP_ASSERT( size1 == size2 )(size1)(size2).error("invalid size");
            for ( int i=0; i <size1; ++i )
                this->operator()( i ) += _e( i );
            return *this;
        }

        Element& operator-=( Element const& _e )
        {
            int size1=super::size();
            int size2=_e.size();
            FEELPP_ASSERT( size1 == size2 )(size1)(size2).error("invalid size");
            for ( int i=0; i <size2; ++i )
                this->operator()( i ) -= _e( i );
            return *this;
        }

        space_element_type expansion( int  N=-1)
        {
            return M_rbspace.expansion( *this, N );
        }

        //evaluate the element to nodes in context
        eigen_vector_type evaluate(  ctxrb_type & context_rb )
        {
           return context_rb.id( *this );
        }

        void updateGlobalValues() const
        {
            //nothing to do
        }

        private:
        rbspace_type M_rbspace;
        std::string M_name;
    };//Element of the reduced basis space

    typedef Element<value_type> element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    //Access to an element of the reduced basis space
    element_type element( std::string const& name="u" )
    {
        element_type u(boost::dynamic_pointer_cast<ReducedBasisSpace<ModelType,MeshType,A1,A2,A3,A4> >(this->shared_from_this() ), name );
        u.setZero();
        return u;
    }

    element_ptrtype elementPtr( std::string const& name="u" )
    {
        element_ptrtype u( new element_type( boost::dynamic_pointer_cast<ReducedBasisSpace<ModelType,MeshType,A1,A2,A3,A4> >(this->shared_from_this() ) , name ) );
        u->setZero();
        return u;
    }

    space_element_type expansion( element_type const& unknown, int  N=-1)
    {
        int number_of_coeff;
        int basis_size = M_basis.size();
        if ( N == -1 )
            number_of_coeff = basis_size;
        else
            number_of_coeff = N;
        FEELPP_ASSERT( number_of_coeff <= basis_size )( number_of_coeff )( basis_size ).error("invalid size");
        return Feel::expansion( M_basis, unknown , number_of_coeff );
    }

private :
    //std::vector< ctxfem_type > M_ctxfem ;
    basis_type M_basis;
    mesh_ptrtype M_mesh;
    model_ptrtype M_model;

};//ReducedBasisSpace

template<int Order, typename ModelType , typename MeshType>
inline
boost::shared_ptr<ReducedBasisSpace<ModelType,MeshType,bases<Lagrange<Order,Scalar,Continuous>>>>
RbSpacePch(  boost::shared_ptr<ModelType> const& model , boost::shared_ptr<MeshType> const& mesh  )
{
    return ReducedBasisSpace<ModelType,MeshType,bases<Lagrange<Order,Scalar,Continuous>>>::New( model , mesh );
}



}//namespace Feel

#endif
