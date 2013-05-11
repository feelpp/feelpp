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

    typedef typename super::Context  ctx_type;
    typedef boost::shared_ptr<ctx_type> ctx_ptrtype;

    typedef ReducedBasisSpace<ModelType, MeshType, A1, A2, A3, A4> this_type;
    typedef boost::shared_ptr<this_type> this_ptrtype;


    ReducedBasisSpace( boost::shared_ptr<ModelType> model , boost::shared_ptr<MeshType> mesh )
        :
        super ( mesh ),
        M_mesh( mesh )
    {
        LOG( INFO ) <<" ReducedBasisSpace constructor " ;
        this->init();
    }

    //copy constructor
    ReducedBasisSpace( ReducedBasisSpace const& rb )
        :
        super( rb.M_mesh ),
        M_ctx( rb.M_ctx ),
        M_basis( rb.M_basis ),
        M_mesh( rb.M_mesh )
    {}

    ReducedBasisSpace( this_ptrtype const& rb )
        :
        super ( rb->M_mesh ),
        M_ctx( rb->M_ctx ),
        M_basis( rb->M_basis ),
        M_mesh( rb->M_mesh )
    {}

    void init()
    {
        M_basis = basis_ptrtype( new basis_type() );
        M_ctx.clear();
    }

    /*
     * Get the mesh
     */
    mesh_ptrtype mesh()
    {
        return M_mesh;
    }
    /*
     * add a context from FunctionSpace
     */
    void addContext( ctx_type ctx)
    {
        M_ctx.push_back( ctx );
        LOG( INFO ) <<"add an existing context from EIM ( for example )";
    }

    /*
     * add a new basis
     */
    void addBasisElement( space_element_type const & e )
    {
        M_basis->push_back( e );
    }

    void addBasisElement( space_element_ptrtype const & e )
    {
        M_basis->push_back( *e );
    }

    /*
     * Get basis of the reduced basis space
     */
    basis_ptrtype basis()
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
        return M_basis->size();
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
        typedef boost::shared_ptr<rbspace_type> rbspace_ptrtype;

        typedef T value_type;

        friend class ReducedBasisSpace<ModelType,MeshType,A1,A2,A3,A4>;

        Element()
        {}

        Element( rbspace_ptrtype const& rbspace , std::string const& name="u")
            :
            M_rbspace( rbspace ),
            M_name( name )
        {}

        void setReducedBasisSpace( rbspace_ptrtype rbspace )
        {
            M_rbspace = rbspace;
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

private :
    std::vector< ctx_type > M_ctx ;
    basis_ptrtype M_basis;
    mesh_ptrtype M_mesh;

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
