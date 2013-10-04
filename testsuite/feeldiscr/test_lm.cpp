/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include <feel/feel.hpp>

#include <boost/noncopyable.hpp>
#include <boost/signals2.hpp>
#include <boost/format.hpp>

#include <feel/feelpoly/im.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>

using namespace Feel;
using namespace Feel::vf;




int
main( int argc, char** argv )
{

    typedef Simplex<2> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    typedef Lagrange<2, Vectorial,Continuous,PointSetFekete> basis_u_type; // velocity
    typedef Lagrange<1, Scalar,Continuous,PointSetFekete> basis_p_type; // pressure
    typedef Lagrange<2, Scalar,Continuous,PointSetFekete> basis_t_type; // temperature
    typedef Lagrange<0, Scalar, Continuous> basis_l_type; // lagrange

    typedef bases< basis_u_type , basis_p_type , basis_l_type, basis_t_type> basis_type;

    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=feel_options() );



    auto mesh = unitSquare();

    space_ptrtype Xh = space_type::New( mesh );

    vector_ptrtype R( backend()->newVector( Xh ) );
    sparse_matrix_ptrtype J( backend()->newMatrix( Xh,Xh ) );
}
