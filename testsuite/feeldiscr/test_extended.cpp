// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>

#if !defined( MESH_DIM )
#define MESH_DIM 2
#endif

#define ORDER_U 2
#define ORDER_P 1
#define MY_DIM 2

int main(int argc, char**argv )
{
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv );

    // Define mesh
    typedef Simplex<MY_DIM> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    // Define spaces (Velocity,Pressure,Lagrange)
    typedef Lagrange<ORDER_U, Vectorial,Continuous,PointSetFekete> basis_u_type;
    typedef Lagrange<ORDER_P, Scalar,Continuous,PointSetFekete> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;
    typedef bases< basis_u_type , basis_p_type, basis_l_type  > basis_type_F;
    typedef FunctionSpace<mesh_type, basis_type_F> space_type_F;
    typedef boost::shared_ptr<space_type_F> space_ptrtype_F;
    // Define backend
    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;


    auto mesh = loadMesh(_mesh=new mesh_type );
    // Function spaces creation with extended doftable
    std::vector<bool> ext_doft( {true, false, false} );
    auto Xh_F = space_type_F::New( mesh,_extended_doftable=ext_doft );

    // Matrix creation with extended pattern
    backend_ptrtype backend( backend_type::build( soption( _name="backend" ) ) );
    auto blockStruct = BlocksStencilPattern(3,3)
        << size_type(Pattern::EXTENDED) << size_type(Pattern::COUPLED) << size_type(Pattern::ZERO) // velocity
        << size_type(Pattern::COUPLED)  << size_type(Pattern::COUPLED) << size_type(Pattern::COUPLED) // pressure
        << size_type(Pattern::ZERO)     << size_type(Pattern::COUPLED) << size_type(Pattern::COUPLED); // lag multiplier
    auto test_matrix = backend->newMatrix( _test=Xh_F, _trial=Xh_F, _pattern_block=blockStruct );
}
