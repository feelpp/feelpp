/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
 */
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>
/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;


template<int Dim, int Order>
class EigenProblem
:
public Simget
{
    typedef Simget super;
public:
    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::vector_type vector_type;
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef bases<Lagrange<Order,Scalar>> basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef Exporter<mesh_type> export_type;

    void run();
private:

}; // EigenProblem

template<int Dim, int Order>
void
EigenProblem<Dim, Order>::run()
{
    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "------------------------------------------------------------\n";
        std::cout << "Execute EigenProblem<" << Dim << ">\n";
    }

    Environment::changeRepository( boost::format( "eigen/%1%/%2%D-P%3%/h_%4%/" )
                                   % this->about().appName()
                                   % Dim
                                   % Order
                                   % doption(_name="gmsh.hsize") );

    auto mesh = loadMesh(_mesh = new mesh_type );

    auto Xh = space_type::New( mesh );
    auto u = Xh->element();
    auto v = Xh->element();

    auto l = form1( _test=Xh );
    auto a = form2( _test=Xh, _trial=Xh);
    a = integrate( elements( mesh ), dxt(u)*dx(v)+dyt(u)*dy(v)+dzt(u)*dz(v));
    a+= on( boundaryfaces(mesh), _element=u, _rhs=l, _expr=cst(0.));

    auto b = form2( _test=Xh, _trial=Xh);
    b = integrate( elements(mesh), trans(idt( u ))*id( v ) );

    auto modes= veigs( _formA=a, _formB=b );

    if ( !modes.empty() )
    {
        int i = 0;
        for( auto const& mode : modes )
        {
            if ( Environment::isMasterRank() )
                std::cout << " -- eigenvalue " << i << " = " << mode.first << "\n";
            ++i;
        }
    }

    auto e =  exporter( _mesh=mesh );

    if ( e->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";
        int i = 0;
        for( auto const& mode: modes )
        {
            e->add( ( boost::format( "mode-%1%" ) % i++ ).str(), mode.second );
        }

        e->save();
        LOG(INFO) << "exportResults done\n";
    }

}

int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="ge_laplacian",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    Application app;

    //app.add( new EigenProblem<2,2>() );
    app.add( new EigenProblem<3,2>() );
    app.run();
}
