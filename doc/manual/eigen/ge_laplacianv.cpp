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
#if defined(USE_LM)
    typedef bases<Lagrange<Order,Vectorial>,Lagrange<Order-1,Scalar>,Lagrange<0,Scalar,Continuous> > basis_type;
#else
    typedef bases<Lagrange<Order,Vectorial>,Lagrange<Order-1,Scalar> > basis_type;
#endif
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
    auto U = Xh->element();
    auto u = U.template element<0>();
    auto p = U.template element<1>();
#if defined(USE_LM)
    auto lambda = U.template element<2>();
#endif
    auto V = Xh->element();
    auto v = U.template element<0>();
    auto q = U.template element<1>();
#if defined(USE_LM)
    auto nu = V.template element<2>();
#endif

    auto l = form1( _test=Xh );
    auto a = form2( _test=Xh, _trial=Xh);
    a = integrate( elements( mesh ), trace(gradt(u)*trans(grad(v))) );
    a += integrate( elements( mesh ), -idt(p)*div(v) );
    a += integrate( elements( mesh ), -divt(u)*id(q) );
#if defined(USE_LM)
    a += integrate( elements( mesh ), idt(lambda)*id(q) );
    a += integrate( elements( mesh ), id(nu)*idt(p) );
#endif
    auto beta = doption(_name="parameters.beta");
    //a += integrate( elements( mesh ), beta*idt(p)*id(q) );
    auto gamma = doption(_name="parameters.gamma");
    //a += integrate( boundaryfaces(mesh), -trans(-idt(p)*N()+gradt(u)*N())*id(v) -trans(-id(q)*N()+grad(u)*N())*idt(v)  + gamma*(trans(idt(u))*N())*(trans(id(u))*N())/hFace() );
    a += integrate( boundaryfaces(mesh), gamma*(trans(idt(u))*N())*(trans(id(u))*N())/hFace() );
    //a+= on( boundaryfaces(mesh), _element=u, _rhs=l, _expr=cst(0.));

    auto b = form2( _test=Xh, _trial=Xh);
    b = integrate( elements(mesh), trans(idt( u ))*id( v ) );
    b += integrate( elements( mesh ), beta*idt(p)*id(q) );
#if defined(USE_LM)
    b += integrate( elements( mesh ), beta*idt(lambda)*id(nu) );
#endif

    int nev = ioption(_name="solvereigen.nev");
    int ncv = ioption(_name="solvereigen.ncv");;

    double eigen_real, eigen_imag;

    SolverEigen<double>::eigenmodes_type modes;

    if ( Environment::worldComm().isMasterRank() )
    {
        std::cout << "nev= " << nev <<std::endl;
        std::cout << "ncv= " << ncv <<std::endl;
    }

    modes=
    eigs( _matrixA=a.matrixPtr(),
          _matrixB=b.matrixPtr(),
          _nev=nev,
          _ncv=ncv,
          _transform=SINVERT,
          _spectrum=SMALLEST_MAGNITUDE,
          _verbose = true );

    auto femodes = std::vector<decltype( Xh->element() )>( modes.size(), Xh->element() );

    if ( !modes.empty() )
    {
        LOG(INFO) << "eigenvalue " << 0 << " = (" << modes.begin()->second.get<0>() << "," <<  modes.begin()->second.get<1>() << ")\n";

        int i = 0;
        for( auto const& mode : modes )
        {
            std::cout << " -- eigenvalue " << i << " = (" << mode.second.get<0>() << "," <<  mode.second.get<1>() << ")\n";
            femodes[i] = *mode.second.get<2>();
            double l2div = normL2(_range=elements(mesh),_expr=divv(femodes[i].template element<0>() ));
            if ( Environment::worldComm().isMasterRank() )
            {
                std::cout << "  - div = " <<  l2div << "\n";
            }
            ++i;
        }
    }

    auto e =  exporter( _mesh=mesh );

    if ( e->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";
        int i = 0;
        for( auto const& mode: femodes )
        {
            e->add( ( boost::format( "mode-u-%1%" ) % i ).str(), mode.template element<0>() );
            e->add( ( boost::format( "mode-p-%1%" ) % i++ ).str(), mode.template element<1>() );
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
                     _about=about(_name="ge_laplacianv",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    Application app;

    app.add( new EigenProblem<2,2>() );
    //app.add( new EigenProblem<3,2>() );
    app.run();
}





