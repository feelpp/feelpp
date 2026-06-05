#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feeldiscr/pdhm.hpp>



using namespace Feel;
using namespace vf;

inline po::options_description makeOptions()
{
    po::options_description options( "SB9 implementation's options" );
    options.add_options()
        ("E", po::value<double>()->default_value(1.0), "Young's modulus")
        ("nu", po::value<double>()->default_value(1.0), "Poisson's ratio")
        ;
    return options;
}



int main(int argc, char **argv)
{
    try
    {
        Environment env(_argc = argc, _argv = argv, _desc = makeOptions());

        
        // ============== lit les paramètres du cfg ==============
        auto mesh_file = Environment::expand( soption(_name = "gmsh.filename") );
        auto mesh = loadMesh(_mesh = new Mesh<Hypercube<3, 1>>(), _filename = mesh_file, _scale = 1, _straighten = false);

        double E_young = doption(_name = "E" );
        double nu_poisson = doption(_name = "nu" );
        double lambda = E_young * nu_poisson / ((1 + nu_poisson) * (1 - 2 * nu_poisson));
        double mu = E_young / (2 * (1 + nu_poisson));

        auto f = expr<3,1>( soption( _name = "functions.f"), "f" );
        auto g = expr<3,1>( soption( _name = "functions.g"), "g" );



        // ============== Construction des espaces, initialisation ==============
        auto Vh = Pchv<1>( mesh );

        auto u = Vh->element();    // 1
        auto v = Vh->element();    // 1

        // auto u = trial( Vh, "u" );    // 2
        // auto v = test( Vh, "v" );

        // auto uExplicit = Vh->element( "u" );
        // auto vExplicit = Vh->element( "v" );
        // auto u = trial( uExplicit );
        // auto v = test( vExplicit );



        auto l = form1( _test = Vh );
        l = integrate( _range = markedfaces(mesh, "ForceApply"), _expr = inner( f, id(v) ));  
        l.vector().printMatlab( "form1.m" );

        auto a = form2( _trial = Vh, _test = Vh );


        // ============== Construction et résolution du système ==============
        // auto epsu = sym(gradt(u));   // 1
        // auto epsv = sym(grad(v));    // 1
        // auto epsu = symm_grad(u);    // 2 mais erreur avec a+= et a.solve
        // auto epsv = symm_grad(v);  
        auto epsu = symm_gradt(u);    // 1   
        auto epsv = symm_grad(v);    // 1
        a = integrate( _range = elements(mesh), _expr = cst( lambda )*trace( epsu )*trace( epsv ) + cst( 2.0*mu )*inner( epsu, epsv ));
        // a.close();
        a += on( _range = markedfaces(mesh,"Dirichlet"), _rhs=l, _element = u, _expr = g ); 
        // a.close();

        a.solve( _rhs = l, _solution = u );
        std::cout << "L2 Norme de la solution attendue : " << u.l2Norm() << std::endl;
        a.matrix().printMatlab( "form2.m" );
        u.printMatlab( "solution.m" );


        // auto epsu_bug = symm_grad( u_bug );
        // auto epsv_bug = symm_grad( v );
        // a = integrate( _range = elements(mesh), _expr = cst( lambda )*trace( epsu_bug )*trace( epsv_bug ) + cst( 2.0*mu )*inner( epsu_bug, epsv_bug )); 
        // a += on( _range = markedfaces(mesh,"Dirichlet"), _rhs=l, _element = u_bug, _expr = g );   

        // a.solve( _rhs = l, _solution = u_bug );
        // std::cout << "L2 Norme de la solution incorrecte : " << u_bug.l2Norm() << std::endl;
        // a.matrix().printMatlab( "form2_bug.m" );
        // u_bug.printMatlab( "solution_bug.m" );



        // ============== Export de la solution ==============
        auto e = exporter( _mesh = mesh, _name = "q1bug" );
        e->add( "displacement", u );
        // e->add( "displacement_bug", u_bug );
        e->save();

    }
    catch (...)
    {
        handleExceptions();
    }
    return 0;
}