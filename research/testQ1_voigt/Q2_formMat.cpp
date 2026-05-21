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
        ("Order", po::value<int>()->default_value(1), "Hexahedron Lagrange's order")
        ("ForceApply_Point", po::value<bool>()->default_value(false), "Apply forces on points")
        ("Dirichlet_Point", po::value<bool>()->default_value(false), "Apply Dirichlet boundary condition on points")
        ;
    return options;
}



int main(int argc, char **argv)
{
    try
    {
        Environment env(_argc = argc, _argv = argv, _desc = makeOptions());

        
        // // ============== lit les paramètres du cfg ==============

        auto mesh_file = Environment::expand( soption(_name = "gmsh.filename") );

        int Order = ioption(_name = "Order" );
        auto mesh = loadMesh(_mesh = new Mesh<Hypercube<3, 2>>(), _filename = mesh_file, _scale = 1, _straighten = false);
        
        double E_young = doption(_name = "E" );
        double nu_poisson = doption(_name = "nu" );
        double lambda = E_young * nu_poisson / ((1 + nu_poisson) * (1 - 2 * nu_poisson));
        double mu = E_young / (2 * (1 + nu_poisson));

        auto f = expr<3,1>( soption(_name="functions.f"), "f" );
        auto g = expr<3,1>( soption(_name="functions.g"), "g" );



        // ============== Construction des espaces, initialisation ==============
        auto Vh = Pchv<1>( mesh );

        auto u = Vh->element();
        auto v = Vh->element();

        auto l = form1( _test = Vh );
        auto a = form2( _trial = Vh, _test = Vh );



        // ============== Construction et résolution du système ==============
        // inner( B(u), D, B(v) )
        auto B_utD = vec( (lambda+2*mu)*gradt(u)(0,0) + lambda* gradt(u)(1,1) + lambda* gradt(u)(2,2),
                          lambda*gradt(u)(0,0) + (lambda+2*mu)* gradt(u)(1,1) + lambda* gradt(u)(2,2),
                          lambda*gradt(u)(0,0) + lambda* gradt(u)(1,1) + (lambda+2*mu)* gradt(u)(2,2),
                          mu*( gradt(u)(1,0) + gradt(u)(0,1) ),
                          mu*( gradt(u)(2,0) + gradt(u)(0,2) ),
                          mu*( gradt(u)(2,1) + gradt(u)(1,2) )
                        );

        auto Eps_v = vec( grad(v)(0,0),
                          grad(v)(1,1),
                          grad(v)(2,2),
                          grad(v)(1,0) + grad(v)(0,1),
                          grad(v)(2,0) + grad(v)(0,2),
                          grad(v)(2,1) + grad(v)(1,2)
                        );

        a = integrate( _range = elements(mesh), _expr  = inner( B_utD, Eps_v ) );
        l = integrate(_range= markedfaces(mesh, "ForceApply"), _expr = inner( f, id(v) ));

        if( boption(_name = "Dirichlet_Point" ) )
            a += on( _range = markedpoints(mesh,"DirichletPoints"), _rhs=l, _element = u, _expr = g );  
        else
            a += on( _range = markedfaces(mesh,"Dirichlet"), _rhs=l, _element = u, _expr = g );  
        // a += on( _range = markedfaces(mesh,"XMoins"), _rhs=l, _element = u, _expr = vec( cst(0.0), cst(0.0), cst(0.0) ) );  
        // a += on( _range = markedfaces(mesh,"XPlus"), _rhs=l, _element = u, _expr = vec( cst(0.0), cst(0.25), cst(0.0) ) );  


        a.solve( _rhs = l, _solution = u );


        
        // ============== Post-traitement ==============
        // espace matriciel symétrique, stocké en vecteur de taille 6 mais considéré comme matrice de taille 3*3
        auto Ah = Pchms<1>( mesh );

        // évaluation du vecteur de dim 6 aux 8 noeuds géométriques (pas pareil que sur Matlab qu'il évalue aux points d'intégration)
        auto Eps_u_mat = mat<3,3>(
            gradv(u)(0,0),                    gradv(u)(0,1) + gradv(u)(1,0),  gradv(u)(0,2) + gradv(u)(2,0),
            gradv(u)(1,0) + gradv(u)(0,1) ,   gradv(u)(1,1),                  gradv(u)(1,2) + gradv(u)(2,1),
            gradv(u)(2,0) + gradv(u)(0,2),    gradv(u)(2,1) + gradv(u)(1,2),  gradv(u)(2,2)
        );
        // projection continue pour visualiser epsilon aux noeuds
        auto epsilon = project( _space = Ah, _range = elements(mesh), _expr = Eps_u_mat );


        auto DEps = mat<3,3>(
            (lambda+2*mu)* gradv(u)(0,0) + lambda* gradv(u)(1,1) + lambda* gradv(u)(2,2),  mu* (gradv(u)(0,1) + gradv(u)(1,0)),  mu* (gradv(u)(0,2) + gradv(u)(2,0)),
            mu* (gradv(u)(1,0) + gradv(u)(0,1)),  lambda* gradv(u)(0,0) + (lambda+2*mu)* gradv(u)(1,1) + lambda* gradv(u)(2,2),  mu* (gradv(u)(2,1) + gradv(u)(1,2)),
            mu* (gradv(u)(0,2) + gradv(u)(2,0)),  mu* (gradv(u)(2,1) + gradv(u)(1,2)),  lambda* gradv(u)(0,0) + lambda* gradv(u)(1,1) + (lambda+2*mu)* gradv(u)(2,2)
        );

        auto sigma = project( _space = Ah, _range = elements(mesh), _expr = DEps );



        // ============== Export de la solution ==============
        auto e = exporter( _mesh = mesh, _name = "q2voigt" );
        e->add( "displacement", u );
        e->add( "epsilon", epsilon );
        e->add( "sigma", sigma );
        e->save();

        l.vector().printMatlab( "form1.m" );
        a.matrix().printMatlab( "form2.m" );
        u.printMatlab( "solution.m" );
        // taille 8*6 (stocké en vecteur 6 mais considéré comme matrice 3*3) -> affichage noeud par noeud
        epsilon.printMatlab( "Epsilon.m" );             
        sigma.printMatlab( "Sigma.m" );             

    }
    catch (...)
    {
        handleExceptions();
    }
    return 0;
}