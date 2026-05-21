// #include <feel/feeldiscr/mesh.hpp>
// #include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/loadmesh.hpp>
// #include <feel/feeldiscr/product.hpp>
// #include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pchm.hpp>
// #include <feel/feeldiscr/pdhv.hpp>
#include <feel/feeldiscr/pdhm.hpp>

// #include <feel/feeldiscr/pdh.hpp>


using namespace Feel;
using namespace vf;


inline po::options_description makeOptions()
{
    po::options_description options( "SB9 implementation's options" );
    options.add_options()
        // ("specs", po::value<std::string>()->default_value("undefined"), "json spec file for the simulation")
        // ("gmsh.filename", po::value<std::string>()->default_value("undefined"), "mesh file")
        ("E", po::value<double>()->default_value(1.0), "Young's modulus")
        ("nu", po::value<double>()->default_value(1.0), "Poisson's ratio")
        ("Order", po::value<int>()->default_value(1), "Hexahedron Lagrange's order")
        ("ForceApply_Point", po::value<bool>()->default_value(false), "Apply forces on points")
        ("Dirichlet_Point", po::value<bool>()->default_value(false), "Apply Dirichlet boundary condition on points")
        ( "moment_x", po::value<bool>()->default_value( false ), "Moment x test" )
        ( "test_partial", po::value<bool>()->default_value( false ), "Tests with partial Dirichlet conditions" )
        ( "force_on_point", po::value<bool>()->default_value( false ), "Tests with partial Dirichlet conditions" )
        ;
    return options; //.add( feel_options() );
}



int main(int argc, char **argv)
{
    try
    {
        Environment env(_argc = argc, _argv = argv, _desc = makeOptions());

        
        // // ============== lit les paramètres du cfg ==============

        auto mesh_file = Environment::expand( soption(_name = "gmsh.filename") );

        // const uint16_type Order = static_cast<uint16_type>( ioption(_name = "Order" ) );
        int Order = ioption(_name = "Order" );
        // if(Order==1)
        //     auto mesh = loadMesh(_mesh = new Mesh<Hypercube<3, 1 >>(), _filename = mesh_file, _scale = 1, _straighten = false);
        // else if(Order==2)
        //     auto mesh = loadMesh(_mesh = new Mesh<Hypercube<3, 2 >>(), _filename = mesh_file, _scale = 1, _straighten = false);

        auto mesh = loadMesh(_mesh = new Mesh<Hypercube<3, 1>>(), _filename = mesh_file, _scale = 1, _straighten = false);   // doit être def à la compilation
        

        double E_young = doption(_name = "E" );
        std::cout << "E = " << E_young << std::endl;

        double nu_poisson = doption(_name = "nu" );
        std::cout << "nu = " << nu_poisson << std::endl;

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


        // // ============== Construction et résolution du système ==============
        // auto C = isotropic_stiffness<3>( lambda, mu );

        // // pour le produit scalaire vaut mieux écrire en Mandel -> std::sqrt(2)*voigt_ij ?
        // auto Eps_u = mandel_vec( gradt(u)(0,0),
        //                std::sqrt(2)* gradt(u)(1,0),
        //                std::sqrt(2)* gradt(u)(2,0),
        //                gradt(u)(1,0) + gradt(u)(0,0),
        //                std::sqrt(2)*( gradt(u)(2,0) + gradt(u)(0,0) ),
        //                gradt(u)(2,0) + gradt(u)(1,0)
        //             );
        // auto Eps_v = mandel_vec( grad(v)(0,0),
        //                std::sqrt(2)* grad(v)(1,0),
        //                std::sqrt(2)* grad(v)(2,0),
        //                grad(v)(1,0) + grad(v)(0,0),
        //                std::sqrt(2)*( grad(v)(2,0) + grad(v)(0,0) ),
        //                grad(v)(2,0) + grad(v)(1,0)
        //             );

        // a += integrate( _range = elements( mesh ), _expr = contract( D, Eps_u, Eps_v ) );
        // l += integrate( _range = elements( mesh ), _expr = cst(1.0) );


        // Sans notation de Voigt:
        // auto B_utD = vec( (lambda+2*mu)*gradt(u)(0,0) + lambda* gradt(u)(1,0) + lambda* gradt(u)(2,0),
        //                   lambda*gradt(u)(0,0) + (lambda+2*mu)* gradt(u)(1,0) + lambda* gradt(u)(2,0),
        //                   lambda*gradt(u)(0,0) + lambda* gradt(u)(1,0) + (lambda+2*mu)* gradt(u)(2,0),
        //                   mu*( gradt(u)(1,0) + gradt(u)(0,0) ),
        //                   mu*( gradt(u)(2,0) + gradt(u)(0,0) ),
        //                   mu*( gradt(u)(2,0) + gradt(u)(1,0) )
        //                 );
        auto B_utD = vec( (lambda+2*mu)*gradt(u)(0,0) + lambda* gradt(u)(1,1) + lambda* gradt(u)(2,2),
                          lambda*gradt(u)(0,0) + (lambda+2*mu)* gradt(u)(1,1) + lambda* gradt(u)(2,2),
                          lambda*gradt(u)(0,0) + lambda* gradt(u)(1,1) + (lambda+2*mu)* gradt(u)(2,2),
                          mu*( gradt(u)(1,0) + gradt(u)(0,1) ),
                          mu*( gradt(u)(2,0) + gradt(u)(0,2) ),
                          mu*( gradt(u)(2,1) + gradt(u)(1,2) )
                        );

        // auto Eps_v = vec( grad(v)(0,0),
        //                grad(v)(1,0),
        //                grad(v)(2,0),
        //                grad(v)(1,0) + grad(v)(0,0),
        //                 grad(v)(2,0) + grad(v)(0,0),
        //                grad(v)(2,0) + grad(v)(1,0)
        //             );
        auto Eps_v = vec( grad(v)(0,0),
                grad(v)(1,1),
                grad(v)(2,2),
                grad(v)(1,0) + grad(v)(0,1),
                grad(v)(2,0) + grad(v)(0,2),
                grad(v)(2,1) + grad(v)(1,2)
            );

        // auto Eps_u = vec( gradt(u)(0,0),
        //         gradt(u)(1,1),
        //         gradt(u)(2,2),
        //         gradt(u)(1,0) + gradt(u)(0,1),
        //         gradt(u)(2,0) + gradt(u)(0,2),
        //         gradt(u)(2,1) + gradt(u)(1,2)
        //     );

        // auto Eps_u_mat = mat<3,3>( gradt(u)(0,0),                  gradt(u)(1,0) + gradt(u)(0,1),  gradt(u)(2,0) + gradt(u)(0,2),
        //                            gradt(u)(1,0) + gradt(u)(0,1),  gradt(u)(1,1),                  gradt(u)(2,1) + gradt(u)(1,2),
        //                            gradt(u)(2,0) + gradt(u)(0,2),  gradt(u)(2,1) + gradt(u)(1,2),  gradt(u)(2,2)
        //                         );



        a = integrate( _range = elements(mesh), _expr  = inner( B_utD, Eps_v ) );

        // Appliquer une force sur un point -> direct changer form1 ddl par ddl
        // if( boption(_name = "ForceApply_Point" ) ) 
        //     // l = integrate(_range= markedpoints(mesh, "ForceApplyPoints"), _expr = inner( f, id(v) ));
        //     l = integrate(_range= markedpoints(mesh, "ForceApplyPoints"), _expr = trans( f )*id(v) );
        // else     
        
        if( boption(_name = "moment_x" ) ) {
            auto force = vec( cst(0.), -6*Pz(), 6*(Py() - 0.5) );
            l = integrate(_range= markedfaces(mesh, "ForceApply"), _expr = inner( force, id(v) ));
        }
        else if( boption(_name = "force_on_point")) {
            auto lvec = l.vectorPtr();
            (*lvec)(0) = 123.0;  // condition de Dirichlet sur le 1e point donc supprime la force appliquée
            (*lvec)(1) = 123.0;
            (*lvec)(2) = 123.0;
            (*lvec)(3) = 123.0;

        }
        else
            l = integrate(_range= markedfaces(mesh, "ForceApply"), _expr = inner( f, id(v) ));



            
        if( boption(_name = "Dirichlet_Point" ) )
            if( boption(_name = "test_partial" ) ) {
                // peut etre avec idv(u)(0,0) pour x 1,1 pour y et 2,2 pour z
                // idv à cette étape vaut encore 0 -> utiliser idt(u) ?
                // a += on( _range = markedpoints(mesh,"DirichletPoint1"), _rhs=l, _element = u, _expr = vec(cst(0.0), cst(0.0), cst(0.0)) ); 
                // a += on( _range = markedpoints(mesh,"DirichletPoint2"), _rhs=l, _element = u, _expr = vec(cst(0.0), idt(u)(1,1), cst(0.0)) ); 
                // a += on( _range = markedpoints(mesh,"DirichletPoint3"), _rhs=l, _element = u, _expr = vec(cst(0.0), idt(u)(1,1), idt(u)(2,2)) ); 
                // a += on( _range = markedpoints(mesh,"DirichletPoint4"), _rhs=l, _element = u, _expr = vec(cst(0.0), cst(0.0), idt(u)(2,2)) ); 

                // .component .comp *Nx()  avec inner(u, Nx())
                // a += on( _range = markedpoints(mesh,"DirichletPoint2"), _rhs=l, _element = u(0), _expr = vec(cst(0.0), cst(0.0), cst(0.0)) ); 

            }
            else
                a += on( _range = markedpoints(mesh,"DirichletPoints"), _rhs=l, _element = u, _expr = g );  
        else
            a += on( _range = markedfaces(mesh,"Dirichlet"), _rhs=l, _element = u, _expr = g );  


        a.solve( _rhs = l, _solution = u );



        
        // ============== Post-traitement ==============

        // auto Eps_u = vec( gradt(u)(0,0),
        //         gradt(u)(1,1),
        //         gradt(u)(2,2),
        //         gradt(u)(1,0) + gradt(u)(0,1),
        //         gradt(u)(2,0) + gradt(u)(0,2),
        //         gradt(u)(2,1) + gradt(u)(1,2)
        //     );
        // auto epsilon = form1( _test = Vh );
        // // epsilon = integrate( _range = elements(mesh), _expr = inner(Eps_u, id(v)) );
        // epsilon = on( _range = elements(mesh), _expr = Eps_u );  // moi je veux évaluer Eps_u sur chaque élément de mon maillage mais le pb c'est que on est défini avec form2 pas form1

        // u appartient à Vh, epsilon est en fonction de u et de taille 6 -> dimension (6,dim Vh) ah mais c'est plutot (6,dim Vh / 3) vu que c'est par rapport à juste une composante
        // produit de fonction space ? mais pas vraiment besoin de faire de projection ?   produit de Pch/Pdh de taille 6 ou Vh de taille 6 qui est lui même de taille 1 (Vh/3) ?


        // Eps de taille 6*1 ; vu que chaque ligne de eps il est en fonction de 1 composante de u
        // donc espace discret vectoriel de dim 6  ; si je fais Pdhv<6>(mesh)
        // auto Ah = Pdhv<6>(mesh);  




        // auto Ah = Pdhms<1>( mesh );    // espace matriciel symétrique discret

        // // mais Eps_u est def sur Vh, donc on doit projeter Eps_u sur Ah
        // auto epsilon = project( _space = Ah, _range = elements(mesh), _expr = Eps_u_mat );   

        auto Ah = Pchms<1>( mesh );

        // évaluation aux points d'intégrations (pas les mêmes que sur Matlab) mais visualisation sur les 8 noeuds géométriques
        auto Eps_u_mat = mat<3,3>(
            gradv(u)(0,0),                    gradv(u)(0,1) + gradv(u)(1,0),  gradv(u)(0,2) + gradv(u)(2,0),
            gradv(u)(1,0) + gradv(u)(0,1) ,   gradv(u)(1,1),                  gradv(u)(1,2) + gradv(u)(2,1),
            gradv(u)(2,0) + gradv(u)(0,2),    gradv(u)(2,1) + gradv(u)(1,2), gradv(u)(2,2)
        );
        // projection continue pour visualiser epsilon aux noeuds
        auto epsilon = project( _space = Ah, _range = elements(mesh), _expr = Eps_u_mat );


        

        // auto DEps = mat<3,3>(
        //     (lambda+2*mu)* gradv(u)(0,0) + lambda* gradv(u)(1,1) + lambda* gradv(u)(2,2),  mu* (gradv(u)(0,1) + gradv(u)(1,0)),  mu* (gradv(u)(0,2) + gradv(u)(2,0)),
        //     mu* (gradv(u)(2,1) + gradv(u)(1,2)),  lambda* gradv(u)(0,0) + (lambda+2*mu)* gradv(u)(1,1) + lambda* gradv(u)(2,2),  mu* (gradv(u)(2,1) + gradv(u)(1,2)),
        //     mu* (gradv(u)(0,2) + gradv(u)(2,0)),  mu* (gradv(u)(2,1) + gradv(u)(1,2)),  lambda* gradv(u)(0,0) + lambda* gradv(u)(1,1) + (lambda+2*mu)* gradv(u)(2,2)
        // );

        auto DEps = mat<3,3>(
            (lambda+2*mu)* gradv(u)(0,0) + lambda* gradv(u)(1,1) + lambda* gradv(u)(2,2),  mu* (gradv(u)(0,1) + gradv(u)(1,0)),  mu* (gradv(u)(0,2) + gradv(u)(2,0)),
            mu* (gradv(u)(1,0) + gradv(u)(0,1)),  lambda* gradv(u)(0,0) + (lambda+2*mu)* gradv(u)(1,1) + lambda* gradv(u)(2,2),  mu* (gradv(u)(2,1) + gradv(u)(1,2)),
            mu* (gradv(u)(0,2) + gradv(u)(2,0)),  mu* (gradv(u)(2,1) + gradv(u)(1,2)),  lambda* gradv(u)(0,0) + lambda* gradv(u)(1,1) + (lambda+2*mu)* gradv(u)(2,2)
        );

        auto sigma = project( _space = Ah, _range = elements(mesh), _expr = DEps );



        // ============== Export de la solution ==============
        auto e = exporter( _mesh = mesh, _name = "q1voigt" );
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