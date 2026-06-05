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

        
        // // ============== lit les paramètres du cfg ==============
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

        auto u = Vh->element();
        auto v = Vh->element();

        auto l = form1( _test = Vh );
        auto a = form2( _trial = Vh, _test = Vh );



        // ============== Construction et résolution du système ==============
        // auto C = isotropic_stiffness<3, SymmetricTensorNotation::Voigt>( lambda, mu );   // par défaut c'est notation Mandel

        // // utiliser voigt_vec ou mandel_vec ? -> Mandel mieux pour simplicité du produit scalaire
        // // si on fait en voigt garder comme ça sinon en mandel ya un facteur 1/sqrt(2) devant les termes extra-diagonaux ?

        // // je crois dans feelpp les vec voigt stockent au format symétrique et pas voigt classique - à voir - non
        
        // -> stocké [11, 22, 33, 12, 13, 31] et pas [11, 22, 33, 2*12, 2*13, 2*31]
        auto Eps_u_voigt = voigt_vec<3>( gradt(u)(0,0),
                          gradt(u)(1,1),
                          gradt(u)(2,2),
                          cst(0.5)* (gradt(u)(1,0) + gradt(u)(0,1)),
                          cst(0.5)* (gradt(u)(2,0) + gradt(u)(0,2)),
                          cst(0.5)* (gradt(u)(2,1) + gradt(u)(1,2))
                        );
        auto Eps_v_voigt = voigt_vec<3>( grad(v)(0,0),
                          grad(v)(1,1),
                          grad(v)(2,2),
                          cst(0.5)* (grad(v)(1,0) + grad(v)(0,1)),
                          cst(0.5)* (grad(v)(2,0) + grad(v)(0,2)),
                          cst(0.5)* (grad(v)(2,1) + grad(v)(1,2))
                        );   


        // constexpr auto sqrt2 = std::numbers::sqrt2_v<double>;
        // avec cst(1/sqrt2) ca me donne 0.63 et sans rien ca me donne 0.61 -> comme sans rien en voigt (au lieu de 0.64)
        // essayer avec 0.5 devant 
        // auto Eps_u_mandel = mandel_vec<3>( gradt(u)(0,0),
        //                   gradt(u)(1,1),
        //                   gradt(u)(2,2),
        //                   cst(1/sqrt2)* (gradt(u)(1,0) + gradt(u)(0,1)),
        //                   cst(1/sqrt2)* (gradt(u)(2,0) + gradt(u)(0,2)),
        //                   cst(1/sqrt2)* (gradt(u)(2,1) + gradt(u)(1,2))
        //                 );
        // auto Eps_v_mandel = mandel_vec<3>( grad(v)(0,0),
        //                   grad(v)(1,1),
        //                   grad(v)(2,2),
        //                   cst(1/sqrt2)* (grad(v)(1,0) + grad(v)(0,1)),
        //                   cst(1/sqrt2)* (grad(v)(2,0) + grad(v)(0,2)),
        //                   cst(1/sqrt2)* (grad(v)(2,1) + grad(v)(1,2))
        //                 );
        auto Eps_u_mandel = mandel_vec<3>( gradt(u)(0,0),
                          gradt(u)(1,1),
                          gradt(u)(2,2),
                          cst(0.5)* (gradt(u)(1,0) + gradt(u)(0,1)),
                          cst(0.5)* (gradt(u)(2,0) + gradt(u)(0,2)),
                          cst(0.5)* (gradt(u)(2,1) + gradt(u)(1,2))
                        );
        auto Eps_v_mandel = mandel_vec<3>( grad(v)(0,0),
                          grad(v)(1,1),
                          grad(v)(2,2),
                          cst(0.5)* (grad(v)(1,0) + grad(v)(0,1)),
                          cst(0.5)* (grad(v)(2,0) + grad(v)(0,2)),
                          cst(0.5)* (grad(v)(2,1) + grad(v)(1,2))
                        );
        
        


        // vecVoigt.vector().printMatlab( "vecVoigt.m" );             
        // std::cout << "Vecteur en notation de Voigt : " << vecVoigt(0) << ", " << vecVoigt(1) << ", " << vecVoigt(2) << ", " << vecVoigt(3) << ", " << vecVoigt(4) << ", " << vecVoigt(5) << std::endl;



        auto vecVoigt = voigt_vec<3>( cst(1.0), cst(2.0), cst(3.0), cst(4.0), cst(5.0), cst(6.0));
        auto devoigt = unvoigt( vecVoigt );

        auto matAttendue = mat<3, 3>( cst( 1.0 ), cst( 4.0 ), cst( 5.0 ),
                                      cst( 4.0 ), cst( 2.0 ), cst( 6.0 ),
                                      cst( 5.0 ), cst( 6.0 ), cst( 3.0 ) );

        auto diff_devoigt = devoigt - matAttendue;

        std::cout << "Test unvoigt_vec : " << integrate( _range = elements( mesh ), _expr = inner( diff_devoigt, diff_devoigt ) ).evaluate()(0,0) << std::endl;



        auto vecMandel = mandel_vec<3>( cst(1.0), cst(2.0), cst(3.0), cst(4.0), cst(5.0), cst(6.0));
        auto deMandel = unmandel( vecMandel );

        auto matAttendue = mat<3, 3>( cst( 1.0 ), cst( 4.0 ), cst( 5.0 ),
                                      cst( 4.0 ), cst( 2.0 ), cst( 6.0 ),
                                      cst( 5.0 ), cst( 6.0 ), cst( 3.0 ) );

        auto diff_deMandel = deMandel - matAttendue;

        std::cout << "Test unMandel_vec : " << integrate( _range = elements( mesh ), _expr = inner( diff_deMandel, diff_deMandel ) ).evaluate()(0,0) << std::endl;





        // auto deft = sym(gradt(u));
        // auto Id = eye<3,3>();
        // auto sigmat = lambda*trace(deft)*Id + 2*mu*deft;
        // a = integrate( _range = elements(mesh), _expr = inner( sigmat, grad(v) ) );  // plutot sym(grad(v)) ?   -> ça fonctionne comme ça
        // a = integrate( _range = elements(mesh), _expr = lambda*trace(sym(gradt(u)))*trace(sym(grad(v))) + 2*mu*inner(sym(gradt(u)), sym(grad(v))));  // ok aussi
        // a = integrate( _range = elements(mesh), _expr = cst(lambda)*trace( sym(gradt(u)) )*trace( sym(grad(v)) ) + cst(2*mu)*inner( sym(gradt(u)), sym(grad(v)) ));

        // auto epsu = sym(gradt(u));    // symm_grad( u );
        // auto epsv = sym(grad(v));     // symm_grad( v );
        // a = integrate( _range = elements(mesh), _expr = cst( lambda )*trace( epsu )*trace( epsv ) + cst( 2.0*mu )*inner( epsu, epsv ));  // ok

        // auto C = isotropic_stiffness<3>( lambda, mu );
        // a = integrate( _range = elements(mesh), _expr = ddot( C, epsu, epsv ));   // ok

        // auto C = isotropic_stiffness<3, SymmetricTensorNotation::Voigt>( lambda, mu );
        // a = integrate( _range = elements(mesh), _expr = ddot<SymmetricTensorNotation::Voigt>( C, voigt(epsu), voigt(epsv) ));  // ok


        // NOTATION VOIGT
        auto C = isotropic_stiffness<3, SymmetricTensorNotation::Voigt>( lambda, mu );
        a = integrate( _range = elements(mesh), _expr = ddot<SymmetricTensorNotation::Voigt>( C, Eps_u_voigt, Eps_v_voigt ));


        // NOTATION MANDEL
        // auto C = isotropic_stiffness<3>( lambda, mu );       // Mandel par défaut
        // a = integrate( _range = elements(mesh), _expr = ddot( C, Eps_u_mandel, Eps_v_mandel ));    // pareil, Mandel par défaut


        // test
        // l = integrate( _range = markedfaces(mesh, "ForceApply"), _expr = inner( f, id(v) ));  
        l = integrate( _range = markedpoints(mesh, "Points"), _expr = inner( f, id(v) ));  

        a += on( _range = markedfaces(mesh,"Dirichlet"), _rhs=l, _element = u, _expr = g );   

        a.solve( _rhs = l, _solution = u );

        

        // ============== Post-traitement ==============
        auto Ah = Pchms<1>( mesh );

        auto Eps_u_mat = mat<3,3>(
            gradv(u)(0,0),                    gradv(u)(0,1) + gradv(u)(1,0),  gradv(u)(0,2) + gradv(u)(2,0),
            gradv(u)(1,0) + gradv(u)(0,1) ,   gradv(u)(1,1),                  gradv(u)(1,2) + gradv(u)(2,1),
            gradv(u)(2,0) + gradv(u)(0,2),    gradv(u)(2,1) + gradv(u)(1,2),  gradv(u)(2,2)
        );
        auto epsilon = project( _space = Ah, _range = elements(mesh), _expr = Eps_u_mat );


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
        epsilon.printMatlab( "Epsilon.m" );             
        sigma.printMatlab( "Sigma.m" );             

    }
    catch (...)
    {
        handleExceptions();
    }
    return 0;
}