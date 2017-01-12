#include <feel/feel.hpp>


using namespace Feel;

int
main( int argc, char** argv )
{
    // Initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _about=about( _name="myintegrals" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    /// [mesh]
    // create the mesh (specify the dimension of geometric entity)
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    /// [mesh]

    /// [expression]
    // our function to integrate
    auto g = expr( soption(_name="functions.g") );
    /// [expression]

    /// [integrals]
    // compute integral of g (global contribution): \(\int_{\partial \Omega} g\)
    auto intf_1 = integrate( _range = elements( mesh ),
                                 _expr = g ).evaluate();

    // compute integral g on boundary: \( \int_{\partial \Omega} g \)
    auto intf_2 = integrate( _range = boundaryfaces( mesh ),
                             _expr = g ).evaluate();

    // compute integral of grad f (global contribution): \( \int_{\Omega} \nabla g \)
    auto grad_g = grad<2>(g);
    auto intgrad_f = integrate( _range = elements( mesh ),
                                _expr = grad_g ).evaluate();

    // only the process with rank 0 prints to the screen to avoid clutter
    if ( Environment::isMasterRank() )
        std::cout << "int_Omega " << g << " = " << intf_1  << std::endl
                  << "int_{boundary of Omega} " << g << " = " << intf_2 << std::endl
                  << "int_Omega grad " << g << " = "
                  << "int_Omega  " << grad_g << " = "
                  << intgrad_f  << std::endl;
    /// [integrals]
}
/// [all]
