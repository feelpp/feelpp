#include <feel/feelmodels/hdg/poroelastic.hpp>

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "mixed-poisson-elasticity-model" ,
                     "mixed-poisson-elasticyt-model" ,
                     "0.1",
                     "Mixed-Poisson-Elasticity-Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Romain Hild", "developer", "", "" );
    about.addAuthor( "Lorenzo Sala", "developer", "", "" );
    return about;
}

template<int nDim, int OrderT>
void
runApplicationPoroelastic()
{
    using namespace Feel;

    typedef FeelModels::MixedPoissonElasticity<nDim,OrderT> mpe_type;

    auto mesh = loadMesh( _mesh=new typename mpe_type::mesh_type );

    decltype( IPtr( _domainSpace=Pdh<OrderT>(mesh), _imageSpace=Pdh<OrderT>(mesh) ) ) Idh_poi ;
    decltype( IPtr( _domainSpace=Pdhv<OrderT>(mesh), _imageSpace=Pdhv<OrderT>(mesh) ) ) Idhv_poi;

    decltype( IPtr( _domainSpace=Pdhv<OrderT>(mesh), _imageSpace=Pdhv<OrderT>(mesh) ) ) Idh_el ;
    decltype( IPtr( _domainSpace=Pdhms<OrderT>(mesh), _imageSpace=Pdhms<OrderT>(mesh) ) ) Idhv_el;

    std::list<std::string> listSubmesh;
    std::string unseparatedList = soption( "hdg.elasticity.gmsh.submesh");

    char help;
    std::string nameSubmesh;
    for ( int i = 0; i < unseparatedList.size(); i++)
    {

        help = unseparatedList[i];
        if ( help == ',' || i == unseparatedList.size()-1 )
        {
            if ( i ==  unseparatedList.size()-1)
                nameSubmesh.push_back(help);
            listSubmesh.push_back(nameSubmesh);
            nameSubmesh.erase();
        }
        else
        {
            nameSubmesh.push_back(help);
        }
    }


    if ( soption( "gmsh.submesh" ).empty() && soption("hdg.elasticity.gmsh.submesh").empty() )
    {
        mpe_type MPE( mesh, mesh, mesh, mesh);
        MPE.run( Idh_el, Idhv_el, Idh_poi, Idhv_poi );
    }
    else if ( soption( "gmsh.submesh" ).empty() )
    {
        Feel::cout << "Using submesh for Elasticity: " << soption("hdg.elasticity.gmsh.submesh") << std::endl;
        auto cmeshElasticity = createSubmesh( _mesh=mesh, _range=markedelements(mesh,listSubmesh ) );

        Idh_el = IPtr( _domainSpace=Pdhv<OrderT>(cmeshElasticity), _imageSpace=Pdhv<OrderT>(mesh) );
        Idhv_el = IPtr( _domainSpace=Pdhms<OrderT>(cmeshElasticity), _imageSpace=Pdhms<OrderT>(mesh) );

        mpe_type MPE( mesh, cmeshElasticity, cmeshElasticity, mesh );
        MPE.run( Idh_el, Idhv_el, Idh_poi, Idhv_poi );

    }
    else if ( soption("hdg.elasticity.gmsh.submesh").empty() )
    {
        Feel::cout << "Using submesh for Poisson: " << soption("gmsh.submesh") << std::endl;
        auto cmeshPoisson = createSubmesh( _mesh=mesh, _range=markedelements(mesh,soption("gmsh.submesh")) );

        Idh_poi = IPtr( _domainSpace=Pdh<OrderT>(cmeshPoisson), _imageSpace=Pdh<OrderT>(mesh) );
        Idhv_poi = IPtr( _domainSpace=Pdhv<OrderT>(cmeshPoisson), _imageSpace=Pdhv<OrderT>(mesh) );

        mpe_type MPE( cmeshPoisson, mesh, cmeshPoisson, mesh );
        MPE.run( Idh_el, Idhv_el, Idh_poi, Idhv_poi );
    }
    else
    {
        Feel::cout << "Using submesh for Poisson: " << soption("gmsh.submesh") << std::endl;
        Feel::cout << "Using submesh for Elasticity: " << soption("hdg.elasticity.gmsh.submesh") << std::endl;

        auto cmeshPoisson = createSubmesh( _mesh=mesh, _range=markedelements(mesh,soption("gmsh.submesh")) );
        auto cmeshElasticity = createSubmesh( _mesh=mesh, _range=markedelements(mesh,listSubmesh ) );

        Idh_poi = IPtr( _domainSpace=Pdh<OrderT>(cmeshPoisson), _imageSpace=Pdh<OrderT>(mesh) );
        Idhv_poi = IPtr( _domainSpace=Pdhv<OrderT>(cmeshPoisson), _imageSpace=Pdhv<OrderT>(mesh) );

        Idh_el = IPtr( _domainSpace=Pdhv<OrderT>(cmeshElasticity), _imageSpace=Pdhv<OrderT>(mesh) );
        Idhv_el = IPtr( _domainSpace=Pdhms<OrderT>(cmeshElasticity), _imageSpace=Pdhms<OrderT>(mesh) );

        auto meshCommon = (soption("gmsh.submesh")<soption("hdg.elasticity.gmsh.submesh")) ? cmeshPoisson : cmeshElasticity ;

        mpe_type MPE( cmeshPoisson, cmeshElasticity, meshCommon, mesh);
        MPE.run( Idh_el, Idhv_el, Idh_poi, Idhv_poi );
    }

}


int main(int argc, char *argv[])
{
    using namespace Feel;

    po::options_description peoptions( "poroelastic options" );
    peoptions.add( FeelModels::makeMixedPoissonElasticityOptions("poroelastic") );
    peoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ;

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=peoptions,
                           _desc_lib=FeelModels::makeMixedPoissonElasticityLibOptions("hdg.poroelasticity").add(feel_options())
                           );

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ),
                                             hana::make_tuple("P3", hana::int_c<3> ) );
#elif FEELPP_INSTANTIATION_ORDER_MAX >= 2
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ) );
#else
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ) );
#endif

    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension]( auto const& d )
                                                                                         {
                                                                                             constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                                                                                             std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                                                                                             constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                                                                                             if ( dimension == _dim && discretization == _discretization )
                                                                                                 runApplicationPoroelastic<_dim,_torder>();
                                                                                         } );

    return 0;
}
