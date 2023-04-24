#include <feel/feelmor/crbplugin.hpp>
#include <feel/feeldiscr/sensors.hpp>
#include <eye2brain.hpp>
#include <iomanip>
namespace Feel
{

po::options_description
makeEye2BrainOptions()
{
    po::options_description eye2brainoptions( "Eye2Brain options" );
    // eye2brainoptions.add_options()
    //     ( "gamma", po::value<double>()->default_value( 10 ), "penalisation term" )
    //     ;
    return eye2brainoptions;
}
AboutData
makeEye2BrainAbout( std::string const& str )
{
    Feel::AboutData about( /*AppName  */ str.c_str(),
                           /*ProgName */ str.c_str(),
                           /*Version  */ "0.1",
                           /*ShortDesc*/ "Eye2Brain 3D Heat Application",
                           /*Licence  */ Feel::AboutData::License_GPL,
                           /*Copyright*/ "Copyright (c) 2023 Feel++ Consortium" );
    return about;
}

Eye2Brain::Eye2Brain()
    :
    super_type( "eye2brain_P1G1" )
{
    this->setPluginName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_NAME) + std::string("P1G1") );
    this->setPluginLibName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_LIBNAME) );
}


void
Eye2Brain::updateSpecificityModelPropertyTree( boost::property_tree::ptree & ptree ) const
{

}

void
Eye2Brain::initBetaQ()
{
    this->M_betaAq.resize( 5 );
    this->M_betaFq.resize( 2 );
    this->M_betaFq[0].resize( 2 );
    this->M_betaFq[1].resize( 1 );
}

Eye2Brain::super_type::betaq_type
Eye2Brain::computeBetaQ( parameter_type const& mu )
{
    //std::cout << "computeBetaQ start \n" << mu << std::endl;
    if ( this->M_betaAq.empty() )
        this->initBetaQ();

    for (int k = 0; k < 5; ++k)
        this->M_betaAq[k] = mu(k);
    for (int k = 0; k < 2; ++k)
        this->M_betaFq[0][k] = mu( 5 + k );
    this->M_betaFq[1][0] = 1.;
    //std::cout << "computeBetaQ finish \n";
    return boost::make_tuple( this->M_betaAq, this->M_betaFq );
}


void
Eye2Brain::initModel()
{

    CHECK( is_linear && !is_time_dependent ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&NonLinear) == NonLinear) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&TimeDependent) == TimeDependent) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";

    // "mu0": "k_lens:k_lens",
    // "mu1": "h_amb:h_amb",
    // "mu2": "h_bl:h_bl",
    // "mu3": "h_r:h_r",
    // "mu4": "1",
    // "mu5": "h_amb * T_amb + h_r * T_amb - E:h_amb:T_amb:E:h_r",
    // "mu6": "h_bl * T_bl:h_bl:T_bl"

    double E_min = 20        , E_max = 320,
           h_amb_min = 8     , h_amb_max = 100,
           h_bl_min = 50     , h_bl_max = 110,
           h_r_min = 6       , h_r_max = 6,
           T_amb_min = 283.15, T_amb_max = 303.15,
           T_bl_min = 308.3  , T_bl_max = 312,
           k_lens_min = 0.21 , k_lens_max = 0.544;

    double k_lens_ref = 0.4, h_amb_ref = 10, h_bl_ref = 65, h_r_ref = 6;

    Dmu->setDimension( 7 );
    auto mu_min = Dmu->element();
    mu_min << k_lens_min, h_amb_min, h_bl_min, h_r_min, 1, h_amb_min*T_amb_min + h_r_min*T_amb_min - E_max, h_bl_min*T_bl_min;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << k_lens_max, h_amb_max, h_bl_max, h_r_max, 1, h_amb_max*T_amb_max + h_r_max*T_amb_max - E_min, h_bl_max*T_bl_max;
    Dmu->setMax( mu_max );
    

    auto mesh = loadMesh( _mesh = new Eye2BrainConfig::mesh_type,
                          _update = size_type(MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES),
                          _savehdf5 = 0 );

    this->setFunctionSpaces( Eye2BrainConfig::space_type::New( mesh ) );
    this->setSymbolicExpressionBuildDir("$repository/crb/eye2brain/symbolicexpr/");

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << this->Xh->nDof() << "\n";
    }

    auto u = Xh->element();
    auto v = Xh->element();


    std::vector<double> muRef = {k_lens_ref, h_amb_ref, h_bl_ref, h_r_ref, 1};

    auto energy = backend()->newMatrix( _test = Xh, _trial = Xh );

    auto a0 = form2( _trial = Xh, _test = Xh );
    a0 = integrate( _range = markedelements(mesh, "Lens"), _expr = gradt( u ) * trans( grad( v ) )  );
    a0.matrixPtr()->close();
    this->addLhs( { a0 , "mu0" } );
    energy->addMatrix(muRef[0], a0.matrixPtr() );

    auto a1 = form2( _trial = Xh, _test = Xh );
    a1 = integrate( _range = markedfaces(mesh, {"BC_Cornea"}), _expr = idt( u )  * id( v ) );
    a1.matrixPtr()->close();
    this->addLhs( { a1 , "mu1" } );
    energy->addMatrix(muRef[1], a1.matrixPtr() );

    auto a2 = form2( _trial = Xh, _test = Xh );
    a2 = integrate( _range = markedfaces(mesh, {"BC_Sclera", "BC_OpticNerve"}), _expr = idt( u )  * id( v ) );
    a2.matrixPtr()->close();
    this->addLhs( { a2 , "mu2" } );
    energy->addMatrix(muRef[2], a2.matrixPtr() );

    auto a3 = form2( _trial = Xh, _test = Xh );
    a3 = integrate( _range = markedfaces(mesh, "BC_Cornea"), _expr = idt( u )  * id( v ));
    a3.matrixPtr()->close();
    this->addLhs( { a3 , "mu3" } );
    energy->addMatrix(muRef[3], a3.matrixPtr() );

    auto a4 = form2( _trial = Xh, _test = Xh );
    std::map < std::string, double > regions = { {"Cornea", 0.58}, {"Sclera", 1.0042}, {"AqueousHumor", 0.28}, {"VitreousHumor", 0.603}, {"Iris", 1.0042}, {"Lamina", 1.0042}, {"Choroid", 0.52}, {"Retina", 0.52}, {"OpticNerve", 1.0042} };
    for (auto const& [key, value] : regions)
    {
        a4 += integrate( _range = markedelements(mesh, key), _expr = value * gradt( u ) * trans( grad( v ) ));
    }
    a4.matrixPtr()->close();
    this->addLhs( { a4 , "mu4" } );
    energy->addMatrix(muRef[4], a4.matrixPtr() );

    auto f0 = form1( _test = Xh );
    f0 = integrate( _range = markedfaces( mesh, "BC_Cornea" ), _expr = id( v ) );
    f0.vectorPtr()->close();
    this->addRhs( { f0, "mu5" } );

    auto f1 = form1( _test = Xh );
    f1 = integrate( _range = markedfaces( mesh, {"BC_Sclera", "BC_OpticNerve" } ), _expr = id( v ) );
    f1.vectorPtr()->close();
    this->addRhs( { f1, "mu6" } );

    /// [energy]
    //a0.matrixPtr()->symmetricPart( energy );
    //energy->addMatrix(1.,a0.matrixPtr() );
    //energy->addMatrix(1./mymu,a0.matrixPtr() );

    energy->close();
    this->addEnergyMatrix( energy );

    /// [output]
#if 1
    auto out1 = form1( _test = Xh );
   
    double meas = integrate( _range = markedelements(mesh, "Cornea"), _expr = cst(1.) ).evaluate()(0,0);
    out1 = integrate( _range = markedelements(mesh, "Cornea"), _expr = id( u )/cst(meas)) ;
#else
    std::vector<double> coord = {-0.013597, 0, 0};
    node_type n(Eye2BrainConfig::space_type::nDim);
    for( int i = 0; i < Eye2BrainConfig::space_type::nDim; ++i )
        n(i) = coord[i];
    auto s = std::make_shared<SensorPointwise<space_type>>(Xh, n, "O");
    auto out1 = form1(_test=Xh,_vector=s->containerPtr());
#endif

    this->addOutput( { out1, "1" } );
}

double
Eye2Brain::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{
    //CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    auto mesh = Xh->mesh();
    double output=0;
    // right hand side (compliant)
    if ( output_index == 0 )
    {
        //output  = integrate( markedfaces( mesh,"BR" ), -mu(0)*expr(soption("functions.f"))*(gradv(u)*N()+doption("gamma")*idv(u)) ).evaluate()(0,0)
        //    + integrate( markedfaces( mesh,"BL" ), -mu(1)*expr(soption("functions.g"))*(gradv(u)*N()+doption("gamma")*idv(u)) ).evaluate()(0,0);
    }
    else if ( output_index == 1 )
    {
#if 1
        output = mean(_range = markedelements(mesh, "Cornea"), _expr = idv(u))(0,0);
#else
        std::vector<double> coord = {-0.013597, 0, 0};
        node_type n(Eye2BrainConfig::space_type::nDim);
        for( int i = 0; i < Eye2BrainConfig::space_type::nDim; ++i )
            n(i) = coord[i];
        auto s = std::make_shared<SensorPointwise<space_type>>(Xh, n, "O");
        auto out1 = form1(_test=Xh,_vector=s->containerPtr());
        out1 = integrate( _range = markedelements(mesh, "Cornea"), _expr = id( u ) );
        out1.vectorPtr()->close();
        output = out1.vectorPtr()->operator()(0);
#endif

    if ( Environment::isMasterRank() )
        std::cout << " Eye2Brain::output " << std::setprecision(16) << output << "\n";

    }
    // else if ( output_index == 2 )
    // {
    //     output = mean(elements(mesh),idv(u)).evaluat()(0,0);
    // }
    else
        throw std::logic_error( "[Eye2Brain::output] error with output_index : only 0 or 1 " );
    return output;

}


FEELPP_CRB_PLUGIN( Eye2Brain, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,P1G1) )

} // namespace Feel
