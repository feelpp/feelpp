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
    eye2brainoptions.add_options()
        ( "measure-index", po::value<int>()->default_value( 1 ), "index of the output to be computed.\n\t          0 : mean over cornea,\n\tFrom 1 to 9 : corresponding point in {O, A, B, B1, C, D, D1, F, G}\n\tWARNING : do not confuse with the option crb.output-index" )
        ;
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

template<int Order, int Dim>
Eye2Brain<Order, Dim>::Eye2Brain()
    :
    super_type( fmt::format("eye2brain_{}D_P{}", Dim, Order) )
{
    this->setPluginName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_NAME) + fmt::format("{}D_P{}", Dim, Order) );
    this->setPluginLibName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_LIBNAME) );
}

template<int Order, int Dim>
void
Eye2Brain<Order, Dim>::updateSpecificityModelPropertyTree( boost::property_tree::ptree & ptree ) const
{

}

template<int Order, int Dim>
void
Eye2Brain<Order, Dim>::initBetaQ()
{
    this->M_betaAq.resize( 5 );
    this->M_betaFq.resize( 2 );
    this->M_betaFq[0].resize( 2 );
    this->M_betaFq[1].resize( 1 );
}

template<int Order, int Dim>
typename Eye2Brain<Order, Dim>::super_type::betaq_type
Eye2Brain<Order, Dim>::computeBetaQ( parameter_type const& mu )
{
    // std::cout << "computeBetaQ start \n" << mu << std::endl;
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

template<int Order, int Dim>
void
Eye2Brain<Order, Dim>::initModel()
{

    CHECK( this->is_linear && !this->is_time_dependent ) << "Invalid model is_linear:" << this->is_linear << " is_time_dependent:" << this->is_time_dependent << "\n";
#if 1  
    LOG_IF( WARNING, ((this->Options & NonLinear) == NonLinear) ) << "Invalid model is_linear:" << this->is_linear << " is_time_dependent:" << this->is_time_dependent << "\n";
    LOG_IF( WARNING, ((this->Options & TimeDependent) == TimeDependent) ) << "Invalid model is_linear:" << this->is_linear << " is_time_dependent:" << this->is_time_dependent << "\n";
#endif

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

    this->Dmu->setDimension( 7 );
    auto mu_min = this->Dmu->element();
    mu_min << k_lens_min, h_amb_min, h_bl_min, h_r_min, 1, h_amb_min*T_amb_min + h_r_min*T_amb_min - E_max, h_bl_min*T_bl_min;
    this->Dmu->setMin( mu_min );
    auto mu_max = this->Dmu->element();
    mu_max << k_lens_max, h_amb_max, h_bl_max, h_r_max, 1, h_amb_max*T_amb_max + h_r_max*T_amb_max - E_min, h_bl_max*T_bl_max;
    this->Dmu->setMax( mu_max );

    LOG(INFO) << "[Eye2brain::initModel] Set parameter space : mu_min " << mu_min << "\n" << std::endl;
    LOG(INFO) << "[Eye2brain::initModel] Set parameter space : mu_max " << mu_max << "\n" << std::endl;


    auto mesh = loadMesh( _mesh = new typename Eye2BrainConfig<Order, Dim>::mesh_type,
                          _update = MESH_UPDATE_FACES | MESH_UPDATE_EDGES | MESH_NO_UPDATE_MEASURES,
                          _savehdf5 = 0 );

    this->setFunctionSpaces( Eye2BrainConfig<Order, Dim>::space_type::New( mesh ) );

    this->setSymbolicExpressionBuildDir("$repository/crb/eye2brain/symbolicexpr/");

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << this->Xh->nDof() << "\n";
    }

    auto u = this->Xh->element();
    auto v = this->Xh->element();

    std::vector<double> muRef = {k_lens_ref, h_amb_ref, h_bl_ref, h_r_ref, 1};
    auto energy = backend()->newMatrix( _test = this->Xh, _trial = this->Xh );


    auto a0 = form2( _trial = this->Xh, _test = this->Xh );
    a0 = integrate( _range = markedelements(mesh, "Lens"), _expr = gradt( u ) * trans( grad( v ) )  );
    a0.matrixPtr()->close();
    this->addLhs( { a0 , "mu0" } );
    energy->addMatrix(muRef[0], a0.matrixPtr() );

    auto a1 = form2( _trial = this->Xh, _test = this->Xh );
    a1 = integrate( _range = markedfaces(mesh, {"BC_Cornea"}), _expr = idt( u )  * id( v ) );
    a1.matrixPtr()->close();
    this->addLhs( { a1 , "mu1" } );
    energy->addMatrix(muRef[1], a1.matrixPtr() );

    auto a2 = form2( _trial = this->Xh, _test = this->Xh );
    a2 = integrate( _range = markedfaces(mesh, {"BC_Sclera", "BC_OpticNerve"}), _expr = idt( u )  * id( v ) );
    a2.matrixPtr()->close();
    this->addLhs( { a2 , "mu2" } );
    energy->addMatrix(muRef[2], a2.matrixPtr() );

    auto a3 = form2( _trial = this->Xh, _test = this->Xh );
    a3 = integrate( _range = markedfaces(mesh, "BC_Cornea"), _expr = idt( u )  * id( v ));
    a3.matrixPtr()->close();
    this->addLhs( { a3 , "mu3" } );
    energy->addMatrix(muRef[3], a3.matrixPtr() );

    auto a4 = form2( _trial = this->Xh, _test = this->Xh );
    std::map < std::string, double > regions = { {"Cornea", 0.58}, {"Sclera", 1.0042}, {"AqueousHumor", 0.28}, {"VitreousHumor", 0.603}, {"Iris", 1.0042}, {"Lamina", 1.0042}, {"Choroid", 0.52}, {"Retina", 0.52}, {"OpticNerve", 1.0042} };
    for (auto const& [key, value] : regions)
    {
        a4 += integrate( _range = markedelements(mesh, key), _expr = value * gradt( u ) * trans( grad( v ) ));
    }
    a4.matrixPtr()->close();
    this->addLhs( { a4 , "mu4" } );
    energy->addMatrix(muRef[4], a4.matrixPtr() );

    auto f0 = form1( _test = this->Xh );
    f0 = integrate( _range = markedfaces( mesh, "BC_Cornea" ), _expr = id( v ) );
    f0.vectorPtr()->close();
    this->addRhs( { f0, "mu5" } );

    auto f1 = form1( _test = this->Xh );
    f1 = integrate( _range = markedfaces( mesh, {"BC_Sclera", "BC_OpticNerve" } ), _expr = id( v ) );
    f1.vectorPtr()->close();
    this->addRhs( { f1, "mu6" } );

    /// [energy]
    energy->close();
    this->addEnergyMatrix( energy );

    /// [output]
    using form1_type = vf::detail::LinearForm<typename Eye2BrainConfig<Order, Dim>::space_type, typename backend_type::vector_type, typename backend_type::vector_type>;
    form1_type out1;
    int measure_index = ioption(_name = "measure-index");
    std::cout << "Measure index = " << measure_index << std::endl;
    if (measure_index >= 1 && measure_index <= 9)    // sensor output
    {
        std::string name = m_outputNames[measure_index-1];
        std::vector<double> coord = m_coordinates[measure_index-1];
        std::cout << "Output " << name << " at coord " << coord << std::endl;
        node_type n(Eye2BrainConfig<Order, Dim>::space_type::nDim);
        for( int i = 0; i < Eye2BrainConfig<Order,Dim>::space_type::nDim; ++i )
            n(i) = coord[i];
        Feel::cout << n << std::endl;
        auto s = std::make_shared<SensorPointwise<space_type>>(this->Xh, n, "O");
        out1 = form1(_test = this->Xh, _vector = s->containerPtr());
        out1.vectorPtr()->close();
    }
    else if ( measure_index == 0 )   // mean over cornea
    {
        out1 = form1( _test = this->Xh );

        double meas = integrate( _range = markedelements(mesh, "Cornea"), _expr = cst(1.) ).evaluate()(0,0);
        out1 = integrate( _range = markedelements(mesh, "Cornea"), _expr = id( u )/cst(meas)) ;
    }
    else
        throw std::logic_error( "[Eye2Brain::output] error with output_index : between 0 and 9 " );


    this->addOutput( { out1, "1" } );

}

template<int Order, int Dim>
double
Eye2Brain<Order, Dim>::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{
    auto mesh = this->Xh->mesh();
    std::string name;
    double output;

    if (output_index == 0) // compliant output
    {
        name = "compliant";
        /*
        output = integrate( _range = markedfaces( mesh, "BC_Cornea" )                    , _expr = mu(5) * id( u ) ).evaluate()(0,0)
               + integrate( _range = markedfaces( mesh, {"BC_Sclera", "BC_OpticNerve" } ), _expr = mu(6) * id( u ) ).evaluate()(0,0);
  */}
    else if ( output_index == 1 )
    {
        int measure_index = ioption(_name = "measure-index");
        if (measure_index == 0)                             // mean over cornea
        {
            name = "mean_cornea";
            output = mean(_range = markedelements(mesh, "Cornea"), _expr = idv(u))(0,0);
        }
        else if (measure_index >= 1 && measure_index <= 9)  // sensor output
        {
            name = m_outputNames[measure_index-1];
            std::vector<double> coord = m_coordinates[measure_index-1];
            node_type n(Eye2BrainConfig<Order, Dim>::space_type::nDim);
            for( int i = 0; i < Eye2BrainConfig<Order,Dim>::space_type::nDim; ++i )
                n(i) = coord[i];
            auto s = std::make_shared<SensorPointwise<space_type>>(this->Xh, n, name);
            auto out1 = form1(_test = this->Xh, _vector = s->containerPtr());
            out1.vectorPtr()->close();
            output = out1(u);
        }
        else
            throw std::logic_error( "[Eye2Brain::output] error with output_index : between 0 and 9 " );
    }
    else
        throw std::logic_error( "[Eye2Brain::output] error with output_index : only 0 or 1 " );

    if ( Environment::isMasterRank() )
        std::cout << " Eye2Brain::output " << std::setprecision(16) << output << "\n";

    return output;

}


template class Eye2Brain<1, 2>;
template class Eye2Brain<1, 3>;
template class Eye2Brain<2, 2>;
template class Eye2Brain<2, 3>;

FEELPP_CRB_PLUGIN_TEMPLATE( Eye2Brain2D_P1, Eye2Brain<1 BOOST_PP_COMMA() 2>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME, 2D_P1) )
FEELPP_CRB_PLUGIN_TEMPLATE( Eye2Brain3D_P1, Eye2Brain<1 BOOST_PP_COMMA() 3>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME, 3D_P1) )
FEELPP_CRB_PLUGIN_TEMPLATE( Eye2Brain2D_P2, Eye2Brain<2 BOOST_PP_COMMA() 2>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME, 2D_P2) )
FEELPP_CRB_PLUGIN_TEMPLATE( Eye2Brain3D_P2, Eye2Brain<2 BOOST_PP_COMMA() 3>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME, 3D_P2) )


} // namespace Feel
