
namespace Feel
{

po::options_description
makeToolboxMorOptions()
{
    po::options_description options( "ToolboxMor options" );
    options.add_options()
        ( "toolboxmor.filename", po::value<std::string>()->default_value("toolbox.json"), "json file defining the parameters" )
        ( "toolboxmor.trainset-deim-size", po::value<int>()->default_value(40), "size of the deim trainset" )
        ( "toolboxmor.trainset-mdeim-size", po::value<int>()->default_value(40), "size of the mdeim trainset" )
        ( "toolboxmor.test-deim", po::value<bool>()->default_value(false), "test deim decomposition" )
        ;
    options.add(crbOptions())
        .add(crbSEROptions())
        .add(eimOptions())
        .add(podOptions())
        .add(deimOptions("vec"))
        .add(deimOptions("mat"))
        .add(backend_options("backend-primal"))
        .add(backend_options("backend-dual"))
        .add(backend_options("backend-l2"))
        .add(bdf_options("OpusHeatTB"));
    return options;
}
AboutData
makeToolboxMorAbout( std::string const& str )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "toolboxmor",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2009-2014 Feel++ Consortium");
    return about;
}

template<typename SpaceType, int Options>
ToolboxMor<SpaceType, Options>::ToolboxMor(std::string const& prefix)
    :
    super_type("ToolboxMor", Environment::worldCommPtr(), prefix),
    M_propertyPath(Environment::expand( soption("toolboxmor.filename"))),
    M_trainsetDeimSize(ioption("toolboxmor.trainset-deim-size")),
    M_trainsetMdeimSize(ioption("toolboxmor.trainset-mdeim-size"))
{
    M_modelProperties = std::make_shared<ModelProperties>(M_propertyPath);

    auto parameters = M_modelProperties->parameters();
    int nbCrbParameters = count_if(parameters.begin(), parameters.end(), [] (auto const& p)
                                   {
                                       return p.second.hasMinMax();
                                   });
    this->Dmu->setDimension(nbCrbParameters);
    auto mu_min = this->Dmu->element();
    auto mu_max = this->Dmu->element();
    int i = 0;
    for( auto const& parameterPair : parameters )
    {
        if( parameterPair.second.hasMinMax() )
        {
            mu_min(i) = parameterPair.second.min();
            mu_max(i) = parameterPair.second.max();
            this->Dmu->setParameterName(i++, parameterPair.first );
        }
    }
    this->Dmu->setMin(mu_min);
    this->Dmu->setMax(mu_max);
    // M_mu = Dmu->element();

}

template<typename SpaceType, int Options>
int ToolboxMor<SpaceType, Options>::Qa()
{
    return 1;
}

template<typename SpaceType, int Options>
int ToolboxMor<SpaceType, Options>::mQA( int q )
{
    return this->mdeim()->size();
}

template<typename SpaceType, int Options>
int ToolboxMor<SpaceType, Options>::Nl()
{
    return 1;
}

template<typename SpaceType, int Options>
int ToolboxMor<SpaceType, Options>::Ql( int l)
{
    return 1;
}
template<typename SpaceType, int Options>
int ToolboxMor<SpaceType, Options>::mLQF( int l, int q )
{
    switch( l )
    {
    case 0:
        return this->deim()->size();
    // case 1:
    //     return 1;
    default:
        return 0;
    }
}
template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::resizeQm( bool resizeMat )
{
    if( resizeMat )
        this->M_Aqm.resize( Qa());
    this->M_betaAqm.resize( Qa() );
    for( int q = 0; q < Qa(); ++q )
    {
        if( resizeMat )
            this->M_Aqm[q].resize(mQA(q), backend()->newMatrix(this->Xh, this->Xh ) );
        this->M_betaAqm[q].resize(mQA(q));
    }

    if( resizeMat )
        this->M_Fqm.resize(Nl());
    this->M_betaFqm.resize(Nl());
    for( int l = 0; l < Nl(); ++l )
    {
        if( resizeMat )
            this->M_Fqm[l].resize(Ql(l));
        this->M_betaFqm[l].resize(Ql(l));
        for( int q = 0; q < Ql(l); ++q )
        {
            if( resizeMat )
                this->M_Fqm[l][q].resize(mLQF(l, q), backend()->newVector(this->Xh) );
            this->M_betaFqm[l][q].resize(mLQF(l, q) );
        }
    }

    // if( resizeMat )
    // {
    //     M_InitialGuess.resize(1);
    //     M_InitialGuess[0].resize(1);
    //     M_InitialGuess[0][0] = Xh->elementPtr();
    // }
}

template<typename SpaceType, int Options>
typename ToolboxMor<SpaceType, Options>::sparse_matrix_ptrtype
ToolboxMor<SpaceType, Options>::assembleForMDEIM( parameter_type const& mu, int const& tag )
{
    return M_assembleForMDEIM(mu);
}

template<typename SpaceType, int Options>
typename ToolboxMor<SpaceType, Options>::vector_ptrtype
ToolboxMor<SpaceType, Options>::assembleForDEIM( parameter_type const& mu, int const& tag )
{
    return M_assembleForDEIM(mu);
}

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    auto ptreeSpecificityOfModel = ptree.get_child_optional( "specifity-of-model" );
    CHECK( ptreeSpecificityOfModel ) << "invalid ptree : section specifity-of-model is missing";
    // M_measureMarkedSurface["IC2"] = ptreeSpecificityOfModel->template get<double>( "surface-measure-IC2" );
    //std::cout << "surface loaded " << M_surface << "\n";
}
template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::updateSpecificityModel( boost::property_tree::ptree & ptree ) const
{
    // ptree.add( "surface-measure-IC2", M_measureMarkedSurface.find("IC2")->second );
}


template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::initModel()
{
    this->addModelFile("property-file", M_propertyPath);

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << this->Xh->nLocalDof() << "\n";
        std::cout << "Number of dof " << this->Xh->nDof() << "\n";
    }

    // for ( std::string const& marker : std::vector<std::string>({"AIR","PCB","IC1","IC2"}) )
    //     M_measureMarkedSurface[marker] = measure(_range=markedelements(M_heatBox->mesh(), marker) );

    auto PsetV = this->Dmu->sampling();
    std::string supersamplingname =(boost::format("DmuDEim-P%1%-Ne%2%-generated-by-master-proc") % this->Dmu->dimension() % M_trainsetDeimSize ).str();
    std::ifstream file ( supersamplingname );
    bool all_proc_same_sampling=true;
    if( ! file )
    {
        PsetV->randomize( M_trainsetDeimSize , all_proc_same_sampling , supersamplingname );
        PsetV->writeOnFile( supersamplingname );
    }
    else
    {
        PsetV->clear();
        PsetV->readFromFile(supersamplingname);
    }
    M_deim = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _sampling=PsetV, _prefix="vec");
    this->addDeim(M_deim);
    this->deim()->run();
    Feel::cout << tc::green << "Electric DEIM construction finished!!" << tc::reset << std::endl;

    auto PsetM = this->Dmu->sampling();
    supersamplingname =(boost::format("DmuMDEim-P%1%-Ne%2%-generated-by-master-proc") % this->Dmu->dimension() % M_trainsetMdeimSize ).str();
    std::ifstream fileM ( supersamplingname );
    if( ! fileM )
    {
        PsetM->randomize( M_trainsetMdeimSize , all_proc_same_sampling , supersamplingname );
        PsetM->writeOnFile( supersamplingname );
    }
    else
    {
        PsetM->clear();
        PsetM->readFromFile(supersamplingname);
    }
    M_mdeim = Feel::mdeim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _sampling=PsetM, _prefix="mat");
    this->addMdeim(M_mdeim);
    this->mdeim()->run();
    Feel::cout << tc::green << "Electric MDEIM construction finished!!" << tc::reset << std::endl;

} // ToolboxMor<SpaceType, Options>::init

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::postInitModel()
{
    this->resizeQm();
    assembleData();
}

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::updateBetaQ_impl( parameter_type const& mu , double time , bool only_terms_time_dependent )
{
    auto betaA = this->mdeim()->beta(mu);
    auto betaF = this->deim()->beta(mu);

    int M_A = this->mdeim()->size();
    for( int i = 0; i < M_A; ++i )
        this->M_betaAqm[0][i] = betaA(i);

    int M_F = this->deim()->size();
    for( int i = 0; i < M_F; ++i )
        this->M_betaFqm[0][0][i] = betaF(i);


    this->M_betaMqm.resize( 1 );
    Feel::cout << this->M_betaMqm.size() << std::endl;

    // for now, only with M independant on mu
    int M_M = 1;
    this->M_betaMqm[0].resize( M_M );
    for ( int i = 0; i < M_M; ++i )
        this->M_betaMqm[0][i] = 1;

    // this->M_betaFqm[1][0][0] = 1./M_measureMarkedSurface["IC2"];
}



template<typename SpaceType, int Options>
typename ToolboxMor<SpaceType, Options>::super_type::betaqm_type
ToolboxMor<SpaceType, Options>::computeBetaQm( parameter_type const& mu , double time , bool only_terms_time_dependent )
{
    return computeBetaQ_impl<super_type>( mu, time, only_terms_time_dependent );
}
template<typename SpaceType, int Options>
typename ToolboxMor<SpaceType, Options>::super_type::betaqm_type
ToolboxMor<SpaceType, Options>::computeBetaQm( parameter_type const& mu )
{
    return computeBetaQ_impl<super_type>( mu );
}
template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::assembleData()
{
    int M_A = this->mdeim()->size();
    auto qa = this->mdeim()->q();
    for( int i = 0; i < M_A; ++i )
        this->M_Aqm[0][i] = qa[i];

    int M_F = this->deim()->size();
    auto qf = this->deim()->q();
    for( int i = 0; i < M_F; ++i )
        this->M_Fqm[0][0][i] = qf[i];

    // output
    // form1(_test=this->Xh, _vector=this->M_Fqm[1][0][0])
    //     += integrate( _range=markedelements( this->Xh->mesh(), "IC2" ), _expr= id( v ) );

    // Energy matrix
    auto mu = this->Dmu->element();
    for( int i = 0; i < this->Dmu->dimension(); ++i )
        mu(i) = 1;
    auto m = this->assembleForMDEIM(mu,0);
    // for( int i = 0; i < Dmu->dimension(); ++i )
    //     M_heatBox->addParameterInModelProperties(Dmu->parameterName(i), 1);
    // M_heatBox->updateParameterValues();
    // M_heatBox->updateFieldVelocityConvection();
    // M_heatBox->algebraicFactory()->assembleLinear(M_heatBox->blockVectorSolution().vectorMonolithic());
    // auto m = M_heatBox->algebraicFactory()->matrix();
    this->M_energy_matrix = backend()->newMatrix(this->Xh, this->Xh );
    m->symmetricPart(this->M_energy_matrix);

    // for now, only with M independant of mu
    auto u = this->Xh->element();

    this->M_Mqm.resize( 1 );
    this->M_Mqm[0].resize( 1 );
    this->M_Mqm[0][0] = backend()->newMatrix(this->Xh, this->Xh);
    form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Mqm[0][0])
        = integrate(_range=elements(this->Xh->mesh()), _expr=inner(id(u),idt(u)));

}

// ToolboxMor<SpaceType, Options>::element_type
// ToolboxMor<SpaceType, Options>::solve( parameter_type const& mu )
// {
//     for( int i = 0; i < mu.size(); ++i )
//         M_heatBox->addParameterInModelProperties(mu.parameterName(i), mu(i));
//     M_heatBox->updateParameterValues();
//     M_heatBox->updateFieldVelocityConvection();
//     M_heatBox->solve();
//     return M_heatBox->fieldTemperature();
// }


template<typename SpaceType, int Options>
double
ToolboxMor<SpaceType, Options>::output( int output_index, parameter_type const& mu, element_type &u, bool need_to_solve )
{
    using namespace vf;

    CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    double s=0;

    //bool dont_use_operators_free = boption(_name="do-not-use-operators-free") ;
    auto fqm = backend()->newVector( this->Xh );
    if ( output_index < Nl() )
    {
        for ( int q=0; q<Ql(output_index); q++ )
        {
            for( int m=0; m<mLQF(output_index, q); m++ )
                s += this->M_betaFqm[output_index][q][m]*dot( *this->M_Fqm[output_index][q][m] , u );
        }
    }
    else
    {
        throw std::logic_error( "[ToolboxMor::output] error with output_index : only 0 or 1 " );
    }

    return s ;
}

}
