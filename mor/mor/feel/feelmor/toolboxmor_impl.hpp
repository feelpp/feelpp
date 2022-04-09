
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
        // .add(deimOptions("vec"))
        // .add(deimOptions("mat"))
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

template<typename ToolboxType>
typename DeimMorModelToolbox<ToolboxType>::deim_function_type
DeimMorModelToolbox<ToolboxType>::deimFunction()
{
    deim_function_type assembleDEIM =
        [&tb=M_tb,&rhs=M_rhs,&mat=M_mat](parameter_type const& mu)
            {
                for( int i = 0; i < mu.size(); ++i )
                    tb->addParameterInModelProperties(mu.parameterName(i), mu(i));
                tb->updateParameterValues();
                rhs->zero();
                tb->algebraicFactory()->applyAssemblyLinear( tb->algebraicBlockVectorSolution()->vectorMonolithic(), mat, rhs, {"ignore-assembly.lhs"} );
                return rhs;
            };
    return assembleDEIM;
}

template<typename ToolboxType>
typename DeimMorModelToolbox<ToolboxType>::mdeim_function_type
DeimMorModelToolbox<ToolboxType>::mdeimFunction()
{
    mdeim_function_type assembleMDEIM =
        [&tb=M_tb,&rhs=M_rhs,&mat=M_mat](parameter_type const& mu)
            {
                for( int i = 0; i < mu.size(); ++i )
                    tb->addParameterInModelProperties(mu.parameterName(i), mu(i));
                tb->updateParameterValues();
                mat->zero();
                tb->algebraicFactory()->applyAssemblyLinear( tb->algebraicBlockVectorSolution()->vectorMonolithic(), mat, rhs, {"ignore-assembly.rhs"} );
                return mat;
            };
    return assembleMDEIM;
}

template<typename ToolboxType>
typename DeimMorModelToolbox<ToolboxType>::deim_function_type
DeimMorModelToolbox<ToolboxType>::deimOnlineFunction(mesh_ptrtype const& mesh)
{
    M_tbDeim = std::make_shared<toolbox_type>(M_prefix);
    M_tbDeim->setMesh(mesh);
    M_tbDeim->init();
    M_rhsDeim = M_tbDeim->algebraicFactory()->rhs()->clone();
    M_matDeim = M_tbDeim->algebraicFactory()->matrix();
    deim_function_type assembleDEIM =
        [&tbDeim=M_tbDeim,&rhs=M_rhsDeim,&mat=M_matDeim](parameter_type const& mu)
            {
                for( int i = 0; i < mu.size(); ++i )
                    tbDeim->addParameterInModelProperties(mu.parameterName(i), mu(i));
                tbDeim->updateParameterValues();
                rhs->zero();
                tbDeim->algebraicFactory()->applyAssemblyLinear( tbDeim->algebraicBlockVectorSolution()->vectorMonolithic(), mat, rhs, {"ignore-assembly.lhs"} );
                return rhs;
            };
    return assembleDEIM;
}

template<typename ToolboxType>
typename DeimMorModelToolbox<ToolboxType>::mdeim_function_type
DeimMorModelToolbox<ToolboxType>::mdeimOnlineFunction(mesh_ptrtype const& mesh)
{
    M_tbMdeim = std::make_shared<toolbox_type>(M_prefix);
    M_tbMdeim->setMesh(mesh);
    M_tbMdeim->init();
    M_rhsMdeim = M_tbMdeim->algebraicFactory()->rhs()->clone();
    M_matMdeim = M_tbMdeim->algebraicFactory()->matrix();
    mdeim_function_type assembleMDEIM =
        [&tbMdeim=M_tbMdeim,&rhs=M_rhsMdeim,&mat=M_matMdeim](parameter_type const& mu)
            {
                for( int i = 0; i < mu.size(); ++i )
                    tbMdeim->addParameterInModelProperties(mu.parameterName(i), mu(i));
                tbMdeim->updateParameterValues();
                mat->zero();
                tbMdeim->algebraicFactory()->applyAssemblyLinear( tbMdeim->algebraicBlockVectorSolution()->vectorMonolithic(), mat, rhs, {"ignore-assembly.rhs"} );
                return mat;
            };
    return assembleMDEIM;
}

template<typename SpaceType, int Options>
ToolboxMor<SpaceType, Options>::ToolboxMor(std::string const& name, std::string const& prefix)
    :
    super_type(name, Environment::worldCommPtr(), prefix)
{
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
    auto outputs = M_modelProperties->outputs().ofTypes({"integrate","mean","sensor"});
    return 1 + outputs.size();
}

template<typename SpaceType, int Options>
int ToolboxMor<SpaceType, Options>::Ql( int l)
{
    if( l == 0 )
        return 1;
    else if( l < Nl() )
        return 1;
    else
        return 0;
}
template<typename SpaceType, int Options>
int ToolboxMor<SpaceType, Options>::mLQF( int l, int q )
{
    if( l == 0 )
        return this->deim()->size();
    else if( l < Nl() )
    {
        auto outputs = M_modelProperties->outputs().ofTypes({"integrate","mean","sensor"});
        auto output = std::next(outputs.begin(), l-1)->second;
        if( M_outputDeim[output.name()] )
            return this->deim(M_outputDeim[output.name()])->size();
        else
            return 1;
    }
    else
        return 0;
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
            this->M_Aqm[q].resize(mQA(q), backend()->newMatrix(_test=this->Xh, _trial=this->Xh ) );
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
    if( tag == 0 )
        return M_assembleForDEIM(mu);
    else
    {
        auto outputs = M_modelProperties->outputs().ofTypes({"integrate","mean","sensor"});
        auto output = outputs[M_outputDeimName[tag-1]];
        if( output.type() == "mean" )
            return assembleOutputMean(mu, output);
        else if( output.type() == "integrate" )
            return assembleOutputIntegrate(mu, output);
        else //if( output.type() == "sensor" )
            return assembleOutputSensor(mu, output);
    }
}

template<typename SpaceType, int Options>
typename ToolboxMor<SpaceType, Options>::vector_ptrtype
ToolboxMor<SpaceType, Options>::assembleOutputMean( parameter_type const& mu, CRBModelOutput& output)
{
    auto f = form1(_test=this->Xh);
    auto u = this->Xh->element();
    if constexpr( is_scalar )
    {
        auto pm = mu.toParameterValues();
        output.setParameterValues(pm);
        auto ex = output.expression();
        auto exU = ex.diff<1>("crb_u");
        auto exGU0 = ex.diff<1>("crb_grad_u_0");
        auto exGU1 = ex.diff<1>("crb_grad_u_1");
        auto exGU2 = ex.diff<1>("crb_grad_u_2");
        auto exDNU = ex.diff<1>("crb_dn_u");
        auto dim = output.dim();
        if( dim == nDim )
        {
            auto range = markedelements(this->Xh->mesh(), output.markers());
            double area = integrate( _range=range, _expr=cst(1.0) ).evaluate()(0,0) ;
            f += integrate( _range=range, _expr=(exU*id(u)+exGU0*dx(u)+exGU1*dy(u)+exGU2*dz(u))/cst(area) );
        }
        else if( dim == nDim-1 )
        {
            auto range = markedfaces(this->Xh->mesh(), output.markers());
            double area = integrate( _range=range, _expr=cst(1.0) ).evaluate()(0,0) ;
            f += integrate( _range=range, _expr=(exU*id(u)+exGU0*dx(u)+exGU1*dy(u)+exGU2*dz(u)+exDNU*dn(u))/cst(area) );
        }
    }
   return f.vectorPtr();
}

template<typename SpaceType, int Options>
typename ToolboxMor<SpaceType, Options>::vector_ptrtype
ToolboxMor<SpaceType, Options>::assembleOutputIntegrate( parameter_type const& mu, CRBModelOutput& output)
{
    auto f = form1(_test=this->Xh);
    auto u = this->Xh->element();
    if constexpr( is_scalar )
    {
        auto pm = mu.toParameterValues();
        output.setParameterValues(pm);
        auto ex = output.expression();
        auto exU = ex.diff<1>("crb_u");
        auto exGU0 = ex.diff<1>("crb_grad_u_0");
        auto exGU1 = ex.diff<1>("crb_grad_u_1");
        auto exGU2 = ex.diff<1>("crb_grad_u_2");
        auto exDNU = ex.diff<1>("crb_dn_u");
        auto dim = output.dim();
        if( dim == nDim )
        {
            auto range = markedelements(this->Xh->mesh(), output.markers());
            f += integrate( _range=range, _expr=exU*id(u)+exGU0*dx(u)+exGU1*dy(u)+exGU2*dz(u) );
        }
        else if( dim == nDim-1 )
        {
            auto range = markedfaces(this->Xh->mesh(), output.markers());
            f += integrate( _range=range, _expr=exU*id(u)+exGU0*dx(u)+exGU1*dy(u)+exGU2*dz(u)+exDNU*dn(u) );
        }
    }
    return f.vectorPtr();    
}

template<typename SpaceType, int Options>
typename ToolboxMor<SpaceType, Options>::vector_ptrtype
ToolboxMor<SpaceType, Options>::assembleOutputSensor( parameter_type const& mu, CRBModelOutput& output)
{
    auto f = form1(_test=this->Xh);
    auto u = this->Xh->element();
    if constexpr( is_scalar )
    {
        auto coord = output.coord();
        auto radius = output.radius();
        if( coord.size() >= nDim )
        {
            if constexpr( nDim == 1 )
            {
                auto phi = exp( -inner(P()-vec(cst(coord[0])))/(2*std::pow(radius,2)));
                auto n = integrate(_range=elements(support(this->Xh)), _expr=phi).evaluate()(0,0);
                f += integrate( _range=elements(support(this->Xh)), _expr=id(u)*phi/n);
            }
            else if constexpr( nDim == 2 )
            {
                auto phi = exp( -inner(P()-vec(cst(coord[0]),cst(coord[1])))/(2*std::pow(radius,2)));
                auto n = integrate(_range=elements(support(this->Xh)), _expr=phi).evaluate()(0,0);
                f += integrate( _range=elements(support(this->Xh)), _expr=id(u)*phi/n);
            }
            else if constexpr( nDim == 3 )
            {
                auto phi = exp( -inner(P()-vec(cst(coord[0]),cst(coord[1]),cst(coord[2])))/(2*std::pow(radius,2)));
                auto n = integrate(_range=elements(support(this->Xh)), _expr=phi).evaluate()(0,0);
                f += integrate( _range=elements(support(this->Xh)), _expr=id(u)*phi/n);
            }
        }
    }
    return f.vectorPtr();    
}

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    Feel::cout << "setupSpecificityModel" << std::endl;

    if( this->hasModelFile("property-file") )
        M_propertyPath = this->additionalModelFiles().find("property-file")->second;
    else
        Feel::cerr << "Warning!! the database does not contain the property file! Expect bugs!"
                   << std::endl;

    M_modelProperties = std::make_shared<CRBModelProperties>(); // TODO : put directoryLibExpr and worldcomm
    M_modelProperties->setup( M_propertyPath );
    auto parameters = M_modelProperties->parameters();
    this->Dmu = parameterspace_type::New(parameters);
    auto parameterNames = std::set<std::string>(this->Dmu->parameterNames().begin(), this->Dmu->parameterNames().end());

    M_deim = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _prefix="vec");
    this->addDeim(M_deim);
    Feel::cout << tc::green << "ToolboxMor DEIM construction finished!!" << tc::reset << std::endl;

    M_mdeim = Feel::mdeim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _prefix="mat");
    this->addMdeim(M_mdeim);
    Feel::cout << tc::green << "ToolboxMor MDEIM construction finished!!" << tc::reset << std::endl;

    // outputs
    int i = 1;
    auto outputs = this->M_modelProperties->outputs().ofTypes({"integrate","mean","sensor"});
    for( auto const& [name,output] : outputs )
    {
        M_outputName.push_back(name);
        if( (output.type() == "mean" || output.type() == "integrate") && output.expression().hasSymbolDependency(parameterNames) )
        {
            Feel::cout << "output " << name << " depends on parameters" << std::endl;
            M_outputDeim[name] = i;
            M_outputDeimName.push_back(name);
            auto dO = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _prefix="output_"+name, _tag=i);
            this->addDeim(dO);
            Feel::cout << tc::green << "ToolboxMor DEIM for output " << name << " construction finished!!" << tc::reset << std::endl;
            i++;
        }
        else 
            M_outputDeim[name] = 0;
    }

    this->resizeQm(false);
}

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::initModel()
{
    M_propertyPath = Environment::expand( soption("toolboxmor.filename"));
    M_trainsetDeimSize = ioption("toolboxmor.trainset-deim-size");
    M_trainsetMdeimSize = ioption("toolboxmor.trainset-mdeim-size");
    this->addModelFile("property-file", M_propertyPath);

    M_modelProperties = std::make_shared<CRBModelProperties>(); // TODO : put directoryLibExpr and worldcomm
    M_modelProperties->setup( M_propertyPath );
    auto parameters = M_modelProperties->parameters();
    this->Dmu = parameterspace_type::New(parameters);
    auto parameterNames = std::set<std::string>(this->Dmu->parameterNames().begin(), this->Dmu->parameterNames().end());

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << this->Xh->nLocalDof() << "\n";
        std::cout << "Number of dof " << this->Xh->nDof() << "\n";
    }

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
    Feel::cout << tc::green << "ToolboxMor DEIM construction finished!!" << tc::reset << std::endl;

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
    Feel::cout << tc::green << "ToolboxMor MDEIM construction finished!!" << tc::reset << std::endl;

    // outputs
    int i = 1;
    auto outputs = this->M_modelProperties->outputs().ofTypes({"integrate","mean","sensor"});
    for( auto const& [name,output] : outputs )
    {
        M_outputName.push_back(name);
        if( (output.type() == "mean" || output.type() == "integrate") && output.expression().hasSymbolDependency(parameterNames) )
        {
            Feel::cout << "output " << name << " depends on parameters" << std::endl;
            M_outputDeim[name] = i;
            M_outputDeimName.push_back(name);
            auto dO = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _sampling=PsetV, _prefix="output_"+name, _tag=i);
            this->addDeim(dO);
            this->deim(i)->run();
            Feel::cout << tc::green << "ToolboxMor DEIM for output " << name << " construction finished!!" << tc::reset << std::endl;
            i++;
        }
        else 
            M_outputDeim[name] = 0;
    }

} // ToolboxMor<SpaceType, Options>::initModel

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::initOnlineModel( std::shared_ptr<super_type> const& model )
{
    auto tbModel = std::dynamic_pointer_cast<self_type>(model);
    this->M_modelProperties = tbModel->modelProperties();
    this->M_outputDeimName = tbModel->outputDeimName();
}

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::initToolbox(std::shared_ptr<DeimMorModelBase<mesh_type>> model )
{
    this->setAssembleDEIM(model->deimFunction());
    this->setAssembleMDEIM(model->mdeimFunction());

    this->initModel();

    this->setOnlineAssembleDEIM(model->deimOnlineFunction(this->getDEIMReducedMesh()));
    this->setOnlineAssembleMDEIM(model->mdeimOnlineFunction(this->getMDEIMReducedMesh()));

    this->postInitModel();
    this->setInitialized(true);
}

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::initOnlineToolbox(std::shared_ptr<DeimMorModelBase<mesh_type>> model )
{
    this->setOnlineAssembleDEIM(model->deimOnlineFunction(this->getDEIMReducedMesh()));
    this->setOnlineAssembleMDEIM(model->mdeimOnlineFunction(this->getMDEIMReducedMesh()));
}

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

    // for now, only with M independant on mu
    int M_M = 1;
    this->M_betaMqm[0].resize( M_M );
    for ( int i = 0; i < M_M; ++i )
        this->M_betaMqm[0][i] = 1;

    int output = 1;
    auto outputs = M_modelProperties->outputs().ofTypes({"integrate","mean","sensor"});
    for( auto const& [name,out] : outputs )
    {
        if( M_outputDeim[name] )
        {
            auto betaO = this->deim(M_outputDeim[name])->beta(mu);
            int M_O = this->deim(M_outputDeim[name])->size();
            for( int i = 0; i < M_O; ++i )
                this->M_betaFqm[output][0][i] = betaO(i);
        }
        else
            this->M_betaFqm[output][0][0] = 1;
        output++;
    }
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

    auto u = this->Xh->element();
    auto mu = this->Dmu->min();

    // output
    auto outputs = M_modelProperties->outputs().ofTypes({"integrate","mean","sensor"});
    int output = 1;
    for( auto const& outp : outputs )
    {
        auto out = outp.second;
        if( M_outputDeim[out.name()] )
        {
            int M_O = this->deim(M_outputDeim[out.name()])->size();
            auto qo = this->deim(M_outputDeim[out.name()])->q();
            for( int i = 0; i < M_O; ++i )
                this->M_Fqm[output][0][i] = qo[i];
        }
        else if( out.type() == "mean" )
        {
            this->M_Fqm[output][0][0] = assembleOutputMean(mu, out);
        }
        else if( out.type() == "integrate")
        {
            this->M_Fqm[output][0][0] = assembleOutputIntegrate(mu, out);
        }
        else if( out.type() == "sensor")
        {
            this->M_Fqm[output][0][0] = assembleOutputSensor(mu, out);
        }
        output++;
    }

    // Energy matrix
    auto m = this->assembleForMDEIM(mu,0);
    this->M_energy_matrix = backend()->newMatrix(_test=this->Xh, _trial=this->Xh );
    m->symmetricPart(this->M_energy_matrix);

    // for now, only with M independant of mu
    this->M_Mqm.resize( 1 );
    this->M_Mqm[0].resize( 1 );
    this->M_Mqm[0][0] = backend()->newMatrix(_test=this->Xh, _trial=this->Xh);
    form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Mqm[0][0])
        = integrate(_range=elements(support(this->Xh)), _expr=inner(id(u),idt(u)));

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
            {
                auto d = dot( *this->M_Fqm[output_index][q][m] , u );
                // Feel::cout << "d["<< output_index << "][" << q << "][" << m << "] = " << d << "\n"
                //            << "beta = " << this->M_betaFqm[output_index][q][m] << std::endl;
                s += this->M_betaFqm[output_index][q][m]*d;
            }
        }
    }
    else
    {
        throw std::logic_error( "[ToolboxMor::output] error with output_index : only 0 or 1 " );
    }

    return s ;
}

}
