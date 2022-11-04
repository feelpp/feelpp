#include <feel/feeldiscr/sensors.hpp>
#include <fmt/ostream.h>
#include <feel/feelcore/enumerate.hpp>
namespace Feel
{

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
    if ( M_toolboxInitFunction )
    {
        M_tbDeim = M_toolboxInitFunction( mesh );
    }
    else
    {
        M_tbDeim = toolbox_type::New(_prefix=M_prefix);
        M_tbDeim->setMesh(mesh);
        M_tbDeim->init();
        //M_tbDeim->printAndSaveInfo();
    }
    M_rhsDeim = M_tbDeim->algebraicFactory()->rhs()->clone();
    M_matDeim = M_tbDeim->algebraicFactory()->matrix();
    deim_function_type assembleDEIM =
        [this](parameter_type const& mu)
            {
                for( int i = 0; i < mu.size(); ++i )
                    M_tbDeim->addParameterInModelProperties(mu.parameterName(i), mu(i));
                M_tbDeim->updateParameterValues();
                M_rhsDeim->zero();
                M_tbDeim->algebraicFactory()->applyAssemblyLinear( M_tbDeim->algebraicBlockVectorSolution()->vectorMonolithic(),
                                                                   M_matDeim, M_rhsDeim, {"ignore-assembly.lhs"} );
                return M_rhsDeim;
            };
    return assembleDEIM;
}

template<typename ToolboxType>
typename DeimMorModelToolbox<ToolboxType>::mdeim_function_type
DeimMorModelToolbox<ToolboxType>::mdeimOnlineFunction(mesh_ptrtype const& mesh)
{
    if ( M_toolboxInitFunction )
    {
        M_tbMdeim = M_toolboxInitFunction( mesh );
    }
    else
    {
        M_tbMdeim = toolbox_type::New(_prefix=M_prefix);
        M_tbMdeim->setMesh(mesh);
        M_tbMdeim->init();
        //M_tbMdeim->printAndSaveInfo();
    }
    M_rhsMdeim = M_tbMdeim->algebraicFactory()->rhs()->clone();
    M_matMdeim = M_tbMdeim->algebraicFactory()->matrix();
    mdeim_function_type assembleMDEIM =
        [this](parameter_type const& mu)
            {
                for( int i = 0; i < mu.size(); ++i )
                    M_tbMdeim->addParameterInModelProperties(mu.parameterName(i), mu(i));
                M_tbMdeim->updateParameterValues();
                M_matMdeim->zero();
                M_tbMdeim->algebraicFactory()->applyAssemblyLinear( M_tbMdeim->algebraicBlockVectorSolution()->vectorMonolithic(),
                                                                    M_matMdeim, M_rhsMdeim, {"ignore-assembly.rhs"} );
                return M_matMdeim;
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
    auto outputs = M_modelProperties->outputs().ofTypes({"integrate","mean","sensor","point"});
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
        auto outputs = M_modelProperties->outputs().ofTypes({"integrate","mean","sensor","point"});
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
        auto outputs = M_modelProperties->outputs().ofTypes({"integrate","mean","sensor","point"});
        auto output = outputs[M_outputDeimName[tag-1]];
        if( output.type() == "mean" )
            return assembleOutputMean(mu, output);
        else if( output.type() == "integrate" )
            return assembleOutputIntegrate(mu, output);
        else if( output.type() == "sensor" )
            return assembleOutputSensor(mu, output);
        else //if( output.type() == "point" )
            return assembleOutputPoint(mu, output);
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
        LOG(INFO) << fmt::format("[ToolboxMor::assembleOutputMean] ex: {}", ex);
        LOG(INFO) << fmt::format("[ToolboxMor::assembleOutputMean] exU: {}", exU);
        LOG(INFO) << fmt::format("[ToolboxMor::assembleOutputMean] exGU0: {}", exGU0);
        LOG(INFO) << fmt::format("[ToolboxMor::assembleOutputMean] exGU1: {}", exGU1);
        LOG(INFO) << fmt::format("[ToolboxMor::assembleOutputMean] exGU2: {}", exGU2);
        LOG(INFO) << fmt::format("[ToolboxMor::assembleOutputMean] exDNU: {}", exDNU);
        auto dim = output.dim();
        LOG(INFO) << "[ToolboxMor::assembleOutputMean] assemble dim = " << dim << std::endl;
        if( dim == nDim )
        {
            auto range = markedelements(this->Xh->mesh(), output.markers());
            double measure = integrate( _range=range, _expr=cst(1.0) ).evaluate()(0,0) ;
            LOG(INFO) << "[ToolboxMor::assembleOutputMean] measure = " << measure << std::endl;
            f += integrate( _range=range, _expr=(exU*id(u)+exGU0*dx(u)+exGU1*dy(u)+exGU2*dz(u))/cst(measure) );
            u.on(_range=range, _expr=cst(1.));
            LOG(INFO) << "[ToolboxMor::assembleOutputMean] f(1) = " << f(u) << std::endl;
        }
        else if( dim == nDim-1 )
        {
            auto range = markedfaces(this->Xh->mesh(), output.markers());
            double measure = integrate( _range=range, _expr=cst(1.0) ).evaluate()(0,0) ;
            LOG(INFO) << "[ToolboxMor::assembleOutputMean] measure = " << measure << std::endl;
            f += integrate( _range=range, _expr=(exU*id(u)+exGU0*dx(u)+exGU1*dy(u)+exGU2*dz(u)+exDNU*dn(u))/cst(measure) );
            u.on(_range=range, _expr=cst(1.));
            LOG(INFO) << "[ToolboxMor::assembleOutputMean] f(1) = " << f(u) << std::endl;
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
    if constexpr( is_scalar )
    {
        auto coord = output.coord();
        node_type n(space_type::nDim);
        for( int i = 0; i < space_type::nDim; ++i )
            n(i) = coord[i];
        auto radius = output.radius();
        auto s = std::make_shared<SensorGaussian<space_type>>(this->Xh, n, radius, output.name());
        return s->containerPtr();
    } else {
        auto f = form1(_test=this->Xh);
        return f.vectorPtr();
    }
}

template<typename SpaceType, int Options>
typename ToolboxMor<SpaceType, Options>::vector_ptrtype
ToolboxMor<SpaceType, Options>::assembleOutputPoint( parameter_type const& mu, CRBModelOutput& output)
{
    if constexpr( is_scalar )
    {
        auto coord = output.coord();
        node_type n(space_type::nDim);
        for( int i = 0; i < space_type::nDim; ++i )
            n(i) = coord[i];
        auto s = std::make_shared<SensorPointwise<space_type>>(this->Xh, n, output.name());
        return s->containerPtr();
    } else {
        auto f = form1(_test=this->Xh);
        return f.vectorPtr();
    }
}

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    Feel::cout << "ToolboxMor::setupSpecificityModel" << std::endl;
    LOG(INFO) << "setupSpecificityModel" << std::endl;

    M_modelProperties = std::make_shared<CRBModelProperties>( Environment::exprRepository(), this->worldCommPtr()/*subWorldCommSeqPtr()*/ );
    if ( this->hasModelData( "crb_properties" ) )
    {
        auto & mdata = this->additionalModelData("crb_properties");
        auto const& jsonData = mdata.template fetch_data<nl::json>( dbDir/*this->crbModelDb().dbRepository()*/ );
        M_modelProperties->setup( jsonData );
    }
    else
        throw std::runtime_error( "the database does not contain the crb_properties file" );


    auto parameters = M_modelProperties->parameters();
    this->Dmu = parameterspace_type::New(parameters);
    auto parameterNames = std::set<std::string>(this->Dmu->parameterNames().begin(), this->Dmu->parameterNames().end());

    //std::cout << tc::green << "ToolboxMor DEIM Parameter names : " << parameterNames << tc::reset << std::endl;

    M_deim = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _prefix="vec");
    this->addDeim(M_deim);
    Feel::cout << tc::green << "ToolboxMor DEIM construction finished!!" << tc::reset << std::endl;
    LOG(INFO) << "ToolboxMor DEIM construction finished!!"  << std::endl;

    M_mdeim = Feel::mdeim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _prefix="mat");
    this->addMdeim(M_mdeim);
    Feel::cout << tc::green << "ToolboxMor MDEIM construction finished!!" << tc::reset << std::endl;
    LOG(INFO) << "ToolboxMor MDEIM construction finished!!" << std::endl;

    // outputs
    int i = 1;
    auto outputs = this->M_modelProperties->outputs().ofTypes({"integrate","mean","sensor"});
    for( auto const& [name,output] : outputs )
    {
        M_outputName.push_back(name);
        if( (output.type() == "mean" || output.type() == "integrate") && output.expression().hasSymbolDependency(parameterNames) )
        {
            LOG(INFO) << "output " << name << " depends on parameters" << std::endl;
            M_outputDeim[name] = i;
            M_outputDeimName.push_back(name);
            auto dO = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _prefix="output_"+name, _tag=i);
            this->addDeim(dO);
            Feel::cout << tc::green << "ToolboxMor DEIM for output " << name << " construction finished!!" << tc::reset << std::endl;
            LOG(INFO) << "ToolboxMor DEIM for output " << name << " construction finished!!" << std::endl;
            i++;
        }
        else
            M_outputDeim[name] = 0;
    }

    this->resizeQm(false);
}

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::initModelImpl()
{
    std::string propertyPath = Environment::expand( soption("toolboxmor.filename"));
    M_trainsetDeimSize = ioption("toolboxmor.trainset-deim-size");
    M_trainsetMdeimSize = ioption("toolboxmor.trainset-mdeim-size");

    M_modelProperties = std::make_shared<CRBModelProperties>( Environment::exprRepository(), this->worldCommPtr() );
    M_modelProperties->setup( propertyPath );
    this->addModelData( "crb_properties", M_modelProperties->jsonData(), "toolboxmor.crb_properties.json" );


    auto parameters = M_modelProperties->parameters();
    this->Dmu = parameterspace_type::New(parameters);
    auto parameterNames = std::set<std::string>(this->Dmu->parameterNames().begin(), this->Dmu->parameterNames().end());

    Feel::cout << "Number of local dof " << this->Xh->nLocalDof() << std::endl;
    Feel::cout << "Number of dof " << this->Xh->nDof() << std::endl;
    LOG(INFO) << "[ToolboxMor] Number of local dof " << this->Xh->nLocalDof() << "\n";
    LOG(INFO) << "[ToolboxMor] Number of dof " << this->Xh->nDof() << "\n";

    auto PsetV = this->Dmu->sampling();
    std::string supersamplingname = fmt::format("DmuDEim-P{}-Ne{}-generated-by-master-proc", this->Dmu->dimension(), M_trainsetDeimSize );
    std::ifstream file ( supersamplingname );
    Feel::cout << tc::blue << "[ToolboxMor] DEIM sampling file \"" << supersamplingname;
    LOG(INFO) << "[ToolboxMor] DEIM sampling file \"" << supersamplingname;
    bool all_proc_same_sampling=true;
    if( ! file )
    {
        this->worldComm().barrier();
        PsetV->randomize( M_trainsetDeimSize , all_proc_same_sampling , supersamplingname );
        PsetV->writeOnFile( supersamplingname );
    }
    else
    {
        Feel::cout << "\" found" << tc::reset << std::endl;
        LOG(INFO) << "\" found.\n";
        PsetV->clear();
        PsetV->readFromFile(supersamplingname);
    }
    M_deim = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _sampling=PsetV, _prefix="vec");
    this->addDeim(M_deim);
    this->deim()->run();
    Feel::cout << tc::green << "[ToolboxMor] DEIM construction finished!!" << tc::reset << std::endl;
    LOG(INFO) << "[ToolboxMor] DEIM construction finished!!" << std::endl;

    auto PsetM = this->Dmu->sampling();
    supersamplingname = fmt::format("DmuMDEim-P{}-Ne{}-generated-by-master-proc", this->Dmu->dimension(), M_trainsetMdeimSize );
    std::ifstream fileM ( supersamplingname );
    Feel::cout << tc::blue << "[ToolboxMor] MDEIM sampling file \"" << supersamplingname;
    LOG(INFO) << "[ToolboxMor] MDEIM sampling file \"" << supersamplingname;
    if( ! fileM )
    {
        this->worldComm().barrier();
        PsetM->randomize( M_trainsetMdeimSize , all_proc_same_sampling , supersamplingname );
        PsetM->writeOnFile( supersamplingname );
    }
    else
    {
        Feel::cout << "\" found" << tc::reset << std::endl;
        LOG(INFO) << "\" found.\n";
        PsetM->clear();
        PsetM->readFromFile(supersamplingname);
    }
    M_mdeim = Feel::mdeim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _sampling=PsetM, _prefix="mat");
    this->addMdeim(M_mdeim);
    this->mdeim()->run();
    LOG(INFO) << "[ToolboxMor] MDEIM construction finished!!" << std::endl;

    // outputs
    int i = 1;
    auto outputs = this->M_modelProperties->outputs().ofTypes({"integrate","mean","sensor"});
    for( auto const& [name,output] : outputs )
    {
        LOG(INFO) << "[ToolboxMor] output " << name << " type " << output.type() << " expression " << output.expression() << "\n";
        M_outputName.push_back(name);
        if( (output.type() == "mean" || output.type() == "integrate") && output.expression().hasSymbolDependency(parameterNames) )
        {
            LOG(INFO) << "[ToolboxMor] output " << name << " depends on parameters" << std::endl;
            M_outputDeim[name] = i;
            M_outputDeimName.push_back(name);
            auto dO = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _sampling=PsetV, _prefix="output_"+name, _tag=i);
            this->addDeim(dO);
            this->deim(i)->run();
            Feel::cout << tc::green << "[ToolboxMor] DEIM for output " << name << " construction finished!!" << tc::reset << std::endl;
            LOG(INFO) << "[ToolboxMor] DEIM for output " << name << " construction finished!!" << std::endl;
            i++;
        }
        else
        {
            Feel::cout << tc::green << "[ToolboxMor] output " << name << " does not depend on parameters" << tc::reset << std::endl;
            LOG(INFO) << "[ToolboxMor] output " << name << " does not depend on parameters" << std::endl;
            M_outputDeim[name] = 0;
        }
    }

} // ToolboxMor<SpaceType, Options>::initModel

// template<typename SpaceType, int Options>
// void
// ToolboxMor<SpaceType, Options>::initOnlineModel( std::shared_ptr<super_type> const& model )
// {
//     auto tbModel = std::dynamic_pointer_cast<self_type>(model);
//     this->M_modelProperties = tbModel->modelProperties();
//     this->M_outputDeimName = tbModel->outputDeimName();
// }

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::initOfflineToolbox(std::shared_ptr<DeimMorModelBase<mesh_type>> model )
{
    M_deimMorOfflineModel = model;
    this->setAssembleDEIM(M_deimMorOfflineModel->deimFunction());
    this->setAssembleMDEIM(M_deimMorOfflineModel->mdeimFunction());

    this->initModelImpl();

    this->initOnlineToolbox( model->createOnline() );

    this->postInitModel();
    //this->setInitialized(true);
}

template<typename SpaceType, int Options>
void
ToolboxMor<SpaceType, Options>::initOnlineToolbox(std::shared_ptr<DeimMorModelBase<mesh_type>> model )
{
    M_deimMorOnlineModel = model;
    this->setOnlineAssembleDEIM(M_deimMorOnlineModel->deimOnlineFunction(this->getDEIMReducedMesh()));
    this->setOnlineAssembleMDEIM(M_deimMorOnlineModel->mdeimOnlineFunction(this->getMDEIMReducedMesh()));
    //this->setInitialized(true);
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
    auto outputs = M_modelProperties->outputs().ofTypes({"integrate","mean","sensor","point"});
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
    auto outputs = M_modelProperties->outputs().ofTypes({"integrate","mean","sensor","point"});
    for( auto [output, outp] : enumerate(outputs) )
    {
        auto out = outp.second;
        if( M_outputDeim[out.name()] )
        {
            int M_O = this->deim(M_outputDeim[out.name()])->size();
            auto qo = this->deim(M_outputDeim[out.name()])->q();
            for( int i = 0; i < M_O; ++i )
                this->M_Fqm[output+1][0][i] = qo[i];
        }
        else if( out.type() == "mean" )
        {
            this->M_Fqm[output+1][0][0] = assembleOutputMean(mu, out);
            auto u = this->Xh->element();
            u.on(_range=elements(this->mesh()), _expr=cst(1.));
            LOG(INFO) << fmt::format("[ToolboxMor] f{}(1) = {}", output+1, dot( *this->M_Fqm[output+1][0][0], u));
        }
        else if( out.type() == "integrate")
        {
            this->M_Fqm[output+1][0][0] = assembleOutputIntegrate(mu, out);
        }
        else if( out.type() == "sensor")
        {
            this->M_Fqm[output+1][0][0] = assembleOutputSensor(mu, out);
        }
        else if( out.type() == "point")
        {
            this->M_Fqm[output+1][0][0] = assembleOutputPoint(mu, out);
        }
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
        throw std::logic_error( "[feelpp.mor.ToolboxMor.output] error with output_index : only 0 or 1 " );
    }

    return s ;
}

}
