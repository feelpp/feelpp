//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!

#include "thermoelectric-nl.hpp"

namespace Feel {

po::options_description
ThermoElectricNLBase::makeOptions()
{
    FEELPP_EXPORT po::options_description options( "Thermoelectric" );
    options.add_options()
        ( "thermoelectric.basename", po::value<std::string>()->default_value("thermoelectric-nl"), "basename for the database" )
        ( "thermoelectric.filename", Feel::po::value<std::string>()->default_value("thermoelectric.json"), "json file containing application parameters")
        ( "thermoelectric.gamma", po::value<double>()->default_value( 1e4 ), "penalisation term" )
        ( "thermoelectric.tolerance", po::value<double>()->default_value(1e-8), "tolerance for picard" )
        ( "thermoelectric.maxit", po::value<int>()->default_value(50), "maximum number of iteration for picard" )
        ( "thermoelectric.trainset-eim-size", po::value<int>()->default_value(10), "size of the trainset" )
        ( "thermoelectric.verbose", po::value<int>()->default_value(0), "level of verbosity" )
        ( "thermoelectric.use-deim", po::value<bool>()->default_value(false), "use deim and mdeim" )
        ( "thermoelectric.test-deim", po::value<bool>()->default_value(false), "test deim" )
        ;
    options.add(backend_options("thermo-electro") );
    options.add(backend_options("electro") );
    options.add(backend_options("thermo") );
    options.add(deimOptions("vec")).add(deimOptions("mat"));
    return options;
}

AboutData
ThermoElectricNLBase::makeAbout( std::string const& str )
{
    AboutData about( str.c_str() );
    return about;
}



THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::ThermoElectricNL( std::string prefix ) :
    super_type(soption("thermoelectric.basename"), Environment::worldCommPtr(), prefix),
    M_propertyPath( Environment::expand( soption("thermoelectric.filename")) ),
    M_gamma(doption("thermoelectric.gamma")),
    M_tolerance(doption("thermoelectric.tolerance")),
    M_maxit(ioption("thermoelectric.maxit")),
    M_trainsetEimSize(ioption("thermoelectric.trainset-eim-size")),
    M_verbose(ioption("thermoelectric.verbose")),
    M_useDEIM(boption("thermoelectric.use-deim"))
{
    Feel::cout << "construct model" << std::endl;
    this->addModelFile("property-file", M_propertyPath);
    M_modelProps = std::make_shared<prop_type>(M_propertyPath);
    M_modelProps->enableBoundaryConditions2();

    M_materials = M_modelProps->materials().materialWithPhysic(std::vector<std::string>({"electric","thermic"}));
    M_elecMaterials = M_modelProps->materials().materialWithPhysic("electric");
    M_therMaterials = M_modelProps->materials().materialWithPhysic("thermic");
    M_nbElecMat = M_elecMaterials.size();
    M_nbTherMat = M_therMaterials.size();
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");
    M_nbPotDir = potentialDirichlet.size();
    M_nbTempRobin = temperatureRobin.size();

    auto parameters = M_modelProps->parameters();
    M_sigmaMax = 0;
    for( auto const& [name,mat] : M_elecMaterials )
    {
        auto sigmaP = parameters[mat.getString("misc.sigmaKey")];
        if( sigmaP.hasMinMax() )
        {
            if( sigmaP.min() > M_sigmaMax )
                M_sigmaMax = sigmaP.min();
        }
        else if( sigmaP.value() > M_sigmaMax )
            M_sigmaMax = sigmaP.value();
    }
    Feel::cout << "sigmaMax=" << M_sigmaMax << std::endl;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::Qa( bool isLinear )
{
    if( M_useDEIM && !isLinear )
        return 1;

    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    return M_nbTherMat + M_nbElecMat + potentialDirichlet.size() + temperatureRobin.size();
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::mMaxA(int q)
{
    if( M_useDEIM )
        return this->mdeim()->size();

    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    if( q < M_nbTherMat )
        return this->scalarContinuousEim()[q]->mMax();
    else if( q < M_nbTherMat + M_nbElecMat )
        return this->scalarContinuousEim()[q]->mMax();
    else if( q < M_nbTherMat + M_nbElecMat + potentialDirichlet.size() )
    {
        auto bd = std::next(potentialDirichlet.begin(), q - (M_nbTherMat + M_nbElecMat) )->second;
        return this->scalarContinuousEim()[M_elecEimIndex.at(bd.material())]->mMax();
    }
    else if( q < M_nbTherMat + M_nbElecMat + potentialDirichlet.size() + temperatureRobin.size() )
        return 1;
    else
        return 0;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::Nl()
{
    if( M_useDEIM )
        return 1;

    auto outputs = M_modelProps->outputs().outputsOfType(std::vector<std::string>({"intensity","averageTemp"}));
    return 1 + outputs.size();
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::Ql( int l)
{
    if( l == 0 )
    {
        if( M_useDEIM )
            return 1;

        auto bdConditions = M_modelProps->boundaryConditions2();
        auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
        auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

        return M_nbElecMat + potentialDirichlet.size() + temperatureRobin.size();
    }
    else if( l < Nl() )
    {
        auto outputs = M_modelProps->outputs().outputsOfType(std::vector<std::string>({"intensity","averageTemp"}));
        auto output = std::next(outputs.begin(), l-1)->second;
        if( output.type() == "intensity" )
            return QIntensity( output );
        else if( output.type() == "averageTemp" )
            return QAverageTemp( output );
        else
            return 0;
    }
    else
        return 0;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::mMaxF(int l, int q)
{
    if( l == 0 )
    {
        if( M_useDEIM )
            return this->deim()->size();

        auto bdConditions = M_modelProps->boundaryConditions2();
        auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
        auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

        if( q < M_nbElecMat )
            return this->scalarContinuousEim()[q+M_nbTherMat]->mMax()*this->scalarDiscontinuousEim()[0]->mMax();
        else if( q < M_nbElecMat + potentialDirichlet.size() )
        {
            auto bd = std::next(potentialDirichlet.begin(), q - M_nbElecMat )->second;
            return this->scalarContinuousEim()[M_elecEimIndex.at(bd.material())]->mMax();
        }
        else if( q < M_nbElecMat + potentialDirichlet.size() + temperatureRobin.size() )
            return 1;
        else
            return 0;
    }
    else if( l < Nl() )
    {
        auto outputs = M_modelProps->outputs().outputsOfType(std::vector<std::string>({"intensity","averageTemp"}));
        auto output = std::next(outputs.begin(), l-1)->second;
        if( output.type() == "intensity" )
            return mMaxIntensity( q, output );
        else if( output.type() == "averageTemp" )
            return mMaxAverageTemp( q, output );
        else
            return 0;
    }
    else
        return 0;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::QIntensity( ModelOutput const& out ) const
{
    return 1;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::QAverageTemp( ModelOutput const& out ) const
{
    return 1;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::mMaxIntensity(int q, ModelOutput const& out ) const
{
    auto mat = out.getString("material");
    auto eimSigma = this->scalarContinuousEim()[M_elecEimIndex.at(mat)];
    return eimSigma->mMax();
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::mMaxAverageTemp(int q, ModelOutput const& out ) const
{
    return 1;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::resize()
{
    if( !M_useDEIM )
        M_eimGradSize = this->scalarDiscontinuousEim()[0]->mMax();
    this->M_Aqm.resize(Qa());
    this->M_betaAqm.resize(Qa());
    for( int q = 0; q < Qa(); ++q)
    {
        this->M_Aqm[q].resize(mMaxA(q), backend()->newMatrix(this->Xh, this->Xh ) );
        this->M_betaAqm[q].resize(mMaxA(q));
    }
    this->M_Fqm.resize(Nl());
    this->M_betaFqm.resize(Nl());
    for( int l = 0; l < Nl(); ++l )
    {
        this->M_Fqm[l].resize(Ql(l));
        this->M_betaFqm[l].resize(Ql(l));
        for( int q = 0; q < Ql(l); ++q )
        {
            this->M_Fqm[l][q].resize(mMaxF(l,q), backend()->newVector(this->Xh) );
            this->M_betaFqm[l][q].resize(mMaxF(l,q));
        }
    }
}


/*************************** mesh support of functionspace **********************/
THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::functionspace_type::mesh_support_vector_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::functionspaceMeshSupport( mesh_ptrtype const& mesh ) const
{
    auto elecDomain = markedelements(mesh, M_materials.markersWithPhysic("electric"));
    auto therDomain = markedelements(mesh, M_materials.markersWithPhysic("thermic"));
    auto suppElec = std::make_shared<MeshSupport<mesh_type>>(mesh,elecDomain);
    auto suppTher = std::make_shared<MeshSupport<mesh_type>>(mesh,therDomain);
    return fusion::make_vector(suppElec,suppTher);
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::initModel()
{
    Feel::cout << "init model" << std::endl;
    auto parameters = M_modelProps->parameters();
    int nbCrbParameters = count_if(parameters.begin(), parameters.end(), [] (auto const& p)
                                   {
                                       return p.second.hasMinMax();
                                   });
    this->Dmu->setDimension(nbCrbParameters);
    auto mu_min = this->Dmu->element();
    auto mu_max = this->Dmu->element();
    int i = 0;
    for( auto const& [key,parameter] : parameters )
    {
        if( parameter.hasMinMax() )
        {
            mu_min(i) = parameter.min();
            mu_max(i) = parameter.max();
            Feel::cout << "Found parameter " << key
                       << " with min/max : " << mu_min(i) << "/" << mu_max(i) << std::endl;
            this->Dmu->setParameterName(i++, key );
        }
    }
    this->Dmu->setMin(mu_min);
    this->Dmu->setMax(mu_max);
    M_mu = this->Dmu->element();

    if( !M_mesh )
        M_mesh = loadMesh( new mesh_type );
    this->setFunctionSpaces(functionspace_type::New( _mesh=M_mesh, _range=this->functionspaceMeshSupport( M_mesh ) ) );
    Feel::cout << "Potential nDof  : " << this->Xh->template functionSpace<0>()->nDof() << std::endl
               << "Temperature nDof: " << this->Xh->template functionSpace<1>()->nDof() << std::endl;

    M_u = this->Xh->element();

    auto mu = this->Dmu->min();
    M_InitialGuess.resize(QInitialGuess());
    this->M_betaInitialGuess.resize(2);
    for( int q = 0; q < QInitialGuess(); ++q )
    {
        M_InitialGuess[q].resize(mMaxInitialGuess(q), this->Xh->elementPtr());
        this->M_betaInitialGuess[q].resize(mMaxInitialGuess(q));
    }
    auto VTInit = this->solveLinear(mu);
    M_InitialGuess[0][0]->template element<0>() = VTInit.template element<0>();
    M_InitialGuess[1][0]->template element<1>() = VTInit.template element<1>();
#if 1
    auto ex = exporter(_mesh=M_mesh, _name="initial-guess");
    ex->add("V_Init",VTInit.template element<0>() );
    ex->add("T_Init",VTInit.template element<1>() );
    ex->save();
#endif
    auto Pset = this->Dmu->sampling();

    std::string supersamplingname = (boost::format("DmuEim-Ne%1%-generated-by-master-proc")
                                     % M_trainsetEimSize ).str();
    std::ifstream file ( supersamplingname );
    bool all_proc_same_sampling=true;
    if( ! file )
    {
        Pset->randomize( M_trainsetEimSize , all_proc_same_sampling , supersamplingname );
        Pset->writeOnFile( supersamplingname );
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(supersamplingname);
    }

    if( M_useDEIM )
    {
        auto d = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _sampling=Pset, _prefix="vec");
        this->addDeim(d);
        this->deim()->run();
        Feel::cout << tc::green << "DEIM construction finished!!" << tc::reset << std::endl;

        auto m = Feel::mdeim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _sampling=Pset, _prefix="mat");
        this->addMdeim(m);
        this->mdeim()->run();
        Feel::cout << tc::green << "MDEIM construction finished!!" << tc::reset << std::endl;
        if( boption("thermoelectric.test-deim") )
        {
            auto musV = this->deim()->mus();
            auto musM = this->mdeim()->mus();
            musV.insert(musV.end(), musM.begin(), musM.end());
            auto qa = this->mdeim()->q();
            auto qf = this->deim()->q();
            for( auto const& mu : musV )
            {
                auto VTFE = this->solve(mu);
                auto A = this->assembleForMDEIMnl(mu, VTFE, 0 );
                auto F = this->assembleForDEIMnl(mu, VTFE, 0 );

                auto betaA = this->mdeim()->beta(mu, VTFE);
                auto betaF = this->deim()->beta(mu, VTFE);
                auto Am = backend()->newMatrix(this->Xh,this->Xh);
                Am->zero();
                for( int m = 0; m < this->mdeim()->size(); ++m )
                {
                    Am->addMatrix( betaA(m), qa[m]);
                    Feel::cout << "betaA(" << m << ") = " << betaA(m) << std::endl;
                }
                auto Fm = backend()->newVector(this->Xh);
                Fm->zero();
                for( int m = 0; m < this->deim()->size(); ++m )
                {
                    Fm->add( betaF(m), qf[m]);
                    Feel::cout << "betaF(" << m << ") = " << betaF(m) << std::endl;
                }

                auto VTEIM = this->Xh->element();
                backend()->solve( _matrix=Am, _rhs=Fm, _solution=VTEIM);

                auto VEIM = VTEIM.template element<0>();
                auto TEIM = VTEIM.template element<1>();
                auto VFE = VTFE.template element<0>();
                auto TFE = VTFE.template element<1>();

                auto e = exporter(_mesh=M_mesh,_name="testdeim");
                e->add("VFE", VFE);
                e->add("VEIM", VEIM);
                e->add("TFE", TFE);
                e->add("TEIM", TEIM);
                e->save();

                auto errV = normL2(elements(M_mesh), idv(VFE)-idv(VEIM) );
                auto normV = normL2(elements(M_mesh), idv(VFE) );
                Feel::cout << "V: err = " << errV << " relative err = " << errV/normV << std::endl;
                auto errT = normL2(elements(M_mesh), idv(TFE)-idv(TEIM) );
                auto normT = normL2(elements(M_mesh), idv(TFE) );
                Feel::cout << "T: err = " << errT << " relative err = " << errT/normT << std::endl;

                Am->addMatrix(-1., A);
                Fm->add(-1., F);
                auto errA = Am->linftyNorm();
                auto normA = A->linftyNorm();
                auto errF = Fm->linftyNorm();
                auto normF = F->linftyNorm();
                Feel::cout << "A: err = " << errA << " relative err = " << errA/normA << std::endl
                           << "F: err = " << errF << " relative err = " << errF/normF << std::endl;
            }
        }
    }
    else
    {
        // eim
        auto T0 = cst(293.0);
        auto L = cst_ref(M_mu.parameterNamed("L"));
        i = 0;
        for( auto const& [key,mat] : M_therMaterials )
        {
            auto name = "eim_k_"+key;
            auto sigma0 = cst_ref(M_mu.parameterNamed(mat.getString("misc.sigmaKey")));
            auto alpha = cst_ref(M_mu.parameterNamed(mat.getString("misc.alphaKey")));
            auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
            auto k = sigma*L*_e1;
            auto Eh = eim_space_type::New( _mesh=M_mesh, _range=markedelements(M_mesh, mat.meshMarkers()) );

            auto eim_k = eim( _model=std::dynamic_pointer_cast<ThermoElectricNL>(this->shared_from_this() ),
                              _element=M_u.template element<1>(),
                              _parameter=M_mu,
                              _expr=k,
                              _space=Eh,
                              _name=name,
                              _sampling=Pset );
            this->addEim( eim_k );
            M_therEimIndex[key] = i++;
            Feel::cout << tc::green << name << " dimension: " << eim_k->mMax() << tc::reset << std::endl;
        }

        for( auto const& [key,mat] : M_elecMaterials )
        {
            auto name = "eim_sigma_"+key;
            auto sigma0 = cst_ref(M_mu.parameterNamed(mat.getString("misc.sigmaKey")));
            auto alpha = cst_ref(M_mu.parameterNamed(mat.getString("misc.alphaKey")));
            auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
            auto Eh = eim_space_type::New( _mesh=M_mesh, _range=markedelements(M_mesh, mat.meshMarkers()) );
            auto eim_sigma = eim( _model=std::dynamic_pointer_cast<ThermoElectricNL>(this->shared_from_this() ),
                                  _element=M_u.template element<1>(),
                                  _parameter=M_mu,
                                  _expr=sigma/M_sigmaMax,
                                  _space=Eh,
                                  _name=name,
                                  _sampling=Pset );
            this->addEim( eim_sigma );
            M_elecEimIndex[key] = i++;
            Feel::cout << tc::green << name << " dimension: " << eim_sigma->mMax() << tc::reset << std::endl;
        }

        auto name = "eim_grad";
        auto gradgrad = _e2v*trans(_e2v);
        auto Jh = eimd_space_type::New( _mesh=M_mesh,
                                        _range=elements(support(this->Xh->template functionSpace<0>())));
        auto eim_grad = eim( _model=std::dynamic_pointer_cast<ThermoElectricNL>(this->shared_from_this() ),
                             _element=M_u.template element<1>(),
                             _element2=M_u.template element<0>(),
                             _parameter=M_mu,
                             _expr=gradgrad,
                             _space=Jh,
                             _name=name,
                             _sampling=Pset );
        this->addEimDiscontinuous( eim_grad );
        M_eimGradSize = eim_grad->mMax();
        Feel::cout << tc::green << name << " dimension: " << eim_grad->mMax() << tc::reset << std::endl;
    }

    auto U = this->Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();

    // Energy matrix
    auto m = form2(_test=this->Xh, _trial=this->Xh);
    m = integrate( elements(support(this->Xh->template functionSpace<0>())),
                   inner(gradt(u1),grad(u1)) );
    m += integrate( elements(support(this->Xh->template functionSpace<1>())),
                    inner(gradt(u2),grad(u2)) );
    this->M_energy_matrix = m.matrixPtr();

    auto mm = form2(_test=this->Xh, _trial=this->Xh);
    mm = integrate( elements(support(this->Xh->template functionSpace<0>())),
                    inner(idt(u1),id(u1)) );
    mm += integrate( elements(support(this->Xh->template functionSpace<1>())),
                     inner(idt(u2),id(u2)) );
    this->M_mass_matrix = mm.matrixPtr();

    this->assemble();
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    Feel::cout << "setup specificity model" << std::endl;
    auto const& ptreeEim = ptree.get_child( "eim" );
    this->clearScalarContinuousEim();
    this->clearScalarDiscontinuousEim();
    auto T0 = cst(293.0);
    auto L = cst_ref(M_mu.parameterNamed("L"));
    int i = 0;
    for( auto const& [key,mat] : M_therMaterials )
    {
        auto name = "eim_k_"+key;
        auto sigma0 = cst_ref(M_mu.parameterNamed(mat.getString("misc.sigmaKey")));
        auto alpha = cst_ref(M_mu.parameterNamed(mat.getString("misc.alphaKey")));
        auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
        auto k = sigma*L*_e1;
        eim_space_ptrtype Eh;
        auto const& ptreeEimK = ptreeEim.get_child( name );
        std::string dbNameEimK = ptreeEimK.template get<std::string>( "database-filename" );
        auto eim_k = eim( _model=std::dynamic_pointer_cast<ThermoElectricNL>(this->shared_from_this() ),
                          _element=M_u.template element<1>(),
                          _parameter=M_mu,
                          _expr=k,
                          _space=Eh,
                          _name=name,
                          _filename=dbNameEimK,
                          _directory=dbDir );
        this->addEim( eim_k );
        M_therEimIndex[key] = i++;
        Feel::cout << tc::green << name << " dimension: " << eim_k->mMax() << tc::reset << std::endl;
    }

    for( auto const& [key,mat] : M_elecMaterials )
    {
        auto name = "eim_sigma_"+key;
        auto sigma0 = cst_ref(M_mu.parameterNamed(mat.getString("misc.sigmaKey")));
        auto alpha = cst_ref(M_mu.parameterNamed(mat.getString("misc.alphaKey")));
        auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
        eim_space_ptrtype Eh;
        auto const& ptreeEimSigma = ptreeEim.get_child( name );
        std::string dbNameEimSigma = ptreeEimSigma.template get<std::string>( "database-filename" );
        auto eim_sigma = eim( _model=std::dynamic_pointer_cast<ThermoElectricNL>(this->shared_from_this() ),
                              _element=M_u.template element<1>(),
                              _parameter=M_mu,
                              _expr=sigma/M_sigmaMax,
                              _space=Eh,
                              _name=name,
                              _filename=dbNameEimSigma,
                              _directory=dbDir );
        this->addEim( eim_sigma );
        M_elecEimIndex[key] = i++;
        Feel::cout << tc::green << name << " dimension: " << eim_sigma->mMax() << tc::reset << std::endl;
    }

    auto name = "eim_grad";
    auto gradgrad = _e2v*trans(_e2v);
    auto const& ptreeEimGrad = ptreeEim.get_child( name );
    std::string dbNameEimGrad = ptreeEimGrad.template get<std::string>( "database-filename" );
    eimd_space_ptrtype Jh;
    auto eim_grad = eim( _model=std::dynamic_pointer_cast<ThermoElectricNL>(this->shared_from_this() ),
                         _element=M_u.template element<1>(),
                         _element2=M_u.template element<0>(),
                         _parameter=M_mu,
                         _expr=gradgrad,
                         _space=Jh,
                         _name=name,
                         _filename=dbNameEimGrad,
                         _directory=dbDir );
    this->addEimDiscontinuous( eim_grad );
    M_eimGradSize = eim_grad->mMax();
    Feel::cout << tc::green << name << " dimension: " << eim_grad->mMax() << tc::reset << std::endl;

    Feel::cout << "eim size = " << this->scalarContinuousEim().size() << std::endl;


    this->resize();

    auto outputs = M_modelProps->outputs().outputsOfType("point");
    bool hasCtxElectro = false, hasCtxThermo = false;
    for( auto const& [name,output] : outputs )
    {
        if( output.getString("field") == "electric-potential" )
            hasCtxElectro = true;
        else if( output.getString("field") == "temperature" )
            hasCtxThermo = true;
    }
    if( hasCtxElectro )
    {
        fs::path dbdir(this->M_crbModelDb.dbRepository());
        std::string filename = ( dbdir / "ctxelectro.crbdb").string();
        if ( this->worldComm().isMasterRank() )
        {
            fs::ifstream ifs( filename );
            if ( ifs )
            {
                boost::archive::binary_iarchive ia( ifs );
                ia >> M_ctxElectro;
                M_ctxElectro.setRbFunctionSpace(this->rBFunctionSpace()->template rbFunctionSpace<0>());
            }
        }
    }
    if( hasCtxThermo )
    {
        fs::path dbdir(this->M_crbModelDb.dbRepository());
        std::string filename = ( dbdir / "ctxthermo.crbdb").string();
        if ( this->worldComm().isMasterRank() )
        {
            fs::ifstream ifs( filename );
            if ( ifs )
            {
                boost::archive::binary_iarchive ia( ifs );
                ia >> M_ctxThermo;
                M_ctxThermo.setRbFunctionSpace(this->rBFunctionSpace()->template rbFunctionSpace<1>());
            }
        }
    }
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::initOutputsPoints()
{
    Feel::cout << "init output points" << std::endl;
    auto outputs = M_modelProps->outputs().outputsOfType("point");
    M_ctxElectro = this->rBFunctionSpace()->template rbFunctionSpace<0>()->context();
    M_ctxThermo = this->rBFunctionSpace()->template rbFunctionSpace<1>()->context();
    bool hasCtxElectro = false, hasCtxThermo = false;
    for( auto const& [name,output] : outputs )
    {
        node_type t(Dim);
        auto coord = expr<Dim,1>(output.getString("coord")).evaluate();
        for( int i = 0; i < Dim; ++i )
            t(i) = coord(i);
        if( output.getString("field") == "electric-potential" )
        {
            hasCtxElectro = true;
            M_ctxElectro.add( t );
        }
        else if( output.getString("field") == "temperature" )
        {
            hasCtxThermo = true;
            M_ctxThermo.add( t );
        }
    }
    if( hasCtxElectro )
    {
        M_ctxElectro.update();

        fs::path dbdir(this->M_crbModelDb.dbRepository());
        std::string filename = ( dbdir / "ctxelectro.crbdb").string();
        this->addModelFile("context-points-electro", filename);
        if ( this->worldComm().isMasterRank() )
        {
            fs::ofstream ofs( filename );
            if ( ofs )
            {
                boost::archive::binary_oarchive oa( ofs );
                oa << M_ctxElectro;
            }
        }
    }
    if( hasCtxThermo )
    {
        M_ctxThermo.update();

        fs::path dbdir(this->M_crbModelDb.dbRepository());
        std::string filename = ( dbdir / "ctxthermo.crbdb").string();
        this->addModelFile("context-points-thermo", filename);
        if ( this->worldComm().isMasterRank() )
        {
            fs::ofstream ofs( filename );
            if ( ofs )
            {
                boost::archive::binary_oarchive oa( ofs );
                oa << M_ctxThermo;
            }
        }
    }
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::vectorN_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::computeOutputsPointsElectro( vectorN_type const& urb )
{
    auto VRbElt = M_ctxElectro.rbFunctionSpace()->element();
    int dimRb = std::min((int)VRbElt.size(),(int)urb.size());
    for ( int k=0; k<dimRb; ++k )
        VRbElt(k) = urb(k);
    auto evaluations = evaluateFromContext( _context=M_ctxElectro , _expr=idv(VRbElt) );
    return evaluations;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::vectorN_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::computeOutputsPointsThermo( vectorN_type const& urb )
{
    auto TRbElt = M_ctxThermo.rbFunctionSpace()->element();
    int dimRb = std::min((int)TRbElt.size(),(int)urb.size());
    for ( int k=0; k<dimRb; ++k )
        TRbElt(k) = urb(k);
    auto evaluations = evaluateFromContext( _context=M_ctxThermo , _expr=idv(TRbElt) );
    return evaluations;
}

// #include <feel/feelcrb/crbplugin.hpp>
// FEELPP_CRB_PLUGIN( ThermoElectricNL, thermoelectricnl )
}
