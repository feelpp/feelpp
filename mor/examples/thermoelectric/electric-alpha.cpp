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

#include <feel/feelmor/crbplugin.hpp>
#include "electric-alpha.hpp"

#include <boost/math/special_functions/binomial.hpp>

namespace Feel {

po::options_description
AlphaElectric::makeOptions( std::string const& prefix )
{
    po::options_description options( "Electric" );
    options.add_options()
        ( "thermoelectric.filename", Feel::po::value<std::string>()->default_value("electric.json"),
          "json file containing application parameters and boundary conditions")
        ( "thermoelectric.penal-dir", po::value<double>()->default_value( 1e5 ), "penalisation term" )
        ( "thermoelectric.trainset-deim-size", po::value<int>()->default_value(40), "size of the deim trainset" )
        ( "thermoelectric.trainset-mdeim-size", po::value<int>()->default_value(40), "size of the mdeim trainset" )
        ( "thermoelectric.test-deim", po::value<bool>()->default_value(false), "test deim interpolation" )
        ( "thermoelectric.db.base", po::value<std::string>()->default_value("alphaelectric"), "database basename" )
        ( "thermoelectric.verbose", po::value<int>()->default_value(0), "verbosity level" )
        ;
    options.add(backend_options("feV") );
    options.add(deimOptions("vec")).add(deimOptions("mat"));
    options.add(crbOptions(prefix));
    options.add(crbSEROptions(prefix));
    options.add(eimOptions());
    options.add(podOptions());
    options.add(backend_options("backend-primal"));
    options.add(backend_options("backend-dual"));
    options.add(backend_options("backend-l2"));

    return options;
}

AlphaElectric::AlphaElectric(std::string const& prefix)
    : AlphaElectric(nullptr, prefix)
{}

AlphaElectric::AlphaElectric( mesh_ptrtype mesh,  std::string const& prefix )
    : super_type( soption("thermoelectric.db.base"), Environment::worldCommPtr(), prefix),
      M_trainsetDeimSize(ioption("thermoelectric.trainset-deim-size")),
      M_trainsetMdeimSize(ioption("thermoelectric.trainset-mdeim-size")),
      M_penalDir(doption("thermoelectric.penal-dir")),
      M_propertyPath(Environment::expand( soption("thermoelectric.filename"))),
      M_testDeim(boption("thermoelectric.test-deim")),
      M_dbBasename(soption("thermoelectric.db.base")),
      M_verbose(ioption("thermoelectric.verbose"))
{
    M_modelProps = std::make_shared<prop_type>(M_propertyPath);

    auto parameters = M_modelProps->parameters();
    int nbCrbParameters = count_if(parameters.begin(), parameters.end(), [] (auto const& p)
                                   {
                                       return p.second.hasMinMax();
                                   });
    Dmu->setDimension(nbCrbParameters);
    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    int i = 0;
    for( auto const& parameterPair : parameters )
    {
        if( parameterPair.second.hasMinMax() )
        {
            mu_min(i) = parameterPair.second.min();
            mu_max(i) = parameterPair.second.max();
            Dmu->setParameterName(i++, parameterPair.first );
        }
    }
    Dmu->setMin(mu_min);
    Dmu->setMax(mu_max);
    M_mu = Dmu->element();

    // keep only materials with physic electric
    M_materials = M_modelProps->materials().materialWithPhysic("electric-geo");
    for( auto const& mp : M_materials )
    {
        if( mp.second.hasProperty("zmin")
            && mp.second.hasProperty("zmax")
            && mp.second.hasProperty("params") )
            M_materialsWithGeo.insert(mp);
        else
            M_materialsWithoutGeo.insert(mp);
    }

    this->M_mesh = mesh;
}

int AlphaElectric::Qa()
{
    return 1;
}

int AlphaElectric::mQA( int q )
{
    return this->mdeim()->size();
}

int AlphaElectric::Nl()
{
    return 1;
}

int AlphaElectric::Ql( int l)
{
    return 1;
}

int AlphaElectric::mLQF( int l, int q )
{
    switch( l )
    {
    case 0:
        return mCompliantQ(q);
    default:
        return 0;
    }
}

int AlphaElectric::mCompliantQ(int q )
{
    return this->deim()->size();
}

void AlphaElectric::resizeQm( bool resizeMat )
{
    if( resizeMat )
        M_Aqm.resize( Qa());
    M_betaAqm.resize( Qa() );
    for( int q = 0; q < Qa(); ++q )
    {
        if( resizeMat )
            M_Aqm[q].resize(mQA(q), backend()->newMatrix(Xh, Xh ) );
        M_betaAqm[q].resize(mQA(q));
    }

    if( resizeMat )
        M_Fqm.resize(Nl());
    M_betaFqm.resize(Nl());
    for( int l = 0; l < Nl(); ++l )
    {
        if( resizeMat )
            M_Fqm[l].resize(Ql(l));
        M_betaFqm[l].resize(Ql(l));
        for( int q = 0; q < Ql(l); ++q )
        {
            if( resizeMat )
                M_Fqm[l][q].resize(mLQF(l, q), backend()->newVector(Xh) );
            M_betaFqm[l][q].resize(mLQF(l, q) );
        }
    }

    if( resizeMat )
    {
        M_InitialGuess.resize(1);
        M_InitialGuess[0].resize(1);
        M_InitialGuess[0][0] = Xh->elementPtr();
    }
}

AlphaElectric::parameter_type AlphaElectric::paramFromVec( std::vector<double> const& x )
{
    auto mu = this->newParameter();
    mu = Eigen::Map<Eigen::VectorXd const>( x.data(), mu.size());
    return mu;
}

AlphaElectric::parameter_type
AlphaElectric::paramFromProperties() const
{
    auto mu = Dmu->element();
    int i = 0;
    auto parameters = M_modelProps->parameters();
    for( auto const& parameterPair : parameters )
    {
        if( parameterPair.second.hasMinMax() )
        {
            mu(i++) = parameterPair.second.value();
        }
    }
    return mu;
}

AlphaElectric::parameter_type AlphaElectric::param0()
{
    std::vector<double> x(Dmu->dimension(), 0.);
    return this->paramFromVec(x);
}

std::string AlphaElectric::alphaRef( parameter_type const& mu, ModelMaterial const& material )
{
    using boost::math::binomial_coefficient;

    double zb = material.getDouble("zmin");
    double zt = material.getDouble("zmax");
    auto paramNames = material.getVecString("params");

    int n = paramNames.size() + 1;

    std::stringstream alphaStream;
    alphaStream << "0 + (z<" << zt << ")*(z>" << zb << ")*(";
    for( int k = 1; k < n; ++k )
    {
        double ckn = binomial_coefficient<double>(n,k);
        alphaStream << " + " << ckn
                    << "*((z-" << zb << ")/" << zt-zb << ")^" << k
                    << "*(1-(z-" << zb << ")/" << zt-zb << ")^" << n-k
                    << "*" << paramNames[k-1];
    }
    alphaStream << "):x:y:z";
    for( auto const& p : paramNames )
        alphaStream << ":" << p;

    return alphaStream.str();
}

std::string AlphaElectric::alpha( ModelMaterial const& material )
{
    using boost::math::binomial_coefficient;

    double zb = material.getDouble("zmin");
    double zt = material.getDouble("zmax");
    auto paramNames = material.getVecString("params");

    std::vector<double> p( 2, 0);
    for( int i = 0; i < paramNames.size(); ++i )
        p.insert( p.begin() +i+1, mu.parameterNamed(paramNames[i]) );
    int n = p.size() - 1;

    std::stringstream alphaStream;
    alphaStream << "0 + (z<" << zt << ")*(z>" << zb << ")*(";
    for( int k = 1; k < n; ++k )
    {
        double ckn = binomial_coefficient<double>(n,k);
        alphaStream << " + " << ckn
                    << "*((z-" << zb << ")/" << zt-zb << ")^" << k
                    << "*(1-(z-" << zb << ")/" << zt-zb << ")^" << n-k
                    << "*" << p[k];
    }
    alphaStream << "):x:y:z";

    return alphaStream.str();
}

std::string AlphaElectric::alphaPrimeRef( ModelMaterial const& material )
{
    using boost::math::binomial_coefficient;

    double zb = material.getDouble("zmin");
    double zt = material.getDouble("zmax");
    auto paramNames = material.getVecString("params");
    paramNames.insert(paramNames.begin(), "0");
    paramNames.push_back("0");

    int n = paramNames.size() - 1;

    std::stringstream alphaPrimeStream;
    alphaPrimeStream << "0 + (z<" << zt << ")*(z>" << zb << ")*(";
    alphaPrimeStream << n << "/" << zt-zb << "*(";
    for( int k = 0; k < n; ++k )
    {
        double ckn1 = binomial_coefficient<double>(n-1,k);
        alphaPrimeStream << " + (" << paramNames[k+1] << "-" << paramNames[k] << ")*" << ckn1
                         << "*((z-" << zb << ")/" << zt-zb << ")^" << k
                         << "*(1-(z-" << zb << ")/" << zt-zb << ")^" << n-1-k;
    }
    alphaPrimeStream << ")):x:y:z";
    for( int k = 1; k < n; ++k )
        alphaPrimeStream << ":" << paramNames[k];

    return alphaPrimeStream.str();
}

std::string AlphaElectric::alphaPrime( parameter_type const& mu, ModelMaterial const& material )
{
    using boost::math::binomial_coefficient;

    double zb = material.getDouble("zmin");
    double zt = material.getDouble("zmax");
    auto paramNames = material.getVecString("params");

    std::vector<double> p( 2, 0);
    for( int i = 0; i < paramNames.size(); ++i )
        p.insert( p.begin() +i+1, mu.parameterNamed(paramNames[i]) );
    int n = p.size() - 1;

    std::stringstream alphaPrimeStream;
    alphaPrimeStream << "0 + (z<" << zt << ")*(z>" << zb << ")*(";
    alphaPrimeStream << n << "/" << zt-zb << "*(";
    for( int k = 0; k < n; ++k )
    {
        double ckn1 = binomial_coefficient<double>(n-1,k);
        alphaPrimeStream << " + " << (p[k+1]-p[k]) << "*" << ckn1
                         << "*((z-" << zb << ")/" << zt-zb << ")^" << k
                         << "*(1-(z-" << zb << ")/" << zt-zb << ")^" << n-1-k;
    }
    alphaPrimeStream << ")):x:y:z";

    return alphaPrimeStream.str();
}

AlphaElectric::sparse_matrix_ptrtype
AlphaElectric::assembleForMDEIM( parameter_type const& mu, int const& tag )
{
    auto mesh = Xh->mesh();
    auto V = Xh->element();
    auto phiV = Xh->element();

    auto bc = M_modelProps->boundaryConditions();
    double sigmaMax = 0;
    for( auto const& material : M_materials )
    {
        auto sigma = material.second.getDouble("sigma");
        if( sigma > sigmaMax )
            sigmaMax = sigma;
    }

    auto a = form2(_test=Xh, _trial=Xh);
    for( auto const& material : M_materialsWithGeo )
    {
        auto sigma = material.second.getDouble("sigma");

        auto alphaExpr = expr(this->alpha(mu, material.second));
        auto alphaPrimeExpr = expr(this->alphaPrime(mu, material.second));

        auto Jinv = mat<3,3>( cos(alphaExpr), -sin(alphaExpr), -alphaPrimeExpr*Py(),
                              sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*Px(),
                              cst(0.), cst(0.), cst(1.) );

        auto gradVJ = gradt(V)*Jinv;
        auto gradPhiVJ = grad(phiV)*Jinv;

        a += integrate( markedelements(mesh, material.first),
                        sigma/sigmaMax*inner(gradVJ,gradPhiVJ) );
    }
    for( auto const& material : M_materialsWithoutGeo )
        a += integrate( markedelements(mesh, material.first),
                        material.second.getDouble("sigma")/sigmaMax*inner(gradt(V),grad(phiV)) );
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        if( isInMaterials(exAtM) )
        {
            auto sigma = M_materials[exAtM.material()].getDouble("sigma");
            a += integrate( markedfaces(mesh, exAtM.marker()),
                            sigma/sigmaMax*(M_penalDir/hFace()*inner(idt(V),id(phiV))
                                            -inner(grad(V)*N(),id(phiV))
                                            -inner(grad(phiV)*N(),idt(V)) ) );
        }
    }

    auto am = a.matrixPtr();
    am->close();

    return am;
}

AlphaElectric::vector_ptrtype
AlphaElectric::assembleForDEIM( parameter_type const& mu, int const& tag )
{
    auto mesh = Xh->mesh();
    auto phiV = Xh->element();

    auto bc = M_modelProps->boundaryConditions();
    double sigmaMax = 0;
    for( auto const& material : M_materials )
    {
        auto sigma = material.second.getDouble("sigma");
        if( sigma > sigmaMax )
            sigmaMax = sigma;
    }

    auto f = form1(_test=Xh);
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        if( isInMaterials(exAtM) )
        {
            auto sigma = M_materials[exAtM.material()].getDouble("sigma");
            auto dif = expr(exAtM.expression());
            f += integrate( markedfaces(mesh, exAtM.marker()),
                            sigma/sigmaMax*dif*(M_penalDir/hFace()*id(phiV) - grad(phiV)*N()) );
        }
    }
    auto fv = f.vectorPtr();
    fv->close();

    return fv;
}

void AlphaElectric::initModel()
{
    this->addModelFile("property-file", M_propertyPath);

    if( !M_mesh )
        M_mesh = loadMesh( new mesh_type );
    std::vector<std::string> range;
    for( auto const& mp : M_materials )
        range.push_back(mp.first);
    auto domain = markedelements(M_mesh, range);
    this->setFunctionSpaces(functionspace_type::New( _mesh=M_mesh, _range=domain ) );

    Feel::cout << "initModel with nDof: " << Xh->nDof() << std::endl;

    if( !pT )
        pT = element_ptrtype( new element_type( Xh ) );

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

    auto d = Feel::deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _sampling=PsetV, _prefix="vec");
    this->addDeim(d);
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

    auto m = Feel::mdeim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()), _sampling=PsetM, _prefix="mat");
    this->addMdeim(m);
    this->mdeim()->run();
    Feel::cout << tc::green << "Electric MDEIM construction finished!!" << tc::reset << std::endl;

    this->resizeQm();
    this->decomposition();

    if( M_testDeim )
    {
        auto VFE = this->solve(M_mu);
        auto A = this->assembleForMDEIM(M_mu, 0 );
        auto F = this->assembleForDEIM(M_mu, 0 );

        auto betaA = this->mdeim()->beta(M_mu);
        auto betaF = deim()->beta(M_mu);
        auto qa = mdeim()->q();
        auto qf = deim()->q();
        auto Am = backend()->newMatrix(Xh,Xh);
        Am->zero();
        for( int m = 0; m < mdeim()->size(); ++m )
        {
            Am->addMatrix( betaA(m), qa[m]);
            Feel::cout << "betaA(" << m << ") = " << betaA(m) << std::endl;
        }
        auto Fm = backend()->newVector(Xh);
        Fm->zero();
        for( int m = 0; m < deim()->size(); ++m )
        {
            Fm->add( betaF(m), qf[m]);
            Feel::cout << "betaF(" << m << ") = " << betaF(m) << std::endl;
        }

        auto VEIM = Xh->element();
        backend()->solve( _matrix=Am, _rhs=Fm, _solution=VEIM);

        auto e = exporter(M_mesh);
        e->add("VFE", VFE);
        e->add("VEIM", VEIM);
        e->save();

        auto errV = normL2(domain, idv(VFE)-idv(VEIM) );
        auto normV = normL2(domain, idv(VFE) );
        Feel::cout << "V: err = " << errV << " relative err = " << errV/normV << std::endl;

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

void AlphaElectric::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    std::string propertyPath;
    if( this->hasModelFile("property-file") )
        propertyPath = this->additionalModelFiles().find("property-file")->second;
    else
        Feel::cerr << "Warning!! the database does not contain the property file! Expect bugs!"
                   << std::endl;
    M_modelProps = std::make_shared<prop_type>(propertyPath);

    auto parameters = M_modelProps->parameters();
    int nbCrbParameters = count_if(parameters.begin(), parameters.end(), [] (auto const& p)
                                   {
                                       return p.second.hasMinMax();
                                   });
    Dmu = parameterspace_type::New( nbCrbParameters );

    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    int i = 0;
    for( auto const& parameterPair : parameters )
    {
        if( parameterPair.second.hasMinMax() )
        {
            mu_min(i) = parameterPair.second.min();
            mu_max(i) = parameterPair.second.max();
            Dmu->setParameterName(i++, parameterPair.first );
        }
    }
    Dmu->setMin(mu_min);
    Dmu->setMax(mu_max);
    M_mu = Dmu->element();

    std::shared_ptr<J_space_type> JspaceEim;
    if ( !pT )
        pT.reset( new element_type );

    // auto const& ptreeEim = ptree.get_child( "eim" );
    // auto const& ptreeEimJoule = ptreeEim.get_child( "eim_joule" );
    // std::string dbnameEimJoule = ptreeEimJoule.template get<std::string>( "database-filename" );

    this->resizeQm( false );
}

void AlphaElectric::decomposition()
{
    int M_A = this->mdeim()->size();
    auto qa = this->mdeim()->q();
    for( int i = 0; i < M_A; ++i )
        M_Aqm[0][i] = qa[i];

    int M_F = this->deim()->size();
    auto qf = this->deim()->q();
    for( int i = 0; i < M_F; ++i )
        M_Fqm[0][0][i] = qf[i];

    // Energy matrix
    auto V = Xh->element();
    auto phiV = Xh->element();
    auto m = form2(_test=Xh, _trial=Xh);
    m = integrate( Xh->dof()->meshSupport()->rangeElements(),
                   gradt(V)*trans(grad(phiV)) );
    M_energy_matrix = m.matrixPtr();
}

AlphaElectric::beta_vector_type
AlphaElectric::computeBetaInitialGuess( parameter_type const& mu )
{
    M_betaInitialGuess.resize( 1 );
    M_betaInitialGuess[0].resize( 1 );
    M_betaInitialGuess[0][0] = 1;
    return this->M_betaInitialGuess;
}

AlphaElectric::beta_type
AlphaElectric::computeBetaQm( element_type const& u, parameter_type const& mu )
{
    auto betaA = this->mdeim()->beta(mu);
    auto betaF = this->deim()->beta(mu);
    this->fillBetaQm(mu, betaA, betaF);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

AlphaElectric::beta_type
AlphaElectric::computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
{
    auto betaA = this->mdeim()->beta(mu);
    auto betaF = this->deim()->beta(mu);
    this->fillBetaQm(mu, betaA, betaF);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

AlphaElectric::beta_type
AlphaElectric::computeBetaQm( parameter_type const& mu )
{
    auto betaA = this->mdeim()->beta(mu);
    auto betaF = this->deim()->beta(mu);
    this->fillBetaQm(mu, betaA, betaF);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

void AlphaElectric::fillBetaQm( parameter_type const& mu, vectorN_type betaA, vectorN_type betaF )
{
    int M_A = this->mdeim()->size();
    for( int i = 0; i < M_A; ++i )
        M_betaAqm[0][i] = betaA(i);

    int M_F = this->deim()->size();
    for( int i = 0; i < M_F; ++i )
        M_betaFqm[0][0][i] = betaF(i);
}

AlphaElectric::beta_vector_type
AlphaElectric::computeBetaLinearDecompositionA( parameter_type const& mu, double time )
{
    beta_vector_type beta;
    return beta;
}

AlphaElectric::affine_decomposition_type
AlphaElectric::computeAffineDecomposition()
{
    return boost::make_tuple( this->M_Aqm, this->M_Fqm);
}

std::vector<std::vector<AlphaElectric::sparse_matrix_ptrtype> >
AlphaElectric::computeLinearDecompositionA()
{
    return this->M_linearAqm;
}

std::vector<std::vector<AlphaElectric::element_ptrtype> >
AlphaElectric::computeInitialGuessAffineDecomposition()
{
    return M_InitialGuess;
}

AlphaElectric::element_type
AlphaElectric::solve( parameter_type const& mu )
{
    auto solution = Xh->element();

    auto V = Xh->element();
    auto phiV = Xh->element();

    auto bc = M_modelProps->boundaryConditions();
    double sigmaMax = 0;
    for( auto const& material : M_materials )
    {
        auto sigma = material.second.getDouble("sigma");
        if( sigma > sigmaMax )
            sigmaMax = sigma;
    }

    /***************************** Electro *****************************/
    tic();
    auto aV = form2(_test=Xh, _trial=Xh);
    for( auto const& material : M_materialsWithGeo )
    {
        auto sigma = material.second.getDouble("sigma");

        auto alphaExpr = expr(this->alpha(mu, material.second));
        auto alphaPrimeExpr = expr(this->alphaPrime(mu, material.second));

        auto Jinv = mat<3,3>( cos(alphaExpr), -sin(alphaExpr), -alphaPrimeExpr*Py(),
                              sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*Px(),
                              cst(0.), cst(0.), cst(1.) );

        auto gradVJ = gradt(V)*Jinv;
        auto gradPhiVJ = grad(phiV)*Jinv;

        aV += integrate( markedelements(M_mesh, material.first),
                         sigma/sigmaMax*inner(gradVJ,gradPhiVJ) );
    }
    for( auto const& material : M_materialsWithoutGeo )
        aV += integrate( markedelements(M_mesh, material.first),
                         material.second.getDouble("sigma")/sigmaMax*inner(gradt(V),grad(phiV)) );

    auto fV = form1(_test=Xh);
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        if( isInMaterials(exAtM) )
        {
            auto sigma = M_materials[exAtM.material()].getDouble("sigma");
            aV += integrate( markedfaces(M_mesh, exAtM.marker()),
                             sigma/sigmaMax*(M_penalDir/hFace()*inner(idt(V),id(phiV))
                                             -inner(grad(V)*N(),id(phiV))
                                             -inner(grad(phiV)*N(),idt(V)) ) );

            auto dif = expr(exAtM.expression());
            fV += integrate( markedfaces(M_mesh, exAtM.marker()),
                             sigma/sigmaMax*dif*(M_penalDir/hFace()*id(phiV) - grad(phiV)*N()) );
        }
    }
    aV.solve(_rhs=fV, _solution=solution, _name="feV");
    toc("solve V", M_verbose > 0);

    return solution;
}

double AlphaElectric::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{
    auto mesh = Xh->mesh();
    double output=0;
    if ( output_index == 0 )
    {
        for ( int q = 0; q < Ql(0); q++ )
        {
            element_ptrtype eltF( new element_type( Xh ) );
            *eltF = *M_Fqm[output_index][q][0];
            output += M_betaFqm[output_index][q][0]*dot( *eltF, u );
            //output += M_betaFqm[output_index][q][m]*dot( M_Fqm[output_index][q][m], U );
        }
    }
    else
        throw std::logic_error( "[Heat2d::output] error with output_index : only 0 or 1 " );
    return output;
}

double AlphaElectric::sigma( std::string mat )
{
    return M_materials[mat].getDouble("sigma");
}

void AlphaElectric::computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu )
{
    auto V = this->solve(mu);
    auto Vh = j.functionSpace();
    for(auto const& mat : M_materials )
    {
        auto sigma = mat.second.getDouble("sigma");
        j += vf::project(Vh, markedelements(M_mesh,mat.first), cst(-1.)*sigma*trans(gradv(V)) );
    }
}

bool AlphaElectric::isInMaterials(ExpressionStringAtMarker const& ex) const
{
    return ex.hasMaterial() && M_materials.find(ex.material()) != M_materials.end();
}

FEELPP_CRB_PLUGIN( AlphaElectric, electric)
}
