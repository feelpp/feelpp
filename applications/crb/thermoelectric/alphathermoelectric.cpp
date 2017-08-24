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

#include <feel/feelcrb/crbplugin.hpp>
#include "alphathermoelectric.hpp"

#include <boost/math/special_functions/binomial.hpp>

namespace Feel {

AlphaThermoelectric::AlphaThermoelectric()
    : super_type( "thermoelectric" )
{}

AlphaThermoelectric::AlphaThermoelectric( mesh_ptrtype mesh )
    : super_type( "thermoelectric" )
{
    this->M_mesh = mesh;
}

int AlphaThermoelectric::Qa()
{
    return 1;
}

int AlphaThermoelectric::mQA( int q )
{
    return this->M_mdeim->size();
}

int AlphaThermoelectric::Nl()
{
    return 1;
}

int AlphaThermoelectric::Ql( int l)
{
    return 1;
}

int AlphaThermoelectric::mLQF( int l, int q )
{
    switch( l )
    {
    case 0:
        return mCompliantQ(q);
    default:
        return 0;
    }
}

int AlphaThermoelectric::mCompliantQ(int q )
{
    return this->M_deim->size();
}

void AlphaThermoelectric::resizeQm( bool resizeMat )
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

std::string AlphaThermoelectric::alpha( parameter_type const& mu )
{
    using boost::math::binomial_coefficient;

    auto parameters = M_modelProps->parameters();
    double zb = parameters["zb"].value();
    double zt = parameters["zt"].value();

    std::vector<double> p( 2, 0);
    for( int i = 0; i < mu.size(); ++i )
        p.insert( p.begin() +i+1, mu(i) );
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

std::string AlphaThermoelectric::alphaPrime( parameter_type const& mu )
{
    using boost::math::binomial_coefficient;

    auto parameters = M_modelProps->parameters();
    double zb = parameters["zb"].value();
    double zt = parameters["zt"].value();

    std::vector<double> p( 2, 0);
    for( int i = 0; i < mu.size(); ++i )
        p.insert( p.begin() +i+1, mu(i) );
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

AlphaThermoelectric::sparse_matrix_ptrtype
AlphaThermoelectric::assembleForMDEIM( parameter_type const& mu )
{
    auto alphaExpr = expr(this->alpha(mu));
    auto alphaPrimeExpr = expr(this->alphaPrime(mu));

    auto Jinv = mat<3,3>( cos(alphaExpr), -sin(alphaExpr), -alphaPrimeExpr*Py(),
                          sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*Px(),
                          cst(0.), cst(0.), cst(1.) );

    auto Vh = Xh->template functionSpace<0>();
    auto Th = Xh->template functionSpace<1>();
    auto V = Vh->element();
    auto phiV = Vh->element();
    auto T = Th->element();
    auto phiT = Th->element();

    auto gradVJ = gradt(V)*Jinv;
    auto gradPhiVJ = grad(phiV)*Jinv;
    auto gradTJ = gradt(T)*Jinv;
    auto gradPhiTJ = grad(phiT)*Jinv;

    auto parameters = M_modelProps->parameters();
    double gamma = doption("thermoelectric.gamma");
    double sigma = parameters["sigma"].value();
    double k = parameters["k"].value();
    double dif = parameters["dif"].value();
    double h = parameters["h"].value();
    double Tw = parameters["Tw"].value();

    // auto bc = M_modelProps->boundaryConditions();

    auto a = form2(_test=Xh, _trial=Xh);
    a = integrate( elements(M_mesh), sigma*inner(gradVJ,gradPhiVJ) );
    a += integrate( markedfaces(M_mesh, "base"),
                    sigma*(gamma/hFace()*inner(idt(V),id(phiV))
                           -inner(gradVJ*N(),id(phiV))
                           -inner(gradPhiVJ*N(),idt(V)) ) );
    a += integrate( markedfaces(M_mesh, "top"),
                    sigma*(gamma/hFace()*inner(idt(V),id(phiV))
                           -inner(gradVJ*N(),id(phiV))
                           -inner(gradPhiVJ*N(),idt(V)) ) );
    a += integrate( elements(M_mesh), k*inner(gradTJ,gradPhiTJ) );
    a += integrate( markedfaces(M_mesh, {"ext","int"}),
                    h*inner(idt(T),id(phiT)) );

    auto am = a.matrixPtr();
    am->close();

    return am;
}

AlphaThermoelectric::vector_ptrtype
AlphaThermoelectric::assembleForDEIMnl( parameter_type const& mu )
{
    auto U = this->solve(mu);
    return this->assembleForDEIMnl(mu, U);
}

AlphaThermoelectric::vector_ptrtype
AlphaThermoelectric::assembleForDEIMnl( parameter_type const& mu, element_type const& U )
{
    auto alphaExpr = expr(this->alpha(mu));
    auto alphaPrimeExpr = expr(this->alphaPrime(mu));

    auto Jinv = mat<3,3>( cos(alphaExpr), -sin(alphaExpr), -alphaPrimeExpr*Py(),
                          sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*Px(),
                          cst(0.), cst(0.), cst(1.) );

    auto Vh = Xh->template functionSpace<0>();
    auto Th = Xh->template functionSpace<1>();
    auto phiV = Vh->element();
    auto phiT = Th->element();
    auto V = U.template element<0>();

    auto gradPhiVJ = grad(phiV)*Jinv;
    auto gradPhiTJ = grad(phiT)*Jinv;
    auto gradvVJ = gradv(V)*Jinv;

    auto parameters = M_modelProps->parameters();
    double gamma = doption("thermoelectric.gamma");
    double sigma = parameters["sigma"].value();
    double k = parameters["k"].value();
    double dif = parameters["dif"].value();
    double h = parameters["h"].value();
    double Tw = parameters["Tw"].value();

    // auto bc = M_modelProps->boundaryConditions();

    auto f = form1(_test=Xh);
    f = integrate( markedfaces(M_mesh, "top"),
                   sigma*dif*(gamma/hFace()*id(phiV) - gradPhiVJ*N()) );
    f += integrate(elements(M_mesh), sigma*inner(gradvVJ,gradvVJ)*id(phiT) );
    f += integrate(markedfaces(M_mesh, {"ext","int"}),
                   h*Tw*id(phiT) );

    auto fv = f.vectorPtr();
    fv->close();

    return fv;
}

void AlphaThermoelectric::initModel()
{
    Feel::cout << "initModel" << std::endl;
    std::string propertyPath = Environment::expand( soption("thermoelectric.filename"));
    M_modelProps = boost::make_shared<prop_type>(propertyPath);
    this->addModelFile("property-file", propertyPath);

    auto parameters = M_modelProps->parameters();
    Feel::cout << "Using parameters:" << std::endl
               << "sigma: " << parameters["sigma"].value() << std::endl
               << "k    : " << parameters["k"].value() << std::endl
               << "dif  : " << parameters["dif"].value() << std::endl
               << "h    : " << parameters["h"].value() << std::endl
               << "Tw   : " << parameters["Tw"].value() << std::endl;

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

    if( !M_mesh )
        M_mesh = loadMesh( new mesh_type );
    this->setFunctionSpaces(functionspace_type::New( M_mesh ) );

    Feel::cout << "Potential nDof  : " << Xh->template functionSpace<0>()->nDof() << std::endl
               << "Temperature nDof: " << Xh->template functionSpace<1>()->nDof() << std::endl;

    if( !pT )
        pT = element_ptrtype( new element_type( Xh ) );
    M_V = pT->template elementPtr<0>();
    M_T = pT->template elementPtr<1>();

    auto PsetV = this->Dmu->sampling();
    int N = ioption(_name="thermoelectric.trainset-deim-size");

    std::string supersamplingname =(boost::format("DmuDEim-Ne%1%-generated-by-master-proc") % N ).str();
    std::ifstream file ( supersamplingname );
    bool all_proc_same_sampling=true;
    if( ! file )
    {
        PsetV->randomize( N , all_proc_same_sampling , supersamplingname );
        PsetV->writeOnFile( supersamplingname );
    }
    else
    {
        PsetV->clear();
        PsetV->readFromFile(supersamplingname);
    }

    this->M_deim = boost::make_shared<deim_type>(boost::dynamic_pointer_cast<self_type>(this->shared_from_this()), PsetV, "vec");
    this->M_deim->assemble = boost::bind( &AlphaThermoelectric::assembleForDEIMnl, boost::ref(*this), _1);
    this->M_deim->run();

    auto PsetM = this->Dmu->sampling();
    N = ioption(_name="thermoelectric.trainset-mdeim-size");

    supersamplingname =(boost::format("DmuMDEim-Ne%1%-generated-by-master-proc") % N ).str();
    std::ifstream fileM ( supersamplingname );
    if( ! fileM )
    {
        PsetM->randomize( N , all_proc_same_sampling , supersamplingname );
        PsetM->writeOnFile( supersamplingname );
    }
    else
    {
        PsetM->clear();
        PsetM->readFromFile(supersamplingname);
    }

    this->M_mdeim = boost::make_shared<mdeim_type>(boost::dynamic_pointer_cast<self_type>(this->shared_from_this()), PsetM, "mat");
    this->M_mdeim->assemble = boost::bind( &AlphaThermoelectric::assembleForMDEIM, boost::ref(*this), _1);
    this->M_mdeim->run();

    this->resizeQm();
    this->decomposition();

    auto s = this->solve(M_mu);
    auto VFE = s.template element<0>();
    auto TFE = s.template element<1>();
    auto A = this->assembleForMDEIM(M_mu);
    auto F = this->assembleForDEIMnl(M_mu,s);

    auto betaA = M_mdeim->beta(A);
    auto betaF = M_deim->beta(F);
    auto qa = M_mdeim->q();
    auto qf = M_deim->q();
    auto Am = backend()->newMatrix(Xh,Xh);
    Am->zero();
    for( int m = 0; m < M_mdeim->size(); ++m )
    {
        Am->addMatrix( betaA(m), qa[m]);
        Feel::cout << "betaA(" << m << ") = " << betaA(m) << std::endl;
    }
    auto Fm = backend()->newVector(Xh);
    Fm->zero();
    for( int m = 0; m < M_deim->size(); ++m )
    {
        Fm->add( betaF(m), qf[m]);
        Feel::cout << "betaF(" << m << ") = " << betaF(m) << std::endl;
    }

    auto VT = Xh->element();
    backend()->solve( _matrix=Am, _rhs=Fm, _solution=VT);
    auto VRB = VT.template element<0>();
    auto TRB = VT.template element<1>();

    auto e = exporter(M_mesh);
    e->add("VFE", VFE);
    e->add("TFE", TFE);
    e->add("VEIM", VRB);
    e->add("TEIM", TRB);
    e->save();

    auto errV = normL2(elements(M_mesh), idv(VFE)-idv(VRB) );
    auto normV = normL2(elements(M_mesh), idv(VFE) );
    auto errT = normL2(elements(M_mesh), idv(TFE)-idv(TRB) );
    auto normT = normL2(elements(M_mesh), idv(TFE) );
    Feel::cout << "V: err = " << errV << " relative err = " << errV/normV << std::endl
               << "T: err = " << errT << " relative err = " << errT/normT << std::endl;

    Am->addMatrix(-1., A);
    Fm->add(-1., F);
    auto errA = Am->linftyNorm();
    auto normA = A->linftyNorm();
    auto errF = Fm->linftyNorm();
    auto normF = F->linftyNorm();
    Feel::cout << "A: err = " << errA << " relative err = " << errA/normA << std::endl
               << "F: err = " << errF << " relative err = " << errF/normF << std::endl;
}

void AlphaThermoelectric::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    std::string propertyPath;
    if( this->hasModelFile("property-file") )
        propertyPath = this->additionalModelFiles().find("property-file")->second;
    else
        Feel::cerr << "Warning!! the database does not contain the property file! Expect bugs!"
                   << std::endl;
    M_modelProps = boost::make_shared<prop_type>(propertyPath);

    auto parameters = M_modelProps->parameters();
    Dmu = parameterspace_type::New( parameters.size(), Environment::worldComm() );

    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    int i = 0;
    for( auto const& parameterPair : parameters )
    {
        mu_min(i) = parameterPair.second.min();
        mu_max(i) = parameterPair.second.max();
        Dmu->setParameterName(i++, parameterPair.first );
    }
    Dmu->setMin(mu_min);
    Dmu->setMax(mu_max);
    M_mu = Dmu->element();

    boost::shared_ptr<J_space_type> JspaceEim;
    if ( !pT )
        pT.reset( new element_type );
    M_V = pT->template elementPtr<0>();
    M_T = pT->template elementPtr<1>();

    // auto const& ptreeEim = ptree.get_child( "eim" );
    // auto const& ptreeEimJoule = ptreeEim.get_child( "eim_joule" );
    // std::string dbnameEimJoule = ptreeEimJoule.template get<std::string>( "database-filename" );

    this->resizeQm( false );
}

void AlphaThermoelectric::decomposition()
{
    int M_A = this->M_mdeim->size();
    auto qa = this->M_mdeim->q();
    for( int i = 0; i < M_A; ++i )
        M_Aqm[0][i] = qa[i];

    int M_F = this->M_deim->size();
    auto qf = this->M_deim->q();
    for( int i = 0; i < M_F; ++i )
        M_Fqm[0][0][i] = qf[i];

    // Energy matrix
    auto U = Xh->element();
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto t = U.template element<1>();
    auto p = U.template element<1>();
    auto m = form2(_test=Xh, _trial=Xh);
    m = integrate( elements(M_mesh),
                   inner(gradt(u),grad(v)) + inner(grad(t), gradt(p)) );
    M_energy_matrix = m.matrixPtr();
}

AlphaThermoelectric::beta_vector_type
AlphaThermoelectric::computeBetaInitialGuess( parameter_type const& mu )
{
    M_betaInitialGuess.resize( 1 );
    M_betaInitialGuess[0].resize( 1 );
    M_betaInitialGuess[0][0] = 1;
    return this->M_betaInitialGuess;
}

AlphaThermoelectric::beta_type
AlphaThermoelectric::computeBetaQm( element_type const& u, parameter_type const& mu )
{
    auto betaA = this->M_mdeim->beta(mu);
    auto T = assembleForDEIMnl(mu, u);
    auto betaF = this->M_deim->beta(T);
    this->fillBetaQm(mu, betaA, betaF);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

AlphaThermoelectric::beta_type
AlphaThermoelectric::computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
{
    auto u = this->solve(mu);
    return this->computeBetaQm( u, mu);
    // auto betaA = this->M_mdeim->beta(mu, urb);
    // auto betaF = this->M_deim->beta(mu, urb);
    // this->fillBetaQm(mu, betaA, betaF);
    // return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

AlphaThermoelectric::beta_type
AlphaThermoelectric::computeBetaQm( parameter_type const& mu )
{
    auto u = this->solve(mu);
    return this->computeBetaQm( u, mu);
}

void AlphaThermoelectric::fillBetaQm( parameter_type const& mu, vectorN_type betaA, vectorN_type betaF )
{
    int M_A = this->M_mdeim->size();
    for( int i = 0; i < M_A; ++i )
        M_betaAqm[0][i] = betaA(i);

    int M_F = this->M_deim->size();
    for( int i = 0; i < M_F; ++i )
        M_betaFqm[0][0][i] = betaF(i);
}

AlphaThermoelectric::beta_vector_type
AlphaThermoelectric::computeBetaLinearDecompositionA( parameter_type const& mu, double time )
{
    beta_vector_type beta;
    return beta;
}

AlphaThermoelectric::affine_decomposition_type
AlphaThermoelectric::computeAffineDecomposition()
{
    return boost::make_tuple( this->M_Aqm, this->M_Fqm);
}

std::vector<std::vector<AlphaThermoelectric::sparse_matrix_ptrtype> >
AlphaThermoelectric::computeLinearDecompositionA()
{
    return this->M_linearAqm;
}

std::vector<std::vector<AlphaThermoelectric::element_ptrtype> >
AlphaThermoelectric::computeInitialGuessAffineDecomposition()
{
    return M_InitialGuess;
}

AlphaThermoelectric::element_type
AlphaThermoelectric::solve( parameter_type const& mu )
{
    Feel::cout << "solve for parameter:" << std::endl << mu << std::endl;

    auto Vh = Xh->template functionSpace<0>();
    auto Th = Xh->template functionSpace<1>();
    auto V = Vh->element();
    auto phiV = Vh->element();
    auto T = Th->element();
    auto phiT = Th->element();

    auto parameters = M_modelProps->parameters();
    double gamma = doption("thermoelectric.gamma");
    double sigma = parameters["sigma"].value();
    double k = parameters["k"].value();
    double dif = parameters["dif"].value();
    double h = parameters["h"].value();
    double Tw = parameters["Tw"].value();

    auto alphaExpr = expr(this->alpha(mu));
    auto alphaPrimeExpr = expr(this->alphaPrime(mu));

    auto Jinv = mat<3,3>( cos(alphaExpr), -sin(alphaExpr), -alphaPrimeExpr*Py(),
                          sin(alphaExpr), cos(alphaExpr), alphaPrimeExpr*Px(),
                          cst(0.), cst(0.), cst(1.) );

    auto gradVJ = gradt(V)*Jinv;
    auto gradPhiVJ = grad(phiV)*Jinv;
    auto gradTJ = gradt(T)*Jinv;
    auto gradPhiTJ = grad(phiT)*Jinv;

    // auto bc = M_modelProps->boundaryConditions();

    /***************************** Electro *****************************/
    tic();
    auto aV = form2(_test=Vh, _trial=Vh);
    aV = integrate( elements(M_mesh), sigma*inner(gradVJ,gradPhiVJ) );
    aV += integrate( markedfaces(M_mesh, "base"),
                     sigma*(gamma/hFace()*inner(idt(V),id(phiV))
                            -inner(gradVJ*N(),id(phiV))
                            -inner(gradPhiVJ*N(),idt(V)) ) );
    aV += integrate( markedfaces(M_mesh, "top"),
                     sigma*(gamma/hFace()*inner(idt(V),id(phiV))
                            -inner(gradVJ*N(),id(phiV))
                            -inner(gradPhiVJ*N(),idt(V)) ) );

    auto fV = form1(_test=Vh);
    fV = integrate( markedfaces(M_mesh, "top"),
                    sigma*dif*(gamma/hFace()*id(phiV) - gradPhiVJ*N()) );
    aV.solve(_rhs=fV, _solution=M_V, _name="feV");
    toc("solve V");

    /**************************** Thermo *****************************/
    auto gradvVJ = gradv(M_V)*Jinv;
    tic();
    auto aT = form2(_test=Th, _trial=Th);
    aT = integrate( elements(M_mesh), k*inner(gradTJ,gradPhiTJ) );
    aT += integrate( markedfaces(M_mesh, {"ext","int"}),
                     h*inner(idt(T),id(phiT)) );
    auto fT = form1(_test=Th);
    fT = integrate(elements(M_mesh), sigma*inner(gradvVJ,gradvVJ)*id(phiT) );
    fT += integrate(markedfaces(M_mesh, {"ext","int"}),
                    h*Tw*id(phiT) );
    aT.solve(_rhs=fT, _solution=M_T, _name="feT");
    toc("solve T");

    auto solution = Xh->element();
    solution.template element<0>() = *M_V;
    solution.template element<1>() = *M_T;

    // auto e = exporter(M_mesh);
    // e->add("sol", solution);
    // e->save();

    return solution;
}

double AlphaThermoelectric::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
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

double AlphaThermoelectric::sigma()
{
    auto parameters = M_modelProps->parameters();
    return parameters["sigma"].value();
}

void AlphaThermoelectric::computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu )
{
    auto VT = this->solve(mu);
    auto V = VT.template element<0>();
    auto parameters = M_modelProps->parameters();
    auto sigma = parameters["sigma"].value();
    auto Vh = j.functionSpace();
    j = vf::project(Vh, elements(M_mesh), cst(-1.)*sigma*trans(gradv(V)) );
}

FEELPP_CRB_PLUGIN( AlphaThermoelectric, "thermoelectric")
}
