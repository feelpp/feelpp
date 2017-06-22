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
#include "nlthermoelectric.hpp"

namespace Feel {

NLThermoelectric::NLThermoelectric()
    : super_type( "thermoelectric" )
{}

NLThermoelectric::NLThermoelectric( mesh_ptrtype mesh )
    : super_type( "thermoelectric" )
{
    this->M_mesh = mesh;
}

int NLThermoelectric::Qa()
{
    auto bc = M_modelProps->boundaryConditions();
    return 2 + bc["potential"]["Dirichlet"].size() + bc["temperature"]["Robin"].size();
}

int NLThermoelectric::mQA( int q )
{
    auto eimSigma = this->scalarContinuousEim()[0];
    // auto eimK = this->scalarContinuousEim()[1];
    auto bc = M_modelProps->boundaryConditions();

    if( q == 0 )
        return eimSigma->mMax();
    else if( q == 1 )
        return 1;//eimK->mMax();
    else if( q < 2 + bc["potential"]["Dirichlet"].size() )
        return eimSigma->mMax();
    else if( q < 2 + bc["potential"]["Dirichlet"].size() + bc["temperature"]["Robin"].size() )
        return 1;
    else
        return 0;
}

int NLThermoelectric::Nl()
{
    return 1;
}

int NLThermoelectric::Ql( int l)
{
    auto bc = M_modelProps->boundaryConditions();
    switch(l) {
    case 0:
        return 1 + bc["potential"]["Dirichlet"].size() + bc["temperature"]["Robin"].size();
    default:
        return 0;
    }
}

int NLThermoelectric::mLQF( int l, int q )
{
    switch( l )
    {
    case 0:
        return mCompliantQ(q);
    default:
        return 0;
    }
}

int NLThermoelectric::mCompliantQ(int q )
{
    auto eimSigma = this->scalarContinuousEim()[0];
    auto eimGradGrad = this->scalarDiscontinuousEim()[0];
    auto bc = M_modelProps->boundaryConditions();
    if( q == 0 )
        return eimGradGrad->mMax();
    else if( q < 1+bc["potential"]["Dirichlet"].size() )
        return eimSigma->mMax();
    else if( q < 1+bc["potential"]["Dirichlet"].size()+bc["temperature"]["Robin"].size() )
        return 1;
    else
        return 0;
}

void NLThermoelectric::resizeQm()
{
    M_Aqm.resize( Qa());
    M_betaAqm.resize( Qa() );
    for( int q = 0; q < Qa(); ++q )
    {
        M_Aqm[q].resize(mQA(q), backend()->newMatrix(Xh, Xh ) );
        M_betaAqm[q].resize(mQA(q));
    }

    M_Fqm.resize(Nl());
    M_betaFqm.resize(Nl());
    for( int l = 0; l < Nl(); ++l )
    {
        M_Fqm[l].resize(Ql(l));
        M_betaFqm[l].resize(Ql(l));
        for( int q = 0; q < Ql(l); ++q )
        {
            M_Fqm[l][q].resize(mLQF(l, q), backend()->newVector(Xh) );
            M_betaFqm[l][q].resize(mLQF(l, q) );
        }
    }

    M_InitialGuess.resize(1);
    M_InitialGuess[0].resize(1);
    M_InitialGuess[0][0] = Xh->elementPtr();
}

void NLThermoelectric::initModel()
{
    Feel::cout << "initModel" << std::endl;
    M_modelProps = boost::make_shared<prop_type>(Environment::expand( soption("thermoelectric.filename")));

    auto parameters = M_modelProps->parameters();
    Dmu->setDimension(parameters.size());
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

    if( !M_mesh )
        M_mesh = loadMesh( new mesh_type );
    this->setFunctionSpaces(functionspace_type::New( M_mesh ) );

    if( !pT )
        pT = element_ptrtype( new element_type( Xh ) );
    M_V = pT->template elementPtr<0>();
    M_T = pT->template elementPtr<1>();

    auto JspaceEim = J_space_type::New( M_mesh );
    auto TspaceEim = T_space_type::New( M_mesh );

    auto Pset = this->Dmu->sampling();
    int Ne = ioption(_name="thermoelectric.trainset-eim-size");
    std::vector<size_type> N(parameterSpace()->dimension(), Ne);

    std::string supersamplingname =(boost::format("DmuEim-Ne%1%-generated-by-master-proc") %Ne ).str();
    std::ifstream file ( supersamplingname );
    bool all_proc_same_sampling=true;
    if( ! file )
    {
        Pset->equidistributeProduct( N , all_proc_same_sampling , supersamplingname );
        Pset->writeOnFile( supersamplingname );
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(supersamplingname);
    }

    auto sigma = M_mu.parameterNamed("sigma")/( cst(1.)+M_mu.parameterNamed("alpha")*(idv(M_T)-cst(293.)));
    // auto k = sigma*M_mu.parameterNamed("L")*idv(M_T);
    auto gradgrad = sigma*inner(gradv(M_V));

    tic();
    auto eim_sigma = eim( _model=boost::dynamic_pointer_cast<NLThermoelectric>(this->shared_from_this() ),
                          _element=*M_T,
                          _parameter=M_mu,
                          _expr=sigma,
                          _space=TspaceEim,
                          _name="eim_sigma",
                          _sampling=Pset );
    this->addEim( eim_sigma );
    toc("eim_sigma");

    // tic();
    // auto eim_k = eim( _model=boost::dynamic_pointer_cast<NLThermoelectric>(this->shared_from_this() ),
    //                   _element=*M_T,
    //                   _parameter=M_mu,
    //                   _expr=k,
    //                   _space=TspaceEim,
    //                   _name="eim_k",
    //                   _sampling=Pset );
    // this->addEim( eim_k );
    // toc("eim_k");

    tic();
    auto eim_gradgrad = eim( _model=boost::dynamic_pointer_cast<NLThermoelectric>(this->shared_from_this() ),
                             _element=*M_T,
                             _element2=*M_V,
                             _parameter=M_mu,
                             _expr=gradgrad,
                             _space=JspaceEim,
                             _name="eim_gradgrad",
                             _sampling=Pset );
    this->addEimDiscontinuous( eim_gradgrad );
    toc("eim_grad_grad");

    this->resizeQm();
    this->decomposition();
}

void NLThermoelectric::decomposition()
{
    auto UT = Xh->element();
    auto VP = Xh->element();
    auto V = UT.template element<0>();
    auto phiV = VP.template element<0>();
    auto T = UT.template element<1>();
    auto phiT = VP.template element<1>();

    auto gamma = doption("thermoelectric.gamma");

    auto eimSigma = this->scalarContinuousEim()[0];
    // auto eimK = this->scalarContinuousEim()[1];
    auto eimGradGrad = this->scalarDiscontinuousEim()[0];
    auto bc = M_modelProps->boundaryConditions();

    /************** Right hand side **************/
    // electro
    for( int m = 0; m < eimSigma->mMax(); ++m )
    {
        auto a0 = form2(_test=Xh, _trial=Xh);
        a0 = integrate( elements(M_mesh),
                        idv(eimSigma->q(m))*inner(gradt(V),grad(phiV)) );
        M_Aqm[0][m] = a0.matrixPtr();
    }

    // thermo
    // for( int m = 0; m < eimK->mMax(); ++m )
    // {
        auto a1 = form2(_test=Xh, _trial=Xh);
        a1 = integrate( elements(M_mesh),
                        inner(gradt(T), grad(phiT)) );
        M_Aqm[1][0] = a1.matrixPtr();
        //                 idv(eimK->q(m))*inner(gradt(T), grad(phiT)) );
        // M_Aqm[1][m] = a1.matrixPtr();
    // }

    int idx = 2;
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        for( int m = 0; m < eimSigma->mMax(); ++m )
        {
            auto aVD = form2(_test=Xh, _trial=Xh);
            aVD = integrate( markedfaces(M_mesh, exAtM.marker() ),
                             idv(eimSigma->q(m))*(gamma/hFace()*inner(idt(V),id(phiV))
                                                  -inner(gradt(V)*N(),id(phiV))
                                                  -inner(grad(phiV)*N(),idt(V)) ) );
            M_Aqm[idx][m] = aVD.matrixPtr();
        }
        idx++;
    }

    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto aTR = form2(_test=Xh, _trial=Xh);
        aTR = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         inner(id(T), idt(phiT)) );
        M_Aqm[idx++][0] = aTR.matrixPtr();
    }

    /************** Left hand side **************/
    for( int m = 0; m < eimGradGrad->mMax(); ++m )
    {
        auto f0 = form1(_test=Xh);
        f0 = integrate(elements(M_mesh),
                       inner(id(phiT), idv(eimGradGrad->q(m))) );
        M_Fqm[0][0][m] = f0.vectorPtr();
    }

    idx = 1;
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        for( int m = 0; m < eimSigma->mMax(); ++m )
        {
            auto fVD = form1(_test=Xh);
            fVD = integrate( markedfaces(M_mesh, exAtM.marker() ),
                             idv(eimSigma->q(m))*(gamma/hFace()*id(phiV) -  grad(phiV)*N() ) );
            M_Fqm[0][idx][m] = fVD.vectorPtr();
        }
        idx++;
    }

    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto fTR = form1(_test=Xh);
        fTR = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         idt(phiT) );
        M_Fqm[0][idx++][0] = fTR.vectorPtr();
    }

    // Energy matrix
    auto m = form2(_test=Xh, _trial=Xh);
    m = integrate( elements(M_mesh),
                   inner(gradt(V),grad(phiV)) + inner(grad(T), gradt(phiT)) );
    M_energy_matrix = m.matrixPtr();
}

NLThermoelectric::beta_vector_type
NLThermoelectric::computeBetaInitialGuess( parameter_type const& mu )
{
    M_betaInitialGuess.resize( 1 );
    M_betaInitialGuess[0].resize( 1 );
    M_betaInitialGuess[0][0] = 1;
    return this->M_betaInitialGuess;
}

NLThermoelectric::beta_type
NLThermoelectric::computeBetaQm( element_type const& T, parameter_type const& mu )
{
    auto eimSigma = this->scalarContinuousEim()[0];
    auto betaEimSigma = eimSigma->beta( mu, T );
    // auto eimK = this->scalarContinuousEim()[1];
    // auto betaEimK = eimK->beta( mu, T );
    auto betaEimK = betaEimSigma;
    auto eimGradGrad = this->scalarDiscontinuousEim()[0];
    auto betaEimGradGrad = eimGradGrad->beta( mu, T );
    this->fillBetaQm(mu, betaEimSigma, betaEimK,  betaEimGradGrad);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

NLThermoelectric::beta_type
NLThermoelectric::computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
{
    auto eimSigma = this->scalarContinuousEim()[0];
    auto betaEimSigma = eimSigma->beta( mu );
    // auto eimK = this->scalarContinuousEim()[1];
    // auto betaEimK = eimK->beta( mu );
    auto betaEimK = betaEimSigma;
    auto eimGradGrad = this->scalarDiscontinuousEim()[0];
    auto betaEimGradGrad = eimGradGrad->beta( mu );
    this->fillBetaQm(mu, betaEimSigma, betaEimK,  betaEimGradGrad);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

NLThermoelectric::beta_type
NLThermoelectric::computeBetaQm( parameter_type const& mu )
{
    auto eimSigma = this->scalarContinuousEim()[0];
    auto betaEimSigma = eimSigma->beta( mu );
    auto eimK = this->scalarContinuousEim()[1];
    // auto betaEimK = eimK->beta( mu );
    auto betaEimK = betaEimSigma;
    auto eimGradGrad = this->scalarDiscontinuousEim()[0];
    auto betaEimGradGrad = eimGradGrad->beta( mu );
    this->fillBetaQm(mu, betaEimSigma, betaEimK,  betaEimGradGrad);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

void NLThermoelectric::fillBetaQm( parameter_type const& mu, vectorN_type betaEimSigma, vectorN_type betaEimK, vectorN_type betaEimGradGrad )
{
    auto bc = M_modelProps->boundaryConditions();

    for( int m = 0; m < betaEimSigma.size(); ++m )
        M_betaAqm[0][m] = betaEimSigma(m);
    // for( int m = 0; m < betaEimK.size(); ++m )
    //     M_betaAqm[1][m] = betaEimK(m);
    M_betaAqm[1][0] = mu.parameterNamed("k");

    int idx = 2;
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        for( int m = 0; m < betaEimSigma.size(); ++m )
            M_betaAqm[idx][m] = betaEimSigma(m);
        idx++;
    }
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression1());
        auto symb = e.expression().symbols();
        if( symb.size() > 0 )
            for( auto const& param : symb )
                e.setParameterValues( { param.get_name(), mu.parameterNamed(param.get_name()) } );
        M_betaAqm[idx++][0] = e.evaluate();
    }

    for( int m = 0; m < betaEimGradGrad.size(); ++m )
        M_betaFqm[0][0][m] = betaEimGradGrad(m);
    idx = 1;
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto e = expr(exAtM.expression());
        auto symb = e.expression().symbols();
        if( symb.size() > 0 )
            for( auto const& param : symb )
                e.setParameterValues( { param.get_name(), mu.parameterNamed(param.get_name()) } );
        for( int m = 0; m < betaEimSigma.size(); ++m )
            M_betaFqm[0][idx][m] = betaEimSigma(m)*e.evaluate();
        idx++;
    }
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression2());
        auto symb = e.expression().symbols();
        if( symb.size() > 0 )
            for( auto const& param : symb )
                e.setParameterValues( { param.get_name(), mu.parameterNamed(param.get_name()) } );
        M_betaFqm[0][idx++][0] = e.evaluate();
    }
}

NLThermoelectric::beta_vector_type
NLThermoelectric::computeBetaLinearDecompositionA( parameter_type const& mu, double time )
{
    beta_vector_type beta;
    return beta;
}

NLThermoelectric::affine_decomposition_type
NLThermoelectric::computeAffineDecomposition()
{
    return boost::make_tuple( this->M_Aqm, this->M_Fqm);
}

std::vector<std::vector<NLThermoelectric::sparse_matrix_ptrtype> >
NLThermoelectric::computeLinearDecompositionA()
{
    return this->M_linearAqm;
}

std::vector<std::vector<NLThermoelectric::element_ptrtype> >
NLThermoelectric::computeInitialGuessAffineDecomposition()
{
    return M_InitialGuess;
}

NLThermoelectric::element_type
NLThermoelectric::solve( parameter_type const& mu )
{
    // Feel::cout << "solve for parameter:" << std::endl << mu << std::endl;
    auto Vh = Xh->template functionSpace<0>();
    auto Th = Xh->template functionSpace<1>();
    auto V = Vh->element();
    auto phiV = Vh->element();
    auto T = Th->element();
    auto phiT = Th->element();
    auto oldV = Vh->element();
    auto oldT = Th->element(cst(293.0));

    auto bc = M_modelProps->boundaryConditions();
    auto gamma = doption("thermoelectric.gamma");
    auto sigma0 = mu.parameterNamed("sigma");
    auto alpha = mu.parameterNamed("alpha");
    // auto L = mu.parameterNamed("L");
    auto k = mu.parameterNamed("k");
    auto current = mu.parameterNamed("current");
    auto h = mu.parameterNamed("h");
    auto Tw = mu.parameterNamed("Tw");

    auto maxIt = ioption("thermoelectric.maxit");
    auto tol = doption("thermoelectric.tol");
    double errV = 0, errT = 0;
    int i = 0;

    /***************************** Electro *****************************/
    tic();
        auto e = exporter(M_mesh);
    do
    {
        auto sigma = cst(sigma0)/(1+alpha*(idv(oldT)-cst(293.0)));
        // auto k = sigma*L*idv(oldT);

        auto a = form2(_test=Vh, _trial=Vh);
        // V
        a = integrate( elements(M_mesh),
                       sigma*inner(gradt(V),grad(phiV)) );

        // V Dirichlet condition
        for( auto const& exAtM : bc["potential"]["Dirichlet"] )
        {
            a += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            sigma*(gamma/hFace()*inner(idt(V),id(phiV))
                                   -inner(gradt(V)*N(),id(phiV))
                                   -inner(grad(phiV)*N(),idt(V)) ) );
        }

        auto f = form1(_test=Vh);
        // V Dirichlet condition
        for( auto const& exAtM : bc["potential"]["Dirichlet"] )
        {
            auto e = expr(exAtM.expression());
            auto symb = e.expression().symbols();
            if( symb.size() > 0 )
                for( auto const& param : symb )
                    e.setParameterValues( { param.get_name(), mu.parameterNamed(param.get_name()) } );

            f += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            sigma*e*(gamma/hFace()*id(phiV) -  grad(phiV)*N()) );
        }
        a.solve( _solution=V, _rhs=f, _name="mono" );

        auto aT = form2(_test=Th, _trial=Th);
        // T
        aT += integrate( elements(M_mesh),
                         k*inner(gradt(T), grad(phiT)) );

        // T Robin condition
        for( auto const& exAtM : bc["temperature"]["Robin"] )
        {
            auto e = expr(exAtM.expression1());
            auto symb = e.expression().symbols();
            if( symb.size() > 0 )
                for( auto const& param : symb )
                    e.setParameterValues( { param.get_name(), mu.parameterNamed(param.get_name()) } );
            aT += integrate( markedfaces(M_mesh, exAtM.marker() ),
                             e*inner(idt(T), id(phiT)) );
        }

        auto fT = form1(_test=Th);
        // T right hand side
        fT = integrate(elements(M_mesh),
                       id(phiT)*sigma*gradv(V)*trans(gradv(V)) );

        // T Robin condition
        for( auto const& exAtM : bc["temperature"]["Robin"] )
        {
            auto e = expr(exAtM.expression2());
            auto symb = e.expression().symbols();
            if( symb.size() > 0 )
                for( auto const& param : symb )
                    e.setParameterValues( { param.get_name(), mu.parameterNamed(param.get_name()) } );
            fT += integrate( markedfaces(M_mesh, exAtM.marker() ),
                             e*id(phiT) );
        }

        aT.solve( _solution=T, _rhs=fT, _name="mono" );

        double errOldV = normL2(elements(M_mesh), idv(oldV));
        double errOldT = normL2(elements(M_mesh), idv(oldT));
        errV = normL2( elements(M_mesh), idv(V) - idv(oldV))/errOldV;
        errT = normL2( elements(M_mesh), idv(T) - idv(oldT))/errOldT;
        Feel::cout << "err(V) = " << errV << " err(T) = " << errT << std::endl;

        oldV = V;
        oldT = T;

        auto s = Vh->element(sigma);
        e->step(i)->add( "V", V);
        e->step(i)->add( "T", T);
        e->step(i)->add("sigma", s);
        e->save();
    } while( ++i < maxIt && ( errV > tol || errT > tol ) );

    auto solution = Xh->element();
    solution.template element<0>() = V;
    solution.template element<1>() = T;
    toc("picard");

    return solution;
}

double NLThermoelectric::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
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

int NLThermoelectric::mMaxSigma()
{
    auto eimSigma = this->scalarContinuousEim()[0];
    return eimSigma->mMax();
}

auto NLThermoelectric::eimSigmaQ(int m)
{
    auto eimSigma = this->scalarContinuousEim()[0];
    return eimSigma->q()[m];
}

NLThermoelectric::vectorN_type NLThermoelectric::eimSigmaBeta( parameter_type const& mu )
{
    auto eimSigma = this->scalarContinuousEim()[0];
    return eimSigma->beta(mu);
}

template<typename vec_space_type>
typename vec_space_type::element_type
NLThermoelectric::computeTruthCurrentDensity( parameter_type const& mu )
{
    auto VT = this->solve(mu);
    auto V = VT.template element<0>();
    auto T = VT.template element<1>();
    auto T0 = 293.0;
    auto sigma = mu.parameterNamed("sigma")/(cst(1.) + mu.parameterNamed("alpha")*(idv(T) - T0) );
    auto Vh = vec_space_type::New(M_mesh);
    auto j = Vh->element();
    j = vf::project(Vh, elements(M_mesh), cst(-1.)*sigma*trans(gradv(V)) );
    return j;
}

FEELPP_CRB_PLUGIN( NLThermoelectric, "nlthermoelectric")
}

