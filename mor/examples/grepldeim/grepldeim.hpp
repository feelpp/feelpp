/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): JB Wahl <wahl@math.unistra.fr>
  Date: 2017-04-19

  Copyright (C) 2011-2014 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file benchmarkgrepl-nonlinear-deim.hpp
   \author JB Wahl <wahl@math.unistra.fr>
   date 2017-04-19
*/

#ifndef FEELPP_GREPLDEIM_HPP
#define FEELPP_GREPLDEIM_HPP 1


#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelmor/modelcrbbase.hpp>
#include <feel/feelmor/deim.hpp>

using namespace Feel;

class ParameterDefinition
{
public :
    static const uint16_type ParameterSpaceDimension = 2;
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
};

template<int Order, int Dim>
class FunctionSpaceDefinition
{
public :
    //static const uint16_type Order = 1;
    typedef double value_type;

    /*mesh*/
    typedef Simplex<Dim,1> entity_type; /*dim,order*/
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;
    typedef bases<Lagrange<Order + 2, Scalar> > basis_type_eimg;
    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef FunctionSpace<mesh_type, basis_type_eimg, value_type> space_type_eimg;

    typedef typename space_type::element_type element_type;

};


template<int Order, int Dim>
class GreplDEIM :
    public ModelCrbBase< ParameterDefinition,
                         FunctionSpaceDefinition<Order,Dim>,
                         NonLinear >
{
public :
    typedef ModelCrbBase<ParameterDefinition,
                         FunctionSpaceDefinition<Order,Dim>,
                         NonLinear> super_type;
    typedef GreplDEIM<Order,Dim> self_type;

    typedef double value_type;

    typedef typename FunctionSpaceDefinition<Order,Dim>::space_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::mesh_ptrtype mesh_ptrtype;
    typedef typename super_type::element_ptrtype element_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::parameter_type parameter_type;

    typedef typename super_type::beta_vector_type beta_vector_type;
    typedef typename super_type::affine_decomposition_type affine_decomposition_type;
    typedef typename super_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super_type::vectorN_type vectorN_type;
    typedef typename super_type::monolithic_type monolithic_type;
    typedef typename super_type::vector_ptrtype vector_ptrtype;
    typedef std::vector< std::vector<sparse_matrix_ptrtype> > vector_sparse_matrix;

    using super_type::computeBetaQm;
    typedef boost::tuple<beta_vector_type,  std::vector<beta_vector_type> > beta_type;

    typedef DEIM<self_type> deim_type;
    typedef boost::shared_ptr<deim_type> deim_ptrtype;

    GreplDEIM() :
        super_type( name() ),
        M_use_newton( boption(_name="crb.use-newton") )
    {}

    static std::string name()
    {
        return "grepldeim"+std::to_string(Order)+"-D"+std::to_string(Dim);
    }

    void initModel();

    void assemble();

    vector_ptrtype assembleForDEIMnl( parameter_type const& mu, element_type const& u, int const& tag )
    {
        vector_ptrtype V = this->M_backend->newVector(Xh);
        auto temp = Xh->element(V,0);
        mesh = Xh->mesh();
        temp.on( elements(mesh), cst(mu(0))/cst(mu(1))*( exp( cst(mu(1))*idv(u) ) - 1  ) );
        return V;
    }

    beta_type computeBetaQm( element_type const& T,parameter_type const& mu );
    beta_type computeBetaQm( parameter_type const& mu );
    beta_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu );
    void fillBetaQm(std::vector<vectorN_type*> betas, parameter_type const& mu);

    affine_decomposition_type computeAffineDecomposition();

    beta_vector_type computeBetaInitialGuess( parameter_type const& mu );
    std::vector< std::vector<element_ptrtype> > computeInitialGuessAffineDecomposition( );

    beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu , double time=1e30 );
    std::vector< std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA();

    void assembleResidualWithAffineDecomposition(std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm);
    void assembleJacobianWithAffineDecomposition(std::vector<std::vector<sparse_matrix_ptrtype> > & Jqm);
    void assembleFunctionalWithAffineDecomposition(std::vector<std::vector<sparse_matrix_ptrtype> > & RF_Aqm,
                                                   std::vector< std::vector<std::vector<vector_ptrtype> > >& RF_Fqm);
    bool updateResidual(element_type const& X, std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm);
    void updateResidualMonolithic(vector_ptrtype const& X, vector_ptrtype & R, parameter_type const& mu);
    void updateJacobianMonolithic(vector_ptrtype const& X, sparse_matrix_ptrtype & J, parameter_type const& mu);
    monolithic_type computeMonolithicFormulationU( parameter_type const& mu , element_type const& solution );

    element_type solve( parameter_type const& mu );
    value_type output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false );
private :
    mesh_ptrtype mesh;

    sparse_matrix_ptrtype M_monoA;
    std::vector<vector_ptrtype> M_monoF;

    std::vector< std::vector< element_ptrtype > > M_InitialGuess;
    using super_type::Xh;
    bool M_use_newton;
};


template<int Order, int Dim>
typename GreplDEIM<Order,Dim>::beta_type
GreplDEIM<Order,Dim>::computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
{
    auto beta_g = this->deim()->beta(mu,urb);
    std::vector<vectorN_type*> betas;
    betas.push_back(&beta_g);

    fillBetaQm(betas, mu);
    if( M_use_newton )
        return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
    else
        return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}
template<int Order, int Dim>
typename GreplDEIM<Order,Dim>::beta_type
GreplDEIM<Order,Dim>::computeBetaQm( element_type const& T,parameter_type const& mu )
{
    auto beta_g = this->deim()->beta(mu,T);

    std::vector<vectorN_type*> betas;
    betas.push_back(&beta_g);

    fillBetaQm(betas, mu);
    if( M_use_newton )
        return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
    else
        return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);

    /*int M = this->deim()->size();

    this->M_betaJqm[0][0] = 1;
    this->M_betaJqm[1].resize(M); //needed if cobuild
    for(int m=0; m<M; m++)
        this->M_betaJqm[1][m] = mu(1)*beta_g(m);

    this->M_betaJqm[2][0] = mu(0);

    this->M_betaRqm[0][0][0] = this->computeBetaInitialGuess( mu )[0][0];
    this->M_betaRqm[0][1].resize(M);
    for(int m=0; m<M; m++)
        this->M_betaRqm[0][1][m] = beta_g(m);
    this->M_betaRqm[0][2][0] = -100;
    //output
    this->M_betaRqm[1][0][0] = 1;*/

    //return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
}

template<int Order, int Dim>
typename GreplDEIM<Order,Dim>::beta_type
GreplDEIM<Order,Dim>::computeBetaQm( parameter_type const& mu )
{
    auto beta_g = this->deim()->beta(mu);
    std::vector<vectorN_type*> betas;
    betas.push_back(&beta_g);

    fillBetaQm(betas, mu);

    if( M_use_newton )
        return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
    else
        return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
    /*int M = this->deim()->size();

    this->M_betaJqm[0][0] = 1;
    this->M_betaJqm[1].resize(M); //needed if cobuild
    for(int m=0; m<M; m++)
        this->M_betaJqm[1][m] = mu(1)*beta_g(m);

    this->M_betaJqm[2][0] = mu(0);

    this->M_betaRqm[0][0][0] = this->computeBetaInitialGuess( mu )[0][0];
    this->M_betaRqm[0][1].resize(M);
    for(int m=0; m<M; m++)
        this->M_betaRqm[0][1][m] = beta_g(m);
    this->M_betaRqm[0][2][0] = -100;
    //output
    this->M_betaRqm[1][0][0] = 1;

    return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);*/
}

template<int Order, int Dim>
void GreplDEIM<Order,Dim>::fillBetaQm(std::vector<vectorN_type*> betas, parameter_type const& mu)
{
    int M = this->deim()->size();

    auto beta_g=*betas[0];

    if( M_use_newton )
    {
        this->M_betaJqm[0][0] = 1;
        this->M_betaJqm[1].resize(M); //needed if cobuild
        for(int m=0; m<M; m++)
        {
            this->M_betaJqm[1][m] = mu(1)*beta_g(m);
        }
        this->M_betaJqm[2][0] = mu(0);

        this->M_betaRqm[0][0][0] = this->computeBetaInitialGuess( mu )[0][0];
        this->M_betaRqm[0][1].resize(M);
        for(int m=0; m<M; m++)
        {
            this->M_betaRqm[0][1][m] = beta_g(m);
        }
        this->M_betaRqm[0][2][0] = -100;
        //output
        this->M_betaRqm[1][0][0] = 1;
    }
    else
    {
        this->M_betaAqm[0][0]=1;
        this->M_betaFqm[0][0].resize(M);
        for(int m=0; m<M; m++)
        {
            this->M_betaFqm[0][0][m]=-beta_g(m) ;
        }
        this->M_betaFqm[0][1][0]=100;
        this->M_betaFqm[1][0][0]=1;
    }

}


template<int Order, int Dim>
void GreplDEIM<Order,Dim>::initModel()
{
    if ( !Xh )
    {
        static mesh_ptrtype static_mesh;
        if ( !static_mesh )
            static_mesh = createGMSHMesh( _mesh=new mesh_type,
                                          _desc=domain( _name = "benchmarkgrepl",
                                                        _shape = "hypercube",
                                                        _dim = Dim,
                                                        _h=doption("gmsh.hsize"),
                                                        _xmin=0,_xmax=1,
                                                        _ymin=0,_ymax=1 ) );
        mesh = static_mesh;
        this->setFunctionSpaces( space_type::New(mesh));
    }



    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << Xh->nDof() << std::endl;
        std::cout << "Number of local dof " << Xh->nLocalDof() << std::endl;
    }

    auto mu_min = this->Dmu->element();
    mu_min <<  0.01, 0.01;
    this->Dmu->setMin( mu_min );
    auto mu_max = this->Dmu->element();
    mu_max << 10, 10;
    this->Dmu->setMax( mu_max );

    auto Pset = this->Dmu->sampling();
    std::vector<size_type> Ne(2);
    Ne[0] = ioption("trainset-deim-size");
    Ne[1] = Ne[0];
    std::string samplingname =(boost::format("DmuEim-Ne%1%-generated-by-master-proc") %Ne[0] ).str();
    std::ifstream file ( samplingname );

    if( ! file )
    {
        Pset->equidistributeProduct( Ne , true , samplingname );
        Pset->writeOnFile( samplingname );
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(samplingname);
    }

    auto d = deim( _model=boost::dynamic_pointer_cast<self_type>(this->shared_from_this()),
                         _sampling=Pset );
    this->addDeim( d );
    this->deim()->run();
    int M = this->deim()->size();

    if ( M_use_newton )
    {
        this->M_betaJqm.resize( 3 );
        this->M_betaJqm[0].resize( 1 );
        this->M_betaJqm[1].resize( M );
        this->M_betaJqm[2].resize( 1 );
        this->M_betaRqm.resize( 2 );
        this->M_betaRqm[0].resize( 3 );
        this->M_betaRqm[0][0].resize( 1 );
        this->M_betaRqm[0][1].resize( M );
        this->M_betaRqm[0][2].resize( 1 );
        this->M_betaRqm[1].resize( 1 );
        this->M_betaRqm[1][0].resize( 1 );

        this->M_Jqm.resize( 3 );
        this->M_Jqm[0].resize( 1 );
        this->M_Jqm[0][0] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        this->M_Jqm[1].resize( M );
        for(int m=0; m<M; m++)
            this->M_Jqm[1][m] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        this->M_Jqm[2].resize( 1 );
        this->M_Jqm[2][0] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );

        this->M_Rqm.resize( 2 );
        this->M_Rqm[0].resize( 3 );
        this->M_Rqm[1].resize( 1 );
        this->M_Rqm[0][0].resize(1);
        this->M_Rqm[0][0][0] = backend()->newVector( this->Xh );
        this->M_Rqm[0][1].resize( M );
        for(int m=0; m<M; m++)
            this->M_Rqm[0][1][m] = backend()->newVector( this->Xh );
        this->M_Rqm[0][2].resize(1);
        this->M_Rqm[0][2][0] = backend()->newVector( this->Xh );
        this->M_Rqm[1][0].resize( 1 );
        this->M_Rqm[1][0][0]= backend()->newVector( this->Xh );
    }
    else
    {
        this->M_betaAqm.resize( 1 );
        this->M_betaAqm[0].resize( 1 );
        this->M_betaFqm.resize(2);
        this->M_betaFqm[0].resize( 2 );
        this->M_betaFqm[0][0].resize( M );
        this->M_betaFqm[0][1].resize( 1 );
        this->M_betaFqm[1].resize( 1 );
        this->M_betaFqm[1][0].resize( 1 );

        this->M_Aqm.resize( 1 );
        this->M_Aqm[0].resize( 1 );
        this->M_Aqm[0][0] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );

        this->M_Fqm.resize( 2 );
        this->M_Fqm[0].resize( 2 );
        this->M_Fqm[0][0].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Fqm[0][0][m] = backend()->newVector( this->Xh );
        }
        this->M_Fqm[0][1].resize(1);
        this->M_Fqm[0][1][0]=backend()->newVector( this->Xh );
        this->M_Fqm[1].resize(1);
        this->M_Fqm[1][0].resize(1);
        this->M_Fqm[1][0][0]=backend()->newVector( this->Xh );
    }


    M_InitialGuess.resize(1);
    M_InitialGuess[0].resize(1);
    M_InitialGuess[0][0] = Xh->elementPtr();

    assemble();
}

template<int Order, int Dim>
void
GreplDEIM<Order,Dim>::assemble()
{
    auto u = Xh->element();
    double gamma = doption(_name="gamma");

    if ( M_use_newton )
    {
        assembleJacobianWithAffineDecomposition( this->M_Jqm );
        assembleResidualWithAffineDecomposition( this->M_Rqm );
    }
    else
    {
        auto v = Xh->element();
        int M = this->deim()->size();
        form2( _test=Xh, _trial=Xh, _matrix=this->M_Aqm[0][0] ) =
            integrate( _range= elements( mesh ), _expr = gradt(u)*trans(grad(v)) );
        form2( _test=Xh, _trial=Xh, _matrix=this->M_Aqm[0][0] )
            += integrate( _range = boundaryfaces( mesh ),
                          _expr = gamma*idt(u)*id(v)/hFace()
                          - (gradt(u)*vf::N())*id(v)
                          - (grad(v)*vf::N())*idt(u) );

        this->M_Fqm[0][0].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Fqm[0][0][m] = backend()->newVector( this->Xh );
            v = *(this->deim()->q(m));
            form1( _test=Xh, _vector=this->M_Fqm[0][0][m] ) =
                integrate( _range= elements( mesh ), _expr=( idv(v)*id(v) ) );
            this->M_Fqm[0][0][m]->close();
        }

        form1( _test=Xh, _vector=this->M_Fqm[0][1][0] ) =
            integrate( _range= elements( mesh ),_expr=sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );

        this->M_Fqm[0][1][0]->close();

        form1( _test=Xh, _vector=this->M_Fqm[1][0][0] ) =
            integrate( _range= elements( mesh ),
                       _expr=id(v) );
        this->M_Fqm[1][0][0]->close();
    }
}

template <int Order, int Dim>
void
GreplDEIM<Order,Dim>::assembleJacobianWithAffineDecomposition( std::vector< std::vector<sparse_matrix_ptrtype> >& Jqm)
{
    auto v = Xh->element();
    auto u = Xh->element();

    int M = this->deim()->size();
    double gamma = doption(_name="gamma");

    form2( _test=Xh, _trial=Xh, _matrix=Jqm[0][0] ) =
        integrate( _range= elements( mesh ), _expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=Jqm[0][0] ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idt(u)*id(v)/hFace()
                   - (gradt(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idt(u) );

    Jqm[0][0]->close();
    Jqm[1].resize(M);
    for(int m=0; m<M; m++)
    {
        this->M_Jqm[1][m] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        v = *(this->deim()->q(m));
        form2( _test=Xh, _trial=Xh, _matrix=Jqm[1][m] ) =
            integrate( _range = elements(mesh),
                       _expr = idv(v)*idt(u)*id(v) );
        Jqm[1][m]->close();
    }

    form2( _test=Xh, _trial=Xh, _matrix=Jqm[2][0] ) =
        integrate( _range= elements( mesh ), _expr = idt(u)*id(v) );
    Jqm[2][0]->close();
}

template <int Order, int Dim>
void
GreplDEIM<Order,Dim>::assembleResidualWithAffineDecomposition(std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm)
{
    double gamma = doption(_name="gamma");

    auto u = Xh->element();
    auto v = Xh->element();
    int M = this->deim()->size();

    u = *this->M_InitialGuess[0][0];

    form1( _test=Xh, _vector=Rqm[0][0][0] ) =
        integrate( _range= elements( mesh ), _expr = gradv(u)*trans(grad(v)) );
    form1( _test=Xh, _vector=Rqm[0][0][0] ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idv(u)*id(v)/hFace()
                   - (gradv(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idv(u) );
    Rqm[0][0][0]->close();

    Rqm[0][1].resize(M);
    for( int m=0; m<M; m++ )
    {
        this->M_Rqm[0][1][m] = backend()->newVector( this->Xh );
        v = *(this->deim()->q(m));
        form1( _test=Xh, _vector=Rqm[0][1][m] ) =
            integrate( _range= elements( mesh ), _expr=( idv(v)*id(v) ) );
        Rqm[0][1][m]->close();
    }

    if( Dim == 2 )
        form1( _test=Xh, _vector=Rqm[0][2][0] ) =
            integrate( _range= elements( mesh ), _expr=sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );
    else if( Dim == 3 )
        form1( _test=Xh, _vector=Rqm[0][2][0] ) =
            integrate( _range= elements( mesh ), _expr=sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * id(v) );

    form1( _test=Xh, _vector=Rqm[1][0][0] ) =
        integrate( _range= elements( mesh ), _expr=id(v) );

    Rqm[0][2][0]->close();
    Rqm[1][0][0]->close();
}

template <int Order, int Dim>
bool
GreplDEIM<Order,Dim>::updateResidual(element_type const& X, std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm)
{
    auto u = Xh->element();
    auto v = Xh->element(); //test
    double gamma = option(_name="gamma").template as<double>();

    u = X;
    Rqm[0][0][0] = backend()->newVector( this->Xh );
    form1( _test=Xh, _vector=Rqm[0][0][0] ) =
        integrate( _range= elements( mesh ), _expr = gradv(u)*trans(grad(v)) );
    form1( _test=Xh, _vector=Rqm[0][0][0] ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idv(u)*id(v)/hFace()
                   - (gradv(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idv(u) );

    Rqm[0][0][0]->close();
    return true;
}

template <int Order, int Dim>
void
GreplDEIM<Order,Dim>::updateJacobianMonolithic( vector_ptrtype const& X,
                                                                  sparse_matrix_ptrtype & J,
                                                                  parameter_type const& mu )
{
    //Here we don't use EIM approximation and so no affine decomposition
    auto u = Xh->element();
    auto v = Xh->element(); //test
    u=*X;

    double gamma = doption(_name="gamma");

    J->zero();
    form2( _test=Xh, _trial=Xh, _matrix=J ) =
        integrate( _range= elements( mesh ),_expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=J ) += integrate( _range = boundaryfaces( mesh ),
                                                          _expr = gamma*idt(u)*id(v)/hFace()
                                                          - (gradt(u)*vf::N())*id(v)
                                                          - (grad(v)*vf::N())*idt(u) );

    form2( _test=Xh, _trial=Xh, _matrix=J ) +=
        integrate( _range = elements(mesh), _expr = mu(0)*exp( mu(1)*idv(u) )*idt(u)*id(v) );

    J->close();
}

template <int Order, int Dim>
void
GreplDEIM<Order,Dim>::updateResidualMonolithic(vector_ptrtype const& X,
                                                                 vector_ptrtype & R,
                                                                 parameter_type const& mu )
{
    double gamma = doption(_name="gamma");

    auto u = Xh->element();
    auto v = Xh->element();
    u=*X;

    R->zero();
    form1( _test=Xh, _vector=R ) =
        integrate( _range= elements( mesh ), _expr = gradv(u)*trans(grad(v)) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idv(u)*id(v)/hFace()
                   - (gradv(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idv(u) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ), _expr= mu(0)/mu(1)*( exp( mu(1)*idv(u) ) * id(v) ) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ), _expr= -mu(0)/mu(1)*( id(v)) );

    if( Dim == 2 )
        form1( _test=Xh, _vector=R ) +=
            integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );
    else if( Dim == 3 )
        form1( _test=Xh, _vector=R ) +=
            integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * id(v) );

    R->close();
}


template<int Order, int Dim>
typename GreplDEIM<Order,Dim>::element_type
GreplDEIM<Order,Dim>::solve( parameter_type const& mu )
{
    sparse_matrix_ptrtype J = backend()->newMatrix( Xh, Xh);
    vector_ptrtype R = backend()->newVector( Xh );
    backend()->nlSolver()->jacobian = std::bind( &self_type::updateJacobianMonolithic,
                                                   std::ref( *this ), std::placeholders::_1, std::placeholders::_2, mu );
    backend()->nlSolver()->residual = std::bind( &self_type::updateResidualMonolithic,
                                                   std::ref( *this ), std::placeholders::_1, std::placeholders::_2, mu );

    auto solution = Xh->element();
    backend()->nlSolve(_jacobian=J, _solution=solution, _residual=R);

    return solution;
}

template<int Order, int Dim>
typename GreplDEIM<Order,Dim>::monolithic_type
GreplDEIM<Order,Dim>::computeMonolithicFormulationU( parameter_type const& mu , element_type const& solution)
{
    auto u=Xh->element();
    auto v=Xh->element();
    double gamma = doption(_name="gamma");

    M_monoA = backend()->newMatrix( Xh, Xh );
    M_monoF.resize(2);
    M_monoF[0] = backend()->newVector( this->functionSpace() );
    M_monoF[1] = backend()->newVector( Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M_monoA ) =
        integrate( _range= elements( mesh ), _expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=M_monoA ) += integrate( _range = boundaryfaces( mesh ),
                                                                _expr = gamma*idt(u)*id(v)/hFace()
                                                                - (gradt(u)*vf::N())*id(v)
                                                                - (grad(v)*vf::N())*idt(u) );

    if( Dim == 2 )
        form1( _test=Xh, _vector=M_monoF[0] ) =
            integrate( _range= elements( mesh ), _expr=-mu(0)/mu(1)*( exp( mu(1)*idv(solution) ) - 1 )*id(v) + 100*sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );
    else if ( Dim == 3 )
        form1( _test=Xh, _vector=M_monoF[0] ) =
            integrate( _range= elements( mesh ), _expr=-mu(0)/mu(1)*( exp( mu(1)*idv(solution) ) - 1 )*id(v) + 100*sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * id(v) );

    form1( _test=Xh, _vector=M_monoF[1] ) =
        integrate( _range= elements( mesh ), _expr=id(v) );

    M_monoF[0]->close();
    M_monoF[1]->close();
    M_monoA->close();

    return boost::make_tuple( M_monoA, M_monoF );

}


template<int Order, int Dim>
typename GreplDEIM<Order,Dim>::beta_vector_type
GreplDEIM<Order,Dim>::computeBetaLinearDecompositionA( parameter_type const& mu , double time )
{
    beta_vector_type beta;
    beta.resize(1);
    beta[0].resize(1);
    beta[0][0]=1;
    return beta;
}

template<int Order, int Dim>
typename GreplDEIM<Order,Dim>::vector_sparse_matrix
GreplDEIM<Order,Dim>::computeLinearDecompositionA()
{
    auto Xh=this->Xh;
    auto u=Xh->element();
    auto v=Xh->element();
    double gamma = doption(_name="gamma");

    this->M_linearAqm.resize(1);
    this->M_linearAqm[0].resize(1);
    this->M_linearAqm[0][0] = backend()->newMatrix( Xh, Xh );
    form2( _test=Xh, _trial=Xh, _matrix=this->M_linearAqm[0][0] ) =
        integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) ) +
        integrate( boundaryfaces( mesh ), gamma*idt( u )*id( v )/h()
                   -gradt( u )*vf::N()*id( v )
                   -grad( v )*vf::N()*idt( u )
                   );

    this->M_linearAqm[0][0]->close();
    return this->M_linearAqm;
}


template<int Order, int Dim>
typename GreplDEIM<Order,Dim>::beta_vector_type
GreplDEIM<Order,Dim>::computeBetaInitialGuess( parameter_type const& mu )
{
    this->M_betaInitialGuess.resize( 1 );
    this->M_betaInitialGuess[0].resize( 1 );

    this->M_betaInitialGuess[0][0] = 1;

    return this->M_betaInitialGuess;
}


template<int Order, int Dim>
typename GreplDEIM<Order,Dim>::affine_decomposition_type
GreplDEIM<Order,Dim>::computeAffineDecomposition()
{
    if( M_use_newton )
        return boost::make_tuple( this->M_Jqm, this->M_Rqm );
    else
        return boost::make_tuple( this->M_Aqm, this->M_Fqm );
}

template<int Order, int Dim>
std::vector< std::vector< typename GreplDEIM<Order,Dim>::element_ptrtype > >
GreplDEIM<Order,Dim>::computeInitialGuessAffineDecomposition( )
{
    return M_InitialGuess;
}


template<int Order, int Dim>
double
GreplDEIM<Order,Dim>::output( int output_index, parameter_type const& mu, element_type &solution, bool need_to_solve )
{
    double s=0;
    if ( output_index==0 )
    {
        if( Dim == 2 )
            s = integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py()) * idv( solution ) ).evaluate()(0,0);
        else if ( Dim == 3 )
            s = integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * idv( solution ) ).evaluate()(0,0);
    }
    if( output_index==1 )
    {
        s = integrate( elements( mesh ), idv( solution ) ).evaluate()( 0,0 );
    }

    return s ;
}


#endif
