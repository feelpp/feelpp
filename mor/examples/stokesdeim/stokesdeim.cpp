#include <feel/feelmor/crbsaddlepointplugin.hpp>

#include "stokesdeim.hpp"

namespace Feel
{
AboutData
    makeStokesDeimAbout( std::string const& str )
{
    Feel::AboutData about( /*AppName  */ str.c_str(),
                           /*ProgName */ str.c_str(),
                           /*Version  */ "0.1",
                           /*ShortDesc*/ "Stokes Application",
                           /*Licence  */ Feel::AboutData::License_GPL,
                           /*Copyright*/ "Copyright (c) 2016 Feel++ Consortium" );
    return about;
}

StokesDeim::StokesDeim() :
    super_type( "StokesDeim" )
{}


void StokesDeim::initModel()
{
    this->setHasDisplacementField(true);

    mesh_ptrtype mesh = loadMesh( _mesh=new mesh_type);
    setFunctionSpaces( space_type::New(mesh));
    cout << "Number of DoF : " << this->functionSpace()->nDof() << std::endl;
    cout << "Number of Local DoF : " << this->functionSpace()->nLocalDof() << std::endl;

    Dmu->setDimension( 2 );
    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    mu_min << 0.1, 0.1;
    mu_max << 0.5, 0.5;
    Dmu->setMin( mu_min );
    Dmu->setMax( mu_max );

    auto Pset = this->Dmu->sampling();
    std::vector<size_type> Ne(2);
    Ne[0] = 4;
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

    auto d = Feel::mdeim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()),
                            _sampling=Pset, _tag=0 );
    d->onlineModel()->setAssembleMDEIM([this,d]( parameter_type const& mu ) {
                                                    return this->assembleForMDEIM( mu, false );
                                                });
    this->addMdeim( d );
    this->mdeim()->run();

    int mMax = this->mdeim()->size();
    this->M_betaAqm.resize(1);
    this->M_betaAqm[0].resize(mMax);
    this->M_Aqm.resize(1);
    this->M_Aqm[0] = this->mdeim()->q();

    this->M_betaFqm.resize(2);
    this->M_betaFqm[0].resize(1);
    this->M_betaFqm[0][0].resize(1);
    this->M_betaFqm[1].resize(1);
    this->M_betaFqm[1][0].resize(1);
    this->M_Fqm.resize(2);
    this->M_Fqm[0].resize(1);
    this->M_Fqm[0][0].resize(1);
    this->M_Fqm[0][0][0]=backend()->newVector(Xh);
    this->M_Fqm[1].resize(1);
    this->M_Fqm[1][0].resize(1);
    this->M_Fqm[1][0][0]=backend()->newVector(Xh);


    auto U = Xh->element();
    auto V = Xh->element();
    auto u = U.template element<0>();
    auto p = U.template element<1>();
    auto v = V.template element<0>();
    auto q = V.template element<1>();

    auto f0 = form1( _test=Xh, _vector=this->M_Fqm[0][0][0] );
    auto f1 = form1( _test=Xh, _vector=this->M_Fqm[1][0][0] );

    auto energy = form2( _trial=Xh, _test=Xh );
    energy += on( _range=markedfaces(mesh,"inlet"), _rhs=f0, _element=u,
                  _expr=4*Py()*(1-Py())*oneX() );

    energy = integrate( _range = elements(mesh), _expr = inner(gradt(u),grad(v)));
    energy += integrate( _range = elements(mesh), _expr = inner(idt(p),id(q)) );
    energy += on( _range=markedfaces(mesh,"inlet"), _rhs=f1, _element=u,
                  _expr=4*Py()*(1-Py())*oneX(), _type="elimination_symmetric" );
    this->addEnergyMatrix( energy );

    f1 = integrate( _range = markedfaces(mesh,"midflux"), _expr = inner(oneX(),id(u)) );
    this->M_Fqm[1][0][0]->close();
}

StokesDeim::beta_type
StokesDeim::computeBetaQm( parameter_type const& mu )
{
    auto Acoeff = this->mdeim()->beta(mu);
    int mMax = this->mdeim()->size();
    for ( int m=0; m<mMax; m++ )
    {
        M_betaAqm[0][m] = Acoeff(m);
    }
    M_betaFqm[0][0][0]=1;
    M_betaFqm[1][0][0]=1./(1.-mu[1]);

    return boost::make_tuple( this->M_betaAqm, M_betaFqm );
}


StokesDeim::sparse_matrix_ptrtype
StokesDeim::assembleForMDEIM( parameter_type const& mu, int const& tag )
{
    tic();
    double nu = 1;
    double mur = 0.3;

    auto Xh=this->functionSpace();
    auto mesh = Xh->mesh();

    element_type U( Xh, "U" );
    element_type V( Xh, "V" );
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = V.template element<1>();
    auto q = V.template element<1>();

    auto J26 = mat<2,2>( cst((2-mur)/(2-mu[0])), cst(0), cst(0), cst(mur/mu[1]) );
    auto detJinv26 = 1./det(J26);
    auto J48 = mat<2,2>( cst((2-mur)/(2-mu[0])),cst(0), cst(0), cst((1-mur)/(1-mu[1])) );
    auto detJinv48 = 1./det(J48);
    auto J3 = mat<2,2>( cst((2-mur)/(2-mu[0])), cst(0), 2*cst((mur-mu[1])/(2-mu[0])), cst(1) );
    auto detJinv3 = 1./det(J3);
    auto J7 = mat<2,2>( cst((2-mur)/(2-mu[0])), cst(0), -2*cst((mur-mu[1])/(2-mu[0])), cst(1) );
    auto detJinv7 = 1./det(J7);
    auto J5 = mat<2,2>( cst(mur/mu[0]), cst(0), cst(0), cst((1-mur)/(1-mu[1])) );
    auto detJinv5 = 1./det(J5);

    auto f2 = form2( _test=Xh, _trial=Xh );
    f2 = integrate( _range = markedelements(mesh, "omega0"),
                    _expr = nu*inner( gradt(u),grad(v) )
                    -idt(p)*div(v) -id(q)*divt(u) );
    f2 += integrate( _range = markedelements( mesh, "omega26"),
                     _expr = detJinv26*( nu*trace(trans(gradt(u)*J26)*(grad(v)*J26))
                                 -idt(p)*trace(grad(v)*J26) -id(q)*trace(gradt(u)*J26) ) );
    f2 += integrate( _range = markedelements( mesh, "omega48"),
                     _expr = detJinv48*( nu*trace(trans(gradt(u)*J48)*(grad(v)*J48))
                                - idt(p)*trace(grad(v)*J48) -id(q)*trace(gradt(u)*J48) ) );
    f2 += integrate( _range = markedelements( mesh, "omega3"),
                     _expr = detJinv3*( nu*trace(trans(gradt(u)*J3)*(grad(v)*J3))
                                 - idt(p)*trace(grad(v)*J3) -id(q)*trace(gradt(u)*J3) ) );
    f2 += integrate( _range = markedelements( mesh, "omega7"),
                     _expr = detJinv7*( nu*trace(trans(gradt(u)*J7)*(grad(v)*J7))
                                - idt(p)*trace(grad(v)*J7) -id(q)*trace(gradt(u)*J7) ) );
    f2 += integrate( _range = markedelements( mesh, "omega5"),
                     _expr = detJinv5*( nu*trace(trans(gradt(u)*J5)*(grad(v)*J5))
                                - idt(p)*trace(grad(v)*J5) -id(q)*trace(gradt(u)*J5) ) );

    auto l = form1( _test = Xh );
    f2 += on( _range=markedfaces(mesh,"inlet"), _rhs=l, _element=u, _expr=zero<2,1>() );
    f2 += on( _range=markedfaces(mesh,"gamma0"), _rhs=l, _element=u, _expr=zero<2,1>() );
    f2 += on( _range=markedfaces(mesh,"gamma26"), _rhs=l, _element=u, _expr=zero<2,1>() );
    f2 += on( _range=markedfaces(mesh,"gamma48"), _rhs=l, _element=u, _expr=zero<2,1>() );
    f2 += on( _range=markedfaces(mesh,"gamma5"), _rhs=l, _element=u, _expr=zero<2,1>() );

    toc("assemblefordeim");
    return f2.matrixPtr();
}


StokesDeim::value_type
StokesDeim::output( int output_index, parameter_type const& mu ,
                    element_type& u, bool need_to_solve )
{
    double output = 0.0;
    auto mesh = this->functionSpace()->mesh();
    auto u0 = u.template element<0>();

    if ( output_index==0 )
    {

    }
    else if ( output_index==1 )
    {
        output = 2./(1.-mu[1])*integrate( _range = markedfaces(mesh,"midflux"), _expr = inner(oneX(),idv(u0)) ).evaluate()(0,0);
    }

    return output;
}


StokesDeim::displacement_field_ptrtype
StokesDeim::meshDisplacementField( parameter_type const& mu )
{
    auto mesh = this->functionSpace()->mesh();
    if ( !M_Dh )
    {
        M_Dh = displacement_space_type::New( mesh );
        mapping = M_Dh->elementPtr();
    }
    else
        mapping->zero();

    double mur = 0.3;
    double mu1 = mu[0];
    double mu2 = mu[1];

    mapping->on( _range = markedelements( mesh, "omega3"),
                _expr = vec( (mur-mu1)/(2-mur)*Px() - cst((mur-mu1)/(2-mur)),
                     -2*(mur-mu2)/(2-mur)*Px() + cst( 2*(mur-mu2)/(2-mur) ) ) );
    mapping->on( _range = markedelements( mesh, "omega5"),
                _expr = vec( (mu1-mur)/mur*Px() + cst(2*(mur-mu1)/mur),
                     (mur-mu2)/(1-mur)*Py() - cst((mur-mu2)/(1-mur)) ));
    mapping->on( _range = markedelements( mesh, "omega7"),
                _expr = vec( (mur-mu1)/(2-mur)*Px() - cst(3*(mur-mu1)/(2-mur)),
                     2*(mur-mu2)/(2-mur)*Px() - cst( 6*(mur-mu2)/(2-mur) ) ) );
    mapping->on( _range = markedelements( mesh, "omega26"),
                _expr = vec( (mur-mu1)/(2-mur)*Px() - cst((mur-mu1)/(2-mur)), (mu2-mur)/mur*Py() )
                + chi(Px()>2)*vec(  - cst(2*(mur-mu1)/(2-mur)), cst(0.) ) );
    mapping->on( _range = markedelements( mesh, "omega48"),
                 _expr = vec( (mur-mu1)/(2-mur)*Px() - cst((mur-mu1)/(2-mur)),
                      (mur-mu2)/(1-mur)*Py() - cst((mur-mu2)/(1-mur)) )
                 + chi(Px()>2)*vec(  - cst(2*(mur-mu1)/(2-mur)), cst(0.) ));

    return mapping;
}





} // namespace Feel
