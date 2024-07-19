/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include "../convection.hpp"

using namespace Feel::vf;

void ConvectionCrb::initModel()
{
    // -- MESH SETTING -- //
    mesh_ptrtype mesh;
    if ( soption( "gmsh.filename" ) == "untitled.geo" )
        throw std::logic_error( "[Room CRB Model:initialization] You did not provide a geometry, this model wont work on default geometry" );
    else
    {
        A_OUT << "Mesh read in file : " << soption( "gmsh.filename" ) << "\n";
        LOG(INFO) << "Mesh read in file : " << soption( "gmsh.filename" ) << "\n";
        mesh  =  loadMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER);
    }

    // -- FUNCTION SPACES AND ELEMENTS SETTING -- //
    this->setFunctionSpaces( space_type::New(mesh) );
    pT = element_ptrtype( new element_type( Xh ) );
    A_OUT << "Number of dof : "<< Xh->nDof() << std::endl;
    LOG(INFO) << "Number of dof : "<< Xh->nDof();
    element_type U( Xh, "u" );
    element_type V( Xh, "v" );

    element_u_type u = U. element<0>(); // velocity
    element_u_type v = V. element<0>();
    element_p_type p = U. element<1>(); // pressure
    element_p_type q = V. element<1>();
    element_l_type xi = U. element<2>(); // lagrange
    element_l_type eta = V. element<2>();
    element_t_type t = U. element<3>(); // temperature
    element_t_type s = V. element<3>();

    // -- PARAMETERS SETTING -- //
    Dmu->setDimension( 4 );
    auto mu_min = Dmu->element();
    mu_min << doption( "TinletMin" ), doption( "Tcpu1Min" ), doption("Tcpu2Min"), doption("UinletMin") ;
    auto mu_max = Dmu->element();
    mu_max << doption( "TinletMax" ), doption( "Tcpu1Max" ), doption("Tcpu2Max"), doption("UinletMax") ;
    Dmu->setMin( mu_min );
    Dmu->setMax( mu_max );

    double penalbc = doption("penalbc");
    double nu = doption("nu");
    double rho = doption("rho");
    double k = doption("k");
    double expansion = 9.81*3.4112e-3;
    double Cp = 1000;

    // -- INITIALIZE MATRIX -- //
    // lhs
    M_Aqm.resize( Qa() );
    M_Aqm[0].resize(1);
    M_Aqm[0][0] = M_backend->newMatrix( _test = Xh, _trial = Xh );

    // rhs
    M_Fqm.resize( Nl() );
    for (int ii=0; ii<Nl(); ii++)
    {
        M_Fqm[ii].resize( Ql(ii) );
        for (int jj=0; jj<Ql(ii); jj++)
        {
            M_Fqm[ii][jj].resize(1);
            M_Fqm[ii][jj][0] = M_backend->newVector( Xh );
        }
    }
    // others
    D = M_backend->newMatrix( _test = Xh, _trial = Xh );
    F = M_backend->newVector( Xh );
    M_A_tril = M_backend->newMatrix( _test = Xh, _trial = Xh );

    // -- BUILD MATRIX -- //
    // fluid diffusion
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0], _init=true )
        = integrate( _range = elements( mesh ),
                     _expr = nu*trace( gradt( u )*trans( grad( v ) ) )  );
    // pressure-velocity terms
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate ( _range = elements( mesh ),
                       _expr = - 1./rho*idt( p )*div( v )
                       + 1./rho*divt( u )*id( q ) );
    // heat diffusion
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate( _range = elements( mesh ),
                      _expr = k*gradt( t )*trans( grad( s ) ) );

    // buyoancy forces c(theta,v)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        +=integrate( _range = elements( mesh ),
                     _expr = -expansion*idt( t )*( trans( oneY() )*id( v ) ) );
    // multipliers for zero-mean pressure (not use in this case)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate ( _range = elements( mesh ),
                       _expr = idt( xi )*id( eta ) );

    // -- BOUNDARY CONDITIONS -- //
    auto SigmaNt = ( -1./rho*idt( p )*N() + nu*gradt( u )*N() );
    auto SigmaN = ( -1./rho*id( q )*N() + nu*grad( v )*N() );

    // weak Dirichlet condition on the walls (u=0)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate (_range =  markedfaces( mesh, "wall" ),
                      _expr =  - trans( SigmaNt )*id(v)
                       - trans( SigmaN )*idt(u)
                       + nu*penalbc*trans( idt(u) )*id(v)/hFace() );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate ( _range = markedfaces( mesh, "cpu1" ),
                       _expr = - trans( SigmaNt )*id(v)
                       - trans( SigmaN )*idt(u)
                       + nu*penalbc*trans( idt(u) )*id(v)/hFace() );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate ( _range = markedfaces( mesh, "cpu2" ),
                       _expr = - trans( SigmaNt )*id(v)
                       - trans( SigmaN )*idt(u)
                       + nu*penalbc*trans( idt(u) )*id(v)/hFace() );

    // weak Dirichlet on inlet (u=U0)
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate ( _range = markedfaces( mesh, "inlet" ),
                       _expr = - trans( SigmaNt )*id(v)
                       - trans( SigmaN )*idt(u)
                       + nu*penalbc*trans( idt(u) )*id(v)/hFace() );
#if CONVECTION_DIM==2
    auto poiseuille = (1.00-Py())*(Py()-0.90)*400*oneX();
#else
    auto poiseuille = (1.00-Py())*(Py()-0.90)*(0.35-Pz())*(Pz()-0.15)*40000*oneX();
#endif
    form1( _test = Xh, _vector=M_Fqm[0][3][0], _init=true )
        = integrate ( _range = markedfaces( mesh, "inlet"),
                      _expr = - trans( SigmaN )*poiseuille
                      + nu*penalbc/hFace()*trans(poiseuille)*id(v) );

    // weak Dirichlet on temperature (inlet + cpu1 )
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate ( _range = markedfaces( mesh, "inlet" ),
                       _expr = - k*gradt( t )*N()*id( s )
                       - k*grad( s )*N()*idt( t )
                       + k*penalbc*idt( t )*id( s )/hFace() );
    form1( _test=Xh, _vector=M_Fqm[0][0][0] )
        = integrate ( _range = markedfaces( mesh, "inlet"),
                      _expr = - k*grad(s)*N()
                      + k*penalbc*id(s)/hFace()  );
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate ( _range = markedfaces( mesh, "cpu1" ),
                       _expr = - k*gradt( t )*N()*id( s )
                       - k*grad( s )*N()*idt( t )
                       + k*penalbc*idt( t )*id( s )/hFace() );
    form1(_test= Xh, _vector=M_Fqm[0][1][0] )
        = integrate ( _range = markedfaces( mesh, "cpu1"),
                      _expr = - k*grad(s)*N()
                      + k*penalbc*id(s)/hFace() );

    // robin condition on temperature (cpu2)
    double k_cpu2 = doption("kCpu2");
    form2( _test=Xh, _trial=Xh, _matrix=M_Aqm[0][0] )
        += integrate (_range =  markedfaces( mesh, "cpu2" ),
                      _expr =  id( s )*k_cpu2/rho/Cp*idt(t) );
    form1(_test =  Xh, _vector=M_Fqm[0][2][0] )
        = integrate ( _range = markedfaces( mesh, "cpu2"),
                      _expr = id(s)*k_cpu2/rho/Cp );

    for (int i=0 ; i<Ql(0) ; i++ )
        M_Fqm[0][i][0]->close();


    M = M_backend->newMatrix(_test = Xh, _trial = Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M, _init=true ) =
        integrate( _range = elements( mesh ),
                   _expr = trans( id( v ) )*idt( u ) + trace( grad( v )*trans( gradt( u ) ) )
                   + id( q )*idt( p ) + grad( q )*trans( gradt( p ) )
                   + id( s )*idt( t ) + grad( s )*trans( gradt( t ) ) );
    M->close();

    // Output 1 : mean temperature on cpu2
    double cpu2_domain = integrate( _range = markedfaces(mesh,"cpu2"), _expr = cst(1.) ).evaluate()(0,0);
    form1( _test=Xh, _vector=M_Fqm[1][0][0] ) =
        integrate( _range = markedfaces( mesh,"cpu2" ),
                   _expr = id(s)/cpu2_domain );
    M_Fqm[1][0][0]->close();

    LOG(INFO) << "Natural Convection : Cavity Model is initilized";
}
