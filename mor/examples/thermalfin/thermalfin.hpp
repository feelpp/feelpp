/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#ifndef FEELPP_THERMAL_FIN_HPP
#define FEELPP_THERMAL_FIN_HPP 1

#include <boost/timer.hpp>

#include <feel/options.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelmor/modelcrbbase.hpp>


namespace Feel
{

po::options_description
makeThermalFinOptions()
{
    po::options_description thermalfinoptions( "Thermal Fin options" );
    return thermalfinoptions;
}
AboutData
makeThermalFinAbout( std::string const& str = "thermalfin" )
{
    Feel::AboutData about( /*AppName  */ str.c_str(),
                           /*ProgName */ str.c_str(),
                           /*Version  */ "0.1",
                           /*ShortDesc*/ "Thermal Fin",
                           /*Licence  */ Feel::AboutData::License_GPL,
                           /*Copyright*/ "Copyright (c) 2009-2014 Feel++ Consortium" );
    return about;
}

class ThermalFin : public ModelCrbBase<ParameterSpaceX, decltype(Pch<3>(Mesh<Simplex<2>>::New()))>
{
public:
    using super = ModelCrbBase<ParameterSpaceX, decltype(Pch<3>(Mesh<Simplex<2>>::New()))>;
    ThermalFin() : super( "thermafin" ) {}
    
    //! initialisation of the model
    void initModel();

    beta_vector_light_type beta;
    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);
};


void

ThermalFin::initModel()
{

    CHECK( is_linear && !is_time_dependent ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&NonLinear) == NonLinear) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&TimeDependent) == TimeDependent) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    /*
     * First we create the mesh
     */
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<3>( mesh );
    this->setFunctionSpaces( Xh );

    Dmu->setDimension( 6 );
    //static const int N = 2;
    auto mu_min = Dmu->element();
    mu_min << 0.1, 0.1, 0.1, 0.1, 0.001, 1;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 10, 10, 10, 10, 1, 4;
    Dmu->setMax( mu_max );

    auto u = Xh->element();
    auto v = Xh->element();
    //lhs
    auto a0 = form2( _trial=Xh, _test=Xh);
    a0 = integrate(_range=markedelements(mesh,"omega0"),
                   _expr=gradt(u)*trans(grad(v)) );
    this->addLhs( { a0 , "1" } );

    auto a1x = form2( _trial=Xh, _test=Xh);
    a1x =  integrate(_range=markedelements(mesh,"omega1"),
                     _expr= 2.5*dxt(u)*trans(dx(v)) );
    this->addLhs( { a1x , "mu0/mu5" } );
    auto a1y = form2( _trial=Xh, _test=Xh);
    a1y =  integrate(_range=markedelements(mesh,"omega1"),
                     _expr= 1./2.5*dyt(u)*trans(dy(v)) );
    this->addLhs( { a1y , "mu0*mu5" } );

    auto a2x = form2( _trial=Xh, _test=Xh);
    a2x =  integrate(_range=markedelements(mesh,"omega2"),
                     _expr= 2.5*dxt(u)*trans(dx(v)) );
    this->addLhs( { a2x , "mu1/mu5" } );
    auto a2y = form2( _trial=Xh, _test=Xh);
    a2y =  integrate(_range=markedelements(mesh,"omega2"),
                     _expr= 1./2.5*dyt(u)*trans(dy(v)) );
    this->addLhs( { a2y , "mu1*mu5" } );

    auto a3x = form2( _trial=Xh, _test=Xh);
    a3x =  integrate(_range=markedelements(mesh,"omega3"),
                     _expr= 2.5*dxt(u)*trans(dx(v)) );
    this->addLhs( { a3x , "mu2/mu5" } );
    auto a3y = form2( _trial=Xh, _test=Xh);
    a3y =  integrate(_range=markedelements(mesh,"omega3"),
                     _expr= 1./2.5*dyt(u)*trans(dy(v)) );
    this->addLhs( { a3y , "mu2*mu5" } );

    auto a4x = form2( _trial=Xh, _test=Xh);
    a4x =  integrate(_range=markedelements(mesh,"omega4"),
                     _expr= 2.5*dxt(u)*trans(dx(v)) );
    this->addLhs( { a4x , "mu3/mu5" } );
    auto a4y = form2( _trial=Xh, _test=Xh);
    a4y =  integrate(_range=markedelements(mesh,"omega4"),
                     _expr= 1./2.5*dyt(u)*trans(dy(v)) );
    this->addLhs( { a4y , "mu3*mu5" } );


    auto a5 = form2( _trial=Xh, _test=Xh  );
    a5 = integrate(_range=markedfaces(mesh,"gamma"),
                   _expr= idt(u)*id(v) );
    this->addLhs( { a5, "mu4" } );

    //rhs
    auto f0 = form1( _test=Xh );
    f0 = integrate(_range=markedfaces(mesh,"root"),
                   _expr=id(v));
    this->addRhs( { f0, "1" } );

    auto out = form1( _test=Xh );
    double meas = integrate( markedfaces(mesh,"root"), cst(1.) ).evaluate()(0,0);
    out = integrate( markedfaces(mesh,"root"), id(v)/meas );
    this->addOutput( { out, "1" } );

    auto energy = form2( _trial=Xh, _test=Xh);
    energy = integrate( elements(mesh), gradt(u)*trans(grad(v)) )
        + integrate( markedfaces(mesh,"gamma"), 0.1*idt(u)*id(v) );
    this->addEnergyMatrix( energy );
}


double
ThermalFin::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{

//CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    auto mesh = Xh->mesh();
    double output=0;
    // right hand side (compliant)
    if ( output_index == 0 )
    {
        output = integrate( markedfaces(mesh,"root"), idv(u) ).evaluate()(0,0);

    }
    else if ( output_index == 1)
    {
        output = mean( markedfaces(mesh,"root"), idv(u) )(0,0);
    }
    else
        throw std::logic_error( "[ThermalFin::output] error with output_index : only 0 or 1 " );
    return output;
}

}

#endif /* FEELPP_THERMAL_FIN_HPP */
