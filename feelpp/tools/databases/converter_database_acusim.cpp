/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4*/

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/converteracusimdatabase.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description converteracusimnoptions( "converter acusim options" );
	converteracusimnoptions.add_options()
        ( "iname", po::value<std::string>(), "input acusim database repository" )
        ( "idir", po::value<std::string>(), "input acusim database repository" )
        ( "acusim.work-dir", po::value<std::string>()->default_value( "ACUSIM.DIR" ), " acusim work repository" )
        ( "odir", po::value<std::string>(), "output feelpp database repository" )
        ( "fields", po::value<std::vector<std::string> >()->multitoken(), "fields name to convert" )
        ( "build-only-scalar-space", po::value<bool>()->default_value(true), "build-only-scalar-space" )
		;

	Environment env( _argc=argc, _argv=argv,
                     _desc=converteracusimnoptions,
                   _about=about(_name="converter_acusim_database",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    std::string idir;
    if ( Environment::vm().count("idir") )
        idir = Environment::vm()["idir"].as<std::string>();
    std::string iname;
    if ( Environment::vm().count("iname") )
        iname = Environment::vm()["iname"].as<std::string>();
    std::string acusimWorkDir = soption(_name="acusim.work-dir");
    std::string odir;
    if ( Environment::vm().count("odir") )
        odir = Environment::vm()["odir"].as<std::string>();

    std::vector<std::string> fieldsToConvert;
    if ( Environment::vm().count("fields"))
        fieldsToConvert = Environment::vm()["fields"].as<std::vector<std::string> >();


    if ( idir.empty() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --idir is missing\n";
        return 0;
    }
    if ( iname.empty() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --iname is missing\n";
        return 0;
    }
    if ( odir.empty() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --odir is missing\n";
        return 0;
    }

    ConverterAcusimDatabase<Mesh<Simplex<3>>> myconverter;
    myconverter.setAcusimProblemName( iname );
    myconverter.setAcusimRepository( idir );
    myconverter.setAcusimWorkRepository( acusimWorkDir );
    myconverter.setOutputRepository( odir );
    myconverter.setBuildOnlyScalarSpace( boption(_name="build-only-scalar-space" ) );
    for (std::string const& fieldName : fieldsToConvert )
        myconverter.addField( fieldName );

    myconverter.run();

    return 0;
}
