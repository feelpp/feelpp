/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feelvf/contractions.hpp>
#include <feel/feelvf/vf.hpp>

#include "qs_shell_benchmark_framework.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numbers>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Feel;

namespace
{
namespace qsb = Feel::Quickstart::ShellBenchmark;
using json = qsb::json;
using qsb::benchmarkSpec;
using qsb::jsonArray;
using qsb::loadBenchmarkSpecs;
using qsb::shellVectorExpression;

enum class FormulationKind
{
    Standard,
    Tensor,
    Mandel,
    Voigt
};

struct BenchmarkRunConfig
{
    qsb::BenchmarkConfig benchmark;
    std::vector<int> orders;
    std::vector<FormulationKind> formulations;
};

node_type
toNode( qsb::shell_vec const& p )
{
    node_type n( 3 );
    n( 0 ) = p( 0 );
    n( 1 ) = p( 1 );
    n( 2 ) = p( 2 );
    return n;
}

template<typename T>
void
pushUnique( std::vector<T>& values, T const& value )
{
    if ( std::find( values.begin(), values.end(), value ) == values.end() )
        values.push_back( value );
}

FormulationKind
formulationKind( std::string name )
{
    std::transform( name.begin(), name.end(), name.begin(),
                    []( unsigned char c ) { return static_cast<char>( std::tolower( c ) ); } );

    if ( name == "standard" || name == "classic" || name == "qs_elasticity" )
        return FormulationKind::Standard;
    if ( name == "tensor" )
        return FormulationKind::Tensor;
    if ( name == "mandel" )
        return FormulationKind::Mandel;
    if ( name == "voigt" )
        return FormulationKind::Voigt;

    throw std::invalid_argument( "unknown formulation '" + name + "'" );
}

std::string
formulationName( FormulationKind kind )
{
    switch ( kind )
    {
    case FormulationKind::Standard: return "standard";
    case FormulationKind::Tensor: return "tensor";
    case FormulationKind::Mandel: return "mandel";
    case FormulationKind::Voigt: return "voigt";
    }
    return "unknown";
}

std::vector<int>
parseOrdersList( std::string const& value )
{
    std::vector<int> orders;
    std::stringstream ss( value );
    std::string item;
    while ( std::getline( ss, item, ',' ) )
    {
        item.erase( std::remove_if( item.begin(), item.end(),
                                    []( unsigned char c ) { return std::isspace( c ) != 0; } ),
                    item.end() );
        if ( item.empty() )
            continue;

        int const order = std::stoi( item );
        if ( order < 1 || order > 3 )
            throw std::invalid_argument( "only Q1, Q2 and Q3 are supported, got Q" + std::to_string( order ) );
        pushUnique( orders, order );
    }
    return orders;
}

std::vector<FormulationKind>
parseFormulationsList( std::string const& value )
{
    std::vector<FormulationKind> kinds;
    std::stringstream ss( value );
    std::string item;
    while ( std::getline( ss, item, ',' ) )
    {
        item.erase( std::remove_if( item.begin(), item.end(),
                                    []( unsigned char c ) { return std::isspace( c ) != 0; } ),
                    item.end() );
        if ( item.empty() )
            continue;

        pushUnique( kinds, formulationKind( item ) );
    }
    return kinds;
}

BenchmarkRunConfig
benchmarkRunConfig( json const& specs, std::string const& name )
{
    BenchmarkRunConfig cfg;
    cfg.benchmark = qsb::benchmarkPreset( specs, name );
    auto const& benchmarkSpecs = benchmarkSpec( specs, name );

    if ( benchmarkSpecs.contains( "orders" ) )
    {
        for ( auto const& entry : jsonArray( benchmarkSpecs, "orders" ) )
        {
            int const order = entry.get<int>();
            if ( order < 1 || order > 3 )
                throw std::invalid_argument( "only Q1, Q2 and Q3 are supported, got Q" + std::to_string( order ) );
            pushUnique( cfg.orders, order );
        }
    }
    if ( cfg.orders.empty() )
        cfg.orders = { 1, 2, 3 };

    if ( benchmarkSpecs.contains( "formulations" ) )
    {
        for ( auto const& entry : jsonArray( benchmarkSpecs, "formulations" ) )
            pushUnique( cfg.formulations, formulationKind( entry.get<std::string>() ) );
    }
    if ( cfg.formulations.empty() )
        cfg.formulations = { FormulationKind::Standard, FormulationKind::Mandel, FormulationKind::Voigt };

    return cfg;
}

qsb::MeshBuildResult
buildMesh( qsb::BenchmarkConfig const& cfg, int order, FormulationKind formulation )
{
    qsb::MeshBuildOptions options;
    std::string const formulationLabel = formulationName( formulation );
    options.repositoryVariantPath = "Q" + std::to_string( order ) + "/" + formulationLabel;
    options.meshVariantTag = "Q" + std::to_string( order ) + "-" + formulationLabel;
    return qsb::buildMesh( cfg, options );
}

std::string
requiredSpecsPath()
{
    if ( !Environment::vm().count( "specs" ) || soption( "specs" ).empty() )
        throw std::invalid_argument( "missing required option 'specs'; set it in qs_elasticity_qk_benchmark.cfg or pass --specs" );
    return Environment::expand( soption( "specs" ) );
}

template<typename PointRangeType>
double
pointLoadScale( PointRangeType const& pointRange, std::string const& quantity, std::string const& label )
{
    if ( quantity == "per-point" || quantity == "per_point" || quantity == "point" )
        return 1.0;

    CHECK( quantity == "total" ) << label << " quantity must be 'per-point' or 'total'";
    size_type const nMarkedPoints = nelements( pointRange, true );
    CHECK( nMarkedPoints > 0 ) << label << " uses quantity=total but selects no marked point";
    return 1.0 / nMarkedPoints;
}

qsb::shell_voigt
toShellVoigtFromSymmetricStorage( Eigen::VectorXd const& values, bool engineeringShear )
{
    double const shearScale = engineeringShear ? 2.0 : 1.0;
    qsb::shell_voigt result;

    if ( values.size() >= 9 )
    {
        result << values( 0 ),                                      // xx
                  values( 4 ),                                      // yy
                  values( 8 ),                                      // zz
                  shearScale * 0.5 * ( values( 1 ) + values( 3 ) ), // xy
                  shearScale * 0.5 * ( values( 2 ) + values( 6 ) ), // xz
                  shearScale * 0.5 * ( values( 5 ) + values( 7 ) ); // yz
        return result;
    }

    CHECK( values.size() >= 6 ) << "expected at least six symmetric tensor components";
    result << values( 0 ),                  // xx
              values( 3 ),                  // yy
              values( 5 ),                  // zz
              shearScale * values( 1 ),     // xy
              shearScale * values( 2 ),     // xz
              shearScale * values( 4 );     // yz
    return result;
}

template <typename SpaceType, typename FieldElementType>
qsb::shell_voigt
evaluateSymmetricFieldReference( SpaceType const& Sh,
                                 FieldElementType const& field,
                                 qsb::shell_vec const& point,
                                 bool engineeringShear )
{
    auto componentValue = [&field,&point]( ComponentType c1, ComponentType c2 )
    {
        auto component = field.comp( c1, c2 );
        auto ctx = component.functionSpace()->context();
        ctx.add( toNode( point ) );
        auto values = component.evaluate( ctx, true );
        CHECK( values.size() >= 1 ) << "expected at least one scalar component value";
        return values( 0 );
    };

    double const shearScale = engineeringShear ? 2.0 : 1.0;
    qsb::shell_voigt result;
    result << componentValue( ComponentType::X, ComponentType::X ),
              componentValue( ComponentType::Y, ComponentType::Y ),
              componentValue( ComponentType::Z, ComponentType::Z ),
              shearScale * componentValue( ComponentType::X, ComponentType::Y ),
              shearScale * componentValue( ComponentType::X, ComponentType::Z ),
              shearScale * componentValue( ComponentType::Y, ComponentType::Z );
    return result;
}

template <int Order, typename SpaceType, typename EpsilonElementType, typename SigmaElementType>
bool
checkFieldReferences( qsb::BenchmarkConfig const& cfg,
                      SpaceType const& Sh,
                      EpsilonElementType const& epsilonh,
                      SigmaElementType const& sigmah,
                      FormulationKind formulation,
                      bool checkReference )
{
    std::string const exactTarget = "qk-q" + std::to_string( Order ) + "-" + formulationName( formulation );
    std::string const fallbackTarget = ( Order == 1 && formulation == FormulationKind::Standard ) ? "hexa8" : "";
    bool ok = true;

    for ( auto const& fieldReference : cfg.fieldReferences )
    {
        std::string target;
        if ( fieldReference.targets.find( exactTarget ) != fieldReference.targets.end() )
            target = exactTarget;
        else if ( !fallbackTarget.empty() && fieldReference.targets.find( fallbackTarget ) != fieldReference.targets.end() )
            target = fallbackTarget;
        else if ( fieldReference.targets.find( "qk" ) != fieldReference.targets.end() )
            target = "qk";
        else
            continue;

        auto const [epsilonValue, sigmaValue] =
            std::make_pair( evaluateSymmetricFieldReference( Sh, epsilonh, fieldReference.point, true ),
                            evaluateSymmetricFieldReference( Sh, sigmah, fieldReference.point, false ) );
        ok = qsb::reportFieldReferenceComparison( fieldReference, target, epsilonValue, sigmaValue,
                                                  checkReference, std::cout ) && ok;
    }

    return ok;
}

po::options_description
makeOptions()
{
    po::options_description options( "qs_elasticity_qk_benchmark options" );
    options.add_options()
        ( "specs", po::value<std::string>(),
          "json benchmark specification file (typically provided by qs_elasticity_qk_benchmark.cfg)" )
        ( "benchmark", po::value<std::string>()->default_value( "square-plate" ),
          "benchmark name in the JSON specification file" )
        ( "orders", po::value<std::string>()->default_value( "" ),
          "comma-separated list of displacement orders to run (allowed: 1,2,3)" )
        ( "formulations", po::value<std::string>()->default_value( "" ),
          "comma-separated list of formulations to run (standard,tensor,mandel,voigt)" )
        ( "E", po::value<double>()->default_value( -1.0 ), "override Young modulus from the benchmark preset" )
        ( "nu", po::value<double>()->default_value( -1.0 ), "override Poisson ratio from the benchmark preset" )
        ( "pressure", po::value<double>()->default_value( -1.0 ),
          "override every pressure load magnitude; positive pressure acts inward as -p*N()" )
        ( "functions.f", po::value<std::string>()->default_value( "" ), "override body force expression" )
        ( "no-solve", po::value<bool>()->default_value( false ), "assemble only" )
        ( "show-timings", po::value<bool>()->default_value( true ),
          "print per-term assembly and solve timings" )
        ( "check-reference", po::value<bool>()->default_value( true ),
          "fail when a benchmark reference value exceeds its configured tolerance" )
        ( "export-results", po::value<bool>()->default_value( false ), "export the displacement field" );
    return options.add( feel_options() );
}

AboutData
makeAbout()
{
    AboutData about( "qs_elasticity_qk_benchmark",
                     "qs_elasticity_qk_benchmark",
                     "0.1",
                     "Qk elasticity benchmark on the SB9 one-layer shell cases using the standard quickstart elasticity formulation plus tensor, Mandel and Voigt variants",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) Feel++ Consortium" );
    about.addAuthor( "Feel++ Consortium", "developer", "feelpp-devel@feelpp.org", "" );
    return about;
}

template<int Order>
bool
runBenchmarkCase( qsb::BenchmarkConfig const& cfg, FormulationKind formulation )
{
    using namespace Feel;
    using namespace vf;

    auto meshData = buildMesh( cfg, Order, formulation );
    auto mesh = meshData.mesh;
    auto Uh = Pchv<Order>( mesh );

    auto u = trial( Uh, "u" );
    auto v = test( Uh, "v" );
    auto uh = Uh->element( "u" );
    auto algebraBackend = backend( _rebuild=true, _worldcomm=Uh->worldCommPtr() );

    double const E = ( doption( "E" ) > 0.0 ) ? doption( "E" ) : cfg.E;
    double const nu = ( doption( "nu" ) >= 0.0 ) ? doption( "nu" ) : cfg.nu;
    double const pressureOverride = doption( "pressure" );
    double const lambda = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
    double const mu = E / ( 2.0 * ( 1.0 + nu ) );
    std::string const bodyForceExpression =
        soption( "functions.f" ).empty() ? cfg.bodyForceExpression : soption( "functions.f" );

    std::cout << "benchmark: " << cfg.name << "\n"
              << "description: " << cfg.description << "\n"
              << "specs: " << requiredSpecsPath() << "\n"
              << "displacement space: Q" << Order << "\n"
              << "formulation: " << formulationName( formulation ) << "\n"
              << "E=" << E << ", nu=" << nu << "\n"
              << "loads: pressures=" << cfg.pressureLoads.size()
              << ", tractions=" << cfg.tractionLoads.size()
              << ", total-forces=" << cfg.totalForceLoads.size()
              << ", point-loads=" << cfg.pointLoads.size()
              << ", point-moments=" << cfg.pointMoments.size() << "\n"
              << "constraints: face-clamps=" << cfg.clampMarkers.size()
              << ", point-constraints=" << cfg.pointConstraints.size() << "\n";

    bool const showTimings = boption( "show-timings" );
    auto ticIf = [showTimings]()
    {
        if ( showTimings )
            tic();
    };
    auto tocIf = [showTimings]( std::string const& label )
    {
        if ( showTimings )
            toc( label, true );
    };

    auto f = expr<3, 1>( bodyForceExpression, "f" );

    ticIf();
    auto l = form1( _test=Uh,
                    _backend=algebraBackend,
                    _vector=algebraBackend->newVector( _test=Uh ) );
    l = integrate( _range=elements( mesh ), _expr=inner( f, v ) );
    tocIf( "elasticity.l.body-force" );

    for ( auto const& pressureLoad : cfg.pressureLoads )
    {
        double const pressureValue = ( pressureOverride >= 0.0 ) ? pressureOverride : pressureLoad.value;
        ticIf();
        l += integrate( _range=markedfaces( mesh, pressureLoad.marker ),
                        _expr=inner( -cst( pressureValue ) * N(), v ) );
        tocIf( "elasticity.l.pressure" );
    }

    for ( auto const& tractionLoad : cfg.tractionLoads )
    {
        auto traction = expr<3, 1>( tractionLoad.expression, "traction" );
        ticIf();
        l += integrate( _range=markedfaces( mesh, tractionLoad.marker ),
                        _expr=inner( traction, v ) );
        tocIf( "elasticity.l.traction" );
    }

    for ( auto const& totalForceLoad : cfg.totalForceLoads )
    {
        double const loadedArea = measure( _range=markedfaces( mesh, totalForceLoad.marker ) );
        CHECK( loadedArea > 0.0 ) << "marked face '" << totalForceLoad.marker
                                  << "' has zero measure, cannot distribute a total force";
        auto traction = expr<3, 1>( shellVectorExpression( totalForceLoad.value / loadedArea ), "traction" );
        ticIf();
        l += integrate( _range=markedfaces( mesh, totalForceLoad.marker ),
                        _expr=inner( traction, v ) );
        tocIf( "elasticity.l.total-force" );
    }

    for ( auto const& pointLoad : cfg.pointLoads )
    {
        auto pointRange = markedpoints( mesh, pointLoad.marker );
        double const scale = pointLoadScale( pointRange, pointLoad.quantity, "point load '" + pointLoad.marker + "'" );
        auto pointLoadExpr = expr<3, 1>( shellVectorExpression( scale * pointLoad.value ), "point_load" );
        ticIf();
        l += integrate( _range=pointRange,
                        _expr=inner( pointLoadExpr, id( v ) ) );
        tocIf( "elasticity.l.point-load" );
    }

    for ( auto const& pointMoment : cfg.pointMoments )
    {
        auto pointRange = markedpoints( mesh, pointMoment.marker );
        double const scale = pointLoadScale( pointRange, pointMoment.quantity, "point moment '" + pointMoment.marker + "'" );
        auto pointMomentExpr = expr<3, 1>( shellVectorExpression( scale * pointMoment.value ), "point_moment" );
        ticIf();
        l += integrate( _range=pointRange,
                        _expr=inner( pointMomentExpr, omega( v ) ) );
        tocIf( "elasticity.l.point-moment" );
    }

    ticIf();
    auto a = form2( _trial=Uh,
                    _test=Uh,
                    _backend=algebraBackend,
                    _matrix=algebraBackend->newMatrix( _test=Uh, _trial=Uh ) );
    auto epsu = symm_grad( u );
    auto epsv = symm_grad( v );
    switch ( formulation )
    {
    case FormulationKind::Standard:
    {
        auto C = isotropic_stiffness<3>( lambda, mu );
        auto sigma = [=]( auto const& w )
        {
            return ddot( C, symm_grad( w ) );
        };
        a = integrate( _range=elements( mesh ), _expr=inner( sigma( u ), epsv ) );
        break;
    }
    case FormulationKind::Tensor:
    {
        auto C = isotropic_stiffness<3>( lambda, mu );
        a = integrate( _range=elements( mesh ), _expr=ddot( C, epsu, epsv ) );
        break;
    }
    case FormulationKind::Mandel:
    {
        auto C = isotropic_stiffness<3, SymmetricTensorNotation::Mandel>( lambda, mu );
        a = integrate( _range=elements( mesh ),
                       _expr=contract( C, mandel( epsu ), mandel( epsv ) ) );
        break;
    }
    case FormulationKind::Voigt:
    {
        auto C = isotropic_stiffness<3, SymmetricTensorNotation::Voigt>( lambda, mu );
        a = integrate( _range=elements( mesh ),
                       _expr=voigt_contract( C, voigt( epsu ), voigt( epsv ) ) );
        break;
    }
    }
    tocIf( "elasticity.a.elastic" );

    ticIf();
    for ( auto const& marker : cfg.clampMarkers )
    {
        a += on( _range=markedfaces( mesh, marker ),
                 _rhs=l,
                 _element=uh,
                 _expr=zero<3, 1>() );
    }
    tocIf( "elasticity.a.face-constraints" );

    auto uhX = uh[ComponentType::X];
    auto uhY = uh[ComponentType::Y];
    auto uhZ = uh[ComponentType::Z];
    auto pointConstraintComponent = [&]( qsb::BenchmarkConfig::PointConstraint const& constraint, int c )
    {
        if ( !constraint.components[c] )
            return;

        switch ( c )
        {
        case 0:
            a += on( _range=markedpoints( mesh, constraint.marker ),
                     _rhs=l,
                     _element=uhX,
                     _expr=cst( constraint.value( 0 ) ) );
            break;
        case 1:
            a += on( _range=markedpoints( mesh, constraint.marker ),
                     _rhs=l,
                     _element=uhY,
                     _expr=cst( constraint.value( 1 ) ) );
            break;
        case 2:
            a += on( _range=markedpoints( mesh, constraint.marker ),
                     _rhs=l,
                     _element=uhZ,
                     _expr=cst( constraint.value( 2 ) ) );
            break;
        }
    };

    ticIf();
    for ( auto const& constraint : cfg.pointConstraints )
        for ( int c = 0; c < 3; ++c )
            pointConstraintComponent( constraint, c );
    tocIf( "elasticity.a.point-constraints" );

    if ( !boption( "no-solve" ) )
    {
        ticIf();
        a.solve( _rhs=l, _solution=uh );
        tocIf( "elasticity.solve" );
    }

    auto dispNormMax = normLinf( _range=elements( mesh ),
                                 _pset=_Q<std::max( 2, 2*Order )>(),
                                 _expr=norm2( idv( uh ) ) );
    auto probeCtx = Uh->context();
    probeCtx.add( toNode( meshData.probe ) );
    auto const probeValue = uh.evaluate( probeCtx, false );
    CHECK( probeValue.size() >= 3 ) << "expected a 3-component probe value";

    double probeDisplacement = 0.0;
    for ( int i = 0; i < 3; ++i )
        probeDisplacement += cfg.probeDirection( i ) * probeValue( i );

    std::cout << "max displacement norm = " << dispNormMax.value()
              << " at " << dispNormMax.arg().transpose() << "\n"
              << cfg.probeLabel << " = " << probeDisplacement << "\n";

    bool referenceOk = true;
    if ( cfg.reference.hasProbeValue )
    {
        double const absoluteError = probeDisplacement - cfg.reference.probeValue;
        double const relativeError = ( std::abs( cfg.reference.probeValue ) > 0.0 ) ?
                                         std::abs( absoluteError / cfg.reference.probeValue ) :
                                         std::abs( absoluteError );
        std::cout << "reference " << cfg.probeLabel << " = " << cfg.reference.probeValue;
        if ( !cfg.reference.description.empty() )
            std::cout << " (" << cfg.reference.description << ")";
        std::cout << "\n"
                  << "absolute error = " << absoluteError << "\n"
                  << "relative error = " << relativeError << "\n";

        if ( boption( "check-reference" ) &&
             ( cfg.reference.hasAbsoluteTolerance || cfg.reference.hasRelativeTolerance ) )
        {
            if ( cfg.reference.hasAbsoluteTolerance )
            {
                std::cout << "absolute tolerance = " << cfg.reference.absoluteTolerance << "\n";
                referenceOk = referenceOk && ( std::abs( absoluteError ) <= cfg.reference.absoluteTolerance );
            }
            if ( cfg.reference.hasRelativeTolerance )
            {
                std::cout << "relative tolerance = " << cfg.reference.relativeTolerance << "\n";
                referenceOk = referenceOk && ( relativeError <= cfg.reference.relativeTolerance );
            }
            std::cout << "reference check = " << ( referenceOk ? "PASS" : "FAIL" ) << "\n";
        }
    }

    auto Sh = Pdhms<Order-1>( mesh );
    auto epsilonTensor = cst( 0.5 ) * ( gradv( uh ) + trans( gradv( uh ) ) );
    auto sigmaTensor = cst( lambda )*trace( epsilonTensor )*eye<3, 3>() + cst( 2.0*mu )*epsilonTensor;
    auto epsilonh = vf::project( _space=Sh, _range=elements( mesh ), _expr=epsilonTensor );
    auto sigmah = vf::project( _space=Sh, _range=elements( mesh ), _expr=sigmaTensor );
    referenceOk = checkFieldReferences<Order>( cfg, Sh, epsilonh, sigmah,
                                               formulation, boption( "check-reference" ) ) && referenceOk;

    if ( boption( "export-results" ) )
    {
        auto e = exporter( _mesh=mesh,
                           _name=( cfg.name + "-Q" + std::to_string( Order ) + "-" + formulationName( formulation ) ) );
        e->addRegions();
        e->add( "u", uh );
        e->add( "epsilon", epsilonh );
        e->add( "sigma", sigmah );
        e->save();
    }

    return referenceOk;
}

template<int Order>
bool
runBenchmarkOrders( qsb::BenchmarkConfig const& cfg, std::vector<FormulationKind> const& formulations )
{
    bool ok = true;
    for ( auto formulation : formulations )
        ok = runBenchmarkCase<Order>( cfg, formulation ) && ok;
    return ok;
}
} // namespace

int
main( int argc, char** argv )
{
    try
    {
        Environment env( _argc=argc, _argv=argv, _desc=makeOptions(), _about=makeAbout() );

        auto const specs = loadBenchmarkSpecs( requiredSpecsPath() );
        auto const cfg = benchmarkRunConfig( specs, soption( "benchmark" ) );

        auto orders = cfg.orders;
        if ( !soption( "orders" ).empty() )
            orders = parseOrdersList( soption( "orders" ) );

        auto formulations = cfg.formulations;
        if ( !soption( "formulations" ).empty() )
            formulations = parseFormulationsList( soption( "formulations" ) );

        bool ok = true;
        for ( int order : orders )
        {
            switch ( order )
            {
            case 1:
                ok = runBenchmarkOrders<1>( cfg.benchmark, formulations ) && ok;
                break;
            case 2:
                ok = runBenchmarkOrders<2>( cfg.benchmark, formulations ) && ok;
                break;
            case 3:
                ok = runBenchmarkOrders<3>( cfg.benchmark, formulations ) && ok;
                break;
            default:
                throw std::invalid_argument( "only Q1, Q2 and Q3 are supported, got Q" + std::to_string( order ) );
            }
        }

        return ok ? 0 : 1;
    }
    catch ( ... )
    {
        handleExceptions();
    }
    return 1;
}
