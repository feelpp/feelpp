#include <feel/feelmodels/maxwell/maxwell.hpp>
#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelfilters/loadmesh.hpp>

using namespace Feel;

void printCvg( std::ostream& out, std::vector<double> h, std::vector<std::vector<double>> err )
{
    auto errA = err[0];
    auto errRelA = err[1];
    auto errB = err[2];
    auto errRelB = err[3];

    boost::format fmter("%1% %|14t|");
    boost::format fmter_endl("%1%\n");

    out << fmter % "h";
    out << fmter % "errA_L2";
    out << fmter % "errA_Curl";
    out << fmter % "errRelA_L2";
    out << fmter % "errRelA_Curl";
    out << fmter % "errB_L2";
    out << fmter % "errB_Div";
    out << fmter % "errRelB_L2";
    out << fmter_endl % "errRelB_Div";

    int maxiter = h.size();
    for( int i = 0; i < maxiter; ++i )
    {
        out << fmter % h[i];
        out << fmter % errA[2*i];
        out << fmter % errA[2*i+1];
        out << fmter % errRelA[2*i];
        out << fmter % errRelA[2*i+1];
        out << fmter % errB[2*i];
        out << fmter % errB[2*i+1];
        out << fmter % errRelB[2*i];
        out << fmter_endl % errRelB[2*i+1];
    }
}

// template <int OrderT>
void
runApplicationMaxwell()
{
    using model_type = FeelModels::Maxwell< Simplex<FEELPP_DIM> >;
    using mesh_type = typename model_type::mesh_type;
    using exporter_type = Exporter<mesh_type>;
    using exporter_ptrtype = std::shared_ptr<exporter_type>;

    std::shared_ptr<model_type> maxwell( new model_type("maxwell", false) );

    int maxiter = ioption("cvg.nb-iter");
    double factor = doption("cvg.factor");
    double size = doption("gmsh.hsize");

    exporter_ptrtype e( exporter_type::New("convergence") );

    std::vector<double> h(maxiter);
    std::vector<double> normA(2*maxiter), errA(2*maxiter), errRelA(2*maxiter);
    std::vector<double> normB(2*maxiter), errB(2*maxiter), errRelB(2*maxiter);

    for( int i = 0; i < maxiter; ++i )
    {
        h[i] = size;
        auto mesh = loadMesh( _mesh=new mesh_type, _h=size );
        maxwell->setMesh(mesh);
        maxwell->init();
        maxwell->printAndSaveInfo();
        maxwell->solve();

        auto functions = maxwell->modelProperties().functions();
        auto params = maxwell->modelProperties().parameters().toParameterValues();
        functions.setParameterValues( params );

        e->step(i)->setMesh(mesh);

        auto A_h = maxwell->fieldMagneticPotential();
        auto AExpr = expr<FEELPP_DIM,1>(functions["A_exact"].expressionString() );
        AExpr.setParameterValues( params );
        auto Ah = maxwell->spaceMagneticPotential();
        auto A_ex = Ah->element(AExpr);
        normA[2*i] = normL2( _range=Ah->template rangeElements<0>(),
                             _expr=idv(A_h) );
        errA[2*i] = normL2( _range=Ah->template rangeElements<0>(),
                            _expr=idv(A_h)-idv(A_ex) );
        errRelA[2*i] = errA[2*i]/normA[2*i];
        normA[2*i+1] = normL2( _range=Ah->template rangeElements<0>(),
                               _expr=idv(A_h) );
        normA[2*i+1] += normL2( _range=Ah->template rangeElements<0>(),
                                _expr=curlv(A_h) );
        errA[2*i+1] = normL2( _range=Ah->template rangeElements<0>(),
                              _expr=idv(A_h)-idv(A_ex) );
        errA[2*i+1] += normL2( _range=Ah->template rangeElements<0>(),
                               _expr=curlv(A_h)-curlv(A_ex) );
        errRelA[2*i+1] = errA[2*i+1]/normA[2*i+1];
        e->step(i)->add("A_h", "A_h", A_h);
        e->step(i)->add("A_ex", "A_ex", A_ex);

        auto B_h = maxwell->fieldMagneticField();
#if FEELPP_DIM==3
        auto BExpr = expr<FEELPP_DIM,1>(functions["B_exact"].expressionString() );
#else
        auto BExpr = expr(functions["B_exact"].expressionString() );
#endif
        BExpr.setParameterValues( params );
        auto Bh = maxwell->spaceMagneticField();
        auto B_ex = Bh->element(BExpr);
        normB[2*i] = normL2( _range=Bh->template rangeElements<0>(),
                             _expr=idv(B_h) );
        errB[2*i] = normL2( _range=Bh->template rangeElements<0>(),
                            _expr=idv(B_h)-idv(B_ex) );
        errRelB[2*i] = errB[2*i]/normB[2*i];
#if FEELPP_DIM==3
        normB[2*i+1] = normL2( _range=Bh->template rangeElements<0>(),
                               _expr=idv(B_h) );
        normB[2*i+1] += normL2( _range=Bh->template rangeElements<0>(),
                               _expr=divv(B_h) );
        errB[2*i+1] = normL2( _range=Bh->template rangeElements<0>(),
                              _expr=idv(B_h)-idv(B_ex) );
        errB[2*i+1] += normL2( _range=Bh->template rangeElements<0>(),
                              _expr=divv(B_h)-divv(B_ex) );
        errRelB[2*i+1] = errB[2*i+1]/normB[2*i+1];
#endif
        e->step(i)->add("B_h", "B_h", B_h);
        e->step(i)->add("B_ex", "N_ex", B_ex);

        e->save();
        size /= factor;
    }

    std::vector<std::vector<double>> err = {errA, errRelA, errB, errRelB};
    Feel::cout << std::endl;
    if( Environment::isMasterRank() )
    {
        printCvg( std::cout, h, err );
        std::ofstream file ( soption("cvg.filename") );
        if( file )
        {
            printCvg( file, h, err );
            file.close();
        }
    }
}

int main(int argc, char** argv) {
    po::options_description maxwelloptions( "application maxwell options" );
    maxwelloptions.add( toolboxes_options("maxwell") );
    maxwelloptions.add_options()
        ("cvg.nb-iter", po::value<int>()->default_value(1), "number of convergence iteration")
        ("cvg.factor", po::value<double>()->default_value(2), "factor of mesh size by iteration")
        ("cvg.filename", po::value<std::string>()->default_value("cvg.dat"), "filename to store convergence results")
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=maxwelloptions,
                     _about=about(_name="cvg_maxwell",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    runApplicationMaxwell();

    return 0;
}
