#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feells/reinit_fms.hpp>
#include <feel/feells/disttocurve.hpp>
#include <feel/feells/curveparametrizations.hpp>
#include <feel/feelvf/vf.hpp>

#undef FEELPP_HAS_GMSH_H

#include <feel/feelmesh/meshadaptation.hpp>

using namespace Feel;


Feel::po::options_description
meshAdaptFromLSOptions()
{
    Feel::po::options_description opts;
    opts.add_options()
        ( "h-min-list", Feel::po::value< std::vector<double> >(), "min mesh sizes" )
        ( "h-max-list", Feel::po::value< std::vector<double> >(), "max mesh sizes" )
        ( "phi-power", Feel::po::value< double >()->default_value(2.), "power of phi for mesh adaptation" );
        return opts;
}


int main( int argc, char** argv )
{
  const int dim = FEELPP_DIM;

  Feel::Environment env( _argc=argc, _argv=argv, _desc=meshAdaptFromLSOptions() );

  auto mesh = unitHypercube<dim>();

  // some ellipse parameters:
  const double x0=0.6;
  const double y0=0.5;
  const double z0 = dim==2 ? 0 : 0.5;
  const double aAxis = 0.1;
  const double bAxis = 0.3;
  auto X0 = Px() - x0;
  auto Y0 = Py() - y0;
  auto Z0 = Pz() - z0;

  const std::vector<double> hmin_list = Environment::vm()["h-min-list"].as< std::vector<double> >();
  const std::vector<double> hmax_list = Environment::vm()["h-max-list"].as< std::vector<double> >();
  CHECK(hmax_list.size()==hmin_list.size())<<"hmin and hmax should have the same size\n";

  const int maxit = hmax_list.size();
  LOG(INFO)<<"number of iteration of mesh adaptation to do = "<< maxit <<std::endl;

  typedef MeshAdaptation<dim, 1, 1> mesh_adapt_type;

  const std::string geofile = Environment::findFile( "hypercube.geo" );
  CHECK(!geofile.empty())<<"File "<< geofile << " can't be found\n";

  auto exp = exporter(_mesh=mesh, _geo="change", _name="meshadaptation");
  mesh_adapt_type mesh_adaptation;

  for (int i=0; i<maxit; ++i)
  {
      const double hmin = hmin_list[i];
      const double hmax = hmax_list[i];

      Feel::cout<< Feel::tc::red << "start making mesh adaptation with:"
                << "hmin="<< hmin << " and hmax="<< hmax << tc::reset << std::endl;

      tic();
      auto Xh1  = mesh_adapt_type::p1_space_type::New( mesh );
      auto Xh0 = mesh_adapt_type::p0_space_type::New( mesh );
      toc("created spaces");


      tic();
      // ellipse function (not exactly a distance function)
      auto phi = vf::project(_space=Xh1, _range=elements(mesh),
                                      _expr=sqrt( (X0/aAxis) * (X0/aAxis)
                                                  + (Y0/bAxis) * (Y0/bAxis)
                                                  + (Z0/bAxis) * (Z0/bAxis) ) - 1 );
      toc("projected phi");

      // ---------- make metric -------------
      tic();
      const double phiMax = phi.max();
      const double phi_power = doption("phi-power");
      const double phiMaxPowered = std::pow(phiMax, phi_power);

      auto linphi = (hmax-hmin)/phiMaxPowered * vf::pow(vf::abs(idv(phi)), phi_power) + hmin;
      auto metric = vf::project(_space=Xh1, _range=elements(mesh),
                                _expr = max(hmin, linphi) );
      toc("projected metric");

      // ------------ perform mesh adaptation --------
      auto metric_list = mesh_adaptation.makeMetricList(metric, "metric" );

      exp->setMesh(mesh);
      exp->step(i)->add("metric", metric);
      exp->step(i)->add("phi", phi);
      exp->step(i)->addRegions();
      exp->save();

      tic();
      mesh = mesh_adaptation.adaptMesh(_mesh=mesh,
                                       _geo=geofile,
                                       _type="isotropic",
                                       _metric=metric_list,
                                       _hmin=hmin,
                                       _hmax=hmax );
      toc("made mesh adaptation");

  }

}
