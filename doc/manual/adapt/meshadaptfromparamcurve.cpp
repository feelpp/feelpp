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
        ( "d-min-list", Feel::po::value< std::vector<double> >(), "h=hmin if d<dmin" )
        ( "d-max-list", Feel::po::value< std::vector<double> >(), "h=hmax if d>dmax" )
        // ( "phi-power", Feel::po::value< double >()->default_value(2.), "power of phi for mesh adaptation" )
        ;
        return opts;
}


int main( int argc, char** argv )
{
  const int dim = 2;

  Feel::Environment env( _argc=argc, _argv=argv, _desc=meshAdaptFromLSOptions() );

  auto mesh = unitHypercube<dim>();

  // ------------- epitrochoid parameters --------------
  const double a_epi=0.1; // 1 / nb_branch
  const double b_epi=0.8;
  auto x_epi = [&](double t) -> double { return ((1+a_epi) * cos(a_epi*t) - a_epi*b_epi * cos( (1+a_epi) * t )) / 4 + 0.5; };
  auto y_epi = [&](double t) -> double { return ((1+a_epi) * sin(a_epi*t) - a_epi*b_epi * sin( (1+a_epi) * t )) / 4 + 0.5; };
  // ----------------------------------------

  const std::vector<double> hmin_list = Environment::vm()["h-min-list"].as< std::vector<double> >();
  const std::vector<double> hmax_list = Environment::vm()["h-max-list"].as< std::vector<double> >();
  const std::vector<double> dmin_list = Environment::vm()["d-min-list"].as< std::vector<double> >();
  const std::vector<double> dmax_list = Environment::vm()["d-max-list"].as< std::vector<double> >();
  CHECK(  (hmax_list.size()==hmin_list.size()) && (hmax_list.size()==dmin_list.size()) && (hmax_list.size()==dmax_list.size()) )
      <<"hmin, hmax, dmin, dmax should have the same size\n";

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
      const double dmin = dmin_list[i];
      const double dmax = dmax_list[i];

      Feel::cout<< Feel::tc::red << "start making mesh adaptation with:\n"
                << "hmin="<< hmin << " and hmax="<< hmax << "\n"
                << "dmin="<< dmin << " and dmax="<< dmax
                << tc::reset << std::endl;

      tic();
      auto Xh1  = mesh_adapt_type::p1_space_type::New( mesh );
      auto Xh0 = mesh_adapt_type::p0_space_type::New( mesh );
      toc("created spaces");

      tic();
      // ------- make level set function ------
      // distance to curve
      auto disttocurve = distToCurve( Xh0, Xh1 );
      // fast marching
      auto fm = fms( Xh1 );
      auto epitro = disttocurve->fromParametrizedCurve( x_epi, y_epi,
                                                        0, 100, hmin/10. );
      toc("made distance to parametrized curve");

      tic();
      *epitro = fm->march( epitro );
      toc("fast marching");

      // ---------- make metric -------------
      tic();
      // const double phiMax = epitro->max();
      // const double phi_power = doption("phi-power");
      // const double phiMaxPowered = std::pow(phiMax, phi_power);

      // auto linphi = (hmax-hmin)/phiMaxPowered * vf::pow(vf::abs(idv(epitro)), phi_power) + hmin;

      const double alpha = (hmin-hmax)/(dmin-dmax);
      const double beta = hmin-alpha*dmin;
      auto phi = abs(idv(epitro));

      auto metric = vf::project(_space=Xh1, _range=elements(mesh),
                                _expr = ( phi <= dmin ) * hmin
                                + ( phi >= dmax ) * hmax
                                + ( phi > dmin )*( phi < dmax )*(alpha*phi+beta) );
      toc("projected metric");

      // ---------- perform mesh adapatation ----------
      auto metric_list = mesh_adaptation.makeMetricList( metric, "metric");

      exp->setMesh(mesh);
      exp->step(i)->add("metric", metric);
      exp->step(i)->add("epitro", *epitro);
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
