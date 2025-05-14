#include <iostream>

#include <chrono>
#include <fmt/chrono.h>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/minmax.hpp>
#include <feel/feeldiscr/sensors.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/measure.hpp>
#include <feel/feelmesh/bvh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feells/distancetorange.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
    using Feel::cout;
    
    try
    {
        po::options_description distop( "Distance options" );
        distop.add_options()
            ("h", po::value<double>()->default_value( 0.1 ), "mesh size" )
            ("M", po::value<int>()->default_value( 100 ), "sampling points")
        ;

        Environment env( _argc=argc, _argv=argv,
                        _desc=distop,
                        _about=about(_name="qs_dist_3d",
                                    _author="Feel++ Consortium",
                                    _email="feelpp-devel@feelpp.org") );

        using bvh_ray_type = BVHRay<FEELPP_DIM>;

        // Load parameters
        double h_ = doption(_name="h");
        int M = ioption(_name="M");
        std::string filename = soption(_name="gmsh.filename");

        // Load the mesh        
        tic();
        auto mesh_ = loadMesh( _mesh = new Mesh<Simplex<FEELPP_DIM,1>>, _filename = filename, _h = h_);
        toc("loadMesh");

        // Define exporter
        auto Xh_ = Pch<1>( mesh_, elements( mesh_) );
        auto e_ = Feel::exporter(_mesh = mesh_, _name = "distance_from_boundary" );
        e_->addRegions();

        // FMM
        tic();
        auto d_FMM = distanceToRange(_space = Xh_, _range = markedfaces(mesh_, "boundary"));
        toc("distance_FMM");
        e_->add( "d_FMM", d_FMM );

        // Ray-Tracing with BVH
        tic();

        // Define distance field
        auto d_BVH = Xh_->element();

        // Preprocessing
        auto bvhThirdParty = boundingVolumeHierarchy(_range=markedfaces(mesh_, "boundary"), _kind="third-party");

        // Directional sampling : on définit toutes les directions
        std::vector<Eigen::Vector3d> directions(M);
        double alpha = (8*M_PI)/((1+std::sqrt(5))*(1+std::sqrt(5)));

        for (int k = 0; k < M; k++)
        {
            double zk = 1 - 2 * (k + 0.5) / M;
            double rk = std::sqrt(1 - zk*zk);
            double thetak = alpha * k;
            directions[k] << rk*std::cos(thetak), rk*std::sin(thetak), zk;
        }

        // Ray Casting
        Eigen::Vector3d origin;
        BVHRaysDistributed<FEELPP_DIM> allrays;
        std::vector<size_type> pointIDs;

        for ( auto const& eltWrap : elements(mesh_) ) // on parcourt tous les éléments du maillage
        {
            auto const& elt = unwrap_ref( eltWrap );
            
            for ( int p=0; p<elt.nVertices(); ++p )
            {
                auto const & point = elt.point(p); 
                
                // pas sure si isGhostCell est la bonne fonction 
                if (!point.isOnBoundary() && !point.isGhostCell()) // on ne considère pas les points qui se trouvent au bord ni les dof ghosts
                {
                    size_type id_p = Xh_->dof()->localToGlobal( elt.id(), p ).index(); 

                    if (std::find(pointIDs.begin(), pointIDs.end(), id_p) == pointIDs.end()) // un nouveau point
                    {
                        pointIDs.push_back(id_p);
                        origin << point.node()[0], point.node()[1], point.node()[2]; // définir l'origin du rayon
                
                        for (int k = 0; k < M; k++) // on stocke tous les rayons pour cet origin
                            allrays.push_back(bvh_ray_type(origin,directions[k]));
                    
                    }
                }
            }
        }

        // Get intersection
        auto multiRayIntersectionResult = bvhThirdParty->intersect(_ray=allrays);
    
        // Get distance 
        std::vector<std::vector<double>> dist(pointIDs.size(), std::vector<double>(M, 0.0));
        
        for (auto const& [fid,rirs] : enumerate(multiRayIntersectionResult))
        {
            if (!rirs.empty()) // on check l'intersection
                dist[static_cast<int>(fid / M)][fid % M] = rirs.front().distance();
        }   
        
        // Get min for each point
        for (int k = 0; k < pointIDs.size(); k++)
            d_BVH[pointIDs[k]] = *(std::min_element(dist[k].begin(), dist[k].end())); 
        

        e_->add( "d_BVH", d_BVH );
        toc("raytracing");

        // Compute exact solution
        tic();
        auto d_exact = Xh_->element();
        auto d_exact_expr = expr("min(min(min(min(min(x, 1-x), y), 1-y), z), 1-z):x:y:z");
        d_exact.on(_range=elements(support(Xh_)), _expr = d_exact_expr);
        e_->add( "d_exact", d_exact );
        toc("exact");
        e_->save();

        // Compute errors
        auto L2_error_FMM = normL2( _range=elements(mesh_), _expr = d_exact_expr - idv(d_FMM) );
        auto L2_error_BVH = normL2( _range=elements(mesh_), _expr = d_exact_expr - idv(d_BVH) );
        //auto Linf_error_FMM = normLinf( _range=elements(mesh_), _expr=d_exact_expr - idv(d_FMM) );
        //auto Linf_error_BVH = normLinf( _range=elements(mesh_), _expr=d_exact_expr - idv(d_BVH) );


        if ( Environment::isMasterRank() )
        {
            std::cout << "L2_error_FMM: " << L2_error_FMM << std::endl;
            std::cout << "L2_error_BVH: " << L2_error_BVH << std::endl;
            //std::cout << "Linf_error_FMM: " << Linf_error_FMM << std::endl;
            //std::cout << "Linf_error_BVH: " << Linf_error_BVH << std::endl;
        }

        return 0;
        
    }
    catch(...)
    {
        handleExceptions();
    }
    return 1;
}