
//#include <fmt/core.h>
//#include <fmt/chrono.h>
#include <feel/feelcore/environment.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/reduce.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/chrono.h>

int main(int argc, char** argv )
{
    using namespace Feel;
    po::options_description viewfactoroptions("view factor options");
    viewfactoroptions.add_options()("nt", po::value<int>()->default_value(1), "number of threads");
    // initialize Feel++ Environment
    Environment env(_argc = argc, _argv = argv,
                    _desc=viewfactoroptions,
                    _about = about(_name = "factor",
                                   _author = "Christophe Prud'homme",
                                   _email = "feelpp-devel@feelpp.org"));

    using mesh_t = Mesh<Simplex<3,1,3>>;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;
    auto beg = std::chrono::high_resolution_clock::now();
    auto mesh = loadMesh(_mesh = new mesh_t);
    auto end = std::chrono::high_resolution_clock::now();
    spdlog::info("load mesh in {}", std::chrono::duration_cast<std::chrono::seconds>(end - beg));

    tf::Executor executor(ioption("nt"));
    tf::Taskflow taskflow;
    //auto [r,meas] = std::tuple{boundaryfaces(mesh),6};
    //auto [r, meas] = std::tuple{markedfaces(mesh,"Walls"), 6};
    //auto [r, meas] = std::tuple{markedfaces(mesh,"1"), 1};
    auto  [r,meas] = std::tuple{elements(mesh),1};
    double area = 0;
    beg = std::chrono::high_resolution_clock::now();
    taskflow.transform_reduce(
        r.get<1>(), r.get<2>(), area, 
        [](double a, double b)
        { 
            return a+b; 
        },
        [](auto const &wf)
        {
            auto const &f = boost::unwrap_ref(wf);
            return f.measure();
        }).name("compute_measure");
    executor.run(taskflow).get();
    end = std::chrono::high_resolution_clock::now();
    spdlog::info("area domain tf({} threads): {} in {}", ioption("nt"), area, std::chrono::duration_cast<std::chrono::milliseconds>(end - beg));
    taskflow.dump(std::cout);
    beg = std::chrono::high_resolution_clock::now();
    double area_seq = 0;
    for (auto const& wf : r)
    {
        auto const &f = boost::unwrap_ref(wf);
        area_seq += f.measure();
    }
    end = std::chrono::high_resolution_clock::now();
    spdlog::info("area domain seq: {} in {}", area_seq, std::chrono::duration_cast<std::chrono::milliseconds>(end - beg));

    if ( std::abs(area-area_seq)>1e-10 || std::abs(area-meas) > 1e-10 )
        spdlog::error("error in computing area: {} {} should be approx {}\n", area, area_seq, meas);
    
}