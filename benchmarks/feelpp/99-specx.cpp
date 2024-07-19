
//#include <fmt/core.h>
//#include <fmt/chrono.h>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/enumerate.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <fmt/chrono.h>
#include <specx/Data/SpDataAccessMode.hpp>
#include <specx/Legacy/SpRuntime.hpp>
#include <specx/Task/SpPriority.hpp>
#include <specx/Utils/SpArrayView.hpp>
#include <specx/Task/SpProbability.hpp>

//#include <spdlog/spdlog.h>
//#include <spdlog/fmt/chrono.h>


template <typename Iterator>
std::vector<std::tuple<Iterator, Iterator, int>>
divide_work( Iterator begin, Iterator end, std::size_t n )
{
    std::vector<std::tuple<Iterator, Iterator,int>> ranges;
    ranges.reserve(n);

    auto dist = std::distance(begin, end);
    auto chunk = dist / n;
    auto remainder = dist % n;

    for (size_t i = 0; i < n-1; ++i) {
        auto next_end = std::next(begin, chunk + (remainder ? 1 : 0));
        ranges.emplace_back(begin, next_end,i);

        begin = next_end;
        if (remainder) remainder -= 1;
    }

    // last chunk
    ranges.emplace_back(begin, end,n-1);
    return ranges;
}

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
//    spdlog::info("load mesh in {}", std::chrono::duration_cast<std::chrono::seconds>(end - beg));
    
  
    const int NumThreads = ioption("nt"); //SpUtils::DefaultNumThreads();
    SpRuntime<SpSpeculativeModel::SP_MODEL_2> runtime( NumThreads );

    // Next we set a predicate that will be called by the runtime each
    // time a speculative task becomes ready to run. It is used to
    // decide if the speculative task should be allowed to run.
    runtime.setSpeculationTest(
        []( [[maybe_unused]] const int nbReadyTasks,
            [[maybe_unused]] const SpProbability& meanProbability ) -> bool
        {
            return true; // Here we always return true, this basically means
                         // that we always allow speculative tasks to run
                         // regardless of runtime conditions.
        } );

    auto [r, meas] = std::tuple{ elements( mesh ), 1 };
    double area = 0;
    auto it = r.get<1>(), en = r.get<2>();
    auto ranges = divide_work( it, en, NumThreads );
    beg = std::chrono::high_resolution_clock::now();
    
    Eigen::VectorXd area_vec(NumThreads);
    area_vec.setZero();
    for( auto const& r : ranges )
    {
        runtime.task( //SpWriteArray( area_vec.data(), SpArrayView( NumThreads ).removeItems( 0, NumThreads ).addItem( std::get<2>( r ) ) ),
                      SpWrite( area_vec( std::get<2>( r ) ) ),
                      [&r]( double& a ) -> double
                      {
                          double buffer = 0;
                          for ( auto it = std::get<0>( r ), en = std::get<1>( r ); it != en; ++it )
                          {
                              auto const& f = boost::unwrap_ref( *it );
                              buffer += std::sqrt(f.measure());
                          }
                          a += buffer;
                          return a;
                      } )
            .setTaskName( fmt::format( "area-{}", std::get<2>( r ) ) );
    }
    runtime.task(SpReadArray( area_vec.data(), SpArrayView( NumThreads )),  
                 [&area](SpArrayAccessor<const double>& a)
                    {
                        Eigen::Map<const Eigen::VectorXd> v(*a.begin(),a.getSize());
                        area = v.sum();
                    }).setTaskName("area");

    // We wait for all tasks to finish
    runtime.waitAllTasks();
#if !defined( USE_FINAL_TASK )
    area = area_vec.sum();
#endif    
    // We make all runtime threads exit
    runtime.stopAllThreads();

    end = std::chrono::high_resolution_clock::now();
    std::cout << fmt::format("area domain specx({} threads): {} in {}", ioption("nt"), area, std::chrono::duration_cast<std::chrono::milliseconds>(end - beg)) << std::endl;
    // We generate the task graph corresponding to the execution
    runtime.generateDot( "/tmp/99-specx.dot", true );

    // We generate an Svg trace of the execution
    runtime.generateTrace( "/tmp/99-specx.svg" );
#if 0
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
    //spdlog::info("area domain tf({} threads): {} in {}", ioption("nt"), area, std::chrono::duration_cast<std::chrono::milliseconds>(end - beg));
    taskflow.dump(std::cout);
#endif
    beg = std::chrono::high_resolution_clock::now();
    double area_seq = 0;
    for (auto const& wf : r)
    {
        auto const &f = boost::unwrap_ref(wf);
        area_seq += std::sqrt(f.measure());
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << fmt::format("area domain standard: {} in {}", area_seq, std::chrono::duration_cast<std::chrono::milliseconds>(end - beg)) << std::endl;
    //spdlog::info("area domain seq: {} in {}", area_seq, std::chrono::duration_cast<std::chrono::milliseconds>(end - beg));

//    if ( std::abs(area-area_seq)>1e-10 || std::abs(area-meas) > 1e-10 )
//        spdlog::error("error in computing area: {} {} should be approx {}\n", area, area_seq, meas);

}