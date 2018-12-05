
#include <feel/feelcore/environment.hpp>
#include <feel/feelalg/datamap.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feelalg/vectorublas.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description benchPerfOptions( "benchmark_perf_io_vector options" );
	benchPerfOptions.add_options()
        ( "ndof", po::value<int>()->default_value( 30000000 ), "coeff" )
		;

	Environment env( _argc=argc, _argv=argv,
                     _desc=benchPerfOptions,
                     _about=about(_name="benchmark_perf_io_vector",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto comm = Environment::worldComm();
    size_type nProc = comm.size();

    size_type nDof = ioption(_name="ndof");
    size_type nLocalDofWithoutGhost = nDof/nProc;
    size_type nDofNotDistributed = nDof%nProc;

    //std::cout<< "nDof=" << nDof << "\n";
    //std::cout<< "nLocalDofWithoutGhost=" << nLocalDofWithoutGhost << "\n";
    //std::cout<< "nDofNotDistributed=" << nDofNotDistributed << "\n";

    std::shared_ptr<DataMap> map( new DataMap( comm ) );
    map->setNDof( nDof );
    for ( rank_type p=0;p<nProc;++p )
        map->setNLocalDofWithoutGhost( p, nLocalDofWithoutGhost );
    for ( rank_type p=0;p<nDofNotDistributed;++p )
        map->setNLocalDofWithoutGhost( p, map->nLocalDofWithoutGhost( p )+1 );
    for ( rank_type p=0;p<nProc;++p )
    {
        map->setNLocalDofWithGhost( p, map->nLocalDofWithoutGhost( p )  );
        map->setLastDof(p,map->nLocalDofWithoutGhost( p ));
        //map->setLastDofGlobalCluster(p,map->nLocalDofWithoutGhost( p ));
    }
    map->setFirstDofGlobalCluster(0, 0 );
    map->setLastDofGlobalCluster(0, (map->nLocalDofWithoutGhost( 0 )>0)? map->nLocalDofWithoutGhost( 0 )-1 : 0 );
    for ( rank_type p=1;p<nProc;++p )
    {
        size_type nLocalDofWithoutGhost = map->nLocalDofWithoutGhost( p );
        map->setFirstDofGlobalCluster(p, (nLocalDofWithoutGhost>0)? map->lastDofGlobalCluster( p-1 )+1 : map->lastDofGlobalCluster( p-1 )  );
        map->setLastDofGlobalCluster(p, map->lastDofGlobalCluster( p-1 )+nLocalDofWithoutGhost );
    }

    size_type nDofCheck=0;
    for ( rank_type p=0;p<nProc;++p )
        nDofCheck += map->nLocalDofWithoutGhost( p );
    CHECK( nDof == nDofCheck ) << "invalid number of dof : " <<nDof << " vs " << nDofCheck;

    comm.barrier();
    if ( comm.isMasterRank() )
        std::cout << "Start benchmarking\n";
    tic();
    VectorUblas<double> vecWrite( map );
    vecWrite.setConstant( 3. );
    toc("init vector write");

#ifdef FEELPP_HAS_HDF5
    tic();
    vecWrite.saveHDF5( "myvec.h5" );
    toc("saveHDF5 vector");
#endif

    tic();
    VectorUblas<double> vecRead( map );
    toc("init vector read");

#ifdef FEELPP_HAS_HDF5
    tic();
    vecRead.loadHDF5( "myvec.h5" );
    toc("loadHDF5 vector");
#endif

    if ( comm.isMasterRank() )
        std::cout << "Finish benchmarking\n";

    return 0;
}
