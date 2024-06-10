#define BOOST_TEST_MODULE test_specx
#include <feel/feelcore/testsuite.hpp>


//#include <feel/feelmesh/ranges.hpp>

///nvme0/lemoinep/feelpp/feelpp/feel/feelcore/range.hpp

//#include <feel/feelcore/enumerate.hpp>
//#include <feel/feelcore/environment.hpp>
//#include <feel/feeldiscr/mesh.hpp>
//#include <feel/feelfilters/loadmesh.hpp>



#include <specx/Data/SpDataAccessMode.hpp>

#include <specx/Legacy/SpRuntime.hpp>

#include <specx/Task/SpPriority.hpp>
#include <specx/Task/SpProbability.hpp>

#include <specx/Utils/SpArrayView.hpp>
#include <specx/Utils/SpTimer.hpp>
#include <specx/Utils/small_vector.hpp>
#include <specx/Utils/SpBufferDataView.hpp>
#include <specx/Utils/SpBufferDataView.hpp>
#include <specx/Utils/SpHeapBuffer.hpp>
#include <specx/Utils/SpUtils.hpp>
#include <specx/Utils/SpConsumerThread.hpp>


#include <specx/Legacy/SpRuntime.hpp>



//#include "UTester.hpp"
//#include "utestUtils.hpp"
//#define USE_GPU_HIP


#include <feel/feeltask/taskHPC.hpp>
//#include "Taskflow_HPC.hpp"

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( testspecx_suite )



BOOST_AUTO_TEST_CASE( test_specx_21 )
{
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T21>\n");
    auto start_time= std::chrono::steady_clock::now();
    int nbThreads=9;
    long int nbN=1000000;
    int sizeBlock=nbN/nbThreads;
    int diffBlock=nbN-sizeBlock*nbThreads;
    double h=1.0/double(nbN);
    double integralValue=0.0;
    std::vector<double> valuesVec(nbThreads,0.0);
    std::vector<double> sumVec(nbThreads*10,0.0);

    auto FC1=[h,sizeBlock](const int& k,double& s) {  
            int vkBegin=k*sizeBlock;
            int vkEnd=(k+1)*sizeBlock;
            double sum=0.0; double x=0.0;
            for(int j=vkBegin;j<vkEnd;j++) { x=h*double(j); sum+=4.0/(1.0+x*x); }
            s=sum;
        return true;
    };

        LEM::Task TsK(nbThreads,3); 
        TsK.setSave(false); 
        TsK.setInfo(false); 
        TsK.setFileName("./PI");
        for(int k=0;k<nbThreads;k++) { 
            auto const& idk = k;
            TsK.add(_param(idk,valuesVec.at(idk)),_task=FC1);
        }

    TsK.run();
    TsK.close();
    integralValue=h*std::reduce(valuesVec.begin(),valuesVec.end());
    //std::cout<<"PI Value= "<<integralValue<<"\n";
    BOOST_MESSAGE("[INFO SPECX] : PI Value= "<<integralValue<<"\n");
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T21> in ms since start :"<<run_time.count()<<"\n");
}


BOOST_AUTO_TEST_CASE( test_specx_22 )
{
    //ADD 2 Vectors
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T22>\n");
    auto start_time= std::chrono::steady_clock::now();

    int nbThreads = 6;
    long int nbN=nbThreads*6;
    int sizeBlock=nbN/nbThreads;
    int diffBlock=nbN-sizeBlock*nbThreads;

    std::vector<double> VecA;
    std::vector<double> VecB;
    std::vector<double> VecR;

    for(int i=0;i<nbN;i++) {  VecA.push_back(i);   VecB.push_back(i);  }

    
    auto FC1=[VecA,VecB,sizeBlock,&VecR](const int& k) {  
            int vkBegin=k*sizeBlock;
            int vkEnd=(k+1)*sizeBlock;
            for(int j=vkBegin;j<vkEnd;j++)
            {
                VecR.push_back(VecA[j]+VecB[j]);    
            }
        return true;
    };

        LEM::Task TsK(nbThreads,3); 
        TsK.setSave(false); 
        TsK.setInfo(false); 
        TsK.setFileName("./AddVectors");
        TsK.getInformation();
        for(int k=0;k<nbThreads;k++) { 
            auto const& idk = k;
            TsK.add(_param(idk),_task=FC1);
        }

        TsK.run();        
        TsK.close();

    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T22> in ms since start :"<<run_time.count()<<"\n");
}


BOOST_AUTO_TEST_CASE( test_specx_23 )
{
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T23>\n");
    auto start_time= std::chrono::steady_clock::now();

    const int nbThreads = 6;
    LEM::Task TsK(nbThreads,3);
    TsK.setSave(false); 
    TsK.setInfo(false); 
    TsK.setFileName("./Test");
    TsK.getInformation();
    const int initVal1 = 100;
    const int initVal2 = 1000;
    int writeVal = 0;
    TsK.add(_param(initVal1,writeVal),
        _task=[](const int& initValParam, int& writeValParam){
            writeValParam += initValParam;
        }
    );
    TsK.run();   

    TsK.add(_param(initVal1,writeVal),
        _task=[](const int& initValParam, int& writeValParam){
            writeValParam += 123;
        }
    );
    TsK.run();   

    TsK.add(_param(initVal2,writeVal),
        _task=[](const int& initValParam, int& writeValParam){
            writeValParam += initValParam;
        }
    );
    TsK.run();   

    TsK.close();


    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T23> in ms since start :"<<run_time.count()<<"\n");
}



BOOST_AUTO_TEST_CASE( test_specx_24 )
{
    //Detach with Sleep
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T24>\n");
    auto start_time= std::chrono::steady_clock::now();

    const int nbThreads = 6;
    LEM::Task TsK(nbThreads,2);
    TsK.setSave(false); 
    TsK.setInfo(false); 
    TsK.setFileName("./Test");
    TsK.getInformation();
    const int initVal1 = 100;
    const int initVal2 = 1000;
    int writeVal = 0;
    TsK.add(_param(initVal1,writeVal),
        _task=[](const int& initValParam, int& writeValParam){
            std::cout <<"I live 1!"<<std::endl; sleep(1);  std::cout <<"YES 1!"<<std::endl;
        }
    );
    TsK.run();   

    TsK.setDetach(true); // <<=== detach mode

    TsK.add(_param(initVal1,writeVal),
        _task=[](const int& initValParam, int& writeValParam){
            std::cout <<"I live 2! (Detach)"<<std::endl; sleep(5);  std::cout <<"YES 2! (Detach)"<<std::endl;
        }
    );
    TsK.run();   

    TsK.add(_param(initVal2,writeVal),
        _task=[](const int& initValParam, int& writeValParam){
            std::cout <<"I live 3!"<<std::endl; sleep(1);  std::cout <<"YES 3!"<<std::endl;
        }
    );
    TsK.run();   
    TsK.close();

    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T24> in ms since start :"<<run_time.count()<<"\n");
}

BOOST_AUTO_TEST_CASE( test_specx_25 )
{
    //for_each
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T25>\n");
    auto start_time= std::chrono::steady_clock::now();

    const int nbThreads = 6;
    int valExternal=0;

     auto FC1=[&](const int &v) {  
        //std::cout <<"I live "<<v<<"\n";
        valExternal=valExternal+v;
        return true;
    };

    double valOutput1=0;
    int           idk=0;
    const int      nb=10;

    LEM::Task TsK(nbThreads,1);
    TsK.setSave(false); 
    TsK.setInfo(false); 
    TsK.setFileName("./ForEach");
    TsK.getInformation();

    std::vector<int> Vec(nb,0); for (int k=0; k<nb; ++k) { Vec[k]=k; }
    TsK.for_each(Vec.begin(),Vec.end(),_param(idk),_task=FC1);

    TsK.run();
    TsK.close();

    //std::cout<<">>> valExternal="<<valExternal<< "\n";
    BOOST_MESSAGE("[INFO SPECX] : Value= "<<valExternal<<"\n"); 
    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T25> in ms since start :"<<run_time.count()<<"\n");
}


//#define EIGEN_DEVICE_FUNC __device__ __host__

template<typename T>
struct matrix_inverse {
  EIGEN_DEVICE_FUNC
  void operator()(int i, const typename T::Scalar* in, typename T::Scalar* out) const
  {
    using namespace Eigen;
    T M(in+i);
    Map<T> res(out+i*T::MaxSizeAtCompileTime);
    res = M.inverse();
  }
};

BOOST_AUTO_TEST_CASE( test_specx_26 )
{
    //for_each
    BOOST_MESSAGE("[INFO SPECX] : Execution of <T26>\n");
    auto start_time= std::chrono::steady_clock::now();

    #ifdef USE_GPU_HIP
        int nthreads = 1;
        Eigen::VectorXf in, out;
        int data_size = nthreads * 9;
        
        in.setRandom(data_size);
        in.setConstant(data_size,0);

        in(0)=1;
        in(4)=2;
        in(8)=1;
        out.setConstant(data_size, -1);
        cfin.setRandom(data_size);
        cfout.setConstant(data_size, -1);
        LEM::Task TsK(1,0);
        TsK.setSave(false); 
        TsK.setInfo(false); 
        TsK.qViewChrono=true;
        TsK.run_gpu_1D(matrix_inverse<Eigen::Matrix3f>(),dim3(128,1,1), nthreads, in, out);
        TsK.run();
        TsK.close();
        TsK.debriefingTasks();
    #endif

    auto stop_time= std::chrono::steady_clock::now();
    auto run_time=std::chrono::duration_cast<std::chrono::microseconds> (stop_time-start_time);
	BOOST_MESSAGE("[INFO SPECX] : Execution Time <T26> in ms since start :"<<run_time.count()<<"\n");
}








BOOST_AUTO_TEST_SUITE_END()



