/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-06-18

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2011-present Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file opusapp.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Cecile Daversin <daversin@math.unistra.fr>
   \author Stephane Veys
   \date 2011-06-18
 */
#ifndef __OpusApp_H
#define __OpusApp_H 1

#include <feel/feelmor/options.hpp>
#include <feel/feelmor/crb.hpp>
#include <feel/feelmor/eim.hpp>
#include <feel/feelmor/ser.hpp>
#include <feel/feelmor/crbmodel.hpp>
#include <boost/serialization/version.hpp>
#include <boost/range/join.hpp>
#include <boost/regex.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>


namespace Feel
{
po::options_description opusapp_options( std::string const& prefix );
std::string _o( std::string const& prefix, std::string const& opt );

enum class SamplingMode
{
    RANDOM = 0, EQUIDISTRIBUTED = 1, LOGEQUIDISTRIBUTED = 2, READFROMCOMMANDLINE = 3
};

#define oprec 4
#define Pdim 7
#define ofill ' '
#define dmanip std::scientific << std::setprecision( oprec )
#define hdrmanip(N) std::setw(N) << std::setfill(ofill) << std::right
#define tabmanip(N) std::setw(N) << std::setfill(ofill) << std::right << dmanip


/**
 * \class OpusApp
 * \brief Certified Reduced Basis application
 *
 * @author Christophe Prud'homme
 */
template<typename ModelType,
         template < typename ReducedMethod > class RM=CRB,
         template < typename ModelInterface > class Model=CRBModel>
class OpusApp   : public Application
{
    typedef Application super;
public:


    //! model type
    typedef ModelType model_type;
    typedef std::shared_ptr<ModelType> model_ptrtype;

    //! function space type
    typedef typename model_type::functionspace_type functionspace_type;
    typedef typename model_type::functionspace_ptrtype functionspace_ptrtype;

    typedef typename model_type::element_type element_type;

    typedef Eigen::VectorXd vectorN_type;

#if 0
    //old
    typedef CRBModel<ModelType> crbmodel_type;
    typedef std::shared_ptr<crbmodel_type> crbmodel_ptrtype;
    typedef CRB<crbmodel_type> crb_type;
    typedef std::shared_ptr<crb_type> crb_ptrtype;
#endif

    typedef Model<ModelType> crbmodel_type;
    typedef std::shared_ptr<crbmodel_type> crbmodel_ptrtype;
    typedef RM<crbmodel_type> crb_type;
    typedef std::shared_ptr<crb_type> crb_ptrtype;

    typedef SER<crb_type> ser_type;
    typedef std::shared_ptr<ser_type> ser_ptrtype;

    typedef CRBModel<ModelType> crbmodelbilinear_type;

    typedef typename ModelType::parameter_type parameter_type;
    typedef typename ModelType::mesh_type mesh_type;
    typedef std::vector< parameter_type > vector_parameter_type;

    typedef typename crb_type::sampling_ptrtype sampling_ptrtype;

    static const int nb_spaces = functionspace_type::nSpaces;

    typedef typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<2> > , fusion::vector< mpl::int_<0>, mpl::int_<1> >  ,
                               typename mpl::if_ < boost::is_same< mpl::int_<nb_spaces> , mpl::int_<3> > ,
                                                   fusion::vector < mpl::int_<0> , mpl::int_<1> , mpl::int_<2> >,
                                                   typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<4> >,
                                                                      fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3> >,
                                                                      fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3>, mpl::int_<4> >
                                                                      >::type >::type >::type index_vector_type;

    OpusApp()
        :
        super(),
        M_modelName( soption(_prefix=this->about().appName(),_name="model-name") ),
        M_mode( ( CRBModelMode )ioption(_name=_o( this->about().appName(),"run.mode" )) ),
        use_newton_( boption(_name="crb.use-newton") && !ModelType::is_linear )
        {
            this->init();
        }

    OpusApp( AboutData const& ad, po::options_description const& od )
        :
        OpusApp(ad, od,  ( CRBModelMode )ioption(_name=_o( this->about().appName(),"run.mode" )) )
        {
        }
    OpusApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        OpusApp(ad, od,  ( CRBModelMode )ioption(_name=_o( this->about().appName(),"run.mode" )) )
        {
        }
    OpusApp( int argc, char** argv, AboutData const& ad, po::options_description const& od, CRBModelMode mode )
        :
        OpusApp(ad,od,mode)
        {}
    OpusApp( AboutData const& ad, po::options_description const& od, CRBModelMode mode )
        :
        super( ad, opusapp_options( ad.appName() )
               .add( od )
               .add( crbOptions() )
               .add( eimOptions() )
               .add( crbSEROptions() )
               .add( podOptions() )),
        M_modelName( soption(_prefix=this->about().appName(),_name="model-name") ),
        M_mode( mode ),
        use_newton_( boption(_name="crb.use-newton") && !ModelType::is_linear )
        {
            this->init();
        }

private:
    void init()
        {
            try
            {
                M_current_path = fs::current_path();
                LOG(INFO) << "[OpusApp] mode:" << ( int )M_mode << "\n";

                crb = newCRB();

                LOG(INFO) << "[OpusApp] get crb done" << "\n";

                //VLOG(1) << "[OpusApp] get crb done" << "\n";
                //crb->setTruthModel( model );
                //VLOG(1) << "[OpusApp] constructor done" << "\n";
            }

            catch ( boost::bad_any_cast const& e )
            {
                std::cout << "[OpusApp] a bad any cast occured, probably a nonexistant or invalid  command line/ config options\n";
                std::cout << "[OpusApp] exception reason: " << e.what() << "\n";
            }

        }
public:
    /**
     * \return a new CRB shared pointer
     */
    crb_ptrtype newCRB( int level=0 )
        {
            model = std::make_shared<crbmodel_type>(M_modelName,crb::stage::offline,level);
            return crb_type::New(model->model()->modelName(), model, crb::stage::offline );
        }
    crb_ptrtype & crbPtr() { return crb; }
    crb_ptrtype const& crbPtr() const { return crb; }

    /* Get parameter space associated to model */
    auto getParameterSpace() const { return model->parameterSpace(); }
    
    /* Returns CRB objects */
    crb_ptrtype getCRB() const { return this->crb; }

    void setMode( std::string const& mode )
        {
            if ( mode == "pfem" ) M_mode = CRBModelMode::PFEM;

            if ( mode == "crb" ) M_mode = CRBModelMode::CRB;

            if ( mode == "scm" ) M_mode = CRBModelMode::SCM;

            if ( mode == "scm_online" ) M_mode = CRBModelMode::SCM_ONLINE;

            if ( mode == "crb_online" ) M_mode = CRBModelMode::CRB_ONLINE;
        }
    void setMode( CRBModelMode mode )
        {
            M_mode = mode;
        }

    void loadDB()
        {
            int proc_number = Environment::worldComm().globalRank();
            int global_size = Environment::worldComm().globalSize();
            std::string pslogfile = ( boost::format("PsLogCrbOffline-%1%_%2%") %global_size %proc_number ).str();

            bool only_master=boption(_name="crb.system-memory-evolution");
            bool all_procs  =boption(_name="crb.system-memory-evolution-on-all-procs");
            bool only_one_proc=only_master * ( Environment::worldComm().globalRank()==Environment::worldComm().masterRank() );
            bool write_memory_evolution = all_procs || only_one_proc ;

            bool crb_use_predefined = boption(_name="crb.use-predefined-WNmu");
            std::string file_name;
            int NlogEquidistributed = ioption(_name="crb.use-logEquidistributed-WNmu");
            int Nequidistributed = ioption(_name="crb.use-equidistributed-WNmu");
            int NlogEquidistributedScm = ioption(_name="crb.scm.use-logEquidistributed-C");
            int NequidistributedScm = ioption(_name="crb.scm.use-equidistributed-C");
            typename crb_type::sampling_ptrtype Sampling( new typename crb_type::sampling_type( model->parameterSpace() ) );
            bool all_proc_have_same_sampling=true;

            if( crb_use_predefined )
            {
                file_name = ( boost::format("SamplingWNmu") ).str();
            }
            if( NlogEquidistributed+Nequidistributed > 0 )
            {
                file_name = ( boost::format("SamplingWNmu") ).str();
                if( NlogEquidistributed > 0 )
                    Sampling->logEquidistribute( NlogEquidistributed,all_proc_have_same_sampling );
                if( Nequidistributed > 0 )
                    Sampling->equidistribute( Nequidistributed,all_proc_have_same_sampling );
                Sampling->writeOnFile(file_name);
            }
            if( NlogEquidistributedScm+NequidistributedScm > 0 )
            {
                file_name = ( boost::format("SamplingC") ).str();
                if( NlogEquidistributedScm > 0 )
                    Sampling->logEquidistribute( NlogEquidistributedScm,all_proc_have_same_sampling );
                if( NequidistributedScm > 0 )
                    Sampling->equidistribute( NequidistributedScm,all_proc_have_same_sampling );
                Sampling->writeOnFile(file_name);
            }

            if ( M_mode == CRBModelMode::PFEM )
                return;

            if ( !crb->scm()->isDBLoaded() || crb->scm()->rebuildDB() )
            {
                if ( M_mode == CRBModelMode::SCM )
                {
                    if( proc_number == Environment::worldComm().masterRank() )
                        std::cout << "No SCM DB available, do scm offline computations first...\n";
                    if( crb->scm()->doScmForMassMatrix() )
                        crb->scm()->setScmForMassMatrix( true );

                    crb->scm()->offline();
                }
            }

            if ( !crb->isDBLoaded() || crb->rebuild() )
            {
                if ( M_mode == CRBModelMode::CRB )
                    //|| M_mode == CRBModelMode::SCM )
                {
                    if( proc_number == Environment::worldComm().masterRank() && !ioption(_name="ser.rb-frequency") )
                        std::cout << "No CRB DB available, do crb offline computations...\n";
                    crb->setOfflineStep( true );
                    do  // SER r-adaptation for RB
                    {
                        crb->setAdaptationSER( false ); //re-init to false
                        crb->offline();
                    }
                    while(crb->adaptationSER());

                    if( write_memory_evolution )
                        this->generateMemoryEvolution(pslogfile);
                }

                else if ( M_mode != CRBModelMode::SCM )
                    throw std::logic_error( "CRB/SCM Database could not be loaded" );
            }

            //if( crb->isDBLoaded() )
            else
            {
                int Nrestart = ioption(_name="crb.restart-from-N");
                bool do_offline = false;
                int current_dimension = crb->dimension();
                int dimension_max = ioption(_name="crb.dimension-max");
                int sampling_size = 0;
                if( crb_use_predefined )
                    sampling_size = Sampling->readFromFile(file_name);

                if( sampling_size > current_dimension )
                    do_offline = true;

                if( current_dimension < dimension_max && !crb_use_predefined )
                    do_offline=true;

                if( Nrestart > 1 )
                    do_offline=true;

                if( ! do_offline )
                {
                    crb->loadSCMDB();
                    crb->setOfflineStep( false );
                }

                if( do_offline )
                {
                    crb->setOfflineStep( true );
                    crb->offline();
                    if( write_memory_evolution )
                        this->generateMemoryEvolution(pslogfile);
                }
            }
        }

    element_type paraFeelRun( std::map<std::string,double> mu_map, int N=-1)
        {
            vectorN_type time_crb;
            double online_tol = doption(_name="crb.online-tolerance");
            bool print_rb_matrix = boption(_name="crb.print-rb-matrix");
            parameter_type mu;
            for( int i=0; i<mu.size(); i++)
            {
                mu(i) = mu_map[crb->Dmu->parameterName(i)];
            }
            auto o = crb->run( mu, time_crb, online_tol, N, print_rb_matrix);
            auto solutions = o.template get<2>();
            auto uN = solutions.template get<0>();
            auto WN = crb->wn();
            auto u_crb = crb->expansion( uN[uN.size()-1], N, false );
            return u_crb;
        }


    FEELPP_DONT_INLINE void run();

    void run( const double * X, unsigned long N,
              double * Y, unsigned long P )
        {

            switch ( M_mode )
            {
            case  CRBModelMode::PFEM:
            {
                model->run( X, N, Y, P );
            }
            break;

            case  CRBModelMode::SCM:
            case  CRBModelMode::SCM_ONLINE:
            {
                crb->scm()->run( X, N, Y, P );
            }
            break;

            case  CRBModelMode::CRB:
            case  CRBModelMode::CRB_ONLINE:
            {
                crb->run( X, N, Y, P );
            }
            break;
            }

            fs::current_path( M_current_path );
        }


    /**
     * @brief Compute the FEM solution for a given parameter \p mu
     * 
     * @param mu parameter
     * @param use_newton use Newton method (default true)
     * @return element_type solution of the FEM problem
     */
    element_type getFEMsolution( parameter_type const& mu, bool use_newton=true )
    {
        element_type u_pfem;
        if( use_newton )
            u_pfem =  model->solveFemUsingAffineDecompositionNewton( mu );
        else
            u_pfem =  model->solveFemUsingAffineDecompositionFixedPoint( mu );
        return u_pfem;
    }

    /**
     * @brief Compute the RB solution for a given parameter mu
     * 
     * @param mu parameter
     * @param N size of the reduced basis (default -1, i.e. use the maximum size)
     * @return auto tuple composed of uN, output, errorBound
     */
    auto getRBsolution( parameter_type const &mu, int N = -1 )
    {
        vectorN_type time_crb;
        double online_tol = doption(_name="crb.online-tolerance");
        bool print_rb_matrix = boption(_name="crb.print-rb-matrix");

        auto o = crb->run( mu, time_crb, online_tol, N, print_rb_matrix);
        auto uN = o.coefficients();
        double errorBound = o.errorbound();
        double output = o.output();

        return std::make_tuple( uN, output, errorBound );
    }

    /**
     * @brief Compute the effectivity of the RB solution for a given parameter \p mu
     * 
     * @param mu parameter
     * @param N size of the reduced basis (default -1, i.e. use the maximum size)
     * @return double effectivity $\eta_N^s(\mu) = \frac{\Delta_N^s(\mu)}{s(\mu) - s_N(\mu)}$
     */
    double computeEffectivity( parameter_type const &mu, int N = -1 )
    {
        element_type u_pfem = getFEMsolution( mu, use_newton_ );
        auto sol_rbm = getRBsolution( mu, N );
        double output_crb = std::get<1>(sol_rbm);
        double output_fem = model->output( 1, mu, u_pfem, false );
        double error_bound = std::get<2>(sol_rbm);

        return error_bound / math::abs( output_crb - output_fem );        
    }


private:
    int printParameterHdr( std::ostream& os, int N, std::vector<std::string> outputhdrs )
        {
            for ( int i = 0; i < N; ++i )
            {
                std::ostringstream s;
                s << "mu" << i;
                os  << hdrmanip( oprec+7 ) << s.str();
            }

            for( auto output : outputhdrs )
            {
                os << hdrmanip( 15 ) << output;
            }
            os << "\n";

            return N*( oprec+7 )+outputhdrs.size()*15;
        }
    void printEntry( std::ostream& os,
                     typename ModelType::parameter_type const& mu,
                     std::vector<double> const& outputs )
        {
            for ( int i = 0; i < mu.size(); ++i )
                os  << std::right <<std::setw( oprec+7 ) << dmanip << mu[i];

            for( auto o : outputs )
            {
                os << tabmanip( 15 ) << o;
            }
            os << "\n";
        }


    /*
     * read the file generated by PsLogger
     * and generate a file to plot memory evolution in function of number of elements
     * ( during offline step )
     */
    void generateMemoryEvolution(std::string pslogname)
    {
        std::ifstream file_read (pslogname);
        std::vector< double > memory;
        std::vector< double > percent;

        // number of elements in the RB
        int N = ioption("crb.dimension-max");
        //int N = crb->dimension();
        memory.resize( N );
        percent.resize( N );
        if( file_read )
        {
            double mem;
            double per;
            std::string str;
            //first line
            for(int i=0; i<5; i++)
                file_read >> str;
            //second line
            for(int i=0; i<3; i++)
                file_read >> str;
            //for each elements in the RB
            for(int i=0; i<N; i++)
            {
                file_read >> str ;
                file_read >> mem;
                file_read >> per;
                file_read >> str;
                memory[i] = mem;
                percent[i] = per;
            }
            file_read.close();
        }//if file
        else
        {
            throw std::logic_error( "[OpusApp::generateMemoryEvaluation] ERROR loading the file generated by ps-log (does not exist)" );
        }//no file

        int proc_number = Environment::worldComm().globalRank();
        int global_size = Environment::worldComm().globalSize();
        std::string file_name = ( boost::format("MemoryEvolution-%1%_%2%") %global_size %proc_number ).str();
        std::ofstream file_write;
        file_write.open( file_name,std::ios::out );
        file_write << "N \t Memory \t Percent \n";
        for(int i=0; i<N; i++)
            file_write << i << "\t"<< memory[i]<< "\t" << percent[i] << "\n";
        file_write.close();
    }

    //double l2Norm( typename ModelType::parameter_type const& mu, int N )
    double l2Norm( element_type const& u )
    {
        static const bool is_composite = functionspace_type::is_composite;
        return l2Norm( u, mpl::bool_< is_composite >() );
    }
    double l2Norm( element_type const& u, mpl::bool_<false> )
    {
        return math::sqrt( integrate( _range=u.functionSpace()->template rangeElements<0>(), _expr=(vf::idv(u))*(vf::idv(u)) ).evaluate()(0,0) );
    }
    double l2Norm( element_type const& u, mpl::bool_<true>)
    {
        index_vector_type index_vector;
        return fusion::fold( index_vector, double(0.), ComputeNormL2InCompositeCase( u ) );
    }
    //double h1Norm( typename ModelType::parameter_type const& mu, int N )
    double h1Norm( element_type const& u )
    {
        static const bool is_composite = functionspace_type::is_composite;
        return h1Norm( u, mpl::bool_< is_composite >() );
    }
    double h1Norm( element_type const& u, mpl::bool_<false> )
    {
        auto mesh = model->functionSpace()->mesh();
        double l22 = integrate( _range=u.functionSpace()->template rangeElements<0>(), _expr=(vf::idv(u))*(vf::idv(u)) ).evaluate()(0,0);
        double semih12 = integrate( _range=u.functionSpace()->template rangeElements<0>(), _expr=(vf::gradv(u))*trans(vf::gradv(u)) ).evaluate()(0,0);
        return math::sqrt( l22+semih12 );
    }
    double h1Norm( element_type const& u, mpl::bool_<true>)
    {
        auto mesh = model->functionSpace()->mesh();
        auto u_femT = u.template element<1>();
        double l22 = integrate( _range=u.functionSpace()->template rangeElements<0>(), _expr=(vf::idv(u_femT))*(vf::idv(u_femT)) ).evaluate()(0,0);
        double semih12 = integrate( _range=u_femT.functionSpace()->template rangeElements<0>(), _expr=(vf::gradv(u_femT))*trans(vf::gradv(u_femT))).evaluate()(0,0);
        return math::sqrt( l22+semih12 );
    }

    struct ComputeNormL2InCompositeCase
    {
        typedef double result_type;

        ComputeNormL2InCompositeCase( element_type const& composite_u )
            :
            M_composite_u( composite_u )
            {}

        template< typename T >
        result_type
        operator()( result_type const& r,const T& t ) const
            {
                int i = T::value;
                if( i == 0 )
                    M_vec.resize( 1 );
                else
                    M_vec.conservativeResize(i+1);

                auto u = M_composite_u.template element< T::value >();
                auto mesh = u.functionSpace()->mesh();

                double norm  = normL2(_range=u.functionSpace()->template rangeElements<0>(),_expr=( idv(u) ) );
                M_vec(i) = norm ;
                return norm + r;
            }

        double norm()
            {
                return M_vec.sum();
            }

    private :
        mutable vectorN_type M_vec;
        element_type const& M_composite_u;
    };


    void initializeConvergenceScmMap( int sampling_size )
    {
        auto N = crb->scm()->KMax();
        M_mapConvSCM["RelativeError"] = std::vector<vectorN_type>(N);

        for(int j=0; j<N; j++)
        {
            M_mapConvSCM["RelativeError"][j].resize(sampling_size);
        }
    }

    void studyEimConvergence( typename ModelType::parameter_type const& mu , element_type & model_solution , int mu_number)
    {
        auto eim_sc_vector = model->scalarContinuousEim();
        auto eim_sd_vector = model->scalarDiscontinuousEim();

        for(int i=0; i<eim_sc_vector.size(); i++)
        {
            auto eim = eim_sc_vector[i];
            std::vector< std::string > all_file_name;
            all_file_name.push_back( eim->name()+"-EimConvergenceL2.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceL2estimated.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceL2ratio.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINF.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINFestimated.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINFratio.dat");
            eim->studyConvergence( mu , model_solution , all_file_name );
        }

        for(int i=0; i<eim_sd_vector.size(); i++)
        {
            std::vector< std::string > all_file_name;
            auto eim = eim_sd_vector[i];
            all_file_name.push_back( eim->name()+"-EimConvergenceL2.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceL2estimated.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceL2ratio.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINF.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINFestimated.dat");
            all_file_name.push_back( eim->name()+"-EimConvergenceLINFratio.dat");
            eim->studyConvergence( mu , model_solution , all_file_name );
        }

        //std::vector<double> l2error;
        //for(int j=0; j<l2error.size(); j++)
        //    M_mapConvEIM[eim->name()][j][mu_number-1] = l2error[j];

    }

    void doTheEimConvergenceStat( int sampling_size )
    {
        auto eim_sc_vector = model->scalarContinuousEim();
        auto eim_sd_vector = model->scalarDiscontinuousEim();

        for(int i=0; i<eim_sc_vector.size(); i++)
        {
            auto eim = eim_sc_vector[i];

            int max=eim->mMax();

            if( max > 1 )
            {
                //load information contained in files
                std::vector< vectorN_type > L2, L2estimated, L2ratio, LINF, LINFestimated, LINFratio;
                std::string filename = eim->name()+"-EimConvergenceL2.dat";
                model->readConvergenceDataFromFile( L2, filename );
                filename = eim->name()+"-EimConvergenceL2estimated.dat";
                model->readConvergenceDataFromFile( L2estimated, filename );
                filename = eim->name()+"-EimConvergenceL2ratio.dat";
                model->readConvergenceDataFromFile( L2ratio, filename );
                filename = eim->name()+"-EimConvergenceLINF.dat";
                model->readConvergenceDataFromFile( LINF, filename );
                filename = eim->name()+"-EimConvergenceLINFestimated.dat";
                model->readConvergenceDataFromFile( LINFestimated, filename );
                filename = eim->name()+"-EimConvergenceLINFratio.dat";
                model->readConvergenceDataFromFile( LINFratio, filename );

                //write files containing statistics
                filename = "cvg-eim-"+eim->name()+"-L2.dat";
                model->writeConvergenceStatistics( L2 , filename);
                filename = "cvg-eim-"+eim->name()+"-L2estimated.dat";
                model->writeConvergenceStatistics( L2estimated , filename);
                filename = "cvg-eim-"+eim->name()+"-L2ratio.dat";
                model->writeConvergenceStatistics( L2ratio , filename);
                filename = "cvg-eim-"+eim->name()+"-LINF.dat";
                model->writeConvergenceStatistics( LINF , filename);
                filename = "cvg-eim-"+eim->name()+"-LINFestimated.dat";
                model->writeConvergenceStatistics( LINFestimated , filename);
                filename = "cvg-eim-"+eim->name()+"-LINFratio.dat";
                model->writeConvergenceStatistics( LINFratio , filename);

            }
        }

        for(int i=0; i<eim_sd_vector.size(); i++)
        {
            auto eim = eim_sd_vector[i];

            int max=eim->mMax();

            if( max > 1 )
            {
                //load information contained in files
                std::vector< vectorN_type > L2, L2estimated, L2ratio, LINF, LINFestimated, LINFratio;
                std::string filename = eim->name()+"-EimConvergenceL2.dat";
                model->readConvergenceDataFromFile( L2, filename );
                filename = eim->name()+"-EimConvergenceL2estimated.dat";
                model->readConvergenceDataFromFile( L2estimated, filename );
                filename = eim->name()+"-EimConvergenceL2ratio.dat";
                model->readConvergenceDataFromFile( L2ratio, filename );
                filename = eim->name()+"-EimConvergenceLINF.dat";
                model->readConvergenceDataFromFile( LINF, filename );
                filename = eim->name()+"-EimConvergenceLINFestimated.dat";
                model->readConvergenceDataFromFile( LINFestimated, filename );
                filename = eim->name()+"-EimConvergenceLINFratio.dat";
                model->readConvergenceDataFromFile( LINFratio, filename );

                //write files containing statistics
                filename = "cvg-eim-"+eim->name()+"-L2.dat";
                model->writeConvergenceStatistics( L2 , filename);
                filename = "cvg-eim-"+eim->name()+"-L2estimated.dat";
                model->writeConvergenceStatistics( L2estimated , filename);
                filename = "cvg-eim-"+eim->name()+"-L2ratio.dat";
                model->writeConvergenceStatistics( L2ratio , filename);
                filename = "cvg-eim-"+eim->name()+"-LINF.dat";
                model->writeConvergenceStatistics( LINF , filename);
                filename = "cvg-eim-"+eim->name()+"-LINFestimated.dat";
                model->writeConvergenceStatistics( LINFestimated , filename);
                filename = "cvg-eim-"+eim->name()+"-LINFratio.dat";
                model->writeConvergenceStatistics( LINFratio , filename);
            }
        }
    }

    void doTheCrbConvergenceStat( int sampling_size )
    {
        int N = ioption("crb.dimension-max");

        std::vector< vectorN_type > L2, H1, OutputError, OutputErrorEstimated, OutputErrorBoundEfficiency, AbsoluteOutputErrorBoundEfficiency;
        std::vector< vectorN_type > PrimalSolutionError, PrimalSolutionErrorEstimated, PrimalSolutionErrorBoundEfficiency , PrimalAbsoluteSolutionErrorBoundEfficiency;
        std::vector< vectorN_type > DualSolutionError, DualSolutionErrorEstimated, DualSolutionErrorBoundEfficiency, DualAbsoluteSolutionErrorBoundEfficiency;
        std::vector< vectorN_type > CrbTimePrediction, CrbTimeErrorEstimation, FemTime;

        //load information contained in files
        std::string filename = "CrbConvergenceL2.dat";
        model->readConvergenceDataFromFile( L2, filename );
        filename = "CrbConvergenceH1.dat";
        model->readConvergenceDataFromFile( H1, filename );
        filename = "CrbConvergenceOutputError.dat";
        model->readConvergenceDataFromFile( OutputError, filename );
        filename = "CrbConvergenceOutputErrorEstimated.dat";
        model->readConvergenceDataFromFile( OutputErrorEstimated, filename );
        filename = "CrbConvergenceOutputErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( OutputErrorBoundEfficiency, filename );
        filename = "CrbConvergenceAbsoluteOutputErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( AbsoluteOutputErrorBoundEfficiency, filename );
        filename = "CrbConvergencePrimalSolutionError.dat";
        model->readConvergenceDataFromFile( PrimalSolutionError , filename );
        filename = "CrbConvergencePrimalSolutionErrorEstimated.dat";
        model->readConvergenceDataFromFile( PrimalSolutionErrorEstimated , filename );
        filename = "CrbConvergencePrimalSolutionErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( PrimalSolutionErrorBoundEfficiency , filename );
        filename = "CrbConvergencePrimalAbsoluteSolutionErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( PrimalAbsoluteSolutionErrorBoundEfficiency , filename );
        filename = "CrbConvergenceDualSolutionError.dat";
        model->readConvergenceDataFromFile( DualSolutionError , filename );
        filename = "CrbConvergenceDualSolutionErrorEstimated.dat";
        model->readConvergenceDataFromFile( DualSolutionErrorEstimated , filename );
        filename = "CrbConvergenceDualSolutionErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( DualSolutionErrorBoundEfficiency , filename );
        filename = "CrbConvergenceDualAbsoluteSolutionErrorBoundEfficiency.dat";
        model->readConvergenceDataFromFile( DualAbsoluteSolutionErrorBoundEfficiency , filename );
        filename = "CrbConvergenceCrbTimePrediction.dat";
        model->readConvergenceDataFromFile( CrbTimePrediction , filename );
        filename = "CrbConvergenceCrbTimeErrorEstimation.dat";
        model->readConvergenceDataFromFile( CrbTimeErrorEstimation , filename );
        filename = "CrbConvergenceFemTime.dat";
        model->readConvergenceDataFromFile( FemTime , filename );

        //write files containing statistics
        filename = "cvg-crb-L2.dat";
        model->writeConvergenceStatistics( L2 , filename);
        filename = "cvg-crb-H1.dat";
        model->writeConvergenceStatistics( H1 , filename);
        filename = "cvg-crb-OutputError.dat";
        model->writeConvergenceStatistics( OutputError , filename);
        filename = "cvg-crb-OutputErrorEstimated.dat";
        model->writeConvergenceStatistics( OutputErrorEstimated , filename);
        filename = "cvg-crb-OutputErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( OutputErrorBoundEfficiency , filename);
        filename = "cvg-crb-AbsoluteOutputErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( AbsoluteOutputErrorBoundEfficiency , filename);
        filename = "cvg-crb-PrimalSolutionError.dat";
        model->writeConvergenceStatistics( PrimalSolutionError , filename);
        filename = "cvg-crb-PrimalSolutionErrorEstimated.dat";
        model->writeConvergenceStatistics( PrimalSolutionErrorEstimated , filename);
        filename = "cvg-crb-PrimalSolutionErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( PrimalSolutionErrorBoundEfficiency , filename);
        filename = "cvg-crb-PrimalAbsoluteSolutionErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( PrimalAbsoluteSolutionErrorBoundEfficiency , filename);
        filename = "cvg-crb-DualSolutionError.dat";
        model->writeConvergenceStatistics( DualSolutionError , filename);
        filename = "cvg-crb-DualSolutionErrorEstimated.dat";
        model->writeConvergenceStatistics( DualSolutionErrorEstimated , filename);
        filename = "cvg-crb-DualSolutionErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( DualSolutionErrorBoundEfficiency , filename);
        filename = "cvg-crb-DualAbsoluteSolutionErrorBoundEfficiency.dat";
        model->writeConvergenceStatistics( DualAbsoluteSolutionErrorBoundEfficiency , filename);
        filename = "cvg-crb-CrbTimePrediction.dat";
        model->writeConvergenceStatistics( CrbTimePrediction , filename, "totaltime");
        filename = "cvg-crb-CrbTimeErrorEstimation.dat";
        model->writeConvergenceStatistics( CrbTimeErrorEstimation , filename, "totaltime");
        filename = "cvg-crb-FemTime.dat";
        model->writeConvergenceStatistics( FemTime , filename , "totaltime" );

        //now we have error and estimated error statistics for all parameters (min max ect ... )
        //but it is interesting to have also the efficiency associated to parameter that generated the max/min
        //error on the output or on primal/dual solution
        //so that what we do here
        filename = "cvg-crb-OutputErrorBoundEfficiency2.dat";
        model->writeVectorsExtremumsRatio( OutputError, OutputErrorEstimated, filename );
        filename = "cvg-crb-PrimalSolutionErrorBoundEfficiency2.dat";
        model->writeVectorsExtremumsRatio( PrimalSolutionError, PrimalSolutionErrorEstimated, filename );
        filename = "cvg-crb-DualSolutionErrorBoundEfficiency2.dat";
        model->writeVectorsExtremumsRatio( DualSolutionError, DualSolutionErrorEstimated, filename );

#if 0
        int Nmax = vector_sampling_for_primal_efficiency_under_1.size();
        for(int N=0; N<Nmax; N++)
        {
            std::string file_name_primal = (boost::format("Sampling-Primal-Problem-Bad-Efficiency-N=%1%") %(N+1) ).str();
            if( vector_sampling_for_primal_efficiency_under_1[N]->size() > 1 )
                vector_sampling_for_primal_efficiency_under_1[N]->writeOnFile(file_name_primal);
        }
        Nmax = vector_sampling_for_dual_efficiency_under_1.size();
        for(int N=0; N<Nmax; N++)
        {
            std::string file_name_dual = (boost::format("Sampling-Dual-Problem-Bad-Efficiency-N=%1%") %(N+1) ).str();
            if( vector_sampling_for_dual_efficiency_under_1[N]->size() > 1 )
                vector_sampling_for_dual_efficiency_under_1[N]->writeOnFile(file_name_dual);
        }
#endif
    }

    void doTheScmConvergenceStat( int sampling_size )
    {
        auto N = crb->scm()->KMax();
        std::list<std::string> list_error_type{"RelativeError"};
        for( auto error_name : list_error_type)
        {
            std::ofstream conv;
            std::string file_name = "cvg-scm-"+ error_name +"-stats.dat";

            if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
            {
                conv.open(file_name, std::ios::app);
                conv << "Nb_basis" << "\t" << "Min" << "\t" << "Max" << "\t" << "Mean" << "\t" << "Variance" << "\n";
            }

            for(int j=0; j<N; j++)
            {
                double mean = M_mapConvSCM[error_name][j].mean();
                double variance = 0.0;
                for( int k=0; k < sampling_size; k++)
                    variance += (M_mapConvSCM[error_name][j](k) - mean)*(M_mapConvSCM[error_name][j](k) - mean)/sampling_size;

                if( Environment::worldComm().globalRank()  == Environment::worldComm().masterRank() )
                {
                    conv << j+1 << "\t"
                         << M_mapConvSCM[error_name][j].minCoeff() << "\t"
                         << M_mapConvSCM[error_name][j].maxCoeff() << "\t"
                         << mean << "\t" << variance << "\n";
                }
            }
            conv.close();
        }
    }


    // Script write current mu in cfg => need to write it in SamplingForTest
    void buildSamplingFromCfg()
    {
        // Size of mu
        int mu_size = model->parameterSpace()->dimension();

        // Clear SamplingForTest is exists, and open a new one
        fs::path input_file ("SamplingForTest");
        if( fs::exists(input_file) )
            std::remove( "SamplingForTest" );

        std::ofstream input( "SamplingForTest" );
        input << "mu= [ ";

        // Check if cfg file is readable
        std::ifstream cfg_file( soption(_name="config-file") );
        if(!cfg_file)
            std::cout << "[Script-mode] Config file cannot be read" << std::endl;

        // OT writes values of mu in config file => read it and copy in SamplingForTest with specific syntax
        for(int i=1; i<=mu_size; i++)
            {
                // convert i into string
                std::ostringstream oss;
                oss << i;
                std::string is = oss.str();

                // Read cfg file, collect line with current mu_i
                std::string cfg_line_mu, tmp_content;
                std::ifstream cfg_file( soption(_name="config-file") );
                while(cfg_file)
                    {
                        std::getline(cfg_file, tmp_content);
                        if(tmp_content.compare(0,2+is.size(),"mu"+is) == 0)
                            cfg_line_mu += tmp_content;
                    }

                //Regular expression : corresponds to one set in xml file (mu<i>=<value>)
                std::string expr_s = "mu"+is+"[[:space:]]*=[[:space:]]*([0-9]+(.?)[0-9]*(e(\\+|-)[0-9]+)?)[[:space:]]*";
                boost::regex expression( expr_s );

                //Match mu<i>=<value> in cfg file and copy to SamplingForTest
                boost::smatch what;
                auto is_match = boost::regex_match(cfg_line_mu, what, expression);
                //std::cout << "is match ?" << is_match << std::endl;
                if(is_match)
                    {
                        // what[0] is the complete string mu<i>=<value>
                        // what[1] is the submatch <value>
                        //std::cout << "what 1 = " << what[1] << std::endl;
                        if( i!=mu_size )
                            input << what[1] << " , ";
                        else
                            input << what[1] << " ]";
                    }
            }
    }

private:
    std::string M_modelName;
    CRBModelMode M_mode;
    crbmodel_ptrtype model;
    bool use_newton_;

    crb_ptrtype crb;

    // For SCM convergence study
    std::map<std::string, std::vector<vectorN_type> > M_mapConvSCM;

    //vector of sampling to stock parameters for which the efficiency is under 1
    std::vector< sampling_ptrtype > vector_sampling_for_primal_efficiency_under_1;
    std::vector< sampling_ptrtype > vector_sampling_for_dual_efficiency_under_1;

    fs::path M_current_path;

    ser_ptrtype M_ser;
}; // OpusApp

} // Feel

#include <feel/feelmor/opusapp_impl.hpp>


#endif /* __OpusApp_H */
