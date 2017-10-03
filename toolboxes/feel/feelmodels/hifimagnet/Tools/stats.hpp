#ifndef __STATS_HPP
#define __STATS_HPP

#include <feel/feelcore/environment.hpp>
#include <feel/feelvf/evaluator.hpp>

/** include modelproperties **/
#include <feel/feelmodels/modelproperties.hpp>

namespace Feel
{
    //! Model properties type
    typedef ModelProperties model_prop_type;
    typedef std::shared_ptr<model_prop_type> model_prop_ptrtype;

    template <typename element_type> void stats_bcs(element_type M_V,  const model_prop_ptrtype& M_modelProps, const std::string Field, const std::string Unit)
    {
        Feel::cout << "\n**** Average Flux(" << Field << ") " << Unit << " ***\n" << std::flush;
        auto M_mesh = M_V.functionSpace()->mesh();

        auto itField = M_modelProps->boundaryConditions().find( Field );
        if ( itField != M_modelProps->boundaryConditions().end() )
            {
                auto mapField = (*itField).second;
                auto itType = mapField.find( "Neumann" );
                if ( itType != mapField.end() )
                    {
                        for ( auto const& exAtMarker : (*itType).second )
                            {
                                std::string marker = exAtMarker.marker();
                                auto g = expr(exAtMarker.expression());
                                auto flux = integrate( markedfaces(M_mesh, marker), -g ).evaluate()(0,0);
                                Feel::cout << marker << "[Neumann]\t" << flux << "\n";
                            }
                    }
                itType = mapField.find( "Robin" );
                if ( itType != mapField.end() )
                    {
                        for ( auto const& exAtMarker : (*itType).second )
                            {
                                std::string marker = exAtMarker.marker();
                                auto g1 = expr(exAtMarker.expression1());
                                auto g2 = expr(exAtMarker.expression2());
                                if ( Field == "Temperature" )
                                    {
                                        auto flux = integrate( markedfaces(M_mesh, marker), g1*(idv(M_V)-g2) ).evaluate()(0,0);
                                        Feel::cout << marker << "[Robin]\t" << flux << "\n";
                                    }
                                else
                                    {
                                        auto flux = integrate( markedfaces(M_mesh, marker), g1*idv(M_V)+g2 ).evaluate()(0,0);
                                        Feel::cout << marker << "[Robin]\t" << flux << "\n";
                                    }

                            }
                    }
            }
        Feel::cout << "***************************************\n" << std::flush;
    }

    template <typename element_type> void stats(element_type M_V, const std::string Field, const std::string Unit, const bool detail = false)
    {
        int proc_rank = Environment::worldComm().globalRank();

        auto M_mesh = M_V.functionSpace()->mesh();

        auto m = minmax( _range=elements(M_mesh), _pset=_Q<2>(), _expr=idv(M_V));
        auto vmax = m.max();
        auto vmin = m.min();
        using mdata_ptrtype = decltype(vmin);

        auto p_vmax = m.argmin();
        auto p_vmin = m.argmax();
        using pdata_ptrtype = decltype(p_vmin);

        auto mean = integrate( elements(M_mesh), idv(M_V) ).evaluate()(0,0);
        auto measure = integrate( elements(M_mesh), cst(1.0) ).evaluate()(0,0);

        auto std_dev = normL2( elements(M_mesh), (idv(M_V)-cst(mean/measure)) );
        using data_ptrtype = decltype(std_dev);

        if( proc_rank == 0 )
            {
                std::cout << "\n**** " << Field << " ****\n" << std::flush;
                Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "{", "}");
                std::cout << "min = " << vmin << " " << Unit << " at " << p_vmin.format(CleanFmt) << "\n";
                std::cout << "max = " << vmax << " " << Unit << " at " << p_vmax.format(CleanFmt) << "\n";
                std::cout << "mean = " << mean / measure << " " << Unit << "\n";
                std::cout << "std_dev = " << math::sqrt(std_dev*std_dev / measure) << " " << Unit << "\n" << std::endl;
            }

        // /* get stats by domain */
        // std::vector<mdata_ptrtype> vmax_;
        // std::vector<mdata_ptrtype> vmin_;
        // std::vector<pdata_ptrtype> p_vmax_;
        // std::vector<pdata_ptrtype> p_vmin_;
        // std::vector<data_ptrtype> mean_;
        // std::vector<data_ptrtype> std_dev_;

        if (detail)
            {
                for( auto marker: M_mesh->markerNames() )
                    {
                        auto name = marker.first;
                        auto data = marker.second;
                        // how to get rid of name that are not in M_mesh (ie traces of createsubmesh)
                        if ( data[1] == M_mesh->dimension() ) //&& nelements( markedelements(M_mesh, name) ) > 0 )
                            {
                                auto m_l = minmax( _range=markedelements(M_mesh, name), _pset=_Q<2>(), _expr=idv(M_V));
                                mdata_ptrtype vmax_l = m_l.max();
                                mdata_ptrtype vmin_l = m_l.min();
                                pdata_ptrtype p_vmax_l = m_l.argmin();
                                pdata_ptrtype p_vmin_l = m_l.argmax();

                                // vmax_.push_back(vmax_l);
                                // vmin_.push_back(vmin_l);
                                // p_vmax_.push_back(p_vmax_l);
                                // p_vmin_.push_back(p_vmin_l);

                                data_ptrtype mean_l = integrate( markedelements(M_mesh, name), idv(M_V) ).evaluate()(0,0);
                                data_ptrtype measure_l =  integrate( markedelements(M_mesh, name), cst(1.0) ).evaluate()(0,0);
                                // mean_.push_back(mean_l/measure_l);

                                data_ptrtype std_dev_l = normL2( markedelements(M_mesh, name), (idv(M_V)-cst(mean_l/measure_l)) );
                                // std_dev_.push_back(math::sqrt(std_dev_l*std_dev_l / measure_l));

                                if( proc_rank == 0 )
                                    {
                                        Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "{", "}");
                                        std::cout << name << ":\n";
                                        std::cout << "\tmin = " << vmin_l << " " << Unit << " at " << p_vmin_l.format(CleanFmt) << "\n";
                                        std::cout << "\tmax = " << vmax_l << " " << Unit << " at " << p_vmax_l.format(CleanFmt) << "\n";
                                        std::cout << "\tmean = " << mean_l/measure_l << " " << Unit << "\n";
                                        std::cout << "\tstd_dev = " << math::sqrt(std_dev_l*std_dev_l / measure_l) << " " << Unit << "\n" << std::endl;
                                    }

                            }
                    }
                if( proc_rank == 0 )
                    std::cout << "*******************\n" << std::endl;
            }

        // if( proc_rank == 0 )
        //     {
        //         std::cout << "\n**** " << Field << " ****\n" << std::flush;
        //         Eigen::IOFormat CleanFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "{", "}");
        //         std::cout << "min = " << vmin << " " << Unit << " at " << p_vmin.format(CleanFmt) << "\n";
        //         std::cout << "max = " << vmax << " " << Unit << " at " << p_vmax.format(CleanFmt) << "\n";
        //         std::cout << "mean = " << mean / measure << " " << Unit << "\n";
        //         std::cout << "std_dev = " << math::sqrt(std_dev*std_dev / measure) << " " << Unit << "\n" << std::flush;

        //         if (detail)
        //             {
        //                 std::cout << "\n";
        //                 int i = 0;
        //                 for( auto marker: M_mesh->markerNames() )
        //                     {
        //                         auto name = marker.first;
        //                         auto data = marker.second;
        //                         if ( data[1] == M_mesh->dimension() ) //&& nelements( markedelements(M_mesh, name) ) > 0 )
        //                             {
        //                                 std::cout << name << ":\n";
        //                                 std::cout << "\tmin = " << vmin_[i] << " " << Unit << " at " << p_vmin_[i].format(CleanFmt) << "\n";
        //                                 std::cout << "\tmax = " << vmax_[i] << " " << Unit << " at " << p_vmax_[i].format(CleanFmt) << "\n";
        //                                 std::cout << "\tmean = " << mean_[i] << " " << Unit << "\n";
        //                                 std::cout << "\tstd_dev = " << std_dev_[i] << " " << Unit << "\n" << std::flush;
        //                                 i++;
        //                             }
        //                     }
        //             }
        //         std::cout << "*******************\n" << std::flush;
        //     }
    }
}

#endif /* __STATS_HPP */
