#ifndef CSSOLVER_HPP
#define CSSOLVER_HPP

#include "feel/feelfmi/solverbase.hpp"

namespace Feel
{
class CSSolver :
        public SolverBase
{
public :
    typedef SolverBase super_type;
    typedef FmuModelBase fmumodel_type;
    typedef boost::shared_ptr<fmumodel_type> fmumodel_ptrtype;

    CSSolver() :
        super_type()
    {}

    CSSolver( fmumodel_ptrtype model ) :
        super_type( model )
    {}

    void initialize( double const& t_init, double const& t_final, double const& tol ) override
    {
        M_model->setupExperiment( t_init, t_final, tol );
        M_model->initialize();
        M_tcur = t_init;
        M_tfinal = t_final;
    }

    void simulate() override
    {
        while( M_tcur<M_tfinal )
        {
            M_model->doStep( M_tcur, M_step, true );
            M_tcur += M_step;
        }
    }

    void doSteps( double t_stop ) override
    {
        while( M_tcur<t_stop && M_tcur<M_tfinal )
        {
            M_model->doStep( M_tcur, M_step, true );
            M_tcur += M_step;
        }
    }

private :

}; // CSSolver

} // namespace Feel

#endif
