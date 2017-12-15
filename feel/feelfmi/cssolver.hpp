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

    void simulate( double const& t_init, double const& t_final, double const& tol ) override
    {
        M_model->setupExperiment( t_init, t_final, tol );
        M_model->initialize();
        double t_cur = t_init;
        while( t_cur<t_final )
        {
            M_model->doStep( t_cur, M_step, true );
            t_cur += M_step;
        }
    }

private :

}; // CSSolver

} // namespace Feel

#endif
