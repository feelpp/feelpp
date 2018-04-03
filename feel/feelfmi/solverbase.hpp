#ifndef SOLVERBASE_HPP
#define SOLVERBASE_HPP

#include "feel/feelfmi/fmumodelbase.hpp"

namespace Feel
{

/**
 * This virtual class provide a generic interface for all FMI solvers
 */
class SolverBase
{
public :
    typedef FmuModelBase fmumodel_type;
    typedef boost::shared_ptr<fmumodel_type> fmumodel_ptrtype;

    SolverBase()
    {}

    SolverBase( fmumodel_ptrtype model ) :
        M_model( model ),
        M_tcur( 0 )
    {}

    virtual ~SolverBase()
    {}

    void setTimeStep( double const& step )
    {
        M_step = step;
    }

    void setRelativeTol( double const& tol )
    {
        M_tol = tol;
    }

    double currentTime()
    {
        return M_tcur;
    }

    virtual void initialize( double const& t_init, double const& t_final, double const& tol )=0;
    virtual void simulate()=0;
    virtual void doSteps( double t_stop )=0;

protected :
    fmumodel_ptrtype M_model;
    double M_step, M_tol, M_tcur, M_tfinal;

}; //classe SovlerBase


} // namespace Feel

#endif
