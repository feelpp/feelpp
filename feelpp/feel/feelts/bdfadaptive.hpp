#ifndef FEELPP_TS_BDF_ADAPTIVE_H
#define FEELPP_TS_BDF_ADAPTIVE_H

#include <vector>
#include <memory>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelts/tsbase.hpp>
#include <feel/feelalg/glas.hpp>

namespace Feel
{

/**
 * \class BdfAdaptive
 * \brief Adaptive BDF for time stepping with variable \Delta t
 */
template<typename SpaceType>
class BdfAdaptive : public TSBase
{
    using super = TSBase;
public:
    using space_type = SpaceType;
    using space_ptrtype = std::shared_ptr<space_type>;
    using element_type = typename space_type::element_type;
    using element_ptrtype = typename space_type::element_ptrtype;
    using value_type = typename element_type::value_type;
    using unknowns_type = std::vector<element_ptrtype>;

    BdfAdaptive(space_ptrtype const& space,
        std::string const& name,
        std::string const& prefix = "",
        po::variables_map const& vm = Environment::vm())
    : super(name, prefix, space->worldComm(), vm),
      M_space(space)
    {
        M_dt_min = doption(_prefix=prefix, _name="bdfadaptive.dt-min", _vm=vm);
        M_dt_max = doption(_prefix=prefix, _name="bdfadaptive.dt-max", _vm=vm);
        M_dt     = doption(_prefix=prefix, _name="bdfadaptive.initial-dt", _vm=vm);
        M_adaptEnabled = boption(_prefix=prefix, _name="bdfadaptive.adapt", _vm=vm);
    }
    bool adaptivityEnabled() const { return M_adaptEnabled; }

    void enableAdaptivity(bool enable = true)
    {
        M_adaptEnabled = enable;
    }
    /**
     * @brief Set the Current Time Step object
     * 
     * @param dt 
     */
    void setCurrentTimeStep(double dt)
    {
        if (dt <= 0.0)
            throw std::invalid_argument("BdfAdaptive: dt must be > 0");
    
        if (dt < M_dt_min)
        {
            LOG(WARNING) << "BdfAdaptive: clamping dt to minimum allowed value " << M_dt_min;
            dt = M_dt_min;
        }
        if (dt > M_dt_max)
        {
            LOG(WARNING) << "BdfAdaptive: clamping dt to maximum allowed value " << M_dt_max;
            dt = M_dt_max;
        }
    
        M_dt = dt;
    }
    void setTimeStepBounds(double dt_min, double dt_max)
    {
        M_dt_min = dt_min;
        M_dt_max = dt_max;
    }
        /// @brief Return the current time step value
    double currentTimeStep() const
    {
        return M_dt;
    }
    void setTimeStamps(std::vector<double> const& times) { M_timeStamps = times; invalidateCoefficients(); }
    void setHistory(std::vector<element_ptrtype> const& history) { M_history = history; invalidateCoefficients(); }

    void setHistory(std::vector<element_type> const& vec)
    {
        M_history.resize(vec.size());
        for (std::size_t i = 0; i < vec.size(); ++i)
        {
            M_history[i] = M_space->elementPtr();
            *M_history[i] = vec[i];
        }
        invalidateCoefficients();
    }
    template<typename... Elements>
    void setHistory(Elements const&... elems)
    {
        std::array<element_type const*, sizeof...(elems)> args = { &elems... };
        M_history.resize(args.size());
        for (std::size_t i = 0; i < args.size(); ++i)
        {
            M_history[i] = M_space->elementPtr();
            *M_history[i] = *args[i];
        }
    }

    int timeOrder() const { return M_timeStamps.size() - 1; }
    unknowns_type const& history() const { return M_history; }
    unknowns_type& history() { return M_history; }

    element_type const& firstDerivative() const
    {
        if (!M_firstDeriv)
            M_firstDeriv = M_space->elementPtr();
        M_firstDeriv->zero();

        const int k = M_timeStamps.size();
        if (M_dirty_d1)
            computeFirstDerivCoefficients();

        for (int i = 0; i < k; ++i)
            M_firstDeriv->add(M_cached_d1[i], *M_history[i]);

        return *M_firstDeriv;
    }

    element_type const& secondDerivative() const
    {
        if (!M_secondDeriv)
            M_secondDeriv = M_space->elementPtr();
        M_secondDeriv->zero();

        const int k = M_timeStamps.size();
        if (M_dirty_d2)
            computeSecondDerivCoefficients();

        for (int i = 0; i < k; ++i)
            M_secondDeriv->add(M_cached_d2[i], *M_history[i]);

        return *M_secondDeriv;
    }
    void setExtrapolationTime(double t_target) const
    {
        if (M_dirty_poly || t_target != M_lastPolyTarget)
        {
            computePolyCoefficients(t_target);
            M_lastPolyTarget = t_target;
        }
    }
    element_type const& extrapolation(double t_target = -1.0) const
    {
        if (!M_extrapolated)
            M_extrapolated = M_space->elementPtr();
        M_extrapolated->zero();

        int k = M_timeStamps.size();
        if (t_target < 0)
            t_target = 2 * M_timeStamps[0] - M_timeStamps[1]; // default: forward extrapolation

        if (M_dirty_poly || t_target != M_lastPolyTarget)
        {
            computePolyCoefficients(t_target);
            M_lastPolyTarget = t_target;
        }

        for (int i = 0; i < k; ++i)
            M_extrapolated->add(M_cached_poly[i], *M_history[i]);

        return *M_extrapolated;
    }
    double polyCoefficient(int i) const
    {
        if (M_dirty_poly)
            computePolyCoefficients(2 * M_timeStamps[0] - M_timeStamps[1]);
        return M_cached_poly[i];
    }

    double polyDerivCoefficient(int i) const
    {
        if (M_dirty_d1)
            computeFirstDerivCoefficients();
        return M_cached_d1[i];
    }

    double polySecondDerivCoefficient(int i) const
    {
        if (M_dirty_d2)
            computeSecondDerivCoefficients();
        return M_cached_d2[i];
    }
    element_type const& poly(double t_target = -1.0) const
    {
        if (!M_poly)
            M_poly = M_space->elementPtr();
        M_poly->zero();

        int k = M_timeStamps.size();
        if (t_target < 0)
            t_target = 2 * M_timeStamps[0] - M_timeStamps[1]; // Default extrapolation time

        computePolyCoefficients(t_target);  // updates M_cached_poly

        for (int i = 0; i < k; ++i)
            M_poly->add(M_cached_poly[i], *M_history[i]);

        return *M_poly;
    }

    element_type const& polyDeriv() const
    {
        if (!M_polyDeriv)
            M_polyDeriv = M_space->elementPtr();
        M_polyDeriv->zero();

        const int k = M_timeStamps.size();
        if (M_dirty_d1)
            computeFirstDerivCoefficients();

        // Use coefficients alpha_1, ..., alpha_k (skip alpha_0 which multiplies u^{n+1})
        for (int i = 1; i < k; ++i)
            M_polyDeriv->add(-M_cached_d1[i], *M_history[i-1]);

        return *M_polyDeriv;
    }

    element_type const& polySecondDeriv() const
    {
        if (!M_polySecondDeriv)
            M_polySecondDeriv = M_space->elementPtr();
        M_polySecondDeriv->zero();

        const int k = M_timeStamps.size();
        if (M_dirty_d2)
            computeSecondDerivCoefficients();

        for (int i = 1; i < k; ++i)
            M_polySecondDeriv->add(-M_cached_d2[i], *M_history[i-1 ]);

        return *M_polySecondDeriv;
    }

    template<typename Element>
    void shiftRight(Element const& newU, double newT)
    {
        for (std::size_t i = M_history.size() - 1; i > 0; --i)
        {
            *M_history[i] = *M_history[i-1];
            M_timeStamps[i] = M_timeStamps[i-1];
        }
        *M_history[0] = newU;
        M_timeStamps[0] = newT;

        invalidateCoefficients();
    }
    std::vector<double> extrapolationCoefficients(double t) const
    {
        int k = M_timeStamps.size();
        std::vector<double> coeffs(k);

        for (int j = 0; j < k; ++j)
        {
            coeffs[j] = 1.0;
            for (int i = 0; i < k; ++i)
                if (i != j)
                    coeffs[j] *= (t - M_timeStamps[i]) / (M_timeStamps[j] - M_timeStamps[i]);
        }
        return coeffs;
    }

    /**
     * @brief Extrapolate trajectory based on target times
     * 
     * @param targetTimes A vector of target times for extrapolation
     * @return std::vector<element_type> A vector of extrapolated elements
     */
    std::vector<element_type> extrapolateTrajectory(std::vector<double> const& targetTimes) const
    {
        std::vector<element_type> results;
        results.reserve(targetTimes.size());

        for (auto const& t : targetTimes)
        {
            auto coeffs = extrapolationCoefficients(t);

            auto result = M_space->element();
            result->zero();
            for (std::size_t i = 0; i < coeffs.size(); ++i)
                result->add(coeffs[i], *M_history[i]);

            results.emplace_back(std::move(*result));
        }

        return results;
    }
    /**
     * @brief Compute extrapolation coefficients for order k-1
     *
     * This function computes the coefficients for the Newton interpolation polynomial
     * of order k-1, which is used for error estimation.
     *
     * @param t The target time for extrapolation.
     * @return A vector of coefficients for the Newton interpolation polynomial.
     */
    std::vector<double> extrapolationCoefficients_kminus1(double t) const
    {
        int k = M_timeStamps.size() - 1;
        std::vector<double> coeffs(k);

        for (int j = 0; j < k; ++j)
        {
            coeffs[j] = 1.0;
            for (int i = 0; i < k; ++i)
            {
                if (i != j)
                    coeffs[j] *= (t - M_timeStamps[i+1]) / (M_timeStamps[j+1] - M_timeStamps[i+1]);
            }
        }
        return coeffs;
    }
    /**
     * @brief Compute L2-norm extrapolation error over multiple target times
     *
     * This function computes the extrapolated solution at each time in `targetTimes`
     * using the cached history and Newton interpolation, then compares it against
     * a user-provided exact solution.
     *
     * @param targetTimes A vector of target times where extrapolation is evaluated.
     * @param exactExpr A callable that returns a Feel++ Expr for the exact solution at time t.
     *                  For example: [](double t) { return expr("1 + 2*t + 3*t^2 : t", {{"t", t}}); }
     * @return A vector of (t_i, error_i) pairs, where error_i is the L2 norm between the extrapolated
     *         and exact solution at time t_i.
     *
     * @note The method internally uses extrapolationCoefficients(t) and projection onto the same
     *       function space for both extrapolated and exact fields.
     * \code
     * auto times = std::vector<double>{1.0, 1.1, 1.2};
     * auto exact_u = "1 + 2*t + 3*t^2 :t"
     * auto errs = bdfAdaptive->extrapolationErrorTrajectory(times, exact_u);
     * for (auto const& [t, e] : errs)
     *     std::cout << "Error at t=" << t << " is " << e << "\n";
     * \endcode
     */
    std::vector<std::pair<double, double>>
    extrapolationErrorTrajectory(std::vector<double> const& targetTimes,
                                 std::string exprs_ref ) const
    {
        std::vector<std::pair<double, double>> errors;

        for (auto const& t : targetTimes)
        {
            auto coeffs = extrapolationCoefficients(t);

            auto u_extrap = M_space->element();
            u_extrap->zero();
            for (std::size_t i = 0; i < coeffs.size(); ++i)
                u_extrap->add(coeffs[i], *M_history[i]);

            auto u_ref_expr = expr(exprs_ref);
            u_ref_expr.setParameterValues({{"t", t}});
            auto u_ref = project(_space=M_space, _expr=u_ref_expr);

            double err = normL2(_range=elements(M_space->mesh()),
                                _expr=idv(u_ref) - idv(u_extrap));

            errors.emplace_back(t, err);
        }

        return errors;
    }
    /**
     * @brief Estimate local error at next time t_{n+1} using extrapolation order reduction
     *
     * This uses the difference between extrapolation of degree k and degree k-1 as an error estimate.
     * @return L2 norm of difference (used for adaptive time-stepping).
     */
    double estimateError(double t_next) const
    {
        auto coeffs_k = extrapolationCoefficients(t_next);        // full order
        LOG(INFO) << fmt::format("[BDFAdaptive] coeffs_k = {}", coeffs_k);
        auto coeffs_km1 = extrapolationCoefficients_kminus1(t_next); // lower order
        LOG(INFO) << fmt::format("[BDFAdaptive] coeffs_km1 = {}", coeffs_km1);
        
        auto u_k = M_space->element();
        auto u_km1 = M_space->element();
        u_k.zero();
        u_km1.zero();

        for (std::size_t i = 0; i < coeffs_k.size(); ++i)
            u_k.add(coeffs_k[i], *M_history[i]);
        for (std::size_t i = 0; i < coeffs_km1.size(); ++i)
            u_km1.add(coeffs_km1[i], *M_history[i]);

        return normL2(_range=elements(M_space->mesh()), _expr=idv(u_k) - idv(u_km1));
    }

    /**
     * @brief Compute new dt using embedded error estimate
     *
     * @param err current error estimate
     * @param tol user-defined tolerance
     * @param safety default = 0.9
     * @param p extrapolation order (usually = BDF order)
     * @return new suggested time step
     */
    double adaptTimeStep(double err, double tol, double safety = 0.9, int p = -1) const
    {
        if (p < 0) p = this->timeOrder();
        if (err < 1e-14) err = 1e-14; // avoid div-by-zero

        return safety * currentTimeStep() * std::pow(tol / err, 1.0 / (p + 1));
    }

    //! Check if the step is accepted based on the error estimate
    /*!
     * \param err The estimated error.
     * \param tol The tolerance.
     * \return true if the step is accepted, false otherwise.
     */
    bool isStepAccepted(double err, double tol) const
    {
        return err <= tol;
    }
    struct TimeAdaptResult {
        bool accepted;
        double new_dt;
        double err;
    };
    
    TimeAdaptResult checkAndAdapt(double t_next, double tol)
    {
        double err = estimateError(t_next);
        double dt_new = adaptTimeStep(err, tol);
        bool ok = isStepAccepted(err, tol);
        return {ok, dt_new, err};
    }

    /**
     * @brief Advance the time step and update the solution
     *
     * @param u The current solution.
     * @param t_final The final time to reach.
     * @param tol The tolerance for adaptivity.
     * @param doStep User-supplied function to perform the time step.
     */
    template<typename StepFunction>
    void advance(element_type& u, double t_final, double tol, StepFunction&& doStep)
    {
        while (this->time() < t_final)
        {
            LOG(INFO) << fmt::format("[BDFAdaptive] t = {:.6f}, dt = {:.3e}", this->time(), this->currentTimeStep());
            double t_next = this->time() + this->currentTimeStep();

            if (adaptivityEnabled())
            {
                auto result = this->checkAndAdapt(t_next, tol);

                if (!result.accepted)
                {
                    this->setCurrentTimeStep(result.new_dt);
                    DVLOG(1) << fmt::format("[BDFAdaptive] Step rejected at t = {:.6f}, new dt = {:.3e}", this->time(), result.new_dt);
                    continue; // retry
                }

                this->setCurrentTimeStep(result.new_dt);
            }

            // User-supplied code to assemble and solve system for u^{n+1}
            doStep(t_next);

            // Shift time and history buffers
            this->shiftRight(u, t_next);

            // Advance internal time state
            this->next();
        }
    }
private:
    void invalidateCoefficients() const
    {
        M_dirty_poly = true;
        M_dirty_d1 = true;
        M_dirty_d2 = true;
    }
    void computePolyCoefficients(double t_target) const
    {
        int k = M_timeStamps.size();
        M_cached_poly.resize(k);
        for (int j = 0; j < k; ++j)
        {
            double prod = 1.0;
            for (int i = 0; i < k; ++i)
                if (i != j)
                    prod *= (t_target - M_timeStamps[i]) / (M_timeStamps[j] - M_timeStamps[i]);
            M_cached_poly[j] = prod;
        }
        M_dirty_poly = false;
    }

    void computeFirstDerivCoefficients() const
    {
        int k = M_timeStamps.size();
        M_cached_d1.resize(k);
        for (int j = 0; j < k; ++j)
        {
            double sum = 0.0;
            for (int m = 0; m < k; ++m)
            {
                if (m == j) continue;
                double prod = 1.0;
                for (int l = 0; l < k; ++l)
                    if (l != j && l != m)
                        prod *= (M_timeStamps[0] - M_timeStamps[l]) / (M_timeStamps[j] - M_timeStamps[l]);
                sum += prod / (M_timeStamps[j] - M_timeStamps[m]);
            }
            M_cached_d1[j] = sum;
        }
        M_dirty_d1 = false;
    }

void computeSecondDerivCoefficients() const
{
    int k = M_timeStamps.size();
    M_cached_d2.resize(k);
    for (int j = 0; j < k; ++j)
    {
        double sum = 0.0;
        for (int m = 0; m < k; ++m)
        {
            if (m == j) continue;
            for (int n = 0; n < k; ++n)
            {
                if (n == j || n == m) continue;
                double prod = 1.0;
                for (int l = 0; l < k; ++l)
                    if (l != j && l != m && l != n)
                        prod *= (M_timeStamps[0] - M_timeStamps[l]) / (M_timeStamps[j] - M_timeStamps[l]);
                sum += prod / ((M_timeStamps[j] - M_timeStamps[m]) * (M_timeStamps[j] - M_timeStamps[n]));
            }
        }
        M_cached_d2[j] = sum;
    }
    M_dirty_d2 = false;
}
private:

    bool M_adaptEnabled = true;
    double M_dt = 0.0;
    double M_dt_min = 1e-10;
    double M_dt_max = 1e3;
    space_ptrtype M_space;
    std::vector<double> M_timeStamps;
    std::vector<element_ptrtype> M_history;
    mutable element_ptrtype M_firstDeriv;
    mutable element_ptrtype M_secondDeriv;
    mutable element_ptrtype M_extrapolated;
    mutable element_ptrtype M_polyDeriv;
    mutable element_ptrtype M_poly, M_polySecondDeriv;
    mutable std::vector<double> M_cached_poly;
    mutable std::vector<double> M_cached_d1;
    mutable std::vector<double> M_cached_d2;
    mutable bool M_dirty_poly = true;
    mutable bool M_dirty_d1 = true;
    mutable bool M_dirty_d2 = true;
    mutable double M_lastPolyTarget = std::numeric_limits<double>::quiet_NaN();
};

} // namespace Feel

#endif // FEELPP_TS_BDF_ADAPTIVE_H