#ifndef __RESULTS_HPP__
#define __RESULTS_HPP__

#include <iostream>
#include <vector>
#include <fstream>
#include <openturns/OT.hxx>
#include <string>
#include <algorithm>
#include <algorithm>
// #include <feel/feel.hpp>


class Results
{
public:
    /**
     * @brief Construct a new Results object
     *
     * @param dim Dimension of the problem
     * @param names vector of names of the variables
     * @param algo algorithm used
     * @param size size of the sample
     */
    Results( size_t dim, const std::vector<std::string>& names, std::string algo, size_t size ) :
        M_dim(dim), M_names(names), M_algo(algo), M_size(size),
        M_firstOrder(dim), M_totalOrder(dim), M_firstOrderMin(dim), M_firstOrderMax(dim), M_totalOrderMin(dim), M_totalOrderMax(dim)
    {
        std::fill(M_firstOrder.begin(), M_firstOrder.end(), 0.0);
        std::fill(M_totalOrder.begin(), M_totalOrder.end(), 0.0);
        std::fill(M_firstOrderMin.begin(), M_firstOrderMin.end(), 1.0);
        std::fill(M_totalOrderMin.begin(), M_totalOrderMin.end(), 1.0);
        std::fill(M_firstOrderMax.begin(), M_firstOrderMax.end(), 0.0);
        std::fill(M_totalOrderMax.begin(), M_totalOrderMax.end(), 0.0);
    };
    ~Results() {};

    // Accessors
    OT::Scalar getFirstOrder( size_t i ) const { return M_firstOrder[i]; };
    OT::Scalar getTotalOrder( size_t i ) const { return M_totalOrder[i]; };
    OT::Scalar getFirstOrderMin( size_t i ) const { return M_firstOrderMin[i]; };
    OT::Scalar getTotalOrderMin( size_t i ) const { return M_totalOrderMin[i]; };
    OT::Scalar getFirstOrderMax( size_t i ) const { return M_firstOrderMax[i]; };
    OT::Scalar getTotalOrderMax( size_t i ) const { return M_totalOrderMax[i]; };

    std::vector<OT::Scalar> getFirstOrder() const { return M_firstOrder; };
    std::vector<OT::Scalar> getTotalOrder() const { return M_totalOrder; };
    std::vector<OT::Scalar> getFirstOrderMin() const { return M_firstOrderMin; };
    std::vector<OT::Scalar> getTotalOrderMin() const { return M_totalOrderMin; };
    std::vector<OT::Scalar> getFirstOrderMax() const { return M_firstOrderMax; };
    std::vector<OT::Scalar> getTotalOrderMax() const { return M_totalOrderMax; };

    // Mutators
    void setSamplingSize( size_t size ) { M_size = size; };

    /**
     * @brief Set the Sobol indice
     *
     * @param indice value obtained with the algorithm
     * @param i index computed
     * @param order order of the indice (1 for order 1, else total order)
     */
    void setIndice( OT::Scalar indice, size_t i, int order )
    {
        if (order == 1)
        {
            M_firstOrder[i] += indice;
            if (indice < M_firstOrderMin[i]) M_firstOrderMin[i] = indice;
            if (indice > M_firstOrderMax[i]) M_firstOrderMax[i] = indice;
        }
        else
        {
            M_totalOrder[i] += indice;
            if (indice < M_totalOrderMin[i]) M_totalOrderMin[i] = indice;
            if (indice > M_totalOrderMax[i]) M_totalOrderMax[i] = indice;
        }
    }

    /**
     * @brief Reset indices and intervals
     *
     */
    void reset()
    {
        std::fill(M_firstOrder.begin(), M_firstOrder.end(), 0.0);
        std::fill(M_totalOrder.begin(), M_totalOrder.end(), 0.0);
        std::fill(M_firstOrderMin.begin(), M_firstOrderMin.end(), 1.0);
        std::fill(M_totalOrderMin.begin(), M_totalOrderMin.end(), 1.0);
        std::fill(M_firstOrderMax.begin(), M_firstOrderMax.end(), 0.0);
        std::fill(M_totalOrderMax.begin(), M_totalOrderMax.end(), 0.0);
    }

    /**
     * @brief Set the Indices object
     *
     * @param P OT::Point containing the indices
     * @param order order of the indice (1 for order 1, else total order)
     */
    void setIndices( OT::Point P, int order )
    {
        size_t n = P.getDimension();
        if ( n != M_dim )
        {
            std::cout << "Error: dimension of the point is not the same as the dimension of the problem" << std::endl;
            return;
        }
        if (order == 1) for (size_t i = 0; i < n; i++)
        {
            M_firstOrder[i] = P[i];
        }
        else for (size_t i = 0; i < n; i++)
        {
            M_totalOrder[i] = P[i];
        }
    }

    void setInterval( const OT::Interval& interval, int order )
    {
        size_t n = interval.getDimension();
        if ( n != M_dim )
        {
            std::cout << "Error: dimension of the interval is not the same as the dimension of the problem" << std::endl;
            return;
        }
        if (order == 1) for (size_t i = 0; i < n; i++)
        {
            M_firstOrderMin[i] = interval.getLowerBound()[i];
            M_firstOrderMax[i] = interval.getUpperBound()[i];
        }
        else for (size_t i = 0; i < n; i++)
        {
            M_totalOrderMin[i] = interval.getLowerBound()[i];
            M_totalOrderMax[i] = interval.getUpperBound()[i];
        }
    }

    /**
     * @brief Normalize indices by the number of execution runned
     *
     * @param nrun number of execution runned
     */
    void normalize( int nrun )
    {
        for (size_t i = 0; i < M_dim; ++i)
        {
            M_firstOrder[i] /= nrun;
            M_totalOrder[i] /= nrun;
        }
    }

    /**
     * @brief Print the results in the console
     */
    void print()
    {
        Feel::cout << "Parameter names: " << M_names << std::endl;
        Feel::cout << "First order indices: " << M_firstOrder << std::endl;
        Feel::cout << "Total order indices: " << M_totalOrder << std::endl;
        Feel::cout << "FirstOrderIntervals" << std::endl;
        for (size_t i=0; i < M_dim; ++i)
            Feel::cout << "\t[" << M_firstOrderMin[i] << ", " << M_firstOrderMax[i] << "]" << std::endl;
        Feel::cout << "TotalOrderIntervals" << std::endl;
        for (size_t i=0; i < M_dim; ++i)
            Feel::cout << "\t[" << M_totalOrderMin[i] << ", " << M_totalOrderMax[i] << "]" << std::endl;
    }

    /**
     * @brief Export values computed in a json file
     *
     * @param filename path to the exported file
     */
    void exportValues( const std::string filename )
    {
        std::ofstream file;
        file.open(filename);
        file << "{\n\t\"N\": " << M_dim << ",\n";
        file << "\t\"sampling-size\": " << M_size << ",\n";
        file << "\t\"algo\": \"" << M_algo << "\",\n";
        file << "\t\"Names\": [";
        for (size_t i = 0; i < M_dim; ++i)
        {
            file << "\"" << M_names[i] << "\"";
            if (i != M_dim - 1)
                file << " ,";
        }
        file << "]," << std::endl;
        file << "\t\"FirstOrder\":\n\t{" << std::endl;
        file << "\t\t\"values\": [";
        for (size_t i = 0; i < M_dim; ++i)
        {
            file << M_firstOrder[i];
            if (i != M_dim - 1)
            {
                file << ", ";
            }
        }
        file << "]," << std::endl;
        file << "\t\t\"intervals\": [";
        for (size_t i = 0; i < M_dim; ++i)
        {
            file << "[" << M_firstOrderMin[i] << ", " << M_firstOrderMax[i] << "]";
            if (i != M_dim - 1)
            {
                file << ", ";
            }
        }
        file << "]" << std::endl;
        file << "\t}," << std::endl;
        file << "\t\"TotalOrder\":\n\t{" << std::endl;
        file << "\t\t\"values\": [";
        for (size_t i = 0; i < M_dim; ++i)
        {
            file << M_totalOrder[i];
            if (i != M_dim - 1)
            {
                file << ", ";
            }
        }
        file << "]," << std::endl;
        file << "\t\t\"intervals\": [";
        for (size_t i = 0; i < M_dim; ++i)
        {
            file << "[" << M_totalOrderMin[i] << ", " << M_totalOrderMax[i] << "]";
            if (i != M_dim - 1)
            {
                file << ", ";
            }
        }
        file << "]" << std::endl;
        file << "\t}" << std::endl;
        file << "}" << std::endl;
        file.close();
        Feel::cout << "Results exported in " << filename << std::endl;
    }

private:
    size_t M_dim;
    std::vector<std::string> M_names;
    std::string M_algo;
    size_t M_size;
    std::vector<OT::Scalar> M_firstOrder, M_totalOrder;
    std::vector<OT::Scalar> M_firstOrderMin, M_firstOrderMax;
    std::vector<OT::Scalar> M_totalOrderMin, M_totalOrderMax;
};



#endif // __RESULTS_HPP__