/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-08-11

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file parameterspace.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-08-11
 */
#ifndef __ParameterSpace_H
#define __ParameterSpace_H 1

#include <vector>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <Eigen/Core>

#if defined(FEELPP_HAS_ANN_H)
#include <ANN/ANN.h>
#endif /* FEELPP_HAS_ANN_H */



#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/enumerate.hpp>
#include <feel/feelcore/commobject.hpp>
#include <feel/feelmesh/kdtree.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feelcore/hashtables.hpp>
#include <feel/feelmor/crbmodelparameters.hpp>

namespace Feel
{
/**
 * \class ParameterSpace
 * \brief Parameter space class
 *
 * @author Christophe Prud'homme
 * @see
 */
template<uint16_type P = invalid_uint16_type_value >
class ParameterSpace: public CommObject, public std::enable_shared_from_this<ParameterSpace<P> >
{
public:


    /** @name Constants
     */
    //@{
    static const int nDimEigenContainer = ( P == invalid_uint16_type_value )? Eigen::Dynamic : P;
    //! dimension of the parameter space
    static const uint16_type Dimension = P;
    static const uint16_type ParameterSpaceDimension = P;
    //@}

    /** @name Typedefs
     */
    //@{
    using super = CommObject;
    typedef ParameterSpace<Dimension> parameterspace_type;
    typedef std::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    using value_type = double;

    //@}

    /**
     * \brief element of a parameter space
     */
    class Element : public Eigen::VectorXd//Eigen::Matrix<double,nDimEigenContainer,1>
    {
        //typedef Eigen::Matrix<double,nDimEigenContainer,1> super;
        using super = Eigen::VectorXd;
    public:
        typedef ParameterSpace<Dimension> parameterspace_type;
        typedef std::shared_ptr<parameterspace_type> parameterspace_ptrtype;
        //typedef typename Eigen::internal::ref_selector<Element>::type Nested;
        typedef typename Eigen::internal::remove_all<Eigen::VectorXd>::type NestedExpression;

        /**
         * default constructor
         */
        Element() : super(), M_space(), param_indices_(){}

        /**
         * @brief Default constructor.
         *
         * Creates an empty Element object.
         */
        Element( Element const& e ) = default;

        /**
         * @brief Constructor with ParameterSpace pointer.
         *
         * Creates an Element object with a size based on the dimension of the ParameterSpace.
         *
         * @param space A shared pointer to a ParameterSpace object.
         */
        Element( parameterspace_ptrtype const& space )
            : super( space->dimension() ), M_space( space ), param_indices_()
        {
            precomputeParamIndices();
        }

        // This constructor allows you to construct Element from Eigen expressions
        template<typename OtherDerived>
        Element(const Eigen::MatrixBase<OtherDerived>& other)
            : super(other), param_indices_()
            { }
        /**
         * destructor
         */
        ~Element()
            {}

        /**
         * copy constructor
         */
        Element& operator=( Element const& e )
            {
                if ( this == &e )
                    return *this;

                this->resize( e.size() );
                super::operator=( e );
                M_space = e.M_space;
                param_indices_ = e.param_indices_;

                return *this;
            }

        template<typename OtherDerived>
        super& operator=( const Eigen::MatrixBase<OtherDerived>& other )
            {
                return super::operator=( other );
            }

        double const& coeff(int i ) const
            {
                return super::operator()(i);
            }
        double& coeff(int i )
            {
                return super::operator()(i);
            }

        /**
         * return name of parameter at index d
         */
        std::string parameterName( int d ) const
            {
                CHECK( M_space ) << "no parameter space";
                return M_space->parameterName(d);
            }

        /**
         * @brief return the list of all parameter names
         * 
         */
        std::vector<std::string> const& parameterNames() const
        {
            CHECK( M_space ) << "no parameter space";
            return M_space->parameterNames();
        }

        /**
         * @brief Accesses an element by its parameter name.
         * 
         * @param name The name of the parameter to access.
         * @return Reference to the value of the specified parameter.
         * @throws std::invalid_argument If the parameter name is not found.
         */
        double const &parameterNamed(std::string name) const
        {
            int index = paramIndexByName(name);
            return this->coeff(index);
        }

        /**
         * @brief Accesses an element by its parameter name (modifiable version).
         * 
         * @param name The name of the parameter to access.
         * @return Reference to the value of the specified parameter.
         * @throws std::invalid_argument If the parameter name is not found.
         */
        double &parameterNamed(std::string name)
        {
            int index = paramIndexByName(name);
            return this->coeff(index);
        }


        void view() const
        {
            for (int i=0; i<this->size(); ++i)
                std::cout << M_space->parameterName( i ) << "\t" << this->operator()( i ) << std::endl;
        }

        /**
         * @brief Sets the value of a parameter with a specific name.
         *
         * Checks whether the given value is within the valid range for the parameter.
         *
         * @param name The name of the parameter to set.
         * @param value The new value for the parameter.
         * @throws std::invalid_argument If the parameter name is not found or the value is out of range.
         */
        void setParameterNamed( std::string name, double value )
        {
            int index = paramIndexByName( name );
            if ( isValueInRange( index, value ) )
                this->coeff( index ) = value;
            else
            {
                LOG( WARNING ) << fmt::format( "{} value not in range [{}, {}] for parameter named {}", value, M_space->min( index ), M_space->max( index ) ) << std::endl;
                throw std::invalid_argument( fmt::format( "Parameter named = {} with index = {} is out of range, cannot set to value = {}", name, index, value ) );
            }
        }
    
        void setParameter( int i, double value)
        {
            element_type min = M_space->min(), max = M_space->max();
            if ( isValueInRange( i, value ) )
                this->operator()(i) = value;
            else
            {
                LOG( WARNING ) << fmt::format("{} value not in range [{}, {}] for parameter number {}", value, min(i), max(i), i) << std::endl;
                throw std::invalid_argument(fmt::format("Parameter with index = {} is out of range, cannot set to value = {}", i, value ) );
            }
        }

        /**
         * set many parameters at once (overloaded function)
         */
        void setParameters( const std::vector<double> &values)
        {
            size_t n = values.size();
            if (n != this->size())
                LOG( WARNING ) << fmt::format("The size of the given vector ({}) is different fom the size ({})", n, this->size()) << std::endl;
            size_t N = (n <= this->size()) ? n : this->size();
            element_type min = M_space->min(), max = M_space->max();
            for (size_t i=0; i<N; ++i)
            {
                if (values[i] >= min(i) && values[i] <= max(i))
                    this->operator()(i) = values[i];
                else
                {
                    LOG( WARNING ) << fmt::format("{} value not in range [{}, {}] for parameter number {}", values[i], min(i), max(i), i) << std::endl; 
                    throw std::invalid_argument("Parameter out of range");
                }
            }

        }

        /**
         * @brief Set the Parameters object
         * 
         * @param values map of parameter values
         */
        void setParameters( const std::map<std::string, double> &values )
        {
            for (auto a : values)
                setParameterNamed( a.first, a.second );
        }

        /**
         * @brief Set the parameters to the map of values
         * 
         * @param values map of values 
         * @return Element& 
         */
        Element& set( const std::map<std::string, double> &values )
        {
            for (auto a : values)
                setParameterNamed( a.first, a.second );
            return *this;
        }

        /**
         * @brief Set the parameters to the vector of values
         * 
         * @param values vector of values ordered by parameter index
         * @return Element& 
         */
        Element& set( const std::vector<double> &values )
        {
            setParameters( values );
            return *this;
        }

        /**
         * get index of named parameter
         */
        int indexOfParameterNamed( std::string name ) const
            {
                auto paramNames = M_space->parameterNames();
                auto it = std::find(paramNames.begin(), paramNames.end(), name);
                if( it != paramNames.end() )
                    return it - paramNames.begin();
                else
                    return -1;
            }

        void setParameterSpace( parameterspace_ptrtype const& space )
            {
                M_space = space;
                precomputeParamIndices();
            }


        /**
         * \brief Retuns the parameter space
         */
        parameterspace_ptrtype const& parameterSpace() const
            {
                return M_space;
            }

        void check() const
            {
#if O// !defined(NDEBUG)
                if ( !M_space->check() )
                {
                    LOG(INFO) << "No need to check element since parameter space is no valid (yet)\n";
                    return;
                }

                if ( M_space->dimension() == 0 )
                    return;
                Element sum( M_space );
                sum.setZero();
                // verify that the element is the same on all processors
                mpi::all_reduce( M_space->worldComm().localComm(), *this, sum,
                                 []( Element const& m1, Element const& m2 ) { return m1+m2; } );
                rank_type proc_number = M_space->worldComm().localRank();
                rank_type nProc = M_space->worldComm().localSize();
                if( (this->array()-sum.array()/nProc).abs().maxCoeff() > 1e-7 )
                {
                    std::cout << "Parameter not identical on all processors: "<< "current parameter on proc "<<proc_number<<" : [";
                    for(int i=0; i<this->size()-1; i++) std::cout <<std::setprecision(15)<< this->operator()(i) <<", ";
                    std::cout<<std::setprecision(15)<<this->operator()(this->size()-1)<<" ]";
                    std::cout <<std::setprecision(15)<< " and test" << (this->array()-sum.array()/nProc).abs().maxCoeff() << "\n";
                }
                M_space->worldComm().barrier();
                CHECK( (this->array()-sum.array()/nProc).abs().maxCoeff() < 1e-7 )
                    << "Parameter not identical on all processors(" << M_space->worldComm().masterRank() << "/" << nProc << ")\n"
                    << "max error: " << (this->array()-sum.array()/nProc).abs().maxCoeff() << "\n"
                    << "error: " << std::setprecision(15) << (this->array()-sum.array()/nProc).abs() << "\n"
                    << "current parameter " << std::setprecision(15) << *this << "\n"
                    << "sum parameter " << std::setprecision(15) << sum << "\n"
                    << "space min : " << M_space->min() << "\n"
                    << "space max : " << M_space->max() << "\n";
#endif
            }


        size_t key() const
        {
            int N = this->size();
            std::vector<std::string> s;
            for ( int i=0; i<N; i++ )
                s.push_back( std::to_string(this->operator[](i)) );
            HashTables::HasherContainers<std::string> h;
            return h(s);
        }

        std::string toString() const
        {
            std::ostringstream mu_str;
            mu_str << "["<<this->operator[](0);
            for ( int i=1; i<this->size(); i++ )
                mu_str <<","<< this->operator[](i);
            mu_str <<"]";
            return mu_str.str();
        }

        std::map<std::string, double> toParameterValues() const
        {
            std::map<std::string, double> pm;
            auto parameterNames = this->parameterNames();
            for( auto const& n : parameterNames )
                pm[n] = this->parameterNamed(n);
            return pm;
        }

    private:

      /**
       * @brief Map from parameter names to their corresponding indices.
       */
      std::unordered_map<std::string, int> param_indices_;

      /**
       * @brief Helper function to precompute the parameter indices during construction.
       */
      void precomputeParamIndices()
      {
            auto paramNames = M_space->parameterNames();
            for ( int i = 0; i < paramNames.size(); ++i )
                param_indices_[paramNames[i]] = i;
      }

      /**
       * @brief Helper function to fetch the parameter index by its name.
       *
       * @param name The name of the parameter to fetch the index for.
       * @return The index of the specified parameter.
       * @throws std::invalid_argument If the parameter name is not found.
       */
      int paramIndexByName( std::string name ) const
      {
            auto it = param_indices_.find( name );
            if ( it != param_indices_.end() )
                return it->second;
            else
                throw std::invalid_argument( fmt::format("Parameter named {} not found", name ) );
      }

      /**
       * @brief Helper function to check if a given value is within the valid range for a parameter.
       *
       * @param index The index of the parameter to check.
       * @param value The value to check.
       * @return true if the value is within the valid range, false otherwise.
       */
      bool isValueInRange( int index, double value ) const
      {
            auto min_val = M_space->min( index );
            auto max_val = M_space->max( index );
            return ( value >= min_val && value <= max_val );
      }

    private:

        friend class boost::serialization::access;
        template<class Archive>
        void save( Archive & ar, const unsigned int version ) const
            {
                ar & boost::serialization::base_object<super>( *this );
            }

        template<class Archive>
        void load( Archive & ar, const unsigned int version )
            {
                ar & boost::serialization::base_object<super>( *this );
            }
        BOOST_SERIALIZATION_SPLIT_MEMBER()

        private:
        parameterspace_ptrtype M_space;
    };

    typedef Element element_type;
    typedef std::shared_ptr<Element> element_ptrtype;
    element_type element( bool broadcast = true, bool apply_log = true )
    {
#if 0
        //first, pick a random element, then
        //look if there is negative elements and if not
        //pick a log-random element
        auto element = parameterspace_type::random( this->shared_from_this(), broadcast );
        auto mu_min = parameterspace_type::min();
        //auto mu_min = min_element.template get<0>();
        int size = mu_min.size();
        bool apply_log=true;
        for(int i=0; i<size; i++)
        {
            if( mu_min(i) < 0 )
                apply_log=false;
        }

        if( apply_log )
        {
            element = parameterspace_type::logRandom( this->shared_from_this(), broadcast );
        }
#else
        auto element = ( !apply_log )?
            parameterspace_type::random( this->shared_from_this(), broadcast ) :
            parameterspace_type::logRandom( this->shared_from_this(), broadcast );

#endif
        return element;
    }
    element_ptrtype elementPtr()  { element_ptrtype e( new element_type( this->shared_from_this() ) ); *e=element(); return e; }
    bool check() const
        {
            if ( this->dimension() == 0 )
                return true;
            return (min()-max()).array().abs().maxCoeff() > 1e-10;
        }
    /**
     * \class Sampling
     * \brief Parameter space sampling class
     */
    class Sampling : public std::vector<Element>
    {
        typedef std::vector<Element> super;
    public:

        typedef Sampling sampling_type;
        typedef std::shared_ptr<sampling_type> sampling_ptrtype;

        typedef ParameterSpace<Dimension> parameterspace_type;
        typedef std::shared_ptr<parameterspace_type> parameterspace_ptrtype;

        typedef typename parameterspace_type::Element element_type;
        typedef std::shared_ptr<element_type> element_ptrtype;

#if defined(FEELPP_HAS_ANN_H)
        typedef ANNkd_tree kdtree_type;
        typedef std::shared_ptr<kdtree_type> kdtree_ptrtype;
#endif /* FEELPP_HAS_ANN_H */

    public:
        Sampling( parameterspace_ptrtype const& space, int N = 0, sampling_ptrtype const& supersampling = sampling_ptrtype() )
            :
            super( N ),
            M_space( space ),
            M_supersampling( supersampling )
#if defined( FEELPP_HAS_ANN_H )
            ,M_kdtree()
#endif
            {}

        /**
         * @brief get the parameter space
         * 
         * @return the parameter space
         */
        auto parameterSpace() const { return M_space; }

        /**
         * \brief return number of elements in the sampling
         */
        int nbElements() const
        {
            return super::size();
        }

        /**
         * \brief return the last element in the sampling
         */
        element_type lastElement()
        {
            return super::back();
        }

        /**
         * \brief create a sampling with elements given by the user
         * \param V : vector of element_type
         */
        Sampling& setElements( std::vector< element_type > const& V )
        {
            super::operator=( V );
            return *this;
        }

        /**
         * \brief create add an element to a sampling
         * \param mu : element_type
         */
        Sampling& addElement( element_type const& mu )
        {
            super::push_back( mu );
            return *this;
        }

        void distributeOnAllProcessors( int N , std::string const& file_name )
        {
            auto const& theworldcomm = (M_space)? M_space->worldComm() : Environment::worldComm();
            rank_type total_proc = theworldcomm.globalSize();
            rank_type proc = theworldcomm.globalRank();
            rank_type master_proc = theworldcomm.masterRank();
            int number_of_elements_per_proc =N/total_proc;
            int extra_elements = N%total_proc;

            LOG( INFO ) <<"total_proc : "<<total_proc;
            LOG( INFO ) <<"proc : "<<proc;
            LOG( INFO ) <<"total_number_of_elements : "<<N;
            LOG( INFO ) <<"number_of_elements_per_proc : "<<number_of_elements_per_proc;
            LOG( INFO ) <<"extra_elements : "<<extra_elements;

            int total_number_of_elements_per_proc=0;
            int proc_rcv_sampling=0;

            std::vector<boost::mpi::request> reqs;
            if ( total_proc > 1 )
                reqs.reserve( theworldcomm.isMasterRank()? total_proc -1 : 1 );

            int tag=0;
            int shift=0;


            if( theworldcomm.isMasterRank() )
            {
                if ( !file_name.empty() )
                {
                    std::ifstream file ( file_name );
                    CHECK( file ) << "The file "<<file_name<<" was not found so we can't distribute the sampling on all processors\n";
                    this->readFromFile( file_name );
                }

                std::vector<Element> temporary_sampling;
                for(int i=0; i<N; i++)
                {
                    temporary_sampling.push_back( super::at( i ) );
                }
                super::clear();


                for(int proc_recv_sampling=0; proc_recv_sampling<total_proc; proc_recv_sampling++)
                {
                    if( proc_recv_sampling < extra_elements )
                    {
                        total_number_of_elements_per_proc=number_of_elements_per_proc+1;
                        shift=proc_recv_sampling*total_number_of_elements_per_proc;
                    }
                    else
                    {
                        int e=number_of_elements_per_proc;
                        total_number_of_elements_per_proc=e;
                        shift=extra_elements*(e+1) + (proc_recv_sampling-extra_elements)*e;
                    }

                    if( proc_recv_sampling != master_proc )
                    {
                        for(int local_element=0; local_element<total_number_of_elements_per_proc;local_element++)
                        {
                            int idx=local_element+shift;
                            super::push_back( temporary_sampling[idx] );
                        }
                    }

                    if( total_proc > 1 )
                    {
                        if( proc_recv_sampling != master_proc )
                        {
                            reqs.push_back( theworldcomm.isend( proc_recv_sampling , tag, /**this*/boost::serialization::base_object<super>( *this ) ) );
                            super::clear();
                        }
                    }//total_proc > 1

                }//proc_recv_sampling

                //the master proc fill its sampling in last
                //because sampling for others are deleted via super::clear
                if( master_proc < extra_elements )
                {
                    total_number_of_elements_per_proc=number_of_elements_per_proc+1;
                    shift=master_proc*total_number_of_elements_per_proc;
                }
                else
                {
                    int e=number_of_elements_per_proc;
                    total_number_of_elements_per_proc=e;
                    shift=extra_elements*(e+1) + (master_proc-extra_elements)*e;
                }
                for(int local_element=0; local_element<total_number_of_elements_per_proc;local_element++)
                {
                    int idx=local_element+shift;
                    super::push_back( temporary_sampling[idx] );
                }

            }//master proc

            if( proc != master_proc )
            {
                reqs.push_back( theworldcomm.irecv( master_proc, tag, /**this*/boost::serialization::base_object<super>( *this )) );
            }//not master proc

            boost::mpi::wait_all(reqs.begin(), reqs.end());

        }

        void sampling( size_type N, std::string const& samplingMode ) { return sample( N, samplingMode ); }

        /**
         * \brief create a sampling with global number of samples
         * \param N the number of samples
         * \param samplingMode : random, log-random, log-equidistribute, equidistribute
         * \param all_procs_have_same_sampling (boolean)
         * \param filename : file name where the sampling is written
         */
        void sample( size_type N, std::string const& samplingMode, bool all_procs_have_same_sampling=true, std::string const& filename="" )
        {
            if( samplingMode == "random" )
                this->randomize( N , all_procs_have_same_sampling, "", false );
            else if( samplingMode == "log-random" )
                this->randomize( N , all_procs_have_same_sampling, "", true );
            else if( samplingMode == "log-equidistribute" )
                this->logEquidistribute( N , all_procs_have_same_sampling );
            else if( samplingMode == "equidistribute" )
                this->equidistribute( N , all_procs_have_same_sampling  );
            else
                CHECK( false ) << "invalid sampling-mode, please select between log-random, log-equidistribute or equidistribute";

            if ( !filename.empty() )
            {
                if ( !all_procs_have_same_sampling || M_space->worldComm().isMasterRank() )
                    this->writeOnFile(filename);
                if ( all_procs_have_same_sampling )
                    M_space->worldComm().barrier();
            }
        }

        /**
         * \brief create a sampling with number of samples in each direction
         * \param N the number of samples in each direction
         * \param samplingMode : random, log-random, log-equidistribute, equidistribute
         * \param all_procs_have_same_sampling (boolean)
         * \param filename : file name where the sampling is written
         */
        void sample( std::vector<size_type> const& N, std::string const& samplingMode, bool all_procs_have_same_sampling=true, std::string const& filename="" )
        {
            if( samplingMode == "log-equidistribute" )
                this->logEquidistributeProduct( N, all_procs_have_same_sampling );
            else if( samplingMode == "equidistribute" )
                this->equidistributeProduct( N, all_procs_have_same_sampling  );
            else
                CHECK( false ) << "invalid sampling-mode, log-equidistribute or equidistribute";

            if ( !filename.empty() )
            {
                if ( !all_procs_have_same_sampling || M_space->worldComm().isMasterRank() )
                    this->writeOnFile(filename);
                if ( all_procs_have_same_sampling )
                    M_space->worldComm().barrier();
            }

        }

        /**
         * \brief create a sampling with random elements
         * \param N the number of samples
         * \param all_procs_have_same_sampling (boolean)
         * \param file_name : file name where the sampling is written
         */
        void randomize( int N , bool all_procs_have_same_sampling=true, std::string const& file_name="", bool apply_log=false )
            {
                CHECK( M_space ) << "Invalid(null pointer) parameter space for parameter generation\n";

                //std::srand(static_cast<unsigned>(std::time(0)));

                this->clear();

                // fill with log Random elements from the parameter space
                //only with one proc and then send a part of the sampling to each other proc
                if( M_space->worldComm().isMasterRank() /*&& generate_the_file*/ )
                {

                    bool already_exist;
                    // first element
                    auto mu = ( apply_log )?parameterspace_type::logRandom( M_space , false ) : parameterspace_type::random( M_space, false );
                    super::push_back( mu );

                    for(int i=1; i<N; i++)
                    {
                        //while mu is already in temporary_sampling
                        //we pick an other parameter
                        do
                        {
                            already_exist=false;
                            if( apply_log )
                            {
                                mu = parameterspace_type::logRandom( M_space, false );
                            }
                            else
                            {
                                mu = parameterspace_type::random( M_space, false );
                            }

                            for( auto const& _mu : *this )
                            {
                                if( mu == _mu )
                                    already_exist=true;
                            }

                        }
                        while( already_exist );

                        super::push_back( mu );

                    }//loop over N

                    this->writeOnFile( file_name );
                }//master proc

                if( all_procs_have_same_sampling )
                {
                    boost::mpi::broadcast( M_space->worldComm() , /**this*/boost::serialization::base_object<super>( *this ) , M_space->worldComm().masterRank() );
                    for(auto& _mu : *this )
                        _mu.setParameterSpace(M_space);
                }
                else
                {
                    this->distributeOnAllProcessors( N , file_name );
                }
            }

        /**
         * \brief create a sampling with log-equidistributed elements
         * \param N the number of samples
         * \param all_procs_have_same_sampling (boolean)
         * \param file_name : file name where the sampling is written
         */
        void logEquidistribute( int N , bool all_procs_have_same_sampling=true,  std::string const& file_name="" )
            {
                this->clear();

                uint16_type nDim = M_space->dimension();

                // define sampling size in each direction (search to be a global sampling size close to N)
                int samplingSizeDirectionEquidistribute = math::floor( math::pow( double(N) , 1./nDim ) );
                std::vector<size_type> samplingSizeDirection( nDim, samplingSizeDirectionEquidistribute );
                for ( uint16_type d = 0; d < nDim;++d )
                {
                    ++samplingSizeDirection[d];
                    size_type newSamplingSize = std::accumulate(samplingSizeDirection.begin(), samplingSizeDirection.end(), 1, std::multiplies<size_type>());
                    if ( newSamplingSize > N )
                    {
                        --samplingSizeDirection[d];
                        break;
                    }
                }
                size_type samplingSize = std::accumulate(samplingSizeDirection.begin(), samplingSizeDirection.end(), 1, std::multiplies<size_type>());

                if( M_space->worldComm().isMasterRank() )
                {
                    this->genericEquidistributeImpl( samplingSizeDirection, 1 );

                    this->writeOnFile( file_name );
                }
                if( all_procs_have_same_sampling )
                {
                    boost::mpi::broadcast( M_space->worldComm() , /**this*/boost::serialization::base_object<super>( *this ) , M_space->worldComm().masterRank() );
                }
                else
                {
                    this->distributeOnAllProcessors( samplingSize , file_name );
                }//if all procs don't have  same sampling

            }

        /**
         * \brief create a sampling with log-equidistributed elements
         * \param N : vector containing the number of samples on each direction
         * if N[direction] < 1 then we take minimum value for this direction
         */
        void logEquidistributeProduct( std::vector<size_type> const&  N, bool all_procs_have_same_sampling=true, std::string const& file_name="" )
        {
            // clear previous sampling
            this->clear();

            size_type samplingSize = std::accumulate(N.begin(), N.end(), 1, std::multiplies<size_type>());

            if( M_space->worldComm().isMasterRank() )
            {
                this->genericEquidistributeImpl( N, 1 );

                this->writeOnFile( file_name );

            }//end of master proc

            if( all_procs_have_same_sampling )
            {
                boost::mpi::broadcast( M_space->worldComm() , /**this*/boost::serialization::base_object<super>( *this ) , M_space->worldComm().masterRank() );
            }
            else
            {
                this->distributeOnAllProcessors( samplingSize , file_name );
            }
        }


        /**
         * \brief create a sampling with log-equidistributed and equidistributed elements
         * \param Nlogequi : vector containing the number of log-equidistributed samples on each direction
         * \param Nequi : vector containing the number of equidistributed samples on each direction
         */
        void mixEquiLogEquidistributeProduct( std::vector<int> const& Nlogequi , std::vector<int> const& Nequi )
        {
            // first empty the set
            this->clear();

            if( M_space->worldComm().isMasterRank() )
            {

                int number_of_directions_equi = Nequi.size();
                int number_of_directions_logequi = Nlogequi.size();
                FEELPP_ASSERT( number_of_directions_equi == number_of_directions_logequi )( number_of_directions_equi )( number_of_directions_logequi )
                    .error( "incompatible number of directions, you don't manipulate the same parameter space for log-equidistributed and equidistributed sampling" );

                //contains values of parameters on each direction
                std::vector< std::vector< double > > components_equi;
                components_equi.resize( number_of_directions_equi );
                std::vector< std::vector< double > > components_logequi;
                components_logequi.resize( number_of_directions_logequi );

                for(int direction=0; direction<number_of_directions_equi; direction++)
                {
                    std::vector<double> coeff_vector_equi = parameterspace_type::equidistributeInDirection( M_space , direction , Nequi[direction] );
                    components_equi[direction]= coeff_vector_equi ;
                    std::vector<double> coeff_vector_logequi = parameterspace_type::logEquidistributeInDirection( M_space , direction , Nlogequi[direction] );
                    components_logequi[direction]= coeff_vector_logequi ;
                }
                std::vector< std::vector< std::vector<double> > > vector_components(2);
                vector_components[0]=components_logequi;
                vector_components[1]=components_equi;
                generateElementsProduct( vector_components );

            }//end of master proc
            boost::mpi::broadcast( M_space->worldComm() , *this , M_space->worldComm().masterRank() );
        }


        /*
         * \param a vector of vector containing values of parameters on each direction
         * example if components has 2 vectors : [1,2,3] and [100,200] then it creates the samping :
         * (1,100) - (2,100) - (3,100) - (1,200) - (2,200) - (3,200)
         */
        void generateElementsProduct( std::vector< std::vector< double > > const& components )
        {
            std::vector< std::vector< std::vector< double > > > vector_components(1);
            vector_components[0]=components;
            generateElementsProduct( vector_components );
        }

        void generateElementsProduct( std::vector< std::vector< std::vector< double > > > const& vector_components )
        {
            int number_of_comp = vector_components.size();
            for(int comp=0; comp<number_of_comp; comp++)
            {
                int number_of_directions=vector_components[comp].size();
                element_type mu( M_space );

                //initialization
                std::vector<int> idx( number_of_directions );
                for(int direction=0; direction<number_of_directions; direction++)
                    idx[direction]=0;

                int last_direction=number_of_directions-1;

                bool stop=false;
                auto min=M_space->min();
                auto max=M_space->max();

                do
                {
                    bool already_exist = true;
                    //construction of the parameter
                    while( already_exist && !stop )
                    {
                        already_exist=false;
                        for(int direction=0; direction<number_of_directions; direction++)
                        {
                            int index = idx[direction];
                            double value = vector_components[comp][direction][index];
                            mu(direction) = value ;
                        }

                        if( mu == min  && comp>0 )
                            already_exist=true;
                        if( mu == max && comp>0 )
                            already_exist=true;

                        //update values or stop
                        for(int direction=0; direction<number_of_directions; direction++)
                        {
                            if( idx[direction] < vector_components[comp][direction].size()-1 )
                            {
                                idx[direction]++;
                                break;
                            }
                            else
                            {
                                idx[direction]=0;
                                if( direction == last_direction )
                                    stop = true;
                            }
                        }//end loop over directions

                    }//while


                    //add the new parameter
                    super::push_back( mu );

                } while( ! stop );

            }//end loop on vector_component
        }


        /**
         * \brief write the sampling in a file
         * \param file_name : name of the file to read
         * in the file we write :
         * mu_0= [ value0 , value1 , ... ]
         * mu_1= [ value0 , value1 , ... ]
         */
        void writeOnFile( std::string const& file_name = "list_of_parameters_taken" )
            {
                if ( file_name.empty() )
                    return;
                //auto const& theworldcomm = (M_space)? M_space->worldComm() : Environment::worldComm();
                //rank_type proc = theworldcomm.globalRank();
                //rank_type total_proc = theworldcomm.globalSize();
                //std::string real_file_name = (boost::format(file_name+"-proc%1%on%2%") %proc %total_proc ).str();
                std::ofstream file;
                file.open( file_name,std::ios::out );
                //element_type mu( M_space );
                //int size = mu.size();
                int size = M_space->dimension();
                int number = 0;
                for( auto const& mu : *this )
                {
                    file << std::setprecision(15) << " mu_"<<number<<"= [ ";
                    for( int i=0; i<size-1; i++ )
                        file << mu[i]<<" , ";
                    if ( size > 0 )
                        file << mu[size-1];
                    file << " ] \n" ;
                    number++;
                }
                file.close();
            }

        /**
         * \brief read the sampling from a file
         * \param file_name : name of the file to read
         * in the file we expect :
         * mu_0= [ value0 , value1 , ... ]
         * mu_1= [ value0 , value1 , ... ]
         * return the size of the sampling
         */
        double readFromFile( std::string const& file_name= "list_of_parameters_taken" )
            {
                std::ifstream file ( file_name );
                double mui;
                std::string str;
                int number=0;
                file>>str;
                while( ! file.eof() )
                {
                    element_type mu( M_space );
                    file>>str;
                    int i=0;
                    while ( str!="]" )
                    {
                        file >> mui;
                        mu[i] = mui;
                        file >> str;
                        i++;
                    }
                    super::push_back( mu );
                    number++;
                    file>>str;
                }
                file.close();
                return number;
            }


        /**
         * build the closest sampling with parameters given from the file
         * look in the supersampling closest parameters
         * \param file_name : give the real parameters we want
         * return the index vector of parameters in the supersampling
         */
        std::vector<int> closestSamplingFromFile( std::string const& file_name="list_of_parameters_taken")
        {
            std::vector<int> index_vector;
            std::ifstream file( file_name );
            double mui;
            std::string str;
            int number=0;
            file>>str;
            while( ! file.eof() )
            {
                element_type mu( M_space );
                file>>str;
                int i=0;
                while ( str!="]" )
                {
                    file >> mui;
                    mu[i] = mui;
                    file >> str;
                    i++;
                }
                //search the closest neighbor of mu in the supersampling
                std::vector<int> local_index_vector;
                auto neighbors = M_supersampling->searchNearestNeighbors( mu, 1 , local_index_vector );
                auto closest_mu = neighbors->at(0);
                int index = local_index_vector[0];
                this->push_back( closest_mu , index );
                number++;
                file>>str;
                index_vector.push_back( index );
            }
            file.close();
            return index_vector;
        }

        /**
         * \brief create a sampling with equidistributed elements
         * \param N the number of samples
         * \param all_procs_have_same_sampling (boolean)
         * \param file_name : file name where the sampling is written
         */
        void equidistribute( int N , bool all_procs_have_same_sampling=true, std::string const& file_name="" )
            {
                this->clear();

                uint16_type nDim = M_space->dimension();

                // define sampling size in each direction (search to be a global sampling size close to N)
                int samplingSizeDirectionEquidistribute = math::floor( math::pow( double(N) , 1./nDim ) );
                std::vector<size_type> samplingSizeDirection( nDim, samplingSizeDirectionEquidistribute );
                for ( uint16_type d = 0; d < nDim;++d )
                {
                    ++samplingSizeDirection[d];
                    size_type newSamplingSize = std::accumulate(samplingSizeDirection.begin(), samplingSizeDirection.end(), 1, std::multiplies<size_type>());
                    if ( newSamplingSize > N )
                    {
                        --samplingSizeDirection[d];
                        break;
                    }
                }
                size_type samplingSize = std::accumulate(samplingSizeDirection.begin(), samplingSizeDirection.end(), 1, std::multiplies<size_type>());

                if( M_space->worldComm().isMasterRank() )
                {
                    this->genericEquidistributeImpl( samplingSizeDirection, 0 );

                    this->writeOnFile( file_name );
                }//master proc

                if( all_procs_have_same_sampling )
                {
                    boost::mpi::broadcast( M_space->worldComm() , /**this*/boost::serialization::base_object<super>( *this ) , M_space->worldComm().masterRank() );
                }
                else
                {
                    this->distributeOnAllProcessors( samplingSize , file_name );
                }//if all procs don't have  same sampling
            }

        /**
         * \brief create a sampling with equidistributed elements
         * \param N : vector containing the number of samples on each direction
         */
        void equidistributeProduct( std::vector<size_type> const&  N, bool all_procs_have_same_sampling=true, std::string const& file_name="" )
        {
            // clear previous sampling
            this->clear();

            size_type samplingSize = std::accumulate(N.begin(), N.end(), 1, std::multiplies<size_type>());

            if( M_space->worldComm().isMasterRank() )
            {
                this->genericEquidistributeImpl( N, 0 );

                this->writeOnFile( file_name );

            }//end of master proc

            if( all_procs_have_same_sampling )
            {
                boost::mpi::broadcast( M_space->worldComm() , /**this*/boost::serialization::base_object<super>( *this ) , M_space->worldComm().masterRank() );
            }
            else
            {
                this->distributeOnAllProcessors( samplingSize , file_name );
            }
        }

        /**
         * \brief create a sampling with random elements
         * \param N : vector of bool for the direction to sample
         */
        void randomInDirections( int N,  std::vector<bool> const& Nd, bool all_procs_have_same_sampling=true, std::string const& file_name="" )
        {
            // clear previous sampling
            this->clear();

            size_type samplingSize = N;

            if( M_space->worldComm().isMasterRank() /*&& generate_the_file*/ )
            {

                bool already_exist;
                // first element
                auto mu = parameterspace_type::random( M_space, Nd );
                super::push_back( mu );

                for(int i=1; i<N; i++)
                {
                    //while mu is already in temporary_sampling
                    //we pick an other parameter
                    do
                    {
                        already_exist=false;
                        mu = parameterspace_type::random( M_space, Nd );

                        for( auto const& _mu : *this )
                        {
                            if( mu == _mu )
                                already_exist=true;
                        }

                    }
                    while( already_exist );

                    super::push_back( mu );

                }//loop over N

                this->writeOnFile( file_name );
            }//master proc

            if( all_procs_have_same_sampling )
            {
                boost::mpi::broadcast( M_space->worldComm() , /**this*/boost::serialization::base_object<super>( *this ) , M_space->worldComm().masterRank() );
            }
            else
            {
                this->distributeOnAllProcessors( N , file_name );
            }
        }


        /**
         * \brief Returns the minimum element in the sampling and its index
         */
        boost::tuple<element_type,size_type>
        min( bool check=true ) const
            {
                element_type mumin( M_space );
                //mumin = M_space->max();

                //element_type mu( M_space );
                int index = 0;
                int i = 0;
                double mumin_norm = 1e30;//mumin.norm();

                for( auto const& mu : *this )
                {
                    double curNorm = mu.norm();
                    if ( curNorm < mumin_norm  )
                    {
                        mumin = mu;
                        index = i;
                        mumin_norm = curNorm;
                    }

                    ++i;
                }

                //do communications to have global min
                rank_type world_size = M_space->worldComm().globalSize();
                std::vector<double> max_world( world_size );
                mpi::all_gather( M_space->worldComm().globalComm(),
                                 mumin_norm,
                                 max_world );
                auto it_min = std::min_element( max_world.begin() , max_world.end() );
                int proc_having_good_mu = it_min - max_world.begin();

                auto tuple = boost::make_tuple( mumin , index );
                boost::mpi::broadcast( M_space->worldComm() , tuple , proc_having_good_mu );

                mumin = tuple.template get<0>();
                index = tuple.template get<1>();
                mumin.setParameterSpace( M_space );

                if( check )
                    mumin.check();

                return boost::make_tuple( mumin, index );
            }

        /**
         * \brief Returns the maximum element in the sampling and its index
         */
        boost::tuple<element_type,size_type>
        max( bool check=true ) const
            {
                element_type mumax( M_space );
                //mumax = M_space->min();

                //element_type mu( M_space );
                int index = 0;
                int i = 0;
                double mumax_norm = -1e30;//mumax.norm();
                for( auto const& mu : *this )
                {
                    if ( mu.norm() > mumax_norm  )
                    {
                        mumax = mu;
                        index = i;
                        mumax_norm = mumax.norm();
                    }

                    ++i;
                }

                //do communications to have global min
                rank_type world_size = M_space->worldComm().globalSize();
                std::vector<double> max_world( world_size );
                mpi::all_gather( M_space->worldComm().globalComm(),
                                 mumax_norm,
                                 max_world );
                auto it_max = std::max_element( max_world.begin() , max_world.end() );
                int proc_having_good_mu = it_max - max_world.begin();
                auto tuple = boost::make_tuple( mumax , index );
                boost::mpi::broadcast( M_space->worldComm() , tuple , proc_having_good_mu );

                mumax = tuple.template get<0>();
                index = tuple.template get<1>();
                mumax.setParameterSpace( M_space );

                if( check )
                    mumax.check();

                return boost::make_tuple( mumax, index );
            }
        /**
         * \breif Returns the supersampling (might be null)
         */
        sampling_ptrtype const& superSampling() const
            {
                return M_supersampling;
            }

        /**
         * set the super sampling
         */
        void setSuperSampling( sampling_ptrtype const& super )
            {
                M_supersampling = super;
            }

        /**
         * \brief Returns the \p M closest points to \p mu in sampling
         * \param mu the point in parameter whom we want to find the neighbors
         * \param M the number of neighbors to find
         * \return vector of neighbors and index of them in the sampling
         */
        sampling_ptrtype searchNearestNeighbors( element_type const& mu, size_type M , std::vector<int>& index_vector, bool broadcast=true  );

        /**
         * \brief if supersampling is != 0, Returns the complement
         */
        sampling_ptrtype complement() const;

        /**
         * \brief add new parameter \p mu in sampling and store \p index in super sampling
         * 
         * \param mu the new parameter
         * \param index the index in the super sampling if the super sampling is not null
         */
        Sampling& push_back( element_type const& mu, size_type index = invalid_v<size_type> )
            {
                if ( M_supersampling && index != invalid_v<size_type> )
                    M_superindices.push_back( index );

                super::push_back( mu );
                return *this;
            }
        
        /**
         * @brief set the sampling from a vector of vector of values
         * 
         * @param data data to set the sampling
         * @return Sampling& 
         */
        Sampling& set( std::vector<std::vector<value_type>> const& data )
        {
            this->clear();
            for (const auto &elem_data : data)
            {
                this->emplace_back(Element(M_space).set(elem_data));
            }
            return *this;
        }

        /**
         * @brief add new parameter \p mu in sampling 
         * 
         * @param data 
         * @return Sampling& 
         */
        Sampling& add( std::vector<std::vector<value_type>> const& data )
        {
            for (const auto &elem_data : data)
            {
                this->emplace_back(Element(M_space).set(elem_data));
            }
            return *this;
        }

        /**
         * @brief set the sampling from a vector of map of pair of name and values
         * 
         * @param data data to set the sampling
         * @return Sampling& 
         */
        Sampling& set( std::vector<std::map<std::string,value_type>> const& data )
        {
            this->clear();
            for (const auto &elem_data : data)
            {
                this->emplace_back(Element(M_space).set(elem_data));
            }
            return *this;
        }

        /**
         * @brief add new parameter \p mu in sampling 
         * 
         * @param data 
         * @return Sampling& 
         */
        Sampling& add( std::vector<std::map<std::string,value_type>> const& data )
        {
            for (const auto &elem_data : data)
            {
                this->emplace_back(Element(M_space).set(elem_data));
            }
            return *this;
        }

        /**
         * \brief given a local index, returns the index in the super sampling
         * \param index index in the local sampling
         * \return the index in the super sampling
         */
        size_type indexInSuperSampling( size_type index ) const
            {
                return M_superindices[ index ];
            }

        /**
         * \brief save sampling data in json
         * \param output filename
         */
        void saveJson( std::string const& filename, bool addParameterSpaceInfo = true )
            {
                if ( M_space->worldComm().isMasterRank() )
                {
                    boost::property_tree::ptree ptree;
                    this->updatePropertyTree( ptree, addParameterSpaceInfo );
                    write_json( filename, ptree/*, std::locale(), false*/ );
                }
                M_space->worldComm().barrier();
            }
        /**
         * \brief update sampling description in property_tree data structure
         * \param ptree to update
         */
        void updatePropertyTree( boost::property_tree::ptree & ptree, bool addParameterSpaceInfo = true )
            {
                // update parameter space ptree
                if ( addParameterSpaceInfo )
                {
                    boost::property_tree::ptree ptreeParameterSpace;
                    M_space->updatePropertyTree( ptreeParameterSpace );
                    ptree.add_child( "parameter_space", ptreeParameterSpace );
                }
                // update sampling ptree
                int nDim = M_space->dimension();
                boost::property_tree::ptree ptreeSampling;
                for ( auto const& mu : *this )
                {
                    boost::property_tree::ptree ptreeMu;
                    for ( int d=0;d<nDim;++d )
                    {
                        // ptreeMu.add( "", mu(d) );
                        // the line commented above fails in debug, need to be replaced by the code below :
                        boost::property_tree::ptree ptreeTmp;
                        ptreeTmp.put( "", mu(d) );
                        ptreeMu.push_back( std::make_pair("", ptreeTmp) );
                    }
                    // ptreeSampling.add_child("",ptreeMu);
                    // the line commented above fails in debug, need to be replaced by the code below :
                    ptreeSampling.push_back( std::make_pair( "", ptreeMu ) );
                }
                ptree.add_child( "sampling", ptreeSampling );
            }

    private:

        Sampling() {}

        void genericEquidistributeImpl( std::vector<size_type> const& samplingSizeDirection, int type )
            {
                this->clear();

                uint16_type nDim = M_space->dimension();
                CHECK( samplingSizeDirection.size() == nDim ) << "invalid size, must be equal " << samplingSizeDirection.size() << " vs " << nDim;
                size_type samplingSize = std::accumulate(samplingSizeDirection.begin(), samplingSizeDirection.end(), 1, std::multiplies<size_type>());
                //std::cout << "samplingSize  : " << samplingSize << "\n";

                // compute sampling in each direction
                typename parameterspace_type::element_type mymu( M_space );

                std::vector<std::vector<double> > distribution( nDim );
                for ( uint16_type d = 0; d < nDim;++d )
                {
                    if ( type == 0 )
                        distribution[d] = parameterspace_type::equidistributeInDirection( M_space,d,samplingSizeDirection[d] );
                    else if ( type == 1 )
                        distribution[d] = parameterspace_type::logEquidistributeInDirection( M_space,d,samplingSizeDirection[d] );
                    mymu(d) = distribution[d][0];
                }

                // store information for next process
                std::vector<size_type> reuseParam( nDim );
                reuseParam[nDim-1] = 1;
                for ( int d = nDim-2; d >= 0;--d )
                {
                    reuseParam[d] = samplingSizeDirection[d+1]*reuseParam[d+1];
                }

                size_type samplingSize2 = samplingSizeDirection[0]*reuseParam[0];
                CHECK( samplingSize == samplingSize2 ) << "invalid data, must be equal : " << samplingSize << " vs " << samplingSize2;

                // update sampling
                std::vector<size_type> currentParam( nDim,0 );
                std::vector<size_type> countReuseParam( nDim,0 );
                super::resize( samplingSize, mymu );
                for ( size_type k=0;k<samplingSize;++k )
                {
                    super::operator[](k) = mymu;

                    for ( uint16_type d = 0; d < (nDim);++d )
                    {
                        ++countReuseParam[d];
                        if ( countReuseParam[d] == reuseParam[d] )
                        {
                            countReuseParam[d] = 0;
                            ++currentParam[d];
                            if ( currentParam[d] >= samplingSizeDirection[d] )
                                currentParam[d] = 0;
                            mymu(d) = distribution[d][currentParam[d]];
                        }
                    }
                }

            }


        friend class boost::serialization::access;
        template<class Archive>
        void save( Archive & ar, const unsigned int version ) const
            {
                ar & boost::serialization::base_object<super>( *this );
                ar & BOOST_SERIALIZATION_NVP( M_space );
                bool hasSupersampling = ( M_supersampling )? true : false;
                ar & hasSupersampling;
                if ( hasSupersampling )
                    ar & BOOST_SERIALIZATION_NVP( M_supersampling );
                ar & BOOST_SERIALIZATION_NVP( M_superindices );
            }

        template<class Archive>
        void load( Archive & ar, const unsigned int version )
            {
                ar & boost::serialization::base_object<super>( *this );
                ar & BOOST_SERIALIZATION_NVP( M_space );
                bool hasSupersampling = false;
                ar & hasSupersampling;
                if ( hasSupersampling )
                    ar & BOOST_SERIALIZATION_NVP( M_supersampling );
                else
                    M_supersampling.reset();
                ar & BOOST_SERIALIZATION_NVP( M_superindices );
            }
        BOOST_SERIALIZATION_SPLIT_MEMBER()
        private:
        parameterspace_ptrtype M_space;

        sampling_ptrtype M_supersampling;
        std::vector<size_type> M_superindices;
#if defined( FEELPP_HAS_ANN_H )
        kdtree_ptrtype M_kdtree;
#endif
    };

    typedef Sampling sampling_type;
    typedef std::shared_ptr<sampling_type> sampling_ptrtype;

    sampling_ptrtype sampling() { return sampling_ptrtype( new sampling_type( this->shared_from_this() ) ); }

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    ParameterSpace( uint16_type dim = 0, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        :
        super( worldComm ),
        M_nDim( (Dimension == invalid_uint16_type_value)? dim : Dimension ),
        M_min(),
        M_max()
        // M_mubar()
        {
            M_min.resize( M_nDim,1 );
            M_min.setZero();
            M_max.resize( M_nDim,1 );
            M_max.setOnes();
            this->updateDefaultParameterNames();
        }
    //! copy constructor
    ParameterSpace( ParameterSpace const & o ) = default;
#if 0
    ParameterSpace( element_type const& emin, element_type const& emax )
        :
        M_min( emin ),
        M_max( emax )
        {
            //M_min.setParameterSpace( this->shared_from_this() );
            //M_max.setParameterSpace( this->shared_from_this() );
        }
#endif
    //! constructor from ModelProperties
    ParameterSpace( CRBModelParameters const& modelParameters, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        :
        super( worldComm ),
        M_nDim(),
        M_min(),
        M_max()
    {
        this->setDimension( modelParameters.size() );

        for( auto const& [index,parameterPair] : enumerate( modelParameters ) )
        {
            setParameterName( index, parameterPair.second.name() );
            M_min( index ) = parameterPair.second.min();
            M_max( index ) = parameterPair.second.max();

            if (!( parameterPair.second.min() < parameterPair.second.max() ))
                throw std::logic_error( fmt::format("invalid range parameter {}\n",parameterPair.second.name()) );
        }

    }

    //! destructor
    ~ParameterSpace() override = default;

    /**
     * generate a shared_ptr out of a parameter space
     */
    static parameterspace_ptrtype create( uint16_type dim)
        {
            return New( dim, Environment::worldCommPtr() );
        }
    static parameterspace_ptrtype New( uint16_type dim = 0, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr())
        {
            return std::make_shared<parameterspace_type>( dim,worldComm );
        }

    static parameterspace_ptrtype New( std::string const& filename, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        {
            auto ps = std::make_shared<parameterspace_type>( 0,worldComm );
            ps->loadJson( filename );
            return ps;
        }
    static parameterspace_ptrtype New( CRBModelParameters const& modelParameters, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr())
        {
	        auto ps = std::make_shared<parameterspace_type>( modelParameters, worldComm );
    	    ps->setSpaces();
	        return ps;
	    }
    void setSpaces()
        {
            M_min.setParameterSpace(this->shared_from_this() );
            M_max.setParameterSpace(this->shared_from_this() );
        }
    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    ParameterSpace& operator=( ParameterSpace const & o ) = default;
    //@}

    /** @name Accessors
     */
    //@{

    //! \return the parameter space dimension
    uint16_type dimension() const
    {
        return M_nDim;//Dimension;
    }

    /**
     * return the minimum element
     */
    element_type const& min() const
    {
        return M_min;
    }

    /**
     * @brief return the minimum element in direction d
     * 
     * @param d direction
     * @return value_type min value
     */
    value_type min( uint16_type d ) const
    {
        return M_min(d);
    }

    /**
     * @brief return the minimum element of the parameter named name
     * 
     * @param name name pof the parameter
     * @return value_type 
     */
    value_type min( std::string const& name ) const
    {
        auto it = std::find(M_parameterNames.begin(), M_parameterNames.end(), name);
        if( it != M_parameterNames.end() )
            return M_min(it - M_parameterNames.begin());
        else
             std::invalid_argument( fmt::format("invalid parameter name {}",name) );
    }
    /**
     * return the maximum element
     */
    element_type const& max() const
    {
        return M_max;
    }

    /**
     * @brief return the maximum element of the parameter named name
     * 
     * @param name name pof the parameter
     * @return value_type 
     */
    value_type max( std::string const& name ) const
    {
        auto it = std::find(M_parameterNames.begin(), M_parameterNames.end(), name);
        if( it != M_parameterNames.end() )
            return M_max(it - M_parameterNames.begin());
        else
            std::invalid_argument( fmt::format("invalid parameter name {}",name) );
    }

    /**
     * @brief return the maximum element in direction d
     * 
     * @param d direction
     * @return value_type max value
     */
    value_type max( uint16_type d ) const
    {
        return M_max(d);
    }

    /**
     * \brief the log-middle point of the parameter space
     */
    element_type logMiddle() const
        {
            return ( ( M_min.array().log() + M_max.array().log() )/2. ).exp();
        }

    /**
     * \brief the middle point of the parameter space
     */
    element_type middle() const
        {
            return ( ( M_min + M_max )/2. );
        }

    /**
     * \brief name of a parameter
     * \param d parameter index
     */
    std::string const& parameterName( uint16_type d ) const
        {
            return M_parameterNames[d];
        }

    /**
     * \brief name of the parameters
     */
    std::vector<std::string> const& parameterNames() const
        {
            return M_parameterNames;
        }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the parameter space dimension
     */
    void setDimension( uint16_type d )
        {
            if ( M_nDim == d )
                return;
            if ( Dimension == invalid_uint16_type_value )
            {
                M_nDim = d;
                M_min.resize( M_nDim,1 );
                M_min.setZero();
                M_max.resize( M_nDim,1 );
                M_max.setOnes();
                this->updateDefaultParameterNames();
            }
        }

    /**
     * set the minimum element
     */
    void setMin( element_type const& min )
        {
            M_min = min;
            M_min.setParameterSpace( this->shared_from_this() );
        }

    /**
     * set the maximum element
     */
    void setMax( element_type const& max )
        {
            M_max = max;
            M_max.setParameterSpace( this->shared_from_this() );
        }



    /**
     * \brief name of a parameter
     * \param d parameter index
     */
    void setParameterName( uint16_type d, std::string const& name )
        {
            M_parameterNames[d] = name;
        }

    /**
     * get index of named parameter
     */
    int indexOfParameterNamed( std::string name ) const
        {
            auto it = std::find(M_parameterNames.begin(), M_parameterNames.end(), name);
            if( it != M_parameterNames.end() )
                return it - M_parameterNames.begin();
            else
                return -1;
        }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * \brief load ParameterSpace from json
     * \param input json filename
     */
    void loadJson( std::string const& filename )
        {
            if ( !fs::exists( filename ) )
            {
                LOG(INFO) << "Could not find " << filename << std::endl;
                return;
            }

            auto json_str_wo_comments = removeComments(readFromFile(filename));
            //LOG(INFO) << "json file without comment:" << json_str_wo_comments;

            boost::property_tree::ptree ptreeParameterSpace;
            std::istringstream istr( json_str_wo_comments );
            boost::property_tree::read_json( istr, ptreeParameterSpace );
            this->setup( ptreeParameterSpace );
        }

    /**
     * \brief setup ParameterSpace from property_tree::ptree
     * \param input property_tree::ptree
     */
    void setup( boost::property_tree::ptree const& ptreeParameterSpace )
        {
            try {
                int nDim  = ptreeParameterSpace.template get<int>( "dimension" );
                if ( Dimension == invalid_uint16_type_value )
                    this->setDimension(nDim);
                else
                    CHECK( Dimension == nDim && M_nDim == nDim ) << "invalid dimension in json " << nDim << " must be " << Dimension;
            }
            catch ( boost::property_tree::ptree_bad_path& e )
            {
                if ( Dimension == invalid_uint16_type_value )
                    CHECK( false ) << "invalid json : dimension is missing";
            }

            std::map<std::string,int> mapParameterNamesToId;
            M_parameterNames.clear();
            if ( M_nDim > 0 )
            {
                for (auto const& item : ptreeParameterSpace.get_child("names"))
                {
                    std::string const& name = item.second.template get_value<std::string>();
                    mapParameterNamesToId[name] = M_parameterNames.size();
                    M_parameterNames.push_back( name );
                }

                element_type mumin( this->shared_from_this() );
                element_type mumax( this->shared_from_this() );
                for ( auto const& ptreeParametersPair : ptreeParameterSpace.get_child("parameters") )
                {
                    std::string const& paramName = ptreeParametersPair.first;
                    auto itFindParameterNames = mapParameterNamesToId.find( paramName );
                    CHECK( itFindParameterNames != mapParameterNamesToId.end() ) << "parameter "<< paramName << " is not registered";
                    int paramId = itFindParameterNames->second;
                    auto const& ptreeParameters = ptreeParametersPair.second;
                    mumin( paramId ) = ptreeParameters.template get<double>("min");
                    mumax( paramId ) = ptreeParameters.template get<double>("max");
                }
                this->setMin( mumin );
                this->setMax( mumax );
            }
        }

        /**
         * \brief save ParameterSpace description in json
         * \param output filename
         */
    void saveJson( std::string const& filename ) const
        {
            if ( this->worldComm().isMasterRank() )
            {
                boost::property_tree::ptree ptree;
                this->updatePropertyTree( ptree );
                write_json( filename, ptree/*, std::locale(), false*/ );
            }
            this->worldComm().barrier();
        }
    /**
     * \brief update ParameterSpace description in property_tree data structure
     * \param ptree to update
     */
    void updatePropertyTree( boost::property_tree::ptree & ptreeParameterSpace ) const
        {
            int nDim = this->dimension();
            ptreeParameterSpace.add( "dimension",nDim);

            if ( nDim > 0 )
            {
                boost::property_tree::ptree ptreeParameterName;
                for ( int d=0;d<nDim;++d )
                {
                    boost::property_tree::ptree ptreeParameterName2;
                    ptreeParameterName2.put( "", this->parameterName(d) );
                    ptreeParameterName.push_back( std::make_pair("", ptreeParameterName2) );
                }
                ptreeParameterSpace.add_child( "names", ptreeParameterName );

                boost::property_tree::ptree ptreeParametersDesc;
                for ( int d=0;d<nDim;++d )
                {
                    boost::property_tree::ptree ptreeParameterDesc;
                    ptreeParameterDesc.add( "min", this->min()(d) );
                    ptreeParameterDesc.add( "max", this->max()(d) );
                    ptreeParametersDesc.add_child( this->parameterName(d), ptreeParameterDesc );
                }
                ptreeParameterSpace.add_child( "parameters", ptreeParametersDesc );
            }

        }

    /**
     * \brief Returns a log random element of the parameter space
     */
    static element_type logRandom( parameterspace_ptrtype const& space, bool broadcast )
        {
            //LOG(INFO) << "call logRandom...\n";
            //LOG(INFO) << "call logRandom broadcast: " << broadcast << "...\n";
            //google::FlushLogFiles(google::GLOG_INFO);
            if ( broadcast )
            {
                element_type mu( space );
                if ( space->dimension() == 0 )
                    return mu;

                if ( !space->check() ) return mu;
                if( space->worldComm().isMasterRank() )
                {
                    //LOG(INFO) << "generate random mu...\n";
                    //google::FlushLogFiles(google::GLOG_INFO);
                    mu = logRandom1( space );
                }
                //LOG(INFO) << "broadcast...\n";
                //google::FlushLogFiles(google::GLOG_INFO);
                boost::mpi::broadcast( space->worldComm() , mu , space->worldComm().masterRank() );
                //Environment::worldComm().barrier();
                //LOG(INFO) << "check...\n";
                mu.check();
                return mu;
            }
            else
            {
                return logRandom1( space );
            }
        }
    static element_type logRandom1( parameterspace_ptrtype const& space )
        {
            element_type mur( space );
            if ( space->dimension() == 0 )
                return mur;
            mur.array() = element_type::Random(space->dimension(),1).array().abs();
            //LOG(INFO) << "random1 generate random mur= " << mur << " \n";
#if 0
            google::FlushLogFiles(google::GLOG_INFO);
            element_type mu( space );
            mu.array() = ( space->min().array().log()+mur.array()*( space->max().array().log()-space->min().array().log() ) ).array().exp();
            //LOG(INFO) << "random1 generate random mu= " << mu << " \n";
#else
            element_type mu( space );
            element_type muShift( space );
            element_type muMin( space );
            element_type muMax( space );
            for (int d=0;d<space->dimension();++d)
            {
                muShift(d) = -space->min()(d) + 1.;
                muMin(d) = space->min()(d) + muShift(d);
                muMax(d) = space->max()(d) + muShift(d);
            }
            mu.array() = ( muMin.array().log()+mur.array()*( muMax.array().log()-muMin.array().log() ) ).exp();
            mu.array() -= muShift.array();

#endif
            return mu;
        }

    /**
     * \brief Returns a log random element of the parameter space
     */
    static element_type random( parameterspace_ptrtype const& space, bool broadcast = true )
        {
            //std::srand( static_cast<unsigned>( std::time( 0 ) ) );
            element_type mur( space );
            if ( space->dimension() == 0 )
                return mur;

            if ( broadcast )
            {
                if( space->worldComm().isMasterRank() )
                {
                    mur.array() = element_type::Random(space->dimension(),1).array().abs();
                }
                boost::mpi::broadcast( space->worldComm() , mur , space->worldComm().masterRank() );

            }
            else
                mur.array() = element_type::Random(space->dimension(),1).array().abs();
            //std::cout << "mur= " << mur << "\n";
            //mur.setRandom()/RAND_MAX;
            //std::cout << "mur= " << mur << "\n";
            element_type mu( space );
            mu.array() = space->min().array()+mur.array()*( space->max().array()-space->min().array() );
            if ( broadcast )
                mu.check();
            return mu;
        }

    /**
     * \brief Returns a random element of the parameter space int the directions N
     */
    static element_type random( parameterspace_ptrtype const& space, std::vector<bool> N )
        {
            element_type mur( space );
            if ( space->dimension() == 0 )
                return mur;
            mur.array() = element_type::Random(space->dimension(),1).array().abs();
            //LOG(INFO) << "random1 generate random mur= " << mur << " \n";

            element_type mu( space );
            element_type muMin = space->min();
            element_type muMax = space->max();
            for (int d=0;d<space->dimension();++d)
                mu(d) = muMin(d) + int(N[d])*mur(d)*( muMax(d)-muMin(d) );

            return mu;
        }

    /**
     * \brief Returns a vector representing the values of log equidistributed element of the parameter space
     * in the given direction
     * \param direction
     * \param N : number of samples in the direction
     */
    static std::vector<double> logEquidistributeInDirection( parameterspace_ptrtype const& space, int direction , int N)
    {
        std::vector<double> result(N);

        CHECK( direction >=0 && direction < space->dimension() ) << "direction invalid " << direction;
        double shift = -space->min()(direction) + 1.;
        double min = space->min()(direction) + shift;
        double max = space->max()(direction) + shift;
        if ( N > 1 )
        {
            for(int i=0; i<N; i++)
            {
                double factor = (double)(i)/(N-1);
                result[i] = math::exp(math::log(min)+factor*( math::log(max)-math::log(min) ) ) - shift;
            }
        }
        else if ( N == 1 )
        {
            result[0] = (min+max)/2. - shift;
        }

        return result;
    }


    /**
     * \brief Returns a log equidistributed element of the parameter space
     * \param factor is a factor in [0,1]
     */
    static element_type logEquidistributed( double factor, parameterspace_ptrtype const& space )
        {
#if 0
            element_type mu( space );
            mu.array() = ( space->min().array().log()+factor*( space->max().array().log()-space->min().array().log() ) ).exp();
            return mu;
#else
            element_type mu( space );
            element_type muShift( space );
            element_type muMin( space );
            element_type muMax( space );
            for (int d=0;d<space->dimension();++d)
            {
                muShift(d) = -space->min()(d) + 1.;
                muMin(d) = space->min()(d) + muShift(d);
                muMax(d) = space->max()(d) + muShift(d);
            }
            mu.array() = ( muMin.array().log()+factor*( muMax.array().log()-muMin.array().log() ) ).exp();
            mu.array() -= muShift.array();
            return mu;
#endif
        }


    /**
     * \brief Returns a vector representing the values of equidistributed element of the parameter space
     * in the given direction
     * \param direction
     * \param N : number of samples in the direction
     */
    static std::vector<double> equidistributeInDirection( parameterspace_ptrtype const& space, int direction , int N)
    {
        std::vector<double> result(N);

        auto min_element = space->min();
        auto max_element = space->max();

        int mu_size = min_element.size();
        if( (mu_size < direction) || (direction < 0)  )
            throw std::logic_error( "[ParameterSpace::equidistributeInDirection] ERROR : bad dimension of vector containing number of samples on each direction" );

        double min = min_element(direction);
        double max = max_element(direction);

        if ( N > 1 )
        {
            for(int i=0; i<N; i++)
            {
                double factor = (double)(i)/(N-1);
                result[i] = min+factor*( max-min );
            }
        }
        else if ( N == 1 )
        {
            result[0] = (min+max)/2.;
        }

        return result;
    }


    /**
     * \brief Returns a equidistributed element of the parameter space
     * \param factor is a factor in [0,1]
     */
    static element_type equidistributed( double factor, parameterspace_ptrtype const& space )
        {
            element_type mu( space );
            mu = space->min()+factor*( space->max()-space->min() );
            return mu;
        }


    //@}



protected:

private:

    void updateDefaultParameterNames()
        {
            M_parameterNames.resize( this->dimension() );
            for ( uint16_type d=0;d<this->dimension();++d )
                M_parameterNames[d] = (boost::format("mu%1%")%d).str();
        }

    friend class boost::serialization::access;
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const
        {
            ar & M_nDim;
            ar & M_min;
            ar & M_max;
            ar & M_parameterNames;
        }

    template<class Archive>
    void load( Archive & ar, const unsigned int version )
        {
            ar & M_nDim;
            ar & M_min;
            ar & M_max;
            ar & M_parameterNames;
        }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    private:

    //! parameter space dimension
    uint16_type M_nDim;

    //! min element
    element_type M_min;

    //! max element
    element_type M_max;

    //! parameter names
    std::vector<std::string> M_parameterNames;
};

template<uint16_type P> const uint16_type ParameterSpace<P>::Dimension;

//!
//! dynamic parameter space type definition
//!
using ParameterSpaceX = ParameterSpace<>;

template<uint16_type P>
//typename ParameterSpace<P>::sampling_ptrtype
std::shared_ptr<typename ParameterSpace<P>::Sampling>
ParameterSpace<P>::Sampling::complement() const
{
    //std::cout << "[ParameterSpace::Sampling::complement] start\n";
    sampling_ptrtype complement;
    bool is_in_sampling;
    int index=0;

    if ( M_supersampling )
    {
        complement = sampling_ptrtype( new sampling_type( M_space, 1, M_supersampling ) );
        complement->clear();
        for( auto const& mu_supersampling : *M_supersampling )
        {
            is_in_sampling=false;
            for( auto const& mu_sampling : *this )
            {
                if( mu_supersampling == mu_sampling )
                {
                    is_in_sampling=true;
                }
            }
            if( ! is_in_sampling )
                complement->push_back( mu_supersampling , index );

            index++;
        }
        return complement;
    }
    return complement;
#if 0
    //this method works if all procs have the same super sampling
    if ( M_supersampling )
    {
        //  std::cout << "[ParameterSpace::Sampling::complement] super sampling available\n";
        //  std::cout << "[ParameterSpace::Sampling::complement] this->size = " << this->size() << std::endl;
        //  std::cout << "[ParameterSpace::Sampling::complement] supersampling->size = " << M_supersampling->size() << std::endl;
        //  std::cout << "[ParameterSpace::Sampling::complement] complement->size = " << M_supersampling->size()-this->size() << std::endl;
        complement = sampling_ptrtype( new sampling_type( M_space, 1, M_supersampling ) );
        complement->clear();

        if ( !M_superindices.empty() )
        {
            for ( size_type i = 0, j = 0; i < M_supersampling->size(); ++i )
            {
                if ( std::find( M_superindices.begin(), M_superindices.end(), i ) == M_superindices.end() )
                {
                    //std::cout << "inserting " << j << "/" << i << " in complement of size (" << complement->size() << " out of " << M_supersampling->size() << std::endl;
                    if ( j < M_supersampling->size()-this->size() )
                    {
                        //complement->at( j ) = M_supersampling->at( i );
                        complement->push_back( M_supersampling->at( i ),i );
                        ++j;
                    }
                }
            }
        }

        return complement;
    }
#endif

}
template<uint16_type P>
std::shared_ptr<typename ParameterSpace<P>::Sampling>
ParameterSpace<P>::Sampling::searchNearestNeighbors( element_type const& mu,
                                                     size_type _M ,
                                                     std::vector<int>& index_vector,
                                                     bool broadcast )
{
    // make sure that the sampling set is big enough
    size_type M = std::min( _M, size_type(this->size()) );
    sampling_ptrtype neighbors( new sampling_type( M_space, M ) );
#if defined(FEELPP_HAS_ANN_H)

    if( M_space->worldComm().isMasterRank() || !broadcast )
    {
        //std::cout << "[ParameterSpace::Sampling::searchNearestNeighbors] start\n";
        //if ( !M_kdtree )
        //{
        //std::cout << "[ParameterSpace::Sampling::searchNearestNeighbors] building data points\n";
        ANNpointArray data_pts;
        data_pts = annAllocPts( this->size(), M_space->dimension() );

        for ( size_type i = 0; i < this->size(); ++i )
        {
            std::copy( this->at( i ).data(), this->at( i ).data()+M_space->dimension(), data_pts[i] );
            FEELPP_ASSERT( data_pts[i] != 0 )( i ) .error( "invalid pointer" );
        }

        //std::cout << "[ParameterSpace::Sampling::searchNearestNeighbors] building tree in R^" <<  M_space->Dimension << "\n";
        M_kdtree = kdtree_ptrtype( new kdtree_type( data_pts, this->size(), M_space->dimension() ) );
        //}


        std::vector<int> nnIdx( M );
        std::vector<double> dists( M );
        double eps = 0;
        M_kdtree->annkSearch( // search
                             const_cast<double*>( mu.data() ), // query point
                             M,       // number of near neighbors
                             nnIdx.data(),   // nearest neighbors (returned)
                             dists.data(),   // distance (returned)
                             eps );   // error bound



        if ( M_supersampling )
        {
            neighbors->setSuperSampling( M_supersampling );

            if ( !M_superindices.empty() )
                neighbors->M_superindices.resize( M );
        }

        //std::cout << "[parameterspace::sampling::searchNearestNeighbors] neighbor size = " << neighbors->size() <<std::endl;
        for ( size_type i = 0; i < M; ++i )
        {
            //std::cout << "[parameterspace::sampling::searchNearestNeighbors] neighbor: " <<i << " distance = " << dists[i] << "\n";
            neighbors->at( i ) = this->at( nnIdx[i] );
            index_vector.push_back( nnIdx[i] );

            //std::cout<<"[parameterspace::sampling::searchNearestNeighbors] M_superindices.size() : "<<M_superindices.size()<<std::endl;
            if ( M_supersampling && !M_superindices.empty() )
                neighbors->M_superindices[i] = M_superindices[ nnIdx[i] ];
            //std::cout << "[parameterspace::sampling::searchNearestNeighbors] " << neighbors->at( i ) << "\n";
        }
        annDeallocPts( data_pts );
    }//end of procMaster

    if ( broadcast )
    {
        boost::mpi::broadcast( M_space->worldComm() , neighbors , M_space->worldComm().masterRank() );
        boost::mpi::broadcast( M_space->worldComm() , index_vector , M_space->worldComm().masterRank() );
    }


#endif /* FEELPP_HAS_ANN_H */

    return neighbors;
}
} // Feel
#endif /* __ParameterSpace_H */
