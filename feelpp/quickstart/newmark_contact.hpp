#pragma once

#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/utility.hpp>

#include <boost/foreach.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <feel/feelcore/parameter.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelts/tsbase.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel
{
namespace ublas = boost::numeric::ublas;

template<typename SpaceType>
class NewmarkContact : public TSBase
{
    friend class boost::serialization::access;
    typedef TSBase super;
public:
    typedef NewmarkContact<SpaceType> NewmarkContact_type;
    typedef std::shared_ptr<NewmarkContact_type> NewmarkContact_ptrtype;
    typedef SpaceType space_type;
    typedef std::shared_ptr<space_type>  space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::element_ptrtype element_ptrtype;
    typedef typename space_type::return_type return_type;
    typedef typename element_type::value_type value_type;
    typedef element_ptrtype unknown_type;
    typedef typename node<value_type>::type node_type;

    typedef typename super::time_iterator time_iterator;

    /**
     * Constructor
     *
     * @param vm option storage
     * @param space approximation space
     * @param name name of the NewmarkContact ts
     * @param prefix prefix of the NewmarkContact ts
     */
    FEELPP_DEPRECATED
    NewmarkContact( po::variables_map const& vm, space_ptrtype const& space, std::string const& name, std::string const& prefix="" )
        : NewmarkContact_type( space,name,prefix ) {}
    NewmarkContact( space_ptrtype const& space, std::string const& name, std::string const& prefix="" );

    /**
     * Constructor
     * @param space approximation space
     * @param name name of the NewmarkContact ts
     */
    NewmarkContact( space_ptrtype const& space, std::string const& name );

    //! copy operator
    NewmarkContact( NewmarkContact const& b )
        :
        super( b ),
        M_space( b.M_space ),
        M_previousUnknown( b.M_previousUnknown ),
        M_previousVel( b.M_previousVel ),
        M_previousAcc( b.M_previousAcc ),
        M_currentVel( b.M_currentVel ),
        M_currentAcc( b.M_currentAcc ),
        M_polyFirstDeriv( b.M_polyFirstDeriv ),
        M_polySecondDeriv( b.M_polySecondDeriv ),
        M_gamma( b.M_gamma ),
        M_beta( b.M_beta )
    {}

    ~NewmarkContact() override;

    /**
       Initialize all the entries of the unknown vector to be derived with the
       vector u0 (duplicated)
    */
    void initialize( element_type const& u0 );

    double start();
    double start( element_type const& u0 );

    /**
       restart the NewmarkContact scheme
    */
    double restart();

    /**
       Update the vectors of the previous time steps by shifting on the right
       the old values.
       @param u_curr current (new) value of the state vector
    */
    template<typename container_type>
    void shiftRight( typename space_type::template Element<value_type, container_type> const& u_curr );

    double next() const override;

    /**
     * shift previousstored solutions with respect to the new one \c u_curr
     * if \c updateVelAcc is true then the velocity and acceleration are updated too
     */
    template<typename container_type>
    double next( typename space_type::template Element<value_type, container_type> const& u_curr, bool updateVelAcc = true );

    /**
     * from the displacement \c u_curr update the velocity and acceleration
     */
    template<typename container_type>
    void updateFromDisp( typename space_type::template Element<value_type, container_type> const& u_curr, int previousTimeStep = 0 );


    //! Coefficients \f$ \gamma \f$ and \f$ \beta \f$  of the time NewmarkContact discretization
    double gamma() const { return M_gamma; }
    double beta() const { return M_beta; }

    void set_gamma(double gamma) { M_gamma = gamma; }
    void set_beta(double beta) { M_beta = beta; }


    double polyDerivCoefficient() const { return this->polySecondDerivCoefficient(); }
    double polyFirstDerivCoefficient() const
    {
        return this->gamma()/(this->beta()*this->timeStep());
    }
    double polySecondDerivCoefficient() const
    {
        return 1.0/(this->beta()*std::pow(this->timeStep(),2));
    }

    //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative formula
    element_type const& polyDeriv() const;
    element_type const& polyFirstDeriv() const;
    element_type const& polySecondDeriv() const;

    //! Return the last unknown
    element_type const& previousUnknown( int k = 0 ) const;
    element_type const& previousVelocity( int k = 0 ) const;
    element_type const& previousAcceleration( int k = 0 ) const;

    element_type & currentVelocity();
    element_type const& currentVelocity() const;
    element_ptrtype & currentVelocityPtr();
    element_ptrtype const& currentVelocityPtr() const;
    element_type & currentAcceleration();
    element_type const& currentAcceleration() const;
    element_ptrtype & currentAccelerationPtr();
    element_ptrtype const& currentAccelerationPtr() const;

    void loadCurrent();

private:
    void init();
    void initPreviousFields();

    void saveCurrent();

    //! save/load ts metadata
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        DVLOG(2) << "[NewmarkContact::serialize] saving/loading archive\n";
        ar & boost::serialization::base_object<TSBase>( *this );
    }

    //! compute first derivative poly of rhs
    void computePolyFirstDeriv();
    //! compute second derivative poly of rhs
    void computePolySecondDeriv();

private:

    //! space
    space_ptrtype M_space;

    //! Last n state vectors
    std::vector<element_ptrtype> M_previousUnknown, M_previousVel, M_previousAcc;
    unknown_type M_currentVel, M_currentAcc;
    unknown_type M_polyFirstDeriv,M_polySecondDeriv;

    //! Coefficients \f$ \gamma \f$ and \f$ \beta \f$  of the time NewmarkContact discretization
    double M_gamma, M_beta;

};

template <typename SpaceType>
NewmarkContact<SpaceType>::NewmarkContact( space_ptrtype const& __space,
                             std::string const& name,
                             std::string const& prefix )
    :
    super( name, prefix, __space->worldComm() ),
    M_space( __space ),
    M_currentVel( unknown_type( new element_type( M_space ) ) ),
    M_currentAcc( unknown_type( new element_type( M_space ) ) ),
    M_polyFirstDeriv( unknown_type( new element_type( M_space ) ) ),
    M_polySecondDeriv( unknown_type( new element_type( M_space ) ) ),
    M_gamma( 0.5 ),
    M_beta( 0.25 )
{
    this->initPreviousFields();
}

template <typename SpaceType>
NewmarkContact<SpaceType>::NewmarkContact( space_ptrtype const& __space,
                             std::string const& name  )
    :
    super( name, __space->worldComm() ),
    M_space( __space ),
    M_currentVel( unknown_type( new element_type( M_space ) ) ),
    M_currentAcc( unknown_type( new element_type( M_space ) ) ),
    M_polyFirstDeriv( unknown_type( new element_type( M_space ) ) ),
    M_polySecondDeriv( unknown_type( new element_type( M_space ) ) ),
    M_gamma( 0.5 ),
    M_beta( 0.25 )
{
    this->initPreviousFields();
}

template <typename SpaceType>
void
NewmarkContact<SpaceType>::initPreviousFields()
{
    int sizePreviousDisp = 4;
    M_previousUnknown.resize( sizePreviousDisp );
    M_previousVel.resize( sizePreviousDisp );
    M_previousAcc.resize( sizePreviousDisp );
    for ( uint8_type k = 0 ; k < sizePreviousDisp ; ++k )
    {
        M_previousUnknown[k] = element_ptrtype( new element_type( M_space ) );
        M_previousVel[k] = element_ptrtype( new element_type( M_space ) );
        M_previousAcc[k] = element_ptrtype( new element_type( M_space ) );
    }
}

template <typename SpaceType>
NewmarkContact<SpaceType>::~NewmarkContact()
{}

template <typename SpaceType>
void
NewmarkContact<SpaceType>::init()
{
    if ( this->path().empty() )
        this->setPathSave( (boost::format("NewmarkContact_dt_%1%")%this->timeStep()).str() );

    super::init();

    if ( this->isRestart() )
    {
        fs::path dirPath = ( this->restartPath().empty() )? this->path() : this->restartPath()/this->path();

        for ( int p = 0; p < std::min( (int)(M_previousUnknown.size()), (int)(M_iteration+1) ); ++p )
        {
            if ( this->fileFormat() == "hdf5")
            {
#ifdef FEELPP_HAS_HDF5
                M_previousUnknown[p]->loadHDF5( ( dirPath / (boost::format("%1%-unknown-%2%.h5")%M_name %(M_iteration-p)).str() ).string() );
                M_previousVel[p]->loadHDF5( ( dirPath / (boost::format("%1%-velocity-%2%.h5")%M_name %(M_iteration-p)).str() ).string() );
                M_previousAcc[p]->loadHDF5( ( dirPath / (boost::format("%1%-acceleration-%2%.h5")%M_name %(M_iteration-p)).str() ).string() );
#else
                CHECK( false ) << "hdf5 not detected";
#endif
            }
            else if ( this->fileFormat() == "binary")
            {
                // create and open a character archive for output
                std::string procsufix = (boost::format("-proc%1%_%2%") %this->worldComm().globalRank() %this->worldComm().globalSize() ).str();

                std::ostringstream ostrUnknown;
                ostrUnknown << M_name << "-unknown-" << M_iteration-p;
                if( M_rankProcInNameOfFiles )
                    ostrUnknown << procsufix;
                DVLOG(2) << "[NewmarkContact::init()] load file: " << ostrUnknown.str() << "\n";
                std::ifstream ifsUnknown;
                ifsUnknown.open( dirPath/ostrUnknown.str() );
                // load data from archive
                boost::archive::binary_iarchive iaUnknown( ifsUnknown );
                iaUnknown >> *(M_previousUnknown[p]);

                std::ostringstream ostrVel;
                ostrVel << M_name << "-velocity-" << M_iteration-p;
                if( M_rankProcInNameOfFiles )
                    ostrVel << procsufix;
                DVLOG(2) << "[NewmarkContact::init()] load file: " << ostrVel.str() << "\n";
                std::ifstream ifsVel;
                ifsVel.open( dirPath/ostrVel.str() );
                // load data from archive
                boost::archive::binary_iarchive iaVel( ifsVel );
                iaVel >> *(M_previousVel[p]);

                std::ostringstream ostrAcc;
                ostrAcc << M_name << "-acceleration-" << M_iteration-p;
                if( M_rankProcInNameOfFiles )
                    ostrAcc << procsufix;
                DVLOG(2) << "[NewmarkContact::init()] load file: " << ostrAcc.str() << "\n";
                std::ifstream ifsAcc;
                ifsAcc.open( dirPath/ostrAcc.str() );
                // load data from archive
                boost::archive::binary_iarchive iaAcc( ifsAcc );
                iaAcc >> *(M_previousAcc[p]);

            } // binary
        } // p
        // current vel,acc are same that previous at this time
        *M_currentVel = *(M_previousVel[0]);
        *M_currentAcc = *(M_previousAcc[0]);

    } // isRestart

    // compute first derivative poly of rhs
    this->computePolyFirstDeriv();
    // compute second derivative poly of rhs
    this->computePolySecondDeriv();
}

template <typename SpaceType>
void
NewmarkContact<SpaceType>::initialize( element_type const& u0 )
{
    M_time_values_map.clear();
    std::ostringstream ostr;

    if( M_rankProcInNameOfFiles )
        ostr << M_name << "-" << 0<<"-proc"<<this->worldComm().globalRank()<<"on"<<this->worldComm().globalSize();
    else
        ostr << M_name << "-" << 0;

    M_time_values_map.push_back( M_Ti );

    std::for_each( M_previousUnknown.begin(), M_previousUnknown.end(), *boost::lambda::_1 = u0 );

    // compute first derivative poly of rhs
    this->computePolyFirstDeriv();
    // compute second derivative poly of rhs
    this->computePolySecondDeriv();

    this->saveCurrent();
}

template <typename SpaceType>
double
NewmarkContact<SpaceType>::start()
{
    if ( this->isRestart() )
        return this->restart();

    this->init();
    return super::start();
}

template <typename SpaceType>
double
NewmarkContact<SpaceType>::start( element_type const& u0 )
{
    if ( this->isRestart() )
        return this->restart();

    this->init();
    this->initialize( u0 );
    auto res = super::start();
    return res;
}

template <typename SpaceType>
double
NewmarkContact<SpaceType>::restart()
{
    this->init();
    return super::restart();
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_type const&
NewmarkContact<SpaceType>::previousUnknown( int k ) const
{
    return *(M_previousUnknown[k]);
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_type const&
NewmarkContact<SpaceType>::previousVelocity( int k ) const
{
    return *(M_previousVel[k]);
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_type const&
NewmarkContact<SpaceType>::previousAcceleration( int k ) const
{
    return *(M_previousAcc[k]);
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_type &
NewmarkContact<SpaceType>::currentVelocity()
{
  return *M_currentVel;
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_type const&
NewmarkContact<SpaceType>::currentVelocity() const
{
  return *M_currentVel;
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_ptrtype &
NewmarkContact<SpaceType>::currentVelocityPtr()
{
  return M_currentVel;
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_ptrtype const&
NewmarkContact<SpaceType>::currentVelocityPtr() const
{
  return M_currentVel;
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_type &
NewmarkContact<SpaceType>::currentAcceleration()
{
  return *M_currentAcc;
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_type const&
NewmarkContact<SpaceType>::currentAcceleration() const
{
  return *M_currentAcc;
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_ptrtype &
NewmarkContact<SpaceType>::currentAccelerationPtr()
{
  return M_currentAcc;
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_ptrtype const&
NewmarkContact<SpaceType>::currentAccelerationPtr() const
{
  return M_currentAcc;
}



template <typename SpaceType>
void
NewmarkContact<SpaceType>::saveCurrent()
{
    if (!this->saveInFile()) return;

    bool doSave=false;
    int sizePreviousDisp = M_previousUnknown.size();
    for ( uint8_type i = 0; i < sizePreviousDisp && !doSave; ++i )
    {
        int iterTranslate = this->iteration() + sizePreviousDisp-(i+1);
        if ( (iterTranslate % this->saveFreq()) == 0 )
            doSave=true;
    }

    if (!doSave) return;

    TSBaseMetadata tssaver( *this );
    tssaver.save();

    if ( this->fileFormat() == "hdf5")
    {
#ifdef FEELPP_HAS_HDF5
        M_previousUnknown[0]->saveHDF5( (M_path_save / (boost::format("%1%-unknown-%2%.h5")%M_name %M_iteration).str() ).string() );
        M_previousVel[0]->saveHDF5( (M_path_save / (boost::format("%1%-velocity-%2%.h5")%M_name %M_iteration).str() ).string() );
        M_previousAcc[0]->saveHDF5( (M_path_save / (boost::format("%1%-acceleration-%2%.h5")%M_name %M_iteration).str() ).string() );
#else
        CHECK( false ) << "hdf5 not detected";
#endif
    }
    else if ( this->fileFormat() == "binary")
    {
        std::string procsufix = (boost::format("-proc%1%_%2%") %this->worldComm().globalRank() %this->worldComm().globalSize() ).str();

        std::ostringstream ostrUnknown;
        ostrUnknown << M_name << "-unknown-" << M_iteration;
        if( M_rankProcInNameOfFiles )
            ostrUnknown << procsufix;
        std::ofstream ofsUnknown( M_path_save / ostrUnknown.str() );
        // save data in archive
        boost::archive::binary_oarchive oaUnknown( ofsUnknown );
        oaUnknown << *(M_previousUnknown[0]);

        std::ostringstream ostrVel;
        ostrVel << M_name << "-velocity-" << M_iteration;
        if( M_rankProcInNameOfFiles )
            ostrVel << procsufix;
        std::ofstream ofsVel( M_path_save / ostrVel.str() );
        // save data in archive
        boost::archive::binary_oarchive oaVel( ofsVel );
        oaVel << *(M_previousVel[0]);

        std::ostringstream ostrAcc;
        ostrAcc << M_name << "-acceleration-" << M_iteration;
        if( M_rankProcInNameOfFiles )
            ostrAcc << procsufix;
        std::ofstream ofsAcc( M_path_save / ostrAcc.str() );
        // save data in archive
        boost::archive::binary_oarchive oaAcc( ofsAcc );
        oaAcc << *(M_previousAcc[0]);
    }
}

template <typename SpaceType>
void
NewmarkContact<SpaceType>::loadCurrent()
{
    if ( this->fileFormat() == "hdf5")
    {
#ifdef FEELPP_HAS_HDF5
        M_previousUnknown[0]->loadHDF5( (M_path_save / (boost::format("%1%-unknown-%2%.h5")%M_name %M_iteration).str() ).string() );
        M_previousVel[0]->loadHDF5( (M_path_save / (boost::format("%1%-velocity-%2%.h5")%M_name %M_iteration).str() ).string() );
        M_previousAcc[0]->loadHDF5( (M_path_save / (boost::format("%1%-acceleration-%2%.h5")%M_name %M_iteration).str() ).string() );
#else
        CHECK( false ) << "hdf5 not detected";
#endif
    }
    else if ( this->fileFormat() == "binary")
    {
        std::string procsufix = (boost::format("-proc%1%_%2%") %this->worldComm().globalRank() %this->worldComm().globalSize() ).str();

        std::ostringstream ostrUnknown;
        ostrUnknown << M_name << "-unknown-" << M_iteration;
        if( M_rankProcInNameOfFiles )
            ostrUnknown << procsufix;
        std::ifstream ifsUnknown( M_path_save / ostrUnknown.str() );
        // load data from archive
        boost::archive::binary_iarchive iaUnknown( ifsUnknown );
        iaUnknown >> *(M_previousUnknown[0]);

        std::ostringstream ostrVel;
        ostrVel << M_name << "-velocity-" << M_iteration;
        if( M_rankProcInNameOfFiles )
            ostrVel << procsufix;
        std::ifstream ifsVel( M_path_save / ostrVel.str() );
        // load data from archive
        boost::archive::binary_iarchive iaVel( ifsVel );
        iaVel >> *(M_previousVel[0]);

        std::ostringstream ostrAcc;
        ostrAcc << M_name << "-acceleration-" << M_iteration;
        if( M_rankProcInNameOfFiles )
            ostrAcc << procsufix;
        std::ifstream ifsAcc( M_path_save / ostrAcc.str() );
        // load data from archive
        boost::archive::binary_iarchive iaAcc( ifsAcc );
        iaAcc >> *(M_previousAcc[0]);
    }
}

template <typename SpaceType>
template<typename container_type>
void
NewmarkContact<SpaceType>::shiftRight( typename space_type::template Element<value_type, container_type> const& __new_unk )
{
    DVLOG(2) << "shiftRight: inserting time " << this->time() << "s\n";
    super::shiftRight();

    // shift all previously stored bdf data
    using namespace boost::lambda;
    typename std::vector<element_ptrtype>::reverse_iterator __itDisp = boost::next( M_previousUnknown.rbegin() );
    std::for_each( M_previousUnknown.rbegin(), boost::prior( M_previousUnknown.rend() ),
                   ( *lambda::_1 = *( *lambda::var( __itDisp ) ), ++lambda::var( __itDisp ) ) );
    typename std::vector<element_ptrtype>::reverse_iterator __itVel = boost::next( M_previousVel.rbegin() );
    std::for_each( M_previousVel.rbegin(), boost::prior( M_previousVel.rend() ),
                   ( *lambda::_1 = *( *lambda::var( __itVel ) ), ++lambda::var( __itVel ) ) );
    typename std::vector<element_ptrtype>::reverse_iterator __itAcc = boost::next( M_previousAcc.rbegin() );
    std::for_each( M_previousAcc.rbegin(), boost::prior( M_previousAcc.rend() ),
                   ( *lambda::_1 = *( *lambda::var( __itAcc ) ), ++lambda::var( __itAcc ) ) );

    // shift all previously stored  data
    *M_previousUnknown[0] = __new_unk;
    // update M_currentVelocity and M_currentAcceleration with new disp and define as previous vel/acc
    this->updateFromDisp(__new_unk,1);
    *(M_previousVel[0]) = *M_currentVel;
    *(M_previousAcc[0]) = *M_currentAcc;

    // save newly stored data
    this->saveCurrent();

    // compute first derivative poly of rhs
    this->computePolyFirstDeriv();
    // compute second derivative poly of rhs
    this->computePolySecondDeriv();
}


template <typename SpaceType>
double
NewmarkContact<SpaceType>::next() const
{
  return super::next();
}

template <typename SpaceType>
template<typename container_type>
double
NewmarkContact<SpaceType>::next( typename space_type::template Element<value_type, container_type> const& u_curr, bool updateVelAcc )
{
    if ( updateVelAcc )
        this->updateFromDisp( u_curr );
    this->shiftRight( u_curr );
    return super::next();
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_type const&
NewmarkContact<SpaceType>::polyDeriv() const
{
    return *M_polySecondDeriv;
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_type const&
NewmarkContact<SpaceType>::polyFirstDeriv() const
{
    return *M_polyFirstDeriv;
}

template <typename SpaceType>
typename NewmarkContact<SpaceType>::element_type const&
NewmarkContact<SpaceType>::polySecondDeriv() const
{
    return *M_polySecondDeriv;
}

template <typename SpaceType>
void
NewmarkContact<SpaceType>::computePolyFirstDeriv()
{
    double deltaT = this->timeStep();
    double cst_vel_displ = M_gamma/(M_beta*deltaT);
    double cst_vel_vel = (M_gamma/M_beta)-1;
    double cst_vel_acc = deltaT*((M_gamma/(2*M_beta))-1 );

    M_polyFirstDeriv->zero();
    M_polyFirstDeriv->add(cst_vel_displ, this->previousUnknown());
    M_polyFirstDeriv->add(cst_vel_vel, this->previousVelocity());
    M_polyFirstDeriv->add(cst_vel_acc, this->previousAcceleration());
}

template <typename SpaceType>
void
NewmarkContact<SpaceType>::computePolySecondDeriv()
{
    double coef_disp = 1.0/(M_beta*std::pow(this->timeStep(),2));
    double coef_vel  = 1.0/(M_beta*this->timeStep());
    double coef_acc  = 1.0/(2*M_beta) - 1.0;
    M_polySecondDeriv->zero();
    M_polySecondDeriv->add(coef_disp, this->previousUnknown());
    M_polySecondDeriv->add(coef_vel, this->previousVelocity());
    M_polySecondDeriv->add(coef_acc, this->previousAcceleration());
}

template <typename SpaceType>
template<typename container_type>
void
NewmarkContact<SpaceType>::updateFromDisp( typename space_type::template Element<value_type, container_type> const& __new_unk, int previousTimeStep )
{
  double deltaT = this->timeStep();

  double cst_vel_displ = M_gamma/(M_beta*deltaT);
  double cst_vel_vel = (M_gamma/M_beta)-1;
  double cst_vel_acc = deltaT*((M_gamma/(2*M_beta))-1 );
  double cst_acc_displ = 1./(M_beta*std::pow(deltaT,2));
  double cst_acc_vel = 1./(M_beta*deltaT);
  double cst_acc_acc = ((1./(2*M_beta))-1);

  M_currentVel->zero();
  M_currentVel->add(cst_vel_displ, __new_unk );
  M_currentVel->add(-cst_vel_displ, this->previousUnknown(previousTimeStep));
  M_currentVel->add(-cst_vel_vel, this->previousVelocity(previousTimeStep));
  M_currentVel->add(-cst_vel_acc, this->previousAcceleration(previousTimeStep));

  M_currentAcc->zero();
  M_currentAcc->add(cst_acc_displ, __new_unk );
  M_currentAcc->add(-cst_acc_displ, this->previousUnknown(previousTimeStep));
  M_currentAcc->add(-cst_acc_vel, this->previousVelocity(previousTimeStep));
  M_currentAcc->add(-cst_acc_acc, this->previousAcceleration(previousTimeStep));
}


template<typename SpaceT>
auto newmarkContact(  SpaceT const& space, bool steady , double initial_time , double final_time, double time_step, double gamma, double beta)
{

    std::string const& prefix = "";
    std::string const& name = "NewmarkContact";
    //double initial_time = args.get_else_invocable( _initial_time, [&prefix,&vm](){ return doption(_prefix=prefix,_name="ts.time-initial",_vm=vm); } );
    //double final_time = args.get_else_invocable( _final_time, [&prefix,&vm](){ return doption(_prefix=prefix,_name="ts.time-final",_vm=vm); } );
    //double time_step = args.get_else_invocable( _time_step, [&prefix,&vm](){ return doption(_prefix=prefix,_name="ts.time-step",_vm=vm); } );
    //bool steady = args.get_else_invocable( _steady, [&prefix,&vm](){ return boption(_prefix=prefix,_name="ts.steady",_vm=vm); } );
    
    bool restart = false;
    std::string const& restart_path = "";
    bool restart_at_last_save = false;
    bool save = true;
    int freq = 1;
    std::string const& format = "hdf5";
    bool rank_proc_in_files_name = false;

    using _space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(space)>>>;
    auto theNewmarkContact = std::make_shared<NewmarkContact<_space_type>>( space,name,prefix );
    
    theNewmarkContact->set_gamma(gamma);
    theNewmarkContact->set_beta(beta);

    std::cout << "Gamma : " << theNewmarkContact->gamma() << std::endl;
    std::cout << "Beta : " << theNewmarkContact->beta() << std::endl;

    theNewmarkContact->setTimeInitial( initial_time );
    theNewmarkContact->setTimeFinal( final_time );
    theNewmarkContact->setTimeStep( time_step );
    theNewmarkContact->setSteady( steady );
    theNewmarkContact->setRestart( restart );
    theNewmarkContact->setRestartPath( restart_path );
    theNewmarkContact->setRestartAtLastSave( restart_at_last_save );
    theNewmarkContact->setSaveInFile( save );
    theNewmarkContact->setSaveFreq( freq );
    theNewmarkContact->setfileFormat( format );
    theNewmarkContact->setRankProcInNameOfFiles( rank_proc_in_files_name );
    
    return theNewmarkContact;
}


}

