/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELMODELS_TIMERTOOL_HPP
#define FEELMODELS_TIMERTOOL_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <unordered_map>
//#include <boost/timer/timer.hpp>

namespace Feel
{
class TimerToolBase
{
public:
    virtual ~TimerToolBase() = default;
    virtual bool isActive() = 0;

    virtual void start() = 0;
    virtual double stop( std::string const& key = "" ) = 0;
    virtual double elapsed( std::string const& key = "" ) = 0;
    virtual void restart() = 0;
    virtual void save() = 0;
    virtual void setDataValue(std::string const& key,double val) = 0;
    virtual void addDataValue(std::string const& key,double val) = 0;
    virtual double accumulateTime() = 0;
    virtual void setAdditionalParameter(std::string const& keyParam,boost::any const& d ) = 0;
};

class TimerTool : public TimerToolBase
{
public:
    typedef boost::mpi::timer timer_type;
    //typedef boost::timer::cpu_timer timer_type;

    TimerTool( std::string fileName = "timers.data", worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    TimerTool( TimerTool const& e ) = default;
    ~TimerTool() = default;
    //TimerTool( TimerTool && e ) = default;
    bool isActive() { return true; }

    void setReinitSaveFile( bool b );
    void setSaveFileMasterRank( bool b );
    void setSaveFileMax( bool b );
    void setSaveFileMin( bool b );
    void setSaveFileMean( bool b );

    void start();
    double stop( std::string const& key = "" );
    double elapsed( std::string const& key = "" );
    double elapsed( std::string const& key, int id );
    void restart();
    void reset();
    void save();

    double accumulateTime();

    void setDataValue(std::string const& key,double val);
    void addDataValue(std::string const& key,double val);
    double dataRegister(std::string const& key) const;
    void setAdditionalParameter(std::string const& keyParam,boost::any const& d );

    int maxIdActiveTimer() const;
    void resize( int newsize );

private :
    void saveImpl( std::string const& filename );

private :

    std::shared_ptr<WorldComm> M_worldComm;
    //boost::timer M_timer;
    std::vector<bool> M_activeTimers;
    std::vector<timer_type> M_timer;
    std::vector<double> M_tElapsedCurrent, M_tElapsedAccumulate;
    std::unordered_map<std::string, double> M_dataRegister;
    std::list<std::string> M_orderingData;
    std::unordered_map<std::string,std::tuple<boost::any> > M_additionalParameters;// time,...
    std::string M_fileName;
    bool M_reinitSaveFile;
    bool M_saveFileMasterRank, M_saveFileMax, M_saveFileMin, M_saveFileMean;

};

class TimerToolNull : public TimerToolBase
{
public:
    virtual ~TimerToolNull() = default;
    bool isActive() { return false; }

    void start() {};
    double stop( std::string const& key = "" ) { return 0; }
    double elapsed( std::string const& key = "" ) { return 0; }
    void restart() {}
    void save() {}
    void setDataValue(std::string const& key,double val) {};
    void addDataValue(std::string const& key,double val) {};
    double accumulateTime() { return 0; }
    void setAdditionalParameter(std::string const& keyParam,boost::any const& d ) {};
};

} // Feel
#endif /* FEELMODELS_TIMERTOOL_HPP */


