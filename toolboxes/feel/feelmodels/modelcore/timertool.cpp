/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/timertool.hpp>

#include <algorithm>

namespace Feel {


TimerTool::TimerTool( std::string fileName, worldcomm_ptr_t const& worldComm )
    :
    TimerToolBase(),
    M_worldComm( worldComm ),
    M_activeTimers(1,false),
    M_timer(1),
    M_tElapsedCurrent(1,0),
    M_tElapsedAccumulate(1,0),
    M_fileName( fileName ),
    M_reinitSaveFile( false ),
    M_saveFileMasterRank( true ),
    M_saveFileMax( false ),
    M_saveFileMin( false ),
    M_saveFileMean( false )
{}

void
TimerTool::setReinitSaveFile( bool b ) { M_reinitSaveFile=b; }

void
TimerTool::setSaveFileMasterRank( bool b ) { M_saveFileMasterRank=b; }
void
TimerTool::setSaveFileMax( bool b ) { M_saveFileMax=b; }
void
TimerTool::setSaveFileMin( bool b ) { M_saveFileMin=b; }
void
TimerTool::setSaveFileMean( bool b ) { M_saveFileMean=b; }


int
TimerTool::maxIdActiveTimer() const
{
    int res = -1;
    for ( int k=0 ; k<M_activeTimers.size() ; ++k)
    {
        if ( M_activeTimers[k] ) res=k;
    }
    return res;
}

void
TimerTool::resize( int newsize )
{
    M_activeTimers.resize(newsize,false);
    M_timer.resize(newsize);
    M_tElapsedCurrent.resize(newsize,0);
    M_tElapsedAccumulate.resize(newsize,0);
}
void
TimerTool::start()
{
    int maxid = this->maxIdActiveTimer();
    if ( maxid == -1 )
    {
        if ( M_activeTimers.size() == 0 )
            this->resize( 1 );
        M_activeTimers[0] = true;
        M_tElapsedCurrent[0] = 0;
        M_tElapsedAccumulate[0] = 0;
        M_timer[0].restart();
    }
    else
    {
        if ( maxid+1 >= M_activeTimers.size() )
            this->resize( maxid+2 );
        M_activeTimers[maxid+1] = true;
        M_tElapsedCurrent[maxid+1] = 0;
        M_tElapsedAccumulate[maxid+1] = 0;
        M_timer[maxid+1].restart();
    }
}

double
TimerTool::stop( std::string const& key )
{
    int maxid = this->maxIdActiveTimer();
    if( maxid < 0 )
        return 0;
    else
    {
        double res = this->elapsed( key, maxid );
        M_activeTimers[maxid] = false;
        return res;
    }
}

double
TimerTool::elapsed( std::string const& key )
{
    int maxid = this->maxIdActiveTimer();
    if( maxid < 0 )
        return 0;
    else
        return this->elapsed( key,maxid );
}
double
TimerTool::elapsed( std::string const& key, int id )
{
    M_tElapsedCurrent[id] = M_timer[id].elapsed();
    M_tElapsedAccumulate[id] += M_tElapsedCurrent[id];
    if ( !key.empty() )
        this->setDataValue(key,M_tElapsedCurrent[id]);
    //M_dataRegister[key] = M_tElapsedCurrent[id];
    return M_tElapsedCurrent[id];
}
void
TimerTool::restart()
{
    int maxid = this->maxIdActiveTimer();
    if( maxid < 0 ) return;
    M_timer[maxid].restart();
}
void
TimerTool::reset()
{
    int maxid = this->maxIdActiveTimer();
    if( maxid < 0 ) return;
    M_tElapsedAccumulate[maxid] = 0;
    M_timer[maxid].restart();
}

double
TimerTool::accumulateTime()
{
    int maxid = this->maxIdActiveTimer();
    if( maxid < 0 )
        return 0;
    else
        return M_tElapsedAccumulate[maxid];
}

void
TimerTool::setDataValue(std::string const& key,double val)
{
    if ( M_dataRegister.find(key) == M_dataRegister.end() )
        M_orderingData.push_back(key);
    M_dataRegister[key] = val;
}
void
TimerTool::addDataValue(std::string const& key,double val)
{
    if ( M_dataRegister.find(key) == M_dataRegister.end() )
    {
        M_orderingData.push_back(key);
        M_dataRegister[key] = val;
    }
    else
        M_dataRegister[key] += val;
}

double
TimerTool::dataRegister(std::string const& key) const
{
    CHECK(M_dataRegister.find( key ) != M_dataRegister.end() ) << "invalid key " << key;
    return M_dataRegister.find( key )->second;
}

void
TimerTool::setAdditionalParameter(std::string const& keyParam, boost::any const& d/*double d*/ )
{
    M_additionalParameters[keyParam] = std::make_tuple(d);
}

void
TimerTool::save()
{
    int numberOfSave = 0;
    if ( M_saveFileMasterRank ) numberOfSave++;
    if ( M_saveFileMax ) numberOfSave++;
    if ( M_saveFileMin ) numberOfSave++;
    if ( M_saveFileMean ) numberOfSave++;
    bool addSaveTypeParameter = numberOfSave > 1;

    if ( M_saveFileMasterRank )
    {
        if ( addSaveTypeParameter )
            this->setAdditionalParameter("saveType","MasterRank");
        this->saveImpl( M_fileName );
    }

    if ( M_saveFileMax || M_saveFileMin || M_saveFileMean )
    {
        std::vector<double> dataSend;
        for (auto const& dataKey : M_orderingData )
            dataSend.push_back( M_dataRegister.find(dataKey)->second );

        std::vector<std::vector<double> > dataRecv;
        mpi::gather( M_worldComm->globalComm(), dataSend, dataRecv, M_worldComm->masterRank() );

        if ( M_worldComm->isMasterRank() )
        {
            int nProc = dataRecv.size();
            if ( M_saveFileMax )
            {
                int k = 0;
                for (auto const& dataKey : M_orderingData )
                {
                    double dataMax = dataRecv[0][k];
                    for (int p=1;p<nProc;++p )
                        if ( dataRecv[p][k] > dataMax )
                            dataMax=dataRecv[p][k];
                    this->setDataValue( dataKey,dataMax );
                    ++k;
                }
                if ( addSaveTypeParameter )
                    this->setAdditionalParameter("saveType","Max");
                this->saveImpl( M_fileName );
                //maybe revert data value?
            }
            if ( M_saveFileMin )
            {
                int k = 0;
                for (auto const& dataKey : M_orderingData )
                {
                    double dataMin = dataRecv[0][k];
                    for (int p=1;p<nProc;++p )
                        if ( dataRecv[p][k] < dataMin )
                            dataMin=dataRecv[p][k];
                    this->setDataValue( dataKey,dataMin );
                    ++k;
                }
                if ( addSaveTypeParameter )
                    this->setAdditionalParameter("saveType","Min");
                this->saveImpl( M_fileName );
                //maybe revert data value?
            }
            if ( M_saveFileMean )
            {
                int k = 0;
                for (auto const& dataKey : M_orderingData )
                {
                    double dataMean = 0;
                    for ( int p=0;p<nProc;++p )
                        dataMean+=dataRecv[p][k];
                    dataMean /= nProc;
                    this->setDataValue( dataKey,dataMean );
                    ++k;
                }
                if ( addSaveTypeParameter )
                    this->setAdditionalParameter("saveType","Mean");
                this->saveImpl( M_fileName );
                //maybe revert data value?
            }
        }
    }
    // synchronisation
    // M_worldComm->barrier();
}

void
TimerTool::saveImpl( std::string const& filename )
{
    if ( M_worldComm->isMasterRank() )
    {
        // if file not exist, force to write the preambule
        if( !M_reinitSaveFile && !fs::exists(filename) )
            M_reinitSaveFile = true;
        std::ofstream file(filename.c_str(), (M_reinitSaveFile)? std::ios::out : (std::ios::out | std::ios::app) );
        if ( M_reinitSaveFile )
        {
            file << "# " << std::setw(10) << std::left << "nProc";
            for (auto const& addParam : M_additionalParameters)
                file << std::setw( std::max( int(addParam.first.size()+2), 20 ) ) << std::left << addParam.first;

            for (auto const& dataKey : M_orderingData )
                file << std::setw( std::max( int(dataKey.size()+2), 20 ) ) << std::right
                     << dataKey;
            file << "\n";
        }

        file << std::setw(12) << std::left << M_worldComm->globalSize();
        for (auto const& addParam : M_additionalParameters)
        {
            file << std::setw( std::max( int(addParam.first.size()+2), 20 ) ) << std::left << std::setprecision( 5 ) << std::fixed;
            boost::any const& param = std::get<0>(addParam.second);
            if ( boost::any_cast<double>( &param ) )
                file << boost::any_cast<double>( param );
            else if ( boost::any_cast<int>( &param ) )
                file << boost::any_cast<int>( param );
            else if ( boost::any_cast<std::string>( &param ) )
                file << boost::any_cast<std::string>( param );
            else if ( boost::any_cast<const char*>( &param ) )
                file << boost::any_cast<const char*>( param );
            else CHECK( false ) << "type unrecognized";
        }

        for (auto const& dataKey : M_orderingData )
            file << std::setw( std::max( int(dataKey.size()+2), 20 ) ) << std::right << std::setprecision(8) << std::scientific
                 << M_dataRegister.find(dataKey)->second;
        file << "\n";
        file.close();
    }
    M_reinitSaveFile = false;
}


} //namespace Feel
