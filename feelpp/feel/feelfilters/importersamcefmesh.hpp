/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_FILTERS_IMPORTERSAMCEFMESH_HPP
#define FEELPP_FILTERS_IMPORTERSAMCEFMESH_HPP

#include <feel/feelfilters/importer.hpp>

namespace Feel
{

template <typename MeshType>
class ImporterSamcefMesh : public Importer<MeshType>
{
    typedef Importer<MeshType> super;

public:
    typedef typename super::mesh_type mesh_type;
    typedef typename super::point_type point_type;
    typedef typename super::node_type node_type;
    typedef typename super::edge_type edge_type;
    typedef typename super::face_type face_type;
    typedef typename super::element_type element_type;
    using index_type = typename mesh_type::index_type;

    ImporterSamcefMesh( std::string const& filename, worldcomm_ptr_t const& _worldcomm = Environment::worldCommPtr() )
        :
        super( SAMCEF, ASCII, _worldcomm ),
        M_filename( filename )
        {}

    ImporterSamcefMesh( ImporterSamcefMesh const& ) = default;

    void visit( mesh_type* mesh ) override;
private:
    std::tuple<std::vector<std::string>,std::string>
    readNextPart( std::ifstream & __is, std::set<std::string> const& stopPartKeyword, bool & commandFinished );

    std::string readSEL( std::ifstream & __is );
private:
    std::string M_filename;

    std::map<index_type, std::tuple<std::string,std::vector<std::tuple<index_type,index_type,index_type>>>> M_selGroupMailles;
};

template <typename MeshType>
void
ImporterSamcefMesh<MeshType>::visit( mesh_type* mesh )
{
    std::ifstream __is( M_filename );
    if ( !__is.is_open() )
    {
        std::ostringstream ostr;
        ostr << "Invalid file name " << M_filename << " (file not found)\n";
        throw std::invalid_argument( ostr.str() );
    }

    rank_type partId = mesh->worldComm().localRank();

    std::string tmpString;

    std::vector<std::string> readedStr;
    std::string currentKeyword, nextKeyword;
    bool commandFinished = false;

    // go to NOE command
    while ( !__is.eof() )
    {
        commandFinished = false;
        std::tie(readedStr,nextKeyword) = this->readNextPart( __is, {}, commandFinished );
        if ( nextKeyword == ".NOE" )
            break;
    }
    if ( nextKeyword != ".NOE" )
        return;


    // load .NOE command
    commandFinished = false;
    __is >> currentKeyword;
    CHECK( currentKeyword == "I" ) << "something wrong " << currentKeyword;
    index_type tmpIndex;
    //double coords[mesh_type::nRealDim];
    node_type coords( mesh_type::nRealDim );
    while ( true )
    {
        CHECK( currentKeyword == "I" ) << "something wrong " << currentKeyword;
        std::tie(readedStr,nextKeyword) = this->readNextPart( __is, {"I"}, commandFinished );
        currentKeyword = nextKeyword;
        CHECK( !readedStr.empty() ) << "should not be empty";

        tmpIndex = std::stoi( readedStr.front() );
        int d=0;
        for ( int k=1;k<readedStr.size();k+=2)
        {
            coords[d++] = std::stod( readedStr[k+1] ); // we guess X x Y y Z z
        }

        point_type pt( tmpIndex, coords, false );
        pt.setProcessIdInPartition( partId );
        pt.setProcessId( partId );
        mesh->addPoint( std::move(pt) );

        if ( commandFinished )
            break;
    }

    if ( nextKeyword != ".MAI" )
        return;

    int markerId = 1;

    // load .MAI command
    index_type idElt = invalid_v<index_type>;
    commandFinished = false;
    __is >> currentKeyword;
    using element_wrapper_type = boost::reference_wrapper<element_type>;
    std::unordered_map<index_type,element_wrapper_type> mapSamcefEltIdToFeelElt;
    element_type e;
    e.setProcessIdInPartition( partId );
    e.setProcessId( partId );
    while ( true )
    {
        std::set<std::string> readedKeywords;

        while ( true ) // read info of current cell
        {
            std::tie(readedStr,nextKeyword) = this->readNextPart( __is, { "I", "ATT", "AT", "N" }, commandFinished );
            if ( currentKeyword == "I" )
            {
                idElt = std::stoi( readedStr.front() );
            }
            if ( currentKeyword == "ATT" || currentKeyword == "AT" )
            {
                markerId = std::stoi( readedStr.front() );
                e.setMarker2( markerId );
            }
            if ( currentKeyword == "N" )
            {
                std::vector<index_type> nodeIdsInElt;
                nodeIdsInElt.reserve( element_type::numPoints );
                for ( auto const& idStr : readedStr )
                {
                    tmpIndex = std::stoi( idStr );
                    if ( tmpIndex == 0 )
                        continue;
                    nodeIdsInElt.push_back( tmpIndex );
                }
                for ( size_type k = 0; k < element_type::numPoints; ++k )
                    e.setPoint( k, mesh->point( nodeIdsInElt[k] ) );
            }

            readedKeywords.insert( currentKeyword );
            currentKeyword = nextKeyword;

            if ( readedKeywords.find( nextKeyword ) != readedKeywords.end() || commandFinished )
            {
                break;
            }
        }

        // add cell in mesh
        auto const& [newEltIt,isAdded] = mesh->addElement( e, true );
        mapSamcefEltIdToFeelElt.emplace( idElt, boost::ref(newEltIt->second) );
        if ( commandFinished )
            break;
    }

    // read file up to the end
    while ( !__is.eof() )
    {
        //std::cout << "currentKeyword : " << currentKeyword << std::endl;
        if ( currentKeyword == ".SEL" )
        {
            // read .SEL command
            currentKeyword = readSEL( __is );
        }
        else
        {
            commandFinished = false;
            std::tie(readedStr,currentKeyword) = this->readNextPart( __is, {}, commandFinished );
        }
    }

    for ( auto const& [_mId,_gData] : M_selGroupMailles )
    {
#if 0
        // case hexa
        if ( _mId == 34 ) continue; // TO REMOVE
#endif
        if ( _mId > 3 ) continue; // TO REMOVE
        std::string const& markerName = std::get<0>( _gData );
        for ( auto const& _gRange : std::get<1>( _gData ) )
        {
            std::cout << "marker " << _mId << " range " << std::get<0>( _gRange ) << " , " << std::get<1>( _gRange ) << std::endl;
            for ( index_type eId = std::get<0>( _gRange ); eId <= std::get<1>( _gRange ); eId+=std::get<2>( _gRange ) )
            {
                auto itFindElt = mapSamcefEltIdToFeelElt.find( eId );
                CHECK( itFindElt !=mapSamcefEltIdToFeelElt.end() ) << "not find samcef elt id " << eId;
                auto & feelElt = unwrap_ref( itFindElt->second );
                feelElt.setMarker( _mId );
            }
        }
        mesh->addMarkerName( markerName, _mId, mesh_type::nDim );
    }
}



template <typename MeshType>
std::tuple<std::vector<std::string>,std::string>
ImporterSamcefMesh<MeshType>::readNextPart( std::ifstream & __is, std::set<std::string> const& stopPartKeyword, bool & commandFinished )
{
    std::string tmpString;
    std::vector<std::string> readedStr;
    std::string nextKeyword;
    bool isNewLine = __is.tellg() == 0;// (__is.tellg() == 0) || (__is.rdbuf()->sgetc() == '\n');
    while ( true )
    {
        if ( __is.eof() )
        {
            commandFinished = true;
            return std::make_tuple( std::move( readedStr ), std::move( nextKeyword ) );
        }


        char peekChar = __is.peek();
        if ( peekChar == '\n' || peekChar == '\r' )
        {
            while ( peekChar == '\n' || peekChar == '\r' )
            {
                if ( peekChar == '\n' )
                    isNewLine = true;
                __is.get();
                peekChar = __is.peek();
            }
        }

        __is >> tmpString;

        if ( isNewLine )
        {
            if ( !tmpString.empty() && tmpString.front() == '.' )
            {
                nextKeyword = tmpString;
                commandFinished = true;
                break;
            }
        }

        isNewLine = false;

        if ( tmpString == "$" )
            continue;

        if ( !tmpString.empty() && tmpString.front() == '!' )
        {
            //if ( __is.peek() != '\n' )
            std::getline(__is, tmpString);
            isNewLine = true;
            continue;
        }

        if ( stopPartKeyword.find( tmpString ) != stopPartKeyword.end() )
        {
            nextKeyword = tmpString;
            break;
        }

        readedStr.push_back( tmpString );
    }
    return std::make_tuple( std::move( readedStr ), std::move( nextKeyword ) );
}

template <typename MeshType>
std::string
ImporterSamcefMesh<MeshType>::readSEL( std::ifstream & __is)
{
    int markerId = 1;
    std::string currentKeyword, nextKeyword, tmpString;
    std::vector<std::string> readedStr;
    bool commandFinished = false;
    __is >> currentKeyword; // read GROUP
    while ( true ) // for each group
    {
        std::tie(readedStr,nextKeyword) = this->readNextPart( __is, {"MAILLES","FACES","NOEUDS","MAILLE","FACE","NOEUD","NOM"}, commandFinished ); // read num or names or nothing of the group
        markerId = std::stoi( readedStr[0] );
        //std::cout << "markerId="<<markerId<<std::endl;

        std::string typeOfGroup;
        std::string groupName;

        if ( nextKeyword == "NOM" )
        {
            std::tie(readedStr,nextKeyword) = this->readNextPart( __is, {"MAILLES","FACES","NOEUDS","MAILLE","FACE","NOEUD","NOM"}, commandFinished );
            groupName = readedStr.front();
        }

        if ( nextKeyword.substr(0,5) == "MAILL" )
        {
            typeOfGroup = "MAILLES";
            std::tie(readedStr,nextKeyword) = this->readNextPart( __is, {"I","NOM"}, commandFinished );
        }
        else if ( nextKeyword.substr(0,4) == "FACE" )
        {
            typeOfGroup = "FACES";
            __is >> tmpString >> tmpString;
        }
        else if ( nextKeyword.substr(0,5) == "NOEUD" )
        {
            typeOfGroup = "NOEUDS";
            std::tie(readedStr,nextKeyword) = this->readNextPart( __is, {"I","NOM"}, commandFinished );
        }
        else CHECK( false ) << "something wrong";


        // name given after the type
        if ( nextKeyword == "NOM" )
        {
            if ( typeOfGroup == "MAILLES" ||  typeOfGroup == "NOEUDS" )
                std::tie(readedStr,nextKeyword) = this->readNextPart( __is, {"I"}, commandFinished );
            else if ( typeOfGroup == "FACES" )
                std::tie(readedStr,nextKeyword) = this->readNextPart( __is, {"MAILL","MAILLE"}, commandFinished );

            groupName = readedStr.front();
            //std::cout << "group name : " << groupName << std::endl;
        }

        // remove quotes in groupname
        if ( groupName.size() > 2 )
        {
            if ( groupName.front() == '"' && groupName.back() == '"' )
                groupName = groupName.substr( 1, groupName.size() - 2 );
        }

        if ( typeOfGroup == "MAILLES" )
        {
            auto & curGroupMailles = M_selGroupMailles[markerId];
            std::get<0>( curGroupMailles ) = groupName;
            while ( true )
            {
                if ( nextKeyword == "I" )
                {
                    std::tie(readedStr,nextKeyword) = this->readNextPart( __is, {"I","GROUPE","GROUP"}, commandFinished );
                    std::tuple<index_type,index_type,index_type> rangeGroupMailles;
                    std::get<0>( rangeGroupMailles ) = std::stoi( readedStr[0] );
                    if (readedStr[1] == "J")
                        std::get<1>( rangeGroupMailles ) = std::stoi( readedStr[2] );
                    if (readedStr[3] == "K")
                        std::get<2>( rangeGroupMailles ) = std::stoi( readedStr[4] );
                    std::get<1>( curGroupMailles ).push_back( std::move( rangeGroupMailles ) );
                }
                else if ( nextKeyword.substr(0,5) == "GROUP" || commandFinished )
                    break;
            }
        }
        else // faces cases
        {
            std::tie(readedStr,nextKeyword) = this->readNextPart( __is, {"GROUPE","GROUP"}, commandFinished );
        }
        if ( commandFinished )
            break;
    }
    return nextKeyword;
}


} // namespace Feel

#endif
