
#include <feel/feelmodels/modelcore/markermanagement.hpp>

namespace Feel {
namespace FeelModels {

MarkerManagementDirichletBC::MarkerManagementDirichletBC()
    :
    M_dirichletBCType(),
    M_dirichletBCMarkersListByType(),
    M_listMarkerEmpty()
{
    /*M_dirichletBCType[ComponentType::NO_COMPONENT].clear();
    M_dirichletBCType[ComponentType::X].clear();
    M_dirichletBCType[ComponentType::Y].clear();
    M_dirichletBCType[ComponentType::Z].clear();
    M_dirichletBCMarkersListByType[ComponentType::NO_COMPONENT].clear();
    M_dirichletBCMarkersListByType[ComponentType::X].clear();
    M_dirichletBCMarkersListByType[ComponentType::Y].clear();
     M_dirichletBCMarkersListByType[ComponentType::Z].clear();*/
}

void
MarkerManagementDirichletBC::clearMarkerDirichletBC()
{
    M_dirichletBCType.clear();
    M_dirichletBCMarkersListByType.clear();
}

void
MarkerManagementDirichletBC::setMarkerDirichletBCByNameId( std::string const& type, std::string const& name,
                                                           std::set<std::string> const& markers,
                                                           ComponentType ct )
{
    if ( markers.empty() ) return;
    //std::cout << "set type " << type<<"\n";
    M_dirichletBCType[ct][type][name] = markers;
    this->updateForUseMarkerDirichletBC();
}
void
MarkerManagementDirichletBC::addMarkerDirichletBC(std::string const& type, std::string const& name,
                                                  std::string const& marker, ComponentType ct)
{
    if ( name.empty() ) return;
    //std::cout << "add type " << type<<"\n";
    M_dirichletBCType[ct][type][name].insert(marker);
    this->updateForUseMarkerDirichletBC();
}
void
MarkerManagementDirichletBC::addMarkerDirichletBC(std::string const& type, std::string const& name,
                                                  std::set<std::string> const& markers,
                                                  ComponentType ct)
{
    if ( name.empty() ) return;
    //std::cout << "add type " << type<<"\n";
    for( auto const& m : markers )
        M_dirichletBCType[ct][type][name].insert(m);
    this->updateForUseMarkerDirichletBC();
}

void
MarkerManagementDirichletBC::updateForUseMarkerDirichletBC()
{
    M_dirichletBCMarkersListByType.clear();
    for ( auto const& bykindCompType : M_dirichletBCType )
    {
        ComponentType ct = bykindCompType.first;
        for ( auto const& bykindDirichlet : bykindCompType.second )
        {
            std::string bcType = bykindDirichlet.first;
            for ( auto const& markId : bykindDirichlet.second )
                for ( auto const& mark : markId.second )
                    M_dirichletBCMarkersListByType[ct][ bcType].insert( mark );
        }
    }
}

bool
MarkerManagementDirichletBC::hasMarkerDirichletBC(std::string const& type,ComponentType ct) const
{
    auto const& itFindCompType = M_dirichletBCType.find(ct);
    if ( itFindCompType == M_dirichletBCType.end() )
        return false;
    auto const& itFindDirichletType = itFindCompType->second.find(type);
    if ( itFindDirichletType == itFindCompType->second.end() )
        return false;
    if ( itFindDirichletType->second.empty() )
        return false;
    return true;
}
bool
MarkerManagementDirichletBC::hasMarkerDirichletBCelimination( ComponentType ct ) const { return this->hasMarkerDirichletBC("elimination",ct); }
bool
MarkerManagementDirichletBC::hasMarkerDirichletBCnitsche( ComponentType ct ) const { return this->hasMarkerDirichletBC("nitsche",ct); }
bool
MarkerManagementDirichletBC::hasMarkerDirichletBClm( ComponentType ct ) const { return this->hasMarkerDirichletBC("lm",ct); }
bool
MarkerManagementDirichletBC::hasMarkerDirichletBCelimination() const
{
    return   this->hasMarkerDirichletBCelimination( ComponentType::NO_COMPONENT ) ||
        /**/ this->hasMarkerDirichletBCelimination( ComponentType::X ) ||
        /**/ this->hasMarkerDirichletBCelimination( ComponentType::Y ) ||
        /**/ this->hasMarkerDirichletBCelimination( ComponentType::Z );
}
bool
MarkerManagementDirichletBC::hasMarkerDirichletBCnitsche() const
{
    return   this->hasMarkerDirichletBCnitsche( ComponentType::NO_COMPONENT ) ||
        /**/ this->hasMarkerDirichletBCnitsche( ComponentType::X ) ||
        /**/ this->hasMarkerDirichletBCnitsche( ComponentType::Y ) ||
        /**/ this->hasMarkerDirichletBCnitsche( ComponentType::Z );
}
bool
MarkerManagementDirichletBC::hasMarkerDirichletBClm() const
{
    return    this->hasMarkerDirichletBClm( ComponentType::NO_COMPONENT ) ||
        /**/  this->hasMarkerDirichletBClm( ComponentType::X ) ||
        /**/  this->hasMarkerDirichletBClm( ComponentType::Y ) ||
        /**/  this->hasMarkerDirichletBClm( ComponentType::Z );
}


std::map<std::string,std::set<std::string> > const&
MarkerManagementDirichletBC::markerDirichletBCByType( ComponentType ct ) const
{
    CHECK( M_dirichletBCMarkersListByType.find(ct) != M_dirichletBCMarkersListByType.end() ) << "invalid comp type";
    return M_dirichletBCMarkersListByType.find(ct)->second;
}

std::set<std::string> const&
MarkerManagementDirichletBC::markerDirichletBCByNameId(std::string const& type,std::string const& markerNameId,ComponentType ct ) const
{
    if ( !this->hasMarkerDirichletBC(type,ct) )
        return M_listMarkerEmpty;
    if ( M_dirichletBCType.find(ct)->second.find(type)->second.find(markerNameId) == M_dirichletBCType.find(ct)->second.find(type)->second.end() )
        return M_listMarkerEmpty;
    return M_dirichletBCType.find(ct)->second.find(type)->second.find(markerNameId)->second;
}
std::set<std::string> const&
MarkerManagementDirichletBC::markerDirichletBCelimination( ComponentType ct ) const
{
    CHECK( hasMarkerDirichletBCelimination(ct) ) << "not has type elimination\n";
    return M_dirichletBCMarkersListByType.find(ct)->second.find("elimination")->second;
}
std::set<std::string> const&
MarkerManagementDirichletBC::markerDirichletBCnitsche( ComponentType ct ) const
{
    CHECK( hasMarkerDirichletBCnitsche(ct) ) << "not has type nitsche\n";
    return M_dirichletBCMarkersListByType.find(ct)->second.find("nitsche")->second;
}
std::set<std::string> const&
MarkerManagementDirichletBC::markerDirichletBClm( ComponentType ct ) const
{
    CHECK( hasMarkerDirichletBClm(ct) ) << "not has type lm\n";
    return M_dirichletBCMarkersListByType.find(ct)->second.find("lm")->second;
}


void
MarkerManagementDirichletBC::updateInformationObjectDirichletBC( pt::ptree & p ) const
{
    for ( auto const& bykindbase : M_dirichletBCType )
    {
        ComponentType ct = bykindbase.first;
        std::string ctStr;
        if ( ct==ComponentType::X ) ctStr = "[X]";
        else if ( ct==ComponentType::Y ) ctStr = "[Y]";
        else if ( ct==ComponentType::Z ) ctStr = "[Z]";
        //for ( auto const& bykind : bykindbase.second/*M_dirichletBCType*/ )
        //ComponentType ct = ComponentType::NO_COMPONENT;
        if ( bykindbase.second.size()>0 )
        {
            for ( auto itBC = this->markerDirichletBCByType(ct).begin(), enBC = this->markerDirichletBCByType(ct).end() ; itBC!=enBC ; ++itBC )
            {
                std::ostringstream _ostr;
                _ostr << "Dirichlet"<< ctStr <<" (" << itBC->first << ")";
                pt::ptree ptTmp;
                for ( auto itMark = itBC->second.begin(), enMark = itBC->second.end() ; itMark!= enMark ; ++itMark)
                    ptTmp.push_back( std::make_pair("", pt::ptree( *itMark ) ) );
                if ( !ptTmp.empty() )
                    p.put_child( _ostr.str(), ptTmp );
            }
        }
    }
}

std::string
MarkerManagementDirichletBC::getInfoDirichletBC() const
{
    std::ostringstream _ostr;

    for ( auto const& bykindbase : M_dirichletBCType )
    {
        ComponentType ct = bykindbase.first;
        std::string ctStr;
        if ( ct==ComponentType::X ) ctStr = "[X]";
        else if ( ct==ComponentType::Y ) ctStr = "[Y]";
        else if ( ct==ComponentType::Z ) ctStr = "[Z]";
        //for ( auto const& bykind : bykindbase.second/*M_dirichletBCType*/ )
        //ComponentType ct = ComponentType::NO_COMPONENT;
        if ( bykindbase.second.size()>0 )
        {
            for ( auto itBC = this->markerDirichletBCByType(ct).begin(), enBC = this->markerDirichletBCByType(ct).end() ; itBC!=enBC ; ++itBC )
            {
                _ostr << "\n       -- Dirichlet"<< ctStr <<" (" << itBC->first << ") : ";
                int cptMark = 0;
                for ( auto itMark = itBC->second.begin(), enMark = itBC->second.end() ; itMark!= enMark ; ++itMark, ++cptMark)
                {
                    if ( cptMark > 0) _ostr << " , ";
                    _ostr << *itMark;
                }
            }
        }
    }
    return _ostr.str();
}

//-------------------------------------------------------------------------//
//-------------------------------------------------------------------------//
//-------------------------------------------------------------------------//

MarkerManagementNeumannBC::MarkerManagementNeumannBC()
    :
    M_containerMarkers(),
    M_listMarkerEmpty()
{}

void
MarkerManagementNeumannBC::clearMarkerNeumannBC()
{
    M_containerMarkers.clear();
}

void
MarkerManagementNeumannBC::setMarkerNeumannBC( NeumannBCShape shape, std::string const& name,std::set<std::string> const& markers )
{
    if ( markers.empty() ) return;
    M_containerMarkers[shape][name] = markers;
}
void
MarkerManagementNeumannBC::addMarkerNeumannBC(NeumannBCShape shape,std::string const& name,std::string const& marker)
{
    if ( name.empty() ) return;
    M_containerMarkers[shape][name].insert(marker);
}
void
MarkerManagementNeumannBC::addMarkerNeumannBC(NeumannBCShape shape,std::string const& name,std::set<std::string> const& markers)
{
    if ( name.empty() ) return;
    M_containerMarkers[shape][name].insert(markers.begin(),markers.end());
}

std::map<std::string,std::set<std::string> > const&
MarkerManagementNeumannBC::markerNeumannBC( NeumannBCShape shape ) const
{
    CHECK( M_containerMarkers.find( shape ) != M_containerMarkers.end() ) << "invalid shape";
    return M_containerMarkers.find( shape )->second;
}
std::set<std::string> const&
MarkerManagementNeumannBC::markerNeumannBC( NeumannBCShape shape,std::string const& markerNameId ) const
{
    CHECK( M_containerMarkers.find( shape ) != M_containerMarkers.end() ) << "invalid shape";
    return M_containerMarkers.find( shape )->second.find(markerNameId)->second;
}

void
MarkerManagementNeumannBC::updateInformationObjectNeumannBC( pt::ptree & p ) const
{
   for ( auto const& markNeumanBase : M_containerMarkers )
    {
        std::string shapeStr = (markNeumanBase.first == NeumannBCShape::SCALAR )? "[scalar]":
            (markNeumanBase.first == NeumannBCShape::VECTORIAL )? "[vectorial]":"[tensor2]";
        std::ostringstream _ostr;
        _ostr << "Neumann" << shapeStr;
        pt::ptree ptTmp;
        for ( auto const& markNeuman : markNeumanBase.second )
        {
            //if ( markNeuman.second.size() == 1 && *markNeuman.second.begin() == markNeuman.first ) continue;
            for ( auto itMark = markNeuman.second.begin(), enMark = markNeuman.second.end() ; itMark!=enMark ; ++itMark )
                ptTmp.push_back( std::make_pair("", pt::ptree( *itMark ) ) );
        }
        if ( !ptTmp.empty() )
            p.put_child( _ostr.str(), ptTmp );
    }
}

std::string
MarkerManagementNeumannBC::getInfoNeumannBC() const
{
    std::ostringstream _ostr;

    for ( auto const& markNeumanBase : M_containerMarkers )
    {
        std::string shapeStr = (markNeumanBase.first == NeumannBCShape::SCALAR )? "[scalar]":
            (markNeumanBase.first == NeumannBCShape::VECTORIAL )? "[vectorial]":"[tensor2]";
        for ( auto const& markNeuman : markNeumanBase.second )
        {
            _ostr << "\n       -- Neumann" << shapeStr << " : " << markNeuman.first;
            if ( markNeuman.second.size() == 1 && *markNeuman.second.begin() == markNeuman.first ) continue;
            _ostr << " -> (";
            int cptMark = 0;
            for ( auto itMark = markNeuman.second.begin(), enMark = markNeuman.second.end() ; itMark!=enMark ; ++itMark,++cptMark )
            {
                if ( cptMark > 0) _ostr << " , ";
                _ostr << *itMark;
            }
            _ostr << ")";
        }
    }
    return _ostr.str();
}

//-------------------------------------------------------------------------//
//-------------------------------------------------------------------------//
//-------------------------------------------------------------------------//

MarkerManagementNeumannEulerianFrameBC::MarkerManagementNeumannEulerianFrameBC()
    :
    M_containerMarkers(),
    M_listMarkerEmpty()
{}

void
MarkerManagementNeumannEulerianFrameBC::clearMarkerNeumannEulerianFrameBC()
{
    M_containerMarkers.clear();
}

void
MarkerManagementNeumannEulerianFrameBC::setMarkerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape, std::string const& name,std::set<std::string> const& markers )
{
    if ( markers.empty() ) return;
    M_containerMarkers[shape][name] = markers;
}
void
MarkerManagementNeumannEulerianFrameBC::addMarkerNeumannEulerianFrameBC(NeumannEulerianFrameBCShape shape,std::string const& name,std::string const& marker)
{
    if ( name.empty() ) return;
    M_containerMarkers[shape][name].insert(marker);
}
void
MarkerManagementNeumannEulerianFrameBC::addMarkerNeumannEulerianFrameBC(NeumannEulerianFrameBCShape shape,std::string const& name,std::set<std::string> const& markers)
{
    if ( name.empty() ) return;
    M_containerMarkers[shape][name].insert(markers.begin(),markers.end());
}

std::map<std::string,std::set<std::string> > const&
MarkerManagementNeumannEulerianFrameBC::markerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape ) const
{
    CHECK( M_containerMarkers.find( shape ) != M_containerMarkers.end() ) << "invalid shape";
    return M_containerMarkers.find( shape )->second;
}
std::set<std::string> const&
MarkerManagementNeumannEulerianFrameBC::markerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape,std::string const& markerNameId ) const
{
    CHECK( M_containerMarkers.find( shape ) != M_containerMarkers.end() ) << "invalid shape";
    return M_containerMarkers.find( shape )->second.find(markerNameId)->second;
}

void
MarkerManagementNeumannEulerianFrameBC::updateInformationObjectNeumannEulerianFrameBC( pt::ptree & p ) const
{
    // TODO
}

std::string
MarkerManagementNeumannEulerianFrameBC::getInfoNeumannEulerianFrameBC() const
{
    std::ostringstream _ostr;

    for ( auto const& markNeumanBase : M_containerMarkers )
    {
        std::string shapeStr = (markNeumanBase.first == NeumannEulerianFrameBCShape::SCALAR )? "[scalar]":
            (markNeumanBase.first == NeumannEulerianFrameBCShape::VECTORIAL )? "[vectorial]":"[tensor2]";
        for ( auto const& markNeuman : markNeumanBase.second )
        {
            _ostr << "\n       -- Neumann (EulerianFrame) " << shapeStr << " : " << markNeuman.first;
            if ( markNeuman.second.size() == 1 && *markNeuman.second.begin() == markNeuman.first ) continue;
            _ostr << " -> (";
            int cptMark = 0;
            for ( auto itMark = markNeuman.second.begin(), enMark = markNeuman.second.end() ; itMark!=enMark ; ++itMark,++cptMark )
            {
                if ( cptMark > 0) _ostr << " , ";
                _ostr << *itMark;
            }
            _ostr << ")";
        }
    }
    return _ostr.str();
}

//-------------------------------------------------------------------------//
//-------------------------------------------------------------------------//
//-------------------------------------------------------------------------//

MarkerManagementALEMeshBC::MarkerManagementALEMeshBC()
    :
    M_containerMarkers(),
    M_listMarkerEmpty()
{}

void
MarkerManagementALEMeshBC::clearMarkerALEMeshBC()
{
    M_containerMarkers.clear();
}

void
MarkerManagementALEMeshBC::setMarkerALEMeshBC( std::string const& type, std::set<std::string> const& markers )
{
    CHECK( type == "fixed" || type == "moving" || type == "free" ) << "error ALE type " << type;
    M_containerMarkers[type] = markers;
}
void
MarkerManagementALEMeshBC::addMarkerALEMeshBC(std::string const& type, std::string const& markerName)
{
    CHECK( type == "fixed" || type == "moving" || type == "free" ) << "error ALE type " << type;
    M_containerMarkers[type].insert(markerName);
}
void
MarkerManagementALEMeshBC::addMarkerALEMeshBC( std::string const& type, std::set<std::string> const& markers )
{
    CHECK( type == "fixed" || type == "moving" || type == "free" ) << "error ALE type " << type;
    M_containerMarkers[type].insert( markers.begin(),markers.end() );
}

std::map<std::string,std::set<std::string> > const&
MarkerManagementALEMeshBC::markerALEMeshBC() const
{
    return M_containerMarkers;
}
std::set<std::string> const&
MarkerManagementALEMeshBC::markerALEMeshBC( std::string const& type ) const
{
    CHECK( type == "fixed" || type == "moving" || type == "free" ) << "error ALE type " << type;
    if ( M_containerMarkers.find(type) == M_containerMarkers.end() )
        return M_listMarkerEmpty;
    return M_containerMarkers.find(type)->second;
}

void
MarkerManagementALEMeshBC::updateInformationObjectALEMeshBC( pt::ptree & p ) const
{
    // TODO
}

std::string
MarkerManagementALEMeshBC::getInfoALEMeshBC() const
{
    std::ostringstream _ostr;
    for ( auto const& markAleMeshBase : M_containerMarkers )
    {
        std::string shapeStr = "["+ markAleMeshBase.first +"]";
        if ( markAleMeshBase.second.empty() ) continue;
        _ostr << "\n       -- MeshALE" << shapeStr << " : ";
        bool doFirstMarker = false;
        for ( auto const& markAleMesh : markAleMeshBase.second )
        {
            if ( !doFirstMarker)
            {
                _ostr << markAleMesh;
                doFirstMarker=true;
            }
            else
                _ostr << " , " << markAleMesh;
        }
    }
    return _ostr.str();

}

//--------------------------------------------------------------//

MarkerManagementSlipBC::MarkerManagementSlipBC()
    :
    M_containerMarkers(),
    M_listMarkerEmpty()
{}
void
MarkerManagementSlipBC::clearMarkerSlipBC()
{
    M_containerMarkers.clear();
}
void
MarkerManagementSlipBC::setMarkerSlipBC( std::set<std::string> const& markers )
{
    M_containerMarkers = markers;
}
void
MarkerManagementSlipBC::addMarkerSlipBC( std::string const& markerName )
{
    M_containerMarkers.insert(markerName);
}
void
MarkerManagementSlipBC::addMarkerSlipBC( std::set<std::string> const& markers )
{
    M_containerMarkers = markers;
}
std::set<std::string> const&
MarkerManagementSlipBC::markerSlipBC() const
{
    return M_containerMarkers;
}

void
MarkerManagementSlipBC::updateInformationObjectSlipBC( pt::ptree & p ) const
{
    // TODO
}

std::string
MarkerManagementSlipBC::getInfoSlipBC() const
{
    std::ostringstream _ostr;

    if( M_containerMarkers.size() > 0 )
    {
        _ostr << "\n       -- slip : ";
        int iMark = 0;
        for ( auto const& markSlip : M_containerMarkers )
        {
            if( iMark > 0 ) _ostr << " , ";
            _ostr << markSlip;
            ++iMark;
        }
    }
    return _ostr.str();
}

//--------------------------------------------------------------//

MarkerManagementPressureBC::MarkerManagementPressureBC()
    :
    M_containerMarkers(),
    M_listMarkers(),
    M_listMarkerEmpty()
{}
void
MarkerManagementPressureBC::clearMarkerPressureBC()
{
    M_containerMarkers.clear();
    M_listMarkers.clear();
}
void
MarkerManagementPressureBC::setMarkerPressureBC( std::string const& name, std::set<std::string> const& markers )
{
    if ( name.empty() ) return;
    M_containerMarkers[name] = markers;
    for ( std::string const& markerName : markers )
    {
        if ( std::find( M_listMarkers.begin(),M_listMarkers.end(),markerName) == M_listMarkers.end() )
            M_listMarkers.insert( markerName );
    }
}
void
MarkerManagementPressureBC::addMarkerPressureBC( std::string const& name,std::string const& marker )
{
    if ( name.empty() ) return;
    M_containerMarkers[name].insert(marker);
    if ( std::find( M_listMarkers.begin(),M_listMarkers.end(),name) == M_listMarkers.end() )
        M_listMarkers.insert( marker );
}
void
MarkerManagementPressureBC::addMarkerPressureBC( std::string const& name,std::set<std::string> const& markers )
{
    if ( name.empty() ) return;
    for(auto const& m : markers )
        M_containerMarkers[name].insert(m);
    if ( std::find( M_listMarkers.begin(),M_listMarkers.end(),name) == M_listMarkers.end() )
        M_listMarkers = markers;
}
std::set<std::string> const&
MarkerManagementPressureBC::markerPressureBC() const
{
    return M_listMarkers;
}
std::set<std::string> const&
MarkerManagementPressureBC::markerPressureBC( std::string const& markerNameId ) const
{
    auto itFind = M_containerMarkers.find( markerNameId );
    if ( itFind != M_containerMarkers.end() )
        return itFind->second;
    else
        return M_listMarkerEmpty;
}
bool
MarkerManagementPressureBC::hasMarkerPressureBC() const
{
    return !M_containerMarkers.empty();
}

void
MarkerManagementPressureBC::updateInformationObjectPressureBC( pt::ptree & p ) const
{
    // TODO
}

std::string
MarkerManagementPressureBC::getInfoPressureBC() const
{
    std::ostringstream _ostr;
    if ( M_containerMarkers.empty() )
            return _ostr.str();

    for ( auto const& markerBase : M_containerMarkers )
    {
        _ostr << "\n       -- Pressure Dirichlet : " << markerBase.first;
        if ( markerBase.second.size() == 1 && *markerBase.second.begin() == markerBase.first ) continue;
        _ostr << " -> (";
        int cptMark = 0;
        for ( auto itMark = markerBase.second.begin(), enMark = markerBase.second.end() ; itMark!=enMark ; ++itMark,++cptMark )
        {
            if ( cptMark > 0) _ostr << " , ";
            _ostr << *itMark;
        }
        _ostr << ")";
    }

    return _ostr.str();
}

//--------------------------------------------------------------//

MarkerManagementRobinBC::MarkerManagementRobinBC()
    :
    M_containerMarkers(),
    M_listMarkerEmpty()
{}
void
MarkerManagementRobinBC::clearMarkerRobinBC()
{
    M_containerMarkers.clear();
}
void
MarkerManagementRobinBC::setMarkerRobinBC( std::string const& name, std::set<std::string> const& markers )
{
    M_containerMarkers[name] = markers;
}
void
MarkerManagementRobinBC::addMarkerRobinBC( std::string const& name, std::string const& marker )
{
    if ( name.empty() ) return;
    M_containerMarkers[name].insert(marker);
}
void
MarkerManagementRobinBC::addMarkerRobinBC( std::string const& name, std::set<std::string> const& markers )
{
    if ( name.empty() ) return;
    for(auto const& m : markers )
        M_containerMarkers[name].insert(m);
}
std::map<std::string,std::set<std::string> > const&
MarkerManagementRobinBC::markerRobinBC() const
{
    return M_containerMarkers;
}
std::set<std::string> const&
MarkerManagementRobinBC::markerRobinBC( std::string const& markerNameId ) const
{
    if ( M_containerMarkers.find( markerNameId ) != M_containerMarkers.end() )
        return M_containerMarkers.find(markerNameId)->second;
    else
        return M_listMarkerEmpty;
}

void
MarkerManagementRobinBC::updateInformationObjectRobinBC( pt::ptree & p ) const
{
    // TODO
}
std::string
MarkerManagementRobinBC::getInfoRobinBC() const
{
    std::ostringstream _ostr;

    for ( auto const& markerBase : M_containerMarkers )
    {
        _ostr << "\n       -- Robin : " << markerBase.first;
        if ( markerBase.second.size() == 1 && *markerBase.second.begin() == markerBase.first ) continue;
        _ostr << " -> (";
        int cptMark = 0;
        for ( auto itMark = markerBase.second.begin(), enMark = markerBase.second.end() ; itMark!=enMark ; ++itMark,++cptMark )
        {
            if ( cptMark > 0) _ostr << " , ";
            _ostr << *itMark;
        }
        _ostr << ")";
    }

    return _ostr.str();
}

//--------------------------------------------------------------//

MarkerManagementFluidStructureInterfaceBC::MarkerManagementFluidStructureInterfaceBC()
    :
    M_containerMarkers(),
    M_listMarkerEmpty()
{}
void
MarkerManagementFluidStructureInterfaceBC::clearMarkerFluidStructureInterfaceBC()
{
    M_containerMarkers.clear();
}
void
MarkerManagementFluidStructureInterfaceBC::setMarkerFluidStructureInterfaceBC( std::set<std::string> const& markers )
{
    M_containerMarkers = markers;
}
void
MarkerManagementFluidStructureInterfaceBC::addMarkerFluidStructureInterfaceBC( std::string const& markerName )
{
    if ( std::find( M_containerMarkers.begin(),M_containerMarkers.end(),markerName) == M_containerMarkers.end() )
        M_containerMarkers.insert(markerName);
}
void
MarkerManagementFluidStructureInterfaceBC::addMarkerFluidStructureInterfaceBC( std::set<std::string> const& markers )
{
    M_containerMarkers = markers;
}
std::set<std::string> const&
MarkerManagementFluidStructureInterfaceBC::markerFluidStructureInterfaceBC() const
{
    return M_containerMarkers;
}
std::string
MarkerManagementFluidStructureInterfaceBC::getInfoFluidStructureInterfaceBC() const
{
    std::ostringstream _ostr;
    return _ostr.str();
}

void
MarkerManagementFluidStructureInterfaceBC::updateInformationObjectFluidStructureInterfaceBC( pt::ptree & p ) const
{
    // TODO
}

//--------------------------------------------------------------//

MarkerManagementIntegralBC::MarkerManagementIntegralBC()
    :
    M_containerMarkers(),
    M_listMarkerEmpty()
{}
void
MarkerManagementIntegralBC::clearMarkerIntegralBC()
{
    M_containerMarkers.clear();
}
void
MarkerManagementIntegralBC::setMarkerIntegralBC( std::string const& name, std::set<std::string> const& markers )
{
    M_containerMarkers[name] = markers;
}
void
MarkerManagementIntegralBC::addMarkerIntegralBC( std::string const& name, std::string const& marker )
{
    if ( name.empty() ) return;
    M_containerMarkers[name].insert(marker);
}
void
MarkerManagementIntegralBC::addMarkerIntegralBC( std::string const& name, std::set<std::string> const& markers )
{
    if ( name.empty() ) return;
    for(auto const& m : markers )
        M_containerMarkers[name].insert(m);
}
std::map<std::string,std::set<std::string> > const&
MarkerManagementIntegralBC::markerIntegralBC() const
{
    return M_containerMarkers;
}
std::set<std::string> const&
MarkerManagementIntegralBC::markerIntegralBC( std::string const& markerNameId ) const
{
    if ( M_containerMarkers.find( markerNameId ) != M_containerMarkers.end() )
        return M_containerMarkers.find(markerNameId)->second;
    else
        return M_listMarkerEmpty;
}

void
MarkerManagementIntegralBC::updateInformationObjectIntegralBC( pt::ptree & p ) const
{
    // TODO
}
std::string
MarkerManagementIntegralBC::getInfoIntegralBC() const
{
    std::ostringstream _ostr;

    for ( auto const& markerBase : M_containerMarkers )
    {
        _ostr << "\n       -- Integral : " << markerBase.first;
        if ( markerBase.second.size() == 1 && *markerBase.second.begin() == markerBase.first ) continue;
        _ostr << " -> (";
        int cptMark = 0;
        for ( auto itMark = markerBase.second.begin(), enMark = markerBase.second.end() ; itMark!=enMark ; ++itMark,++cptMark )
        {
            if ( cptMark > 0) _ostr << " , ";
            _ostr << *itMark;
        }
        _ostr << ")";
    }

    return _ostr.str();
}



} // namespace FeelModels
} // namespace Feel
