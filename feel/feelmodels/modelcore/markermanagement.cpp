
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
MarkerManagementDirichletBC::setMarkerDirichletBCByNameId( std::string type,std::string markerNameId,std::list<std::string> const& markers,
                                                           ComponentType ct )
{
    if ( markers.empty() ) return;
    //std::cout << "set type " << type<<"\n";
    M_dirichletBCType[ct][type][markerNameId] = markers;
    this->updateForUseMarkerDirichletBC();
}
void
MarkerManagementDirichletBC::addMarkerDirichletBC(std::string type,std::string markerNameId, ComponentType ct)
{
    if ( markerNameId.empty() ) return;
    //std::cout << "add type " << type<<"\n";
    M_dirichletBCType[ct][type][markerNameId].push_back(markerNameId);
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
                    M_dirichletBCMarkersListByType[ct][ bcType].push_back( mark );
        }
    }
}

bool
MarkerManagementDirichletBC::hasMarkerDirichletBC(std::string type,ComponentType ct) const
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


std::map<std::string,std::list<std::string> > const&
MarkerManagementDirichletBC::markerDirichletBCByType( ComponentType ct ) const
{
    CHECK( M_dirichletBCMarkersListByType.find(ct) != M_dirichletBCMarkersListByType.end() ) << "invalid comp type";
    return M_dirichletBCMarkersListByType.find(ct)->second;
}

std::list<std::string> const&
MarkerManagementDirichletBC::markerDirichletBCByNameId(std::string type,std::string markerNameId,ComponentType ct ) const
{
    if ( !this->hasMarkerDirichletBC(type,ct) )
        return M_listMarkerEmpty;
    if ( M_dirichletBCType.find(ct)->second.find(type)->second.find(markerNameId) == M_dirichletBCType.find(ct)->second.find(type)->second.end() )
        return M_listMarkerEmpty;
    return M_dirichletBCType.find(ct)->second.find(type)->second.find(markerNameId)->second;
}
std::list<std::string> const&
MarkerManagementDirichletBC::markerDirichletBCelimination( ComponentType ct ) const
{
    CHECK( hasMarkerDirichletBCelimination(ct) ) << "not has type elimination\n";
    return M_dirichletBCMarkersListByType.find(ct)->second.find("elimination")->second;
}
std::list<std::string> const&
MarkerManagementDirichletBC::markerDirichletBCnitsche( ComponentType ct ) const
{
    CHECK( hasMarkerDirichletBCnitsche(ct) ) << "not has type nitsche\n";
    return M_dirichletBCMarkersListByType.find(ct)->second.find("nitsche")->second;
}
std::list<std::string> const&
MarkerManagementDirichletBC::markerDirichletBClm( ComponentType ct ) const
{
    CHECK( hasMarkerDirichletBClm(ct) ) << "not has type lm\n";
    return M_dirichletBCMarkersListByType.find(ct)->second.find("lm")->second;
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
MarkerManagementNeumannBC::setMarkerNeumannBC( NeumannBCShape shape, std::string markerNameId,std::list<std::string> const& markers )
{
    if ( markers.empty() ) return;
    M_containerMarkers[shape][markerNameId] = markers;
}
void
MarkerManagementNeumannBC::addMarkerNeumannBC(NeumannBCShape shape,std::string markerNameId)
{
    if ( markerNameId.empty() ) return;
    M_containerMarkers[shape][markerNameId].push_back(markerNameId);
}

std::map<std::string,std::list<std::string> > const&
MarkerManagementNeumannBC::markerNeumannBC( NeumannBCShape shape ) const
{
    CHECK( M_containerMarkers.find( shape ) != M_containerMarkers.end() ) << "invalid shape";
    return M_containerMarkers.find( shape )->second;
}
std::list<std::string> const&
MarkerManagementNeumannBC::markerNeumannBC( NeumannBCShape shape,std::string markerNameId ) const
{
    CHECK( M_containerMarkers.find( shape ) != M_containerMarkers.end() ) << "invalid shape";
    return M_containerMarkers.find( shape )->second.find(markerNameId)->second;
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
            if ( markNeuman.second.size() == 1 && markNeuman.second.front() == markNeuman.first ) continue;
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
MarkerManagementNeumannEulerianFrameBC::setMarkerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape, std::string markerNameId,std::list<std::string> const& markers )
{
    if ( markers.empty() ) return;
    M_containerMarkers[shape][markerNameId] = markers;
}
void
MarkerManagementNeumannEulerianFrameBC::addMarkerNeumannEulerianFrameBC(NeumannEulerianFrameBCShape shape,std::string markerNameId)
{
    if ( markerNameId.empty() ) return;
    M_containerMarkers[shape][markerNameId].push_back(markerNameId);
}

std::map<std::string,std::list<std::string> > const&
MarkerManagementNeumannEulerianFrameBC::markerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape ) const
{
    CHECK( M_containerMarkers.find( shape ) != M_containerMarkers.end() ) << "invalid shape";
    return M_containerMarkers.find( shape )->second;
}
std::list<std::string> const&
MarkerManagementNeumannEulerianFrameBC::markerNeumannEulerianFrameBC( NeumannEulerianFrameBCShape shape,std::string markerNameId ) const
{
    CHECK( M_containerMarkers.find( shape ) != M_containerMarkers.end() ) << "invalid shape";
    return M_containerMarkers.find( shape )->second.find(markerNameId)->second;
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
            if ( markNeuman.second.size() == 1 && markNeuman.second.front() == markNeuman.first ) continue;
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
MarkerManagementALEMeshBC::setMarkerALEMeshBC( std::string type, std::list<std::string> const& markers )
{
    CHECK( type == "fixed" || type == "moving" || type == "free" ) << "error ALE type " << type;
    M_containerMarkers[type] = markers;
}
void
MarkerManagementALEMeshBC::addMarkerALEMeshBC(std::string type, std::string markerName)
{
    CHECK( type == "fixed" || type == "moving" || type == "free" ) << "error ALE type " << type;
    if ( std::find( M_containerMarkers[type].begin(),M_containerMarkers[type].end(),markerName) == M_containerMarkers[type].end() )
        M_containerMarkers[type].push_back(markerName);
}

std::map<std::string,std::list<std::string> > const&
MarkerManagementALEMeshBC::markerALEMeshBC() const
{
    return M_containerMarkers;
}
std::list<std::string> const&
MarkerManagementALEMeshBC::markerALEMeshBC( std::string type ) const
{
    CHECK( type == "fixed" || type == "moving" || type == "free" ) << "error ALE type " << type;
    if ( M_containerMarkers.find(type) == M_containerMarkers.end() )
        return M_listMarkerEmpty;
    return M_containerMarkers.find(type)->second;
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
MarkerManagementSlipBC::setMarkerSlipBC( std::list<std::string> const& markers )
{
    M_containerMarkers = markers;
}
void
MarkerManagementSlipBC::addMarkerSlipBC( std::string markerName )
{
    if ( std::find( M_containerMarkers.begin(),M_containerMarkers.end(),markerName) == M_containerMarkers.end() )
        M_containerMarkers.push_back(markerName);
}
std::list<std::string> const&
MarkerManagementSlipBC::markerSlipBC() const
{
    return M_containerMarkers;
}
std::string
MarkerManagementSlipBC::getInfoSlipBC() const
{
    std::ostringstream _ostr;
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
MarkerManagementPressureBC::setMarkerPressureBC( std::string const& markerNameId, std::list<std::string> const& markers )
{
    if ( markerNameId.empty() ) return;
    M_containerMarkers[markerNameId] = markers;
    for ( std::string const& markerName : markers )
    {
        if ( std::find( M_listMarkers.begin(),M_listMarkers.end(),markerName) == M_listMarkers.end() )
            M_listMarkers.push_back( markerName );
    }
}
void
MarkerManagementPressureBC::addMarkerPressureBC( std::string const& markerName )
{
    if ( markerName.empty() ) return;
    M_containerMarkers[markerName].push_back(markerName);
    if ( std::find( M_listMarkers.begin(),M_listMarkers.end(),markerName) == M_listMarkers.end() )
        M_listMarkers.push_back( markerName );
}
std::list<std::string> const&
MarkerManagementPressureBC::markerPressureBC() const
{
    return M_listMarkers;
}
std::list<std::string> const&
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
std::string
MarkerManagementPressureBC::getInfoPressureBC() const
{
    std::ostringstream _ostr;
    if ( M_containerMarkers.empty() )
            return _ostr.str();

    for ( auto const& markerBase : M_containerMarkers )
    {
        _ostr << "\n       -- Pressure Dirichlet : " << markerBase.first;
        if ( markerBase.second.size() == 1 && markerBase.second.front() == markerBase.first ) continue;
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
MarkerManagementRobinBC::setMarkerRobinBC( std::string const& markerNameId, std::list<std::string> const& markers )
{
    M_containerMarkers[markerNameId] = markers;
}
void
MarkerManagementRobinBC::addMarkerRobinBC( std::string const& markerNameId )
{
    if ( markerNameId.empty() ) return;
    M_containerMarkers[markerNameId].push_back(markerNameId);
}
std::map<std::string,std::list<std::string> > const&
MarkerManagementRobinBC::markerRobinBC() const
{
    return M_containerMarkers;
}
std::list<std::string> const&
MarkerManagementRobinBC::markerRobinBC( std::string const& markerNameId ) const
{
    if ( M_containerMarkers.find( markerNameId ) != M_containerMarkers.end() )
        return M_containerMarkers.find(markerNameId)->second;
    else
        return M_listMarkerEmpty;
}
std::string
MarkerManagementRobinBC::getInfoRobinBC() const
{
    std::ostringstream _ostr;
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
MarkerManagementFluidStructureInterfaceBC::setMarkerFluidStructureInterfaceBC( std::list<std::string> const& markers )
{
    M_containerMarkers = markers;
}
void
MarkerManagementFluidStructureInterfaceBC::addMarkerFluidStructureInterfaceBC( std::string markerName )
{
    if ( std::find( M_containerMarkers.begin(),M_containerMarkers.end(),markerName) == M_containerMarkers.end() )
        M_containerMarkers.push_back(markerName);
}
std::list<std::string> const&
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

} // namespace FeelModels
} // namespace Feel
