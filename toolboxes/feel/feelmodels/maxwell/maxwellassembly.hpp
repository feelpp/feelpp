#ifndef FEELPP_TOOLBOXES_MAXWELL_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_MAXWELL_ASSEMBLY_HPP 1

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType>
template< typename SymbolsExpr>
void
Maxwell<ConvexType>::updateLinearPDE( DataUpdateLinear & data, SymbolsExpr const& symbolsExpr ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    // bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    std::string sc=(buildCstPart)?" (build cst part)":" (build non cst part)";
    this->log("Maxwell","updateLinearPDE", "start"+sc);
    boost::mpi::timer thetimer;

    auto mesh = this->mesh();
    auto XhV = this->spaceMagneticPotential();
    auto const& v = this->fieldMagneticPotential();

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix() ,
                                              _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=XhV, _vector=F,
                               _rowstart=this->rowStartInVector() );

    for ( auto const& rangeData : M_maxwellProperties->rangeMeshElementsByMaterial() )
    {
        std::string const& matName = rangeData.first;
        auto const& range = rangeData.second;
        auto const& magneticPermeability = M_maxwellProperties->magneticPermeability( matName );
        auto mu = expr( magneticPermeability.expr(), symbolsExpr);
        bool buildDiffusion = mu.expression().isConstant()? buildCstPart : buildNonCstPart;
        if( buildDiffusion )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=range,
                           _expr=inner(curlt(v),curl(v))/mu + M_epsilon*inner(idt(v),id(v)),
                           _geomap=this->geomap() );
        }
    }

    // update source term
    if ( buildNonCstPart )
    {
        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto rangeEltUsed = (markers(d).empty())? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            myLinearForm +=
                integrate( _range=rangeEltUsed,
                           _expr= inner(expression(d, symbolsExpr),id(v)),
                           _geomap=this->geomap() );
        }
    }

    // update bc
    if ( buildNonCstPart )
        this->updateLinearPDEWeakBC( data, symbolsExpr, hana::int_<nDim>() );
}

template< typename ConvexType>
template< typename SymbolsExpr>
void
Maxwell<ConvexType>::updateLinearPDEWeakBC( DataUpdateLinear & data, SymbolsExpr const& symbolsExpr, hana::int_<3> ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto mesh = this->mesh();
    auto XhV = this->spaceMagneticPotential();
    auto const& u = this->fieldMagneticPotential();

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix() ,
                                              _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=XhV, _vector=F,
                               _rowstart=this->rowStartInVector() );

    auto bcs = this->modelProperties().boundaryConditions2().boundaryConditions("magnetic-potential", "Dirichlet");
    for( auto const& bcpair : bcs )
    {
        auto bc = bcpair.second;
        auto markers = bc.markers();
        auto const& magneticPermeability = M_maxwellProperties->magneticPermeability( bc.material() );
        auto mu = expr( magneticPermeability.expr(), symbolsExpr);
        // Feel::cout << "Dirichlet on " << bcpair.first
        //            << " with " << markers.size() << " markers" << std::endl;
        auto Ad = bc.template expr<3,1>();

        bilinearForm_PatternCoupled +=
            integrate(_range=markedfaces(M_mesh, markers),
                      _expr=1e5/hFace()*inner(cross(id(u),N()),cross(idt(u),N()))/mu );
        bilinearForm_PatternCoupled +=
            integrate(_range=markedfaces(M_mesh, markers),
                      _expr=inner(curl(u),cross(idt(u),N()))/mu );
        bilinearForm_PatternCoupled +=
            integrate(_range=markedfaces(M_mesh, markers),
                      _expr=inner(curlt(u),cross(id(u),N()))/mu );
        myLinearForm += integrate(_range=markedfaces(M_mesh, markers), _expr=inner(curl(u),Ad)/mu );
        myLinearForm += integrate(_range=markedfaces(M_mesh, markers),
                                  _expr=1e5/hFace()*inner(cross(id(u),N()),Ad)/mu );
    }
}

template< typename ConvexType>
template< typename SymbolsExpr>
void
Maxwell<ConvexType>::updateLinearPDEWeakBC( DataUpdateLinear & data, SymbolsExpr const& symbolsExpr, hana::int_<2> ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto mesh = this->mesh();
    auto XhV = this->spaceMagneticPotential();
    auto const& u = this->fieldMagneticPotential();

    auto bilinearForm_PatternCoupled = form2( _test=XhV,_trial=XhV,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=this->rowStartInMatrix() ,
                                              _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=XhV, _vector=F,
                               _rowstart=this->rowStartInVector() );

    auto bcs = this->modelProperties().boundaryConditions2().boundaryConditions("magnetic-potential", "Dirichlet");
    for( auto const& bcpair : bcs )
    {
        auto bc = bcpair.second;
        auto markers = bc.markers();
        auto const& magneticPermeability = M_maxwellProperties->magneticPermeability( bc.material() );
        auto mu = expr( magneticPermeability.expr(), symbolsExpr);
        // Feel::cout << "Dirichlet on " << bcpair.first
        //            << " with " << markers.size() << " markers" << std::endl;
        auto Ad = bc.expr();

        bilinearForm_PatternCoupled += integrate(_range=markedfaces(M_mesh, markers),
                       _expr=1e5/hFace()*inner(cross(id(u),N()),cross(idt(u),N()))/mu );
        bilinearForm_PatternCoupled += integrate(_range=markedfaces(M_mesh, markers),
                       _expr=inner(curl(u),cross(idt(u),N()))/mu );
        bilinearForm_PatternCoupled += integrate(_range=markedfaces(M_mesh, markers),
                       _expr=inner(curlt(u),cross(id(u),N()))/mu );
        myLinearForm += integrate(_range=markedfaces(M_mesh, markers), _expr=inner(curl(u),Ad)/mu );
        myLinearForm += integrate(_range=markedfaces(M_mesh, markers),
                       _expr=1e5/hFace()*inner(cross(id(u),N()),Ad)/mu );
    }
}


} // namespace FeelModels
} // namespace Feel

#endif
