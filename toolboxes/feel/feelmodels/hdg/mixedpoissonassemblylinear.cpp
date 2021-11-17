#include <feel/feelmodels/hdg/mixedpoisson.cpp>

namespace Feel {
namespace FeelModels {

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    auto A = std::dynamic_pointer_cast<condensed_matrix_t<value_type>>(data.matrix());
    auto F = std::dynamic_pointer_cast<condensed_vector_t<value_type>>(data.rhs());
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string scst = buildCstPart ? "(build cst part)" : "build non cst part)";
    this->log("MixedPoisson", "updateLinearPDE", "start"+scst);

    this->updateLinearPDE( data, this->modelContext() );
}

MIXEDPOISSON_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDPOISSON_CLASS_TEMPLATE_TYPE::updatePostPDE( DataUpdateLinear & data ) const
{
    this->updatePostPDE( data, this->modelContext() );
}

}
}
