#include <ginac/parse_context.h>

namespace Feel
{
/**
 * Default prototype table for feelpp
 *
 * It supports all defined GiNaC functions and "pow", "sqrt", and "power".
 */
extern const GiNaC::prototype_table& get_default_reader();
}