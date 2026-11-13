# Clang 19 and newer reject the old ExprTk template disambiguator without an
# explicit argument list. These calls can deduce their template arguments.
function(feelpp_patch_openturns_exprtk source_dir)
    set(header "${source_dir}/lib/src/Base/Func/openturns/exprtk.hpp")
    file(READ "${header}" original)
    string(REPLACE "::ndb_t::template collect(" "::ndb_t::collect(" patched "${original}")
    if(NOT original STREQUAL patched)
        file(WRITE "${header}" "${patched}")
        message(STATUS "[mor] Applied bundled ExprTk compatibility fix for Clang 19+")
    endif()
endfunction()
