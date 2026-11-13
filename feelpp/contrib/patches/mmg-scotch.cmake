# Backport feelpp/mmg commit 3e4f752e585964e4ab38d38c6679d5692d18a0ac
# (MmgTools/mmg#329). Scotch's private symbol prefix depends on its build options;
# sorting locally avoids that dependency (MmgTools/mmg#328).
function(feelpp_patch_mmg_scotch source_dir)
    set(source "${source_dir}/src/common/librnbg.c")
    file(READ "${source}" original)
    set(old_call "_SCOTCHintSort2asc1(sortPartTb, vertNbr);")
    string(FIND "${original}" "${old_call}" call_position)
    if(NOT call_position EQUAL -1)
        set(comparator [=[/**
 * \param lhs pointer to the first {box, vertex} pair.
 * \param rhs pointer to the second {box, vertex} pair.
 * \return a negative value, zero or a positive value if the first box index
 * is respectively lower than, equal to or greater than the second one.
 *
 * Compare the box index of two Scotch partition pairs.
 *
 */
static int
MMG5_compareScotchPairsFirst(const void *lhs, const void *rhs) {
  const SCOTCH_Num *a = (const SCOTCH_Num *)lhs;
  const SCOTCH_Num *b = (const SCOTCH_Num *)rhs;

  return (a[0] > b[0]) - (a[0] < b[0]);
}
]=])
        string(REPLACE "#include \"librnbg.h\""
                       "#include \"librnbg.h\"\n\n${comparator}" patched "${original}")
        string(REPLACE "${old_call}"
                       "qsort(sortPartTb, (size_t)vertNbr, 2 * sizeof(*sortPartTb),\n        MMG5_compareScotchPairsFirst);"
                       patched "${patched}")
        file(WRITE "${source}" "${patched}")
        message(STATUS "[feelpp] Removed bundled MMG's dependency on a private Scotch sort symbol")
    endif()

    set(header "${source_dir}/src/common/librnbg.h")
    file(READ "${header}" original)
    string(REPLACE "int    _SCOTCHintSort2asc1(SCOTCH_Num * sortPartTb, int vertNbr);\n"
                   "" patched "${original}")
    if(NOT original STREQUAL patched)
        file(WRITE "${header}" "${patched}")
    endif()
endfunction()
