# NLopt's build-time cmake -P scripts do not inherit the parent's policy settings.
# CMake 4 requires these scripts to declare compatibility with at least CMake 3.5.
function(feelpp_patch_nlopt_cmake source_dir)
    foreach(generator IN ITEMS generate-cpp generate-fortran)
        set(script "${source_dir}/cmake/${generator}.cmake")
        file(READ "${script}" original)
        string(REPLACE "cmake_minimum_required (VERSION 3.2)"
                       "cmake_minimum_required (VERSION 3.5)" patched "${original}")
        if(NOT original STREQUAL patched)
            file(WRITE "${script}" "${patched}")
            message(STATUS "[feelpp] Applied bundled NLopt ${generator} compatibility fix for CMake 4")
        endif()
    endforeach()
endfunction()
