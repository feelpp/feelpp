# Check if symbols are prefixed by an underscore in assembly language.
set(ASM_UNDERSCORE)

set(_conftest_c ${CMAKE_CURRENT_BINARY_DIR}/conftest.c)
set(_conftest_s ${CMAKE_CURRENT_BINARY_DIR}/conftest.s)
set(_cc_ret)

file(WRITE ${_conftest_c} "int foo() { return 0; }")

# CC='ccache gcc' => "${CMAKE_C_COMPILER_ARG1}" == " gcc"
# (notice the leading whitespace). Grrr!
string(STRIP "${CMAKE_C_COMPILER_ARG1}" _c_compiler_arg1)
execute_process(
	COMMAND ${CMAKE_C_COMPILER} ${_c_compiler_arg1} -S ${_conftest_c} -o ${_conftest_s}
	RESULT_VARIABLE _cc_ret
)
if ("${_cc_ret}" EQUAL "0")
	file(STRINGS ${_conftest_s} _asm_underscore REGEX "_foo")
	if (_asm_underscore)
		set(ASM_UNDERSCORE true)
	endif()
endif()
file(REMOVE ${_conftest_s} ${_conftest_c})

