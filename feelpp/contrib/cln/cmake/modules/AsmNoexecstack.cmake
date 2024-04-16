# Check whether the stack can be marked nonexecutable by passing
# an option to the C compiler when acting on .s files. 
# 
# ASM_NOEXECSTACK_FLAG -- compiler option(s) for marking the stack
#                         nonexecutable

# CC='ccache gcc' => "${CMAKE_C_COMPILER_ARG1}" == " gcc"
# (notice the leading whitespace). Grrr!
string(STRIP "${CMAKE_C_COMPILER_ARG1}" _c_compiler_arg1)

set(_conftest_c "${CMAKE_CURRENT_BINARY_DIR}/conftest.c")
set(_conftest_s "${CMAKE_CURRENT_BINARY_DIR}/conftest.s")
set(_conftest_o "${CMAKE_CURRENT_BINARY_DIR}/conftest.o")
set(_need_noexecstack)
set(_cc_ret)

file(WRITE ${_conftest_c} "void foo() { }")

execute_process(
	COMMAND ${CMAKE_C_COMPILER} ${_c_compiler_arg1} -S ${_conftest_c} -o ${_conftest_s}
	RESULT_VARIABLE _cc_ret
)

if ("${_cc_ret}" EQUAL "0")
	file(STRINGS ${_conftest_s} _need_noexecstack REGEX "\\.note\\.GNU-stack")
endif()

if (_need_noexecstack)
	execute_process(COMMAND ${CMAKE_C_COMPILER} ${_c_compiler_arg1} -Wa,--noexecstack -c ${_conftest_s}
			RESULT_VARIABLE _cc_ret)
	if ("${_cc_ret}" EQUAL "0")
		set(ASM_NOEXECSTACK_FLAG "-Wa,--noexecstack")
	endif()
endif()
file(REMOVE ${_conftest_o} ${_conftest_s} ${_conftest_c})

