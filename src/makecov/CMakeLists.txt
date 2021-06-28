# Add executable that is built from the listed source files.
# The extensions are automatically found.
add_executable( makecov 
	makecov.cpp 
)


file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/${LOCALSTATEDIR}" COTHERSTATEDIRNAT)
string(REGEX REPLACE "([\\])" "\\\\\\1" COTHERSTATEDIRNAT "${COTHERSTATEDIRNAT}")


## ONE executable
## Specify include directories
target_include_directories( makecov PUBLIC 
	${CMAKE_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions( makecov PUBLIC
	COTHER_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}
	COTHERSTATEDIR=${COTHERSTATEDIRNAT}
)

# Link the executable to the libraries 
target_link_libraries( makecov PUBLIC pro HDP seg lib extsp )


## install targets
install(TARGETS 
	makecov
	DESTINATION ${LOCALBINDIR}
)

