# Add executable that is built from the listed source files.
# The extensions are automatically found.
add_executable( cother 
	cother.cpp 
)


file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/${LOCALSTATEDIR}" COTHERSTATEDIRNAT)
string(REGEX REPLACE "([\\])" "\\\\\\1" COTHERSTATEDIRNAT "${COTHERSTATEDIRNAT}")


## ONE executable
## Specify include directories
target_include_directories( cother PUBLIC 
	${CMAKE_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions( cother PUBLIC
	COTHER_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}
	COTHERSTATEDIR=${COTHERSTATEDIRNAT}
)

set_target_properties( cother PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties( mycu PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)

# Link the executable to the libraries 
target_link_libraries( cother PUBLIC mycu pro HDP seg lib extsp Threads::Threads)


## install targets
install(TARGETS 
	cother 
	DESTINATION ${LOCALBINDIR}
)
#install(FILES 
#	cother.sh 
#	DESTINATION ${LOCALBINDIR} 
#	PERMISSIONS OWNER_EXECUTE OWNER_READ GROUP_READ WORLD_READ
#)

