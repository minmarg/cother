# Add executable that is built from the listed source files.
# The extensions are automatically found.
add_executable( adddist 
	adddist.cpp 
)


file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/${LOCALSTATEDIR}" COTHERSTATEDIRNAT)
string(REGEX REPLACE "([\\])" "\\\\\\1" COTHERSTATEDIRNAT "${COTHERSTATEDIRNAT}")


## ONE executable
## Specify include directories
target_include_directories( adddist PUBLIC 
	${CMAKE_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions( adddist PUBLIC
	COTHER_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}
	COTHERSTATEDIR=${COTHERSTATEDIRNAT}
)

# Link the executable to the libraries 
target_link_libraries( adddist PUBLIC pro HDP lib extsp )


## install targets
install(TARGETS 
	adddist 
	DESTINATION ${LOCALBINDIR}
)
install(FILES 
	batchadddist.py 
	batchadddist_ray.py 
	DESTINATION ${LOCALBINDIR} 
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

