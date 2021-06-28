# Add executable that is built from the listed source files.
# The extensions are automatically found.
add_executable( makepro 
	makepro.cpp 
)
add_executable( select 
	select.cpp 
)


file(TO_NATIVE_PATH "${CMAKE_INSTALL_PREFIX}/${LOCALSTATEDIR}" COTHERSTATEDIRNAT)
string(REGEX REPLACE "([\\])" "\\\\\\1" COTHERSTATEDIRNAT "${COTHERSTATEDIRNAT}")


## ONE executable
## Specify include directories
target_include_directories( makepro PUBLIC 
	${CMAKE_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions( makepro PUBLIC
	COTHER_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}
	COTHERSTATEDIR=${COTHERSTATEDIRNAT}
)

# Link the executable to the libraries 
target_link_libraries( makepro PUBLIC pro HDP seg lib extsp )


## ANOTHER executable
##
## Specify include directories
target_include_directories( select PUBLIC 
	${CMAKE_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions( select PUBLIC
	COTHER_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}
	COTHERSTATEDIR=${COTHERSTATEDIRNAT}
)

# Link the executable to the libraries 
target_link_libraries( select PUBLIC pro HDP seg lib extsp )

set(MAKEPROSH "makepro.sh")
if(WIN32)
  set(MAKEPROSH "makepro.cmd")
endif()

## install targets
install(TARGETS 
	makepro 
	select 
	DESTINATION ${LOCALBINDIR}
)
install(FILES 
	ssp2.pl 
	inssp2.pl 
	${MAKEPROSH} 
	DESTINATION ${LOCALBINDIR} 
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

