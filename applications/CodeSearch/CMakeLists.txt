###
### SharedObjects
###

# compile all the objects
file(GLOB files "${PROJECT_SOURCE_DIR}/applications/CodeSearch/sharedLibraries/source/*.cc")
foreach(file ${files})
    get_filename_component(fileName "${file}" NAME_WE)
    add_library(${fileName} MODULE ${file})
    target_link_libraries(${fileName} pdb-shared-common)
    target_link_libraries(${fileName} ${GSL_LIBRARIES})
    add_dependencies(shared-libraries ${fileName})
endforeach()


# create the target
add_pdb_application(CodeSearchLoadData)
add_pdb_application(CodeSearch)
add_pdb_application(CountDB)

add_dependencies(TopProgram SearchProgramData)


add_dependencies(CodeSearchLoadData ProgramResultWriter)
add_dependencies(CodeSearchLoadData SearchProgramData)
add_dependencies(CodeSearchLoadData TopProgram)

add_dependencies(CodeSearch ProgramResultWriter)
add_dependencies(CodeSearch SearchProgramData)
add_dependencies(CodeSearch TopProgram)
add_dependencies(CodeSearch CodeSearchLoadData)

# add a build dependency to build-tests target
add_dependencies(build-ml-tests CodeSearchLoadData)
add_dependencies(build-ml-tests CodeSearch)
