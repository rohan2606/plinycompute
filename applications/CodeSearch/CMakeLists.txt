
# create the target
add_pdb_application(CodeSearchLoadData)

# add a build dependency to build-tests target
add_dependencies(build-ml-tests CodeSearchLoadData)
