# create the target
add_pdb_integration_test(TestMultiSelectionWithJoin)

# add a build dependency to build-tests target
add_dependencies(build-integration-tests TestMultiSelectionWithJoin)

# add build dependencies to shared libraries it uses
add_dependencies(TestMultiSelectionWithJoin OptimizedMethodJoin)
add_dependencies(TestMultiSelectionWithJoin ScanStringIntPairSet)
add_dependencies(TestMultiSelectionWithJoin StringSelectionOfStringIntPair)
add_dependencies(TestMultiSelectionWithJoin StringIntPairMultiSelection)
add_dependencies(TestMultiSelectionWithJoin WriteStringSet)