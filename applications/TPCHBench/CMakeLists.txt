###
### SharedObjects
###

# compile all the objects
file(GLOB files "${PROJECT_SOURCE_DIR}/applications/TPCHBench/sharedLibraries/source/*.cc")
foreach(file ${files})
    get_filename_component(fileName "${file}" NAME_WE)
    add_library(${fileName} MODULE ${file})
    target_link_libraries(${fileName} pdb-shared-common)
    target_link_libraries(${fileName} ${GSL_LIBRARIES})
    add_dependencies(shared-libraries ${fileName})
endforeach()


###
### tpchJaccard
###

# create the target
add_pdb_application(tpchJaccard)

# add a build dependency to build-tpch-tests target
add_dependencies(build-tpch-tests tpchJaccard)

# add build dependencies to shared libraries it uses
add_dependencies(tpchJaccard Customer)
add_dependencies(tpchJaccard LineItem)
add_dependencies(tpchJaccard Order)
add_dependencies(tpchJaccard Part)
add_dependencies(tpchJaccard Supplier)
add_dependencies(tpchJaccard ScanCustomerSet)
add_dependencies(tpchJaccard TopJaccard)
add_dependencies(tpchJaccard AllParts)
add_dependencies(tpchJaccard JaccardResultWriter)

# add build dependencies to shared libraries it uses
#add_dependencies(tpchDataGeneratorAll Part)
#add_dependencies(tpchDataGeneratorAll Supplier)
#add_dependencies(tpchDataGeneratorAll LineItem)
#add_dependencies(tpchDataGeneratorAll Order)
#add_dependencies(tpchDataGeneratorAll Customer)
#add_dependencies(tpchDataGeneratorAll SumResultWriteSet)
#add_dependencies(tpchDataGeneratorAll CustomerWriteSet)
#add_dependencies(tpchDataGeneratorAll ScanCustomerSet)
#add_dependencies(tpchDataGeneratorAll CustomerMultiSelection)
#add_dependencies(tpchDataGeneratorAll CustomerSupplierPartGroupBy)
#add_dependencies(tpchDataGeneratorAll SupplierData)
#add_dependencies(tpchDataGeneratorAll CustomerSupplierPartFlat)
#add_dependencies(tpchDataGeneratorAll CountAggregation)
#add_dependencies(tpchDataGeneratorAll CountCustomer)

###
### tcphDataGeneratorNew
###

# create the target
add_pdb_application(tpchDataGenerator)

# add a build dependency to build-tpch-tests target
add_dependencies(build-tpch-tests tpchDataGenerator)

###
### tpchGetCustomerCount
###

# create the target
add_pdb_application(tpchGetCustomerCount)

# add a build dependency to build-tpch-tests target
add_dependencies(build-tpch-tests tpchGetCustomerCount)

###
### tpchQuery
###

# create the target
add_pdb_application(tpchQuery)

# add a build dependency to build-tpch-tests target
add_dependencies(build-tpch-tests tpchQuery)

###
### tpchRegisterAndCreateSets
###

# create the target
add_pdb_application(tpchRegisterAndCreateSets)

# add a build dependency to build-tpch-tests target
add_dependencies(build-tpch-tests tpchRegisterAndCreateSets)

# add build dependencies to shared libraries it uses
add_dependencies(tpchRegisterAndCreateSets Part)
add_dependencies(tpchRegisterAndCreateSets Supplier)
add_dependencies(tpchRegisterAndCreateSets LineItem)
add_dependencies(tpchRegisterAndCreateSets Order)
add_dependencies(tpchRegisterAndCreateSets Customer)
add_dependencies(tpchRegisterAndCreateSets SumResultWriteSet)
add_dependencies(tpchRegisterAndCreateSets CustomerWriteSet)
add_dependencies(tpchRegisterAndCreateSets ScanCustomerSet)
add_dependencies(tpchRegisterAndCreateSets CustomerMultiSelection)
add_dependencies(tpchRegisterAndCreateSets CustomerSupplierPartGroupBy)
add_dependencies(tpchRegisterAndCreateSets SupplierInfo)
add_dependencies(tpchRegisterAndCreateSets CustomerSupplierPartFlat)
add_dependencies(tpchRegisterAndCreateSets CountAggregation)
add_dependencies(tpchRegisterAndCreateSets CountCustomer)
add_dependencies(tpchRegisterAndCreateSets SupplierInfoWriteSet)