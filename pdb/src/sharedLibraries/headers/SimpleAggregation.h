/*****************************************************************************
 *                                                                           *
 *  Copyright 2018 Rice University                                           *
 *                                                                           *
 *  Licensed under the Apache License, Version 2.0 (the "License");          *
 *  you may not use this file except in compliance with the License.         *
 *  You may obtain a copy of the License at                                  *
 *                                                                           *
 *      http://www.apache.org/licenses/LICENSE-2.0                           *
 *                                                                           *
 *  Unless required by applicable law or agreed to in writing, software      *
 *  distributed under the License is distributed on an "AS IS" BASIS,        *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 *  See the License for the specific language governing permissions and      *
 *  limitations under the License.                                           *
 *                                                                           *
 *****************************************************************************/
#ifndef SIMPLE_AGG_H
#define SIMPLE_AGG_H

// by Jia, Mar 2017

#include "AggregateComp.h"
#include "DepartmentTotal.h"
#include "Employee.h"
#include "LambdaCreationFunctions.h"


using namespace pdb;

class SimpleAggregation : public AggregateComp<DepartmentTotal, Employee, String, double> {

public:
    ENABLE_DEEP_COPY

    SimpleAggregation() {}

    // the below constructor is NOT REQUIRED
    // user can also set output later by invoking the setOutput (std :: string dbName, std :: string
    // setName)  method
    SimpleAggregation(std::string dbName, std::string setName) {
        this->setOutput(dbName, setName);
    }


    // the key type must have == and size_t hash () defined
    Lambda<String> getKeyProjection(Handle<Employee> aggMe) override {
        return makeLambdaFromMember(aggMe, department);
    }

    // the value type must have + defined
    Lambda<double> getValueProjection(Handle<Employee> aggMe) override {
        return makeLambdaFromMethod(aggMe, getSalary);
    }
};


#endif
