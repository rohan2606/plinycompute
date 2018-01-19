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
//
// Created by barnett on 10/25/16.
//

#ifndef PDB_TCAPINTERMEDIARYREPTESTS_TCAPIRTESTSRUNNER_H
#define PDB_TCAPINTERMEDIARYREPTESTS_TCAPIRTESTSRUNNER_H

#include "qunit.h"

using QUnit::UnitTest;

namespace pdb_tests {
void runBuildTcapIrTests(UnitTest& qunit);
}


#endif  // PDB_TCAPIRTESTSRUNNER_H