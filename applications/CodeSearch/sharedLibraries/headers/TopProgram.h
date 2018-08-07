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
#ifndef TOP_PROGRAM_H
#define TOP_PROGRAM_H

#include "Object.h"
#include "PDBVector.h"
#include "PDBString.h"
#include "Handle.h"
#include "TopKComp.h"
#include "SearchProgramData.h"
#include "AllParts.h"

namespace pdb {

// For each customer, constructs a list of unique part IDs, and then checks to compute the
// Jaccard similarity to the query part list.  It then returns the k part lists that are most
// simiar to the query part list
class TopProgram : public TopKComp<SearchProgramData, double , Handle<SearchProgramData>> {

    Vector<double> query;
    int k = 1;

public:
    ENABLE_DEEP_COPY

    TopProgram() {}

    TopProgram(int kIn, Vector<double> query)
        : k(kIn), query(query) {}

    Lambda<TopKQueue<double, Handle<SearchProgramData>>> getValueProjection(
        Handle<SearchProgramData> Prog) override {
        return makeLambda(Prog, [&](Handle<SearchProgramData>& Prog) {
            //Handle<SearchProgramData> myGuy = makeObject<SearchProgramData>();
            return TopKQueue<double, Handle<SearchProgramData>>(
                k, Prog->getSquaredDistance(query) , Prog);
        });
    }
};
}
#endif
