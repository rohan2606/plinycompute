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

#ifndef WRITE_USER_SET_H
#define WRITE_USER_SET_H

#include "Computation.h"
#include "VectorSink.h"
namespace pdb {

template <class OutputClass>
class WriteUserSet : public Computation {

public:

        ENABLE_DEEP_COPY

        ComputeSinkPtr getComputeSink (TupleSpec &consumeMe, TupleSpec &projection, ComputePlan &plan) override {
        
             return std :: make_shared<VectorSink <OutputClass>> (consumeMe, projection);

        }



	std :: string getComputationType () override {
		return std :: string ("WriteUserSet");
	}
};

}

#endif
