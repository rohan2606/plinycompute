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
// Created by barnett on 8/26/16.
//

#ifndef QUERYALGO_H
#define QUERYALGO_H


#include "Object.h"

namespace pdb {
class QueryBase;

class QueryAlgo : public Object {
public:
    virtual void forSelection(QueryBase& selection) = 0;

    virtual void forSet(QueryBase& set) = 0;

    virtual void forQueryOutput(QueryBase& output) = 0;
};
}

#endif  // QUERYALGO_H