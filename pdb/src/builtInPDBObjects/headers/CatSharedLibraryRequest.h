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

#ifndef CAT_SHARED_LIB_REQ_H
#define CAT_SHARED_LIB_REQ_H

#include "Object.h"
#include "Handle.h"

// PRELOAD %CatSharedLibraryRequest%

namespace pdb {

// encapsulates a request to obtain a shared library from the catalog
class CatSharedLibraryRequest : public Object {

public:
    CatSharedLibraryRequest(){};
    CatSharedLibraryRequest(int16_t searchForMe) : findMe(searchForMe) {}
    ~CatSharedLibraryRequest(){};

    int16_t getTypeID() {
        return findMe;
    }

    ENABLE_DEEP_COPY

private:
    int16_t findMe;
};
}

#endif
