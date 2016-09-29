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
// Created by barnett on 9/28/16.
//

#ifndef PDB_QUERYINTERMEDIARYREP_MATERIALIZATIONMODEALGO_H
#define PDB_QUERYINTERMEDIARYREP_MATERIALIZATIONMODEALGO_H

#include "MaterializationMode.h"

namespace pdb_detail
{
    class MaterializationModeNone;

    class MaterializationModeNamedSet;

    class MaterializationModeAlgo
    {

    public:

        virtual void forNone(MaterializationModeNone &mode) = 0;

        virtual void forNamedSet(MaterializationModeNamedSet &mode) = 0;
    };
}

#endif //PDB_QUERYINTERMEDIARYREP_MATERIALIZATIONMODEALGO_H
