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
#include "ApplyBase.h"

using std::invalid_argument;

namespace pdb_detail {
ApplyBase::ApplyBase(const string& executorId,
                     const string& functionId,
                     const string& outputTableId,
                     const string& outputColumnId,
                     TableColumns inputColumns,
                     shared_ptr<vector<TableColumn>> columnsToCopyToOutputTable,
                     InstructionType type)

    : Instruction(type),
      executorId(executorId),
      functionId(functionId),
      outputTableId(outputTableId),
      outputColumnId(outputColumnId),
      inputColumns(inputColumns),
      columnsToCopyToOutputTable(columnsToCopyToOutputTable) {
    if (columnsToCopyToOutputTable == nullptr)
        throw invalid_argument("columnsToCopyToOutputTable may not be null");
}

uint ApplyBase::getOutputTableColumnCount() {
    return columnsToCopyToOutputTable->size() + 1;
}
}