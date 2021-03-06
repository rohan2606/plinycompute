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
#ifndef COMBINER_PROCESSOR_H
#define COMBINER_PROCESSOR_H


#include "UseTemporaryAllocationBlock.h"
#include "InterfaceFunctions.h"
#include "PDBMap.h"
#include "AggregationMap.h"
#include "PDBVector.h"
#include "Handle.h"
#include "SimpleSingleTableQueryProcessor.h"
#include "DataTypes.h"
#include <vector>

namespace pdb {

//this class is used for combining key-value pairs as pre-aggregation

template <class KeyType, class ValueType>
class CombinerProcessor : public SimpleSingleTableQueryProcessor {

public:
    ~CombinerProcessor() {
        blockPtr = nullptr;
        begin = nullptr;
        end = nullptr;
    };
    CombinerProcessor(std::vector<HashPartitionID>& partitions);
    void initialize() override;
    void loadInputPage(void* pageToProcess) override;
    void loadOutputPage(void* pageToWriteTo, size_t numBytesInPage) override;
    bool fillNextOutputPage() override;
    void finalize() override;
    void clearOutputPage() override;
    void clearInputPage() override;
    void addNodePartition(HashPartitionID partitionId);

private:
    UseTemporaryAllocationBlockPtr blockPtr;
    Handle<Vector<Handle<Map<KeyType, ValueType>>>> inputData;
    Handle<Vector<Handle<AggregationMap<KeyType, ValueType>>>> outputData;
    bool finalized;
    int numNodePartitions;
    HashPartitionID curPartId;
    int curPartPos;
    Handle<Map<KeyType, ValueType>> curMap;
    Handle<AggregationMap<KeyType, ValueType>> curOutputMap;

    // the iterators for current map partition
    PDBMapIterator<KeyType, ValueType>* begin;
    PDBMapIterator<KeyType, ValueType>* end;

    // partitions on this node
    std::vector<HashPartitionID> nodePartitionIds;
    int count;
};
}


#include "CombinerProcessor.cc"


#endif
