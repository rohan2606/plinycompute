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
#ifndef PIPELINE_STAGE_H
#define PIPELINE_STAGE_H

//by Jia, Mar 2017

#include "DataTypes.h"
#include "PDBLogger.h"
#include "Configuration.h"
#include "SharedMem.h"
#include "TupleSetJobStage.h"
#include "HermesExecutionServer.h"
#include <vector>
#include <memory>
#include <unordered_map>

namespace pdb {

class PipelineStage;
typedef std::shared_ptr<PipelineStage> PipelineStagePtr;


// this class encapsulates the pipeline stage
class PipelineStage {

private:

    // Job stage 
    Handle<TupleSetJobStage> jobStage;

    // batch size
    size_t batchSize;

    // number of threads
    int numThreads;

    // node id
    NodeID nodeId;

    // logger
    PDBLoggerPtr logger;

    // configuration
    ConfigurationPtr conf;

    // shared memory
    SharedMemPtr shm;

    // operator id
    OperatorID id;

public:

    //destructor
    ~PipelineStage();

    //constructor
    PipelineStage(Handle<TupleSetJobStage> stage, SharedMemPtr shm, PDBLoggerPtr logger, ConfigurationPtr conf, NodeID nodeId, size_t batchSize, int numThreads);

    //return the root job stage corresponding to the pipeline 
    Handle<TupleSetJobStage> & getJobStage(); 

    //return the number of threads that are required to run the pipeline network
    int getNumThreads();

    //run a Map pipeline with UserSet source and UserSet sink 
    void runMapPipeline (HermesExecutionServer * server);

    //run a Map pipeline with UserSet source and ShuffleSet sink
    void runMapPipelineWithShuffleSink  (HermesExecutionServer * server);

    //run a Reduce pipeline with ShuffleSet source and UserSet sink
    void runReducePipeline  (HermesExecutionServer * server);

    //run a Reduce pipeline with ShuffleSet source and ShuffleSet sink
    void runReducePipelineWithShuffleSink  (HermesExecutionServer * server);

};





}





#endif