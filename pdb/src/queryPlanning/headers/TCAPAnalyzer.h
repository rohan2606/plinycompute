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
#ifndef PDB_TCAPANALYZERNEW_H
#define PDB_TCAPANALYZERNEW_H

#include <set>
#include "AggregationJobStage.h"
#include "AtomicComputationList.h"
#include "BroadcastJoinBuildHTJobStage.h"
#include "Computation.h"
#include "ComputePlan.h"
#include "Configuration.h"
#include "HashPartitionedJoinBuildHTJobStage.h"
#include "LogicalPlan.h"
#include "PDBLogger.h"
#include "Statistics.h"
#include "TupleSetJobStage.h"
#include "AbstractTCAPAnalyzerNodeFactory.h"

namespace pdb {

/**
 * This class takes in the computations and the TCAP string that are sent from the client and iteratively creates a
 * series of partial physical plans, that we can execute to get a result.
 *
 * This is accomplished by iteratively calling the method getNextStagesOptimized to generate a sequence of JobStages.
 * As an input to the getNextStagesOptimized we have to provide the storage statistics, so
 * it can determine the best starting source.
 *
 * The statistics should be iteratively updated after the execution of each sequence of JobStages, to reflect the
 * current state.
 *
 * There are four types of JobStages that can be generated
 * -- TupleSetJobStage: a pipeline
 * -- AggregationJobStage: This stage performs the aggregation on shuffled data
 * -- BroadcastJoinBuildHTJobStage: Builds a hash table for the broadcast join
 * -- HashPartitionedJoinBuildHTJobStage: Builds the hash table for the partitioned join
 *
 * The TCAP computations (AtomicComputation) are :
 * -- ApplyLambda - is a computation that applies a lambda to a tuple set
 * -- HashLeft - is a computation that applies a hash to a particular attribute in a tuple set (if left child of equals lambda)
 * -- HashRight - is a computation that applies a hash to a particular attribute in a tuple set (if right child of equals lambda)
 * -- HashOne -  is a computation that adds 1 to each tuple of a tuple set
 * -- Flatten - is a computation that flatten each tuple of a tuple set
 * -- ApplyFilter - is a computation that performs a filter over a tuple set
 * -- ApplyAgg -  is a computation that aggregates a tuple set
 * -- ScanSet - is a computation that produces a tuple set by scanning a set stored in the database
 * -- WriteSet - is a computation that writes out a tuple set
 * -- ApplyJoin - is a computation that performs the join of two sets
 */
class TCAPAnalyzer {
 public:

  /**
   * The constructor for the TCAPAnalyzer from a TCAP string and a list of computations associated with it
   * @param jobId the id of the current job, can any legal database name, generated by the QuerySchedulingServer
   * @param computations a list of computations that we want to execute
   * @param TCAPString the TCAP string extracted from these computations
   * @param logger an instance of the PDBLogger
   * @param conf the configuration of the nodes, used to fill the necessary information in the JobStages
   */
  TCAPAnalyzer(std::string &jobId,
                  PDBLoggerPtr logger,
                  ConfigurationPtr &conf,
                  std::string TCAPString,
                  pdb::Handle<pdb::Vector<pdb::Handle<pdb::Computation>>> computations);

  /**
     * Returns a sequence of job stages that, make up a partial physical plan. After the execution we gather the
     * statistics about the newly created sets and use them to generate the next partial plan.
     * @param physicalPlanToOutput a list where we want to put the sequence of job stages
     * @param interGlobalSets a list of intermediates sets that need to be created
     * @param stats the statistics about
     * @param jobStageId the id of the current job stage
     * @return true if we succeeded in creating the partial physical plan.
     */
  bool getNextStagesOptimized(std::vector<pdb::Handle<AbstractJobStage>> &physicalPlanToOutput,
                              std::vector<pdb::Handle<SetIdentifier>> &interGlobalSets,
                              StatisticsPtr &stats,
                              int &jobStageId);

  /**
   * Check if we still have some sources to process
   * @return true if we do, false otherwise
   */
  bool hasSources();

  /**
   * Returns true if the the provided source still has any consumers that we need to process
   * @param name source set name in the form of "databaseName:setName"
   * @return true if the the provided source has any consumers, false otherwise
   */
  bool hasConsumers(std::string &name);


  /**
   * Returns the best source node based on heuristics
   * @return the node
   */
  AbstractTCAPAnalyzerNodePtr getBestNode(StatisticsPtr &ptr);

 private:

  /**
   * This is the factor applied to the cost of the source if penalized
   */
  static constexpr double SOURCE_PENALIZE_FACTOR = 1000.00;

  /**
    * Hash map where the key is the name of the source set in the form of "databaseName:setName"
    * and the AbstractTCAPAnalyzerNodePtr associated with it.
    */
  std::map<std::string, AbstractTCAPAnalyzerNodePtr> sourceNodes;

  /**
   * The starting sources as AtomicComputations
   */
  std::vector<AtomicComputationPtr> sourcesComputations;

  /**
   * Penalized source sets in the form databaseName:setName
   */
  std::set<std::string> penalizedSets;

  /**
   * An instance of the PDBLogger
   */
  PDBLoggerPtr logger;

  /**
   * A configuration object for this node
   */
  ConfigurationPtr conf;

  /**
 * The ComputePlan generated from input computations and the input TCAP string
 */
  Handle<ComputePlan> computePlan;

  /**
   * Logical plan generated from the compute plan
   */
  LogicalPlanPtr logicalPlan;

  /**
   * The computation graph generated from the logical plan
   */
  AtomicComputationList computationGraph;

};

}


#endif //PDB_TCAPANALYZERNEW_H
