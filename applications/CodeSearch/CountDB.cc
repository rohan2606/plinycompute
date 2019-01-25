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
#ifndef CODE_SEARCH_COUNTDB
#define CODE_SEARCH_COUNTDB

// By Rohan, August 2018
// CodeSearch

#include "Lambda.h"
#include "PDBClient.h"
#include "PDBDebug.h"
#include "PDBString.h"
#include "Query.h"

#include "DataTypes.h"
#include "SearchProgramData.h"
#include "TopProgram.h"
#include "ProgramResultWriter.h"
#include "Set.h"

#include <chrono>
#include <cstddef>
#include <ctime>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <random>
#include <sys/stat.h>
#include <sys/types.h>
#include <thread>
#include <time.h>
#include <unistd.h>
#include <cstring>

#include <string>
#include <streambuf>

#include <vector>


// for convenience
#define BLOCKSIZE DEFAULT_NET_PAGE_SIZE

#define COUT std::cout
//#define COUT term

using namespace pdb;
int main(int argc, char *argv[]) {
  bool printResult = true;

  freopen("/dev/tty", "w", stdout);

  const std::string red("\033[0;31m");
  const std::string green("\033[1;32m");
  const std::string yellow("\033[1;33m");
  const std::string blue("\033[1;34m");
  const std::string cyan("\033[0;36m");
  const std::string magenta("\033[0;35m");
  const std::string reset("\033[0m");

  //***********************************************************************************
  //**** INPUT PARAMETERS
  //***************************************************************
  //***********************************************************************************


  COUT << magenta << std::endl;
  COUT << "Usage: ./bin/CountDB blocksSize[MB] #managerIp"
          "./bin/CountDB 8 localhost"  << std::endl;
  COUT << reset << std::endl;


  int blocksize = 64; // by default we add 64MB data
  if (argc > 1) {
    blocksize = atoi(argv[1]);
  }

  std::string managerIp = "localhost";
  if (argc > 2) {
    managerIp = argv[2];
  }




  //***********************************************************************************
  //*******************************  REGISTER TYPES ***********************************
  //***********************************************************************************


  PDBClient pdbClient(8108, managerIp);
  pdbClient.registerType("libraries/libSearchProgramData.so");
  pdbClient.registerType("libraries/libTopProgram.so");
  pdbClient.registerType("libraries/libProgramResultWriter.so");


   //***********************************************************************************
   //*********************************LOAD QUERY********************************
   //***********************************************************************************
   //***********************************************************************************
   // Use the JAVA Eclipse DOM driver to get the Query as a JSON and parse it similar to before.

   string errMsg;
   // Step 2. Load the QUERY
   // Load the Query and Get Top-k
   pdb::makeObjectAllocatorBlock(blocksize * 1024 * 1024, true);


   // now iterate through the result0
   SetIterator<SearchProgramData> inputData =
           pdbClient.getSetIterator<SearchProgramData>("code_search_db", "code_search_input_set");

   int i = 0;
   for (auto& a : inputData) {
  	   i += 1;
   }
   std::cout << "Got back " << i << " items from the top-k query.\n";

}

#endif
