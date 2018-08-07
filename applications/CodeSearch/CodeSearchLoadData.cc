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
#ifndef CODE_SEARCH_LOAD_DATA
#define CODE_SEARCH_LOAD_DATA

// By Tania, September 2017
// Gaussian Mixture Model (based on EM) - Load DATA!!

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

#define BLOCKSIZE DEFAULT_NET_PAGE_SIZE

#define COUT std::cout
//#define COUT term

using namespace pdb;
int main(int argc, char *argv[]) {
  bool printResult = true;
  bool clusterMode = false;

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

  COUT << "Usage: #printResult[Y/N] #clusterMode[Y/N] blocksSize[MB] #managerIp "
          "#randomData[Y/N] #addData[Y/N] "
          "#niter, #clusters #nDatapoints #nDimensions "
          "#pathToInputFile(randomData == N)"
          "./bin/CodeSearchLoadData Y N 256 localhost Y"
       << std::endl;
  if (argc > 1) {
    if (strcmp(argv[1], "N") == 0) {
      printResult = false;
      COUT << "You successfully disabled printing result." << std::endl;
    } else {
      printResult = true;
      COUT << "Will print result." << std::endl;
    }
  } else {
    COUT << "Will print result. If you don't want to print result, you can add "
            "N as the first "
            "parameter to disable result printing."
         << std::endl;
  }

  if (argc > 2) {
    if (strcmp(argv[2], "Y") == 0) {
      clusterMode = true;
      COUT << "You successfully set the test to run on cluster." << std::endl;
    } else {
      clusterMode = false;
    }
  } else {
    COUT << "Will run on local node. If you want to run on cluster, you can "
            "add any character "
            "as the second parameter to run on the cluster configured by "
            "$PDB_HOME/conf/serverlist."
         << std::endl;
  }

  int blocksize = 256; // by default we add 64MB data
  if (argc > 3) {
    blocksize = atoi(argv[3]);
  }
  // numOfMb = 10; //Force it to be 64 by now.

  COUT << "To add data with size: " << blocksize << "MB" << std::endl;

  std::string managerIp = "localhost";
  if (argc > 4) {
    managerIp = argv[4];
  }
  COUT << "Manager IP Address is " << managerIp << std::endl;


  bool whetherToAddData = true;
  if (argc > 5) {
    if (strcmp(argv[5], "N") == 0) {
      whetherToAddData = false;
    }
  }

  COUT << blue << std::endl;
  COUT << "*****************************************" << std::endl;
  COUT << "Code Search starts : " << std::endl;
  COUT << "*****************************************" << std::endl;
  COUT << reset << std::endl;

  COUT << "The Code Search paramers are: " << std::endl;
  COUT << std::endl;

  int iter = 1;
  int k = 1;
  int dim = 3;
  int numData = 10;
  double convergenceTol = 0.001; // Convergence threshold


  if (argc > 6) {
    numData = std::stoi(argv[6]);
  }
  COUT << "The number of data points: " << numData << std::endl;

  if (argc >7) {
    dim = std::stoi(argv[7]);
  }
  COUT << "The dimension of each data point: " << dim << std::endl;

  std::string fileName = "/home/ubuntu/code_search_data.txt";
  if (argc > 8) {
    fileName = argv[8];
  }

  COUT << "Input file: " << fileName << std::endl;
  COUT << std::endl;

  COUT << std::endl;

  //***********************************************************************************
  //**** LOAD DATA
  //********************************************************************
  //***********************************************************************************

  PDBClient pdbClient(8108, managerIp);

  pdbClient.registerType("libraries/libSearchProgramData.so");
  pdbClient.registerType("libraries/libTopProgram.so");
  pdbClient.registerType("libraries/libProgramResultWriter.so");
  string errMsg;

  //    srand(time(0));
  // For the random number generator
  std::random_device rd;
  std::mt19937 randomGen(rd());

  //***********************************************************************************
  //****READ INPUT DATA
  //***************************************************************
  //***********************************************************************************

  // Step 1. Create Database and Set
  // now, register a type for user data
  // TODO: once sharedLibrary is supported, add this line back!!!

  if (whetherToAddData == true) {
    // now, create a new database
    pdbClient.createDatabase("code_search_db04");

    // now, create a new set in that database
    pdbClient.createSet<SearchProgramData>("code_search_db04", "code_search_input_set");
  }

  // Step 2. Add data

  auto begin = std::chrono::high_resolution_clock::now();

  if (whetherToAddData == true) {
      numData = 0;
      std::ifstream inFile(fileName.c_str());
      std::string line;
      bool rollback = false;
      bool end = false;

      numData = 0;
      while (!end) {
        pdb::makeObjectAllocatorBlock(blocksize * 1024 * 1024, true);
        pdb::Handle<pdb::Vector<pdb::Handle<SearchProgramData>>> storeMe =
            pdb::makeObject<pdb::Vector<pdb::Handle<SearchProgramData>>>();
        try {

          while (1) {
                if (!rollback) {
                    //      std::istringstream iss(line);
                    if (!std::getline(inFile, line)) {
                      end = true;
                      break;
                    }
                    else
                    {
                          pdb::Handle<SearchProgramData> myData =
                              pdb::makeObject<SearchProgramData>(dim);
                          std::stringstream lineStream(line);
                          double value;
                          int index = 0;
                          while (lineStream >> value) {
                              myData->setDouble(index, value);
                              index++;
                              if (index == dim){
                                break;
                              }
                          }
                          //Read ProbY
                          while (lineStream >> value)
                          {
                            myData->setProbY(value);
                            break;
                          }
                          // Read the Program from here
                          std::string temp;
                          std::string prog = "";
                          while (lineStream >> temp){
                            prog = prog + temp;
                          }
                          myData->setProg(prog);

                          storeMe->push_back(myData);
                          // myData->print();
                    }
                }
                else
                {
                  rollback = false;
                  pdb::Handle<SearchProgramData> myData =
                      pdb::makeObject<SearchProgramData>(dim);
                  std::stringstream lineStream(line);
                  double value;
                  int index = 0;
                  while (lineStream >> value) {
                    //(*myData)[index] = value;
                    myData->setDouble(index, value);
                    index++;
                    if (index == dim){
                      break;
                    }
                  }
                  //Read ProbY
                  while (lineStream >> value){
                    myData->setProbY(value);
                    break;
                  }
                  // Read the Program from here
                  std::string temp;
                  std::string prog = "";
                  while (lineStream >> temp){
                    prog = prog + temp;
                  }
                  myData->setProg(prog);

                  storeMe->push_back(myData);
            }
          }

          end = true;

          // send the rest of data at the end, it can happen that the exception
          // never
          // happens.
          pdbClient.sendData<SearchProgramData>(
                  std::pair<std::string, std::string>("code_search_input_set",
                                                      "code_search_db04"),
                  storeMe);

          numData += storeMe->size();
          COUT << "Added " << storeMe->size() << " Total: " << numData
               << std::endl;

          pdbClient.flushData();
        } catch (pdb::NotEnoughSpace &n) {
          pdbClient.sendData<SearchProgramData>(
                  std::pair<std::string, std::string>("code_search_input_set",
                                                      "code_search_db04"),
                  storeMe);

          numData += storeMe->size();
          COUT << "Added " << storeMe->size() << " Total: " << numData
               << std::endl;

          rollback = true;
        }
        PDB_COUT << blocksize << "MB data sent to dispatcher server~~"
                 << std::endl;

      } // while not end
      inFile.close();


  } // End if - whetherToAddData = true


  //pdbClient.removeSet("code_search_db04", "result");

  auto end = std::chrono::high_resolution_clock::now();

  COUT << "Time loading data: "
       << std::chrono::duration_cast<std::chrono::duration<float>>(end - begin)
              .count()
       << " secs." << std::endl;




   // this is the object allocation block where all of this stuff will reside
   // pdb::makeObjectAllocatorBlock(blocksize * 1024 * 1024, true);

   pdbClient.createSet<TopKQueue<double, SearchProgramData>>("code_search_db04", "result");
   Handle<Vector<double>> myQuery = makeObject<Vector<double>>();


   int i = 0;
   while ( i < dim) {
       double val = 0.0;
       myQuery->push_back(val);
       i+=1;
   }

   Handle<Computation> myScanSet =
    makeObject<ScanUserSet<SearchProgramData>>("code_search_db04", "code_search_input_set");
   Handle<Computation> myTopK = makeObject<TopProgram>(k, *myQuery);
   myTopK->setInput(myScanSet);

   Handle<Computation> myWriter = makeObject<ProgramResultWriter>("code_search_db04", "result");
   myWriter->setInput(myTopK);

   std::cout << "Ready to start computations" << std::endl;
   begin = std::chrono::high_resolution_clock::now();
   pdbClient.executeComputations(myWriter);
   end = std::chrono::high_resolution_clock::now();

   std::cout << "The query is executed successfully!" << std::endl;
   float timeDifference =
       (float(std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count())) /
       (float)1000000000;
   std::cout << "#TimeDuration: " << timeDifference << " Second " << std::endl;



}




#endif
