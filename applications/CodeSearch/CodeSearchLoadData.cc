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
#include <cstring>

#include "rapidjson/filereadstream.h"
#include "rapidjson/writer.h"
#include "rapidjson/document.h"
#include <vector>


// for convenience
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


  // read a JSON file
  // std::ifstream inp("/home/ubuntu/code_search_data.json");
  // document holds a json document retrieved from a http GET request
  // I did not include all of that in this example.  I am only showing
  // the part of iterating through a nested object and retrieving members.

  //***********************************************************************************
  //**** INPUT PARAMETERS
  //***************************************************************
  //***********************************************************************************

  COUT << "Usage: #printResult[Y/N] #clusterMode[Y/N] blocksSize[MB] #managerIp "
          "#addData[Y/N]" "#nDimensions" "#topKResults" "#pathToInputFile(addData == Y)"
          "./bin/CodeSearchLoadData Y Y 8 localhost Y 64 1 /home/ubuntu/Program_output_json.json"
       << std::endl;


  if (argc > 1) {
    if (strcmp(argv[1], "N") == 0) {
      printResult = false;
      COUT << "You successfully disabled printing result." << std::endl;
    } else {
      printResult = true;
      COUT << "Will print result." << std::endl;
    }
  }
  else
  {
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
  }
  else {
    COUT << "Will run on local node. If you want to run on cluster, you can "
            "add any character "
            "as the second parameter to run on the cluster configured by "
            "$PDB_HOME/conf/serverlist."
         << std::endl;
  }

  int blocksize = 64; // by default we add 64MB data
  if (argc > 3) {
    blocksize = atoi(argv[3]);
  }
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

  int dim = 64;
  if (argc > 6) {
    dim = std::stoi(argv[6]);
  }
  COUT << "The dimension of each data point: " << dim << std::endl;


  int Tk = 3;
  if (argc > 7) {
    Tk = std::stoi(argv[7]);
  }
  COUT << "The number of TopK points: " << Tk << std::endl;

  std::string fileName = "/home/ubuntu/Program_output_json.json";
  if (argc > 8) {
    fileName = argv[8];
  }

  COUT << blue << std::endl;
  COUT << "*****************************************" << std::endl;
  COUT << "Code Search starts : " << std::endl;
  COUT << "*****************************************" << std::endl;
  COUT << reset << std::endl;

  COUT << "The Code Search parameters are: " << std::endl;
  COUT << std::endl;

  int numData;

  COUT << "Input file: " << fileName << std::endl;
  COUT << std::endl;

  COUT << std::endl;

  //***********************************************************************************
  //*********************** LOAD DATA****************
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

  auto begin = std::chrono::high_resolution_clock::now();
  FILE* fp = fopen("/home/ubuntu/Program_output_json.json", "r");


  int k;
  if (whetherToAddData == true) {
      // now, create a new database
      pdbClient.createDatabase("code_search_db13");
      // now, create a new set in that database

      pdbClient.createSet<SearchProgramData>("code_search_db13", "code_search_input_set");

      char readBuffer[256*256];
      rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
      rapidjson::Document document;
      document.ParseStream(is);

      COUT << "Document Parsed" << std::endl;
      assert(document.IsObject());


      rapidjson::Value& program = document["programs"];

      COUT << program.Size() << endl;
      fflush(stdout);
      assert(program.IsArray());

      bool rollback = false;
      bool end = false;

      rapidjson::SizeType i = 0;
      numData = 0;
      while (!end) {
        pdb::makeObjectAllocatorBlock(blocksize * 1024 * 1024, true);
        pdb::Handle<pdb::Vector<pdb::Handle<SearchProgramData>>> storeMe =
            pdb::makeObject<pdb::Vector<pdb::Handle<SearchProgramData>>>();
          try {

                while (1) {
                    //      std::istringstream iss(line);
                          if (i >= program.Size()) {
                                    end = true;
                                    break;
                          }
                          else
                          {
                                    pdb::Handle<SearchProgramData> myData =
                                        pdb::makeObject<SearchProgramData>(dim);
                                    //std::stringstream lineStream(line);
                                    // First Set A1, A2 and ProbY
                                    myData->setDoubleA1(program[i]["a1"].GetDouble());
                                    myData->setDoubleA2(program[i]["a2"].GetDouble());
                                    myData->setProbY(program[i]["ProbY"].GetDouble());
                                    // Now set B1
                                    k = 0;
                                    for (auto& v : program[i]["b1"].GetArray()){
                                      myData->setDoubleArrB1(k, v.GetDouble());
                                      k+=1;
                                    }
                                    // Now set B2
                                    k = 0;
                                    for (auto& v : program[i]["b2"].GetArray()){
                                      myData->setDoubleArrB2(k, v.GetDouble());
                                      k+=1;
                                    }
                                    // Setting Program Handles
                                    myData->setFilePtr(program[i]["file"].GetString());
                                    myData->setMethod(program[i]["method"].GetString());
                                    // //Setting Program components
                                    if (!program[i]["javadoc"].IsNull())
                                      myData->setJavaDoc(program[i]["javadoc"].GetString());

                                    // std::cout << program[i]["returnType"].GetString() << std::endl;
                                    if (program[i].HasMember("returnType"))
                                        myData->setReturnType(program[i]["returnType"].GetString());

                                    if (program[i].HasMember("body"))
                                        myData->setBody(program[i]["body"].GetString());

                                    if (program[i].HasMember("formalParam")){
                                          const rapidjson::Value& a = program[i]["formalParam"].GetArray();
                                          for(int k = 0; k < a.Size(); k++){
                                            myData->setFormalParams(a[k].GetString());
                                          }
                                    }


                                    storeMe->push_back(myData);
                                    i++;
                                    //std::cout << i << endl;
                                    //fflush(stdout);
                          }

                  } // end while (1)

                  end = true;

                  // send the rest of data at the end, it can happen that the exception
                  // never happens.
                  pdbClient.sendData<SearchProgramData>(
                          std::pair<std::string, std::string>("code_search_input_set",
                                                              "code_search_db13"), storeMe);

                  numData += storeMe->size();
                  COUT << "Added " << storeMe->size() << " Total: " << numData
                       << std::endl;

                  fflush(stdout);
                  pdbClient.flushData();
              } catch (pdb::NotEnoughSpace &n) {
                  pdbClient.sendData<SearchProgramData>(
                          std::pair<std::string, std::string>("code_search_input_set",
                                                          "code_search_db13"), storeMe);

                  numData += storeMe->size();
                  COUT << "Added " << storeMe->size() << " Total: " << numData
                       << std::endl;
                  fflush(stdout);
                  rollback = true;
              }

              PDB_COUT << blocksize << "MB data sent to dispatcher server~~"
                       << std::endl;

      } // while not end
      fclose(fp);


  } // End if - whetherToAddData = true

  pdb::makeObjectAllocatorBlock(blocksize * 1024 * 1024, true);
  //pdbClient.removeSet("code_search_db13", "result");

  auto end = std::chrono::high_resolution_clock::now();

  COUT << "Time loading data: "
       << std::chrono::duration_cast<std::chrono::duration<float>>(end - begin)
              .count()
       << " secs." << std::endl;




   // this is the object allocation block where all of this stuff will reside
   // pdb::makeObjectAllocatorBlock(blocksize * 1024 * 1024, true);

   pdbClient.createSet<TopKQueue<double, SearchProgramData>>("code_search_db13", "result");
   Handle<SearchProgramData> myQuery = makeObject<SearchProgramData>(dim);

   myQuery->setDoubleA1(-0.01);
   for(int i=0;i < dim; i++){
     myQuery->setDoubleArrB1(i,0.01);
   }

   Handle<Computation> myScanSet =
    makeObject<ScanUserSet<SearchProgramData>>("code_search_db13", "code_search_input_set");
   Handle<Computation> myTopK = makeObject<TopProgram>(Tk, myQuery);
   myTopK->setInput(myScanSet);

   Handle<Computation> myWriter = makeObject<ProgramResultWriter>("code_search_db13", "result");
   myWriter->setInput(myTopK);

   std::cout << "Ready to start computations" << std::endl;
   begin = std::chrono::high_resolution_clock::now();
   pdbClient.executeComputations(myWriter);
   end = std::chrono::high_resolution_clock::now();

   std::cout << "The query is executed successfully!" << std::endl;
   float timeDifference =
       (float(std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count())) / (float)1000000000;
   std::cout << "#TimeDuration: " << timeDifference << " Second " << std::endl;




   // now iterate through the result
   SetIterator<TopKQueue<double, Handle<SearchProgramData>>> result =
           pdbClient.getSetIterator<TopKQueue<double, Handle<SearchProgramData>>>("code_search_db13", "result");
   //std::map<int, int> resultMap;

   for (auto& a : result) {
       std::cout << "Got back " << a->size() << " items from the top-k query.\n";
       std::cout << "These items are:\n";

       for (int i = 0; i < a->size(); i++) {
           std::cout << "score: " << (*a)[i].getScore() << "\n";
           std::cout << "data: ";
           (*a)[i].getValue()->print();
           std::cout << "\n\n";
       }
   }

   // for (auto& a : resultMap) {
   //     std::cout << a.first << "  => " << a.second << '\n';
   // }


}




#endif
