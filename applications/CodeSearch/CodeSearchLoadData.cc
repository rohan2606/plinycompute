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
  COUT << "Usage: ./bin/CodeSearchLoadData blocksSize[MB] #managerIp #createDB[Y/N]" "#nDimensions" "#pathToInputFile(addData == Y)"
          "./bin/CodeSearchLoadData 8 localhost Y 128 /home/ubuntu/DATABASE/DatabaseOfCodes-01.json"  << std::endl;
  COUT << reset << std::endl;


  int blocksize = 64; // by default we add 64MB data
  if (argc > 1) {
    blocksize = atoi(argv[1]);
  }

  std::string managerIp = "localhost";
  if (argc > 2) {
    managerIp = argv[2];
  }


  bool createDatabase = true;
  if (argc > 3) {
    if (strcmp(argv[3], "N") == 0) {
      createDatabase = false;
    }
  }

  int dim = 128;
  if (argc > 4) {
    dim = std::stoi(argv[4]);
  }

  std::string fileName = "/home/ubuntu/DatabaseOfCodes.json";
  if (argc > 5) {
    fileName = argv[5];
  }

  COUT << blue << std::endl;
  COUT << "*****************************************" << std::endl;
  COUT << "Code Search starts : " << std::endl;
  COUT << "*****************************************" << std::endl;
  COUT << reset << std::endl;

  COUT << green << std::endl;
  COUT << "*****************************************" << std::endl;
  COUT << "The Code Search parameters are: " << std::endl;
  COUT << "Manager IP Address is " << managerIp << std::endl;
  COUT << "To add data with size: " << blocksize << "MB" << std::endl;
  COUT << "The Dimension of each data point: " << dim << std::endl;
  COUT << "Input file: " << fileName << std::endl;
  COUT << "*****************************************" << std::endl;
  COUT << reset << std::endl;


  //***********************************************************************************
  //*******************************  REGISTER TYPES ***********************************
  //***********************************************************************************

  PDBClient pdbClient(8108, managerIp);
  pdbClient.registerType("libraries/libSearchProgramData.so");
  pdbClient.registerType("libraries/libTopProgram.so");
  pdbClient.registerType("libraries/libProgramResultWriter.so");
  string errMsg;

  //***********************************************************************************
  //************************************READ INPUT DATA********************************
  //**********************************************************************************
  //***********************************************************************************

  // Step 1. Create Database and Set
  // now, register a type for user data



	  if (createDatabase == true) {
		  // now, create a new database
		  // and, create a new set in that database
		  pdbClient.createDatabase("code_search_db");
		  pdbClient.createSet<SearchProgramData>("code_search_db", "code_search_input_set");
		}
      // read a JSON file
      // document holds a json document retrieved from a Tensorflow Test Code with Programs + latent vactors
	  
	  
	  auto tbegin = std::chrono::high_resolution_clock::now();
	  int numData;

      char readBuffer[256*256];
      FILE* fp = fopen(fileName.c_str(), "r");
      rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
      rapidjson::Document document;
      document.ParseStream(is);

      COUT << "Document Parsed" << std::endl;
	  fclose(fp);
      assert(document.IsObject());

      rapidjson::Value& program = document["programs"];
      bool end = false;

      COUT << " Total Number of Programs in JSON :: "<< program.Size() << endl;
      assert(program.IsArray());

      rapidjson::SizeType i = 0;
      numData = 0;
      while (!end) {
        pdb::makeObjectAllocatorBlock(blocksize * 1024 * 1024, true);
        pdb::Handle<pdb::Vector<pdb::Handle<SearchProgramData>>> storeMe =
            pdb::makeObject<pdb::Vector<pdb::Handle<SearchProgramData>>>();
          try {

                while (1) {
                          if (i >= program.Size()) {
                                    end = true;
                                    break;
                          }
                          else
                          {
                                    pdb::Handle<SearchProgramData> myData =
                                        pdb::makeObject<SearchProgramData>(dim);
                                    // First Set A1, A2 and ProbY
                                    myData->setDoubleA1(program[i]["a1"].GetDouble());
                                    myData->setDoubleA2(program[i]["a2"].GetDouble());
                                    myData->setProbY(program[i]["ProbY"].GetDouble());
                                    // Now set B1
                                    int k = 0;
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
                          }

                  } // end while (1)

                  end = true;
                  // send the rest of data at the end, it can happen that the exception
                  // never happens.
                  pdbClient.sendData<SearchProgramData>(
                          std::pair<std::string, std::string>("code_search_input_set",
                                                              "code_search_db"), storeMe);

                  numData += storeMe->size();
                  COUT << "Added " << storeMe->size() << " Total: " << numData
                       << std::endl;

                  pdbClient.flushData();
              } catch (pdb::NotEnoughSpace &n) {
                  pdbClient.sendData<SearchProgramData>(
                          std::pair<std::string, std::string>("code_search_input_set",
                                                          "code_search_db"), storeMe);

                  numData += storeMe->size();
                  COUT << "Added " << storeMe->size() << " Total: " << numData
                       << std::endl;
              } // end try-catch

              PDB_COUT << blocksize << "MB data sent to dispatcher server~~"
                       << std::endl;

      } // while not end


    auto tend = std::chrono::high_resolution_clock::now();

    COUT << "Time loading data: "
         << std::chrono::duration_cast<std::chrono::duration<float>>(tend - tbegin).count()
         << " secs." << std::endl;


  

}

#endif
