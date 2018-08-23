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

#include "rapidjson/filereadstream.h"
#include "rapidjson/istreamwrapper.h"
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
  COUT << "Usage: ./bin/CodeSearchLoadData blocksSize[MB] #managerIp #addData[Y/N]" "#nDimensions" "#topKResults" "#pathToQuery" "pathToOutput" "#pathToInputFile(addData == Y)"
          "./bin/CodeSearchLoadData 8 localhost Y 64 3 /home/ubuntu/QueryProgWEncoding.json /home/ubuntu/TopPrograms.txt /home/ubuntu/DatabaseOfCodes.json"  << std::endl;
  COUT << reset << std::endl;


  int blocksize = 64; // by default we add 64MB data
  if (argc > 1) {
    blocksize = atoi(argv[1]);
  }

  std::string managerIp = "localhost";
  if (argc > 2) {
    managerIp = argv[2];
  }


  bool whetherToAddData = true;
  if (argc > 3) {
    if (strcmp(argv[3], "N") == 0) {
      whetherToAddData = false;
    }
  }

  int dim = 64;
  if (argc > 4) {
    dim = std::stoi(argv[4]);
  }


  int Tk = 3;
  if (argc > 5) {
    Tk = std::stoi(argv[5]);
  }

  std::string queryFileName = "/home/ubuntu/QueryProgWEncoding.json";
  if (argc > 6) {
    queryFileName = argv[6];
  }

  std::string outputFileName = "/home/ubuntu/TopPrograms.json";
  if (argc > 7) {
    outputFileName = argv[7];
  }

  std::string fileName = "/home/ubuntu/DatabaseOfCodes.json.json";
  if (argc > 8) {
    fileName = argv[8];
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
  COUT << "The number of TopK points: " << Tk << std::endl;
  COUT << "Latent Space Dimensions: " << dim << std::endl;
  COUT << "Query File: " << queryFileName << std::endl;
  COUT << "Input file: " << fileName << std::endl;
  COUT << "Output file: " << outputFileName << std::endl;
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



  if (whetherToAddData == true) {
    int numData;
    auto tbegin = std::chrono::high_resolution_clock::now();
      // now, create a new database
      // and, create a new set in that database
      pdbClient.createDatabase("code_search_db0");
      pdbClient.createSet<SearchProgramData>("code_search_db0", "code_search_input_set");

      // read a JSON file
      // document holds a json document retrieved from a Tensorflow Test Code with Programs + latent vactors
      char readBuffer[256*256];
      FILE* fp = fopen(fileName.c_str(), "r");
      rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
      rapidjson::Document document;
      document.ParseStream(is);

      COUT << "Document Parsed" << std::endl;
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

                                    // Setting Program components
                                    if (!program[i]["javadoc"].IsNull())
                                      myData->setJavaDoc(program[i]["javadoc"].GetString());

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
                                                              "code_search_db0"), storeMe);

                  numData += storeMe->size();
                  COUT << "Added " << storeMe->size() << " Total: " << numData
                       << std::endl;

                  pdbClient.flushData();
              } catch (pdb::NotEnoughSpace &n) {
                  pdbClient.sendData<SearchProgramData>(
                          std::pair<std::string, std::string>("code_search_input_set",
                                                          "code_search_db0"), storeMe);

                  numData += storeMe->size();
                  COUT << "Added " << storeMe->size() << " Total: " << numData
                       << std::endl;
              } // end try-catch

              PDB_COUT << blocksize << "MB data sent to dispatcher server~~"
                       << std::endl;

      } // while not end
      fclose(fp);


    auto tend = std::chrono::high_resolution_clock::now();

    COUT << "Time loading data: "
         << std::chrono::duration_cast<std::chrono::duration<float>>(tend - tbegin).count()
         << " secs." << std::endl;

   } // End if - whetherToAddData = true

   //***********************************************************************************
   //*********************************LOAD QUERY********************************
   //***********************************************************************************
   //***********************************************************************************
   // Use the JAVA Eclipse DOM driver to get the Query as a JSON and parse it similar to before.

   // Step 2. Load the QUERY
   // Load the Query and Get Top-k
   pdb::makeObjectAllocatorBlock(blocksize * 1024 * 1024, true);
   pdb::Handle<SearchProgramData> myQuery = pdb::makeObject<SearchProgramData>(dim);

   /*
   myQuery->setDoubleA1(-0.01);
   for(int i=0;i < dim; i++){
     myQuery->setDoubleArrB1(i,0.01);
   }*/

	ifstream ifs(queryFileName.c_str());
	rapidjson::IStreamWrapper isw(ifs);


   rapidjson::Document queryDocument;
   queryDocument.ParseStream<0>(isw);
   COUT << "Query Document Parsed" << std::endl;
   assert(queryDocument.IsObject());
   rapidjson::Value& queryProgram = queryDocument["programs"];

   assert(queryProgram.Size() == 1);
   assert(queryProgram.IsArray());


	// First Set A1, A2 and ProbY
	myQuery->setDoubleA1(queryProgram[0]["a1"].GetDouble());
	// Now set B1
	int k = 0;
	for (auto& v : queryProgram[0]["b1"].GetArray()){
	  myQuery->setDoubleArrB1(k, v.GetDouble());
	  k+=1;
	}
	/*
	// Setting queryProgram Handles

	myQuery->setFilePtr(queryProgram["file"].GetString());
	myQuery->setMethod(queryProgram["method"].GetString());

	// Setting queryProgram components
	if (!queryProgram["javadoc"].IsNull())
	  myQuery->setJavaDoc(queryProgram["javadoc"].GetString());

	if (queryProgram[i].HasMember("returnType"))
		myQuery->setReturnType(queryProgram["returnType"].GetString());

	if (queryProgram[i].HasMember("body"))
		myQuery->setBody(queryProgram["body"].GetString());

	if (queryProgram[i].HasMember("formalParam")){
		  const rapidjson::Value& a = queryProgram["formalParam"].GetArray();
		  for(int k = 0; k < a.Size(); k++){
			myQuery->setFormalParams(a[k].GetString());
		  }
	}*/
   //Need to call the destructor for myQuery


   //***********************************************************************************
   //*********************************CREATE COMPUTATION********************************
   //***********************************************************************************
   //***********************************************************************************

   // Step 3. Create the Computation Graph
   // Load the Query and Get Top-k


   pdbClient.createSet<TopKQueue<double, SearchProgramData>>("code_search_db0", "result");

   Handle<Computation> myScanSet =
        makeObject<ScanUserSet<SearchProgramData>>("code_search_db0", "code_search_input_set");
   Handle<Computation> myTopK = makeObject<TopProgram>(Tk, myQuery);
   myTopK->setInput(myScanSet);

   Handle<Computation> myWriter = makeObject<ProgramResultWriter>("code_search_db0", "result");
   myWriter->setInput(myTopK);

   std::cout << "Ready to start computations" << std::endl;
   tbegin = std::chrono::high_resolution_clock::now();
   pdbClient.executeComputations(myWriter);
   tend = std::chrono::high_resolution_clock::now();

   std::cout << "The query is executed successfully!" << std::endl;
   float timeDifference =
       (float(std::chrono::duration_cast<std::chrono::nanoseconds>(tend - tbegin).count())) / (float)1000000000;
   std::cout << "#TimeDuration: " << timeDifference << " Second " << std::endl;

   // now iterate through the result
   SetIterator<TopKQueue<double, Handle<SearchProgramData>>> result =
           pdbClient.getSetIterator<TopKQueue<double, Handle<SearchProgramData>>>("code_search_db0", "result");


   FILE* fout = fopen(outputFileName.c_str(), "w");
   for (auto& a : result) {
       std::cout << "Got back " << a->size() << " items from the top-k query.\n";
       std::cout << "These items are:\n";

       for (int i = 0; i < a->size(); i++) {
           std::cout << "score: " << (*a)[i].getScore() << "\n";
           std::cout << "data: ";
           (*a)[i].getValue()->print();
           (*a)[i].getValue()->fprint(fout);
           std::cout << "\n\n";
       }
   }
   pdbClient.removeSet("code_search_db0", "result");

}

#endif
