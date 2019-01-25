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
#ifndef CODE_SEARCH
#define CODE_SEARCH

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
#include <vector>

#include <string>
#include <streambuf>
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/writer.h"
#include "rapidjson/document.h"
#include "rapidjson/stringbuffer.h"


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
  COUT << "Usage: ./bin/CodeSearch blocksSize[MB] #managerIp" "#nDimensions" "#topKResults" "#pathToQuery" "pathToOutput"
          "./bin/CodeSearch 8 localhost 128 10 /home/ubuntu/QueryProgWEncoding.json /home/ubuntu/TopPrograms.txt"  << std::endl;
  COUT << reset << std::endl;


  int blocksize = 64; // by default we add 64MB data
  if (argc > 1) {
    blocksize = atoi(argv[1]);
  }

  std::string managerIp = "localhost";
  if (argc > 2) {
    managerIp = argv[2];
  }

  int dim = 128;
  if (argc > 3) {
    dim = std::stoi(argv[3]);
  }

  int Tk = 10;
  if (argc > 4) {
    Tk = std::stoi(argv[4]);
  }

  std::string queryFileName = "/home/ubuntu/QueryProgWEncoding.json";
  if (argc > 5) {
    queryFileName = argv[5];
  }

  std::string outputFileName = "/home/ubuntu/TopPrograms.json";
  if (argc > 6) {
    outputFileName = argv[6];
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
  COUT << "Input Set: code_search_input_set" << std::endl;
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


   //***********************************************************************************
   //*********************************LOAD QUERY********************************
   //***********************************************************************************
   //***********************************************************************************
   // Use the JAVA Eclipse DOM driver to get the Query as a JSON and parse it similar to before.

   string errMsg;
   // Step 2. Load the QUERY
   // Load the Query and Get Top-k
   pdb::makeObjectAllocatorBlock(blocksize * 1024 * 1024, true);


   pdb::Handle<SearchProgramData> myQuery = pdb::makeObject<SearchProgramData>(dim);
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

   //Need to call the destructor for myQuery


   //***********************************************************************************
   //*********************************CREATE COMPUTATION********************************
   //***********************************************************************************
   //***********************************************************************************

   // Step 3. Create the Computation Graph
   // Load the Query and Get Top-k



   pdbClient.createSet<TopKQueue<double, SearchProgramData>>("code_search_db", "result0");

   Handle<Computation> myScanSet =
        makeObject<ScanUserSet<SearchProgramData>>("code_search_db", "code_search_input_set");
   Handle<Computation> myTopK = makeObject<TopProgram>(Tk, myQuery);
   myTopK->setInput(myScanSet);

   Handle<Computation> myWriter = makeObject<ProgramResultWriter>("code_search_db", "result0");
   myWriter->setInput(myTopK);

   std::cout << "Ready to start computations" << std::endl;
   auto tbegin = std::chrono::high_resolution_clock::now();
   pdbClient.executeComputations(myWriter);
   auto tend = std::chrono::high_resolution_clock::now();

   std::cout << "The query is executed successfully!" << std::endl;
   float timeDifference =
       (float(std::chrono::duration_cast<std::chrono::nanoseconds>(tend - tbegin).count())) / (float)1000000000;
   std::cout << "#TimeDuration: " << timeDifference << " Second " << std::endl;

   // now iterate through the result0
   SetIterator<TopKQueue<double, Handle<SearchProgramData>>> result0 =
           pdbClient.getSetIterator<TopKQueue<double, Handle<SearchProgramData>>>("code_search_db", "result0");


   // FILE* fout = fopen(outputFileName.c_str(), "w");
   for (auto& a : result0) {
       std::cout << "Got back " << a->size() << " items from the top-k query.\n";
	   std::vector<SearchProgramData> TopProgs;
	   for (int i = 0; i < a->size(); i++) {
		   TopProgs.push_back( *((*a)[i].getValue()) );
		   TopProgs.at(i).setPostProb((*a)[i].getScore());
	   }
	   std::sort(TopProgs.begin(), TopProgs.end());
	   std::reverse(TopProgs.begin(), TopProgs.end());
     // for (auto& prog : TopProgs) {
     //     prog.fprint(fout);
     // }
     rapidjson::StringBuffer s;
     rapidjson::Writer<rapidjson::StringBuffer> writer(s);

     writer.StartObject();
     writer.Key("TopPrograms");
     writer.StartArray();                // Between StartArray()/EndArray(),
     int count = 0;
     for (auto& prog : TopProgs){
        count += 1;
        writer.Key("Rank");
        writer.Int(count);
        writer.Key("Prob");
        writer.Double(prog->PostProb);
        writer.Key("File");
        writer.Double(prog->filePtr.c_str());
        writer.Key("Method");
        writer.Double(prog->method.c_str());
        writer.Key("Body");
        writer.Double(prog->body.c_str());
     }
                         // all values are elements of the array.
     writer.EndArray();
     writer.EndObject();
   }

   std::ofstream of(outputFileName.c_str());
   of << s.GetString();
   if (!of.good())
        throw std::runtime_error("Can't write the JSON string to the file!");

   pdbClient.removeSet("code_search_db", "result0");

}

#endif
