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
// By Rohan, August 2018
// CodeSearch


#include <chrono>
#include <cstddef>
#include <ctime>

#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <random>

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


  std::string queryFileName = "/home/ubuntu/query.json";
  if (argc > 1) {
    queryFileName = argv[1];
  }
  
  std::string fileName = "/home/ubuntu/Program_output_json.json";
  if (argc > 2) {
    fileName = argv[2];
  }





  // For the random number generator
  //    srand(time(0));
  std::random_device rd;
  std::mt19937 randomGen(rd());

  //***********************************************************************************
  //************************************READ INPUT DATA********************************
  //**********************************************************************************
  //***********************************************************************************

  // Step 1. Create Database and Set
  // now, register a type for user data

  auto tbegin = std::chrono::high_resolution_clock::now();


  int numData;



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

      COUT << " Total Number of Programs in JSON :: "<< program.Size() << endl;
      assert(program.IsArray());


      rapidjson::SizeType i = rand() % program.Size();
	  rapidjson::Document d;
	  d.parse(programs[i])
	
		FILE* fp = fopen("/home/ubuntu/query.json", "w"); 
		char writeBuffer[65536];
		FileWriteStream os(fp, writeBuffer, sizeof(writeBuffer));
		Writer<FileWriteStream> writer(os);
		d.Accept(writer);
		fclose(fp);

	  /*
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
                                    
		*/

      fclose(fp);

}

