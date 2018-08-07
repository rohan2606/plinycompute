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
#ifndef TEST_LDA1_CC
#define TEST_LDA1_CC

#include "PDBDebug.h"
#include "PDBVector.h"
#include "Query.h"
#include "Lambda.h"
#include "PDBClient.h"
#include "Set.h"
#include "DataTypes.h"
#include "TopKQueue.h"
#include "sharedLibraries/headers/ScanCustomerSet.h"
#include "sharedLibraries/headers/TopJaccard.h"
#include "sharedLibraries/headers/JaccardResultWriter.h"

#include <ctime>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>
#include <fcntl.h>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <sstream>

using namespace pdb;
int main(int argc, char* argv[]) {
    std::ofstream term("/dev/tty", std::ios_base::out);

    if (argc != 4) {
        std::cout << "Usage: #managerIp #k #queryFile\n";
        return (-1);
    }

    std::string managerIp;
    managerIp = argv[1];
    std::cout << "Manager IP Address is " << managerIp << std::endl;

    int k = std::stoi(argv[2]);
    term << "The number of query results: " << k << std::endl;

    std::string queryFile;
    queryFile = argv[3];
    std::ifstream read(queryFile);
    int next;
    pdb::makeObjectAllocatorBlock(24 * 1024 * 1024, true);
    Handle<Vector<int>> myQuery = makeObject<Vector<int>>();
    while (read >> next) {
        myQuery->push_back(next);
    }

    PDBClient pdbClient(8108, managerIp);

    string errMsg;
    std::vector<std::string> v = {"libraries/libCustomer.so",
                                  "libraries/libLineItem.so",
                                  "libraries/libOrder.so",
                                  "libraries/libPart.so",
                                  "libraries/libSupplier.so",
                                  "libraries/libScanCustomerSet.so",
                                  "libraries/libTopJaccard.so",
                                  "libraries/libAllParts.so",
                                  "libraries/libJaccardResultWriter.so"};

    for (auto& a : v) {
        pdbClient.registerType(a);
    }

    pdbClient.removeSet("TPCH_db", "result");
    pdbClient.createSet<TopKQueue<double, AllParts>>("TPCH_db", "result");

    // Create computations
    Handle<Computation> myScanSet = makeObject<ScanCustomerSet>("TPCH_db", "tpch_bench_set1");
    Handle<Computation> myTopK = makeObject<TopJaccard>(k, *myQuery);
    myTopK->setInput(myScanSet);

    Handle<Computation> myWriter = makeObject<JaccardResultWriter>("TPCH_db", "result");
    myWriter->setInput(myTopK);

    auto begin = std::chrono::high_resolution_clock::now();
    pdbClient.executeComputations(myWriter);
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "The query is executed successfully!" << std::endl;
    float timeDifference =
        (float(std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count())) /
        (float)1000000000;
    std::cout << "#TimeDuration: " << timeDifference << " Second " << std::endl;

    // now iterate through the result
    SetIterator<TopKQueue<double, Handle<AllParts>>> result =
            pdbClient.getSetIterator<TopKQueue<double, Handle<AllParts>>>("TPCH_db", "result");
    std::map<int, int> resultMap;

    for (auto& a : result) {
        std::cout << "Got back " << a->size() << " items from the top-k query.\n";
        std::cout << "These items are:\n";

        for (int i = 0; i < a->size(); i++) {
            std::cout << "score: " << (*a)[i].getScore() << "\n";
            std::cout << "data: ";
            (*a)[i].getValue()->print();
            if (resultMap.count((*a)[i].getValue()->custKey) == 0) {
                resultMap[(*a)[i].getValue()->custKey] = 1;
            } else {
                resultMap[(*a)[i].getValue()->custKey] = resultMap[(*a)[i].getValue()->custKey] + 1;
            }
            std::cout << "\n\n";
        }
    }

    for (auto& a : resultMap) {
        std::cout << a.first << "  => " << a.second << '\n';
    }
}

#endif
