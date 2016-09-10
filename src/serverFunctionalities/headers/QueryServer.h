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

#ifndef QUERY_SERVER_H
#define QUERY_SERVER_H

#include "ServerFunctionality.h"
#include "PDBServer.h"
#include "Record.h"
#include "MyDB_Catalog.h"
#include <vector>
#include "PDBVector.h"
#include "QueryBase.h"
#include "MyDB_BufferManager.h"

namespace pdb {

class QueryServer : public ServerFunctionality {

public:

	// creates a query server... the param is the number of threads to use
	// to answer queries
	QueryServer (int numThreads);

	// from the ServerFunctionality interface... registers the QueryServer's 
	// single handler, which accepts a query and executes it
	void registerHandlers (PDBServer &forMe) override;

	// this recursively traverses a simple query graph, where each node can only have one input,
	// makes sure that each node has been computed... setPrefix is the string that we'll use
	// to create each set name, whichNode is the counter that we use to name each set, and
	// computeMe is the node that we are wirred about computing
	void computeQuery (std :: string setPrefix, int &whichNode, Handle <QueryBase> &computeMe);

	// this actually computes a selection query.  The params are like the above function
	void doSelection (std :: string setPrefix, int whichNode, Handle <QueryBase> &computeMe);

	// destructor
	~QueryServer ();

private:

	// the number of threads to use
	int numThreadsToUse;
		
	// used to count up temporary file names: tempSet0, tempSet1, tempSet2, ...
	int tempSetName;
};

}

#endif
