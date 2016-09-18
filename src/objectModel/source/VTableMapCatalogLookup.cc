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

#ifndef VTABLEMAP_CAT_LOOKUP_CC
#define VTABLEMAP_CAT_LOOKUP_CC

/*
 * VTableMap.cc
 *
 *  Created on: Dec 11, 2015
 *      Author: carlos
 */

// This header is auto-generated by SConstruct... it includes all of the header files
// in BuiltInPDBObjects/headers

#include <dlfcn.h>
#include <unistd.h>
#include <unistd.h>
#include "PDBLogger.h"
#include <cctype>
#include "CatalogClient.h"

namespace pdb {

// note: this should only be called while protected by a lock on the vTableMap
void *VTableMap :: getVTablePtrUsingCatalog (int16_t objectTypeID) {

	// in ths case, we do not have the vTable pointer for this guy, so we will try to load it
	if (theVTable->catalog == nullptr) { 
		if (theVTable->logger != nullptr)
			theVTable->logger->writeLn ("unable to obtain shared library file.");
		return nullptr;
	}

	std :: string sharedLibraryFile = "/var/tmp/objectFile.";
	sharedLibraryFile += to_string (getpid ()) + "." + to_string (objectTypeID) + ".so";
	unlink (sharedLibraryFile.c_str ());
	theVTable->catalog->getSharedLibrary (objectTypeID, sharedLibraryFile);
	
	// open up the shared library
	void *so_handle = dlopen(sharedLibraryFile.c_str(), RTLD_LOCAL | RTLD_LAZY );
	theVTable->so_handles.push_back (so_handle);

	if (!so_handle) {
		const char* dlsym_error = dlerror();
		if (theVTable->logger != nullptr)
			theVTable->logger->error ("Cannot load Stored Data Type library: " + sharedLibraryFile + 
				" error " + (std::string)dlsym_error + '\n');

	// if we were able to open it
	} else {
		const char* dlsym_error = dlerror();

		// first we need to correctly set all of the global variables in the shared library
		typedef void setGlobalVars (Allocator *, VTableMap *, void *, void *);
		std::string getInstance = "setAllGlobalVariables";
		setGlobalVars *setGlobalVarsFunc = (setGlobalVars *) dlsym (so_handle, getInstance.c_str());

		// see if we were able to get the function
		if ((dlsym_error = dlerror())) {
			if (theVTable->logger != nullptr)
				theVTable->logger->error ("Error, can't set global variables in .so file; error is " +
					(std::string)dlsym_error + "\n");
			return nullptr;
		// if we were able to, then run it
		} else {
			setGlobalVarsFunc (mainAllocatorPtr, theVTable, stackBase, stackEnd);
		}

		// get the function that will give us access to the vTable
		typedef void *getObjectVTable ();
		getInstance = "getObjectVTable";
		getObjectVTable *getObjectFunc = (getObjectVTable *) dlsym (so_handle, getInstance.c_str());

		// see if we were able to get the function
		if ((dlsym_error = dlerror())) {
			if (theVTable->logger != nullptr)
				theVTable->logger->error ("Error, can't load function getInstance (); error is " +
					(std::string)dlsym_error + "\n");
			return nullptr;
		// if we were able to, then run it
		} else {
			theVTable->allVTables[objectTypeID] = getObjectFunc ();
		}

	}
	return theVTable->allVTables[objectTypeID];
}

int16_t VTableMap :: lookupTypeNameInCatalog (std :: string objectTypeName) {
        return theVTable->catalog->searchForObjectTypeName (objectTypeName);
}


} /* namespace pdb */

#endif
