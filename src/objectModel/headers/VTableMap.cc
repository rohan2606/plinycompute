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

#include "Object.h"
#include "Nothing.h"

#ifndef PDBCATALOG_VTABLEMAP_C_
#define PDBCATALOG_VTABLEMAP_C_

/*
 * VTableMap.cc
 *
 *  Created on: Dec 11, 2015
 *      Author: carlos
 */

// This header is auto-generated by SConstruct... it includes all of the header files
// in BuiltInPDBObjects/headers

#include <dlfcn.h>
#include "LockGuard.h"
#include <unistd.h>
#include "PDBLogger.h"
#include <cctype>

namespace pdb {

inline int16_t VTableMap :: lookupBuiltInType (std :: string objectTypeName) {
	if (theVTable->objectTypeNamesList.count (objectTypeName) != 0)
		if (theVTable->objectTypeNamesList[objectTypeName] >= 0 && theVTable->objectTypeNamesList[objectTypeName] <= 8191)
			return theVTable->objectTypeNamesList[objectTypeName];
	
	// otherwise, return -1 to signify that this was not a built-in type
	return -1;
}

inline std :: string VTableMap :: lookupBuiltInType (int16_t objectType) {

	// loop through every type in the mape
	for (auto &a : theVTable->objectTypeNamesList) {
		if (a.second == objectType)
			return a.first;
	}
	
	// otherwise, return "" to signify that this was not a built-in type
	return "";
}

// returns the number of registered built-in objects
inline int VTableMap :: totalBuiltInObjects () {
        int count = 0;
        for (int i = 0; i < theVTable->allVTables.size();i++){
            if (theVTable->allVTables[i] != nullptr){
                count = count + 1;
            }
        }
        return count;
}

inline void VTableMap :: listVtableEntries () {
        for (int i=0;i< theVTable->allVTables.size();i++){
            if (theVTable->allVTables[i] != nullptr) std::cout << "vtpr " << i << ": " << theVTable->allVTables[i] << std::endl;
        }
}    

inline void VTableMap :: listVtableLabels () {
        std::map<std::string, int16_t> iterator;

        for(auto &iterator : theVTable->objectTypeNamesList){
            std::cout << "label= " << iterator.first << " id= " <<  iterator.second << std::endl;
        }
}

inline int16_t VTableMap :: getIDByName (std::string objectTypeName) {

	// one important issue is that we might need to lookup soething nasty like:
	//
	// pdb::PairArray<pdb::Handle<pdb::String>,pdb::Handle<pdb::Employee>>
	//
	// this needs to be normalized to:
	//
	// pdb::PairArray<<pdb::Nothing>,<pdb::Nothing>?
	//
	// basically, we have a stack that pushes and pops the depth of the template by counting instances of '>' and '<'
	// Any time that we see a "," at depth 1, it means we've hit an additional template arguement, and so we add
	// a "<pdb::Nothing>" to the replacement string
	std :: string replacementString ("");
	std :: string prefix ("");
	bool isTemplate = false;
	int depth = 0;
	int length = objectTypeName.size ();
	for (int loc = 0; loc < length; loc++) {
		if (objectTypeName[loc] == '<') {
			depth++;
			if (depth == 1) {
				isTemplate = true;
				prefix = objectTypeName.substr (0, loc);
				replacementString = std :: string ("pdb::Nothing");
			}
		} else if (objectTypeName[loc] == '>') {
			depth--;
		} else if (objectTypeName[loc] == ',' && depth == 1) {
			std :: string nextReplacement (",pdb::Nothing");
			replacementString += nextReplacement;	
		}	
	}

	// if this was a template, do the normalization
	if (isTemplate)
		objectTypeName = prefix + std :: string ("<") + replacementString + std :: string (">");

	// now, check to make sure that we have seen the given object type before
	if (theVTable->objectTypeNamesList.count (objectTypeName) == 0 && theVTable->catalog != nullptr) {
		
		// make sure no one is modifying the map
		const LockGuard guard {theVTable->myLock};
	
		// in this case, we do not have this object type, and we have never looked for it before
		// so, go to the catalog and ask for it...
		int16_t identifier = lookupTypeNameInCatalog (objectTypeName); 

		// if the identifier is -1, then it means the catalog has never seen this type before
		// so let the caller know, and remember that we have not seen it
		if (identifier == -1) {
			theVTable->objectTypeNamesList[objectTypeName] = -1;
			return -1;

		// otherwise, return the ID
		} else {
			theVTable->objectTypeNamesList[objectTypeName] = identifier;
			return identifier;
		}

	// in the easy case, we have seen it before, so just return the typeID
	} else {
		return theVTable->objectTypeNamesList[objectTypeName];	
	}
}

inline void VTableMap :: setCatalogClient (CatalogClient *catalogIn) {
	const LockGuard guard {theVTable->myLock};
	theVTable->catalog = catalogIn;
}

inline CatalogClient * VTableMap :: getCatalogClient () {
	const LockGuard guard {theVTable->myLock};
	return theVTable->catalog;
}

extern bool inSharedLibrary;

inline VTableMap :: ~VTableMap () {
	const LockGuard guard {theVTable->myLock};
	if (!inSharedLibrary)
		for (void *v : theVTable->so_handles) {
			int res = dlclose (v);
			if (res != 0)
				std :: cout << dlerror () << "\n";
                        
		} 
	theVTable->so_handles.clear ();
}

inline void *VTableMap :: getVTablePtr (int16_t objectTypeID) {

	// OK, first, we check to see if we have the v table pointer for this guy...
	// this is done without a lock, so we can be very fast...
	void *returnVal = theVTable->allVTables[objectTypeID];
	if (returnVal != nullptr) {
                //std :: cout << "VTablePtr for objectTypeID("<<objectTypeID<<") exists" << std :: endl;
		return returnVal;}

	// we do not, so get the lock...
	const LockGuard guard {theVTable->myLock};

	// before we go out to the network for the v table pointer, just verify
	// that another thread has not since gotten it for us
	returnVal = theVTable->allVTables[objectTypeID];
	if (returnVal != nullptr) {
                //std :: cout << "VTablePtr for objectTypeID=" << objectTypeID << " exists" << std :: endl;
		return returnVal;}
	else {
                //std :: cout << "VTablePtr for objectTypeID=" << objectTypeID << " doesn't exist" << std :: endl;
		// if they have not gotten it for us, then go and get it
		return getVTablePtrUsingCatalog (objectTypeID);
        }
}

} /* namespace pdb */

#endif
