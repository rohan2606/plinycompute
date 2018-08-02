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
#ifndef SEARCHPROGRAMDATA_H
#define SEARCHPROGRAMDATA_H

#include "Handle.h"
#include "LambdaCreationFunctions.h"
#include "Object.h"
#include "PDBVector.h"
#include <vector>

using namespace pdb;

// SearchProgramData
class SearchProgramData : public Object {

private:
  Handle<Vector<double>> DoubleVec = nullptr;
  size_t size = 0;

public:
  ENABLE_DEEP_COPY

  SearchProgramData() {}

  SearchProgramData(size_t dim) {
    this->size = dim;
    DoubleVec = makeObject<Vector<double>>(size, size);
     //this->DoubleVec.resize(dim);
  }

  void setDouble(int i,double value){

    if (i < this->size) {
        (*DoubleVec)[i] = value;
    } else {
        std::cout << "Error in SearchProgramData: Cannot assign the value " << value << "to the pos "
                  << i << std::endl;
        exit(-1);
    }
  }

  size_t getSize() {
      return this->size;
  }

  void print() {
      double* rawData = DoubleVec->c_ptr();
      for (int i = 0; i < this->getSize(); i++) {
          std::cout << i << ": " << rawData[i] << "; ";
      }
      std::cout << std::endl;
  }

  ~SearchProgramData() {}
};

#endif
