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

namespace pdb{

// SearchProgramData
class SearchProgramData : public Object {

private:
  Handle<Vector<double>> LatentDoubleVec = nullptr;
  double ProbY = -1000.00;
  size_t size = 0;
  String Prog = "";

public:
  ENABLE_DEEP_COPY

  SearchProgramData() {}

  SearchProgramData(size_t dim) {
    this->size = dim;
    LatentDoubleVec = makeObject<Vector<double>>(size, size);
     //this->LatentDoubleVec.resize(dim);
  }

  void setDouble(int i,double value){

    if (i < this->size) {
        (*LatentDoubleVec)[i] = value;
    } else {
        std::cout << "Error in SearchProgramData: Cannot assign the value " << value << "to the pos "
                  << i << std::endl;
        exit(-1);
    }
  }

  void setProbY(double value){
     this->ProbY = value;
  }


  void setProg(String inpProg){
    this->Prog = inpProg;
  }

  String getProg(){
    return this->Prog;
  }

  size_t getSize() {
      return this->size;
  }

  void print() {
      double* rawData = LatentDoubleVec->c_ptr();
      for (int i = 0; i < this->getSize(); i++) {
          std::cout << i << ": " << rawData[i] << "; ";
      }
      std::cout << std::endl;
  }

  // double* getRawData() {
  //     if (data == nullptr) {
  //         return nullptr;
  //     }
  //     return data->c_ptr();
  // }
  //


  // to get squared distance following SparkMLLib

  double getSquaredDistance(Vector<double> queryProg) {

      return 0.0;
      // # a1 = tf.placeholder(tf.float32,[self.model.config.latent_size])
      // a1_in = tf.placeholder(tf.float32,[])
      // a1 = tf.tile(tf.reshape(a1_in,[1]),[self.model.config.latent_size])
      //
      // b1_in = tf.placeholder(tf.float32,[self.model.config.latent_size])
      // b1 = b1_in
      // # a2 = tf.placeholder(tf.float32,[self.model.config.batch_size,self.model.config.latent_size])
      // a2_in = tf.placeholder(tf.float32,[self.model.config.batch_size])
      // a2 = tf.tile(tf.expand_dims(a2_in,axis=1),[1,self.model.config.latent_size])
      //
      // b2_in = tf.placeholder(tf.float32,[self.model.config.batch_size,self.model.config.latent_size])
      // b2 = b2_in
      // t1 = tf.reduce_sum(tf.square(b1)/(4*a1), axis=0) + 0.5 * tf.reduce_sum(tf.log(-a1/np.pi), axis=0)
      // t2 = tf.reduce_sum(tf.square(b2)/(4*a2), axis=1) + 0.5 * tf.reduce_sum(tf.log(-a2/np.pi), axis=1)
      // t3 = 0.5 * self.model.config.latent_size * tf.log(2*np.pi)
      // c = t1 + t2 - t3
      //
      // b_star = b1 + b2
      // a_star = a1 + a2 + 0.5
      // c_star = tf.reduce_sum(tf.square(b_star)/(4*a_star), axis=1) + 0.5 * tf.reduce_sum(tf.log(-a_star/np.pi), axis=1)
      // prob = (c - c_star)
      //
      // _prob = self.sess.run(prob, feed_dict={a1_in:_a1, b1_in:_b1, a2_in:_a2, b2_in:_b2})
      //
      // _prob += np.array(prob_Y)
      // return _prob
  }



  ~SearchProgramData() {}

};
}
#endif
