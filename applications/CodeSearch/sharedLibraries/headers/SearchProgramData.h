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
#include <cmath>

namespace pdb{

// SearchProgramData
class SearchProgramData : public Object {

private:

  // Search parameters
  double a1, a2, ProbY;
  size_t dim = 0;
  Handle<Vector<double>> b1 = nullptr;
  Handle<Vector<double>> b2 = nullptr;

  //Program Handles
  String filePtr = "";
  String method = "";

  // Program components
  String javadoc = "";
  Handle<Vector<String>> FormalParams =  makeObject<Vector<String>>();
  String returnType = "";
  String body = "";

  double PostProb = -100000.0;
  
public:
  ENABLE_DEEP_COPY

  SearchProgramData() {}

  SearchProgramData(size_t dim) {
    this->dim = dim;
    b1 = makeObject<Vector<double>>(dim, dim);
    b2 = makeObject<Vector<double>>(dim, dim);
  }


  size_t getSize() {
    return this->dim;
  }

  void setDoubleArrB1(int i, String inpValue){
    double value = atof(inpValue.c_str());
    if (i < this->dim) {
        (*b1)[i] = value;
    } else {
        std::cout << "Error in SearchProgramData: Cannot assign the value " << value << "to the pos "
                  << i << std::endl;
        exit(-1);
    }
    return;
  }

  void setDoubleArrB2(int i, String inpValue){
    double value = atof(inpValue.c_str());
    if (i < this->dim) {
        (*b2)[i] = value;
    } else {
        std::cout << "Error in SearchProgramData: Cannot assign the value " << value << "to the pos "
                  << i << std::endl;
        exit(-1);
    }
    return;
  }

  void setDoubleA1(String inpValue){
     double value = atof(inpValue.c_str());
     this->a1 = value;
     return;
  }

  void setDoubleA2(String inpValue){
     double value = atof(inpValue.c_str());
     this->a1 = value;
     this->a2 = value;
     return;
  }

  void setProbY(String inpValue){
     double value = atof(inpValue.c_str());
     this->a1 = value;
     this->ProbY = value;
     return;
  }

  void setFilePtr(String inpFilePtr){
    this->filePtr = inpFilePtr;
    return;
  }

  void setMethod(String inpMethod){
    this->method = inpMethod;
    return;
  }

  void setJavaDoc(String inpJavaDoc){
    this->javadoc = inpJavaDoc;
    return;
  }

  void setReturnType(String inpRetType){
    this->returnType = inpRetType;
    return;
  }

  void setFormalParams(String inpFormalParamI){
      (*FormalParams).push_back(inpFormalParamI);
      return;
  }

  void setPostProb(double val){
	this->PostProb = val;  
  }
  
  void setBody(String inpBody){
    this->body = inpBody;
    return;
  }

  double getDoubleA1(){
    return this->a1;
  }

  double getDoubleA2(){
    return this->a2;
  }

  Handle<Vector<double>> getB1(){
    return this->b1;
  }

  Handle<Vector<double>> getB2(){
    return this->b2;
  }

  double getProbY(){
    return this->ProbY;
  }



  String getFilePtr(){
    return this->filePtr;
  }

  String getMethod(){
    return this->method;
  }

  String getJavaDoc(){
    return this->javadoc;
  }

  String getReturnType(){
    return this->returnType;
  }

  String getFormalParams(){
      std::string temp = "";
      for(int i=0;  i < (*FormalParams).size(); i++)
        temp = temp + (*FormalParams)[i].c_str() + " ";
      return temp;
  }

  String getBody(){
    return this->body;
  }
  // to get squared distance following SparkMLLib

  double getSquaredDistance(Handle<SearchProgramData> queryProg) {
      double a1_in = queryProg->getDoubleA1();
      Handle<Vector<double>> b1_in = queryProg->getB1();

      double a2 = this->getDoubleA2();
      Handle<Vector<double>> b2 = this->getB2();

      // t1 = tf.reduce_sum(tf.square(b1)/(4*a1), axis=0) + 0.5 * tf.reduce_sum(tf.log(-a1/np.pi), axis=0)
      double t1 = 0.5 * this->getSize() * log(-1 * a1_in / M_PI);
      for (int i=0; i<this->getSize(); i++){
          t1 += pow( (*b1_in)[i], 2) / (4 * a1_in);
      }

      // t2 = tf.reduce_sum(tf.square(b2)/(4*a2), axis=1) + 0.5 * tf.reduce_sum(tf.log(-a2/np.pi), axis=1)
      double t2 = 0.5 * this->getSize() * log(-1 * a2 / M_PI);
      for (int i=0; i<this->getSize(); i++){
          t2 += pow( (*b2)[i], 2) / (4 * a2);
      }
      // t3 = 0.5 * self.model.config.latent_size * tf.log(2*np.pi)
      double t3 = 0.5 * this->getSize() * log(2 * M_PI);

      double c = t1 + t2 - t3;
      //%%###############################################################%%//

      // b_star = b1 + b2
      Handle<Vector<double>> b_star = makeObject<Vector<double>>(this->getSize(), this->getSize());
      for(int i=0;i<this->getSize();i++){
          (*b_star)[i] = (*b1_in)[i] + (*b2)[i];
      }

      // a_star = a1 + a2 + 0.5
      double a_star = a1_in + a2 + 0.5;
      //c_star = tf.reduce_sum(tf.square(b_star)/(4*a_star), axis=1) + 0.5 * tf.reduce_sum(tf.log(-a_star/np.pi), axis=1)
      double c_star = 0.5 * this->getSize() * log(-1 * a_star / M_PI);
      for (int i=0; i<this->getSize(); i++){
          c_star += pow( (*b_star)[i], 2) / (4 * a_star);
      }


      double prob = (c - c_star) + this->getProbY();
      return prob;
  }

  void print(){
    std::cout << this->body << std::endl;
    return;
  }

  void fprint(FILE* fout, int rank){
	  
	fprintf(fout, "Prob :: %lf\n", this->PostProb);
        fprintf(fout, "Rank :: %d\n", rank);
	fprintf(fout, "File :: %s \n", this->filePtr.c_str());
	fprintf(fout, "Method :: %s \n ", this->method.c_str());
    fprintf(fout, "%s\n\n", this->body.c_str());
    return;
  }

  bool operator< (const SearchProgramData &other) const {
    return this->PostProb < other.PostProb;
  }

  ~SearchProgramData() {}

};
}
#endif
