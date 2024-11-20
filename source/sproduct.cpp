
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2023 Neil Christensen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

//File:  SPINAS/source/sproduct.cpp

#include <iostream>
#include <sstream>
#include <cmath>

#include "utilities.h"
#include "cmatrix.h"
#include "cvector.h"
#include "particle.h"
#include "sproduct.h"


namespace spinas {
  //Constructors
  sproduct::sproduct(){}

  //0 Internal Momenta
  sproduct::sproduct(const bool& as,  particle* partL,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(true), N(0) {
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(true), N(0) {
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(iRU), N(0) {
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(iRU), N(0) {
    update();
  }
  
  //1 Internal Momenta
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(true), N(1) {
    p[0] = p0;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(true), N(1) {
    p[0] = p0;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(iRU), N(1) {
    p[0] = p0;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(iRU), N(1) {
    p[0] = p0;
    update();
  }

  //2 Internal Momenta
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(true), N(2) {
    p[0] = p0;
    p[1] = p1;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(true), N(2) {
    p[0] = p0;
    p[1] = p1;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(iRU), N(2) {
    p[0] = p0;
    p[1] = p1;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(iRU), N(2) {
    p[0] = p0;
    p[1] = p1;
    update();
  }

  //3 Internal Momenta
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(true), N(3) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(true), N(3) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(iRU), N(3) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(iRU), N(3) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    update();
  }

  //4 Internal Momenta
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(true), N(4) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(true), N(4) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(iRU), N(4) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(iRU), N(4) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    update();
  }

  //5 Internal Momenta
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(true), N(5) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(true), N(5) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(iRU), N(5) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(iRU), N(5) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    update();
  }


  //6 Internal Momenta
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(true), N(6) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    p[5] = p5;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(true), N(6) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    p[5] = p5;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(true), pR(partR), isRightUpper(iRU), N(6) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    p[5] = p5;
    update();
  }
  sproduct::sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const bool& iRU):
    isLeftAngle(as), pL(partL), isLeftUpper(iLU), pR(partR), isRightUpper(iRU), N(6) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    p[5] = p5;
    update();
  }



  //If masses or momenta change, this function must be run.
  void sproduct::update(){
    isLeftMassive = false;
    isRightMassive = false;
    if(pL->mass()!=0) isLeftMassive = true;
    if(pR->mass()!=0) isRightMassive = true;
    //Check whether massless spinors have lower indices
    if(!isLeftMassive && !isLeftUpper)
      throw std::runtime_error("Incorrect usage of sproduct.  The left particle can only be specified to have a lower index if it is massive.");
    if(!isRightMassive && !isRightUpper)
      throw std::runtime_error("Incorrect usage of sproduct.  The right particle can only be specified to have a lower index if it is massive.");
    //Update isCalculated
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	isCalculated[i][j] = false;
    //Calculate pMat as the product of the momenta in the middle.
    if(N>0){
      if(isLeftAngle){
	pMat = p[0]->lmat();
	for(int i=1;i<N;i++){
	  if(i%2==0) pMat *= p[i]->lmat();
	  else pMat *= p[i]->umat();
	}
      }
      else {
	pMat = p[0]->umat();
	for(int i=1;i<N;i++){
	  if(i%2==0) pMat *= p[i]->umat();
	  else pMat *= p[i]->lmat();
	}
      }
    }
  }


  //Products
  //Both massless
  cdouble sproduct::v(){
    if(isLeftMassive || isRightMassive){
      throw std::runtime_error("Incorrect usage of sproduct.val().  Both particles must be massless.");
      return cdouble(0,0);
    }
    //Check whether it is already calculted.
    if(isCalculated[0][0]) return product[0][0];
    //Calcualte it.
    cvector vec;
    if(isLeftAngle){
      vec = pL->langle();
      if(N>0)
	vec = vec * pMat;
      if(N%2==0)
	product[0][0] = vec * pR->rangle();
      else
	product[0][0] = vec * pR->rsquare();
    }
    else{
      vec = pL->lsquare();
      if(N>0)
	vec = vec * pMat;
      if(N%2==0)
	product[0][0] = vec * pR->rsquare();
      else
	product[0][0] = vec * pR->rangle();
    }
    isCalculated[0][0] = true;
    return product[0][0];
  }
  //One massive spinor
  cdouble sproduct::v(const int& spin){
    if(isLeftMassive){
      if(isRightMassive){
	throw std::runtime_error("Incorrect usage of sproduct.val(j).  Only one particle must be massive.");
	return cdouble(0,0);
      }
      int jL=0;
      if(spin>0) jL=1;
      //Check whether it is already calculted.
      if(isCalculated[jL][0]) return product[jL][0];
      //Calcualte it.
      cvector vec;
      if(isLeftAngle){
	vec = pL->langle(spin,isLeftUpper);
	if(N>0)
	  vec = vec * pMat;
	if(N%2==0)
	  product[jL][0] = vec * pR->rangle();
	else
	  product[jL][0] = vec * pR->rsquare();
      }
      else{
	vec = pL->lsquare(spin,isLeftUpper);
	if(N>0)
	  vec = vec * pMat;
	if(N%2==0)
	  product[jL][0] = vec * pR->rsquare();
	else
	  product[jL][0] = vec * pR->rangle();
      }
      isCalculated[jL][0] = true;
      return product[jL][0];
    }
    else if(isRightMassive){
      if(isLeftMassive){
	throw std::runtime_error("Incorrect usage of sproduct.val(j).  Only one particle must be massive.");
	return cdouble(0,0);
      }
      int jR=0;
      if(spin>0) jR=1;
      //Check whether it is already calculted.
      if(isCalculated[0][jR]) return product[0][jR];
      //Calcualte it.
      cvector vec;
      if(isLeftAngle){
	vec = pL->langle();
	if(N>0)
	  vec = vec * pMat;
	if(N%2==0)
	  product[0][jR] = vec * pR->rangle(spin,isRightUpper);
	else
	  product[0][jR] = vec * pR->rsquare(spin,isRightUpper);
      }
      else{
	vec = pL->lsquare();
	if(N>0)
	  vec = vec * pMat;
	if(N%2==0)
	  product[0][jR] = vec * pR->rsquare(spin,isRightUpper);
	else
	  product[0][jR] = vec * pR->rangle(spin,isRightUpper);
      }
      isCalculated[0][jR] = true;
      return product[0][jR];
    }
    throw std::runtime_error("Incorrect usage of sproduct.val(j).  One particle must be massive.");
    return cdouble(0,0);
  }
  //Both massive
  cdouble sproduct::v(const int& spinL, const int& spinR){
    if(!isLeftMassive || !isRightMassive){
      throw std::runtime_error("Incorrect usage of sproduct.val(jL,jR).  Both particles must be massive.");
      return cdouble(0,0);
    }
    int jL=0;
    if(spinL>0) jL=1;
    int jR=0;
    if(spinR>0) jR=1;
    //Check whether it is already calculted.
    if(isCalculated[jL][jR]) return product[jL][jR];
    //Calcualte it.
    cvector vec;
    if(isLeftAngle){
      vec = pL->langle(spinL,isLeftUpper);
      if(N>0)
	vec = vec * pMat;
      if(N%2==0)
	product[jL][jR] = vec * pR->rangle(spinR,isRightUpper);
      else
	product[jL][jR] = vec * pR->rsquare(spinR,isRightUpper);
    }
    else{
      vec = pL->lsquare(spinL,isLeftUpper);
      if(N>0)
	vec = vec * pMat;
      if(N%2==0)
	product[jL][jR] = vec * pR->rsquare(spinR,isRightUpper);
      else
	product[jL][jR] = vec * pR->rangle(spinR,isRightUpper);
    }
    isCalculated[jL][jR] = true;
    return product[jL][jR];
  }


 
  
}
