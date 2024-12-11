
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2024 Neil Christensen

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

//File:  SPINAS/source/schain.cpp

#include <iostream>
#include <sstream>
#include <cmath>

#include "utilities.h"
#include "cmatrix.h"
#include "cvector.h"
#include "particle.h"
#include "schain.h"


namespace spinas {

  //Constructors
  schain::schain(){}

  //0 Internal Momenta
  schain::schain(particle* partR, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(true), N(0), dimension(dim) {
    update();
  }
  schain::schain(particle* partR, const bool& iRU, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(iRU), N(0), dimension(dim) {
    update();
  }
  
  //1 Internal Momenta
  schain::schain(particle* p0,  particle* partR, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(true), N(1), dimension(dim) {
    p[0] = p0;
    update();
  }
  schain::schain(particle* p0,  particle* partR, const bool& iRU, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(iRU), N(1), dimension(dim) {
    p[0] = p0;
    update();
  }

  //2 Internal Momenta
  schain::schain(particle* p0,  particle* p1,  particle* partR, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(true), N(2), dimension(dim) {
    p[0] = p0;
    p[1] = p1;
    update();
  }
  schain::schain(particle* p0,  particle* p1,  particle* partR, const bool& iRU, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(iRU), N(2), dimension(dim) {
    p[0] = p0;
    p[1] = p1;
    update();
  }

  //3 Internal Momenta
  schain::schain(particle* p0,  particle* p1,  particle* p2,  particle* partR, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(true), N(3), dimension(dim) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    update();
  }
  schain::schain(particle* p0,  particle* p1,  particle* p2,  particle* partR, const bool& iRU, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(iRU), N(3), dimension(dim) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    update();
  }

  //4 Internal Momenta
  schain::schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(true), N(4), dimension(dim) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    update();
  }
  schain::schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const bool& iRU, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(iRU), N(4), dimension(dim) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    update();
  }

  //5 Internal Momenta
  schain::schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(true), N(5), dimension(dim) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    update();
  }
  schain::schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const bool& iRU, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(iRU), N(5), dimension(dim) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    update();
  }


  //6 Internal Momenta
  schain::schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(true), N(6), dimension(dim) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    p[5] = p5;
    update();
  }
  schain::schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const bool& iRU, const bool& as, const int& dim):
    isRightAngle(as), pR(partR), isRightUpper(iRU), N(6), dimension(dim) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    p[4] = p4;
    p[5] = p5;
    update();
  }


  //If masses or momenta change, this function must be run.
  void schain::update(){
    isRightMassive = false;
    if(pR->mass()!=0) isRightMassive = true;
    //Check whether massless spinors have lower indices
    if(!isRightMassive && !isRightUpper)
        throw std::runtime_error("Incorrect usage of schain.  The right particle can only be specified to have a lower index if it is massive.");
    //Update isCalculated
    for(int i=0;i<dimension;i++)
	    isCalculated[i] = false;
    //Calculate pMat as the product of the momenta in the middle.
    if(N>0){
      if((isRightAngle&&N%2==0)||(!isRightAngle&&N%2==1)){
    	pMat = p[0]->lmat(dimension);
	    for(int i=1;i<N;i++){
	        if(i%2==0) pMat *= p[i]->lmat(dimension);
	        else pMat *= p[i]->umat(dimension);
	    }
      }
      else {
	    pMat = p[0]->umat(dimension);
	    for(int i=1;i<N;i++){
	        if(i%2==0) pMat *= p[i]->umat(dimension);
	        else pMat *= p[i]->lmat(dimension);
	    }
      }
    }
  }


  //Products
  //Massless
  cvector schain::v(){
    if(isRightMassive){
      throw std::runtime_error("Incorrect usage of schain.v().  The particle at the right must be massless.");
      return cvector(0,0);
    }
    //Check whether it is already calculted.
    if(isCalculated[0]) return product[0];
    //Calcualte it.
    if(isRightAngle){
      if(N>0)
	    product[0] = pMat * pR->rangle(dimension);
      else
	    product[0] = pR->rangle(dimension);
    }
    else{
      if(N>0)
	    product[0] = pMat * pR->rsquare(dimension);
      else
	    product[0] = pR->rsquare(dimension);
    }
    isCalculated[0] = true;
    return product[0];
  }

  //Massive spinor
  cvector schain::v(const int& spin){
    if(!isRightMassive){
      throw std::runtime_error("Incorrect usage of schain.v(j).  The particle must be massive.");
	  return cvector(0,0);
    }
    int jR=0;
    if(dimension==2)
        jR = (spin+1)/2;
    else if(dimension==3)
        jR = (spin+2)/2;
    //Check whether it is already calculted.
    if(isCalculated[jR]) return product[jR];
    //Calcualte it.
    if(isRightAngle){
	    if(N>0)
	        product[jR] = pMat * pR->rangle(spin,isRightUpper,dimension);
	    else
	        product[jR] = pR->rangle(spin,isRightUpper,dimension);
    }
    else{
	    if(N>0)
	        product[jR] = pMat * pR->rsquare(spin,isRightUpper,dimension);
	    else
	        product[jR] = pR->rsquare(spin,isRightUpper,dimension);
    }
    isCalculated[jR] = true;
    return product[jR];
    throw std::runtime_error("Incorrect usage of schain.v(j).  Particle must be massive.");
    return cvector(0,0);
  }





}