
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2025 Neil Christensen

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

//File:  SPINAS/source/psproduct.cpp

#include <iostream>
#include <sstream>
#include <cmath>

#include "utilities.h"
#include "cmatrix.h"
#include "cvector.h"
#include "particle.h"
#include "psproduct.h"


namespace spinas {
  //Constructors
  psproduct::psproduct(){}

  //0 Internal Momenta
  psproduct::psproduct(particle* partmu,  particle* partR):
    pmu(partmu), pR(partR), isRightUpper(true), dimension(3) {
    update();
  }
  psproduct::psproduct(particle* partmu,  particle* partR, const bool& iRU):
    pmu(partmu), pR(partR), isRightUpper(iRU), dimension(3) {
    update();
  }



  //If masses or momenta change, this function must be run.
  void psproduct::update(){
    if(pR->mass()==0) 
      throw std::runtime_error("Incorrect usage of psproduct.  The right particle must be massive.");
    //Update isCalculated
    for(int i=0;i<dimension;i++)
	    isCalculated[i] = false;
  }


  //Product
  cdouble psproduct::v(const int& spin){
    int jR = (spin+2)/2;
    //Check whether it is already calculted.
    if(isCalculated[jR]) return product[jR];
    //Calcualte it.
    product[jR] = pmu->pmu() * pR->rcurly(spin,isRightUpper,dimension);
	  isCalculated[jR] = true;
    return product[jR];
  }

 
  
}
