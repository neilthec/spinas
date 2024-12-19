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

//File:  SPINAS/source/supleproduct.cpp

#include <iostream>
#include <sstream>
#include <cmath>

#include "utilities.h"
#include "cmatrix.h"
#include "cvector.h"
#include "particle.h"
#include "sproduct.h"
#include "schain.h"
#include "supleproduct.h"

namespace spinas {
  //Constructors
  
  //3 chains for a triple product
  supleproduct::supleproduct(schain *s0, schain *s1, schain *s2):
    N(3) {
    s[0] = s0;
    s[1] = s1;
    s[2] = s2;

    //Test that the chains all have the same Lorentz index exposed.
    check();
    //Update
    update();
  }

  //Check that the suple product makes sense
  void supleproduct::check(){
    isDotted = s[0]->is_Lindex_dotted();
    //Check that the chains all have the same Lorentz index exposed.
    for(int i=1;i<N;i++)
      if(s[i]->is_Lindex_dotted()!=isDotted)
        throw std::runtime_error("Incorrect usage of supleproduct.  All chains must have the same Lorentz index exposed.");
  
    //Determine which chains are massive
    for(int i=0;i<N;i++)
      isRightMassive[i] = s[i]->is_right_massive();
  }

  //Update
  //Must be run after masses or momenta of particles is updated.
  void supleproduct::update(){
    //Update isCalculated
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
        isCalculated[i][j] = false;
  } 


  //products
  cdouble supleproduct::v(){
    //Check that there are no spin indices.
    int nmassive=0;
    for(int i=0;i<N;i++)
      if(isRightMassive[i]) nmassive++;
    if(nmassive>0){
      throw std::runtime_error("Incorrect usage of supleproduct.v().  All the particles at the ends must be massless.");
      return cdouble(0,0);
    }

    //Check whether it is already calculted.
    if(isCalculated[0][0]) return product[0][0];

    //Calcualte it.
    cvector vec[3];
    for(int i=0;i<N;i++)
      vec[i] = s[i]->v();
    product[0][0] = v(vec);
    isCalculated[0][0] = true;

    return product[0][0];
  }

   cdouble supleproduct::v(const int& ds){
    //Check that there is 1 spin index.
    int nmassive=0, imassive=0;
    for(int i=0;i<N;i++)
      if(isRightMassive[i]){ 
        imassive=i;
        nmassive++;
      }
    if(nmassive>1){
      throw std::runtime_error("Incorrect usage of supleproduct.v(int j).  Only 1 particle at the end must be massive.");
      return cdouble(0,0);
    }
    
    //Get index from spin
    int jmassive=(ds+2)/2;
    //Check whether it is already calculted.
    if(isCalculated[imassive][jmassive]) return product[imassive][jmassive];

    //Calcualte it.
    cvector vec[3];
    for(int i=0;i<N;i++){
      if(i==imassive)
        vec[i] = s[i]->v(ds);
      else
        vec[i] = s[i]->v();
    }
    product[imassive][jmassive] = v(vec);
    isCalculated[imassive][jmassive] = true;

    return product[imassive][jmassive];
  }

  cdouble supleproduct::v(cvector vec[3]){
    if(N!=3)
      throw std::runtime_error("Only triple products are implemented.");
    return vec[0].get(0) * vec[1].get(1) * vec[2].get(2)
        + vec[0].get(1) * vec[1].get(2) * vec[2].get(0)
        + vec[0].get(2) * vec[1].get(0) * vec[2].get(1)
        - vec[0].get(2) * vec[1].get(1) * vec[2].get(0)
        - vec[0].get(0) * vec[1].get(2) * vec[2].get(1)
        - vec[0].get(1) * vec[1].get(0) * vec[2].get(2);
  }

}