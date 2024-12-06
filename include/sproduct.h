
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

//File:  SPINAS/include/sproduct.h

#pragma once

#include "types.h"


namespace spinas {
  
  class sproduct{
  private:
    int dimension;//Whether spinors are for spin-1/2 (2) or spin-1 (3).
    particle *pL, *pR;//Particles in the left and right spinors.
    particle *p[6];//Particles in the middle.
    cmatrix pMat;//Product of all momenta in between the spinors
    int N;//Number of particles/momenta in product

    //Whether the ends are massive
    bool isLeftMassive = false, isRightMassive = false;

    //Whether the left spinor is angle or square
    bool isLeftAngle = false;
    
    //Whether the ends are upper
    bool isLeftUpper = true, isRightUpper = true;

    //Whether the spinor product is calculated
    bool isCalculated[2][2] = {{false,false},{false,false}};
    cdouble product[2][2];
    

  public:
    sproduct();

    //0 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL,  particle* partR, const bool& iRU, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* partR, const bool& iRU, const int& dim);

    //1 Internal Momentum
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* partR, const bool& iRU, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* partR, const bool& iRU, const int& dim);

    //2 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* partR, const bool& iRU, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* partR, const bool& iRU, const int& dim);

    //3 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* partR, const bool& iRU, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* partR, const bool& iRU, const int& dim);

    //4 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const bool& iRU, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const bool& iRU, const int& dim);

    //5 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const bool& iRU, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const bool& iRU, const int& dim);

    //6 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const int& dim);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const bool& iRU, const int& dim);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const bool& iRU, const int& dim);

    //Update
    //Must be run after masses or momenta of particles is updated.
    void update();
    

    //products
    cdouble v();
    cdouble v(const int& spin);
    cdouble v(const int& jL, const int& jR);

  };
  
  int test_sproduct_sub(sproduct* sp1, sproduct* sp2, const int& ni, const int& np, particle* p1, particle* p2, particle* p3, particle* p4, particle* p5, const ldouble& expected, const char* spstring, const char* resstring);
  int test_sproduct();
}
