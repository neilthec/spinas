
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

//File:  SPINAS/include/schain.h

#pragma once

#include "types.h"


namespace spinas {
  
  class schain{
  private:
    int dimension;//Whether spinors are for spin-1/2 (2) or spin-1 (3).
    particle *pR;//Particles in the right spinors.
    particle *p[6];//Particles to the left.
    cmatrix pMat;//Product of all momenta in between the spinors
    int N;//Number of particles/momenta in product

    //Whether the end is massive
    bool isRightMassive = false;

    //Whether the right spinor is angle or square
    bool isRightAngle = false;
    
    //Whether the end is upper
    bool isRightUpper = true;

    //Whether the spinor product is calculated
    bool isCalculated[3] = {false,false,false};
    cvector product[3];
    

  public:
    schain();

    //0 Internal Momenta
    schain(particle* partR, const bool& as, const int& dim);
    schain(particle* partR, const bool& iRU, const bool& as, const int& dim);

    //1 Internal Momentum
    schain(particle* p0,  particle* partR, const bool& as, const int& dim);
    schain(particle* p0,  particle* partR, const bool& iRU, const bool& as, const int& dim);

    //2 Internal Momenta
    schain(particle* p0,  particle* p1,  particle* partR, const bool& as, const int& dim);
    schain(particle* p0,  particle* p1,  particle* partR, const bool& iRU, const bool& as, const int& dim);

    //3 Internal Momenta
    schain(particle* p0,  particle* p1,  particle* p2,  particle* partR, const bool& as, const int& dim);
    schain(particle* p0,  particle* p1,  particle* p2,  particle* partR, const bool& iRU, const bool& as, const int& dim);

    //4 Internal Momenta
    schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const bool& as, const int& dim);
    schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const bool& iRU, const bool& as, const int& dim);

    //5 Internal Momenta
    schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const bool& as, const int& dim);
    schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const bool& iRU, const bool& as, const int& dim);

    //6 Internal Momenta
    schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const bool& as, const int& dim);
    schain(particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const bool& iRU, const bool& as, const int& dim);

    //Update
    //Must be run after masses or momenta of particles is updated.
    void update();
    

    //products
    cvector v();
    cvector v(const int& jR);

  };
  
  int test_schain();
}
