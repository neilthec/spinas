
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

#define SPINAS_HAS_MOMENTUM_CHAIN 1

#include "types.h"
#include <initializer_list>


namespace spinas {

  // One explicitly updated middle-momentum product, shared by brackets with
  // different endpoints. The pointed-to particles must outlive this object.
  // There is no global cache; each process or worker owns its own chains.
  class momentum_chain {
    friend class sproduct;
  private:
    particle* particles[6] = {};
    int length = 0;
    bool leftAngle = false;
    bool ready = false;
    cmatrix product;
  public:
    momentum_chain() = default;
    momentum_chain(const bool& angle, std::initializer_list<particle*> momenta);
    // Call after any inserted particle momentum or mass changes.
    void update();
    const cmatrix& matrix() const;
  };

  class sproduct{
  private:
    particle *pL = nullptr, *pR = nullptr;//Particles in the left and right spinors.
    particle *p[6] = {};//Particles in the middle.
    cmatrix pMat;//Product of all momenta in between the spinors
    int N = 0;//Number of particles/momenta in product

    //Whether the ends are massive
    bool isLeftMassive = false, isRightMassive = false;

    //Whether the left spinor is angle or square
    bool isLeftAngle = false;
    
    //Whether the ends are upper
    bool isLeftUpper = true, isRightUpper = true;

    //Whether the spinor product is calculated
    bool isCalculated[2][2] = {{false,false},{false,false}};
    cdouble product[2][2];
    void update_endpoints();
    

  public:
    sproduct();

    //0 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* partR);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* partR);
    sproduct(const bool& as,  particle* partL,  particle* partR, const bool& iRU);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* partR, const bool& iRU);

    //1 Internal Momentum
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* partR);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* partR);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* partR, const bool& iRU);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* partR, const bool& iRU);

    //2 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* partR);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* partR);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* partR, const bool& iRU);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* partR, const bool& iRU);

    //3 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* partR);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* partR);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* partR, const bool& iRU);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* partR, const bool& iRU);

    //4 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const bool& iRU);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* partR, const bool& iRU);

    //5 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const bool& iRU);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4,  particle* partR, const bool& iRU);

    //6 Internal Momenta
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR);
    sproduct(const bool& as,  particle* partL,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const bool& iRU);
    sproduct(const bool& as,  particle* partL, const bool& iLU,  particle* p0,  particle* p1,  particle* p2,  particle* p3,  particle* p4, particle* p5,  particle* partR, const bool& iRU);

    //Update
    //Must be run after masses or momenta of particles is updated.
    void update();
    // Use a ready shared chain with the exact same inserted particles and
    // angle/square orientation. Endpoint spin caches are invalidated normally.
    void update(const momentum_chain& chain);

    //products
    cdouble v();
    cdouble v(const int& spin);
    cdouble v(const int& jL, const int& jR);

  };
  
  int test_sproduct_sub(sproduct* sp1, sproduct* sp2, const int& ni, const int& np, particle* p1, particle* p2, particle* p3, particle* p4, particle* p5, const ldouble& expected, const char* spstring, const char* resstring);
  int test_sproduct();
}
