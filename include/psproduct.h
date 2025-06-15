
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

//File:  SPINAS/include/psproduct.h

#pragma once

#include "types.h"


namespace spinas {
  
  class psproduct{
  private:
    int dimension;//Whether spinors are for spin-1/2 (2) or spin-1 (3).//Must be spin-1 curly.
    particle *pR;//Curly spinor particle on the right.
    particle *pmu;//Particle whose 4 momentum is multiplying the curly spinor.

    //Whether the spin index is upper
    bool isRightUpper = true;

    //Whether pj |i}} is calculated
    bool isCalculated[3] = {false,false,false};
    cdouble product[3];
    

  public:
    psproduct();

    //0 Internal Momenta
    psproduct(particle* partmu,  particle* partR);
    psproduct(particle* partmu,  particle* partR, const bool& iRU);

    //Update
    //Must be run after masses or momenta of particles is updated.
    void update();
    

    //products
    cdouble v(const int& spin);

  };
  
  
}
