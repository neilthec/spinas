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

//File:  SPINAS/include/supleproduct.h

#pragma once

#include "types.h"



namespace spinas {
  
  class supleproduct{
  private:
    schain *s[3]; //The three chains for spin-1.  Needs to be increased for higher spin.
    int N;//Number of chains

    //Whether the Lorentz indices are dotted.
    bool isDotted= false;

    //Whether the ends are massive
    bool isRightMassive[3] = {false,false,false};

    //Whether the spinor uple product is calculated
    bool isCalculated[3][3] = {{false,false,false},{false,false,false},{false,false,false}};
    cdouble product[3][3];

    //Check that it makes sense
    void check();

    //products
    cdouble v(cvector vec[3]);
    

  public:
    //3 chains for a triple product
    supleproduct(schain *s0, schain *s1, schain *s2);

    //Update
    //Must be run after masses or momenta of particles is updated.
    void update();
    

    //products
    cdouble v();
    cdouble v(const int& ds);
    cdouble v(const int& ds1, const int& ds2);
    cdouble v(const int& ds1, const int& ds2, const int& ds3);

  };
  
  int test_supleproduct();
}