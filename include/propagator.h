
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

//File:  SPINAS/include/propagator.h

#pragma once

#include "types.h"

namespace spinas {

  class propagator{
  private:
    ldouble m;//Mass
    ldouble w;//Width

    
  public:
    propagator();
    propagator(const ldouble& mass, const ldouble& width);

    void set_mass(const ldouble& mass);
    void set_width(const ldouble& width);

    //Denominator
    cdouble den(const ldouble p[4]);

  };

  //Tests
  int test_propagator();
}
