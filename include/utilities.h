
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

//File:  SPINAS/include/utilities.h

#pragma once

#include <string>
#include "types.h"

namespace spinas {
  
  //Double to Mathematica String
  std::string double_to_Mathematica_string(const ldouble& x);

  //Rotate a Momentum by 3 angles
  void rotate_momentum(ldouble mom[4], const ldouble u[3], const ldouble& theta);
  //Boost a momentum by v[3]={vx,vy,vz};
  void boost_momentum(ldouble mom[4], const ldouble v[3]);

  //Random Number Generation
  ldouble choose_random_momentum(ldouble p[4], ldouble begin, ldouble end); 
  void choose_random_massless_momentum(ldouble p[4], ldouble begin, ldouble end);  
  cdouble choose_random_cdouble(ldouble begin, ldouble end);  
  ldouble choose_random_ldouble(ldouble begin, ldouble end);
  int choose_random_int(int begin, int end);
}

