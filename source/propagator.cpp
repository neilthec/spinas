
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

//File:  SPINAS/source/propagator.cpp

#include <iostream>
#include <sstream>
#include <cmath>

#include "propagator.h"

namespace spinas {
  //Constructors
  propagator::propagator():
    m(0), w(0) {}

  propagator::propagator(const ldouble& mass, const ldouble& width):
    m(mass), w(width) {}

  void propagator::set_mass(const ldouble& mass){
    m = mass;
  }
  void propagator::set_width(const ldouble& width){
    w = width;
  }


  //Denominator
  cdouble propagator::den(const ldouble p[4]){
    return cdouble(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3] - m*m, -m*w);
  }


}
