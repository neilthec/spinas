
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

//File:  SPINAS/source/cvector.h

#pragma once

#include "types.h"

namespace spinas{
  
  class cvector{
  private:
    cdouble vec[3];
    int sizeN;
    
  public:
    cvector();
    cvector(const cdouble& v0, const cdouble& v1);
    
    //Get Element
    cdouble get(const int& i) const;
    
    //Conjugate
    cvector get_conjugate() const;
    
    //Algebra
    //*
    const cdouble operator*(const cvector& v2) const;
    const cvector operator*(const cmatrix& m) const;
    friend cvector operator*(const cmatrix &m, const cvector& vec2);
    const cvector operator*(const cdouble& d) const;
    friend cvector operator*(const cdouble &d, const cvector &v2);
    const cvector operator*(const ldouble& d) const;
    friend cvector operator*(const ldouble& d, const cvector& v2);
    const cvector operator*(const int& d) const;
    friend cvector operator*(const int& d, const cvector& v2);
    // -
    friend cvector operator-(const cvector &v);
    
    //Outer Product: For example |i>[i| = p_i
    friend cmatrix outer(const cvector& v1, const cvector& v2);

    // ==
    bool operator==(const cvector &v2) const;
    // !=
    bool operator!=(const cvector &v2) const;

    //Screen Output
    std::string to_string(const int& which) const;
    std::string to_string() const;
    friend std::ostream& operator<<(std::ostream& o, const cvector & v);
    std::string to_Mathematica_string() const;
    
  };


  int test_cvector();
}

