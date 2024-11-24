
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

//File:  SPINAS/include/cmatrix.h

#pragma once

#include "types.h"


namespace spinas{
  
  class cmatrix{
  private:
    cdouble mat[3][3];
    int dimension;
    
  public:
    cmatrix();
    explicit cmatrix(const int& dim);
    cmatrix(const ldouble p[4], const bool& upp);
    cmatrix(const cdouble& m00, const cdouble& m01, const cdouble& m10, const cdouble& m11);
    
    //Get
    int get_dimension() const;
    cdouble get(const int& i, const int& j) const;

    //Set Elements
    void set(const int& i, const int& j, const cdouble& val);
    
    //Determinant
    ldouble get_det() const;
    friend ldouble det(const cmatrix& z);
    
    
    //Algebra
    // +
    const cmatrix operator+(const cmatrix& m) const;
    cmatrix & operator+=(const cmatrix& m);
    // -
    const cmatrix operator-(const cmatrix& m) const;
    cmatrix & operator-=(const cmatrix& m);
    friend cmatrix operator-(const cmatrix& m);
    // *
    cmatrix operator*(const cmatrix& m);
    cmatrix & operator*=(const cmatrix& m);
    cmatrix & operator*=(const cdouble& d);
    const cmatrix operator*(const cdouble& d) const;
    friend cmatrix operator*(const cdouble &d, const cmatrix &m);
    cmatrix & operator*=(const ldouble& d);
    const cmatrix operator*(const ldouble& d) const;
    friend cmatrix operator*(const ldouble &d, const cmatrix &m);
    cmatrix & operator*=(const int& d);
    const cmatrix operator*(const int& d) const;
    friend cmatrix operator*(const int &d, const cmatrix &m);
    // /
    cmatrix & operator/=(const cdouble &d);
    const cmatrix operator/(const cdouble &d) const;
    cmatrix & operator/=(const ldouble &d);
    const cmatrix operator/(const ldouble &d) const;
    cmatrix & operator/=(const int &d);
    const cmatrix operator/(const int &d) const;
    
    //Comparison
    bool operator==(const cmatrix &m) const;
    bool operator!=(const cmatrix &m) const;
    
    
    
    //Output
    std::string to_string() const;
    friend std::ostream& operator<<(std::ostream& o, const cmatrix& m);
    std::string to_Mathematica_string() const;
    
  };//end class cmatrix
  

  int test_cmatrix();
  
}//end namespace spinas
