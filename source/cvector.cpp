
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

//File:  SPINAS/source/cvector.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>

#include "types.h"
#include "aliases.h"
#include "cmatrix.h"
#include "cvector.h"

namespace spinas {
  //Constructors
  cvector::cvector():
    vec({cdouble(),cdouble()}){}

  cvector::cvector(const cdouble& v0, const cdouble& v1):
    vec({v0,v1}){}

  cvector::cvector(const cdouble vnew[2]):
    vec({vnew[0],vnew[1]}){}


  //Get Element
  cdouble cvector::get(const int& i) const{
    if (i < 0 || i >= 2) 
      throw std::out_of_range("cvector index out of bounds");
    return vec[i];
  }

  //Conjugate
  cvector cvector::get_conjugate() const{
    return cvector(std::conj(vec[0]), std::conj(vec[1]));
  }


  //Algebra
  //*
  const cdouble cvector::operator*(const cvector& v2) const {
    return vec[0]*v2.vec[0]+vec[1]*v2.vec[1];
  }
  const cvector cvector::operator*(const cmatrix& m) const {
    constexpr ldouble zero = 0;
    cdouble vnew[2] = {zero, zero};
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	vnew[i] += vec[j]*m.get(j,i);
    return cvector(vnew);
  }
  cvector operator*(const cmatrix &m, const cvector& v2){
    constexpr ldouble zero = 0;
    cdouble vnew[2] = {zero, zero};
    for(int i=0;i<2;i++)
      for(int j=0;j<2;j++)
	vnew[i] += m.get(i,j)*v2.vec[j];
    return cvector(vnew);
  }
  const cvector cvector::operator*(const cdouble& d) const {
    return cvector(d*vec[0],d*vec[1]);
  }
  cvector operator*(const cdouble &d, const cvector &v2) {
    return cvector(d*v2.vec[0],d*v2.vec[1]);
  }
  const cvector cvector::operator*(const ldouble& d) const {
    return cvector(d*vec[0],d*vec[1]);
  }
  cvector operator*(const ldouble& d, const cvector& v2) {
    return cvector(d*v2.vec[0],d*v2.vec[1]);
  }
  const cvector cvector::operator*(const int& d) const {
    ldouble dd = static_cast<ldouble>(d);
    return cvector(dd*vec[0],dd*vec[1]);
  }
  cvector operator*(const int& d, const cvector& v2) {
    ldouble dd = static_cast<ldouble>(d);
    return cvector(dd*v2.vec[0],dd*v2.vec[1]);
  }
  //-
  cvector operator-(const cvector& v){
    return cvector(-v.vec[0],-v.vec[1]);
  }

  //Outer Product: For example |i>[i| = p_i
  cmatrix outer(const cvector& v1, const cvector& v2){
    return cmatrix(v1.vec[0]*v2.vec[0],v1.vec[0]*v2.vec[1],
		   v1.vec[1]*v2.vec[0],v1.vec[1]*v2.vec[1]);
  }



  // ==
  bool cvector::operator==(const cvector &v2) const {
    constexpr ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 10000000;
    for(int i=0;i<2;i++)
      if(std::abs(vec[i]-v2.vec[i]) > epsilon) return false;
    return true;
  }
  // !=
  bool cvector::operator!=(const cvector &v2) const { 
    return !(*this == v2);
  }



  //Screen Output
  std::string cvector::to_string(const int& which) const {
    std::stringstream tmpStr;
    if(which==1) tmpStr<<"{ "<<vec[0]<<" , "<<vec[1]<<" }";
    else tmpStr<<"{ "<<std::abs(vec[0])<<" exp(i * "<<std::arg(vec[0])<<" ) , "<<std::abs(vec[1])<<" exp(i * "<<std::arg(vec[1])<<" ) }";
    return tmpStr.str();
  }
  std::string cvector::to_string() const {
    std::stringstream tmpStr;
    tmpStr<<"{ "<<vec[0]<<" , "<<vec[1]<<" }";
    return tmpStr.str();
  }
  
  std::ostream& operator<<(std::ostream& o, const cvector & v){
    o << v.to_string();
    return o;
  }
  
  
  std::string cvector::to_Mathematica_string() const {
    std::stringstream tmpStr;
    tmpStr<<"{ "
	  <<std::real(vec[0])<<" + I * "<<std::imag(vec[0])<<
      " , "<<std::real(vec[1])<<" + I * "<<std::imag(vec[1])<<
      " }";
    return tmpStr.str();
  }


}
