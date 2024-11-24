
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
//#include "aliases.h"
#include "cmatrix.h"
#include "cvector.h"

namespace spinas {
  //Constructors
  cvector::cvector():
    vec{0,0,0},
    dimension(2){}

  cvector::cvector(const int& dim):
    vec{0,0,0},
    dimension(dim){}

  cvector::cvector(const cdouble& v0, const cdouble& v1):
    vec{v0,v1,0},
    dimension(2){}

  cvector::cvector(const cdouble& v0, const cdouble& v1, const cdouble& v2):
    vec{v0,v1,v2},
    dimension(3){}

  //Get Element
  cdouble cvector::get(const int& i) const{
    if (i < 0 || i >= dimension) 
      throw std::out_of_range("cvector index out of bounds");
    return vec[i];
  }

  //Conjugate
  cvector cvector::get_conjugate() const{
    cvector nvec(dimension);
    for(int i=0;i<dimension;i++)
      nvec.vec[i] = std::conj(vec[i]);
    return nvec;
    //return cvector(std::conj(vec[0]), std::conj(vec[1]));
  }



  //Algebra
  //*
  const cdouble cvector::operator*(const cvector& v2) const {
    if(dimension != v2.dimension)
      throw std::invalid_argument("cvector dimensions do not match");
    cdouble prod = 0;
    for(int i=0;i<dimension;i++)
      prod += vec[i]*v2.vec[i];
    return prod;
  }
  const cvector cvector::operator*(const cmatrix& m) const {
    if(dimension != m.get_dimension())
      throw std::invalid_argument("cvector and cmatrix dimensions do not match");
    cvector vnew(dimension);
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
	      vnew.vec[i] += vec[j]*m.get(j,i);
    return vnew;
  }
  cvector operator*(const cmatrix &m, const cvector& v2){
    if(m.get_dimension() != v2.dimension)
      throw std::invalid_argument("cvector and cmatrix dimensions do not match");
    cvector vnew(v2.dimension);
    for(int i=0;i<v2.dimension;i++)
      for(int j=0;j<v2.dimension;j++)
	      vnew.vec[i] += m.get(i,j)*v2.vec[j];
    return vnew;
  }
  const cvector cvector::operator*(const cdouble& d) const {
    cvector vnew(dimension);
    for(int i=0;i<dimension;i++)
      vnew.vec[i] = vec[i]*d;
    return vnew;
  }
  cvector operator*(const cdouble &d, const cvector &v2) {
    cvector vnew(v2.dimension);
    for(int i=0;i<v2.dimension;i++)
      vnew.vec[i] = d*v2.vec[i];
    return vnew;
  }
  const cvector cvector::operator*(const ldouble& d) const {
    cvector vnew(dimension);
    for(int i=0;i<dimension;i++)
      vnew.vec[i] = d*vec[i];
    return vnew;
  }
  cvector operator*(const ldouble& d, const cvector& v2) {
    cvector vnew(v2.dimension);
    for(int i=0;i<v2.dimension;i++)
      vnew.vec[i] = d*v2.vec[i];
    return vnew;
  }
  const cvector cvector::operator*(const int& d) const {
    ldouble dd = static_cast<ldouble>(d);
    cvector vnew(dimension);
    for(int i=0;i<dimension;i++)
      vnew.vec[i] = dd*vec[i];
    return vnew;
  }
  cvector operator*(const int& d, const cvector& v2) {
    ldouble dd = static_cast<ldouble>(d);
    cvector vnew(v2.dimension);
    for(int i=0;i<v2.dimension;i++)
      vnew.vec[i] = dd*v2.vec[i];
    return vnew;
  }
  //-
  cvector operator-(const cvector& v){
    cvector nvec(v.dimension);
    for(int i=0;i<v.dimension;i++)
      nvec.vec[i] = -v.vec[i];
    return nvec;
  }

  //Outer Product: For example |i>[i| = p_i
  cmatrix outer(const cvector& v1, const cvector& v2){
    if(v1.dimension != v2.dimension)
      throw std::invalid_argument("cvector dimensions do not match");
    cmatrix m(v1.dimension);
    for(int i=0;i<v1.dimension;i++)
      for(int j=0;j<v2.dimension;j++)
        m.set(i,j,v1.vec[i]*v2.vec[j]);
    return m;
  }



  // ==
  bool cvector::operator==(const cvector &v2) const {
    constexpr ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 10000000;
    if(dimension != v2.dimension) return false;
    for(int i=0;i<dimension;i++)
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
    tmpStr<<"{ ";
    if(which==1){
      tmpStr<<vec[0];
      for(int i=1;i<dimension;i++)
        tmpStr<<" , "<<vec[i];
    } 
    else{
      tmpStr<<std::abs(vec[0])<<" exp(i * "<<std::arg(vec[0])<<" )";
      for(int i=1;i<dimension;i++)
        tmpStr<<" , "<<std::abs(vec[i])<<" exp(i * "<<std::arg(vec[i])<<" )";
    }
    tmpStr<<" }";
    return tmpStr.str();
  }
  std::string cvector::to_string() const {
    std::stringstream tmpStr;
    tmpStr<<"{ "<<vec[0];
    for(int i=1;i<dimension;i++)
      tmpStr<<" , "<<vec[i];
    tmpStr<<" }";
    return tmpStr.str();
  }
  
  std::ostream& operator<<(std::ostream& o, const cvector & v){
    o << v.to_string();
    return o;
  }
  
  
  std::string cvector::to_Mathematica_string() const {
    std::stringstream tmpStr;
    tmpStr<<"{ "<<std::real(vec[0])<<" + I * "<<std::imag(vec[0]);
    for(int i=1;i<dimension;i++)
      tmpStr<<" , "<<std::real(vec[i])<<" + I * "<<std::imag(vec[i]);
    tmpStr<<" }";
    return tmpStr.str();
  }


}
