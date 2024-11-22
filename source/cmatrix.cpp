
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

//File:  SPINAS/source/cmatrix.cpp

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>
#include <string> 

#include "types.h"
//#include "aliases.h"
#include "cmatrix.h"

namespace spinas {
  //Constructors
  cmatrix::cmatrix():
    mat{{cdouble(),cdouble()},{cdouble(),cdouble()}},
    sizeN(2){}
  
  cmatrix::cmatrix(const ldouble p[4], const bool& upp):
  sizeN(2)
  {
    if(upp){//Upper Lorentz indices
      mat[0][0] = cdouble(p[0]-p[3],0);
      mat[0][1] = cdouble(-p[1],p[2]);
      mat[1][0] = cdouble(-p[1],-p[2]);
      mat[1][1] = cdouble(p[0]+p[3],0);
    }
    else{//Lower Lorentz indices
      mat[0][0] = cdouble(p[0]+p[3],0);
      mat[0][1] = cdouble(p[1],-p[2]);
      mat[1][0] = cdouble(p[1],p[2]);
      mat[1][1] = cdouble(p[0]-p[3],0);
    }
  }
  
  cmatrix::cmatrix(const cdouble& m00, const cdouble& m01, const cdouble& m10, const cdouble& m11):
    sizeN(2),
    mat{{m00,m01},{m10,m11}}{}


  //Get Elements
  cdouble cmatrix::get(const int& i, const int& j) const{
    if (i < 0 || i >= sizeN || j < 0 || j >= sizeN) 
      throw std::out_of_range("cmatrix index out of bounds");
    return mat[i][j];
  }
  
  
  
  //Determinant
  ldouble cmatrix::get_det() const{
    constexpr ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
    cdouble det = mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
    if (std::abs(std::imag(det)) > epsilon)
      throw std::domain_error("Determinant has a non-zero imaginary part.");
    return std::real(det);
  }
  ldouble det(const cmatrix& m){
    return m.get_det();
  }
  
  
  //Algebra
  // +
  const cmatrix cmatrix::operator+(const cmatrix& m) const{
    return cmatrix(*this) += m;
  }
  cmatrix & cmatrix::operator+=(const cmatrix& m){
    for(int i=0;i<sizeN;i++)
      for(int j=0;j<sizeN;j++)
	mat[i][j]+=m.mat[i][j];
    return *this;		 
  }
  // -
  const cmatrix cmatrix::operator-(const cmatrix& m) const{
    return cmatrix(*this) -= m;
  }
  cmatrix & cmatrix::operator-=(const cmatrix& m){
    for(int i=0;i<sizeN;i++)
      for(int j=0;j<sizeN;j++)
	mat[i][j]-=m.mat[i][j];
    return *this;		 
  }
  cmatrix operator-(const cmatrix& m){
    return cmatrix(-m.mat[0][0],-m.mat[0][1],-m.mat[1][0],-m.mat[1][1]);
  }
  // *
  cmatrix cmatrix::operator*(const cmatrix& m){
    return cmatrix(*this) *= m;
  }
  cmatrix & cmatrix::operator*=(const cmatrix& m){
    cdouble mnew[3][3];
    for(int i=0;i<sizeN;i++)
      for(int j=0;j<sizeN;j++)
	mnew[i][j] = mat[i][0]*m.mat[0][j]+mat[i][1]*m.mat[1][j];
    for(int i=0;i<sizeN;i++)
      for(int j=0;j<sizeN;j++)
	mat[i][j]=mnew[i][j];
    return *this;		 
  }
  cmatrix & cmatrix::operator*=(const cdouble& d){
    for(int i=0;i<sizeN;i++)
      for(int j=0;j<sizeN;j++)
	mat[i][j]=mat[i][j]*d;
    return *this;
  }
  const cmatrix cmatrix::operator*(const cdouble& d) const {
    return cmatrix(*this) *= d;
  }
  cmatrix operator*(const cdouble &d, const cmatrix &m){
    return m*d;
  }
  cmatrix & cmatrix::operator*=(const ldouble& d){
    for(int i=0;i<sizeN;i++)
      for(int j=0;j<sizeN;j++)
	mat[i][j]=mat[i][j]*d;
    return *this;
  }
  const cmatrix cmatrix::operator*(const ldouble& d) const {
    return cmatrix(*this) *= d;
  }
  cmatrix operator*(const ldouble &d, const cmatrix &m){
    return m*d;
  }
  cmatrix & cmatrix::operator*=(const int& d){
    for(int i=0;i<sizeN;i++)
      for(int j=0;j<sizeN;j++)
	mat[i][j]=mat[i][j] * static_cast<ldouble>(d);
    return *this;
  }
  const cmatrix cmatrix::operator*(const int& d) const {
    return cmatrix(*this) *= d;
  }
  cmatrix operator*(const int &d, const cmatrix &m){
    return m*d;
  }
  
  // /
  cmatrix & cmatrix::operator/=(const cdouble &d){
    for(int i=0;i<sizeN;i++)
      for(int j=0;j<sizeN;j++)
	mat[i][j]=mat[i][j]/d;
    return *this;
  }
  const cmatrix cmatrix::operator/(const cdouble &d) const{
    return cmatrix(*this) /= d;
  }
  cmatrix & cmatrix::operator/=(const ldouble &d){
    for(int i=0;i<sizeN;i++)
      for(int j=0;j<sizeN;j++)
	mat[i][j]=mat[i][j]/d;
    return *this;
  }
  const cmatrix cmatrix::operator/(const ldouble &d) const{
    return cmatrix(*this) /= d;
  }
  cmatrix & cmatrix::operator/=(const int &d){
    for(int i=0;i<sizeN;i++)
      for(int j=0;j<sizeN;j++)
	mat[i][j]=mat[i][j]/static_cast<ldouble>(d);
    return *this;
  }
  const cmatrix cmatrix::operator/(const int &d) const{
    return cmatrix(*this) /= d;
  }
  
  //Comparison
  //==
  bool cmatrix::operator==(const cmatrix &m) const{
    constexpr ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
    for(int i=0;i<sizeN;i++)
      for(int j=0;j<sizeN;j++)
	if(std::abs(mat[i][j]-m.mat[i][j]) > epsilon) return false;
    return true;
  }
  bool cmatrix::operator!=(const cmatrix &m) const{
    return !(*this == m);
  }
  
  
  
  
  
  
  //Screen Output
  std::string cmatrix::to_string() const {
    std::stringstream tmpStr;
    tmpStr<<"{ { "<<std::real(mat[0][0])<<" + I* "<<std::imag(mat[0][0])<<
      " , "<<std::real(mat[0][1])<<" + I* "<<std::imag(mat[0][1])<<
      " } , { "<<std::real(mat[1][0])<<" + I* "<<std::imag(mat[1][0])<<
      " , "<<std::real(mat[1][1])<<" + I* "<<std::imag(mat[1][1])<<" } }";
    return tmpStr.str();
  }
  
  std::ostream& operator<<(std::ostream& o, const cmatrix & m){
    o << m.to_string();
    return o;
  }
  
  
  std::string cmatrix::to_Mathematica_string() const {
    std::stringstream tmpStr;
    tmpStr<<"{ { "<<std::real(mat[0][0])<<" + I* "<<std::imag(mat[0][0])<<
      " , "<<std::real(mat[0][1])<<" + I* "<<std::imag(mat[0][1])<<
      " } , { "<<std::real(mat[1][0])<<" + I* "<<std::imag(mat[1][0])<<
      " , "<<std::real(mat[1][1])<<" + I* "<<std::imag(mat[1][1])<<" } }";
    return tmpStr.str();
  }
  
  
  
  
}//end namespace spinas
