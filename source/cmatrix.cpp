
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
    mat{{0,0,0},{0,0,0},{0,0,0}},
    dimension(2){}

  cmatrix::cmatrix(const int& dim):
    mat{{0,0,0},{0,0,0},{0,0,0}},
    dimension(dim){}

  cmatrix::cmatrix(const ldouble p[4], const bool& upp):
    cmatrix(p,upp,2){}

  
  cmatrix::cmatrix(const ldouble p[4], const bool& upp, const int& dim):
  dimension(dim)
  {
    if(dim==2){
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
    else if(dim==3){
      cdouble eppz = cdouble(p[0]+p[3],0), empz = cdouble(p[0]-p[3],0);
      cdouble pxppy = cdouble(p[1],p[2]), pxmpy = cdouble(p[1],-p[2]);
      cdouble sqrt2 = std::sqrt(2), two = 2.0;
      cdouble m2 = cdouble(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3],0);
      if(upp){//Upper Lorentz indices
        mat[0][0] = empz*empz;
        mat[0][1] = -sqrt2*empz*pxmpy;
        mat[0][2] = pxmpy*pxmpy;
        mat[1][0] = -sqrt2*empz*pxppy;
        mat[1][1] = m2+two*pxppy*pxmpy;
        mat[1][2] = -sqrt2*eppz*pxmpy;
        mat[2][0] = pxppy*pxppy;
        mat[2][1] = -sqrt2*eppz*pxppy;
        mat[2][2] = eppz*eppz;
      }
      else{//Lower Lorentz indices
        mat[0][0] = eppz*eppz;
        mat[0][1] = sqrt2*eppz*pxmpy;
        mat[0][2] = pxmpy*pxmpy;
        mat[1][0] = sqrt2*eppz*pxppy;
        mat[1][1] = m2+two*pxppy*pxmpy;
        mat[1][2] = sqrt2*empz*pxmpy;
        mat[2][0] = pxppy*pxppy;
        mat[2][1] = sqrt2*empz*pxppy;
        mat[2][2] = empz*empz;
      }
    }
    //We don't use this in practice.  This was just for an early test.
  }
  
  cmatrix::cmatrix(const cdouble& m00, const cdouble& m01, const cdouble& m10, const cdouble& m11):
    dimension(2),
    mat{{m00,m01},{m10,m11}}{}

  cmatrix::cmatrix(const cdouble& m00, const cdouble& m01, const cdouble& m02,
    const cdouble& m10, const cdouble& m11, const cdouble& m12,
    const cdouble& m20, const cdouble& m21, const cdouble& m22):
    dimension(3),
    mat{{m00,m01,m02},{m10,m11,m12},{m20,m21,m22}}{}


  //Get
  int cmatrix::get_dimension() const{
    return dimension;
  }
  cdouble cmatrix::get(const int& i, const int& j) const{
    if (i < 0 || i >= dimension || j < 0 || j >= dimension) 
      throw std::out_of_range("cmatrix index out of bounds");
    return mat[i][j];
  }
  
  //Set Elements
  void cmatrix::set(const int& i, const int& j, const cdouble& val){
    if (i < 0 || i >= dimension || j < 0 || j >= dimension) 
      throw std::out_of_range("cmatrix index out of bounds");
    mat[i][j] = val;
  }
  
  //Determinant
  ldouble cmatrix::get_det() const{
    ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
    cdouble det = 0;
    if(dimension==2)
      det = mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
    else if (dimension==3)
    {
      epsilon *= 1000000;
      det = mat[0][0]*mat[1][1]*mat[2][2]
          + mat[0][1]*mat[1][2]*mat[2][0]
          + mat[0][2]*mat[1][0]*mat[2][1]
          - mat[0][2]*mat[1][1]*mat[2][0]
          - mat[0][1]*mat[1][0]*mat[2][2]
          - mat[0][0]*mat[1][2]*mat[2][1];
    }
    
    if (std::abs(std::imag(det)) > epsilon)
      throw std::domain_error("Determinant has a non-zero imaginary part: "+std::to_string(std::imag(det))+".");
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
    if(dimension != m.dimension)
      throw std::invalid_argument("cmatrix dimensions do not match");
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
	mat[i][j]+=m.mat[i][j];
    return *this;		 
  }
  // -
  const cmatrix cmatrix::operator-(const cmatrix& m) const{
    return cmatrix(*this) -= m;
  }
  cmatrix & cmatrix::operator-=(const cmatrix& m){
    if(dimension != m.dimension)
      throw std::invalid_argument("cmatrix dimensions do not match");
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
	mat[i][j]-=m.mat[i][j];
    return *this;		 
  }
  cmatrix operator-(const cmatrix& m){
    cmatrix mnew(m.dimension);
    for(int i=0;i<m.dimension;i++)
      for(int j=0;j<m.dimension;j++)
        mnew.set(i,j,-m.mat[i][j]);
    return mnew;
  }
  // *
  cmatrix cmatrix::operator*(const cmatrix& m){
    return cmatrix(*this) *= m;
  }
  cmatrix & cmatrix::operator*=(const cmatrix& m){
    if(dimension != m.dimension)
      throw std::invalid_argument("cmatrix dimensions do not match");
    cdouble mnew[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
        for(int k=0;k<dimension;k++)
  	      mnew[i][j] += mat[i][k]*m.mat[k][j];
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
	      mat[i][j]=mnew[i][j];
    return *this;		 
  }
  cmatrix & cmatrix::operator*=(const cdouble& d){
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
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
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
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
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
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
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
	      mat[i][j]=mat[i][j]/d;
    return *this;
  }
  const cmatrix cmatrix::operator/(const cdouble &d) const{
    return cmatrix(*this) /= d;
  }
  cmatrix & cmatrix::operator/=(const ldouble &d){
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
	      mat[i][j]=mat[i][j]/d;
    return *this;
  }
  const cmatrix cmatrix::operator/(const ldouble &d) const{
    return cmatrix(*this) /= d;
  }
  cmatrix & cmatrix::operator/=(const int &d){
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
	      mat[i][j]=mat[i][j]/static_cast<ldouble>(d);
    return *this;
  }
  const cmatrix cmatrix::operator/(const int &d) const{
    return cmatrix(*this) /= d;
  }
  
  //Comparison
  //==
  bool cmatrix::operator==(const cmatrix &m) const{
    ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
    if(dimension==3) epsilon *= 1000;
    if(dimension != m.dimension) return false;
    for(int i=0;i<dimension;i++)
      for(int j=0;j<dimension;j++)
	      if(std::abs(mat[i][j]-m.mat[i][j]) > epsilon) return false;
    return true;
  }
  bool cmatrix::operator!=(const cmatrix &m) const{
    return !(*this == m);
  }
  
  
  
  
  
  
  //Screen Output
  std::string cmatrix::to_string() const {
    std::stringstream tmpStr;
    tmpStr<<"{ ";
    for(int i=0;i<dimension;i++){
      tmpStr<<"{ ";
      for(int j=0;j<dimension;j++){
        tmpStr<<std::real(mat[i][j])<<" + I* "<<std::imag(mat[i][j]);
        if(j<dimension-1) tmpStr<<" , ";
      }
      tmpStr<<" }";
      if(i<dimension-1) tmpStr<<" , ";
    }
    tmpStr<<" }";
    /*tmpStr<<"{ "<<std::real(mat[0][0])<<" + I* "<<std::imag(mat[0][0])<<
      " , "<<std::real(mat[0][1])<<" + I* "<<std::imag(mat[0][1])<<
      " } , { "<<std::real(mat[1][0])<<" + I* "<<std::imag(mat[1][0])<<
      " , "<<std::real(mat[1][1])<<" + I* "<<std::imag(mat[1][1])<<" } }";*/
    return tmpStr.str();
  }
  
  std::ostream& operator<<(std::ostream& o, const cmatrix & m){
    o << m.to_string();
    return o;
  }
  
  
  std::string cmatrix::to_Mathematica_string() const {
    std::stringstream tmpStr;
    tmpStr<<"{ ";
    for(int i=0;i<dimension;i++){
      tmpStr<<"{ ";
      for(int j=0;j<dimension;j++){
        tmpStr<<std::real(mat[i][j])<<" + I* "<<std::imag(mat[i][j]);
        if(j<dimension-1) tmpStr<<" , ";
      }
      tmpStr<<" }";
      if(i<dimension-1) tmpStr<<" , ";
    }
    tmpStr<<" }";
    /*tmpStr<<"{ { "<<std::real(mat[0][0])<<" + I* "<<std::imag(mat[0][0])<<
      " , "<<std::real(mat[0][1])<<" + I* "<<std::imag(mat[0][1])<<
      " } , { "<<std::real(mat[1][0])<<" + I* "<<std::imag(mat[1][0])<<
      " , "<<std::real(mat[1][1])<<" + I* "<<std::imag(mat[1][1])<<" } }";*/
    return tmpStr.str();
  }
  
  
  
  
}//end namespace spinas
