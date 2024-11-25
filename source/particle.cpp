
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

//File:  SPINAS/source/particle.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>

#include "types.h"
//#include "aliases.h"
#include "cmatrix.h"
#include "cvector.h"
#include "particle.h"

namespace spinas {
  //Constructors
  particle::particle()
  {}

  particle::particle(const ldouble& mass):
    m(mass), p{m,0,0,0}{}

  particle::particle(const ldouble momentum[4], const ldouble& mass):
    m(mass),p{momentum[0],momentum[1],momentum[2],momentum[3]}
  {
    update();
  }

  //Update momentum and mass
  void particle::set_mass(const ldouble& mass){
    m=mass;
  }
  void particle::set_momentum(const ldouble momentum[4]){
    for(int j=0;j<4;j++) p[j]=momentum[j];
    update();
  }
  //Get momentum and mass
  const ldouble particle::mass() const{
    return m;
  }
  const ldouble particle::get_mass() const {
    return m;
  }
  const ldouble particle::get_momentum(const int& mu) const {
    return p[mu];
  }

  //Calc particle properties
  void particle::update(){
    constexpr ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000000;
    pxymag = std::sqrt(p[1]*p[1]+p[2]*p[2]);
    pmag = std::sqrt(pxymag*pxymag+p[3]*p[3]);
    if(std::abs(p[0]*p[0]-pmag*pmag-m*m) > epsilon){
      throw std::runtime_error("Momentum and mass don't match:  p^2 = " + std::to_string(p[0]*p[0]-pmag*pmag) + " != " + std::to_string(m*m));
    }
    theta = atan2(pxymag,p[3]);
    phi = atan2(p[2],p[1]);
    //if(phi<0) phi+=2.*PI;
    c = cdouble(std::cos(theta/2.),0.);
    s = std::polar(std::sin(theta/2.),phi);
    sc = std::conj(s);
    sqrtEpP = cdouble(std::sqrt(p[0]+pmag),0);
    sqrtEmP = cdouble(std::sqrt(p[0]-pmag),0);
    
    //Reset all the calcs.
    //The matrices and spinors
    upMat2dimCalculated = false, loMat2dimCalculated = false;
    //Helicity Spinors
    m0rangle2dimCalculated = false, m0langle2dimCalculated = false;
    m0rsquare2dimCalculated = false, m0lsquare2dimCalculated = false;
    //Spin Spinors
    //rangle
    rangleUpperP12dimCalculated = false, rangleUpperM12dimCalculated = false;
    rangleLowerP12dimCalculated = false, rangleLowerM12dimCalculated = false;
    //langle
    langleUpperP12dimCalculated = false, langleUpperM12dimCalculated = false;
    langleLowerP12dimCalculated = false, langleLowerM12dimCalculated = false;
    //lsquare
    lsquareUpperP12dimCalculated = false, lsquareUpperM12dimCalculated = false;
    lsquareLowerP12dimCalculated = false, lsquareLowerM12dimCalculated = false;
    //rsquare
    rsquareUpperP12dimCalculated = false, rsquareUpperM12dimCalculated = false;
    rsquareLowerP12dimCalculated = false, rsquareLowerM12dimCalculated = false;
  }

  //Test angles
  bool particle::test_angles() const {
    constexpr ldouble epsilon = std::numeric_limits<ldouble>::epsilon() * 1000000;
    ldouble p0 = std::real(sqrtEpP*std::conj(sqrtEpP)+sqrtEmP*std::conj(sqrtEmP))/2;
    ldouble p1 = pmag*std::sin(theta)*std::cos(phi);
    ldouble p2 = pmag*std::sin(theta)*std::sin(phi);
    ldouble p3 = pmag*std::cos(theta);
    if(std::abs(p[0]-p0) > epsilon ||
       std::abs(p[1]-p1) > epsilon ||
       std::abs(p[2]-p2) > epsilon ||
       std::abs(p[3]-p3) > epsilon
       ){
      std::cout<<"{ "<<p[0]<<" , "<<p[1]<<" , "<<p[2]<<" , "<<p[3]<<" } != ";
      std::cout<<"{ "<<p0<<" , "<<p1<<" , "<<p2<<" , "<<p3<<" }"<<std::endl;
      return false;
    }
    return true;
  }


  //dot: p1.p2
  ldouble particle::dot(const particle& p2) const{
    return p[0]*p2.p[0]-p[1]*p2.p[1]-p[2]*p2.p[2]-p[3]*p2.p[3];
  }

  


  //Matrices
  //Upper indices
  cmatrix particle::umat(const int& dim) {
    if(dim==2){
      if(!upMat2dimCalculated) {
        upMat2dim = cmatrix(cdouble(p[0]-p[3],0),-std::polar(pxymag,-phi),
		        -std::polar(pxymag,+phi),cdouble(p[0]+p[3],0));
        upMat2dimCalculated = true;
      }
    return upMat2dim;
    }
    return upMat2dim;
  }
  //Lower indices
  cmatrix particle::lmat(const int& dim) {
    if(dim==2){
      if(!loMat2dimCalculated) {
        loMat2dim = cmatrix(cdouble(p[0]+p[3],0),std::polar(pxymag,-phi),
		        std::polar(pxymag,phi),cdouble(p[0]-p[3],0));
        loMat2dimCalculated = true;
      }
      return loMat2dim;
    }
    return loMat2dim;
  }

  //Helicity Spinors
  cvector particle::rangle(const int& dim) {
    if(!m0rangle2dimCalculated){
      if(m!=0) {
	      usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m!=0");
      }
      m0rangle2dim = cvector(sqrtEpP*c,sqrtEpP*s);
      m0rangle2dimCalculated = true;
    }
    return m0rangle2dim;
  }
  cvector particle::lsquare(const int& dim) {
    if(!m0lsquare2dimCalculated){
      if(m!=0){
      	usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m!=0");
      }
      m0lsquare2dim = cvector(sqrtEpP*c,sqrtEpP*sc);
      m0lsquare2dimCalculated = true;
    }
    return m0lsquare2dim;
  }
  cvector particle::langle(const int& dim) {
    if(!m0langle2dimCalculated){
      if(m!=0){
	      usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m!=0");
      }
      m0langle2dim = cvector(sqrtEpP*s,-sqrtEpP*c);
      m0langle2dimCalculated = true;
    }
    return m0langle2dim;
  }
  cvector particle::rsquare(const int& dim) {
    if(!m0rsquare2dimCalculated){
      if(m!=0){
	      usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m!=0");
      }
      m0rsquare2dim = cvector(sqrtEpP*sc,-sqrtEpP*c);
      m0rsquare2dimCalculated = true;
    }
    return m0rsquare2dim;
  }


  //Spin Spinors
  //rangle
  cvector particle::rangle(const int& spin2, const int& dim) {
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!rangleUpperM12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	      }
	      rangleUpperM12dim = cvector(sqrtEpP*c,sqrtEpP*s);
	      rangleUpperM12dimCalculated = true;
      }
      return rangleUpperM12dim;
    }
    if(!rangleUpperP12dimCalculated){
	    if(m==0){
	      usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m==0");
	    }
	    rangleUpperP12dim = cvector(-sqrtEmP*sc,sqrtEmP*c);
	    rangleUpperP12dimCalculated = true;
    }
    return rangleUpperP12dim;
  }
  cvector particle::rangle(const int& spin2, const bool& upper, const int& dim) {
    if(upper) return rangle(spin2, dim);
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!rangleLowerM12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	      }
	      rangleLowerM12dim = cvector(-sqrtEmP*sc,sqrtEmP*c);
	      rangleLowerM12dimCalculated = true;
      }
      return rangleLowerM12dim;
    }
    if(!rangleLowerP12dimCalculated){
	    if(m==0){
	      usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m==0");
	    }
	    rangleLowerP12dim = cvector(-sqrtEpP*c,-sqrtEpP*s);
	    rangleLowerP12dimCalculated = true;
    }
    return rangleLowerP12dim;
  }
  //langle
  cvector particle::langle(const int& spin2, const int& dim) {
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!langleUpperM12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	      }
	      langleUpperM12dim = cvector(sqrtEpP*s,-sqrtEpP*c);
	      langleUpperM12dimCalculated = true;
      }
      return langleUpperM12dim;
    }
    if(!langleUpperP12dimCalculated){
	    if(m==0){
	      usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m==0");
	    }
	    langleUpperP12dim = cvector(sqrtEmP*c,sqrtEmP*sc);
	    langleUpperP12dimCalculated = true;
    }
    return langleUpperP12dim;
  }
  cvector particle::langle(const int& spin2, const bool& upper, const int& dim) {
    if(upper) return langle(spin2, dim);
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!langleLowerM12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	      }
	      langleLowerM12dim = cvector(sqrtEmP*c,sqrtEmP*sc);
	      langleLowerM12dimCalculated = true;
      }
      return langleLowerM12dim;
    }
    if(!langleLowerP12dimCalculated){
	    if(m==0){
	      usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m==0");
	    }
	    langleLowerP12dim = cvector(-sqrtEpP*s,sqrtEpP*c);
	    langleLowerP12dimCalculated = true;
    }
    return langleLowerP12dim;
  }
  //lsquare
  cvector particle::lsquare(const int& spin2, const int& dim) {
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!lsquareUpperM12dimCalculated){
	  if(m==0){
	    usage("Incorrect usage:");
	    throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	  }
	  lsquareUpperM12dim = cvector(sqrtEmP*s,-sqrtEmP*c);
	  lsquareUpperM12dimCalculated = true;
      }
      return lsquareUpperM12dim;
    }
    if(!lsquareUpperP12dimCalculated){
	    if(m==0){
	      usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m==0");
	    }
	    lsquareUpperP12dim = cvector(sqrtEpP*c,sqrtEpP*sc);
	    lsquareUpperP12dimCalculated = true;
    }
    return lsquareUpperP12dim;
  }
  cvector particle::lsquare(const int& spin2, const bool& upper, const int& dim) {
    if(upper) return lsquare(spin2,dim);
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!lsquareLowerM12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	      }
	      lsquareLowerM12dim = cvector(sqrtEpP*c,sqrtEpP*sc);
	      lsquareLowerM12dimCalculated = true;
      }
      return lsquareLowerM12dim;
    }
    if(!lsquareLowerP12dimCalculated){
	  if(m==0){
	    usage("Incorrect usage:");
	    throw std::runtime_error("Incorrect usage: m==0");
	  }
	  lsquareLowerP12dim = cvector(-sqrtEmP*s,sqrtEmP*c);
	  lsquareLowerP12dimCalculated = true;
    }
    return lsquareLowerP12dim;
  }
  //rsquare
  cvector particle::rsquare(const int& spin2, const int& dim) {
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!rsquareUpperM12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	      }
	      rsquareUpperM12dim = cvector(-sqrtEmP*c,-sqrtEmP*s);
	      rsquareUpperM12dimCalculated = true;
      }
      return rsquareUpperM12dim;
    }
    if(!rsquareUpperP12dimCalculated){
	    if(m==0){
	      usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m==0");
	    }
	    rsquareUpperP12dim = cvector(sqrtEpP*sc,-sqrtEpP*c);
	    rsquareUpperP12dimCalculated = true;
    }
    return rsquareUpperP12dim;
  }
  cvector particle::rsquare(const int& spin2, const bool& upper, const int& dim) {
    if(upper) return rsquare(spin2, dim);
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!rsquareLowerM12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	      }
	      rsquareLowerM12dim = cvector(sqrtEpP*sc,-sqrtEpP*c);
	      rsquareLowerM12dimCalculated = true;
      }
      return rsquareLowerM12dim;
    }
    if(!rsquareLowerP12dimCalculated){
	    if(m==0){
	      usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m==0");
	    }
	    rsquareLowerP12dim = cvector(sqrtEmP*c,sqrtEmP*s);
	    rsquareLowerP12dimCalculated = true;
    }
  return rsquareLowerP12dim;
  }

  

  //Massive cmatrix form with both spins
  cmatrix particle::rangle_matrix(const int& dim){
    return rangle_matrix(UPPER, dim);
  }
  cmatrix particle::rangle_matrix(const bool& upper, const int& dim){
    cvector ram=rangle(-1,upper,2), rap=rangle(+1,upper,2);
    return cmatrix(ram.get(0), rap.get(0),
		   ram.get(1), rap.get(1));
  }
  cmatrix particle::langle_matrix(const int& dim){
    return langle_matrix(UPPER,dim);
  }
  cmatrix particle::langle_matrix(const bool& upper, const int& dim){
    cvector lam=langle(-1,upper,2), lap=langle(+1,upper,2);
    return cmatrix(lam.get(0), lap.get(0),
		   lam.get(1), lap.get(1)
		   );
  }
  
  cmatrix particle::lsquare_matrix(const int& dim){
    return lsquare_matrix(UPPER, dim);
  }
  cmatrix particle::lsquare_matrix(const bool& upper, const int& dim){
    cvector lsm=lsquare(-1,upper,2), lsp=lsquare(+1,upper,2);
    return cmatrix(lsm.get(0), lsp.get(0),
		   lsm.get(1), lsp.get(1)
		   );
  }
  cmatrix particle::rsquare_matrix(const int& dim){
    return rsquare_matrix(UPPER, dim);
  }
  cmatrix particle::rsquare_matrix(const bool& upper, const int& dim){
    cvector rsm=rsquare(-1,upper,2), rsp=rsquare(+1,upper,2);
    return cmatrix(rsm.get(0), rsp.get(0),
		   rsm.get(1), rsp.get(1)
		   );
  }



  //Usage Message

  void particle::usage(const char* msg) const {
    std::cout<<msg<<std::endl;
    std::cout<<"* umat() and lmat() return a momentum matrix with upper and lower indices, respectively."<<std::endl;
    std::cout<<"* langle(int dim) and rangle(int dim) return a spinor <i| and |i> of dimension dim, respectively, for a massless particle."<<std::endl;
    std::cout<<"* lsquare(int dim) and rsquare(int dim) return a spinor [i| and |i] of dimension dim, respectively, for a massless particle."<<std::endl;
    std::cout<<"* langle(int I, int dim) and rangle(int I, int dim) return a spinor <i|^I and |i>^I with upper indices of dimension dim, respectively, for a massive particle.  I is twice the spin and can only take the values 1 and -1."<<std::endl;
    std::cout<<"* lsquare(int I, int dim) and rsquare(int I, int dim) return a spinor [i|^I and |i]^I with upper indices of dimension dim, respectively, for a massive particle.  I is twice the spin and can only take the values 1 and -1."<<std::endl;
    std::cout<<"* langle(int I, bool upper, int dim) and rangle(int I, bool upper, int dim) return a spinor <i|^I and |i>^I with upper indices with dimension dim, respectively, for a massive particle if upper is true and <i|_I and |i>_I with lower indices, respectively, if false.  I is twice the spin and can only take the values 1 and -1."<<std::endl;
    std::cout<<"* lsquare(int I, bool upper, int dim) and rsquare(int I, bool upper, int dim) return a spinor [i|^I and |i]^I with upper indices of dimension dim, respectively, for a massive particle if upper is true and [i|_I and |i]_I with lower indices, respectively, if false.  I is twice the spin and can only take the values 1 and -1."<<std::endl;
  }



  

  //Generators of Lorentz Rotations with no dot (acting on left chiral angle brackets)
  //LU acts on |j> (from the left)
  cmatrix particle::lorentz_j3_lu(const int& dim) const{
    constexpr ldouble half=0.5;
    constexpr cdouble one=cdouble(1,0);
    
    return -half*cmatrix(std::cos(theta)*one,std::polar(std::sin(theta),-phi),
		      std::polar(std::sin(theta),phi),-std::cos(theta)*one);
  }
  cmatrix particle::lorentz_jp_lu(const int& dim) const{
    constexpr ldouble half=0.5, one=1;
    ldouble pre=0;
    if(m!=0) pre = half*std::sqrt((p[0]-pmag)/(p[0]+pmag));

    return pre*cmatrix(
		       std::sin(theta), -std::polar(-one+std::cos(theta),-phi) ,
		       -std::polar(one+std::cos(theta),phi), -std::sin(theta));
  }
  cmatrix particle::lorentz_jm_lu(const int& dim) const{
    constexpr ldouble half=0.5, one=1;
    ldouble pre=1;
    if(m!=0) pre=half*std::sqrt((p[0]+pmag)/(p[0]-pmag));
    return pre*cmatrix(
		       std::sin(theta), -std::polar(one+std::cos(theta),-phi) ,
		       -std::polar(-one+std::cos(theta),phi), -std::sin(theta));
  }
  //UL acts on <j| (from the left)
  cmatrix particle::lorentz_j3_ul(const int& dim) const{
    constexpr ldouble half=0.5;
    constexpr cdouble one=cdouble(1,0);
    
    return half*cmatrix(std::cos(theta)*one,std::polar(std::sin(theta),phi),
			std::polar(std::sin(theta),-phi),-std::cos(theta)*one);
  }
  cmatrix particle::lorentz_jp_ul(const int& dim) const{
    constexpr ldouble half=0.5, one=1;
    ldouble pre=0;
    if(m!=0) pre = half*std::sqrt((p[0]-pmag)/(p[0]+pmag));

    return pre*cmatrix(
		       -std::sin(theta), std::polar(one+std::cos(theta),phi) ,
		       std::polar(-one+std::cos(theta),-phi), std::sin(theta));
  }
  cmatrix particle::lorentz_jm_ul(const int& dim) const{
    constexpr ldouble half=0.5, one=1;
    ldouble pre=1;
    if(m!=0) pre=half*std::sqrt((p[0]+pmag)/(p[0]-pmag));
    return pre*cmatrix(
		       -std::sin(theta), std::polar(-one+std::cos(theta),phi) ,
		       std::polar(one+std::cos(theta),-phi), std::sin(theta));
  }

  
  //Generators of Lorentz Rotations with a dot (acting on right chiral square brackets)
  //LU acts on [j| (from the left)
  cmatrix particle::lorentz_j3_lu_dot(const int& dim) const{
    constexpr ldouble half=0.5;
    constexpr cdouble one=cdouble(1,0);
    
    return half*cmatrix(std::cos(theta)*one,std::polar(std::sin(theta),phi),
			std::polar(std::sin(theta),-phi),-std::cos(theta)*one);
  }
  cmatrix particle::lorentz_jp_lu_dot(const int& dim) const{
    constexpr ldouble half=0.5, one=1;
    ldouble pre=1;
    if(m!=0) pre = half*std::sqrt((p[0]+pmag)/(p[0]-pmag));

    return pre*cmatrix(
		       -std::sin(theta), std::polar(one+std::cos(theta),phi) ,
		       std::polar(-one+std::cos(theta),-phi), std::sin(theta));
  }
  cmatrix particle::lorentz_jm_lu_dot(const int& dim) const{
    constexpr ldouble half=0.5, one=1;
    ldouble pre=0;
    if(m!=0) pre=half*std::sqrt((p[0]-pmag)/(p[0]+pmag));
    return pre*cmatrix(
		       -std::sin(theta), std::polar(-one+std::cos(theta),phi) ,
		       std::polar(one+std::cos(theta),-phi), std::sin(theta));
  }
  //UL acts on |j] (from the right)
  cmatrix particle::lorentz_j3_ul_dot(const int& dim) const{
    constexpr ldouble half=0.5;
    constexpr cdouble one=cdouble(1,0);
    
    return half*cmatrix(-std::cos(theta)*one,-std::polar(std::sin(theta),-phi),
			-std::polar(std::sin(theta),phi),std::cos(theta)*one);
  }
  cmatrix particle::lorentz_jp_ul_dot(const int& dim) const{
    constexpr ldouble half=0.5, one=1;
    ldouble pre=1;
    if(m!=0) pre=half*std::sqrt((p[0]+pmag)/(p[0]-pmag));

    return pre*cmatrix(
		       std::sin(theta), -std::polar(-one+std::cos(theta),-phi) ,
		       -std::polar(one+std::cos(theta),phi), -std::sin(theta));
  }
  cmatrix particle::lorentz_jm_ul_dot(const int& dim) const{
    constexpr ldouble half=0.5, one=1;
    ldouble pre=0;
    if(m!=0) pre = half*std::sqrt((p[0]-pmag)/(p[0]+pmag));
    
    return pre*cmatrix(
		       std::sin(theta), -std::polar(one+std::cos(theta),-phi) ,
		       -std::polar(-one+std::cos(theta),phi), -std::sin(theta));
  }
  
  
  //Generators of spin
  //Acts on spinors with upper indices (the default) from the right
  cmatrix particle::spin_j3_lu(const int& dim) const{
    constexpr ldouble half=0.5;
    constexpr cdouble one=cdouble(1,0), zero=cdouble(0,0);
    return half*cmatrix(-one,zero,
		       zero,one);
  }
  cmatrix particle::spin_jp_lu(const int& dim) const{
    constexpr cdouble zero=cdouble(0,0);
    constexpr ldouble one=1;
    return cmatrix(zero,zero,
		   -std::polar(one,phi),zero);
  }
  cmatrix particle::spin_jm_lu(const int& dim) const{
    constexpr cdouble zero=cdouble(0,0);
    constexpr ldouble one=1;
    return cmatrix(zero,-std::polar(one,-phi),
		   zero,zero);
  }

  
  //Acts on spinors with lower indices (not the default) from the right
  cmatrix particle::spin_j3_ul(const int& dim) const{
    constexpr ldouble half=0.5;
    constexpr cdouble one=cdouble(1,0), zero=cdouble(0,0);
    return half*cmatrix(one,zero,
		       zero,-one);
  }
  cmatrix particle::spin_jp_ul(const int& dim) const{
    constexpr cdouble zero=cdouble(0,0);
    constexpr ldouble one=1;
    return cmatrix(zero,std::polar(one,phi),
		   zero,zero);
  }
  cmatrix particle::spin_jm_ul(const int& dim) const{
    constexpr cdouble zero=cdouble(0,0);
    constexpr ldouble one=1;
    return cmatrix(zero,zero,
		   std::polar(one,-phi),zero);
  }

}
