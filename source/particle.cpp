
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
#include "aliases.h"
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
    upMatCalculated = false, loMatCalculated = false;
    //Helicity Spinors
    m0rangleCalculated = false, m0langleCalculated = false;
    m0rsquareCalculated = false, m0lsquareCalculated = false;
    //Spin Spinors
    //rangle
    rangleUpperP1Calculated = false, rangleUpperM1Calculated = false;
    rangleLowerP1Calculated = false, rangleLowerM1Calculated = false;
    //langle
    langleUpperP1Calculated = false, langleUpperM1Calculated = false;
    langleLowerP1Calculated = false, langleLowerM1Calculated = false;
    //lsquare
    lsquareUpperP1Calculated = false, lsquareUpperM1Calculated = false;
    lsquareLowerP1Calculated = false, lsquareLowerM1Calculated = false;
    //rsquare
    rsquareUpperP1Calculated = false, rsquareUpperM1Calculated = false;
    rsquareLowerP1Calculated = false, rsquareLowerM1Calculated = false;
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
  cmatrix particle::umat() {
    if(!upMatCalculated) {
      upMat = cmatrix(cdouble(p[0]-p[3],0),-std::polar(pxymag,-phi),
		      -std::polar(pxymag,+phi),cdouble(p[0]+p[3],0));
      upMatCalculated = true;
    }
    return upMat;
  }
  //Lower indices
  cmatrix particle::lmat() {
    if(!loMatCalculated) {
      loMat = cmatrix(cdouble(p[0]+p[3],0),std::polar(pxymag,-phi),
		      std::polar(pxymag,phi),cdouble(p[0]-p[3],0));
      loMatCalculated = true;
    }
    return loMat;
  }

  //Helicity Spinors
  cvector particle::rangle() {
    if(!m0rangleCalculated){
      if(m!=0) {
	usage("Incorrect usage:");
	throw std::runtime_error("Incorrect usage: m!=0");
      }
      m0rangle = cvector(sqrtEpP*c,sqrtEpP*s);
      m0rangleCalculated = true;
    }
    return m0rangle;
  }
  cvector particle::lsquare() {
    if(!m0lsquareCalculated){
      if(m!=0){
	usage("Incorrect usage:");
	throw std::runtime_error("Incorrect usage: m!=0");
      }
      m0lsquare = cvector(sqrtEpP*c,sqrtEpP*sc);
      m0lsquareCalculated = true;
    }
    return m0lsquare;
  }
  cvector particle::langle() {
    if(!m0langleCalculated){
      if(m!=0){
	usage("Incorrect usage:");
	throw std::runtime_error("Incorrect usage: m!=0");
      }
      m0langle = cvector(sqrtEpP*s,-sqrtEpP*c);
      m0langleCalculated = true;
    }
    return m0langle;
  }
  cvector particle::rsquare() {
    if(!m0rsquareCalculated){
      if(m!=0){
	usage("Incorrect usage:");
	throw std::runtime_error("Incorrect usage: m!=0");
      }
      m0rsquare = cvector(sqrtEpP*sc,-sqrtEpP*c);
      m0rsquareCalculated = true;
    }
    return m0rsquare;
  }


  //Spin Spinors
  //rangle
  cvector particle::rangle(const int& spin2) {
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!rangleUpperM1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	}
	rangleUpperM1 = cvector(sqrtEpP*c,sqrtEpP*s);
	rangleUpperM1Calculated = true;
      }
      return rangleUpperM1;
    }
    if(!rangleUpperP1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: m==0");
	}
	rangleUpperP1 = cvector(-sqrtEmP*sc,sqrtEmP*c);
	rangleUpperP1Calculated = true;
      }
      return rangleUpperP1;
  }
  cvector particle::rangle(const int& spin2, const bool& upper) {
    if(upper) return rangle(spin2);
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!rangleLowerM1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	}
	rangleLowerM1 = cvector(-sqrtEmP*sc,sqrtEmP*c);
	rangleLowerM1Calculated = true;
      }
      return rangleLowerM1;
    }
    if(!rangleLowerP1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: m==0");
	}
	rangleLowerP1 = cvector(-sqrtEpP*c,-sqrtEpP*s);
	rangleLowerP1Calculated = true;
      }
      return rangleLowerP1;
  }
  //langle
  cvector particle::langle(const int& spin2) {
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!langleUpperM1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	}
	langleUpperM1 = cvector(sqrtEpP*s,-sqrtEpP*c);
	langleUpperM1Calculated = true;
      }
      return langleUpperM1;
    }
    if(!langleUpperP1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: m==0");
	}
	langleUpperP1 = cvector(sqrtEmP*c,sqrtEmP*sc);
	langleUpperP1Calculated = true;
      }
      return langleUpperP1;
  }
  cvector particle::langle(const int& spin2, const bool& upper) {
    if(upper) return langle(spin2);
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!langleLowerM1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	}
	langleLowerM1 = cvector(sqrtEmP*c,sqrtEmP*sc);
	langleLowerM1Calculated = true;
      }
      return langleLowerM1;
    }
    if(!langleLowerP1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: m==0");
	}
	langleLowerP1 = cvector(-sqrtEpP*s,sqrtEpP*c);
	langleLowerP1Calculated = true;
      }
      return langleLowerP1;
  }
  //lsquare
  cvector particle::lsquare(const int& spin2) {
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!lsquareUpperM1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	}
	lsquareUpperM1 = cvector(sqrtEmP*s,-sqrtEmP*c);
	lsquareUpperM1Calculated = true;
      }
      return lsquareUpperM1;
    }
    if(!lsquareUpperP1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: m==0");
	}
	lsquareUpperP1 = cvector(sqrtEpP*c,sqrtEpP*sc);
	lsquareUpperP1Calculated = true;
      }
      return lsquareUpperP1;
  }
  cvector particle::lsquare(const int& spin2, const bool& upper) {
    if(upper) return lsquare(spin2);
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!lsquareLowerM1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	}
	lsquareLowerM1 = cvector(sqrtEpP*c,sqrtEpP*sc);
	lsquareLowerM1Calculated = true;
      }
      return lsquareLowerM1;
    }
    if(!lsquareLowerP1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: m==0");
	}
	lsquareLowerP1 = cvector(-sqrtEmP*s,sqrtEmP*c);
	lsquareLowerP1Calculated = true;
      }
    return lsquareLowerP1;
  }
  //rsquare
  cvector particle::rsquare(const int& spin2) {
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!rsquareUpperM1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	}
	rsquareUpperM1 = cvector(-sqrtEmP*c,-sqrtEmP*s);
	rsquareUpperM1Calculated = true;
      }
      return rsquareUpperM1;
    }
    if(!rsquareUpperP1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: m==0");
	}
	rsquareUpperP1 = cvector(sqrtEpP*sc,-sqrtEpP*c);
	rsquareUpperP1Calculated = true;
      }
      return rsquareUpperP1;
  }
  cvector particle::rsquare(const int& spin2, const bool& upper) {
    if(upper) return rsquare(spin2);
    if(spin2!=1&&spin2!=-1){
      usage("Incorrect usage:");
      throw std::runtime_error("Incorrect usage: spin2 != 1 or -1");
    }
    if(spin2==-1){
      if(!rsquareLowerM1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: spin2==-1 and m==0");
	}
	rsquareLowerM1 = cvector(sqrtEpP*sc,-sqrtEpP*c);
	rsquareLowerM1Calculated = true;
      }
      return rsquareLowerM1;
    }
    if(!rsquareLowerP1Calculated){
	if(m==0){
	  usage("Incorrect usage:");
	  throw std::runtime_error("Incorrect usage: m==0");
	}
	rsquareLowerP1 = cvector(sqrtEmP*c,sqrtEmP*s);
	rsquareLowerP1Calculated = true;
      }
    return rsquareLowerP1;
  }




  //Usage Message

  void particle::usage(const char* msg) const {
    std::cout<<msg<<std::endl;
    std::cout<<"* umat() and lmat() return a momentum matrix with upper and lower indices, respectively."<<std::endl;
    std::cout<<"* langle() and rangle() return a spinor <i| and |i>, respectively, for a massless particle."<<std::endl;
    std::cout<<"* lsquare() and rsquare() return a spinor [i| and |i], respectively, for a massless particle."<<std::endl;
    std::cout<<"* langle(int I) and rangle(int I) return a spinor <i|^I and |i>^I with upper indices, respectively, for a massive particle.  I is twice the spin and can only take the values 1 and -1."<<std::endl;
    std::cout<<"* lsquare(int I) and rsquare(int I) return a spinor [i|^I and |i]^I with upper indices, respectively, for a massive particle.  I is twice the spin and can only take the values 1 and -1."<<std::endl;
    std::cout<<"* langle(int I, bool upper) and rangle(int I, bool upper) return a spinor <i|^I and |i>^I with upper indices, respectively, for a massive particle if upper is true and <i|_I and |i>_I with lower indices, respectively, if false.  I is twice the spin and can only take the values 1 and -1."<<std::endl;
    std::cout<<"* lsquare(int I, bool upper) and rsquare(int I, bool upper) return a spinor [i|^I and |i]^I with upper indices, respectively, for a massive particle if upper is true and [i|_I and |i]_I with lower indices, respectively, if false.  I is twice the spin and can only take the values 1 and -1."<<std::endl;
  }




}
