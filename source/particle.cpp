
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
  particle::particle():
  sqrt2(std::sqrt(2))
  {}

  particle::particle(const ldouble& mass):
  sqrt2(std::sqrt(2)),
  m(mass), p{m,0,0,0}{}

  particle::particle(const ldouble momentum[4], const ldouble& mass):
  sqrt2(std::sqrt(2)),
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
    epP = p[0]+pmag;
    emP = p[0]-pmag;
    sqrtEpP = cdouble(std::sqrt(p[0]+pmag),0);
    sqrtEmP = cdouble(std::sqrt(p[0]-pmag),0);

    eppz = cdouble(p[0]+p[3],0);
    empz = cdouble(p[0]-p[3],0);
    pxppy = cdouble(p[1],p[2]);
    pxmpy = cdouble(p[1],-p[2]);
    
    //Reset all the calcs.
    //2-dim
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

    //3-dim
    //The matrices and spinors
    upMat3dimCalculated = false, loMat3dimCalculated = false;
    //Helicity Spinors
    m0rangle3dimCalculated = false, m0langle3dimCalculated = false;
    m0rsquare3dimCalculated = false, m0lsquare3dimCalculated = false;
    //Spin Spinors
    //rangle
    rangleUpperP23dimCalculated = false, rangleUpper03dimCalculated = false, rangleUpperM23dimCalculated = false;
    rangleLowerP23dimCalculated = false, rangleLower03dimCalculated = false, rangleLowerM23dimCalculated = false;
    //langle
    langleUpperP23dimCalculated = false, langleUpper03dimCalculated = false, langleUpperM23dimCalculated = false;
    langleLowerP23dimCalculated = false, langleLower03dimCalculated = false, langleLowerM23dimCalculated = false;
    //lsquare
    lsquareUpperP23dimCalculated = false, lsquareUpper03dimCalculated = false, lsquareUpperM23dimCalculated = false;
    lsquareLowerP23dimCalculated = false, lsquareLower03dimCalculated = false, lsquareLowerM23dimCalculated = false;
    //rsquare
    rsquareUpperP23dimCalculated = false, rsquareUpper03dimCalculated = false, rsquareUpperM23dimCalculated = false;
    rsquareLowerP23dimCalculated = false, rsquareLower03dimCalculated = false, rsquareLowerM23dimCalculated = false;
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
    constexpr ldouble two=2;
    if(dim==2){
      if(!upMat2dimCalculated) {
        upMat2dim = cmatrix(cdouble(p[0]-p[3],0),-std::polar(pxymag,-phi),
		        -std::polar(pxymag,+phi),cdouble(p[0]+p[3],0));
        upMat2dimCalculated = true;
      }
    return upMat2dim;
    }
    else if (dim==3){
      if(!upMat3dimCalculated) {
        upMat3dim = cmatrix(
          empz*empz, -sqrt2*empz*pxmpy, pxmpy*pxmpy,
          -sqrt2*empz*pxppy, m*m+two*pxppy*pxmpy, -sqrt2*eppz*pxmpy,
          pxppy*pxppy, -sqrt2*eppz*pxppy, eppz*eppz);
        upMat3dimCalculated = true;
      }
      return upMat3dim;
    }
    
    return upMat2dim;
  }
  //Lower indices
  cmatrix particle::lmat(const int& dim) {
    constexpr ldouble two=2;
    if(dim==2){
      if(!loMat2dimCalculated) {
        loMat2dim = cmatrix(cdouble(p[0]+p[3],0),std::polar(pxymag,-phi),
		        std::polar(pxymag,phi),cdouble(p[0]-p[3],0));
        loMat2dimCalculated = true;
      }
      return loMat2dim;
    }
    if(dim==3){
      if(!loMat3dimCalculated) {
        loMat3dim = cmatrix(
          eppz*eppz, sqrt2*eppz*pxmpy, pxmpy*pxmpy,
          sqrt2*eppz*pxppy, m*m+two*pxppy*pxmpy, sqrt2*empz*pxmpy,
          pxppy*pxppy, sqrt2*empz*pxppy, empz*empz);
        loMat3dimCalculated = true;
      }
      return loMat3dim;      
    }
    return loMat2dim;
  }

  //Helicity Spinors
  cvector particle::rangle(const int& dim) {
    constexpr ldouble two=2;
    if(dim==2){
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
    else if (dim==3){
      if(!m0rangle3dimCalculated){
        if(m!=0) {
          usage("Incorrect usage:");
          throw std::runtime_error("Incorrect usage: m!=0");
        }
        m0rangle3dim = two*p[0]*cvector(c*c,sqrt2*c*s,s*s);
        m0rangle3dimCalculated = true;
      }
      return m0rangle3dim;
    }
    
    return m0rangle2dim;
  }
  cvector particle::lsquare(const int& dim) {
    constexpr ldouble two=2;
    if(dim==2){
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
    else if (dim==3){
      if(!m0lsquare3dimCalculated){
        if(m!=0){
          usage("Incorrect usage:");
          throw std::runtime_error("Incorrect usage: m!=0");
        }
        m0lsquare3dim = two*p[0]*cvector(c*c,sqrt2*c*sc,sc*sc);
        m0lsquare3dimCalculated = true;
      }
      return m0lsquare3dim;
    }
    
    return m0lsquare2dim;
  }
  cvector particle::langle(const int& dim) {
    constexpr ldouble two=2;
    if(dim==2){
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
    else if(dim==3){
      if(!m0langle3dimCalculated){
        if(m!=0){
          usage("Incorrect usage:");
          throw std::runtime_error("Incorrect usage: m!=0");
        }
        m0langle3dim = two*p[0]*cvector(s*s,-sqrt2*c*s,c*c);
        m0langle3dimCalculated = true;
      }
      return m0langle3dim;
    }
    return m0langle2dim;
  }
  cvector particle::rsquare(const int& dim) {
    constexpr ldouble two=2;
    if(dim==2){
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
    else if(dim==3){
      if(!m0rsquare3dimCalculated){
        if(m!=0){
          usage("Incorrect usage:");
          throw std::runtime_error("Incorrect usage: m!=0");
        }
        m0rsquare3dim = two*p[0]*cvector(sc*sc,-sqrt2*c*sc,c*c);
        m0rsquare3dimCalculated = true;
      }
      return m0rsquare3dim;
    }
    return m0rsquare2dim;
  }


  //Spin Spinors
  /********************************************************************
                            |j>^I
  *********************************************************************/
  cvector particle::rangle(const int& spin2, const int& dim) {
    if(dim==2){
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
	        rangleUpperM12dim = sqrtEpP*cvector(c,s);
	        rangleUpperM12dimCalculated = true;
        }
        return rangleUpperM12dim;
      }
      if(!rangleUpperP12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: m==0");
	      }
	      rangleUpperP12dim = sqrtEmP*cvector(-sc,c);
	      rangleUpperP12dimCalculated = true;
      }
      return rangleUpperP12dim;
    }
    if(dim==3){
      if(spin2!=2&&spin2!=0&&spin2!=-2){
        usage("Incorrect usage:");
        throw std::runtime_error("Incorrect usage: spin2 != 2, 0 or -2");
      }
      if(spin2==-2){
        if(!rangleUpperM23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==-2 and m==0");
          }
          rangleUpperM23dim = epP*cvector(c*c,sqrt2*c*s,s*s);//empz*cvector(c*c,sqrt2*c*sc,sc*sc)
          rangleUpperM23dimCalculated = true;
        }
        return rangleUpperM23dim;
      }
      else if(spin2==0){
        if(!rangleUpper03dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==0 and m==0");
          }
          rangleUpper03dim = m*cvector(-sqrt2*c*sc,c*c-s*sc,sqrt2*c*s);
          rangleUpper03dimCalculated = true;
        }
        return rangleUpper03dim;
      }
      else if(spin2==2){
        if(!rangleUpperP23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==2 and m==0");
          }
          rangleUpperP23dim = emP*cvector(sc*sc,-sqrt2*c*sc,c*c);
          rangleUpperP23dimCalculated = true;
        }
        return rangleUpperP23dim;
      }
    }
    return rangleUpperP12dim;
  }
  /********************************************************************
                            |j>_I
  *********************************************************************/
  cvector particle::rangle(const int& spin2, const bool& upper, const int& dim) {
    if(upper) return rangle(spin2, dim);
    if(dim==2){
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
	        rangleLowerM12dim = -sqrtEpP*cvector(c,s);
	        rangleLowerM12dimCalculated = true;
        }
        return rangleLowerM12dim;
      }
      if(!rangleLowerP12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: m==0");
	      }
	      rangleLowerP12dim = sqrtEmP*cvector(-sc,c);
	      rangleLowerP12dimCalculated = true;
      }
      return rangleLowerP12dim;
    }
    else if(dim==3){
      if(spin2!=2&&spin2!=0&&spin2!=-2){
        usage("Incorrect usage:");
        throw std::runtime_error("Incorrect usage: spin2 != 2, 0 or -2");
      }
      if(spin2==2){
        if(!rangleLowerM23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==-2 and m==0");
          }
          rangleLowerM23dim = emP*cvector(sc*sc,-sqrt2*c*sc,c*c);
          rangleLowerM23dimCalculated = true;
        }
        return rangleLowerM23dim;
      }
      else if(spin2==0){
        if(!rangleLower03dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==0 and m==0");
          }
          rangleLower03dim = -m*cvector(-sqrt2*c*sc,c*c-s*sc,sqrt2*c*s);
          rangleLower03dimCalculated = true;
        }
        return rangleLower03dim;
      }
      else if(spin2==-2){
        if(!rangleLowerP23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==2 and m==0");
          }
          rangleLowerP23dim = epP*cvector(c*c,sqrt2*c*s,s*s);
          rangleLowerP23dimCalculated = true;
        }
        return rangleLowerP23dim;
      }
    }
    return rangleLowerP12dim;
  }
/********************************************************************
                            <j|^I
*********************************************************************/
  cvector particle::langle(const int& spin2, const int& dim) {
    if(dim==2){
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
	        langleUpperM12dim = sqrtEpP*cvector(s,-c);
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
    else if(dim==3){
      if(spin2!=2&&spin2!=0&&spin2!=-2){
        usage("Incorrect usage:");
        throw std::runtime_error("Incorrect usage: spin2 != 2, 0 or -2");
      }
      if(spin2==-2){
        if(!langleUpperM23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==-2 and m==0");
          }
          langleUpperM23dim = epP*cvector(s*s,-sqrt2*c*s,c*c);
          langleUpperM23dimCalculated = true;
        }
        return langleUpperM23dim;
      }
      else if(spin2==0){
        if(!langleUpper03dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==0 and m==0");
          }
          langleUpper03dim = m*cvector(sqrt2*c*s,-c*c+s*sc,-sqrt2*c*sc);
          langleUpper03dimCalculated = true;
        }
        return langleUpper03dim;
      }
      else if(spin2==2){
        if(!langleUpperP23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==2 and m==0");
          }
          langleUpperP23dim = emP*cvector(c*c,sqrt2*c*sc,sc*sc);
          langleUpperP23dimCalculated = true;
        }
        return langleUpperP23dim;
      }
    }
    return langleUpperM12dim;
  }
/********************************************************************
                            <j|_I
*********************************************************************/
cvector particle::langle(const int& spin2, const bool& upper, const int& dim) {
    if(upper) return langle(spin2, dim);
    if(dim==2){
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
	        langleLowerM12dim = -sqrtEpP*cvector(s,-c);
	        langleLowerM12dimCalculated = true;
        }
        return langleLowerM12dim;
      }
      if(!langleLowerP12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: m==0");
	      }
	      langleLowerP12dim = sqrtEmP*cvector(c,sc);
	      langleLowerP12dimCalculated = true;
      }
      return langleLowerP12dim;
    }
    else if(dim==3){
      if(spin2!=2&&spin2!=0&&spin2!=-2){
        usage("Incorrect usage:");
        throw std::runtime_error("Incorrect usage: spin2 != 2, 0 or -2");
      }
      if(spin2==2){
        if(!langleLowerM23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==-2 and m==0");
          }
          langleLowerM23dim = emP*cvector(c*c,sqrt2*c*sc,sc*sc);
          langleLowerM23dimCalculated = true;
        }
        return langleLowerM23dim;
      }
      else if(spin2==0){
        if(!langleLower03dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==0 and m==0");
          }
          langleLower03dim = -m*cvector(sqrt2*c*s,-c*c+s*sc,-sqrt2*c*sc);
          langleLower03dimCalculated = true;
        }
        return langleLower03dim;
      }
      else if(spin2==-2){
        if(!langleLowerP23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==2 and m==0");
          }
          langleLowerP23dim = epP*cvector(s*s,-sqrt2*c*s,c*c);
          langleLowerP23dimCalculated = true;
        }
        return langleLowerP23dim;
      }
    }
    return langleLowerP12dim;
  }
/********************************************************************
                            [j|^I
*********************************************************************/
  cvector particle::lsquare(const int& spin2, const int& dim) {
    if(dim==2){
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
	        lsquareUpperM12dim = -sqrtEmP*cvector(-s,c);
	        lsquareUpperM12dimCalculated = true;
        }
        return lsquareUpperM12dim;
      }
      if(!lsquareUpperP12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: m==0");
	      }
	      lsquareUpperP12dim = sqrtEpP*cvector(c,sc);
	      lsquareUpperP12dimCalculated = true;
      }
      return lsquareUpperP12dim;
    }
    else if(dim==3){
      if(spin2!=2&&spin2!=0&&spin2!=-2){
        usage("Incorrect usage:");
        throw std::runtime_error("Incorrect usage: spin2 != 2, 0 or -2");
      }
      if(spin2==-2){
        if(!lsquareUpperM23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==-2 and m==0");
          }
          lsquareUpperM23dim = emP*cvector(s*s,-sqrt2*c*s,c*c);
          lsquareUpperM23dimCalculated = true;
        }
        return lsquareUpperM23dim;
      }
      else if(spin2==0){
        if(!lsquareUpper03dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==0 and m==0");
          }
          lsquareUpper03dim = -m*cvector(-sqrt2*c*s,c*c-s*sc,sqrt2*c*sc);
          lsquareUpper03dimCalculated = true;
        }
        return lsquareUpper03dim;
      }
      else if(spin2==2){
        if(!lsquareUpperP23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==2 and m==0");
          }
          lsquareUpperP23dim = epP*cvector(c*c,sqrt2*c*sc,sc*sc);
          lsquareUpperP23dimCalculated = true;
        }
        return lsquareUpperP23dim;
      }
    }
    return lsquareUpperP12dim;
  }
/********************************************************************
                            [j|_I
*********************************************************************/
  cvector particle::lsquare(const int& spin2, const bool& upper, const int& dim) {
    if(upper) return lsquare(spin2,dim);
    if(dim==2){
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
	        lsquareLowerM12dim = sqrtEmP*cvector(-s,c);
	        lsquareLowerM12dimCalculated = true;
        }
        return lsquareLowerM12dim;
      }
      if(!lsquareLowerP12dimCalculated){
	    if(m==0){
	      usage("Incorrect usage:");
	      throw std::runtime_error("Incorrect usage: m==0");
	    }
	    lsquareLowerP12dim = sqrtEpP*cvector(c,sc);
	    lsquareLowerP12dimCalculated = true;
      }
      return lsquareLowerP12dim;
    }
    else if(dim==3){
      if(spin2!=2&&spin2!=0&&spin2!=-2){
        usage("Incorrect usage:");
        throw std::runtime_error("Incorrect usage: spin2 != 2, 0 or -2");
      }
      if(spin2==-2){
        if(!lsquareLowerM23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==-2 and m==0");
          }
          lsquareLowerM23dim = emP*cvector(s*s,-sqrt2*c*s,c*c);
          lsquareLowerM23dimCalculated = true;
        }
        return lsquareLowerM23dim;
      }
      else if(spin2==0){
        if(!lsquareLower03dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==0 and m==0");
          }
          lsquareLower03dim = m*cvector(-sqrt2*c*s,c*c-s*sc,sqrt2*c*sc);
          lsquareLower03dimCalculated = true;
        }
        return lsquareLower03dim;
      }
      else if(spin2==2){
        if(!lsquareLowerP23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==2 and m==0");
          }
          lsquareLowerP23dim = epP*cvector(c*c,sqrt2*c*sc,sc*sc);
          lsquareLowerP23dimCalculated = true;
        }
        return lsquareLowerP23dim;
      }
    }
    return lsquareLowerP12dim;
  }
/********************************************************************
                            |j]^I
*********************************************************************/
  cvector particle::rsquare(const int& spin2, const int& dim) {
    if(dim==2){
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
	        rsquareUpperM12dim = -sqrtEmP*cvector(c,s);
	        rsquareUpperM12dimCalculated = true;
        }
        return rsquareUpperM12dim;
      }
      if(!rsquareUpperP12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: m==0");
	      }
	      rsquareUpperP12dim = sqrtEpP*cvector(sc,-c);
	      rsquareUpperP12dimCalculated = true;
      }
      return rsquareUpperP12dim;
    }
    else if(dim==3){
      if(spin2!=2&&spin2!=0&&spin2!=-2){
        usage("Incorrect usage:");
        throw std::runtime_error("Incorrect usage: spin2 != 2, 0 or -2");
      }
      if(spin2==-2){
        if(!rsquareUpperM23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==-2 and m==0");
          }
          rsquareUpperM23dim = emP*cvector(c*c,sqrt2*c*s,s*s);
          rsquareUpperM23dimCalculated = true;
        }
        return rsquareUpperM23dim;
      }
      else if(spin2==0){
        if(!rsquareUpper03dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==0 and m==0");
          }
          rsquareUpper03dim = -m*cvector(sqrt2*c*sc,-c*c+s*sc,-sqrt2*c*s);
          rsquareUpper03dimCalculated = true;
        }
        return rsquareUpper03dim;
      }
      else if(spin2==2){
        if(!rsquareUpperP23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==2 and m==0");
          }
          rsquareUpperP23dim = epP*cvector(sc*sc,-sqrt2*c*sc,c*c);
          rsquareUpperP23dimCalculated = true;
        }
        return rsquareUpperP23dim;
      }
    }
    return rsquareUpperP12dim;
  }
/********************************************************************
                            |j]_I
*********************************************************************/
  cvector particle::rsquare(const int& spin2, const bool& upper, const int& dim) {
    if(upper) return rsquare(spin2, dim);
    if(dim==2){
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
	        rsquareLowerM12dim = sqrtEmP*cvector(c,s);
	        rsquareLowerM12dimCalculated = true;
        }
        return rsquareLowerM12dim;
      }
      if(!rsquareLowerP12dimCalculated){
	      if(m==0){
	        usage("Incorrect usage:");
	        throw std::runtime_error("Incorrect usage: m==0");
	      }
	      rsquareLowerP12dim = sqrtEpP*cvector(sc,-c);
	      rsquareLowerP12dimCalculated = true;
      }
      return rsquareLowerP12dim;
    }
    else if(dim==3){
      if(spin2!=2&&spin2!=0&&spin2!=-2){
        usage("Incorrect usage:");
        throw std::runtime_error("Incorrect usage: spin2 != 2, 0 or -2");
      }
      if(spin2==-2){
        if(!rsquareLowerM23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==-2 and m==0");
          }
          rsquareLowerM23dim = emP*cvector(c*c,sqrt2*c*s,s*s);
          rsquareLowerM23dimCalculated = true;
        }
        return rsquareLowerM23dim;
      }
      else if(spin2==0){
        if(!rsquareLower03dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==0 and m==0");
          }
          rsquareLower03dim = m*cvector(sqrt2*c*sc,-c*c+s*sc,-sqrt2*c*s);
          rsquareLower03dimCalculated = true;
        }
        return rsquareLower03dim;
      }
      else if(spin2==2){
        if(!rsquareLowerP23dimCalculated){
          if(m==0){
            usage("Incorrect usage:");
            throw std::runtime_error("Incorrect usage: spin2==2 and m==0");
          }
          rsquareLowerP23dim = epP*cvector(sc*sc,-sqrt2*c*sc,c*c);
          rsquareLowerP23dimCalculated = true;
        }
        return rsquareLowerP23dim;
      }
    }
    return rsquareLowerP12dim;
  }

  

  //Massive cmatrix form with both spins
  cmatrix particle::rangle_matrix(const int& dim){
    return rangle_matrix(UPPER, dim);
  }
  cmatrix particle::rangle_matrix(const bool& upper, const int& dim){
    int flip=1;
    if(!upper) flip=-1;
    if(dim==2){
      cvector ral=rangle(-flip,upper,2), rar=rangle(+flip,upper,2);
      return cmatrix(ral.get(0), rar.get(0),
		    ral.get(1), rar.get(1));
    }
    else if(dim==3){
      cvector ral=rangle(-2*flip,upper,3), ra0=rangle(0,upper,3), rar=rangle(+2*flip,upper,3);
      return cmatrix(ral.get(0), ra0.get(0), rar.get(0),
        ral.get(1), ra0.get(1), rar.get(1),
        ral.get(2), ra0.get(2), rar.get(2));
    }
    return cmatrix(0,0,0,0);
  }
  cmatrix particle::langle_matrix(const int& dim){
    return langle_matrix(UPPER,dim);
  }
  cmatrix particle::langle_matrix(const bool& upper, const int& dim){
    int flip=1;
    if(!upper) flip=-1;
    if(dim==2){
      cvector lal=langle(-flip,upper,2), lar=langle(+flip,upper,2);
      return cmatrix(lal.get(0), lar.get(0),
		     lal.get(1), lar.get(1)
		     );
    }
    else if(dim==3){
      cvector lal=langle(-2*flip,upper,3), la0=langle(0,upper,3), lar=langle(+2*flip,upper,3);
      return cmatrix(lal.get(0), la0.get(0), lar.get(0),
         lal.get(1), la0.get(1), lar.get(1),
         lal.get(2), la0.get(2), lar.get(2)
         );
    }
    return cmatrix(0,0,0,0);
  }
  
  cmatrix particle::lsquare_matrix(const int& dim){
    return lsquare_matrix(UPPER, dim);
  }
  cmatrix particle::lsquare_matrix(const bool& upper, const int& dim){
    int flip=1;
    if(!upper) flip=-1;
    if(dim==2){
      cvector lsl=lsquare(-flip,upper,2), lsr=lsquare(+flip,upper,2);
      return cmatrix(lsl.get(0), lsr.get(0),
		     lsl.get(1), lsr.get(1)
		     );
    }
    else if(dim==3){
      cvector lsl=lsquare(-2*flip,upper,3), ls0=lsquare(0,upper,3), lsr=lsquare(+2*flip,upper,3);
      return cmatrix(lsl.get(0), ls0.get(0), lsr.get(0),
         lsl.get(1), ls0.get(1), lsr.get(1),
         lsl.get(2), ls0.get(2), lsr.get(2)
         );
    }
    return cmatrix(0,0,0,0);
  }
  cmatrix particle::rsquare_matrix(const int& dim){
    return rsquare_matrix(UPPER, dim);
  }
  cmatrix particle::rsquare_matrix(const bool& upper, const int& dim){
    int flip=1;
    if(!upper) flip=-1;
    if(dim==2){
      cvector rsl=rsquare(-flip,upper,2), rsr=rsquare(+flip,upper,2);
      return cmatrix(rsl.get(0), rsr.get(0),
		     rsl.get(1), rsr.get(1)
		     );
    }
    else if(dim==3){
      cvector rsl=rsquare(-2*flip,upper,3), rs0=rsquare(0,upper,3), rsr=rsquare(+2*flip,upper,3);
      return cmatrix(rsl.get(0), rs0.get(0), rsr.get(0),
         rsl.get(1), rs0.get(1), rsr.get(1),
         rsl.get(2), rs0.get(2), rsr.get(2)
         );
    }
    return cmatrix(0,0,0,0);
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
  //Todo: Still need to construct these for 3 dimensions.
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
