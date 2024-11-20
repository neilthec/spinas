
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

//File:  SPINAS/source/utilities.cpp

#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <random>

#include "types.h"
#include "aliases.h"
#include "utilities.h"


namespace spinas {


  // double to Mathematica String
  std::string double_to_Mathematica_string(const ldouble& x) {
    std::stringstream tmpStr;
    std::string num;
    size_t epos;
    
    std::stringstream xStr;
    xStr << std::setprecision(std::numeric_limits<ldouble>::digits10 + 1) << x;
    num=xStr.str();
    epos=num.find("e");
    if(epos != std::string::npos)
      num.replace(epos,1," 10^");
    tmpStr<<num;
    
    return tmpStr.str();
  }
  
  

  //Rotate a Momentum by theta around axis given by normalized u[3]={ux,uy,uz}.
  void rotate_momentum(ldouble mom[4], const ldouble u[3], const ldouble& theta){
    ldouble cost = std::cos(theta), ux=u[0], uy=u[1], uz=u[2];
    ldouble sint = std::sqrt(1.0-cost*cost);
    ldouble R[4][4] = {{1,0,0,0},
		       {0, cost + ux*ux*(1 - cost), ux*uy*(1 - cost) - uz*sint, ux*uz*(1 - cost) + uy*sint},
		       {0, ux*uy*(1 - cost) + uz*sint, cost + uy*uy*(1 - cost), uy*uz*(1 - cost) - ux*sint},
		       {0, uz*ux*(1 - cost) - uy*sint, uz*uy*(1 - cost) + ux*sint,    cost + uz*uz*(1 - cost)}};
    ldouble nMom[4] = {0,0,0,0};
    for(int i=0;i<4;i++)
      for(int j=0;j<4;j++)
	nMom[i] += R[i][j]*mom[j];

    for(int i=0;i<4;i++) mom[i]=nMom[i];
  }
  //Boost a momentum by v[3]={vx,vy,vz};
  void boost_momentum(ldouble mom[4], const ldouble v[3]){
    ldouble vx=v[0], vy=v[1], vz=v[2];
    ldouble v2=vx*vx+vy*vy+vz*vz;
    ldouble gamma = 1./std::sqrt(1-v2);
    ldouble B[4][4] = {{gamma, -gamma*vx, -gamma*vy, -gamma*vz},
		       {-gamma*vx, 1 + (gamma - 1)*vx*vx/v2, (gamma - 1)*vx*vy/v2, (gamma - 1)*vx*vz/v2},
		       {-gamma*vy, (gamma - 1)*vx*vy/v2, 1 + (gamma - 1)*vy*vy/v2, (gamma - 1)*vy*vz/v2},
		       {-gamma*vz, (gamma - 1)*vx*vz/v2, (gamma - 1)*vz*vy/v2, 1 + (gamma - 1)*vz*vz/v2}};
    ldouble nMom[4] = {0,0,0,0};
    for(int i=0;i<4;i++)
      for(int j=0;j<4;j++)
	nMom[i] += B[i][j]*mom[j];

    for(int i=0;i<4;i++) mom[i]=nMom[i];
  }


  
  //Random Generation Used in Tests
  std::random_device rd;
  std::mt19937 gen(rd());  
  
  int choose_random_int(int begin, int end){
    std::uniform_int_distribution<> disInt(begin, end);
    return disInt(gen);
  }


  ldouble choose_random_ldouble(ldouble begin, ldouble end){
    std::uniform_real_distribution<ldouble> dis(begin, end);
    ldouble ld = dis(gen);
    return ld;
  }
  
  ldouble choose_random_momentum(ldouble p[4], ldouble begin, ldouble end){
    ldouble p2 = 0;
    for(int j=1;j<4;j++){
      p[j] = choose_random_ldouble(begin,end);
      p2 += p[j]*p[j];
    }
    // Ensure the range for p[0] is positive
    ldouble mass = choose_random_ldouble(0,end);
    p[0] = std::sqrt(mass*mass + p2);
    return mass;
  }
  
  void choose_random_massless_momentum(ldouble p[4], ldouble begin, ldouble end){
    ldouble p2 = 0;
    for(int j=1;j<4;j++){
      p[j] = choose_random_ldouble(begin,end);
      p2 += p[j]*p[j];
    }
    p[0] = std::sqrt(p2);
  }
  
  cdouble choose_random_cdouble(ldouble begin, ldouble end){
    ldouble ld1 = choose_random_ldouble(begin,end);
    ldouble ld2 = choose_random_ldouble(begin,end);
    return cdouble(ld1, ld2);
  }
  
}
