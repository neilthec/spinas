
/*
SPINAS - Spinor Amplitudes
Copyright (C) 2023-2026 Neil Christensen, Nick Majestic

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

//File:  SPINAS/SM/eeAAA.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>
#include <random>

#include "spinas.h"
#include "include/eeAAA.h"

namespace spinas {

  eeAAA::eeAAA(const ldouble& echarge, const ldouble& masse):
    e(echarge), me(masse), prop(masse,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(me);
    p2=particle(me);
    p3=particle(0);
    p4=particle(0);
    p5=particle(0);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s25s = sproduct(SQUARE,&p2,&p5);
    a25a = sproduct(ANGLE,&p2,&p5);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s15s = sproduct(SQUARE,&p1,&p5);
    a15a = sproduct(ANGLE,&p1,&p5);
    s35s = sproduct(SQUARE,&p3,&p5);
    a35a = sproduct(ANGLE,&p3,&p5);
    s45s = sproduct(SQUARE,&p4,&p5);
    a45a = sproduct(ANGLE,&p4,&p5);
    s134a = sproduct(SQUARE,&p1,&p3,&p4);
    s135a = sproduct(SQUARE,&p1,&p3,&p5);
    s143a = sproduct(SQUARE,&p1,&p4,&p3);
    s145a = sproduct(SQUARE,&p1,&p4,&p5);
    s153a = sproduct(SQUARE,&p1,&p5,&p3);
    s234a = sproduct(SQUARE,&p2,&p3,&p4);
    s243a = sproduct(SQUARE,&p2,&p4,&p3);
    s245a = sproduct(SQUARE,&p2,&p4,&p5);
    s253a = sproduct(SQUARE,&p2,&p5,&p3);
    s254a = sproduct(SQUARE,&p2,&p5,&p4);
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    s315a = sproduct(SQUARE,&p3,&p1,&p5);
    s324a = sproduct(SQUARE,&p3,&p2,&p4);
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
    s414a = sproduct(SQUARE,&p4,&p1,&p4);
    s415a = sproduct(SQUARE,&p4,&p1,&p5);
    s423a = sproduct(SQUARE,&p4,&p2,&p3);
    s425a = sproduct(SQUARE,&p4,&p2,&p5);
    s453a = sproduct(SQUARE,&p4,&p5,&p3);
    s434a = sproduct(SQUARE,&p4,&p3,&p4);
    s513a = sproduct(SQUARE,&p5,&p1,&p3);
    s523a = sproduct(SQUARE,&p5,&p2,&p3);
    s524a = sproduct(SQUARE,&p5,&p2,&p4);
    s525a = sproduct(SQUARE,&p5,&p2,&p5);
    s3123s = sproduct(SQUARE,&p3,&p1,&p2,&p3);
    a3123a = sproduct(ANGLE,&p3,&p1,&p2,&p3);
    s4124s = sproduct(SQUARE,&p4,&p1,&p2,&p4);
    a4124a = sproduct(ANGLE,&p4,&p1,&p2,&p4);
    s5125s = sproduct(SQUARE,&p5,&p1,&p2,&p5);
    a5125a = sproduct(ANGLE,&p5,&p1,&p2,&p5);
    s3145s = sproduct(SQUARE,&p3,&p1,&p4,&p5);
    s3154s = sproduct(SQUARE,&p3,&p1,&p5,&p4);
    s4135s = sproduct(SQUARE,&p4,&p1,&p3,&p5);
    s4153s = sproduct(SQUARE,&p4,&p1,&p5,&p3);
    s5134s = sproduct(SQUARE,&p5,&p1,&p3,&p4);
    s5143s = sproduct(SQUARE,&p5,&p1,&p4,&p3);
  }
  void eeAAA::set_masses(const ldouble& masse){
    me=masse;
    p1.set_mass(me);
    p2.set_mass(me);
    prop.set_mass(me);
  }
  void eeAAA::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4], const ldouble mom5[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    p5.set_momentum(mom5);
    s34s.update();
    a34a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    s15s.update();
    a15a.update();
    s24s.update();
    a24a.update();
    s23s.update();
    a23a.update();
    s25s.update();
    a25a.update();
    s14s.update();
    a14a.update();
    s35s.update();
    a35a.update();
    s45s.update();
    a45a.update();
    s134a.update();
    s135a.update();
    s143a.update();
    s145a.update();
    s153a.update();
    s234a.update();
    s243a.update();
    s245a.update();
    s253a.update();
    s254a.update();
    s314a.update();
    s315a.update();
    s324a.update();
    s413a.update();
    s414a.update();
    s415a.update();
    s423a.update();
    s425a.update();
    s453a.update();
    s434a.update();
    s513a.update();
    s523a.update();
    s524a.update();
    s525a.update();
    s3123s.update();
    a3123a.update();
    s4124s.update();
    a4124a.update();
    s5125s.update();
    a5125a.update();
    s3145s.update();
    s3154s.update();
    s4135s.update();
    s4153s.update();
    s5134s.update();
    s5143s.update();
    //Propagator Momentum
    ldouble propS13P[4], propS14P[4], propS15P[4], propS23P[4], propS24P[4], propS25P[4];
    for(int j=0;j<4;j++){
      propS13P[j] = mom1[j]-mom3[j];
      propS14P[j] = mom1[j]-mom4[j];
      propS15P[j] = mom1[j]-mom5[j];
      propS23P[j] = mom2[j]-mom3[j];
      propS24P[j] = mom2[j]-mom4[j];
      propS25P[j] = mom2[j]-mom5[j];
    }
    pDenS13 = prop.denominator(propS13P);
    pDenS14 = prop.denominator(propS14P);
    pDenS15 = prop.denominator(propS15P);
    pDenS23 = prop.denominator(propS23P);
    pDenS24 = prop.denominator(propS24P);
    pDenS25 = prop.denominator(propS25P);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eeAAA::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5){
    
    cdouble one(1,0);
    cdouble two(2,0);

    if(ds3>0&&ds4>0&&ds5>0){

      //<12>([45]^2[3|p_1p_2|3](*2 denominators) + [35]^2)
      return sqrt(2)*sqrt(2)*sqrt(2)*e*e*e*me*a12a.v(ds1,ds2)*(
        s45s.v()*s45s.v()*s3123s.v()/pDenS13/pDenS23*((one/pDenS24/pDenS25)+(one/pDenS14/pDenS15)) +
        s35s.v()*s35s.v()*s4124s.v()/pDenS14/pDenS24*((one/pDenS23/pDenS25)+(one/pDenS13/pDenS15)) +
        s34s.v()*s34s.v()*s5125s.v()/pDenS15/pDenS25*((one/pDenS23/pDenS24)+(one/pDenS13/pDenS14))
      )/two;
    }

    return cdouble(0,0);    
  }

  cdouble eeAAA::amp_permutation(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5) {
    
    cdouble one(1,0);
    
    if(ds3>0&&ds4>0&&ds5>0){

      return -sqrt(2)*sqrt(2)*sqrt(2)*e*e*e*me*a12a.v(ds1,ds2)*(
        s3145s.v() / (pDenS13 * pDenS25 * a34a.v() * a45a.v())
        - s3154s.v() / (pDenS13 * pDenS24 * a35a.v() * a45a.v())
        - s4135s.v() / (pDenS14 * pDenS25 * a34a.v() * a35a.v())
        - s4153s.v() / (pDenS14 * pDenS23 * a35a.v() * a45a.v())
        - s5134s.v() / (pDenS15 * pDenS24 * a35a.v() * a34a.v())
        + s5143s.v() / (pDenS15 * pDenS23 * a45a.v() * a34a.v())
      );
    }

    return cdouble(0,0);  
  }

  cdouble eeAAA::amp_feynman(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5){
    
    cdouble one(1,0);
    
    if(ds3>0&&ds4>0&&ds5>0){
      
      return -sqrt(2)*sqrt(2)*sqrt(2)*e*e*e*(
      - a23a.v(ds2)*s525a.v()*s414a.v()*s13s.v(ds1) + a23a.v(ds2)*s525a.v()*s434a.v()*s13s.v(ds1)
      - me*a23a.v(ds2)*s525a.v()*s34s.v()*a14a.v(ds1) + me*me*a23a.v(ds2)*s45s.v()*a45a.v()*s13s.v(ds1) 
      + me*a23a.v(ds2)*s45s.v()*s315a.v()*a14a.v(ds1) - me*s25s.v(ds2)*a35a.v()*s414a.v()*s13s.v(ds1) 
      + me*s25s.v(ds2)*a35a.v()*s434a.v()*s13s.v(ds1) - me*me*s25s.v(ds2)*a35a.v()*s34s.v()*a14a.v(ds1) 
      - me*s25s.v(ds2)*s423a.v()*a45a.v()*s13s.v(ds1) + me*s25s.v(ds2)*s453a.v()*a45a.v()*s13s.v(ds1) 
      - s25s.v(ds2)*s423a.v()*s315a.v()*a14a.v(ds1) + s25s.v(ds2)*s453a.v()*s315a.v()*a14a.v(ds1)) / (
        pDenS13 * pDenS25 * a34a.v() * a35a.v() * a45a.v()
      );
    }

    return cdouble(0,0);    
  }

  cdouble eeAAA::amp_feynman_r(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5){
    
    cdouble one(1,0);
    
    if(ds3>0&&ds4>0&&ds5>0){
      
      return -sqrt(2)*sqrt(2)*sqrt(2)*e*e*e / (a34a.v() * a35a.v() * a45a.v()) * ((
      (-s134a.v(ds1)*s34s.v() - s14s.v(ds1)*s314a.v()) * a25a.v(ds2)*s523a.v() + 
      a15a.v(ds1)*s314a.v()*(s45s.v()*s253a.v(ds2) - s24s.v(ds2)*s523a.v())) / (
        pDenS13 * pDenS25) 
      
      + ((-s135a.v(ds1) * s35s.v() - s15s.v(ds1) * s315a.v()) * a24a.v(ds2) * s423a.v() -
      a14a.v(ds1) * s315a.v() * (s45s.v() * s243a.v(ds2) + s25s.v(ds2) * s423a.v())) / (
        pDenS13 * pDenS24)

      + ((s143a.v(ds1)*s34s.v() - s13s.v(ds1)*s413a.v()) * a25a.v(ds2)*s524a.v() + 
      a15a.v(ds1)*s413a.v()*(s35s.v()*s254a.v(ds2) - s23s.v(ds2)*s524a.v())) / (
        pDenS14 * pDenS25) 

      - ((-s145a.v(ds1)*s35s.v() - s15s.v(ds1)*s415a.v()) * a23a.v(ds2)*s324a.v() - 
      a13a.v(ds1)*s415a.v()*(s35s.v()*s234a.v(ds2) + s25s.v(ds2)*s324a.v())) / (
        pDenS14 * pDenS23)

      - ((s153a.v(ds1)*s35s.v() - s13s.v(ds1)*s513a.v()) * a24a.v(ds2)*s425a.v() + 
      a14a.v(ds1)*s513a.v()*(s34s.v()*s245a.v(ds2) - s23s.v(ds2)*s425a.v())) / (
        pDenS15 * pDenS24)
      );
    }

    return cdouble(0,0);    
  }

  //set_momenta(...) must be called before amp2().
  ldouble eeAAA::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-2;j4<=2;j4+=4)
    for(int j5=-2;j5<=2;j5+=4){
	    M = amp(j1,j2,j3,j4,j5);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Symmetry factor for identical photons 1/6
    return amp2/24.0;
  }

    //set_momenta(...) must be called before amp2().
  ldouble eeAAA::amp2_feynman(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-2;j4<=2;j4+=4)
    for(int j5=-2;j5<=2;j5+=4){
	    M = amp_feynman(j1,j2,j3,j4,j5);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Symmetry factor for identical photons 1/6
    return amp2/24.0;
  }

  //set_momenta(...) must be called before amp2().
  ldouble eeAAA::amp2_permutation(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-2;j4<=2;j4+=4)
    for(int j5=-2;j5<=2;j5+=4){
	    M = amp_permutation(j1,j2,j3,j4,j5);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Symmetry factor for identical photons 1/6
    return amp2/24.0;
  }
  

  bool check_phase_space(
      const ldouble p1[4],
      const ldouble p2[4],
      const ldouble p3[4],
      const ldouble p4[4],
      const ldouble p5[4],
      ldouble m1,
      ldouble m2,
      ldouble m3,
      ldouble m4,
      ldouble m5,
      ldouble tol
      ){
      // Momentum conservation
      for(int j=0; j<4; j++){
          ldouble diff = p1[j] + p2[j] - p3[j] - p4[j] - p5[j];

          if(std::abs(diff) > tol)
            return false;
      }

      auto p_squared = [](const ldouble p[4]){
        return p[0]*p[0]
          - p[1]*p[1]
          - p[2]*p[2]
          - p[3]*p[3];
      };

      // On-shell conditions
      if(std::abs(p_squared(p1) - m1*m1) > tol)
        return false;

      if(std::abs(p_squared(p2) - m2*m2) > tol)
        return false;

      if(std::abs(p_squared(p3) - m3*m3) > tol)
        return false;

      if(std::abs(p_squared(p4) - m4*m4) > tol)
        return false;

      if(std::abs(p_squared(p5) - m5*m5) > tol)
        return false;

      return true;
  }

  void make_phase_space(
    ldouble energy,
    ldouble me,
    ldouble E3,
    ldouble theta3,
    ldouble phi3,
    ldouble theta4,
    ldouble phi4,
    ldouble p1[4],
    ldouble p2[4],
    ldouble p3[4],
    ldouble p4[4],
    ldouble p5[4]
    ){
    // Incoming e- and e+ in CM frame
    p1[0] = energy/2.0;
    p1[1] = 0;
    p1[2] = 0;
    p1[3] = std::sqrt(energy*energy/4.0 - me*me);

    p2[0] = energy/2.0;
    p2[1] = 0;
    p2[2] = 0;
    p2[3] = -p1[3];

    // Photon 3 direction
    ldouble n3[3] = {
        std::sin(theta3)*std::cos(phi3),
        std::sin(theta3)*std::sin(phi3),
        std::cos(theta3)
    };

    // Photon 4 direction
    ldouble n4[3] = {
        std::sin(theta4)*std::cos(phi4),
        std::sin(theta4)*std::sin(phi4),
        std::cos(theta4)
    };

    // n3 . n4
    ldouble c =
        n3[0]*n4[0] +
        n3[1]*n4[1] +
        n3[2]*n4[2];

    // Determine E4 from energy-momentum conservation
    ldouble E4 =
        (energy*energy - 2.0*energy*E3) /
        (2.0*(energy - E3*(1.0-c)));

    // Photon 3
    p3[0] = E3;
    p3[1] = E3*n3[0];
    p3[2] = E3*n3[1];
    p3[3] = E3*n3[2];

    // Photon 4
    p4[0] = E4;
    p4[1] = E4*n4[0];
    p4[2] = E4*n4[1];
    p4[3] = E4*n4[2];

    // Photon 5 from momentum conservation
    p5[1] = -(p3[1] + p4[1]);
    p5[2] = -(p3[2] + p4[2]);
    p5[3] = -(p3[3] + p4[3]);

    // Since photon 5 is massless, E5 = |p5|
    p5[0] = std::sqrt(
        p5[1]*p5[1] +
        p5[2]*p5[2] +
        p5[3]*p5[3]
    );
  }

  bool make_random_phase_space(
    ldouble energy,
    ldouble me,
    ldouble p1[4],
    ldouble p2[4],
    ldouble p3[4],
    ldouble p4[4],
    ldouble p5[4],
    std::mt19937& rng
    ){
    const ldouble pi = 3.14159265358979323846;

    std::uniform_real_distribution<ldouble> U(0.0, 1.0);

    // ------------------------------------------------------------
    // Incoming e- and e+ in the CM frame
    // ------------------------------------------------------------

    p1[0] = energy/2.0;
    p1[1] = 0.0;
    p1[2] = 0.0;
    p1[3] = std::sqrt(energy*energy/4.0 - me*me);

    p2[0] = energy/2.0;
    p2[1] = 0.0;
    p2[2] = 0.0;
    p2[3] = -p1[3];

    // ------------------------------------------------------------
    // Random photon 3 energy
    //
    // Keep away from the soft endpoint E3 = 0 and
    // the endpoint E3 = energy/2.
    // ------------------------------------------------------------

    ldouble E3 =
        energy * (0.05 + 0.40*U(rng));

    // ------------------------------------------------------------
    // Random direction for photon 3
    //
    // Uniform in solid angle:
    // cos(theta), not theta, is sampled uniformly.
    // ------------------------------------------------------------

    ldouble cos_theta3 = 2.0*U(rng) - 1.0;
    ldouble theta3 = std::acos(cos_theta3);
    ldouble phi3 = 2.0*pi*U(rng);

    // ------------------------------------------------------------
    // Random direction for photon 4
    // ------------------------------------------------------------

    ldouble cos_theta4 = 2.0*U(rng) - 1.0;
    ldouble theta4 = std::acos(cos_theta4);
    ldouble phi4 = 2.0*pi*U(rng);

    // ------------------------------------------------------------
    // Construct photon directions
    // ------------------------------------------------------------

    ldouble n3[3] = {
        std::sin(theta3)*std::cos(phi3),
        std::sin(theta3)*std::sin(phi3),
        std::cos(theta3)
    };

    ldouble n4[3] = {
        std::sin(theta4)*std::cos(phi4),
        std::sin(theta4)*std::sin(phi4),
        std::cos(theta4)
    };

    // Dot product n3 . n4
    ldouble c =
        n3[0]*n4[0] +
        n3[1]*n4[1] +
        n3[2]*n4[2];

    // ------------------------------------------------------------
    // Solve energy conservation for E4
    // ------------------------------------------------------------

    ldouble denominator =
        2.0*(energy - E3*(1.0-c));

    if(std::abs(denominator) < 1e-14)
        return false;

    ldouble E4 =
        (energy*energy - 2.0*energy*E3)
        / denominator;

    // Reject unphysical energies
    if(!std::isfinite(E4) || E4 <= 0.0)
        return false;

    // ------------------------------------------------------------
    // Photon 3
    // ------------------------------------------------------------

    p3[0] = E3;
    p3[1] = E3*n3[0];
    p3[2] = E3*n3[1];
    p3[3] = E3*n3[2];

    // ------------------------------------------------------------
    // Photon 4
    // ------------------------------------------------------------

    p4[0] = E4;
    p4[1] = E4*n4[0];
    p4[2] = E4*n4[1];
    p4[3] = E4*n4[2];

    // ------------------------------------------------------------
    // Photon 5 from spatial momentum conservation
    // ------------------------------------------------------------

    p5[1] = -(p3[1] + p4[1]);
    p5[2] = -(p3[2] + p4[2]);
    p5[3] = -(p3[3] + p4[3]);

    p5[0] = std::sqrt(
        p5[1]*p5[1] +
        p5[2]*p5[2] +
        p5[3]*p5[3]
    );

    // ------------------------------------------------------------
    // Check all components are finite
    // ------------------------------------------------------------

    for(int i=0; i<4; i++){
        if(!std::isfinite(p1[i])) return false;
        if(!std::isfinite(p2[i])) return false;
        if(!std::isfinite(p3[i])) return false;
        if(!std::isfinite(p4[i])) return false;
        if(!std::isfinite(p5[i])) return false;
    }

    return true;
}

  //  Tests
  int test_eeAAA(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> A , A , A   :";
    if (1 == 2) {//amp^2
      int i=0;
      // me=0.0005, pspatial=250
      ldouble me=0.0005;
      ldouble EE=0.31333;
      eeAAA eeAAAAmp = eeAAA(EE,me);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.761473098865852E-01,1.196559098890515E-01,6.884618493553914E-02,4.748300512537967E-02,3.599742458224401E-02,2.906647080620063E-02,2.465909507168581E-02,2.184718935306326E-02,2.016436086672569E-02,1.937355800660445E-02,1.937355800660445E-02,2.016436086672569E-02,2.184718935306325E-02,2.465909507168581E-02,2.906647080620063E-02,3.599742458224402E-02,4.748300512537965E-02,6.884618493553911E-02,1.196559098890515E-01,3.761473098865843E-01};
      i += eeAAAAmp.test_2to2_amp2([&]() { return eeAAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eeAAAAmp.test_2to2_amp2_rotations([&]() { return eeAAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eeAAAAmp.test_2to2_amp2_boosts([&]() { return eeAAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eeAAAAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      //Close to threshold
      pspatial = 0.0001;
      ldouble dataCH2[20] = {2.081251342993773E-02,2.079750983522615E-02,2.078115831044515E-02,2.076461266895772E-02,2.074882438042710E-02,2.073456777596373E-02,2.072246020291247E-02,2.071297785405228E-02,2.070646782725441E-02,2.070315683100089E-02,2.070315683100089E-02,2.070646782725441E-02,2.071297785405228E-02,2.072246020291247E-02,2.073456777596374E-02,2.074882438042710E-02,2.076461266895772E-02,2.078115831044515E-02,2.079750983522615E-02,2.081251342993773E-02};
      i += eeAAAAmp.test_2to2_amp2([&]() { return eeAAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eeAAAAmp.test_2to2_amp2_rotations([&]() { return eeAAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eeAAAAmp.test_2to2_amp2_boosts([&]() { return eeAAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eeAAAAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeAAAAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }

    {
      ldouble me = 0.0005;
      ldouble EE = 0.31333;

      eeAAA eeAAAAmp = eeAAA(EE, me);

      ldouble energy = 300.0;

      ldouble p1[4], p2[4], p3[4], p4[4], p5[4];

      // Fixed seed makes the test reproducible.
      std::mt19937 rng(12345);

      const int Npoints = 10;

      int point = 0;
      int attempts = 0;

      while(point < Npoints){

        attempts++;

        bool generated = make_random_phase_space(energy, me, p1, p2, p3, p4, p5, rng);

        if(!generated)
          continue;

        // Check phase space
        bool phase_space_ok = check_phase_space(p1, p2, p3, p4, p5, me, me, 0, 0, 0, 1e-10);

        if(!phase_space_ok)
          continue;

        point++;

        std::cout << "\nPhase-space point " << point << "\n";

        std::cout << "  Phase space: PASS\n";

        eeAAAAmp.set_momenta(p1, p2, p3, p4, p5);

        cdouble amp_x = eeAAAAmp.amp(1,1,2,2,2);
        // cdouble amp_f = eeAAAAmp.amp_feynman(1,1,2,2,2);
        cdouble amp_p = eeAAAAmp.amp_permutation(1,1,2,2,2);
        cdouble amp_fr = eeAAAAmp.amp_feynman_r(1,1,2,2,2);


        // std::cout << "  Feynman     = " << amp_f << "\n";
        std::cout << "  Reduced     = " << amp_fr << "\n";
        std::cout << "  x-factor    = " << amp_x << "\n";
        std::cout << "  Permutation = " << amp_p << "\n";
        // std::cout << "  |Feynman - Reduced| = " << std::abs(amp_f - amp_fr) << "\n";
      }

      std::cout << "\nGenerated " << Npoints << " valid phase-space points after " << attempts << " attempts.\n";
    }

    return n;
  }


}
