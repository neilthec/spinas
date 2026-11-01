
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
#include <algorithm>
#include <iomanip>
#include <limits>

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
    s154a = sproduct(SQUARE,&p1,&p5,&p4);
    s234a = sproduct(SQUARE,&p2,&p3,&p4);
    s235a = sproduct(SQUARE,&p2,&p3,&p5);
    s243a = sproduct(SQUARE,&p2,&p4,&p3);
    s245a = sproduct(SQUARE,&p2,&p4,&p5);
    s253a = sproduct(SQUARE,&p2,&p5,&p3);
    s254a = sproduct(SQUARE,&p2,&p5,&p4);
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s313a = sproduct(SQUARE,&p3,&p1,&p3);
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    s315a = sproduct(SQUARE,&p3,&p1,&p5);
    s321a = sproduct(SQUARE,&p3,&p2,&p1);
    s323a = sproduct(SQUARE,&p3,&p2,&p3);
    s324a = sproduct(SQUARE,&p3,&p2,&p4);
    s325a = sproduct(SQUARE,&p3,&p2,&p5);
    s343a = sproduct(SQUARE,&p3,&p4,&p3);
    s345a = sproduct(SQUARE,&p3,&p4,&p5);
    s353a = sproduct(SQUARE,&p3,&p5,&p3);
    s354a = sproduct(SQUARE,&p3,&p5,&p4);
    s412a = sproduct(SQUARE,&p4,&p1,&p2);
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
    s414a = sproduct(SQUARE,&p4,&p1,&p4);
    s415a = sproduct(SQUARE,&p4,&p1,&p5);
    s421a = sproduct(SQUARE,&p4,&p2,&p1);
    s423a = sproduct(SQUARE,&p4,&p2,&p3);
    s424a = sproduct(SQUARE,&p4,&p2,&p4);
    s425a = sproduct(SQUARE,&p4,&p2,&p5);
    s434a = sproduct(SQUARE,&p4,&p3,&p4);
    s435a = sproduct(SQUARE,&p4,&p3,&p5);
    s453a = sproduct(SQUARE,&p4,&p5,&p3);
    s454a = sproduct(SQUARE,&p4,&p5,&p4);
    s512a = sproduct(SQUARE,&p5,&p1,&p2);
    s513a = sproduct(SQUARE,&p5,&p1,&p3);
    s514a = sproduct(SQUARE,&p5,&p1,&p4);
    s515a = sproduct(SQUARE,&p5,&p1,&p5);
    s521a = sproduct(SQUARE,&p5,&p2,&p1);
    s523a = sproduct(SQUARE,&p5,&p2,&p3);
    s524a = sproduct(SQUARE,&p5,&p2,&p4);
    s525a = sproduct(SQUARE,&p5,&p2,&p5);
    s534a = sproduct(SQUARE,&p5,&p3,&p4);
    s535a = sproduct(SQUARE,&p5,&p3,&p5);
    s543a = sproduct(SQUARE,&p5,&p4,&p3);
    s545a = sproduct(SQUARE,&p5,&p4,&p5);
    s3123s = sproduct(SQUARE,&p3,&p1,&p2,&p3);
    a3123a = sproduct(ANGLE,&p3,&p1,&p2,&p3);
    s3124s = sproduct(SQUARE,&p3,&p1,&p2,&p4);
    s3125s = sproduct(SQUARE,&p3,&p1,&p2,&p5);
    s3134s = sproduct(SQUARE,&p3,&p1,&p3,&p4);
    s3135s = sproduct(SQUARE,&p3,&p1,&p3,&p5);
    s3214s = sproduct(SQUARE,&p3,&p2,&p1,&p4);
    s3215s = sproduct(SQUARE,&p3,&p2,&p1,&p5);
    s3424s = sproduct(SQUARE,&p3,&p4,&p2,&p4);
    s3525s = sproduct(SQUARE,&p3,&p5,&p2,&p5);
    s4123s = sproduct(SQUARE,&p4,&p1,&p2,&p3);
    s4124s = sproduct(SQUARE,&p4,&p1,&p2,&p4);
    a4124a = sproduct(ANGLE,&p4,&p1,&p2,&p4);
    s4125s = sproduct(SQUARE,&p4,&p1,&p2,&p5);
    s4143s = sproduct(SQUARE,&p4,&p1,&p4,&p3);
    s4145s = sproduct(SQUARE,&p4,&p1,&p4,&p5);
    s4213s = sproduct(SQUARE,&p4,&p2,&p1,&p3);
    s4215s = sproduct(SQUARE,&p4,&p2,&p1,&p5);
    s4323s = sproduct(SQUARE,&p4,&p3,&p2,&p3);
    s4525s = sproduct(SQUARE,&p4,&p5,&p2,&p5);
    s5124s = sproduct(SQUARE,&p5,&p1,&p2,&p4);
    s5125s = sproduct(SQUARE,&p5,&p1,&p2,&p5);
    a5125a = sproduct(ANGLE,&p5,&p1,&p2,&p5);
    s5154s = sproduct(SQUARE,&p5,&p1,&p5,&p4);
    s3145s = sproduct(SQUARE,&p3,&p1,&p4,&p5);
    s3154s = sproduct(SQUARE,&p3,&p1,&p5,&p4);
    s4135s = sproduct(SQUARE,&p4,&p1,&p3,&p5);
    s4153s = sproduct(SQUARE,&p4,&p1,&p5,&p3);
    s5123s = sproduct(SQUARE,&p5,&p1,&p2,&p3);
    s5134s = sproduct(SQUARE,&p5,&p1,&p3,&p4);
    s5153s = sproduct(SQUARE,&p5,&p1,&p5,&p3);
    s5213s = sproduct(SQUARE,&p5,&p2,&p1,&p3);
    s5214s = sproduct(SQUARE,&p5,&p2,&p1,&p4);
    s5143s = sproduct(SQUARE,&p5,&p1,&p4,&p3);
    s5323s = sproduct(SQUARE,&p5,&p3,&p2,&p3);
    s5424s = sproduct(SQUARE,&p5,&p4,&p2,&p4);
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
    s154a.update();
    s234a.update();
    s235a.update();
    s243a.update();
    s245a.update();
    s253a.update();
    s254a.update();
    s312a.update();
    s313a.update();
    s314a.update();
    s315a.update();
    s321a.update();
    s323a.update();
    s324a.update();
    s325a.update();
    s343a.update();
    s345a.update();
    s353a.update();
    s354a.update();
    s412a.update();
    s413a.update();
    s414a.update();
    s415a.update();
    s421a.update();
    s423a.update();
    s424a.update();
    s425a.update();
    s434a.update();
    s435a.update();
    s453a.update();
    s454a.update();
    s512a.update();
    s513a.update();
    s514a.update();
    s515a.update();
    s521a.update();
    s523a.update();
    s524a.update();
    s525a.update();
    s534a.update();
    s535a.update();
    s543a.update();
    s545a.update();
    s3123s.update();
    s3124s.update();
    a3123a.update();
    s3125s.update();
    s3134s.update();
    s3135s.update();
    s3214s.update();
    s3215s.update();
    s3424s.update();
    s3525s.update();
    s4123s.update();
    s4124s.update();
    s4125s.update();
    a4124a.update();
    s4143s.update();
    s4145s.update();
    s4213s.update();
    s4215s.update();
    s4323s.update();
    s5123s.update();
    s5125s.update();
    a5125a.update();
    s3145s.update();
    s3154s.update();
    s4135s.update();
    s4153s.update();
    s4525s.update();
    s5124s.update();
    s5134s.update();
    s5143s.update();
    s5153s.update();
    s5213s.update();
    s5214s.update();
    s5154s.update();
    s5323s.update();
    s5424s.update();
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

  // H in the PDF, with only electron propagator poles. Photons are outgoing.
  cdouble eeAAA::amp_nonspurious(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5){
    if(ds3>0&&ds4>0&&ds5>0){
      return 2*sqrt2*e*e*e*me*a12a.v(ds1,ds2)
        / (pDenS13*pDenS14*pDenS15) * (
          s45s.v()*s45s.v()*s3123s.v()/pDenS23
          + s35s.v()*s35s.v()*s4124s.v()/pDenS24
          + s34s.v()*s34s.v()*s5125s.v()/pDenS25
        );
    }
    return cdouble(0,0);
  }

  // S2 in the PDF, equivalent to the full permutation sum below.
  // The signs include crossing each photon angle spinor to an outgoing photon.
  cdouble eeAAA::amp_spurious(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5){
    if(ds3>0&&ds4>0&&ds5>0){
      return 2*sqrt2*e*e*e*me*a12a.v(ds1,ds2) * (
          - s34s.v()*s523a.v()/(a34a.v()*a35a.v()*pDenS13*pDenS25)
          - s35s.v()*s423a.v()/(a35a.v()*a34a.v()*pDenS13*pDenS24)
          - s34s.v()*s524a.v()/(a34a.v()*a45a.v()*pDenS14*pDenS25)
          + s45s.v()*s324a.v()/(a45a.v()*a34a.v()*pDenS14*pDenS23)
          + s35s.v()*s425a.v()/(a35a.v()*a45a.v()*pDenS15*pDenS24)
          + s45s.v()*s325a.v()/(a45a.v()*a35a.v()*pDenS15*pDenS23)
        );
    }
    return cdouble(0,0);
  }

  // Full total-photon-momentum numerator [a|p1 (k3+k4+k5)|c].
  // The kc insertion vanishes at the massless endpoint; the ka insertion must remain off shell.
  cdouble eeAAA::amp_permutation(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5) {
    
    if(ds3>0&&ds4>0&&ds5>0){

      return -sqrt(2)*sqrt(2)*sqrt(2)*e*e*e*me*a12a.v(ds1,ds2)*(
        (s3145s.v() + s3135s.v()) / (pDenS13 * pDenS25 * a34a.v() * a45a.v())
        - (s3154s.v() + s3134s.v()) / (pDenS13 * pDenS24 * a35a.v() * a45a.v())
        - (s4135s.v() + s4145s.v()) / (pDenS14 * pDenS25 * a34a.v() * a35a.v())
        - (s4153s.v() + s4143s.v()) / (pDenS14 * pDenS23 * a35a.v() * a45a.v())
        - (s5134s.v() + s5154s.v()) / (pDenS15 * pDenS24 * a35a.v() * a34a.v())
        + (s5143s.v() + s5153s.v()) / (pDenS15 * pDenS23 * a45a.v() * a34a.v())
      );
    }

    return cdouble(0,0);  
  }

  // q3=|4>, q4=|5>, q5=|3>; fixed for all six orders.
  cdouble eeAAA::amp_feynman(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5) {
    if(ds3<=0 || ds4<=0 || ds5<=0) return cdouble(0,0);
    sproduct f_s31s(SQUARE,&p3,&p1);
    sproduct f_a41a(ANGLE,&p4,&p1);
    sproduct f_a513s(ANGLE,&p5,&p1,&p3);
    sproduct f_a324s(ANGLE,&p3,&p2,&p4);
    sproduct f_a354s(ANGLE,&p3,&p5,&p4);
    sproduct f_s54s(SQUARE,&p5,&p4);
    sproduct f_s43s(SQUARE,&p4,&p3);
    sproduct f_a54a(ANGLE,&p5,&p4);
    sproduct f_a313s(ANGLE,&p3,&p1,&p3);
    sproduct f_a525s(ANGLE,&p5,&p2,&p5);
    sproduct f_a545s(ANGLE,&p5,&p4,&p5);
    sproduct f_s53s(SQUARE,&p5,&p3);
    sproduct f_a53a(ANGLE,&p5,&p3);
    sproduct f_s41s(SQUARE,&p4,&p1);
    sproduct f_a51a(ANGLE,&p5,&p1);
    sproduct f_a414s(ANGLE,&p4,&p1,&p4);
    sproduct f_a323s(ANGLE,&p3,&p2,&p3);
    sproduct f_a353s(ANGLE,&p3,&p5,&p3);
    sproduct f_a314s(ANGLE,&p3,&p1,&p4);
    sproduct f_a425s(ANGLE,&p4,&p2,&p5);
    sproduct f_a435s(ANGLE,&p4,&p3,&p5);
    sproduct f_a43a(ANGLE,&p4,&p3);
    sproduct f_s51s(SQUARE,&p5,&p1);
    sproduct f_a31a(ANGLE,&p3,&p1);
    sproduct f_a415s(ANGLE,&p4,&p1,&p5);
    sproduct f_a523s(ANGLE,&p5,&p2,&p3);
    sproduct f_a543s(ANGLE,&p5,&p4,&p3);
    sproduct f_a515s(ANGLE,&p5,&p1,&p5);
    sproduct f_a424s(ANGLE,&p4,&p2,&p4);
    sproduct f_a434s(ANGLE,&p4,&p3,&p4);
    const cdouble reference_denominator=f_a43a.v()*f_a54a.v()*a35a.v();
    // 345: eight terms from the two off-shell fermion numerators.
    const cdouble M345=(
      - a23a.v(ds2)*s525a.v()*(s414a.v() - s434a.v())*f_s31s.v(ds1)
      + me*a23a.v(ds2)*s525a.v()*f_s43s.v()*f_a41a.v(ds1)
      + me*me*a23a.v(ds2)*f_s54s.v()*f_a54a.v()*f_s31s.v(ds1)
      - me*a23a.v(ds2)*f_s54s.v()*f_a513s.v()*f_a41a.v(ds1)
      - me*s25s.v(ds2)*a35a.v()*(s414a.v() - s434a.v())*f_s31s.v(ds1)
      + me*me*s25s.v(ds2)*a35a.v()*f_s43s.v()*f_a41a.v(ds1)
      + me*s25s.v(ds2)*(f_a324s.v() - f_a354s.v())*f_a54a.v()*f_s31s.v(ds1)
      - s25s.v(ds2)*(f_a324s.v() - f_a354s.v())*f_a513s.v()*f_a41a.v(ds1)
    )/(pDenS13*pDenS25);

    // 354: eight terms from the two off-shell fermion numerators.
    const cdouble M354=(
      - a25a.v(ds2)*s423a.v()*(s514a.v() - s534a.v())*f_s31s.v(ds1)
      + me*a25a.v(ds2)*s423a.v()*f_s53s.v()*f_a41a.v(ds1)
      + me*me*a25a.v(ds2)*s45s.v()*a34a.v()*f_s31s.v(ds1)
      - me*a25a.v(ds2)*s45s.v()*f_a313s.v()*f_a41a.v(ds1)
      - me*s24s.v(ds2)*f_a53a.v()*(s514a.v() - s534a.v())*f_s31s.v(ds1)
      + me*me*s24s.v(ds2)*f_a53a.v()*f_s53s.v()*f_a41a.v(ds1)
      + me*s24s.v(ds2)*(f_a525s.v() - f_a545s.v())*a34a.v()*f_s31s.v(ds1)
      - s24s.v(ds2)*(f_a525s.v() - f_a545s.v())*f_a313s.v()*f_a41a.v(ds1)
    )/(pDenS13*pDenS24);

    // 435: eight terms from the two off-shell fermion numerators.
    const cdouble M435=(
      - a23a.v(ds2)*s524a.v()*(s315a.v() - s345a.v())*f_s41s.v(ds1)
      + me*a23a.v(ds2)*s524a.v()*s34s.v()*f_a51a.v(ds1)
      + me*me*a23a.v(ds2)*f_s53s.v()*a45a.v()*f_s41s.v(ds1)
      - me*a23a.v(ds2)*f_s53s.v()*f_a414s.v()*f_a51a.v(ds1)
      - me*s25s.v(ds2)*a34a.v()*(s315a.v() - s345a.v())*f_s41s.v(ds1)
      + me*me*s25s.v(ds2)*a34a.v()*s34s.v()*f_a51a.v(ds1)
      + me*s25s.v(ds2)*(f_a323s.v() - f_a353s.v())*a45a.v()*f_s41s.v(ds1)
      - s25s.v(ds2)*(f_a323s.v() - f_a353s.v())*f_a414s.v()*f_a51a.v(ds1)
    )/(pDenS14*pDenS25);

    // 453: eight terms from the two off-shell fermion numerators.
    const cdouble M453=(
      - a24a.v(ds2)*s323a.v()*(s515a.v() - s545a.v())*f_s41s.v(ds1)
      + me*a24a.v(ds2)*s323a.v()*f_s54s.v()*f_a51a.v(ds1)
      + me*me*a24a.v(ds2)*s35s.v()*a35a.v()*f_s41s.v(ds1)
      - me*a24a.v(ds2)*s35s.v()*f_a314s.v()*f_a51a.v(ds1)
      - me*s23s.v(ds2)*f_a43a.v()*(s515a.v() - s545a.v())*f_s41s.v(ds1)
      + me*me*s23s.v(ds2)*f_a43a.v()*f_s54s.v()*f_a51a.v(ds1)
      + me*s23s.v(ds2)*(f_a425s.v() - f_a435s.v())*a35a.v()*f_s41s.v(ds1)
      - s23s.v(ds2)*(f_a425s.v() - f_a435s.v())*f_a314s.v()*f_a51a.v(ds1)
    )/(pDenS14*pDenS23);

    // 534: eight terms from the two off-shell fermion numerators.
    const cdouble M534=(
      - a25a.v(ds2)*s424a.v()*(s313a.v() - s353a.v())*f_s51s.v(ds1)
      + me*a25a.v(ds2)*s424a.v()*s35s.v()*f_a31a.v(ds1)
      + me*me*a25a.v(ds2)*f_s43s.v()*f_a43a.v()*f_s51s.v(ds1)
      - me*a25a.v(ds2)*f_s43s.v()*f_a415s.v()*f_a31a.v(ds1)
      - me*s24s.v(ds2)*f_a54a.v()*(s313a.v() - s353a.v())*f_s51s.v(ds1)
      + me*me*s24s.v(ds2)*f_a54a.v()*s35s.v()*f_a31a.v(ds1)
      + me*s24s.v(ds2)*(f_a523s.v() - f_a543s.v())*f_a43a.v()*f_s51s.v(ds1)
      - s24s.v(ds2)*(f_a523s.v() - f_a543s.v())*f_a415s.v()*f_a31a.v(ds1)
    )/(pDenS15*pDenS24);

    // 543: eight terms from the two off-shell fermion numerators.
    const cdouble M543=(
      - a24a.v(ds2)*s325a.v()*(s413a.v() - s453a.v())*f_s51s.v(ds1)
      + me*a24a.v(ds2)*s325a.v()*s45s.v()*f_a31a.v(ds1)
      + me*me*a24a.v(ds2)*s34s.v()*f_a53a.v()*f_s51s.v(ds1)
      - me*a24a.v(ds2)*s34s.v()*f_a515s.v()*f_a31a.v(ds1)
      - me*s23s.v(ds2)*a45a.v()*(s413a.v() - s453a.v())*f_s51s.v(ds1)
      + me*me*s23s.v(ds2)*a45a.v()*s45s.v()*f_a31a.v(ds1)
      + me*s23s.v(ds2)*(f_a424s.v() - f_a434s.v())*f_a53a.v()*f_s51s.v(ds1)
      - s23s.v(ds2)*(f_a424s.v() - f_a434s.v())*f_a515s.v()*f_a31a.v(ds1)
    )/(pDenS15*pDenS23);
    return -2*sqrt2*e*e*e*(M345+M354+M435+M453+M534+M543)/reference_denominator;
  }

  // Compatibility name for the same explicit cyclic eight-term expression.
  cdouble eeAAA::amp_feynman_r(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5) {
    return amp_feynman(ds1,ds2,ds3,ds4,ds5);
  }

  // <q3|=[3|p2, <q4|=[4|p1, <q5|=[5|p1; fixed for all orders.
  cdouble eeAAA::amp_feynman_x(const int& ds1, const int& ds2, const int& ds3, const int& ds4, const int& ds5) {
    if(ds3<=0 || ds4<=0 || ds5<=0) return cdouble(0,0);
    sproduct f_s31s(SQUARE,&p3,&p1);
    sproduct f_s4113s(SQUARE,&p4,&p1,&p1,&p3);
    sproduct f_s54s(SQUARE,&p5,&p4);
    sproduct f_s43s(SQUARE,&p4,&p3);
    sproduct f_s5114s(SQUARE,&p5,&p1,&p1,&p4);
    sproduct f_s5113s(SQUARE,&p5,&p1,&p1,&p3);
    sproduct f_s53s(SQUARE,&p5,&p3);
    sproduct f_s4115s(SQUARE,&p4,&p1,&p1,&p5);
    sproduct f_s41s(SQUARE,&p4,&p1);
    sproduct f_s411a(SQUARE,&p4,&p1,&p1);
    sproduct f_s5223s(SQUARE,&p5,&p2,&p2,&p3);
    sproduct f_s3114s(SQUARE,&p3,&p1,&p1,&p4);
    sproduct f_s3414s(SQUARE,&p3,&p4,&p1,&p4);
    sproduct f_s322a(SQUARE,&p3,&p2,&p2);
    sproduct f_s5414s(SQUARE,&p5,&p4,&p1,&p4);
    sproduct f_s3225s(SQUARE,&p3,&p2,&p2,&p5);
    sproduct f_s3235s(SQUARE,&p3,&p2,&p3,&p5);
    sproduct f_s51s(SQUARE,&p5,&p1);
    sproduct f_s511a(SQUARE,&p5,&p1,&p1);
    sproduct f_s4223s(SQUARE,&p4,&p2,&p2,&p3);
    sproduct f_s3115s(SQUARE,&p3,&p1,&p1,&p5);
    sproduct f_s3515s(SQUARE,&p3,&p5,&p1,&p5);
    sproduct f_s4515s(SQUARE,&p4,&p5,&p1,&p5);
    sproduct f_s3224s(SQUARE,&p3,&p2,&p2,&p4);
    sproduct f_s3234s(SQUARE,&p3,&p2,&p3,&p4);
    const cdouble reference_denominator=s323a.v()*s414a.v()*s515a.v();
    // 345: eight terms from the two off-shell fermion numerators.
    const cdouble M345=(
      + s512a.v(ds2)*s5214s.v()*(s4123s.v() - s4323s.v())*f_s31s.v(ds1)
      + me*s512a.v(ds2)*s5214s.v()*f_s43s.v()*s321a.v(ds1)
      + me*me*s512a.v(ds2)*f_s54s.v()*s4123s.v()*f_s31s.v(ds1)
      + me*s512a.v(ds2)*f_s54s.v()*f_s4113s.v()*s321a.v(ds1)
      - me*s25s.v(ds2)*f_s5114s.v()*(s4123s.v() - s4323s.v())*f_s31s.v(ds1)
      - me*me*s25s.v(ds2)*f_s5114s.v()*f_s43s.v()*s321a.v(ds1)
      - me*s25s.v(ds2)*(s5124s.v() - s5154s.v())*s4123s.v()*f_s31s.v(ds1)
      - s25s.v(ds2)*(s5124s.v() - s5154s.v())*f_s4113s.v()*s321a.v(ds1)
    )/(pDenS13*pDenS25);

    // 354: eight terms from the two off-shell fermion numerators.
    const cdouble M354=(
      + s412a.v(ds2)*s4215s.v()*(s5123s.v() - s5323s.v())*f_s31s.v(ds1)
      + me*s412a.v(ds2)*s4215s.v()*f_s53s.v()*s321a.v(ds1)
      + me*me*s412a.v(ds2)*s45s.v()*s5123s.v()*f_s31s.v(ds1)
      + me*s412a.v(ds2)*s45s.v()*f_s5113s.v()*s321a.v(ds1)
      - me*s24s.v(ds2)*f_s4115s.v()*(s5123s.v() - s5323s.v())*f_s31s.v(ds1)
      - me*me*s24s.v(ds2)*f_s4115s.v()*f_s53s.v()*s321a.v(ds1)
      - me*s24s.v(ds2)*(s4125s.v() - s4145s.v())*s5123s.v()*f_s31s.v(ds1)
      - s24s.v(ds2)*(s4125s.v() - s4145s.v())*f_s5113s.v()*s321a.v(ds1)
    )/(pDenS13*pDenS24);

    // 435: eight terms from the two off-shell fermion numerators.
    const cdouble M435=(
      + s512a.v(ds2)*f_s5223s.v()*(f_s3114s.v() - f_s3414s.v())*f_s41s.v(ds1)
      + me*s512a.v(ds2)*f_s5223s.v()*s34s.v()*f_s411a.v(ds1)
      + me*me*s512a.v(ds2)*f_s53s.v()*s3214s.v()*f_s41s.v(ds1)
      + me*s512a.v(ds2)*f_s53s.v()*s3214s.v()*f_s411a.v(ds1)
      - me*s25s.v(ds2)*s5123s.v()*(f_s3114s.v() - f_s3414s.v())*f_s41s.v(ds1)
      - me*me*s25s.v(ds2)*s5123s.v()*s34s.v()*f_s411a.v(ds1)
      - me*s25s.v(ds2)*(s5123s.v() - s5153s.v())*s3214s.v()*f_s41s.v(ds1)
      - s25s.v(ds2)*(s5123s.v() - s5153s.v())*s3214s.v()*f_s411a.v(ds1)
    )/(pDenS14*pDenS25);

    // 453: eight terms from the two off-shell fermion numerators.
    const cdouble M453=(
      + f_s322a.v(ds2)*s3215s.v()*(f_s5114s.v() - f_s5414s.v())*f_s41s.v(ds1)
      + me*f_s322a.v(ds2)*s3215s.v()*f_s54s.v()*f_s411a.v(ds1)
      + me*me*f_s322a.v(ds2)*s35s.v()*f_s5114s.v()*f_s41s.v(ds1)
      + me*f_s322a.v(ds2)*s35s.v()*f_s5114s.v()*f_s411a.v(ds1)
      - me*s23s.v(ds2)*s3215s.v()*(f_s5114s.v() - f_s5414s.v())*f_s41s.v(ds1)
      - me*me*s23s.v(ds2)*s3215s.v()*f_s54s.v()*f_s411a.v(ds1)
      - me*s23s.v(ds2)*(f_s3225s.v() - f_s3235s.v())*f_s5114s.v()*f_s41s.v(ds1)
      - s23s.v(ds2)*(f_s3225s.v() - f_s3235s.v())*f_s5114s.v()*f_s411a.v(ds1)
    )/(pDenS14*pDenS23);

    // 534: eight terms from the two off-shell fermion numerators.
    const cdouble M534=(
      + s412a.v(ds2)*f_s4223s.v()*(f_s3115s.v() - f_s3515s.v())*f_s51s.v(ds1)
      + me*s412a.v(ds2)*f_s4223s.v()*s35s.v()*f_s511a.v(ds1)
      + me*me*s412a.v(ds2)*f_s43s.v()*s3215s.v()*f_s51s.v(ds1)
      + me*s412a.v(ds2)*f_s43s.v()*s3215s.v()*f_s511a.v(ds1)
      - me*s24s.v(ds2)*s4123s.v()*(f_s3115s.v() - f_s3515s.v())*f_s51s.v(ds1)
      - me*me*s24s.v(ds2)*s4123s.v()*s35s.v()*f_s511a.v(ds1)
      - me*s24s.v(ds2)*(s4123s.v() - s4143s.v())*s3215s.v()*f_s51s.v(ds1)
      - s24s.v(ds2)*(s4123s.v() - s4143s.v())*s3215s.v()*f_s511a.v(ds1)
    )/(pDenS15*pDenS24);

    // 543: eight terms from the two off-shell fermion numerators.
    const cdouble M543=(
      + f_s322a.v(ds2)*s3214s.v()*(f_s4115s.v() - f_s4515s.v())*f_s51s.v(ds1)
      + me*f_s322a.v(ds2)*s3214s.v()*s45s.v()*f_s511a.v(ds1)
      + me*me*f_s322a.v(ds2)*s34s.v()*f_s4115s.v()*f_s51s.v(ds1)
      + me*f_s322a.v(ds2)*s34s.v()*f_s4115s.v()*f_s511a.v(ds1)
      - me*s23s.v(ds2)*s3214s.v()*(f_s4115s.v() - f_s4515s.v())*f_s51s.v(ds1)
      - me*me*s23s.v(ds2)*s3214s.v()*s45s.v()*f_s511a.v(ds1)
      - me*s23s.v(ds2)*(f_s3224s.v() - f_s3234s.v())*f_s4115s.v()*f_s51s.v(ds1)
      - s23s.v(ds2)*(f_s3224s.v() - f_s3234s.v())*f_s4115s.v()*f_s511a.v(ds1)
    )/(pDenS15*pDenS23);
    return -2*sqrt2*e*e*e*(M345+M354+M435+M453+M534+M543)/reference_denominator;
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

  ldouble eeAAA::amp2(const ldouble mom1[4], const ldouble mom2[4],
                      const ldouble mom3[4], const ldouble mom4[4],
                      const ldouble mom5[4]) {
    set_momenta(mom1,mom2,mom3,mom4,mom5);
    return amp2();
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
  

  bool eeAAA::check_phase_space(
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

  void eeAAA::make_phase_space(
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

  bool eeAAA::make_random_phase_space(
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
    const ldouble relative_tolerance=1e-9L;
    auto relative_error=[](const cdouble& a, const cdouble& b){
      const ldouble scale=std::max(std::abs(a),std::abs(b));
      return scale>0 ? std::abs(a-b)/scale : ldouble(0);
    };
    auto agrees=[&](const cdouble& a, const cdouble& b){
      return std::isfinite(a.real()) && std::isfinite(a.imag())
        && std::isfinite(b.real()) && std::isfinite(b.imag())
        && relative_error(a,b)<=relative_tolerance;
    };
    auto print_reference_comparison=[&](eeAAA& amplitude){
      const cdouble h=amplitude.amp_nonspurious(1,1,2,2,2);
      const cdouble s2=amplitude.amp_spurious(1,1,2,2,2);
      const bool pass=agrees(h,s2);
      const std::streamsize precision=std::cout.precision();
      std::cout << std::setprecision(12)
        << "  H  (nonspurious) = " << h << "\n"
        << "  S2 (spurious)    = " << s2 << "\n"
        << "  H = S2: " << (pass ? "PASS" : "FAIL")
        << "; relative complex difference = " << relative_error(h,s2) << "\n";
      std::cout.precision(precision);
      if(!pass) n++;
    };
    auto print_ratio=[](const char* label,const cdouble& a,const cdouble& b){
      std::cout << "  " << label << " = ";
      if(std::abs(b)==0) std::cout << "undefined (zero denominator)";
      else std::cout << a/b << "; magnitude = " << std::abs(a/b);
      std::cout << "\n";
    };
    std::cout<<"\t* e , E  -> A , A , A   :";
    {
      eeAAA amplitude(0.31333L,1.L);
      ldouble p[5][4];
      const ldouble pi=std::acos(-1.L);
      amplitude.make_phase_space(6.L,1.L,2.L,pi/2,pi,pi/2,pi/3,
        p[0],p[1],p[2],p[3],p[4]);
      bool pass=amplitude.check_phase_space(p[0],p[1],p[2],p[3],p[4],1,1,0,0,0,1e-10L);
      std::mt19937 rng(12345);
      pass=eeAAA::make_random_phase_space(6.L,1.L,p[0],p[1],p[2],p[3],p[4],rng) && pass;
      const ldouble squared=amplitude.amp2(p[0],p[1],p[2],p[3],p[4]);
      pass=std::isfinite(squared) && squared==amplitude.amp2() && pass;
      if(!pass) n++;
      std::cout << "\n  Declared phase-space helpers and momentum amp2 overload: "
                << (pass ? "PASS" : "FAIL") << "\n";
    }

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

      const ldouble pi = 3.14159265358979323846;

      ldouble p1[4], p2[4], p3[4], p4[4], p5[4];

      p1[0] = energy/2.0;
      p1[1] = 0.0;
      p1[2] = 0.0;
      p1[3] = std::sqrt(energy*energy/4.0 - me*me);

      p2[0] = energy/2.0;
      p2[1] = 0.0;
      p2[2] = 0.0;
      p2[3] = -p1[3];

      p3[0] = energy/3.0;
      p3[1] = -energy/3.0;
      p3[2] = 0;
      p3[3] = 0;

      p4[0] = energy/3.0;
      p4[1] = energy/3.0 * std::cos(pi/3.0);
      p4[2] = energy/3.0 * std::sin(pi/3.0);
      p4[3] = 0;

      p5[0] = energy/3.0;
      p5[1] = energy/3.0 * std::cos(pi/3.0);
      p5[2] = -energy/3.0 * std::sin(pi/3.0);
      p5[3] = 0;

      bool phase_space_ok = eeAAA::check_phase_space(p1, p2, p3, p4, p5, me, me, 0, 0, 0, 1e-10);
      if (phase_space_ok)
        std::cout << "\nPhase-space point " << "PASS";
      eeAAAAmp.set_momenta(p1, p2, p3, p4, p5);

      cdouble amp_x = eeAAAAmp.amp(1,1,2,2,2);
      cdouble amp_f = eeAAAAmp.amp_feynman(1,1,2,2,2);
      cdouble amp_fx = eeAAAAmp.amp_feynman_x(1,1,2,2,2);
      cdouble amp_p = eeAAAAmp.amp_permutation(1,1,2,2,2);


      std::cout << "\n " << " Feynman_s   = " << amp_f << "\n";
      std::cout << "  Feynman_ns  = " << amp_fx << "\n";
      std::cout << "  x-factor    = " << amp_x << "\n";
      std::cout << "  Permutation = " << amp_p << "\n";
      print_ratio("Feynman_ns / Feynman_s",amp_fx,amp_f);
      print_ratio("x-factor / Permutation",amp_x,amp_p);
      print_ratio("Feynman_s / x-factor",amp_f,amp_x);
      print_ratio("Feynman_ns / x-factor",amp_fx,amp_x);
      print_ratio("Feynman_s / Permutation",amp_f,amp_p);
      print_reference_comparison(eeAAAAmp);

      energy = 600.0;

      p1[0] = energy/2.0;
      p1[1] = 0.0;
      p1[2] = 0.0;
      p1[3] = std::sqrt(energy*energy/4.0 - me*me);

      p2[0] = energy/2.0;
      p2[1] = 0.0;
      p2[2] = 0.0;
      p2[3] = -p1[3];

      p3[0] = energy/3.0;
      p3[2] = -energy/3.0;
      p3[1] = 0;
      p3[3] = 0;

      p4[0] = energy/3.0;
      p4[2] = energy/3.0 * std::cos(pi/3.0);
      p4[1] = energy/3.0 * std::sin(pi/3.0);
      p4[3] = 0;

      p5[0] = energy/3.0;
      p5[2] = energy/3.0 * std::cos(pi/3.0);
      p5[1] = -energy/3.0 * std::sin(pi/3.0);
      p5[3] = 0;

      phase_space_ok = eeAAA::check_phase_space(p1, p2, p3, p4, p5, me, me, 0, 0, 0, 1e-10);
      if (phase_space_ok)
        std::cout << "\nPhase-space point " << "PASS";
      eeAAAAmp.set_momenta(p1, p2, p3, p4, p5);

      amp_x = eeAAAAmp.amp(1,1,2,2,2);
      amp_f = eeAAAAmp.amp_feynman(1,1,2,2,2);
      amp_fx = eeAAAAmp.amp_feynman_x(1,1,2,2,2);
      amp_p = eeAAAAmp.amp_permutation(1,1,2,2,2);


      std::cout << "\n " << " Feynman_s   = " << amp_f << "\n";
      std::cout << "  Feynman_ns  = " << amp_fx << "\n";
      std::cout << "  x-factor    = " << amp_x << "\n";
      std::cout << "  Permutation = " << amp_p << "\n";
      print_ratio("Feynman_ns / Feynman_s",amp_fx,amp_f);
      print_ratio("x-factor / Permutation",amp_x,amp_p);
      print_ratio("Feynman_s / x-factor",amp_f,amp_x);
      print_ratio("Feynman_ns / x-factor",amp_fx,amp_x);
      print_ratio("Feynman_s / Permutation",amp_f,amp_p);
      print_reference_comparison(eeAAAAmp);
      if (2 == 2) {
        // Fixed seed makes the test reproducible.
        std::mt19937 rng(12345);

        const int Npoints = 10;

        int point = 0;
        int attempts = 0;

        while(point < Npoints && attempts < 100*Npoints){

          attempts++;

          bool generated = eeAAA::make_random_phase_space(energy, me, p1, p2, p3, p4, p5, rng);

          if(!generated)
            continue;

          // Check phase space
          bool phase_space_ok = eeAAA::check_phase_space(p1, p2, p3, p4, p5, me, me, 0, 0, 0, 1e-10);

          if(!phase_space_ok)
            continue;

          point++;

          std::cout << "\nPhase-space point " << point << "\n";

          std::cout << "  Phase space: PASS\n";

          eeAAAAmp.set_momenta(p1, p2, p3, p4, p5);

          cdouble amp_x = eeAAAAmp.amp(1,1,2,2,2);
          cdouble amp_f = eeAAAAmp.amp_feynman(1,1,2,2,2);
          cdouble amp_fx = eeAAAAmp.amp_feynman_x(1,1,2,2,2);
          cdouble amp_p = eeAAAAmp.amp_permutation(1,1,2,2,2);


          std::cout << "\n " << " Feynman_s   = " << amp_f << "\n";
          std::cout << "  Feynman_ns  = " << amp_fx << "\n";
          std::cout << "  x-factor    = " << amp_x << "\n";
          std::cout << "  Permutation = " << amp_p << "\n";
      print_ratio("Feynman_ns / Feynman_s",amp_fx,amp_f);
      print_ratio("x-factor / Permutation",amp_x,amp_p);
      print_ratio("Feynman_s / x-factor",amp_f,amp_x);
      print_ratio("Feynman_ns / x-factor",amp_fx,amp_x);
      print_ratio("Feynman_s / Permutation",amp_f,amp_p);
          print_reference_comparison(eeAAAAmp);
        }

        std::cout << "\nGenerated " << point << " valid phase-space points after " << attempts << " attempts.\n";
        if(point != Npoints) n++;
      }
    }

    {
      // Exercise the complex identity beyond one spin, photon order, or CM frame.
      std::mt19937 rng(20261015);
      const ldouble mass_energy[][2]={{0.0005L,600.0L},{1.0L,6.0L},{1.0L,2.1L}};
      const int Npoints=40;
      int checks=0, permutation_checks=0, permutation_sum_checks=0, failures=0, reported_failures=0;
      ldouble max_relative_error=0, max_permutation_error=0, max_permutation_sum_error=0;
      for(const auto& configuration : mass_energy){
        const ldouble mass=configuration[0], energy=configuration[1];
        const ldouble phase_tolerance=128*std::numeric_limits<ldouble>::epsilon()*energy*energy;
        eeAAA amplitude(0.31333L,mass);
        int point=0, attempts=0;
        while(point<Npoints && attempts<100*Npoints){
          attempts++;
          ldouble p1[4],p2[4],p3[4],p4[4],p5[4];
          if(!eeAAA::make_random_phase_space(energy,mass,p1,p2,p3,p4,p5,rng)
              || !eeAAA::check_phase_space(p1,p2,p3,p4,p5,mass,mass,0,0,0,phase_tolerance))
            continue;
          point++;
          for(int frame=0;frame<3;frame++){
            ldouble* momenta[]={p1,p2,p3,p4,p5};
            if(frame==1){
              const ldouble axis[]={1.0L/3,2.0L/3,2.0L/3};
              for(auto momentum : momenta) rotate_momentum(momentum,axis,0.7L);
            }
            if(frame==2){
              const ldouble velocity[]={0.21L,-0.13L,0.17L};
              for(auto momentum : momenta) boost_momentum(momentum,velocity);
            }
            if(!eeAAA::check_phase_space(p1,p2,p3,p4,p5,mass,mass,0,0,0,phase_tolerance)){
              failures++;
              continue;
            }
            amplitude.set_momenta(p1,p2,p3,p4,p5);
            cdouble reference[2][2];
            for(int i=0;i<2;i++)
              for(int j=0;j<2;j++)
                reference[i][j]=amplitude.amp_nonspurious(2*i-1,2*j-1,2,2,2);
            const ldouble* photons[]={p3,p4,p5};
            int order[]={0,1,2};
            do{
              amplitude.set_momenta(p1,p2,photons[order[0]],photons[order[1]],photons[order[2]]);
              for(int i=0;i<2;i++){
                for(int j=0;j<2;j++){
                  const cdouble h=amplitude.amp_nonspurious(2*i-1,2*j-1,2,2,2);
                  const cdouble s2=amplitude.amp_spurious(2*i-1,2*j-1,2,2,2);
                  checks++;
                  max_relative_error=std::max(max_relative_error,relative_error(h,s2));
                  if(!agrees(h,s2)){
                    failures++;
                    if(reported_failures++<5)
                      std::cout << "  H/S2 FAIL: m=" << mass << ", energy=" << energy
                        << ", point=" << point << ", frame=" << frame
                        << ", order=" << order[0]+3 << order[1]+3 << order[2]+3
                        << ", spins=" << 2*i-1 << "," << 2*j-1
                        << ", H=" << h << ", S2=" << s2 << "\n";
                  }
                  const cdouble permutation_sum=amplitude.amp_permutation(2*i-1,2*j-1,2,2,2);
                  permutation_sum_checks++;
                  max_permutation_sum_error=std::max(max_permutation_sum_error,relative_error(h,permutation_sum));
                  if(!agrees(h,permutation_sum)){
                    failures++;
                    if(reported_failures++<5)
                      std::cout << "  H/permutation FAIL: m=" << mass << ", energy=" << energy
                        << ", point=" << point << ", frame=" << frame
                        << ", order=" << order[0]+3 << order[1]+3 << order[2]+3
                        << ", spins=" << 2*i-1 << "," << 2*j-1
                        << ", H=" << h << ", permutation=" << permutation_sum << "\n";
                  }
                  permutation_checks++;
                  max_permutation_error=std::max(max_permutation_error,relative_error(h,reference[i][j]));
                  if(!agrees(h,reference[i][j])) failures++;
                }
              }
            }while(std::next_permutation(order,order+3));
            for(int ds3=-2;ds3<=2;ds3+=4)
              for(int ds4=-2;ds4<=2;ds4+=4)
                for(int ds5=-2;ds5<=2;ds5+=4)
                  if(ds3<0 || ds4<0 || ds5<0)
                    if(amplitude.amp_nonspurious(1,1,ds3,ds4,ds5)!=cdouble(0,0)
                        || amplitude.amp_spurious(1,1,ds3,ds4,ds5)!=cdouble(0,0)
                        || amplitude.amp_permutation(1,1,ds3,ds4,ds5)!=cdouble(0,0))
                      failures++;
          }
        }
        if(point!=Npoints) failures++;
      }
      n+=failures;
      std::cout << "\n  H = S2 = permutation regression: " << (failures==0 ? "PASS" : "FAIL")
        << " (" << checks << " complex amplitude pairs; all fermion spins,\n"
        << "    all 6 photon orders, CM/rotated/boosted frames, 3 mass/energy choices)\n"
        << "    Maximum relative H/S2 difference: " << max_relative_error << "\n"
        << "    Photon permutation checks: " << permutation_checks
        << "; maximum relative difference: " << max_permutation_error << "\n"
        << "    Permutation-sum/H checks: " << permutation_sum_checks
        << "; maximum relative difference: " << max_permutation_sum_error << "\n"
        << "    Failures: " << failures << "\n";
    }
    {
      // Numerical agreement of the raw eight-term sums with the constructive H.
      // Relative tolerance for nonzero amplitudes; absolute tolerance for zero spin components.
      std::mt19937 rng(20261031);
      int comparisons=0, failures=0;
      ldouble max_relative_error=0;
      for(ldouble energy : {2.1L,6.L}){
        eeAAA amplitude(0.31333L,1.L);
        for(int point=0;point<12;point++){
          ldouble p1[4],p2[4],p3[4],p4[4],p5[4];
          if(!eeAAA::make_random_phase_space(energy,1.L,p1,p2,p3,p4,p5,rng)){
            failures++;
            continue;
          }
          for(int frame=0;frame<3;frame++){
            ldouble* momenta[]={p1,p2,p3,p4,p5};
            if(frame==1){
              const ldouble axis[]={1.L/3,2.L/3,2.L/3};
              for(auto p:momenta) rotate_momentum(p,axis,0.7L);
            }
            if(frame==2){
              const ldouble velocity[]={0.21L,-0.13L,0.17L};
              for(auto p:momenta) boost_momentum(p,velocity);
            }
            const ldouble* photons[]={p3,p4,p5};
            int order[]={0,1,2};
            do{
              amplitude.set_momenta(p1,p2,photons[order[0]],photons[order[1]],photons[order[2]]);
              for(int s1:{-1,1}) for(int s2:{-1,1}){
                const cdouble h=amplitude.amp_nonspurious(s1,s2,2,2,2);
                const cdouble values[]={amplitude.amp_feynman(s1,s2,2,2,2),
                                        amplitude.amp_feynman_x(s1,s2,2,2,2)};
                for(const auto& value:values){
                  comparisons++;
                  const ldouble scale=0.31333L*0.31333L*0.31333L/energy;
                  if(std::abs(h)>1e-8L*scale)
                    max_relative_error=std::max(max_relative_error,relative_error(h,value));
                  if(!std::isfinite(std::abs(value))
                      || std::abs(value-h)>1e-9L*std::abs(h)+1e-12L*scale) failures++;
                }
              }
            }while(std::next_permutation(order,order+3));
          }
        }
      }
      eeAAA guard(0.31333L,1.L);
      for(int a:{-2,2}) for(int b:{-2,2}) for(int c:{-2,2})
        if(a<0 || b<0 || c<0)
          if(guard.amp_feynman(1,1,a,b,c)!=cdouble(0,0)
              || guard.amp_feynman_r(1,1,a,b,c)!=cdouble(0,0)
              || guard.amp_feynman_x(1,1,a,b,c)!=cdouble(0,0)) failures++;
      n+=failures;
      std::cout << "\n  Eight-term Feynman/H comparison: " << (failures==0?"PASS":"FAIL")
        << " (" << comparisons << " comparisons, fixed cyclic and momentum-inserted references)\n"
        << "    Maximum relative difference (nonzero components): " << max_relative_error
        << "; failures: " << failures << "\n";
    }
    return n;
  }


}
