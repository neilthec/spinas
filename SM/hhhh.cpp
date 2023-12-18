
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

//File:  SPINAS/SM/hhhh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/hhhh.h"

namespace spinas {
  //Constructors
  hhhh::hhhh(const ldouble& echarge, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW):
    e(echarge), mh(massh), wh(widthh), MW(massW), SW(sinW),
    proph(mh,wh),
    p1(particle(mh)), p2(particle(mh)),
    p3(particle(mh)), p4(particle(mh))
  {
    constexpr ldouble sqrt2 = std::sqrt(2);
    preh = 3*e*e*mh*mh/(4*MW*MW*SW*SW);
    prehSTU = 3*mh*mh*preh;
  }
  void hhhh::set_masses(const ldouble& massh, const ldouble& massW){
    constexpr ldouble sqrt2 = std::sqrt(2);
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    p1.set_mass(mh);
    p2.set_mass(mh);
    p3.set_mass(mh);
    p4.set_mass(mh);
    preh = 3*e*e*mh*mh/(4*MW*MW*SW*SW);
    prehSTU = 3*mh*mh*preh;
  }
  void hhhh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //Propagator Momentum
    ldouble propPS[4], propPT[4], propPU[4];
    for(int j=0;j<4;j++){
      propPS[j] = mom1[j]+mom2[j];
      propPT[j] = mom1[j]-mom3[j];
      propPU[j] = mom1[j]-mom4[j];
    }
    pDenS = proph.denominator(propPS);
    pDenT = proph.denominator(propPT);
    pDenU = proph.denominator(propPU);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble hhhh::amp(){
    constexpr ldouble one = 1;
    cdouble amplitude(0,0);
    
    amplitude += prehSTU*(one/pDenS+one/pDenT+one/pDenU) + preh;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble hhhh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    M = amp();
    amp2 += std::pow(std::abs(M),2);

    //Symmetry factor 1/2
    return amp2/2.0;
  }

  



  //  Tests
  int test_hhhh(){
    int n=0;//Number of fails
    std::cout<<"\t* h , h  -> h , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125, wh=0, MW=80.385, SW=0.474;//Set width to 0 for comparison with Feynman diagrams.
      hhhh hhhhAmp = hhhh(0.31333,mh,wh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.266676389938904E-01,4.914424463497754E-02,5.566725917090273E-04,7.754876562737883E-03,2.494210057543500E-02,4.165927854682620E-02,5.533296766464278E-02,6.546231198419260E-02,7.209269658467868E-02,7.536282644302122E-02,7.536282644302118E-02,7.209269658467871E-02,6.546231198419264E-02,5.533296766464280E-02,4.165927854682616E-02,2.494210057543499E-02,7.754876562737869E-03,5.566725917089327E-04,4.914424463497737E-02,4.266676389938894E-01};
      i += hhhhAmp.test_2to2_amp2([&]() { return hhhhAmp.amp2(); }, mh,mh,mh,mh,pspatial,dataCH);
      i += hhhhAmp.test_2to2_amp2_rotations([&]() { return hhhhAmp.amp2(); }, mh,mh,mh,mh,pspatial,dataCH);
      i += hhhhAmp.test_2to2_amp2_boosts([&]() { return hhhhAmp.amp2(); }, mh,mh,mh,mh,pspatial,dataCH);
      i += hhhhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhhhAmp.amp2(); }, mh,mh,mh,mh,pspatial,dataCH);
      //std::cout<<"\n########### mh=125, MW=80.385, pspatial=0.008\n";
      pspatial = 0.008;
      ldouble dataCH2[20] = {5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00,5.023970030423047E+00};
      i += hhhhAmp.test_2to2_amp2([&]() { return hhhhAmp.amp2(); }, mh,mh,mh,mh,pspatial,dataCH2);
      i += hhhhAmp.test_2to2_amp2_rotations([&]() { return hhhhAmp.amp2(); }, mh,mh,mh,mh,pspatial,dataCH2);
      i += hhhhAmp.test_2to2_amp2_boosts([&]() { return hhhhAmp.amp2(); }, mh,mh,mh,mh,pspatial,dataCH2);
      i += hhhhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return hhhhAmp.amp2(); }, mh,mh,mh,mh,pspatial,dataCH2);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
