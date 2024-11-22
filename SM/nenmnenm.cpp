
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

//File:  SPINAS/SM/nenmnenm.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/nenmnenm.h"

namespace spinas {
  //Constructors
  nenmnenm::nenmnenm(const ldouble& echarge, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    p1(particle(0)), p2(particle(0)),
    p3(particle(0)), p4(particle(0)),
    a12a(sproduct(ANGLE,&p1,&p2)),
    s34s(sproduct(SQUARE,&p3,&p4))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    preZ = e*e/(4.0*CW*CW*SW*SW);
  }
  void nenmnenm::set_masses(const ldouble& massW){
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
  }
  void nenmnenm::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a12a.update();
    s34s.update();
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]-mom3[j];
    pDenTZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble nenmnenm::amp(){
    constexpr ldouble two = 2;
    
    //Z Boson
    //Defined above:
    //gL=1.0;
    //gR=0;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //neNeNmnm all in:
    //+ preZ [23] <14> /(s-MZ^2)
    //nenmNeNm: 4->2->3->4
    //+ preZ [34] <12> /(t-MZ^2)
    //34 out:
    //+ preZ [34] <12> /(t-MZ^2)
    return + two*preZ*s34s.v()*a12a.v()/pDenTZ;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble nenmnenm::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    M = amp();
    amp2 = std::pow(std::abs(M),2);

    return amp2;
  }

  



  //  Tests
  int test_nenmnenm(){
    int n=0;//Number of fails
    std::cout<<"\t* ne, nm -> ne, nm      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### MW=80.385, pspatial=250\n";
      ldouble MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      nenmnenm nenmnenmAmp = nenmnenm(0.31333,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.333348441203347E+01,6.765724637646326E+00,3.167401797998023E+00,1.829511734007181E+00,1.189857096389561E+00,8.352505216849498E-01,6.184217313224624E-01,4.762451742915185E-01,3.780055923582760E-01,3.073014584801713E-01,2.547271581483846E-01,2.145752401186570E-01,1.832194372137332E-01,1.582662161998009E-01,1.380840716310245E-01,1.215299059602808E-01,1.077835321923638E-01,9.624409005627466E-02,8.646325281239219E-02,7.810101460838127E-02};
      i += nenmnenmAmp.test_2to2_amp2([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += nenmnenmAmp.test_2to2_amp2_rotations([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += nenmnenmAmp.test_2to2_amp2_boosts([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += nenmnenmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      //std::cout<<"########### MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {1.829171706102774E-20,1.829171706014983E-20,1.829171705927193E-20,1.829171705839402E-20,1.829171705751612E-20,1.829171705663821E-20,1.829171705576031E-20,1.829171705488240E-20,1.829171705400450E-20,1.829171705312659E-20,1.829171705224868E-20,1.829171705137078E-20,1.829171705049287E-20,1.829171704961497E-20,1.829171704873706E-20,1.829171704785916E-20,1.829171704698125E-20,1.829171704610335E-20,1.829171704522544E-20,1.829171704434754E-20};
      i += nenmnenmAmp.test_2to2_amp2([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += nenmnenmAmp.test_2to2_amp2_rotations([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += nenmnenmAmp.test_2to2_amp2_boosts([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += nenmnenmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      //std::cout<<"########### MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      nenmnenmAmp.set_masses(MW);
      ldouble dataCH3[20] = {1.507614303092817E-05,1.505539137780902E-05,1.503468254083861E-05,1.501401640230984E-05,1.499339284491983E-05,1.497281175176825E-05,1.495227300635565E-05,1.493177649258184E-05,1.491132209474421E-05,1.489090969753614E-05,1.487053918604533E-05,1.485021044575223E-05,1.482992336252840E-05,1.480967782263491E-05,1.478947371272074E-05,1.476931091982123E-05,1.474918933135647E-05,1.472910883512971E-05,1.470906931932584E-05,1.468907067250981E-05};
      i += nenmnenmAmp.test_2to2_amp2([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += nenmnenmAmp.test_2to2_amp2_rotations([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += nenmnenmAmp.test_2to2_amp2_boosts([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += nenmnenmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      //std::cout<<"########### MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      nenmnenmAmp.set_masses(MW);
      ldouble dataCH4[20] = {3.989643605062781E+00,2.174880248567598E+00,1.366293251210003E+00,9.372206902943574E-01,6.825811827506973E-01,5.192012997320203E-01,4.081618210057981E-01,3.292780014961992E-01,2.712341263686462E-01,2.272871237440712E-01,1.932159449402035E-01,1.662688220254061E-01,1.445896898207890E-01,1.268898077200900E-01,1.122517110825307E-01,1.000078794456197E-01,8.966326504434634E-02,8.084445472334295E-02,7.326549275664614E-02,6.670440441781883E-02};
      i += nenmnenmAmp.test_2to2_amp2([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += nenmnenmAmp.test_2to2_amp2_rotations([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += nenmnenmAmp.test_2to2_amp2_boosts([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += nenmnenmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nenmnenmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
