
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

//File:  SPINAS/SM/nenenmnm.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/nenenmnm.h"

namespace spinas {
  //Constructors
  nenenmnm::nenenmnm(const ldouble& echarge, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propZ(MZ,WZ),
    p1(particle(0)), p2(particle(0)),
    p3(particle(0)), p4(particle(0)),
    a14a(sproduct(ANGLE,&p1,&p4)),
    s23s(sproduct(SQUARE,&p2,&p3))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    preZ = e*e/(4.0*CW*CW*SW*SW);
  }
  void nenenmnm::set_masses(const ldouble& massW){
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
  }
  void nenenmnm::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a14a.update();
    s23s.update();
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]+mom2[j];
    pDenSZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble nenenmnm::amp(){
    constexpr ldouble two = 2;
    
    //Z Boson
    //Defined above:
    //gL=1.0;
    //gR=0;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //all in:
    //- EE^2 [23] <14> /(4 CW^2 SW^2 (s-MZ^2))
    //34 out:
    //+ EE^2 [23] <14>/(4 CW^2 SW^2 (s-MZ^2))
    //=+preZ [23] <14> /(s-MZ^2)
    return + two*preZ*s23s.v()*a14a.v()/pDenSZ;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble nenenmnm::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    M = amp();
    amp2 = std::pow(std::abs(M),2);

    return amp2;
  }

  



  //  Tests
  int test_nenenmnm(){
    int n=0;//Number of fails
    std::cout<<"\t* ne, Ne -> nm, Nm      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### MW=80.385, pspatial=250\n";
      ldouble MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      nenenmnm nenenmnmAmp = nenenmnm(0.31333,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {8.078433614990031E-02,7.271121379961443E-02,6.506299262565936E-02,5.783967262803514E-02,5.104125380674175E-02,4.466773616177921E-02,3.871911969314749E-02,3.319540440084662E-02,2.809659028487658E-02,2.342267734523737E-02,1.917366558192901E-02,1.534955499495148E-02,1.195034558430478E-02,8.976037349988926E-03,6.426630292003905E-03,4.302124410349723E-03,2.602519705026376E-03,1.327816176033866E-03,4.780138233721919E-04,5.311264704135482E-05};
      i += nenenmnmAmp.test_2to2_amp2([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += nenenmnmAmp.test_2to2_amp2_rotations([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += nenenmnmAmp.test_2to2_amp2_boosts([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      i += nenenmnmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH);
      //std::cout<<"########### MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {1.738856354824795E-20,1.565085042574059E-20,1.400459588862836E-20,1.244979993691125E-20,1.098646257058927E-20,9.614583789662410E-21,8.334163594130673E-21,7.145201983994062E-21,6.047698959252574E-21,5.041654519906210E-21,4.127068665954971E-21,3.303941397398855E-21,2.572272714237863E-21,1.932062616471995E-21,1.383311104101250E-21,9.260181771256306E-22,5.601838355451346E-22,2.858080793597627E-22,1.028909085695147E-22,1.143232317439056E-23};
      i += nenenmnmAmp.test_2to2_amp2([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += nenenmnmAmp.test_2to2_amp2_rotations([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += nenenmnmAmp.test_2to2_amp2_boosts([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      i += nenenmnmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH2);
      //std::cout<<"########### MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      nenenmnmAmp.set_masses(MW);
      ldouble dataCH3[20] = {1.474532126945223E-05,1.327175859163714E-05,1.187575184423338E-05,1.055730102724094E-05,9.316406140659822E-06,8.153067184490022E-06,7.067284158731541E-06,6.059057063384381E-06,5.128385898448540E-06,4.275270663924019E-06,3.499711359810818E-06,2.801707986108938E-06,2.181260542818377E-06,1.638369029939137E-06,1.173033447471216E-06,7.852537954146159E-07,4.750300737693354E-07,2.423622825353754E-07,8.725042171273526E-08,9.694491301415056E-09};
      i += nenenmnmAmp.test_2to2_amp2([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += nenenmnmAmp.test_2to2_amp2_rotations([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += nenenmnmAmp.test_2to2_amp2_boosts([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      i += nenenmnmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH3);
      //std::cout<<"########### MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      nenenmnmAmp.set_masses(MW);
      ldouble dataCH4[20] = {9.661663465289290E-02,8.696132336608177E-02,7.781418635752387E-02,6.917522362721919E-02,6.104443517516771E-02,5.342182100136945E-02,4.630738110582441E-02,3.970111548853259E-02,3.360302414949399E-02,2.801310708870860E-02,2.293136430617643E-02,1.835779580189747E-02,1.429240157587174E-02,1.073518162809921E-02,7.686135958579909E-03,5.145264567313825E-03,3.112567454300955E-03,1.588044619541305E-03,5.716960630348705E-04,6.352178478165247E-05};
      i += nenenmnmAmp.test_2to2_amp2([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += nenenmnmAmp.test_2to2_amp2_rotations([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += nenenmnmAmp.test_2to2_amp2_boosts([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      i += nenenmnmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return nenenmnmAmp.amp2(); }, 0,0,0,0,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
