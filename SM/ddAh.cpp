
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

//File:  SPINAS/SM/ddAh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ddAh.h"

namespace spinas {

  ddAh::ddAh(const ldouble& echarge, const ldouble& massd, const ldouble& massh, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qd(-1.0/3.0), md(massd), mh(massh), prop(massd,0), MW(massW), SW(sinW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(md);
    p2=particle(md);
    p3=particle(0);
    p4=particle(mh);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s342a = sproduct(SQUARE,&p3,&p4,&p2,2);
    a342s = sproduct(ANGLE,&p3,&p4,&p2,2);
    s341a = sproduct(SQUARE,&p3,&p4,&p1,2);
    a341s = sproduct(ANGLE,&p3,&p4,&p1,2);
    s3243s = sproduct(SQUARE,&p3,&p2,&p4,&p3,2);
    a3243a = sproduct(ANGLE,&p3,&p2,&p4,&p3,2);
    //2.*MW*SW/EE
    pre = sqrt2*e*e*Qd*md/(2.0*MW*SW);
  }
  void ddAh::set_masses(const ldouble& massd, const ldouble& massh, const ldouble& massW){
    md=massd;
    mh=massh;
    MW=massW;
    p1.set_mass(md);
    p2.set_mass(md);
    p4.set_mass(mh);
    prop.set_mass(md);
    pre = sqrt2*e*e*Qd*md/(2.0*MW*SW);
  }
  void ddAh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s13s.update();
    a13a.update();
    s23s.update();
    a23a.update();
    s12s.update();
    a12a.update();
    s342a.update();
    a342s.update();
    s341a.update();
    a341s.update();
    s3243s.update();
    a3243a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT=prop.denominator(propTP);
    pDenU=prop.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ddAh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    /*- mh^2[13][23] - md[13][342> - md[23][341> - <12>[3243]
      - mh^2<13><23> - md<13><342] - md<23><341] - [12]<3243>
      Becomes after a sign change due to p3 and p4 being outgoing:*/
    //pre = sqrt2*e*e*Qd*md/(2.0*MW*SW);
    if(ds3>0){
      //- mh^2[13][23] + md[13][342> + md[23][341> + <12>[3243]
      return pre*(- mh*mh*s13s.v(ds1)*s23s.v(ds2) + md*s13s.v(ds1)*s342a.v(ds2) + md*s23s.v(ds2)*s341a.v(ds1) + a12a.v(ds1,ds2)*s3243s.v())/pDenT/pDenU;
    }
    else if(ds3<0){
      //- mh^2<13><23> + md<13><342] + md<23><341] + [12]<3243>
      return pre*(- mh*mh*a13a.v(ds1)*a23a.v(ds2) + md*a13a.v(ds1)*a342s.v(ds2) + md*a23a.v(ds2)*a341s.v(ds1) + s12s.v(ds1,ds2)*a3243a.v())/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ddAh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=4){
	  M = amp(j1,j2,j3);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2^2=1/4
    //Average over colors 1/3^2=1/9
    return amp2/36.0;
  }
  



  //  Tests
  int test_ddAh(){
    int n=0;//Number of fails
    std::cout<<"\t* d , D  -> A , h       :";
    {//amp^2
      int i=0;
      // md=0.0042, pspatial=250
      ldouble md=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ddAh ddAhAmp = ddAh(EE,md,mh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.620339228621996E-10,5.693083814571288E-11,3.611041738630147E-11,2.735637682314550E-11,2.264990339803551E-11,1.980979012288921E-11,1.800376937302786E-11,1.685152813439796E-11,1.616195153626187E-11,1.583790238298139E-11,1.583790238298139E-11,1.616195153626187E-11,1.685152813439796E-11,1.800376937302786E-11,1.980979012288921E-11,2.264990339803551E-11,2.735637682314549E-11,3.611041738630146E-11,5.693083814571287E-11,1.620339228621993E-10};
      i += ddAhAmp.test_2to2_amp2([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH);
      i += ddAhAmp.test_2to2_amp2_rotations([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH);
      i += ddAhAmp.test_2to2_amp2_boosts([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH);
      i += ddAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH);
      //Close to Higgs threshold
      pspatial = 125.1;
      ldouble dataCH2[20] = {2.676193854471343E-10,9.402843602963509E-11,5.964089412549134E-11,4.518249568188764E-11,3.740916315033510E-11,3.271835902462888E-11,2.973548870416051E-11,2.783241743781897E-11,2.669349499127624E-11,2.615828707448638E-11,2.615828707448638E-11,2.669349499127623E-11,2.783241743781897E-11,2.973548870416051E-11,3.271835902462888E-11,3.740916315033511E-11,4.518249568188762E-11,5.964089412549132E-11,9.402843602963506E-11,2.676193854471340E-10};
      i += ddAhAmp.test_2to2_amp2([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH2);
      i += ddAhAmp.test_2to2_amp2_rotations([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH2);
      i += ddAhAmp.test_2to2_amp2_boosts([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH2);
      i += ddAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH2);
      //md~mh and close to threshold
      md = 125.1;
      mh = 125;
      pspatial = 95;
      ddAhAmp.set_masses(md,mh,MW);
      ldouble dataCH4[20] = {5.333452747432132E-03,4.260843781373609E-03,3.533443303206279E-03,3.023679723241644E-03,2.659229261240801E-03,2.396884547712880E-03,2.209882066183146E-03,2.081343610410291E-03,2.000722947337476E-03,1.961842261483735E-03,1.961842261483735E-03,2.000722947337475E-03,2.081343610410290E-03,2.209882066183146E-03,2.396884547712880E-03,2.659229261240801E-03,3.023679723241644E-03,3.533443303206278E-03,4.260843781373609E-03,5.333452747432132E-03};
      i += ddAhAmp.test_2to2_amp2([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH4);
      i += ddAhAmp.test_2to2_amp2_rotations([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH4);
      i += ddAhAmp.test_2to2_amp2_boosts([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH4);
      i += ddAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH4);
      //Invert md and mh and close to threshold
      md = 125;
      mh = 0.0005;
      pspatial = 125.1;
      ddAhAmp.set_masses(md,mh,MW);
      ldouble dataCH3[20] = {6.383167573991343E-03,4.709856154597882E-03,3.719189552460563E-03,3.088077391351115E-03,2.666777622953964E-03,2.378251617968889E-03,2.179869947529311E-03,2.046939443631887E-03,1.964965061897607E-03,1.925813199045049E-03,1.925813199045049E-03,1.964965061897607E-03,2.046939443631886E-03,2.179869947529311E-03,2.378251617968888E-03,2.666777622953964E-03,3.088077391351115E-03,3.719189552460562E-03,4.709856154597881E-03,6.383167573991337E-03};
      i += ddAhAmp.test_2to2_amp2([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH3);
      i += ddAhAmp.test_2to2_amp2_rotations([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH3);
      i += ddAhAmp.test_2to2_amp2_boosts([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH3);
      i += ddAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddAhAmp.amp2(); }, md,md,0,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
