
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

//File:  SPINAS/SM/uuAh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uuAh.h"

namespace spinas {

  uuAh::uuAh(const ldouble& echarge, const ldouble& massu, const ldouble& massh, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qu(2.0/3.0), mu(massu), mh(massh), prop(massu,0), MW(massW), SW(sinW) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(mu);
    p2=particle(mu);
    p3=particle(0);
    p4=particle(mh);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s342a = sproduct(SQUARE,&p3,&p4,&p2);
    a342s = sproduct(ANGLE,&p3,&p4,&p2);
    s341a = sproduct(SQUARE,&p3,&p4,&p1);
    a341s = sproduct(ANGLE,&p3,&p4,&p1);
    s3243s = sproduct(SQUARE,&p3,&p2,&p4,&p3);
    a3243a = sproduct(ANGLE,&p3,&p2,&p4,&p3);
    //2.*MW*SW/EE
    pre = e*e*Qu*mu/(sqrt2*MW*SW);
  }
  void uuAh::set_masses(const ldouble& massu, const ldouble& massh, const ldouble& massW){
    mu=massu;
    mh=massh;
    MW=massW;
    p1.set_mass(mu);
    p2.set_mass(mu);
    p4.set_mass(mh);
    prop.set_mass(mu);
    pre = sqrt2*e*e*Qu*mu/(2.0*MW*SW);
  }
  void uuAh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
  cdouble uuAh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    /*- mh^2[13][23] - mu[13][342> - mu[23][341> - <12>[3243]
      - mh^2<13><23> - mu<13><342] - mu<23><341] - [12]<3243>
      Becomes after a sign change due to p3 and p4 being outgoing:*/
    //pre = e*e*Qu*mu/(sqrt2*MW*SW);
    if(ds3>0){
      //- mh^2[13][23] + mu[13][342> + mu[23][341> + <12>[3243]
      return pre*(- mh*mh*s13s.v(ds1)*s23s.v(ds2) + mu*s13s.v(ds1)*s342a.v(ds2) + mu*s23s.v(ds2)*s341a.v(ds1) + a12a.v(ds1,ds2)*s3243s.v())/pDenT/pDenU;
    }
    else if(ds3<0){
      //- mh^2<13><23> + mu<13><342] + mu<23><341] + [12]<3243>
      return pre*(- mh*mh*a13a.v(ds1)*a23a.v(ds2) + mu*a13a.v(ds1)*a342s.v(ds2) + mu*a23a.v(ds2)*a341s.v(ds1) + s12s.v(ds1,ds2)*a3243a.v())/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uuAh::amp2(){
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
  int test_uuAh(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> A , h       :";
    {//amp^2
      int i=0;
      // mu=0.0042, mmu=0.105, pspatial=250
      ldouble mu=0.0042, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      uuAh uuAhAmp = uuAh(EE,mu,mh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.032553544288259E-10,7.141404359875552E-11,4.529690767301439E-11,3.431583915225801E-11,2.841203887068495E-11,2.484940076919558E-11,2.258392833521652E-11,2.113855692225198E-11,2.027355203568877E-11,1.986706477695718E-11,1.986706477695718E-11,2.027355203568876E-11,2.113855692225198E-11,2.258392833521652E-11,2.484940076919558E-11,2.841203887068495E-11,3.431583915225798E-11,4.529690767301434E-11,7.141404359875545E-11,2.032553544288248E-10};
      i += uuAhAmp.test_2to2_amp2([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH);
      i += uuAhAmp.test_2to2_amp2_rotations([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH);
      i += uuAhAmp.test_2to2_amp2_boosts([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH);
      i += uuAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH);
      //Close to Higgs threshold
      pspatial = 125.1;
      ldouble dataCH2[20] = {3.357017709702232E-10,1.179492721815369E-10,7.481353851850339E-11,5.667692317211186E-11,4.692605469254505E-11,4.104190991573570E-11,3.730019733789016E-11,3.491298471248409E-11,3.348432037884490E-11,3.281295556033407E-11,3.281295556033407E-11,3.348432037884490E-11,3.491298471248408E-11,3.730019733789016E-11,4.104190991573570E-11,4.692605469254505E-11,5.667692317211186E-11,7.481353851850337E-11,1.179492721815368E-10,3.357017709702223E-10};
      i += uuAhAmp.test_2to2_amp2([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH2);
      i += uuAhAmp.test_2to2_amp2_rotations([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH2);
      i += uuAhAmp.test_2to2_amp2_boosts([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH2);
      i += uuAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH2);
      //mu~mh and close to threshold
      mu = 125.1;
      mh = 125;
      pspatial = 95;
      uuAhAmp.set_masses(mu,mh,MW);
      ldouble dataCH4[20] = {2.133381098972853E-02,1.704337512549444E-02,1.413377321282512E-02,1.209471889296658E-02,1.063691704496320E-02,9.587538190851522E-03,8.839528264732584E-03,8.325374441641165E-03,8.002891789349903E-03,7.847369045934941E-03,7.847369045934941E-03,8.002891789349900E-03,8.325374441641161E-03,8.839528264732584E-03,9.587538190851520E-03,1.063691704496320E-02,1.209471889296658E-02,1.413377321282511E-02,1.704337512549444E-02,2.133381098972853E-02};
      i += uuAhAmp.test_2to2_amp2([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH4);
      i += uuAhAmp.test_2to2_amp2_rotations([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH4);
      i += uuAhAmp.test_2to2_amp2_boosts([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH4);
      i += uuAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH4);
      //Invert mu and mh and close to threshold
      mu = 125;
      mh = 0.0005;
      pspatial = 125.1;
      uuAhAmp.set_masses(mu,mh,MW);
      ldouble dataCH3[20] = {2.553267029596537E-02,1.883942461839153E-02,1.487675820984225E-02,1.235230956540446E-02,1.066711049181586E-02,9.513006471875555E-03,8.719479790117245E-03,8.187757774527548E-03,7.859860247590429E-03,7.703252796180197E-03,7.703252796180197E-03,7.859860247590427E-03,8.187757774527545E-03,8.719479790117245E-03,9.513006471875553E-03,1.066711049181586E-02,1.235230956540446E-02,1.487675820984225E-02,1.883942461839152E-02,2.553267029596535E-02};
      i += uuAhAmp.test_2to2_amp2([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH3);
      i += uuAhAmp.test_2to2_amp2_rotations([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH3);
      i += uuAhAmp.test_2to2_amp2_boosts([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH3);
      i += uuAhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuAhAmp.amp2(); }, mu,mu,0,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
