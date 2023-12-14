
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

//File:  SPINAS/SM/uuAA.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uuAA.h"

namespace spinas {

  uuAA::uuAA(const ldouble& echarge, const ldouble& massu):
    e(echarge), Qu(2.0/3.0), mu(massu), prop(massu,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(mu);
    p2=particle(mu);
    p3=particle(0);
    p4=particle(0);
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
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
  }
  void uuAA::set_masses(const ldouble& massu){
    mu=massu;
    p1.set_mass(mu);
    p2.set_mass(mu);
    prop.set_mass(mu);
  }
  void uuAA::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s34s.update();
    a34a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    s24s.update();
    a24a.update();
    s23s.update();
    a23a.update();
    s14s.update();
    a14a.update();
    s314a.update();
    s413a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT = prop.den(propTP);
    pDenU = prop.den(propUP);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uuAA::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    //No sign changes due to p3 and p4 being outgoing.
    if(ds3>0&&ds4>0){
      //mu[34]^2<12>
      return 2.0*e*e*Qu*Qu*mu*s34s.v()*s34s.v()*a12a.v(ds1,ds2)/pDenT/pDenU;
    }
    else if(ds3<0&&ds4<0){
      //mu<34>^2[12]
      return 2.0*e*e*Qu*Qu*mu*a34a.v()*a34a.v()*s12s.v(ds1,ds2)/pDenT/pDenU;
    }
    else if(ds3>0&&ds4<0){
      //([13]<24>+[23]<14>)[314>
      return 2.0*e*e*Qu*Qu*(s13s.v(ds1)*a24a.v(ds2)+s23s.v(ds2)*a14a.v(ds1))*s314a.v()/pDenT/pDenU;
    }
    else if(ds3<0&&ds4>0){
      //(<13>[24]+<23>[14])*[413>
      return 2.0*e*e*Qu*Qu*(a13a.v(ds1)*s24s.v(ds2)+a23a.v(ds2)*s14s.v(ds1))*s413a.v()/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }
  //set_momenta(...) must be called before amp2().
  ldouble uuAA::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-2;j4<=2;j4+=4){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/3^2=1/9
    //Symmetry factor for identical photons 1/2
    return amp2/72.0;
  }


  



  //  Tests
  int test_uuAA(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> A , A       :";
    {//amp^2
      int i=0;
      // mu=0.0042, pspatial=250
      ldouble mu=0.0042;
      ldouble EE=0.31333;
      uuAA uuAAAmp = uuAA(EE,mu);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.476690100267933E-02,7.878578424827100E-03,4.533082135212090E-03,3.126453012698946E-03,2.370200796112170E-03,1.913841699834484E-03,1.623644120638945E-03,1.438498065102407E-03,1.327694543371012E-03,1.275625219248172E-03,1.275625219248172E-03,1.327694543371011E-03,1.438498065102407E-03,1.623644120638945E-03,1.913841699834484E-03,2.370200796112169E-03,3.126453012698945E-03,4.533082135212088E-03,7.878578424827096E-03,2.476690100267927E-02};
      i += uuAAAmp.test_2to2_amp2([&]() { return uuAAAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH);
      i += uuAAAmp.test_2to2_amp2_rotations([&]() { return uuAAAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH);
      i += uuAAAmp.test_2to2_amp2_boosts([&]() { return uuAAAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH);
      i += uuAAAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuAAAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH);
      //Close to threshold
      pspatial = 0.001;
      ldouble dataCH2[20] = {1.412302631327627E-03,1.410332670637965E-03,1.408200785119260E-03,1.406057476526193E-03,1.404023853515509E-03,1.402196361860851E-03,1.400650465121158E-03,1.399443468207727E-03,1.398616626873141E-03,1.398196647154393E-03,1.398196647154393E-03,1.398616626873141E-03,1.399443468207727E-03,1.400650465121159E-03,1.402196361860851E-03,1.404023853515509E-03,1.406057476526193E-03,1.408200785119260E-03,1.410332670637965E-03,1.412302631327628E-03};
      i += uuAAAmp.test_2to2_amp2([&]() { return uuAAAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH2);
      i += uuAAAmp.test_2to2_amp2_rotations([&]() { return uuAAAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH2);
      i += uuAAAmp.test_2to2_amp2_boosts([&]() { return uuAAAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH2);
      i += uuAAAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuAAAmp.amp2(); }, mu,mu,0,0,pspatial,dataCH2);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
