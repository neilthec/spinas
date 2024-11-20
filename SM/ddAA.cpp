
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

//File:  SPINAS/SM/ddAA.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ddAA.h"

namespace spinas {

  ddAA::ddAA(const ldouble& echarge, const ldouble& massd):
    e(echarge), Qd(-1.0/3.0), md(massd), prop(massd,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(md);
    p2=particle(md);
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
  void ddAA::set_masses(const ldouble& massd){
    md=massd;
    p1.set_mass(md);
    p2.set_mass(md);
    prop.set_mass(md);
  }
  void ddAA::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenT = prop.denominator(propTP);
    pDenU = prop.denominator(propUP);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ddAA::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    //No sign changes due to p3 and p4 being outgoing.
    if(ds3>0&&ds4>0){
      //md[34]^2<12>
      return 2.0*e*e*Qd*Qd*md*s34s.v()*s34s.v()*a12a.v(ds1,ds2)/pDenT/pDenU;
    }
    else if(ds3<0&&ds4<0){
      //md<34>^2[12]
      return 2.0*e*e*Qd*Qd*md*a34a.v()*a34a.v()*s12s.v(ds1,ds2)/pDenT/pDenU;
    }
    else if(ds3>0&&ds4<0){
      //([13]<24>+[23]<14>)[314>
      return 2.0*e*e*Qd*Qd*(s13s.v(ds1)*a24a.v(ds2)+s23s.v(ds2)*a14a.v(ds1))*s314a.v()/pDenT/pDenU;
    }
    else if(ds3<0&&ds4>0){
      //(<13>[24]+<23>[14])*[413>
      return 2.0*e*e*Qd*Qd*(a13a.v(ds1)*s24s.v(ds2)+a23a.v(ds2)*s14s.v(ds1))*s413a.v()/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }
  //set_momenta(...) must be called before amp2().
  ldouble ddAA::amp2(){
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
  int test_ddAA(){
    int n=0;//Number of fails
    std::cout<<"\t* d , D  -> A , A       :";
    {//amp^2
      int i=0;
      // md=0.0075, pspatial=250
      ldouble md=0.0075;
      ldouble EE=0.31333;
      ddAA ddAAAmp = ddAA(EE,md);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.547931304367646E-03,4.924111509853038E-04,2.833176333867474E-04,1.954033133392363E-04,1.481375498365880E-04,1.196151063313484E-04,1.014777576360364E-04,8.990612916647732E-05,8.298090905864879E-05,7.972657630102263E-05,7.972657630102263E-05,8.298090905864876E-05,8.990612916647729E-05,1.014777576360364E-04,1.196151063313484E-04,1.481375498365880E-04,1.954033133392362E-04,2.833176333867474E-04,4.924111509853036E-04,1.547931304367643E-03};
      i += ddAAAmp.test_2to2_amp2([&]() { return ddAAAmp.amp2(); }, md,md,0,0,pspatial,dataCH);
      i += ddAAAmp.test_2to2_amp2_rotations([&]() { return ddAAAmp.amp2(); }, md,md,0,0,pspatial,dataCH);
      i += ddAAAmp.test_2to2_amp2_boosts([&]() { return ddAAAmp.amp2(); }, md,md,0,0,pspatial,dataCH);
      i += ddAAAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddAAAmp.amp2(); }, md,md,0,0,pspatial,dataCH);
      //Close to threshold
      pspatial = 0.001;
      ldouble dataCH2[20] = {8.214416623040461E-05,8.213185167966661E-05,8.211830299994840E-05,8.210447301716053E-05,8.209117364563712E-05,8.207908494462669E-05,8.206876263910799E-05,8.206064422083176E-05,8.205505372256967E-05,8.205220523766732E-05,8.205220523766732E-05,8.205505372256966E-05,8.206064422083175E-05,8.206876263910801E-05,8.207908494462670E-05,8.209117364563713E-05,8.210447301716053E-05,8.211830299994840E-05,8.213185167966663E-05,8.214416623040461E-05};
      i += ddAAAmp.test_2to2_amp2([&]() { return ddAAAmp.amp2(); }, md,md,0,0,pspatial,dataCH2);
      i += ddAAAmp.test_2to2_amp2_rotations([&]() { return ddAAAmp.amp2(); }, md,md,0,0,pspatial,dataCH2);
      i += ddAAAmp.test_2to2_amp2_boosts([&]() { return ddAAAmp.amp2(); }, md,md,0,0,pspatial,dataCH2);
      i += ddAAAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddAAAmp.amp2(); }, md,md,0,0,pspatial,dataCH2);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
