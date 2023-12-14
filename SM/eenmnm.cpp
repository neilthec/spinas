
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

//File:  SPINAS/SM/eenmnm.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eenmnm.h"

namespace spinas {
  //Constructors
  eenmnm::eenmnm(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propZ(MZ,WZ),
    p1(particle(me)), p2(particle(me)),
    p3(particle(0)), p4(particle(0)),
    a13a(sproduct(ANGLE,&p1,&p3)),
    s13s(sproduct(SQUARE,&p1,&p3)),
    a14a(sproduct(ANGLE,&p1,&p4)),
    s14s(sproduct(SQUARE,&p1,&p4)),
    a23a(sproduct(ANGLE,&p2,&p3)),
    s23s(sproduct(SQUARE,&p2,&p3)),
    a24a(sproduct(ANGLE,&p2,&p4)),
    s24s(sproduct(SQUARE,&p2,&p4)),
    s12s(sproduct(SQUARE,&p1,&p2)),
    a12a(sproduct(ANGLE,&p1,&p2)),
    s34s(sproduct(SQUARE,&p3,&p4)),
    a34a(sproduct(ANGLE,&p3,&p4))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
  }
  void eenmnm::set_masses(const ldouble& masse, const ldouble& massW){
    me=masse;
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(me);
    p2.set_mass(me);
  }
  void eenmnm::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a13a.update();
    s13s.update();
    a14a.update();
    s14s.update();
    a23a.update();
    s23s.update();
    a24a.update();
    s24s.update();
    s12s.update();
    a12a.update();
    s34s.update();
    a34a.update();
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]+mom2[j];
    pDenSZ = propZ.den(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eenmnm::amp(const int& ds1, const int& ds2){
    constexpr ldouble two = 2;
    
    //Z Boson
    //Defined above:
    //gL=2.0*SW*SW-1.0;
    //gR=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //all in:
    //+(EE^2 (gL [23] <14> + gR [13] <24>)/(4 CW^2 SW^2 (s-MZ^2))
    //34 out:
    //-(EE^2 (gL [23] <14> + gR [13] <24>)/(4 CW^2 SW^2 (s-MZ^2))
    //=-preZ (gL [23] <14> + gR [13] <24>)/(s-MZ^2)
    return - two*preZ*( gL*s23s.v(ds2)*a14a.v(ds1) + gR*s13s.v(ds1)*a24a.v(ds2) )/pDenSZ;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eenmnm::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2){
	    M = amp(j1,j2);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    return amp2/4.0;
  }

  



  //  Tests
  int test_eenmnm(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> nm, Nm      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### me=0.0005, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      eenmnm eenmnmAmp = eenmnm(0.31333,me,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {6.126400759938745E-03,5.535880120378437E-03,4.999017110336632E-03,4.515811729813331E-03,4.086263978808533E-03,3.710373857322239E-03,3.388141365354448E-03,3.119566502905161E-03,2.904649269974378E-03,2.743389666562098E-03,2.635787692668322E-03,2.581843348293049E-03,2.581556633436280E-03,2.634927548098014E-03,2.741956092278253E-03,2.902642265976994E-03,3.116986069194239E-03,3.384987501929987E-03,3.706646564184240E-03,4.081963255956995E-03};
      i += eenmnmAmp.test_2to2_amp2([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eenmnmAmp.test_2to2_amp2_rotations([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eenmnmAmp.test_2to2_amp2_boosts([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eenmnmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      //std::cout<<"########### me=0.0005, MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {1.684489065507020E-21,1.522187298106318E-21,1.374322566791849E-21,1.240894871563611E-21,1.121904212421606E-21,1.017350589365833E-21,9.272340023962916E-22,8.515544515129827E-22,7.903119367159059E-22,7.435064580050613E-22,7.111380153804488E-22,6.932066088420685E-22,6.897122383899204E-22,7.006549040240044E-22,7.260346057443205E-22,7.658513435508689E-22,8.201051174436494E-22,8.887959274226620E-22,9.719237734879069E-22,1.069488655639384E-21};
      i += eenmnmAmp.test_2to2_amp2([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eenmnmAmp.test_2to2_amp2_rotations([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eenmnmAmp.test_2to2_amp2_boosts([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eenmnmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      //std::cout<<"########### me=0.105, MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      eenmnmAmp.set_masses(me,MW);
      ldouble dataCH3[20] = {1.438463782351513E-06,1.299866732956433E-06,1.173598142059417E-06,1.059658009660466E-06,9.580463357595792E-07,8.687631203567569E-07,7.918083634519993E-07,7.271820650453062E-07,6.748842251366776E-07,6.349148437261135E-07,6.072739208136139E-07,5.919614563991788E-07,5.889774504828084E-07,5.983219030645024E-07,6.199948141442610E-07,6.539961837220840E-07,7.003260117979715E-07,7.589842983719237E-07,8.299710434439402E-07,9.132862470140215E-07};
      i += eenmnmAmp.test_2to2_amp2([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH3);
      i += eenmnmAmp.test_2to2_amp2_rotations([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH3);
      i += eenmnmAmp.test_2to2_amp2_boosts([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH3);
      i += eenmnmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH3);
      //std::cout<<"########### me=0.105, MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      eenmnmAmp.set_masses(me,MW);
      ldouble dataCH4[20] = {5.687447183702413E-03,5.139457440809541E-03,4.640212377778950E-03,4.189711994610638E-03,3.787956291304607E-03,3.434945267860856E-03,3.130678924279385E-03,2.875157260560194E-03,2.668380276703284E-03,2.510347972708653E-03,2.401060348576303E-03,2.340517404306232E-03,2.328719139898443E-03,2.365665555352933E-03,2.451356650669703E-03,2.585792425848754E-03,2.768972880890084E-03,3.000898015793695E-03,3.281567830559586E-03,3.610982325187757E-03};
      i += eenmnmAmp.test_2to2_amp2([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH4);
      i += eenmnmAmp.test_2to2_amp2_rotations([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH4);
      i += eenmnmAmp.test_2to2_amp2_boosts([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH4);
      i += eenmnmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eenmnmAmp.amp2(); }, me,me,0,0,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
