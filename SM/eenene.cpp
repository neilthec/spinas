
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

//File:  SPINAS/SM/eenene.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eenene.h"

namespace spinas {
  //Constructors
  eenene::eenene(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& widthW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), MW(massW), WW(widthW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    propW(MW,WW),
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
  void eenene::set_masses(const ldouble& masse, const ldouble& massW){
    me=masse;
    MW=massW;
    propW.set_mass(MW);
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(me);
    p2.set_mass(me);
  }
  void eenene::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    ldouble propPS[4], propPT[4];
    for(int j=0;j<4;j++){
      propPS[j] = mom1[j]+mom2[j];
      propPT[j] = mom1[j]-mom3[j];
    }
    pDenSZ = propZ.denominator(propPS);
    pDenTW = propW.denominator(propPT);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eenene::amp(const int& ds1, const int& ds2){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Z Boson
    //Defined above:
    //gL=2.0*SW*SW-1.0;
    //gR=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //all in:
    //-2(EE^2 (gL [23] <14> + gR [13] <24>)/(4 CW^2 SW^2 (s-MZ^2))
    //34 out:
    //+2(EE^2 (gL [23] <14> + gR [13] <24>)/(4 CW^2 SW^2 (s-MZ^2))
    //=2preZ (gL [23] <14> + gR [13] <24>)/(s-MZ^2)
    amplitude = two*preZ*( gL*s23s.v(ds2)*a14a.v(ds1) + gR*s13s.v(ds1)*a24a.v(ds2) )/pDenSZ;

    //W Boson
    //All ingoing: -e^2  ( 2 MW^2 [23] <14> + Me^2 [13] <24>)/(4 MW^2 SW^2 (t-MW^2))
    //34 outgoing:  e^2  ( 2 MW^2 [23] <14> + Me^2 [13] <24>)/(4 MW^2 SW^2 (t-MW^2))
    amplitude +=  e*e/SW/SW*(   2.0*MW*MW*a14a.v(ds1)*s23s.v(ds2)
				+ me*me*s13s.v(ds1)*a24a.v(ds2)
				)/(2.0*MW*MW*pDenTW);

    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eenene::amp2(){
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
  int test_eenene(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> ne , Ne     :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### me=0.0005, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, MW=80.385, WW=0, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      eenene eeneneAmp = eenene(0.31333,me,MW,WW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.821326959543661E+01,4.319078646843732E+00,1.789120643723311E+00,9.287643824528204E-01,5.436984967075855E-01,3.422015887833890E-01,2.256208863673576E-01,1.533443036625346E-01,1.062685251435159E-01,7.449338182225153E-02,5.249521558422088E-02,3.700761158381562E-02,2.600591261130697E-02,1.818542475222659E-02,1.267751013203036E-02,8.887552565195949E-03,6.398653936260428E-03,4.912367890965814E-03,4.211105740919303E-03,4.133647927063595E-03};
      i += eeneneAmp.test_2to2_amp2([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eeneneAmp.test_2to2_amp2_rotations([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eeneneAmp.test_2to2_amp2_boosts([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      i += eeneneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH);
      //std::cout<<"########### me=0.0005, MW=80.385, pspatial=0.001\n";
      pspatial = 0.001;
      ldouble dataCH2[20] = {1.331179459023632E-20,1.211280105931882E-20,1.097961593744878E-20,9.912239224608997E-21,8.910670920782284E-21,7.974911025951425E-21,7.104959540099225E-21,6.300816463208482E-21,5.562481795261996E-21,4.889955536242563E-21,4.283237686132981E-21,3.742328244916051E-21,3.267227212574570E-21,2.857934589091337E-21,2.514450374449151E-21,2.236774568630810E-21,2.024907171619113E-21,1.878848183396858E-21,1.798597603946844E-21,1.784155433251870E-21};
      i += eeneneAmp.test_2to2_amp2([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eeneneAmp.test_2to2_amp2_rotations([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eeneneAmp.test_2to2_amp2_boosts([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      i += eeneneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH2);
      //std::cout<<"########### me=0.105, MW=0.015, pspatial=0.001\n";
      MW = 0.015;
      eeneneAmp.set_masses(me,MW);
      ldouble dataCH3[20] = {1.082747045738139E-05,9.827217934180198E-06,8.885938012579608E-06,8.003202416927408E-06,7.178585711675719E-06,6.411665281691866E-06,5.702021312719340E-06,5.049236771987604E-06,4.452897388969493E-06,3.912591636284931E-06,3.427910710749766E-06,2.998448514568506E-06,2.623801636669696E-06,2.303569334182817E-06,2.037353514055442E-06,1.824758714809525E-06,1.665392088435636E-06,1.558863382423991E-06,1.504784921931135E-06,1.502771592081145E-06};
      i += eeneneAmp.test_2to2_amp2([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH3);
      i += eeneneAmp.test_2to2_amp2_rotations([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH3);
      i += eeneneAmp.test_2to2_amp2_boosts([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH3);
      i += eeneneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH3);
      //std::cout<<"########### me=0.105, MW=0.0006, pspatial=0.001\n";
      MW = 0.0006;
      eeneneAmp.set_masses(me,MW);
      ldouble dataCH4[20] = {4.830513485750748E+00,2.125621077672104E+00,1.147160338796207E+00,6.928642821388733E-01,4.487284045920820E-01,3.043114378732191E-01,2.129181690018390E-01,1.521440372241083E-01,1.101949988660430E-01,8.040528220728960E-02,5.878939437342853E-02,4.285292033678506E-02,3.097438839754686E-02,2.206559577486331E-02,1.537605361368256E-02,1.037585555627793E-02,6.683085749100616E-03,4.017474519856446E-03,2.170029315639724E-03,9.826690619561341E-04};
      i += eeneneAmp.test_2to2_amp2([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH4);
      i += eeneneAmp.test_2to2_amp2_rotations([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH4);
      i += eeneneAmp.test_2to2_amp2_boosts([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH4);
      i += eeneneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeneneAmp.amp2(); }, me,me,0,0,pspatial,dataCH4);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
