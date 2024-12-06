
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

//File:  SPINAS/SM/eheh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eheh.h"

namespace spinas {
  //Constructors
  eheh::eheh(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW):
    e(echarge), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW),
    prope(me,0), proph(mh,wh),
    p1(particle(me)), p2(particle(mh)),
    p3(particle(me)), p4(particle(mh)),
    s13s(sproduct(SQUARE,&p1,&p3,2)),
    a13a(sproduct(ANGLE,&p1,&p3,2)),
    s123a(sproduct(SQUARE,&p1,&p2,&p3,2)),
    s321a(sproduct(SQUARE,&p3,&p2,&p1,2))
  {
    prehT = 3*e*e*me*mh*mh/(4*MW*MW*SW*SW);
    prehSU = e*e*me*me/(4.0*MW*MW*SW*SW);
  }
  void eheh::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    prope.set_mass(me);
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    p1.set_mass(me);
    p2.set_mass(mh);
    p3.set_mass(me);
    p4.set_mass(mh);
    prehT = 3.0*e*e*me*mh*mh/(4.0*MW*MW*SW*SW);
    prehSU = e*e*me*me/(4.0*MW*MW*SW*SW);
  }
  void eheh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s13s.update();
    a13a.update();
    s123a.update();
    s321a.update();
    //Propagator Momentum
    ldouble propPS[4], propPT[4], propPU[4];
    for(int j=0;j<4;j++){
      propPS[j] = mom1[j]+mom2[j];
      propPT[j] = mom1[j]-mom3[j];
      propPU[j] = mom1[j]-mom4[j];
    }
    pDenTh = proph.denominator(propPT);
    pDenSe = prope.denominator(propPS);
    pDenUe = prope.denominator(propPU);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eheh::amp(const int& ds1, const int& ds3){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Higgs
    //T Channel
    //prehT = 3*e*e*me*mh*mh/(4*MW*MW*SW*SW);
    //eEhh all ingoing:
    //+ prehT ([12]+<12>)/(s-mh^2)
    //ehEh: 2<->3
    //+ prehT ([13]+<13>)/(t-mh^2)
    //34 out:
    //+ prehT ([13]-<13>)/(t-mh^2)
    amplitude += prehT*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3))/pDenTh;

    //S Channel
    //prehSU = e*e*me*me/(4*MW*MW*SW*SW);
    //eEhh all ingoing:
    //+ prehSU ( 2me([12]+<12>) + [132> + [231> )/(t-me^2)
    //ehEh: 2<->3
    //+ prehSU ( 2me([13]+<13>) + [123> + [321> )/(s-me^2)
    //34 out:
    //+ prehSU ( 2me([13]-<13>) - [123> + [321> )/(s-me^2)
    amplitude += prehSU*( two*me*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3)) - s123a.v(ds1,ds3) + s321a.v(ds3,ds1) )/pDenSe;

    //U Channel
    //eEhh all ingoing:
    //+ prehSU * ( 2me([12]+<12>) - [132> - [231> )/(u-me^2)
    //ehEh: 2<->3
    //+ prehSU * ( 2me([13]+<13>) - [123> - [321> )/(u-me^2)
    //34 out:
    //+ prehSU * ( 2me([13]-<13>) + [123> - [321> )/(u-me^2)
    amplitude += prehSU*( two*me*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3)) + s123a.v(ds1,ds3) - s321a.v(ds3,ds1) )/pDenUe;

    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eheh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-1;j3<=1;j3+=2){
	M = amp(j1,j3);
	amp2 += std::pow(std::abs(M),2);
      }
    //Average over initial spins 1/2
    return amp2/2.0;
  }

  



  //  Tests
  int test_eheh(){
    int n=0;//Number of fails
    std::cout<<"\t* e , h  -> e , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### me=0.0005, mh=125, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, mh=125, wh=0, MW=80.385, SW=0.474;//Set width to 0 for comparison with Feynman diagrams.
      eheh ehehAmp = eheh(0.31333,me,mh,wh,MW,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.050600057306506E-12,2.491224862905813E-12,2.232875617856252E-12,1.948354070987034E-12,1.709479442793550E-12,1.516150110908555E-12,1.359242701842188E-12,1.230360034361087E-12,1.123043949848151E-12,1.032508736263662E-12,9.552138791268814E-13,8.885145192774486E-13,8.304082877035313E-13,7.793576202196993E-13,7.341654526822179E-13,6.938876113735407E-13,6.577703418584389E-13,6.252051732212782E-13,5.956958670234038E-13,5.688338973200450E-13};
      i += ehehAmp.test_2to2_amp2([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH);
      i += ehehAmp.test_2to2_amp2_rotations([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH);
      i += ehehAmp.test_2to2_amp2_boosts([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH);
      i += ehehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH);
      //std::cout<<"\n########### me=0.0005, mh=125, MW=80.385, pspatial=126\n";
      pspatial = 126;
      ldouble dataCH2[20] = {8.412889966649463E-13,1.798948319264609E-12,2.244641557806150E-12,2.440460970614084E-12,2.506969868336024E-12,2.504229124031296E-12,2.463967028436724E-12,2.403685431630976E-12,2.333342025877302E-12,2.258717188519032E-12,2.183197544060331E-12,2.108760585012413E-12,2.036537149446771E-12,1.967142008349460E-12,1.900872787535833E-12,1.837831984578111E-12,1.778002962340040E-12,1.721297821649854E-12,1.667587785878193E-12,1.616722539768436E-12};
      i += ehehAmp.test_2to2_amp2([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH2);
      i += ehehAmp.test_2to2_amp2_rotations([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH2);
      i += ehehAmp.test_2to2_amp2_boosts([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH2);
      i += ehehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH2);
      //std::cout<<"\n########### me=0.1, mh=0.15, MW=80.385, pspatial=0.2\n";
      me=0.10;
      mh=0.15;
      pspatial = 0.2;
      ehehAmp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {3.529499456952278E-13,2.415383046226727E-13,1.902207104895969E-13,1.635010706027575E-13,1.486629047580430E-13,1.403097060040832E-13,1.358839101205679E-13,1.340913589125564E-13,1.342692504853799E-13,1.361108047071021E-13,1.395443696229677E-13,1.446919502533162E-13,1.518828630816734E-13,1.617274427557318E-13,1.752893438926153E-13,1.944660273991694E-13,2.228904216766091E-13,2.683853012883027E-13,3.512814823054526E-13,5.454659082625271E-13};
      i += ehehAmp.test_2to2_amp2([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH4);
      i += ehehAmp.test_2to2_amp2_rotations([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH4);
      i += ehehAmp.test_2to2_amp2_boosts([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH4);
      i += ehehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH4);
      //std::cout<<"\n########### me=0.1, mh=0.15, MW=0.11, pspatial=0.2\n";
      me=0.10;
      mh=0.15;
      MW=0.11;
      pspatial = 0.2;
      ehehAmp.set_masses(me,mh,MW);
      ldouble dataCH5[20] = {1.006566432386016E-01,6.888352088841179E-02,5.424842368123804E-02,4.662833677555567E-02,4.239668867937969E-02,4.001446718555772E-02,3.875228890014571E-02,3.824107707072800E-02,3.829180938787165E-02,3.881699622685900E-02,3.979620340053371E-02,4.126422512250993E-02,4.331497808607276E-02,4.612252163770529E-02,4.999019596666304E-02,5.545913175702459E-02,6.356539200427269E-02,7.653992825823670E-02,1.001808196091335E-01,1.555596423698312E-01};
      i += ehehAmp.test_2to2_amp2([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH5);
      i += ehehAmp.test_2to2_amp2_rotations([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH5);
      i += ehehAmp.test_2to2_amp2_boosts([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH5);
      i += ehehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH5);
      //std::cout<<"\n########### me=0.1, mh=0.15, MW=0.006, pspatial=0.2\n";
      me=0.10;
      mh=0.15;
      MW=0.006;
      pspatial = 0.2;
      ehehAmp.set_masses(me,mh,MW);
      ldouble dataCH6[20] = {1.137124933376825E+04,7.781818127524977E+03,6.128481258618874E+03,5.267634866750854E+03,4.789582707984553E+03,4.520461528269680E+03,4.377872390332047E+03,4.320120442843585E+03,4.325851707159173E+03,4.385182420967922E+03,4.495804120271712E+03,4.661647531008240E+03,4.893322485788513E+03,5.210492587173172E+03,5.647426382314148E+03,6.265255772026211E+03,7.181025496408613E+03,8.646767666889225E+03,1.131749521525712E+04,1.757367842543748E+04};
      i += ehehAmp.test_2to2_amp2([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH6);
      i += ehehAmp.test_2to2_amp2_rotations([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH6);
      i += ehehAmp.test_2to2_amp2_boosts([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH6);
      i += ehehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ehehAmp.amp2(); }, me,mh,me,mh,pspatial,dataCH6);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
