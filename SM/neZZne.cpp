
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

//File:  SPINAS/SM/neZZne.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/neZZne.h"

namespace spinas {

  neZZne::neZZne(const ldouble& echarge, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propne(0,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    p1=particle(0);
    p2=particle(MZ);
    p3=particle(MZ);
    p4=particle(0);
    //[23], [24], <13>, <14>, <34>
    s34s = sproduct(SQUARE,&p3,&p4);
    s24s = sproduct(SQUARE,&p2,&p4);
    a13a = sproduct(ANGLE,&p1,&p3);
    a12a = sproduct(ANGLE,&p1,&p2);
    a23a = sproduct(ANGLE,&p2,&p3);
    //[314>, [413>
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
  }
  void neZZne::set_masses(const ldouble& massW){
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p2.set_mass(MZ);
    p3.set_mass(MZ);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
  }
  void neZZne::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //[23], [24], <13>, <14>, <34>
    s34s.update();
    s24s.update();
    a13a.update();
    a12a.update();
    a23a.update();
    //[314>, [413>
    s312a.update();
    s213a.update();
    //Propagator Momentum
    ldouble propTP[4], propSP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propSP[j] = mom1[j]+mom2[j];
    }
    pDenT=propne.den(propTP);
    pDenS=propne.den(propSP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble neZZne::amp(const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    int ds2a, ds2b, ds3a, ds3b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds2,ds3);
    ldouble normFactor=get_spin_normalization(ds2,ds3);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, ds3,ds3a,ds3b, i);
      
      //T-Channel e
      //preTU = e*e/(4.0*MW*MW*SW*SW);
      //nNZZ all ingoing:
      //+preTU [24] <13> (MZ <34>+[314>))/(Me^2-t)
      //nZZN: 2<->4
      //+preTU [24] <13> (-MZ <23>+[312>))/(t-Me^2)
      //34 out:
      //-preTU [24] <13> (MZ <23>+[312>))/(t-Me^2)
      amplitude += - normFactor*preTU*s24s.v(ds2a)*a13a.v(ds3a)*(MZ*a23a.v(ds2b,ds3b)+s312a.v(ds3b,ds2b))/pDenT;

      //S-Channel e
      //preTU = e*e/(4.0*MW*MW*SW*SW);
      //nNZZ all ingoing:
      //+preTU [23] <14> (MZ <34>-[413>)/(u-Me^2)
      //nZZN: 2<->4
      //-preTU [34] <12> (-MZ <23>-[213>)/(s-Me^2)
      //34 out:
      //-preTU [34] <12> (MZ <23>+[213>)/(s-Me^2)
      amplitude += -normFactor*preTU*s34s.v(ds3a)*a12a.v(ds2a)*(MZ*a23a.v(ds2b,ds3b)+s213a.v(ds2b,ds3b))/pDenS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble neZZne::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-2;j3<=2;j3+=2){
	M = amp(j2,j3);
	amp2 += std::pow(std::abs(M),2);
      }
    //Average over initial spins 1/3
    return amp2/3.0;
  }
  



  //  Tests
  int test_neZZne(){
    int n=0;//Number of fails
    std::cout<<"\t* ne, Z  -> Z , ne      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      neZZne neZZneAmp = neZZne(EE,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.094461419151703E+00,3.472077700279796E-01,2.091913086041192E-01,1.517583625031387E-01,1.206605890070925E-01,1.014177925884506E-01,8.852113409645064E-02,7.941714150148020E-02,7.276023548775576E-02,6.777403917338271E-02,6.397918728135561E-02,6.106369346642204E-02,5.881576144609052E-02,5.708655573495644E-02,5.576842463338999E-02,5.478157711851412E-02,5.406561751369150E-02,5.357398990895435E-02,5.327022870954217E-02,5.312536548573234E-02};
      i += neZZneAmp.test_2to2_amp2([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH);
      i += neZZneAmp.test_2to2_amp2_rotations([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH);
      i += neZZneAmp.test_2to2_amp2_boosts([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH);
      i += neZZneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH);
      //std::cout<<"\n# MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {6.034465206784677E+00,5.018888867307031E-01,2.509646168386984E-01,1.701472681000258E-01,1.311214308073457E-01,1.084788397214498E-01,9.390196125158647E-02,8.387986530575652E-02,7.667833951482834E-02,7.134336332197493E-02,6.730712227757267E-02,6.421069824614506E-02,6.181616554132763E-02,5.995960315492051E-02,5.852441442940950E-02,5.742540741632385E-02,5.659889360608122E-02,5.599630894208470E-02,5.557997759648878E-02,5.532022356277697E-02};
      i += neZZneAmp.test_2to2_amp2([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH2);
      i += neZZneAmp.test_2to2_amp2_rotations([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH2);
      i += neZZneAmp.test_2to2_amp2_boosts([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH2);
      i += neZZneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH2);
      //std::cout<<"\n# MW=80.385, pspatial=95\n";
      pspatial = 95;
      ldouble dataCH4[20] = {7.456798714137383E+00,1.581221827536240E+00,4.452740563750762E-01,2.489067218942878E-01,1.743583610950463E-01,1.364220981311065E-01,1.139144844349525E-01,9.925546888889708E-02,8.910210672601479E-02,8.176301573807961E-02,7.629466672937971E-02,7.213079171579097E-02,6.891165745496576E-02,6.639809884625641E-02,6.442507661265942E-02,6.287511224151682E-02,6.166235575189644E-02,6.072264124559407E-02,6.000706441321789E-02,5.947771043468283E-02};
      i += neZZneAmp.test_2to2_amp2([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH4);
      i += neZZneAmp.test_2to2_amp2_rotations([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH4);
      i += neZZneAmp.test_2to2_amp2_boosts([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH4);
      i += neZZneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH4);
      //std::cout<<"\n# MW=8.0, pspatial=125.1\n";
      pspatial = 125.1;
      MW=8;
      MZ=MW/CW;
      neZZneAmp.set_masses(MW);
      ldouble dataCH3[20] = {1.056949431485677E+00,3.541366604687709E-01,2.146457769190266E-01,1.556198540533428E-01,1.234147556566583E-01,1.034007561099862E-01,8.995113419990319E-02,8.044011452158384E-02,7.347759079573257E-02,6.825876893471727E-02,6.428544510510059E-02,6.123270481243348E-02,5.887960677879667E-02,5.707067414447375E-02,5.569332042551545E-02,5.466401384505959E-02,5.391947283650067E-02,5.341088025197682E-02,5.309997404159650E-02,5.295634080669820E-02};
      i += neZZneAmp.test_2to2_amp2([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH3);
      i += neZZneAmp.test_2to2_amp2_rotations([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH3);
      i += neZZneAmp.test_2to2_amp2_boosts([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH3);
      i += neZZneAmp.test_2to2_amp2_boosts_and_rotations([&]() { return neZZneAmp.amp2(); }, 0,MZ,MZ,0,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
