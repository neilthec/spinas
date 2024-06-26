
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

//File:  SPINAS/SM/ddZZ.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ddZZ.h"

namespace spinas {

  ddZZ::ddZZ(const ldouble& echarge, const ldouble& massf, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), Qf(-1.0/3.0), mf(massf), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propf(massf,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    proph = propagator(mh,wh);  
    p1=particle(mf);
    p2=particle(mf);
    p3=particle(MZ);
    p4=particle(MZ);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    //[314>,[413>
    s314a = sproduct(SQUARE,&p3,&p1,&p4);
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
    gL=-2.0*Qf*SW*SW-1.0;
    gR=-2.0*Qf*SW*SW;
    preh = e*e*mf/(2.0*MW*MW*SW*SW);
  }
  void ddZZ::set_masses(const ldouble& massf, const ldouble& massh, const ldouble& massW){
    mf=massf;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(mf);
    p2.set_mass(mf);
    p3.set_mass(MZ);
    p4.set_mass(MZ);
    propf.set_mass(mf);
    proph.set_mass(mh);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
    preh = e*e*mf/(2.0*MW*MW*SW*SW);
  }
  void ddZZ::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s34s.update();
    a34a.update();
    s24s.update();
    a24a.update();
    s14s.update();
    a14a.update();
    //[314>,[413>
    s314a.update();
    s413a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=proph.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=propf.denominator(propTP);
    pDenU=propf.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ddZZ::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);
      
      //S-Channel h
      //preh = e*e*mf/(4.0*MW*MW*SW*SW);
      //all ingoing: 
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //34 outgoing:
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      amplitude += normFactor*preh*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))/pDenS;
      
      //T-Channel e
      //preTU = e*e/(4.0*MW*MW*SW*SW);
      //all ingoing:
      //+preTU gLe^2 [24] <13> (MZ <34>+[314>))/(Mf^2-t)
      //+preTU gRe^2 [13] <24> (MZ [34]+[413>))/(Mf^2-t)
      //+preTU gLe gRe Mf ([13] [24] <34>+[34] <13> <24>))/(Mf^2-t)
      //34 outgoing:
      //+preTU gLe^2 [24] <13> (MZ <34>-[314>))/(t-Mf^2)
      //+preTU gRe^2 [13] <24> (MZ [34]-[413>))/(t-Mf^2)
      //-preTU gLe gRe Mf ([13] [24] <34>+[34] <13> <24>))/(t-Mf^2)
      amplitude += normFactor*preTU*gL*gL*s24s.v(ds2,ds4a)*a13a.v(ds1,ds3a)*(MZ*a34a.v(ds3b,ds4b)-s314a.v(ds3b,ds4b))/pDenT
      	+          normFactor*preTU*gR*gR*a24a.v(ds2,ds4a)*s13s.v(ds1,ds3a)*(MZ*s34s.v(ds3b,ds4b)-s413a.v(ds4b,ds3b))/pDenT
      	-          normFactor*preTU*gL*gR*mf*(s13s.v(ds1,ds3a)*s24s.v(ds2,ds4a)*a34a.v(ds3b,ds4b)+a13a.v(ds1,ds3a)*a24a.v(ds2,ds4a)*s34s.v(ds3b,ds4b))/pDenT;

      //U-Channel e
      //preTU = e*e/(4.0*MW*MW*SW*SW);
      //all ingoing:
      //+preTU gLe^2 [23] <14> (MZ <34>-[413>)/(u-Mf^2)
      //+preTU gRe^2 [14] <23> (MZ [34]-[314>)/(u-Mf^2)
      //+preTU gLe gRe Mf ([14] [23] <34>+[34] <14> <23>)/(u-Mf^2)
      //34 outgoing:
      //-preTU gLe^2 [23] <14> (MZ <34>+[413>)/(u-Mf^2)
      //-preTU gRe^2 [14] <23> (MZ [34]+[314>)/(u-Mf^2)
      //+preTU gLe gRe Mf ([14] [23] <34>+[34] <14> <23>)/(u-Mf^2)
      amplitude += -normFactor*preTU*gL*gL*s23s.v(ds2,ds3a)*a14a.v(ds1,ds4a)*(MZ*a34a.v(ds3b,ds4b)+s413a.v(ds4b,ds3b))/pDenU
      	-           normFactor*preTU*gR*gR*s14s.v(ds1,ds4a)*a23a.v(ds2,ds3a)*(MZ*s34s.v(ds3b,ds4b)+s314a.v(ds3b,ds4b))/pDenU
      	+           normFactor*preTU*gL*gR*mf*(s14s.v(ds1,ds4a)*s23s.v(ds2,ds3a)*a34a.v(ds3b,ds4b)+a14a.v(ds1,ds4a)*a23a.v(ds2,ds3a)*s34s.v(ds3b,ds4b))/pDenU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ddZZ::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3^2=1/9
    //Symmetry factor 1/2
    return amp2/72.0;
  }
  



  //  Tests
  int test_ddZZ(){
    int n=0;//Number of fails
    std::cout<<"\t* d , D  -> Z , Z       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mf=0.0075, mh=125, MW=80.385, pspatial=250\n";
      ldouble mf=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      ddZZ ddZZAmp = ddZZ(EE,mf,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {7.149729867530870E-02,2.454165267978785E-02,1.452242665603786E-02,1.023440475684747E-02,7.910268695675113E-03,6.501381043899090E-03,5.602959572285331E-03,5.028742857465567E-03,4.684711841611840E-03,4.522944201974127E-03,4.522944201974161E-03,4.684711841611831E-03,5.028742857465593E-03,5.602959572285334E-03,6.501381043899070E-03,7.910268695675090E-03,1.023440475684746E-02,1.452242665603786E-02,2.454165267978785E-02,7.149729867530859E-02};
      i += ddZZAmp.test_2to2_amp2([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH);
      i += ddZZAmp.test_2to2_amp2_rotations([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH);
      i += ddZZAmp.test_2to2_amp2_boosts([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH);
      i += ddZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH);
      //std::cout<<"\n# mf=0.0075, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {4.281460910017963E-02,2.814829131284285E-02,2.086280221229616E-02,1.673301562583663E-02,1.415524442257640E-02,1.245731531657764E-02,1.131614856461435E-02,1.056152147858391E-02,1.009960618588336E-02,9.879813823343669E-03,9.879813823343669E-03,1.009960618588337E-02,1.056152147858391E-02,1.131614856461434E-02,1.245731531657765E-02,1.415524442257640E-02,1.673301562583662E-02,2.086280221229615E-02,2.814829131284285E-02,4.281460910017960E-02};
      i += ddZZAmp.test_2to2_amp2([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH2);
      i += ddZZAmp.test_2to2_amp2_rotations([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH2);
      i += ddZZAmp.test_2to2_amp2_boosts([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH2);
      i += ddZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH2);
      //std::cout<<"\n# mf=125.1, mh=125, MW=80.385, pspatial=95\n";
      mf = 125.1;
      mh = 125;
      pspatial = 95;
      ddZZAmp.set_masses(mf,mh,MW);
      ldouble dataCH4[20] = {2.591865116447664E-02,2.127971852415857E-02,1.764087623876426E-02,1.477550219515841E-02,1.252585566304810E-02,1.078097100956944E-02,9.462502452798826E-03,8.515655192601375E-03,7.903380963831030E-03,7.602733582785445E-03,7.602733582785447E-03,7.903380963831029E-03,8.515655192601371E-03,9.462502452798829E-03,1.078097100956944E-02,1.252585566304810E-02,1.477550219515840E-02,1.764087623876426E-02,2.127971852415857E-02,2.591865116447663E-02};
      i += ddZZAmp.test_2to2_amp2([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH4);
      i += ddZZAmp.test_2to2_amp2_rotations([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH4);
      i += ddZZAmp.test_2to2_amp2_boosts([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH4);
      i += ddZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH4);
      //std::cout<<"\n# mf=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      mf = 125;
      mh = 0.0005;
      pspatial = 125.1;
      ddZZAmp.set_masses(mf,mh,MW);
      ldouble dataCH3[20] = {3.585146370646807E-02,2.802690355605272E-02,2.235783513028667E-02,1.816094152210952E-02,1.501965562470507E-02,1.267115678569556E-02,1.094532052518677E-02,9.730938257367582E-03,8.956525911213341E-03,8.579337054721426E-03,8.579337054721416E-03,8.956525911213343E-03,9.730938257367577E-03,1.094532052518676E-02,1.267115678569556E-02,1.501965562470507E-02,1.816094152210951E-02,2.235783513028667E-02,2.802690355605271E-02,3.585146370646806E-02};
      i += ddZZAmp.test_2to2_amp2([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH3);
      i += ddZZAmp.test_2to2_amp2_rotations([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH3);
      i += ddZZAmp.test_2to2_amp2_boosts([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH3);
      i += ddZZAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddZZAmp.amp2(); }, mf,mf,MZ,MZ,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
