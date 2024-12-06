
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

//File:  SPINAS/SM/dZZd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dZZd.h"

namespace spinas {

  dZZd::dZZd(const ldouble& echarge, const ldouble& massf, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), Qf(-1.0/3.0), mf(massf), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propf(massf,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    proph = propagator(mh,wh);  
    p1=particle(mf);
    p2=particle(MZ);
    p3=particle(MZ);
    p4=particle(mf);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<24>,[24],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s34s = sproduct(SQUARE,&p3,&p4,2);
    a34a = sproduct(ANGLE,&p3,&p4,2);
    s24s = sproduct(SQUARE,&p2,&p4,2);
    a24a = sproduct(ANGLE,&p2,&p4,2);
    s14s = sproduct(SQUARE,&p1,&p4,2);
    a14a = sproduct(ANGLE,&p1,&p4,2);
    //[312>,[213>
    s312a = sproduct(SQUARE,&p3,&p1,&p2,2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3,2);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
    gL=-2.0*Qf*SW*SW-1.0;
    gR=-2.0*Qf*SW*SW;
    preh = e*e*mf/(2.0*MW*MW*SW*SW);
  }
  void dZZd::set_masses(const ldouble& massf, const ldouble& massh, const ldouble& massW){
    mf=massf;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(mf);
    p2.set_mass(MZ);
    p3.set_mass(MZ);
    p4.set_mass(mf);
    propf.set_mass(mf);
    proph.set_mass(mh);
    //Couplings
    preTU = e*e/(2.0*MW*MW*SW*SW);
    preh = e*e*mf/(2.0*MW*MW*SW*SW);
  }
  void dZZd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    //[312>,[213>
    s312a.update();
    s213a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propf.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=propf.denominator(propTP);
    pDenU=proph.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble dZZd::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds2a, ds2b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds2);
    ldouble normFactor=get_spin_normalization(ds3,ds2);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds2,ds2a,ds2b, i);
      
      //S-Channel h
      //preh = e*e*mf/(4.0*MW*MW*SW*SW);
      //dDZZ all ingoing: 
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //dZZD: 2<->4
      //preh [23] <23> ([14]+<14>)/(u-Mh^2)
      //34 out:
      //- preh [23] <23> ([14]-<14>)/(u-Mh^2)
      amplitude += - normFactor*preh*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(s14s.v(ds1,ds4)-a14a.v(ds1,ds4))/pDenU;
      
      //T-Channel e
      //preTU = e*e/(4.0*MW*MW*SW*SW);
      //dDZZ all ingoing:
      //+preTU gLe^2 [24] <13> (MZ <34>+[314>))/(Mf^2-t)
      //+preTU gRe^2 [13] <24> (MZ [34]+[413>))/(Mf^2-t)
      //+preTU gLe gRe Mf ([13] [24] <34>+[34] <13> <24>))/(Mf^2-t)
      //dZZD: 2<->4
      //+preTU gLe^2 [24] <13> (-MZ <23>+[312>))/(t-mf^2)
      //+preTU gRe^2 [13] <24> (-MZ [23]+[213>))/(t-mf^2)
      //-preTU gLe gRe Mf ([13] [24] <23>+[23] <13> <24>))/(t-mf^2)
      //34 out:
      //-preTU gLe^2 [24] <13> (MZ <23>+[312>))/(t-mf^2)
      //+preTU gRe^2 [13] <24> (MZ [23]+[213>))/(t-mf^2)
      //+preTU gLe gRe Mf ([13] [24] <23>-[23] <13> <24>))/(t-mf^2)
      amplitude +=
	- normFactor*preTU*gL*gL*s24s.v(ds2a,ds4)*a13a.v(ds1,ds3a)*(MZ*a23a.v(ds2b,ds3b)+s312a.v(ds3b,ds2b))/pDenT
      	+ normFactor*preTU*gR*gR*a24a.v(ds2a,ds4)*s13s.v(ds1,ds3a)*(MZ*s23s.v(ds2b,ds3b)+s213a.v(ds2b,ds3b))/pDenT
      	+ normFactor*preTU*gL*gR*mf*(s13s.v(ds1,ds3a)*s24s.v(ds2a,ds4)*a23a.v(ds2b,ds3b)-a13a.v(ds1,ds3a)*a24a.v(ds2a,ds4)*s23s.v(ds2b,ds3b))/pDenT;

      //U-Channel e
      //preTU = e*e/(4.0*MW*MW*SW*SW);
      //dDZZ all ingoing:
      //+preTU gLe^2 [23] <14> (MZ <34>-[413>)/(u-Mf^2)
      //+preTU gRe^2 [14] <23> (MZ [34]-[314>)/(u-Mf^2)
      //+preTU gLe gRe Mf ([14] [23] <34>+[34] <14> <23>)/(u-Mf^2)
      //dZZD: 2<->4
      //+preTU gLe^2 [34] <12> (MZ <23>+[213>)/(s-Mf^2)
      //+preTU gRe^2 [12] <34> (MZ [23]+[312>)/(s-Mf^2)
      //+preTU gLe gRe Mf ([12] [34] <23>+[23] <12> <34>)/(s-Mf^2)
      //34 out:
      //-preTU gLe^2 [34] <12> (MZ <23>+[213>)/(s-Mf^2)
      //+preTU gRe^2 [12] <34> (MZ [23]+[312>)/(s-Mf^2)
      //+preTU gLe gRe Mf (-[12] [34] <23>+[23] <12> <34>)/(s-Mf^2)
      amplitude +=
	- normFactor*preTU*gL*gL*s34s.v(ds3a,ds4)*a12a.v(ds1,ds2a)*(MZ*a23a.v(ds2b,ds3b)+s213a.v(ds2b,ds3b))/pDenS
      	+ normFactor*preTU*gR*gR*s12s.v(ds1,ds2a)*a34a.v(ds3a,ds4)*(MZ*s23s.v(ds2b,ds3b)+s312a.v(ds3b,ds2b))/pDenS
      	+ normFactor*preTU*gL*gR*mf*(-s12s.v(ds1,ds2a)*s34s.v(ds3a,ds4)*a23a.v(ds2b,ds3b)+a12a.v(ds1,ds2a)*a34a.v(ds3a,ds4)*s23s.v(ds2b,ds3b))/pDenS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble dZZd::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2*1/3=1/6
    //Average over initial colors 1/3
    return amp2/18.0;
  }
  



  //  Tests
  int test_dZZd(){
    int n=0;//Number of fails
    std::cout<<"\t* d , Z  -> Z , d       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mf=0.0075, mh=125, MW=80.385, pspatial=250\n";
      ldouble mf=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      dZZd dZZdAmp = dZZd(EE,mf,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.862237645398273E-01,9.080184529280307E-02,5.470775292749400E-02,3.968787741355533E-02,3.155518147449662E-02,2.652280150423533E-02,2.315006479194621E-02,2.076918682332629E-02,1.902827142703494E-02,1.772428036571513E-02,1.673184993618968E-02,1.596938940692239E-02,1.538150977157856E-02,1.492928754399172E-02,1.458456962974838E-02,1.432648906390158E-02,1.413925115566979E-02,1.401068059106300E-02,1.393124087149531E-02,1.389335612977252E-02};
      i += dZZdAmp.test_2to2_amp2([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH);
      i += dZZdAmp.test_2to2_amp2_rotations([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH);
      i += dZZdAmp.test_2to2_amp2_boosts([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH);
      i += dZZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH);
      //std::cout<<"\n# mf=0.0075, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {1.578134199098067E+00,1.312540771346248E-01,6.563231565130595E-02,4.449694715376095E-02,3.429090261406786E-02,2.836940776437320E-02,2.455725963992628E-02,2.193627912881013E-02,2.005293467667007E-02,1.865773064180079E-02,1.760217207827179E-02,1.679239464474251E-02,1.616617596866774E-02,1.568064739528200E-02,1.530531654309293E-02,1.501790399596761E-02,1.480175392079488E-02,1.464416583339183E-02,1.453528673935436E-02,1.446735575494091E-02};
      i += dZZdAmp.test_2to2_amp2([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH2);
      i += dZZdAmp.test_2to2_amp2_rotations([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH2);
      i += dZZdAmp.test_2to2_amp2_boosts([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH2);
      i += dZZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH2);
      //std::cout<<"\n# mf=125.1, mh=125, MW=80.385, pspatial=95\n";
      mf = 125.1;
      mh = 125;
      pspatial = 95;
      dZZdAmp.set_masses(mf,mh,MW);
      ldouble dataCH4[20] = {2.709868693140395E-01,2.559572000737879E-01,2.439743006358444E-01,2.343809321330275E-01,2.267166808968272E-01,2.206557710556637E-01,2.159685067155196E-01,2.124971885629665E-01,2.101413967983937E-01,2.088498161595859E-01,2.086171712178791E-01,2.094858135706063E-01,2.115523520460445E-01,2.149806880049101E-01,2.200242088894815E-01,2.270621919490070E-01,2.366596086978057E-01,2.496674418635549E-01,2.673966515719223E-01,2.919333098443554E-01};
      i += dZZdAmp.test_2to2_amp2([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH4);
      i += dZZdAmp.test_2to2_amp2_rotations([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH4);
      i += dZZdAmp.test_2to2_amp2_boosts([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH4);
      i += dZZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH4);
      //std::cout<<"\n# mf=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      mf = 125;
      mh = 0.0005;
      pspatial = 125.1;
      dZZdAmp.set_masses(mf,mh,MW);
      ldouble dataCH3[20] = {4.305955448543033E-01,3.940334449717797E-01,3.691759007086457E-01,3.524930894016898E-01,3.419845093681995E-01,3.365261181244348E-01,3.355623514053744E-01,3.389704645200742E-01,3.470309992246364E-01,3.604894614736421E-01,3.807333220621787E-01,4.101630824960047E-01,4.529500917337088E-01,5.166620256008919E-01,6.160737256144868E-01,7.833099475291655E-01,1.100161974668830E+00,1.833378534107305E+00,4.351055699357554E+00,3.383224666128863E+01};
      i += dZZdAmp.test_2to2_amp2([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH3);
      i += dZZdAmp.test_2to2_amp2_rotations([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH3);
      i += dZZdAmp.test_2to2_amp2_boosts([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH3);
      i += dZZdAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dZZdAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
