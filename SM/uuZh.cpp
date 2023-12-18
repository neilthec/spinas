
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

//File:  SPINAS/SM/uuZh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uuZh.h"

namespace spinas {

  uuZh::uuZh(const ldouble& echarge, const ldouble& massu, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), mu(massu), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propu(massu,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);  
    p1=particle(mu);
    p2=particle(mu);
    p3=particle(MZ);
    p4=particle(mh);
    //<12>,[12],<23>,[23],<13>,[13]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //[312>,[213>,[343>,[143>,[341>
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
    s343a = sproduct(SQUARE,&p3,&p4,&p3);
    s143a = sproduct(SQUARE,&p1,&p4,&p3);
    s341a = sproduct(SQUARE,&p3,&p4,&p1);
    //Couplings
    preTU = sqrt2*e*e*mu/(4.0*MW*MW*SW*SW);
    gL=-4.0/3.0*SW*SW+1.0;
    gR=-4.0/3.0*SW*SW;
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void uuZh::set_masses(const ldouble& massu, const ldouble& massh, const ldouble& massW){
    mu=massu;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(mu);
    p2.set_mass(mu);
    p3.set_mass(MZ);
    p4.set_mass(mh);
    propu.set_mass(mu);
    propZ.set_mass(MZ);
    //Couplings
    preTU = sqrt2*e*e*mu/(4.0*MW*MW*SW*SW);
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void uuZh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    //[312>,[213>,[343>,[143>,[341>
    s312a.update();
    s213a.update();
    s343a.update();
    s143a.update();
    s341a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propZ.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=propu.denominator(propTP);
    pDenU=propu.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uuZh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;

    //Symmetrize the Z-Boson Spin indices
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      //S-Channel Z
      //preZ = e*e/(4.0*MW*MW*SW*SW);
      //all ingoing: 
      //preZ (-Mu (gLe-gRe) (<12>-[12]) (MZ [33]-[343>) + 2 MZ^2(gLe[23] <13> + gRe[13] <23>)))/(s-MZ^2)
      //34 outgoing:
      //preZ (-Mu (gLe-gRe) (<12>-[12]) (MZ [33]-[343>) - 2 MZ^2(gLe[23] <13> + gRe[13] <23>)))/(s-MZ^2)
      //[33]=0 since it is symmetrized
      //preZ ( Mu (gLe-gRe) (<12>-[12]) [343> - 2 MZ^2(gLe[23] <13> + gRe[13] <23>)))/(s-MZ^2)
      amplitude += normFactor*preZ*( mu*(gL-gR)*(a12a.v(ds1,ds2)-s12s.v(ds1,ds2))*s343a.v(ds3a,ds3b) - 2.0*MZ*MZ*(gL*s23s.v(ds2,ds3a)*a13a.v(ds1,ds3b) + gR*a23a.v(ds2,ds3a)*s13s.v(ds1,ds3b)))/pDenS;
      
      //T-Channel e
      //preTU = e*e*mu/(4.0*MW*MW*SW*SW);
      //all ingoing:
      //preh*(gLe <13> (Mu [23]-[312>+MZ <23>)+gRe [13] (MZ [23]-[213>+Mu <23>)))/(t-Mu^2)
      //34 outgoing:
      //- preh*( gLe <13> (Mu [23]-[312>-MZ <23>) + gRe [13] (Mu <23>-[213>-MZ [23])) )/(t-Mu^2)
      amplitude += -normFactor*preTU*( gL*a13a.v(ds1,ds3a)*(mu*s23s.v(ds2,ds3b)-s312a.v(ds3b,ds2)-MZ*a23a.v(ds2,ds3b)) + gR*s13s.v(ds1,ds3a)*(mu*a23a.v(ds2,ds3b)-s213a.v(ds2,ds3b)-MZ*s23s.v(ds2,ds3b)) )/pDenT;

      //U-Channel e
      //preTU = e*e*mu/(4.0*MW*MW*SW*SW);
      //all ingoing:
      //preh (gLe [23] ([143>+2 Mu <13>)+gRe <23> (2 Mu [13]+[341>) )/(u-Mu^2)
      //34 outgoing:
      //preh (gLe [23] ([143>-2 Mu <13>)+gRe <23> ([341> - 2 Mu [13]) )/(u-Mu^2)
      amplitude += normFactor*preTU*( gL*s23s.v(ds2,ds3a)*(s143a.v(ds1,ds3b)-2.0*mu*a13a.v(ds1,ds3b)) + gR*a23a.v(ds2,ds3a)*(s341a.v(ds3b,ds1)-2.0*mu*s13s.v(ds1,ds3b)) )/pDenU;
    }

    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uuZh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3^2=1/9
    return amp2/36.0;
  }
  



  //  Tests
  int test_uuZh(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> Z , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, mh=125, MW=80.385, pspatial=250\n";
      ldouble mu=0.0042, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      uuZh uuZhAmp = uuZh(EE,mu,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.550798794235204E-04,5.047307213059481E-04,6.377537450600018E-04,7.541488990868183E-04,8.539161765148959E-04,9.370555753786010E-04,1.003567094908711E-03,1.053450734748344E-03,1.086706494716404E-03,1.103334374720276E-03,1.103334374720276E-03,1.086706494716404E-03,1.053450734748344E-03,1.003567094908711E-03,9.370555753786010E-04,8.539161765148959E-04,7.541488990868183E-04,6.377537450600016E-04,5.047307213059481E-04,3.550798794235204E-04};
      i += uuZhAmp.test_2to2_amp2([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH);
      i += uuZhAmp.test_2to2_amp2_rotations([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH);
      i += uuZhAmp.test_2to2_amp2_boosts([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH);
      i += uuZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {1.391701223395004E-03,1.448756846497878E-03,1.499472952217119E-03,1.543849543528141E-03,1.581886621340945E-03,1.613584185995968E-03,1.638942237638352E-03,1.657960776335624E-03,1.670639802120840E-03,1.676979315010122E-03,1.676979315010122E-03,1.670639802120840E-03,1.657960776335624E-03,1.638942237638352E-03,1.613584185995968E-03,1.581886621340945E-03,1.543849543528141E-03,1.499472952217119E-03,1.448756846497878E-03,1.391701223395004E-03};
      i += uuZhAmp.test_2to2_amp2([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH2);
      i += uuZhAmp.test_2to2_amp2_rotations([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH2);
      i += uuZhAmp.test_2to2_amp2_boosts([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH2);
      i += uuZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, mh=125, MW=80.385, pspatial=95\n";
      mu = 125.1;
      mh = 125;
      pspatial = 95;
      uuZhAmp.set_masses(mu,mh,MW);
      ldouble dataCH4[20] = {5.480644350327243E-02,5.144691238463379E-02,4.853967172994656E-02,4.609148619214275E-02,4.407532395348427E-02,4.245634462683409E-02,4.120178883420263E-02,4.028455190577993E-02,3.968426757823585E-02,3.938746654980071E-02,3.938746654980072E-02,3.968426757823584E-02,4.028455190577993E-02,4.120178883420263E-02,4.245634462683410E-02,4.407532395348427E-02,4.609148619214275E-02,4.853967172994657E-02,5.144691238463380E-02,5.480644350327244E-02};
      i += uuZhAmp.test_2to2_amp2([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH4);
      i += uuZhAmp.test_2to2_amp2_rotations([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH4);
      i += uuZhAmp.test_2to2_amp2_boosts([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH4);
      i += uuZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      mu = 125;
      mh = 0.0005;
      pspatial = 125.1;
      uuZhAmp.set_masses(mu,mh,MW);
      ldouble dataCH3[20] = {9.007262096417734E-02,7.799808579949234E-02,6.900526479855870E-02,6.224639672285202E-02,5.714189493582707E-02,5.330427357977629E-02,5.047410140773295E-02,4.847816665652709E-02,4.720360620441878E-02,4.658235237429811E-02,4.658235237429811E-02,4.720360620441878E-02,4.847816665652709E-02,5.047410140773295E-02,5.330427357977629E-02,5.714189493582708E-02,6.224639672285202E-02,6.900526479855870E-02,7.799808579949236E-02,9.007262096417736E-02};
      i += uuZhAmp.test_2to2_amp2([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH3);
      i += uuZhAmp.test_2to2_amp2_rotations([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH3);
      i += uuZhAmp.test_2to2_amp2_boosts([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH3);
      i += uuZhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuZhAmp.amp2(); }, mu,mu,MZ,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
