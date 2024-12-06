
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

//File:  SPINAS/SM/ZZWW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/ZZWW.h"

namespace spinas {

  ZZWW::ZZWW(const ldouble& echarge, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), mh(massh), wh(widthh), MW(massW), WW(widthW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW) {
    CW = std::sqrt(1.0-sinW*sinW);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    proph = propagator(mh,wh);
    propW = propagator(MW,WW);
    p1=particle(MZ);
    p2=particle(MZ);
    p3=particle(MW);
    p4=particle(MW);
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s34s = sproduct(SQUARE,&p3,&p4,2);
    a34a = sproduct(ANGLE,&p3,&p4,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s24s = sproduct(SQUARE,&p2,&p4,2);
    a24a = sproduct(ANGLE,&p2,&p4,2);
    s14s = sproduct(SQUARE,&p1,&p4,2);
    a14a = sproduct(ANGLE,&p1,&p4,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    //
    s431a = sproduct(SQUARE,&p4,&p3,&p1,2);
    s134a = sproduct(SQUARE,&p1,&p3,&p4,2);
    s341a = sproduct(SQUARE,&p3,&p4,&p1,2);
    s143a = sproduct(SQUARE,&p1,&p4,&p3,2);
    //Couplings
    preh = e*e/(MW*MW*SW*SW);
    preW = e*e/(2.0*MW*MW*MW*SW*SW);
  }
  void ZZWW::set_masses(const ldouble& massh, const ldouble& massW){
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(MZ);
    p2.set_mass(MZ);
    p3.set_mass(MW);
    p4.set_mass(MW);
    proph.set_mass(mh);
    propW.set_mass(MW);
    //Couplings
    preh = e*e/(MW*MW*SW*SW);
    preW = e*e/(2.0*MW*MW*MW*SW*SW);
  }
  void ZZWW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<24>,[24],<34>,[34],<14>,[14],<13>,[13]
    s12s.update();
    a12a.update();
    s34s.update();
    a34a.update();
    s23s.update();
    a23a.update();
    s24s.update();
    a24a.update();
    s14s.update();
    a14a.update();
    s13s.update();
    a13a.update();
    //
    s431a.update();
    s134a.update();
    s341a.update();
    s143a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=proph.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=propW.denominator(propTP);
    pDenU=propW.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble ZZWW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    constexpr ldouble one=1, two=2, three=3, four=4, six=6;
    int ds1a, ds1b, ds2a, ds2b, ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds1,ds2,ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds1,ds2,ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds1,ds1a,ds1b, ds2,ds2a,ds2b, ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);

      
      //S-Channel h
      //preh = e*e/(MW*MW*SW*SW);
      //all ingoing = 34 outgoing:
      //preh [12]<12>[34]<34>/(s-Mh^2)
      amplitude += normFactor*preh*s12s.v(ds1a,ds2a)*a12a.v(ds1b,ds2b)*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)/pDenS;

      //T-Channel W
      //preW = e*e/(2.0*MW*MW*MW*SW*SW);
      //No Mandelstahm s
      amplitude += normFactor*preW*(
        two*MW*a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s13s.v(ds1b,ds3b)*s24s.v(ds2b,ds4b)
        +(four*MW*a23a.v(ds2a,ds3a)*a24a.v(ds2b,ds4a)*s13s.v(ds1a,ds3b)*s14s.v(ds1b,ds4b)
        +two*MW*a12a.v(ds1a,ds2b)*a34a.v(ds3a,ds4a)*s13s.v(ds1b,ds3b)*s24s.v(ds2a,ds4b)
        +three*MW*a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)
        +three*MW*a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s13s.v(ds1b,ds3b)*s24s.v(ds2b,ds4b)
        -six*MW*a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s13s.v(ds1b,ds3b)*s24s.v(ds2b,ds4b)
        +four*MW*a13a.v(ds1a,ds3a)*a14a.v(ds1b,ds4a)*s23s.v(ds2a,ds3b)*s24s.v(ds2b,ds4b))*CW*CW
        +(four*MW*a24a.v(ds2a,ds4a)*a34a.v(ds3a,ds4b)*s12s.v(ds1a,ds2b)*s13s.v(ds1b,ds3b)
        -three*MW*a13a.v(ds1a,ds3a)*a23a.v(ds2a,ds3b)*a24a.v(ds2b,ds4a)*s14s.v(ds1b,ds4b)
        -MW*a24a.v(ds2a,ds4a)*s13s.v(ds1a,ds3a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)
        -three*MW*a13a.v(ds1a,ds3a)*a14a.v(ds1b,ds4a)*a23a.v(ds2a,ds3b)*s24s.v(ds2b,ds4b)
        -MW*a14a.v(ds1a,ds4a)*s13s.v(ds1b,ds3a)*s23s.v(ds2a,ds3b)*s24s.v(ds2b,ds4b)
        +four*MW*a12a.v(ds1a,ds2b)*a13a.v(ds1b,ds3b)*s24s.v(ds2a,ds4a)*s34s.v(ds3a,ds4b)
        -a23a.v(ds2a,ds3a)*s13s.v(ds1a,ds3b)*s24s.v(ds2b,ds4a)*s134a.v(ds1b,ds4b)
        -a13a.v(ds1a,ds3a)*s23s.v(ds2a,ds3b)*s24s.v(ds2b,ds4a)*s134a.v(ds1b,ds4b)
        -three*a23a.v(ds2a,ds3a)*a24a.v(ds2b,ds4a)*s13s.v(ds1a,ds3b)*s431a.v(ds4b,ds1b)
        +-three*a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s23s.v(ds2b,ds3b)*s431a.v(ds4b,ds1b))*CW*CW*CW
        +(-two*MW*a12a.v(ds1a,ds2b)*a34a.v(ds3a,ds4a)*s13s.v(ds1b,ds3b)*s24s.v(ds2a,ds4b)
        +three*MW*a14a.v(ds1a,ds4a)*a24a.v(ds2a,ds4b)*s13s.v(ds1b,ds3a)*s23s.v(ds2b,ds3b)
        +MW*a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)
        +MW*a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s13s.v(ds1b,ds3b)*s24s.v(ds2b,ds4b)
        +two*MW*a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s13s.v(ds1b,ds3b)*s24s.v(ds2b,ds4b)
        +MW*a13a.v(ds1a,ds3a)*a23a.v(ds2a,ds3b)*s14s.v(ds1b,ds4a)*s24s.v(ds2b,ds4b)
        +three*a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s23s.v(ds2b,ds3b)*s134a.v(ds1b,ds4b)
        +three*a13a.v(ds1a,ds3a)*a23a.v(ds2a,ds3b)*s24s.v(ds2b,ds4a)*s134a.v(ds1b,ds4b)
        +a24a.v(ds2a,ds4a)*s13s.v(ds1a,ds3a)*s23s.v(ds2b,ds3b)*s431a.v(ds4b,ds1b)
        +a23a.v(ds2a,ds3a)*s13s.v(ds1a,ds3b)*s24s.v(ds2b,ds4a)*s431a.v(ds4b,ds1b))*CW*CW*CW*CW
        )/pDenT;

      //U-Channel W
      //No Mandelstahm s
      amplitude += normFactor*preW*(
        two*MW*a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)
        +(four*MW*a23a.v(ds2a,ds3a)*a24a.v(ds2b,ds4a)*s13s.v(ds1a,ds3b)*s14s.v(ds1b,ds4b)
        -two*MW*a12a.v(ds1a,ds2b)*a34a.v(ds3a,ds4a)*s14s.v(ds1b,ds4b)*s23s.v(ds2a,ds3b)
        -six*MW*a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)
        +three*MW*a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)
        +three*MW*a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s13s.v(ds1b,ds3b)*s24s.v(ds2b,ds4b)
        +four*MW*a13a.v(ds1a,ds3a)*a14a.v(ds1b,ds4a)*s23s.v(ds2a,ds3b)*s24s.v(ds2b,ds4b))*CW*CW
        +(-three*MW*a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*a24a.v(ds2b,ds4b)*s13s.v(ds1b,ds3b)
        -four*MW*a23a.v(ds2a,ds3b)*a34a.v(ds3a,ds4a)*s12s.v(ds1a,ds2b)*s14s.v(ds1b,ds4b)
        -three*MW*a13a.v(ds1a,ds3a)*a14a.v(ds1b,ds4a)*a24a.v(ds2a,ds4b)*s23s.v(ds2b,ds3b)
        -MW*a23a.v(ds2a,ds3a)*s13s.v(ds1a,ds3b)*s14s.v(ds1b,ds4a)*s24s.v(ds2b,ds4b)
        -MW*a13a.v(ds1a,ds3a)*s14s.v(ds1b,ds4a)*s23s.v(ds2a,ds3b)*s24s.v(ds2b,ds4b)
        -four*MW*a12a.v(ds1a,ds2b)*a14a.v(ds1b,ds4a)*s23s.v(ds2a,ds3b)*s34s.v(ds3a,ds4b)
        -a24a.v(ds2a,ds4a)*s14s.v(ds1a,ds4b)*s23s.v(ds2b,ds3a)*s143a.v(ds1b,ds3b)
        -a14a.v(ds1a,ds4a)*s23s.v(ds2a,ds3a)*s24s.v(ds2b,ds4b)*s143a.v(ds1b,ds3b)
        -three*a23a.v(ds2a,ds3a)*a24a.v(ds2b,ds4a)*s14s.v(ds1a,ds4b)*s341a.v(ds3b,ds1b)
        -three*a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s24s.v(ds2b,ds4b)*s341a.v(ds3b,ds1b))*CW*CW*CW
        +(two*MW*a12a.v(ds1a,ds2b)*a34a.v(ds3a,ds4a)*s14s.v(ds1b,ds4b)*s23s.v(ds2a,ds3b)
        +MW*a14a.v(ds1a,ds4a)*a24a.v(ds2a,ds4b)*s13s.v(ds1b,ds3a)*s23s.v(ds2b,ds3b)
        +two*MW*a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)
        +MW*a13a.v(ds1a,ds3a)*a24a.v(ds2a,ds4a)*s14s.v(ds1b,ds4b)*s23s.v(ds2b,ds3b)
        +MW*a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s13s.v(ds1b,ds3b)*s24s.v(ds2b,ds4b)
        +three*MW*a13a.v(ds1a,ds3a)*a23a.v(ds2a,ds3b)*s14s.v(ds1b,ds4a)*s24s.v(ds2b,ds4b)
        +three*a14a.v(ds1a,ds4a)*a24a.v(ds2a,ds4b)*s23s.v(ds2b,ds3a)*s143a.v(ds1b,ds3b)
        +three*a14a.v(ds1a,ds4a)*a23a.v(ds2a,ds3a)*s24s.v(ds2b,ds4b)*s143a.v(ds1b,ds3b)
        +a24a.v(ds2a,ds4a)*s14s.v(ds1a,ds4b)*s23s.v(ds2b,ds3a)*s341a.v(ds3b,ds1b)
        +a23a.v(ds2a,ds3a)*s14s.v(ds1a,ds4a)*s24s.v(ds2b,ds4b)*s341a.v(ds3b,ds1b))*CW*CW*CW*CW
        )/pDenU;


      
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble ZZWW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over spins 1/3^2=1/9
    return amp2/9.0;
  }
  



  //  Tests
  int test_ZZWW(){
    int n=0;//Number of fails
    std::cout<<"\t* Z , Z  -> W+, W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=250\n";
      ldouble mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      ZZWW ZZWWAmp = ZZWW(EE,mh,0,MW,0,SW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.218784675064574E+02,3.035340606909353E+01,1.347840875626388E+01,7.654802496655186E+00,5.018015978860242E+00,3.634216607622210E+00,2.846024916348151E+00,2.383362924509608E+00,2.122184053542665E+00,2.003628011521073E+00,2.003628011520370E+00,2.122184053542223E+00,2.383362924509694E+00,2.846024916348093E+00,3.634216607622649E+00,5.018015978860762E+00,7.654802496655284E+00,1.347840875626420E+01,3.035340606909304E+01,1.218784675064576E+02};
      i += ZZWWAmp.test_2to2_amp2([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH);
      i += ZZWWAmp.test_2to2_amp2_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH);
      i += ZZWWAmp.test_2to2_amp2_boosts([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH);
      i += ZZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH2[20] = {2.897427960518776E+01,1.437235517314631E+01,8.543721261435447E+00,5.684265799442969E+00,4.100104442565673E+00,3.154892264080381E+00,2.568401066136697E+00,2.203747929416787E+00,1.990103424164924E+00,1.891073747438559E+00,1.891073747438550E+00,1.990103424164917E+00,2.203747929416789E+00,2.568401066136700E+00,3.154892264080387E+00,4.100104442565669E+00,5.684265799442965E+00,8.543721261435456E+00,1.437235517314630E+01,2.897427960518773E+01};
      i += ZZWWAmp.test_2to2_amp2([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH2);
      i += ZZWWAmp.test_2to2_amp2_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH2);
      i += ZZWWAmp.test_2to2_amp2_boosts([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH2);
      i += ZZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH2);
      //std::cout<<"\n# mh=125, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH4[20] = {5.299651142764252E+00,5.299277865346616E+00,5.298946084430754E+00,5.298655792521403E+00,5.298406983060965E+00,5.298199650429202E+00,5.298033789942983E+00,5.297909397856059E+00,5.297826471358894E+00,5.297785008578538E+00,5.297785008578538E+00,5.297826471358893E+00,5.297909397856059E+00,5.298033789942983E+00,5.298199650429202E+00,5.298406983060965E+00,5.298655792521402E+00,5.298946084430754E+00,5.299277865346616E+00,5.299651142764252E+00};
      i += ZZWWAmp.test_2to2_amp2([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH4);
      i += ZZWWAmp.test_2to2_amp2_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH4);
      i += ZZWWAmp.test_2to2_amp2_boosts([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH4);
      i += ZZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH4);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=250\n";
      mh=1;
      ZZWWAmp.set_masses(mh,MW);
      pspatial=250;
      ldouble dataCH5[20] = {1.221350687402042E+02,3.047757047129016E+01,1.355867225887735E+01,7.713545181063446E+00,5.064290916243988E+00,3.672635298320071E+00,2.879310415187293E+00,2.413315692057453E+00,2.150120306798409E+00,2.030610893868498E+00,2.030610893867794E+00,2.150120306797949E+00,2.413315692057539E+00,2.879310415187236E+00,3.672635298320510E+00,5.064290916244522E+00,7.713545181063554E+00,1.355867225887767E+01,3.047757047128967E+01,1.221350687402043E+02};
      i += ZZWWAmp.test_2to2_amp2([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH5);
      i += ZZWWAmp.test_2to2_amp2_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH5);
      i += ZZWWAmp.test_2to2_amp2_boosts([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH5);
      i += ZZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH5);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=125\n";
      pspatial = 125;
      ldouble dataCH6[20] = {2.911022131522191E+01,1.446513995183542E+01,8.613130909593968E+00,5.739295560711144E+00,4.145643790136731E+00,3.193937136850835E+00,2.602964086350918E+00,2.235292982022554E+00,2.019779160450637E+00,1.919854033476574E+00,1.919854033476566E+00,2.019779160450629E+00,2.235292982022556E+00,2.602964086350922E+00,3.193937136850840E+00,4.145643790136726E+00,5.739295560711140E+00,8.613130909593977E+00,1.446513995183541E+01,2.911022131522189E+01};
      i += ZZWWAmp.test_2to2_amp2([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH6);
      i += ZZWWAmp.test_2to2_amp2_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH6);
      i += ZZWWAmp.test_2to2_amp2_boosts([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH6);
      i += ZZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH6);
      //std::cout<<"\n# mh=1, MW=80.385, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH7[20] = {5.400456960104258E+00,5.400079783557148E+00,5.399744536875505E+00,5.399451212519208E+00,5.399199803891410E+00,5.398990305338236E+00,5.398822712148521E+00,5.398697020553595E+00,5.398613227727105E+00,5.398571331784889E+00,5.398571331784889E+00,5.398613227727104E+00,5.398697020553595E+00,5.398822712148521E+00,5.398990305338236E+00,5.399199803891411E+00,5.399451212519208E+00,5.399744536875505E+00,5.400079783557148E+00,5.400456960104258E+00};
      i += ZZWWAmp.test_2to2_amp2([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH7);
      i += ZZWWAmp.test_2to2_amp2_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH7);
      i += ZZWWAmp.test_2to2_amp2_boosts([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH7);
      i += ZZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH7);
      //std::cout<<"\n# mh=125, MW=50, pspatial=125\n";
      mh=125;
      MW=50;
      MZ=MW/CW;
      ZZWWAmp.set_masses(mh,MW);
      pspatial = 125;
      ldouble dataCH9[20] = {8.037053939959351E+01,2.492965811824679E+01,1.195271913619398E+01,7.049291955728941E+00,4.722182876960156E+00,3.466434569025781E+00,2.738545103648279E+00,2.306462970277772E+00,2.060841706754497E+00,1.948921210776397E+00,1.948921210776458E+00,2.060841706754377E+00,2.306462970277771E+00,2.738545103648250E+00,3.466434569025708E+00,4.722182876960120E+00,7.049291955728955E+00,1.195271913619403E+01,2.492965811824680E+01,8.037053939959331E+01};
      i += ZZWWAmp.test_2to2_amp2([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH9);
      i += ZZWWAmp.test_2to2_amp2_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH9);
      i += ZZWWAmp.test_2to2_amp2_boosts([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH9);
      i += ZZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH9);
      //std::cout<<"\n# mh=125, MW=50, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH10[20] = {6.735644317458705E+00,6.734603771230734E+00,6.733678989509639E+00,6.732869919947022E+00,6.732176516751526E+00,6.731598740683431E+00,6.731136559050035E+00,6.730789945701794E+00,6.730558881029247E+00,6.730443351960704E+00,6.730443351960704E+00,6.730558881029247E+00,6.730789945701794E+00,6.731136559050035E+00,6.731598740683431E+00,6.732176516751526E+00,6.732869919947022E+00,6.733678989509638E+00,6.734603771230734E+00,6.735644317458705E+00};
      i += ZZWWAmp.test_2to2_amp2([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH10);
      i += ZZWWAmp.test_2to2_amp2_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH10);
      i += ZZWWAmp.test_2to2_amp2_boosts([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH10);
      i += ZZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH10);
      //std::cout<<"\n# mh=1, MW=50, pspatial=125\n";
      mh=1;
      MW=50;
      MZ=MW/CW;
      ZZWWAmp.set_masses(mh,MW);
      pspatial = 125;
      ldouble dataCH11[20] = {8.095572928488492E+01,2.522635050632390E+01,1.213814497633243E+01,7.176797932512432E+00,4.815042184981791E+00,3.537034173768330E+00,2.794417325328955E+00,2.352694702877831E+00,2.101209928688888E+00,1.986509243437368E+00,1.986509243437429E+00,2.101209928688769E+00,2.352694702877820E+00,2.794417325328927E+00,3.537034173768262E+00,4.815042184981752E+00,7.176797932512439E+00,1.213814497633249E+01,2.522635050632392E+01,8.095572928488473E+01};
      i += ZZWWAmp.test_2to2_amp2([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH11);
      i += ZZWWAmp.test_2to2_amp2_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH11);
      i += ZZWWAmp.test_2to2_amp2_boosts([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH11);
      i += ZZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH11);
      //std::cout<<"\n# mh=1, MW=50, pspatial=1\n";
      pspatial = 1;
      ldouble dataCH12[20] = {5.400154835027087E+00,5.399179883366282E+00,5.398313402400197E+00,5.397555341731683E+00,5.396905657276364E+00,5.396364311257388E+00,5.395931272200915E+00,5.395606514932374E+00,5.395390020573461E+00,5.395281776539889E+00,5.395281776539889E+00,5.395390020573461E+00,5.395606514932374E+00,5.395931272200915E+00,5.396364311257388E+00,5.396905657276365E+00,5.397555341731683E+00,5.398313402400197E+00,5.399179883366282E+00,5.400154835027087E+00};
      i += ZZWWAmp.test_2to2_amp2([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH12);
      i += ZZWWAmp.test_2to2_amp2_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH12);
      i += ZZWWAmp.test_2to2_amp2_boosts([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH12);
      i += ZZWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ZZWWAmp.amp2(); }, MZ,MZ,MW,MW,pspatial,dataCH12);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
  

}
