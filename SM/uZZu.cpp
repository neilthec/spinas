
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

//File:  SPINAS/SM/uZZu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uZZu.h"

namespace spinas {

  uZZu::uZZu(const ldouble& echarge, const ldouble& massf, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), Qf(2.0/3.0), mf(massf), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propf(massf,0) {
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
    gL=-2.0*Qf*SW*SW+1.0;
    gR=-2.0*Qf*SW*SW;
    preh = e*e*mf/(2.0*MW*MW*SW*SW);
  }
  void uZZu::set_masses(const ldouble& massf, const ldouble& massh, const ldouble& massW){
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
  void uZZu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
  cdouble uZZu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
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
      //uUZZ all ingoing: 
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //uZZU: 2<->4
      //preh [23] <23> ([14]+<14>)/(u-Mh^2)
      //34 out:
      //- preh [23] <23> ([14]-<14>)/(u-Mh^2)
      amplitude += - normFactor*preh*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(s14s.v(ds1,ds4)-a14a.v(ds1,ds4))/pDenU;
      
      //T-Channel e
      //preTU = e*e/(4.0*MW*MW*SW*SW);
      //uUZZ all ingoing:
      //+preTU gLe^2 [24] <13> (MZ <34>+[314>))/(Mf^2-t)
      //+preTU gRe^2 [13] <24> (MZ [34]+[413>))/(Mf^2-t)
      //+preTU gLe gRe Mf ([13] [24] <34>+[34] <13> <24>))/(Mf^2-t)
      //uZZU: 2<->4
      //+preTU gLe^2 [24] <13> (-MZ <23>+[312>))/(t-Mf^2)
      //+preTU gRe^2 [13] <24> (-MZ [23]+[213>))/(t-Mf^2)
      //-preTU gLe gRe Mf ([13] [24] <23>+[23] <13> <24>))/(t-Mf^2)
      //34 out:
      //-preTU gLe^2 [24] <13> (MZ <23>+[312>))/(t-Mf^2)
      //+preTU gRe^2 [13] <24> (MZ [23]+[213>))/(t-Mf^2)
      //+preTU gLe gRe Mf ([13] [24] <23>-[23] <13> <24>))/(t-Mf^2)
      amplitude +=
	-          normFactor*preTU*gL*gL*s24s.v(ds2a,ds4)*a13a.v(ds1,ds3a)*(MZ*a23a.v(ds2b,ds3b)+s312a.v(ds3b,ds2b))/pDenT
      	+          normFactor*preTU*gR*gR*a24a.v(ds2a,ds4)*s13s.v(ds1,ds3a)*(MZ*s23s.v(ds2b,ds3b)+s213a.v(ds2b,ds3b))/pDenT
      	+          normFactor*preTU*gL*gR*mf*(s13s.v(ds1,ds3a)*s24s.v(ds2a,ds4)*a23a.v(ds2b,ds3b)-a13a.v(ds1,ds3a)*a24a.v(ds2a,ds4)*s23s.v(ds2b,ds3b))/pDenT;

      //U-Channel e
      //preTU = e*e/(4.0*MW*MW*SW*SW);
      //uUZZ all ingoing:
      //+preTU gLe^2 [23] <14> (MZ <34>-[413>)/(u-Mf^2)
      //+preTU gRe^2 [14] <23> (MZ [34]-[314>)/(u-Mf^2)
      //+preTU gLe gRe Mf ([14] [23] <34>+[34] <14> <23>)/(u-Mf^2)
      //uZZU: 2<->4
      //+preTU gLe^2 [34] <12> (MZ <23>+[213>)/(s-Mf^2)
      //+preTU gRe^2 [12] <34> (MZ [23]+[312>)/(s-Mf^2)
      //+preTU gLe gRe Mf ([12] [34] <23>+[23] <12> <34>)/(s-Mf^2)
      //34 out:/
      //-preTU gLe^2 [34] <12> (MZ <23>+[213>)/(s-Mf^2)
      //+preTU gRe^2 [12] <34> (MZ [23]+[312>)/(s-Mf^2)
      //+preTU gLe gRe Mf (-[12] [34] <23>+[23] <12> <34>)/(s-Mf^2)
      amplitude +=
	-           normFactor*preTU*gL*gL*s34s.v(ds3a,ds4)*a12a.v(ds1,ds2a)*(MZ*a23a.v(ds2b,ds3b)+s213a.v(ds2b,ds3b))/pDenS
      	+           normFactor*preTU*gR*gR*s12s.v(ds1,ds2a)*a34a.v(ds3a,ds4)*(MZ*s23s.v(ds2b,ds3b)+s312a.v(ds3b,ds2b))/pDenS
      	+           normFactor*preTU*gL*gR*mf*(-s12s.v(ds1,ds2a)*s34s.v(ds3a,ds4)*a23a.v(ds2b,ds3b)+a12a.v(ds1,ds2a)*a34a.v(ds3a,ds4)*s23s.v(ds2b,ds3b))/pDenS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uZZu::amp2(){
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
  int test_uZZu(){
    int n=0;//Number of fails
    std::cout<<"\t* u , Z  -> Z , u       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mf=0.0042, mh=125, MW=80.385, pspatial=250\n";
      ldouble mf=0.0042, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      uZZu uZZuAmp = uZZu(EE,mf,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.361218317809021E-01,4.318339326177314E-02,2.601782372444087E-02,1.887469586162279E-02,1.500696161625270E-02,1.261367057507442E-02,1.100966996281555E-02,9.877375910975025E-03,9.049434213333573E-03,8.429284276763873E-03,7.957305836206456E-03,7.594696104937499E-03,7.315113269897005E-03,7.100046163216358E-03,6.936105782369646E-03,6.813368250443147E-03,6.724321954536023E-03,6.663176579439148E-03,6.625396774719249E-03,6.607379629665790E-03};
      i += uZZuAmp.test_2to2_amp2([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH);
      i += uZZuAmp.test_2to2_amp2_rotations([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH);
      i += uZZuAmp.test_2to2_amp2_boosts([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH);
      i += uZZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH);
      //std::cout<<"\n# mf=0.0042, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {7.505265966559196E-01,6.242160109411553E-02,3.121330979491984E-02,2.116178555748002E-02,1.630801149932568E-02,1.349187662396934E-02,1.167890143340509E-02,1.043241979714247E-02,9.536741916125836E-03,8.873213059161875E-03,8.371212241730278E-03,7.986099613360579E-03,7.688283559604103E-03,7.457376675616048E-03,7.278877444468423E-03,7.142190257976118E-03,7.039393962507049E-03,6.964448483016264E-03,6.912667949376896E-03,6.880361454826674E-03};
      i += uZZuAmp.test_2to2_amp2([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH2);
      i += uZZuAmp.test_2to2_amp2_rotations([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH2);
      i += uZZuAmp.test_2to2_amp2_boosts([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH2);
      i += uZZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH2);
      //std::cout<<"\n# mf=125.1, mh=125, MW=80.385, pspatial=95\n";
      mf = 125.1;
      mh = 125;
      pspatial = 95;
      uZZuAmp.set_masses(mf,mh,MW);
      ldouble dataCH4[20] = {2.798028969980267E-01,2.627771461529683E-01,2.495564743773021E-01,2.392083720174515E-01,2.311093784120436E-01,2.248339110612301E-01,2.200898314145681E-01,2.166800408025495E-01,2.144795462897187E-01,2.134224549075689E-01,2.134960505103904E-01,2.147407377325697E-01,2.172558380162700E-01,2.212123967113974E-01,2.268756858339422E-01,2.346424982528375E-01,2.451025949284506E-01,2.591417754361685E-01,2.781203932388273E-01,3.041961483541079E-01};
      i += uZZuAmp.test_2to2_amp2([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH4);
      i += uZZuAmp.test_2to2_amp2_rotations([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH4);
      i += uZZuAmp.test_2to2_amp2_boosts([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH4);
      i += uZZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH4);
      //std::cout<<"\n# mf=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      mf = 125;
      mh = 0.0005;
      pspatial = 125.1;
      uZZuAmp.set_masses(mf,mh,MW);
      ldouble dataCH3[20] = {4.330581810832074E-01,3.940470871296521E-01,3.683342768970890E-01,3.515455932839719E-01,3.413110861395466E-01,3.363357057417961E-01,3.359847322995867E-01,3.401032008759065E-01,3.489666648773297E-01,3.633348417532424E-01,3.846271420149551E-01,4.152975078370743E-01,4.596026062750684E-01,5.252494152552792E-01,6.272536409107764E-01,7.981943811761409E-01,1.120835140695923E+00,1.864450196361229E+00,4.406392181807410E+00,3.400901240112917E+01};
      i += uZZuAmp.test_2to2_amp2([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH3);
      i += uZZuAmp.test_2to2_amp2_rotations([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH3);
      i += uZZuAmp.test_2to2_amp2_boosts([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH3);
      i += uZZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uZZuAmp.amp2(); }, mf,MZ,MZ,mf,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
