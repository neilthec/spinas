
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

//File:  SPINAS/SM/dZdh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dZdh.h"

namespace spinas {

  dZdh::dZdh(const ldouble& echarge, const ldouble& massd, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), md(massd), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propd(massd,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);  
    p1=particle(md);
    p2=particle(MZ);
    p3=particle(md);
    p4=particle(mh);
    //<12>,[12],<23>,[23],<13>,[13]
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    //[312>,[213>,[242>,[142>,[241>
    s312a = sproduct(SQUARE,&p3,&p1,&p2,2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3,2);
    s242a = sproduct(SQUARE,&p2,&p4,&p2,2);
    s142a = sproduct(SQUARE,&p1,&p4,&p2,2);
    s241a = sproduct(SQUARE,&p2,&p4,&p1,2);
    //Couplings
    preTU = sqrt2*e*e*md/(4.0*MW*MW*SW*SW);
    gL=2.0/3.0*SW*SW-1.0;
    gR=2.0/3.0*SW*SW;
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void dZdh::set_masses(const ldouble& massd, const ldouble& massh, const ldouble& massW){
    md=massd;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(md);
    p2.set_mass(MZ);
    p3.set_mass(md);
    p4.set_mass(mh);
    propd.set_mass(md);
    propZ.set_mass(MZ);
    //Couplings
    preTU = sqrt2*e*e*md/(4.0*MW*MW*SW*SW);
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void dZdh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    //[312>,[213>,[242>,[142>,[241>
    s312a.update();
    s213a.update();
    s242a.update();
    s142a.update();
    s241a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propd.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenT=propZ.denominator(propTP);
    pDenU=propd.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble dZdh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    int ds2a, ds2b;

    //Symmetrize the Z-Boson Spin indices
    int nCombs=get_num_spin_loops(ds2);
    ldouble normFactor=get_spin_normalization(ds2);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, i);
      
      //S-Channel Z
      //preZ = e*e/(4.0*MW*MW*SW*SW);
      //dDZh all in: 
      //preZ (Md (gLe-gRe) (<12>-[12]) [343> + 2 MZ^2(gLe[23] <13> + gRe[13] <23>)))/(s-MZ^2)
      //dZDh: 2<->3: 
      //preZ (Md (gLe-gRe) (<13>-[13]) [242> - 2 MZ^2(gLe[23] <12> + gRe[12] <23>)))/(t-MZ^2)
      //34 out:
      //preZ (Md (gLe-gRe) (<13>+[13]) [242> - 2 MZ^2(gLe[23] <12> - gRe[12] <23>)))/(t-MZ^2)
      amplitude += normFactor*preZ*(md*(gL-gR)*(a13a.v(ds1,ds3)+s13s.v(ds1,ds3))*s242a.v(ds2a,ds2b) - 2.0*MZ*MZ*(gL*s23s.v(ds2a,ds3)*a12a.v(ds1,ds2b) - gR*a23a.v(ds2a,ds3)*s12s.v(ds1,ds2b)))/pDenT;
      
      //T-Channel e
      //preTU = e*e*md/(4.0*MW*MW*SW*SW);
      //dDZh all in:
      //preh*(gLe <13> (Md [23]-[312>+MZ <23>)+gRe [13] (MZ [23]-[213>+Md <23>)))/(t-Md^2)
      //dZDh: 2<->3:
      //-preh*(gLe <12> (Md [23]+[213>+MZ <23>)+gRe [12] (MZ [23]+[312>+Md <23>)))/(s-Md^2)
      //34 out:
      //-preh*(gLe <12> (Md [23]-[213>-MZ <23>)+gRe [12] (MZ [23]+[312>-Md <23>)))/(s-Md^2)
      amplitude += -normFactor*preTU*( gL*a12a.v(ds1,ds2a)*(md*s23s.v(ds2b,ds3)-s213a.v(ds2b,ds3)-MZ*a23a.v(ds2b,ds3)) + gR*s12s.v(ds1,ds2a)*(-md*a23a.v(ds2b,ds3)+s312a.v(ds3,ds2b)+MZ*s23s.v(ds2b,ds3)) )/pDenS;

      //U-Channel e
      //preTU = e*e*md/(4.0*MW*MW*SW*SW);
      //dDZh all in:
      //preh (gLe [23] ([143>+2 Md <13>)+gRe <23> (2 Md [13]+[341>) )/(u-Md^2)
      //dZDh: 2<->3:
      //- preh (gLe [23] ([142>+2 Md <12>)+gRe <23> (2 Md [12]+[241>) )/(u-Md^2)
      //34 out:
      //- preh (gLe [23] (-[142>+2 Md <12>)-gRe <23> (2 Md [12]-[241>) )/(u-Md^2)
      amplitude += -normFactor*preTU*( gL*s23s.v(ds2a,ds3)*(-s142a.v(ds1,ds2b)+2.0*md*a12a.v(ds1,ds2b)) - gR*a23a.v(ds2a,ds3)*(-s241a.v(ds2b,ds1)+2.0*md*s12s.v(ds1,ds2b)) )/pDenU;
    }

    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble dZdh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-1;j3<=1;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2*1/3=1/6
    //Average over initial colors 1/3
    return amp2/18.0;
  }
  



  //  Tests
  int test_dZdh(){
    int n=0;//Number of fails
    std::cout<<"\t* d , Z  -> d , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# md=0.0075, mh=125, MW=80.385, pspatial=250\n";
      ldouble md=0.0075, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      dZdh dZdhAmp = dZdh(EE,md,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.003327946296211E+00,8.416894528240376E-01,3.762784774205344E-01,2.064038841199361E-01,1.269392911045942E-01,8.391688539883109E-02,5.825945681761747E-02,4.186880349807850E-02,3.084600886537151E-02,2.313224027335849E-02,1.756048786472848E-02,1.343103812487691E-02,1.030476313699596E-02,7.895442489438379E-03,6.010472128173930E-03,4.516663248163005E-03,3.319718318341954E-03,2.351488016772902E-03,1.561821120199909E-03,9.132149585316279E-04};
      i += dZdhAmp.test_2to2_amp2([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH);
      i += dZdhAmp.test_2to2_amp2_rotations([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH);
      i += dZdhAmp.test_2to2_amp2_boosts([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH);
      i += dZdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH);
      //std::cout<<"\n# md=0.0042, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {4.484180750001232E-01,2.594423197072280E-01,1.660322797705285E-01,1.135059939376067E-01,8.126932570483239E-02,6.019185225129742E-02,4.573382078040696E-02,3.543688405613662E-02,2.787873900484719E-02,2.219181634736965E-02,1.782352729276281E-02,1.440898210047854E-02,1.169976309312878E-02,9.522258605435846E-03,7.752348679971561E-03,6.299513576084354E-03,5.096568343320697E-03,4.092862571336806E-03,3.249674167484711E-03,2.537075819187221E-03};
      i += dZdhAmp.test_2to2_amp2([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH2);
      i += dZdhAmp.test_2to2_amp2_rotations([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH2);
      i += dZdhAmp.test_2to2_amp2_boosts([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH2);
      i += dZdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH2);
      //std::cout<<"\n# md=125.1, mh=125, MW=80.385, pspatial=95\n";
      md = 125.1;
      mh = 125;
      pspatial = 95;
      dZdhAmp.set_masses(md,mh,MW);
      ldouble dataCH4[20] = {3.944772955561968E-01,3.063043892806311E-01,2.477461333796074E-01,2.072748545964575E-01,1.784810114521640E-01,1.575818981002063E-01,1.422378707368183E-01,1.309477222723602E-01,1.227213242075518E-01,1.168938298003504E-01,1.130162477209700E-01,1.107893218722422E-01,1.100232169939400E-01,1.106134923075564E-01,1.125282013081167E-01,1.158035484083707E-01,1.205473253192786E-01,1.269509355178881E-01,1.353127086091445E-01,1.460780919840986E-01};
      i += dZdhAmp.test_2to2_amp2([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH4);
      i += dZdhAmp.test_2to2_amp2_rotations([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH4);
      i += dZdhAmp.test_2to2_amp2_boosts([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH4);
      i += dZdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH4);
      //std::cout<<"\n# md=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      md = 125;
      mh = 0.0005;
      pspatial = 125.1;
      dZdhAmp.set_masses(md,mh,MW);
      ldouble dataCH3[20] = {7.879419676024333E-01,4.309705171728301E-01,2.748693114970894E-01,1.941853129472658E-01,1.478742883324591E-01,1.193988966254154E-01,1.011081983724801E-01,8.910185571226936E-02,8.125158608905460E-02,7.634637730028920E-02,7.369142461327297E-02,7.291049647843524E-02,7.384971009952579E-02,7.654003471338860E-02,8.120403880595076E-02,8.831167718195147E-02,9.871409043923031E-02,1.139342862118837E-01,1.368267501121754E-01,1.732416675835024E-01};
      i += dZdhAmp.test_2to2_amp2([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH3);
      i += dZdhAmp.test_2to2_amp2_rotations([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH3);
      i += dZdhAmp.test_2to2_amp2_boosts([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH3);
      i += dZdhAmp.test_2to2_amp2_boosts_and_rotations([&]() { return dZdhAmp.amp2(); }, md,MZ,md,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
