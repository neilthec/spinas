
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

//File:  SPINAS/SM/eZeh.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eZeh.h"

namespace spinas {

  eZeh::eZeh(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), prope(masse,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);  
    p1=particle(me);
    p2=particle(MZ);
    p3=particle(me);
    p4=particle(mh);
    //<12>,[12],<23>,[23],<13>,[13]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    //[312>,[213>,[242>,[142>,[241>
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
    s242a = sproduct(SQUARE,&p2,&p4,&p2);
    s142a = sproduct(SQUARE,&p1,&p4,&p2);
    s241a = sproduct(SQUARE,&p2,&p4,&p1);
    //Couplings
    preSU = sqrt2*e*e*me/(4.0*MW*MW*SW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void eZeh::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(me);
    p2.set_mass(MZ);
    p3.set_mass(me);
    p4.set_mass(mh);
    prope.set_mass(me);
    propZ.set_mass(MZ);
    //Couplings
    preSU = sqrt2*e*e*me/(4.0*MW*MW*SW*SW);
    preZ = sqrt2*e*e/(4.0*MW*MW*SW*SW);
  }
  void eZeh::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenT=propZ.den(propTP);
    pDenS=prope.den(propSP);
    pDenU=prope.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eZeh::amp(const int& ds1, const int& ds2, const int& ds3){//Double Spin
    cdouble amplitude(0,0);
    int ds2a, ds2b;

    //Symmetrize the Z-Boson Spin indices
    int nCombs=get_num_spin_loops(ds2);
    ldouble normFactor=get_spin_normalization(ds2);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, i);
      
      //T-Channel Z
      //preZ = e*e/(4.0*MW*MW*SW*SW);
      //eEZh all ingoing: 
      //+ preZ ( Me (gLe-gRe) (<12>-[12]) [343> + 2 MZ^2(gLe[23] <13> + gRe[13] <23>) )/(s-MZ^2)
      //eZEh: 2<->3
      //+ preZ ( Me (gLe-gRe) (<13>-[13]) [242> - 2 MZ^2(gLe[23] <12> + gRe[12] <23>) )/(t-MZ^2)
      //34 out:
      //+ preZ ( Me (gLe-gRe) (<13>+[13]) [242> - 2 MZ^2(gLe[23] <12> - gRe[12] <23>) )/(t-MZ^2)
      amplitude += normFactor*preZ*(
				    + me*(gL-gR)*(a13a.v(ds1,ds3)+s13s.v(ds1,ds3))*s242a.v(ds2a,ds2b)
				    - 2.0*MZ*MZ*(gL*s23s.v(ds2a,ds3)*a12a.v(ds1,ds2b) - gR*s12s.v(ds1,ds2a)*a23a.v(ds2b,ds3))
				    )/pDenT;
      
      //S-Channel e
      //preSU = e*e*me/(4.0*MW*MW*SW*SW);
      //eEZh all ingoing:
      //preh*(gLe <13> (Me [23]-[312>+MZ <23>)+gRe [13] (MZ [23]-[213>+Me <23>)))/(t-Me^2)
      //eZEh: 2<->3
      //- preh*(gLe <12> (Me [23]+[213>+MZ <23>)+gRe [12] (MZ [23]+[312>+Me <23>)))/(s-Me^2)
      //34 out:
      //- preh*(gLe <12> (Me [23]-[213>-MZ <23>)+gRe [12] (MZ [23]+[312>-Me <23>)))/(s-Me^2)
      amplitude += -normFactor*preSU*(
				      + gL*a12a.v(ds1,ds2a)*(me*s23s.v(ds2b,ds3)-s213a.v(ds2b,ds3)-MZ*a23a.v(ds2b,ds3))
				      + gR*s12s.v(ds1,ds2a)*(MZ*s23s.v(ds2b,ds3)+s312a.v(ds3,ds2b)-me*a23a.v(ds2b,ds3))
				      )/pDenS;
      
      //U-Channel e
      //preSU = e*e*me/(4.0*MW*MW*SW*SW);
      //eEZh all ingoing:
      //+ preh (gLe [23] ([143>+2 Me <13>)+gRe <23> (2 Me [13]+[341>) )/(u-Me^2)
      //eZEh: 2<->3
      //- preh (gLe [23] ([142>+2 Me <12>)+gRe <23> (2 Me [12]+[241>) )/(u-Me^2)
      //34 out:
      //- preh (gLe [23] (-[142>+2 Me <12>)-gRe <23> (2 Me [12]-[241>) )/(u-Me^2)
      amplitude += - normFactor*preSU*(
				       + gL*s23s.v(ds2a,ds3)*(-s142a.v(ds1,ds2b)+2.0*me*a12a.v(ds1,ds2b))
				       - gR*a23a.v(ds2a,ds3)*(-s241a.v(ds2b,ds1)+2.0*me*s12s.v(ds1,ds2b))
				       )/pDenU;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eZeh::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-1;j3<=1;j3+=2){
	  M = amp(j1,j2,j3);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2*1/3=1/6
    return amp2/6.0;
  }
  



  //  Tests
  int test_eZeh(){
    int n=0;//Number of fails
    std::cout<<"\t* e , Z  -> e , h       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      eZeh eZehAmp = eZeh(EE,me,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.035512264110678E+00,5.704569178494980E-01,2.550235834686065E-01,1.398906961420196E-01,8.603338966404971E-02,5.687485752219632E-02,3.948547764611723E-02,2.837667552612985E-02,2.090595169247571E-02,1.567792769278886E-02,1.190165997373524E-02,9.102916128186197E-03,6.984076223299059E-03,5.351153731143494E-03,4.073610819946204E-03,3.061178475213540E-03,2.249946327745205E-03,1.593726006532460E-03,1.058527350403071E-03,6.189322476925001E-04};
      i += eZehAmp.test_2to2_amp2([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH);
      i += eZehAmp.test_2to2_amp2_rotations([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH);
      i += eZehAmp.test_2to2_amp2_boosts([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH);
      i += eZehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=35\n";
      pspatial = 35;
      ldouble dataCH2[20] = {1.778226578610318E-03,1.751187583202719E-03,1.724682735195740E-03,1.698698463595645E-03,1.673221620922363E-03,1.648239467533086E-03,1.623739656616632E-03,1.599710219826065E-03,1.576139553518763E-03,1.553016405574835E-03,1.530329862766312E-03,1.508069338651021E-03,1.486224561966428E-03,1.464785565500037E-03,1.443742675414150E-03,1.423086501003970E-03,1.402807924869073E-03,1.382898093479355E-03,1.363348408117479E-03,1.344150516180798E-03};
      i += eZehAmp.test_2to2_amp2([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH2);
      i += eZehAmp.test_2to2_amp2_rotations([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH2);
      i += eZehAmp.test_2to2_amp2_boosts([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH2);
      i += eZehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, mh=125, MW=80.385, pspatial=65\n";
      me = 125.1;
      mh = 125;
      pspatial = 65;
      eZehAmp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {1.551236140094245E-01,1.511544130707358E-01,1.475713141049073E-01,1.443425784069567E-01,1.414402363478566E-01,1.388396123197874E-01,1.365189211698484E-01,1.344589244346613E-01,1.326426368240397E-01,1.310550751170403E-01,1.296830430180666E-01,1.285149466436164E-01,1.275406362256420E-01,1.267512703674330E-01,1.261391998054927E-01,1.256978681421845E-01,1.254217274397159E-01,1.253061669228558E-01,1.253474533389861E-01,1.255426817804170E-01};
      i += eZehAmp.test_2to2_amp2([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH4);
      i += eZehAmp.test_2to2_amp2_rotations([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH4);
      i += eZehAmp.test_2to2_amp2_boosts([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH4);
      i += eZehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH4);
      //std::cout<<"\n# me=125, mh=0.0005, MW=80.385, pspatial=1\n";
      me = 125;
      mh = 0.0005;
      pspatial = 1;
      eZehAmp.set_masses(me,mh,MW);
      ldouble dataCH3[20] = {8.078324914831650E-02,8.074087324000968E-02,8.069902919537844E-02,8.065771581143902E-02,8.061693189821847E-02,8.057667627869479E-02,8.053694778873793E-02,8.049774527705110E-02,8.045906760511308E-02,8.042091364712083E-02,8.038328228993309E-02,8.034617243301417E-02,8.030958298837873E-02,8.027351288053687E-02,8.023796104643997E-02,8.020292643542744E-02,8.016840800917319E-02,8.013440474163387E-02,8.010091561899654E-02,8.006793963962794E-02};
      i += eZehAmp.test_2to2_amp2([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH3);
      i += eZehAmp.test_2to2_amp2_rotations([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH3);
      i += eZehAmp.test_2to2_amp2_boosts([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH3);
      i += eZehAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eZehAmp.amp2(); }, me,MZ,me,mh,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
