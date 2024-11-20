
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

//File:  SPINAS/SM/eZZe.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eZZe.h"

namespace spinas {

  eZZe::eZZe(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), prope(masse,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    proph = propagator(mh,wh);  
    p1=particle(me);
    p2=particle(MZ);
    p3=particle(MZ);
    p4=particle(me);
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
    s312a = sproduct(SQUARE,&p3,&p1,&p2);
    s213a = sproduct(SQUARE,&p2,&p1,&p3);
    //Couplings
    preTS = e*e/(2.0*MW*MW*SW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
    preh = e*e*me/(2.0*MW*MW*SW*SW);
  }
  void eZZe::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(me);
    p2.set_mass(MZ);
    p3.set_mass(MZ);
    p4.set_mass(me);
    prope.set_mass(me);
    proph.set_mass(mh);
    //Couplings
    preTS = e*e/(2.0*MW*MW*SW*SW);
    preh = e*e*me/(2.0*MW*MW*SW*SW);
  }
  void eZZe::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    s312a.update();
    s213a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenU=proph.denominator(propUP);
    pDenT=prope.denominator(propTP);
    pDenS=prope.denominator(propSP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eZZe::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds2a, ds2b;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds2,ds3);
    ldouble normFactor=get_spin_normalization(ds2,ds3);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, ds3,ds3a,ds3b, i);
      
      //U-Channel h
      //preh = e*e*me/(2.0*MW*MW*SW*SW);
      //eEZZ all ingoing: 
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //eZZE: 2<->4
      //preh [23] <23> ([14]+<14>)/(u-Mh^2)
      //34 out:
      //- preh [23] <23> ([14]-<14>)/(u-Mh^2)
      amplitude += - normFactor*preh*s23s.v(ds2a,ds3a)*a23a.v(ds2b,ds3b)*(s14s.v(ds1,ds4)-a14a.v(ds1,ds4))/pDenU;
      
      //T-Channel e
      //preTS = e*e/(2.0*MW*MW*SW*SW);
      //eEZZ all ingoing:
      //+preTS gLe^2 [24] <13> (MZ <34>+[314>))/(Me^2-t)
      //+preTS gRe^2 [13] <24> (MZ [34]+[413>))/(Me^2-t)
      //+preTS gLe gRe Me ([13] [24] <34>+[34] <13> <24>))/(Me^2-t)
      //eZZE: 2<->4
      //+preTS gLe^2 [24] <13> (- MZ <23>+[312>))/(t-Me^2)
      //+preTS gRe^2 [13] <24> (- MZ [23]+[213>))/(t-Me^2)
      //-preTS gLe gRe Me ([13] [24] <23>+[23] <13> <24>))/(t-Me^2)
      //34 out:
      //-preTS gLe^2 [24] <13> (MZ <23>+[312>))/(t-Me^2)
      //+preTS gRe^2 [13] <24> (MZ [23]+[213>))/(t-Me^2)
      //+preTS gLe gRe Me ([13] [24] <23>-[23] <13> <24>))/(t-Me^)
      amplitude += - normFactor*preTS*gL*gL*s24s.v(ds2a,ds4)*a13a.v(ds1,ds3a)*(MZ*a23a.v(ds2b,ds3b)+s312a.v(ds3b,ds2b))/pDenT
      	+            normFactor*preTS*gR*gR*a24a.v(ds2a,ds4)*s13s.v(ds1,ds3a)*(MZ*s23s.v(ds2b,ds3b)+s213a.v(ds2b,ds3b))/pDenT
      	+            normFactor*preTS*gL*gR*me*(s13s.v(ds1,ds3a)*s24s.v(ds2a,ds4)*a23a.v(ds2b,ds3b)-a13a.v(ds1,ds3a)*a24a.v(ds2a,ds4)*s23s.v(ds2b,ds3b))/pDenT;

      //S-Channel e
      //preTS = e*e/(2.0*MW*MW*SW*SW);
      //eEZZ all ingoing:
      //+preTS gLe^2 [23] <14> (MZ <34>-[413>)/(u-Me^2)
      //+preTS gRe^2 [14] <23> (MZ [34]-[314>)/(u-Me^2)
      //+preTS gLe gRe Me ([14] [23] <34>+[34] <14> <23>)/(u-Me^2)
      //eZZE: 2<->4
      //-preTS gLe^2 [34] <12> (-MZ <23>-[213>)/(s-Me^2)
      //-preTS gRe^2 [12] <34> (-MZ [23]-[312>)/(s-Me^2)
      //+preTS gLe gRe Me ([12] [34] <23>+[23] <12> <34>)/(s-Me^2)
      //34 out:
      //-preTS gLe^2 [34] <12> (MZ <23>+[213>)/(s-Me^2)
      //+preTS gRe^2 [12] <34> (MZ [23]+[312>)/(s-Me^2)
      //+preTS gLe gRe Me (-[12] [34] <23>+[23] <12> <34>)/(s-Me^2)
      amplitude += -normFactor*preTS*gL*gL*s34s.v(ds3a,ds4)*a12a.v(ds1,ds2a)*(MZ*a23a.v(ds2b,ds3b)+s213a.v(ds2b,ds3b))/pDenS
	+           normFactor*preTS*gR*gR*s12s.v(ds1,ds2a)*a34a.v(ds3a,ds4)*(MZ*s23s.v(ds2b,ds3b)+s312a.v(ds3b,ds2b))/pDenS
      	+           normFactor*preTS*gL*gR*me*(-s12s.v(ds1,ds2a)*s34s.v(ds3a,ds4)*a23a.v(ds2b,ds3b)+a12a.v(ds1,ds2a)*a34a.v(ds3a,ds4)*s23s.v(ds2b,ds3b))/pDenS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eZZe::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/3=1/6
    return amp2/6.0;
  }
  



  //  Tests
  int test_eZZe(){
    int n=0;//Number of fails
    std::cout<<"\t* e , Z  -> Z , e       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      eZZe eZZeAmp = eZZe(EE,me,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {7.262234196914057E-02,2.303876680301348E-02,1.388076590524249E-02,1.006983663974295E-02,8.006362220369032E-03,6.729517813090797E-03,5.873767645086135E-03,5.269677585826917E-03,4.827962513415049E-03,4.497106397746591E-03,4.245301238587764E-03,4.051845366002976E-03,3.902685162581659E-03,3.787944737469834E-03,3.700480925628678E-03,3.634999240792869E-03,3.587492163455797E-03,3.554870503660180E-03,3.534714608343737E-03,3.525102294305462E-03};
      i += eZZeAmp.test_2to2_amp2([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH);
      i += eZZeAmp.test_2to2_amp2_rotations([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH);
      i += eZZeAmp.test_2to2_amp2_boosts([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH);
      i += eZZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {4.004133793743134E-01,3.330254106787883E-02,1.665260913448643E-02,1.129002162448603E-02,8.700485207996217E-03,7.198049430833038E-03,6.230809257222156E-03,5.565798991629806E-03,5.087945995125090E-03,4.733946796331096E-03,4.466124402333261E-03,4.260663011937129E-03,4.101775206548599E-03,3.978584104345683E-03,3.883353002884104E-03,3.810429040664828E-03,3.755586204191805E-03,3.715602054122860E-03,3.687976633105390E-03,3.670740807196489E-03};
      i += eZZeAmp.test_2to2_amp2([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH2);
      i += eZZeAmp.test_2to2_amp2_rotations([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH2);
      i += eZZeAmp.test_2to2_amp2_boosts([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH2);
      i += eZZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, mh=125, MW=80.385, pspatial=95\n";
      me = 125.1;
      mh = 125;
      pspatial = 95;
      eZZeAmp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {2.844354978263298E-01,2.663967334165577E-01,2.525398170680189E-01,2.417962040400655E-01,2.334616148726172E-01,2.270608095504110E-01,2.222704276196023E-01,2.188735295765546E-01,2.167326061693552E-01,2.157741803365995E-01,2.159814609990380E-01,2.173934600897569E-01,2.201103566081738E-01,2.243061646216646E-01,2.302513527581811E-01,2.383505273367076E-01,2.492046164485062E-01,2.637151871679697E-01,2.832650285644218E-01,3.100444314319180E-01};
      i += eZZeAmp.test_2to2_amp2([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH4);
      i += eZZeAmp.test_2to2_amp2_rotations([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH4);
      i += eZZeAmp.test_2to2_amp2_boosts([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH4);
      i += eZZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH4);
      //std::cout<<"\n# me=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      me = 125;
      mh = 0.0005;
      pspatial = 125.1;
      eZZeAmp.set_masses(me,mh,MW);
      ldouble dataCH3[20] = {4.348995442070926E-01,3.945989980091484E-01,3.683925992577726E-01,3.514896623442730E-01,3.413362951485038E-01,3.365511667834200E-01,3.364586931999741E-01,3.408863237203111E-01,3.501055020023897E-01,3.648812180655990E-01,3.866466263854679E-01,4.178794651136264E-01,4.628748961928034E-01,5.294031184765088E-01,6.325894454942563E-01,8.052199766742321E-01,1.130501395610619E+00,1.878859460899439E+00,4.431867460545117E+00,3.408985823469851E+01};
      i += eZZeAmp.test_2to2_amp2([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH3);
      i += eZZeAmp.test_2to2_amp2_rotations([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH3);
      i += eZZeAmp.test_2to2_amp2_boosts([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH3);
      i += eZZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eZZeAmp.amp2(); }, me,MZ,MZ,me,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
