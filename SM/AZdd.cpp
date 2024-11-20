
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

//File:  SPINAS/SM/AZdd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AZdd.h"

namespace spinas {

  AZdd::AZdd(const ldouble& echarge, const ldouble& massd, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qd(-1.0/3.0), md(massd), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), prope(massd,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    p1=particle(0);
    p2=particle(MZ);
    p3=particle(md);
    p4=particle(md);
    //<12>,[12],<23>,[23],<13>,[13],<24>,[24],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    //[1341], <1341>
    s1341s = sproduct(SQUARE,&p1,&p3,&p4,&p1);
    a1341a = sproduct(ANGLE,&p1,&p3,&p4,&p1);
    //Couplings
    preTU = 2.0*e*e*Qd/(2.0*MW*SW);
    gL=-2.0*Qd*SW*SW-1.0;
    gR=-2.0*Qd*SW*SW;
  }
  void AZdd::set_masses(const ldouble& massd, const ldouble& massW){
    md=massd;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(0);
    p2.set_mass(MZ);
    p3.set_mass(md);
    p4.set_mass(md);
    prope.set_mass(md);
    //Couplings
    preTU = 2.0*e*e*Qd/(2.0*MW*SW);
  }
  void AZdd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13],<24>,[24],<14>,[14]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s24s.update();
    a24a.update();
    s14s.update();
    a14a.update();
    //[1341], <1341>
    s1341s.update();
    a1341a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT=prope.denominator(propTP);
    pDenU=prope.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AZdd::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds2a, ds2b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds2);
    ldouble normFactor=get_spin_normalization(ds2);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, i);
      
      if(ds1>0){
	//all ingoing:
	//preTU = 2.0*e*e/(2.0*MW*SW);
	//- preTU [1341] (gRe [24] <23> + gLe <24> [23])/((t-md^2) (u-md^2))
	//- preTU [12] ( gRe [14]<23>/(u-md^2) - gLe [13]<24>/(t-md^2) )
	//34 outgoing:
	//+ preTU [1341] (gRe [24] <23> + gLe <24> [23])/((t-md^2) (u-md^2))
	//+ preTU [12] ( gRe [14]<23>/(u-md^2) - gLe [13]<24>/(t-md^2) )
	
	amplitude += normFactor*preTU*s1341s.v()*(gR*s24s.v(ds2a,ds4)*a23a.v(ds2b,ds3) + gL*s23s.v(ds2a,ds3)*a24a.v(ds2b,ds4))/pDenT/pDenU;
	amplitude += normFactor*preTU*s12s.v(ds2a)*(gR*s14s.v(ds4)*a23a.v(ds2b,ds3)/pDenU - gL*s13s.v(ds3)*a24a.v(ds2b,ds4)/pDenT);

	
      }
      else if(ds1<0){
	//all ingoing:
	//preTU = 2.0*e*e/(2.0*MW*SW);
	//- preTU <1341> (gRe [24] <23> + gLe <24> [23])/((t-md^2) (u-md^2))
	//- preTU <12> ( gLe <14>[23]/(u-md^2) - gRe <13>[24]/(t-md^2) )
	//34 outgoing:
	//+ preTU <1341> (gRe [24] <23> + gLe <24> [23])/((t-md^2) (u-md^2))
	//+ preTU <12> ( gLe <14>[23]/(u-md^2) - gRe <13>[24]/(t-md^2) )
	
	amplitude += normFactor*preTU*a1341a.v()*(gR*s24s.v(ds2a,ds4)*a23a.v(ds2b,ds3) + gL*s23s.v(ds2a,ds3)*a24a.v(ds2b,ds4))/pDenT/pDenU;
	amplitude += normFactor*preTU*a12a.v(ds2a)*(gL*a14a.v(ds4)*s23s.v(ds2b,ds3)/pDenU - gR*a13a.v(ds3)*s24s.v(ds2b,ds4)/pDenT);

      }


      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AZdd::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2*1/3=1/6
    return amp2/6.0;
  }
  //set_momenta(...) must be called before amp2_Aplus().
  ldouble AZdd::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }  
  //set_momenta(...) must be called before amp2_Aminus().
  ldouble AZdd::amp2_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(-2,j2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }



  //  Tests
  int test_AZdd(){
    int n=0;//Number of fails
    std::cout<<"\t* A , Z  -> d , D       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# md=0.0075, MW=80.385, pspatial=250\n";
      ldouble md=0.0075;
      ldouble EE=0.31333, MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      AZdd AZddAmp = AZdd(EE,md,MW,SW);
      ldouble pspatial=250;
      ldouble dataCHp[20] = {6.339981599570700E-03,2.559433006521695E-03,2.043377866941699E-03,2.024854743603380E-03,2.202681068279126E-03,2.501601387982142E-03,2.900655041640943E-03,3.399747329199391E-03,4.011091634733001E-03,4.757956545732457E-03,5.677114651871997E-03,6.824984864894147E-03,8.289668453082899E-03,1.021463962732194E-02,1.284861316568451E-02,1.666212523863122E-02,2.266413218047140E-02,3.347858499677635E-02,5.872744556225561E-02,1.850101724346488E-01};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCHp);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCHp);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCHp);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCHp);
      ldouble dataCHm[20] = {1.850101724346493E-01,5.872744556225563E-02,3.347858499677635E-02,2.266413218047141E-02,1.666212523863123E-02,1.284861316568449E-02,1.021463962732194E-02,8.289668453082903E-03,6.824984864894149E-03,5.677114651872000E-03,4.757956545732458E-03,4.011091634733003E-03,3.399747329199392E-03,2.900655041640942E-03,2.501601387982145E-03,2.202681068279126E-03,2.024854743603379E-03,2.043377866941698E-03,2.559433006521693E-03,6.339981599570637E-03};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCHm);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCHm);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCHm);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCHm);
      ldouble dataCH[20] = {9.567507701710998E-02,3.064343928438866E-02,1.776098143185903E-02,1.234449346203740E-02,9.432403153455176E-03,7.675107276833317E-03,6.557647334481442E-03,5.844707891141147E-03,5.418038249813574E-03,5.217535598802228E-03,5.217535598802227E-03,5.418038249813575E-03,5.844707891141145E-03,6.557647334481438E-03,7.675107276833328E-03,9.432403153455169E-03,1.234449346203739E-02,1.776098143185902E-02,3.064343928438865E-02,9.567507701710973E-02};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH);
      //std::cout<<"\n# md=0.0042, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2p[20] = {1.055173214873646E-02,4.577297846612861E-03,3.664565241206136E-03,3.511575926766089E-03,3.647689667753213E-03,3.952642787838753E-03,4.389563914908140E-03,4.952632189989589E-03,5.653180205541048E-03,6.516807271380238E-03,7.585651613525952E-03,8.925320045277380E-03,1.063888408726429E-02,1.289464153403732E-02,1.598468048334876E-02,2.046187533576176E-02,2.751198965027456E-02,4.021904136452575E-02,6.989235208789678E-02,2.183186130064260E-01};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH2p);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH2p);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH2p);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH2p);
      ldouble dataCH2m[20] = {2.183186130064270E-01,6.989235208789689E-02,4.021904136452579E-02,2.751198965027458E-02,2.046187533576177E-02,1.598468048334876E-02,1.289464153403733E-02,1.063888408726430E-02,8.925320045277386E-03,7.585651613525954E-03,6.516807271380241E-03,5.653180205541048E-03,4.952632189989591E-03,4.389563914908140E-03,3.952642787838753E-03,3.647689667753212E-03,3.511575926766087E-03,3.664565241206134E-03,4.577297846612855E-03,1.055173214873641E-02};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH2m);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH2m);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH2m);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH2m);
      ldouble dataCH2[20] = {1.144351725775817E-01,3.723482496725487E-02,2.194180330286596E-02,1.551178278852033E-02,1.205478250175749E-02,9.968661635593758E-03,8.642102724472734E-03,7.795758138626943E-03,7.289250125409217E-03,7.051229442453096E-03,7.051229442453096E-03,7.289250125409213E-03,7.795758138626941E-03,8.642102724472732E-03,9.968661635593757E-03,1.205478250175749E-02,1.551178278852032E-02,2.194180330286594E-02,3.723482496725482E-02,1.144351725775812E-01};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH2);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH2);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH2);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH2);
      //std::cout<<"\n# md=125.1, MW=80.385, pspatial=95\n";
      md = 125;
      pspatial = 250;
      AZddAmp.set_masses(md,MW);
      ldouble dataCH4p[20] = {1.098806628940977E-01,6.594509890286443E-02,4.810200256225918E-02,3.885704935821192E-02,3.345357253752718E-02,3.011060189937588E-02,2.802700336208176E-02,2.680032265547784E-02,2.621730617656152E-02,2.616742170723430E-02,2.660474380121407E-02,2.753232227765699E-02,2.899978757776802E-02,3.111226830877193E-02,3.405401061227468E-02,3.813759737691318E-02,4.390556132627795E-02,5.234995494617890E-02,6.540105622172149E-02,8.665178495065155E-02};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH4p);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH4p);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH4p);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH4p);
      ldouble dataCH4m[20] = {8.665178495065158E-02,6.540105622172150E-02,5.234995494617891E-02,4.390556132627796E-02,3.813759737691318E-02,3.405401061227469E-02,3.111226830877193E-02,2.899978757776803E-02,2.753232227765700E-02,2.660474380121406E-02,2.616742170723430E-02,2.621730617656151E-02,2.680032265547783E-02,2.802700336208176E-02,3.011060189937587E-02,3.345357253752718E-02,3.885704935821192E-02,4.810200256225917E-02,6.594509890286442E-02,1.098806628940976E-01};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH4m);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH4m);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH4m);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH4m);
      ldouble dataCH4[20] = {9.826622392237462E-02,6.567307756229296E-02,5.022597875421905E-02,4.138130534224494E-02,3.579558495722018E-02,3.208230625582528E-02,2.956963583542685E-02,2.790005511662293E-02,2.687481422710926E-02,2.638608275422418E-02,2.638608275422418E-02,2.687481422710926E-02,2.790005511662293E-02,2.956963583542684E-02,3.208230625582528E-02,3.579558495722018E-02,4.138130534224493E-02,5.022597875421903E-02,6.567307756229296E-02,9.826622392237457E-02};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH4);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH4);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH4);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH4);
      //std::cout<<"\n# md=125, MW=80.385, pspatial=125.1\n";
      md = 125;
      pspatial = 125.1;
      AZddAmp.set_masses(md,MW);
      ldouble dataCH3p[20] = {3.890524371349773E-02,3.592394088850055E-02,3.355444578603105E-02,3.166169732839630E-02,3.014912063165381E-02,2.894640884025673E-02,2.800168108324577E-02,2.727631698263540E-02,2.674147978101941E-02,2.637573521489460E-02,2.616340013760091E-02,2.609338901895682E-02,2.615840751133831E-02,2.635439165396818E-02,2.668012043476368E-02,2.713694381866091E-02,2.772856897614826E-02,2.846083036969282E-02,2.934132273626563E-02,3.037867219219167E-02};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH3p);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH3p);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH3p);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2_Aplus(); }, 0,MZ,md,md,pspatial,dataCH3p);
      ldouble dataCH3m[20] = {3.037867219219166E-02,2.934132273626562E-02,2.846083036969281E-02,2.772856897614825E-02,2.713694381866090E-02,2.668012043476367E-02,2.635439165396817E-02,2.615840751133830E-02,2.609338901895682E-02,2.616340013760090E-02,2.637573521489459E-02,2.674147978101940E-02,2.727631698263539E-02,2.800168108324576E-02,2.894640884025672E-02,3.014912063165380E-02,3.166169732839628E-02,3.355444578603103E-02,3.592394088850054E-02,3.890524371349773E-02};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH3m);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH3m);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH3m);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2_Aminus(); }, 0,MZ,md,md,pspatial,dataCH3m);
      ldouble dataCH3[20] = {3.464195795284470E-02,3.263263181238309E-02,3.100763807786193E-02,2.969513315227227E-02,2.864303222515736E-02,2.781326463751020E-02,2.717803636860697E-02,2.671736224698685E-02,2.641743439998812E-02,2.626956767624775E-02,2.626956767624775E-02,2.641743439998812E-02,2.671736224698685E-02,2.717803636860697E-02,2.781326463751020E-02,2.864303222515736E-02,2.969513315227227E-02,3.100763807786193E-02,3.263263181238309E-02,3.464195795284470E-02};
      i += AZddAmp.test_2to2_amp2([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH3);
      i += AZddAmp.test_2to2_amp2_rotations([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH3);
      i += AZddAmp.test_2to2_amp2_boosts([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH3);
      i += AZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
