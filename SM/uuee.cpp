
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

//File:  SPINAS/SM/uuee.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uuee.h"

namespace spinas {
  //Constructors
  uuee::uuee(const ldouble& echarge, const ldouble& massu, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), mu(massu), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    propA(0,0), proph(mh,wh),
    p1(particle(mu)), p2(particle(mu)),
    p3(particle(me)), p4(particle(me)),
    a13a(sproduct(ANGLE,&p1,&p3)),
    s13s(sproduct(SQUARE,&p1,&p3)),
    a14a(sproduct(ANGLE,&p1,&p4)),
    s14s(sproduct(SQUARE,&p1,&p4)),
    a23a(sproduct(ANGLE,&p2,&p3)),
    s23s(sproduct(SQUARE,&p2,&p3)),
    a24a(sproduct(ANGLE,&p2,&p4)),
    s24s(sproduct(SQUARE,&p2,&p4)),
    s12s(sproduct(SQUARE,&p1,&p2)),
    a12a(sproduct(ANGLE,&p1,&p2)),
    s34s(sproduct(SQUARE,&p3,&p4)),
    a34a(sproduct(ANGLE,&p3,&p4))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    preh = e*e*mu*me/(4.0*MW*MW*SW*SW);
    gLu=1.0-4.0/3.0*SW*SW;
    gRu=-4.0/3.0*SW*SW;
    gLe=-1.0+2.0*SW*SW;
    gRe=2.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLe-gRe)*mu*me/MZ/MZ;//=-preh!
  }
  void uuee::set_masses(const ldouble& massu, const ldouble& masse, const ldouble& massh, const ldouble& massW){
    mu=massu;
    me=masse;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(mu);
    p2.set_mass(mu);
    p3.set_mass(me);
    p4.set_mass(me);
    preh = e*e*mu*me/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLe-gRe)*mu*me/MZ/MZ;//=-preh
  }
  void uuee::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a13a.update();
    s13s.update();
    a14a.update();
    s14s.update();
    a23a.update();
    s23s.update();
    a24a.update();
    s24s.update();
    s12s.update();
    a12a.update();
    s34s.update();
    a34a.update();
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]+mom2[j];
    pDenSA = propA.denominator(propP);
    pDenSh = proph.denominator(propP);
    pDenSZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uuee::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //Sign changes due to p3 and p4 being outgoing.
    // (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    amplitude = -two*e*e*2.0/3.0*(
			  a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenSA;
    
    //Higgs
    //EE^2 Mu Me / (4 MW^2 SW^2) * ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //preh = e*e*mu*me/(4*MW*MW*SW*SW);
    amplitude += preh*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))*(s34s.v(ds3,ds4)+a34a.v(ds3,ds4))/pDenSh;
    
    //Z Boson
    //Defined above:
    //gLu=1.0-4.0/3.0*SW*SW;
    //gRu=-4.0/3.0*SW*SW;
    //gLe=-1.0+2.0*SW*SW;
    //gRe=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gLu-gRu)*(gLe-gRe)*mu*me/MZ/MZ; // = -preh
    //all in:
    //
    //+(EE^2 Mu Me (gLu-gRu)(gLe-gRe) (<12>-[12]) (<34>-[34]))/(8 CW^2 MZ^2 SW^2 (s-MZ^2))
    //+(EE^2 ( gLu gLe [23] <14> + gLe gRu [13] <24> + gLu gRe [24] <13> + gRu gRe [14] <23>)/(4 CW^2 SW^2 (s-MZ^2))
    //= - preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //  - preZ 2( gLu gLe [23] <14> + gLe gRu [13] <24> + gLu gRe [24] <13> + gRu gRe [14] <23> )/(s-MZ^2)
    //34 out:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //+ preZ 2( gLu gLe [23] <14> + gLe gRu [13] <24> + gLu gRe [24] <13> + gRu gRe [14] <23> )/(s-MZ^2)
    amplitude += 
      - preZ0*(a12a.v(ds1,ds2)-s12s.v(ds1,ds2))*(a34a.v(ds3,ds4)-s34s.v(ds3,ds4))/pDenSZ
      + two*preZ*(
	        gLu*gLe*s23s.v(ds2,ds3)*a14a.v(ds1,ds4)
	      + gLe*gRu*s13s.v(ds1,ds3)*a24a.v(ds2,ds4)
	      + gLu*gRe*s24s.v(ds2,ds4)*a13a.v(ds1,ds3)
	      + gRu*gRe*s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
	      )/pDenSZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uuee::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble two=2, three = 3;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += three*std::pow(std::abs(M),2);// Color factor 3
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over colors 1/3*1/3 = 1/9
    return amp2/36.0;
  }

  



  //  Tests
  int test_uuee(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> e , E       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrame\n";
      //std::cout<<"########### mu=0.0042, me=0.0005, pspatial=250\n";
      ldouble mu=0.0042, me=0.0005, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrame.
      uuee uueeAmp = uuee(0.31333,mu,me,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {6.984390188069010E-03,6.289950220814459E-03,5.635739915716351E-03,5.021759272774685E-03,4.448008291989460E-03,3.914486973360679E-03,3.421195316888339E-03,2.968133322572441E-03,2.555300990412986E-03,2.182698320409974E-03,1.850325312563403E-03,1.558181966873275E-03,1.306268283339589E-03,1.094584261962345E-03,9.231299027415437E-04,7.919052056771847E-04,7.009101707692680E-04,6.501447980177935E-04,6.396090874227613E-04,6.693030389841714E-04};
      i += uueeAmp.test_2to2_amp2([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH);
      i += uueeAmp.test_2to2_amp2_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH);
      i += uueeAmp.test_2to2_amp2_boosts([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH);
      i += uueeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH);
      //std::cout<<"########### mu=0.0042, me=1.23, pspatial=1.25\n";
      me=1.23;
      pspatial = 1.25;
      uueeAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH2[20] = {2.850952372403282E-03,2.842832003681638E-03,2.835618465727526E-03,2.829311758540945E-03,2.823911882121895E-03,2.819418836470377E-03,2.815832621586389E-03,2.813153237469932E-03,2.811380684121006E-03,2.810514961539612E-03,2.810556069725749E-03,2.811504008679417E-03,2.813358778400616E-03,2.816120378889346E-03,2.819788810145607E-03,2.824364072169399E-03,2.829846164960723E-03,2.836235088519577E-03,2.843530842845963E-03,2.851733427939879E-03};
      i += uueeAmp.test_2to2_amp2([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH2);
      i += uueeAmp.test_2to2_amp2_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH2);
      i += uueeAmp.test_2to2_amp2_boosts([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH2);
      i += uueeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH2);
      //std::cout<<"########### mu=1.23, me=0.0042, pspatial=0.005\n";
      mu=1.23;
      me=0.0042;
      pspatial = 0.005;
      uueeAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH3[20] = {2.855756837659720E-03,2.855753498636772E-03,2.855750631500775E-03,2.855748236251729E-03,2.855746312889634E-03,2.855744861414490E-03,2.855743881826296E-03,2.855743374125053E-03,2.855743338310761E-03,2.855743774383420E-03,2.855744682343029E-03,2.855746062189590E-03,2.855747913923100E-03,2.855750237543563E-03,2.855753033050975E-03,2.855756300445339E-03,2.855760039726653E-03,2.855764250894918E-03,2.855768933950134E-03,2.855774088892301E-03};
      i += uueeAmp.test_2to2_amp2([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH3);
      i += uueeAmp.test_2to2_amp2_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH3);
      i += uueeAmp.test_2to2_amp2_boosts([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH3);
      i += uueeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH3);
      //std::cout<<"########### mu=1.2, me=1.23, pspatial=0.3\n";
      mu=1.2;
      me=1.23;
      pspatial = 0.3;
      uueeAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH4[20] = {4.184466191800508E-03,4.184303009063705E-03,4.184158601157469E-03,4.184032968081799E-03,4.183926109836696E-03,4.183838026422160E-03,4.183768717838191E-03,4.183718184084789E-03,4.183686425161953E-03,4.183673441069684E-03,4.183679231807982E-03,4.183703797376846E-03,4.183747137776278E-03,4.183809253006276E-03,4.183890143066840E-03,4.183989807957972E-03,4.184108247679671E-03,4.184245462231936E-03,4.184401451614768E-03,4.184576215828167E-03};
      i += uueeAmp.test_2to2_amp2([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH4);
      i += uueeAmp.test_2to2_amp2_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH4);
      i += uueeAmp.test_2to2_amp2_boosts([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH4);
      i += uueeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH4);
      //std::cout<<"########### mu=1.2, me=1.23, MW=2.11, pspatial=0.3\n";
      mu=1.2;
      me=1.23;
      MW=2.11;
      pspatial = 0.3;
      uueeAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH5[20] = {1.083716798842535E-02,1.060370353324234E-02,1.037194487465056E-02,1.014189201265000E-02,9.913544947240682E-03,9.686903678422591E-03,9.461968206195728E-03,9.238738530560097E-03,9.017214651515696E-03,8.797396569062524E-03,8.579284283200582E-03,8.362877793929871E-03,8.148177101250390E-03,7.935182205162138E-03,7.723893105665116E-03,7.514309802759325E-03,7.306432296444764E-03,7.100260586721433E-03,6.895794673589332E-03,6.693034557048462E-03};
      i += uueeAmp.test_2to2_amp2([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH5);
      i += uueeAmp.test_2to2_amp2_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH5);
      i += uueeAmp.test_2to2_amp2_boosts([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH5);
      i += uueeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH5);
      //std::cout<<"########### mu=1.2, me=1.23, MW=0.006, pspatial=0.3\n";
      mu=1.2;
      me=1.23;
      MW=0.006;
      pspatial = 0.3;
      uueeAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH6[20] = {6.686841430748373E+06,6.686841430546358E+06,6.686841430344368E+06,6.686841430142405E+06,6.686841429940467E+06,6.686841429738555E+06,6.686841429536670E+06,6.686841429334810E+06,6.686841429132976E+06,6.686841428931168E+06,6.686841428729386E+06,6.686841428527630E+06,6.686841428325901E+06,6.686841428124197E+06,6.686841427922519E+06,6.686841427720866E+06,6.686841427519240E+06,6.686841427317641E+06,6.686841427116066E+06,6.686841426914518E+06};
      i += uueeAmp.test_2to2_amp2([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH6);
      i += uueeAmp.test_2to2_amp2_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH6);
      i += uueeAmp.test_2to2_amp2_boosts([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH6);
      i += uueeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH6);
      //std::cout<<"########### mu=1.2, me=1.23, MW=2.11, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      me=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      uueeAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH7[20] = {1.092236755079614E-02,1.068001998884545E-02,1.043937822348598E-02,1.020044225471774E-02,9.963212082540739E-03,9.727687706954964E-03,9.493869127960419E-03,9.261756345557103E-03,9.031349359745018E-03,8.802648170524164E-03,8.575652777894538E-03,8.350363181856144E-03,8.126779382408978E-03,7.904901379553044E-03,7.684729173288339E-03,7.466262763614864E-03,7.249502150532620E-03,7.034447334041605E-03,6.821098314141821E-03,6.609455090833266E-03};
      i += uueeAmp.test_2to2_amp2([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH7);
      i += uueeAmp.test_2to2_amp2_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH7);
      i += uueeAmp.test_2to2_amp2_boosts([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH7);
      i += uueeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH7);
      //std::cout<<"########### mu=1.2, me=1.23, MW=0.006, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      me=1.23;
      MW=0.006;
      pspatial = 0.3;
      uueeAmp.set_masses(mu,me,mh,MW);
      ldouble dataCH8[20] = {6.699238288563340E+06,6.699237459308634E+06,6.699236630053954E+06,6.699235800799300E+06,6.699234971544672E+06,6.699234142290070E+06,6.699233313035494E+06,6.699232483780944E+06,6.699231654526420E+06,6.699230825271922E+06,6.699229996017450E+06,6.699229166763004E+06,6.699228337508583E+06,6.699227508254189E+06,6.699226678999821E+06,6.699225849745478E+06,6.699225020491162E+06,6.699224191236872E+06,6.699223361982608E+06,6.699222532728369E+06};
      i += uueeAmp.test_2to2_amp2([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH8);
      i += uueeAmp.test_2to2_amp2_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH8);
      i += uueeAmp.test_2to2_amp2_boosts([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH8);
      i += uueeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uueeAmp.amp2(); }, mu,mu,me,me,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
