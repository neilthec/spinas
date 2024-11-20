
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

//File:  SPINAS/SM/udud.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/udud.h"

namespace spinas {
  //Constructors
  udud::udud(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& widthW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), gs(gscharge), mu(massu), md(massd), mh(massh), wh(widthh), MW(massW), WW(widthW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propAG(0,0), proph(mh,wh), propZ(MZ,WZ), propW(MW,WW),
    p1(particle(mu)), p2(particle(md)),
    p3(particle(mu)), p4(particle(md)),
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
    preh = e*e*mu*md/(4.0*MW*MW*SW*SW);
    gLu=1.0-4.0/3.0*SW*SW;
    gRu=-4.0/3.0*SW*SW;
    gLd=-1.0+2.0/3.0*SW*SW;
    gRd=2.0/3.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*md/MZ/MZ;//=-preh!
    preW = e*e/(2.0*MW*MW*SW*SW);
  }
  void udud::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& massW){
    mu=massu;
    md=massd;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    propW.set_mass(MW);
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(mu);
    p2.set_mass(md);
    p3.set_mass(mu);
    p4.set_mass(md);
    preh = e*e*mu*md/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*md/MZ/MZ;//=-preh
    preW = e*e/(2.0*MW*MW*SW*SW);
  }
  void udud::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    ldouble propPT[4], propPU[4];
    for(int j=0;j<4;j++){
      propPT[j] = mom1[j]-mom3[j];
      propPU[j] = mom1[j]-mom4[j];
    }
    pDenTAG = propAG.denominator(propPT);
    pDenTh = proph.denominator(propPT);
    pDenTZ = propZ.denominator(propPT);
    pDenUW = propW.denominator(propPU);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  //We have to separate the gluon so we can separate the color factor between the gluon^2 diagram
  // And the rest^2.
  // The cross term vanishes due to the trace of the adjoint rep.
  cdouble udud::amp_gluon(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Gluon
    //uUDd:
    //all in:
    // - (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //udUD: 4->2->3->4
    // - (- <14>[23] + <12>[34] - [14]<23> + [12]<34>)
    //34 out:
    // - (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude = - two*gs*gs*(
			     + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
			     + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
			     )/pDenTAG;

    return amplitude;
  }
  cdouble udud::amp_rest(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //uUDd:
    //all in:
    // - (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //udUD: 4->2->3->4
    // - (- <14>[23] + <12>[34] - [14]<23> + [12]<34>)
    //34 out:
    // - (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude = two*e*e*2.0/9.0*(
				 + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
				 + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
				 )/pDenTAG;
    
    //Higgs
    //preh = e*e*mu*md/(4*MW*MW*SW*SW);
    //uUDd:
    //all in:
    //preh ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //udUD: 4->2->3->4
    //- preh ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //34 out:
    //- preh ([13]-<13>) ([24]-<24>)/(t-Mh^2)
    amplitude += - preh*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3))*(s24s.v(ds2,ds4)-a24a.v(ds2,ds4))/pDenTh;
    
    //Z Boson
    //Defined above:
    //gLu=1.0-4.0/3.0*SW*SW;
    //gRu=-4.0/3.0*SW*SW;
    //gLd=-1.0+2.0/3.0*SW*SW;
    //gRd=2.0/3.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*md/MZ/MZ; // = -preh
    //uUDd:
    //all in:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ ( gLu gLd [23] <14> + gLd gRu [13] <24> + gLu gRd [24] <13> + gRu gRd [14] <23> )/(s-MZ^2)
    //udUD: 4->2->3->4
    //+ preZ0 (<13>-[13]) (<24>-[24]) / (t-MZ^2)
    //- preZ ( gLu gLd [34] <12> - gLd gRu [14] <23> - gLu gRd [23] <14> + gRu gRd [12] <34> )/(t-MZ^2)
    //34 out:
    //+ preZ0 (<13>+[13]) (<24>+[24]) / (t-MZ^2)
    //- preZ ( gLu gLd [34] <12> + gLd gRu [14] <23> + gLu gRd [23] <14> + gRu gRd [12] <34> )/(t-MZ^2)
    amplitude += 
      + preZ0*(a13a.v(ds1,ds3)+s13s.v(ds1,ds3))*(a24a.v(ds2,ds4)+s24s.v(ds2,ds4))/pDenTZ
      - two*preZ*(
	        gLu*gLd*s34s.v(ds3,ds4)*a12a.v(ds1,ds2)
	      + gLd*gRu*s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
	      + gLu*gRd*s23s.v(ds2,ds3)*a14a.v(ds1,ds4)
	      + gRu*gRd*s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
	      )/pDenTZ;

    return amplitude;
  }
  
  
  cdouble udud::amp_W(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble sqrt2=std::sqrt(2.0);
    //W Boson
    //preW = e*e/(2.0*MW*MW*SW*SW);
    //uUDd:
    //all in:
    //- preW ( 2 MW^2 [23] <14> + (Md <13>-Mu [13]) (Md [24]-Mu <24>) )/(t-MW^2)
    //udUD: 4->2->3->4
    //- preW ( 2 MW^2 [34] <12> - (Md <14>-Mu [14]) (Md [23]-Mu <23>) )/(u-MW^2)
    //34 out:
    //- preW ( 2 MW^2 [34] <12> + (Md <14>+Mu [14]) (Md [23]+Mu <23>) )/(u-MW^2)
    return - preW*(
		   2.0*MW*MW*a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
		   + (mu*s14s.v(ds1,ds4)+md*a14a.v(ds1,ds4))*(mu*a23a.v(ds2,ds3)+md*s23s.v(ds2,ds3))
		   )/pDenUW;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble udud::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble two=2, three = 3, four = 4, nine = 9;
    cdouble M_rest, M_W, M_gluon;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M_rest = amp_rest(j1,j2,j3,j4);
	    M_W = amp_W(j1,j2,j3,j4);
	    M_gluon = amp_gluon(j1,j2,j3,j4);
	    amp2 += nine*std::pow(std::abs(M_rest),2);// Color factor 3*3=9
	    amp2 += three*two*std::real(M_rest*std::conj(M_W));//Color factor 3
	    amp2 += nine*std::pow(std::abs(M_W),2);//Color factor 9
	    //Cross term with gluon and rest color factor 0 (Trace(Ta)*Trace(Ta)=0)
	    amp2 += four*two*std::real(M_gluon*std::conj(M_W));//Color factor 4
	    amp2 += two*std::pow(std::abs(M_gluon),2);//Color factor C^2*delta^ab*delta^ab = 1/4*8=2
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over colors 1/3*1/3 = 1/9
    return amp2/36.0;
  }

  



  //  Tests
  int test_udud(){
    int n=0;//Number of fails
    std::cout<<"\t* u , d  -> u , d       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mu=0.0042, md=0.0075, pspatial=250\n";
      ldouble mu=0.0042, md=0.0075, mh=125, wh=0, MW=80.385, WW=0, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      udud ududAmp = udud(0.31333,1.238,mu,md,mh,wh,MW,WW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.276658721524500E+03,3.497459109556606E+02,1.211094768689896E+02,5.957247332932436E+01,3.485262826414292E+01,2.265391685151779E+01,1.582751670901319E+01,1.167268844625483E+01,8.991469319824606E+00,7.189772539576206E+00,5.948619696406961E+00,5.087585354569023E+00,4.503065508771460E+00,4.139914155486395E+00,3.981645243366304E+00,4.057966502563196E+00,4.487023223277483E+00,5.639271925828444E+00,8.964205785029229E+00,2.527542147094941E+01};
      i += ududAmp.test_2to2_amp2([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH);
      i += ududAmp.test_2to2_amp2_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH);
      i += ududAmp.test_2to2_amp2_boosts([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH);
      i += ududAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH);
      //std::cout<<"########### mu=0.0042, md=0.0075, pspatial=0.004\n";
      pspatial = 0.004;
      ldouble dataCH2[20] = {1.365820699914168E+04,1.459151287370682E+03,5.045950312263190E+02,2.470561139205552E+02,1.432713949467134E+02,9.183947112717115E+01,6.289070194804111E+01,4.512373715189995E+01,3.351388476215678E+01,2.555829401794013E+01,1.990007980003788E+01,1.575365945233315E+01,1.263967378540635E+01,1.025276541646211E+01,8.391314469024904E+00,6.918090068081318E+00,5.737220614811263E+00,4.780197184293779E+00,3.997085068748989E+00,3.350846137147871E+00};
      i += ududAmp.test_2to2_amp2([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH2);
      i += ududAmp.test_2to2_amp2_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH2);
      i += ududAmp.test_2to2_amp2_boosts([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH2);
      i += ududAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH2);
      //std::cout<<"########### mu=1.23, md=0.0042, pspatial=0.005\n";
      mu=1.23;
      md=0.0042;
      pspatial = 0.005;
      ududAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH3[20] = {8.554579574060710E+07,9.221042637535851E+06,3.217319203763663E+06,1.589318057334964E+06,9.298789862196008E+05,6.013538786546729E+05,4.154286929567790E+05,3.006714989387955E+05,2.252412933280807E+05,1.732365336328379E+05,1.360137230615578E+05,1.085552428315538E+05,8.779104544750448E+04,7.176007093352733E+04,5.916383873418027E+04,4.911601810551775E+04,4.099559851461209E+04,3.435749267304111E+04,2.887630286748811E+04,2.430992576315572E+04};
      i += ududAmp.test_2to2_amp2([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH3);
      i += ududAmp.test_2to2_amp2_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH3);
      i += ududAmp.test_2to2_amp2_boosts([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH3);
      i += ududAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH3);
      //std::cout<<"########### mu=1.2, md=1.23, pspatial=0.3\n";
      mu=1.2;
      md=1.23;
      pspatial = 0.3;
      ududAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH4[20] = {2.815780078120495E+05,3.096409918367171E+04,1.103136544136370E+04,5.569382407246357E+03,3.333625906368681E+03,2.207900970383047E+03,1.563885768972831E+03,1.161981026826327E+03,8.948205039957601E+02,7.085010493430543E+02,5.735668719047803E+02,4.728269579322347E+02,3.957060632338486E+02,3.354121767282559E+02,2.874233325651161E+02,2.486357201645457E+02,2.168624334319344E+02,1.905271930407249E+02,1.684709945073363E+02,1.498265154662878E+02};
      i += ududAmp.test_2to2_amp2([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH4);
      i += ududAmp.test_2to2_amp2_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH4);
      i += ududAmp.test_2to2_amp2_boosts([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH4);
      i += ududAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH4);
      //std::cout<<"########### mu=1.2, md=1.23, MW=2.11, pspatial=0.3\n";
      mu=1.2;
      md=1.23;
      MW=2.11;
      pspatial = 0.3;
      ududAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH5[20] = {2.816555270254476E+05,3.098988038701960E+04,1.104679908101133E+04,5.580381405254339E+03,3.342161237668185E+03,2.214868533427101E+03,1.569767971271569E+03,1.167067313459030E+03,8.992981634690246E+02,7.124982268922743E+02,5.771751002713415E+02,4.761138921508365E+02,3.987231215050977E+02,3.381993513386374E+02,2.900123418692255E+02,2.510521437990330E+02,2.191272004950878E+02,1.926576450780505E+02,1.704816612038954E+02,1.517296914437174E+02};
      i += ududAmp.test_2to2_amp2([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH5);
      i += ududAmp.test_2to2_amp2_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH5);
      i += ududAmp.test_2to2_amp2_boosts([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH5);
      i += ududAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH5);
      //std::cout<<"########### mu=1.2, md=1.23, MW=0.006, pspatial=0.3\n";
      mu=1.2;
      md=1.23;
      MW=0.006;
      pspatial = 0.3;
      ududAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH6[20] = {9.252336266767460E+07,8.934738445405926E+07,8.878361349932419E+07,8.857978902990051E+07,8.849959534417202E+07,8.848073244231939E+07,8.850081221367596E+07,8.855112695517008E+07,8.862922615419436E+07,8.873653376447183E+07,8.887786250371973E+07,8.906207118680136E+07,8.930399592970949E+07,8.962860103949787E+07,9.007998600831628E+07,9.074281760713050E+07,9.180160001585545E+07,9.374758543339592E+07,9.846531358779052E+07,1.260249887079488E+08};
      i += ududAmp.test_2to2_amp2([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH6);
      i += ududAmp.test_2to2_amp2_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH6);
      i += ududAmp.test_2to2_amp2_boosts([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH6);
      i += ududAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH6);
      //std::cout<<"########### mu=1.2, md=1.23, MW=2.11, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      md=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      ududAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH7[20] = {2.816562279978520E+05,3.099011287459264E+04,1.104693787402431E+04,5.580480045195970E+03,3.342237571717878E+03,2.214930673952569E+03,1.569820286482112E+03,1.167112424289975E+03,8.993377658282558E+02,7.125334811454417E+02,5.772068352579203E+02,4.761427204975562E+02,3.987495088140022E+02,3.382236597499391E+02,2.900348586024879E+02,2.510731004875215E+02,2.191467866612286E+02,1.926760177587598E+02,1.704989519744698E+02,1.517460116407280E+02};
      i += ududAmp.test_2to2_amp2([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH7);
      i += ududAmp.test_2to2_amp2_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH7);
      i += ududAmp.test_2to2_amp2_boosts([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH7);
      i += ududAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH7);
      //std::cout<<"########### mu=1.2, md=1.23, MW=0.006, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      md=1.23;
      MW=0.006;
      pspatial = 0.3;
      ududAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH8[20] = {9.186010357237682E+07,8.860584571516001E+07,8.803072997313763E+07,8.782486487536305E+07,8.774546900009432E+07,8.772841970683087E+07,8.775055976164722E+07,8.780273995233406E+07,8.788217139863750E+07,8.798994379831356E+07,8.813047158538750E+07,8.831207578994839E+07,8.854879949673867E+07,8.886434678560141E+07,8.930065161417821E+07,8.993828641409777E+07,9.095293471301453E+07,9.281266973114601E+07,9.731514185941888E+07,1.236796329530243E+08};
      i += ududAmp.test_2to2_amp2([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH8);
      i += ududAmp.test_2to2_amp2_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH8);
      i += ududAmp.test_2to2_amp2_boosts([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH8);
      i += ududAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ududAmp.amp2(); }, mu,md,mu,md,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
