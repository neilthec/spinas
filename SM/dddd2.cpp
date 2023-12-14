
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

//File:  SPINAS/SM/dddd2.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dddd2.h"

namespace spinas {
  //Constructors
  dddd2::dddd2(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), gs(gscharge), md(massd), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propAG(0,0), proph(mh,wh), propZ(MZ,WZ),
    p1(particle(md)), p2(particle(md)),
    p3(particle(md)), p4(particle(md)),
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
    preh = e*e*md*md/(4.0*MW*MW*SW*SW);
    gL=-1.0+2.0/3.0*SW*SW;
    gR=2.0/3.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*md*md/MZ/MZ;//=preh!
  }
  void dddd2::set_masses(const ldouble& massd, const ldouble& massh, const ldouble& massW){
    md=massd;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(md);
    p2.set_mass(md);
    p3.set_mass(md);
    p4.set_mass(md);
    preh = e*e*md*md/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*md*md/MZ/MZ;//=preh
  }
  void dddd2::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenTAG = propAG.den(propPT);
    pDenTh = proph.den(propPT);
    pDenTZ = propZ.den(propPT);
    pDenUAG = propAG.den(propPU);
    pDenUh = proph.den(propPU);
    pDenUZ = propZ.den(propPU);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  //We have to separate the gluon and the S and T channels because each has its own color factor.
  cdouble dddd2::amp_gluon_T(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Gluon
    //dDDd all in:
    //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //ddDD: 4->2->3->4
    //- (- <14>[23] + <12>[34] - [14]<23> + [12]<34>)
    //34 out:
    //- (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude = - two*gs*gs*(
			     a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
			     + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
			     )/pDenTAG;

    return amplitude;
  }
  cdouble dddd2::amp_gluon_U(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Gluon
    //dDDd all in:
    //(<12>[34] - <14>[23] + [12]<34> - [14]<23>)
    //ddDD: 4->2->3->4
    //(- <13>[24] - <12>[34] - [13]<24> - [12]<34>)
    //34 out:
    //(<13>[24] - <12>[34] + [13]<24> - [12]<34>)
    amplitude = two*gs*gs*(
			  a13a.v(ds1,ds3)*s24s.v(ds2,ds4) - a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
			  + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) - s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
			  )/pDenUAG;

    return amplitude;
  }
  cdouble dddd2::amp_rest_T(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //dDDd all in:
    //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //ddDD: 4->2->3->4
    //- (- <14>[23] + <12>[34] - [14]<23> + [12]<34>)
    //34 out:
    //- (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude = - two*e*e*1.0/9.0*(
				   a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
				   + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
				   )/pDenTAG;
    
    //Higgs
    //preh = e*e*me*mm/(4*MW*MW*SW*SW);
    //dDDd all in:
    //preh ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //ddDD: 4->2->3->4
    //- preh ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //34 out:
    //- preh ([13]-<13>) ([24]-<24>)/(t-Mh^2)
    amplitude += - preh*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3))*(s24s.v(ds2,ds4)-a24a.v(ds2,ds4))/pDenTh;
    
    //Z Boson
    //Defined above:
    //gL=1.0-4.0/3.0*SW*SW;
    //gR=-4.0/3.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gL-gR)*(gL-gR)*md*md/MZ/MZ; // = preh
    //dDDd all in:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    //ddDD: 4->2->3->4
    //+ preZ0 (<13>-[13]) (<24>-[24]) / (t-MZ^2)
    //- preZ (gL^2 [34] <12> - gLgR( [14] <23>+ [23] <14> ) + gR^2 [12] <34>))/(t-MZ^2)
    //34 out:
    //+ preZ0 (<13>+[13]) (<24>+[24]) / (t-MZ^2)
    //- preZ (gL^2 [34] <12> + gLgR( [14] <23>+ [23] <14> ) + gR^2 [12] <34>))/(t-MZ^2)
    amplitude += 
      + preZ0*(a13a.v(ds1,ds3)+s13s.v(ds1,ds3))*(a24a.v(ds2,ds4)+s24s.v(ds2,ds4))/pDenTZ
      - two*preZ*(
	      gL*gL*s34s.v(ds3,ds4)*a12a.v(ds1,ds2)
	      + gL*gR*(s14s.v(ds1,ds4)*a23a.v(ds2,ds3)+s23s.v(ds2,ds3)*a14a.v(ds1,ds4))
	      + gR*gR*s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
	      )/pDenTZ;
    
    return amplitude;
  }
  cdouble dddd2::amp_rest_U(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //dDDd all in:
    //(<12>[34] - <14>[23] + [12]<34> - [14]<23>)
    //ddDD: 4->2->3->4
    //(- <13>[24] - <12>[34] - [13]<24> - [12]<34>)
    //34 out:
    //(<13>[24] - <12>[34] + [13]<24> - [12]<34>)
    amplitude += two*e*e*1.0/9.0*(
				  a13a.v(ds1,ds3)*s24s.v(ds2,ds4) - a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
				  + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) - s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
				  )/pDenUAG;
    
    //Higgs
    //T Channel
    //dDDd all in:
    //- preh ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //ddDD: 4->2->3->4
    //+ preh ([14]+<14>) ([23]+<23>)/(u-Mh^2)
    //34 out:
    //+ preh ([14]-<14>) ([23]-<23>)/(u-Mh^2)
    amplitude += preh*(s14s.v(ds1,ds4)-a14a.v(ds1,ds4))*(s23s.v(ds2,ds3)-a23a.v(ds2,ds3))/pDenUh;
    
    //Z Boson
    //T Channel
    //dDDd all in:
    //+ preZ0 (<13>-[13]) (<24>-[24]) / (t-MZ^2)
    //+ preZ (- gL^2 [23] <14> + gLgR( [12] <34>+ [34] <12> ) - gR^2 [14] <23>))/(t-MZ^2)
    //ddDD: 4->2->3->4
    //- preZ0 (<14>-[14]) (<23>-[23]) / (u-MZ^2)
    //- preZ (gL^2 [34] <12> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [12] <34>))/(u-MZ^2)
    //34 out:
    //- preZ0 (<14>+[14]) (<23>+[23]) / (u-MZ^2)
    //- preZ (gL^2 [34] <12> - gLgR( [13] <24>+ [24] <13> ) + gR^2 [12] <34>))/(u-MZ^2)
    amplitude += 
      - preZ0*(a14a.v(ds1,ds4)+s14s.v(ds1,ds4))*(a23a.v(ds2,ds3)+s23s.v(ds2,ds3))/pDenUZ
      - two*preZ*(
	      gL*gL*s34s.v(ds3,ds4)*a12a.v(ds1,ds2)
	      - gL*gR*(s13s.v(ds1,ds3)*a24a.v(ds2,ds4)+s24s.v(ds2,ds4)*a13a.v(ds1,ds3))
	      + gR*gR*s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
	      )/pDenUZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble dddd2::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble two=2, three = 3, four = 4, nine = 9;
    cdouble M_rest_S, M_gluon_S, M_rest_T, M_gluon_T;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M_rest_S = amp_rest_T(j1,j2,j3,j4);
	    M_rest_T = amp_rest_U(j1,j2,j3,j4);
	    M_gluon_S = amp_gluon_T(j1,j2,j3,j4);
	    M_gluon_T = amp_gluon_U(j1,j2,j3,j4);
	    amp2 += nine*std::pow(std::abs(M_rest_S),2);// Color factor Tr(1)^2 = 3*3=9
	    amp2 += two*std::pow(std::abs(M_gluon_S),2);// Color factor Tr(Ta,Tb)^2 = C^2*delta^ab*delta^ab = 1/4*8=2
	    amp2 += nine*std::pow(std::abs(M_rest_T),2);// Color factor Tr(1)^2 = 3*3=9
	    amp2 += two*std::pow(std::abs(M_gluon_T),2);// Color factor Tr(Ta,Tb)^2 = 2
	    //Cross terms
	    //No cross term (color factor Tr(Ta)^2=0) for rest_S*gluon_S and rest_T*gluon_T
	    amp2 += three*two*std::real(M_rest_S*std::conj(M_rest_T));// Color factor Tr(1) = 3 -- 2: 2real(w*conj(z)) = w*conj(z)+z*conj(w)
	    amp2 += (-two/three)*two*std::real(M_gluon_S*std::conj(M_gluon_T));// Color factor Tr(Ta,Tb,Ta,Tb) = -2/3 -- 2: 2real(w*conj(z)) = w*conj(z)+z*conj(w)
	    amp2 += four*two*std::real(M_gluon_S*std::conj(M_rest_T));// Color factor Tr(Ta,Ta) = 4 -- 2: 2real(w*conj(z)) = w*conj(z)+z*conj(w)
	    amp2 += four*two*std::real(M_gluon_T*std::conj(M_rest_S));// Color factor Tr(Ta,Ta) = 4 -- 2: 2real(w*conj(z)) = w*conj(z)+z*conj(w)
	    
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over colors 1/9
    //Symmetry factor 1/2
    return amp2/72.0;
  }

  



  //  Tests
  int test_dddd2(){
    int n=0;//Number of fails
    std::cout<<"\t* d , d  -> d , d       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### md=0.0042, pspatial=250\n";
      ldouble md=0.0042, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      dddd2 dddd2Amp = dddd2(0.31333,1.238,md,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.623354962635813E+03,1.705920817854226E+02,5.822661999130911E+01,2.831316023489040E+01,1.647375217004471E+01,1.076536061414886E+01,7.697783989510825E+00,5.969164759039909E+00,5.019795495498731E+00,4.595686416790152E+00,4.595686416790152E+00,5.019795495498729E+00,5.969164759039906E+00,7.697783989510822E+00,1.076536061414886E+01,1.647375217004470E+01,2.831316023489037E+01,5.822661999130906E+01,1.705920817854225E+02,1.623354962635806E+03};
      i += dddd2Amp.test_2to2_amp2([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH);
      i += dddd2Amp.test_2to2_amp2_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH);
      i += dddd2Amp.test_2to2_amp2_boosts([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH);
      i += dddd2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH);
      //std::cout<<"\n########### md=0.0042, pspatial=0.005\n";
      pspatial = 0.005;
      ldouble dataCH2[20] = {2.976195422500911E+03,3.133268457200552E+02,1.070148084077215E+02,5.198685680318425E+01,3.016301990474965E+01,1.961920072705204E+01,1.394238106267905E+01,1.073815535523031E+01,8.976131414741829E+00,8.188352242148122E+00,8.188352242148122E+00,8.976131414741825E+00,1.073815535523031E+01,1.394238106267905E+01,1.961920072705204E+01,3.016301990474965E+01,5.198685680318420E+01,1.070148084077214E+02,3.133268457200548E+02,2.976195422500899E+03};
      i += dddd2Amp.test_2to2_amp2([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH2);
      i += dddd2Amp.test_2to2_amp2_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH2);
      i += dddd2Amp.test_2to2_amp2_boosts([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH2);
      i += dddd2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH2);
      //std::cout<<"\n########### md=1.2, pspatial=0.3\n";
      md=1.2;
      pspatial = 0.3;
      dddd2Amp.set_masses(md,mh,MW);
      ldouble dataCH4[20] = {1.355172754614375E+05,1.516999058521940E+04,5.545826054672679E+03,2.903112292821218E+03,1.825809068091800E+03,1.291929157459660E+03,9.979348237351013E+02,8.288490346946886E+02,7.345437429301810E+02,6.920083772909296E+02,6.920083772909295E+02,7.345437429301807E+02,8.288490346946882E+02,9.979348237351010E+02,1.291929157459660E+03,1.825809068091800E+03,2.903112292821216E+03,5.545826054672675E+03,1.516999058521939E+04,1.355172754614369E+05};
      i += dddd2Amp.test_2to2_amp2([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH4);
      i += dddd2Amp.test_2to2_amp2_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH4);
      i += dddd2Amp.test_2to2_amp2_boosts([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH4);
      i += dddd2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH4);
      //std::cout<<"\n########### md=1.2, MW=2.11, pspatial=0.3\n";
      md=1.2;
      MW=2.11;
      pspatial = 0.3;
      dddd2Amp.set_masses(md,mh,MW);
      ldouble dataCH5[20] = {1.355405620378523E+05,1.517814158111313E+04,5.550978763580448E+03,2.907004358010571E+03,1.829023367640441E+03,1.294734460377288E+03,1.000480046670243E+03,8.312283265608744E+02,7.368237307623535E+02,6.942416996881998E+02,6.942416996881997E+02,7.368237307623532E+02,8.312283265608739E+02,1.000480046670243E+03,1.294734460377288E+03,1.829023367640441E+03,2.907004358010569E+03,5.550978763580444E+03,1.517814158111312E+04,1.355405620378517E+05};
      i += dddd2Amp.test_2to2_amp2([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH5);
      i += dddd2Amp.test_2to2_amp2_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH5);
      i += dddd2Amp.test_2to2_amp2_boosts([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH5);
      i += dddd2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH5);
      //std::cout<<"\n########### md=1.2, md=1.23, MW=0.006, pspatial=0.3\n";
      md=1.2;
      MW=0.006;
      pspatial = 0.3;
      dddd2Amp.set_masses(md,mh,MW);
      ldouble dataCH6[20] = {2.351708507798129E+07,2.267414581897025E+07,2.252146030678531E+07,2.245867483361864E+07,2.242526468369552E+07,2.240522027039191E+07,2.239251983924446E+07,2.238443555296431E+07,2.237960431304163E+07,2.237733578479601E+07,2.237733578479602E+07,2.237960431304162E+07,2.238443555296430E+07,2.239251983924446E+07,2.240522027039191E+07,2.242526468369552E+07,2.245867483361864E+07,2.252146030678531E+07,2.267414581897024E+07,2.351708507798124E+07};
      i += dddd2Amp.test_2to2_amp2([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH6);
      i += dddd2Amp.test_2to2_amp2_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH6);
      i += dddd2Amp.test_2to2_amp2_boosts([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH6);
      i += dddd2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH6);
      //std::cout<<"\n########### md=1.2, MW=2.11, Mh=3.125, pspatial=0.3\n";
      md=1.2;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      dddd2Amp.set_masses(md,mh,MW);
      ldouble dataCH7[20] = {1.355457272267546E+05,1.517991524999498E+04,5.552080589446201E+03,2.907823704981770E+03,1.829690844579736E+03,1.295310291818244E+03,1.000997601099032E+03,8.317087003242707E+02,7.372818532124655E+02,6.946893657032772E+02,6.946893657032772E+02,7.372818532124652E+02,8.317087003242702E+02,1.000997601099031E+03,1.295310291818243E+03,1.829690844579736E+03,2.907823704981768E+03,5.552080589446197E+03,1.517991524999497E+04,1.355457272267540E+05};
      i += dddd2Amp.test_2to2_amp2([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH7);
      i += dddd2Amp.test_2to2_amp2_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH7);
      i += dddd2Amp.test_2to2_amp2_boosts([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH7);
      i += dddd2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH7);
      //std::cout<<"\n########### md=1.2, MW=0.006, Mh=3.125, pspatial=0.3\n";
      md=1.2;
      MW=0.006;
      pspatial = 0.3;
      dddd2Amp.set_masses(md,mh,MW);
      ldouble dataCH8[20] = {2.598412631568862E+07,2.475370420227984E+07,2.452460523359960E+07,2.442988074426384E+07,2.437941460869066E+07,2.434914886295770E+07,2.432999094817749E+07,2.431780983598817E+07,2.431053731068075E+07,2.430712465513479E+07,2.430712465513479E+07,2.431053731068074E+07,2.431780983598816E+07,2.432999094817750E+07,2.434914886295770E+07,2.437941460869066E+07,2.442988074426385E+07,2.452460523359960E+07,2.475370420227983E+07,2.598412631568857E+07};
      i += dddd2Amp.test_2to2_amp2([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH8);
      i += dddd2Amp.test_2to2_amp2_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH8);
      i += dddd2Amp.test_2to2_amp2_boosts([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH8);
      i += dddd2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return dddd2Amp.amp2(); }, md,md,md,md,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
