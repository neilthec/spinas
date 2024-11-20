
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

//File:  SPINAS/SM/dddd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/dddd.h"

namespace spinas {
  //Constructors
  dddd::dddd(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
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
  void dddd::set_masses(const ldouble& massd, const ldouble& massh, const ldouble& massW){
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
  void dddd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    ldouble propPS[4], propPT[4];
    for(int j=0;j<4;j++){
      propPS[j] = mom1[j]+mom2[j];
      propPT[j] = mom1[j]-mom3[j];
    }
    pDenSAG = propAG.denominator(propPS);
    pDenSh = proph.denominator(propPS);
    pDenSZ = propZ.denominator(propPS);
    pDenTAG = propAG.denominator(propPT);
    pDenTh = proph.denominator(propPT);
    pDenTZ = propZ.denominator(propPT);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  //We have to separate the gluon and the S and T channels because each has its own color factor.
  cdouble dddd::amp_gluon_S(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Gluon
    //Sign changes due to p3 and p4 being outgoing.
    // (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    amplitude = two*gs*gs*(
			  a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenSAG;

    return amplitude;
  }
  cdouble dddd::amp_gluon_T(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Gluon
    // 2 <-> 3 with a minus sign for fermions
    // all ingoing:  (<12>[34] - <14>[23] + [12]<34> - [14]<23>)
    // 34 outgoing:  (<12>[34] + <14>[23] + [12]<34> + [14]<23>)/pDenT
    amplitude = two*gs*gs*(
			  a12a.v(ds1,ds2)*s34s.v(ds3,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s12s.v(ds1,ds2)*a34a.v(ds3,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenTAG;

    return amplitude;
  }
  cdouble dddd::amp_rest_S(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //Sign changes due to p3 and p4 being outgoing.
    // (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    amplitude = two*e*e*1.0/9.0*(
			  a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenSAG;
    
    //Higgs
    //EE^2 Me Mm / (4 MW^2 SW^2) * ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //preh = e*e*me*mm/(4*MW*MW*SW*SW);
    amplitude += preh*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))*(s34s.v(ds3,ds4)+a34a.v(ds3,ds4))/pDenSh;
    
    //Z Boson
    //Defined above:
    //gL=1.0-4.0/3.0*SW*SW;
    //gR=-4.0/3.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gL-gR)*(gL-gR)*md*md/MZ/MZ; // = preh
    //all in:
    //+(EE^2 Md Md (gL-gR)^2 (<12>-[12]) (<34>-[34]))/(8 CW^2 MZ^2 SW^2 (s-MZ^2))
    //+(EE^2 (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>)/(4 CW^2 SW^2 (s-MZ^2))
    //= + preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //  + preZ (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    //34 out:
    //+ preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    amplitude += 
      - preZ0*(a12a.v(ds1,ds2)-s12s.v(ds1,ds2))*(a34a.v(ds3,ds4)-s34s.v(ds3,ds4))/pDenSZ
      + two*preZ*(
	      gL*gL*s23s.v(ds2,ds3)*a14a.v(ds1,ds4)
	      + gL*gR*(s13s.v(ds1,ds3)*a24a.v(ds2,ds4)+s24s.v(ds2,ds4)*a13a.v(ds1,ds3))
	      + gR*gR*s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
	      )/pDenSZ;
    
    return amplitude;
  }
  cdouble dddd::amp_rest_T(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //T Channel
    // 2 <-> 3 with a minus sign for fermions
    // all ingoing:  (<12>[34] - <14>[23] + [12]<34> - [14]<23>)
    // 34 outgoing:  (<12>[34] + <14>[23] + [12]<34> + [14]<23>)/pDenT
    amplitude += two*e*e*1.0/9.0*(
			  a12a.v(ds1,ds2)*s34s.v(ds3,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s12s.v(ds1,ds2)*a34a.v(ds3,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenTAG;
    
    //Higgs
    //T Channel
    // 2 <-> 3 with a minus sign for fermions
    //all ingoing:  - preh * ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //34 outgoing:  - preh * ([13]-<13>) ([24]-<24>)/(t-Mh^2)
    amplitude += - preh*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3))*(s24s.v(ds2,ds4)-a24a.v(ds2,ds4))/pDenTh;
    
    //Z Boson
    //T Channel
    // 2 <-> 3 with a minus sign for fermions
    //all in:
    //= - preZ0 (<13>-[13]) (<24>-[24]) / (t-MZ^2)
    //  - preZ (- gL^2 [23] <14> + gLgR( [12] <34>+ [34] <12> ) - gR^2 [14] <23>))/(t-MZ^2)
    //34 outgoing:
    //= - preZ0 (<13>+[13]) (<24>+[24]) / (t-MZ^2)
    //  - preZ ( gL^2 [23] <14> + gLgR( [12] <34>+ [34] <12> ) + gR^2 [14] <23>))/(t-MZ^2)
    amplitude += 
      + preZ0*(a13a.v(ds1,ds3)+s13s.v(ds1,ds3))*(a24a.v(ds2,ds4)+s24s.v(ds2,ds4))/pDenTZ
      + two*preZ*(
	      gL*gL*s23s.v(ds2,ds3)*a14a.v(ds1,ds4)
	      + gL*gR*(s12s.v(ds1,ds2)*a34a.v(ds3,ds4)+s34s.v(ds3,ds4)*a12a.v(ds1,ds2))
	      + gR*gR*s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
	      )/pDenTZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble dddd::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble two=2, three = 3, four = 4, nine = 9;
    cdouble M_rest_S, M_gluon_S, M_rest_T, M_gluon_T;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M_rest_S = amp_rest_S(j1,j2,j3,j4);
	    M_rest_T = amp_rest_T(j1,j2,j3,j4);
	    M_gluon_S = amp_gluon_S(j1,j2,j3,j4);
	    M_gluon_T = amp_gluon_T(j1,j2,j3,j4);
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
    return amp2/36.0;
  }

  



  //  Tests
  int test_dddd(){
    int n=0;//Number of fails
    std::cout<<"\t* d , D  -> d , D       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### md=0.0042, pspatial=250\n";
      ldouble md=0.0042, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      dddd ddddAmp = dddd(0.31333,1.238,md,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.280196742439907E+03,3.507201178134367E+02,1.214732809133873E+02,5.965919534190422E+01,3.476799603629121E+01,2.244465485655904E+01,1.551861646502431E+01,1.127875172563933E+01,8.520660660475707E+00,6.644971821112491E+00,5.327068598620898E+00,4.379210419040074E+00,3.686966890852712E+00,3.177529766739191E+00,2.802884084520010E+00,2.530390123359564E+00,2.337275063417056E+00,2.207287629062720E+00,2.128599103244099E+00,2.092448351230577E+00};
      i += ddddAmp.test_2to2_amp2([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH);
      i += ddddAmp.test_2to2_amp2_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH);
      i += ddddAmp.test_2to2_amp2_boosts([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH);
      i += ddddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH);
      //std::cout<<"\n########### md=0.0042, pspatial=0.005\n";
      pspatial = 0.005;
      ldouble dataCH2[20] = {6.015250840294934E+03,6.467171692746253E+02,2.252994432274650E+02,1.112505379021992E+02,6.514726541425375E+01,4.222907397105967E+01,2.929048792344175E+01,2.132743142453092E+01,1.611158868995814E+01,1.253126949716057E+01,9.982826505170776E+00,8.116206010577109E+00,6.717454893524716E+00,5.650053941204408E+00,4.823708884288839E+00,4.176854633077796E+00,3.666404466281325E+00,3.261513228290784E+00,2.939656956945574E+00,2.684097384143378E+00};
      i += ddddAmp.test_2to2_amp2([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH2);
      i += ddddAmp.test_2to2_amp2_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH2);
      i += ddddAmp.test_2to2_amp2_boosts([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH2);
      i += ddddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH2);
      //std::cout<<"\n########### md=1.2, pspatial=0.3\n";
      md=1.2;
      pspatial = 0.3;
      ddddAmp.set_masses(md,mh,MW);
      ldouble dataCH4[20] = {2.695842535632465E+05,2.971422449060211E+04,1.061128269575677E+04,5.370335392465041E+03,3.222483486325252E+03,2.139721446935554E+03,1.519537735066558E+03,1.132035357789601E+03,8.741334889685027E+02,6.940492958218922E+02,5.634686163585163E+02,4.658567561161001E+02,3.910360698339728E+02,3.324659476974399E+02,2.857891551369601E+02,2.480131269364518E+02,2.170280483897213E+02,1.913122715213296E+02,1.697462706237379E+02,1.514917473950011E+02};
      i += ddddAmp.test_2to2_amp2([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH4);
      i += ddddAmp.test_2to2_amp2_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH4);
      i += ddddAmp.test_2to2_amp2_boosts([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH4);
      i += ddddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH4);
      //std::cout<<"\n########### md=1.2, MW=2.11, pspatial=0.3\n";
      md=1.2;
      MW=2.11;
      pspatial = 0.3;
      ddddAmp.set_masses(md,mh,MW);
      ldouble dataCH5[20] = {2.691314255004546E+05,2.956554247446141E+04,1.052341873036909E+04,5.308528369465182E+03,3.175146396049303E+03,2.101587488526045E+03,1.487770957265093E+03,1.104934205866996E+03,8.505969711071427E+02,6.733241145594455E+02,5.450167136293753E+02,4.492804474047923E+02,3.760331246668007E+02,3.188013069887004E+02,2.732764011740592E+02,2.365019307888662E+02,2.063954211331385E+02,1.814563151117678E+02,1.605816177243597E+02,1.429461704778715E+02};
      i += ddddAmp.test_2to2_amp2([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH5);
      i += ddddAmp.test_2to2_amp2_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH5);
      i += ddddAmp.test_2to2_amp2_boosts([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH5);
      i += ddddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH5);
      //std::cout<<"\n########### md=1.2, md=1.23, MW=0.006, pspatial=0.3\n";
      md=1.2;
      MW=0.006;
      pspatial = 0.3;
      ddddAmp.set_masses(md,mh,MW);
      ldouble dataCH6[20] = {4.271546050230096E+07,4.387626109858434E+07,4.413741601507217E+07,4.425228620189720E+07,4.431685958714175E+07,4.435822821461102E+07,4.438699048281719E+07,4.440814384391049E+07,4.442435273094194E+07,4.443716764570987E+07,4.444755218338761E+07,4.445613671640427E+07,4.446335096933167E+07,4.446949796551257E+07,4.447479751662409E+07,4.447941295018253E+07,4.448346815722834E+07,4.448705881730759E+07,4.449025999598312E+07,4.449313141259501E+07};
      i += ddddAmp.test_2to2_amp2([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH6);
      i += ddddAmp.test_2to2_amp2_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH6);
      i += ddddAmp.test_2to2_amp2_boosts([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH6);
      i += ddddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH6);
      //std::cout<<"\n########### md=1.2, MW=2.11, Mh=3.125, pspatial=0.3\n";
      md=1.2;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      ddddAmp.set_masses(md,mh,MW);
      ldouble dataCH7[20] = {2.691333465568053E+05,2.956609895718539E+04,1.052370243422359E+04,5.308695264852588E+03,3.175248472889325E+03,2.101648378425918E+03,1.487803384789618E+03,1.104945805628550E+03,8.505926830380373E+02,6.733073185146125E+02,5.449898237020637E+02,4.492452478441311E+02,3.759909714283352E+02,3.187532546365795E+02,2.732232860238464E+02,2.364444272300873E+02,2.063340808843887E+02,1.813915952600731E+02,1.605139011863305E+02,1.428757812144056E+02};
      i += ddddAmp.test_2to2_amp2([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH7);
      i += ddddAmp.test_2to2_amp2_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH7);
      i += ddddAmp.test_2to2_amp2_boosts([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH7);
      i += ddddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH7);
      //std::cout<<"\n########### md=1.2, MW=0.006, Mh=3.125, pspatial=0.3\n";
      md=1.2;
      MW=0.006;
      pspatial = 0.3;
      ddddAmp.set_masses(md,mh,MW);
      ldouble dataCH8[20] = {4.708144302667494E+07,4.800353911914250E+07,4.819577978407289E+07,4.826611458860529E+07,4.829435884724329E+07,4.830319221123844E+07,4.830151532342476E+07,4.829353289096233E+07,4.828148929014210E+07,4.826669151667900E+07,4.824995214570561E+07,4.823180280877794E+07,4.821260586034089E+07,4.819261670817567E+07,4.817202049327806E+07,4.815095464357247E+07,4.812952326547867E+07,4.810780662339706E+07,4.808586755779622E+07,4.806375593642240E+07};
      i += ddddAmp.test_2to2_amp2([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH8);
      i += ddddAmp.test_2to2_amp2_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH8);
      i += ddddAmp.test_2to2_amp2_boosts([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH8);
      i += ddddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ddddAmp.amp2(); }, md,md,md,md,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
