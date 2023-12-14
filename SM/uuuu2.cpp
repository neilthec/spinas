
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

//File:  SPINAS/SM/uuuu2.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uuuu2.h"

namespace spinas {
  //Constructors
  uuuu2::uuuu2(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), gs(gscharge), mu(massu), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propAG(0,0), proph(mh,wh), propZ(MZ,WZ),
    p1(particle(mu)), p2(particle(mu)),
    p3(particle(mu)), p4(particle(mu)),
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
    preh = e*e*mu*mu/(4.0*MW*MW*SW*SW);
    gL=1.0-4.0/3.0*SW*SW;
    gR=-4.0/3.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*mu*mu/MZ/MZ;//=preh!
  }
  void uuuu2::set_masses(const ldouble& massu, const ldouble& massh, const ldouble& massW){
    mu=massu;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(mu);
    p2.set_mass(mu);
    p3.set_mass(mu);
    p4.set_mass(mu);
    preh = e*e*mu*mu/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*mu*mu/MZ/MZ;//=preh
  }
  void uuuu2::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    ldouble propPU[4], propPT[4];
    for(int j=0;j<4;j++){
      propPU[j] = mom1[j]-mom4[j];
      propPT[j] = mom1[j]-mom3[j];
    }
    pDenUAG = propAG.den(propPU);
    pDenUh = proph.den(propPU);
    pDenUZ = propZ.den(propPU);
    pDenTAG = propAG.den(propPT);
    pDenTh = proph.den(propPT);
    pDenTZ = propZ.den(propPT);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  //We have to separate the gluon and the U and T channels because each has its own color factor.
  cdouble uuuu2::amp_gluon_U(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Gluon
    //uUUu
    //all in:
    // - (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //uuUU: 2<->4
    // + (<13>[24] + <12>[34] + [13]<24> + [12]<34>)
    //34 out:
    // + (- <13>[24] + <12>[34] - [13]<24> + [12]<34>)
    amplitude = two*gs*gs*(
			  - a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
			  - s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
			  )/pDenUAG;

    return amplitude;
  }
  cdouble uuuu2::amp_gluon_T(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Gluon
    //uUUu
    //all in:
    //  (<12>[34] - <14>[23] + [12]<34> - [14]<23>)
    //uuUU: 2<->4
    //  (- <14>[23] + <12>[34] - [14]<23> + [12]<34>)
    //34 outgoing:
    //  (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude = two*gs*gs*(
			  + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
			  + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
			  )/pDenTAG;

    return amplitude;
  }
  cdouble uuuu2::amp_rest_U(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //uUUu
    //all in: 
    // - (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //uuUU: 2<->4
    // + (<13>[24] + <12>[34] + [13]<24> + [12]<34>)
    //34 out:
    // + (- <13>[24] + <12>[34] - [13]<24> + [12]<34>)
    amplitude = two*e*e*4.0/9.0*(
				 - a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
				 - s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
				 )/pDenUAG;
    
    //Higgs
    //preh = e*e*me*mm/(4*MW*MW*SW*SW);
    //uUUu
    //all in:
    //preh ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //uuUU: 2<->4
    //- preh ([14]+<14>) ([23]+<23>)/(u-Mh^2)
    //34 out:
    //- preh ([14]-<14>) ([23]-<23>)/(u-Mh^2)    
    amplitude += - preh*(s14s.v(ds1,ds4)-a14a.v(ds1,ds4))*(s23s.v(ds2,ds3)-a23a.v(ds2,ds3))/pDenUh;
    
    //Z Boson
    //Defined above:
    //gL=1.0-4.0/3.0*SW*SW;
    //gR=-4.0/3.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gL-gR)*(gL-gR)*mu*mu/MZ/MZ; // = preh
    //uUUu:
    //all in:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- 2 preZ (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    //uuUU: 2<->4
    //+ preZ0 (<14>-[14]) (<23>-[23]) / (u-MZ^2)
    //+ 2 preZ (gL^2 [34] <12> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [12] <34>))/(u-MZ^2)
    //34 out:
    //+ preZ0 (<14>+[14]) (<23>+[23]) / (u-MZ^2)
    //+ 2 preZ (gL^2 [34] <12> - gLgR( [13] <24>+ [24] <13> ) + gR^2 [12] <34>))/(u-MZ^2)
    amplitude += 
      + preZ0*(a14a.v(ds1,ds4)+s14s.v(ds1,ds4))*(a23a.v(ds2,ds3)+s23s.v(ds2,ds3))/pDenUZ
      + two*preZ*(
	      gL*gL*s34s.v(ds3,ds4)*a12a.v(ds1,ds2)
	      - gL*gR*(s13s.v(ds1,ds3)*a24a.v(ds2,ds4)+s24s.v(ds2,ds4)*a13a.v(ds1,ds3))
	      + gR*gR*s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
	      )/pDenUZ;
    
    return amplitude;
  }
  cdouble uuuu2::amp_rest_T(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //uUUu
    //all in:
    //  (<12>[34] - <14>[23] + [12]<34> - [14]<23>)
    //uuUU: 2<->4
    //  (- <14>[23] + <12>[34] - [14]<23> + [12]<34>)
    //34 outgoing:
    //  (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude += two*e*e*4.0/9.0*(
				  + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
				  + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
				  )/pDenTAG;
    
    //Higgs
    //uUUu:
    //all in:
    //- preh * ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //uuUU: 2<->4
    //+ preh * ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //34 out:
    //+ preh * ([13]-<13>) ([24]-<24>)/(t-Mh^2)
    amplitude += preh*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3))*(s24s.v(ds2,ds4)-a24a.v(ds2,ds4))/pDenTh;
    
    //Z Boson
    //T Channel
    //uUUu:
    //all in:
    // + preZ0 (<13>-[13]) (<24>-[24]) / (t-MZ^2)
    // + 2 preZ (- gL^2 [23] <14> + gLgR( [12] <34>+ [34] <12> ) - gR^2 [14] <23>))/(t-MZ^2)
    //uuUU: 2<->4
    // - preZ0 (<13>-[13]) (<24>-[24]) / (t-MZ^2)
    // + 2 preZ (+ gL^2 [34] <12> - gLgR( [14] <23>+ [23] <14> ) + gR^2 [12] <34>))/(t-MZ^2)
    //34 out:
    // - preZ0 (<13>+[13]) (<24>+[24]) / (t-MZ^2)
    // + 2 preZ (+ gL^2 [34] <12> + gLgR( [14] <23>+ [23] <14> ) + gR^2 [12] <34>))/(t-MZ^2)
    amplitude += 
      - preZ0*(a13a.v(ds1,ds3)+s13s.v(ds1,ds3))*(a24a.v(ds2,ds4)+s24s.v(ds2,ds4))/pDenTZ
      + two*preZ*(
	      gL*gL*s34s.v(ds3,ds4)*a12a.v(ds1,ds2)
	      + gL*gR*(s14s.v(ds1,ds4)*a23a.v(ds2,ds3)+s23s.v(ds2,ds3)*a14a.v(ds1,ds4))
	      + gR*gR*s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
	      )/pDenTZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uuuu2::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble two=2, three = 3, four = 4, nine = 9;
    cdouble M_rest_U, M_gluon_U, M_rest_T, M_gluon_T;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M_rest_U = amp_rest_U(j1,j2,j3,j4);
	    M_rest_T = amp_rest_T(j1,j2,j3,j4);
	    M_gluon_U = amp_gluon_U(j1,j2,j3,j4);
	    M_gluon_T = amp_gluon_T(j1,j2,j3,j4);
	    amp2 += nine*std::pow(std::abs(M_rest_U),2);// Color factor Tr(1)^2 = 3*3=9
	    amp2 += two*std::pow(std::abs(M_gluon_U),2);// Color factor Tr(Ta,Tb)^2 = C^2*delta^ab*delta^ab = 1/4*8=2
	    amp2 += nine*std::pow(std::abs(M_rest_T),2);// Color factor Tr(1)^2 = 3*3=9
	    amp2 += two*std::pow(std::abs(M_gluon_T),2);// Color factor Tr(Ta,Tb)^2 = 2
	    //Cross terms
	    //No cross term (color factor Tr(Ta)^2=0) for rest_U*gluon_U and rest_T*gluon_T
	    amp2 += three*two*std::real(M_rest_U*std::conj(M_rest_T));// Color factor Tr(1) = 3 -- 2: 2real(w*conj(z)) = w*conj(z)+z*conj(w)
	    amp2 += (-two/three)*two*std::real(M_gluon_U*std::conj(M_gluon_T));// Color factor Tr(Ta,Tb,Ta,Tb) = -2/3 -- 2: 2real(w*conj(z)) = w*conj(z)+z*conj(w)
	    amp2 += four*two*std::real(M_gluon_U*std::conj(M_rest_T));// Color factor Tr(Ta,Ta) = 4 -- 2: 2real(w*conj(z)) = w*conj(z)+z*conj(w)
	    amp2 += four*two*std::real(M_gluon_T*std::conj(M_rest_U));// Color factor Tr(Ta,Ta) = 4 -- 2: 2real(w*conj(z)) = w*conj(z)+z*conj(w)
	    
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over colors 1/9
    //Symmetry Factor 1/2
    return amp2/72.0;
  }

  



  //  Tests
  int test_uuuu2(){
    int n=0;//Number of fails
    std::cout<<"\t* u , u  -> u , u       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mu=0.0042, pspatial=250\n";
      ldouble mu=0.0042, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      uuuu2 uuuu2Amp = uuuu2(0.31333,1.238,mu,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.631331854019200E+03,1.719846893556283E+02,5.894911559839137E+01,2.881324077060831E+01,1.686793663070660E+01,1.110006095756243E+01,7.996329285396084E+00,6.245304817999785E+00,5.282769007490297E+00,4.852534582451461E+00,4.852534582451461E+00,5.282769007490295E+00,6.245304817999781E+00,7.996329285396081E+00,1.110006095756243E+01,1.686793663070659E+01,2.881324077060829E+01,5.894911559839132E+01,1.719846893556283E+02,1.631331854019193E+03};
      i += uuuu2Amp.test_2to2_amp2([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH);
      i += uuuu2Amp.test_2to2_amp2_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH);
      i += uuuu2Amp.test_2to2_amp2_boosts([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH);
      i += uuuu2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH);
      //std::cout<<"\n########### mu=0.0042, pspatial=0.005\n";
      pspatial = 0.005;
      ldouble dataCH2[20] = {2.989654377010874E+03,3.155520283014097E+02,1.081126447922916E+02,5.271912189294680E+01,3.072527313233418E+01,2.008783687011044E+01,1.435500052705444E+01,1.111646372502473E+01,9.334437690583947E+00,8.537396599824246E+00,8.537396599824246E+00,9.334437690583941E+00,1.111646372502472E+01,1.435500052705443E+01,2.008783687011044E+01,3.072527313233418E+01,5.271912189294675E+01,1.081126447922916E+02,3.155520283014093E+02,2.989654377010861E+03};
      i += uuuu2Amp.test_2to2_amp2([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH2);
      i += uuuu2Amp.test_2to2_amp2_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH2);
      i += uuuu2Amp.test_2to2_amp2_boosts([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH2);
      i += uuuu2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH2);
      //std::cout<<"\n########### mu=1.2, pspatial=0.3\n";
      mu=1.2;
      pspatial = 0.3;
      uuuu2Amp.set_masses(mu,mh,MW);
      ldouble dataCH4[20] = {1.357442054688230E+05,1.513880986687428E+04,5.512113495054864E+03,2.873134238543552E+03,1.799013129190489E+03,1.267450756647455E+03,9.750856795930729E+02,8.071025404467632E+02,7.134806893238718E+02,6.712725779870998E+02,6.712725779870998E+02,7.134806893238714E+02,8.071025404467628E+02,9.750856795930727E+02,1.267450756647454E+03,1.799013129190488E+03,2.873134238543550E+03,5.512113495054860E+03,1.513880986687427E+04,1.357442054688224E+05};
      i += uuuu2Amp.test_2to2_amp2([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH4);
      i += uuuu2Amp.test_2to2_amp2_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH4);
      i += uuuu2Amp.test_2to2_amp2_boosts([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH4);
      i += uuuu2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH4);
      //std::cout<<"\n########### mu=1.2, MW=2.11, pspatial=0.3\n";
      mu=1.2;
      MW=2.11;
      pspatial = 0.3;
      uuuu2Amp.set_masses(mu,mh,MW);
      ldouble dataCH5[20] = {1.357703831145274E+05,1.514798296105951E+04,5.517918066174580E+03,2.877522497435610E+03,1.802639929039143E+03,1.270618054821265E+03,9.779607815161086E+02,8.097912209904414E+02,7.160578030303037E+02,6.737972636455107E+02,6.737972636455107E+02,7.160578030303034E+02,8.097912209904409E+02,9.779607815161083E+02,1.270618054821264E+03,1.802639929039143E+03,2.877522497435608E+03,5.517918066174576E+03,1.514798296105950E+04,1.357703831145269E+05};
      i += uuuu2Amp.test_2to2_amp2([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH5);
      i += uuuu2Amp.test_2to2_amp2_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH5);
      i += uuuu2Amp.test_2to2_amp2_boosts([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH5);
      i += uuuu2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH5);
      //std::cout<<"\n########### mu=1.2, mu=1.23, MW=0.006, pspatial=0.3\n";
      mu=1.2;
      MW=0.006;
      pspatial = 0.3;
      uuuu2Amp.set_masses(mu,mh,MW);
      ldouble dataCH6[20] = {2.352861705823458E+07,2.267805904847016E+07,2.252390755396980E+07,2.246050996586506E+07,2.242677176481464E+07,2.240652974239203E+07,2.239370379020315E+07,2.238553947759005E+07,2.238066036556913E+07,2.237836934644146E+07,2.237836934644146E+07,2.238066036556912E+07,2.238553947759004E+07,2.239370379020316E+07,2.240652974239204E+07,2.242677176481464E+07,2.246050996586507E+07,2.252390755396980E+07,2.267805904847015E+07,2.352861705823452E+07};
      i += uuuu2Amp.test_2to2_amp2([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH6);
      i += uuuu2Amp.test_2to2_amp2_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH6);
      i += uuuu2Amp.test_2to2_amp2_boosts([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH6);
      i += uuuu2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH6);
      //std::cout<<"\n########### mu=1.2, MW=2.11, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      uuuu2Amp.set_masses(mu,mh,MW);
      ldouble dataCH7[20] = {1.357751325895691E+05,1.514961163790949E+04,5.518928528476250E+03,2.878273029037016E+03,1.803250713208439E+03,1.271144508801571E+03,9.784336102394273E+02,8.102298368510983E+02,7.164759436973450E+02,6.742057825179271E+02,6.742057825179271E+02,7.164759436973446E+02,8.102298368510978E+02,9.784336102394270E+02,1.271144508801570E+03,1.803250713208438E+03,2.878273029037015E+03,5.518928528476246E+03,1.514961163790948E+04,1.357751325895685E+05};
      i += uuuu2Amp.test_2to2_amp2([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH7);
      i += uuuu2Amp.test_2to2_amp2_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH7);
      i += uuuu2Amp.test_2to2_amp2_boosts([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH7);
      i += uuuu2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH7);
      //std::cout<<"\n########### mu=1.2, MW=0.006, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      MW=0.006;
      pspatial = 0.3;
      uuuu2Amp.set_masses(mu,mh,MW);
      ldouble dataCH8[20] = {2.596237799811850E+07,2.474600953338573E+07,2.451972816655959E+07,2.442619165794419E+07,2.437636502882333E+07,2.434648546446338E+07,2.432757322771700E+07,2.431554890116436E+07,2.430837022264750E+07,2.430500167216714E+07,2.430500167216715E+07,2.430837022264749E+07,2.431554890116435E+07,2.432757322771700E+07,2.434648546446338E+07,2.437636502882333E+07,2.442619165794419E+07,2.451972816655959E+07,2.474600953338572E+07,2.596237799811845E+07};
      i += uuuu2Amp.test_2to2_amp2([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH8);
      i += uuuu2Amp.test_2to2_amp2_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH8);
      i += uuuu2Amp.test_2to2_amp2_boosts([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH8);
      i += uuuu2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuu2Amp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
