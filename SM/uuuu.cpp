
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

//File:  SPINAS/SM/uuuu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uuuu.h"

namespace spinas {
  //Constructors
  uuuu::uuuu(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
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
  void uuuu::set_masses(const ldouble& massu, const ldouble& massh, const ldouble& massW){
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
  void uuuu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
  cdouble uuuu::amp_gluon_S(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
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
  cdouble uuuu::amp_gluon_T(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
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
  cdouble uuuu::amp_rest_S(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //Sign changes due to p3 and p4 being outgoing.
    // (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    amplitude = two*e*e*4.0/9.0*(
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
    //preZ0 = preZ*(gL-gR)*(gL-gR)*mu*mu/MZ/MZ; // = preh
    //all in:
    //+(EE^2 Mu Mu (gL-gR)^2 (<12>-[12]) (<34>-[34]))/(8 CW^2 MZ^2 SW^2 (s-MZ^2))
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
  cdouble uuuu::amp_rest_T(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //T Channel
    // 2 <-> 3 with a minus sign for fermions
    // all ingoing:  (<12>[34] - <14>[23] + [12]<34> - [14]<23>)
    // 34 outgoing:  (<12>[34] + <14>[23] + [12]<34> + [14]<23>)/pDenT
    amplitude += two*e*e*4.0/9.0*(
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
  ldouble uuuu::amp2(){
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
  int test_uuuu(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> u , U       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mu=0.0042, pspatial=250\n";
      ldouble mu=0.0042, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      uuuu uuuuAmp = uuuu(0.31333,1.238,mu,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.285162064580300E+03,3.501286817780895E+02,1.209697791284904E+02,5.929591088812482E+01,3.450399746419477E+01,2.224968893233480E+01,1.537295347790775E+01,1.116927590358787E+01,8.438331883081862E+00,6.583337365131206E+00,5.281379292315666E+00,4.345877979344921E+00,3.663218708994452E+00,3.161181906287536E+00,2.792187504831790E+00,2.523923719686284E+00,2.333869509029420E+00,2.205970010713291E+00,2.128551920666219E+00,2.092978646667690E+00};
      i += uuuuAmp.test_2to2_amp2([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH);
      i += uuuuAmp.test_2to2_amp2_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH);
      i += uuuuAmp.test_2to2_amp2_boosts([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH);
      i += uuuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH);
      //std::cout<<"\n########### mu=0.0042, pspatial=0.005\n";
      pspatial = 0.005;
      ldouble dataCH2[20] = {6.024234006489967E+03,6.453045831291670E+02,2.240333163251409E+02,1.102716514584603E+02,6.438456929037260E+01,4.162395202411605E+01,2.880281521629490E+01,2.092967259820590E+01,1.578431233732639E+01,1.226034580316262E+01,9.757696605322916E+00,7.928789071390195E+00,6.561437473487648E+00,5.520414317212088E+00,4.716389993815790E+00,4.088533683597716E+00,3.594327078318527E+00,3.203371693675780E+00,2.893499308620652E+00,2.648258407613977E+00};
      i += uuuuAmp.test_2to2_amp2([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH2);
      i += uuuuAmp.test_2to2_amp2_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH2);
      i += uuuuAmp.test_2to2_amp2_boosts([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH2);
      i += uuuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH2);
      //std::cout<<"\n########### mu=1.2, pspatial=0.3\n";
      mu=1.2;
      pspatial = 0.3;
      uuuuAmp.set_masses(mu,mh,MW);
      ldouble dataCH4[20] = {2.704093774429977E+05,2.978389127063369E+04,1.062855552046007E+04,5.375224573967216E+03,3.223103684797658E+03,2.138595879988551E+03,1.517645830349706E+03,1.129811395266634E+03,8.717868268216185E+02,6.916860658555493E+02,5.611438402218322E+02,4.635987242329334E+02,3.888585394145855E+02,3.303745733042836E+02,2.837850197739428E+02,2.460947430219357E+02,2.151925186027655E+02,1.895559731347899E+02,1.680652623058193E+02,1.498820108979344E+02};
      i += uuuuAmp.test_2to2_amp2([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH4);
      i += uuuuAmp.test_2to2_amp2_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH4);
      i += uuuuAmp.test_2to2_amp2_boosts([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH4);
      i += uuuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH4);
      //std::cout<<"\n########### mu=1.2, MW=2.11, pspatial=0.3\n";
      mu=1.2;
      MW=2.11;
      pspatial = 0.3;
      uuuuAmp.set_masses(mu,mh,MW);
      ldouble dataCH5[20] = {2.702255339144453E+05,2.972405703445178E+04,1.059351659734705E+04,5.350807558976173E+03,3.204584125742109E+03,2.123826311842330E+03,1.505469937030984E+03,1.119535396898797E+03,8.629618156638037E+02,6.840049022693900E+02,5.543871457965738E+02,4.576043561602148E+02,3.835032751130502E+02,3.255625785163225E+02,2.794402931805873E+02,2.421559980785551E+02,2.116096163334628E+02,1.862873741282294E+02,1.650761673606742E+02,1.471429784257655E+02};
      i += uuuuAmp.test_2to2_amp2([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH5);
      i += uuuuAmp.test_2to2_amp2_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH5);
      i += uuuuAmp.test_2to2_amp2_boosts([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH5);
      i += uuuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH5);
      //std::cout<<"\n########### mu=1.2, mu=1.23, MW=0.006, pspatial=0.3\n";
      mu=1.2;
      MW=0.006;
      pspatial = 0.3;
      uuuuAmp.set_masses(mu,mh,MW);
      ldouble dataCH6[20] = {4.269398026575736E+07,4.386892625528717E+07,4.413298470491174E+07,4.424910675018226E+07,4.431437749838014E+07,4.435619056754918E+07,4.438526080351583E+07,4.440664013253763E+07,4.442302187485366E+07,4.443597327674414E+07,4.444646830864288E+07,4.445514411442180E+07,4.446243502597629E+07,4.446864731043839E+07,4.447400312941124E+07,4.447866755452256E+07,4.448276579851384E+07,4.448639456080157E+07,4.448962970649454E+07,4.449253159060198E+07};
      i += uuuuAmp.test_2to2_amp2([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH6);
      i += uuuuAmp.test_2to2_amp2_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH6);
      i += uuuuAmp.test_2to2_amp2_boosts([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH6);
      i += uuuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH6);
      //std::cout<<"\n########### mu=1.2, MW=2.11, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      uuuuAmp.set_masses(mu,mh,MW);
      ldouble dataCH7[20] = {2.702284639527313E+05,2.972494982627628E+04,1.059400200476158E+04,5.351118419111275E+03,3.204798050738656E+03,2.123978579754362E+03,1.505579543572516E+03,1.119613739359406E+03,8.630162686385811E+02,6.840405114257449E+02,5.544075150854190E+02,4.576121492831510E+02,3.835005164026732E+02,3.255508423793764E+02,2.794208282158887E+02,2.421298111568743E+02,2.115775311986503E+02,1.862500732551606E+02,1.650342225002335E+02,1.470968733177810E+02};
      i += uuuuAmp.test_2to2_amp2([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH7);
      i += uuuuAmp.test_2to2_amp2_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH7);
      i += uuuuAmp.test_2to2_amp2_boosts([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH7);
      i += uuuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH7);
      //std::cout<<"\n########### mu=1.2, MW=0.006, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      MW=0.006;
      pspatial = 0.3;
      uuuuAmp.set_masses(mu,mh,MW);
      ldouble dataCH8[20] = {4.714030326683572E+07,4.802270433607113E+07,4.820710622696318E+07,4.827409265449134E+07,4.830048011194506E+07,4.830813349797463E+07,4.830564069567200E+07,4.829706063287691E+07,4.828456057225213E+07,4.826940289176049E+07,4.825237255880871E+07,4.823398320651425E+07,4.821458495402800E+07,4.819442460123484E+07,4.817368105712346E+07,4.815248712669383E+07,4.813094341468939E+07,4.810912748458108E+07,4.808710006014419E+07,4.806490932680301E+07};
      i += uuuuAmp.test_2to2_amp2([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH8);
      i += uuuuAmp.test_2to2_amp2_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH8);
      i += uuuuAmp.test_2to2_amp2_boosts([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH8);
      i += uuuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuuuAmp.amp2(); }, mu,mu,mu,mu,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
