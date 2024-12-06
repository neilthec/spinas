
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

//File:  SPINAS/SM/usus.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/usus.h"

namespace spinas {
  //Constructors
  usus::usus(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& masss, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), gs(gscharge), mu(massu), ms(masss), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    propAG(0,0), proph(mh,wh),
    p1(particle(mu)), p2(particle(ms)),
    p3(particle(mu)), p4(particle(ms)),
    a13a(sproduct(ANGLE,&p1,&p3,2)),
    s13s(sproduct(SQUARE,&p1,&p3,2)),
    a14a(sproduct(ANGLE,&p1,&p4,2)),
    s14s(sproduct(SQUARE,&p1,&p4,2)),
    a23a(sproduct(ANGLE,&p2,&p3,2)),
    s23s(sproduct(SQUARE,&p2,&p3,2)),
    a24a(sproduct(ANGLE,&p2,&p4,2)),
    s24s(sproduct(SQUARE,&p2,&p4,2)),
    s12s(sproduct(SQUARE,&p1,&p2,2)),
    a12a(sproduct(ANGLE,&p1,&p2,2)),
    s34s(sproduct(SQUARE,&p3,&p4,2)),
    a34a(sproduct(ANGLE,&p3,&p4,2))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    preh = e*e*mu*ms/(4.0*MW*MW*SW*SW);
    gLu=1.0-4.0/3.0*SW*SW;
    gRu=-4.0/3.0*SW*SW;
    gLd=-1.0+2.0/3.0*SW*SW;
    gRd=2.0/3.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*ms/MZ/MZ;//=-preh!
  }
  void usus::set_masses(const ldouble& massu, const ldouble& masss, const ldouble& massh, const ldouble& massW){
    mu=massu;
    ms=masss;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(mu);
    p2.set_mass(ms);
    p3.set_mass(mu);
    p4.set_mass(ms);
    preh = e*e*mu*ms/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*ms/MZ/MZ;//=-preh
  }
  void usus::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
      propP[j] = mom1[j]-mom3[j];
    pDenTAG = propAG.denominator(propP);
    pDenTh = proph.denominator(propP);
    pDenTZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  //We have to separate the gluon so we can separate the color factor between the gluon^2 diagram
  // And the rest^2.
  // The cross term vanishes due to the trace of the adjoint rep.
  cdouble usus::amp_gluon(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Gluon
    //uUSs:
    //all in:
    // - (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //usUS: 4->2->3->4
    // - (- <14>[23] + <12>[34] - [14]<23> + [12]<34>)
    //34 out:
    // - (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude = - two*gs*gs*(
			     + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
			     + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
			     )/pDenTAG;

    return amplitude;
  }
  cdouble usus::amp_rest(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //uUSs:
    //all in:
    // - (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //usUS: 4->2->3->4
    // - (- <14>[23] + <12>[34] - [14]<23> + [12]<34>)
    //34 out:
    // - (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude = two*e*e*2.0/9.0*(
				  + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
				  + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
				  )/pDenTAG;
    
    //Higgs
    //preh = e*e*mu*ms/(4*MW*MW*SW*SW);
    //uUSs:
    //all in:
    //preh ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //usUS: 4->2->3->4
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
    //preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*ms/MZ/MZ; // = -preh
    //uUSs:
    //all in:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ ( gLu gLd [23] <14> + gLd gRu [13] <24> + gLu gRd [24] <13> + gRu gRd [14] <23> )/(s-MZ^2)
    //usUS: 4->2->3->4
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
  
  //set_momenta(...) must be called before amp2().
  ldouble usus::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble two=2, three = 3, nine = 9;
    cdouble M_rest, M_gluon;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M_rest = amp_rest(j1,j2,j3,j4);
	    M_gluon = amp_gluon(j1,j2,j3,j4);
	    amp2 += nine*std::pow(std::abs(M_rest),2);// Color factor 3*3=9
	    //Cross term color factor 0 (Trace(Ta)*Trace(Ta)=0)
	    amp2 += two*std::pow(std::abs(M_gluon),2);//Color factor C^2*delta^ab*delta^ab = 1/4*8=2
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over colors 1/3*1/3 = 1/9
    return amp2/36.0;
  }

  



  //  Tests
  int test_usus(){
    int n=0;//Number of fails
    std::cout<<"\t* u , s  -> u , s       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mu=0.0042, ms=1.23, pspatial=250\n";
      ldouble mu=0.0042, ms=1.23, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      usus ususAmp = usus(0.31333,1.238,mu,ms,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.265090964917122E+03,3.456866621254358E+02,1.185226888484794E+02,5.759849056453773E+01,3.320343519529108E+01,2.119614560843673E+01,1.448634139000884E+01,1.039958049003045E+01,7.750540684092007E+00,5.950728004042229E+00,4.682171855313013E+00,3.761417235251554E+00,3.076919358275862E+00,2.557883100633553E+00,2.157740176341411E+00,1.844916927633297E+00,1.597443254265064E+00,1.399683735470590E+00,1.240289962539918E+00,1.110880866523161E+00};
      i += ususAmp.test_2to2_amp2([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH);
      i += ususAmp.test_2to2_amp2_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH);
      i += ususAmp.test_2to2_amp2_boosts([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH);
      i += ususAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH);
      //std::cout<<"########### mu=0.0042, ms=1.23, pspatial=1.25\n";
      pspatial = 1.25;
      ldouble dataCH2[20] = {4.707362289895174E+03,4.971525489313572E+02,1.699897289776446E+02,8.231577836894006E+01,4.722916895431791E+01,2.996704014969722E+01,2.032438258299762E+01,1.445328479354512E+01,1.064884817365128E+01,8.064833200688373E+00,6.244108993984576E+00,4.922984806978280E+00,3.941160756947429E+00,3.196915768121126E+00,2.623349972066220E+00,2.175112035116045E+00,1.820649780247469E+00,1.537511966816183E+00,1.309405631659035E+00,1.124299842128789E+00};
      i += ususAmp.test_2to2_amp2([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH2);
      i += ususAmp.test_2to2_amp2_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH2);
      i += ususAmp.test_2to2_amp2_boosts([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH2);
      i += ususAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH2);
      //std::cout<<"########### mu=1.23, ms=0.0042, pspatial=0.005\n";
      mu=1.23;
      ms=0.0042;
      pspatial = 0.005;
      ususAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH3[20] = {8.554579573594129E+07,9.221042636026561E+06,3.217319202885678E+06,1.589318056727537E+06,9.298789857624845E+05,6.013538782932088E+05,4.154286926615356E+05,3.006714986921138E+05,2.252412931185346E+05,1.732365334526094E+05,1.360137229050625E+05,1.085552426946642E+05,8.779104532708366E+04,7.176007082713548E+04,5.916383863988237E+04,4.911601802175328E+04,4.099559844010428E+04,3.435749260673205E+04,2.887630280849144E+04,2.430992571072147E+04};
      i += ususAmp.test_2to2_amp2([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH3);
      i += ususAmp.test_2to2_amp2_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH3);
      i += ususAmp.test_2to2_amp2_boosts([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH3);
      i += ususAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH3);
      //std::cout<<"########### mu=1.2, ms=1.23, pspatial=0.3\n";
      mu=1.2;
      ms=1.23;
      pspatial = 0.3;
      ususAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH4[20] = {2.815779524031042E+05,3.096408080369855E+04,1.103135446718514E+04,5.569374606979876E+03,3.333619869386775E+03,2.207896055491394E+03,1.563881630912950E+03,1.161977458443131E+03,8.948173712477198E+02,7.084982605175460E+02,5.735643614927866E+02,4.728246775140469E+02,3.957039760104885E+02,3.354102540782957E+02,2.874215517888030E+02,2.486340629556271E+02,2.168608848125675E+02,1.905257406006970E+02,1.684696278489018E+02,1.498252257913540E+02};
      i += ususAmp.test_2to2_amp2([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH4);
      i += ususAmp.test_2to2_amp2_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH4);
      i += ususAmp.test_2to2_amp2_boosts([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH4);
      i += ususAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH4);
      //std::cout<<"########### mu=1.2, ms=1.23, MW=2.11, pspatial=0.3\n";
      mu=1.2;
      ms=1.23;
      MW=2.11;
      pspatial = 0.3;
      ususAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH5[20] = {2.815783246905214E+05,3.096420866922901E+04,1.103143342540494E+04,5.569432588763860E+03,3.333666185195739E+03,2.207934937489164E+03,1.563915358070353E+03,1.162007398239350E+03,8.948444086061735E+02,7.085230009860622E+02,5.735872376713361E+02,4.728460092095651E+02,3.957240063297860E+02,3.354291721532391E+02,2.874395076652390E+02,2.486511776794415E+02,2.168772574593731E+02,1.905414532970501E+02,1.684847494364429E+02,1.498398145514343E+02};
      i += ususAmp.test_2to2_amp2([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH5);
      i += ususAmp.test_2to2_amp2_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH5);
      i += ususAmp.test_2to2_amp2_boosts([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH5);
      i += ususAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH5);
      //std::cout<<"########### mu=1.2, ms=1.23, MW=0.006, pspatial=0.3\n";
      mu=1.2;
      ms=1.23;
      MW=0.006;
      pspatial = 0.3;
      ususAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH6[20] = {2.034375311590081E+07,2.009169324188472E+07,2.007163769897423E+07,2.006614020939152E+07,2.006388918178958E+07,2.006275540088154E+07,2.006210656874272E+07,2.006170152901981E+07,2.006143219965163E+07,2.006124430908073E+07,2.006110819527241E+07,2.006100654396894E+07,2.006092870233887E+07,2.006086782739885E+07,2.006081936248441E+07,2.006078017944296E+07,2.006074807384570E+07,2.006072145651417E+07,2.006069915876072E+07,2.006068030586727E+07};
      i += ususAmp.test_2to2_amp2([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH6);
      i += ususAmp.test_2to2_amp2_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH6);
      i += ususAmp.test_2to2_amp2_boosts([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH6);
      i += ususAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH6);
      //std::cout<<"########### mu=1.2, ms=1.23, MW=2.11, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      ms=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      ususAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH7[20] = {2.815790267982625E+05,3.096444228945029E+04,1.103157334837095E+04,5.569532355956384E+03,3.333743643787115E+03,2.207998199839874E+03,1.563968792382313E+03,1.162053625440340E+03,8.948851245963727E+02,7.085593661232518E+02,5.736200807869576E+02,4.728759429221481E+02,3.957514962330052E+02,3.354545803787665E+02,2.874631214239247E+02,2.486732285956621E+02,2.168979350464882E+02,1.905609145827292E+02,1.685031259865534E+02,1.498572176929332E+02};
      i += ususAmp.test_2to2_amp2([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH7);
      i += ususAmp.test_2to2_amp2_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH7);
      i += ususAmp.test_2to2_amp2_boosts([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH7);
      i += ususAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH7);
      //std::cout<<"########### mu=1.2, ms=1.23, MW=0.006, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      ms=1.23;
      MW=0.006;
      pspatial = 0.3;
      ususAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH8[20] = {2.781061532582761E+07,2.749218163872934E+07,2.747290556130367E+07,2.747775662827629E+07,2.748901602358511E+07,2.750280651855364E+07,2.751782344652037E+07,2.753351190728379E+07,2.754959776387193E+07,2.756593065424604E+07,2.758242147641550E+07,2.759901419273806E+07,2.761567185741873E+07,2.763236916129292E+07,2.764908820950472E+07,2.766581600954505E+07,2.768254291378829E+07,2.769926161924270E+07,2.771596650536307E+07,2.773265318391837E+07};
      i += ususAmp.test_2to2_amp2([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH8);
      i += ususAmp.test_2to2_amp2_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH8);
      i += ususAmp.test_2to2_amp2_boosts([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH8);
      i += ususAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ususAmp.amp2(); }, mu,ms,mu,ms,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
