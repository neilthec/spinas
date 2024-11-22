
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

//File:  SPINAS/SM/emem.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/emem.h"

namespace spinas {
  //Constructors
  emem::emem(const ldouble& echarge, const ldouble& masse, const ldouble& massmu, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), mm(massmu), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    propA(0,0), proph(mh,wh),
    p1(particle(me)), p2(particle(mm)),
    p3(particle(me)), p4(particle(mm)),
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
    preh = e*e*me*mm/(4.0*MW*MW*SW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*me*mm/MZ/MZ;//=preh!
    preZLL = preZ*2.0*gL*gL;
    preZLR = preZ*2.0*gL*gR;
    preZRR = preZ*2.0*gR*gR;
  }
  void emem::set_masses(const ldouble& masse, const ldouble& massmu, const ldouble& massh, const ldouble& massW){
    me=masse;
    mm=massmu;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(me);
    p2.set_mass(mm);
    p3.set_mass(me);
    p4.set_mass(mm);
    preh = e*e*me*mm/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*me*mm/MZ/MZ;
  }
  void emem::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenTA = propA.denominator(propP);
    pDenTh = proph.denominator(propP);
    pDenTZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble emem::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble one=1, two = 2, four=4, eight=8;
    cdouble amplitude(0,0);
    
    //Photon
    //eEMm all in:
    //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //emEM: 4->2->3->4
    //- (- <14>[23] + <12>[34] - [14]<23> + [12]<34>)
    //34 out:
    //- (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude += - two*e*e*(
			    + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
			    + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
			    )/pDenTA;
    
    //Higgs
    //preh = e*e*me*mm/(4*MW*MW*SW*SW);
    //eEMm all in:
    //preh ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //emEM: 4->2->3->4
    //- preh ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //34 out:
    //- preh ([13]-<13>) ([24]-<24>)/(t-Mh^2)
    amplitude += - preh*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3))*(s24s.v(ds2,ds4)-a24a.v(ds2,ds4))/pDenTh;
    
    //Z Boson
    //Defined above:
    //gL=2.0*SW*SW-1.0;
    //gR=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gL-gR)*(gL-gR)*me*mm/MZ/MZ; // = preh
    //preZLL = preZ*2.0*gL*gL;
    //preZLR = preZ*2.0*gL*gR;
    //preZRR = preZ*2.0*gR*gR;
    //eEMm all in:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    //emEM: 4->2->3->4
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
  
  //set_momenta(...) must be called before amp2().
  ldouble emem::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    //M = amp(j1,j2,j3,j4);
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    return amp2/4.0;
  }

  



  //  Tests
  int test_emem(){
    int n=0;//Number of fails
    std::cout<<"\t* e , m  -> e , m       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### me=0.0005, mmu=0.105, pspatial=250\n";
      ldouble me=0.0005, mmu=0.105, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      emem ememAmp = emem(0.31333,me,mmu,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {6.227461063835386E+01,7.039128047362604E+00,2.534081160434054E+00,1.284402725315226E+00,7.697889557731630E-01,5.099764561510306E-01,3.612330517626481E-01,2.684599244914039E-01,2.068794462265394E-01,1.640227450376409E-01,1.330658803385414E-01,1.100228075596124E-01,9.244053731321415E-02,7.874363009235799E-02,6.788352225076523E-02,5.914090193990398E-02,5.200934395631676E-02,4.612417464004259E-02,4.121755342001841E-02,3.708929385759174E-02};
      i += ememAmp.test_2to2_amp2([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH);
      i += ememAmp.test_2to2_amp2_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH);
      i += ememAmp.test_2to2_amp2_boosts([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH);
      i += ememAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH);
      //std::cout<<"########### me=0.0005, mmu=0.105, pspatial=0.11\n";
      pspatial = 0.11;
      ldouble dataCH2[20] = {8.536542023055135E+01,9.015880770871000E+00,3.082986018745986E+00,1.493073280154556E+00,8.567993533234547E-01,5.437613309761956E-01,3.688976977451064E-01,2.624289650497791E-01,1.934378533029140E-01,1.465783330430451E-01,1.135606188498470E-01,8.960284227195248E-02,7.179806412859423E-02,5.830163562811980E-02,4.790036901788118E-02,3.977184588801153E-02,3.334388699771718E-02,2.820935295921192E-02,2.407278093976095E-02,2.071599739618701E-02};
      i += ememAmp.test_2to2_amp2([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH2);
      i += ememAmp.test_2to2_amp2_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH2);
      i += ememAmp.test_2to2_amp2_boosts([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH2);
      i += ememAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH2);
      //std::cout<<"########### me=0.105, mmu=0.0005, pspatial=0.005\n";
      me=0.105;
      mmu=0.0005;
      pspatial = 0.005;
      ememAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH3[20] = {7.364510816177090E+03,7.767207156357254E+02,2.646646512298400E+02,1.274061293112309E+02,7.246097045968226E+01,4.542091812973273E+01,3.031170293278287E+01,2.110923872472001E+01,1.514406175073547E+01,1.109097837767800E+01,8.234028930157422E+00,6.160159806305013E+00,4.618234451758954E+00,3.448865216243085E+00,2.547209474923595E+00,1.842182944492978E+00,1.284322007697004E+00,8.384242763946134E-01,4.789386605271361E-01,1.869944309116256E-01};
      i += ememAmp.test_2to2_amp2([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH3);
      i += ememAmp.test_2to2_amp2_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH3);
      i += ememAmp.test_2to2_amp2_boosts([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH3);
      i += ememAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH3);
      //std::cout<<"########### me=0.10, mmu=0.105, pspatial=0.05\n";
      me=0.10;
      mmu=0.105;
      pspatial = 0.05;
      ememAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH4[20] = {5.849835524942619E+02,6.323224092108515E+01,2.213404001276569E+01,1.097481731760567E+01,6.448584111119253E+00,4.190568542259298E+00,2.910866806988827E+00,2.119859570151809E+00,1.599150630805529E+00,1.239597998385934E+00,9.818354616961827E-01,7.913765895182775E-01,6.471059346225061E-01,5.355259512208452E-01,4.476930365496384E-01,3.774979032404815E-01,3.206568799337504E-01,2.740974819753368E-01,2.355717344851645E-01,2.034057919127872E-01};
      i += ememAmp.test_2to2_amp2([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH4);
      i += ememAmp.test_2to2_amp2_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH4);
      i += ememAmp.test_2to2_amp2_boosts([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH4);
      i += ememAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH4);
      //std::cout<<"########### me=0.10, mmu=0.105, MW=0.11, pspatial=0.05\n";
      me=0.10;
      mmu=0.105;
      MW=0.11;
      pspatial = 0.05;
      ememAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH5[20] = {5.851421621203569E+02,6.334046872824040E+01,2.222946145512348E+01,1.106301114540145E+01,6.531551937040927E+00,4.269321272349898E+00,2.986023474890327E+00,2.191848352804633E+00,1.668292646554371E+00,1.306147766223919E+00,1.046002673857736E+00,8.533388627596810E-01,7.070168067680098E-01,5.935201875357728E-01,5.038903215472664E-01,4.320055136401404E-01,3.735716862546871E-01,3.255074907298343E-01,2.855573054353559E-01,2.520406187058514E-01};
      i += ememAmp.test_2to2_amp2([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH5);
      i += ememAmp.test_2to2_amp2_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH5);
      i += ememAmp.test_2to2_amp2_boosts([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH5);
      i += ememAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH5);
      //std::cout<<"########### me=0.10, mmu=0.105, MW=0.006, pspatial=0.05\n";
      me=0.10;
      mmu=0.105;
      MW=0.006;
      pspatial = 0.05;
      ememAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH6[20] = {1.686366628124198E+03,1.090928976714863E+03,1.042271017081930E+03,1.028859211403645E+03,1.023358964618994E+03,1.020588475285801E+03,1.019004086534166E+03,1.018016121410862E+03,1.017360093631845E+03,1.016903157795580E+03,1.016572709350690E+03,1.016326378019481E+03,1.016138105345759E+03,1.015991161026698E+03,1.015874411349312E+03,1.015780218244371E+03,1.015703203122808E+03,1.015639491625625E+03,1.015586236918553E+03,1.015541310115033E+03};
      i += ememAmp.test_2to2_amp2([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH6);
      i += ememAmp.test_2to2_amp2_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH6);
      i += ememAmp.test_2to2_amp2_boosts([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH6);
      i += ememAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH6);
      //std::cout<<"########### me=0.10, mmu=0.105, MW=0.11, Mh=0.125, pspatial=0.05\n";
      me=0.10;
      mmu=0.105;
      MW=0.11;
      mh=0.125;
      pspatial = 0.05;
      ememAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH7[20] = {5.730397516380955E+02,5.950159005842288E+01,2.003605193701078E+01,9.570109943850820E+00,5.424498808974794E+00,3.405325273582793E+00,2.288374775062998E+00,1.614658156665682E+00,1.181976305164199E+00,8.905498439507217E-01,6.867969363079703E-01,5.399943776486121E-01,4.315783058404629E-01,3.498407045768514E-01,2.871288864380274E-01,2.382872438227948E-01,1.997511772034313E-01,1.690000527054805E-01,1.442160133247172E-01,1.240650550505993E-01};
      i += ememAmp.test_2to2_amp2([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH7);
      i += ememAmp.test_2to2_amp2_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH7);
      i += ememAmp.test_2to2_amp2_boosts([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH7);
      i += ememAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH7);
      //std::cout<<"########### me=0.10, mmu=0.105, MW=0.006, Mh=0.125, pspatial=0.05\n";
      me=0.10;
      mmu=0.105;
      MW=0.006;
      pspatial = 0.05;
      ememAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH8[20] = {4.775609210787762E+03,6.695243805649394E+03,6.951737319920967E+03,6.943413154987129E+03,6.853563000692598E+03,6.733305767513593E+03,6.601931268384440E+03,6.467968019963037E+03,6.335524423291958E+03,6.206651333407994E+03,6.082357440457563E+03,5.963092531248137E+03,5.848995559140284E+03,5.740029465419238E+03,5.636057686397045E+03,5.536888983726409E+03,5.442304310539764E+03,5.352073129631125E+03,5.265963363125084E+03,5.183747411755379E+03};
      i += ememAmp.test_2to2_amp2([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH8);
      i += ememAmp.test_2to2_amp2_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH8);
      i += ememAmp.test_2to2_amp2_boosts([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH8);
      i += ememAmp.test_2to2_amp2_boosts_and_rotations([&]() { return ememAmp.amp2(); }, me,mmu,me,mmu,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
