
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

//File:  SPINAS/SM/eeWW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eeWW.h"

namespace spinas {

  eeWW::eeWW(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ), propne(0,0), propA(0,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    propZ = propagator(MZ,WZ);
    proph = propagator(mh,wh);  
    p1=particle(me);
    p2=particle(me);
    p3=particle(MW);
    p4=particle(MW);
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
    //[413>,[231>,[132>
    s413a = sproduct(SQUARE,&p4,&p1,&p3);
    s231a = sproduct(SQUARE,&p2,&p3,&p1);
    s132a = sproduct(SQUARE,&p1,&p3,&p2);
    //Couplings
    prene = e*e/(MW*MW*SW*SW);
    preh = e*e*me/(2.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
  }
  void eeWW::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    mh=massh;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(me);
    p2.set_mass(me);
    p3.set_mass(MW);
    p4.set_mass(MW);
    proph.set_mass(mh);
    propZ.set_mass(MZ);
    //Couplings
    prene = e*e/(MW*MW*SW*SW);
    preh = e*e*me/(2.0*MW*MW*SW*SW);
    preZ = e*e/(2.0*MW*MW*SW*SW);
  }
  void eeWW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    //[413>,[231>,[132>
    s413a.update();
    s231a.update();
    s132a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenhS=proph.denominator(propSP);//std::cout<<"pDenS="<<pDenS<<"\n";
    pDenneT=propne.denominator(propTP);
    pDenneU=propne.denominator(propUP);
    pDenZS=propZ.denominator(propSP);
    pDenAS=propA.denominator(propSP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eeWW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds4a, ds4b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);
      
      //S-Channel h
      //preh = e*e*me/(2.0*MW*MW*SW*SW);
      //all ingoing: 
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      //34 outgoing:
      //preh [34] <34> ([12]+<12>)/(s-Mh^2)
      amplitude += normFactor*preh*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))/pDenhS;
      
      //U-Channel ne
      //prene = e*e/(MW*MW*SW*SW);
      //all ingoing:
      // - prene [23] <14> ([413>-MW <34>))/u
      //34 outgoing:
      // - prene [23] <14> (MW <34> + [413>)/u
      amplitude += - normFactor*prene*s23s.v(ds2,ds3a)*a14a.v(ds1,ds4a)*(MW*a34a.v(ds3b,ds4b)+s413a.v(ds4b,ds3b))/pDenneU;

      //S-Channel A
      //all in:
      //+ 2e^2/MW*( <13>[24] + [13]<24> + <14>[23] + [14]<23> )( <34> + [34] )/s
      //+ 2e^2/MW/MW*[34]<34>([231>+[132>)/s
      amplitude +=
	- normFactor*e*e*two/MW*(
				a13a.v(ds1,ds3a)*s24s.v(ds2,ds4a) + s13s.v(ds1,ds3a)*a24a.v(ds2,ds4a)
				+ a14a.v(ds1,ds4a)*s23s.v(ds2,ds3a) + s14s.v(ds1,ds4a)*a23a.v(ds2,ds3a)
				)*(s34s.v(ds3b,ds4b)+a34a.v(ds3b,ds4b))/pDenAS
	- normFactor*e*e*two/MW/MW*(
				s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s231a.v(ds2,ds1)+s132a.v(ds1,ds2))
				)/pDenAS;

      //S-Channel Z
      //preZ = e*e/(2.0*MW*MW*SW*SW); //=prene
      //all ingoing
      //- preZ gLe ( Me[34]<34>( [12]-<12>) + 2[34]<34>[231> + 2MW([23]<14>+[24]<13>)([34]+<34>) )/(s-MZ^2)
      //- preZ gRe ( Me[34]<34>(-[12]+<12>) + 2[34]<34>[132> + 2MW([13]<24>+[14]<23>)([34]+<34>) )/(s-MZ^2)
      // = - preZ ( (gL-gR)Me[34]<34>([12]-<12>)
      //          + 2[34]<34>(gL[231>+gR[132>)
      //          + 2MW([34]+<34>)( gL([23]<14>+[24]<13>) + gR([13]<24>+[14]<23>) )
      //        )/(s-MZ^2)
      //34 outgoing
      //- preZ ( (gL-gR)Me[34]<34>([12]-<12>)
      //        - 2[34]<34>(gL[231>+gR[132>)
      //        - 2MW([34]+<34>)( gL([23]<14>+[24]<13>) + gR([13]<24>+[14]<23>) )
      //      )/(s-MZ^2)
      amplitude += - normFactor*preZ*( (gL-gR)*me*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(s12s.v(ds1,ds2)-a12a.v(ds1,ds2))
				       -two*s34s.v(ds3a,ds4a)*a34a.v(ds3b,ds4b)*(gL*s231a.v(ds2,ds1)+gR*s132a.v(ds1,ds2))
				       -two*MW*(s34s.v(ds3a,ds4a)+a34a.v(ds3a,ds4a))*( gL*(s23s.v(ds2,ds3b)*a14a.v(ds1,ds4b)+s24s.v(ds2,ds4b)*a13a.v(ds1,ds3b))
										       + gR*(s13s.v(ds1,ds3b)*a24a.v(ds2,ds4b)+s14s.v(ds1,ds4b)*a23a.v(ds2,ds3b)) )
				       )/pDenZS;
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eeWW::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2^2=1/4
    return amp2/4.0;
  }

  
  

  //  Tests
  int test_eeWW(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> W+, W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=250\n";
      ldouble me=0.0005, mh=125;
      ldouble EE=0.31333,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      eeWW eeWWAmp = eeWW(EE,me,mh,0,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {4.427057986735773E-03,7.735122730595210E-03,1.082595495325175E-02,1.376043913921246E-02,1.661514324546994E-02,1.948771624985716E-02,2.250467958382392E-02,2.583294999057593E-02,2.969735952269897E-02,3.440815999832256E-02,4.040584697251608E-02,4.833751887277057E-02,5.919409307422187E-02,7.457363186517232E-02,9.723045616930849E-02,1.323518114173552E-01,1.910054365163739E-01,3.017679146640585E-01,5.679949904526829E-01,1.874672973859846E+00};
      i += eeWWAmp.test_2to2_amp2([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH);
      i += eeWWAmp.test_2to2_amp2_rotations([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH);
      i += eeWWAmp.test_2to2_amp2_boosts([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH);
      i += eeWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=2500\n";
      pspatial = 2500;
      ldouble dataCH1[20] = {1.770562268328381E-03,5.013830738654090E-03,7.963637936796804E-03,1.067503121327328E-02,1.321726367870653E-02,1.567869276146775E-02,1.817385584750952E-02,2.085394228715550E-02,2.392272900634522E-02,2.766162223609481E-02,3.247051500073769E-02,3.893748521636568E-02,4.796425222114187E-02,6.100743937799832E-02,8.058303195035148E-02,1.114434452107251E-01,1.637724220467787E-01,2.640476765361803E-01,5.093540344037539E-01,1.773450774162739E+00};
      i += eeWWAmp.test_2to2_amp2([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH1);
      i += eeWWAmp.test_2to2_amp2_rotations([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH1);
      //std::cout<<"\n# me=0.0005, mh=125, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2[20] = {1.639183344855925E-02,2.017916398967428E-02,2.401656342566092E-02,2.798387132183932E-02,3.218077415649641E-02,3.673337939295766E-02,4.180358722607697E-02,4.760275452548505E-02,5.441213288371305E-02,6.261434957100698E-02,7.274357648115047E-02,8.556873936123811E-02,1.022382538918203E-01,1.245467402516593E-01,1.554629661960497E-01,2.002746508857516E-01,2.693889134018722E-01,3.864305683179654E-01,6.182238534697527E-01,1.229539589641644E+00};
      i += eeWWAmp.test_2to2_amp2([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH2);
      i += eeWWAmp.test_2to2_amp2_rotations([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH2);
      i += eeWWAmp.test_2to2_amp2_boosts([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH2);
      i += eeWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, mh=125, MW=80.385, pspatial=95\n";
      me = 125.1;
      mh = 125;
      pspatial = 95;
      eeWWAmp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {1.058333088017236E-01,1.179546365511695E-01,1.316057537892103E-01,1.470766821403135E-01,1.647270803637689E-01,1.850084518480779E-01,2.084954154374346E-01,2.359306609567241E-01,2.682911081850997E-01,3.068878954929486E-01,3.535221793552363E-01,4.107366275340256E-01,4.822385149070731E-01,5.736471893472228E-01,6.938941458432443E-01,8.580370143322010E-01,1.093411485907179E+00,1.454361977227176E+00,2.058113932129806E+00,3.031957493793098E+00};
      i += eeWWAmp.test_2to2_amp2([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH4);
      i += eeWWAmp.test_2to2_amp2_rotations([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH4);
      i += eeWWAmp.test_2to2_amp2_boosts([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH4);
      i += eeWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH4);
      //std::cout<<"\n# me=125, mh=0.0005, MW=80.385, pspatial=125.1\n";
      me = 125;
      mh = 0.0005;
      pspatial = 125.1;
      eeWWAmp.set_masses(me,mh,MW);
      ldouble dataCH3[20] = {8.570021274904190E-02,9.772641133345955E-02,1.111654897160825E-01,1.263115059933083E-01,1.435320203821205E-01,1.632925852618677E-01,1.861917228698985E-01,2.130120311065706E-01,2.447968116660867E-01,2.829684321922070E-01,3.295175686439002E-01,3.873182930320367E-01,4.606786916431793E-01,5.563613038607091E-01,6.856174585757228E-01,8.686394739689423E-01,1.145588175619176E+00,1.609052378672457E+00,2.527452392828179E+00,4.983827317946191E+00};
      i += eeWWAmp.test_2to2_amp2([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH3);
      i += eeWWAmp.test_2to2_amp2_rotations([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH3);
      i += eeWWAmp.test_2to2_amp2_boosts([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH3);
      i += eeWWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eeWWAmp.amp2(); }, me,me,MW,MW,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }


  
    
  

}
