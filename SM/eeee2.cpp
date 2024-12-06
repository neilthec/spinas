
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

//File:  SPINAS/SM/eeee2.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eeee2.h"

namespace spinas {
  //Constructors
  eeee2::eeee2(const ldouble& echarge, const ldouble& masse, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), WZ(widthZ),
    propA(0,0), proph(mh,wh),
    p1(particle(me)), p2(particle(me)),
    p3(particle(me)), p4(particle(me)),
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
    preh = e*e*me*me/(4.0*MW*MW*SW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*me*me/MZ/MZ;//=preh!
  }
  void eeee2::set_masses(const ldouble& masse, const ldouble& massh, const ldouble& massW){
    me=masse;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(me);
    p2.set_mass(me);
    p3.set_mass(me);
    p4.set_mass(me);
    preh = e*e*me*me/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*me*me/MZ/MZ;
  }
  void eeee2::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenUA = propA.denominator(propPU);
    pDenUh = proph.denominator(propPU);
    pDenUZ = propZ.denominator(propPU);
    pDenTA = propA.denominator(propPT);
    pDenTh = proph.denominator(propPT);
    pDenTZ = propZ.denominator(propPT);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eeee2::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble one=1, two = 2, four=4, eight=8;
    cdouble amplitude(0,0);
    
    //Photon
    //U Channel
    //eEEe all ingoing:
    //- (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    //eeEE: 2<->4
    //- (- <13>[24] - <12>[34] - [13]<24> - [12]<34>)
    //34 out:
    //- (<13>[24] - <12>[34] + [13]<24> - [12]<34>)
    amplitude += - two*e*e*(
			    + a13a.v(ds1,ds3)*s24s.v(ds2,ds4) - a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
			    + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) - s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
			    )/pDenUA;

    //T Channel
    // 2 <-> 3 with a minus sign for fermions
    //eEEe all ingoing:
    //  (<12>[34] - <14>[23] + [12]<34> - [14]<23>)
    //eeEE: 2<->4
    //- (<14>[23] - <12>[34] + [14]<23> - [12]<34>)
    //34 out:
    //+ (<14>[23] + <12>[34] + [14]<23> + [12]<34>)
    amplitude += two*e*e*(
			  + a14a.v(ds1,ds4)*s23s.v(ds2,ds3) + a12a.v(ds1,ds2)*s34s.v(ds3,ds4)
			  + s14s.v(ds1,ds4)*a23a.v(ds2,ds3) + s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
			  )/pDenTA;
    
    //Higgs
    //U Channel
    //preh = e*e*me*me/(4*MW*MW*SW*SW);
    //eEEe all in:
    //preh ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //eeEE: 2<->4
    //- preh ([14]+<14>) ([23]+<23>)/(u-Mh^2)
    //34 out:
    //- preh ([14]-<14>) ([23]-<23>)/(u-Mh^2)
    amplitude += - preh*(s14s.v(ds1,ds4)-a14a.v(ds1,ds4))*(s23s.v(ds2,ds3)-a23a.v(ds2,ds3))/pDenUh;

    //T Channel
    //eEEe all ingoing:
    //- preh * ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //eeEE: 2<->4
    //+ preh * ([13]+<13>) ([24]+<24>)/(t-Mh^2)
    //34 out:
    //+ preh * ([13]-<13>) ([24]-<24>)/(t-Mh^2)
    amplitude += + preh*(s13s.v(ds1,ds3)-a13a.v(ds1,ds3))*(s24s.v(ds2,ds4)-a24a.v(ds2,ds4))/pDenTh;
    
    
    //Z Boson
    //Defined above:
    //gL=2.0*SW*SW-1.0;
    //gR=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gL-gR)*(gL-gR)*me*me/MZ/MZ; // = preh
    //U Channel
    //eEEe all in:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    //eeEE: 2<->4
    //+ preZ0 (<14>-[14]) (<23>-[23]) / (u-MZ^2)
    //+ preZ (gL^2 [34] <12> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [12] <34>))/(u-MZ^2)
    //34 out:
    //+ preZ0 (<14>+[14]) (<23>+[23]) / (u-MZ^2)
    //+ preZ (gL^2 [34] <12> - gLgR( [13] <24>+ [24] <13> ) + gR^2 [12] <34>))/(u-MZ^2)
    amplitude += 
      + preZ0*(a14a.v(ds1,ds4)+s14s.v(ds1,ds4))*(a23a.v(ds2,ds3)+s23s.v(ds2,ds3))/pDenUZ
      + two*preZ*(
	      gL*gL*s34s.v(ds3,ds4)*a12a.v(ds1,ds2)
	      - gL*gR*(s13s.v(ds1,ds3)*a24a.v(ds2,ds4)+s24s.v(ds2,ds4)*a13a.v(ds1,ds3))
	      + gR*gR*s12s.v(ds1,ds2)*a34a.v(ds3,ds4)
	      )/pDenUZ;

    //T Channel
    //eEEe all in:
    //+ preZ0 (<13>-[13]) (<24>-[24]) / (t-MZ^2)
    //+ preZ (- gL^2 [23] <14> + gLgR( [12] <34>+ [34] <12> ) - gR^2 [14] <23>))/(t-MZ^2)
    //eeEE: 2<->4
    //- preZ0 (<13>-[13]) (<24>-[24]) / (t-MZ^2)
    //+ preZ (+ gL^2 [34] <12> - gLgR( [14] <23>+ [23] <14> ) + gR^2 [12] <34>))/(t-MZ^2)
    //34 out:
    //- preZ0 (<13>+[13]) (<24>+[24]) / (t-MZ^2)
    //+ preZ (+ gL^2 [34] <12> + gLgR( [14] <23>+ [23] <14> ) + gR^2 [12] <34>))/(t-MZ^2)
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
  ldouble eeee2::amp2(){
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
    //Symmetry factor 1/2
    return amp2/8.0;
  }

  



  //  Tests
  int test_eeee2(){
    int n=0;//Number of fails
    std::cout<<"\t* e , e  -> e , e       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### me=0.0005, pspatial=250\n";
      ldouble me=0.0005, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      eeee2 eeee2Amp = eeee2(0.31333,me,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.239151231292354E+01,4.010353993958522E+00,1.596640108811971E+00,9.036288843569612E-01,6.109477957934305E-01,4.616466545511128E-01,3.774798480779697E-01,3.281804602926848E-01,3.003213277005353E-01,2.876559449425260E-01,2.876559449425260E-01,3.003213277005352E-01,3.281804602926847E-01,3.774798480779697E-01,4.616466545511127E-01,6.109477957934304E-01,9.036288843569608E-01,1.596640108811970E+00,4.010353993958518E+00,3.239151231292340E+01};
      i += eeee2Amp.test_2to2_amp2([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH);
      i += eeee2Amp.test_2to2_amp2_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH);
      i += eeee2Amp.test_2to2_amp2_boosts([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH);
      i += eeee2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH);
      //std::cout<<"########### me=0.0005, pspatial=0.11\n";
      pspatial = 0.005;
      ldouble dataCH2[20] = {3.118412297136729E+01,3.500612710324413E+00,1.289032220945836E+00,6.823058037131237E-01,4.351814273586435E-01,3.128065781619428E-01,2.454620207821479E-01,2.067506971417634E-01,1.851685576636600E-01,1.754365498610951E-01,1.754365498610951E-01,1.851685576636599E-01,2.067506971417634E-01,2.454620207821479E-01,3.128065781619428E-01,4.351814273586434E-01,6.823058037131232E-01,1.289032220945835E+00,3.500612710324409E+00,3.118412297136714E+01};
      i += eeee2Amp.test_2to2_amp2([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH2);
      i += eeee2Amp.test_2to2_amp2_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH2);
      i += eeee2Amp.test_2to2_amp2_boosts([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH2);
      i += eeee2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH2);
      //std::cout<<"########### me=0.10, pspatial=0.05\n";
      me=0.10;
      pspatial = 0.05;
      eeee2Amp.set_masses(me,mh,MW);
      ldouble dataCH4[20] = {2.714628174468743E+02,2.884221171993742E+01,9.939423699160043E+00,4.871078261033936E+00,2.850751442644789E+00,1.869966338400375E+00,1.339611690477008E+00,1.039155908219962E+00,8.734691469486325E-01,7.992619792629213E-01,7.992619792629212E-01,8.734691469486321E-01,1.039155908219962E+00,1.339611690477008E+00,1.869966338400374E+00,2.850751442644788E+00,4.871078261033932E+00,9.939423699160031E+00,2.884221171993739E+01,2.714628174468731E+02};
      i += eeee2Amp.test_2to2_amp2([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH4);
      i += eeee2Amp.test_2to2_amp2_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH4);
      i += eeee2Amp.test_2to2_amp2_boosts([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH4);
      i += eeee2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH4);
      //std::cout<<"########### me=0.10, MW=0.11, pspatial=0.05\n";
      me=0.10;
      MW=0.11;
      pspatial = 0.05;
      eeee2Amp.set_masses(me,mh,MW);
      ldouble dataCH5[20] = {2.743314972284846E+02,2.990404285701647E+01,1.064231860221773E+01,5.422489272565070E+00,3.320381827157619E+00,2.290043884415155E+00,1.728061052870422E+00,1.407363428106136E+00,1.229534219030809E+00,1.149612740084618E+00,1.149612740084618E+00,1.229534219030809E+00,1.407363428106136E+00,1.728061052870422E+00,2.290043884415155E+00,3.320381827157618E+00,5.422489272565064E+00,1.064231860221772E+01,2.990404285701645E+01,2.743314972284833E+02};
      i += eeee2Amp.test_2to2_amp2([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH5);
      i += eeee2Amp.test_2to2_amp2_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH5);
      i += eeee2Amp.test_2to2_amp2_boosts([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH5);
      i += eeee2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH5);
      //std::cout<<"########### me=0.10, MW=0.006, pspatial=0.05\n";
      me=0.10;
      MW=0.006;
      pspatial = 0.05;
      eeee2Amp.set_masses(me,mh,MW);
      ldouble dataCH6[20] = {2.028069761124929E+03,1.529963393073678E+03,1.464264596649519E+03,1.440017290079247E+03,1.427816510757349E+03,1.420738197117430E+03,1.416348757213876E+03,1.413593740433501E+03,1.411961913761497E+03,1.411199466445355E+03,1.411199466445355E+03,1.411961913761497E+03,1.413593740433501E+03,1.416348757213876E+03,1.420738197117431E+03,1.427816510757349E+03,1.440017290079247E+03,1.464264596649519E+03,1.529963393073678E+03,2.028069761124927E+03};
      i += eeee2Amp.test_2to2_amp2([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH6);
      i += eeee2Amp.test_2to2_amp2_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH6);
      i += eeee2Amp.test_2to2_amp2_boosts([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH6);
      i += eeee2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH6);
      //std::cout<<"########### me=0.10, MW=0.11, Mh=0.125, pspatial=0.05\n";
      me=0.10;
      MW=0.11;
      mh=0.125;
      pspatial = 0.05;
      eeee2Amp.set_masses(me,mh,MW);
      ldouble dataCH7[20] = {2.709677106717437E+02,2.878113133527047E+01,9.962394448288221E+00,4.928301995174162E+00,2.925743316598655E+00,1.955291736034680E+00,1.431282082155070E+00,1.134753151970281E+00,9.713618302713348E-01,8.982174228927258E-01,8.982174228927257E-01,9.713618302713344E-01,1.134753151970281E+00,1.431282082155070E+00,1.955291736034680E+00,2.925743316598654E+00,4.928301995174158E+00,9.962394448288210E+00,2.878113133527044E+01,2.709677106717425E+02};
      i += eeee2Amp.test_2to2_amp2([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH7);
      i += eeee2Amp.test_2to2_amp2_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH7);
      i += eeee2Amp.test_2to2_amp2_boosts([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH7);
      i += eeee2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH7);
      //std::cout<<"########### me=0.10, MW=0.006, Mh=0.125, pspatial=0.05\n";
      me=0.10;
      MW=0.006;
      mh=0.125;
      pspatial = 0.05;
      eeee2Amp.set_masses(me,mh,MW);
      ldouble dataCH8[20] = {9.143086050966065E+02,1.361318339859029E+03,1.468939003124573E+03,1.507165554463950E+03,1.522250748919556E+03,1.527879441236381E+03,1.529406427266916E+03,1.529276899582259E+03,1.528702585560780E+03,1.528288217896787E+03,1.528288217896787E+03,1.528702585560780E+03,1.529276899582259E+03,1.529406427266916E+03,1.527879441236382E+03,1.522250748919557E+03,1.507165554463950E+03,1.468939003124574E+03,1.361318339859030E+03,9.143086050966081E+02};
      i += eeee2Amp.test_2to2_amp2([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH8);
      i += eeee2Amp.test_2to2_amp2_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH8);
      i += eeee2Amp.test_2to2_amp2_boosts([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH8);
      i += eeee2Amp.test_2to2_amp2_boosts_and_rotations([&]() { return eeee2Amp.amp2(); }, me,me,me,me,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
