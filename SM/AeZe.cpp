
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

//File:  SPINAS/SM/AeZe.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AeZe.h"

namespace spinas {

  AeZe::AeZe(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW):
    e(echarge), me(masse), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), prope(masse,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    p1=particle(0);
    p2=particle(me);
    p3=particle(MZ);
    p4=particle(me);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    //[1421], <1421>
    s1421s = sproduct(SQUARE,&p1,&p4,&p2,&p1);
    a1421a = sproduct(ANGLE,&p1,&p4,&p2,&p1);
    //Couplings
    preTU = 2.0*e*e/(2.0*MW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
  }
  void AeZe::set_masses(const ldouble& masse, const ldouble& massW){
    me=masse;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p2.set_mass(me);
    p3.set_mass(MZ);
    p4.set_mass(me);
    prope.set_mass(me);
    //Couplings
    preTU = 2.0*e*e/(2.0*MW*SW);
  }
  void AeZe::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13],<34>,[34],<14>,[14]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s34s.update();
    a34a.update();
    s14s.update();
    a14a.update();
    //[1421], <1421>
    s1421s.update();
    a1421a.update();
    //Propagator Momentum
    ldouble propUP[4], propSP[4];
    for(int j=0;j<4;j++){
      propUP[j] = mom1[j]-mom4[j];
      propSP[j] = mom1[j]+mom2[j];
    }
    pDenU=prope.den(propUP);
    pDenS=prope.den(propSP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AeZe::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      //preTU = 2.0*e*e/(2.0*MW*SW);
      
      if(ds1>0){
	//AZEe all ingoing:
	//- preTU [1341] (gRe [24] <23> + gLe <24> [23])/((t-me^2) (u-me^2))
	//- preTU [12] ( gRe [14]<23>/(u-me^2) - gLe [13]<24>/(t-me^2) )
	//AeZE: 4->2->3->4
	//+ preTU [1421] (gRe [23] <34> + gLe <23> [34])/((u-me^2) (s-me^2))
	//- preTU [13] ( gRe [12]<34>/(s-me^2) + gLe [14]<23>/(u-me^2) )
	//AeZE: HEPCAT: Same up to overall sign
	//+ preTU [1241] (gRe <34> [23] + gLe <23> [34])/((s-me^2)(u-me^2))
	//+ preTU [13] ( gRe [12]<34>/(s-me^2) + gLe [14]<23>/(u-me^2))
	//34 out:
	//- preTU [1421] (gRe [23] <34> - gLe <23> [34])/((u-me^2) (s-me^2))
	//- preTU [13] ( gRe [12]<34>/(s-me^2) - gLe [14]<23>/(u-me^2) )
	amplitude += - normFactor*preTU*s1421s.v()*(
						    + gR*s23s.v(ds2,ds3a)*a34a.v(ds3b,ds4)
						    - gL*a23a.v(ds2,ds3a)*s34s.v(ds3b,ds4)
						    )/pDenU/pDenS;
	amplitude += - normFactor*preTU*s13s.v(ds3a)*(
						      + gR*s12s.v(ds2)*a34a.v(ds3b,ds4)/pDenS
						      - gL*s14s.v(ds4)*a23a.v(ds2,ds3b)/pDenU
						      );

	
      }
      else if(ds1<0){
	//AZEe all ingoing:
	//- preTU <1341> (gRe [24] <23> + gLe <24> [23])/((t-me^2) (u-me^2))
	//- preTU <12> ( gLe <14>[23]/(u-me^2) - gRe <13>[24]/(t-me^2) )
	//AeZE: 4->2->3->4	
	//+ preTU <1421> (gRe [23] <34> + gLe <23> [34])/((u-me^2) (s-me^2))
	//- preTU <13> ( gLe <12>[34]/(s-me^2) + gRe <14>[23]/(u-me^2) )
	//34 out:
	//- preTU <1421> (gRe [23] <34> - gLe <23> [34])/((u-me^2) (s-me^2))
	//+ preTU <13> ( gLe <12>[34]/(s-me^2) - gRe <14>[23]/(u-me^2) )
	amplitude += - normFactor*preTU*a1421a.v()*(
						    + gR*s23s.v(ds2,ds3a)*a34a.v(ds3b,ds4)
						    - gL*a23a.v(ds2,ds3a)*s34s.v(ds3b,ds4)
						    )/pDenU/pDenS;
	amplitude += + normFactor*preTU*a13a.v(ds3a)*(
						      + gL*a12a.v(ds2)*s34s.v(ds3b,ds4)/pDenS
						      - gR*a14a.v(ds4)*s23s.v(ds2,ds3b)/pDenU
						      );

      }


      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AeZe::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2=1/4
    return amp2/4.0;
  }
  //set_momenta(...) must be called before amp2_Aplus().
  ldouble AeZe::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2
    return amp2/2.0;
  }  
  //set_momenta(...) must be called before amp2_Aminus().
  ldouble AeZe::amp2_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(-2,j2,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/2
    return amp2/2.0;
  }  



  //  Tests
  int test_AeZe(){
    int n=0;//Number of fails
    std::cout<<"\t* A , e  -> Z , e       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=250\n";
      ldouble me=0.0005;
      ldouble EE=0.31333, MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      AeZe AeZeAmp = AeZe(EE,me,MW,SW);
      ldouble pspatial=250;
      ldouble dataCHp[20] = {1.401372179722165E-02,1.390818340270411E-02,1.383691997220342E-02,1.380616331735902E-02,1.382375345926568E-02,1.389969318345807E-02,1.404694906377748E-02,1.428263697769437E-02,1.462982615136459E-02,1.512037408757479E-02,1.579955197449542E-02,1.673393508884779E-02,1.802560552160443E-02,1.983946511202866E-02,2.246034430975693E-02,2.642625618632420E-02,3.289008589501827E-02,4.484933483438120E-02,7.329484146568851E-02,2.171441519389301E-01};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCHp);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCHp);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCHp);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCHp);
      ldouble dataCHm[20] = {1.396075528994397E-02,1.414067300732706E-02,1.437200937321001E-02,1.466411323277459E-02,1.502874603641079E-02,1.548091377254718E-02,1.604006864840596E-02,1.673188757752050E-02,1.759097857378348E-02,1.866513365638007E-02,2.002226779994171E-02,2.176225596716165E-02,2.403823976750264E-02,2.709760180086565E-02,3.136763928526057E-02,3.765546918803702E-02,4.769062786779412E-02,6.596990985685916E-02,1.089815617735448E-01,3.251198031753933E-01};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCHm);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCHm);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCHm);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCHm);
      ldouble dataCH[20] = {1.398723854358281E-02,1.402442820501558E-02,1.410446467270672E-02,1.423513827506681E-02,1.442624974783823E-02,1.469030347800263E-02,1.504350885609172E-02,1.550726227760743E-02,1.611040236257403E-02,1.689275387197743E-02,1.791090988721857E-02,1.924809552800472E-02,2.103192264455353E-02,2.346853345644715E-02,2.691399179750876E-02,3.204086268718061E-02,4.029035688140620E-02,5.540962234562018E-02,9.113820161961668E-02,2.711319775571617E-01};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=92\n";
      pspatial = 92;
      ldouble dataCH2p[20] = {1.530603876632118E-02,1.526070486688893E-02,1.524632451316218E-02,1.526852562254193E-02,1.533438847821013E-02,1.545294654491616E-02,1.563590986955738E-02,1.589873571403832E-02,1.626225777918686E-02,1.675524641376537E-02,1.741858579759416E-02,1.831239972452677E-02,1.952887801180564E-02,2.121694267189988E-02,2.363382263543949E-02,2.726539477139581E-02,3.315282379900462E-02,4.400313798832712E-02,6.974210766326845E-02,1.997016657504402E-01};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH2p);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH2p);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH2p);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH2p);
      ldouble dataCH2m[20] = {1.380224449842647E-02,1.396716757293270E-02,1.417500190234123E-02,1.443354953299782E-02,1.475262594256323E-02,1.514475432667357E-02,1.562616845744964E-02,1.621829691409083E-02,1.695002169464831E-02,1.786122746309727E-02,1.900859242613888E-02,2.047546688774469E-02,2.238965464733640E-02,2.495760799628461E-02,2.853592634695036E-02,3.379818636711009E-02,4.218776663051866E-02,5.745745142872134E-02,9.336744982914678E-02,2.737596433868220E-01};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH2m);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH2m);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH2m);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH2m);
      ldouble dataCH2[20] = {1.455414163237382E-02,1.461393621991082E-02,1.471066320775171E-02,1.485103757776988E-02,1.504350721038668E-02,1.529885043579486E-02,1.563103916350351E-02,1.605851631406457E-02,1.660613973691758E-02,1.730823693843132E-02,1.821358911186652E-02,1.939393330613573E-02,2.095926632957102E-02,2.308727533409224E-02,2.608487449119493E-02,3.053179056925294E-02,3.767029521476164E-02,5.073029470852423E-02,8.155477874620762E-02,2.367306545686311E-01};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH2);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH2);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH2);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, MW=80.385, pspatial=250\n";
      me = 125;
      pspatial = 250;
      AeZeAmp.set_masses(me,MW);
      ldouble dataCH4p[20] = {1.417061666608737E-02,1.444541937738647E-02,1.491693888537775E-02,1.561839755456245E-02,1.659094498924402E-02,1.788617243875278E-02,1.956964906605391E-02,2.172599978524518E-02,2.446637011449536E-02,2.793969844632593E-02,3.235027184010787E-02,3.798607027746888E-02,4.526651670980986E-02,5.482712989721544E-02,6.767927237769514E-02,8.553619724623671E-02,1.115497635306112E-01,1.522216216939633E-01,2.234636658170986E-01,3.772763813551800E-01};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH4p);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH4p);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH4p);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH4p);
      ldouble dataCH4m[20] = {1.411864508129129E-02,1.465597431140994E-02,1.540058491818810E-02,1.638708968992425E-02,1.765826086082050E-02,1.926758458360794E-02,2.128283897504735E-02,2.379120835826977E-02,2.690676225424271E-02,3.078168084856972E-02,3.562361496304542E-02,4.172348147152304E-02,4.950182243972506E-02,5.958999112301089E-02,7.298101945891747E-02,9.133146305800563E-02,1.176252743725313E-01,1.578302902664077E-01,2.258347407451065E-01,3.627146533938100E-01};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH4m);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH4m);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH4m);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH4m);
      ldouble dataCH4[20] = {1.414463087368933E-02,1.455069684439821E-02,1.515876190178293E-02,1.600274362224335E-02,1.712460292503226E-02,1.857687851118036E-02,2.042624402055063E-02,2.275860407175748E-02,2.568656618436903E-02,2.936068964744783E-02,3.398694340157664E-02,3.985477587449596E-02,4.738416957476746E-02,5.720856051011317E-02,7.033014591830630E-02,8.843383015212117E-02,1.145875189515713E-01,1.550259559801855E-01,2.246492032811025E-01,3.699955173744950E-01};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH4);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH4);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH4);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH4);
      //std::cout<<"\n# me=125, MW=80.385, pspatial=92\n";
      pspatial = 92;
      ldouble dataCH3p[20] = {1.667583780235594E-02,1.757081037940531E-02,1.856244905049405E-02,1.965878233847982E-02,2.086870967035667E-02,2.220211262901593E-02,2.366998059826334E-02,2.528455175075007E-02,2.705946964566245E-02,2.900995437332637E-02,3.115298474617320E-02,3.350748370261304E-02,3.609449150470716E-02,3.893729810071552E-02,4.206148298329345E-02,4.549477040770977E-02,4.926653594073952E-02,5.340667056625567E-02,5.794326975514632E-02,6.289816463950243E-02};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH3p);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH3p);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH3p);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2_Aplus(); }, 0,me,MZ,me,pspatial,dataCH3p);
      ldouble dataCH3m[20] = {1.544592128016146E-02,1.640380022507265E-02,1.744833428886892E-02,1.858584523281800E-02,1.982318845923605E-02,2.116778497718239E-02,2.262764569573469E-02,2.421138179472063E-02,2.592819151150301E-02,2.778780843297382E-02,2.980038823785982E-02,3.197629807460216E-02,3.432575253570469E-02,3.685820770296550E-02,3.958137175374496E-02,4.249960269722415E-02,4.561131503803231E-02,4.890475981590993E-02,5.235108606717637E-02,5.589275914121893E-02};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH3m);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH3m);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH3m);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2_Aminus(); }, 0,me,MZ,me,pspatial,dataCH3m);
      ldouble dataCH3[20] = {1.606087954125870E-02,1.698730530223898E-02,1.800539166968149E-02,1.912231378564891E-02,2.034594906479636E-02,2.168494880309916E-02,2.314881314699901E-02,2.474796677273535E-02,2.649383057858273E-02,2.839888140315010E-02,3.047668649201651E-02,3.274189088860760E-02,3.521012202020592E-02,3.789775290184051E-02,4.082142736851920E-02,4.399718655246696E-02,4.743892548938592E-02,5.115571519108280E-02,5.514717791116135E-02,5.939546189036068E-02};
      i += AeZeAmp.test_2to2_amp2([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH3);
      i += AeZeAmp.test_2to2_amp2_rotations([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH3);
      i += AeZeAmp.test_2to2_amp2_boosts([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH3);
      i += AeZeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AeZeAmp.amp2(); }, 0,me,MZ,me,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
