
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

//File:  SPINAS/SM/uuss.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uuss.h"

namespace spinas {
  //Constructors
  uuss::uuss(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& masss, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), gs(gscharge), mu(massu), ms(masss), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propAG(0,0), proph(mh,wh), propZ(MZ,WZ),
    p1(particle(mu)), p2(particle(mu)),
    p3(particle(ms)), p4(particle(ms)),
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
    preh = e*e*mu*ms/(4.0*MW*MW*SW*SW);
    gLu=1.0-4.0/3.0*SW*SW;
    gRu=-4.0/3.0*SW*SW;
    gLd=-1.0+2.0/3.0*SW*SW;
    gRd=2.0/3.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*ms/MZ/MZ;//=-preh!
  }
  void uuss::set_masses(const ldouble& massu, const ldouble& masss, const ldouble& massh, const ldouble& massW){
    mu=massu;
    ms=masss;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(mu);
    p2.set_mass(mu);
    p3.set_mass(ms);
    p4.set_mass(ms);
    preh = e*e*mu*ms/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*ms/MZ/MZ;//=-preh
  }
  void uuss::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
      propP[j] = mom1[j]+mom2[j];
    pDenSAG = propAG.denominator(propP);
    pDenSh = proph.denominator(propP);
    pDenSZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  //We have to separate the gluon so we can separate the color factor between the gluon^2 diagram
  // And the rest^2.
  // The cross term vanishes due to the trace of the adjoint rep.
  cdouble uuss::amp_gluon(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
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
  cdouble uuss::amp_rest(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble two = 2;
    cdouble amplitude(0,0);
    
    //Photon
    //Sign changes due to p3 and p4 being outgoing.
    // (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    amplitude = -two*e*e*2.0/9.0*(
			  a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenSAG;
    
    //Higgs
    //EE^2 Mu Ms / (4 MW^2 SW^2) * ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //preh = e*e*mu*ms/(4*MW*MW*SW*SW);
    amplitude += preh*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))*(s34s.v(ds3,ds4)+a34a.v(ds3,ds4))/pDenSh;
    
    //Z Boson
    //Defined above:
    //gLu=1.0-4.0/3.0*SW*SW;
    //gRu=-4.0/3.0*SW*SW;
    //gLd=-1.0+2.0/3.0*SW*SW;
    //gRd=2.0/3.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*ms/MZ/MZ; // = -preh
    //all in:
    //
    //+(EE^2 Mu Ms (gLu-gRu)(gLd-gRd) (<12>-[12]) (<34>-[34]))/(8 CW^2 MZ^2 SW^2 (s-MZ^2))
    //+(EE^2 ( gLu gLd [23] <14> + gLd gRu [13] <24> + gLu gRd [24] <13> + gRu gRd [14] <23>)/(4 CW^2 SW^2 (s-MZ^2))
    //= + preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //  + preZ ( gLu gLd [23] <14> + gLd gRu [13] <24> + gLu gRd [24] <13> + gRu gRd [14] <23> )/(s-MZ^2)
    //34 out:
    //+ preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //- preZ ( gLu gLd [23] <14> + gLd gRu [13] <24> + gLu gRd [24] <13> + gRu gRd [14] <23> )/(s-MZ^2)
    amplitude += 
      - preZ0*(a12a.v(ds1,ds2)-s12s.v(ds1,ds2))*(a34a.v(ds3,ds4)-s34s.v(ds3,ds4))/pDenSZ
      + two*preZ*(
	        gLu*gLd*s23s.v(ds2,ds3)*a14a.v(ds1,ds4)
	      + gLd*gRu*s13s.v(ds1,ds3)*a24a.v(ds2,ds4)
	      + gLu*gRd*s24s.v(ds2,ds4)*a13a.v(ds1,ds3)
	      + gRu*gRd*s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
	      )/pDenSZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uuss::amp2(){
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
  int test_uuss(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> s , S       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mu=0.0042, ms=1.23, pspatial=250\n";
      ldouble mu=0.0042, ms=1.23, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      uuss uussAmp = uuss(0.31333,1.238,mu,ms,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.005084973580743E+00,9.099315655191025E-01,8.252822867023432E-01,7.511371371304647E-01,6.874961168034672E-01,6.343592257213504E-01,5.917264638841144E-01,5.595978312917594E-01,5.379733279442851E-01,5.268529538416917E-01,5.262367089839791E-01,5.361245933711475E-01,5.565166070031966E-01,5.874127498801266E-01,6.288130220019373E-01,6.807174233686291E-01,7.431259539802014E-01,8.160386138366549E-01,8.994554029379890E-01,9.933763212842039E-01};
      i += uussAmp.test_2to2_amp2([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH);
      i += uussAmp.test_2to2_amp2_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH);
      i += uussAmp.test_2to2_amp2_boosts([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH);
      i += uussAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH);
      //std::cout<<"########### mu=0.0042, ms=1.23, pspatial=1.25\n";
      pspatial = 1.25;
      ldouble dataCH2[20] = {1.043337948243725E+00,1.040351613956156E+00,1.037697099154046E+00,1.035374403837396E+00,1.033383528006204E+00,1.031724471660471E+00,1.030397234800197E+00,1.029401817425383E+00,1.028738219536027E+00,1.028406441132130E+00,1.028406482213692E+00,1.028738342780713E+00,1.029402022833193E+00,1.030397522371132E+00,1.031724841394529E+00,1.033383979903386E+00,1.035374937897702E+00,1.037697715377477E+00,1.040352312342711E+00,1.043338728793403E+00};
      i += uussAmp.test_2to2_amp2([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH2);
      i += uussAmp.test_2to2_amp2_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH2);
      i += uussAmp.test_2to2_amp2_boosts([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH2);
      i += uussAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH2);
      //std::cout<<"########### mu=1.23, ms=0.0042, pspatial=0.005\n";
      mu=1.23;
      ms=0.0042;
      pspatial = 0.005;
      uussAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH3[20] = {1.044955508723670E+00,1.044953955614376E+00,1.044952575173602E+00,1.044951367401348E+00,1.044950332297615E+00,1.044949469862403E+00,1.044948780095711E+00,1.044948262997540E+00,1.044947918567889E+00,1.044947746806759E+00,1.044947747714149E+00,1.044947921290060E+00,1.044948267534491E+00,1.044948786447443E+00,1.044949478028915E+00,1.044950342278908E+00,1.044951379197421E+00,1.044952588784456E+00,1.044953971040010E+00,1.044955525964085E+00};
      i += uussAmp.test_2to2_amp2([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH3);
      i += uussAmp.test_2to2_amp2_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH3);
      i += uussAmp.test_2to2_amp2_boosts([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH3);
      i += uussAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH3);
      //std::cout<<"########### mu=1.2, ms=1.23, pspatial=0.3\n";
      mu=1.2;
      ms=1.23;
      pspatial = 0.3;
      uussAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH4[20] = {1.531162081650840E+00,1.531100258200932E+00,1.531045304666243E+00,1.530997221046773E+00,1.530956007342523E+00,1.530921663553492E+00,1.530894189679680E+00,1.530873585721088E+00,1.530859851677715E+00,1.530852987549561E+00,1.530852993336627E+00,1.530859869038912E+00,1.530873614656416E+00,1.530894230189140E+00,1.530921715637083E+00,1.530956071000245E+00,1.530997296278627E+00,1.531045391472228E+00,1.531100356581048E+00,1.531162191605088E+00};
      i += uussAmp.test_2to2_amp2([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH4);
      i += uussAmp.test_2to2_amp2_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH4);
      i += uussAmp.test_2to2_amp2_boosts([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH4);
      i += uussAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH4);
      //std::cout<<"########### mu=1.2, ms=1.23, MW=2.11, pspatial=0.3\n";
      mu=1.2;
      ms=1.23;
      MW=2.11;
      pspatial = 0.3;
      uussAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH5[20] = {1.667596351160644E+00,1.665463668697492E+00,1.663345356659372E+00,1.661241415046284E+00,1.659151843858227E+00,1.657076643095202E+00,1.655015812757209E+00,1.652969352844247E+00,1.650937263356317E+00,1.648919544293418E+00,1.646916195655551E+00,1.644927217442715E+00,1.642952609654911E+00,1.640992372292139E+00,1.639046505354399E+00,1.637115008841689E+00,1.635197882754012E+00,1.633295127091366E+00,1.631406741853752E+00,1.629532727041169E+00};
      i += uussAmp.test_2to2_amp2([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH5);
      i += uussAmp.test_2to2_amp2_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH5);
      i += uussAmp.test_2to2_amp2_boosts([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH5);
      i += uussAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH5);
      //std::cout<<"########### mu=1.2, ms=1.23, MW=0.006, pspatial=0.3\n";
      mu=1.2;
      ms=1.23;
      MW=0.006;
      pspatial = 0.3;
      uussAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH6[20] = {2.006052580892816E+07,2.006052580857571E+07,2.006052580823017E+07,2.006052580789153E+07,2.006052580755979E+07,2.006052580723496E+07,2.006052580691703E+07,2.006052580660601E+07,2.006052580630189E+07,2.006052580600467E+07,2.006052580571436E+07,2.006052580543095E+07,2.006052580515444E+07,2.006052580488484E+07,2.006052580462214E+07,2.006052580436635E+07,2.006052580411746E+07,2.006052580387547E+07,2.006052580364039E+07,2.006052580341221E+07};
      i += uussAmp.test_2to2_amp2([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH6);
      i += uussAmp.test_2to2_amp2_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH6);
      i += uussAmp.test_2to2_amp2_boosts([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH6);
      i += uussAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH6);
      //std::cout<<"########### mu=1.2, ms=1.23, MW=2.11, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      ms=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      uussAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH7[20] = {1.668118836858382E+00,1.665931411705388E+00,1.663758356977426E+00,1.661599672674495E+00,1.659455358796595E+00,1.657325415343728E+00,1.655209842315891E+00,1.653108639713087E+00,1.651021807535314E+00,1.648949345782573E+00,1.646891254454863E+00,1.644847533552185E+00,1.642818183074538E+00,1.640803203021923E+00,1.638802593394340E+00,1.636816354191788E+00,1.634844485414268E+00,1.632886987061780E+00,1.630943859134323E+00,1.629015101631897E+00};
      i += uussAmp.test_2to2_amp2([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH7);
      i += uussAmp.test_2to2_amp2_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH7);
      i += uussAmp.test_2to2_amp2_boosts([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH7);
      i += uussAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH7);
      //std::cout<<"########### mu=1.2, ms=1.23, MW=0.006, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      ms=1.23;
      MW=0.006;
      pspatial = 0.3;
      uussAmp.set_masses(mu,ms,mh,MW);
      ldouble dataCH8[20] = {2.009770395603885E+07,2.009770277656351E+07,2.009770159709508E+07,2.009770041763355E+07,2.009769923817892E+07,2.009769805873120E+07,2.009769687929038E+07,2.009769569985646E+07,2.009769452042945E+07,2.009769334100934E+07,2.009769216159614E+07,2.009769098218984E+07,2.009768980279044E+07,2.009768862339795E+07,2.009768744401236E+07,2.009768626463368E+07,2.009768508526190E+07,2.009768390589702E+07,2.009768272653905E+07,2.009768154718798E+07};
      i += uussAmp.test_2to2_amp2([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH8);
      i += uussAmp.test_2to2_amp2_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH8);
      i += uussAmp.test_2to2_amp2_boosts([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH8);
      i += uussAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uussAmp.amp2(); }, mu,mu,ms,ms,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
