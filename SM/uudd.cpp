
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

//File:  SPINAS/SM/uudd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uudd.h"

namespace spinas {
  //Constructors
  uudd::uudd(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& widthW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), gs(gscharge), mu(massu), md(massd), mh(massh), wh(widthh), MW(massW), WW(widthW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propAG(0,0), proph(mh,wh), propZ(MZ,WZ), propW(MW,WW),
    p1(particle(mu)), p2(particle(mu)),
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
    preh = e*e*mu*md/(4.0*MW*MW*SW*SW);
    gLu=1.0-4.0/3.0*SW*SW;
    gRu=-4.0/3.0*SW*SW;
    gLd=-1.0+2.0/3.0*SW*SW;
    gRd=2.0/3.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*md/MZ/MZ;//=-preh!
  }
  void uudd::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massh, const ldouble& massW){
    mu=massu;
    md=massd;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    propW.set_mass(MW);
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(mu);
    p2.set_mass(mu);
    p3.set_mass(md);
    p4.set_mass(md);
    preh = e*e*mu*md/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*md/MZ/MZ;//=-preh
  }
  void uudd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenTW = propW.denominator(propPT);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  //We have to separate the gluon so we can separate the color factor between the gluon^2 diagram
  // And the rest^2.
  // The cross term vanishes due to the trace of the adjoint rep.
  cdouble uudd::amp_gluon(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
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
  cdouble uudd::amp_rest(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
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
    //EE^2 Mu Md / (4 MW^2 SW^2) * ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //preh = e*e*mu*md/(4*MW*MW*SW*SW);
    amplitude += preh*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))*(s34s.v(ds3,ds4)+a34a.v(ds3,ds4))/pDenSh;
    
    //Z Boson
    //Defined above:
    //gLu=1.0-4.0/3.0*SW*SW;
    //gRu=-4.0/3.0*SW*SW;
    //gLd=-1.0+2.0/3.0*SW*SW;
    //gRd=2.0/3.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gLu-gRu)*(gLd-gRd)*mu*md/MZ/MZ; // = -preh
    //all in:
    //
    //+(EE^2 Mu Md (gLu-gRu)(gLd-gRd) (<12>-[12]) (<34>-[34]))/(8 CW^2 MZ^2 SW^2 (s-MZ^2))
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
  
  
  cdouble uudd::amp_W(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble sqrt2=std::sqrt(2.0);
    //W Boson
    //All ingoing uDTb: 2MW^2 [23]⟨14⟩ + mdmt⟨12⟩⟨34⟩ − mtmu[12]⟨34⟩ − mbmd⟨12⟩[34] + mumb[12][34]
    //All ingoing uTDb: 2MW^2 [32]⟨14⟩ + mdmt⟨13⟩⟨24⟩ − mtmu[13]⟨24⟩ − mbmd⟨13⟩[24] + mumb[13][24]
    //All ingoing uUDd: 2MW^2 [32]⟨14⟩ + mdmu⟨13⟩⟨24⟩ − mumu[13]⟨24⟩ − mdmd⟨13⟩[24] + mumd[13][24]
    //All ingoing uUDd: 2MW^2 [32]⟨14⟩ − mu^2[13]⟨24⟩ + mdmu( ⟨13⟩⟨24⟩ + [13][24] ) − md^2⟨13⟩[24]
    //Sign changes due to p3 and p4 being outgoing.
    // + 2MW^2 [23]⟨14⟩ + mu^2[13]⟨24⟩ + mdmu( ⟨13⟩⟨24⟩ + [13][24] ) + md^2⟨13⟩[24]
    //All ingoing:  e^2  ( 2 MW^2 [23] <14> + (Md <13>-Mu [13]) (Md [24]-Mu <24>) )/(4 MW^2 SW^2 (t-MW^2))
    //34 outgoing: -e^2  ( 2 MW^2 [23] <14> + (Md <13>+Mu [13]) (Md [24]+Mu <24>) )/(4 MW^2 SW^2 (t-MW^2))
    return + e*e/SW/SW*(   2.0*MW*MW*a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			   + (mu*s13s.v(ds1,ds3)+md*a13a.v(ds1,ds3))*(mu*a24a.v(ds2,ds4)+md*s24s.v(ds2,ds4))
			   )/(2.0*MW*MW*pDenTW);
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble uudd::amp2(){
    ldouble amp2 = 0;
    constexpr ldouble two=2, three = 3, four = 4, nine = 9;
    cdouble M_rest, M_W, M_gluon;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M_rest = amp_rest(j1,j2,j3,j4);
	    M_W = amp_W(j1,j2,j3,j4);
	    M_gluon = amp_gluon(j1,j2,j3,j4);
	    amp2 += nine*std::pow(std::abs(M_rest),2);// Color factor 3*3=9
	    amp2 += three*two*std::real(M_rest*std::conj(M_W));//Color factor 3
	    amp2 += nine*std::pow(std::abs(M_W),2);//Color factor 9
	    //Cross term with gluon and rest color factor 0 (Trace(Ta)*Trace(Ta)=0)
	    amp2 += four*two*std::real(M_gluon*std::conj(M_W));//Color factor 4
	    amp2 += two*std::pow(std::abs(M_gluon),2);//Color factor C^2*delta^ab*delta^ab = 1/4*8=2
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    //Average over colors 1/3*1/3 = 1/9
    return amp2/36.0;
  }

  



  //  Tests
  int test_uudd(){
    int n=0;//Number of fails
    std::cout<<"\t* u , U  -> d , D       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### mu=0.0042, md=1.23, pspatial=250\n";
      ldouble mu=0.0042, md=1.23, mh=125, wh=0, MW=80.385, WW=0, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      uudd uuddAmp = uudd(0.31333,1.238,mu,md,mh,wh,MW,WW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.328751719410853E+01,2.534780154017806E+00,1.000973453318062E+00,6.014755889223842E-01,4.683158035863678E-01,4.191650687429382E-01,4.024039968969574E-01,4.009869194169851E-01,4.089226713039439E-01,4.241872542957835E-01,4.462735288184110E-01,4.752617864778155E-01,5.114456629419802E-01,5.551771896352140E-01,6.068031598486040E-01,6.666408527286304E-01,7.349709158806232E-01,8.120375845373611E-01,8.980518048143988E-01,9.931952576741613E-01};
      i += uuddAmp.test_2to2_amp2([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH);
      i += uuddAmp.test_2to2_amp2_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH);
      i += uuddAmp.test_2to2_amp2_boosts([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH);
      i += uuddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH);
      //std::cout<<"########### mu=0.0042, md=1.23, pspatial=1.25\n";
      pspatial = 1.25;
      ldouble dataCH2[20] = {1.043171672809875E+00,1.040188283592273E+00,1.037536668595511E+00,1.035216827820757E+00,1.033228761269182E+00,1.031572468941956E+00,1.030247950840247E+00,1.029255206965227E+00,1.028594237318064E+00,1.028265041899928E+00,1.028267620711989E+00,1.028601973755416E+00,1.029268101031380E+00,1.030266002541048E+00,1.031595678285591E+00,1.033257128266178E+00,1.035250352483979E+00,1.037575350940162E+00,1.040232123635897E+00,1.043220670572354E+00};
      i += uuddAmp.test_2to2_amp2([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH2);
      i += uuddAmp.test_2to2_amp2_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH2);
      i += uuddAmp.test_2to2_amp2_boosts([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH2);
      i += uuddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH2);
      //std::cout<<"########### mu=1.23, md=0.0042, pspatial=0.005\n";
      mu=1.23;
      md=0.0042;
      pspatial = 0.005;
      uuddAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH3[20] = {1.044817108036189E+00,1.044815611185082E+00,1.044814286979695E+00,1.044813135420026E+00,1.044812156506077E+00,1.044811350237846E+00,1.044810716615335E+00,1.044810255638543E+00,1.044809967307470E+00,1.044809851622115E+00,1.044809908582480E+00,1.044810138188564E+00,1.044810540440367E+00,1.044811115337890E+00,1.044811862881131E+00,1.044812783070092E+00,1.044813875904771E+00,1.044815141385170E+00,1.044816579511288E+00,1.044818190283125E+00};
      i += uuddAmp.test_2to2_amp2([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH3);
      i += uuddAmp.test_2to2_amp2_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH3);
      i += uuddAmp.test_2to2_amp2_boosts([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH3);
      i += uuddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH3);
      //std::cout<<"########### mu=1.2, md=1.23, pspatial=0.3\n";
      mu=1.2;
      md=1.23;
      pspatial = 0.3;
      uuddAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH4[20] = {1.530954360792078E+00,1.530892903313619E+00,1.530838314832607E+00,1.530790595349045E+00,1.530749744862937E+00,1.530715763374286E+00,1.530688650883095E+00,1.530668407389368E+00,1.530655032893108E+00,1.530648527394319E+00,1.530648890893003E+00,1.530656123389164E+00,1.530670224882806E+00,1.530691195373931E+00,1.530719034862544E+00,1.530753743348647E+00,1.530795320832244E+00,1.530843767313338E+00,1.530899082791932E+00,1.530961267268030E+00};
      i += uuddAmp.test_2to2_amp2([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH4);
      i += uuddAmp.test_2to2_amp2_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH4);
      i += uuddAmp.test_2to2_amp2_boosts([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH4);
      i += uuddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH4);
      //std::cout<<"########### mu=1.2, md=1.23, MW=2.11, pspatial=0.3\n";
      mu=1.2;
      md=1.23;
      MW=2.11;
      pspatial = 0.3;
      uuddAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH5[20] = {1.418851612066934E+00,1.417107565483293E+00,1.415376383543563E+00,1.413658069914767E+00,1.411952628267551E+00,1.410260062275891E+00,1.408580375616822E+00,1.406913571970164E+00,1.405259655018254E+00,1.403618628445683E+00,1.401990495939045E+00,1.400375261186679E+00,1.398772927878427E+00,1.397183499705391E+00,1.395606980359695E+00,1.394043373534257E+00,1.392492682922556E+00,1.390954912218415E+00,1.389430065115778E+00,1.387918145308500E+00};
      i += uuddAmp.test_2to2_amp2([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH5);
      i += uuddAmp.test_2to2_amp2_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH5);
      i += uuddAmp.test_2to2_amp2_boosts([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH5);
      i += uuddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH5);
      //std::cout<<"########### mu=1.2, md=1.23, MW=0.006, pspatial=0.3\n";
      mu=1.2;
      md=1.23;
      MW=0.006;
      pspatial = 0.3;
      uuddAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH6[20] = {9.609735836452958E+07,9.428398968880223E+07,9.307160626655106E+07,9.220398222347404E+07,9.155239456529854E+07,9.104509966299106E+07,9.063895093355531E+07,9.030644303301983E+07,9.002921375933129E+07,8.979453672059040E+07,8.959331336021218E+07,8.941886632731687E+07,8.926618437819876E+07,8.913143326291876E+07,8.901162935371499E+07,8.890441616328084E+07,8.880790780082829E+07,8.872057708792524E+07,8.864117414495783E+07,8.856866618738356E+07};
      i += uuddAmp.test_2to2_amp2([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH6);
      i += uuddAmp.test_2to2_amp2_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH6);
      i += uuddAmp.test_2to2_amp2_boosts([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH6);
      i += uuddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH6);
      //std::cout<<"########### mu=1.2, md=1.23, MW=2.11, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      md=1.23;
      MW=2.11;
      mh=3.125;
      pspatial = 0.3;
      uuddAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH7[20] = {1.419444879910220E+00,1.417638598456531E+00,1.415845207770816E+00,1.414064711383704E+00,1.412297112830388E+00,1.410542415650335E+00,1.408800623387001E+00,1.407071739587550E+00,1.405355767802578E+00,1.403652711585846E+00,1.401962574494013E+00,1.400285360086380E+00,1.398621071924633E+00,1.396969713572593E+00,1.395331288595976E+00,1.393705800562153E+00,1.392093253039909E+00,1.390493649599222E+00,1.388906993811030E+00,1.387333289247017E+00};
      i += uuddAmp.test_2to2_amp2([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH7);
      i += uuddAmp.test_2to2_amp2_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH7);
      i += uuddAmp.test_2to2_amp2_boosts([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH7);
      i += uuddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH7);
      //std::cout<<"########### mu=1.2, md=1.23, MW=0.006, Mh=3.125, pspatial=0.3\n";
      mu=1.2;
      md=1.23;
      MW=0.006;
      pspatial = 0.3;
      uuddAmp.set_masses(mu,md,mh,MW);
      ldouble dataCH8[20] = {9.572466267949896E+07,9.414310242646402E+07,9.308720128151496E+07,9.233231743461339E+07,9.176581974529694E+07,9.132502545534734E+07,9.097227976250891E+07,9.068359981008488E+07,9.044298667699768E+07,9.023935864722821E+07,9.006479707249156E+07,8.991349308243747E+07,8.978108887459169E+07,8.966425127698977E+07,8.956038733878280E+07,8.946744966659765E+07,8.938380011967003E+07,8.930811242448087E+07,8.923930133343299E+07,8.917647025379072E+07};
      i += uuddAmp.test_2to2_amp2([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH8);
      i += uuddAmp.test_2to2_amp2_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH8);
      i += uuddAmp.test_2to2_amp2_boosts([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH8);
      i += uuddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uuddAmp.amp2(); }, mu,mu,md,md,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
