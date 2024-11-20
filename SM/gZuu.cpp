
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

//File:  SPINAS/SM/gZuu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/gZuu.h"

namespace spinas {

  gZuu::gZuu(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qu(2.0/3.0), gs(gscharge), mu(massu), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), prope(massu,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    p1=particle(0);
    p2=particle(MZ);
    p3=particle(mu);
    p4=particle(mu);
    //<12>,[12],<23>,[23],<13>,[13],<24>,[24],<14>,[14]
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    //[1341], <1341>
    s1341s = sproduct(SQUARE,&p1,&p3,&p4,&p1);
    a1341a = sproduct(ANGLE,&p1,&p3,&p4,&p1);
    //Couplings
    preTU = 2.0*e*gs/(2.0*MW*SW);
    gL=-2.0*Qu*SW*SW+1.0;
    gR=-2.0*Qu*SW*SW;
  }
  void gZuu::set_masses(const ldouble& massu, const ldouble& massW){
    mu=massu;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(0);
    p2.set_mass(MZ);
    p3.set_mass(mu);
    p4.set_mass(mu);
    prope.set_mass(mu);
    //Couplings
    preTU = 2.0*e*gs/(2.0*MW*SW);
  }
  void gZuu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //<12>,[12],<23>,[23],<13>,[13],<24>,[24],<14>,[14]
    s12s.update();
    a12a.update();
    s23s.update();
    a23a.update();
    s13s.update();
    a13a.update();
    s24s.update();
    a24a.update();
    s14s.update();
    a14a.update();
    //[1341], <1341>
    s1341s.update();
    a1341a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT=prope.denominator(propTP);
    pDenU=prope.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble gZuu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds2a, ds2b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds2);
    ldouble normFactor=get_spin_normalization(ds2);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds2,ds2a,ds2b, i);
      
      if(ds1>0){
	//all ingoing:
	//preTU = 2.0*e*gs/(2.0*MW*SW);
	//- preTU [1341] (gRe [24] <23> + gLe <24> [23])/((t-mu^2) (u-mu^2))
	//- preTU [12] ( gRe [14]<23>/(u-mu^2) - gLe [13]<24>/(t-mu^2) )
	//34 outgoing:
	//+ preTU [1341] (gRe [24] <23> + gLe <24> [23])/((t-mu^2) (u-mu^2))
	//+ preTU [12] ( gRe [14]<23>/(u-mu^2) - gLe [13]<24>/(t-mu^2) )
	
	amplitude += normFactor*preTU*s1341s.v()*(gR*s24s.v(ds2a,ds4)*a23a.v(ds2b,ds3) + gL*s23s.v(ds2a,ds3)*a24a.v(ds2b,ds4))/pDenT/pDenU;
	amplitude += normFactor*preTU*s12s.v(ds2a)*(gR*s14s.v(ds4)*a23a.v(ds2b,ds3)/pDenU - gL*s13s.v(ds3)*a24a.v(ds2b,ds4)/pDenT);

	
      }
      else if(ds1<0){
	//all ingoing:
	//preTU = 2.0*e*gs/(2.0*MW*SW);
	//- preTU <1341> (gRe [24] <23> + gLe <24> [23])/((t-mu^2) (u-mu^2))
	//- preTU <12> ( gLe <14>[23]/(u-mu^2) - gRe <13>[24]/(t-mu^2) )
	//34 outgoing:
	//+ preTU <1341> (gRe [24] <23> + gLe <24> [23])/((t-mu^2) (u-mu^2))
	//+ preTU <12> ( gLe <14>[23]/(u-mu^2) - gRe <13>[24]/(t-mu^2) )
	
	amplitude += normFactor*preTU*a1341a.v()*(gR*s24s.v(ds2a,ds4)*a23a.v(ds2b,ds3) + gL*s23s.v(ds2a,ds3)*a24a.v(ds2b,ds4))/pDenT/pDenU;
	amplitude += normFactor*preTU*a12a.v(ds2a)*(gL*a14a.v(ds4)*s23s.v(ds2b,ds3)/pDenU - gR*a13a.v(ds3)*s24s.v(ds2b,ds4)/pDenT);

      }


      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble gZuu::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4 // It's the same for both diagrams
	  }
    //Average over initial spins 1/2*1/3=1/6
    //Average over initial colors 1/8
    return amp2/48.0;
  }
  //set_momenta(...) must be called before amp2_gplus().
  ldouble gZuu::amp2_gplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
	}
    //Average over initial spins 1/3
    //Average over initial colors 1/8
    return amp2/24.0;
  }  
  //set_momenta(...) must be called before amp2_gminus().
  ldouble gZuu::amp2_gminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(-2,j2,j3,j4);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
	}
    //Average over initial spins 1/3
    //Average over initial colors 1/8
    return amp2/24.0;
  }



  //  Tests
  int test_gZuu(){
    int n=0;//Number of fails
    std::cout<<"\t* g , Z  -> u , U       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, MW=80.385, pspatial=250\n";
      ldouble mu=0.0042;
      ldouble EE=0.31333, gs=1.238, MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      gZuu gZuuAmp = gZuu(EE,gs,mu,MW,SW);
      ldouble pspatial=250;
      ldouble dataCHp[20] = {5.473019908480842E-01,1.823805756361486E-01,1.132120111108584E-01,8.678940132563741E-02,7.510033394915083E-02,7.061450262891711E-02,7.056253904744168E-02,7.380604739305831E-02,7.993727887099061E-02,8.897763409500940E-02,1.012959089429656E-01,1.176482123563313E-01,1.393392695809126E-01,1.685823158185813E-01,2.092819790250086E-01,2.688813672449640E-01,3.633906636636570E-01,5.344970104901092E-01,9.351270452261717E-01,2.941785092119764E+00};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCHp);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCHp);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCHp);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCHp);
      ldouble dataCHm[20] = {2.941785092119745E+00,9.351270452261694E-01,5.344970104901087E-01,3.633906636636567E-01,2.688813672449638E-01,2.092819790250088E-01,1.685823158185812E-01,1.393392695809126E-01,1.176482123563315E-01,1.012959089429657E-01,8.897763409500928E-02,7.993727887099050E-02,7.380604739305834E-02,7.056253904744172E-02,7.061450262891708E-02,7.510033394915090E-02,8.678940132563748E-02,1.132120111108585E-01,1.823805756361490E-01,5.473019908480882E-01};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCHm);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCHm);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCHm);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCHm);
      ldouble dataCH[20] = {1.744543541483915E+00,5.587538104311590E-01,3.238545108004836E-01,2.250900324946471E-01,1.719908505970573E-01,1.399482408269629E-01,1.195724274330115E-01,1.065726584869855E-01,9.879274561366105E-02,9.513677151898757E-02,9.513677151898742E-02,9.879274561366087E-02,1.065726584869854E-01,1.195724274330115E-01,1.399482408269629E-01,1.719908505970575E-01,2.250900324946472E-01,3.238545108004839E-01,5.587538104311605E-01,1.744543541483927E+00};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2p[20] = {6.944019229078121E-01,2.412752666173151E-01,1.551418638469000E-01,1.220194350239670E-01,1.071381524010151E-01,1.011441732097194E-01,1.005892583751070E-01,1.040454928684880E-01,1.109863981100301E-01,1.214102636464466E-01,1.357345875355204E-01,1.548386121620309E-01,1.802508959460022E-01,2.145716988546501E-01,2.623938504917775E-01,3.324767096792454E-01,4.436656163040619E-01,6.450337392643064E-01,1.116607593651550E+00,3.478830047602526E+00};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH2p);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH2p);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH2p);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH2p);
      ldouble dataCH2m[20] = {3.478830047602567E+00,1.116607593651554E+00,6.450337392643077E-01,4.436656163040626E-01,3.324767096792458E-01,2.623938504917777E-01,2.145716988546501E-01,1.802508959460023E-01,1.548386121620311E-01,1.357345875355204E-01,1.214102636464466E-01,1.109863981100301E-01,1.040454928684880E-01,1.005892583751069E-01,1.011441732097194E-01,1.071381524010150E-01,1.220194350239668E-01,1.551418638468997E-01,2.412752666173142E-01,6.944019229078040E-01};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH2m);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH2m);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH2m);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH2m);
      ldouble dataCH2[20] = {2.086615985255189E+00,6.789414301344346E-01,4.000878015556039E-01,2.828425256640148E-01,2.198074310401304E-01,1.817690118507486E-01,1.575804786148785E-01,1.421481944072451E-01,1.329125051360306E-01,1.285724255909835E-01,1.285724255909835E-01,1.329125051360305E-01,1.421481944072451E-01,1.575804786148785E-01,1.817690118507485E-01,2.198074310401302E-01,2.828425256640144E-01,4.000878015556030E-01,6.789414301344319E-01,2.086615985255165E+00};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH2);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH2);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH2);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, MW=80.385, pspatial=95\n";
      mu = 125;
      pspatial = 250;
      gZuuAmp.set_masses(mu,MW);
      ldouble dataCH4p[20] = {2.315155231542889E+00,1.443008096564127E+00,1.072482426547620E+00,8.742214984628480E-01,7.551782004013885E-01,6.795253338329851E-01,6.307919586490055E-01,6.005550345469810E-01,5.842809860790007E-01,5.795824634805944E-01,5.854433199612180E-01,6.019044282611521E-01,6.300315883757593E-01,6.721397325350570E-01,7.323736870922876E-01,8.179520636779771E-01,9.418801063696202E-01,1.129412180615404E+00,1.435717003529024E+00,2.003848938751837E+00};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH4p);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH4p);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH4p);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH4p);
      ldouble dataCH4m[20] = {2.003848938751838E+00,1.435717003529025E+00,1.129412180615404E+00,9.418801063696203E-01,8.179520636779772E-01,7.323736870922877E-01,6.721397325350571E-01,6.300315883757596E-01,6.019044282611522E-01,5.854433199612180E-01,5.795824634805945E-01,5.842809860790006E-01,6.005550345469809E-01,6.307919586490055E-01,6.795253338329851E-01,7.551782004013885E-01,8.742214984628479E-01,1.072482426547620E+00,1.443008096564127E+00,2.315155231542887E+00};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH4m);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH4m);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH4m);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH4m);
      ldouble dataCH4[20] = {2.159502085147364E+00,1.439362550046576E+00,1.100947303581512E+00,9.080508024162341E-01,7.865651320396828E-01,7.059495104626363E-01,6.514658455920312E-01,6.152933114613702E-01,5.930927071700764E-01,5.825128917209063E-01,5.825128917209063E-01,5.930927071700763E-01,6.152933114613701E-01,6.514658455920312E-01,7.059495104626363E-01,7.865651320396828E-01,9.080508024162340E-01,1.100947303581512E+00,1.439362550046576E+00,2.159502085147362E+00};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH4);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH4);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH4);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, MW=80.385, pspatial=125.1\n";
      mu = 125;
      pspatial = 125.1;
      gZuuAmp.set_masses(mu,MW);
      ldouble dataCH3p[20] = {8.321865603170451E-01,7.717338321980346E-01,7.235923759148859E-01,6.851208151172424E-01,6.544306019437387E-01,6.301504110024373E-01,6.112749349090042E-01,5.970655874906741E-01,5.869842362310244E-01,5.806486848222924E-01,5.778030352862824E-01,5.782987227498725E-01,5.820837078651155E-01,5.891984675638831E-01,5.997783093612370E-01,6.140623327932887E-01,6.324102360961334E-01,6.553293054375197E-01,6.835156040588876E-01,7.179160657637574E-01};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH3p);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH3p);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH3p);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2_gplus(); }, 0,MZ,mu,mu,pspatial,dataCH3p);
      ldouble dataCH3m[20] = {7.179160657637573E-01,6.835156040588874E-01,6.553293054375194E-01,6.324102360961333E-01,6.140623327932886E-01,5.997783093612369E-01,5.891984675638829E-01,5.820837078651154E-01,5.782987227498725E-01,5.778030352862823E-01,5.806486848222923E-01,5.869842362310242E-01,5.970655874906740E-01,6.112749349090041E-01,6.301504110024372E-01,6.544306019437386E-01,6.851208151172423E-01,7.235923759148858E-01,7.717338321980345E-01,8.321865603170450E-01};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH3m);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH3m);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH3m);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2_gminus(); }, 0,MZ,mu,mu,pspatial,dataCH3m);
      ldouble dataCH3[20] = {7.750513130404012E-01,7.276247181284611E-01,6.894608406762027E-01,6.587655256066879E-01,6.342464673685136E-01,6.149643601818371E-01,6.002367012364436E-01,5.895746476778947E-01,5.826414794904484E-01,5.792258600542873E-01,5.792258600542873E-01,5.826414794904484E-01,5.895746476778947E-01,6.002367012364436E-01,6.149643601818371E-01,6.342464673685136E-01,6.587655256066879E-01,6.894608406762027E-01,7.276247181284611E-01,7.750513130404012E-01};
      i += gZuuAmp.test_2to2_amp2([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH3);
      i += gZuuAmp.test_2to2_amp2_rotations([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH3);
      i += gZuuAmp.test_2to2_amp2_boosts([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH3);
      i += gZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
