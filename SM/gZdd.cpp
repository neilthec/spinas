
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

//File:  SPINAS/SM/gZdd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/gZdd.h"

namespace spinas {

  gZdd::gZdd(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qd(-1.0/3.0), gs(gscharge), md(massd), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), propd(massd,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    p1=particle(0);
    p2=particle(MZ);
    p3=particle(md);
    p4=particle(md);
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
    gL=-2.0*Qd*SW*SW-1.0;
    gR=-2.0*Qd*SW*SW;
  }
  void gZdd::set_masses(const ldouble& massd, const ldouble& massW){
    md=massd;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(0);
    p2.set_mass(MZ);
    p3.set_mass(md);
    p4.set_mass(md);
    propd.set_mass(md);
    //Couplings
    preTU = 2.0*e*gs/(2.0*MW*SW);
  }
  void gZdd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenT=propd.den(propTP);
    pDenU=propd.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble gZdd::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
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
	//- preTU [1341] (gRe [24] <23> + gLe <24> [23])/((t-md^2) (u-md^2))
	//- preTU [12] ( gRe [14]<23>/(u-md^2) - gLe [13]<24>/(t-md^2) )
	//34 outgoing:
	//+ preTU [1341] (gRe [24] <23> + gLe <24> [23])/((t-md^2) (u-md^2))
	//+ preTU [12] ( gRe [14]<23>/(u-md^2) - gLe [13]<24>/(t-md^2) )
	
	amplitude += normFactor*preTU*s1341s.v()*(gR*s24s.v(ds2a,ds4)*a23a.v(ds2b,ds3) + gL*s23s.v(ds2a,ds3)*a24a.v(ds2b,ds4))/pDenT/pDenU;
	amplitude += normFactor*preTU*s12s.v(ds2a)*(gR*s14s.v(ds4)*a23a.v(ds2b,ds3)/pDenU - gL*s13s.v(ds3)*a24a.v(ds2b,ds4)/pDenT);

	
      }
      else if(ds1<0){
	//all ingoing:
	//preTU = 2.0*e*gs/(2.0*MW*SW);
	//- preTU <1341> (gRe [24] <23> + gLe <24> [23])/((t-md^2) (u-md^2))
	//- preTU <12> ( gLe <14>[23]/(u-md^2) - gRe <13>[24]/(t-md^2) )
	//34 outgoing:
	//+ preTU <1341> (gRe [24] <23> + gLe <24> [23])/((t-md^2) (u-md^2))
	//+ preTU <12> ( gLe <14>[23]/(u-md^2) - gRe <13>[24]/(t-md^2) )
	
	amplitude += normFactor*preTU*a1341a.v()*(gR*s24s.v(ds2a,ds4)*a23a.v(ds2b,ds3) + gL*s23s.v(ds2a,ds3)*a24a.v(ds2b,ds4))/pDenT/pDenU;
	amplitude += normFactor*preTU*a12a.v(ds2a)*(gL*a14a.v(ds4)*s23s.v(ds2b,ds3)/pDenU - gR*a13a.v(ds3)*s24s.v(ds2b,ds4)/pDenT);

      }


      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble gZdd::amp2(){
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
  ldouble gZdd::amp2_gplus(){
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
  ldouble gZdd::amp2_gminus(){
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
  int test_gZdd(){
    int n=0;//Number of fails
    std::cout<<"\t* g , Z  -> d , D       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# md=0.0075, MW=80.385, pspatial=250\n";
      ldouble md=0.0075;
      ldouble EE=0.31333, gs=1.238, MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      gZdd gZddAmp = gZdd(EE,gs,md,MW,SW);
      ldouble pspatial=250;
      ldouble dataCHp[20] = {1.484624381183095E-01,5.993387494601176E-02,4.784948589499014E-02,4.741573257738442E-02,5.157986607025876E-02,5.857964024455908E-02,6.792422231183459E-02,7.961139469447824E-02,9.392715645243195E-02,1.114163847554765E-01,1.329401780521202E-01,1.598196898917430E-01,1.941179750143844E-01,2.391947484003425E-01,3.008741282696622E-01,3.901745925141662E-01,5.307227266363457E-01,7.839632142952950E-01,1.375212154120943E+00,4.332356623589680E+00};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCHp);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCHp);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCHp);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCHp);
      ldouble dataCHm[20] = {4.332356623589690E+00,1.375212154120943E+00,7.839632142952954E-01,5.307227266363459E-01,3.901745925141665E-01,3.008741282696618E-01,2.391947484003427E-01,1.941179750143845E-01,1.598196898917431E-01,1.329401780521202E-01,1.114163847554765E-01,9.392715645243201E-02,7.961139469447827E-02,6.792422231183456E-02,5.857964024455916E-02,5.157986607025874E-02,4.741573257738440E-02,4.784948589499011E-02,5.993387494601172E-02,1.484624381183080E-01};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCHm);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCHm);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCHm);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCHm);
      ldouble dataCH[20] = {2.240409530854000E+00,7.175730145334774E-01,4.159063500951428E-01,2.890692296068652E-01,2.208772292922126E-01,1.797268842571104E-01,1.535594853560886E-01,1.368646848544314E-01,1.268734231720875E-01,1.221782814037984E-01,1.221782814037983E-01,1.268734231720875E-01,1.368646848544313E-01,1.535594853560885E-01,1.797268842571107E-01,2.208772292922125E-01,2.890692296068650E-01,4.159063500951426E-01,7.175730145334772E-01,2.240409530853994E+00};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH);
      //std::cout<<"\n# md=0.0042, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2p[20] = {2.470883955371162E-01,1.071859259572375E-01,8.581259768796699E-02,8.223006889492504E-02,8.541742480928935E-02,9.255846820328760E-02,1.027897874542004E-01,1.159750763438911E-01,1.323797085615549E-01,1.526031394668007E-01,1.776321138434232E-01,2.090029370113353E-01,2.491292199586657E-01,3.019519679575139E-01,3.743109660224891E-01,4.791527945484563E-01,6.442443084128243E-01,9.418035125153911E-01,1.636659165934345E+00,5.112335590333008E+00};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH2p);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH2p);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH2p);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH2p);
      ldouble dataCH2m[20] = {5.112335590333030E+00,1.636659165934347E+00,9.418035125153922E-01,6.442443084128249E-01,4.791527945484566E-01,3.743109660224892E-01,3.019519679575140E-01,2.491292199586659E-01,2.090029370113355E-01,1.776321138434233E-01,1.526031394668008E-01,1.323797085615549E-01,1.159750763438912E-01,1.027897874542004E-01,9.255846820328761E-02,8.541742480928934E-02,8.223006889492500E-02,8.581259768796692E-02,1.071859259572373E-01,2.470883955371151E-01};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH2m);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH2m);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH2m);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH2m);
      ldouble dataCH2[20] = {2.679711992935073E+00,8.719225459457923E-01,5.138080551016796E-01,3.632371886538749E-01,2.822851096788729E-01,2.334347171128884E-01,2.023708777058572E-01,1.825521481512785E-01,1.706913227864452E-01,1.651176266551120E-01,1.651176266551120E-01,1.706913227864451E-01,1.825521481512785E-01,2.023708777058572E-01,2.334347171128884E-01,2.822851096788728E-01,3.632371886538747E-01,5.138080551016790E-01,8.719225459457910E-01,2.679711992935061E+00};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH2);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH2);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH2);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH2);
      //std::cout<<"\n# md=125.1, MW=80.385, pspatial=95\n";
      md = 125;
      pspatial = 250;
      gZddAmp.set_masses(md,MW);
      ldouble dataCH4p[20] = {2.573059694119367E+00,1.544226905285537E+00,1.126397682171464E+00,9.099099414096501E-01,7.833775011311249E-01,7.050956380526152E-01,6.563043104993342E-01,6.275793046558726E-01,6.139268915433598E-01,6.127587540930778E-01,6.229994627285157E-01,6.447204346266215E-01,6.790838586976395E-01,7.285515164309091E-01,7.974378732498885E-01,8.930626377424133E-01,1.028130220739688E+00,1.225871372752823E+00,1.531487135995296E+00,2.029112397100374E+00};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH4p);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH4p);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH4p);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH4p);
      ldouble dataCH4m[20] = {2.029112397100374E+00,1.531487135995296E+00,1.225871372752823E+00,1.028130220739688E+00,8.930626377424136E-01,7.974378732498886E-01,7.285515164309092E-01,6.790838586976399E-01,6.447204346266217E-01,6.229994627285156E-01,6.127587540930778E-01,6.139268915433597E-01,6.275793046558724E-01,6.563043104993342E-01,7.050956380526152E-01,7.833775011311249E-01,9.099099414096499E-01,1.126397682171464E+00,1.544226905285536E+00,2.573059694119365E+00};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH4m);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH4m);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH4m);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH4m);
      ldouble dataCH4[20] = {2.301086045609871E+00,1.537857020640417E+00,1.176134527462144E+00,9.690200810746690E-01,8.382200694367692E-01,7.512667556512519E-01,6.924279134651217E-01,6.533315816767562E-01,6.293236630849907E-01,6.178791084107966E-01,6.178791084107967E-01,6.293236630849907E-01,6.533315816767560E-01,6.924279134651217E-01,7.512667556512519E-01,8.382200694367692E-01,9.690200810746689E-01,1.176134527462144E+00,1.537857020640416E+00,2.301086045609869E+00};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH4);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH4);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH4);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH4);
      //std::cout<<"\n# md=125, MW=80.385, pspatial=125.1\n";
      md = 125;
      pspatial = 125.1;
      gZddAmp.set_masses(md,MW);
      ldouble dataCH3p[20] = {9.110385017022786E-01,8.412257618359584E-01,7.857396354967533E-01,7.414174168354006E-01,7.059976256715692E-01,6.778338964713864E-01,6.557113423344177E-01,6.387255954187795E-01,6.262013895198760E-01,6.176368015691649E-01,6.126645815748365E-01,6.110251412695215E-01,6.125476699630410E-01,6.171370024698305E-01,6.247645465253966E-01,6.354619140646702E-01,6.493159153662903E-01,6.664631751994888E-01,6.870815384286206E-01,7.113730015801715E-01};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH3p);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH3p);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH3p);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2_gplus(); }, 0,MZ,md,md,pspatial,dataCH3p);
      ldouble dataCH3m[20] = {7.113730015801712E-01,6.870815384286204E-01,6.664631751994886E-01,6.493159153662902E-01,6.354619140646700E-01,6.247645465253964E-01,6.171370024698303E-01,6.125476699630408E-01,6.110251412695213E-01,6.126645815748363E-01,6.176368015691648E-01,6.262013895198759E-01,6.387255954187794E-01,6.557113423344175E-01,6.778338964713863E-01,7.059976256715690E-01,7.414174168354004E-01,7.857396354967531E-01,8.412257618359581E-01,9.110385017022784E-01};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH3m);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH3m);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH3m);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2_gminus(); }, 0,MZ,md,md,pspatial,dataCH3m);
      ldouble dataCH3[20] = {8.112057516412249E-01,7.641536501322894E-01,7.261014053481210E-01,6.953666661008454E-01,6.707297698681196E-01,6.512992214983914E-01,6.364241724021240E-01,6.256366326909102E-01,6.186132653946987E-01,6.151506915720006E-01,6.151506915720006E-01,6.186132653946986E-01,6.256366326909102E-01,6.364241724021240E-01,6.512992214983914E-01,6.707297698681196E-01,6.953666661008454E-01,7.261014053481210E-01,7.641536501322894E-01,8.112057516412249E-01};
      i += gZddAmp.test_2to2_amp2([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH3);
      i += gZddAmp.test_2to2_amp2_rotations([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH3);
      i += gZddAmp.test_2to2_amp2_boosts([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH3);
      i += gZddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gZddAmp.amp2(); }, 0,MZ,md,md,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
