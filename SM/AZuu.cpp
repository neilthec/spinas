
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

//File:  SPINAS/SM/AZuu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AZuu.h"

namespace spinas {

  AZuu::AZuu(const ldouble& echarge, const ldouble& massu, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qu(2.0/3.0), mu(massu), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), prope(massu,0) {
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
    preTU = 2.0*e*e*Qu/(2.0*MW*SW);
    gL=-2.0*Qu*SW*SW+1.0;
    gR=-2.0*Qu*SW*SW;
  }
  void AZuu::set_masses(const ldouble& massu, const ldouble& massW){
    mu=massu;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(0);
    p2.set_mass(MZ);
    p3.set_mass(mu);
    p4.set_mass(mu);
    prope.set_mass(mu);
    //Couplings
    preTU = 2.0*e*e*Qu/(2.0*MW*SW);
  }
  void AZuu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    pDenT=prope.den(propTP);
    pDenU=prope.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AZuu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
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
	//preTU = 2.0*e*e/(2.0*MW*SW);
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
	//preTU = 2.0*e*e/(2.0*MW*SW);
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
  ldouble AZuu::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2*1/3=1/6
    return amp2/6.0;
  }
  //set_momenta(...) must be called before amp2_Aplus().
  ldouble AZuu::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }  
  //set_momenta(...) must be called before amp2_Aminus().
  ldouble AZuu::amp2_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(-2,j2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }



  //  Tests
  int test_AZuu(){
    int n=0;//Number of fails
    std::cout<<"\t* A , Z  -> U , u       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, MW=80.385, pspatial=250\n";
      ldouble mu=0.0042;
      ldouble EE=0.31333, MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      AZuu AZuuAmp = AZuu(EE,mu,MW,SW);
      ldouble pspatial=250;
      ldouble dataCHp[20] = {9.348855091872113E-02,3.115372503126614E-02,1.933854990906865E-02,1.482511575092981E-02,1.282842290330164E-02,1.206216610758889E-02,1.205328983815461E-02,1.260733659880059E-02,1.365465591392116E-02,1.519890337476879E-02,1.730307563178806E-02,2.009632903824818E-02,2.380153300558408E-02,2.879674600119634E-02,3.574894533479064E-02,4.592954130099921E-02,6.207334731354482E-02,9.130124102723922E-02,1.597357105309098E-01,5.025072628557495E-01};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCHp);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCHp);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCHp);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCHp);
      ldouble dataCHm[20] = {5.025072628557463E-01,1.597357105309094E-01,9.130124102723911E-02,6.207334731354479E-02,4.592954130099917E-02,3.574894533479066E-02,2.879674600119633E-02,2.380153300558409E-02,2.009632903824822E-02,1.730307563178809E-02,1.519890337476877E-02,1.365465591392114E-02,1.260733659880059E-02,1.205328983815462E-02,1.206216610758888E-02,1.282842290330165E-02,1.482511575092982E-02,1.933854990906867E-02,3.115372503126621E-02,9.348855091872182E-02};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCHm);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCHm);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCHm);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCHm);
      ldouble dataCH[20] = {2.979979068872337E-01,9.544471778108775E-02,5.531989546815388E-02,3.844923153223730E-02,2.937898210215041E-02,2.390555572118978E-02,2.042501791967547E-02,1.820443480219234E-02,1.687549247608469E-02,1.625098950327844E-02,1.625098950327841E-02,1.687549247608466E-02,1.820443480219234E-02,2.042501791967548E-02,2.390555572118976E-02,2.937898210215043E-02,3.844923153223732E-02,5.531989546815394E-02,9.544471778108798E-02,2.979979068872358E-01};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2p[20] = {1.186157379534263E-01,4.121394664329283E-02,2.650088667757490E-02,2.084300871377262E-02,1.830103084507182E-02,1.727715657053672E-02,1.718236761555625E-02,1.777275164452862E-02,1.895837710167417E-02,2.073895181228521E-02,2.318579159301541E-02,2.644908609753656E-02,3.078993927590227E-02,3.665252005205577E-02,4.482136236056360E-02,5.679271466557612E-02,7.578568368891815E-02,1.101828069071335E-01,1.907356955043479E-01,5.942437365135418E-01};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH2p);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH2p);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH2p);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH2p);
      ldouble dataCH2m[20] = {5.942437365135488E-01,1.907356955043487E-01,1.101828069071337E-01,7.578568368891828E-02,5.679271466557618E-02,4.482136236056363E-02,3.665252005205578E-02,3.078993927590228E-02,2.644908609753658E-02,2.318579159301542E-02,2.073895181228522E-02,1.895837710167416E-02,1.777275164452862E-02,1.718236761555625E-02,1.727715657053671E-02,1.830103084507180E-02,2.084300871377259E-02,2.650088667757484E-02,4.121394664329267E-02,1.186157379534249E-01};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH2m);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH2m);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH2m);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH2m);
      ldouble dataCH2[20] = {3.564297372334875E-01,1.159748210738207E-01,6.834184679235432E-02,4.831434620134545E-02,3.754687275532400E-02,3.104925946555017E-02,2.691744383380601E-02,2.428134546021545E-02,2.270373159960537E-02,2.196237170265031E-02,2.196237170265031E-02,2.270373159960536E-02,2.428134546021544E-02,2.691744383380601E-02,3.104925946555016E-02,3.754687275532396E-02,4.831434620134537E-02,6.834184679235417E-02,1.159748210738203E-01,3.564297372334833E-01};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH2);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH2);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH2);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, MW=80.385, pspatial=95\n";
      mu = 125;
      pspatial = 250;
      AZuuAmp.set_masses(mu,MW);
      ldouble dataCH4p[20] = {3.954681535388729E-01,2.464904900176182E-01,1.831983614537253E-01,1.493320003215018E-01,1.289973667582184E-01,1.160745618231079E-01,1.077500669308662E-01,1.025850825788866E-01,9.980519645697623E-02,9.900260834926040E-02,1.000037429854912E-01,1.028155821295299E-01,1.076201826691676E-01,1.148129746654643E-01,1.251019654270979E-01,1.397202174173295E-01,1.608892489998421E-01,1.929229381973783E-01,2.452450464894317E-01,3.422917085567956E-01};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH4p);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH4p);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH4p);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH4p);
      ldouble dataCH4m[20] = {3.422917085567957E-01,2.452450464894317E-01,1.929229381973784E-01,1.608892489998421E-01,1.397202174173295E-01,1.251019654270979E-01,1.148129746654643E-01,1.076201826691677E-01,1.028155821295300E-01,1.000037429854912E-01,9.900260834926042E-02,9.980519645697622E-02,1.025850825788866E-01,1.077500669308662E-01,1.160745618231079E-01,1.289973667582184E-01,1.493320003215017E-01,1.831983614537253E-01,2.464904900176182E-01,3.954681535388727E-01};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH4m);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH4m);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH4m);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH4m);
      ldouble dataCH4[20] = {3.688799310478343E-01,2.458677682535249E-01,1.880606498255519E-01,1.551106246606719E-01,1.343587920877739E-01,1.205882636251029E-01,1.112815207981653E-01,1.051026326240271E-01,1.013103892932531E-01,9.950317566737582E-02,9.950317566737582E-02,1.013103892932531E-01,1.051026326240271E-01,1.112815207981653E-01,1.205882636251029E-01,1.343587920877739E-01,1.551106246606719E-01,1.880606498255518E-01,2.458677682535249E-01,3.688799310478342E-01};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH4);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH4);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH4);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, MW=80.385, pspatial=125.1\n";
      mu = 125;
      pspatial = 125.1;
      AZuuAmp.set_masses(mu,MW);
      ldouble dataCH3p[20] = {1.421517131657403E-01,1.318253521339229E-01,1.236019671247538E-01,1.170303658320564E-01,1.117879519454715E-01,1.076404796082803E-01,1.044162251064072E-01,1.019890252755024E-01,1.002669578678783E-01,9.918473721021977E-02,9.869865154637657E-02,9.878332345229784E-02,9.942986371633868E-02,1.006451864915071E-01,1.024524046181053E-01,1.048923603907534E-01,1.080264964269146E-01,1.119414658897731E-01,1.167561682989325E-01,1.226323561613953E-01};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH3p);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH3p);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH3p);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2_Aplus(); }, 0,MZ,mu,mu,pspatial,dataCH3p);
      ldouble dataCH3m[20] = {1.226323561613953E-01,1.167561682989325E-01,1.119414658897731E-01,1.080264964269146E-01,1.048923603907534E-01,1.024524046181053E-01,1.006451864915071E-01,9.942986371633866E-02,9.878332345229783E-02,9.869865154637654E-02,9.918473721021977E-02,1.002669578678782E-01,1.019890252755024E-01,1.044162251064072E-01,1.076404796082803E-01,1.117879519454715E-01,1.170303658320563E-01,1.236019671247538E-01,1.318253521339229E-01,1.421517131657402E-01};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH3m);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH3m);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH3m);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2_Aminus(); }, 0,MZ,mu,mu,pspatial,dataCH3m);
      ldouble dataCH3[20] = {1.323920346635678E-01,1.242907602164277E-01,1.177717165072634E-01,1.125284311294855E-01,1.083401561681124E-01,1.050464421131928E-01,1.025307057989572E-01,1.007094444959205E-01,9.952514066008805E-02,9.894169437829815E-02,9.894169437829817E-02,9.952514066008805E-02,1.007094444959205E-01,1.025307057989572E-01,1.050464421131928E-01,1.083401561681124E-01,1.125284311294855E-01,1.177717165072634E-01,1.242907602164277E-01,1.323920346635678E-01};
      i += AZuuAmp.test_2to2_amp2([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH3);
      i += AZuuAmp.test_2to2_amp2_rotations([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH3);
      i += AZuuAmp.test_2to2_amp2_boosts([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH3);
      i += AZuuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZuuAmp.amp2(); }, 0,MZ,mu,mu,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
