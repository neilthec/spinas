
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

//File:  SPINAS/SM/AuZu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AuZu.h"

namespace spinas {

  AuZu::AuZu(const ldouble& echarge, const ldouble& massu, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qu(2.0/3.0), mu(massu), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), prope(massu,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    p1=particle(0);
    p2=particle(mu);
    p3=particle(MZ);
    p4=particle(mu);
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
    //[1341], <1341>
    s1421s = sproduct(SQUARE,&p1,&p4,&p2,&p1);
    a1421a = sproduct(ANGLE,&p1,&p4,&p2,&p1);
    //Couplings
    preTU = 2.0*e*e*Qu/(2.0*MW*SW);
    gL=-2.0*Qu*SW*SW+1.0;
    gR=-2.0*Qu*SW*SW;
  }
  void AuZu::set_masses(const ldouble& massu, const ldouble& massW){
    mu=massu;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(0);
    p2.set_mass(mu);
    p3.set_mass(MZ);
    p4.set_mass(mu);
    prope.set_mass(mu);
    //Couplings
    preTU = 2.0*e*e*Qu/(2.0*MW*SW);
  }
  void AuZu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    //[1341], <1341>
    s1421s.update();
    a1421a.update();
    //Propagator Momentum
    ldouble propSP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=prope.den(propSP);
    pDenU=prope.den(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble AuZu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      if(ds1>0){
	//preTU = 2.0*e*e/(2.0*MW*SW);
	//AZUu all in:
	//- preTU [1341] (gRe [24] <23> + gLe <24> [23])/((t-mu^2) (u-mu^2))
	//- preTU [12] ( gRe [14]<23>/(u-mu^2) - gLe [13]<24>/(t-mu^2) )
	//AuZU 4->2->3->4:
	//  preTU [1421] (gRe [23] <34> + gLe <23> [34])/((u-mu^2) (s-mu^2))
	//- preTU [13] ( gRe [12]<34>/(s-mu^2) + gLe [14]<23>/(u-mu^2) )
	//34 out:
	//- preTU [1421] (gRe [23] <34> - gLe <23> [34])/((u-mu^2) (s-mu^2))
	//- preTU [13] ( gRe [12]<34>/(s-mu^2) - gLe [14]<23>/(u-mu^2) )
	
	amplitude += - normFactor*preTU*s1421s.v()*(gR*s23s.v(ds2,ds3a)*a34a.v(ds3b,ds4) - gL*s34s.v(ds3a,ds4)*a23a.v(ds2,ds3b))/pDenU/pDenS;
	amplitude += - normFactor*preTU*s13s.v(ds3a)*(gR*s12s.v(ds2)*a34a.v(ds3b,ds4)/pDenS - gL*s14s.v(ds4)*a23a.v(ds2,ds3b)/pDenU);

	
      }
      else if(ds1<0){
	//preTU = 2.0*e*e/(2.0*MW*SW);
	//AZUu all in:
	//- preTU <1341> (gRe [24] <23> + gLe <24> [23])/((t-mu^2) (u-mu^2))
	//- preTU <12> ( gLe <14>[23]/(u-mu^2) - gRe <13>[24]/(t-mu^2) )
	//AuZU: 4->2->3->4
	//  preTU <1421> (gRe [23] <34> + gLe <23> [34])/((u-mu^2) (s-mu^2))
	//- preTU <13> ( gLe <12>[34]/(s-mu^2) + gRe <14>[23]/(u-mu^2) )
	//34 out:
	//- preTU <1421> (gRe [23] <34> - gLe <23> [34])/((u-mu^2) (s-mu^2))
	//  preTU <13> ( gLe <12>[34]/(s-mu^2) - gRe <14>[23]/(u-mu^2) )
	
	amplitude += - normFactor*preTU*a1421a.v()*(gR*s23s.v(ds2,ds3a)*a34a.v(ds3b,ds4) - gL*s34s.v(ds3a,ds4)*a23a.v(ds2,ds3b))/pDenU/pDenS;
	amplitude += + normFactor*preTU*a13a.v(ds3a)*(gL*a12a.v(ds2)*s34s.v(ds3b,ds4)/pDenS - gR*a14a.v(ds4)*s23s.v(ds2,ds3b)/pDenU);

      }


      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AuZu::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over colors 1/3
    return amp2/12.0;
  }
  //set_momenta(...) must be called before amp2_Aplus().
  ldouble AuZu::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/3
    return amp2/6.0;
  }  
  //set_momenta(...) must be called before amp2_Aminus().
  ldouble AuZu::amp2_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(-2,j2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/3
    return amp2/6.0;
  }



  //  Tests
  int test_AuZu(){
    int n=0;//Number of fails
    std::cout<<"\t* A , u  -> Z , u       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, MW=80.385, pspatial=250\n";
      ldouble mu=0.0042;
      ldouble EE=0.31333, MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      AuZu AuZuAmp = AuZu(EE,mu,MW,SW);
      ldouble pspatial=250;
      ldouble dataCHp[20] = {7.188794020568123E-03,6.956750492693353E-03,6.731509200135351E-03,6.514306912952145E-03,6.306699567668470E-03,6.110672324678114E-03,5.928798540047365E-03,5.764475043784086E-03,5.622280171986434E-03,5.508536388790908E-03,5.432228248498373E-03,5.406568330965885E-03,5.451815925313021E-03,5.600697575513200E-03,5.909740948454557E-03,6.485720526101144E-03,7.557438770005292E-03,9.719782081046223E-03,1.515400060187050E-02,4.349134581725004E-02};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCHp);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCHp);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCHp);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCHp);
      ldouble dataCHm[20] = {7.095635286767396E-03,7.365658655596717E-03,7.672636973606860E-03,8.023289322465960E-03,8.426068740401369E-03,8.891760138254413E-03,9.434343956003996E-03,1.007227437442124E-02,1.083042441976139E-02,1.174314255685062E-02,1.285923976201682E-02,1.425049688279350E-02,1.602697788538488E-02,1.836647840747699E-02,2.157610000500913E-02,2.623592339458992E-02,3.358897868272197E-02,4.686714414414430E-02,7.792064010391045E-02,2.334014489708591E-01};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCHm);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCHm);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCHm);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCHm);
      ldouble dataCH[20] = {7.142214653667760E-03,7.161204574145034E-03,7.202073086871106E-03,7.268798117709053E-03,7.366384154034919E-03,7.501216231466263E-03,7.681571248025681E-03,7.918374709102663E-03,8.226352295873911E-03,8.625839472820767E-03,9.145734005257597E-03,9.828532606879693E-03,1.073939690534895E-02,1.198358799149510E-02,1.374292047673184E-02,1.636082196034553E-02,2.057320872636363E-02,2.829346311259526E-02,4.653732035289048E-02,1.384463973940546E-01};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2p[20] = {7.812312120973011E-03,7.610718125981472E-03,7.415966212833532E-03,7.229300396412359E-03,7.052285727700692E-03,6.886918995884398E-03,6.735788631391063E-03,6.602311362387306E-03,6.491092344231581E-03,6.408491077186372E-03,6.363544745844597E-03,6.369543327565465E-03,6.446864787664596E-03,6.628427377869981E-03,6.971089896765957E-03,7.582252297608401E-03,8.692057338543046E-03,1.089887619565803E-02,1.639673655006546E-02,4.493188870180785E-02};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH2p);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH2p);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH2p);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH2p);
      ldouble dataCH2m[20] = {6.612894366557194E-03,6.856228386137994E-03,7.132837883763960E-03,7.448772946347565E-03,7.811644973551552E-03,8.231165061495451E-03,8.719921668109174E-03,9.294531557855952E-03,9.977391239353523E-03,1.079942922449185E-02,1.180459655646493E-02,1.305752712229067E-02,1.465732621533596E-02,1.676408700168514E-02,1.965433404065941E-02,2.385039149914546E-02,3.047152544679507E-02,4.242782765526033E-02,7.038963455240241E-02,2.103899717384897E-01};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH2m);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH2m);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH2m);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH2m);
      ldouble dataCH2[20] = {7.212603243765102E-03,7.233473256059733E-03,7.274402048298746E-03,7.339036671379961E-03,7.431965350626122E-03,7.559042028689925E-03,7.727855149750118E-03,7.948421460121629E-03,8.234241791792551E-03,8.603960150839111E-03,9.084070651154762E-03,9.713535224928068E-03,1.055209550150028E-02,1.169625718977756E-02,1.331271196871268E-02,1.571632189837693E-02,1.958179139266906E-02,2.666335192545918E-02,4.339318555123393E-02,1.276609302201488E-01};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH2);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH2);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH2);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, MW=80.385, pspatial=95\n";
      mu = 125;
      pspatial = 250;
      AuZuAmp.set_masses(mu,MW);
      ldouble dataCH4p[20] = {7.255170115361143E-03,7.201363293763426E-03,7.229978950419989E-03,7.355205441183829E-03,7.594684518186338E-03,7.970629597637306E-03,8.511408433353365E-03,9.253831807858904E-03,1.024654495487835E-02,1.155519523473241E-02,1.327056425755426E-02,1.552185605414214E-02,1.849940090725320E-02,2.249558799415393E-02,2.798370354868042E-02,3.578298364642562E-02,4.744393838084317E-02,6.629277912016067E-02,1.009590090412518E-01,1.825836726442554E-01};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH4p);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH4p);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH4p);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH4p);
      ldouble dataCH4m[20] = {7.163761274526369E-03,7.571692244220536E-03,8.080626922903213E-03,8.707199176374310E-03,9.471904676242978E-03,1.040028963921062E-02,1.152460697877688E-02,1.288616902603107E-02,1.453876381931825E-02,1.655374046057461E-02,1.902779682408420E-02,2.209530241549160E-02,2.594855584203938E-02,3.087261989370141E-02,3.730854051478683E-02,4.597583291707331E-02,5.812968917558978E-02,7.615743679502490E-02,1.051293111698474E-01,1.569721643802168E-01};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH4m);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH4m);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH4m);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH4m);
      ldouble dataCH4[20] = {7.209465694943756E-03,7.386527768991981E-03,7.655302936661601E-03,8.031202308779069E-03,8.533294597214658E-03,9.185459618423961E-03,1.001800770606512E-02,1.107000041694499E-02,1.239265438709830E-02,1.405446784765351E-02,1.614918054081923E-02,1.880857923481687E-02,2.222397837464629E-02,2.668410394392767E-02,3.264612203173362E-02,4.087940828174946E-02,5.278681377821647E-02,7.122510795759279E-02,1.030441601055496E-01,1.697779185122361E-01};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH4);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH4);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH4);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, MW=80.385, pspatial=125.1\n";
      mu = 125;
      pspatial = 125.1;
      AuZuAmp.set_masses(mu,MW);
      ldouble dataCH3p[20] = {8.024871352569568E-03,8.227834832297677E-03,8.505329069866079E-03,8.867640883958090E-03,9.327040938612453E-03,9.898286068404849E-03,1.059928232904104E-02,1.145197213007191E-02,1.248353904179510E-02,1.372807135972146E-02,1.522890193595946E-02,1.704196814760806E-02,1.924075123610089E-02,2.192373414242533E-02,2.522601388384114E-02,2.933804221014313E-02,3.453717510721791E-02,4.124354236125905E-02,5.012529063305870E-02,6.231276593130378E-02};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH3p);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH3p);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH3p);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2_Aplus(); }, 0,mu,MZ,mu,pspatial,dataCH3p);
      ldouble dataCH3m[20] = {6.985043707197275E-03,7.481774548899924E-03,8.047005789471354E-03,8.688705415348372E-03,9.416087074191797E-03,1.023985060381684E-02,1.117247513756522E-02,1.222857571406322E-02,1.342533476047939E-02,1.478301772398990E-02,1.632557384589000E-02,1.808130045593516E-02,2.008349380161393E-02,2.237087879453053E-02,2.498730271208365E-02,2.797944145848228E-02,3.138943957335385E-02,3.523463378083674E-02,3.945316976262692E-02,4.375351692850396E-02};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH3m);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH3m);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH3m);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2_Aminus(); }, 0,mu,MZ,mu,pspatial,dataCH3m);
      ldouble dataCH3[20] = {7.504957529883422E-03,7.854804690598800E-03,8.276167429668716E-03,8.778173149653231E-03,9.371564006402124E-03,1.006906833611084E-02,1.088587873330313E-02,1.184027392206756E-02,1.295443690113725E-02,1.425554454185568E-02,1.577723789092473E-02,1.756163430177161E-02,1.966212251885741E-02,2.214730646847793E-02,2.510665829796239E-02,2.865874183431270E-02,3.296330734028588E-02,3.823908807104789E-02,4.478923019784281E-02,5.303314142990388E-02};
      i += AuZuAmp.test_2to2_amp2([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH3);
      i += AuZuAmp.test_2to2_amp2_rotations([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH3);
      i += AuZuAmp.test_2to2_amp2_boosts([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH3);
      i += AuZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AuZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
