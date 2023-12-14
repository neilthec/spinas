
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

//File:  SPINAS/SM/guZu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/guZu.h"

namespace spinas {

  guZu::guZu(const ldouble& echarge, const ldouble& gscharge, const ldouble& massu, const ldouble& massW, const ldouble& sinW):
    e(echarge), Qu(2.0/3.0), gs(gscharge), mu(massu), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), prope(massu,0) {
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
    //[1421], <1421>
    s1421s = sproduct(SQUARE,&p1,&p4,&p2,&p1);
    a1421a = sproduct(ANGLE,&p1,&p4,&p2,&p1);
    //Couplings
    preTU = 2.0*e*gs/(2.0*MW*SW);
    gL=-2.0*Qu*SW*SW+1.0;
    gR=-2.0*Qu*SW*SW;
  }
  void guZu::set_masses(const ldouble& massu, const ldouble& massW){
    mu=massu;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(0);
    p2.set_mass(mu);
    p3.set_mass(MZ);
    p4.set_mass(mu);
    prope.set_mass(mu);
    //Couplings
    preTU = 2.0*e*gs/(2.0*MW*SW);
  }
  void guZu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
  cdouble guZu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    cdouble amplitude(0,0);
    int ds3a, ds3b;
    constexpr ldouble two=2;

    //Symmetrize the Z-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3);
    ldouble normFactor=get_spin_normalization(ds3);
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, i);
      
      if(ds1>0){
	//preTU = 2.0*e*gs/(2.0*MW*SW);
	//gZUu all in:
	//- preTU [1341] (gRe [24] <23> + gLe <24> [23])/((t-mu^2) (u-mu^2))
	//- preTU [12] ( gRe [14]<23>/(u-mu^2) - gLe [13]<24>/(t-mu^2) )
	//guZU: 4->2->3->4
	//+ preTU [1421] (gRe [23] <34> + gLe <23> [34])/((u-mu^2) (s-mu^2))
	//- preTU [13] ( gRe [12]<34>/(s-mu^2) + gLe [14]<23>/(u-mu^2) )
	//34 out:
	//- preTU [1421] (gRe [23] <34> - gLe <23> [34])/((u-mu^2) (s-mu^2))
	//- preTU [13] ( gRe [12]<34>/(s-mu^2) - gLe [14]<23>/(u-mu^2) )
	amplitude += - normFactor*preTU*s1421s.v()*(gR*s23s.v(ds2,ds3a)*a34a.v(ds3b,ds4) - gL*s34s.v(ds3a,ds4)*a23a.v(ds2,ds3b))/pDenU/pDenS;
	amplitude += - normFactor*preTU*s13s.v(ds3a)*(gR*s12s.v(ds2)*a34a.v(ds3b,ds4)/pDenS - gL*s14s.v(ds4)*a23a.v(ds2,ds3b)/pDenU);

	
      }
      else if(ds1<0){
	//preTU = 2.0*e*gs/(2.0*MW*SW);
	//gZUu all in:
	//- preTU <1341> (gRe [24] <23> + gLe <24> [23])/((t-mu^2) (u-mu^2))
	//- preTU <12> ( gLe <14>[23]/(u-mu^2) - gRe <13>[24]/(t-mu^2) )
	//guZU: 4->2->3->4
	//+ preTU <1421> (gRe [23] <34> + gLe <23> [34])/((u-mu^2) (s-mu^2))
	//- preTU <13> ( gLe <12>[34]/(s-mu^2) + gRe <14>[23]/(u-mu^2) )
	//34 out:
	//- preTU <1421> (gRe [23] <34> - gLe <23> [34])/((u-mu^2) (s-mu^2))
	//+ preTU <13> ( gLe <12>[34]/(s-mu^2) - gRe <14>[23]/(u-mu^2) )
	amplitude += - normFactor*preTU*a1421a.v()*(gR*s23s.v(ds2,ds3a)*a34a.v(ds3b,ds4) - gL*s34s.v(ds3a,ds4)*a23a.v(ds2,ds3b))/pDenU/pDenS;
	amplitude += + normFactor*preTU*a13a.v(ds3a)*(gL*a12a.v(ds2)*s34s.v(ds3b,ds4)/pDenS - gR*a14a.v(ds4)*s23s.v(ds2,ds3b)/pDenU);

      }


      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble guZu::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4 // It's the same for both diagrams
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/8*1/3=1/24
    return amp2/96.0;
  }
  //set_momenta(...) must be called before amp2_gplus().
  ldouble guZu::amp2_gplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/8*1/3=1/24
    return amp2/48.0;
  }  
  //set_momenta(...) must be called before amp2_gminus().
  ldouble guZu::amp2_gminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-1;j2<=1;j2+=2)
      for(int j3=-2;j3<=2;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(-2,j2,j3,j4);
	  amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
	}
    //Average over initial spins 1/2
    //Average over initial colors 1/8*1/3
    return amp2/48.0;
  }



  //  Tests
  int test_guZu(){
    int n=0;//Number of fails
    std::cout<<"\t* g , u  -> Z , u       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, MW=80.385, pspatial=250\n";
      ldouble mu=0.0042;
      ldouble EE=0.31333, gs=1.238, MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      guZu guZuAmp = guZu(EE,gs,mu,MW,SW);
      ldouble pspatial=250;
      ldouble dataCHp[20] = {4.208473915350701E-02,4.072630666609273E-02,3.940769592042651E-02,3.813614723319232E-02,3.692076938989657E-02,3.577318394446211E-02,3.470845587419792E-02,3.374647094915400E-02,3.291403172896647E-02,3.224815127218597E-02,3.180142708307828E-02,3.165120843263906E-02,3.191609753642680E-02,3.278768258886972E-02,3.459688865321723E-02,3.796879640487209E-02,4.424286443515851E-02,5.690168508645553E-02,8.871476364980165E-02,2.546076489209643E-01};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCHp);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCHp);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCHp);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCHp);
      ldouble dataCHm[20] = {4.153936798267400E-02,4.312014258965531E-02,4.491725938035426E-02,4.697005329725849E-02,4.932801023575392E-02,5.205426737270950E-02,5.523066919665557E-02,5.896525043244615E-02,6.340362310052799E-02,6.874687046721771E-02,7.528074243725721E-02,8.342545945875263E-02,9.382535955205232E-02,1.075212964421677E-01,1.263111083810175E-01,1.535907120657006E-01,1.966370718441239E-01,2.743702950081536E-01,4.561641082029795E-01,1.366382048210766E+00};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCHm);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCHm);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCHm);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCHm);
      ldouble dataCH[20] = {4.181205356809051E-02,4.192322462787402E-02,4.216247765039038E-02,4.255310026522541E-02,4.312438981282524E-02,4.391372565858580E-02,4.496956253542675E-02,4.635586069080008E-02,4.815882741474723E-02,5.049751086970184E-02,5.354108476016774E-02,5.753833394569585E-02,6.287072854423956E-02,7.015448951551873E-02,8.045399851711735E-02,9.577975423528637E-02,1.204399681396412E-01,1.656359900473046E-01,2.724394359263905E-01,8.104948485658653E-01};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2p[20] = {4.573494759430417E-02,4.455477319093994E-02,4.341465379942206E-02,4.232187214919185E-02,4.128559038374404E-02,4.031749813443980E-02,3.943274862712622E-02,3.865134411993439E-02,3.800024417819160E-02,3.751667867785861E-02,3.725355323219781E-02,3.728867021375153E-02,3.774132659393573E-02,3.880423286541419E-02,4.081025260727365E-02,4.438812813775834E-02,5.088516443137724E-02,6.380435445059343E-02,9.599000593200263E-02,2.630408953545942E-01};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH2p);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH2p);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH2p);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH2p);
      ldouble dataCH2m[20] = {3.871329929192512E-02,4.013782873482367E-02,4.175715994972075E-02,4.360671144058780E-02,4.573104197047788E-02,4.818700237499925E-02,5.104828818177581E-02,5.441217748724620E-02,5.840978424744114E-02,6.322217059205462E-02,6.910663533686176E-02,7.644155909310887E-02,8.580712546435219E-02,9.814055411715089E-02,1.150606790193825E-01,1.396252976617118E-01,1.783868332252586E-01,2.483815859256558E-01,4.120759894876539E-01,1.231667934259203E+00};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH2m);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH2m);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH2m);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH2m);
      ldouble dataCH2[20] = {4.222412344311464E-02,4.234630096288181E-02,4.258590687457141E-02,4.296429179488982E-02,4.350831617711096E-02,4.425225025471953E-02,4.524051840445101E-02,4.653176080359030E-02,4.820501421281637E-02,5.036942463495662E-02,5.318009428452979E-02,5.686511465343020E-02,6.177422602914397E-02,6.847239349128253E-02,7.793546581332805E-02,9.200671289973504E-02,1.146359988283179E-01,1.560929701881246E-01,2.540329977098283E-01,7.473544148068987E-01};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH2);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH2);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH2);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH2);
      //std::cout<<"\n# mu=125.1, MW=80.385, pspatial=95\n";
      mu = 125;
      pspatial = 250;
      guZuAmp.set_masses(mu,MW);
      ldouble dataCH4p[20] = {4.247331902203578E-02,4.215832236958748E-02,4.232584454961548E-02,4.305894723468710E-02,4.446090901847691E-02,4.666177199486360E-02,4.982760706757987E-02,5.417391243787379E-02,5.998546771780501E-02,6.764658660829580E-02,7.768872408934147E-02,9.086826821687556E-02,1.082994532115381E-01,1.316939879109116E-01,1.638225912476781E-01,2.094812465881472E-01,2.777469719483004E-01,3.880921207438537E-01,5.910356519523505E-01,1.068883906666107E+00};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH4p);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH4p);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH4p);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH4p);
      ldouble dataCH4m[20] = {4.193819209923769E-02,4.432630732455881E-02,4.730571988489277E-02,5.097380799156415E-02,5.545055769382946E-02,6.088552201312201E-02,6.746751640968311E-02,7.543839211373882E-02,8.511303581655540E-02,9.690914043508526E-02,1.113927960262002E-01,1.293506304413950E-01,1.519084008942205E-01,1.807349259826359E-01,2.184121830825203E-01,2.691523682754426E-01,3.403036446766632E-01,4.458419385710855E-01,6.154494933798710E-01,9.189486544477142E-01};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH4m);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH4m);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH4m);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH4m);
      ldouble dataCH4[20] = {4.220575556063673E-02,4.324231484707314E-02,4.481578221725412E-02,4.701637761312562E-02,4.995573335615319E-02,5.377364700399280E-02,5.864756173863149E-02,6.480615227580630E-02,7.254925176718020E-02,8.227786352169053E-02,9.454076005777085E-02,1.101094493291353E-01,1.301039270528793E-01,1.562144569467737E-01,1.911173871650992E-01,2.393168074317950E-01,3.090253083124818E-01,4.169670296574696E-01,6.032425726661107E-01,9.939162805569105E-01};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH4);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH4);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH4);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH4);
      //std::cout<<"\n# mu=125, MW=80.385, pspatial=125.1\n";
      mu = 125;
      pspatial = 125.1;
      guZuAmp.set_masses(mu,MW);
      ldouble dataCH3p[20] = {4.697931484016164E-02,4.816750650034925E-02,4.979201717227201E-02,5.191306809467572E-02,5.460249436171289E-02,5.794668564209184E-02,6.205046781919286E-02,6.704229645589706E-02,7.308130999211027E-02,8.036706860728074E-02,8.915325209422845E-02,9.976732983340929E-02,1.126394767428362E-01,1.283462330420046E-01,1.476784855911480E-01,1.717511797049480E-01,2.021880167003858E-01,2.414485262860587E-01,2.934441820492293E-01,3.647922715060921E-01};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH3p);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH3p);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH3p);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2_gplus(); }, 0,mu,MZ,mu,pspatial,dataCH3p);
      ldouble dataCH3m[20] = {4.089194120073165E-02,4.379990988682093E-02,4.710889457226360E-02,5.086554155542519E-02,5.512378950214430E-02,5.994627755358279E-02,6.540606122782133E-02,7.158870008972519E-02,7.859478374630550E-02,8.654295046370222E-02,9.557340422744504E-02,1.058518036637363E-01,1.175730721333425E-01,1.309639135589164E-01,1.462810192890536E-01,1.637976400592417E-01,1.837604990266460E-01,2.062710890857704E-01,2.309673018156620E-01,2.561424547311097E-01};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH3m);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH3m);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH3m);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2_gminus(); }, 0,mu,MZ,mu,pspatial,dataCH3m);
      ldouble dataCH3[20] = {4.393562802044664E-02,4.598370819358508E-02,4.845045587226781E-02,5.138930482505046E-02,5.486314193192859E-02,5.894648159783731E-02,6.372826452350709E-02,6.931549827281112E-02,7.583804686920788E-02,8.345500953549148E-02,9.236332816083674E-02,1.028095667485728E-01,1.151062744380894E-01,1.296550733004605E-01,1.469797524401008E-01,1.677744098820949E-01,1.929742578635159E-01,2.238598076859145E-01,2.622057419324457E-01,3.104673631186009E-01};
      i += guZuAmp.test_2to2_amp2([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH3);
      i += guZuAmp.test_2to2_amp2_rotations([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH3);
      i += guZuAmp.test_2to2_amp2_boosts([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH3);
      i += guZuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return guZuAmp.amp2(); }, 0,mu,MZ,mu,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
