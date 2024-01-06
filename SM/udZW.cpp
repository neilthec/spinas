
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

//File:  SPINAS/SM/udZW.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/udZW.h"

namespace spinas {

  udZW::udZW(const ldouble& echarge, const ldouble& massu, const ldouble& massd, const ldouble& massW, const ldouble& sinW, const ldouble& widthW):
    e(echarge), Qu(2.0/3.0), Qd(-1.0/3.0), mu(massu), md(massd), MW(massW), SW(sinW), WW(widthW) {
    constexpr ldouble sqrt2 = std::sqrt(2);
    CW=std::sqrt(1.0-SW*SW);
    MZ=MW/CW;
    propW = propagator(MW,WW);
    propu = propagator(mu,0);
    propd = propagator(md,0);
    p1=particle(mu);
    p2=particle(md);
    p3=particle(MZ);
    p4=particle(MW);
    //Spinor Products
    a12a= sproduct(ANGLE,&p1,&p2);
    a13a= sproduct(ANGLE,&p1,&p3);
    a23a= sproduct(ANGLE,&p2,&p3);
    a24a= sproduct(ANGLE,&p2,&p4);
    a34a= sproduct(ANGLE,&p3,&p4);
    s12s= sproduct(SQUARE,&p1,&p2);
    s13s= sproduct(SQUARE,&p1,&p3);
    s14s= sproduct(SQUARE,&p1,&p4);
    s23s= sproduct(SQUARE,&p2,&p3);
    s34s= sproduct(SQUARE,&p3,&p4);
    s132a= sproduct(SQUARE,&p1,&p3,&p2);
    s314a= sproduct(SQUARE,&p3,&p1,&p4);
    s413a= sproduct(SQUARE,&p4,&p1,&p3);
    //prefactor
    preud = e*e/(sqrt2*MW*MW*SW*SW);
    preW = preud/(MZ*MZ);
    gLu=-2.0*Qu*SW*SW+1.0;
    gRu=-2.0*Qu*SW*SW;
    gLd=-2.0*Qd*SW*SW-1.0;
    gRd=-2.0*Qd*SW*SW;
  }
  void udZW::set_masses(const ldouble& massu, const ldouble& massd, const ldouble& massW){
    constexpr ldouble sqrt2 = std::sqrt(2);
    mu=massu;
    md=massd;
    MW=massW;
    MZ=MW/CW;
    p1.set_mass(mu);
    p2.set_mass(md);
    p3.set_mass(MZ);
    p4.set_mass(MW);
    propW.set_mass(MW);
    propu.set_mass(mu);
    propd.set_mass(md);
    preud = e*e/(sqrt2*MW*MW*SW*SW);
    preW = preud/(MZ*MZ);
  }
  void udZW::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    constexpr ldouble one=1, two=2, three=3;
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    //Spinor Products
    a12a.update();
    a13a.update();
    a23a.update();
    a24a.update();
    a34a.update();
    s12s.update();
    s13s.update();
    s14s.update();
    s23s.update();
    s34s.update();
    s132a.update();
    s314a.update();
    s413a.update();
    //Propagator Momentum
    ldouble propSP[4], propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propSP[j] = mom1[j]+mom2[j];
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenS=propW.denominator(propSP);
    pDenT=propu.denominator(propTP);
    pDenU=propd.denominator(propUP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble udZW::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
    constexpr ldouble two=2;
    cdouble amplitude(0,0);
    int ds3a, ds3b, ds4a, ds4b;

    //Symmetrize the W-Boson Spin indices for the spin-0 cases
    int nCombs=get_num_spin_loops(ds3,ds4);
    ldouble normFactor=get_spin_normalization(ds3,ds4);
    //Start the loop
    for(int i=0;i<nCombs;i++){
      get_spinor_spins(ds3,ds3a,ds3b, ds4,ds4a,ds4b, i);
      
      
      //W Diagram
      //preW = e*e/(sqrt2*MW*MW*MZ*MZ*SW*SW);
      //UdZW all ingoing:
      // - preW (2(MZ<34>+MW[34])MW^2 <23>[14]+2(MW<34>+MZ[34])MW^2 <24>[13]+<34>[34]((2MW^2 -MZ^2)(Mu<12>-Md[12])+2MW^2 [132>)))/(s-MW^2)
      amplitude += + normFactor*preW*(two*(MW*a34a.v(ds3a,ds4b)+MZ*s34s.v(ds3a,ds4b))*MW*MW*a24a.v(ds2,ds4a)*s13s.v(ds1,ds3b)
				      +two*(MZ*a34a.v(ds3a,ds4b)+MW*s34s.v(ds3a,ds4b))*MW*MW*a23a.v(ds2,ds3b)*s14s.v(ds1,ds4a)
				      +a34a.v(ds3a,ds4a)*s34s.v(ds3b,ds4b)*(
									    (MZ*MZ-two*MW*MW)*mu*a12a.v(ds1,ds2)
									    +(two*MW*MW-MZ*MZ)*md*s12s.v(ds1,ds2)
									    +two*MW*MW*s132a.v(ds1,ds2)
									    )
				      )/pDenS;
      

      //u Diagram
      //preu = e*e/(sqrt2*MW*MW*SW*SW);
      //UdZW all in:
      // - preu<24>(gRu*mu<13>[34]+(MZ[34]+[413>)gLu[13]))/(t-mu^2)
      //34 out:
      // - preu<24>(gRu*mu<13>[34]-(MZ[34]-[413>)gLu[13]))/(t-mu^2)
      amplitude += - normFactor*preud*a24a.v(ds2,ds4a)*(
						       +gRu*mu*a13a.v(ds1,ds3b)*s34s.v(ds3a,ds4b)
						       -gLu*s13s.v(ds1,ds3a)*(MZ*s34s.v(ds3b,ds4b)-s413a.v(ds4b,ds3b))
						       )/pDenT;
      
      //d Diagram
      //pred = e*e/(sqrt2*MW*MW*SW*SW);
      //UdZW all in:
      // + pred[14](gRd*md<34>[23]+(MW[34]-[314>)gLd<23>))/(u-md^2)
      //34 out:
      // + pred[14](gRd*md<34>[23]-(MW[34]+[314>)gLd<23>))/(u-md^2)
      amplitude += + normFactor*preud*s14s.v(ds1,ds4a)*(
						       +gRd*md*a34a.v(ds3a,ds4b)*s23s.v(ds2,ds3b)
						       -gLd*a23a.v(ds2,ds3a)*(MW*s34s.v(ds3b,ds4b)+s314a.v(ds3b,ds4b))
						       )/pDenU;
      
      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble udZW::amp2(){
    constexpr ldouble three=3;
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-2;j3<=2;j3+=2)
	  for(int j4=-2;j4<=2;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    //Color factor 3
	    amp2 += three*std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2^2=1/4
    //Average over initial colors 1/3^2=1/9
    return amp2/36.0;
  }



  



  //  Tests
  int test_udZW(){
    int n=0;//Number of fails
    std::cout<<"\t* U , d  -> Z , W-      :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=250\n";
      ldouble EE=0.31333,mu=0.0042,md=0.0075,MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;
      udZW udZWAmp = udZW(EE,mu,md,MW,SW,0);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.853981861815233E-01,5.035120659362118E-02,2.359292682800378E-02,1.312639968155064E-02,8.057074895575613E-03,5.372274896139220E-03,3.916673057722036E-03,3.165002369040677E-03,2.859519260411860E-03,2.878844726634189E-03,3.186668293407384E-03,3.815737728599906E-03,4.874790501200399E-03,6.581716804302270E-03,9.344820465262320E-03,1.395969701164108E-02,2.214610310409112E-02,3.835944240376929E-02,7.867694993353912E-02,2.787903598858085E-01};
      i += udZWAmp.test_2to2_amp2([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH);
      i += udZWAmp.test_2to2_amp2_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH);
      i += udZWAmp.test_2to2_amp2_boosts([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH);
      i += udZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH);
      //std::cout<<"\n# mu=0.0042, md=0.0075, MW=80.385, pspatial=90\n";
      pspatial = 90;
      ldouble dataCH2[20] = {7.853396605194745E-03,6.679539409602917E-03,5.703713116222286E-03,4.943371836087850E-03,4.396239956068659E-03,4.053820522737164E-03,3.907240819317174E-03,3.949822039178119E-03,4.178231800802172E-03,4.593045325411996E-03,5.199096312374271E-03,6.005797411666049E-03,7.027511623392201E-03,8.283995059901130E-03,9.800868438984544E-03,1.160996613141707E-02,1.374917889826977E-02,1.626087040227040E-02,1.918665942548638E-02,2.255311582430134E-02};      
      i += udZWAmp.test_2to2_amp2([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH2);
      i += udZWAmp.test_2to2_amp2_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH2);
      i += udZWAmp.test_2to2_amp2_boosts([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH2);
      i += udZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH2);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=250\n";
      mu=80;
      udZWAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH3[20] = {1.808134234672233E-01,8.235658415491745E-02,4.825229128219211E-02,3.188024663413759E-02,2.281554109889982E-02,1.743669828494671E-02,1.417008297855480E-02,1.223881959658271E-02,1.123586099472753E-02,1.095431797862807E-02,1.131651361024739E-02,1.235036625277856E-02,1.419845019155493E-02,1.716545596783116E-02,2.183649724778160E-02,2.936633043993690E-02,4.227428865187535E-02,6.715498588422031E-02,1.282343285419850E-01,4.433175714517522E-01};
      i += udZWAmp.test_2to2_amp2([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH3);
      i += udZWAmp.test_2to2_amp2_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH3);
      i += udZWAmp.test_2to2_amp2_boosts([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH3);
      i += udZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH3);
      //std::cout<<"\n# mu=80, md=0.0075, MW=80.385, pspatial=70\n";
      pspatial = 70;
      ldouble dataCH4[20] = {2.718453708750528E-02,2.867643854932377E-02,3.036482216549734E-02,3.226851868319656E-02,3.441049365977095E-02,3.681871978035253E-02,3.952731559709746E-02,4.257804489945773E-02,4.602231181637014E-02,4.992384866631025E-02,5.436238924217993E-02,5.943877141594949E-02,6.528215795930405E-02,7.206047277919701E-02,7.999585196644807E-02,8.938816065090731E-02,1.006519507332056E-01,1.143767596695453E-01,1.314299650361124E-01,1.531418964279022E-01};
      i += udZWAmp.test_2to2_amp2([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH4);
      i += udZWAmp.test_2to2_amp2_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH4);
      i += udZWAmp.test_2to2_amp2_boosts([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH4);
      i += udZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH4);
      //std::cout<<"\n# mu=80, md=40, MW=30, pspatial=250\n";
      mu=80;
      md=40;
      MW=30;
      MZ=MW/CW;
      udZWAmp.set_masses(mu,md,MW);
      pspatial=250;
      ldouble dataCH5[20] = {3.387626811235170E+00,1.628149637487752E+00,1.054006735031317E+00,7.726772606228829E-01,6.081373230622453E-01,5.021421165846151E-01,4.299242456873747E-01,3.792646566227967E-01,3.435605606552616E-01,3.190863011739691E-01,3.037922556712566E-01,2.967740946037333E-01,2.981090061440526E-01,3.089936785127473E-01,3.323095663868519E-01,3.741056876004766E-01,4.476459431651361E-01,5.866531482780066E-01,9.060649219011626E-01,2.198352014545872E+00};
      i += udZWAmp.test_2to2_amp2([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH5);
      i += udZWAmp.test_2to2_amp2_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH5);
      i += udZWAmp.test_2to2_amp2_boosts([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH5);
      i += udZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH5);
      //std::cout<<"\n# mu=80, md=40, MW=1, pspatial=1\n";
      mu=80;
      md=40;
      MW=1;
      MZ=MW/CW;
      udZWAmp.set_masses(mu,md,MW);
      pspatial = 1;
      ldouble dataCH6[20] = {1.833943887570910E+05,1.833816697704371E+05,1.833703395366404E+05,1.833604027996034E+05,1.833518643777542E+05,1.833447291647934E+05,1.833390021304511E+05,1.833346883212588E+05,1.833317928613326E+05,1.833303209531701E+05,1.833302778784603E+05,1.833316689989070E+05,1.833344997570661E+05,1.833387756771966E+05,1.833445023661257E+05,1.833516855141285E+05,1.833603308958225E+05,1.833704443710762E+05,1.833820318859343E+05,1.833950994735564E+05};
      i += udZWAmp.test_2to2_amp2([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH6);
      i += udZWAmp.test_2to2_amp2_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH6);
      i += udZWAmp.test_2to2_amp2_boosts([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH6);
      i += udZWAmp.test_2to2_amp2_boosts_and_rotations([&]() { return udZWAmp.amp2(); }, mu,md,MZ,MW,pspatial,dataCH6);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }



}
