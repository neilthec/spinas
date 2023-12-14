
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

//File:  SPINAS/SM/AZee.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AZee.h"

namespace spinas {

  AZee::AZee(const ldouble& echarge, const ldouble& masse, const ldouble& massW, const ldouble& sinW):
    e(echarge), me(masse), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), prope(masse,0) {
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    MZ = MW/CW;//std::cout<<"constructor: MZ="<<MZ<<"\n";
    p1=particle(0);
    p2=particle(MZ);
    p3=particle(me);
    p4=particle(me);
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
    preTU = 2.0*e*e/(2.0*MW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
  }
  void AZee::set_masses(const ldouble& masse, const ldouble& massW){
    me=masse;
    MW=massW;
    MZ=MW/CW;//std::cout<<"set_masses: MZ="<<MZ<<"\n";
    p1.set_mass(0);
    p2.set_mass(MZ);
    p3.set_mass(me);
    p4.set_mass(me);
    prope.set_mass(me);
    //Couplings
    preTU = 2.0*e*e/(2.0*MW*SW);
  }
  void AZee::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
  cdouble AZee::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){//Double Spin
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
	//- preTU [1341] (gRe [24] <23> + gLe <24> [23])/((t-me^2) (u-me^2))
	//- preTU [12] ( gRe [14]<23>/(u-me^2) - gLe [13]<24>/(t-me^2) )
	//34 outgoing:
	//+ preTU [1341] (gRe [24] <23> + gLe <24> [23])/((t-me^2) (u-me^2))
	//+ preTU [12] ( gRe [14]<23>/(u-me^2) - gLe [13]<24>/(t-me^2) )
	
	amplitude += normFactor*preTU*s1341s.v()*(gR*s24s.v(ds2a,ds4)*a23a.v(ds2b,ds3) + gL*s23s.v(ds2a,ds3)*a24a.v(ds2b,ds4))/pDenT/pDenU;
	amplitude += normFactor*preTU*s12s.v(ds2a)*(gR*s14s.v(ds4)*a23a.v(ds2b,ds3)/pDenU - gL*s13s.v(ds3)*a24a.v(ds2b,ds4)/pDenT);

	
      }
      else if(ds1<0){
	//all ingoing:
	//preTU = 2.0*e*e/(2.0*MW*SW);
	//- preTU <1341> (gRe [24] <23> + gLe <24> [23])/((t-me^2) (u-me^2))
	//- preTU <12> ( gLe <14>[23]/(u-me^2) - gRe <13>[24]/(t-me^2) )
	//34 outgoing:
	//+ preTU <1341> (gRe [24] <23> + gLe <24> [23])/((t-me^2) (u-me^2))
	//+ preTU <12> ( gLe <14>[23]/(u-me^2) - gRe <13>[24]/(t-me^2) )
	
	amplitude += normFactor*preTU*a1341a.v()*(gR*s24s.v(ds2a,ds4)*a23a.v(ds2b,ds3) + gL*s23s.v(ds2a,ds3)*a24a.v(ds2b,ds4))/pDenT/pDenU;
	amplitude += normFactor*preTU*a12a.v(ds2a)*(gL*a14a.v(ds4)*s23s.v(ds2b,ds3)/pDenU - gR*a13a.v(ds3)*s24s.v(ds2b,ds4)/pDenT);

      }


      
    }
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble AZee::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/3=1/6
    return amp2/6.0;
  }
  //set_momenta(...) must be called before amp2_Aplus().
  ldouble AZee::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(2,j2,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }  
  //set_momenta(...) must be called before amp2_Aminus().
  ldouble AZee::amp2_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j2=-2;j2<=2;j2+=2)
      for(int j3=-1;j3<=1;j3+=2)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(-2,j2,j3,j4);
	  amp2 += std::pow(std::abs(M),2);
	}
    //Average over initial spins 1/3
    return amp2/3.0;
  }



  //  Tests
  int test_AZee(){
    int n=0;//Number of fails
    std::cout<<"\t* A , Z  -> E , e       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=250\n";
      ldouble me=0.0005;
      ldouble EE=0.31333, MW=80.385, SW=0.474;
      ldouble CW=std::sqrt(1-SW*SW);
      ldouble MZ=MW/CW;//std::cout<<"MZ="<<MZ<<"\n";
      AZee AZeeAmp = AZee(EE,me,MW,SW);
      ldouble pspatial=250;
      ldouble dataCHp[20] = {1.557732521509090E-01,5.012150932248519E-02,2.929342746741558E-02,2.062226688292527E-02,1.604182460314358E-02,1.336089670772296E-02,1.174676526044491E-02,1.082302595795997E-02,1.040584765529472E-02,1.040919725155873E-02,1.080798131501265E-02,1.162667763983240E-02,1.294455717553334E-02,1.491999571407773E-02,1.785003016182664E-02,2.231516960131901E-02,2.957678176222770E-02,4.293184243734666E-02,7.449044898311881E-02,2.332907389389991E-01};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCHp);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCHp);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCHp);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCHp);
      ldouble dataCHm[20] = {2.332907389389995E-01,7.449044898311885E-02,4.293184243734668E-02,2.957678176222772E-02,2.231516960131901E-02,1.785003016182662E-02,1.491999571407774E-02,1.294455717553335E-02,1.162667763983241E-02,1.080798131501265E-02,1.040919725155873E-02,1.040584765529471E-02,1.082302595795997E-02,1.174676526044490E-02,1.336089670772297E-02,1.604182460314358E-02,2.062226688292527E-02,2.929342746741557E-02,5.012150932248516E-02,1.557732521509087E-01};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCHm);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCHm);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCHm);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCHm);
      ldouble dataCH[20] = {1.945319955449543E-01,6.230597915280202E-02,3.611263495238113E-02,2.509952432257649E-02,1.917849710223129E-02,1.560546343477479E-02,1.333338048726133E-02,1.188379156674666E-02,1.101626264756356E-02,1.060858928328569E-02,1.060858928328569E-02,1.101626264756356E-02,1.188379156674666E-02,1.333338048726132E-02,1.560546343477480E-02,1.917849710223130E-02,2.509952432257648E-02,3.611263495238112E-02,6.230597915280199E-02,1.945319955449539E-01};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH);
      //std::cout<<"\n# me=0.0005, MW=80.385, pspatial=125.1\n";
      pspatial = 125.1;
      ldouble dataCH2p[20] = {1.876054316499420E-01,6.153923081208371E-02,3.668359415244473E-02,2.633305153385477E-02,2.086298352368973E-02,1.765875069923413E-02,1.572661570986701E-02,1.461726487049104E-02,1.411109491015079E-02,1.410509640217563E-02,1.456882303316146E-02,1.553073884183784E-02,1.708429054733000E-02,1.941661138034301E-02,2.287894609689913E-02,2.815795192284883E-02,3.674582224753741E-02,5.254304493143119E-02,8.987668790833418E-02,2.777467544278416E-01};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH2p);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH2p);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH2p);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH2p);
      ldouble dataCH2m[20] = {2.777467544278464E-01,8.987668790833454E-02,5.254304493143130E-02,3.674582224753749E-02,2.815795192284885E-02,2.287894609689915E-02,1.941661138034302E-02,1.708429054733001E-02,1.553073884183785E-02,1.456882303316147E-02,1.410509640217563E-02,1.411109491015079E-02,1.461726487049103E-02,1.572661570986700E-02,1.765875069923412E-02,2.086298352368971E-02,2.633305153385471E-02,3.668359415244465E-02,6.153923081208346E-02,1.876054316499388E-01};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH2m);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH2m);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH2m);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH2m);
      ldouble dataCH2[20] = {2.326760930388942E-01,7.570795936020913E-02,4.461331954193802E-02,3.153943689069613E-02,2.451046772326929E-02,2.026884839806664E-02,1.757161354510501E-02,1.585077770891052E-02,1.482091687599432E-02,1.433695971766855E-02,1.433695971766855E-02,1.482091687599431E-02,1.585077770891052E-02,1.757161354510501E-02,2.026884839806662E-02,2.451046772326927E-02,3.153943689069606E-02,4.461331954193792E-02,7.570795936020883E-02,2.326760930388903E-01};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH2);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH2);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH2);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH2);
      //std::cout<<"\n# me=125.1, MW=80.385, pspatial=95\n";
      me = 125;
      pspatial = 250;
      AZeeAmp.set_masses(me,MW);
      ldouble dataCH4p[20] = {2.734282173851775E-01,1.787652063682291E-01,1.357318576224262E-01,1.116762275620454E-01,9.673552613714616E-02,8.693851258371067E-02,8.039902311103191E-02,7.612780880918668E-02,7.358106469784126E-02,7.246656622384685E-02,7.265630189625377E-02,7.415159490338541E-02,7.708206417773766E-02,8.173758987285859E-02,8.864939185099913E-02,9.876772760765483E-02,1.138665646753578E-01,1.375748655861465E-01,1.785291691193781E-01,2.633501837423999E-01};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH4p);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH4p);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH4p);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH4p);
      ldouble dataCH4m[20] = {2.633501837424000E-01,1.785291691193782E-01,1.375748655861466E-01,1.138665646753578E-01,9.876772760765484E-02,8.864939185099914E-02,8.173758987285859E-02,7.708206417773769E-02,7.415159490338542E-02,7.265630189625377E-02,7.246656622384685E-02,7.358106469784126E-02,7.612780880918665E-02,8.039902311103191E-02,8.693851258371066E-02,9.673552613714616E-02,1.116762275620454E-01,1.357318576224261E-01,1.787652063682291E-01,2.734282173851774E-01};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH4m);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH4m);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH4m);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH4m);
      ldouble dataCH4[20] = {2.683892005637887E-01,1.786471877438036E-01,1.366533616042864E-01,1.127713961187016E-01,9.775162687240049E-02,8.779395221735491E-02,8.106830649194526E-02,7.660493649346219E-02,7.386632980061335E-02,7.256143406005031E-02,7.256143406005032E-02,7.386632980061333E-02,7.660493649346216E-02,8.106830649194524E-02,8.779395221735489E-02,9.775162687240049E-02,1.127713961187016E-01,1.366533616042863E-01,1.786471877438036E-01,2.683892005637886E-01};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH4);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH4);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH4);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH4);
      //std::cout<<"\n# me=125, MW=80.385, pspatial=125.1\n";
      me = 125;
      pspatial = 125.1;
      AZeeAmp.set_masses(me,MW);
      ldouble dataCH3p[20] = {9.903169412376710E-02,9.251216211368650E-02,8.729335036820035E-02,8.311144605556794E-02,7.977734236248796E-02,7.715392332736466E-02,7.514145219583854E-02,7.366799874428896E-02,7.268312263139391E-02,7.215375276391060E-02,7.206162950002686E-02,7.240194329421103E-02,7.318298481943469E-02,7.442676256132111E-02,7.617067605858130E-02,7.847048546938944E-02,8.140502700459946E-02,8.508344476474393E-02,8.965624081790566E-02,9.533237318116787E-02};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH3p);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH3p);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH3p);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2_Aplus(); }, 0,MZ,me,me,pspatial,dataCH3p);
      ldouble dataCH3m[20] = {9.533237318116787E-02,8.965624081790566E-02,8.508344476474393E-02,8.140502700459946E-02,7.847048546938944E-02,7.617067605858130E-02,7.442676256132111E-02,7.318298481943469E-02,7.240194329421103E-02,7.206162950002684E-02,7.215375276391060E-02,7.268312263139391E-02,7.366799874428896E-02,7.514145219583854E-02,7.715392332736466E-02,7.977734236248796E-02,8.311144605556794E-02,8.729335036820035E-02,9.251216211368650E-02,9.903169412376708E-02};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH3m);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH3m);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH3m);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2_Aminus(); }, 0,MZ,me,me,pspatial,dataCH3m);
      ldouble dataCH3[20] = {9.718203365246748E-02,9.108420146579609E-02,8.618839756647213E-02,8.225823653008370E-02,7.912391391593870E-02,7.666229969297297E-02,7.478410737857982E-02,7.342549178186182E-02,7.254253296280247E-02,7.210769113196873E-02,7.210769113196873E-02,7.254253296280247E-02,7.342549178186182E-02,7.478410737857982E-02,7.666229969297297E-02,7.912391391593870E-02,8.225823653008370E-02,8.618839756647215E-02,9.108420146579607E-02,9.718203365246748E-02};
      i += AZeeAmp.test_2to2_amp2([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH3);
      i += AZeeAmp.test_2to2_amp2_rotations([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH3);
      i += AZeeAmp.test_2to2_amp2_boosts([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH3);
      i += AZeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AZeeAmp.amp2(); }, 0,MZ,me,me,pspatial,dataCH3);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
