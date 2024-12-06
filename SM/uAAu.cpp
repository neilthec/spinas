
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

//File:  SPINAS/SM/uAAu.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/uAAu.h"

namespace spinas {

  uAAu::uAAu(const ldouble& echarge, const ldouble& massu):
    e(echarge), Qu(2.0/3.0), mu(massu), prop(massu,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(mu);
    p2=particle(0);
    p3=particle(0);
    p4=particle(mu);
    s34s = sproduct(SQUARE,&p3,&p4,2);
    a34a = sproduct(ANGLE,&p3,&p4,2);
    s12s = sproduct(SQUARE,&p1,&p2,2);
    a12a = sproduct(ANGLE,&p1,&p2,2);
    s13s = sproduct(SQUARE,&p1,&p3,2);
    a13a = sproduct(ANGLE,&p1,&p3,2);
    s24s = sproduct(SQUARE,&p2,&p4,2);
    a24a = sproduct(ANGLE,&p2,&p4,2);
    s23s = sproduct(SQUARE,&p2,&p3,2);
    a23a = sproduct(ANGLE,&p2,&p3,2);
    s14s = sproduct(SQUARE,&p1,&p4,2);
    a14a = sproduct(ANGLE,&p1,&p4,2);
    s342a = sproduct(SQUARE,&p3,&p4,&p2,2);
    s243a = sproduct(SQUARE,&p2,&p4,&p3,2);
  }
  void uAAu::set_masses(const ldouble& massu){
    mu=massu;
    p1.set_mass(mu);
    p4.set_mass(mu);
    prop.set_mass(mu);
  }
  void uAAu::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s34s.update();
    a34a.update();
    s12s.update();
    a12a.update();
    s13s.update();
    a13a.update();
    s24s.update();
    a24a.update();
    s23s.update();
    a23a.update();
    s14s.update();
    a14a.update();
    s342a.update();
    s243a.update();
    //Propagator Momentum
    ldouble propTP[4], propSP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propSP[j] = mom1[j]+mom2[j];
    }
    pDenT = prop.denominator(propTP);
    pDenS = prop.denominator(propSP);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble uAAu::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds2>0&&ds3>0){
      //AAUu all in:        mu[12]^2<34>/ut
      //uAAU 4->1->3->4:  - mu[23]^2<14>/ts
      //34 out:             mu[23]^2<14>/ts
      return 2.0*e*e*Qu*Qu*mu*s23s.v()*s23s.v()*a14a.v(ds1,ds4)/pDenT/pDenS;
    }
    else if(ds2<0&&ds3<0){
      //AAUu:              mu<12>^2[34]/ut
      //uAAU 4->1->3->4: - mu<23>^2[14]/ts
      //34 out:          - mu<23>^2[14]/ts
      return - 2.0*e*e*Qu*Qu*mu*a23a.v()*a23a.v()*s14s.v(ds1,ds4)/pDenT/pDenS;
    }
    else if(ds3>0&&ds2<0){
      //AAUu:    ([31]<42>+[41]<32>)[132>/ut
      //uAAU 4->1->3->4: - ([34]<12>+[13]<24>)[342>/ts
      //34 out:            ([34]<12>-[13]<24>)[342>/ts
      return 2.0*e*e*Qu*Qu*(s34s.v(ds4)*a12a.v(ds1)-s13s.v(ds1)*a24a.v(ds4))*s342a.v()/pDenT/pDenS;
    }
    else if(ds3<0&&ds2>0){
      //AAUu:   (<31>[42]+<41>[32])[231>
      //uAAU 4->1->3->4: - (<34>[12]+<13>[24])[243>
      //34 out:          - (<34>[12]-<13>[24])[243>
      return -2.0*e*e*Qu*Qu*(a34a.v(ds4)*s12s.v(ds1)-a13a.v(ds1)*s24s.v(ds4))*s243a.v()/pDenT/pDenS;
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble uAAu::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-2;j3<=2;j3+=4)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial color 1/3
    return amp2/12.0;
  }

  //u, A+ -> A, u
  ldouble uAAu::amp2_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spin 1/2
    //Average over initial color 1/3
    return amp2/6.0;
  }

  //u, A- -> A, u
  ldouble uAAu::amp2_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j3=-2;j3<=2;j3+=4)
	for(int j4=-1;j4<=1;j4+=2){
	  M = amp(j1,-2,j3,j4);
	  amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	}
    //Average over initial spin 1/2
    //Average over initial color 1/3
    return amp2/6.0;
  }



  //  Tests
  int test_uAAu(){
    int n=0;//Number of fails
    std::cout<<"\t* u , A  -> A , u       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n#  mu=0.0042, pspatial=250\n";
      ldouble mu=0.0042;
      ldouble EE=0.31333;
      uAAu uAAuAmp = uAAu(EE,mu);
      ldouble pspatial=250;
      ldouble dataCH[20] = {1.524067540229081E-01,5.105610403697496E-02,3.093828555368211E-02,2.242515731312624E-02,1.778025915590884E-02,1.489364740075031E-02,1.295380524687658E-02,1.158202485384500E-02,1.057781383550926E-02,9.825097652508974E-03,9.252020637615101E-03,8.811720399153455E-03,8.472330510783903E-03,8.211426450266096E-03,8.012769770459878E-03,7.864312591228618E-03,7.756927730650376E-03,7.683574219813775E-03,7.638733461463834E-03,7.618018878295231E-03};
      i += uAAuAmp.test_2to2_amp2([&]() { return uAAuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH);
      i += uAAuAmp.test_2to2_amp2_rotations([&]() { return uAAuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH);
      i += uAAuAmp.test_2to2_amp2_boosts([&]() { return uAAuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH);
      i += uAAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAAuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH);
      //std::cout<<"\n#  u , A+ -> A , u\n";
      ldouble dataCHpp[20] = {1.524067540229081E-01,5.105610403697496E-02,3.093828555368211E-02,2.242515731312624E-02,1.778025915590884E-02,1.489364740075031E-02,1.295380524687658E-02,1.158202485384500E-02,1.057781383550926E-02,9.825097652508974E-03,9.252020637615101E-03,8.811720399153455E-03,8.472330510783903E-03,8.211426450266096E-03,8.012769770459878E-03,7.864312591228618E-03,7.756927730650376E-03,7.683574219813775E-03,7.638733461463834E-03,7.618018878295231E-03};
      i += uAAuAmp.test_2to2_amp2([&]() { return uAAuAmp.amp2_Aplus(); }, mu,0,0,mu,pspatial,dataCHpp);
      i += uAAuAmp.test_2to2_amp2_rotations([&]() { return uAAuAmp.amp2_Aplus(); }, mu,0,0,mu,pspatial,dataCHpp);
      i += uAAuAmp.test_2to2_amp2_boosts([&]() { return uAAuAmp.amp2_Aplus(); }, mu,0,0,mu,pspatial,dataCHpp);
      i += uAAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAAuAmp.amp2_Aplus(); }, mu,0,0,mu,pspatial,dataCHpp);
      //std::cout<<"\n#  u , A- -> A , u\n";
      ldouble dataCHpm[20] = {1.524067540229081E-01,5.105610403697496E-02,3.093828555368211E-02,2.242515731312624E-02,1.778025915590884E-02,1.489364740075031E-02,1.295380524687658E-02,1.158202485384500E-02,1.057781383550926E-02,9.825097652508974E-03,9.252020637615101E-03,8.811720399153455E-03,8.472330510783903E-03,8.211426450266096E-03,8.012769770459878E-03,7.864312591228618E-03,7.756927730650376E-03,7.683574219813775E-03,7.638733461463834E-03,7.618018878295231E-03};
      i += uAAuAmp.test_2to2_amp2([&]() { return uAAuAmp.amp2_Aminus(); }, mu,0,0,mu,pspatial,dataCHpm);
      i += uAAuAmp.test_2to2_amp2_rotations([&]() { return uAAuAmp.amp2_Aminus(); }, mu,0,0,mu,pspatial,dataCHpm);
      i += uAAuAmp.test_2to2_amp2_boosts([&]() { return uAAuAmp.amp2_Aminus(); }, mu,0,0,mu,pspatial,dataCHpm);
      i += uAAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAAuAmp.amp2_Aminus(); }, mu,0,0,mu,pspatial,dataCHpm);
      //Close to threshold
      //std::cout<<"\n#  mu=0.0042, pspatial=0.006\n";
      pspatial = 0.006;
      ldouble dataCH2[20] = {2.922261414639316E-02,1.967993469579533E-02,1.508799181348436E-02,1.248736690279479E-02,1.086478720704472E-02,9.788243737713162E-03,9.045054711568853E-03,8.519252071632857E-03,8.142576093589773E-03,7.872386772156488E-03,7.680750106496581E-03,7.548639630628491E-03,7.462668590241912E-03,7.413156077591419E-03,7.392934468742411E-03,7.396587654583180E-03,7.419949410722770E-03,7.459764172083856E-03,7.513452197439041E-03,7.578943584940903E-03};
      i += uAAuAmp.test_2to2_amp2([&]() { return uAAuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH2);
      i += uAAuAmp.test_2to2_amp2_rotations([&]() { return uAAuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH2);
      i += uAAuAmp.test_2to2_amp2_boosts([&]() { return uAAuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH2);
      i += uAAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAAuAmp.amp2(); }, mu,0,0,mu,pspatial,dataCH2);
      //std::cout<<"\n#  u , A+ -> A , u\n";
      ldouble dataCH2pp[20] = {2.922261414639316E-02,1.967993469579533E-02,1.508799181348436E-02,1.248736690279479E-02,1.086478720704472E-02,9.788243737713162E-03,9.045054711568853E-03,8.519252071632857E-03,8.142576093589773E-03,7.872386772156488E-03,7.680750106496581E-03,7.548639630628491E-03,7.462668590241912E-03,7.413156077591419E-03,7.392934468742411E-03,7.396587654583180E-03,7.419949410722770E-03,7.459764172083856E-03,7.513452197439041E-03,7.578943584940903E-03};
      i += uAAuAmp.test_2to2_amp2([&]() { return uAAuAmp.amp2_Aplus(); }, mu,0,0,mu,pspatial,dataCH2pp);
      i += uAAuAmp.test_2to2_amp2_rotations([&]() { return uAAuAmp.amp2_Aplus(); }, mu,0,0,mu,pspatial,dataCH2pp);
      i += uAAuAmp.test_2to2_amp2_boosts([&]() { return uAAuAmp.amp2_Aplus(); }, mu,0,0,mu,pspatial,dataCH2pp);
      i += uAAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAAuAmp.amp2_Aplus(); }, mu,0,0,mu,pspatial,dataCH2pp);
      //std::cout<<"\n#  u , A- -> A , u\n";
      ldouble dataCH2pm[20] = {2.922261414639316E-02,1.967993469579533E-02,1.508799181348436E-02,1.248736690279479E-02,1.086478720704472E-02,9.788243737713162E-03,9.045054711568853E-03,8.519252071632857E-03,8.142576093589773E-03,7.872386772156488E-03,7.680750106496581E-03,7.548639630628491E-03,7.462668590241912E-03,7.413156077591419E-03,7.392934468742411E-03,7.396587654583180E-03,7.419949410722770E-03,7.459764172083856E-03,7.513452197439041E-03,7.578943584940903E-03};
      i += uAAuAmp.test_2to2_amp2([&]() { return uAAuAmp.amp2_Aminus(); }, mu,0,0,mu,pspatial,dataCH2pm);
      i += uAAuAmp.test_2to2_amp2_rotations([&]() { return uAAuAmp.amp2_Aminus(); }, mu,0,0,mu,pspatial,dataCH2pm);
      i += uAAuAmp.test_2to2_amp2_boosts([&]() { return uAAuAmp.amp2_Aminus(); }, mu,0,0,mu,pspatial,dataCH2pm);
      i += uAAuAmp.test_2to2_amp2_boosts_and_rotations([&]() { return uAAuAmp.amp2_Aminus(); }, mu,0,0,mu,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
