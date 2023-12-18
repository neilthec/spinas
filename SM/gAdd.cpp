
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

//File:  SPINAS/SM/gAdd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/gAdd.h"

namespace spinas {

  gAdd::gAdd(const ldouble& echarge, const ldouble& gscharge, const ldouble& massd):
    e(echarge), Qd(-1.0/3.0), gs(gscharge), md(massd), prop(massd,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(0);
    p2=particle(0);
    p3=particle(md);
    p4=particle(md);
    s34s = sproduct(SQUARE,&p3,&p4);
    a34a = sproduct(ANGLE,&p3,&p4);
    s12s = sproduct(SQUARE,&p1,&p2);
    a12a = sproduct(ANGLE,&p1,&p2);
    s13s = sproduct(SQUARE,&p1,&p3);
    a13a = sproduct(ANGLE,&p1,&p3);
    s24s = sproduct(SQUARE,&p2,&p4);
    a24a = sproduct(ANGLE,&p2,&p4);
    s23s = sproduct(SQUARE,&p2,&p3);
    a23a = sproduct(ANGLE,&p2,&p3);
    s14s = sproduct(SQUARE,&p1,&p4);
    a14a = sproduct(ANGLE,&p1,&p4);
    s132a = sproduct(SQUARE,&p1,&p3,&p2);
    s231a = sproduct(SQUARE,&p2,&p3,&p1);
  }
  void gAdd::set_masses(const ldouble& massd){
    md=massd;
    p3.set_mass(md);
    p4.set_mass(md);
    prop.set_mass(md);
  }
  void gAdd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
    s132a.update();
    s231a.update();
    //Propagator Momentum
    ldouble propTP[4], propUP[4];
    for(int j=0;j<4;j++){
      propTP[j] = mom1[j]-mom3[j];
      propUP[j] = mom1[j]-mom4[j];
    }
    pDenT = prop.denominator(propTP);
    pDenU = prop.denominator(propUP);

  }



  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble gAdd::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds1>0&&ds2>0){
      //gAdd:   md[12]^2<34>
      //34 out: md[12]^2<34>
      return 2.0*e*Qd*gs*md*s12s.v()*s12s.v()*a34a.v(ds3,ds4)/pDenT/pDenU;
    }
    else if(ds1<0&&ds2<0){
      //gAdd:   md<12>^2[34]
      //34 out: md<12>^2[34]
      return 2.0*e*Qd*gs*md*a12a.v()*a12a.v()*s34s.v(ds3,ds4)/pDenT/pDenU;
    }
    else if(ds1>0&&ds2<0){
      //gAdd:   ([31]<42>+[41]<32>)[132>
      //34 out: -([13]<24>+[14]<23>)[132>
      return -2.0*e*Qd*gs*(s13s.v(ds3)*a24a.v(ds4)+s14s.v(ds4)*a23a.v(ds3))*s132a.v()/pDenT/pDenU;
    }
    else if(ds1<0&&ds2>0){
      //gAdd:   (<31>[42]+<41>[32])[231>
      //34 out: -(<13>[24]+<14>[23])[231>
      return -2.0*e*Qd*gs*(a13a.v(ds3)*s24s.v(ds4)+a14a.v(ds4)*s23s.v(ds3))*s231a.v()/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) mdst be called before amp2().
  ldouble gAdd::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4 //The color factor is the same for both diagrams
	  }
    //Average over initial spins 1/2*1/2=1/4
    //Average over initial colors 1/8
    return amp2/32.0;
  }

  //A+, A+ -> e, E
  ldouble gAdd::amp2_gplus_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(2,2,j3,j4);
	amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
      }
    //Average over initial colors 1/8
    return amp2/8.0;
  }

  //A+, A- -> u, U
  ldouble gAdd::amp2_gplus_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(2,-2,j3,j4);
	amp2 += 4.0*std::pow(std::abs(M),2);//Color factor Tr(Ta,Ta) = 4
      }
    //Average over initial colors 1/8
    return amp2/8.0;
  }



  //  Tests
  int test_gAdd(){
    int n=0;//Number of fails
    std::cout<<"\t* g , A  -> d , D       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n#  md=0.0075, pspatial=250\n";
      ldouble md=0.0075;
      ldouble EE=0.31333, gs=1.238;
      gAdd gAddAmp = gAdd(EE,gs,md);
      ldouble pspatial=250;
      ldouble dataCH[20] = {6.524583287480624E-01,2.075529810155417E-01,1.194193495943857E-01,8.236316359351274E-02,6.244048293260498E-02,5.041817563206104E-02,4.277322124405065E-02,3.789574034367803E-02,3.497673642857076E-02,3.360502406236970E-02,3.360502406236970E-02,3.497673642857075E-02,3.789574034367802E-02,4.277322124405064E-02,5.041817563206104E-02,6.244048293260498E-02,8.236316359351270E-02,1.194193495943856E-01,2.075529810155417E-01,6.524583287480609E-01};
      i += gAddAmp.test_2to2_amp2([&]() { return gAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      i += gAddAmp.test_2to2_amp2_rotations([&]() { return gAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      i += gAddAmp.test_2to2_amp2_boosts([&]() { return gAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      i += gAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      //std::cout<<"\n#  g+ , A+ -> e , E\n";
      ldouble dataCHpp[20] = {1.266269058434672E-08,1.563180960053788E-09,6.288964131219419E-10,3.609368779587473E-10,2.474268405992016E-10,1.892665469090980E-10,1.563295173167902E-10,1.369596681769972E-10,1.259800261377601E-10,1.209788473504567E-10,1.209788468422370E-10,1.259800252280467E-10,1.369596659069489E-10,1.563295177081195E-10,1.892665422097592E-10,2.474268381970162E-10,3.609368770735979E-10,6.288964140693482E-10,1.563180960800024E-09,1.266269059290546E-08};
      i += gAddAmp.test_2to2_amp2([&]() { return gAddAmp.amp2_gplus_Aplus(); }, 0,0,md,md,pspatial,dataCHpp);
      i += gAddAmp.test_2to2_amp2_rotations([&]() { return gAddAmp.amp2_gplus_Aplus(); }, 0,0,md,md,pspatial,dataCHpp);
      i += gAddAmp.test_2to2_amp2_boosts([&]() { return gAddAmp.amp2_gplus_Aplus(); }, 0,0,md,md,pspatial,dataCHpp);
      i += gAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAddAmp.amp2_gplus_Aplus(); }, 0,0,md,md,pspatial,dataCHpp);
      //std::cout<<"\n#  g+ , A- -> e , E\n";
      ldouble dataCHpm[20] = {1.304916644833434E+00,4.151059604679025E-01,2.388386985598749E-01,1.647263268260886E-01,1.248809656177831E-01,1.008363510748555E-01,8.554644233177178E-02,7.579148055039640E-02,6.995347273116150E-02,6.721004800376057E-02,6.721004800376057E-02,6.995347273116148E-02,7.579148055039638E-02,8.554644233177176E-02,1.008363510748555E-01,1.248809656177831E-01,1.647263268260885E-01,2.388386985598748E-01,4.151059604679024E-01,1.304916644833431E+00};
      i += gAddAmp.test_2to2_amp2([&]() { return gAddAmp.amp2_gplus_Aminus(); }, 0,0,md,md,pspatial,dataCHpm);
      i += gAddAmp.test_2to2_amp2_rotations([&]() { return gAddAmp.amp2_gplus_Aminus(); }, 0,0,md,md,pspatial,dataCHpm);
      i += gAddAmp.test_2to2_amp2_boosts([&]() { return gAddAmp.amp2_gplus_Aminus(); }, 0,0,md,md,pspatial,dataCHpm);
      i += gAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAddAmp.amp2_gplus_Aminus(); }, 0,0,md,md,pspatial,dataCHpm);
      //Close to threshold
      //std::cout<<"\n#  md=0.0075, pspatial=0.008\n";
      pspatial = 0.008;
      ldouble dataCH2[20] = {4.251740181333088E-02,4.222129054559189E-02,4.191132593764320E-02,4.160893285203752E-02,4.132940872801938E-02,4.108367360944004E-02,4.087947522546553E-02,4.072221991104089E-02,4.061553991455593E-02,4.056166857297065E-02,4.056166857297065E-02,4.061553991455593E-02,4.072221991104089E-02,4.087947522546553E-02,4.108367360944004E-02,4.132940872801938E-02,4.160893285203752E-02,4.191132593764320E-02,4.222129054559189E-02,4.251740181333088E-02};
      i += gAddAmp.test_2to2_amp2([&]() { return gAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      i += gAddAmp.test_2to2_amp2_rotations([&]() { return gAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      i += gAddAmp.test_2to2_amp2_boosts([&]() { return gAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      i += gAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      //std::cout<<"\n#  g+ , A+ -> e , E\n";
      ldouble dataCH2pp[20] = {8.305613887776056E-02,7.913564947032858E-02,7.587920794559862E-02,7.319192577479136E-02,7.100055634711140E-02,6.924872719627594E-02,6.789352677079241E-02,6.690306416743368E-02,6.625474828816369E-02,6.593411851897575E-02,6.593411851897575E-02,6.625474828816369E-02,6.690306416743366E-02,6.789352677079241E-02,6.924872719627592E-02,7.100055634711140E-02,7.319192577479136E-02,7.587920794559862E-02,7.913564947032858E-02,8.305613887776055E-02};
      i += gAddAmp.test_2to2_amp2([&]() { return gAddAmp.amp2_gplus_Aplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      i += gAddAmp.test_2to2_amp2_rotations([&]() { return gAddAmp.amp2_gplus_Aplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      i += gAddAmp.test_2to2_amp2_boosts([&]() { return gAddAmp.amp2_gplus_Aplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      i += gAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAddAmp.amp2_gplus_Aplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      //std::cout<<"\n#  g+ , A- -> e , E\n";
      ldouble dataCH2pm[20] = {1.978664748901198E-03,5.306931620855186E-03,7.943443929687780E-03,1.002593992928367E-02,1.165826110892737E-02,1.291862002260415E-02,1.386542368013864E-02,1.454137565464809E-02,1.497633154094817E-02,1.518921862696554E-02,1.518921862696554E-02,1.497633154094817E-02,1.454137565464810E-02,1.386542368013864E-02,1.291862002260415E-02,1.165826110892737E-02,1.002593992928368E-02,7.943443929687780E-03,5.306931620855190E-03,1.978664748901204E-03};
      i += gAddAmp.test_2to2_amp2([&]() { return gAddAmp.amp2_gplus_Aminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      i += gAddAmp.test_2to2_amp2_rotations([&]() { return gAddAmp.amp2_gplus_Aminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      i += gAddAmp.test_2to2_amp2_boosts([&]() { return gAddAmp.amp2_gplus_Aminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      i += gAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return gAddAmp.amp2_gplus_Aminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
