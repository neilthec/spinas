
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

//File:  SPINAS/SM/AAee.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AAee.h"

namespace spinas {

  AAee::AAee(const ldouble& echarge, const ldouble& masse):
    e(echarge), me(masse), prop(masse,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(0);
    p2=particle(0);
    p3=particle(me);
    p4=particle(me);
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
    s132a = sproduct(SQUARE,&p1,&p3,&p2,2);
    s231a = sproduct(SQUARE,&p2,&p3,&p1,2);
  }
  void AAee::set_masses(const ldouble& masse){
    me=masse;
    p3.set_mass(me);
    p4.set_mass(me);
    prop.set_mass(me);
  }
  void AAee::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
  cdouble AAee::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds1>0&&ds2>0){
      //eeAA:   me[34]^2<12>
      //AAee:   me[12]^2<34>
      //34 out: me[12]^2<34>
      return 2.0*e*e*me*s12s.v()*s12s.v()*a34a.v(ds3,ds4)/pDenT/pDenU;
    }
    else if(ds1<0&&ds2<0){
      //eeAA:   me<34>^2[12]
      //AAee:   me<12>^2[34]
      //34 out: me<12>^2[34]
      return 2.0*e*e*me*a12a.v()*a12a.v()*s34s.v(ds3,ds4)/pDenT/pDenU;
    }
    else if(ds1>0&&ds2<0){
      //eeAA:   ([13]<24>+[23]<14>)[314>
      //AAee:   ([31]<42>+[41]<32>)[132>
      //34 out: -([13]<24>+[14]<23>)[132>
      return -2.0*e*e*(s13s.v(ds3)*a24a.v(ds4)+s14s.v(ds4)*a23a.v(ds3))*s132a.v()/pDenT/pDenU;
    }
    else if(ds1<0&&ds2>0){
      //eeAA:   (<13>[24]+<23>[14])[413>
      //AAee:   (<31>[42]+<41>[32])[231>
      //34 out: -(<13>[24]+<14>[23])[231>
      return -2.0*e*e*(a13a.v(ds3)*s24s.v(ds4)+a14a.v(ds4)*s23s.v(ds3))*s231a.v()/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble AAee::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2=1/4
    return amp2/4.0;
  }

  //A+, A+ -> e, E
  ldouble AAee::amp2_Aplus_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(2,2,j3,j4);
	amp2 += std::pow(std::abs(M),2);
      }
    return amp2;
  }

  //A+, A- -> e, E
  ldouble AAee::amp2_Aplus_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(2,-2,j3,j4);
	amp2 += std::pow(std::abs(M),2);
      }
    return amp2;
  }



  //  Tests
  int test_AAee(){
    int n=0;//Number of fails
    std::cout<<"\t* A , A  -> e , E       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n#  me=0.0005, pspatial=250\n";
      ldouble me=0.0005;
      ldouble EE=0.31333;
      AAee AAeeAmp = AAee(EE,me);
      ldouble pspatial=250;
      ldouble dataCH[20] = {7.522946197731721E-01,2.393118197781031E-01,1.376923698710783E-01,9.496601025075935E-02,7.199484916448802E-02,5.813294161240127E-02,4.931819014337163E-02,4.369437870612651E-02,4.032872173345139E-02,3.874711601320890E-02,3.874711601320890E-02,4.032872173345137E-02,4.369437870612649E-02,4.931819014337162E-02,5.813294161240125E-02,7.199484916448799E-02,9.496601025075925E-02,1.376923698710782E-01,2.393118197781027E-01,7.522946197731676E-01};
      i += AAeeAmp.test_2to2_amp2([&]() { return AAeeAmp.amp2(); }, 0,0,me,me,pspatial,dataCH);
      i += AAeeAmp.test_2to2_amp2_rotations([&]() { return AAeeAmp.amp2(); }, 0,0,me,me,pspatial,dataCH);
      i += AAeeAmp.test_2to2_amp2_boosts([&]() { return AAeeAmp.amp2(); }, 0,0,me,me,pspatial,dataCH);
      i += AAeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAeeAmp.amp2(); }, 0,0,me,me,pspatial,dataCH);
      //std::cout<<"\n#  A+ , A+ -> e , E\n";
      ldouble dataCHpp[20] = {6.489013194551735E-11,8.010541909278626E-12,3.222786637845131E-12,1.849625790731140E-12,1.267944374440302E-12,9.698999212819949E-13,8.011098447617939E-13,7.018475613917085E-13,6.455920003369539E-13,6.199559688177666E-13,6.199553911412966E-13,6.455907484222579E-13,7.018546358108840E-13,8.011162872943907E-13,9.698994266147537E-13,1.267939541270305E-12,1.849623589292510E-12,3.222786596340517E-12,8.010543715999903E-12,6.489013576214194E-11};
      i += AAeeAmp.test_2to2_amp2([&]() { return AAeeAmp.amp2_Aplus_Aplus(); }, 0,0,me,me,pspatial,dataCHpp);
      i += AAeeAmp.test_2to2_amp2_rotations([&]() { return AAeeAmp.amp2_Aplus_Aplus(); }, 0,0,me,me,pspatial,dataCHpp);
      i += AAeeAmp.test_2to2_amp2_boosts([&]() { return AAeeAmp.amp2_Aplus_Aplus(); }, 0,0,me,me,pspatial,dataCHpp);
      i += AAeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAeeAmp.amp2_Aplus_Aplus(); }, 0,0,me,me,pspatial,dataCHpp);
      //std::cout<<"\n#  A+ , A- -> e , E\n";
      ldouble dataCHpm[20] = {1.504589239481454E+00,4.786236395481956E-01,2.753847397389338E-01,1.899320204996691E-01,1.439896983277081E-01,1.162658832238326E-01,9.863638028594214E-02,8.738875741155117E-02,8.065744346625718E-02,7.749423202579785E-02,7.749423202579785E-02,8.065744346625715E-02,8.738875741155112E-02,9.863638028594213E-02,1.162658832238326E-01,1.439896983277081E-01,1.899320204996689E-01,2.753847397389337E-01,4.786236395481949E-01,1.504589239481445E+00};
      i += AAeeAmp.test_2to2_amp2([&]() { return AAeeAmp.amp2_Aplus_Aminus(); }, 0,0,me,me,pspatial,dataCHpm);
      i += AAeeAmp.test_2to2_amp2_rotations([&]() { return AAeeAmp.amp2_Aplus_Aminus(); }, 0,0,me,me,pspatial,dataCHpm);
      i += AAeeAmp.test_2to2_amp2_boosts([&]() { return AAeeAmp.amp2_Aplus_Aminus(); }, 0,0,me,me,pspatial,dataCHpm);
      i += AAeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAeeAmp.amp2_Aplus_Aminus(); }, 0,0,me,me,pspatial,dataCHpm);
      //Close to threshold
      //std::cout<<"\n#  me=0.0005, pspatial=0.0006\n";
      pspatial = 0.0006;
      ldouble dataCH2[20] = {7.095517491415575E-02,6.787650295313721E-02,6.499522353873709E-02,6.244154484907705E-02,6.026194050136540E-02,5.846595419023571E-02,5.704788745810282E-02,5.599715439187751E-02,5.530330527724128E-02,5.495846266515966E-02,5.495846266515966E-02,5.530330527724127E-02,5.599715439187750E-02,5.704788745810282E-02,5.846595419023570E-02,6.026194050136541E-02,6.244154484907704E-02,6.499522353873709E-02,6.787650295313721E-02,7.095517491415573E-02};
      i += AAeeAmp.test_2to2_amp2([&]() { return AAeeAmp.amp2(); }, 0,0,me,me,pspatial,dataCH2);
      i += AAeeAmp.test_2to2_amp2_rotations([&]() { return AAeeAmp.amp2(); }, 0,0,me,me,pspatial,dataCH2);
      i += AAeeAmp.test_2to2_amp2_boosts([&]() { return AAeeAmp.amp2(); }, 0,0,me,me,pspatial,dataCH2);
      i += AAeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAeeAmp.amp2(); }, 0,0,me,me,pspatial,dataCH2);
      //std::cout<<"\n#  A+ , A+ -> e , E\n";
      ldouble dataCH2pp[20] = {1.332816670543560E-01,1.151310926439474E-01,1.019386729607098E-01,9.217038792958972E-02,8.487328956706151E-02,7.943453319223830E-02,7.545128212165271E-02,7.265723022970423E-02,7.087987414918964E-02,7.001555341187911E-02,7.001555341187911E-02,7.087987414918964E-02,7.265723022970422E-02,7.545128212165271E-02,7.943453319223830E-02,8.487328956706151E-02,9.217038792958969E-02,1.019386729607097E-01,1.151310926439474E-01,1.332816670543560E-01};
      i += AAeeAmp.test_2to2_amp2([&]() { return AAeeAmp.amp2_Aplus_Aplus(); }, 0,0,me,me,pspatial,dataCH2pp);
      i += AAeeAmp.test_2to2_amp2_rotations([&]() { return AAeeAmp.amp2_Aplus_Aplus(); }, 0,0,me,me,pspatial,dataCH2pp);
      i += AAeeAmp.test_2to2_amp2_boosts([&]() { return AAeeAmp.amp2_Aplus_Aplus(); }, 0,0,me,me,pspatial,dataCH2pp);
      i += AAeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAeeAmp.amp2_Aplus_Aplus(); }, 0,0,me,me,pspatial,dataCH2pp);
      //std::cout<<"\n#  A+ , A- -> e , E\n";
      ldouble dataCH2pm[20] = {8.628682773955452E-03,2.062191326232702E-02,2.805177411676443E-02,3.271270176856437E-02,3.565059143566930E-02,3.749737518823312E-02,3.864449279455293E-02,3.933707855405078E-02,3.972673640529290E-02,3.990137191844022E-02,3.990137191844022E-02,3.972673640529290E-02,3.933707855405078E-02,3.864449279455293E-02,3.749737518823312E-02,3.565059143566930E-02,3.271270176856438E-02,2.805177411676444E-02,2.062191326232704E-02,8.628682773955450E-03};
      i += AAeeAmp.test_2to2_amp2([&]() { return AAeeAmp.amp2_Aplus_Aminus(); }, 0,0,me,me,pspatial,dataCH2pm);
      i += AAeeAmp.test_2to2_amp2_rotations([&]() { return AAeeAmp.amp2_Aplus_Aminus(); }, 0,0,me,me,pspatial,dataCH2pm);
      i += AAeeAmp.test_2to2_amp2_boosts([&]() { return AAeeAmp.amp2_Aplus_Aminus(); }, 0,0,me,me,pspatial,dataCH2pm);
      i += AAeeAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAeeAmp.amp2_Aplus_Aminus(); }, 0,0,me,me,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
