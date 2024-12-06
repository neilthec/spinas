
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

//File:  SPINAS/SM/AAdd.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/AAdd.h"

namespace spinas {

  AAdd::AAdd(const ldouble& echarge, const ldouble& massd):
    e(echarge), Qd(-1.0/3.0), md(massd), prop(massd,0){
    constexpr ldouble two = 2;
    sqrt2 = std::sqrt(two);
    p1=particle(0);
    p2=particle(0);
    p3=particle(md);
    p4=particle(md);
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
  void AAdd::set_masses(const ldouble& massd){
    md=massd;
    p3.set_mass(md);
    p4.set_mass(md);
    prop.set_mass(md);
  }
  void AAdd::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
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
  cdouble AAdd::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    if(ds1>0&&ds2>0){
      //eeAA:   md[34]^2<12>
      //AAdd:   md[12]^2<34>
      //34 out: md[12]^2<34>
      return 2.0*e*e*Qd*Qd*md*s12s.v()*s12s.v()*a34a.v(ds3,ds4)/pDenT/pDenU;
    }
    else if(ds1<0&&ds2<0){
      //eeAA:   md<34>^2[12]
      //AAdd:   md<12>^2[34]
      //34 out: md<12>^2[34]
      return 2.0*e*e*Qd*Qd*md*a12a.v()*a12a.v()*s34s.v(ds3,ds4)/pDenT/pDenU;
    }
    else if(ds1>0&&ds2<0){
      //eeAA:   ([13]<24>+[23]<14>)[314>
      //AAdd:   ([31]<42>+[41]<32>)[132>
      //34 out: -([13]<24>+[14]<23>)[132>
      return -2.0*e*e*Qd*Qd*(s13s.v(ds3)*a24a.v(ds4)+s14s.v(ds4)*a23a.v(ds3))*s132a.v()/pDenT/pDenU;
    }
    else if(ds1<0&&ds2>0){
      //eeAA:   (<13>[24]+<23>[14])[413>
      //AAdd:   (<31>[42]+<41>[32])[231>
      //34 out: -(<13>[24]+<14>[23])[231>
      return -2.0*e*e*Qd*Qd*(a13a.v(ds3)*s24s.v(ds4)+a14a.v(ds4)*s23s.v(ds3))*s231a.v()/pDenT/pDenU;
    }
    return cdouble(0,0);    
  }

 
  //set_momenta(...) must be called before amp2().
  ldouble AAdd::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-2;j1<=2;j1+=4)
      for(int j2=-2;j2<=2;j2+=4)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    M = amp(j1,j2,j3,j4);
	    amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
	  }
    //Average over initial spins 1/2*1/2=1/4
    return amp2/4.0;
  }

  //A+, A+ -> e, E
  ldouble AAdd::amp2_Aplus_Aplus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(2,2,j3,j4);
	amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
      }
    return amp2;
  }

  //A+, A- -> u, U
  ldouble AAdd::amp2_Aplus_Aminus(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j3=-1;j3<=1;j3+=2)
      for(int j4=-1;j4<=1;j4+=2){
	M = amp(2,-2,j3,j4);
	amp2 += 3.0*std::pow(std::abs(M),2);//Color factor 3
      }
    return amp2;
  }



  //  Tests
  int test_AAdd(){
    int n=0;//Number of fails
    std::cout<<"\t* A , A  -> d , D       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n#  md=0.0075, pspatial=250\n";
      ldouble md=0.0075;
      ldouble EE=0.31333;
      AAdd AAddAmp = AAdd(EE,md);
      ldouble pspatial=250;
      ldouble dataCH[20] = {2.786276347861765E-02,8.863400717735465E-03,5.099717400961453E-03,3.517259640106254E-03,2.666475897058585E-03,2.153071913964272E-03,1.826599637448656E-03,1.618310324996591E-03,1.493656363055678E-03,1.435078373418408E-03,1.435078373418408E-03,1.493656363055678E-03,1.618310324996591E-03,1.826599637448656E-03,2.153071913964272E-03,2.666475897058585E-03,3.517259640106252E-03,5.099717400961451E-03,8.863400717735464E-03,2.786276347861759E-02};
      i += AAddAmp.test_2to2_amp2([&]() { return AAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      i += AAddAmp.test_2to2_amp2_rotations([&]() { return AAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      i += AAddAmp.test_2to2_amp2_boosts([&]() { return AAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      i += AAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH);
      //std::cout<<"\n#  A+ , A+ -> e , E\n";
      ldouble dataCHpp[20] = {5.407510904666671E-10,6.675451817407725E-11,2.685656882493397E-11,1.541354967536712E-11,1.056618520168913E-11,8.082491706426418E-12,6.675939556282049E-12,5.848764085536649E-12,5.379886372258350E-12,5.166314630599044E-12,5.166314608735007E-12,5.379886333400714E-12,5.848763988498437E-12,6.675939572905071E-12,8.082491505838428E-12,1.056618509888051E-11,1.541354963756828E-11,2.685656886527392E-11,6.675451820598657E-11,5.407510908321618E-10};
      i += AAddAmp.test_2to2_amp2([&]() { return AAddAmp.amp2_Aplus_Aplus(); }, 0,0,md,md,pspatial,dataCHpp);
      i += AAddAmp.test_2to2_amp2_rotations([&]() { return AAddAmp.amp2_Aplus_Aplus(); }, 0,0,md,md,pspatial,dataCHpp);
      i += AAddAmp.test_2to2_amp2_boosts([&]() { return AAddAmp.amp2_Aplus_Aplus(); }, 0,0,md,md,pspatial,dataCHpp);
      i += AAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAddAmp.amp2_Aplus_Aplus(); }, 0,0,md,md,pspatial,dataCHpp);
      //std::cout<<"\n#  A+ , A- -> e , E\n";
      ldouble dataCHpm[20] = {5.572552641648421E-02,1.772680136871642E-02,1.019943477506634E-02,7.034519264798958E-03,5.332951783550984E-03,4.306143819846052E-03,3.653199268221373E-03,3.236620644144419E-03,2.987312720731470E-03,2.870156741670500E-03,2.870156741670500E-03,2.987312720731469E-03,3.236620644144418E-03,3.653199268221373E-03,4.306143819846052E-03,5.332951783550984E-03,7.034519264798955E-03,1.019943477506633E-02,1.772680136871641E-02,5.572552641648408E-02};
      i += AAddAmp.test_2to2_amp2([&]() { return AAddAmp.amp2_Aplus_Aminus(); }, 0,0,md,md,pspatial,dataCHpm);
      i += AAddAmp.test_2to2_amp2_rotations([&]() { return AAddAmp.amp2_Aplus_Aminus(); }, 0,0,md,md,pspatial,dataCHpm);
      i += AAddAmp.test_2to2_amp2_boosts([&]() { return AAddAmp.amp2_Aplus_Aminus(); }, 0,0,md,md,pspatial,dataCHpm);
      i += AAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAddAmp.amp2_Aplus_Aminus(); }, 0,0,md,md,pspatial,dataCHpm);
      //Close to threshold
      //std::cout<<"\n#  md=0.0075, pspatial=0.008\n";
      pspatial = 0.008;
      ldouble dataCH2[20] = {1.815675052724516E-03,1.803029834090865E-03,1.789793018531183E-03,1.776879539385369E-03,1.764942662789696E-03,1.754448721359744E-03,1.745728576294978E-03,1.739013101239312E-03,1.734457408746791E-03,1.732156871865383E-03,1.732156871865383E-03,1.734457408746791E-03,1.739013101239312E-03,1.745728576294978E-03,1.754448721359744E-03,1.764942662789696E-03,1.776879539385369E-03,1.789793018531183E-03,1.803029834090865E-03,1.815675052724515E-03};
      i += AAddAmp.test_2to2_amp2([&]() { return AAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      i += AAddAmp.test_2to2_amp2_rotations([&]() { return AAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      i += AAddAmp.test_2to2_amp2_boosts([&]() { return AAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      i += AAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAddAmp.amp2(); }, 0,0,md,md,pspatial,dataCH2);
      //std::cout<<"\n#  A+ , A+ -> e , E\n";
      ldouble dataCH2pp[20] = {3.546852651017117E-03,3.379430971705727E-03,3.240366979435754E-03,3.125608527858963E-03,3.032027673163035E-03,2.957217069735940E-03,2.899344210646499E-03,2.857047217818231E-03,2.829361354664038E-03,2.815669091067328E-03,2.815669091067328E-03,2.829361354664039E-03,2.857047217818230E-03,2.899344210646498E-03,2.957217069735940E-03,3.032027673163035E-03,3.125608527858962E-03,3.240366979435754E-03,3.379430971705727E-03,3.546852651017116E-03};
      i += AAddAmp.test_2to2_amp2([&]() { return AAddAmp.amp2_Aplus_Aplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      i += AAddAmp.test_2to2_amp2_rotations([&]() { return AAddAmp.amp2_Aplus_Aplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      i += AAddAmp.test_2to2_amp2_boosts([&]() { return AAddAmp.amp2_Aplus_Aplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      i += AAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAddAmp.amp2_Aplus_Aplus(); }, 0,0,md,md,pspatial,dataCH2pp);
      //std::cout<<"\n#  A+ , A- -> e , E\n";
      ldouble dataCH2pm[20] = {8.449745443191446E-05,2.266286964760033E-04,3.392190576266117E-04,4.281505509117760E-04,4.978576524163568E-04,5.516803729835484E-04,5.921129419434580E-04,6.209789846603932E-04,6.395534628295437E-04,6.486446526634381E-04,6.486446526634381E-04,6.395534628295437E-04,6.209789846603932E-04,5.921129419434581E-04,5.516803729835484E-04,4.978576524163569E-04,4.281505509117761E-04,3.392190576266117E-04,2.266286964760035E-04,8.449745443191471E-05};
      i += AAddAmp.test_2to2_amp2([&]() { return AAddAmp.amp2_Aplus_Aminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      i += AAddAmp.test_2to2_amp2_rotations([&]() { return AAddAmp.amp2_Aplus_Aminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      i += AAddAmp.test_2to2_amp2_boosts([&]() { return AAddAmp.amp2_Aplus_Aminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      i += AAddAmp.test_2to2_amp2_boosts_and_rotations([&]() { return AAddAmp.amp2_Aplus_Aminus(); }, 0,0,md,md,pspatial,dataCH2pm);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

  
  

}
