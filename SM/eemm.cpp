
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

//File:  SPINAS/SM/eemm.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/eemm.h"

namespace spinas {
  //Constructors
  eemm::eemm(const ldouble& echarge, const ldouble& masse, const ldouble& massmu, const ldouble& massh, const ldouble& widthh, const ldouble& massW, const ldouble& sinW, const ldouble& widthZ):
    e(echarge), me(masse), mm(massmu), mh(massh), wh(widthh), MW(massW), SW(sinW), CW(std::sqrt(1.0-sinW*sinW)), MZ(massW/CW), WZ(widthZ),
    propA(0,0), proph(mh,wh), propZ(MZ,WZ),
    p1(particle(me)), p2(particle(me)),
    p3(particle(mm)), p4(particle(mm)),
    a13a(sproduct(ANGLE,&p1,&p3)),
    s13s(sproduct(SQUARE,&p1,&p3)),
    a14a(sproduct(ANGLE,&p1,&p4)),
    s14s(sproduct(SQUARE,&p1,&p4)),
    a23a(sproduct(ANGLE,&p2,&p3)),
    s23s(sproduct(SQUARE,&p2,&p3)),
    a24a(sproduct(ANGLE,&p2,&p4)),
    s24s(sproduct(SQUARE,&p2,&p4)),
    s12s(sproduct(SQUARE,&p1,&p2)),
    a12a(sproduct(ANGLE,&p1,&p2)),
    s34s(sproduct(SQUARE,&p3,&p4)),
    a34a(sproduct(ANGLE,&p3,&p4))
  {
    //For some reason, MZ doesn't get set correctly above.  Redo it here.
    MZ=MW/CW;
    propZ.set_mass(MZ);
    preh = e*e*me*mm/(4.0*MW*MW*SW*SW);
    gL=2.0*SW*SW-1.0;
    gR=2.0*SW*SW;
    preZ = e*e/(4.0*CW*CW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*me*mm/MZ/MZ;//=preh!
    preZLL = preZ*2.0*gL*gL;
    preZLR = preZ*2.0*gL*gR;
    preZRR = preZ*2.0*gR*gR;
  }
  void eemm::set_masses(const ldouble& masse, const ldouble& massmu, const ldouble& massh, const ldouble& massW){
    me=masse;
    mm=massmu;
    mh=massh;
    proph.set_mass(mh);
    MW=massW;
    MZ=MW/CW;//std::cout<<"MZ="<<MZ<<std::endl;
    propZ.set_mass(MZ);
    p1.set_mass(me);
    p2.set_mass(me);
    p3.set_mass(mm);
    p4.set_mass(mm);
    preh = e*e*me*mm/(4*MW*MW*SW*SW);
    preZ0 = preZ*(gL-gR)*(gL-gR)*me*mm/MZ/MZ;
  }
  void eemm::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    a13a.update();
    s13s.update();
    a14a.update();
    s14s.update();
    a23a.update();
    s23s.update();
    a24a.update();
    s24s.update();
    s12s.update();
    a12a.update();
    s34s.update();
    a34a.update();
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]+mom2[j];
    pDenSA = propA.denominator(propP);
    pDenSh = proph.denominator(propP);
    pDenSZ = propZ.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble eemm::amp(const int& ds1, const int& ds2, const int& ds3, const int& ds4){
    constexpr ldouble one=1, two = 2, four=4, eight=8;
    cdouble amplitude(0,0);
    
    //Photon
    //Sign changes due to p3 and p4 being outgoing.
    // (<13>[24] + <14>[23] + [13]<24> + [14]<23>)
    amplitude += two*e*e*(
			  a13a.v(ds1,ds3)*s24s.v(ds2,ds4) + a14a.v(ds1,ds4)*s23s.v(ds2,ds3)
			  + s13s.v(ds1,ds3)*a24a.v(ds2,ds4) + s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
			  )/pDenSA;
    
    //Higgs
    //EE^2 Me Mm / (4 MW^2 SW^2) * ([12]+<12>) ([34]+<34>)/(s-Mh^2)
    //preh = e*e*me*mm/(4*MW*MW*SW*SW);
    amplitude += preh*(s12s.v(ds1,ds2)+a12a.v(ds1,ds2))*(s34s.v(ds3,ds4)+a34a.v(ds3,ds4))/pDenSh;
    
    //Z Boson
    //Defined above:
    //gL=2.0*SW*SW-1.0;
    //gR=2.0*SW*SW;
    //preZ = e*e/(4.0*CW*CW*SW*SW);
    //preZ0 = preZ*(gL-gR)*(gL-gR)*me*mm/MZ/MZ; // = preh
    //preZLL = preZ*2.0*gL*gL;
    //preZLR = preZ*2.0*gL*gR;
    //preZRR = preZ*2.0*gR*gR;
    //all in:
    //+(EE^2 Me Mm (gL-gR)^2 (<12>-[12]) (<34>-[34]))/(8 CW^2 MZ^2 SW^2 (s-MZ^2))
    //+(EE^2 (gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>)/(4 CW^2 SW^2 (s-MZ^2))
    //= - preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //  - preZ 2(gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    //34 out:
    //- preZ0 (<12>-[12]) (<34>-[34]) / (s-MZ^2)
    //+ preZ 2(gL^2 [23] <14> + gLgR( [13] <24>+ [24] <13> ) + gR^2 [14] <23>))/(s-MZ^2)
    amplitude += 
      - preZ0*(a12a.v(ds1,ds2)-s12s.v(ds1,ds2))*(a34a.v(ds3,ds4)-s34s.v(ds3,ds4))/pDenSZ
      + two*preZ*(
	      gL*gL*s23s.v(ds2,ds3)*a14a.v(ds1,ds4)
	      + gL*gR*(s13s.v(ds1,ds3)*a24a.v(ds2,ds4)+s24s.v(ds2,ds4)*a13a.v(ds1,ds3))
	      + gR*gR*s14s.v(ds1,ds4)*a23a.v(ds2,ds3)
	      )/pDenSZ;
    
    return amplitude;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble eemm::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2)
	for(int j3=-1;j3<=1;j3+=2)
	  for(int j4=-1;j4<=1;j4+=2){
	    //M = amp(j1,j2,j3,j4);
	    M = amp(j1,j2,j3,j4);
	    amp2 += std::pow(std::abs(M),2);
	  }
    //Average over initial spins 1/2*1/2 = 1/4
    return amp2/4.0;
  }

  



  //  Tests
  int test_eemm(){
    int n=0;//Number of fails
    std::cout<<"\t* e , E  -> m , M       :";
    {//amp^2
      int i=0;
      //std::cout<<"\n###########  All Diagrams\n";
      //std::cout<<"########### me=0.0005, mmu=0.105, pspatial=250\n";
      ldouble me=0.0005, mmu=0.105, mh=125, wh=0, MW=80.385, SW=0.474, WZ=0;//Set width to 0 for comparison with Feynman diagrams.
      eemm eemmAmp = eemm(0.31333,me,mmu,mh,wh,MW,SW,WZ);
      ldouble pspatial=250;
      ldouble dataCH[20] = {3.475296719289669E-02,3.131901456177322E-02,2.810640362061734E-02,2.511513436942906E-02,2.234520680820838E-02,1.979662093695529E-02,1.746937675566979E-02,1.536347426435189E-02,1.347891346300158E-02,1.181569435161887E-02,1.037381693020375E-02,9.153281198756229E-03,8.154087157276302E-03,7.376234805763969E-03,6.819724144219231E-03,6.484555172642088E-03,6.370727891032542E-03,6.478242299390588E-03,6.807098397716230E-03,7.357296186009467E-03};
      i += eemmAmp.test_2to2_amp2([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
      i += eemmAmp.test_2to2_amp2_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
      i += eemmAmp.test_2to2_amp2_boosts([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
      i += eemmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH);
      //std::cout<<"########### me=0.0005, mmu=0.105, pspatial=0.11\n";
      pspatial = 0.11;
      ldouble dataCH2[20] = {1.919358346433044E-02,1.903942059310635E-02,1.890238719587382E-02,1.878248327263286E-02,1.867970882338346E-02,1.859406384812563E-02,1.852554834685936E-02,1.847416231958466E-02,1.843990576630153E-02,1.842277868700996E-02,1.842278108170996E-02,1.843991295040152E-02,1.847417429308465E-02,1.852556510975934E-02,1.859408540042560E-02,1.867973516508343E-02,1.878251440373282E-02,1.890242311637377E-02,1.903946130300629E-02,1.919362896363038E-02};
      i += eemmAmp.test_2to2_amp2([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH2);
      i += eemmAmp.test_2to2_amp2_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH2);
      i += eemmAmp.test_2to2_amp2_boosts([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH2);
      i += eemmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH2);
      //std::cout<<"########### me=0.105, mmu=0.0005, pspatial=0.005\n";
      me=0.105;
      mmu=0.0005;
      pspatial = 0.005;
      eemmAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH3[20] = {1.927501920178369E-02,1.927109447257587E-02,1.926760586316216E-02,1.926455337354255E-02,1.926193700371706E-02,1.925975675368568E-02,1.925801262344841E-02,1.925670461300524E-02,1.925583272235619E-02,1.925539695150124E-02,1.925539730044040E-02,1.925583376917368E-02,1.925670635770106E-02,1.925801506602255E-02,1.925975989413815E-02,1.926194084204786E-02,1.926455790975169E-02,1.926761109724962E-02,1.927110040454166E-02,1.927502583162780E-02};
      i += eemmAmp.test_2to2_amp2([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH3);
      i += eemmAmp.test_2to2_amp2_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH3);
      i += eemmAmp.test_2to2_amp2_boosts([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH3);
      i += eemmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH3);
      //std::cout<<"########### me=0.10, mmu=0.105, pspatial=0.05\n";
      me=0.10;
      mmu=0.105;
      pspatial = 0.05;
      eemmAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH4[20] = {2.605564194232666E-02,2.601469901591307E-02,2.597830544519909E-02,2.594646123018471E-02,2.591916637086995E-02,2.589642086725479E-02,2.587822471933924E-02,2.586457792712330E-02,2.585548049060696E-02,2.585093240979024E-02,2.585093368467312E-02,2.585548431525561E-02,2.586458430153771E-02,2.587823364351942E-02,2.589643234120074E-02,2.591918039458166E-02,2.594647780366219E-02,2.597832456844234E-02,2.601472068892209E-02,2.605566616510144E-02};
      i += eemmAmp.test_2to2_amp2([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH4);
      i += eemmAmp.test_2to2_amp2_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH4);
      i += eemmAmp.test_2to2_amp2_boosts([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH4);
      i += eemmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH4);
      //std::cout<<"########### me=0.10, mmu=0.105, MW=0.11, pspatial=0.05\n";
      me=0.10;
      mmu=0.105;
      MW=0.11;
      pspatial = 0.05;
      eemmAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH5[20] = {3.842257265146085E-02,3.805758797327145E-02,3.769846457436183E-02,3.734520245473198E-02,3.699780161438190E-02,3.665626205331159E-02,3.632058377152105E-02,3.599076676901027E-02,3.566681104577928E-02,3.534871660182805E-02,3.503648343715659E-02,3.473011155176490E-02,3.442960094565298E-02,3.413495161882084E-02,3.384616357126846E-02,3.356323680299585E-02,3.328617131400301E-02,3.301496710428995E-02,3.274962417385666E-02,3.249014252270313E-02};
      i += eemmAmp.test_2to2_amp2([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH5);
      i += eemmAmp.test_2to2_amp2_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH5);
      i += eemmAmp.test_2to2_amp2_boosts([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH5);
      i += eemmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH5);
      //std::cout<<"########### me=0.10, mmu=0.105, MW=0.006, pspatial=0.05\n";
      me=0.10;
      mmu=0.105;
      MW=0.006;
      pspatial = 0.05;
      eemmAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH6[20] = {1.015220287293079E+03,1.015220026649366E+03,1.015219771187368E+03,1.015219520907085E+03,1.015219275808516E+03,1.015219035891662E+03,1.015218801156523E+03,1.015218571603099E+03,1.015218347231390E+03,1.015218128041395E+03,1.015217914033115E+03,1.015217705206550E+03,1.015217501561700E+03,1.015217303098565E+03,1.015217109817144E+03,1.015216921717439E+03,1.015216738799448E+03,1.015216561063172E+03,1.015216388508611E+03,1.015216221135764E+03};
      i += eemmAmp.test_2to2_amp2([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH6);
      i += eemmAmp.test_2to2_amp2_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH6);
      i += eemmAmp.test_2to2_amp2_boosts([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH6);
      i += eemmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH6);
      //std::cout<<"########### me=0.10, mmu=0.105, MW=0.11, Mh=0.125, pspatial=0.05\n";
      me=0.10;
      mmu=0.105;
      MW=0.11;
      mh=0.125;
      pspatial = 0.05;
      eemmAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH7[20] = {4.220806754664559E-02,4.149184053397863E-02,4.078147480059144E-02,4.007697034648402E-02,3.937832717165637E-02,3.868554527610849E-02,3.799862465984038E-02,3.731756532285204E-02,3.664236726514347E-02,3.597303048671467E-02,3.530955498756565E-02,3.465194076769638E-02,3.400018782710690E-02,3.335429616579719E-02,3.271426578376724E-02,3.208009668101706E-02,3.145178885754666E-02,3.082934231335603E-02,3.021275704844516E-02,2.960203306281407E-02};
      i += eemmAmp.test_2to2_amp2([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH7);
      i += eemmAmp.test_2to2_amp2_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH7);
      i += eemmAmp.test_2to2_amp2_boosts([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH7);
      i += eemmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH7);
      //std::cout<<"########### me=0.10, mmu=0.105, MW=0.006, Mh=0.125, pspatial=0.05\n";
      me=0.10;
      mmu=0.105;
      MW=0.006;
      pspatial = 0.05;
      eemmAmp.set_masses(me,mmu,mh,MW);
      ldouble dataCH8[20] = {1.067029084815383E+03,1.066910963483066E+03,1.066792847332465E+03,1.066674736363578E+03,1.066556630576406E+03,1.066438529970949E+03,1.066320434547206E+03,1.066202344305179E+03,1.066084259244866E+03,1.065966179366268E+03,1.065848104669385E+03,1.065730035154217E+03,1.065611970820763E+03,1.065493911669025E+03,1.065375857699001E+03,1.065257808910692E+03,1.065139765304097E+03,1.065021726879218E+03,1.064903693636053E+03,1.064785665574603E+03};
      i += eemmAmp.test_2to2_amp2([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH8);
      i += eemmAmp.test_2to2_amp2_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH8);
      i += eemmAmp.test_2to2_amp2_boosts([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH8);
      i += eemmAmp.test_2to2_amp2_boosts_and_rotations([&]() { return eemmAmp.amp2(); }, me,me,mmu,mmu,pspatial,dataCH8);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
