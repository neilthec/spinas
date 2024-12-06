
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

//File:  SPINAS/SM/mnen.cpp

#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>

#include "spinas.h"
#include "include/mnen.h"

namespace spinas {
  //Constructors
  mnen::mnen(const ldouble& echarge, const ldouble& masse, const ldouble& massmu, const ldouble& massW, const ldouble& widthW, const ldouble& sinW):
    e(echarge), me(masse), mm(massmu), MW(massW), WW(widthW), SW(sinW), prop(massW,widthW),
    p1(particle(mm)), p2(particle(0)),
    p3(particle(me)), p4(particle(0)),
    //[24], ⟨13⟩, [14], ⟨23⟩
    s34s(sproduct(SQUARE,&p3,&p4,2)),
    a12a(sproduct(ANGLE,&p1,&p2,2)),
    s14s(sproduct(SQUARE,&p1,&p4,2)),
    a23a(sproduct(ANGLE,&p2,&p3,2))
  {
    preW = e*e/(4.0*MW*MW*SW*SW);
  }
  void mnen::set_masses(const ldouble& masse, const ldouble& massmu, const ldouble& massW){
    me=masse;
    mm=massmu;
    MW=massW;
    p1.set_mass(mm);
    p3.set_mass(me);
    prop.set_mass(massW);
    preW = e*e/(4.0*MW*MW*SW*SW);
  }
  void mnen::set_momenta(const ldouble mom1[4], const ldouble mom2[4], const ldouble mom3[4], const ldouble mom4[4]){
    //Particles
    p1.set_momentum(mom1);
    p2.set_momentum(mom2);
    p3.set_momentum(mom3);
    p4.set_momentum(mom4);
    s34s.update();
    a12a.update();
    s14s.update();
    a23a.update();
    //Propagator Momentum
    ldouble propP[4];
    for(int j=0;j<4;j++)
      propP[j] = mom1[j]-mom4[j];
    pDenU = prop.denominator(propP);
  }

  
  //Amplitude
  //set_momenta(...) must be called before amp(...).
  cdouble mnen::amp(const int& ds1, const int& ds3){
    //preW = e*e/(4.0*MW*MW*SW*SW);
    //mEneNm all ingoing:
    //- preW ( 2MW^2 [24]⟨13⟩ + mμme[14]⟨23⟩ )/(u−MW^2)
    //mneENm: 2<->3
    //- preW ( 2MW^2 [34]⟨12⟩ - mμme[14]⟨23⟩ )/(u−MW^2)
    //34 out:
    //- preW ( 2MW^2 [34]⟨12⟩ + mμme[14]⟨23⟩ )/(u−MW^2)
    return - preW*2.0*(
		       + 2.0*MW*MW*s34s.v(ds3)*a12a.v(ds1)
		       + mm*me*s14s.v(ds1)*a23a.v(ds3)
		       )/pDenU;
  }
  
  //set_momenta(...) must be called before amp2().
  ldouble mnen::amp2(){
    ldouble amp2 = 0;
    cdouble M;

    //Sum over spins
    for(int j1=-1;j1<=1;j1+=2)
      for(int j2=-1;j2<=1;j2+=2){
	M = amp(j1,j2);
	amp2 += std::pow(std::abs(M),2);
      }
    //Average over initial spins 1/2
    return amp2/2.0;
  }

  



  //  Tests
  int test_mnen(){
    int n=0;//Number of fails
    std::cout<<"\t* m , ne -> e , nm      :";
    {//amp^2
      int i=0;
      // me=0.0005, mmu=0.106, pspatial=250
      ldouble me=0.0005, mmu=0.106, MW=80.385, WW=0.0;//WW=2.10;
      ldouble sinW=0.474;
      mnen mnenAmp = mnen(0.31333,me,mmu,MW,WW,sinW);
      ldouble pspatial=250;
      ldouble dataCH[20] = {9.530796478311548E-02,1.055949863344363E-01,1.176420291236528E-01,1.318747461218597E-01,1.488557018847690E-01,1.693408290933995E-01,1.943650367858616E-01,2.253755337542719E-01,2.644458705578949E-01,3.146317949579756E-01,3.805874233571693E-01,4.696844905185881E-01,5.941653407188850E-01,7.755843867442443E-01,1.054807539970082E+00,1.517213777044642E+00,2.366649312143242E+00,4.195572832580421E+00,9.387256319365109E+00,3.692615125342182E+01};
      i += mnenAmp.test_2to2_amp2([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH);
      i += mnenAmp.test_2to2_amp2_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH);
      i += mnenAmp.test_2to2_amp2_boosts([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH);
      i += mnenAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH);
      // me=0.0005, mmu=0.106, pspatial=0.05
      pspatial = 0.05;
      ldouble dataCH2[20] = {1.068741198879758E-12,1.068741475417470E-12,1.068741751955289E-12,1.068742028493215E-12,1.068742305031249E-12,1.068742581569390E-12,1.068742858107639E-12,1.068743134645995E-12,1.068743411184458E-12,1.068743687723029E-12,1.068743964261707E-12,1.068744240800492E-12,1.068744517339385E-12,1.068744793878384E-12,1.068745070417492E-12,1.068745346956706E-12,1.068745623496028E-12,1.068745900035457E-12,1.068746176574994E-12,1.068746453114638E-12};
      i += mnenAmp.test_2to2_amp2([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH2);
      i += mnenAmp.test_2to2_amp2_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH2);
      i += mnenAmp.test_2to2_amp2_boosts([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH2);
      i += mnenAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH2);
      // me=0.106, mmu=0.0005, pspatial=0.005
      me=0.106;
      mmu=0.0005;
      pspatial = 250;
      mnenAmp.set_masses(me,mmu,MW);
      ldouble dataCH3[20] = {9.530796434063593E-02,1.055949858184189E-01,1.176420285168561E-01,1.318747454016783E-01,1.488557010210992E-01,1.693408280454459E-01,1.943650354972361E-01,2.253755321452584E-01,2.644458685128411E-01,3.146317923039599E-01,3.805874198263032E-01,4.696844856778796E-01,5.941653338314011E-01,7.755843764725229E-01,1.054807523678641E+00,1.517213748940532E+00,2.366649257391169E+00,4.195572703343499E+00,9.387255886843567E+00,3.692614787898892E+01};
      i += mnenAmp.test_2to2_amp2([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH3);
      i += mnenAmp.test_2to2_amp2_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH3);
      i += mnenAmp.test_2to2_amp2_boosts([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH3);
      i += mnenAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH3);
      // me=0.106, mmu=0.0005, pspatial=0.05
      me=0.106;
      mmu=0.0005;
      pspatial = 0.11;
      mnenAmp.set_masses(me,mmu,MW);
      ldouble dataCH4[20] = {4.112780355804128E-12,4.112782721227851E-12,4.112785086653616E-12,4.112787452081422E-12,4.112789817511267E-12,4.112792182943154E-12,4.112794548377081E-12,4.112796913813049E-12,4.112799279251058E-12,4.112801644691107E-12,4.112804010133197E-12,4.112806375577329E-12,4.112808741023500E-12,4.112811106471712E-12,4.112813471921965E-12,4.112815837374259E-12,4.112818202828593E-12,4.112820568284968E-12,4.112822933743383E-12,4.112825299203840E-12};
      i += mnenAmp.test_2to2_amp2([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH4);
      i += mnenAmp.test_2to2_amp2_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH4);
      i += mnenAmp.test_2to2_amp2_boosts([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH4);
      i += mnenAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH4);
      // me=0.105, mmu=0.106, pspatial=250
      me=0.105;
      mmu=0.106;
      pspatial=250;
      mnenAmp.set_masses(me,mmu,MW);
      ldouble dataCH5[20] = {9.530796876908343E-02,1.055949907379726E-01,1.176420340139137E-01,1.318747515841406E-01,1.488557080254952E-01,1.693408360470866E-01,1.943650447250136E-01,2.253755429037683E-01,2.644458812164259E-01,3.146318075308895E-01,3.805874384084647E-01,4.696845088563475E-01,5.941653635419489E-01,7.755844159073542E-01,1.054807578493301E+00,1.517213830164071E+00,2.366649389648119E+00,4.195572954195883E+00,9.387256521133901E+00,3.692615122629651E+01};
      i += mnenAmp.test_2to2_amp2([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH5);
      i += mnenAmp.test_2to2_amp2_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH5);
      i += mnenAmp.test_2to2_amp2_boosts([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH5);
      i += mnenAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH5);
      // me=0.105, mmu=0.106, pspatial=0.05
      pspatial=0.05;
      ldouble dataCH6[20] = {6.472710093476381E-13,6.472710885880502E-13,6.472711678284801E-13,6.472712470689276E-13,6.472713263093931E-13,6.472714055498763E-13,6.472714847903770E-13,6.472715640308958E-13,6.472716432714321E-13,6.472717225119863E-13,6.472718017525582E-13,6.472718809931479E-13,6.472719602337553E-13,6.472720394743805E-13,6.472721187150234E-13,6.472721979556842E-13,6.472722771963627E-13,6.472723564370588E-13,6.472724356777727E-13,6.472725149185045E-13};
      i += mnenAmp.test_2to2_amp2([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH6);
      i += mnenAmp.test_2to2_amp2_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH6);
      i += mnenAmp.test_2to2_amp2_boosts([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH6);
      i += mnenAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH6);
      // me=105, mmu=106, MW=0.1, pspatial=500
      me=105;
      mmu=106;
      MW=0.1;
      pspatial=500;
      mnenAmp.set_masses(me,mmu,MW);
      ldouble dataCH7[20] = {3.024500072291648E+10,3.028191918183252E+10,3.032308399299144E+10,3.036927244781155E+10,3.042146386552105E+10,3.048090978759284E+10,3.054923562264361E+10,3.062859150709253E+10,3.072188262110615E+10,3.083313247451267E+10,3.096807824940524E+10,3.113519171130177E+10,3.134752852426555E+10,3.162631261546312E+10,3.200850597022922E+10,3.256470656915260E+10,3.344863607544114E+10,3.507023307932272E+10,3.900690448558545E+10,6.192865231312740E+10};
      i += mnenAmp.test_2to2_amp2([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH7);
      i += mnenAmp.test_2to2_amp2_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH7);
      i += mnenAmp.test_2to2_amp2_boosts([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH7);
      i += mnenAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH7);
      // me=105, mmu=106, MW=0.1, pspatial=5
      pspatial=5;
      ldouble dataCH8[20] = {4.395259672966754E+08,4.346730067052537E+08,4.298525910599265E+08,4.250646138993171E+08,4.203089691317419E+08,4.155855510337607E+08,4.108942542487346E+08,4.062349737853901E+08,4.016076050163892E+08,3.970120436769064E+08,3.924481858632112E+08,3.879159280312577E+08,3.834151669952803E+08,3.789457999263946E+08,3.745077243512059E+08,3.701008381504233E+08,3.657250395574793E+08,3.613802271571561E+08,3.570662998842183E+08,3.527831570220510E+08};
      i += mnenAmp.test_2to2_amp2([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH8);
      i += mnenAmp.test_2to2_amp2_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH8);
      i += mnenAmp.test_2to2_amp2_boosts([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH8);
      i += mnenAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH8);
      // me=0.105, mmu=0.106, MW=0.11, pspatial=250
      me=0.105;
      mmu=0.106;
      MW=0.11;
      pspatial=250;
      mnenAmp.set_masses(me,mmu,MW);
      ldouble dataCH9[20] = {1.206221995504843E-01,1.317727273147437E-01,1.448889097521171E-01,1.604614488123332E-01,1.791443013562088E-01,2.018244731974548E-01,2.297292768016694E-01,2.645959296713681E-01,3.089486231844368E-01,3.665685704336158E-01,4.433274504455260E-01,5.487447461294915E-01,6.990880237509389E-01,9.240469248385320E-01,1.282600047224204E+00,1.906010248427453E+00,3.137563272896713E+00,6.130236330752228E+00,1.699252486347331E+01,1.527707765023219E+02};
      i += mnenAmp.test_2to2_amp2([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH9);
      i += mnenAmp.test_2to2_amp2_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH9);
      i += mnenAmp.test_2to2_amp2_boosts([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH9);
      i += mnenAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH9);
      // me=0.105, mmu=0.106, MW=0.11, pspatial=0.005
      pspatial=0.005;
      ldouble dataCH10[20] = {3.501121442715958E-02,3.391477572095502E-02,3.280608096105405E-02,3.168497748712708E-02,3.055131038580769E-02,2.940492245166537E-02,2.824565414739786E-02,2.707334356322519E-02,2.588782637546731E-02,2.468893580428666E-02,2.347650257057645E-02,2.225035485197514E-02,2.101031823798678E-02,1.975621568418682E-02,1.848786746549192E-02,1.720509112847215E-02,1.590770144268326E-02,1.459551035099606E-02,1.326832691889939E-02,1.192595728275270E-02};
      i += mnenAmp.test_2to2_amp2([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH10);
      i += mnenAmp.test_2to2_amp2_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH10);
      i += mnenAmp.test_2to2_amp2_boosts([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH10);
      i += mnenAmp.test_2to2_amp2_boosts_and_rotations([&]() { return mnenAmp.amp2(); }, mmu,0,me,0,pspatial,dataCH10);
      // Done
      if(i==0) std::cout<<"                                         Pass"<<std::endl;
      else std::cout<<"                                         Fail!"<<std::endl;
      n+=i;
    }
    
    return n;
  }

 
  
  

}
